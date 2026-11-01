"""A file that fails in the BED/modkit merge is reported, and its batch goes on (ERR-1 layer 1).

Unit tests drive ``merge_modkit_files``, the parallel ``_process_bams`` and the real
``bed_conversion_handler``; the Coordinator tests run the real Coordinator on both schedulers
(``failure_fixtures``) with a synthetic handler that leaves ``ctx.metadata["file_failures"]``.
"""

from types import SimpleNamespace
from unittest import mock

import polars as pl
import pytest

from failure_fixtures import SCHEDULERS, Harness, init_ray, shutdown_ray
from robin import workflow_ray as wr
from robin.analysis import bed_conversion
from robin.analysis.temp_utilities import merge_modkit_files


def _good(path, pos=100):
    pl.DataFrame(
        {
            "chrom": ["chr1"], "chromStart": [pos], "mod_code": ["m"], "strand": ["+"],
            "valid_cov": [10], "percent_modified": [50.0], "n_mod": [5], "n_canonical": [5],
        },
        schema={
            "chrom": pl.Utf8, "chromStart": pl.UInt32, "mod_code": pl.Utf8, "strand": pl.Utf8,
            "valid_cov": pl.UInt32, "percent_modified": pl.Float32, "n_mod": pl.UInt32,
            "n_canonical": pl.UInt32,
        },
    ).write_parquet(path)


def _merge(tmp_path, files):
    panel = tmp_path / "panel.bed"
    panel.write_text("chr1\t100\t101\tcg0001\n")
    out = tmp_path / "out.parquet"
    ret = merge_modkit_files(
        new_files=[str(f) for f in files], existing_file=str(out), output_file=str(out),
        filter_bed_file=str(panel), sample_id="s", output_dir=str(tmp_path),
        mnpflex_config=None, num_bam_files_seen=len(files),
    )
    return ret, out


def test_merge_returns_the_broken_file_and_keeps_the_good_one(tmp_path):
    good, broken = tmp_path / "good.parquet", tmp_path / "broken.parquet"
    _good(good)
    pl.DataFrame({"chrom": ["chr1"], "chromStart": [100]}).write_parquet(broken)
    failures, out = _merge(tmp_path, [good, broken])
    assert [f["file"] for f in failures] == [str(broken)]
    assert failures[0]["error"]
    assert pl.read_parquet(out).select("chrom", "chromStart").rows() == [("chr1", 100)]


def test_merge_with_no_broken_file_returns_empty_list(tmp_path):
    good = tmp_path / "good.parquet"
    _good(good)
    failures, _ = _merge(tmp_path, [good])
    assert failures == []


def test_merge_returns_failures_when_nothing_is_left_to_merge(tmp_path):
    broken = tmp_path / "broken.parquet"
    pl.DataFrame({"chrom": ["chr1"], "chromStart": [100]}).write_parquet(broken)
    failures, out = _merge(tmp_path, [broken])
    assert [f["file"] for f in failures] == [str(broken)]
    assert not out.exists()


def _analysis(tmp_path, threads):
    return bed_conversion.BedConversionAnalysis(work_dir=str(tmp_path), threads=threads)


def test_parallel_process_bams_returns_the_bam_that_raised(tmp_path):
    bams = []
    for n in ("a", "bad", "c"):
        b = tmp_path / f"{n}.bam"
        b.write_text("x")
        bams.append(str(b))

    def fake_matkit(bam, out):
        if "bad" in bam:
            raise RuntimeError("matkit blew up")
        _good(out)

    with mock.patch.object(bed_conversion, "run_matkit", fake_matkit):
        processed, failures = _analysis(tmp_path, 3)._process_bams(bams, str(tmp_path))
    assert len(processed) == 2
    assert [f["file"] for f in failures] == [bams[1]]
    assert "matkit blew up" in failures[0]["error"]


def test_sequential_process_bams_still_raises(tmp_path):
    b = tmp_path / "bad.bam"
    b.write_text("x")
    boom = mock.Mock(side_effect=RuntimeError("boom"))
    with mock.patch.object(bed_conversion, "run_matkit", boom):
        with pytest.raises(RuntimeError):
            _analysis(tmp_path, 1)._process_bams([str(b)], str(tmp_path))


def _batched_job(tmp_path, names, fail_only=False):
    paths = []
    for n in names:
        p = tmp_path / n
        p.write_text("x")
        paths.append(str(p))
    inner = [
        wr.WorkflowContext(filepath=p, metadata={"bam_metadata": {"sample_id": "s"}}) for p in paths
    ]
    batched = SimpleNamespace(
        get_file_count=lambda: len(paths), get_sample_id=lambda: "s", batch_id="b1",
        get_filepaths=lambda: paths, contexts=inner,
    )
    meta = {"_batched_job": batched}
    if fail_only:
        meta["fail_only_bam_submission"] = True
    ctx = wr.WorkflowContext(filepath=paths[0], metadata=meta)
    job = wr.Job(7, "bed_conversion", "bed_conversion", ["bed_conversion:bed_conversion"], 0, ctx)
    return job, paths


def _run_handler_with_failing_merge(tmp_path, names, fail_only=False):
    job, paths = _batched_job(tmp_path, names, fail_only)

    def fake_merge(data, state, *a, **k):
        _good(state)
        return [{"file": data[1], "error": "synthetic merge failure"}]

    with mock.patch.object(bed_conversion, "run_matkit", lambda bam, out: _good(out)), \
            mock.patch.object(bed_conversion, "merge_modkit_files", fake_merge):
        bed_conversion.bed_conversion_handler(job, work_dir=str(tmp_path / "work"))
    return job, paths


def test_real_handler_puts_a_merge_failure_in_file_failures_and_files_no_error(tmp_path):
    job, paths = _run_handler_with_failing_merge(tmp_path, ["a_pass_0.bam", "b_pass_0.bam"])
    ff = job.context.metadata["file_failures"]
    assert len(ff) == 1
    assert ff[0]["file"] == paths[1] and ff[0]["step"] == "merge"
    assert "synthetic merge failure" in ff[0]["error"]
    assert job.context.errors == []
    assert job.context.results["bed_conversion"]["status"] == "success"


def test_real_handler_with_no_failure_sets_no_file_failures(tmp_path):
    job, _ = _batched_job(tmp_path, ["a_pass_0.bam", "b_pass_0.bam"])
    with mock.patch.object(bed_conversion, "run_matkit", lambda bam, out: _good(out)), \
            mock.patch.object(
                bed_conversion, "merge_modkit_files", lambda data, state, *a, **k: _good(state) or []
            ):
        bed_conversion.bed_conversion_handler(job, work_dir=str(tmp_path / "work"))
    assert not job.context.metadata.get("file_failures")
    assert job.context.errors == []
    assert job.context.results["bed_conversion"]["status"] == "success"


def test_real_handler_failed_bam_in_the_parallel_step_is_a_file_failure(tmp_path):
    job, paths = _batched_job(tmp_path, ["a_pass_0.bam", "bad_pass_0.bam"])

    def fake_matkit(bam, out):
        if "bad" in bam:
            raise RuntimeError("matkit blew up")
        _good(out)

    def fake_merge(data, state, *a, **k):
        _good(state)
        return []

    real_process = bed_conversion.process_multiple_files

    def threaded(*a, **k):
        k["threads"] = 2
        return real_process(*a, **k)

    with mock.patch.object(bed_conversion, "run_matkit", fake_matkit), \
            mock.patch.object(bed_conversion, "merge_modkit_files", fake_merge), \
            mock.patch.object(bed_conversion, "process_multiple_files", threaded):
        bed_conversion.bed_conversion_handler(job, work_dir=str(tmp_path / "work"))
    ff = job.context.metadata["file_failures"]
    assert [(f["file"], f["step"]) for f in ff] == [(paths[1], "bam")]
    assert job.context.errors == []


def test_fail_only_submission_file_failure_records_nothing(tmp_path):
    job, _ = _run_handler_with_failing_merge(
        tmp_path, ["a_fail_0.bam", "b_fail_0.bam"], fail_only=True
    )
    assert not job.context.metadata.get("file_failures")
    assert job.context.errors == []


# --- the real Coordinator, both schedulers --------------------------------------------------------


@pytest.fixture(scope="module")
def ray_cluster():
    init_ray()
    yield
    shutdown_ray()


@pytest.fixture(scope="module", params=SCHEDULERS)
def fcampaign(request, ray_cluster, tmp_path_factory):
    harness = Harness(request.param, str(tmp_path_factory.mktemp(f"filefail_{request.param}")))
    c = SimpleNamespace(scheduler=request.param)
    try:
        c.before = harness.stats()
        c.ff = harness.run_one("bed_conversion", "file_failure", tag="ff", downstream=True)
        c.ok = harness.run_one("bed_conversion", "ok", tag="ok")
        c.sturgeon_after = harness.finished_count("sturgeon", timeout=30)
        c.stats = harness.stats()
    finally:
        c.rows = harness.close()
    return c


def _publish(stats):
    sent = []
    wr.publish_failures(stats, {"last_seq": 0}, lambda t, d, p=0: sent.append((t.value, d)))
    return [d for t, d in sent if t == "error_update"]


def test_F4_job_with_a_file_failure_completes_and_triggers_its_next_jobs(fcampaign):
    o = fcampaign.ff
    assert o.d_failed == 0 and o.d_failed_by_type == {}
    assert o.d_completed_by_type.get("bed_conversion") == 1
    assert fcampaign.sturgeon_after >= 1, "the classifiers did not run after the file failure"
    assert fcampaign.stats["failed"] == fcampaign.before["failed"]


def test_F4_one_file_failed_row_per_failed_file(fcampaign):
    rows = [r for r in fcampaign.rows if r["sample_id"] == fcampaign.ff.sample_id]
    assert len(rows) == 1
    r = rows[0]
    assert r["failure_kind"] == "file_failed" and r["job_type"] == "bed_conversion"
    assert r["filepath"].endswith("_bad.bam") and "synthetic merge failure" in r["error"]
    assert r["is_oom"] == "no"
    assert [x for x in fcampaign.rows if x["sample_id"] == fcampaign.ok.sample_id] == []


def test_F4_error_update_lists_the_failed_file(fcampaign):
    errs = _publish(fcampaign.stats)
    assert len(errs) == 1
    mine = [r for r in errs[0]["recent"] if r["sample_id"] == fcampaign.ff.sample_id]
    assert len(mine) == 1
    assert mine[0]["job_type"] == "bed_conversion"
    assert "_bad.bam" in mine[0]["error"] and "synthetic merge failure" in mine[0]["error"]
    assert errs[0]["preprocessing_errors"] == 1
    assert errs[0]["analysis_errors"] == 0 and errs[0]["classification_errors"] == 0


def test_F4_file_failure_is_one_event_and_job_counts_unchanged(fcampaign):
    events = fcampaign.stats["recent_events"]
    ev = [e for e in events if e["kind"] == "file_failed"]
    assert len(ev) == 1 and ev[0]["sample_id"] == fcampaign.ff.sample_id
    assert not [e for e in events if e["kind"] == "failure"]
    assert fcampaign.stats["failed_by_type"].get("bed_conversion", 0) == 0
