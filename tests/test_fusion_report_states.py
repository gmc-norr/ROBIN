"""REP-6: the fusion report tells "not run", "unreadable" and "none found" apart, and its
methods text matches the filters (batch 7, Task 7).

Rows from the locked plan:
  R1  No fusion output at all -> "Fusion analysis has not been run for this sample."
  R2  A processed fusion pickle that can't be read -> says the fusion results could not be read
  R3  Fusion's preprocessing raises -> the fusion handler files an error; judged failed
  C5  Fusion ran and no read passed the filters (new run: an empty processed pickle; old folder:
      ``sv_count.txt``, no pickle, no candidate rows) -> "No fusion candidates were identified."
  Also: the methods sentence states MAPQ >= 50, mapping span >= 250 bp and >= 3 supporting reads,
  built from ``FUSION_DETECTION_THRESHOLDS``; the filter is checked at the boundaries.

The report section is the real ``FusionSection`` on a temp sample folder; the writer is the real
``accumulate_fusion_candidates`` on synthetic staged rows; the R3 handler path is the real
``fusion_handler`` with only the preprocessing replaced by a function that raises. "Judged failed"
is the shared outcome helper applied to the job's context.

Which tests can fail on the old code by assertion: every report-text test (the old section says
"No fusion candidates were identified." for all three cases and ">50"/">200bp"/"4 or more" in the
methods), the empty-pickle tests (no pickle was written), and R3 (the preprocessing failure was
only a warning). ``test_methods_sentence_follows_the_config`` and the sentence helper tests can only
fail by ImportError there.
"""

import os
import pickle
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pysam
import pytest

from robin import classification_config as cc
from robin import workflow_ray as wr
from robin.analysis import fusion_analysis, fusion_work
from robin.job_outcome import job_outcome

SAMPLE = "S1"

NOT_RUN = "Fusion analysis has not been run for this sample."
UNREADABLE = "The fusion results could not be read."
NONE_FOUND = "No fusion candidates were identified."
NONE_FOUND_SUMMARY = "No fusion networks identified within the target panel."


# ------------------------------------------------------------------ synthetic data
def _rows(n_reads, gene_a="GENEA", gene_b="GENEB", prefix="r", mapq=60, span=400):
    """n_reads reads, each aligned to gene_a (primary) and gene_b (supplementary)."""
    rows = []
    for i in range(n_reads):
        common = dict(
            col2=1000, col3=5000, read_id=f"{prefix}{i}", mapping_quality=mapq, strand="+",
            is_secondary=False, mapping_span=span,
        )
        rows.append(dict(common, col1="chr1", col4=gene_a, reference_id="chr1",
                         reference_start=2000 + i, reference_end=2000 + i + span,
                         read_start=0, read_end=span, is_supplementary=False))
        rows.append(dict(common, col1="chr2", col4=gene_b, reference_id="chr2",
                         reference_start=3000 + i, reference_end=3000 + i + span,
                         read_start=span, read_end=2 * span, is_supplementary=True))
    return pd.DataFrame(rows)


def _sample_dir(tmp_path):
    d = tmp_path / SAMPLE
    d.mkdir()
    return d


def _report_text(sample_dir):
    from robin.reporting.sections.fusion import FusionSection
    from robin.reporting.styling.styles import ReportStyles

    report = SimpleNamespace(
        output=str(sample_dir), styles=ReportStyles(tempfile.mkdtemp())
    )
    section = FusionSection(report)
    section.add_content()
    body = "\n".join(
        el.getPlainText() for el in section.elements if hasattr(el, "getPlainText")
    )
    summary = "\n".join(
        el.getPlainText() for el in section.summary_elements if hasattr(el, "getPlainText")
    )
    return body, summary


def _write_pickles(sample_dir, df, names=("master", "all")):
    for n in names:
        fusion_work.preprocess_fusion_data_standalone(
            df, str(sample_dir / f"fusion_candidates_{n}_processed.pkl")
        )


def _stage(work_dir, df, kinds=("target", "genome"), n=1):
    staging = Path(fusion_work._get_staging_dir(str(work_dir), SAMPLE))
    for kind in kinds:
        df.to_parquet(staging / f"{kind}_{n}.parquet", index=False)


def _accumulate(work_dir):
    return fusion_work.accumulate_fusion_candidates(
        str(work_dir), SAMPLE, "rCNS2", force=True, batch_size=1
    )


def _load(path):
    with open(path, "rb") as f:
        return pickle.load(f)


# ------------------------------------------------------------------ R1, R2
def test_r1_no_output_says_not_run(tmp_path):
    body, summary = _report_text(_sample_dir(tmp_path))
    assert NOT_RUN in body
    assert NOT_RUN in summary
    assert NONE_FOUND not in body
    assert NONE_FOUND_SUMMARY not in summary


@pytest.mark.parametrize("content", [b"", b"not a pickle at all", b"\x80\x04\x95"])
def test_r2_unreadable_pickle_says_could_not_be_read(tmp_path, content):
    d = _sample_dir(tmp_path)
    (d / "fusion_candidates_master_processed.pkl").write_bytes(content)
    body, summary = _report_text(d)
    assert UNREADABLE in body
    assert UNREADABLE in summary
    assert NONE_FOUND not in body
    assert NOT_RUN not in body
    assert NONE_FOUND_SUMMARY not in summary


def test_r2_unreadable_beats_sv_count(tmp_path):
    d = _sample_dir(tmp_path)
    (d / "sv_count.txt").write_text("0")
    (d / "fusion_candidates_all_processed.pkl").write_bytes(b"garbage")
    body, _ = _report_text(d)
    assert UNREADABLE in body
    assert NONE_FOUND not in body


# ------------------------------------------------------------------ C5
def test_c5_new_run_where_nothing_passes_writes_empty_pickles(tmp_path):
    # Two reads for the gene pair: below the 3-read support.
    _stage(tmp_path, _rows(2))
    result = _accumulate(tmp_path)
    assert result["status"] == "success"
    d = tmp_path / SAMPLE
    for name in ("master", "all"):
        pkl = d / f"fusion_candidates_{name}_processed.pkl"
        assert pkl.exists(), f"no processed pickle for {name}"
        data = _load(pkl)
        assert data["annotated_data"].empty
        assert len(data["goodpairs"]) == 0
        assert data["gene_pairs"] == [] and data["gene_groups"] == []
        assert data["candidate_count"] == 0
    assert not list(d.glob("*.tmp")), "temp files left behind"
    body, summary = _report_text(d)
    assert NONE_FOUND in body
    assert NONE_FOUND_SUMMARY in summary
    assert NOT_RUN not in body and UNREADABLE not in body


def test_c5_empty_pickle_replaces_a_stale_populated_one(tmp_path):
    d = _sample_dir(tmp_path)
    _write_pickles(d, _rows(4, prefix="old"))
    assert not _load(d / "fusion_candidates_master_processed.pkl")["annotated_data"].empty
    _stage(tmp_path, _rows(2))
    _accumulate(tmp_path)
    for name in ("master", "all"):
        data = _load(d / f"fusion_candidates_{name}_processed.pkl")
        assert data["annotated_data"].empty
    body, _ = _report_text(d)
    assert NONE_FOUND in body


def test_c5_old_folder_with_sv_count_and_no_pickle_says_none_found(tmp_path):
    d = _sample_dir(tmp_path)
    (d / "sv_count.txt").write_text("0")
    body, summary = _report_text(d)
    assert NONE_FOUND in body
    assert NONE_FOUND_SUMMARY in summary
    assert NOT_RUN not in body


def test_c5_old_folder_with_header_only_csv_says_none_found(tmp_path):
    d = _sample_dir(tmp_path)
    (d / "sv_count.txt").write_text("0")
    _rows(1).iloc[0:0].to_csv(d / "fusion_candidates_master.csv", index=False)
    body, _ = _report_text(d)
    assert NONE_FOUND in body
    assert UNREADABLE not in body


def test_c5_old_folder_with_zero_byte_csv_says_none_found(tmp_path):
    d = _sample_dir(tmp_path)
    (d / "sv_count.txt").write_text("0")
    (d / "fusion_candidates_master.csv").write_text("")
    body, summary = _report_text(d)
    assert NONE_FOUND in body and NONE_FOUND_SUMMARY in summary
    assert UNREADABLE not in body


def test_old_folder_with_candidate_rows_but_no_pickle_is_unreadable(tmp_path):
    d = _sample_dir(tmp_path)
    (d / "sv_count.txt").write_text("1")
    _rows(4).to_csv(d / "fusion_candidates_master.csv", index=False)
    body, summary = _report_text(d)
    assert UNREADABLE in body
    assert UNREADABLE in summary
    assert NONE_FOUND not in body


# ------------------------------------------------------------------ control: candidates
def test_candidates_are_reported_with_the_methods_sentence(tmp_path):
    d = _sample_dir(tmp_path)
    _write_pickles(d, _rows(4))
    body, summary = _report_text(d)
    assert "GENEA" in summary and "GENEB" in summary
    assert NONE_FOUND not in body and NOT_RUN not in body and UNREADABLE not in body
    assert fusion_report_sentence() in body
    assert "(>50)" not in body and "(>200bp)" not in body


# ------------------------------------------------------------------ methods sentence
def fusion_report_sentence():
    from robin.reporting.sections.fusion import fusion_methods_sentence

    return fusion_methods_sentence()


def test_methods_sentence_states_the_locked_values():
    s = fusion_report_sentence()
    assert "mapping quality of at least 50" in s
    assert "mapping span of at least 250 bp" in s
    assert "at least 3 supporting reads" in s


def test_methods_sentence_follows_the_config(monkeypatch):
    t = cc.FUSION_DETECTION_THRESHOLDS
    monkeypatch.setitem(t, "mapping_quality", dict(t["mapping_quality"], min_threshold=59))
    monkeypatch.setitem(t, "mapping_span", dict(t["mapping_span"], min_threshold=399))
    monkeypatch.setitem(t, "read_support", dict(t["read_support"], min_threshold=5))
    s = fusion_report_sentence()
    assert "mapping quality of at least 60" in s
    assert "mapping span of at least 400 bp" in s
    assert "at least 5 supporting reads" in s


# ------------------------------------------------------------------ the filters at the boundaries
GENES = {
    "chr1": [fusion_work.GeneRegion(0, 100000, "GENEA")],
    "chr2": [fusion_work.GeneRegion(0, 100000, "GENEB")],
}


def _bam(tmp_path, mapq, span):
    path = str(tmp_path / "t.bam")
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 200000}, {"SN": "chr2", "LN": 200000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for flag, ref, pos, other_ref, other_pos in (
            (0, 0, 2000, "chr2", 3000),
            (2048, 1, 3000, "chr1", 2000),
        ):
            a = pysam.AlignedSegment(out.header)
            a.query_name = "read1"
            a.query_sequence = "A" * span
            a.flag = flag
            a.reference_id = ref
            a.reference_start = pos
            a.mapping_quality = mapq
            a.cigarstring = f"{span}M"
            a.set_tag("SA", f"{other_ref},{other_pos + 1},+,{span}M,{mapq},0;")
            out.write(a)
    return path


@pytest.mark.parametrize(
    "mapq,span,kept",
    [(49, 400, False), (50, 400, True), (60, 249, False), (60, 250, True)],
)
def test_filter_boundaries_match_the_stated_values(tmp_path, mapq, span, kept):
    df = fusion_work._process_reads_for_fusions(_bam(tmp_path, mapq, span), {"read1"}, GENES)
    assert (df is not None and not df.empty) is kept


def test_master_bed_filter_boundaries_match_the_stated_values(tmp_path):
    for mapq, span, kept in [(49, 400, False), (50, 400, True), (60, 249, False), (60, 250, True)]:
        df = fusion_work._process_reads_for_master_bed_fusions(
            _bam(tmp_path, mapq, span), {"read1"}
        )
        assert (df is not None and not df.empty) is kept, (mapq, span)


@pytest.mark.parametrize("n_reads,good", [(2, False), (3, True)])
def test_read_support_boundary(n_reads, good):
    _, goodpairs = fusion_work._annotate_results(_rows(n_reads))
    assert bool(goodpairs.any()) is good


# ------------------------------------------------------------------ R3 through the real handler
def _batched_job(tmp_path):
    bam = str(tmp_path / "in_0.bam")
    ctx = wr.WorkflowContext(
        filepath=bam,
        metadata={
            "bam_metadata": {"sample_id": SAMPLE, "has_supplementary_reads": True},
            "target_panel": "rCNS2",
        },
    )
    ctx.metadata["_batched_job"] = wr.BatchedJob(
        job_id=1, job_type="fusion", origin="fusion", workflow=["fusion:fusion"], step=0,
        contexts=[ctx], batch_id="b1", sample_id=SAMPLE,
    )
    return wr.Job(1, "fusion", "fusion", ["fusion:fusion"], 0, ctx)


def _run_handler(tmp_path, monkeypatch, preprocess=None):
    def staged(bam_path, temp_dir, metadata, fusion_metadata, target_panel, **kw):
        _stage(kw["work_dir"], _rows(4), kinds=("target", "genome"))
        return {}, False

    monkeypatch.setattr(fusion_analysis, "process_bam_with_staging", staged)
    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    if preprocess is not None:
        monkeypatch.setattr(fusion_work, "preprocess_fusion_data_standalone", preprocess)
    job = _batched_job(tmp_path)
    before = len(job.context.errors)
    fusion_analysis.fusion_handler(job, work_dir=str(tmp_path / "work"), target_panel="rCNS2")
    return job, job.context.errors[before:], job_outcome("fusion", job.context, before)


def test_r3_a_raising_preprocessing_is_filed_and_the_job_is_judged_failed(tmp_path, monkeypatch):
    def boom(df, output_file):
        raise RuntimeError("synthetic preprocessing crash")

    _, filed, outcome = _run_handler(tmp_path, monkeypatch, preprocess=boom)
    assert [e["job_type"] for e in filed] == ["fusion_analysis"]
    assert "synthetic preprocessing crash" in filed[0]["error"]
    assert outcome["status"] == "failed"


def test_r3_control_a_working_preprocessing_files_nothing(tmp_path, monkeypatch):
    job, filed, outcome = _run_handler(tmp_path, monkeypatch)
    assert filed == []
    assert outcome["status"] == "success"
    assert (tmp_path / "work" / SAMPLE / "fusion_candidates_master_processed.pkl").exists()


# ------------------------------------------------------------------ a failed build (Ruling J)
def _boom(df, output_file):
    raise RuntimeError("synthetic preprocessing crash")


def test_a_raising_preprocessing_still_finishes_the_steps_other_outputs(tmp_path, monkeypatch):
    _, filed, outcome = _run_handler(tmp_path, monkeypatch, preprocess=_boom)
    assert outcome["status"] == "failed" and filed
    d = tmp_path / "work" / SAMPLE
    assert (d / "sv_count.txt").exists()
    assert list((d / "bed_files").glob("fusion_breakpoints_*.bed")), "breakpoint BED not written"
    assert not list((d / "_fusion_staging").glob("*.parquet")), "staging not cleaned"
    assert not (d / "fusion_candidates_master_processed.pkl").exists()
    # The report does not call it "none found": results exist but were not processed.
    body, _ = _report_text(d)
    assert UNREADABLE in body


def test_a_failed_build_is_rebuilt_on_the_next_run_without_double_counting(tmp_path, monkeypatch):
    work = tmp_path / "work"
    d = work / SAMPLE
    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    _stage(work, _rows(4))
    real = fusion_work.preprocess_fusion_data_standalone
    monkeypatch.setattr(fusion_work, "preprocess_fusion_data_standalone", _boom)
    assert _accumulate(work)["status"] == "error"
    assert not (d / "fusion_candidates_master_processed.pkl").exists()
    csv_rows = len(pd.read_csv(d / "fusion_candidates_master.csv"))
    assert csv_rows == 8

    # Second run: preprocessing works again; the only new rows are one read that does not pass.
    monkeypatch.setattr(fusion_work, "preprocess_fusion_data_standalone", real)
    _stage(work, _rows(1, "GENEC", "GENED", prefix="x"), n=2)
    assert _accumulate(work)["status"] == "success"
    for name in ("master", "all"):
        data = _load(d / f"fusion_candidates_{name}_processed.pkl")
        assert not data["annotated_data"].empty
        assert data["annotated_data"]["read_id"].nunique() == 4
    assert len(pd.read_csv(d / "fusion_candidates_master.csv")) == csv_rows
    body, summary = _report_text(d)
    assert "GENEA" in summary


# ------------------------------------------------------------------ one candidate type at a time (I1)
# The target-panel (master) and genome-wide (all) results each have their own state: one of them
# missing or unreadable must not read as "none found" because the other holds candidates.
TARGET_UNREADABLE = "The target-panel fusion results could not be read."
GENOME_WIDE_UNREADABLE = "The genome-wide fusion results could not be read."
ONE_GENOME_WIDE = "Additional 1 genome-wide fusion networks identified."


def _target_unreadable_genome_wide_shown(d):
    body, summary = _report_text(d)
    assert NONE_FOUND_SUMMARY not in summary, summary
    assert TARGET_UNREADABLE in summary, summary
    # the genome-wide candidates still show
    assert ONE_GENOME_WIDE in summary, summary
    assert "Genome-wide Fusions" in body
    assert NONE_FOUND not in body and NOT_RUN not in body


def test_i1_target_pickle_missing_while_its_csv_has_rows(tmp_path):
    # P4 A: the genome-wide pickle holds candidates (writing it also writes sv_count.txt); the
    # target panel's candidates are in its CSV only.
    d = _sample_dir(tmp_path)
    _write_pickles(d, _rows(4), names=("all",))
    _rows(4).to_csv(d / "fusion_candidates_master.csv", index=False)
    assert (d / "sv_count.txt").exists()
    _target_unreadable_genome_wide_shown(d)


def test_i1_target_pickle_unreadable(tmp_path):
    # P4 B: a truncated target-panel pickle next to a genome-wide one with candidates.
    d = _sample_dir(tmp_path)
    _write_pickles(d, _rows(4), names=("all",))
    (d / "fusion_candidates_master_processed.pkl").write_bytes(b"\x80\x04\x95")
    _target_unreadable_genome_wide_shown(d)


def test_i1_only_the_target_panel_build_raised(tmp_path, monkeypatch):
    # P4 C, through the real writer: both builds run (Ruling J); the preprocessing raises only
    # for the target-panel pickle.
    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    real = fusion_work.preprocess_fusion_data_standalone

    def target_panel_only_raises(df, output_file):
        if "master" in os.path.basename(output_file):
            raise RuntimeError("synthetic target-panel crash")
        return real(df, output_file)

    monkeypatch.setattr(
        fusion_work, "preprocess_fusion_data_standalone", target_panel_only_raises
    )
    _stage(tmp_path, _rows(4))
    result = _accumulate(tmp_path)
    assert result["status"] == "error"
    d = tmp_path / SAMPLE
    assert not (d / "fusion_candidates_master_processed.pkl").exists()
    assert len(pd.read_csv(d / "fusion_candidates_master.csv")) == 8
    assert not _load(d / "fusion_candidates_all_processed.pkl")["annotated_data"].empty
    _target_unreadable_genome_wide_shown(d)


def test_i1_genome_wide_pickle_unreadable_is_said_beside_the_target_panel_candidates(tmp_path):
    d = _sample_dir(tmp_path)
    _write_pickles(d, _rows(4), names=("master",))
    (d / "fusion_candidates_all_processed.pkl").write_bytes(b"garbage")
    _, summary = _report_text(d)
    assert "Found 1 fusion networks within the target panel" in summary
    assert GENOME_WIDE_UNREADABLE in summary


def test_i1_target_panel_with_no_trace_of_a_run_is_not_called_none_found(tmp_path):
    # No target-panel pickle or CSV and no sv_count.txt (a pickle write makes one, so it is
    # removed by hand): the target panel shows no sign of a run.
    d = _sample_dir(tmp_path)
    _write_pickles(d, _rows(4), names=("all",))
    (d / "sv_count.txt").unlink()
    _, summary = _report_text(d)
    assert "Fusion analysis has not been run for the target panel." in summary
    assert NONE_FOUND_SUMMARY not in summary
    assert ONE_GENOME_WIDE in summary


def test_i1_control_target_panel_ran_and_found_nothing_is_none_found(tmp_path):
    # The target-panel build ran and nothing passed (an empty pickle): "none found" is right.
    d = _sample_dir(tmp_path)
    _write_pickles(d, _rows(4), names=("all",))
    fusion_work.write_empty_fusion_pickle(str(d / "fusion_candidates_master_processed.pkl"))
    _, summary = _report_text(d)
    assert NONE_FOUND_SUMMARY in summary
    assert ONE_GENOME_WIDE in summary
    assert TARGET_UNREADABLE not in summary and UNREADABLE not in summary


# ------------------------------------------------------------------ I1 follow-ups (a), (b)
def test_i1a_a_target_panel_filter_step_that_raises_reads_as_could_not_be_read(tmp_path, monkeypatch):
    # (a) The target-panel build fails before its candidate CSV is written: only sv_count.txt
    # (and the genome-wide outputs) are left. Through the real writer.
    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    real_iter = fusion_work._iter_fusion_candidates_parquet_batches

    def target_filter_raises(candidate_type, *args, **kwargs):
        if candidate_type == "target_candidates":
            raise RuntimeError("synthetic target-panel filter crash")
        return real_iter(candidate_type, *args, **kwargs)

    monkeypatch.setattr(fusion_work, "_iter_fusion_candidates_parquet_batches", target_filter_raises)
    _stage(tmp_path, _rows(4))
    assert _accumulate(tmp_path)["status"] == "error"
    d = tmp_path / SAMPLE
    assert (d / "sv_count.txt").exists()
    assert not (d / "fusion_candidates_master.csv").exists()
    assert not (d / "fusion_candidates_master_processed.pkl").exists()
    _target_unreadable_genome_wide_shown(d)

    # The next good build removes the record of the failure.
    monkeypatch.setattr(fusion_work, "_iter_fusion_candidates_parquet_batches", real_iter)
    _stage(tmp_path, _rows(1, "GENEC", "GENED", prefix="x"), n=2)
    assert _accumulate(tmp_path)["status"] == "success"
    _, summary = _report_text(d)
    assert "Found 1 fusion networks within the target panel" in summary
    assert TARGET_UNREADABLE not in summary


def test_i1a_a_good_build_that_finds_nothing_after_a_failed_one_reads_as_none_found(
    tmp_path, monkeypatch
):
    # The record of a failed build must go with the next good build, or a run that worked and
    # found nothing would still read as "could not be read".
    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    real_iter = fusion_work._iter_fusion_candidates_parquet_batches

    def target_filter_raises(candidate_type, *args, **kwargs):
        if candidate_type == "target_candidates":
            raise RuntimeError("synthetic target-panel filter crash")
        return real_iter(candidate_type, *args, **kwargs)

    monkeypatch.setattr(fusion_work, "_iter_fusion_candidates_parquet_batches", target_filter_raises)
    _stage(tmp_path, _rows(2))  # below the read support: nothing passes
    assert _accumulate(tmp_path)["status"] == "error"
    d = tmp_path / SAMPLE
    body, _ = _report_text(d)
    assert UNREADABLE in body

    monkeypatch.setattr(fusion_work, "_iter_fusion_candidates_parquet_batches", real_iter)
    _stage(tmp_path, _rows(1, "GENEC", "GENED", prefix="x"), n=2)
    assert _accumulate(tmp_path)["status"] == "success"
    body, summary = _report_text(d)
    assert NONE_FOUND in body and NONE_FOUND_SUMMARY in summary
    assert UNREADABLE not in body and "could not be read" not in summary


TARGET_EARLIER = "The latest target-panel fusion build failed; these results are from an earlier run."
GENOME_EARLIER = "The latest genome-wide fusion build failed; these results are from an earlier run."


def _target_filter_raises(real_iter):
    def target_filter_raises(candidate_type, *args, **kwargs):
        if candidate_type == "target_candidates":
            raise RuntimeError("synthetic target-panel filter crash")
        return real_iter(candidate_type, *args, **kwargs)
    return target_filter_raises


def test_earlier_candidates_after_a_failed_build_say_so(tmp_path, monkeypatch):
    # A good build with candidates, then a build that fails for the target panel, then a good
    # one: through the real writer.
    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    real_iter = fusion_work._iter_fusion_candidates_parquet_batches
    d = tmp_path / SAMPLE
    found = "Found 1 fusion networks within the target panel"

    _stage(tmp_path, _rows(4))
    assert _accumulate(tmp_path)["status"] == "success"
    _, summary = _report_text(d)
    assert found in summary and TARGET_EARLIER not in summary

    monkeypatch.setattr(fusion_work, "_iter_fusion_candidates_parquet_batches",
                        _target_filter_raises(real_iter))
    _stage(tmp_path, _rows(1, "GENEC", "GENED", prefix="x"), n=2)
    assert _accumulate(tmp_path)["status"] == "error"
    _, summary = _report_text(d)
    assert found in summary  # the earlier candidates still show
    assert TARGET_EARLIER in summary
    assert GENOME_EARLIER not in summary  # the genome-wide build worked

    monkeypatch.setattr(fusion_work, "_iter_fusion_candidates_parquet_batches", real_iter)
    _stage(tmp_path, _rows(1, "GENEE", "GENEF", prefix="y"), n=3)
    assert _accumulate(tmp_path)["status"] == "success"
    _, summary = _report_text(d)
    assert found in summary and TARGET_EARLIER not in summary


def test_earlier_result_note_from_the_shared_helper(tmp_path):
    from robin import fusion_status as fs

    d = _sample_dir(tmp_path)
    _write_pickles(d, _rows(4))
    states = fs.fusion_states(str(d))
    assert (states["master_note"], states["all_note"]) == (None, None)
    Path(fs.failed_build_marker(str(d / "fusion_candidates_all_processed.pkl"))).write_text("x\n")
    states = fs.fusion_states(str(d))
    assert (states["master_state"], states["all_state"]) == (fs.STATE_FOUND, fs.STATE_FOUND)
    assert (states["master_note"], states["all_note"]) == (None, GENOME_EARLIER)
    assert fs.earlier_result_note(str(d), "all", fs.STATE_FOUND) == GENOME_EARLIER
    # no note where the record already makes the type "could not be read"
    assert fs.earlier_result_note(str(d), "all", fs.STATE_UNREADABLE) is None
    _, summary = _report_text(d)
    assert GENOME_EARLIER in summary and TARGET_EARLIER not in summary


def test_earlier_result_note_in_the_gui_panel():
    import ast

    tree = _gui_fusion_tree()
    fns = {n.name: n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    load_calls = [
        ast.unparse(c) for c in ast.walk(fns["_load_fusion_data"])
        if isinstance(c, ast.Call) and getattr(c.func, "id", None) == "earlier_result_note"
    ]
    assert any("'master'" in c for c in load_calls) and any("'all'" in c for c in load_calls), load_calls
    shown = [
        ast.unparse(c) for c in ast.walk(fns["_update_fusion_ui"])
        if isinstance(c, ast.Call) and getattr(c.func, "id", None) == "_show_earlier_result_note"
    ]
    # each section shows its note wherever it redraws its candidates (two paths each)
    assert sum("'target'" in c for c in shown) == 2, shown
    assert sum("'genome'" in c for c in shown) == 2, shown


def test_i1b_body_says_which_results_could_not_be_read(tmp_path):
    # (b) The target panel's results can't be read; the genome-wide ones hold candidates
    # (3 reads) but none at a validated breakpoint (4 reads), so no table is drawn.
    d = _sample_dir(tmp_path)
    _write_pickles(d, _rows(3), names=("all",))
    assert not _load(d / "fusion_candidates_all_processed.pkl")["annotated_data"].empty
    (d / "fusion_candidates_master_processed.pkl").write_bytes(b"garbage")
    body, summary = _report_text(d)
    assert "Detailed Analysis" not in body
    assert NONE_FOUND not in body
    assert "Target panel: the results could not be read." in body
    assert "Genome-wide: no fusion candidates were identified." in body
    assert TARGET_UNREADABLE in summary


# ------------------------------------------------------------------ M8: the GUI panel's state
def _type_state(d, kind):
    from robin.fusion_status import candidate_type_state

    return candidate_type_state(str(d), kind)


def test_m8_shared_state_per_type(tmp_path):
    from robin import fusion_status as fs

    d = _sample_dir(tmp_path)
    assert (_type_state(d, "master"), _type_state(d, "all")) == (fs.STATE_NOT_RUN, fs.STATE_NOT_RUN)
    _write_pickles(d, _rows(4), names=("all",))  # (also writes sv_count.txt)
    assert (_type_state(d, "master"), _type_state(d, "all")) == (fs.STATE_NONE, fs.STATE_FOUND)
    (d / "fusion_candidates_master_processed.pkl").write_bytes(b"garbage")
    assert _type_state(d, "master") == fs.STATE_UNREADABLE
    fusion_work.write_empty_fusion_pickle(str(d / "fusion_candidates_master_processed.pkl"))
    assert _type_state(d, "master") == fs.STATE_NONE
    # a failed build's record makes an empty (earlier) pickle unreadable, not "none"
    marker = fs.failed_build_marker(str(d / "fusion_candidates_master_processed.pkl"))
    Path(marker).write_text("RuntimeError: x\n")
    assert _type_state(d, "master") == fs.STATE_UNREADABLE


def test_m8_shared_state_unreadable_shapes(tmp_path):
    from robin import fusion_status as fs

    d = _sample_dir(tmp_path)
    # a pickle that is not the processed structure
    (d / "fusion_candidates_all_processed.pkl").write_bytes(pickle.dumps({"annotated_data": "x"}))
    assert _type_state(d, "all") == fs.STATE_UNREADABLE
    # a candidate CSV that can't be parsed, no pickle
    (d / "fusion_candidates_master.csv").write_text('a,b\n"1,2\n')
    assert _type_state(d, "master") == fs.STATE_UNREADABLE


def test_m8_gui_has_candidates():
    from robin.gui.components.fusion import _has_candidates

    assert _has_candidates(None) is False
    assert _has_candidates({"annotated_data": pd.DataFrame()}) is False
    assert _has_candidates({"annotated_data": _rows(1)}) is True


def test_m8_shared_messages_for_the_panel():
    from robin import fusion_status as fs

    assert fs.status_message("master", fs.STATE_NOT_RUN) == (
        "Target panel fusion analysis has not been run yet")
    assert fs.status_message("all", fs.STATE_NOT_RUN) == (
        "Genome-wide fusion analysis has not been run yet")
    assert fs.status_message("master", fs.STATE_UNREADABLE) == (
        "The target panel fusion results could not be read")
    assert fs.status_message("all", fs.STATE_UNREADABLE) == (
        "The genome-wide fusion results could not be read")
    for kind in ("master", "all"):
        for state in (fs.STATE_NONE, fs.STATE_FOUND):
            assert fs.status_message(kind, state) == "No fusion candidates were identified"


def _gui_fusion_tree():
    import ast

    import robin

    return ast.parse((Path(robin.__file__).parent / "gui" / "components" / "fusion.py").read_text())


def _called_names(node):
    import ast

    names = set()
    for c in ast.walk(node):
        if isinstance(c, ast.Call):
            f = c.func
            names.add(f.id if isinstance(f, ast.Name) else getattr(f, "attr", None))
    return names


def test_m8_gui_panel_uses_the_shared_state():
    import ast

    tree = _gui_fusion_tree()
    fns = {n.name: n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef)}
    # the loader decides each type's state with the shared helper, for both types
    loader = _called_names(fns["_load_fusion_data"])
    assert "candidate_type_state" in loader
    load_src = ast.unparse(fns["_load_fusion_data"])
    assert '"master"' in load_src or "'master'" in load_src
    assert '"all"' in load_src or "'all'" in load_src
    # the panel shows the shared message for both sections
    ui_fn = fns["_update_fusion_ui"]
    msg_calls = [
        ast.unparse(c) for c in ast.walk(ui_fn)
        if isinstance(c, ast.Call) and getattr(c.func, "id", getattr(c.func, "attr", "")) == "status_message"
    ]
    assert any("'master'" in m or '"master"' in m for m in msg_calls), msg_calls
    assert any("'all'" in m or '"all"' in m for m in msg_calls), msg_calls
    # and no longer says "no data" whatever the reason
    strings = {n.value for n in ast.walk(tree) if isinstance(n, ast.Constant) and isinstance(n.value, str)}
    for old in ("No fusion data available", "No fusion candidates found yet",
                "(Fusion data file not found or could not be loaded)"):
        assert old not in strings, old


def test_a_step_inside_the_real_preprocessing_that_raises_fails_the_job(tmp_path, monkeypatch):
    # Pins preprocess_fusion_data_standalone's own re-raise: the heavy step it calls raises,
    # the function itself is real.
    def boom(df):
        raise RuntimeError("synthetic annotate crash")

    monkeypatch.setattr(fusion_work, "_annotate_results", boom)
    _, filed, outcome = _run_handler(tmp_path, monkeypatch)
    assert [e["job_type"] for e in filed] == ["fusion_analysis"]
    assert "synthetic annotate crash" in filed[0]["error"]
    assert outcome["status"] == "failed"
    body, _ = _report_text(tmp_path / "work" / SAMPLE)
    assert UNREADABLE in body


# ------------------------------------------------------------------ second review: C1 and m-1
_TYPE_PICKLE = {
    "master": "fusion_candidates_master_processed.pkl",
    "all": "fusion_candidates_all_processed.pkl",
}
_TYPE_CSV = {"master": "fusion_candidates_master.csv", "all": "fusion_candidates_all.csv"}


def _good_pair_tags(pkl):
    data = _load(pkl)
    ad = data["annotated_data"]
    return sorted(set(ad[data["goodpairs"]]["tag"])) if not ad.empty else []


def _csv_tags(csv):
    return sorted(set(pd.read_csv(csv)["tag"]))


def _good_then_failed_build(tmp_path, monkeypatch, kind):
    """Build 1: GENEA-GENEB, good. Build 2: new valid reads GENEC-GENED, and that type's
    preprocessing raises after its CSV was appended. Returns (real preprocessing, sample dir)."""
    from robin import fusion_status as fs

    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    real = fusion_work.preprocess_fusion_data_standalone

    def raises_for_kind(df, output_file):
        if os.path.basename(output_file) == _TYPE_PICKLE[kind]:
            raise RuntimeError(f"synthetic {kind} preprocess crash")
        return real(df, output_file)

    _stage(tmp_path, _rows(4))
    assert _accumulate(tmp_path)["status"] == "success"
    monkeypatch.setattr(fusion_work, "preprocess_fusion_data_standalone", raises_for_kind)
    _stage(tmp_path, _rows(4, "GENEC", "GENED", prefix="c"), n=2)
    assert _accumulate(tmp_path)["status"] == "error"
    d = tmp_path / SAMPLE
    pkl = d / _TYPE_PICKLE[kind]
    assert os.path.exists(fs.failed_build_marker(str(pkl)))
    assert _good_pair_tags(pkl) == ["GENEA,GENEB"]  # still the earlier build's
    assert _csv_tags(d / _TYPE_CSV[kind]) == ["GENEA,GENEB", "GENEC,GENED"]
    return real, d


@pytest.mark.parametrize("kind", ["master", "all"])
def test_a_quiet_build_after_a_failed_one_rebuilds_from_the_csv(tmp_path, monkeypatch, kind):
    # Second review C1: the next build brings no new valid read (one read, below support).
    # It must not clear the failure while leaving the earlier build's pickle in place.
    from robin import fusion_status as fs

    real, d = _good_then_failed_build(tmp_path, monkeypatch, kind)
    monkeypatch.setattr(fusion_work, "preprocess_fusion_data_standalone", real)
    _stage(tmp_path, _rows(1, "GENEE", "GENEF", prefix="e"), n=3)
    assert _accumulate(tmp_path)["status"] == "success"
    pkl = d / _TYPE_PICKLE[kind]
    assert _good_pair_tags(pkl) == _csv_tags(d / _TYPE_CSV[kind]) == ["GENEA,GENEB", "GENEC,GENED"]
    assert not os.path.exists(fs.failed_build_marker(str(pkl)))
    assert fs.fusion_states(str(d))[f"{kind}_note"] is None
    if kind == "master":
        _, summary = _report_text(d)
        assert "Found 2 fusion networks within the target panel" in summary, summary


@pytest.mark.parametrize("kind", ["master", "all"])
def test_a_quiet_build_while_the_build_still_fails_keeps_the_failure_shown(
    tmp_path, monkeypatch, kind
):
    # Second review C1, the other way: the preprocessing still raises on the quiet build.
    from robin import fusion_status as fs

    _good_then_failed_build(tmp_path, monkeypatch, kind)
    d = tmp_path / SAMPLE
    _stage(tmp_path, _rows(1, "GENEE", "GENEF", prefix="e"), n=3)
    assert _accumulate(tmp_path)["status"] == "error"
    assert os.path.exists(fs.failed_build_marker(str(d / _TYPE_PICKLE[kind])))
    note = TARGET_EARLIER if kind == "master" else GENOME_EARLIER
    assert fs.fusion_states(str(d))[f"{kind}_note"] == note


def test_a_quiet_build_with_no_failure_record_does_not_rebuild(tmp_path, monkeypatch):
    # Control for C1: with no failed build behind it, a quiet build leaves the pickles alone.
    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    real = fusion_work.preprocess_fusion_data_standalone
    _stage(tmp_path, _rows(4))
    assert _accumulate(tmp_path)["status"] == "success"
    d = tmp_path / SAMPLE
    before = {k: (d / n).read_bytes() for k, n in _TYPE_PICKLE.items()}
    calls = []

    def counting(df, output_file):
        calls.append(os.path.basename(output_file))
        return real(df, output_file)

    monkeypatch.setattr(fusion_work, "preprocess_fusion_data_standalone", counting)
    _stage(tmp_path, _rows(1, "GENEE", "GENEF", prefix="e"), n=2)
    assert _accumulate(tmp_path)["status"] == "success"
    assert calls == []
    assert {k: (d / n).read_bytes() for k, n in _TYPE_PICKLE.items()} == before


M1_SENTENCE = (
    "The latest target-panel fusion build failed; "
    "an earlier run identified no fusion networks within the target panel."
)


def test_a_failed_build_over_an_earlier_result_with_no_network_says_so(tmp_path):
    # Second review m-1: candidates (3 reads) but no validated network, latest build failed.
    from robin import fusion_status as fs

    d = _sample_dir(tmp_path)
    _write_pickles(d, _rows(3), names=("master",))
    marker = fs.failed_build_marker(str(d / "fusion_candidates_master_processed.pkl"))
    Path(marker).write_text("RuntimeError: x\n")
    _, summary = _report_text(d)
    assert M1_SENTENCE in summary, summary
    assert NONE_FOUND_SUMMARY not in summary
    assert TARGET_EARLIER not in summary


def test_an_earlier_result_with_no_network_and_no_failed_build_says_none_found(tmp_path):
    # Control for m-1: the same pickle with no failed build still reads "none found".
    d = _sample_dir(tmp_path)
    _write_pickles(d, _rows(3), names=("master",))
    _, summary = _report_text(d)
    assert NONE_FOUND_SUMMARY in summary
    assert M1_SENTENCE not in summary and TARGET_EARLIER not in summary


# ------------------------------------------------------------------ external review, findings 1 and 3
def _spoil_csv(csv, how):
    if how == "missing":
        csv.unlink()
    elif how == "zero_bytes":
        csv.write_text("")
    else:  # header only
        pd.read_csv(csv).iloc[:0].to_csv(csv, index=False)


@pytest.mark.parametrize("how", ["missing", "zero_bytes", "header_only"])
@pytest.mark.parametrize("kind", ["master", "all"])
def test_a_failed_build_whose_csv_cannot_rebuild_keeps_its_record(tmp_path, monkeypatch, kind, how):
    # External review finding 1: after the failed build, the CSV is gone or empty; the quiet
    # build can't rebuild, so it must not clear the record over the stale pickle.
    from robin import fusion_status as fs

    real, d = _good_then_failed_build(tmp_path, monkeypatch, kind)
    monkeypatch.setattr(fusion_work, "preprocess_fusion_data_standalone", real)
    _spoil_csv(d / _TYPE_CSV[kind], how)
    _stage(tmp_path, _rows(1, "GENEE", "GENEF", prefix="e"), n=3)
    assert _accumulate(tmp_path)["status"] == "error"
    pkl = d / _TYPE_PICKLE[kind]
    assert os.path.exists(fs.failed_build_marker(str(pkl)))
    assert _good_pair_tags(pkl) == ["GENEA,GENEB"]
    note = TARGET_EARLIER if kind == "master" else GENOME_EARLIER
    assert fs.fusion_states(str(d))[f"{kind}_note"] == note


def test_a_quiet_build_with_no_record_and_no_csv_still_succeeds(tmp_path, monkeypatch):
    # Control for finding 1: without a failed build behind it, a missing CSV changes nothing.
    from robin import fusion_status as fs

    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    _stage(tmp_path, _rows(4))
    assert _accumulate(tmp_path)["status"] == "success"
    d = tmp_path / SAMPLE
    (d / _TYPE_CSV["master"]).unlink()
    _stage(tmp_path, _rows(1, "GENEE", "GENEF", prefix="e"), n=2)
    assert _accumulate(tmp_path)["status"] == "success"
    assert not os.path.exists(fs.failed_build_marker(str(d / _TYPE_PICKLE["master"])))


def test_the_body_says_a_types_latest_build_failed(tmp_path):
    # External review finding 3: 3-read candidates for both types (none validated) and the
    # target panel's latest build failed.
    from robin import fusion_status as fs

    d = _sample_dir(tmp_path)
    _write_pickles(d, _rows(3))
    Path(fs.failed_build_marker(str(d / _TYPE_PICKLE["master"]))).write_text("RuntimeError: x\n")
    body, _ = _report_text(d)
    assert "Target panel: the latest build failed; an earlier run identified no fusion candidates." in body
    assert "Genome-wide: no fusion candidates were identified." in body
    assert NONE_FOUND not in body


def test_the_body_without_a_failed_build_says_none_found(tmp_path):
    # Control for finding 3.
    d = _sample_dir(tmp_path)
    _write_pickles(d, _rows(3))
    body, _ = _report_text(d)
    assert NONE_FOUND in body
    assert "latest build failed" not in body


# ------------------------------------------------------------------ batch 7b: unreadable candidate files
def _staged(tmp_path, name):
    return Path(fusion_work._get_staging_dir(str(tmp_path), SAMPLE)) / name


def _csv_rows_per_tag(csv):
    return pd.read_csv(csv)["tag"].value_counts().to_dict()


def test_a_damaged_staged_file_fails_and_stays_staged(tmp_path, monkeypatch):
    # 7b test 1: a damaged staged target file next to good genome-wide ones, first build.
    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    _stage(tmp_path, _rows(4), kinds=("genome",))
    _staged(tmp_path, "target_1.parquet").write_bytes(b"PAR1 half written PAR1")
    assert _accumulate(tmp_path)["status"] == "error"
    assert _staged(tmp_path, "target_1.parquet").exists()
    _, summary = _report_text(tmp_path / SAMPLE)
    assert NONE_FOUND_SUMMARY not in summary, summary


def test_a_damaged_staged_file_after_a_good_build_says_the_latest_build_failed(tmp_path, monkeypatch):
    # 7b test 2 (F5): the accumulation fails before its builds; earlier results say so.
    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    _stage(tmp_path, _rows(4))
    assert _accumulate(tmp_path)["status"] == "success"
    _stage(tmp_path, _rows(4, "GENEC", "GENED", prefix="c"), kinds=("genome",), n=2)
    _staged(tmp_path, "target_2.parquet").write_bytes(b"PAR1 half written PAR1")
    assert _accumulate(tmp_path)["status"] == "error"
    _, summary = _report_text(tmp_path / SAMPLE)
    assert TARGET_EARLIER in summary, summary


def _dataset_read_fails(monkeypatch, which_pass, part=None, only_once=False):
    """pq.ParquetFile for the target dataset raises in one pass of the build (first: the
    read_id/col4 pass; second: the full-row pass), optionally for one part only, optionally once."""
    real_pf = fusion_work.pq.ParquetFile
    fired = []

    class Faulty:
        def __init__(self, path, *args, **kwargs):
            self._pf = real_pf(path, *args, **kwargs)
            self._path = str(path)

        def __getattr__(self, name):
            return getattr(self._pf, name)

        def iter_batches(self, columns=None, **kwargs):
            hit = (
                "target_candidates_dataset" in self._path
                and (columns is None) == (which_pass == "second")
                and (part is None or self._path.endswith(part))
                and not (only_once and fired)
            )
            if hit:
                fired.append(self._path)
                raise OSError("synthetic dataset read failure")
            return self._pf.iter_batches(columns=columns, **kwargs)

    monkeypatch.setattr(fusion_work.pq, "ParquetFile", Faulty)
    return fired


@pytest.mark.parametrize("which_pass", ["first", "second"])
def test_a_dataset_read_error_fails_the_types_build(tmp_path, monkeypatch, which_pass):
    # 7b test 3 (F3): the reviewer's case, first build.
    from robin import fusion_status as fs

    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    _stage(tmp_path, _rows(4))
    fired = _dataset_read_fails(monkeypatch, which_pass)
    assert _accumulate(tmp_path)["status"] == "error"
    assert fired
    d = tmp_path / SAMPLE
    assert os.path.exists(fs.failed_build_marker(str(d / _TYPE_PICKLE["master"])))
    assert fs.fusion_states(str(d))["master_state"] != fs.STATE_NONE


def test_a_read_error_part_way_through_the_second_pass_appends_nothing_twice(tmp_path, monkeypatch):
    # 7b test 4 (F4): build 2 reads part 2, then fails on part 3 in the second pass; build 3
    # retries. Each read's rows must reach the CSV once.
    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    _stage(tmp_path, _rows(4))
    assert _accumulate(tmp_path)["status"] == "success"
    _stage(tmp_path, _rows(4, "GENEC", "GENED", prefix="c"), n=2)
    _stage(tmp_path, _rows(4, "GENEG", "GENEH", prefix="g"), n=3)
    fired = _dataset_read_fails(monkeypatch, "second", part="part_000003.parquet", only_once=True)
    _accumulate(tmp_path)
    assert fired
    _stage(tmp_path, _rows(1, "GENEE", "GENEF", prefix="e"), n=4)
    assert _accumulate(tmp_path)["status"] == "success"
    rows = _csv_rows_per_tag(tmp_path / SAMPLE / _TYPE_CSV["master"])
    assert rows == {"GENEA,GENEB": 8, "GENEC,GENED": 8, "GENEG,GENEH": 8}, rows


def test_after_a_damaged_staged_file_is_repaired_each_read_counts_once(tmp_path, monkeypatch):
    # 7b test 5 (F2): the genome-wide staged file is damaged, the target one is good. The failed
    # accumulation must append nothing, so the retry appends each read once.
    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    _stage(tmp_path, _rows(4), kinds=("target",))
    _staged(tmp_path, "genome_1.parquet").write_bytes(b"PAR1 half written PAR1")
    assert _accumulate(tmp_path)["status"] == "error"
    _rows(4).to_parquet(_staged(tmp_path, "genome_1.parquet"), index=False)
    assert _accumulate(tmp_path)["status"] == "success"
    d = tmp_path / SAMPLE
    for kind in ("master", "all"):
        assert _csv_rows_per_tag(d / _TYPE_CSV[kind]) == {"GENEA,GENEB": 8}, kind


def test_parquet_writes_go_through_a_temporary_name(tmp_path, monkeypatch):
    # 7b test 6 (F1): a reader never sees half a staged file or dataset part.
    import inspect

    calls = []
    real_replace = os.replace

    def recording_replace(a, b):
        calls.append((str(a), str(b)))
        return real_replace(a, b)

    monkeypatch.setattr(fusion_work.os, "replace", recording_replace)
    path = tmp_path / "x.parquet"
    fusion_work._write_parquet_atomically(_rows(1), str(path), index=False)
    assert calls == [(str(path) + ".tmp", str(path))]
    assert len(pd.read_parquet(path)) == 2
    for fn in (fusion_work._append_fusion_candidates_parquet, fusion_work.process_bam_with_staging):
        src = inspect.getsource(fn)
        assert "_write_parquet_atomically(" in src and ".to_parquet(" not in src, fn.__name__
