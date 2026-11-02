"""Batch 10e: a BAM whose coverage the target step could not add shows as a failed file.

Before, an unreadable staged set, a failed commit on the one-BAM path, a batch BAM that
could not be written and a BAM that could not be staged in a batch were only logged: the
step said "success" and the coverage quietly missed that BAM. Synthetic BAMs only.
"""

import glob
import os
from types import SimpleNamespace

import pytest

from target_fixtures import SAMPLE, write_bams, write_panel
from robin import workflow_ray as wr
from robin.analysis import target_analysis as ta_mod
from robin.analysis.target_analysis import TargetAnalysis, target_handler

PANEL_NAME = "TESTPANEL"
STEP = "target_coverage"


@pytest.fixture
def panel(tmp_path, monkeypatch):
    path = write_panel(str(tmp_path / "panel.bed"))
    monkeypatch.setattr(TargetAnalysis, "_find_target_bed", lambda self, p: path)
    return path


def _batch_size_one(monkeypatch):
    """The one-BAM path adds staged BAMs every 10; make it every one."""
    real_init = TargetAnalysis.__init__

    def init(self, *args, **kwargs):
        kwargs["batch_size"] = 1
        real_init(self, *args, **kwargs)

    monkeypatch.setattr(TargetAnalysis, "__init__", init)


def _one_bam_job(bam):
    ctx = wr.WorkflowContext(filepath=bam, metadata={"bam_metadata": {"sample_id": SAMPLE}})
    return wr.Job(1, "target", "target", ["target:target"], 0, ctx)


def _batched_job(bams):
    inner = [wr.WorkflowContext(filepath=b, metadata={"bam_metadata": {"sample_id": SAMPLE}}) for b in bams]
    batched = SimpleNamespace(
        get_file_count=lambda: len(bams), get_sample_id=lambda: SAMPLE, batch_id="b1",
        get_filepaths=lambda: list(bams), contexts=inner,
    )
    ctx = wr.WorkflowContext(filepath=bams[0], metadata={"_batched_job": batched})
    return wr.Job(2, "target", "target", ["target:target"], 0, ctx)


def _failures(job):
    return [(f["file"], f["step"]) for f in job.context.metadata.get("file_failures") or []]


def _stage(work, bam):
    """Stage a BAM without adding it (the one-BAM path's batch size is 10)."""
    ta = TargetAnalysis(work_dir=str(work), target_panel=PANEL_NAME, batch_size=10, use_staging=True)
    _, should_accumulate = ta.process_file_with_staging(bam, {"sample_id": SAMPLE})
    assert not should_accumulate
    staged = glob.glob(os.path.join(str(work), SAMPLE, "_staging", "coverage_*.parquet"))
    assert len(staged) == 1
    return staged[0]


def test_one_bam_path_an_unreadable_staged_set_is_a_failed_file(tmp_path, panel, monkeypatch):
    bams, _ = write_bams(str(tmp_path), 2)
    work = tmp_path / "work"
    with open(_stage(work, bams[0]), "wb") as f:
        f.write(b"not a parquet file")
    _batch_size_one(monkeypatch)
    job = _one_bam_job(bams[1])
    target_handler(job, work_dir=str(work), target_panel=PANEL_NAME)
    assert _failures(job) == [(bams[0], STEP)]
    assert "coverage" in job.context.metadata["file_failures"][0]["error"]
    assert job.context.errors == []
    assert job.context.results["target_analysis"]["status"] == "success"


def test_one_bam_path_a_failed_commit_is_a_failed_file(tmp_path, panel, monkeypatch):
    bams, _ = write_bams(str(tmp_path), 1)
    work = tmp_path / "work"

    def failing_commit(self, *args, **kwargs):
        raise TimeoutError("lock held")

    monkeypatch.setattr(TargetAnalysis, "_commit_batch", failing_commit)
    _batch_size_one(monkeypatch)
    job = _one_bam_job(bams[0])
    target_handler(job, work_dir=str(work), target_panel=PANEL_NAME)
    assert _failures(job) == [(bams[0], STEP)]
    assert "lock held" in job.context.metadata["file_failures"][0]["error"]
    assert job.context.errors == []


def test_all_staged_sets_unreadable_names_each_bam(tmp_path, panel):
    bams, _ = write_bams(str(tmp_path), 1)
    work = tmp_path / "work"
    with open(_stage(work, bams[0]), "wb") as f:
        f.write(b"not a parquet file")
    ta = TargetAnalysis(work_dir=str(work), target_panel=PANEL_NAME, batch_size=1, use_staging=True)
    result = ta.accumulate_staged_files(SAMPLE, force=True)
    assert result["status"] == "load_failed"
    assert [(f["file"], f["step"]) for f in result.get("file_failures") or []] == [(bams[0], STEP)]


def test_batched_path_a_batch_bam_that_cannot_be_written_is_a_failed_file(tmp_path, panel, monkeypatch):
    bams, _ = write_bams(str(tmp_path), 2)
    work = tmp_path / "work"

    def failing_write(self, *args, **kwargs):
        raise RuntimeError("disk full")

    monkeypatch.setattr(TargetAnalysis, "_write_batch_bam", failing_write)
    job = _batched_job(bams)
    target_handler(job, work_dir=str(work), target_panel=PANEL_NAME)
    assert sorted(_failures(job)) == sorted((b, STEP) for b in bams)
    for f in job.context.metadata["file_failures"]:
        assert "target.bam" in f["error"] and "disk full" in f["error"]
    assert job.context.errors == []


def test_one_bam_path_a_batch_bam_that_cannot_be_moved_into_place_is_a_failed_file(tmp_path, panel, monkeypatch):
    bams, _ = write_bams(str(tmp_path), 1)
    work = tmp_path / "work"
    real_replace = os.replace

    def replace(src, dst, *args, **kwargs):
        name = os.path.basename(str(dst))
        if name.startswith("batch_") and name.endswith(".bam"):
            raise OSError("rename failed")
        return real_replace(src, dst, *args, **kwargs)

    monkeypatch.setattr(ta_mod.os, "replace", replace)
    _batch_size_one(monkeypatch)
    job = _one_bam_job(bams[0])
    target_handler(job, work_dir=str(work), target_panel=PANEL_NAME)
    assert _failures(job) == [(bams[0], STEP)]
    assert "target.bam" in job.context.metadata["file_failures"][0]["error"]
    assert job.context.errors == []
    claims = glob.glob(os.path.join(str(work), SAMPLE, "_staging", "_batch_*"))
    assert len(claims) == 1, "the batch's folder must stay, listed as reads not in target.bam"


def test_batched_path_a_bam_that_cannot_be_staged_is_a_failed_file(tmp_path, panel):
    bams, _ = write_bams(str(tmp_path), 1)
    bad = tmp_path / "bad.bam"
    bad.write_text("not a BAM")
    work = tmp_path / "work"
    job = _batched_job([bams[0], str(bad)])
    target_handler(job, work_dir=str(work), target_panel=PANEL_NAME)
    assert _failures(job) == [(str(bad), STEP)]
    assert job.context.errors == []
    assert job.context.results["target_analysis"]["status"] == "success"


def test_finalize_shows_a_bam_its_accumulation_could_not_add(tmp_path, panel):
    """Found after the lock: finalize adds the last staged sets too, and only kept the result."""
    bams, _ = write_bams(str(tmp_path), 1)
    work = tmp_path / "work"
    with open(_stage(work, bams[0]), "wb") as f:
        f.write(b"not a parquet file")
    ctx = wr.WorkflowContext(
        filepath=str(work / SAMPLE),
        metadata={"sample_id": SAMPLE, "work_dir": str(work), "target_panel": PANEL_NAME},
    )
    job = wr.Job(3, "target_bam_finalize", "target_bam_finalize", ["x:target_bam_finalize"], 0, ctx)
    ta_mod.target_bam_finalize_handler(job, work_dir=str(work))
    assert _failures(job) == [(bams[0], STEP)]


# --- round 2 (the external review): whose failure, and when a set is taken -----------------------


def _ledger_with_ended_target_steps(work, bams):
    """A work folder whose ledger says each BAM's target step ended well (as after its job)."""
    from robin.bam_identity import bam_key
    from robin.file_ledger import FileLedger

    records = []
    for bam in bams:
        k = bam_key(bam)
        records += [
            {"e": "take", "key": k, "file": bam},
            {"e": "start", "key": k, "step": "preprocessing", "sample": SAMPLE, "file": bam},
            {"e": "end", "key": k, "step": "preprocessing", "status": "success", "outcome": "success",
             "downstream": True},
            {"e": "start", "key": k, "step": "target", "sample": SAMPLE, "file": bam},
            {"e": "end", "key": k, "step": "target", "status": "success", "outcome": "success"},
        ]
    os.makedirs(str(work), exist_ok=True)
    FileLedger(str(work)).append(records)


def _rows_at_next_start(work, bam):
    from robin.file_ledger import BamIntake

    return [
        (r["job_type"], r["failure_kind"])
        for r in BamIntake(str(work), ["fast:preprocessing", "fast:target"]).decide(bam).cut_short
    ]


@pytest.mark.parametrize("found_by", ["the_next_target_job", "finalize"])
def test_an_earlier_jobs_bam_left_out_is_reported_at_every_start(tmp_path, panel, monkeypatch, found_by):
    bams, _ = write_bams(str(tmp_path), 2)
    work = tmp_path / "work"
    _ledger_with_ended_target_steps(work, bams)
    with open(_stage(work, bams[0]), "wb") as f:  # BAM 0's job staged it and ended well
        f.write(b"not a parquet file")
    if found_by == "finalize":
        ta_mod.finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
    else:
        _batch_size_one(monkeypatch)
        job = _one_bam_job(bams[1])
        target_handler(job, work_dir=str(work), target_panel=PANEL_NAME)
        assert _failures(job) == [(bams[0], STEP)]
    for _ in range(2):  # at every start
        assert _rows_at_next_start(work, bams[0]) == [("target", "failed_before")]
    assert _rows_at_next_start(work, bams[1]) == []


@pytest.mark.parametrize("case", ["no_ledger", "bam_not_in_ledger"])
def test_no_ledger_line_is_written_where_the_ledger_has_no_such_bam(tmp_path, panel, case):
    """A failure line alone would make the next start skip a target step that never ran."""
    from robin.bam_identity import bam_key

    bams, _ = write_bams(str(tmp_path), 2)
    work = tmp_path / "work"
    if case == "bam_not_in_ledger":
        _ledger_with_ended_target_steps(work, bams[1:])
    with open(_stage(work, bams[0]), "wb") as f:
        f.write(b"not a parquet file")
    ta = TargetAnalysis(work_dir=str(work), target_panel=PANEL_NAME, batch_size=1, use_staging=True)
    assert ta.accumulate_staged_files(SAMPLE, force=True)["status"] == "load_failed"
    ledger = work / "robin_bam_ledger.jsonl"
    if case == "no_ledger":
        assert not ledger.exists()
    else:
        assert bam_key(bams[0]) not in ledger.read_text(), "a line for a BAM the ledger never took in"


def test_a_bam_renamed_after_staging_is_reported_under_its_new_name(tmp_path, panel):
    """The plan check's case: the staged set names the old path, which is gone by the time
    its coverage fails. The key saved at staging still ties the failure to the BAM."""
    bams, _ = write_bams(str(tmp_path), 1)
    work = tmp_path / "work"
    _ledger_with_ended_target_steps(work, bams)
    with open(_stage(work, bams[0]), "wb") as f:  # BAM 0's job staged it and ended well
        f.write(b"not a parquet file")
    moved = str(tmp_path / "renamed.bam")
    os.rename(bams[0], moved)  # the same bytes, so the same key
    ta_mod.finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
    assert _rows_at_next_start(work, moved) == [("target", "failed_before")]


def test_a_failed_first_move_leaves_the_set_staged_and_reports_nothing(tmp_path, panel, monkeypatch):
    bams, _ = write_bams(str(tmp_path), 1)
    work = tmp_path / "work"
    staged = _stage(work, bams[0])
    real_move = ta_mod.shutil.move

    def move(src, dst, *args, **kwargs):
        raise OSError("synthetic: the first move fails")

    monkeypatch.setattr(ta_mod.shutil, "move", move)
    ta = TargetAnalysis(work_dir=str(work), target_panel=PANEL_NAME, batch_size=1, use_staging=True)
    result = ta.accumulate_staged_files(SAMPLE, force=True)
    assert result["status"] == "error"
    assert not result.get("file_failures")
    assert os.path.exists(staged)
    assert glob.glob(os.path.join(str(work), SAMPLE, "_staging", "_batch_*")) == [], "an empty claim folder was kept"
    monkeypatch.setattr(ta_mod.shutil, "move", real_move)
    again = ta.accumulate_staged_files(SAMPLE, force=True)
    assert again["status"] == "success" and again["files_processed"] == 1


# --- guards: no failed file where nothing failed -------------------------------------------------


def test_clean_runs_report_no_failed_file(tmp_path, panel, monkeypatch):
    bams, _ = write_bams(str(tmp_path), 3)
    batched = _batched_job(bams[:2])
    target_handler(batched, work_dir=str(tmp_path / "work_a"), target_panel=PANEL_NAME)
    assert _failures(batched) == [] and batched.context.errors == []
    _batch_size_one(monkeypatch)
    single = _one_bam_job(bams[2])
    target_handler(single, work_dir=str(tmp_path / "work_b"), target_panel=PANEL_NAME)
    assert _failures(single) == [] and single.context.errors == []


def test_a_lock_timeout_before_the_claim_leaves_the_files_staged_and_reports_nothing(tmp_path, panel, monkeypatch):
    bams, _ = write_bams(str(tmp_path), 1)
    work = tmp_path / "work"
    staged = _stage(work, bams[0])

    def lock_held(staging_dir):
        raise TimeoutError("lock held")

    monkeypatch.setattr(ta_mod, "_complete_staged_counters", lock_held)
    ta = TargetAnalysis(work_dir=str(work), target_panel=PANEL_NAME, batch_size=1, use_staging=True)
    result = ta.accumulate_staged_files(SAMPLE, force=True)
    assert result["status"] == "lock_timeout"
    assert not result.get("file_failures")
    assert os.path.exists(staged), "the staged set left the staging folder"


def test_a_claim_that_fails_before_any_file_moves_reports_nothing(tmp_path, panel, monkeypatch):
    """A set counts as taken only once its files begin to move: before that it is still
    staged, and the next accumulation takes it."""
    bams, _ = write_bams(str(tmp_path), 1)
    work = tmp_path / "work"
    staged = _stage(work, bams[0])
    real_makedirs = os.makedirs

    def makedirs(path, *args, **kwargs):
        if os.path.basename(str(path)).startswith("_batch_"):
            raise OSError("no space left")
        return real_makedirs(path, *args, **kwargs)

    monkeypatch.setattr(ta_mod.os, "makedirs", makedirs)
    ta = TargetAnalysis(work_dir=str(work), target_panel=PANEL_NAME, batch_size=1, use_staging=True)
    result = ta.accumulate_staged_files(SAMPLE, force=True)
    assert result["status"] == "error"
    assert not result.get("file_failures")
    assert os.path.exists(staged), "the staged set left the staging folder"
