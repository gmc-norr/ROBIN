"""Each BAM is taken in once, and after a restart only its unstarted steps run (batch 10d-2:
JOB-2 with the ledger, JOB-8).

Synthetic BAMs (``target_fixtures``). The start-up scan and the watcher run as written, with a
stand-in Coordinator that records what is sent; the steps run through the real worker wrappers
(``_wrap_real_handler`` and ``WorkflowManager.worker``), the preprocessing step itself for real.

- A BAM taken in before (by its bytes, so a copy too) is not taken in again: not by the next
  session's scan, not by the watcher after the scan.
- After a restart, a BAM's steps that never started run, through its preprocessing (which counts
  it once, batch 10d-1) told which steps to start. A step that was running when ROBIN stopped is
  not run again: it may have added its data. It is reported once as a failed file.
- A BAM still being written is taken in once it is complete (it ends with the BAM end-of-file
  block).
- The watcher's BAMs get the reference genome, as the scan's do (JOB-8).
"""

import ast
import asyncio
import inspect
import os
import queue
import shutil
import threading
from itertools import count

import pytest

pytest.importorskip("pysam")

from target_fixtures import fixture_reads, write_bam  # noqa: E402
from robin import workflow_ray as wr  # noqa: E402
from robin import workflow_simple as ws  # noqa: E402
from robin.analysis.bam_preprocessor import bam_preprocessing_handler  # noqa: E402

_ids = count(70_000_000)
STEPS = ["bed_conversion", "mgmt", "cnv", "target", "fusion"]  # the per-BAM steps after preprocessing
PLAN = [f"{wr.job_queue_of(t)}:{t}" for t in ["preprocessing"] + STEPS]
REFERENCE = "/ref/genome.fa"
CoordinatorClass = wr.Coordinator.__ray_metadata__.modified_class


# A stand-in Coordinator --------------------------------------------------------------------


class _Done:
    """What a Ray ``.remote()`` call gives back: awaitable (the scan) or for ``ray.get`` (the watcher)."""

    def __init__(self, value=None):
        self.value = value

    def __await__(self):
        if False:
            yield
        return self.value


class _Remote:
    def __init__(self, fn):
        self.fn = fn

    def remote(self, *args, **kwargs):
        return _Done(self.fn(*args, **kwargs))


class FakeCoord:
    def __init__(self):
        self.jobs, self.cut_short = [], []
        self.get_reference = _Remote(lambda: REFERENCE)
        self.get_target_panel = _Remote(lambda: "rCNS2")
        self.can_accept_jobs = _Remote(lambda n: True)
        self.submit_jobs = _Remote(self.jobs.extend)
        self.record_cut_short_steps = _Remote(self.cut_short.extend)


@pytest.fixture(autouse=True)
def ray_get(monkeypatch):
    monkeypatch.setattr(wr.ray, "get", lambda x, *a, **k: x.value if isinstance(x, _Done) else x)


# BAMs, sessions and steps ------------------------------------------------------------------


def _bam(directory, name="S_pass_0.bam", part=0):
    os.makedirs(directory, exist_ok=True)
    return write_bam(os.path.join(str(directory), name), fixture_reads()[part::2])


def _scan(coord, folder, work, intake=None):
    """One session's start-up scan of ``folder``; the jobs it sent."""
    before = len(coord.jobs)
    kwargs = {"intake": intake} if intake is not None else {}
    asyncio.run(wr.submit_existing_paths(
        coord, [str(folder)], PLAN, patterns=["*.bam"], work_dir=str(work), **kwargs
    ))
    return coord.jobs[before:]


def _noop(job):
    pass


def _step_job(job, step):
    q = wr.job_queue_of(step)
    return wr.Job(next(_ids), step, q, [f"{q}:{step}"], 0, job.context)


def _work_through(job, steps=STEPS, preprocessing=bam_preprocessing_handler):
    """What the workers do for a BAM: its preprocessing (for real), then each of ``steps``."""
    wr._wrap_real_handler(preprocessing, "preprocessing")(job)
    for step in steps:
        wr._wrap_real_handler(_noop, step)(_step_job(job, step))


def _types(jobs):
    return [j.job_type for j in jobs]


# The Ray runner ----------------------------------------------------------------------------


def test_a_bam_scanned_again_after_a_restart_is_not_taken_again(tmp_path):
    _bam(tmp_path / "run")
    first = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    assert _types(first) == ["preprocessing"]
    _work_through(first[0])
    assert _scan(FakeCoord(), tmp_path / "run", tmp_path / "work") == [], "a finished BAM was taken in again"


def test_different_bams_in_one_scan_are_each_taken_in(tmp_path):
    # Control: the intake takes in every new BAM, with the reference (as the scan did before).
    a = _bam(tmp_path / "run", "S_pass_0.bam", 0)
    b = _bam(tmp_path / "run", "S_pass_1.bam", 1)
    jobs = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    assert sorted(j.context.filepath for j in jobs) == sorted([a, b])
    assert {j.context.metadata.get("reference") for j in jobs} == {REFERENCE}


def test_a_copy_of_a_finished_bam_is_not_taken_again(tmp_path):
    bam = _bam(tmp_path / "run")
    (tmp_path / "copy").mkdir()
    shutil.copyfile(bam, tmp_path / "copy" / "S_pass_0.bam")
    first = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    _work_through(first[0])
    assert _scan(FakeCoord(), tmp_path / "copy", tmp_path / "work") == [], "a copy of a finished BAM was taken in"


def test_a_bam_seen_by_the_scan_and_then_the_watcher_is_taken_once(tmp_path):
    bam = _bam(tmp_path / "run")
    coord = FakeCoord()
    watcher = wr.RayFileWatcher(coord, PLAN, "rCNS2", patterns=["*.bam"], work_dir=str(tmp_path / "work"))
    _scan(coord, tmp_path / "run", tmp_path / "work", intake=watcher.intake)
    watcher._handle(bam)
    watcher.flush_remaining()
    assert _types(coord.jobs) == ["preprocessing"], "the scan and the watcher both took the BAM in"


def test_after_a_restart_only_the_steps_a_bam_never_started_run(tmp_path):
    _bam(tmp_path / "run")
    first = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    _work_through(first[0], steps=["cnv", "target"])  # stopped before the other three started
    again = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    assert _types(again) == ["preprocessing"]
    assert again[0].context.metadata.get("only_steps") == ["bed_conversion", "mgmt", "fusion"]


def test_a_step_cut_short_by_a_stop_is_not_run_again_and_is_reported_at_each_start(tmp_path):
    # Round 2 ruling (the user's choice): reported at every start, not once.
    from robin.file_ledger import record_step_start

    bam = _bam(tmp_path / "run")
    first = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    _work_through(first[0], steps=["bed_conversion", "mgmt", "target", "fusion"])
    record_step_start(_step_job(first[0], "cnv"), "cnv")  # ROBIN stopped while CNV ran
    second = FakeCoord()
    assert _scan(second, tmp_path / "run", tmp_path / "work") == [], "a step cut short was run again"
    assert [(r["job_type"], r["filepath"]) for r in second.cut_short] == [("cnv", bam)]
    assert "may be incomplete" in second.cut_short[0]["error"]
    third = FakeCoord()
    assert _scan(third, tmp_path / "run", tmp_path / "work") == [], "a step cut short was run again"
    assert [(r["job_type"], r["filepath"]) for r in third.cut_short] == [("cnv", bam)], "not reported at the next start"


def test_a_bam_whose_preprocessing_was_cut_short_runs_again(tmp_path):
    from robin.file_ledger import record_step_start

    _bam(tmp_path / "run")
    first = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    record_step_start(first[0], "preprocessing")  # stopped during preprocessing (it counts a BAM once)
    again = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    assert _types(again) == ["preprocessing"]
    assert again[0].context.metadata.get("only_steps") == STEPS


def test_a_bam_its_preprocessing_skipped_runs_nothing_more(tmp_path):
    def skipping(job):  # as the real step does for a BAM over 50,000 reads
        job.context.add_metadata("skip_downstream", True)
        job.context.add_result("preprocessing", {"status": "skipped", "reason": "too_many_reads"})

    _bam(tmp_path / "run")
    first = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    _work_through(first[0], steps=[], preprocessing=skipping)
    assert _scan(FakeCoord(), tmp_path / "run", tmp_path / "work") == [], "a skipped BAM was taken in again"


def test_a_bam_still_being_written_is_taken_in_once_it_is_complete(tmp_path):
    whole = _bam(tmp_path / "made")
    (tmp_path / "run").mkdir()
    path = str(tmp_path / "run" / "S_pass_0.bam")
    with open(whole, "rb") as f:
        data = f.read()
    with open(path, "wb") as f:
        f.write(data[: len(data) - 100])  # the copy is not finished: no end-of-file block yet
    coord = FakeCoord()
    watcher = wr.RayFileWatcher(coord, PLAN, "rCNS2", patterns=["*.bam"], work_dir=str(tmp_path / "work"))
    watcher._handle(path)
    watcher.flush_remaining()
    assert coord.jobs == [], "a BAM still being written was taken in"
    with open(path, "wb") as f:
        f.write(data)
    watcher._handle(path)
    watcher.flush_remaining()
    assert _types(coord.jobs) == ["preprocessing"]


def test_the_watcher_gives_its_bams_the_reference(tmp_path):
    bam = _bam(tmp_path / "run")
    coord = FakeCoord()
    watcher = wr.RayFileWatcher(coord, PLAN, "rCNS2", patterns=["*.bam"], work_dir=str(tmp_path / "work"))
    watcher._handle(bam)
    watcher.flush_remaining()
    assert coord.jobs[0].context.metadata.get("reference") == REFERENCE


def test_the_worker_records_each_bam_of_a_batch_and_a_failed_step(tmp_path):
    from robin.file_ledger import FileLedger

    work = str(tmp_path / "work")
    contexts = [
        wr.WorkflowContext(filepath=f"/data/S_pass_{i}.bam", metadata={
            "bam_key": f"k{i}", "work_dir": work, "bam_metadata": {"sample_id": "S"},
        })
        for i in range(2)
    ]
    q = wr.job_queue_of("cnv")
    batch = wr.BatchedJob(next(_ids), "cnv", q, [f"{q}:cnv"], 0, contexts, "b1", "S")
    job = wr.Job(batch.job_id, "cnv", q, [f"{q}:cnv"], 0, contexts[0])
    job.context.metadata["_batched_job"] = batch
    wr._wrap_real_handler(_noop, "cnv")(job)

    def broken(job):
        raise RuntimeError("synthetic CNV failure")

    fusion = wr.Job(next(_ids), "fusion", wr.job_queue_of("fusion"), ["x:fusion"], 0, contexts[1])
    with pytest.raises(RuntimeError):
        wr._wrap_real_handler(broken, "fusion")(fusion)
    state = FileLedger(work).load()
    assert {k: (state[k].steps["cnv"].ended, state[k].steps["cnv"].status) for k in ("k0", "k1")} == {
        "k0": (True, "success"), "k1": (True, "success"),
    }
    assert (state["k1"].steps["fusion"].ended, state["k1"].steps["fusion"].status) == (True, "failed")


def test_a_step_whose_start_cannot_be_recorded_does_not_run(tmp_path, monkeypatch):
    from robin import file_ledger

    ran = []

    def disk_full(self, records):
        raise OSError("disk full")

    monkeypatch.setattr(file_ledger.FileLedger, "append", disk_full)
    ctx = wr.WorkflowContext(filepath="/data/S_pass_0.bam", metadata={"bam_key": "k0", "work_dir": str(tmp_path)})
    job = wr.Job(next(_ids), "cnv", wr.job_queue_of("cnv"), ["x:cnv"], 0, ctx)
    with pytest.raises(Exception):
        wr._wrap_real_handler(lambda j: ran.append(j), "cnv")(job)
    assert ran == [], "the step ran though the ledger could not record its start"


def test_a_preprocessing_told_its_steps_starts_only_those(tmp_path):
    c = CoordinatorClass(target_panel="rCNS2", preset="standard", enable_batching=False)
    c.using_pools = True
    c._spool_enabled = False
    c.set_work_dir(str(tmp_path))
    sent = []
    proc = type("P", (), {"enqueue": _Remote(sent.append)})()
    proc.enqueue = type("E", (), {"remote": staticmethod(sent.append)})()
    c.processors = {t: proc for t in STEPS}
    ctx = wr.WorkflowContext(filepath=str(tmp_path / "S_pass_0.bam"), metadata={
        "bam_metadata": {"sample_id": "S"}, "original_workflow": PLAN, "only_steps": ["cnv"],
        "work_dir": str(tmp_path),
    })
    ctx.add_result("preprocessing", {"status": "success", "sample_id": "S"})
    job = wr.Job(next(_ids), "preprocessing", wr.job_queue_of("preprocessing"), PLAN, 0, ctx)
    asyncio.run(c._on_finish(job, True, ctx, None))
    assert _types(sent) == ["cnv"], "steps the BAM had finished were started again"


def test_the_coordinator_reports_a_cut_short_step_as_a_failed_file(tmp_path):
    c = CoordinatorClass(target_panel="rCNS2", preset="standard", enable_batching=False)
    c.set_work_dir(str(tmp_path))
    c.record_cut_short_steps([{
        "sample_id": "S", "filepath": "/data/S_pass_0.bam", "job_type": "cnv",
        "error": "cnv was running when ROBIN stopped: it may be incomplete for this BAM",
    }])
    assert [(r["job_type"], r["failure_kind"]) for r in c._failed_jobs_log] == [("cnv", "cut_short")]
    assert [e["kind"] for e in c._events] == ["file_failed"]
    assert c.file_failed_by_type == {"cnv": 1}


def test_the_ray_run_watches_before_it_scans():
    tree = ast.parse(inspect.getsource(wr.run))
    lines = {}
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            name = getattr(node.func, "attr", getattr(node.func, "id", ""))
            lines.setdefault(name, []).append(node.lineno)
    assert min(lines["start"]) < min(lines["submit_existing_paths"]), "the scan runs before watching starts"


def test_a_ledger_that_cannot_be_read_stops_the_intake_and_a_torn_last_line_is_ignored(tmp_path):
    from robin.file_ledger import LEDGER_NAME, BamIntake, LedgerError

    good = '{"v": 1, "e": "take", "key": "k0", "file": "/data/a.bam"}\n'
    (tmp_path / LEDGER_NAME).write_text(good + '{"v": 1, "e": "start", "key": "k0", "st')  # torn at a stop
    BamIntake(str(tmp_path), PLAN)  # reads
    (tmp_path / LEDGER_NAME).write_text(good + "not json\n" + good)
    with pytest.raises(LedgerError):
        BamIntake(str(tmp_path), PLAN)


# Round 2: the external review's findings ---------------------------------------------------


def test_a_torn_last_line_is_cut_off_before_the_next_record(tmp_path):
    # B1: the next record must not be glued to the torn one.
    from robin.file_ledger import FileLedger

    led = FileLedger(str(tmp_path))
    led.append([{"e": "take", "key": "k0", "file": "/data/a.bam"}])
    with open(led.path, "ab") as f:
        f.write(b'{"e": "start", "key": "k0", "st')  # a stop part-way through a record
    led.load()
    led.append([{"e": "take", "key": "k1", "file": "/data/b.bam"}])  # the next session takes a BAM in
    assert sorted(led.load()) == ["k0", "k1"]


def _ctx(i, work=None, module=wr):
    metadata = {"bam_key": f"k{i}", "bam_metadata": {"sample_id": "S"}}
    if work is not None:
        metadata["work_dir"] = str(work)
    return module.WorkflowContext(filepath=f"/data/S_pass_{i}.bam", metadata=metadata)


def _ray_jobs(step, contexts):
    q = wr.job_queue_of(step)
    return [wr.Job(next(_ids), step, q, [f"{q}:{step}"], 0, c) for c in contexts]


def _simple_jobs(step, contexts):
    return [ws.Job(next(_ids), step, c, "fast", [f"fast:{step}"]) for c in contexts]


def test_two_batches_that_share_a_bam_keep_their_own_lists(tmp_path):
    # Finding 1: a BAM's steps share its context; a later batch must not change an earlier one's list.
    from robin.file_ledger import FileLedger

    a, b, c = (_ctx(i, tmp_path) for i in range(3))
    batcher = wr.SampleJobBatcher()
    coord = CoordinatorClass(target_panel="rCNS2", preset="standard", enable_batching=False)
    (cnv_job,) = coord._convert_batched_jobs([batcher._create_batched_job(_ray_jobs("cnv", [a, b]), "S", "cnv")])
    coord._convert_batched_jobs([batcher._create_batched_job(_ray_jobs("mgmt", [a, c]), "S", "mgmt")])
    listed = [x.filepath for x in cnv_job.context.metadata["_batched_job"].contexts]
    assert listed == [a.filepath, b.filepath], "the MGMT batch changed the CNV batch's BAMs"
    wr._wrap_real_handler(_noop, "cnv")(cnv_job)
    state = FileLedger(str(tmp_path)).load()
    assert sorted(k for k in state if "cnv" in state[k].steps) == ["k0", "k1"]


def test_two_non_ray_batches_that_share_a_bam_keep_their_own_lists():
    a, b, c = (_ctx(i, module=ws) for i in range(3))
    batcher = ws.SampleJobBatcher()
    mgr = ws.WorkflowManager(enable_batching=False)
    (cnv_job,) = mgr._convert_batched_jobs([batcher._create_batched_job(_simple_jobs("cnv", [a, b]), "S", "cnv")])
    mgr._convert_batched_jobs([batcher._create_batched_job(_simple_jobs("mgmt", [a, c]), "S", "mgmt")])
    listed = [x.filepath for x in cnv_job.context.metadata["_batched_job"].contexts]
    assert listed == [a.filepath, b.filepath], "the MGMT batch changed the CNV batch's BAMs"


def test_a_converted_batch_sees_its_batch_and_its_results_reach_the_bam():
    # Control: the batch's job still gets its batch, and (non-Ray, one process) what it writes
    # reaches the first BAM's own context, as before.
    a, b = (_ctx(i, module=ws) for i in range(2))
    batcher = ws.SampleJobBatcher()
    mgr = ws.WorkflowManager(enable_batching=False)
    (job,) = mgr._convert_batched_jobs([batcher._create_batched_job(_simple_jobs("cnv", [a, b]), "S", "cnv")])
    assert [x.filepath for x in job.context.metadata["_batched_job"].contexts] == [a.filepath, b.filepath]
    job.context.add_result("cnv", {"status": "success"})
    assert a.results.get("cnv") == {"status": "success"}


_PER_BAM_REMOTES = {
    "preprocessing": "preprocessing_handler_remote",
    "bed_conversion": "bed_conversion_handler_remote",
    "mgmt": "mgmt_handler_remote",
    "cnv": "cnv_handler_remote",
    "target": "target_handler_remote",
    "fusion": "fusion_handler_remote",
}


def test_a_per_bam_step_is_not_run_again_when_its_worker_dies():
    # Finding 2: Ray retries a task whose worker died (default 3), which would add its data twice.
    got = {t: getattr(wr, name)._default_options.get("max_retries") for t, name in _PER_BAM_REMOTES.items()}
    assert got == {t: 0 for t in _PER_BAM_REMOTES}


def test_a_classifier_keeps_rays_retries():
    # Control: the per-sample classifiers are left as they were (they add nothing up).
    assert wr.sturgeon_handler_remote._default_options.get("max_retries") is None


def _broken(text):
    def handler(job):
        raise RuntimeError(text)

    return handler


def _rows(coord):
    return [(r["job_type"], r["filepath"], r.get("failure_kind")) for r in coord.cut_short]


def test_a_step_that_failed_is_reported_at_each_start_and_not_run(tmp_path):
    # Finding 5: an ended failure must not disappear at the next start.
    bam = _bam(tmp_path / "run")
    first = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    _work_through(first[0], steps=["bed_conversion", "mgmt", "target", "fusion"])
    with pytest.raises(RuntimeError):
        wr._wrap_real_handler(_broken("synthetic CNV failure"), "cnv")(_step_job(first[0], "cnv"))
    for _ in range(2):  # every start
        coord = FakeCoord()
        assert _scan(coord, tmp_path / "run", tmp_path / "work") == [], "a failed step was run again"
        assert _rows(coord) == [("cnv", bam, "failed_before")]
        assert "synthetic CNV failure" in coord.cut_short[0]["error"]


def test_a_preprocessing_that_failed_is_reported_at_the_next_start_and_not_run(tmp_path):
    bam = _bam(tmp_path / "run")
    first = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    with pytest.raises(RuntimeError):
        _work_through(first[0], steps=[], preprocessing=_broken("synthetic preprocessing failure"))
    coord = FakeCoord()
    assert _scan(coord, tmp_path / "run", tmp_path / "work") == []
    assert _rows(coord) == [("preprocessing", bam, "failed_before")]


def test_a_bam_that_failed_inside_a_batch_that_went_on_is_reported_at_the_next_start(tmp_path):
    a = _bam(tmp_path / "run", "S_pass_0.bam", 0)
    b = _bam(tmp_path / "run", "S_pass_1.bam", 1)
    first = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    by_path = {j.context.filepath: j for j in first}
    for j in first:
        _work_through(j, steps=["bed_conversion", "mgmt", "target", "fusion"])
    q = wr.job_queue_of("cnv")
    contexts = [by_path[a].context, by_path[b].context]
    batch = wr.BatchedJob(next(_ids), "cnv", q, [f"{q}:cnv"], 0, contexts, "b1", "S")
    job = wr.Job(batch.job_id, "cnv", q, [f"{q}:cnv"], 0, contexts[0])
    job.context.metadata["_batched_job"] = batch

    def one_file_fails(job):  # the batch went on; one of its BAMs failed (batch 7's file failures)
        job.context.add_metadata("file_failures", [{"file": b, "error": "synthetic: unreadable BAM", "step": "bam"}])

    wr._wrap_real_handler(one_file_fails, "cnv")(job)
    coord = FakeCoord()
    assert _scan(coord, tmp_path / "run", tmp_path / "work") == []
    assert _rows(coord) == [("cnv", b, "failed_before")]
    assert "synthetic: unreadable BAM" in coord.cut_short[0]["error"]


def test_an_expected_failure_of_a_fail_only_bam_is_not_reported_at_the_next_start(tmp_path):
    # Control: a fail-only BAM's expected failure (batch 7) is not a failure.
    _bam(tmp_path / "run", "S_fail_0.bam", 0)
    first = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    assert first[0].context.metadata.get("fail_only_bam_submission")
    _work_through(first[0], steps=["mgmt", "cnv", "target", "fusion"])
    wr._wrap_real_handler(_broken("synthetic: fail reads not supported"), "bed_conversion")(
        _step_job(first[0], "bed_conversion")
    )
    coord = FakeCoord()
    assert _scan(coord, tmp_path / "run", tmp_path / "work") == [] and coord.cut_short == []


def test_a_bam_whose_steps_all_ended_well_reports_nothing_at_the_next_start(tmp_path):
    # Control.
    _bam(tmp_path / "run")
    first = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    _work_through(first[0])
    coord = FakeCoord()
    assert _scan(coord, tmp_path / "run", tmp_path / "work") == [] and coord.cut_short == []


def test_a_cut_short_step_whose_report_was_lost_is_reported_at_the_next_start(tmp_path):
    # Finding 4: a second stop after the intake decided, before the report.
    from robin.file_ledger import BamIntake, record_step_start

    bam = _bam(tmp_path / "run")
    first = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    _work_through(first[0], steps=["bed_conversion", "mgmt", "target", "fusion"])
    record_step_start(_step_job(first[0], "cnv"), "cnv")  # session 1 stopped while CNV ran
    BamIntake(str(tmp_path / "work"), PLAN).decide(bam)  # session 2 decided, then stopped before reporting
    coord = FakeCoord()
    assert _scan(coord, tmp_path / "run", tmp_path / "work") == []
    assert _rows(coord) == [("cnv", bam, "cut_short")]


def test_the_coordinator_keeps_the_kind_of_an_earlier_failure(tmp_path):
    c = CoordinatorClass(target_panel="rCNS2", preset="standard", enable_batching=False)
    c.set_work_dir(str(tmp_path))
    c.record_cut_short_steps([{
        "sample_id": "S", "filepath": "/data/S_pass_0.bam", "job_type": "cnv",
        "error": "cnv failed for this BAM in an earlier session", "failure_kind": "failed_before",
    }])
    assert [(r["job_type"], r["failure_kind"]) for r in c._failed_jobs_log] == [("cnv", "failed_before")]


# Round 3: the second external review's findings ---------------------------------------------


def test_a_steps_own_failed_file_counts_when_it_replaces_the_list(tmp_path):
    # Finding 2: the BED step replaces the list of failed files; an earlier step's entry was in it.
    from robin.file_ledger import FileLedger

    work = str(tmp_path)
    ctx = wr.WorkflowContext(filepath="/data/S_pass_1.bam", metadata={
        "bam_key": "k1", "work_dir": work, "bam_metadata": {"sample_id": "S"},
        "file_failures": [{"file": "/data/S_pass_0.bam", "error": "synthetic: an earlier step's failure", "step": "master_csv"}],
    })

    def replaces(job):  # as bed_conversion does: its own failures, in place of the list
        job.context.add_metadata("file_failures", [
            {"file": "/data/S_pass_1.bam", "error": "synthetic merge failure", "step": "merge"},
        ])

    q = wr.job_queue_of("bed_conversion")
    wr._wrap_real_handler(replaces, "bed_conversion")(wr.Job(next(_ids), "bed_conversion", q, [f"{q}:bed_conversion"], 0, ctx))
    state = FileLedger(work).load()
    assert state["k1"].steps["bed_conversion"].file_failed == "synthetic merge failure", "the step's own failed file was lost"


def test_a_non_ray_step_is_not_failed_by_its_siblings_failed_file(tmp_path):
    # Finding 2: with batching off, sibling steps run at once on one shared context.
    from robin.file_ledger import FileLedger

    work = str(tmp_path)
    ctx = ws.WorkflowContext(filepath="/data/S_pass_0.bam", metadata={
        "bam_key": "k0", "work_dir": work, "bam_metadata": {"sample_id": "S"},
    })
    bed_done = threading.Event()

    def cnv(job):  # starts first, ends after BED has reported its failed file
        if not bed_done.wait(5):
            raise RuntimeError("synthetic: BED never ran")

    def bed(job):
        job.context.add_metadata("file_failures", [
            {"file": "/data/S_pass_0.bam", "error": "synthetic BED-only merge failure", "step": "merge"},
        ])
        bed_done.set()

    mgr = ws.WorkflowManager(enable_batching=False)
    mgr.enqueue_jobs = lambda jobs: None
    queues, threads = [], []
    for step, handler in [("cnv", cnv), ("bed_conversion", bed)]:
        q = queue.Queue()
        q.put(ws.Job(next(_ids), step, ctx, "fast", [f"fast:{step}"]))
        t = threading.Thread(target=mgr.worker, args=(step, q, {step: handler}), daemon=True)
        t.start()
        queues.append(q)
        threads.append(t)
    for q in queues:
        q.join()
    mgr.running = False
    for t in threads:
        t.join(timeout=5)
    steps = FileLedger(work).load()["k0"].steps
    assert steps["bed_conversion"].file_failed == "synthetic BED-only merge failure"
    assert steps["cnv"].file_failed is None, "a step that worked was marked failed by its sibling's failed file"


def test_the_coordinator_records_a_repeated_row_once(tmp_path):
    # Finding 3: a row sent again (its first answer lost) is not a second failure.
    c = CoordinatorClass(target_panel="rCNS2", preset="standard", enable_batching=False)
    c.set_work_dir(str(tmp_path))
    row = {
        "sample_id": "unknown", "filepath": "/data/S_pass_0.bam", "job_type": "preprocessing",
        "error": "S_pass_0.bam has no BAM end-of-file block", "failure_kind": "incomplete",
    }
    c.record_cut_short_steps([row])
    c.record_cut_short_steps([dict(row)])
    assert len(c._failed_jobs_log) == 1 and c.file_failed_by_type == {"preprocessing": 1}


# The non-Ray runner ------------------------------------------------------------------------


class _Manager:
    def __init__(self):
        self.jobs = []

    def enqueue_jobs(self, jobs):
        self.jobs.extend(jobs)


def _simple_watcher(folder, work, manager):
    def build(path):
        jobs = ws.default_file_classifier(path, PLAN, "rCNS2")
        for j in jobs:
            j.context.add_metadata("work_dir", str(work))
        return jobs

    return ws.FileWatcher(str(folder), build, manager, "rCNS2", patterns=["*.bam"], show_progress=False)


def _simple_steps(job):
    """The non-Ray worker loop runs the BAM's preprocessing and its steps."""
    mgr = ws.WorkflowManager(enable_batching=False)
    mgr.enqueue_jobs = lambda jobs: None

    def prep(j):
        j.context.add_result("preprocessing", {"status": "success", "sample_id": "S"})

    for step, handler in [("preprocessing", prep)] + [(s, _noop) for s in STEPS]:
        q = queue.Queue()
        q.put(ws.Job(next(_ids), step, job.context, "fast", [f"fast:{step}"]))
        t = threading.Thread(target=mgr.worker, args=("t", q, {step: handler}), daemon=True)
        t.start()
        q.join()
        mgr.running = False
        t.join(timeout=5)
        mgr.running = True


def test_the_non_ray_watcher_does_not_take_a_finished_bam_in_again(tmp_path):
    bam = _bam(tmp_path / "run")
    first = _Manager()
    _simple_watcher(tmp_path / "run", tmp_path / "work", first).handle_file(bam)
    assert _types(first.jobs) == ["preprocessing"]
    _simple_steps(first.jobs[0])
    again = _Manager()
    _simple_watcher(tmp_path / "run", tmp_path / "work", again).handle_file(bam)
    assert again.jobs == [], "a finished BAM was taken in again"


def test_a_non_ray_preprocessing_starts_every_step_of_the_plan():
    # Control: with no list, every step of the plan starts, as before.
    ctx = ws.WorkflowContext(filepath="/data/S_pass_0.bam", metadata={"original_workflow": PLAN})
    job = ws.Job(1, "preprocessing", ctx, "fast", PLAN)
    assert sorted(_types(job.can_trigger_jobs())) == sorted(STEPS)


def test_a_non_ray_preprocessing_told_its_steps_starts_only_those():
    ctx = ws.WorkflowContext(filepath="/data/S_pass_0.bam", metadata={"original_workflow": PLAN, "only_steps": ["cnv"]})
    job = ws.Job(1, "preprocessing", ctx, "fast", PLAN)
    assert _types(job.can_trigger_jobs()) == ["cnv"]


# Batch 10e part 1, round 2: a failure found by another job is recorded for its own BAM ----------


def _ledger_lines(work, *records):
    from robin.file_ledger import FileLedger

    FileLedger(str(work)).append(list(records))


def test_a_failure_line_survives_a_later_end_of_the_same_step(tmp_path):
    # The race: another job records the BAM's failure before the BAM's own step has ended.
    from robin.file_ledger import FileLedger

    _ledger_lines(
        tmp_path,
        {"e": "take", "key": "k1", "file": "/data/A.bam"},
        {"e": "start", "key": "k1", "step": "target", "sample": "S", "file": "/data/A.bam"},
        {"e": "file_failed", "key": "k1", "step": "target", "error": "staged coverage unreadable"},
        {"e": "end", "key": "k1", "step": "target", "status": "success", "outcome": "success"},
    )
    step = FileLedger(str(tmp_path)).load()["k1"].steps["target"]
    assert step.ended and step.file_failed == "staged coverage unreadable"


def test_a_step_with_only_a_failure_line_is_reported_at_the_next_start(tmp_path):
    from robin.bam_identity import bam_key
    from robin.file_ledger import BamIntake

    bam = _bam(tmp_path / "run")
    key = bam_key(str(bam))
    _ledger_lines(
        tmp_path / "work",
        {"e": "take", "key": key, "file": str(bam)},
        {"e": "file_failed", "key": key, "step": "target", "error": "staged coverage unreadable"},
    )
    rows = BamIntake(str(tmp_path / "work"), PLAN).decide(str(bam)).cut_short
    assert [(r["job_type"], r["failure_kind"]) for r in rows] == [("target", "failed_before")]
    assert "staged coverage unreadable" in rows[0]["error"]
