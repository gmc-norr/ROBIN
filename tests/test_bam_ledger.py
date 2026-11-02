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


def test_a_step_cut_short_by_a_stop_is_not_run_again_and_is_reported_once(tmp_path):
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
    assert _scan(third, tmp_path / "run", tmp_path / "work") == [] and third.cut_short == [], "reported twice"


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
