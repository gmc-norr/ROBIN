"""A step that fails or never runs must be counted as failed, on both schedulers (batch 7, ERR-1).

Every case goes through a real local Ray ``Coordinator`` (see ``failure_fixtures.py``), once on
``standard`` (pool actors) and once on ``high`` (one processor per job type). Only the job
handlers are synthetic; each is wrapped exactly like the real ones.

Rows from the locked plan:
  F1  the handler FILES an error under the name the real handler uses  -> counted failed
  F2  the handler RAISES                                               -> counted failed
  C1  the handler succeeds                                             -> counted completed
  C2  a classifier / bed_conversion handler files an error (or raises) on a fail-only BAM
      submission                                                       -> not failed; result
                                                                          is expected_failure
  C3  preprocessing marks skip_downstream                              -> counted skipped

Each scheduler's whole matrix runs once, in a module-scoped fixture, then the Coordinator is shut
down (which is what writes ``failed_jobs.csv``) and the tests assert on what was recorded.
"""

import os
import time

import pytest

from robin import workflow_ray as wr
from failure_fixtures import (
    FAIL_ONLY_TYPES, JOB_TYPES, SCHEDULERS, Harness, init_ray, shutdown_ray,
)


@pytest.fixture(scope="module")
def ray_cluster():
    init_ray()
    yield
    shutdown_ray()


class Campaign:
    def __init__(self):
        self.outcomes = {}  # (case, job_type) -> Outcome
        self.csv_rows = []  # failed_jobs.csv rows after shutdown

    def rows_for(self, outcome):
        return [r for r in self.csv_rows if r["sample_id"] == outcome.sample_id]


@pytest.fixture(scope="module", params=SCHEDULERS)
def campaign(request, ray_cluster, tmp_path_factory):
    scheduler = request.param
    work = str(tmp_path_factory.mktemp(f"failures_{scheduler}"))
    harness = Harness(scheduler, work)
    c = Campaign()
    c.scheduler = scheduler
    try:
        for jt in JOB_TYPES:
            c.outcomes[("F1", jt)] = harness.run_one(jt, "file_error", tag="F1")
            c.outcomes[("F2", jt)] = harness.run_one(jt, "raise", tag="F2")
            c.outcomes[("C1", jt)] = harness.run_one(jt, "ok", tag="C1")
        for jt in FAIL_ONLY_TYPES:
            c.outcomes[("C2", jt)] = harness.run_one(jt, "file_error", tag="C2", fail_only=True)
            c.outcomes[("C2raise", jt)] = harness.run_one(jt, "raise", tag="C2", fail_only=True)
        c.outcomes[("C3", "preprocessing")] = harness.run_one("preprocessing", "skip", tag="C3")
        # F5: two jobs in a row (same worker process) each send a worker warning.
        c.outcomes[("F5a", "mgmt")] = harness.run_one("mgmt", "warn", tag="F5a")
        c.outcomes[("F5b", "mgmt")] = harness.run_one("mgmt", "warn", tag="F5b")
        c.stats = harness.stats()
    finally:
        c.csv_rows = harness.close()
    return c


def _counted_failed(c, o, text):
    assert o.d_failed == 1, f"{c.scheduler}/{o.job_type}: failed count moved by {o.d_failed}, expected 1"
    assert o.d_completed == 0, f"{c.scheduler}/{o.job_type}: also counted completed"
    assert o.d_failed_by_type == {o.job_type: 1}, (
        f"{c.scheduler}/{o.job_type}: failed_by_type moved by {o.d_failed_by_type}"
    )
    rows = c.rows_for(o)
    assert len(rows) == 1, f"{c.scheduler}/{o.job_type}: {len(rows)} failed_jobs.csv rows, expected 1"
    assert rows[0]["job_type"] == o.job_type
    assert text in rows[0]["error"], f"csv error text {rows[0]['error']!r} lacks {text!r}"


@pytest.mark.parametrize("job_type", JOB_TYPES)
def test_F1_filed_error_counts_as_failed(campaign, job_type):
    o = campaign.outcomes[("F1", job_type)]
    _counted_failed(campaign, o, f"synthetic {job_type} failure")


@pytest.mark.parametrize("job_type", JOB_TYPES)
def test_F2_raised_exception_counts_as_failed(campaign, job_type):
    o = campaign.outcomes[("F2", job_type)]
    _counted_failed(campaign, o, f"synthetic {job_type} crash")


@pytest.mark.parametrize("job_type", JOB_TYPES)
def test_C1_success_counts_as_completed(campaign, job_type):
    o = campaign.outcomes[("C1", job_type)]
    assert o.d_completed == 1, f"{campaign.scheduler}/{job_type}: completed moved by {o.d_completed}"
    assert o.d_failed == 0 and o.d_failed_by_type == {}
    assert o.d_completed_by_type == {job_type: 1}
    assert campaign.rows_for(o) == []


@pytest.mark.parametrize("job_type", FAIL_ONLY_TYPES)
def test_C2_filed_error_on_fail_only_bam_is_not_failed(campaign, job_type):
    o = campaign.outcomes[("C2", job_type)]
    assert o.d_failed == 0 and o.d_failed_by_type == {}, (
        f"{campaign.scheduler}/{job_type}: a fail-only BAM's filed error counted as failed"
    )
    assert o.d_completed == 1
    assert campaign.rows_for(o) == []


@pytest.mark.parametrize("job_type", FAIL_ONLY_TYPES)
def test_C2_filed_error_on_fail_only_bam_result_is_expected_failure(campaign, job_type):
    o = campaign.outcomes[("C2", job_type)]
    assert o.result.get("status") == "expected_failure", (
        f"{campaign.scheduler}/{job_type}: result is {o.result!r}"
    )
    assert o.result.get("reason") == "fail_only_bam_submission"


@pytest.mark.parametrize("job_type", FAIL_ONLY_TYPES)
def test_C2_raised_exception_on_fail_only_bam_is_expected_failure(campaign, job_type):
    """The behaviour the filed-error case must match: it works today for an exception."""
    o = campaign.outcomes[("C2raise", job_type)]
    assert o.d_failed == 0 and o.d_completed == 1
    assert o.result.get("status") == "expected_failure"
    assert campaign.rows_for(o) == []


def test_C3_deliberate_skip_is_skipped_not_failed(campaign):
    o = campaign.outcomes[("C3", "preprocessing")]
    assert o.d_skipped == 1, f"{campaign.scheduler}: skipped moved by {o.d_skipped}"
    assert o.d_failed == 0 and o.d_completed == 0 and o.d_failed_by_type == {}
    assert campaign.rows_for(o) == []


# ---------------------------------------------------------------------------------------------
# Task 2: one rule for a job's outcome (a job's own filed errors decide it, whatever their name)
# ---------------------------------------------------------------------------------------------


@pytest.fixture(scope="module", params=SCHEDULERS)
def rule_campaign(request, ray_cluster, tmp_path_factory):
    scheduler = request.param
    work = str(tmp_path_factory.mktemp(f"rule_{scheduler}"))
    harness = Harness(scheduler, work)
    c = Campaign()
    c.scheduler = scheduler
    try:
        # Review Focus 1: an error an earlier job left in the same context, own handler succeeds.
        c.outcomes[("stale", "mgmt")] = harness.run_one("mgmt", "ok", tag="stale", pre_error=True)
        # ... and the job's own error still fails it when an old one is there too.
        c.outcomes[("stale_fail", "mgmt")] = harness.run_one(
            "mgmt", "file_error", tag="stale", pre_error=True
        )
        # Review Focus 3: a successful target triggers igv_bam (control); a failed one must not.
        c.outcomes[("trig_ok", "target")] = harness.run_one(
            "target", "ok", tag="trig", downstream=True
        )
        c.igv_after_ok = harness.finished_count("igv_bam", timeout=30)
        c.outcomes[("trig_fail", "target")] = harness.run_one(
            "target", "file_error", tag="trig", downstream=True
        )
        time.sleep(3)  # long enough for a wrongly triggered igv_bam to have run
        c.igv_after_fail = harness.finished_count("igv_bam") - c.igv_after_ok
    finally:
        c.csv_rows = harness.close()
    return c


def test_error_from_an_earlier_job_does_not_fail_a_later_one(rule_campaign):
    o = rule_campaign.outcomes[("stale", "mgmt")]
    assert o.d_completed == 1 and o.d_failed == 0 and o.d_failed_by_type == {}, (
        f"{rule_campaign.scheduler}: stale error counted against the later job "
        f"(completed {o.d_completed}, failed {o.d_failed})"
    )
    assert o.result.get("status") == "success"
    assert rule_campaign.rows_for(o) == []


def test_own_error_still_fails_when_an_old_one_is_present(rule_campaign):
    o = rule_campaign.outcomes[("stale_fail", "mgmt")]
    _counted_failed(rule_campaign, o, "synthetic mgmt failure")


def test_filed_error_fallback_result_says_failed(rule_campaign):
    o = rule_campaign.outcomes[("stale_fail", "mgmt")]
    assert o.result == {"status": "failed"}


def test_failed_target_submits_no_igv_bam(rule_campaign):
    # Control: a successful target does trigger igv_bam (so the check below can fail).
    assert rule_campaign.igv_after_ok == 1, (
        f"{rule_campaign.scheduler}: control broken, igv_bam finished "
        f"{rule_campaign.igv_after_ok}x after a successful target"
    )
    o = rule_campaign.outcomes[("trig_fail", "target")]
    assert o.d_failed == 1 and o.d_failed_by_type == {"target": 1}
    assert rule_campaign.igv_after_fail == 0, (
        f"{rule_campaign.scheduler}: igv_bam ran after target filed target_analysis"
    )


def test_cnv_handler_files_its_batch_error(tmp_path):
    """F3: process_multiple_bams returns an error_message -> the handler files cnv_analysis."""
    from types import SimpleNamespace
    from unittest import mock

    from robin.analysis import cnv_analysis
    from robin.workflow_ray import Job, WorkflowContext

    bam = str(tmp_path / "s_pass_0.bam")
    inner = WorkflowContext(filepath=bam, metadata={"bam_metadata": {"sample_id": "s"}})
    batched = SimpleNamespace(
        get_file_count=lambda: 1,
        get_sample_id=lambda: "s",
        batch_id="b1",
        get_filepaths=lambda: [bam],
        contexts=[inner],
    )
    ctx = WorkflowContext(filepath=bam, metadata={"_batched_job": batched})
    job = Job(7, "cnv", "cnv", ["cnv:cnv"], 0, ctx)
    with mock.patch.object(
        cnv_analysis,
        "process_multiple_bams",
        return_value={"error_message": "synthetic batch failure", "files_processed": 0},
    ):
        cnv_analysis.cnv_handler(job, work_dir=str(tmp_path), target_panel="rCNS2")
    assert [(e["job_type"], e["error"]) for e in ctx.errors] == [
        ("cnv_analysis", "synthetic batch failure")
    ]


def _run_simple_worker(
    job_type, handler, metadata=None, filepath="x_pass_0.bam", old_error=False, enqueued=None
):
    """Drive the real WorkflowManager.worker loop on one job; return (manager, context).

    ``enqueued``: a list that receives every job the worker hands on (triggered or next).
    """
    import queue
    import threading

    from robin.workflow_simple import Job, WorkflowContext, WorkflowManager

    mgr = WorkflowManager(enable_batching=False)
    if enqueued is not None:
        mgr.enqueue_jobs = enqueued.extend
    ctx = WorkflowContext(
        filepath=filepath,
        metadata={"bam_metadata": {"sample_id": "simple_s"}, **(metadata or {})},
    )
    if old_error:
        ctx.add_error("earlier_job", "left over from an earlier job")
    q = queue.Queue()
    q.put(Job(1, job_type, ctx, "fast", [f"fast:{job_type}"]))
    t = threading.Thread(target=mgr.worker, args=("t", q, {job_type: handler}), daemon=True)
    t.start()
    q.join()
    mgr.running = False
    t.join(timeout=5)
    return mgr, ctx


def test_non_ray_runner_counts_a_filed_error_as_failed():
    def handler(job):
        job.context.add_error("mgmt_analysis", "synthetic mgmt failure")

    mgr, ctx = _run_simple_worker("mgmt", handler)
    assert mgr.failed_jobs_by_type == {"mgmt": 1}
    assert mgr.completed_jobs_by_type == {}
    assert [e["job_type"] for e in ctx.errors] == ["mgmt_analysis"]  # filed once, not doubled


def test_non_ray_runner_ignores_an_error_from_an_earlier_job():
    mgr, _ = _run_simple_worker("mgmt", lambda job: None, old_error=True)
    assert mgr.completed_jobs_by_type == {"mgmt": 1} and mgr.failed_jobs_by_type == {}


def test_non_ray_runner_counts_a_raised_exception_as_failed():
    def handler(job):
        raise RuntimeError("synthetic crash")

    mgr, _ = _run_simple_worker("mgmt", handler)
    assert mgr.failed_jobs_by_type == {"mgmt": 1}


def test_non_ray_runner_handler_that_files_an_error_and_raises_is_failed_once():
    """M10a: the filed error stays once, the exception is recorded once, the job fails once."""

    def handler(job):
        job.context.add_error("mgmt_analysis", "synthetic mgmt failure")
        raise RuntimeError("synthetic crash")

    mgr, ctx = _run_simple_worker("mgmt", handler)
    assert mgr.failed_jobs_by_type == {"mgmt": 1}
    assert mgr.completed_jobs_by_type == {}
    assert [(e["job_type"], e["error"]) for e in ctx.errors] == [
        ("mgmt_analysis", "synthetic mgmt failure"),
        ("mgmt", "Job mgmt failed: synthetic crash"),
    ]


PREPROCESSING_PLAN = {"original_workflow": ["preprocessing:preprocessing", "mgmt:mgmt", "cnv:cnv"]}


@pytest.mark.parametrize("how", ["files_an_error", "raises"])
def test_non_ray_runner_failed_job_triggers_nothing(how):
    """M10b: a failed preprocessing job hands on none of the jobs it would trigger."""

    def handler(job):
        if how == "raises":
            raise RuntimeError("synthetic crash")
        job.context.add_error("preprocessing", "synthetic preprocessing failure")

    enqueued = []
    mgr, _ = _run_simple_worker("preprocessing", handler, metadata=PREPROCESSING_PLAN, enqueued=enqueued)
    assert mgr.failed_jobs_by_type == {"preprocessing": 1}
    assert enqueued == []


def test_non_ray_runner_control_successful_job_triggers_its_next_jobs():
    enqueued = []
    mgr, _ = _run_simple_worker(
        "preprocessing", lambda job: None, metadata=PREPROCESSING_PLAN, enqueued=enqueued
    )
    assert mgr.completed_jobs_by_type == {"preprocessing": 1}
    assert sorted(j.job_type for j in enqueued) == ["cnv", "mgmt"]


def test_non_ray_runner_fail_only_filed_error_is_expected_not_failed():
    def handler(job):
        job.context.add_error("sturgeon_analysis", "synthetic sturgeon failure")

    mgr, ctx = _run_simple_worker(
        "sturgeon", handler, metadata={"fail_only_bam_submission": True}, filepath="x_fail_0.bam"
    )
    assert mgr.failed_jobs_by_type == {} and mgr.completed_jobs_by_type == {"sturgeon": 1}
    assert ctx.results["sturgeon"]["status"] == "expected_failure"


# --no-use-ray fail-only submissions: the non-Ray FileWatcher marks them with the Ray runner's
# two rules, and the real worker then judges a fail-only BAM's filed classifier error expected.


def _sturgeon_files_an_error(job):
    job.context.add_error("sturgeon_analysis", "synthetic sturgeon failure")


def _non_ray_watch(tmp_path, names, existing):
    """Run names through the non-Ray FileWatcher, then each job through the real worker.

    existing: the files are already there when the watcher starts (process_existing);
    otherwise each arrives as a watch event, in order. Returns {name: (status, context)}.
    """
    import itertools
    import queue
    import threading

    from watchdog.events import FileCreatedEvent

    from robin.workflow_simple import FileWatcher, Job, WorkflowContext, WorkflowManager

    ids = itertools.count(1)
    jobs = []

    def classify(fp):
        ctx = WorkflowContext(filepath=fp, metadata={"bam_metadata": {"sample_id": "s"}})
        return [Job(next(ids), "sturgeon", ctx, "classification", ["classification:sturgeon"])]

    mgr = WorkflowManager(enable_batching=False)
    mgr.enqueue_jobs = jobs.extend
    watcher = FileWatcher(str(tmp_path), classify, mgr, "rCNS2", patterns=["*.bam"],
                          show_progress=False)
    paths = [tmp_path / n for n in names]
    if existing:
        for p in paths:
            p.write_bytes(b"")
        watcher._process_existing_files()
    else:
        for p in paths:
            p.write_bytes(b"")
            watcher.on_created(FileCreatedEvent(str(p)))
    assert len(jobs) == len(names)

    q = queue.Queue()
    for j in jobs:
        q.put(j)
    t = threading.Thread(target=mgr.worker, args=("t", q, {"sturgeon": _sturgeon_files_an_error}),
                         daemon=True)
    t.start()
    q.join()
    mgr.running = False
    t.join(timeout=5)
    out = {}
    for j in jobs:
        name = os.path.basename(j.context.filepath)
        if j.job_id in mgr.failed_jobs:
            out[name] = ("failed", j.context)
        else:
            out[name] = (j.context.results["sturgeon"]["status"], j.context)
    return out


def test_non_ray_existing_fail_only_bams_errors_are_expected(tmp_path):
    out = _non_ray_watch(tmp_path, ["x_fail_0.bam", "x_fail_1.bam"], existing=True)
    assert {n: s for n, (s, _) in out.items()} == {
        "x_fail_0.bam": "expected_failure", "x_fail_1.bam": "expected_failure"}


def test_non_ray_control_existing_fail_bam_beside_a_pass_bam_is_failed(tmp_path):
    out = _non_ray_watch(tmp_path, ["x_fail_0.bam", "x_pass_0.bam"], existing=True)
    assert {n: s for n, (s, _) in out.items()} == {"x_fail_0.bam": "failed", "x_pass_0.bam": "failed"}


def test_non_ray_watched_fail_bam_before_any_pass_bam_is_expected_after_one_failed(tmp_path):
    # The watch-session rule; the "fail" BAM after a "pass" BAM is the control.
    out = _non_ray_watch(tmp_path, ["a_fail_0.bam", "b_pass_0.bam", "c_fail_1.bam"], existing=False)
    assert {n: s for n, (s, _) in out.items()} == {
        "a_fail_0.bam": "expected_failure", "b_pass_0.bam": "failed", "c_fail_1.bam": "failed"}


def test_non_ray_runner_sibling_error_on_a_shared_context_fails_only_its_filer():
    """can_trigger_jobs gives sibling jobs one context object and each runs on its own thread:
    an error one sibling files must not fail another that succeeds."""
    import queue
    import threading

    from robin.workflow_simple import Job, WorkflowContext, WorkflowManager

    mgr = WorkflowManager(enable_batching=False)
    ctx = WorkflowContext(filepath="x_pass_0.bam", metadata={"bam_metadata": {"sample_id": "s"}})
    started = threading.Event()

    def mgmt(job):  # succeeds, just takes a moment
        started.set()
        time.sleep(0.5)

    def target(job):  # files its own error while mgmt is still running
        started.wait(5)
        job.context.add_error("target_analysis", "target broke")

    qa, qb = queue.Queue(), queue.Queue()
    qa.put(Job(1, "mgmt", ctx, "mgmt", ["mgmt:mgmt"]))
    qb.put(Job(2, "target", ctx, "target", ["target:target"]))
    ta = threading.Thread(target=mgr.worker, args=("a", qa, {"mgmt": mgmt}), daemon=True)
    tb = threading.Thread(target=mgr.worker, args=("b", qb, {"target": target}), daemon=True)
    ta.start()
    tb.start()
    qa.join()
    qb.join()
    mgr.running = False
    ta.join(timeout=5)
    tb.join(timeout=5)
    assert mgr.failed_jobs_by_type == {"target": 1}
    assert mgr.completed_jobs_by_type == {"mgmt": 1}


# ---------------------------------------------------------------------------------------------
# Task 3 (ERR-1 layers 3 and 4): failures and worker messages reach the GUI
# ---------------------------------------------------------------------------------------------


def _publish(stats, state=None):
    """Run the real publisher against recorded stats; return (state, [(update_type_value, data)])."""
    state = {"last_seq": 0} if state is None else state
    sent = []
    wr.publish_failures(stats, state, lambda t, d, p=0: sent.append((t.value, d)))
    return state, sent


def _error_updates(sent):
    return [d for t, d in sent if t == "error_update"]


@pytest.mark.parametrize("job_type", JOB_TYPES)
def test_F1_error_update_names_job_type_sample_and_error(campaign, job_type):
    o = campaign.outcomes[("F1", job_type)]
    _, sent = _publish(campaign.stats)
    errs = _error_updates(sent)
    assert len(errs) == 1
    mine = [r for r in errs[0]["recent"] if r["sample_id"] == o.sample_id]
    assert len(mine) == 1, f"{campaign.scheduler}/{job_type}: no ERROR_UPDATE entry for {o.sample_id}"
    assert mine[0]["job_type"] == job_type
    assert f"synthetic {job_type} failure" in mine[0]["error"]
    assert mine[0]["time"]


def test_error_update_counts_by_category(campaign):
    # Each of the 13 types failed twice (filed error, raised exception).
    _, sent = _publish(campaign.stats)
    payload = _error_updates(sent)[0]
    assert payload["preprocessing_errors"] == 2 * 2  # preprocessing, bed_conversion
    assert payload["classification_errors"] == 4 * 2  # sturgeon, nanodx, pannanodx, random_forest
    assert payload["analysis_errors"] == 7 * 2  # the other seven


def test_error_update_names_raised_jobs_too(campaign):
    o = campaign.outcomes[("F2", "target")]
    _, sent = _publish(campaign.stats)
    mine = [r for r in _error_updates(sent)[0]["recent"] if r["sample_id"] == o.sample_id]
    assert len(mine) == 1 and mine[0]["job_type"] == "target"


def test_error_update_is_newest_first(campaign):
    _, sent = _publish(campaign.stats)
    times = [r["time"] for r in _error_updates(sent)[0]["recent"]]
    assert times == sorted(times, reverse=True)


@pytest.mark.parametrize("case", ["C1", "C2", "C2raise", "C3"])
def test_succeeded_expected_and_skipped_jobs_send_no_failure(campaign, case):
    _, sent = _publish(campaign.stats)
    seen = {r["sample_id"] for r in _error_updates(sent)[0]["recent"]}
    clean = [o.sample_id for (c, _), o in campaign.outcomes.items() if c == case]
    assert clean and not (seen & set(clean)), f"{case} jobs appear in the error panel"


def test_a_second_publish_sends_nothing_new(campaign):
    state, first = _publish(campaign.stats)
    assert first
    _, second = _publish(campaign.stats, state)
    assert second == []


@pytest.mark.parametrize("case", ["F5a", "F5b"])
def test_F5_worker_warning_reaches_the_gui_once(campaign, case):
    o = campaign.outcomes[(case, "mgmt")]
    expected = f"alignment warning for {o.sample_id}"
    state, sent = _publish(campaign.stats)
    warns = [d for t, d in sent if t == "warning_notification" and d["message"] == expected]
    assert len(warns) == 1, f"{campaign.scheduler}: {len(warns)} copies of {expected!r}"
    assert warns[0]["sample_id"] == o.sample_id and warns[0]["title"] == "BAM Alignment Warning"
    _, again = _publish(campaign.stats, state)
    assert [d for t, d in again if t == "warning_notification"] == []


def test_F5_two_jobs_in_one_worker_each_send_their_own_message(campaign):
    _, sent = _publish(campaign.stats)
    msgs = [d["message"] for t, d in sent if t == "warning_notification"]
    assert sorted(msgs) == sorted(
        f"alignment warning for {campaign.outcomes[(c, 'mgmt')].sample_id}" for c in ("F5a", "F5b")
    )


def test_F5_too_many_reads_skip_reaches_the_gui_once(campaign):
    o = campaign.outcomes[("C3", "preprocessing")]
    _, sent = _publish(campaign.stats)
    logs = [d for t, d in sent if t == "log_message" and d["message"].startswith("Skipped (too many reads)")]
    assert [d["message"] for d in logs] == [f"Skipped (too many reads): {o.sample_id}_pass_0.bam"]
    assert logs[0]["level"] == "WARNING"


# ---------------------------------------------------------------------------------------------
# JOB-3: the wrapper calls a handler once, with the arguments it takes
# ---------------------------------------------------------------------------------------------


def _wrapped_job(job_type, metadata, job_id=1):
    ctx = wr.WorkflowContext(filepath="a_pass_0.bam", metadata=dict(metadata))
    return wr.Job(job_id, job_type, "analysis", [f"analysis:{job_type}"], 0, ctx)


def test_a_raising_handler_runs_once_and_writes_nothing_outside_work_dir(tmp_path, monkeypatch):
    work = tmp_path / "work"
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)  # where a call without work_dir would write
    calls = []

    def handler(job, work_dir=None):
        calls.append(work_dir)
        with open(os.path.join(work_dir or ".", f"out_{len(calls)}.txt"), "w") as f:
            f.write("x")
        raise RuntimeError(f"boom from call {len(calls)}")

    job = _wrapped_job("mgmt", {"target_panel": "rCNS2", "work_dir": str(work)})
    with pytest.raises(RuntimeError) as raised:
        wr._wrap_real_handler(handler, "mgmt")(job)
    assert calls == [str(work)]
    assert str(raised.value) == "boom from call 1"
    assert [(e["job_type"], e["error"]) for e in job.context.errors] == [("mgmt", "boom from call 1")]
    # the handler wrote only into work_dir (robin.log here is the wrapper's logging setup)
    assert [p.name for p in elsewhere.glob("out_*")] == []
    assert [p.name for p in work.glob("out_*")] == ["out_1.txt"]


def test_control_a_handler_that_takes_only_the_job_still_runs_once(tmp_path):
    calls = []
    job = _wrapped_job("mgmt", {"target_panel": "rCNS2", "work_dir": str(tmp_path / "work")})
    out = wr._wrap_real_handler(lambda j: calls.append(j.job_id), "mgmt")(job)
    assert calls == [1]
    assert out.errors == []


def _wd_ref_tp(seen):
    def handler(job, work_dir=None, reference=None, target_panel=None):
        seen.append(("wd" if work_dir else None, reference, target_panel))
    return handler


def _wd_ref(seen):
    def handler(job, work_dir=None, reference=None):
        seen.append(("wd" if work_dir else None, reference, None))
    return handler


def _tp_only(seen):
    def handler(job, target_panel=None):
        seen.append((None, None, target_panel))
    return handler


def _wd_positional_only(seen):
    def handler(job, work_dir=None, /):  # the keyword call can't bind: called with the job alone
        seen.append(("wd" if work_dir else None, None, None))
    return handler


FULL_META = {"target_panel": "rCNS2", "reference": "/ref.fa"}


@pytest.mark.parametrize(
    "job_type,make,meta,expected",
    [
        ("fusion", _wd_ref_tp, FULL_META, ("wd", "/ref.fa", "rCNS2")),
        ("mgmt", _wd_ref, FULL_META, ("wd", "/ref.fa", None)),
        ("cnv", _tp_only, FULL_META, (None, None, "rCNS2")),
        ("mgmt", _wd_positional_only, FULL_META, (None, None, None)),
        ("mgmt", _wd_ref, {"reference": "/ref.fa"}, (None, None, None)),  # no target_panel
    ],
    ids=["all_three", "work_dir_and_reference", "target_panel_only", "cannot_bind", "no_panel"],
)
def test_control_a_handler_gets_the_arguments_it_takes_once(tmp_path, job_type, make, meta, expected):
    seen = []
    job = _wrapped_job(job_type, {**meta, "work_dir": str(tmp_path / "work")})
    out = wr._wrap_real_handler(make(seen), job_type)(job)
    assert seen == [expected]
    assert out.errors == []


def test_the_ray_watcher_keeps_its_watch_session_rule():
    """The rule moved to robin.bam_submission; the Ray watcher still applies it."""
    w = wr.RayFileWatcher(coord=None, plan=["preprocessing:preprocessing"], target_panel="rCNS2")
    w._flush_interval_s = 1e9
    w._batch_size = 10**6
    for name in ("a_fail_0.bam", "notes.txt", "b_pass_0.bam", "c_fail_1.bam"):
        w._handle(f"/nonexistent/{name}")
    assert [
        (os.path.basename(j.context.filepath),
         bool(j.context.metadata.get("fail_only_bam_submission", False)))
        for j in w._pending_jobs
    ] == [("a_fail_0.bam", True), ("notes.txt", False), ("b_pass_0.bam", False),
          ("c_fail_1.bam", False)]


def test_existing_files_rule_over_a_folder_with_patterns_and_ignores(tmp_path):
    from robin.bam_submission import detect_fail_only_bam_submission

    sub = tmp_path / "run" / "sub"
    sub.mkdir(parents=True)
    for name in ("x_fail_0.bam", "x_fail_1.bam", "skip_pass_0.bam", ".hidden_pass.bam", "y_pass.txt"):
        (sub / name).write_bytes(b"")
    run = str(tmp_path / "run")
    # the pass BAMs are ignored (pattern, dotfile) and the .txt is not a BAM
    assert detect_fail_only_bam_submission([run], ["*.bam", "*.txt"], ["skip_*"], True) is True
    assert detect_fail_only_bam_submission([run], None, None, True) is False  # skip_pass counts
    assert detect_fail_only_bam_submission([run], None, None, False) is False  # nothing at top
    assert detect_fail_only_bam_submission([str(tmp_path / "missing")], None, None, True) is False


def test_a_message_from_a_crashed_job_rides_with_the_next_job_once():
    """Worker messages are kept when the handler raises, and go out with the next normal job."""
    from robin.gui.app import drain_gui_outbox, send_gui_update
    from robin.gui_launcher import UpdateType

    drain_gui_outbox()

    calls = []

    def crashing(job):
        # (the wrapper's legacy fallback calls a failing handler a second time; send only once)
        if not calls:
            send_gui_update(UpdateType.LOG_MESSAGE, {"message": "from the crash"}, 0)
        calls.append(1)
        raise RuntimeError("boom")

    def fine(job):
        pass

    ctx1 = wr.WorkflowContext(filepath="a_pass_0.bam", metadata={"target_panel": "rCNS2"})
    with pytest.raises(RuntimeError):
        wr._wrap_real_handler(crashing, "mgmt")(wr.Job(1, "mgmt", "analysis", ["analysis:mgmt"], 0, ctx1))
    assert "gui_messages" not in ctx1.metadata
    ctx2 = wr.WorkflowContext(filepath="b_pass_0.bam", metadata={"target_panel": "rCNS2"})
    out = wr._wrap_real_handler(fine, "mgmt")(wr.Job(2, "mgmt", "analysis", ["analysis:mgmt"], 0, ctx2))
    assert [m["data"]["message"] for m in out.metadata["gui_messages"]] == ["from the crash"]
    ctx3 = wr.WorkflowContext(filepath="c_pass_0.bam", metadata={"target_panel": "rCNS2"})
    out3 = wr._wrap_real_handler(fine, "mgmt")(wr.Job(3, "mgmt", "analysis", ["analysis:mgmt"], 0, ctx3))
    assert "gui_messages" not in out3.metadata


def test_event_log_is_bounded_and_numbered():
    class FakeJob:
        job_type = "mgmt"
        context = wr.WorkflowContext(filepath="x_pass_0.bam", metadata={"bam_metadata": {"sample_id": "s"}})

    coord_cls = wr.Coordinator.__ray_metadata__.modified_class
    obj = object.__new__(coord_cls)
    obj._events = wr.deque(maxlen=wr.EVENT_LOG_MAX)
    obj._event_seq = 0
    for i in range(wr.EVENT_LOG_MAX + 50):
        obj._record_events(FakeJob, False, FakeJob.context, "e" * 1000, "returned_error", False)
    assert len(obj._events) == wr.EVENT_LOG_MAX
    seqs = [e["seq"] for e in obj._events]
    assert seqs == sorted(seqs) and seqs[-1] == wr.EVENT_LOG_MAX + 50
    assert all(len(e["error"]) <= 300 for e in obj._events)
