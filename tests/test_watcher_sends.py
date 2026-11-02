"""The Ray watcher sends every BAM it takes in, and says so when it can't (batch 10d part 3, JOB-4).

Synthetic BAMs (``target_fixtures``), the real ``RayFileWatcher`` and intake, and a stand-in
Coordinator that records what is sent.

- Waiting jobs are sent on a timer, not only when the next file arrives.
- A send that fails is logged with the BAM's name and tried again; nothing is dropped.
- At shutdown, ROBIN stops watching and sends what waits before it stops the Coordinator; what
  can't be sent is logged by name (the BAM ledger runs it at the next start).
- A half-written BAM is checked again on the timer: taken in once it completes, with no new
  file event; shown once as a failed file when it stays unchanged and incomplete too long.
"""

import ast
import inspect
import logging
import time

import pytest

pytest.importorskip("pysam")

from test_bam_ledger import PLAN, FakeCoord, _bam, _Done, _Remote, _scan, _step_job, _work_through  # noqa: E402
from robin import workflow_ray as wr  # noqa: E402


@pytest.fixture(autouse=True)
def ray_get(monkeypatch):
    monkeypatch.setattr(wr.ray, "get", lambda x, *a, **k: x.value if isinstance(x, _Done) else x)


def _until(check, seconds=5.0):
    end = time.monotonic() + seconds
    while time.monotonic() < end:
        if check():
            return True
        time.sleep(0.02)
    return check()


def _watcher(coord, work, intake=None, tick=0.05):
    w = wr.RayFileWatcher(coord, PLAN, "rCNS2", patterns=["*.bam"], work_dir=str(work), intake=intake)
    w._flush_interval_s = 1e9  # the watcher's own events never send: only the timer, or a stop
    w._tick_s = tick
    return w


def _half(path_whole, path):
    with open(path_whole, "rb") as f:
        data = f.read()
    with open(path, "wb") as f:
        f.write(data[: len(data) - 100])  # no BAM end-of-file block yet
    return data


class _SendFails(FakeCoord):
    """A Coordinator whose first ``fails`` sends raise."""

    def __init__(self, fails):
        super().__init__()
        self.fails = fails

        def submit(jobs):
            if self.fails > 0:
                self.fails -= 1
                raise RuntimeError("synthetic: the Coordinator could not take the jobs")
            self.jobs.extend(jobs)

        self.submit_jobs = _Remote(submit)


class _ReportFails(FakeCoord):
    """A Coordinator whose first ``fails`` failed-file reports raise (it was briefly away)."""

    def __init__(self, fails):
        super().__init__()
        self.fails = fails

        def record(rows):
            if self.fails > 0:
                self.fails -= 1
                raise RuntimeError("synthetic: the Coordinator was briefly unavailable")
            self.cut_short.extend(rows)

        self.record_cut_short_steps = _Remote(record)


def test_the_last_bam_is_sent_without_another_file_event(tmp_path):
    bam = _bam(tmp_path / "run")
    coord = FakeCoord()
    w = _watcher(coord, tmp_path / "work")
    w.start()
    try:
        w._handle(bam)  # the run's last BAM: no event comes after it
        assert _until(lambda: [j.context.filepath for j in coord.jobs] == [bam]), "the last BAM was never sent"
    finally:
        w.stop()


def test_a_failed_send_is_logged_and_tried_again(tmp_path, caplog):
    bam = _bam(tmp_path / "run")
    coord = _SendFails(1)
    w = _watcher(coord, tmp_path / "work")
    w._handle(bam)
    with caplog.at_level(logging.ERROR):
        w.flush_remaining()  # the Coordinator refuses
    assert any("S_pass_0.bam" in r.getMessage() for r in caplog.records if r.levelno >= logging.ERROR), \
        "the failed send was not logged with the BAM's name"
    w.flush_remaining()  # tried again
    assert [j.context.filepath for j in coord.jobs] == [bam], "the BAM's jobs were dropped"


def test_what_cannot_be_sent_at_the_end_is_logged_and_runs_at_the_next_start(tmp_path, caplog):
    bam = _bam(tmp_path / "run")
    w = _watcher(_SendFails(10**6), tmp_path / "work")
    w._handle(bam)
    with caplog.at_level(logging.ERROR):
        w.flush_remaining()
    text = " ".join(r.getMessage() for r in caplog.records if r.levelno >= logging.ERROR)
    assert "S_pass_0.bam" in text and "next start" in text, "the unsent BAM was not named"
    again = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    assert [j.job_type for j in again] == ["preprocessing"], "the ledger did not run it at the next start"


def test_the_ray_run_stops_watching_before_it_stops_the_coordinator():
    tree = ast.parse(inspect.getsource(wr.run))
    shutdown_blocks = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Try):
            blocks = [h.body for h in node.handlers] + [node.finalbody]
            shutdown_blocks += [b for b in blocks if b]
    checked = 0
    for body in shutdown_blocks:
        calls = []  # (line, the call's full name, e.g. "coord.shutdown.remote")
        for stmt in body:
            for n in ast.walk(stmt):
                if isinstance(n, ast.Call):
                    calls.append((n.lineno, ast.unparse(n.func)))
        calls.sort()
        stops = [line for line, name in calls if "shutdown" in name]
        if stops:
            checked += 1
            watching = [line for line, name in calls if name == "_stop_watching"]
            assert watching and min(watching) < min(stops), \
                "the Coordinator is stopped before watching stops and waiting jobs are sent"
    assert checked >= 2, "expected the Ctrl-C path and the final cleanup to stop the Coordinator"


def test_stop_sends_what_waits_and_ends_the_timer(tmp_path):
    bam = _bam(tmp_path / "run")
    coord = FakeCoord()
    w = _watcher(coord, tmp_path / "work", tick=60.0)  # the timer won't send in time
    w.start()
    w._handle(bam)
    w.stop()
    assert [j.context.filepath for j in coord.jobs] == [bam], "stop did not send what waited"
    assert not w._timer.is_alive(), "the timer still runs after stop"


def test_a_bam_the_scan_found_half_written_is_taken_in_once_it_completes(tmp_path):
    whole = _bam(tmp_path / "made")
    (tmp_path / "run").mkdir()
    path = str(tmp_path / "run" / "S_pass_0.bam")
    data = _half(whole, path)
    coord = FakeCoord()
    w = _watcher(coord, tmp_path / "work")
    assert _scan(coord, tmp_path / "run", tmp_path / "work", intake=w.intake) == []
    w.start()
    try:
        with open(path, "wb") as f:
            f.write(data)  # completes; no observer runs here, so no file event
        assert _until(lambda: [j.job_type for j in coord.jobs] == ["preprocessing"]), \
            "a BAM that completed with no new file event was never taken in"
    finally:
        w.stop()


def test_a_bam_still_incomplete_and_unchanged_after_the_set_time_shows_once_as_a_failed_file(tmp_path):
    whole = _bam(tmp_path / "made")
    (tmp_path / "run").mkdir()
    path = str(tmp_path / "run" / "S_pass_0.bam")
    _half(whole, path)
    coord = FakeCoord()
    w = _watcher(coord, tmp_path / "work")
    w._incomplete_after_s = 0.2
    w.start()
    try:
        w._handle(path)
        assert _until(lambda: coord.cut_short), "a BAM that never completed was not shown"
        time.sleep(0.3)  # more ticks
    finally:
        w.stop()
    assert [(r["filepath"], r["failure_kind"]) for r in coord.cut_short] == [(path, "incomplete")], "shown more than once"
    assert "not analysed" in coord.cut_short[0]["error"]
    assert coord.jobs == []


def test_a_bam_still_growing_is_not_shown_as_incomplete(tmp_path):
    # Control: a copy still being written is not a failure, however long it takes.
    whole = _bam(tmp_path / "made")
    (tmp_path / "run").mkdir()
    path = str(tmp_path / "run" / "S_pass_0.bam")
    _half(whole, path)
    coord = FakeCoord()
    w = _watcher(coord, tmp_path / "work")
    w._incomplete_after_s = 0.3
    w.start()
    try:
        w._handle(path)
        for _ in range(10):
            with open(path, "ab") as f:
                f.write(b"\0" * 10)  # still being written (no end-of-file block)
            time.sleep(0.1)
    finally:
        w.stop()
    assert coord.cut_short == [], "a BAM still being written was shown as incomplete"


def test_a_failed_report_of_an_incomplete_bam_is_tried_again(tmp_path):
    # Round 3, finding 3: the report is kept until the Coordinator has it.
    whole = _bam(tmp_path / "made")
    (tmp_path / "run").mkdir()
    path = str(tmp_path / "run" / "S_pass_0.bam")
    _half(whole, path)
    coord = _ReportFails(1)
    w = _watcher(coord, tmp_path / "work")
    w._incomplete_after_s = 0.2
    w.start()
    try:
        w._handle(path)
        assert _until(lambda: coord.cut_short), "the report was lost when the Coordinator was briefly away"
        time.sleep(0.3)
    finally:
        w.stop()
    assert [(r["filepath"], r["failure_kind"]) for r in coord.cut_short] == [(path, "incomplete")]


def test_a_failed_report_of_an_earlier_failure_is_tried_again(tmp_path):
    from robin.file_ledger import record_step_start

    bam = _bam(tmp_path / "run")
    first = _scan(FakeCoord(), tmp_path / "run", tmp_path / "work")
    _work_through(first[0], steps=["bed_conversion", "mgmt", "target", "fusion"])
    record_step_start(_step_job(first[0], "cnv"), "cnv")  # ROBIN stopped while CNV ran
    coord = _ReportFails(1)
    w = _watcher(coord, tmp_path / "work")  # the next session: the watcher sees the BAM
    w.start()
    try:
        w._handle(bam)
        assert _until(lambda: coord.cut_short), "the report was lost when the Coordinator was briefly away"
        time.sleep(0.3)
    finally:
        w.stop()
    assert [(r["job_type"], r["failure_kind"]) for r in coord.cut_short] == [("cnv", "cut_short")]


def test_a_bam_taken_in_is_sent_once_however_many_ticks_pass(tmp_path):
    # Control: the timer sends what waits once.
    bam = _bam(tmp_path / "run")
    coord = FakeCoord()
    w = _watcher(coord, tmp_path / "work")
    w.start()
    try:
        w._handle(bam)
        assert _until(lambda: coord.jobs)
        time.sleep(0.3)
    finally:
        w.stop()
    assert [j.job_type for j in coord.jobs] == ["preprocessing"]
