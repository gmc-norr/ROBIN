"""The error panel and the GUI-side plumbing for failures and worker messages (batch 7, ERR-1).

``GUILauncher._update_errors`` is tested on a stand-in object that holds fake labels. The call
from the GUI's publisher to ``publish_failures`` is checked on the source, so reverting the call
site turns a test red.
"""

import ast
import inspect
import logging
import queue
import time
from pathlib import Path

from robin import gui_launcher, workflow_ray
from robin.gui import app as gui_app
from robin.gui_launcher import GUILauncher, UpdateType


class Label:
    def __init__(self, text=""):
        self.text = text
        self.calls = []

    def set_text(self, t):
        self.text = t
        self.calls.append(t)


class Panel:
    """Stand-in for the launcher: only the labels _update_errors touches."""

    def __init__(self):
        self.preprocessing_errors = Label("0")
        self.analysis_errors = Label("0")
        self.classification_errors = Label("0")
        self.error_summary_label = Label("No errors detected.")


def _payload(recent):
    return {
        "preprocessing_errors": 1,
        "analysis_errors": 2,
        "classification_errors": 3,
        "recent": recent,
    }


def test_update_errors_sets_the_three_counts():
    p = Panel()
    GUILauncher._update_errors(p, _payload([]))
    assert (p.preprocessing_errors.text, p.analysis_errors.text, p.classification_errors.text) == ("1", "2", "3")


def test_update_errors_replaces_no_errors_with_the_newest_failures():
    p = Panel()
    t = time.mktime((2026, 11, 1, 12, 34, 56, 0, 0, -1))
    recent = [
        {"time": t, "sample_id": "S1", "job_type": "mgmt", "error": "boom"},
        {"time": t - 60, "sample_id": "S2", "job_type": "cnv", "error": "bust"},
    ]
    GUILauncher._update_errors(p, _payload(recent))
    lines = p.error_summary_label.text.splitlines()
    assert lines == ["12:34:56 S1 mgmt: boom", "12:33:56 S2 cnv: bust"]
    assert "No errors detected." not in p.error_summary_label.text


def test_update_errors_shows_at_most_five():
    p = Panel()
    recent = [{"time": 0, "sample_id": f"S{i}", "job_type": "mgmt", "error": "e"} for i in range(8)]
    GUILauncher._update_errors(p, _payload(recent))
    assert len(p.error_summary_label.text.splitlines()) == 5
    assert "S0 " in p.error_summary_label.text and "S5 " not in p.error_summary_label.text


def test_update_errors_without_recent_leaves_the_summary_alone():
    # Changed by the external review's finding 4 (round 5): a count above 0 with no recent
    # failures no longer leaves "No errors detected."; details already shown are kept
    # (test_shown_failure_details_are_kept_when_an_update_has_none).
    p = Panel()
    GUILauncher._update_errors(p, {"preprocessing_errors": 4})
    assert p.error_summary_label.text == "4 failures recorded; the details are no longer available."
    assert p.preprocessing_errors.text == "4"


def test_handle_update_routes_error_update_to_update_errors():
    src = inspect.getsource(GUILauncher._handle_update)
    assert "UpdateType.ERROR_UPDATE" in src and "self._update_errors(update.data)" in src


def test_publish_gui_calls_publish_failures():
    """_publish_gui is nested in run(); if the call goes, nothing reaches the GUI."""
    tree = ast.parse(Path(workflow_ray.__file__).read_text())
    fns = [n for n in ast.walk(tree) if isinstance(n, ast.AsyncFunctionDef) and n.name == "_publish_gui"]
    assert len(fns) == 1
    called = {
        c.func.id
        for c in ast.walk(fns[0])
        if isinstance(c, ast.Call) and isinstance(c.func, ast.Name)
    }
    assert "publish_failures" in called


def test_send_gui_update_without_a_launcher_keeps_the_update_in_the_outbox():
    assert gui_app.get_gui_launcher() is None
    gui_app.drain_gui_outbox()
    gui_app.send_gui_update(UpdateType.WARNING_NOTIFICATION, {"message": "m"}, priority=5)
    assert gui_app.drain_gui_outbox() == [("warning_notification", {"message": "m"}, 5)]
    assert gui_app.drain_gui_outbox() == []  # taken once


def test_the_outbox_is_bounded_and_keeps_the_newest():
    gui_app.drain_gui_outbox()
    for i in range(gui_app.OUTBOX_MAX + 25):
        gui_app.send_gui_update(UpdateType.LOG_MESSAGE, {"i": i}, 0)
    items = gui_app.drain_gui_outbox()
    assert len(items) == gui_app.OUTBOX_MAX
    assert items[-1][1] == {"i": gui_app.OUTBOX_MAX + 24}
    assert items[0][1] == {"i": 25}


def test_send_gui_update_with_a_launcher_sends_and_does_not_keep(monkeypatch):
    got = []

    class L:
        def send_update(self, t, d, p):
            got.append((t, d, p))

    gui_app.drain_gui_outbox()
    monkeypatch.setattr(gui_app, "_current_gui_launcher", L())
    gui_app.send_gui_update(UpdateType.LOG_MESSAGE, {"x": 1}, 2)
    assert got == [(UpdateType.LOG_MESSAGE, {"x": 1}, 2)]
    assert gui_app.drain_gui_outbox() == []


def test_failures_by_category_mapping():
    f = workflow_ray.failures_by_category(
        {"preprocessing": 1, "bed_conversion": 2, "sturgeon": 3, "nanodx": 1, "pannanodx": 1,
         "random_forest": 1, "mgmt": 5, "cnv": 1, "igv_bam": 1}
    )
    assert f == {"preprocessing": 3, "classification": 6, "analysis": 7}


# ---------------------------------------------------------------- publish_failures (M1, M2, M3)


def _failure(seq):
    return {"seq": seq, "kind": "failure", "sample_id": "S", "job_type": "mgmt", "error": "boom",
            "time": float(seq)}


def _msg(seq, message, update_type="warning_notification", priority=5):
    return {"seq": seq, "kind": "gui_message", "update_type": update_type,
            "data": {"message": message}, "priority": priority}


def _stats(events):
    return {"recent_events": events,
            "failures_by_category": {"preprocessing": 0, "classification": 0, "analysis": 1}}


def test_an_unknown_update_type_does_not_block_the_events_behind_it(caplog):
    # P1: one event's update_type is not an UpdateType; the failure behind it still goes out.
    events = [_msg(1, "w1"), _msg(2, "", update_type="not_a_real_type", priority=0), _failure(3)]
    state = {"last_seq": 0}
    sent = []
    with caplog.at_level(logging.WARNING):
        workflow_ray.publish_failures(_stats(events), state, lambda t, d, p=0: sent.append(t.value))
    assert sent == ["warning_notification", "error_update"]
    assert state["last_seq"] == 3
    assert "not_a_real_type" in caplog.text
    again = []
    workflow_ray.publish_failures(_stats(events), state, lambda t, d, p=0: again.append(t.value))
    assert again == []  # nothing is sent twice


def test_a_send_that_raises_does_not_block_the_events_behind_it(caplog):
    # P1b: send itself raises on the second event.
    events = [_msg(1, "w1"), _msg(2, "w2"), _failure(3)]
    state = {"last_seq": 0}
    sent = []

    def flaky(t, d, p=0):
        if d.get("message") == "w2":
            raise RuntimeError("queue closed")
        sent.append(t.value)

    with caplog.at_level(logging.WARNING):
        workflow_ray.publish_failures(_stats(events), state, flaky)
    assert sent == ["warning_notification", "error_update"]
    assert state["last_seq"] == 3
    assert "queue closed" in caplog.text
    sent.clear()
    workflow_ray.publish_failures(_stats(events), state, flaky)
    assert sent == []  # w1 is not re-sent, the failure not announced twice


def _launcher():
    """A real GUILauncher, NiceGUI not started: only the attributes send_update reads."""
    g = object.__new__(GUILauncher)
    g.update_queue = queue.PriorityQueue()
    g._last_update_times = {}
    g._last_queue_size_logged = 0
    g._last_log_time = 0.0
    g._update_seq = 0
    g.total_updates_enqueued = 0
    return g


def _enqueued(g):
    out = []
    while not g.update_queue.empty():
        _, _, update = g.update_queue.get_nowait()
        if update is not None:
            out.append(update)
    return out


def test_error_update_is_enqueued_when_the_launcher_queue_is_full():
    # M2: the launcher drops updates below priority 5 once its queue holds more than 1000.
    g = _launcher()
    for i in range(1001):
        g.update_queue.put((0, -1 - i, None))
    g.send_update(UpdateType.LOG_MESSAGE, {"message": "control"}, 0)
    assert g.update_queue.qsize() == 1001, "the pre-fill is not past the launcher's drop threshold"
    workflow_ray.publish_failures(_stats([_failure(1)]), {"last_seq": 0}, g.send_update)
    assert [u.update_type for u in _enqueued(g)] == [UpdateType.ERROR_UPDATE]


def test_two_forwarded_messages_of_one_type_within_half_a_second_are_both_enqueued():
    # M2: the launcher merges same-type updates below priority 5 that arrive within 0.5 s.
    g = _launcher()
    skipped = [
        _msg(i, f"Skipped (too many reads): {name}", update_type="log_message", priority=0)
        for i, name in ((1, "a_pass_0.bam"), (2, "b_pass_0.bam"))
    ]
    workflow_ray.publish_failures({"recent_events": skipped}, {"last_seq": 0}, g.send_update)
    assert sorted(u.data["message"] for u in _enqueued(g)) == [
        "Skipped (too many reads): a_pass_0.bam",
        "Skipped (too many reads): b_pass_0.bam",
    ]


def test_events_dropped_from_the_log_between_two_polls_are_reported(caplog):
    # M3 / P2: 230 worker messages then one failure, all between two polls; the log keeps 200.
    class FakeJob:
        job_type = "mgmt"
        context = workflow_ray.WorkflowContext(
            filepath="x_pass_0.bam", metadata={"bam_metadata": {"sample_id": "s"}}
        )

    coord = object.__new__(workflow_ray.Coordinator.__ray_metadata__.modified_class)
    coord._events = workflow_ray.deque(maxlen=workflow_ray.EVENT_LOG_MAX)
    coord._event_seq = 0
    for i in range(230):
        ctx = workflow_ray.WorkflowContext(filepath="x_pass_0.bam", metadata={
            "bam_metadata": {"sample_id": f"s{i}"},
            "gui_messages": [{"update_type": "warning_notification",
                              "data": {"message": f"warn {i}"}, "priority": 5}],
        })
        coord._record_events(FakeJob, True, ctx, None, None, False)
    coord._record_events(FakeJob, False, FakeJob.context, "the failure", "returned_error", False)
    stats = {"recent_events": list(coord._events), "failures_by_category": {"analysis": 1}}
    assert len(stats["recent_events"]) == 200 and stats["recent_events"][0]["seq"] == 32

    sent = []
    with caplog.at_level(logging.WARNING):
        workflow_ray.publish_failures(
            stats, {"last_seq": 0}, lambda t, d, p=0: sent.append((t.value, d, p))
        )
    lost = [(d, p) for t, d, p in sent if t == "log_message"]
    assert len(lost) == 1
    assert lost[0][0]["message"].startswith("31 earlier messages were not shown")
    assert lost[0][1] >= 5  # the launcher neither drops nor merges it
    assert "31" in caplog.text
    # what the log still holds goes out as before
    assert len([d for t, d, _ in sent if t == "warning_notification"]) == 199
    assert [t for t, _, _ in sent].count("error_update") == 1


# ------------------------------------------------------------------ second review, I-1
def _coordinator_stats_failure_first(n_messages=230):
    """A real Coordinator event log: one failure, then n_messages worker messages."""

    class FakeJob:
        job_type = "mgmt"
        context = workflow_ray.WorkflowContext(
            filepath="x_pass_0.bam", metadata={"bam_metadata": {"sample_id": "s"}}
        )

    coord = object.__new__(workflow_ray.Coordinator.__ray_metadata__.modified_class)
    coord._events = workflow_ray.deque(maxlen=workflow_ray.EVENT_LOG_MAX)
    coord._event_seq = 0
    coord._record_events(FakeJob, False, FakeJob.context, "the failure", "returned_error", False)
    for i in range(n_messages):
        ctx = workflow_ray.WorkflowContext(filepath="x_pass_0.bam", metadata={
            "bam_metadata": {"sample_id": f"s{i}"},
            "gui_messages": [{"update_type": "warning_notification",
                              "data": {"message": f"warn {i}"}, "priority": 5}],
        })
        coord._record_events(FakeJob, True, ctx, None, None, False)
    return {"recent_events": list(coord._events),
            "failures_by_category": {"preprocessing": 0, "classification": 0, "analysis": 1}}


def test_a_failure_that_left_the_log_still_refreshes_the_failure_counts():
    # Second review I-1: the failure first, then 230 worker messages, all between two polls.
    # The failure leaves the 200-event log; its count must still reach the GUI.
    stats = _coordinator_stats_failure_first()
    kinds = {e.get("kind") for e in stats["recent_events"]}
    assert not kinds & set(workflow_ray._FAILURE_KINDS)  # the failure has left the log
    sent = []
    workflow_ray.publish_failures(
        stats, {"last_seq": 0}, lambda t, d, p=0: sent.append((t.value, d, p))
    )
    updates = [d for t, d, _ in sent if t == "error_update"]
    assert len(updates) == 1, [t for t, _, _ in sent if t != "warning_notification"]
    assert updates[0]["analysis_errors"] == 1


def test_a_poll_with_no_lost_event_and_no_failure_sends_no_error_update():
    # Control for I-1: nothing lost, nothing failed -> no ERROR_UPDATE.
    msgs = [_msg(i, f"note {i}", update_type="log_message", priority=0) for i in (1, 2, 3)]
    sent = []
    workflow_ray.publish_failures(_stats(msgs), {"last_seq": 0}, lambda t, d, p=0: sent.append(t.value))
    assert sent == ["log_message"] * 3


# ------------------------------------------------------------------ external review, finding 4
def test_an_evicted_failure_does_not_leave_no_errors_detected():
    # The I-1 setup through the real _update_errors: the count is 1, and the summary must not
    # still say "No errors detected." though the failure's details have left the log.
    stats = _coordinator_stats_failure_first()
    sent = []
    workflow_ray.publish_failures(stats, {"last_seq": 0}, lambda t, d, p=0: sent.append((t.value, d)))
    update = next(d for t, d in sent if t == "error_update")
    p = Panel()
    GUILauncher._update_errors(p, update)
    assert p.analysis_errors.text == "1"
    assert p.error_summary_label.text == "1 failure recorded; the details are no longer available."


def test_no_failures_keeps_no_errors_detected():
    # Control for finding 4.
    p = Panel()
    GUILauncher._update_errors(p, {"preprocessing_errors": 0, "analysis_errors": 0,
                                   "classification_errors": 0, "recent": []})
    assert p.error_summary_label.text == "No errors detected."


def test_shown_failure_details_are_kept_when_an_update_has_none():
    # Control for finding 4: details already on screen are still true; keep them.
    p = Panel()
    p.error_summary_label = Label("12:00:00 S1 mgmt: boom")
    GUILauncher._update_errors(p, _payload([]))
    assert p.error_summary_label.text == "12:00:00 S1 mgmt: boom"
