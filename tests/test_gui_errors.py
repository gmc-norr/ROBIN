"""The error panel and the GUI-side plumbing for failures and worker messages (batch 7, ERR-1).

``GUILauncher._update_errors`` is tested on a stand-in object that holds fake labels. The call
from the GUI's publisher to ``publish_failures`` is checked on the source, so reverting the call
site turns a test red.
"""

import ast
import inspect
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
    p = Panel()
    GUILauncher._update_errors(p, {"preprocessing_errors": 4})
    assert p.error_summary_label.text == "No errors detected."
    assert p.error_summary_label.calls == []
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
