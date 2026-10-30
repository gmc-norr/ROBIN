"""GUI-1: the classification panel must keep refreshing when a classifier is disabled.

The first five tests cover ``classification_files_to_refresh`` directly. They can
only fail by ``ImportError`` on the old code, so they cannot tell a half fix from a
whole one. ``test_panel_keeps_refreshing_with_a_classifier_disabled`` is the test
that reddens for GUI-1's own reason: it drives ``add_classification_section`` with
``pannanodx``/``random_forest`` absent from ``workflow_steps`` and watches the panel
either advance or freeze.

That end-to-end test runs in a subprocess because the flow needs NiceGUI stubbed,
and stubbing NiceGUI in-process would poison the real-NiceGUI imports other tests
in the same pytest session rely on.
"""

import os
import subprocess
import sys
from pathlib import Path

from robin.gui.components.classification import classification_files_to_refresh

_REPO_ROOT = Path(__file__).resolve().parents[1]
_STUBS = _REPO_ROOT / "reviews" / "2026-10-16-repro" / "gui" / "stubs2"


def test_only_charts_tools_are_considered(tmp_path):
    """A tool absent from charts must not be looked up at all (GUI-1: KeyError)."""
    (tmp_path / "sturgeon_scores.csv").write_text("x\n")
    charts = {"Sturgeon": {"file": "sturgeon_scores.csv", "last_mtime": None}}
    mtimes, changed, fresh = classification_files_to_refresh(charts, tmp_path)
    assert set(mtimes) == {"Sturgeon"}
    assert set(changed) == {"Sturgeon"}
    assert changed["Sturgeon"] is True
    assert fresh is True


def test_unchanged_file_is_not_flagged(tmp_path):
    p = tmp_path / "sturgeon_scores.csv"
    p.write_text("x\n")
    charts = {"Sturgeon": {"file": "sturgeon_scores.csv", "last_mtime": p.stat().st_mtime}}
    mtimes, changed, fresh = classification_files_to_refresh(charts, tmp_path)
    assert changed["Sturgeon"] is False
    assert fresh is False


def test_missing_file_reports_zero_and_no_change(tmp_path):
    charts = {"NanoDX": {"file": "NanoDX_scores.csv", "last_mtime": 1.0}}
    mtimes, changed, fresh = classification_files_to_refresh(charts, tmp_path)
    assert mtimes["NanoDX"] == 0.0
    assert changed["NanoDX"] is False


def test_no_enabled_classifiers_is_not_a_fresh_visit(tmp_path):
    """Review Focus 5: charts empty. Must not crash and must not claim a fresh visit."""
    assert classification_files_to_refresh({}, tmp_path) == ({}, {}, False)


def test_missing_sample_dir_is_handled(tmp_path):
    charts = {"Sturgeon": {"file": "sturgeon_scores.csv", "last_mtime": None}}
    mtimes, changed, fresh = classification_files_to_refresh(charts, None)
    assert mtimes == {"Sturgeon": 0.0}
    assert changed == {"Sturgeon": False}
    assert fresh is True


_DRIVER_SCRIPT = '''\
"""Drive add_classification_section with pannanodx/random_forest disabled.

Prints one "tick N: <label>" line per refresh tick. argv[1] is the sample dir.
"""
import logging
import os
import pathlib
import sys
import time
import types

logging.basicConfig(level=logging.ERROR, format="%(levelname)s %(message)s")

theme = types.ModuleType("robin.gui.theme")
theme.register_theme_sync_callback = lambda cb, **k: (lambda: None)
theme.styled_table = lambda **k: (None, None)
sys.modules["robin.gui.theme"] = theme

import nicegui

from robin.gui.components import classification as c

d = pathlib.Path(sys.argv[1])


def write(rows):
    with open(d / "sturgeon_scores.csv", "w") as f:
        f.write("timestamp,Glioblastoma,Medulloblastoma,number_probes\\n")
        for r in rows:
            f.write(",".join(map(str, r)) + "\\n")


write([(1000, 0.40, 0.60, 5000)])
(d / "NanoDX_scores.csv").write_text("timestamp,GBM,MB,number_probes\\n1000,0.3,0.7,4000\\n")


class Launcher:
    # pannanodx / random_forest were not run, so charts has no entry for them.
    workflow_steps = ["sturgeon", "nanodx", "cnv"]


c.add_classification_section(d, Launcher())

refresh = [t for t in nicegui.TIMERS if t.interval == 30.0][0].cb
charts = [
    x
    for x in (cell.cell_contents for cell in refresh.__closure__)
    if isinstance(x, dict) and "Sturgeon" in x
][0]
label = charts["Sturgeon"]["summary"]["class"]

refresh()
print("tick 1: %s" % label.text)

time.sleep(0.05)
# A new batch lands: Glioblastoma is now the top class.
write([(1000, 0.40, 0.60, 5000), (2000, 0.97, 0.03, 9000)])
os.utime(d / "sturgeon_scores.csv", None)

refresh()
print("tick 2: %s" % label.text)
refresh()
print("tick 3: %s" % label.text)
'''


def _stubbed_env():
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(_STUBS), str(_REPO_ROOT / "src"), env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    return env


def test_panel_keeps_refreshing_with_a_classifier_disabled(tmp_path):
    """GUI-1 end to end: workflow_steps omits pannanodx/random_forest, so charts
    lacks them; before the fix the second tick raised KeyError and the label froze."""
    script = tmp_path / "drive.py"
    script.write_text(_DRIVER_SCRIPT)
    r = subprocess.run(
        [sys.executable, str(script), str(tmp_path)],
        capture_output=True,
        text=True,
        env=_stubbed_env(),
    )
    assert r.returncode == 0, r.stderr
    assert "KeyError" not in r.stderr, r.stderr
    ticks = [
        line.split(": ", 1)[1] for line in r.stdout.splitlines() if line.startswith("tick")
    ]
    assert len(ticks) == 3, r.stdout
    assert ticks[0] != ticks[1], f"panel froze: {ticks}"
    assert "Glioblastoma" in ticks[1]
