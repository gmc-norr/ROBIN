"""The LGA status timer must stop itself when the analysis ends.

``check_lga_status`` polls every 0.5 s while the variant analysis runs in a
thread. It returned ``False`` to "stop checking", but NiceGUI's timer only
looks at a callback's return value to see whether to await it, so the timer
never stopped: after completion it re-ran ``_check_lga_results`` every tick,
and that rebuilds the results panel, so an opened results view was replaced
by the summary half a second later (fix-series review, finding 3).

These tests compile the real ``check_lga_status`` out of ``coverage.py`` (it
is nested deep inside ``add_coverage_section``) and run it in NiceGUI's real
``Timer`` loop with stand-in widgets, bounded at ``MAX_TICKS``. They lean on
NiceGUI 3.18's private ``_skip_registration``, ``_should_stop`` and
``_run_in_loop`` (the version is pinned). No browser is involved.
"""

import ast
import asyncio
from pathlib import Path

from nicegui.timer import Timer

COVERAGE = Path(__file__).resolve().parents[1] / "src/robin/gui/components/coverage.py"
MAX_TICKS = 5


class _Widget:
    """Accepts any call, like a NiceGUI element would."""

    def __getattr__(self, name):
        return lambda *args, **kwargs: self


def _compiled_check_lga_status():
    tree = ast.parse(COVERAGE.read_text())
    found = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "check_lga_status"]
    assert len(found) == 1
    return compile(ast.Module(body=found, type_ignores=[]), "check_lga_status", "exec")


def _drive(status):
    """Run check_lga_status in a real Timer loop; return what happened."""
    view = {"value": "nothing"}
    refreshes = []
    errors = []

    def refresh():  # stands in for _check_lga_results: it redraws the summary
        refreshes.append(view["value"])
        view["value"] = "summary"

    ns = {
        "lga_status": {"status": status, "message": "done", "error": "boom"},
        "lga_status_label": _Widget(),
        "lga_analysis_button": _Widget(),
        "_check_lga_results": refresh,
    }
    exec(_compiled_check_lga_status(), ns)
    ticks = []

    class ProbeTimer(Timer):
        def _skip_registration(self):
            return True

        def _should_stop(self):
            return self._is_canceled or len(ticks) >= MAX_TICKS

    def tick():
        ticks.append(1)
        try:
            ns["check_lga_status"]()
        except Exception as e:  # NiceGUI would swallow it; keep it visible
            errors.append(repr(e))
        # The user opens the results once, after the first refresh drew them.
        if len(refreshes) == 1 and view["value"] == "summary":
            view["value"] = "results opened by the user"

    timer = ProbeTimer(0.001, tick)
    ns["status_timer"] = timer
    asyncio.run(timer._run_in_loop())
    return {"ticks": len(ticks), "refreshes": len(refreshes), "view": view["value"],
            "canceled": timer._is_canceled, "errors": errors}


def test_completion_refreshes_once_and_stops_the_timer():
    r = _drive("completed")
    assert r["errors"] == []
    assert (r["ticks"], r["refreshes"], r["canceled"]) == (1, 1, True)


def test_an_opened_results_view_survives_after_completion():
    assert _drive("completed")["view"] == "results opened by the user"


def test_an_error_stops_the_timer():
    r = _drive("error")
    assert r["errors"] == []
    assert (r["ticks"], r["refreshes"], r["canceled"]) == (1, 0, True)


def test_a_running_analysis_keeps_polling():
    r = _drive("running")
    assert r["errors"] == []
    assert (r["ticks"], r["refreshes"], r["canceled"]) == (MAX_TICKS, 0, False)
