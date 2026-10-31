"""The LGA panel's "Hide Results" must offer the file that was shown.

``_hide_lga_results`` rebuilt the summary view itself, with "View Analysis
Results" pointing at a hard-coded ``<sample>/<sample>/lightweight_gene_analysis_
results.json``. Since batch 2b the panel can show the legacy file instead (the
first one ``pick_lga_results`` trusts), so after "Hide" the button could point at
a file that does not exist, or at one the trust check had turned down.

Hide now re-runs ``_check_lga_results``, which picks the trusted file again and
draws the same view through ``_display_lga_results``. ``pick_lga_results`` is
tested in ``tests/test_lga_sample_inputs.py``.

Lint-shaped: GUI code has no harness here (see ``test_gui_timer_cleanup.py``),
so this reads the source. It pins the wiring, not what a browser shows.
"""

import ast
from pathlib import Path

COVERAGE = Path(__file__).resolve().parents[1] / "src/robin/gui/components/coverage.py"


def _function(name):
    tree = ast.parse(COVERAGE.read_text())
    found = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == name]
    assert len(found) == 1, f"expected one {name}, found {len(found)}"
    return found[0]


def test_hide_re_runs_the_results_check():
    hide = _function("_hide_lga_results")
    calls = {
        n.func.id
        for n in ast.walk(hide)
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name)
    }
    assert "_check_lga_results" in calls


def test_results_from_either_location_are_not_called_legacy():
    """Both result locations are current (a run given the monitored root writes
    the flat one), so neither may be captioned "(legacy location)" (batch 3c
    final review, Important 1)."""
    check = _function("_check_lga_results")
    strings = [
        n.value for n in ast.walk(check) if isinstance(n, ast.Constant) and isinstance(n.value, str)
    ]
    assert not [s for s in strings if "legacy" in s.lower()]


def test_hide_does_not_name_a_results_file_itself():
    hide = _function("_hide_lga_results")
    strings = [
        n.value for n in ast.walk(hide) if isinstance(n, ast.Constant) and isinstance(n.value, str)
    ]
    assert not [s for s in strings if "lightweight_gene_analysis_results" in s]
