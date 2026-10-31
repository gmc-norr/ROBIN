"""Both CNV views must draw breakpoint lines from the array ROBIN saves.

``save_cnv_files`` saves ``cnv_data_array.npy`` with fields ``name, start,
end``. The absolute-value view read ``r["end_pos"]``, which does not exist, so
its dashed breakpoint lines never rendered: the ValueError was swallowed by an
``except Exception: pass``. The difference view drew the midpoint of ``start``
and ``end``. Both now use ``breakpoint_midpoints``.

The call-site checks read the source (GUI code has no harness here).
"""

import ast
import logging
from pathlib import Path

import numpy as np

from robin.analysis.cnv_analysis import detect_breakpoints_from_cnv, save_cnv_files

CNV_GUI = Path(__file__).resolve().parents[1] / "src/robin/gui/components/cnv.py"
LOG = logging.getLogger("test_cnv_breakpoint_lines")


def _saved_array(tmp_path):
    """The array exactly as save_cnv_files writes it: chr1 and chr2, one shift each."""
    signal = {c: np.r_[np.zeros(100), np.ones(100) * 5] for c in ("chr1", "chr2")}
    breakpoints = detect_breakpoints_from_cnv(signal, 1000, LOG)
    save_cnv_files(str(tmp_path), 1, signal, signal, signal, 1000, 0.1, breakpoints, "Unknown", {}, LOG)
    return np.load(tmp_path / "cnv_data_array.npy", allow_pickle=True)


def test_breakpoint_midpoints_reads_the_saved_array(tmp_path):
    from robin.gui.components.cnv import breakpoint_midpoints

    arr = _saved_array(tmp_path)
    assert breakpoint_midpoints(arr, "chr1") == [100000]
    assert breakpoint_midpoints(arr, "chr3") == []


def _function(tree, name):
    found = [n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == name]
    assert len(found) == 1, name
    return found[0]


def _calls(node, name):
    return [
        n for n in ast.walk(node)
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == name
    ]


def test_both_views_use_the_helper():
    tree = ast.parse(CNV_GUI.read_text())
    assert _calls(_function(tree, "_render_cnv_from_state"), "breakpoint_midpoints")
    assert len(_calls(_function(tree, "_refresh_cnv_sync"), "breakpoint_midpoints")) == 2


def test_no_view_reads_end_pos_from_a_breakpoint_row():
    # Cytoband and gene frames legitimately have an "end_pos" column (e.g.
    # _build_cyto_rows); the breakpoint rows ("r" over bp_array) in the two
    # drawing functions never do.
    tree = ast.parse(CNV_GUI.read_text())
    bad = [
        n.lineno
        for fn in ("_render_cnv_from_state", "_refresh_cnv_sync")
        for n in ast.walk(_function(tree, fn))
        if isinstance(n, ast.Subscript)
        and isinstance(n.value, ast.Name) and n.value.id == "r"
        and isinstance(n.slice, ast.Constant) and n.slice.value == "end_pos"
    ]
    assert bad == []
