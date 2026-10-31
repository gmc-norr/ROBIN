"""A results file is replaced whole or not at all.

``export_json`` and ``save_cnv_files`` opened their outputs with "w" and wrote
in place, so a reader (the GUI's refresh, a second run) could see a half
written file, and a write that failed part-way destroyed the previous good
one. They now write a temporary file in the same folder and rename it over
the old one. These tests make a write fail part-way and check that the
previous file is still intact. Synthetic data only.
"""

import json
import logging

import numpy as np
import pytest

from lga_fixtures import write_clinvar
from robin.analysis.cnv_analysis import save_cnv_files
from robin.analysis.lightweight_gene_analysis import LightweightGeneAnalysis

LOG = logging.getLogger("test_atomic_writes")


def test_a_failed_results_write_leaves_the_previous_file(tmp_path, monkeypatch):
    vcf = write_clinvar(tmp_path / "clinvar.vcf", [(11, "C", "T")])
    analysis = LightweightGeneAnalysis(str(tmp_path / "work"), str(vcf))
    out = tmp_path / "lightweight_gene_analysis_results.json"
    analysis.export_json({}, str(out))
    before = out.read_text()

    def half_then_fail(obj, f, **kwargs):
        f.write('{"metadata": {"trunc')
        raise OSError("disk full (synthetic)")

    monkeypatch.setattr(json, "dump", half_then_fail)
    with pytest.raises(OSError):
        analysis.export_json({}, str(out))
    assert out.read_text() == before
    assert [p.name for p in tmp_path.iterdir() if p.name.endswith(".tmp")] == []


def test_a_failed_cnv_write_leaves_the_previous_arrays(tmp_path, monkeypatch):
    signal = {"chr1": np.r_[np.zeros(100), np.ones(100) * 5]}
    save_cnv_files(str(tmp_path), 1, signal, signal, signal, 1000, 0.1, [], "Unknown", {}, LOG)
    before = np.load(tmp_path / "CNV.npy", allow_pickle=True).item()

    import numpy.lib.format as fmt

    def half_then_fail(fp, array, *args, **kwargs):
        fp.write(b"\x93NUMPY trunc")
        raise OSError("disk full (synthetic)")

    monkeypatch.setattr(fmt, "write_array", half_then_fail)
    save_cnv_files(str(tmp_path), 2, signal, signal, signal, 1000, 0.1, [], "Unknown", {}, LOG)
    monkeypatch.undo()
    after = np.load(tmp_path / "CNV.npy", allow_pickle=True).item()
    assert after.keys() == before.keys()
    assert [p.name for p in tmp_path.iterdir() if p.name.endswith(".tmp")] == []
