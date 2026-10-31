"""A CNV run with no breakpoints must replace the last run's breakpoints.

Since CNV-5 dropped the fake end-of-chromosome breakpoint, a flat signal gives
an empty breakpoint list in a normal, successful run. ``save_cnv_files`` wrote
``cnv_data_array.npy`` and ``breakpoints_NNN.bed`` only when the list was not
empty, so after a run that found a breakpoint, a flat run left the old array
for the GUI and the old BED as the newest one for the master BED.

Each test saves two runs of a synthetic 200-bin chr1 signal, as the
2026-10-31 fix-series review's ``cnv_probe`` does.
"""

import logging

import numpy as np

from robin.analysis.cnv_analysis import detect_breakpoints_from_cnv, save_cnv_files
from robin.analysis.master_bed_generator import _build_master_bed_data, _load_bed_file

BIN = 1000
LOG = logging.getLogger("test_cnv_breakpoint_files")
SHIFT = {"chr1": np.r_[np.zeros(100), np.ones(100) * 5]}
FLAT = {"chr1": np.zeros(200)}


def _save(sample, counter, signal):
    breakpoints = detect_breakpoints_from_cnv(signal, BIN, LOG)
    save_cnv_files(str(sample), counter, signal, signal, signal, BIN, 0.1, breakpoints, "Unknown", {}, LOG)
    return breakpoints


def _two_runs(tmp_path, first, second):
    sample = tmp_path / "root" / "A"
    sample.mkdir(parents=True)
    assert len(_save(sample, 1, first)) == (1 if first is SHIFT else 0)
    assert len(_save(sample, 2, second)) == (1 if second is SHIFT else 0)
    return sample


def test_a_flat_run_clears_the_saved_breakpoint_array(tmp_path):
    sample = _two_runs(tmp_path, SHIFT, FLAT)
    arr = np.load(sample / "cnv_data_array.npy", allow_pickle=True)
    assert len(arr) == 0
    assert arr.dtype.names == ("name", "start", "end")


def test_a_flat_run_writes_an_empty_breakpoint_bed_for_its_counter(tmp_path):
    sample = _two_runs(tmp_path, SHIFT, FLAT)
    bed = sample / "bed_files" / "breakpoints_002.bed"
    assert bed.exists()
    assert bed.read_text() == ""


def test_the_master_bed_for_a_flat_run_has_no_old_breakpoints(tmp_path):
    sample = _two_runs(tmp_path, SHIFT, FLAT)
    built = _build_master_bed_data("A", str(sample.parent), 2, None, str(tmp_path / "missing.fa"), LOG)
    regions = set() if built is None else set(zip(built[0]["start"], built[0]["end"]))
    assert (99000, 101000) not in regions


def test_breakpoints_after_a_flat_run_are_shown(tmp_path):
    sample = _two_runs(tmp_path, FLAT, SHIFT)
    arr = np.load(sample / "cnv_data_array.npy", allow_pickle=True)
    assert [(str(r["name"]), int(r["start"]), int(r["end"])) for r in arr] == [("chr1", 99000, 101000)]


def test_an_empty_bed_loads_as_empty_without_a_warning(tmp_path, caplog):
    bed = tmp_path / "breakpoints_002.bed"
    bed.write_text("")
    with caplog.at_level(logging.WARNING, logger="robin.analysis.master_bed"):
        df = _load_bed_file(str(bed))
    assert df.empty
    assert "Error loading BED file" not in caplog.text
