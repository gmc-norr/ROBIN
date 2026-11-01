"""Batch 8: the panel's summary depth counts each base once.

"Targets Estimated Coverage" (and the enrichment and quality label made from it) is samtools
bedcov -j over the panel with its overlapping rows merged, divided by the merged length, on all
of a sample's reads. Each row keeps its own depth. A folder, or a staged BAM, from before this
change has no merged-panel counts: its number stays the row sum (the old one). The merged panels
are written by hand in tests/target_fixtures.py, not by the code under test. Synthetic data only.
"""

import os
import types
from pathlib import Path

import pandas as pd
import pysam
import pytest

from target_fixtures import (
    MERGED_PANEL, MERGED_PANEL_OVERLAPS, PANEL, PANEL_OVERLAPS, SAMPLE,
    samtools_merged_depth, write_bams, write_panel,
)
from robin.analysis.target_analysis import TargetAnalysis

import robin

SRC = Path(robin.__file__).resolve().parent
UNION_FILE = "target_union_coverage.csv"


def _analysis(tmp_path, monkeypatch, rows):
    panel = write_panel(str(tmp_path / "panel.bed"), rows)
    monkeypatch.setattr(TargetAnalysis, "_find_target_bed", lambda self, p: panel)
    ta = TargetAnalysis(work_dir=str(tmp_path / "work"), target_panel="TESTPANEL",
                        batch_size=1, use_staging=True)
    return ta, panel


def _add(ta, bam, drop_union_part=False):
    ta.process_file_with_staging(bam, {"sample_id": SAMPLE})
    if drop_union_part:  # a BAM staged before this change: no merged-panel counts
        staging = Path(ta._get_staging_dir(SAMPLE))
        for p in staging.glob("union_*.parquet"):
            p.unlink()
    ta.accumulate_staged_files(SAMPLE, force=True)


def _row_sum(sample_dir):
    bed = pd.read_csv(sample_dir / "bed_coverage_main.csv")
    return float(bed["bases"].sum() / (bed["endpos"] - bed["startpos"]).sum())


def _shown(sample_dir):
    """The report's number and the GUI summary card's text."""
    from robin.gui.components.summary import _extract_coverage_data
    from robin.reporting.sections.coverage import CoverageSection

    section = CoverageSection.__new__(CoverageSection)
    section.report = types.SimpleNamespace(output=str(sample_dir))
    section._initialize_data()
    return section.coverage_data["target_coverage"], _extract_coverage_data(Path(sample_dir))["target_coverage"]


def test_every_overlap_shape_counts_each_base_once(tmp_path, monkeypatch):
    # U1 (every overlap shape) and U2 (each row unchanged).
    from robin.analysis.target_coverage import panel_depth

    ta, panel = _analysis(tmp_path, monkeypatch, PANEL_OVERLAPS)
    bams, merged = write_bams(str(tmp_path), 2)
    for bam in bams:
        _add(ta, bam)
    sample_dir = tmp_path / "work" / SAMPLE
    want = samtools_merged_depth(merged, MERGED_PANEL_OVERLAPS, str(tmp_path))
    assert abs(want - _row_sum(sample_dir)) > 0.01  # the fixture tells the two apart
    report, summary = _shown(sample_dir)
    assert report == pytest.approx(want, rel=1e-9)
    assert summary == f"{want:.2f}x"
    bed = pd.read_csv(sample_dir / "bed_coverage_main.csv")
    assert panel_depth(str(sample_dir), bed) == pytest.approx(want, rel=1e-9)
    # U2: each row's own bases are samtools bedcov -j for that row
    ref = [int(line.split("\t")[-1]) for line in pysam.bedcov("-j", panel, merged).splitlines()]
    assert list(bed["bases"]) == ref


def test_the_gui_coverage_panel_uses_the_shared_helper():
    src = (SRC / "gui" / "components" / "coverage.py").read_text()
    assert "panel_depth(" in src
    assert "target_depth(state" not in src


def test_a_bam_staged_before_this_change_falls_back_to_the_row_sum(tmp_path, monkeypatch):
    # U3 (guard): the first BAM has no merged-panel part; the second has one.
    ta, _ = _analysis(tmp_path, monkeypatch, PANEL)
    bams, _ = write_bams(str(tmp_path), 2)
    _add(ta, bams[0], drop_union_part=True)
    _add(ta, bams[1])
    sample_dir = tmp_path / "work" / SAMPLE
    assert not (sample_dir / UNION_FILE).exists()
    report, summary = _shown(sample_dir)
    assert report == pytest.approx(_row_sum(sample_dir), rel=1e-9)
    assert summary == f"{_row_sum(sample_dir):.2f}x"


def test_a_bam_staged_before_this_change_after_a_new_one_removes_the_union_file(tmp_path, monkeypatch):
    # U3b (guard, added after the lock): a new BAM writes the merged-panel file; a later batch
    # without the part makes it unknown again, so the file goes and the row sum shows.
    ta, _ = _analysis(tmp_path, monkeypatch, PANEL)
    bams, _ = write_bams(str(tmp_path), 2)
    _add(ta, bams[0])
    sample_dir = tmp_path / "work" / SAMPLE
    _add(ta, bams[1], drop_union_part=True)
    assert not (sample_dir / UNION_FILE).exists()
    report, _ = _shown(sample_dir)
    assert report == pytest.approx(_row_sum(sample_dir), rel=1e-9)


def test_a_folder_written_before_this_change_falls_back_to_the_row_sum(tmp_path, monkeypatch):
    # U4 (guard): totals exist without the merged-panel file; then a new BAM arrives.
    ta, _ = _analysis(tmp_path, monkeypatch, PANEL)
    bams, _ = write_bams(str(tmp_path), 2)
    _add(ta, bams[0])
    sample_dir = tmp_path / "work" / SAMPLE
    if (sample_dir / UNION_FILE).exists():
        (sample_dir / UNION_FILE).unlink()  # as a folder written before this change
    _add(ta, bams[1])
    assert not (sample_dir / UNION_FILE).exists()
    report, _ = _shown(sample_dir)
    assert report == pytest.approx(_row_sum(sample_dir), rel=1e-9)


def test_a_panel_without_overlaps_shows_the_row_sum(tmp_path, monkeypatch):
    # K1 (control): no overlaps, so merging changes nothing.
    rows = [r for r in PANEL if r[3] != "GENE_B"]
    ta, _ = _analysis(tmp_path, monkeypatch, rows)
    bams, merged = write_bams(str(tmp_path), 2)
    for bam in bams:
        _add(ta, bam)
    sample_dir = tmp_path / "work" / SAMPLE
    want = samtools_merged_depth(merged, [(c, s, e) for c, s, e, _ in rows], str(tmp_path))
    assert want == pytest.approx(_row_sum(sample_dir), rel=1e-12)
    report, _ = _shown(sample_dir)
    assert report == pytest.approx(want, rel=1e-9)


def test_an_old_folder_with_no_new_bams_shows_the_row_sum(tmp_path, monkeypatch):
    # K2 (control): an old folder (no merged-panel file) read as it is.
    ta, _ = _analysis(tmp_path, monkeypatch, PANEL)
    bams, _ = write_bams(str(tmp_path), 2)
    for bam in bams:
        _add(ta, bam)
    sample_dir = tmp_path / "work" / SAMPLE
    if (sample_dir / UNION_FILE).exists():
        (sample_dir / UNION_FILE).unlink()
    report, summary = _shown(sample_dir)
    assert report == pytest.approx(_row_sum(sample_dir), rel=1e-9)
    assert summary == f"{_row_sum(sample_dir):.2f}x"


# ------------------------------------------------------------------ the external review's findings
def _union_write_fails(monkeypatch):
    """Writing the union file raises (a disk full, say); every other file is written."""
    from robin.analysis import target_analysis as ta_mod

    real = ta_mod.atomic_open

    def failing(path, *args, **kwargs):
        if str(path).endswith(UNION_FILE):
            raise OSError("synthetic disk full")
        return real(path, *args, **kwargs)

    monkeypatch.setattr(ta_mod, "atomic_open", failing)
    return lambda: monkeypatch.setattr(ta_mod, "atomic_open", real)


def test_a_commit_that_fails_part_way_never_shows_a_partial_total(tmp_path, monkeypatch):
    # R1: the second of four commits writes bed_coverage_main.csv, then fails on the union file.
    ta, _ = _analysis(tmp_path, monkeypatch, PANEL)
    bams, _ = write_bams(str(tmp_path), 4)
    _add(ta, bams[0])
    restore = _union_write_fails(monkeypatch)
    _add(ta, bams[1])
    restore()
    sample_dir = tmp_path / "work" / SAMPLE
    # read right after the failed commit: the union file is a generation behind the rows
    report, _ = _shown(sample_dir)
    assert report == pytest.approx(_row_sum(sample_dir), rel=1e-9), report
    _add(ta, bams[2])
    _add(ta, bams[3])
    report, summary = _shown(sample_dir)
    assert report == pytest.approx(_row_sum(sample_dir), rel=1e-9), report
    assert summary == f"{_row_sum(sample_dir):.2f}x"
    assert not (sample_dir / UNION_FILE).exists()


def test_a_first_commit_that_fails_part_way_falls_back_to_the_row_sum(tmp_path, monkeypatch):
    # R1 guard: the first commit fails on the union file; the sample stays on the row sum.
    ta, _ = _analysis(tmp_path, monkeypatch, PANEL)
    bams, _ = write_bams(str(tmp_path), 2)
    restore = _union_write_fails(monkeypatch)
    _add(ta, bams[0])
    restore()
    _add(ta, bams[1])
    sample_dir = tmp_path / "work" / SAMPLE
    report, _ = _shown(sample_dir)
    assert report == pytest.approx(_row_sum(sample_dir), rel=1e-9)
    assert not (sample_dir / UNION_FILE).exists()


@pytest.mark.parametrize("when", ["between_commits", "within_one_batch"])
def test_a_panel_changed_mid_run_is_not_added_to_the_old_one(tmp_path, monkeypatch, when):
    # R2: BAM 1 with PANEL, then the same BED file holds PANEL_OVERLAPS for BAM 2.
    ta, panel = _analysis(tmp_path, monkeypatch, PANEL)
    bams, _ = write_bams(str(tmp_path), 2)
    if when == "between_commits":
        _add(ta, bams[0])
        write_panel(panel, PANEL_OVERLAPS)
        _add(ta, bams[1])
    else:  # both BAMs staged, then committed in one batch
        ta.process_file_with_staging(bams[0], {"sample_id": SAMPLE})
        write_panel(panel, PANEL_OVERLAPS)
        ta.process_file_with_staging(bams[1], {"sample_id": SAMPLE})
        ta.accumulate_staged_files(SAMPLE, force=True)
    sample_dir = tmp_path / "work" / SAMPLE
    assert not (sample_dir / UNION_FILE).exists()
    report, _ = _shown(sample_dir)
    assert report == pytest.approx(_row_sum(sample_dir), rel=1e-9), report


def test_an_unreadable_union_part_keeps_the_bams_other_data(tmp_path, monkeypatch):
    # R3: two BAMs staged; the first one's union part is damaged.
    ta, panel = _analysis(tmp_path, monkeypatch, PANEL)
    bams, merged = write_bams(str(tmp_path), 2)
    for bam in bams:
        ta.process_file_with_staging(bam, {"sample_id": SAMPLE})
    staging = Path(ta._get_staging_dir(SAMPLE))
    sorted(staging.glob("union_*.parquet"))[0].write_bytes(b"PAR1 damaged PAR1")
    result = ta.accumulate_staged_files(SAMPLE, force=True)
    assert result["files_processed"] == 2, result
    sample_dir = tmp_path / "work" / SAMPLE
    assert not (sample_dir / UNION_FILE).exists()
    bed = pd.read_csv(sample_dir / "bed_coverage_main.csv")
    ref = [int(line.split("\t")[-1]) for line in pysam.bedcov("-j", panel, merged).splitlines()]
    assert list(bed["bases"]) == ref  # both BAMs counted
    report, _ = _shown(sample_dir)
    assert report == pytest.approx(_row_sum(sample_dir), rel=1e-9)
