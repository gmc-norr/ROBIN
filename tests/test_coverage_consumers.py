"""Locked rule 5 (plan 7278fb9): the report and the GUI show samtools' numbers.

The global number is mean depth (it was covered bases over genome length, shown
with an "x"); a target's length is end - start for a 0-based BED row (it was one
base longer); covered bases are unknown in a folder written before step 6.
Synthetic data only.
"""

import re
import types
from pathlib import Path

import pandas as pd
import pysam
import pytest

from target_fixtures import (
    CONTIGS, MERGED_PANEL, SAMPLE, samtools_merged_depth, samtools_reference, write_bams, write_panel,
)
from robin.analysis.target_analysis import TargetAnalysis

import robin

# The code the tests import (PYTHONPATH decides), not this checkout's src/.
SRC = Path(robin.__file__).resolve().parent


@pytest.fixture
def sample(tmp_path, monkeypatch):
    panel = write_panel(str(tmp_path / "panel.bed"))
    monkeypatch.setattr(TargetAnalysis, "_find_target_bed", lambda self, p: panel)
    bams, merged = write_bams(str(tmp_path), 2)
    work = tmp_path / "work"
    ta = TargetAnalysis(work_dir=str(work), target_panel="TESTPANEL", batch_size=1, use_staging=True)
    for bam in bams:
        ta.process_file_with_staging(bam, {"sample_id": SAMPLE})
        ta.accumulate_staged_files(SAMPLE, force=True)
    ref = samtools_reference(merged, panel)
    # Batch 8: the target depth is samtools over the panel with its overlapping rows merged
    # (each base once); it was the row sum, which counts GENE_A/GENE_B's overlap twice.
    return work / SAMPLE, ref, samtools_merged_depth(merged, MERGED_PANEL, str(tmp_path))


def _report_section(sample_dir):
    from robin.reporting.sections.coverage import CoverageSection

    section = CoverageSection.__new__(CoverageSection)
    section.report = types.SimpleNamespace(output=str(sample_dir))
    section._initialize_data()
    return section


def test_the_report_shows_samtools_depth(sample):
    sample_dir, ref, target_depth = sample
    section = _report_section(sample_dir)
    assert section.coverage_data["global_coverage"] == pytest.approx(ref["global"], rel=1e-9)
    assert section.coverage_data["target_coverage"] == pytest.approx(target_depth, rel=1e-9)
    covered = {d["name"]: d["covered_bases"] for d in section.chromosome_data}
    assert covered == {name: ref["contigs"][name][1] for name, _ in CONTIGS}


def test_the_gui_summary_shows_samtools_depth(sample):
    from robin.gui.components.summary import _extract_coverage_data

    sample_dir, ref, target_depth = sample
    data = _extract_coverage_data(Path(sample_dir))
    assert data["global_coverage"] == f"{ref['global']:.2f}x"
    assert data["target_coverage"] == f"{target_depth:.2f}x"


def test_an_old_folder_shows_depth_and_no_covered_bases(tmp_path):
    """A folder written before step 6: summed per-BAM covered bases can't be trusted."""
    bams, merged = write_bams(str(tmp_path), 1)
    sample_dir = tmp_path / SAMPLE
    sample_dir.mkdir()
    old = pd.read_csv(pd.io.common.StringIO(pysam.coverage(bams[0])), sep="\t")
    old[["#rname", "startpos", "endpos", "numreads", "covbases", "meandepth"]].to_csv(
        sample_dir / "coverage_main.csv", index=False)
    panel = write_panel(str(tmp_path / "p.bed"))
    bed = pd.read_csv(pd.io.common.StringIO(pysam.bedcov(panel, bams[0])), sep="\t",
                      names=["chrom", "startpos", "endpos", "name", "bases"])
    bed["length"] = bed["endpos"] - bed["startpos"] + 1
    bed["coverage"] = bed["bases"] / bed["length"]
    for name in ("bed_coverage_main.csv", "target_coverage.csv"):
        bed[["chrom", "startpos", "endpos", "name", "length", "coverage", "bases"]].to_csv(sample_dir / name, index=False)
    section = _report_section(sample_dir)
    assert section.coverage_data["global_coverage"] == pytest.approx(
        samtools_reference(merged, panel)["global"], rel=1e-4)
    assert all(d["covered_bases"] is None for d in section.chromosome_data)


OLD_FORMULAS = [
    (re.compile(r'\["startpos"\]\s*\+\s*1'), "a target length of end - start + 1"),
    (re.compile(r'\["covbases"\]\.sum\(\)'), "global coverage from covered bases"),
    (re.compile(r'\["covbases"\]\s*/'), "depth from covered bases"),
]


def test_no_screen_or_report_computes_coverage_the_old_way():
    places = list((SRC / "gui").rglob("*.py")) + list((SRC / "reporting").rglob("*.py")) + [SRC / "gui_launcher.py"]
    found = []
    for path in places:
        for n, line in enumerate(path.read_text().splitlines(), 1):
            for pattern, what in OLD_FORMULAS:
                if pattern.search(line):
                    found.append(f"{path.relative_to(SRC)}:{n}: {what}")
    assert not found, "\n".join(found)
