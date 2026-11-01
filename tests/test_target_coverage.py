"""Per-BAM totals that add up, covered stretches, and the accumulation's edge cases
(plan 7278fb9, Tasks 3-4). Synthetic BAMs only.
"""

import math
import os
import shutil
from collections import Counter

import pandas as pd
import pysam
import pytest

from target_fixtures import (
    CONTIGS, PANEL, SAMPLE, bam_records, samtools_reference, write_bams, write_panel,
)
from robin.analysis.target_analysis import TargetAnalysis, finalize_accumulation_for_sample

PANEL_NAME = "TESTPANEL"


@pytest.fixture
def panel(tmp_path, monkeypatch):
    path = write_panel(str(tmp_path / "panel.bed"))
    monkeypatch.setattr(TargetAnalysis, "_find_target_bed", lambda self, p: path)
    return path


def _analysis(work):
    return TargetAnalysis(work_dir=str(work), target_panel=PANEL_NAME, batch_size=1, use_staging=True)


def test_contig_totals_match_samtools(tmp_path, panel):
    from robin.analysis.target_coverage import contig_totals

    _, merged = write_bams(str(tmp_path), 1)
    ref = samtools_reference(merged, panel)
    got = contig_totals(merged).set_index("#rname")
    for name, length in CONTIGS:
        numreads, _, depth = ref["contigs"][name]
        assert int(got.loc[name, "numreads"]) == numreads
        assert int(got.loc[name, "bases"]) == round(depth * length)


def test_target_bases_match_bedcov_j_one_row_per_target(tmp_path, panel):
    from robin.analysis.target_coverage import target_bases

    _, merged = write_bams(str(tmp_path), 1)
    ref = samtools_reference(merged, panel)
    rows = PANEL + [PANEL[0], ("chr9", 10, 20, "NOT_IN_BAM")]  # a duplicate row and a missing contig
    got = target_bases(merged, write_panel(str(tmp_path / "wide.bed"), rows))
    keyed = {(r.chrom, r.startpos, r.endpos, r.name): r.bases for r in got.itertuples(index=False)}
    assert len(got) == len(PANEL) + 1
    for key, (bases, _, _) in ref["targets"].items():
        assert keyed[key] == bases, key
    assert keyed[("chr9", 10, 20, "NOT_IN_BAM")] == 0


def test_covered_stretches_add_up_to_samtools_covered_bases(tmp_path, panel):
    from robin.analysis.target_coverage import covered_bases, covered_stretches, union_stretches

    bams, merged = write_bams(str(tmp_path), 4)
    ref = samtools_reference(merged, panel)
    union = union_stretches([covered_stretches(b) for b in bams])
    got = covered_bases(union)
    assert {name: got.get(name, 0) for name, _ in CONTIGS} == {n: v[1] for n, v in ref["contigs"].items()}


def test_a_half_written_staging_set_is_not_claimed(tmp_path, panel):
    """A set whose coverage file isn't there yet stays staged; nothing is paired with it."""
    bams, merged = write_bams(str(tmp_path), 1)
    work = tmp_path / "work"
    ta = _analysis(work)
    ta.process_file_with_staging(bams[0], {"sample_id": SAMPLE})  # counter 0, complete
    staging = work / SAMPLE / "_staging"
    # Another worker is half-way through staging counter 1: no coverage file yet.
    bed = pd.read_parquet(staging / "bedcov_000000.parquet")
    bed["bases"] = bed["bases"] * 1000
    bed.to_parquet(staging / "bedcov_000001.parquet", index=False)
    (staging / "timestamp_000001.txt").write_text("1")
    (staging / "source_bam_000001.txt").write_text(bams[0])
    # Rename counter 0's files to 2 so the half-written set sorts first.
    for f in os.listdir(staging):
        if "_000000." in f:
            os.rename(staging / f, staging / f.replace("_000000.", "_000002."))
    ta.accumulate_staged_files(SAMPLE, force=True)
    ref = samtools_reference(merged, panel)
    out = pd.read_csv(work / SAMPLE / "bed_coverage_main.csv")
    got = {r.name: int(r.bases) for r in out.itertuples(index=False)}
    assert got == {key[3]: v[0] for key, v in ref["targets"].items()}
    assert (staging / "bedcov_000001.parquet").exists()


def test_an_old_coverage_folder_carries_on(tmp_path, panel):
    """A sample folder written by the old code: depth carries on, covered bases become unknown."""
    bams, merged = write_bams(str(tmp_path), 2)
    work = tmp_path / "work"
    sd = work / SAMPLE
    sd.mkdir(parents=True)
    old = pd.read_csv(pd.io.common.StringIO(pysam.coverage(bams[0])), sep="\t")
    old[["#rname", "startpos", "endpos", "numreads", "covbases", "meandepth"]].to_csv(sd / "coverage_main.csv", index=False)
    bed = pd.read_csv(pd.io.common.StringIO(pysam.bedcov(panel, bams[0])), sep="\t",
                      names=["chrom", "startpos", "endpos", "name", "bases"])
    bed["length"] = bed["endpos"] - bed["startpos"] + 1
    bed["coverage"] = bed["bases"] / bed["length"]
    bed[["chrom", "startpos", "endpos", "name", "length", "coverage", "bases"]].to_csv(sd / "bed_coverage_main.csv", index=False)

    ta = _analysis(work)
    ta.process_file_with_staging(bams[1], {"sample_id": SAMPLE})
    ta.accumulate_staged_files(SAMPLE, force=True)
    ref = samtools_reference(merged, panel)
    cov = pd.read_csv(sd / "coverage_main.csv").set_index("#rname")
    for name, _ in CONTIGS:
        assert abs(cov.loc[name, "meandepth"] - ref["contigs"][name][2]) <= 1e-4
        assert math.isnan(cov.loc[name, "covbases"]), "covered bases can't be known for an old folder"


def test_a_bam_before_the_threshold_still_reaches_target_bam(tmp_path, panel):
    """TGT-1: every BAM's panel reads are kept when it arrives, whatever the coverage so far."""
    bams, merged = write_bams(str(tmp_path), 2)
    work = tmp_path / "work"
    ta = _analysis(work)
    ta.callthreshold = 1e9  # no target ever passes
    for bam in bams:
        ta.process_file_with_staging(bam, {"sample_id": SAMPLE})
        ta.accumulate_staged_files(SAMPLE, force=True)
    finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
    target_bam = work / SAMPLE / "target.bam"
    assert target_bam.exists()
    assert bam_records(str(target_bam)) == samtools_reference(merged, panel)["target_bam"]


def test_a_batch_bam_left_by_a_crash_is_not_merged_twice(tmp_path, panel):
    """A guard for Task 5's design (target.bam is an input of every merge). On db7269c's
    code it fails only on which records target.bam holds (measured); its crash case
    can't fail there, since that code deletes every leftover batch BAM (read, not run)."""
    bams, merged = write_bams(str(tmp_path), 2)
    work = tmp_path / "work"
    ta = _analysis(work)
    for bam in bams:
        ta.process_file_with_staging(bam, {"sample_id": SAMPLE})
        ta.accumulate_staged_files(SAMPLE, force=True)
    sd = work / SAMPLE
    kept = tmp_path / "kept"
    kept.mkdir()
    for f in os.listdir(sd):
        if f.startswith("batch_"):
            shutil.copy2(sd / f, kept / f)
    finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
    for f in os.listdir(kept):  # the crash: the merged batch BAMs were never deleted
        shutil.copy2(kept / f, sd / f)
    finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
    assert bam_records(str(sd / "target.bam")) == samtools_reference(merged, panel)["target_bam"]
    assert not [f for f in os.listdir(sd) if f.startswith("batch_")]
