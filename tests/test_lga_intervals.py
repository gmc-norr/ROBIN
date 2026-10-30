"""VAR-10: every BED interval of a gene is analysed, not only the first."""

import pytest

from lga_fixtures import write_bam, write_bed, write_clinvar
from robin.analysis.lightweight_gene_analysis import LightweightGeneAnalysis

# Reads carry the alternate T everywhere, so every ClinVar C>T site is supported.
EXON1_READS = [(0, "T" * 30)] * 10    # covers 1-based 1..30 at 10x
EXON2_READS = [(100, "T" * 30)] * 30  # covers 1-based 101..130 at 30x


def _report(tmp_path, bed_rows, sites, reads):
    bam = write_bam(tmp_path / "s.bam", reads)
    bed = write_bed(tmp_path / "s.bed", bed_rows)
    vcf = write_clinvar(tmp_path / "clinvar.vcf", sites)
    analysis = LightweightGeneAnalysis(str(tmp_path / "work"), str(vcf))
    return analysis.analyze_bed_regions(str(bam), str(bed))["TEST"]


def test_every_interval_with_the_same_name_is_analysed(tmp_path):
    report = _report(
        tmp_path,
        [(0, 30, "TEST"), (100, 130, "TEST")],
        [(11, "C", "T"), (111, "C", "T")],
        EXON1_READS + EXON2_READS,
    )
    assert sorted(v["position"] for v in report.pathogenic_snps_details) == [11, 111]


def test_coverage_covers_every_interval(tmp_path):
    report = _report(
        tmp_path,
        [(0, 30, "TEST"), (100, 130, "TEST")],
        [(11, "C", "T"), (111, "C", "T")],
        EXON1_READS + EXON2_READS,
    )
    # 30 positions at 10x and 30 at 30x
    assert report.mean_coverage == pytest.approx(20.0)


def test_site_in_overlapping_intervals_counts_once(tmp_path):
    report = _report(
        tmp_path,
        [(0, 30, "TEST"), (5, 30, "TEST")],
        [(11, "C", "T")],
        EXON1_READS,
    )
    assert [v["position"] for v in report.pathogenic_snps_details] == [11]


def test_overlapping_intervals_are_piled_up_once(tmp_path):
    # Both intervals hold a site. Positions 1..20 have 20x, 21..40 have 10x.
    report = _report(
        tmp_path,
        [(0, 30, "TEST"), (5, 40, "TEST")],
        [(11, "C", "T"), (35, "C", "T")],
        [(0, "T" * 40)] * 10 + [(0, "T" * 20)] * 10,
    )
    assert sorted(v["position"] for v in report.pathogenic_snps_details) == [11, 35]
    # 20 positions at 20x and 20 at 10x; counting the overlap twice gives ~15.4
    assert report.mean_coverage == pytest.approx(15.0)
