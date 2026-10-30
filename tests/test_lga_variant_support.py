"""VAR-8: a pathogenic variant counts as found only when reads support it.

"Supported" means 3 or more alternate reads (LGA calls 1-2 reads "low_support").
"""

import json

from lga_fixtures import write_bam, write_bed, write_clinvar
from robin.analysis.lightweight_gene_analysis import LightweightGeneAnalysis

REF_READ = (0, "C" * 50)
ALT_READ = (0, "C" * 10 + "T" + "C" * 39)  # ClinVar POS 11 (1-based) is the T
SNV = [(11, "C", "T")]


def _analyse(tmp_path, reads, sites=SNV):
    bam = write_bam(tmp_path / "s.bam", reads)
    bed = write_bed(tmp_path / "s.bed", [(0, 50, "TEST")])
    vcf = write_clinvar(tmp_path / "clinvar.vcf", sites)
    analysis = LightweightGeneAnalysis(str(tmp_path / "work"), str(vcf))
    reports = analysis.analyze_bed_regions(str(bam), str(bed))
    out = tmp_path / "results.json"
    analysis.export_json(reports, str(out))
    return reports["TEST"], json.loads(out.read_text())


def test_reference_only_site_is_not_reported_as_present(tmp_path):
    report, data = _analyse(tmp_path, [REF_READ] * 20)
    assert report.pathogenic_snps_found == 0
    gene = data["genes"]["TEST"]
    assert gene["summary"]["pathogenic_status"] == "none_detected"
    assert gene["summary"]["total_variants"] == 0
    assert gene["summary"]["high_coverage_variants"] == 0
    assert gene["summary"]["sites_checked"] == 1
    # The checked site is still listed, with its evidence.
    assert gene["variants"][0]["supported"] is False
    assert gene["variants"][0]["evidence_analysis"]["evidence_level"] == "no_support"
    meta = data["metadata"]
    assert meta["total_variants_found"] == 0
    assert meta["genes_with_pathogenic_variants"] == 0
    assert meta["genes_with_good_coverage_variants"] == 0
    assert meta["total_sites_checked"] == 1
    assert meta["genes_with_good_coverage"] == 1


def test_supported_site_is_reported(tmp_path):
    report, data = _analyse(tmp_path, [ALT_READ] * 20)
    assert report.pathogenic_snps_found == 1
    gene = data["genes"]["TEST"]
    assert gene["summary"]["pathogenic_status"] == "present_with_good_coverage"
    assert gene["variants"][0]["supported"] is True
    assert data["metadata"]["total_variants_found"] == 1
    assert data["metadata"]["genes_with_pathogenic_variants"] == 1
    assert data["metadata"]["genes_with_good_coverage_variants"] == 1


def test_supported_site_at_low_coverage(tmp_path):
    _, data = _analyse(tmp_path, [ALT_READ] * 5)
    summary = data["genes"]["TEST"]["summary"]
    assert summary["pathogenic_status"] == "present_but_low_coverage"
    assert summary["low_coverage_variants"] == 1
    assert data["metadata"]["genes_with_good_coverage"] == 0


def test_two_alternate_reads_are_insufficient_evidence(tmp_path):
    report, data = _analyse(tmp_path, [REF_READ] * 18 + [ALT_READ] * 2)
    assert report.pathogenic_snps_found == 0
    summary = data["genes"]["TEST"]["summary"]
    assert summary["pathogenic_status"] == "insufficient_evidence"
    assert summary["total_variants"] == 0
    assert summary["insufficient_evidence_sites"] == 1


def test_three_alternate_reads_count_as_found(tmp_path):
    report, data = _analyse(tmp_path, [REF_READ] * 17 + [ALT_READ] * 3)
    assert report.pathogenic_snps_found == 1
    assert data["genes"]["TEST"]["summary"]["pathogenic_status"] == "present_with_good_coverage"


def test_indel_site_is_insufficient_evidence(tmp_path):
    # LGA can't call indels from a pileup, so it must not say "none detected".
    _, data = _analyse(tmp_path, [REF_READ] * 20, sites=[(11, "CC", "C")])
    summary = data["genes"]["TEST"]["summary"]
    assert summary["pathogenic_status"] == "insufficient_evidence"
    assert summary["total_variants"] == 0


def test_multi_base_substitution_is_insufficient_evidence(tmp_path):
    # CC>TT is not an indel, but a one-base pileup can't call it either.
    # The reads carry TT here, so "none detected" would be wrong.
    mnv_read = (0, "C" * 10 + "TT" + "C" * 38)
    _, data = _analyse(tmp_path, [mnv_read] * 20, sites=[(11, "CC", "TT")])
    summary = data["genes"]["TEST"]["summary"]
    assert summary["pathogenic_status"] == "insufficient_evidence"
    assert data["genes"]["TEST"]["variants"][0]["needs_review"] is True


def test_site_where_reads_show_neither_allele_is_insufficient_evidence(tmp_path):
    # ClinVar C>T, but every read has G: wrong build or position, not "none detected".
    g_read = (0, "C" * 10 + "G" + "C" * 39)
    _, data = _analyse(tmp_path, [g_read] * 20)
    summary = data["genes"]["TEST"]["summary"]
    assert summary["pathogenic_status"] == "insufficient_evidence"


def test_sites_needing_review_are_flagged_per_variant(tmp_path):
    # The GUI lists these sites, so the rule lives in one place.
    _, data = _analyse(tmp_path, [REF_READ] * 18 + [ALT_READ] * 2)
    variant = data["genes"]["TEST"]["variants"][0]
    assert variant["supported"] is False
    assert variant["needs_review"] is True


def test_reference_only_site_needs_no_review(tmp_path):
    _, data = _analyse(tmp_path, [REF_READ] * 20)
    assert data["genes"]["TEST"]["variants"][0]["needs_review"] is False


def _pileup(tmp_path, reads):
    """Pile up ``reads`` over chr1:0-50 and index the entries by 1-based position."""
    bam = write_bam(tmp_path / "s.bam", reads)
    vcf = write_clinvar(tmp_path / "clinvar.vcf", SNV)
    analysis = LightweightGeneAnalysis(str(tmp_path / "work"), str(vcf))
    entries = analysis._perform_targeted_pileup(
        str(bam), {"chromosome": "chr1", "start": 0, "end": 50}
    )
    return {entry["position"]: entry for entry in entries}


def test_pileup_reference_is_the_most_common_base(tmp_path):
    # Every read carries C here, so the pileup's reference must be C, not the
    # letter that happens to sort highest.
    entry = _pileup(tmp_path, [REF_READ] * 20)[11]
    assert entry["reference"] == "C"
    assert entry["coverage"] == 20
    assert entry["variant_alleles"] == []
    assert entry["has_variants"] is False


def test_insertion_site_on_reference_reads_needs_review(tmp_path):
    # C>CA on reads that all show C: the reference base matches, so the site is
    # one LGA can't call, not a position mismatch.
    _, data = _analyse(tmp_path, [REF_READ] * 20, sites=[(11, "C", "CA")])
    variant = data["genes"]["TEST"]["variants"][0]
    assert variant["evidence_analysis"]["evidence_level"] == "insertion_requires_review"
    assert variant["needs_review"] is True
    assert data["genes"]["TEST"]["summary"]["pathogenic_status"] == "insufficient_evidence"


def test_column_where_every_read_is_deleted_is_piled_up(tmp_path):
    # 10M1D39M deletes 1-based position 11 in every read, so that column has no
    # base counts at all. It must still be reported, not abort the pileup.
    by_position = _pileup(tmp_path, [(0, "C" * 49, "10M1D39M")] * 20)
    assert 11 in by_position
    entry = by_position[11]
    assert entry["base_counts"]["indel"] == 20
    assert entry["coverage"] == 0
    assert entry["reference"] in ("A", "T", "G", "C", "N")
    assert entry["variant_alleles"] == []
    assert entry["has_variants"] is False
