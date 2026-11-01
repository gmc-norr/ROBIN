"""VAR-1 and VAR-2: which ClinVar records count as pathogenic, and that each one is checked.

The rule (user's choice, 2026-11-01): a record counts when ClinVar's own aggregate label says so in
any of its three classifications. CLNSIG (germline) must be Pathogenic or Likely_pathogenic, ONC
(oncogenicity) Oncogenic or Likely_oncogenic, and SCI (somatic clinical impact) Tier I or Tier II.
Labels are matched whole. "Conflicting_classifications_of_pathogenicity" contains the word
"pathogenic" but is not a pathogenic call, and never counts by itself.

The INFO values below are copied from ClinVar's GRCh38 VCF (fileDate 2026-10-25), keeping only the
classification keys and GENEINFO. They are public reference data, not patient data. The reads are
synthetic.

On the old code, the helper test fails by ImportError. Every other fault test runs on the old
code and fails by assertion.
"""

import json
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pytest

from lga_fixtures import write_bam, write_bed

# ClinVar 2026-10-25, chr7:140753336 A>T (BRAF V600E). Germline label: conflicting.
BRAF_V600E = (
    "CLNSIG=Conflicting_classifications_of_pathogenicity;"
    "CLNSIGCONF=Pathogenic_(4)|Likely_pathogenic_(1)|Uncertain_significance_(1);"
    "ONC=Oncogenic;SCI=Tier_I_-_Strong;GENEINFO=BRAF:673"
)
# chr5:1295113 G>A (TERT C228T). Germline label: conflicting.
TERT_C228T = (
    "CLNSIG=Conflicting_classifications_of_pathogenicity;"
    "CLNSIGCONF=Likely_pathogenic_(4)|Uncertain_significance_(1);"
    "ONC=Oncogenic;SCI=Tier_I_-_Strong;GENEINFO=TERT:7015|LOC110806263:110806263"
)
# chr5:1295113 G>T: a somatic Tier I label only.
TERT_SCI_ONLY = "SCI=Tier_I_-_Strong;GENEINFO=TERT:7015|LOC110806263:110806263"
# chr1:226064454 G>A in H3-3A: germline label "other", oncogenicity and somatic labels.
H3_3A = "CLNSIG=other;ONC=Likely_oncogenic;SCI=Tier_I_-_Strong;GENEINFO=H3-3A:3020"
# chr12:25245350 CC>AA in KRAS: an oncogenicity label only.
KRAS_ONC_ONLY = "ONC=Oncogenic;GENEINFO=KRAS:3845"
# chr17:7668477 G>A in TP53: conflicting, between "uncertain" and "likely benign".
TP53_CONFLICTING = (
    "CLNSIG=Conflicting_classifications_of_pathogenicity;"
    "CLNSIGCONF=Uncertain_significance_(1)|Likely_benign_(1);GENEINFO=TP53:7157"
)
# chr17:31095037 A>T in NF1: conflicting, with more pathogenic than benign submissions.
NF1_CONFLICTING = (
    "CLNSIG=Conflicting_classifications_of_pathogenicity;"
    "CLNSIGCONF=Likely_pathogenic_(1)|Uncertain_significance_(1);"
    "GENEINFO=NF1:4763|MIR4733HG:105371722|LOC111811965:111811965"
)
PLAIN_PATHOGENIC = "CLNSIG=Pathogenic;GENEINFO=TESTG:1"
BENIGN = "CLNSIG=Benign;GENEINFO=TESTG:1"
NOT_IN_CLINVAR = "DP=10"

# (name, INFO, counts as pathogenic)
RECORDS = [
    ("braf_v600e_conflicting_germline_oncogenic", BRAF_V600E, True),
    ("tert_c228t_conflicting_germline_oncogenic", TERT_C228T, True),
    ("tert_somatic_tier_only", TERT_SCI_ONLY, True),
    ("h3_3a_other_germline_oncogenic", H3_3A, True),
    ("kras_oncogenicity_only", KRAS_ONC_ONLY, True),
    ("tp53_conflicting_only", TP53_CONFLICTING, False),
    ("nf1_conflicting_pathogenic_majority", NF1_CONFLICTING, False),
    ("plain_pathogenic", PLAIN_PATHOGENIC, True),
    ("benign", BENIGN, False),
    ("not_in_clinvar", NOT_IN_CLINVAR, False),
]


def _info_dict(info):
    out = {}
    for part in info.split(";"):
        key, sep, value = part.partition("=")
        out[key] = value if sep else True
    return out


# ------------------------------------------------------------ the shared rule


@pytest.mark.parametrize(
    "field,value,expected",
    [
        # every CLNSIG value seen with "athogenic" in it, most common first, and the rest
        ("CLNSIG", "Pathogenic", True),
        ("CLNSIG", "Conflicting_classifications_of_pathogenicity", False),
        ("CLNSIG", "Likely_pathogenic", True),
        ("CLNSIG", "Pathogenic/Likely_pathogenic", True),
        ("CLNSIG", "Pathogenic,_low_penetrance", True),
        ("CLNSIG", "Pathogenic|other", True),
        ("CLNSIG", "Likely_pathogenic,_low_penetrance", True),
        ("CLNSIG", "Conflicting_classifications_of_pathogenicity|other", False),
        ("CLNSIG", "Likely_pathogenic/Likely_risk_allele", True),
        ("CLNSIG", "Conflicting_classifications_of_pathogenicity|association|risk_factor", False),
        ("CLNSIG", "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance/Established_risk_allele|risk_factor", True),
        ("CLNSIG", "Uncertain_significance", False),
        ("CLNSIG", "Likely_benign", False),
        ("CLNSIG", "Benign/Likely_benign", False),
        ("CLNSIG", "Likely_risk_allele", False),
        ("CLNSIG", "risk_factor", False),
        ("CLNSIG", "other", False),
        ("CLNSIG", "not_provided", False),
        ("CLNSIG", "no_classification_for_the_single_variant", False),
        ("ONC", "Oncogenic", True),
        ("ONC", "Likely_oncogenic", True),
        ("ONC", "Oncogenic/Likely_oncogenic", True),
        ("ONC", "Conflicting_classifications_of_oncogenicity", False),
        ("ONC", "Uncertain_significance", False),
        ("ONC", "Likely_benign", False),
        ("SCI", "Tier_I_-_Strong", True),
        ("SCI", "Tier_II_-_Potential", True),
        ("SCI", "Tier_III_-_Unknown", False),
        ("SCI", "Tier_IV_-_Benign/Likely_benign", False),
        ("SCI", "no_classification_for_the_single_variant", False),
        # a label under another key does not count
        ("CLNSIGCONF", "Pathogenic_(4)", False),
        ("CLNSIGINCL", "12345:Pathogenic", False),
    ],
)
def test_the_rule_on_clinvar_labels(field, value, expected):
    from robin.analysis.clinvar_significance import is_clinvar_pathogenic

    assert is_clinvar_pathogenic(f"{field}={value};GENEINFO=X:1") is expected
    assert is_clinvar_pathogenic({field: value, "GENEINFO": "X:1"}) is expected


# ------------------------------------------------------------ the three checks agree


def _lga_check(info):
    from robin.analysis.lightweight_gene_analysis import LightweightGeneAnalysis

    record = SimpleNamespace(
        chrom="1", pos=11, id="1", ref="C", alts=["T"], info=_info_dict(info)
    )
    return LightweightGeneAnalysis.__new__(LightweightGeneAnalysis)._is_pathogenic(record)


def _snp_table_check(info):
    from robin.analysis import snp_processing

    return snp_processing._is_pathogenic(info)


def _report_check(info):
    from robin.reporting.sections.variants import VariantAnalysis

    return VariantAnalysis.__new__(VariantAnalysis)._is_pathogenic(info)


CHECKS = {"variant_check": _lga_check, "snp_table": _snp_table_check, "report": _report_check}


@pytest.mark.parametrize("check", list(CHECKS))
@pytest.mark.parametrize("name,info,expected", RECORDS, ids=[r[0] for r in RECORDS])
def test_every_check_follows_the_rule(check, name, info, expected):
    assert bool(CHECKS[check](info)) is expected


# ------------------------------------------------------------ the report


def _annotated_vcf(path, rows):
    """rows: (chrom, pos, ref, alt, info). A SnpSift-annotated file with the ClinVar status line."""
    head = (
        "##fileformat=VCFv4.2\n##ROBIN_ClinVar=annotated\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    )
    body = "".join(
        f"{c}\t{p}\t.\t{r}\t{a}\t50\tPASS\t{i}\tGT\t0/1\n" for c, p, r, a, i in rows
    )
    path.write_text(head + body)
    return path


REPORT_ROWS = [
    ("chr7", 140753336, "A", "T", BRAF_V600E),
    ("chr5", 1295113, "G", "A", TERT_C228T),
    ("chr5", 1295113, "G", "T", TERT_SCI_ONLY),
    ("chr1", 226064454, "G", "A", H3_3A),
    ("chr17", 7668477, "G", "A", TP53_CONFLICTING),
    ("chr17", 31095037, "A", "T", NF1_CONFLICTING),
    ("chr9", 300, "G", "A", PLAIN_PATHOGENIC),
]
COUNTED = {
    ("chr7", 140753336, "T"),
    ("chr5", 1295113, "A"),
    ("chr5", 1295113, "T"),
    ("chr1", 226064454, "A"),
    ("chr9", 300, "A"),
}


def _report_section(output):
    from robin.reporting.sections.variants import VariantsSection
    from robin.reporting.styling.styles import ReportStyles

    report = SimpleNamespace(output=str(output), styles=ReportStyles(tempfile.mkdtemp()))
    section = VariantsSection(report)
    section.add_content()
    text = "\n".join(
        el.getPlainText()
        for el in section.elements + section.summary_elements
        if hasattr(el, "getPlainText")
    )
    return section, text


def _clair(tmp_path, snp_rows):
    out = tmp_path / "sample"
    (out / "clair3").mkdir(parents=True)
    _annotated_vcf(out / "clair3" / "snpsift_output.vcf", snp_rows)
    _annotated_vcf(out / "clair3" / "snpsift_indel_output.vcf", [])
    return out


def test_the_report_lists_the_records_clinvar_calls_pathogenic(tmp_path):
    section, text = _report_section(_clair(tmp_path, REPORT_ROWS))
    listed = {
        (v["chromosome"], v["position"], v["alternate"]) for v in section.variant_result.snp_data
    }
    assert listed == COUNTED
    assert section.variant_result.affected_genes == {"BRAF", "TERT", "H3-3A", "TESTG"}
    # The note under the table states the rule.
    assert "Conflicting classifications" in text
    assert "Tier I or II" in text


def test_control_report_with_plain_pathogenic_rows(tmp_path):
    # Passes before and after: a plain Pathogenic row is listed, a Benign one is not.
    rows = [("chr9", 300, "G", "A", PLAIN_PATHOGENIC), ("chr9", 400, "G", "A", BENIGN)]
    section, _ = _report_section(_clair(tmp_path, rows))
    assert [(v["position"], v["alternate"]) for v in section.variant_result.snp_data] == [(300, "A")]


# ------------------------------------------------------------ the GUI's SNP table

ANN = "ANN=A|missense_variant|MODERATE|{g}|{g}|transcript|T1|protein_coding|1/1|c.1G>A|p.X1Y|1/1|1/1|1/1||"


def test_the_snp_table_flags_the_records_clinvar_calls_pathogenic(tmp_path):
    from robin.analysis.snp_processing import build_snp_display_data

    rows = [(c, p, r, a, i + ";" + ANN.format(g="G" + str(n))) for n, (c, p, r, a, i) in enumerate(REPORT_ROWS)]
    data = build_snp_display_data(_annotated_vcf(tmp_path / "snpsift_output.vcf", rows))
    flagged = {(r["CHROM"], int(r["POS"]), r["ALT"]) for r in data["rows_pathogenic"]}
    assert flagged == COUNTED
    assert data["summary"]["total_variants"] == len(REPORT_ROWS)
    assert data["summary"]["pathogenic_variants"] == len(COUNTED)


def test_a_record_without_a_gene_annotation_stays_in_the_snp_table(tmp_path):
    # VAR-2's second part: rows without ANN had no "Allele" and groupby dropped them.
    from robin.analysis.snp_processing import build_snp_display_data, parse_vcf

    rows = [
        ("chr9", 300, "G", "A", PLAIN_PATHOGENIC + ";" + ANN.format(g="G1")),
        ("chr9", 400, "G", "A", PLAIN_PATHOGENIC),  # no ANN
        ("chr9", 500, "C", "T", "DP=10;" + ANN.format(g="G2")),
    ]
    vcf = _annotated_vcf(tmp_path / "snpsift_output.vcf", rows)
    assert sorted(int(p) for p in parse_vcf(vcf)["POS"]) == [300, 400, 500]
    data = build_snp_display_data(vcf)
    assert sorted(int(r["POS"]) for r in data["rows_all"]) == [300, 400, 500]
    assert sorted(int(r["POS"]) for r in data["rows_pathogenic"]) == [300, 400]


def test_control_every_record_annotated_keeps_its_rows(tmp_path):
    # Passes before and after: with ANN on every row nothing was dropped.
    from robin.analysis.snp_processing import build_snp_display_data

    rows = [
        ("chr9", 300, "G", "A", PLAIN_PATHOGENIC + ";" + ANN.format(g="G1")),
        ("chr9", 500, "C", "T", "DP=10;" + ANN.format(g="G2")),
    ]
    data = build_snp_display_data(_annotated_vcf(tmp_path / "snpsift_output.vcf", rows))
    assert sorted(int(r["POS"]) for r in data["rows_all"]) == [300, 500]
    assert [int(r["POS"]) for r in data["rows_pathogenic"]] == [300]


def test_the_gui_snp_table_shows_the_oncogenicity_and_somatic_labels():
    # Source check: the table's columns are picked inside a nested GUI function. A row can be
    # flagged by ONC or SCI while CLNSIG says "Conflicting", so both must be visible.
    import ast

    import robin.gui.components.snp as snp

    tree = ast.parse(Path(snp.__file__).read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and any(
            getattr(t, "id", None) == "preferred_display_fields" for t in node.targets
        ):
            fields = [e.value for e in node.value.elts]
            break
    else:
        raise AssertionError("preferred_display_fields not found")
    at = fields.index("CLNSIG")
    assert fields[at + 1 : at + 3] == ["ONC", "SCI"]


# ------------------------------------------------------------ the variant check (LGA)


def _clinvar(path, records):
    """records: (pos, ref, alt, info) on contig 1, in file order."""
    lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    lines += [f"1\t{p}\t{n}\t{r}\t{a}\t.\t.\t{i}" for n, (p, r, a, i) in enumerate(records, 1)]
    Path(path).write_text("\n".join(lines) + "\n")
    return path


def _variant_check(tmp_path, records, reads, gene):
    from robin.analysis.lightweight_gene_analysis import LightweightGeneAnalysis

    bam = write_bam(tmp_path / "s.bam", reads)
    bed = write_bed(tmp_path / "s.bed", [(50, 150, gene)])
    vcf = _clinvar(tmp_path / "clinvar.vcf", records)
    analysis = LightweightGeneAnalysis(str(tmp_path / "work"), str(vcf))
    reports = analysis.analyze_bed_regions(str(bam), str(bed))
    out = tmp_path / "results.json"
    analysis.export_json(reports, str(out))
    return json.loads(out.read_text())


def _by_allele(data, gene):
    return {
        (v["position"], v["reference"], v["alternate"]): v for v in data["genes"][gene]["variants"]
    }


def _reads(base_at_101, n=20):
    return [(0, "G" * 100 + base_at_101 + "G" * 49)] * n


# ClinVar 2026-10-25, chr2:208248389 (IDH1 c.394), in the file's order: G>A (R132C), G>C, G>T.
IDH1_C394 = [
    (101, "G", "A", "CLNSIG=Pathogenic/Likely_pathogenic;ONC=Oncogenic;SCI=Tier_I_-_Strong;GENEINFO=IDH1:3417"),
    (101, "G", "C", "ONC=Oncogenic;SCI=Tier_I_-_Strong;GENEINFO=IDH1:3417"),
    (101, "G", "T", "CLNSIG=Pathogenic;ONC=Oncogenic;SCI=Tier_II_-_Potential;GENEINFO=IDH1:3417"),
]


def test_every_clinvar_allele_at_one_position_is_checked(tmp_path):
    # VAR-1. Old code keeps only the last record (G>T): an R132C-type sample, every read A,
    # read "reference_mismatch" and the gene "insufficient_evidence".
    data = _variant_check(tmp_path, IDH1_C394, _reads("A"), "IDH1")
    alleles = _by_allele(data, "IDH1")
    assert set(alleles) == {(101, "G", "A"), (101, "G", "C"), (101, "G", "T")}
    r132c = alleles[(101, "G", "A")]
    assert r132c["evidence_analysis"]["alternate_support"] == 20
    assert r132c["evidence_analysis"]["evidence_level"] == "high_frequency"
    assert r132c["supported"] is True
    for other in [(101, "G", "C"), (101, "G", "T")]:
        assert alleles[other]["supported"] is False
    assert data["genes"]["IDH1"]["summary"]["pathogenic_status"] == "present_with_good_coverage"


# chr12:25245350 (KRAS c.35), in the file's order: C>A, C>G, C>T, CC>AA (moved to 101/102).
KRAS_C35 = [
    (101, "C", "A", "CLNSIG=Pathogenic;ONC=Oncogenic;SCI=Tier_I_-_Strong;GENEINFO=KRAS:3845"),
    (101, "C", "G", "CLNSIG=Pathogenic/Likely_pathogenic;ONC=Oncogenic;SCI=Tier_I_-_Strong;GENEINFO=KRAS:3845"),
    (101, "C", "T", "CLNSIG=Pathogenic/Likely_pathogenic;ONC=Oncogenic;SCI=Tier_I_-_Strong;GENEINFO=KRAS:3845"),
    (101, "CC", "AA", KRAS_ONC_ONLY),
]


def test_a_snv_next_to_a_multi_base_record_at_the_same_position_is_checked(tmp_path):
    # VAR-1. Old code keeps only CC>AA, which a one-base pileup can't call: a G12D-type sample
    # (every read T) is never checked for C>T.
    reads = [(0, "C" * 100 + "T" + "C" * 49)] * 20
    data = _variant_check(tmp_path, KRAS_C35, reads, "KRAS")
    alleles = _by_allele(data, "KRAS")
    assert len(alleles) == 4
    assert alleles[(101, "C", "T")]["supported"] is True
    assert alleles[(101, "CC", "AA")]["evidence_analysis"]["evidence_level"] == "mnv_requires_review"


def test_the_variant_check_uses_the_rule_and_shows_why(tmp_path):
    # VAR-2 in the variant check. Old code checks the TP53-type conflicting record and skips the
    # H3-3A-type one (its CLNSIG "other" hides its oncogenicity label).
    records = [
        (101, "G", "A", TP53_CONFLICTING.replace("TP53:7157", "TESTG:1")),
        (121, "G", "A", H3_3A.replace("H3-3A:3020", "TESTG:1")),
    ]
    reads = [(0, "G" * 100 + "A" + "G" * 19 + "A" + "G" * 29)] * 20
    data = _variant_check(tmp_path, records, reads, "TESTG")
    alleles = _by_allele(data, "TESTG")
    assert set(alleles) == {(121, "G", "A")}
    h3 = alleles[(121, "G", "A")]
    assert h3["supported"] is True
    assert h3["clinical_significance"] == (
        "other; oncogenicity: Likely_oncogenic; somatic impact: Tier_I_-_Strong"
    )


def test_control_a_plain_pathogenic_record_keeps_its_label(tmp_path):
    # Passes before and after: one record per position, CLNSIG only, label shown as written.
    data = _variant_check(tmp_path, [(101, "G", "A", PLAIN_PATHOGENIC)], _reads("A"), "TESTG")
    v = _by_allele(data, "TESTG")[(101, "G", "A")]
    assert v["supported"] is True
    assert v["clinical_significance"] == "Pathogenic"


# ------------------------------------------------------------ external review, round 2

# ClinVar 2026-10-25, chr17:43067696 in BRCA1, in the file's order: C>. (ID 125742, HGVS
# "g.43067696="), then C>A, C>G, C>T. ALT "." means ClinVar gives no alternate allele: the record
# is the reference base itself, which a pileup can't check. 13 selected records are like this.
BRCA1_43067696 = [
    (101, "C", ".", "CLNSIG=Pathogenic;GENEINFO=BRCA1:672"),
    (101, "C", "A", "CLNSIG=Pathogenic;GENEINFO=BRCA1:672"),
    (101, "C", "G", "CLNSIG=Pathogenic;GENEINFO=BRCA1:672"),
    (101, "C", "T", "CLNSIG=Pathogenic;GENEINFO=BRCA1:672"),
]
BRCA1_NEIGHBOUR = (121, "C", "T", "CLNSIG=Pathogenic;GENEINFO=BRCA1:672")
# chr10:87864144 in PTEN: G>. (ID 1705890) is alone among the selected records there (G>C is Benign).
PTEN_87864144 = [
    (101, "G", ".", "CLNSIG=Pathogenic;GENEINFO=PTEN:5728"),
    (101, "G", "C", "CLNSIG=Benign;GENEINFO=PTEN:5728"),
]
PTEN_NEIGHBOUR = (121, "G", "A", "CLNSIG=Pathogenic;GENEINFO=PTEN:5728")


def _neighbour_reads(background, alt_at_121, covered):
    """Reads carrying ``alt_at_121`` at 121. ``covered``: the reads also cover 101."""
    if covered:
        return [(0, background * 120 + alt_at_121 + background * 29)] * 20
    return [(110, background * 10 + alt_at_121 + background * 29)] * 20  # 111..150


@pytest.mark.parametrize("covered", [True, False], ids=["site_covered", "site_uncovered"])
def test_a_record_without_an_alternate_allele_does_not_drop_its_region(tmp_path, covered):
    # Review finding 1. On cf4f84e the position key hid C>. behind C>T (and C>A, C>G with it);
    # on 502fd72 C>. was kept, raised in the evidence code, and the whole BRCA1 region vanished.
    data = _variant_check(
        tmp_path, BRCA1_43067696 + [BRCA1_NEIGHBOUR], _neighbour_reads("C", "T", covered), "BRCA1"
    )
    assert "BRCA1" in data["genes"]
    alleles = _by_allele(data, "BRCA1")
    assert set(alleles) == {(101, "C", "A"), (101, "C", "G"), (101, "C", "T"), (121, "C", "T")}
    assert alleles[(121, "C", "T")]["supported"] is True
    assert data["genes"]["BRCA1"]["summary"]["pathogenic_status"] == "present_with_good_coverage"
    assert data["metadata"]["clinvar_records_not_checkable"] == 1


@pytest.mark.parametrize("covered", [True, False], ids=["site_covered", "site_uncovered"])
def test_a_lone_record_without_an_alternate_allele_does_not_drop_its_region(tmp_path, covered):
    # Review finding 1, the part that predates batch 9: a lone G>. dropped its region on cf4f84e too.
    data = _variant_check(
        tmp_path, PTEN_87864144 + [PTEN_NEIGHBOUR], _neighbour_reads("G", "A", covered), "PTEN"
    )
    assert "PTEN" in data["genes"]
    assert set(_by_allele(data, "PTEN")) == {(121, "G", "A")}
    assert data["genes"]["PTEN"]["summary"]["pathogenic_status"] == "present_with_good_coverage"
    assert data["metadata"]["clinvar_records_not_checkable"] == 1


def _variant_check_with_html(tmp_path, records, reads, gene):
    from robin.analysis.lightweight_gene_analysis import LightweightGeneAnalysis

    bam = write_bam(tmp_path / "s.bam", reads)
    bed = write_bed(tmp_path / "s.bed", [(50, 150, gene)])
    vcf = _clinvar(tmp_path / "clinvar.vcf", records)
    analysis = LightweightGeneAnalysis(str(tmp_path / "work"), str(vcf))
    reports = analysis.analyze_bed_regions(str(bam), str(bed))
    data = json.loads(Path(analysis.export_json(reports, str(tmp_path / "r.json"))).read_text())
    html = Path(analysis.generate_report(reports, str(tmp_path / "r.html"))).read_text()
    return data, html


@pytest.mark.parametrize(
    "reads,status,found,phrase",
    [
        (_reads("A"), "present_with_good_coverage", 1, "PATHOGENIC SNPS PRESENT (1 with good coverage)"),
        (_reads("G"), "none_detected", 0, "NO PATHOGENIC SNPS DETECTED"),
        (_reads("G", 18) + _reads("A", 2), "insufficient_evidence", 0, "INSUFFICIENT EVIDENCE"),
    ],
    ids=["all_alt", "reference_only", "two_alt_reads"],
)
def test_the_html_export_counts_what_the_json_counts(tmp_path, reads, status, found, phrase):
    # Review finding 2 (predates batch 9; VAR-1 made it larger): the HTML export counted every
    # checked record with reads as a found pathogenic variant (3 here), whatever the reads showed.
    data, html = _variant_check_with_html(tmp_path, IDH1_C394, reads, "IDH1")
    summary = data["genes"]["IDH1"]["summary"]
    assert (summary["pathogenic_status"], summary["total_variants"]) == (status, found)
    assert phrase in html
    assert f"<strong>Total Pathogenic Variants:</strong> {found} " in html
