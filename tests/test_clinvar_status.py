"""VAR-4: the ClinVar step records whether it ran, and the report says when it didn't.

Before: with the ClinVar database missing, or SnpSift failing, ROBIN copied the
unannotated snpEff file over ``snpsift_output.vcf`` and the report then said
"No pathogenic variants were identified." as if ClinVar had been checked.

Real SnpSift is not installed here (and must not be). Tests put a stand-in
executable named ``SnpSift`` first on ``PATH``.

Which tests can fail on the old code by assertion: the two report-text tests
(V1, V5) and the display-data test read VCF files written by hand, so on the old
code they run and fail on what the report/summary says. The tests that call
``annotate_with_clinvar`` / ``clinvar_status`` can only fail by ImportError there.
"""

import ast
import os
import stat
import sys
import tempfile
import textwrap
import logging
from pathlib import Path

import pytest

SNPEFF_VCF = textwrap.dedent(
    """\
    ##fileformat=VCFv4.2
    ##source=snpEff-stand-in
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
    chr1\t100\t.\tA\tG\t50\tPASS\tDP=10
    chr2\t200\t.\tC\tT\t50\tPASS\tDP=12
    """
)

PATHOGENIC_ROW = "chr9\t300\t.\tG\tA\t50\tPASS\tGENEINFO=BRCA9:1;CLNSIG=Pathogenic;CLNDN=Some_disease\n"

NOT_ANNOTATED_SENTENCE = "No pathogenic variants were identified."

KINDS = [
    ("SNP", "snpsift_output.vcf"),
    ("INDEL", "snpsift_indel_output.vcf"),
]


def _stand_in(bindir: Path, body: str) -> None:
    bindir.mkdir(parents=True, exist_ok=True)
    exe = bindir / "SnpSift"
    exe.write_text(f"#!{sys.executable}\n" + textwrap.dedent(body))
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)


def _use_path(monkeypatch, bindir: Path) -> None:
    # Only the stand-in's folder: nothing else named SnpSift can be found.
    monkeypatch.setenv("PATH", str(bindir))


PASSTHROUGH = """\
    import sys
    sys.stdout.write(open(sys.argv[-1]).read())
"""


def _annotate(snpeff: Path, final: Path, db="/fake/clinvar.vcf.gz", label="SNP"):
    from robin.analysis.target_analysis import annotate_with_clinvar

    return annotate_with_clinvar(
        str(snpeff), str(final), db, logging.getLogger("test"), label=label
    )


@pytest.fixture
def folder(tmp_path):
    d = tmp_path / "clair3"
    d.mkdir()
    (d / "snpeff.vcf").write_text(SNPEFF_VCF)
    return d


def _lines(path: Path):
    return path.read_text().splitlines()


def _no_leftovers(d: Path, final_names):
    assert sorted(p.name for p in d.iterdir()) == sorted(
        ["snpeff.vcf", *final_names]
    )


# ---------------------------------------------------------------- V1 .. V4, C4


@pytest.mark.parametrize("label,final_name", KINDS)
def test_v1_missing_database_is_recorded_unavailable(folder, label, final_name):
    final = folder / final_name
    status = _annotate(folder / "snpeff.vcf", final, db=None, label=label)
    assert status == "unavailable"
    lines = _lines(final)
    assert lines[0] == "##fileformat=VCFv4.2"
    assert lines[1] == "##ROBIN_ClinVar=unavailable;reason=ClinVar database not found"
    # the snpEff records are all still there
    assert lines[2:] == SNPEFF_VCF.splitlines()[1:]
    _no_leftovers(folder, [final_name])


@pytest.mark.parametrize("label,final_name", KINDS)
def test_v2_snpsift_nonzero_exit_is_recorded_unavailable(
    folder, tmp_path, monkeypatch, label, final_name
):
    _stand_in(
        tmp_path / "bin",
        """\
        import sys
        sys.stdout.write("##half-written\\n")
        sys.stderr.write("java.lang.OutOfMemoryError\\n")
        sys.exit(1)
        """,
    )
    _use_path(monkeypatch, tmp_path / "bin")
    final = folder / final_name
    status = _annotate(folder / "snpeff.vcf", final, label=label)
    assert status == "unavailable"
    lines = _lines(final)
    assert lines[1].startswith("##ROBIN_ClinVar=unavailable;reason=")
    assert "exited with code 1" in lines[1]
    assert "\n" not in lines[1] and "\t" not in lines[1]
    assert "##half-written" not in final.read_text()
    assert lines[2:] == SNPEFF_VCF.splitlines()[1:]
    _no_leftovers(folder, [final_name])


@pytest.mark.parametrize("label,final_name", KINDS)
def test_v3_snpsift_not_found_is_recorded_unavailable(
    folder, tmp_path, monkeypatch, label, final_name
):
    empty = tmp_path / "emptybin"
    empty.mkdir()
    _use_path(monkeypatch, empty)
    final = folder / final_name
    status = _annotate(folder / "snpeff.vcf", final, label=label)
    assert status == "unavailable"
    assert _lines(final)[1] == "##ROBIN_ClinVar=unavailable;reason=SnpSift not found"
    _no_leftovers(folder, [final_name])


def test_v4_while_snpsift_runs_the_previous_file_is_read_whole(
    folder, tmp_path, monkeypatch
):
    final = folder / "snpsift_output.vcf"
    previous = "##fileformat=VCFv4.2\n##ROBIN_ClinVar=annotated\n#CHROM\tPOS\nchr1\t1\n"
    final.write_text(previous)
    seen = tmp_path / "seen.txt"
    _stand_in(
        tmp_path / "bin",
        f"""\
        import sys
        open({str(seen)!r}, "w").write(open({str(final)!r}).read())
        sys.stdout.write(open(sys.argv[-1]).read())
        """,
    )
    _use_path(monkeypatch, tmp_path / "bin")
    assert _annotate(folder / "snpeff.vcf", final) == "annotated"
    assert seen.read_text() == previous  # whole, not truncated
    assert "chr2\t200" in final.read_text()  # and then replaced by the new one


def test_v4_failed_run_replaces_the_previous_file_only_by_a_whole_one(
    folder, tmp_path, monkeypatch
):
    final = folder / "snpsift_output.vcf"
    final.write_text("OLD\n")
    seen = tmp_path / "seen.txt"
    _stand_in(
        tmp_path / "bin",
        f"""\
        import sys
        open({str(seen)!r}, "w").write(open({str(final)!r}).read())
        sys.exit(2)
        """,
    )
    _use_path(monkeypatch, tmp_path / "bin")
    assert _annotate(folder / "snpeff.vcf", final) == "unavailable"
    assert seen.read_text() == "OLD\n"
    assert _lines(final)[0] == "##fileformat=VCFv4.2"


@pytest.mark.parametrize("label,final_name", KINDS)
def test_c4_success_without_pathogenic_record_is_annotated(
    folder, tmp_path, monkeypatch, label, final_name
):
    _stand_in(tmp_path / "bin", PASSTHROUGH)
    _use_path(monkeypatch, tmp_path / "bin")
    final = folder / final_name
    assert _annotate(folder / "snpeff.vcf", final, label=label) == "annotated"
    lines = _lines(final)
    assert lines[0] == "##fileformat=VCFv4.2"
    assert lines[1] == "##ROBIN_ClinVar=annotated"
    assert lines[2:] == SNPEFF_VCF.splitlines()[1:]
    _no_leftovers(folder, [final_name])


def test_stale_status_line_in_the_input_is_not_carried_over(
    folder, tmp_path, monkeypatch
):
    (folder / "snpeff.vcf").write_text(
        "##fileformat=VCFv4.2\n##ROBIN_ClinVar=annotated\n#CHROM\tPOS\nchr1\t1\n"
    )
    final = folder / "snpsift_output.vcf"
    _annotate(folder / "snpeff.vcf", final, db=None)
    text = final.read_text()
    assert text.count("##ROBIN_ClinVar=") == 1
    assert "##ROBIN_ClinVar=unavailable;reason=ClinVar database not found" in text


def test_status_line_goes_right_after_a_fileformat_line_that_is_not_first(folder):
    # M11: strict readers require ##fileformat= first; the status line must not go above it.
    (folder / "snpeff.vcf").write_text(
        "##source=x\n##fileformat=VCFv4.2\n##INFO=<ID=DP>\n#CHROM\tPOS\nchr1\t1\n"
    )
    final = folder / "snpsift_output.vcf"
    _annotate(folder / "snpeff.vcf", final, db=None)
    assert _lines(final) == [
        "##source=x",
        "##fileformat=VCFv4.2",
        "##ROBIN_ClinVar=unavailable;reason=ClinVar database not found",
        "##INFO=<ID=DP>",
        "#CHROM\tPOS",
        "chr1\t1",
    ]


def test_status_line_goes_first_without_a_fileformat_line(folder):
    (folder / "snpeff.vcf").write_text(
        "##source=x\n##ROBIN_ClinVar=annotated\n#CHROM\tPOS\nchr1\t1\n"
    )
    final = folder / "snpsift_output.vcf"
    _annotate(folder / "snpeff.vcf", final, db=None)
    assert _lines(final) == [
        "##ROBIN_ClinVar=unavailable;reason=ClinVar database not found",
        "##source=x",
        "#CHROM\tPOS",
        "chr1\t1",
    ]


def test_a_status_line_below_the_column_header_is_not_carried_over(folder):
    (folder / "snpeff.vcf").write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\n##ROBIN_ClinVar=annotated\nchr1\t1\n"
    )
    final = folder / "snpsift_output.vcf"
    _annotate(folder / "snpeff.vcf", final, db=None)
    assert _lines(final) == [
        "##fileformat=VCFv4.2",
        "##ROBIN_ClinVar=unavailable;reason=ClinVar database not found",
        "#CHROM\tPOS",
        "chr1\t1",
    ]


def test_empty_snpsift_output_is_unavailable(folder, tmp_path, monkeypatch):
    _stand_in(tmp_path / "bin", "pass\n")
    _use_path(monkeypatch, tmp_path / "bin")
    final = folder / "snpsift_output.vcf"
    assert _annotate(folder / "snpeff.vcf", final) == "unavailable"
    assert "no output" in _lines(final)[1]


# ---------------------------------------------------------------- the reader


def test_clinvar_status_reader(tmp_path):
    from robin.analysis.clinvar_status import clinvar_status

    def write(name, text):
        p = tmp_path / name
        p.write_text(text)
        return p

    assert clinvar_status(
        write("a.vcf", "##fileformat=VCFv4.2\n##ROBIN_ClinVar=annotated\n#CHROM\n")
    ) == ("annotated", "")
    assert clinvar_status(
        write(
            "u.vcf",
            "##fileformat=VCFv4.2\n##ROBIN_ClinVar=unavailable;reason=SnpSift not found\n",
        )
    ) == ("unavailable", "SnpSift not found")
    assert clinvar_status(write("n.vcf", "##fileformat=VCFv4.2\n#CHROM\n"))[0] == "not_recorded"
    assert clinvar_status(tmp_path / "missing.vcf")[0] == "not_recorded"
    # a line in the body does not count
    assert (
        clinvar_status(
            write("b.vcf", "##fileformat=VCFv4.2\n#CHROM\n##ROBIN_ClinVar=annotated\n")
        )[0]
        == "not_recorded"
    )
    assert (
        clinvar_status(write("x.vcf", "##ROBIN_ClinVar=garbage\n"))[0] == "not_recorded"
    )


# ---------------------------------------------------------------- the report


def _write_vcf(path: Path, status_line, rows=""):
    head = "##fileformat=VCFv4.2\n"
    if status_line:
        head += status_line + "\n"
    head += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
    path.write_text(head + "chr1\t100\t.\tA\tG\t50\tPASS\tDP=10\tGT\t0/1\n" + rows)


def _report_text(output: Path):
    from types import SimpleNamespace

    from robin.reporting.sections.variants import VariantsSection
    from robin.reporting.styling.styles import ReportStyles

    report = SimpleNamespace(
        output=str(output), styles=ReportStyles(tempfile.mkdtemp())
    )
    section = VariantsSection(report)
    section.add_content()
    texts = []
    for el in section.elements + section.summary_elements:
        if hasattr(el, "getPlainText"):
            texts.append(el.getPlainText())
    return "\n".join(texts), section


def _clair(tmp_path, snp_line, snp_rows="", indel_line="##ROBIN_ClinVar=annotated"):
    out = tmp_path / "sample"
    (out / "clair3").mkdir(parents=True)
    _write_vcf(out / "clair3" / "snpsift_output.vcf", snp_line, snp_rows)
    _write_vcf(out / "clair3" / "snpsift_indel_output.vcf", indel_line)
    return out


def test_v1_report_says_unavailable_not_none_found(tmp_path):
    out = _clair(
        tmp_path, "##ROBIN_ClinVar=unavailable;reason=ClinVar database not found"
    )
    text, _ = _report_text(out)
    assert NOT_ANNOTATED_SENTENCE not in text
    assert "ClinVar annotation was unavailable for the SNP file" in text
    assert "ClinVar database not found" in text
    assert "Found 0 candidate pathogenic SNPs" not in text


def test_v1_report_both_files_unavailable(tmp_path):
    line = "##ROBIN_ClinVar=unavailable;reason=SnpSift not found"
    out = _clair(tmp_path, line, indel_line=line)
    text, _ = _report_text(out)
    assert NOT_ANNOTATED_SENTENCE not in text
    assert "unavailable for the SNP file" in text
    assert "unavailable for the INDEL file" in text


def test_v5_report_says_not_recorded_for_a_file_without_the_line(tmp_path):
    out = _clair(tmp_path, None)
    text, _ = _report_text(out)
    assert NOT_ANNOTATED_SENTENCE not in text
    assert "not recorded" in text


def test_c4_report_says_none_found_when_annotated(tmp_path):
    out = _clair(tmp_path, "##ROBIN_ClinVar=annotated")
    text, _ = _report_text(out)
    assert NOT_ANNOTATED_SENTENCE in text
    assert "unavailable" not in text
    assert "not recorded" not in text
    assert "Found 0 candidate pathogenic SNPs and 0 candidate pathogenic indels" in text


def test_pathogenic_rows_are_still_listed_when_the_other_file_was_not_annotated(
    tmp_path,
):
    out = _clair(
        tmp_path,
        "##ROBIN_ClinVar=annotated",
        snp_rows=PATHOGENIC_ROW,
        indel_line="##ROBIN_ClinVar=unavailable;reason=SnpSift not found",
    )
    text, section = _report_text(out)
    assert len(section.variant_result.snp_data) == 1
    assert "Found 1 candidate pathogenic SNPs" in text
    assert "BRCA9" in section.variant_result.affected_genes
    assert "unavailable for the INDEL file" in text
    assert NOT_ANNOTATED_SENTENCE not in text
    assert "Detailed Analysis" in text


def test_not_recorded_file_with_pathogenic_rows_keeps_its_count_and_says_they_are_annotated(
    tmp_path,
):
    # I2 / P5: a pre-batch SNP file (no status line) that lists a CLNSIG=Pathogenic row. Such a
    # row can only come from ClinVar, so "could not be assessed" would contradict the list.
    out = _clair(tmp_path, None, snp_rows=PATHOGENIC_ROW)
    text, section = _report_text(out)
    assert len(section.variant_result.snp_data) == 1
    assert "Found 1 candidate pathogenic SNPs" in text
    assert "BRCA9" in text
    assert "ClinVar annotation status was not recorded for the SNP file" in text
    assert "the SNP variants listed carry ClinVar annotations" in text
    assert "could not be assessed" not in text
    assert NOT_ANNOTATED_SENTENCE not in text


def test_not_recorded_files_old_folder_rows_in_one_file_only(tmp_path):
    # A typical old folder: neither file has a status line; only the SNP file has a row.
    out = _clair(tmp_path, None, snp_rows=PATHOGENIC_ROW, indel_line=None)
    text, _ = _report_text(out)
    assert "Found 1 candidate pathogenic SNPs" in text
    assert "the SNP variants listed carry ClinVar annotations" in text
    assert "pathogenic SNP variants could not be assessed" not in text
    # the INDEL file has no rows: its sentence stays, and no count is claimed for it
    assert "pathogenic INDEL variants could not be assessed" in text
    assert "0 candidate pathogenic indels" not in text


def test_report_does_not_claim_zero_for_the_unannotated_file(tmp_path):
    out = _clair(
        tmp_path,
        "##ROBIN_ClinVar=annotated",
        indel_line="##ROBIN_ClinVar=unavailable;reason=SnpSift not found",
    )
    text, _ = _report_text(out)
    assert "0 candidate pathogenic indels" not in text
    assert "0 candidate pathogenic SNPs" in text


# ---------------------------------------------------------------- GUI data


@pytest.mark.parametrize(
    "line,expected",
    [
        ("##ROBIN_ClinVar=annotated", "annotated"),
        ("##ROBIN_ClinVar=unavailable;reason=SnpSift not found", "unavailable"),
        (None, "not_recorded"),
    ],
)
def test_display_data_summary_carries_the_clinvar_status(tmp_path, line, expected):
    from robin.analysis.snp_processing import build_snp_display_data

    vcf = tmp_path / "snpsift_output.vcf"
    _write_vcf(vcf, line)
    data = build_snp_display_data(vcf)
    assert data["summary"]["clinvar"] == expected
    assert data["summary"]["total_variants"] == 1
    if expected == "unavailable":
        assert data["summary"]["clinvar_reason"] == "SnpSift not found"


# ---------------------------------------------------------------- run_snp_analysis source


def _run_snp_analysis_tree():
    import robin.analysis.target_analysis as ta

    tree = ast.parse(Path(ta.__file__).read_text())
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == "run_snp_analysis":
            return node
    raise AssertionError("run_snp_analysis not found")


def test_run_snp_analysis_annotates_both_files_with_the_one_function():
    fn = _run_snp_analysis_tree()
    calls = [
        n
        for n in ast.walk(fn)
        if isinstance(n, ast.Call)
        and isinstance(n.func, ast.Name)
        and n.func.id == "annotate_with_clinvar"
    ]
    assert len(calls) == 2
    first_args = {ast.unparse(c.args[0]) for c in calls}
    second_args = {ast.unparse(c.args[1]) for c in calls}
    assert first_args == {"snpeff_out", "snpeff_indel_out"}
    assert second_args == {"snpsift_out", "snpsift_indel_out"}


def test_run_snp_analysis_never_copies_snpeff_output_to_the_final_file():
    fn = _run_snp_analysis_tree()
    for n in ast.walk(fn):
        if (
            isinstance(n, ast.Call)
            and isinstance(n.func, ast.Attribute)
            and n.func.attr in {"copy", "copy2", "copyfile", "move"}
        ):
            assert not any("snpeff" in ast.unparse(a) for a in n.args), ast.unparse(n)
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "open":
            if len(n.args) >= 2 and ast.unparse(n.args[0]) in {
                "snpsift_out",
                "snpsift_indel_out",
            }:
                raise AssertionError("final file opened directly: " + ast.unparse(n))


def test_final_file_is_replaced_only_by_rename():
    import robin.analysis.target_analysis as ta

    tree = ast.parse(Path(ta.__file__).read_text())
    fn = next(
        n
        for n in ast.walk(tree)
        if isinstance(n, ast.FunctionDef) and n.name == "annotate_with_clinvar"
    )
    replaces = [
        ast.unparse(n)
        for n in ast.walk(fn)
        if isinstance(n, ast.Call)
        and isinstance(n.func, ast.Attribute)
        and n.func.attr == "replace"
    ]
    assert replaces == ["os.replace(out_path, final_vcf)"]
    for n in ast.walk(fn):
        if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == "open":
            assert ast.unparse(n.args[0]) != "final_vcf"
