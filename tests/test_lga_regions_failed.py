"""Batch 10e: the variant check names the regions it could not check.

Before, an error inside one region dropped that region: its genes vanished from the
results, the job reported success and the GUI said "completed successfully!". A BED that
could not be read was read as no regions. Synthetic reads and ClinVar records only.
"""

import ast
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

from lga_fixtures import write_bam, write_bed
from robin.analysis import lightweight_gene_analysis as lga
from robin.analysis.lightweight_gene_analysis import LightweightGeneAnalysis

RESULTS = "lightweight_gene_analysis_results.json"
READS = [(0, "C" * 50)] * 20 + [(100, "C" * 50)] * 20
REGIONS = [(0, 50, "GENE1"), (100, 150, "GENE2")]


def _clinvar(path):
    lines = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        "1\t11\t1\tC\tT\t.\t.\tCLNSIG=Pathogenic;GENEINFO=GENE1:1;CLNDN=test",
        "1\t111\t2\tC\tT\t.\t.\tCLNSIG=Pathogenic;GENEINFO=GENE2:2;CLNDN=test",
    ]
    Path(path).write_text("\n".join(lines) + "\n")
    return path


def _fail_region(monkeypatch, name):
    real = LightweightGeneAnalysis._analyze_single_region

    def analyse(self, bam_path, region_id, variants):
        if region_id == name:
            raise RuntimeError("synthetic region failure")
        return real(self, bam_path, region_id, variants)

    monkeypatch.setattr(LightweightGeneAnalysis, "_analyze_single_region", analyse)


def _check(tmp_path, bed_rows=REGIONS):
    bam = write_bam(tmp_path / "s.bam", READS)
    bed = write_bed(tmp_path / "s.bed", bed_rows)
    analysis = LightweightGeneAnalysis(str(tmp_path / "work"), str(_clinvar(tmp_path / "clinvar.vcf")))
    reports = analysis.analyze_bed_regions(str(bam), str(bed))
    out = tmp_path / RESULTS
    analysis.export_json(reports, str(out))
    return json.loads(out.read_text())


def test_a_region_that_fails_is_listed_and_the_others_are_kept(tmp_path, monkeypatch):
    _fail_region(monkeypatch, "GENE2")
    data = _check(tmp_path)
    assert "GENE1" in data["genes"]
    failed = data["metadata"]["regions_failed"]
    assert [f["region"] for f in failed] == ["GENE2"]
    assert failed[0]["genes"] == ["GENE2"]
    assert failed[0]["intervals"] == ["chr1:100-150"]
    assert "synthetic region failure" in failed[0]["error"]


def test_a_clean_check_lists_no_failed_region(tmp_path):
    data = _check(tmp_path)
    assert data["metadata"]["regions_failed"] == []
    assert set(data["genes"]) == {"GENE1", "GENE2"}


def test_a_bed_line_that_is_not_a_bed_row_stops_the_check(tmp_path):
    bam = write_bam(tmp_path / "s.bam", READS)
    bed = tmp_path / "bad.bed"
    bed.write_text("chr1\t0\t50\tGENE1\nchr1\tabc\t150\tGENE2\n")
    analysis = LightweightGeneAnalysis(str(tmp_path / "work"), str(_clinvar(tmp_path / "clinvar.vcf")))
    with pytest.raises(ValueError, match="line 2"):
        analysis.analyze_bed_regions(str(bam), str(bed))


def _handler_layout(tmp_path):
    root = tmp_path / "root"
    write_bam(root / "A" / "igv" / "igv_ready.bam", READS)
    write_bed(root / "A" / "targets_exceeding_threshold.bed", REGIONS)
    vcf = _clinvar(tmp_path / "clinvar.vcf")
    calls = {"results": [], "errors": []}
    context = SimpleNamespace(
        filepath=str(root / "A"),
        metadata={"clinvar_path": str(vcf)},
        get_sample_id=lambda: "A",
        add_result=lambda *args: calls["results"].append(args),
        add_error=lambda *args: calls["errors"].append(args),
    )
    job = SimpleNamespace(job_id=1, job_type="lightweight_gene_analysis", context=context)
    return root, job, calls


def test_the_handler_fails_the_job_and_still_writes_the_other_genes(tmp_path, monkeypatch):
    _fail_region(monkeypatch, "GENE2")
    root, job, calls = _handler_layout(tmp_path)
    lga.lightweight_gene_analysis_handler(job, work_dir=str(root))
    assert len(calls["errors"]) == 1
    job_type, text = calls["errors"][0]
    assert job_type == "lightweight_gene_analysis"
    assert "GENE2" in text and "synthetic region failure" in text
    data = json.loads((root / "A" / RESULTS).read_text())
    assert "GENE1" in data["genes"]
    assert [f["region"] for f in data["metadata"]["regions_failed"]] == ["GENE2"]


def test_the_handler_reports_no_error_for_a_clean_check(tmp_path):
    root, job, calls = _handler_layout(tmp_path)
    lga.lightweight_gene_analysis_handler(job, work_dir=str(root))
    assert calls["errors"] == []


def test_the_gui_status_names_the_regions_not_checked(tmp_path, monkeypatch):
    _fail_region(monkeypatch, "GENE2")
    _check(tmp_path)
    text = lga.lga_status_text(str(tmp_path / RESULTS))
    assert "GENE2" in text and "could not be checked" in text
    assert "successfully" not in text


def test_the_gui_status_for_a_clean_or_older_results_file(tmp_path):
    _check(tmp_path)
    assert "successfully" in lga.lga_status_text(str(tmp_path / RESULTS))
    old = tmp_path / "old.json"
    old.write_text(json.dumps({"metadata": {"analysis_type": "lightweight_gene_analysis"}, "genes": {}}))
    assert "successfully" in lga.lga_status_text(str(old))


def test_the_gui_uses_the_status_text():
    """The coverage view's results check sets its status from ``lga_status_text``."""
    source = Path(lga.__file__).resolve().parents[1] / "gui" / "components" / "coverage.py"
    tree = ast.parse(source.read_text())
    check = next(
        n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == "_check_lga_results"
    )
    called = {
        n.func.id if isinstance(n.func, ast.Name) else getattr(n.func, "attr", "")
        for n in ast.walk(check) if isinstance(n, ast.Call)
    }
    assert "lga_status_text" in called


def test_the_gui_run_button_keeps_the_jobs_errors():
    """Found after the lock: the GUI's own run of the check threw its job errors away
    (``add_error`` was a no-op) and then said "completed successfully!"."""
    source = Path(lga.__file__).resolve().parents[1] / "gui" / "components" / "coverage.py"
    tree = ast.parse(source.read_text())
    # The innermost function that builds the job and defines run_lga_analysis_direct
    def defines_run(n):
        return any(
            isinstance(c, ast.FunctionDef) and c.name == "run_lga_analysis_direct" for c in ast.walk(n)
            if c is not n
        )

    starter = min(
        (
            n for n in ast.walk(tree)
            if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef))
            and n.name != "run_lga_analysis_direct" and defines_run(n)
        ),
        key=lambda n: n.end_lineno - n.lineno,
    )
    error_hooks = [
        node.value for node in ast.walk(starter)
        if isinstance(node, ast.Assign)
        and any(isinstance(t, ast.Attribute) and t.attr == "add_error" for t in node.targets)
    ]
    assert error_hooks, "the GUI's run sets no add_error"
    for hook in error_hooks:
        assert not (
            isinstance(hook, ast.Lambda) and isinstance(hook.body, ast.Constant) and hook.body.value is None
        ), "a job error from the GUI's run is thrown away"
    run = next(
        n for n in ast.walk(starter)
        if isinstance(n, ast.FunctionDef) and n.name == "run_lga_analysis_direct"
    )
    assert "lga_errors" in {n.id for n in ast.walk(run) if isinstance(n, ast.Name)}
