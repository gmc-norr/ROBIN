"""VAR-9: the lightweight gene analysis uses only the sample's own files."""

import json
from types import SimpleNamespace

import pytest

from lga_fixtures import write_bam, write_bed, write_clinvar
from robin.analysis.lightweight_gene_analysis import (
    lga_results_belong_to_sample,
    lightweight_gene_analysis_handler,
    pick_lga_results,
    resolve_sample_folder,
)

RESULTS = "lightweight_gene_analysis_results.json"
REF_ONLY = [(0, "C" * 50)] * 20
ALT_AT_11 = [(0, "C" * 10 + "T" + "C" * 39)] * 20  # ClinVar POS 11 (1-based) is the T


def _layout(tmp_path, *, a_bam, a_bed):
    """Monitored root with sample A (reference-only) and sample B (20 alternate reads)."""
    root = tmp_path / "root"
    write_bam(root / "A" / "target.bam", REF_ONLY)
    if a_bam:
        write_bam(root / "A" / "igv" / "igv_ready.bam", REF_ONLY)
    if a_bed:
        write_bed(root / "A" / "targets_exceeding_threshold.bed", [(0, 50, "TEST")])
    write_bam(root / "B" / "igv" / "igv_ready.bam", ALT_AT_11)
    write_bed(root / "B" / "targets_exceeding_threshold.bed", [(0, 50, "TEST")])
    vcf = write_clinvar(tmp_path / "clinvar.vcf", [(11, "C", "T")])
    return root, vcf


def _job(vcf, filepath, sample_id="A", metadata=None):
    calls = {"results": [], "errors": []}
    context = SimpleNamespace(
        filepath=str(filepath),
        metadata={"clinvar_path": str(vcf), **(metadata or {})},
        get_sample_id=lambda: sample_id,
        add_result=lambda *args: calls["results"].append(args),
        add_error=lambda *args: calls["errors"].append(args),
    )
    return SimpleNamespace(job_id=1, job_type="lightweight_gene_analysis", context=context), calls


def test_missing_bam_does_not_use_another_samples_bam(tmp_path):
    root, vcf = _layout(tmp_path, a_bam=False, a_bed=True)
    job, calls = _job(vcf, root / "A")
    with pytest.raises(FileNotFoundError, match=r"A/igv/igv_ready\.bam"):
        lightweight_gene_analysis_handler(job, work_dir=str(root))
    assert list(root.rglob(RESULTS)) == []
    assert calls["errors"]


@pytest.mark.parametrize("give_bed_path", [False, True])
def test_missing_bed_does_not_use_another_samples_bed(tmp_path, give_bed_path):
    # The GUI passes A's (missing) BED path; other callers pass none.
    root, vcf = _layout(tmp_path, a_bam=True, a_bed=False)
    metadata = {"bed_path": str(root / "A" / "targets_exceeding_threshold.bed")} if give_bed_path else {}
    job, _ = _job(vcf, root / "A", metadata=metadata)
    with pytest.raises(FileNotFoundError, match=r"A/targets_exceeding_threshold\.bed"):
        lightweight_gene_analysis_handler(job, work_dir=str(root))
    assert list(root.rglob(RESULTS)) == []


def test_results_record_the_inputs_used(tmp_path):
    root, vcf = _layout(tmp_path, a_bam=True, a_bed=True)
    job, calls = _job(vcf, root / "A")
    lightweight_gene_analysis_handler(job, work_dir=str(root))
    meta = json.loads((root / "A" / RESULTS).read_text())["metadata"]
    assert meta["input_bam"] == str(root / "A" / "igv" / "igv_ready.bam")
    assert meta["input_bed"] == str(root / "A" / "targets_exceeding_threshold.bed")
    assert meta["sample_id"] == "A"
    result = calls["results"][0][1]
    assert result["input_bam"] == meta["input_bam"]
    assert result["input_bed"] == meta["input_bed"]


def test_sample_folder_passed_directly_still_works(tmp_path):
    # The GUI passes the sample folder itself when no monitored root is set.
    root, vcf = _layout(tmp_path, a_bam=True, a_bed=True)
    job, _ = _job(vcf, root / "A")
    lightweight_gene_analysis_handler(job, work_dir=str(root / "A"))
    # The output stays where the GUI looks first.
    meta = json.loads((root / "A" / "A" / RESULTS).read_text())["metadata"]
    assert meta["input_bam"] == str(root / "A" / "igv" / "igv_ready.bam")


def test_unknown_sample_id_looks_only_in_the_given_folder(tmp_path):
    root, vcf = _layout(tmp_path, a_bam=False, a_bed=True)
    job, _ = _job(vcf, root, sample_id="unknown")
    with pytest.raises(FileNotFoundError):
        lightweight_gene_analysis_handler(job, work_dir=str(root))
    assert list(root.rglob(RESULTS)) == []


@pytest.mark.parametrize("sample_id", ["../B", "A/../B", "/abs/B", "..", "."])
def test_sample_id_cannot_leave_the_base_folder(tmp_path, sample_id):
    with pytest.raises(ValueError):
        resolve_sample_folder(str(tmp_path / "root"), sample_id)


def test_trailing_slash_resolves_to_the_same_folder(tmp_path):
    folder = str(tmp_path / "A") + "/"
    assert resolve_sample_folder(folder, "A") == folder


@pytest.mark.parametrize("sample_id", [None, "", "unknown"])
def test_no_sample_id_means_the_given_folder(tmp_path, sample_id):
    assert resolve_sample_folder(str(tmp_path), sample_id) == str(tmp_path)


def test_sample_id_names_a_subfolder_of_the_root(tmp_path):
    assert resolve_sample_folder(str(tmp_path / "root"), "A") == str(tmp_path / "root" / "A")


def test_root_named_like_the_sample_resolves_to_the_sample_folder(tmp_path):
    root = tmp_path / "A"                      # monitored root happens to be called "A"
    (root / "A" / "igv").mkdir(parents=True)   # sample A lives at /x/A/A
    assert resolve_sample_folder(str(root), "A") == str(root / "A")


def test_output_folder_does_not_change_the_resolved_folder(tmp_path):
    sample = tmp_path / "root" / "A"
    (sample / "igv").mkdir(parents=True)
    (sample / "A").mkdir()                     # the handler's output folder in layout (b)
    assert resolve_sample_folder(str(sample), "A") == str(sample)


# Result files written before the fix can't prove which sample's BAM they used.


def _results_file(path, metadata):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"metadata": metadata, "genes": {}}))
    return path


def test_old_results_without_inputs_are_not_trusted(tmp_path):
    sample_a = tmp_path / "root" / "A"
    results = _results_file(sample_a / RESULTS, {"analysis_type": "lightweight_gene_analysis"})
    assert lga_results_belong_to_sample(str(results), str(sample_a)) is False


def test_results_from_another_samples_bam_are_not_trusted(tmp_path):
    root = tmp_path / "root"
    results = _results_file(root / "A" / RESULTS, {"input_bam": str(root / "B" / "igv" / "igv_ready.bam")})
    assert lga_results_belong_to_sample(str(results), str(root / "A")) is False


def test_results_from_a_sibling_with_a_longer_name_are_not_trusted(tmp_path):
    # "A2" starts with "A", but it is another sample.
    root = tmp_path / "root"
    results = _results_file(root / "A" / RESULTS, {"input_bam": str(root / "A2" / "igv" / "igv_ready.bam")})
    assert lga_results_belong_to_sample(str(results), str(root / "A")) is False


def test_results_from_this_samples_bam_are_trusted(tmp_path):
    root = tmp_path / "root"
    results = _results_file(root / "A" / "A" / RESULTS, {"input_bam": str(root / "A" / "igv" / "igv_ready.bam")})
    assert lga_results_belong_to_sample(str(results), str(root / "A")) is True


def test_unreadable_results_are_not_trusted(tmp_path):
    results = tmp_path / "A" / RESULTS
    results.parent.mkdir(parents=True)
    results.write_text("{not json")
    assert lga_results_belong_to_sample(str(results), str(tmp_path / "A")) is False


# The GUI must show the first results file it can trust, not the first that exists.


def test_pick_lga_results_skips_an_untrusted_candidate_for_a_trusted_one(tmp_path):
    root = tmp_path / "root"
    sample = root / "A"
    untrusted = _results_file(
        sample / "A" / RESULTS, {"input_bam": str(root / "B" / "igv" / "igv_ready.bam")}
    )
    trusted = _results_file(
        sample / RESULTS, {"input_bam": str(sample / "igv" / "igv_ready.bam")}
    )
    assert pick_lga_results([str(untrusted), str(trusted)], str(sample)) == (
        str(trusted),
        False,
    )


def test_pick_lga_results_reports_untrusted_when_no_candidate_is_trusted(tmp_path):
    sample = tmp_path / "root" / "A"
    untrusted = _results_file(
        sample / "A" / RESULTS, {"analysis_type": "lightweight_gene_analysis"}
    )
    assert pick_lga_results([str(untrusted), str(sample / RESULTS)], str(sample)) == (
        None,
        True,
    )


def test_pick_lga_results_finds_nothing_when_no_candidate_exists(tmp_path):
    sample = tmp_path / "root" / "A"
    candidates = [str(sample / "A" / RESULTS), str(sample / RESULTS)]
    assert pick_lga_results(candidates, str(sample)) == (None, False)


def test_pick_lga_results_keeps_the_first_of_two_trusted_candidates(tmp_path):
    sample = tmp_path / "root" / "A"
    bam = str(sample / "igv" / "igv_ready.bam")
    first = _results_file(sample / "A" / RESULTS, {"input_bam": bam})
    _results_file(sample / RESULTS, {"input_bam": bam})
    assert pick_lga_results([str(first), str(sample / RESULTS)], str(sample)) == (
        str(first),
        False,
    )
