"""An SNP run calls on every read the sample has, and earlier calls are reused only when
they came from those reads (batch 10b: JOB-6). Synthetic BAMs; Clair3 is never run.

- The SNP job first merges what is not yet in target.bam: new batch BAMs and staged files.
- A finished call records which batch BAMs its target.bam held, in
  ``clair3/target_bam_batches.txt`` (one name per line).
- Earlier calls are reused only when that record matches target.bam and nothing waits to
  be merged; otherwise the variants are called again.
"""

import os

import pytest

from target_fixtures import SAMPLE, bam_records, samtools_reference, write_bams, write_panel
from robin.analysis import target_analysis as ta_mod
from robin.analysis.target_analysis import (
    TargetAnalysis, finalize_accumulation_for_sample, run_snp_analysis, snp_analysis_handler,
)
from robin.gui_launcher import GUILauncher

PANEL_NAME = "TESTPANEL"
STAMP = "target_bam_batches.txt"


@pytest.fixture
def panel(tmp_path, monkeypatch):
    path = write_panel(str(tmp_path / "panel.bed"))
    monkeypatch.setattr(TargetAnalysis, "_find_target_bed", lambda self, p: path)
    return path


class _Context:
    def __init__(self, sample_dir, metadata):
        self.filepath = str(sample_dir)
        self.metadata = metadata
        self.errors, self.results = {}, {}

    def get_sample_id(self):
        return SAMPLE

    def add_error(self, key, value):
        self.errors[key] = value

    def add_result(self, key, value):
        self.results[key] = value


class _Job:
    job_id = 1
    job_type = "snp_analysis"

    def __init__(self, context):
        self.context = context


@pytest.fixture
def sample(tmp_path, panel):
    """A sample finalized after its first BAM, with two more BAMs on hand."""
    bams, merged = write_bams(str(tmp_path), 3)
    work = tmp_path / "work"
    ta = TargetAnalysis(work_dir=str(work), target_panel=PANEL_NAME, batch_size=1, use_staging=True)
    ta.process_file_with_staging(bams[0], {"sample_id": SAMPLE})
    ta.accumulate_staged_files(SAMPLE, force=True)
    finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
    sd = work / SAMPLE
    (sd / "targets_exceeding_threshold.bed").write_text("chr1\t100\t400\tGENE_A\n")
    return {"work": work, "dir": sd, "ta": ta, "bams": bams, "merged": merged, "panel": panel}


def _batch(sample, i):
    sample["ta"].process_file_with_staging(sample["bams"][i], {"sample_id": SAMPLE})
    sample["ta"].accumulate_staged_files(SAMPLE, force=True)


def _stage(sample, i):
    sample["ta"].process_file_with_staging(sample["bams"][i], {"sample_id": SAMPLE})


def _record_snp_runs(monkeypatch, sample):
    """Stand in for run_snp_analysis: note target.bam's reads when it is called."""
    seen = []

    def fake_run(sample_dir, **kwargs):
        seen.append(bam_records(os.path.join(sample_dir, "target.bam")))
        clair = os.path.join(sample_dir, "clair3")
        os.makedirs(clair, exist_ok=True)
        return clair

    monkeypatch.setattr(ta_mod, "run_snp_analysis", fake_run)
    return seen


def _run_handler(sample, **metadata):
    meta = {"work_dir": str(sample["work"]), "target_panel": PANEL_NAME, "reference": "/ref.fa"}
    meta.update(metadata)
    job = _Job(_Context(sample["dir"], meta))
    snp_analysis_handler(job)
    return job


def _batch_names(sample):
    return ta_mod._merged_batch_names(str(sample["dir"] / "target.bam"))


def _write_calls(sample, stamp_names=None):
    clair = sample["dir"] / "clair3"
    clair.mkdir(exist_ok=True)
    for name in ("snpsift_output.vcf", "snpsift_indel_output.vcf"):
        (clair / name).write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    if stamp_names is not None:
        (clair / STAMP).write_text("".join(f"{n}\n" for n in sorted(stamp_names)))


def _run_until_input_checks(sample):
    """run_snp_analysis on the sample. Without a targets BED, a run that does not reuse
    earlier calls stops at its input checks and returns "": Clair3 is never reached."""
    os.remove(sample["dir"] / "targets_exceeding_threshold.bed")
    return run_snp_analysis(str(sample["dir"]))


def _stamp(sample):
    path = sample["dir"] / "clair3" / STAMP
    if not path.exists():
        return None
    return {line.strip() for line in path.read_text().splitlines() if line.strip()}


# The SNP job merges first -------------------------------------------------------------


def test_an_snp_run_first_merges_reads_not_yet_in_target_bam(sample, monkeypatch):
    _batch(sample, 1)  # a batch BAM not yet in target.bam
    _stage(sample, 2)  # a file still staged
    seen = _record_snp_runs(monkeypatch, sample)
    job = _run_handler(sample)
    assert job.context.errors == {}
    assert len(seen) == 1
    assert seen[0] == samtools_reference(sample["merged"], sample["panel"])["target_bam"]


def test_an_snp_run_fails_when_its_merge_fails(sample, monkeypatch):
    _batch(sample, 1)
    seen = _record_snp_runs(monkeypatch, sample)
    monkeypatch.setattr(ta_mod, "finalize_accumulation_for_sample",
                        lambda *a, **k: {"status": "error", "error": "disk full", "sample_id": SAMPLE})
    job = _run_handler(sample)
    assert seen == [], "variants were called on a target.bam missing reads"
    assert "disk full" in job.context.errors.get("snp_analysis", "")


def test_a_finished_call_records_the_batches_it_read(sample, monkeypatch):
    _batch(sample, 1)
    _record_snp_runs(monkeypatch, sample)
    job = _run_handler(sample)
    assert job.context.errors == {}
    assert _stamp(sample) == _batch_names(sample)
    assert len(_stamp(sample)) == 2


def test_an_annotation_only_run_records_nothing(sample, monkeypatch):
    _record_snp_runs(monkeypatch, sample)
    _run_handler(sample, annotation_only=True)
    assert _stamp(sample) is None


# Earlier calls are reused only when they came from the same reads -----------------------


def test_calls_from_the_current_target_bam_are_reused(sample):
    _write_calls(sample, _batch_names(sample))
    assert _run_until_input_checks(sample) == str(sample["dir"] / "clair3")


def test_calls_from_an_older_target_bam_are_not_reused(sample):
    _write_calls(sample, _batch_names(sample))
    _batch(sample, 1)
    finalize_accumulation_for_sample(SAMPLE, str(sample["work"]), PANEL_NAME)
    assert _run_until_input_checks(sample) == ""
    assert _stamp(sample) is None, "a call that may be redone kept its old record"


def test_calls_without_a_record_are_not_reused(sample):
    _write_calls(sample, None)
    assert _run_until_input_checks(sample) == ""


@pytest.mark.parametrize("waiting", ["batch", "staged"])
def test_calls_are_not_reused_while_reads_wait_to_be_merged(sample, waiting):
    _write_calls(sample, _batch_names(sample))
    (_batch if waiting == "batch" else _stage)(sample, 1)
    assert _run_until_input_checks(sample) == ""


# The GUI's bulk SNP list --------------------------------------------------------------


def _gui(sample):
    g = GUILauncher()
    g.monitored_directory = str(sample["work"])
    return g


def test_bulk_snp_lists_a_sample_whose_calls_are_older_than_target_bam(sample):
    _write_calls(sample, set())  # made before the first batch was merged
    assert _gui(sample)._list_samples_needing_snp_calling() == [SAMPLE]


def test_bulk_snp_leaves_out_a_sample_whose_calls_are_current(sample):
    _write_calls(sample, _batch_names(sample))
    assert _gui(sample)._list_samples_needing_snp_calling() == []
