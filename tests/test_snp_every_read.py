"""Batch 10b, round 2: fixes from the external review. Synthetic BAMs; Clair3 and the
annotation tools are never run (a stand-in container writes the VCFs a run would).

1. Reads a worker is still adding (moved into ``_staging/_batch_<id>/``) count as waiting:
   finalize waits for that worker, and a claim folder left behind counts as reads not in
   target.bam.
2. An SNP run goes ahead only when its merge got every read in: an accumulation that
   timed out, could not read its files or could not write its batch BAM fails the run.
3. A Clair3 run counts only when every region ran and wrote both VCFs during this run;
   a failed region, a missing VCF, or a file left from an earlier run fails it.
4. Calls are current only for target.bam as it is now: a target.bam rewritten without a
   new batch is not the one they came from.
"""

import gzip
import os
import shutil
import threading
import time
from types import SimpleNamespace

import pytest

from target_fixtures import SAMPLE, bam_records, samtools_reference, write_bams, write_panel
from robin.analysis import target_analysis as ta_mod
from robin.analysis.target_analysis import (
    TargetAnalysis, finalize_accumulation_for_sample, run_snp_analysis, snp_analysis_handler,
)
from robin.gui_launcher import GUILauncher

PANEL_NAME = "TESTPANEL"


@pytest.fixture
def panel(tmp_path, monkeypatch):
    path = write_panel(str(tmp_path / "panel.bed"))
    monkeypatch.setattr(TargetAnalysis, "_find_target_bed", lambda self, p: path)
    return path


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


def _stage(sample, i):
    sample["ta"].process_file_with_staging(sample["bams"][i], {"sample_id": SAMPLE})


def _claim_folders(sample):
    staging = sample["dir"] / "_staging"
    return sorted(p.name for p in staging.glob("_batch_*")) if staging.exists() else []


# The SNP job, with run_snp_analysis replaced -------------------------------------------


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


def _record_snp_runs(monkeypatch):
    """Stand in for run_snp_analysis: note target.bam's reads when it is called."""
    seen = []

    def fake_run(sample_dir, **kwargs):
        seen.append(bam_records(os.path.join(sample_dir, "target.bam")))
        clair = os.path.join(sample_dir, "clair3")
        os.makedirs(clair, exist_ok=True)
        return clair

    monkeypatch.setattr(ta_mod, "run_snp_analysis", fake_run)
    return seen


def _run_handler(sample):
    meta = {"work_dir": str(sample["work"]), "target_panel": PANEL_NAME, "reference": "/ref.fa"}
    job = _Job(_Context(sample["dir"], meta))
    snp_analysis_handler(job)
    return job


def _current_calls(sample, monkeypatch):
    """Calls made from target.bam as it is now: the SNP job writes their record."""
    clair = sample["dir"] / "clair3"
    clair.mkdir(exist_ok=True)
    for name in ("snpsift_output.vcf", "snpsift_indel_output.vcf"):
        (clair / name).write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    with monkeypatch.context() as m:
        _record_snp_runs(m)
        assert _run_handler(sample).context.errors == {}
    assert ta_mod.snp_calls_are_current(str(sample["dir"]))


class _PausedWorker:
    """A worker that has claimed staged files and is still adding them: it stops inside
    its batch BAM write until released."""

    def __init__(self, sample, monkeypatch):
        self.started, self.release = threading.Event(), threading.Event()
        orig, calls = TargetAnalysis._write_batch_bam, []

        def paused(ta_self, *a, **k):
            if not calls:
                calls.append(1)
                self.started.set()
                self.release.wait(30)
            return orig(ta_self, *a, **k)

        monkeypatch.setattr(TargetAnalysis, "_write_batch_bam", paused)
        self.thread = threading.Thread(
            target=lambda: sample["ta"].accumulate_staged_files(SAMPLE, force=True), daemon=True
        )

    def __enter__(self):
        self.thread.start()
        assert self.started.wait(10), "the worker never claimed the staged files"
        return self

    def __exit__(self, *exc):
        self.release.set()
        self.thread.join(30)


def _left_by_a_stopped_worker(sample, i):
    """Staged files a worker claimed, then stopped (crashed) before adding."""
    _stage(sample, i)
    staging = sample["dir"] / "_staging"
    dead = staging / "_batch_dead"
    dead.mkdir()
    for f in list(staging.iterdir()):
        if f.is_file():
            shutil.move(str(f), str(dead / f.name))
    return dead


# 1. Reads a worker is still adding ------------------------------------------------------


def test_finalize_waits_for_a_worker_still_adding_reads(sample, monkeypatch):
    _stage(sample, 1)
    _stage(sample, 2)
    with _PausedWorker(sample, monkeypatch) as worker:
        threading.Timer(0.5, worker.release.set).start()
        result = finalize_accumulation_for_sample(SAMPLE, str(sample["work"]), PANEL_NAME)
    assert bam_records(str(sample["dir"] / "target.bam")) == samtools_reference(
        sample["merged"], sample["panel"]
    )["target_bam"], "finalize merged before the worker had added its reads"
    assert not result.get("claims_left")


def test_calls_are_not_current_while_a_worker_adds_reads(sample, monkeypatch):
    _current_calls(sample, monkeypatch)
    _stage(sample, 1)
    with _PausedWorker(sample, monkeypatch):
        assert not ta_mod.snp_calls_are_current(str(sample["dir"]))


def test_the_gui_does_not_call_finalize_redundant_while_a_worker_adds_reads(sample, monkeypatch):
    _stage(sample, 1)
    g = GUILauncher()
    g.monitored_directory = str(sample["work"])
    with _PausedWorker(sample, monkeypatch):
        assert not g._is_target_bam_finalize_redundant(SAMPLE)


def test_an_snp_run_fails_while_a_worker_still_adds_reads_after_the_wait(sample, monkeypatch):
    monkeypatch.setattr(ta_mod, "_CLAIM_WAIT_SECONDS", 0.3, raising=False)
    _stage(sample, 1)
    seen = _record_snp_runs(monkeypatch)
    with _PausedWorker(sample, monkeypatch):
        job = _run_handler(sample)
    assert seen == [], "variants were called while a worker was still adding reads"
    assert "_batch_" in job.context.errors.get("snp_analysis", "")


def test_finalize_reports_a_claim_folder_a_stopped_worker_left(sample):
    dead = _left_by_a_stopped_worker(sample, 1)
    result = finalize_accumulation_for_sample(SAMPLE, str(sample["work"]), PANEL_NAME)
    assert result.get("claims_left") == ["_batch_dead"]
    assert dead.exists(), "the reads were deleted"


def test_an_snp_run_fails_on_a_claim_folder_a_stopped_worker_left(sample, monkeypatch):
    _left_by_a_stopped_worker(sample, 1)
    seen = _record_snp_runs(monkeypatch)
    job = _run_handler(sample)
    assert seen == []
    assert "_batch_dead" in job.context.errors.get("snp_analysis", "")
    assert not ta_mod.snp_calls_are_current(str(sample["dir"]))


def test_a_finished_claim_leaves_no_folder_and_no_lock(sample):
    _stage(sample, 1)
    sample["ta"].accumulate_staged_files(SAMPLE, force=True)
    assert _claim_folders(sample) == []
    assert list((sample["dir"] / "_locks").glob("claim*")) == []


# 2. An SNP run goes ahead only when every read got in ----------------------------------


def _batch_bam_fails(monkeypatch):
    def fail(self, *a, **k):
        raise OSError("disk full")

    monkeypatch.setattr(TargetAnalysis, "_write_batch_bam", fail)


def test_a_batch_bam_that_cannot_be_written_keeps_its_staged_files(sample, monkeypatch):
    _batch_bam_fails(monkeypatch)
    _stage(sample, 1)
    result = sample["ta"].accumulate_staged_files(SAMPLE, force=True)
    assert "disk full" in result.get("batch_bam_error", "")
    folders = _claim_folders(sample)
    assert len(folders) == 1, "the reads that never reached target.bam were deleted"
    assert list((sample["dir"] / "_staging" / folders[0]).glob("coverage_*.parquet"))


def test_a_staged_file_that_cannot_be_read_keeps_its_folder(sample):
    _stage(sample, 1)
    _stage(sample, 2)
    first = sorted((sample["dir"] / "_staging").glob("coverage_*.parquet"))[0]
    first.write_bytes(b"not a parquet file")
    sample["ta"].accumulate_staged_files(SAMPLE, force=True)
    folders = _claim_folders(sample)
    assert len(folders) == 1, "a staged file that could not be read was deleted"
    assert (sample["dir"] / "_staging" / folders[0] / first.name).exists()


def test_an_snp_run_fails_again_after_a_batch_bam_could_not_be_written(sample, monkeypatch):
    _batch_bam_fails(monkeypatch)
    _stage(sample, 1)
    seen = _record_snp_runs(monkeypatch)
    first, second = _run_handler(sample), _run_handler(sample)
    assert seen == [], "variants were called on a target.bam missing reads"
    assert first.context.errors.get("snp_analysis")
    assert second.context.errors.get("snp_analysis"), "the second run called the lost reads complete"


@pytest.mark.parametrize("merge", [
    {"status": "lock_timeout", "error": "Could not acquire lock"},
    {"status": "load_failed", "files_attempted": 1},
    {"status": "success", "files_processed": 1, "batch_bam_error": "disk full"},
], ids=["lock_timeout", "load_failed", "batch_bam_error"])
def test_an_snp_run_fails_when_its_merge_missed_reads(sample, monkeypatch, merge):
    seen = _record_snp_runs(monkeypatch)
    monkeypatch.setattr(ta_mod, "finalize_accumulation_for_sample", lambda *a, **k: dict(merge))
    job = _run_handler(sample)
    assert seen == []
    assert job.context.errors.get("snp_analysis")


# 3. A Clair3 run counts only when every region ran and wrote both VCFs now -----------------


_HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _vcf_gz(path, ids, old=False):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with gzip.open(path, "wt") as f:
        f.write(_HEADER)
        for k, vid in enumerate(ids):
            f.write(f"chr1\t{100 + k}\t{vid}\tA\tG\t50\tPASS\t.\n")
    if old:
        an_hour_ago = time.time() - 3600
        os.utime(path, (an_hour_ago, an_hour_ago))


class _FakeDocker:
    """Stands in for the docker module. Each started container acts as the next region's
    run: its exit code, and which of snv.vcf.gz / indel.vcf.gz it writes."""

    def __init__(self, clair_dir, regions):
        self.clair_dir, self.regions, self.started = clair_dir, list(regions), 0
        self.api = self
        self.images = SimpleNamespace(get=lambda name: True, pull=lambda name: True)

    def from_env(self):
        return self

    def ping(self):
        return True

    def create_host_config(self, binds):
        return {}

    def create_container(self, image, command, volumes, host_config):
        out = command.split("--output_dir ")[1].split()[0]
        self._out = out.replace("/data/output", self.clair_dir)
        return {"Id": str(self.started)}

    def start(self, container):
        r = self.regions[self.started]
        self.started += 1
        if r.get("snv", True):
            _vcf_gz(os.path.join(self._out, "snv.vcf.gz"), [f"NEW_SNV_{self.started}"])
        if r.get("indel", True):
            _vcf_gz(os.path.join(self._out, "indel.vcf.gz"), [f"NEW_INDEL_{self.started}"])
        self._exit = r.get("exit", 0)

    def logs(self, container, stream=False, follow=False, stdout=True, stderr=True):
        return [] if stream else b""

    def wait(self, container):
        return {"StatusCode": self._exit}

    def inspect_container(self, container):
        return {"State": {}}

    def remove_container(self, container):
        pass


@pytest.fixture
def caller(sample, tmp_path, monkeypatch):
    """run_snp_analysis on the sample with a stand-in Clair3 container and stand-in
    snpEff/SnpSift that pass the calls through unchanged."""
    sd = sample["dir"]
    shutil.rmtree(sd / "bed_files", ignore_errors=True)  # use the BED below
    ref = tmp_path / "ref" / "ref.fa"
    ref.parent.mkdir()
    ref.write_text(">chr1\nACGT\n")

    def fake_subprocess_run(cmd, stdout=None, **kwargs):
        if cmd[:2] == ["snpEff", "-version"]:
            return SimpleNamespace(returncode=0, stdout="stub", stderr="")
        if not os.path.exists(cmd[-1]):
            return SimpleNamespace(returncode=1, stderr="no input")
        with gzip.open(cmd[-1], "rt") as f:
            stdout.write(f.read())
        return SimpleNamespace(returncode=0, stderr="")

    def fake_clinvar(src, dst, db, logger, **kwargs):
        shutil.copy(src, dst)
        return "stub"

    monkeypatch.setattr(ta_mod.subprocess, "run", fake_subprocess_run)
    monkeypatch.setattr(ta_mod, "annotate_with_clinvar", fake_clinvar)
    monkeypatch.setattr(ta_mod, "_resolve_clinvar_db_for_snpsift", lambda logger: None)
    monkeypatch.delenv("ROBIN_CLAIRS_SPLIT_REGIONS", raising=False)

    def run(regions, bed="chr1\t100\t400\tGENE_A\n", split=True, docker_cls=None):
        (sd / "targets_exceeding_threshold.bed").write_text(bed)
        monkeypatch.setenv("ROBIN_CLAIRS_SPLIT_REGIONS", "1" if split else "0")
        docker = (docker_cls or _FakeDocker)(str(sd / "clair3"), regions)
        monkeypatch.setattr(ta_mod, "docker", docker)
        out = run_snp_analysis(str(sd), threads=1, force_regenerate=True, reference=str(ref))
        return out, docker

    return run


_TWO_CHUNKS = "chr1\t100\t200000000\tA\nchr2\t100\t100000000\tB\n"


def _ids(path):
    with open(path) as f:
        return {line.split("\t")[2] for line in f if not line.startswith("#")}


def test_a_complete_run_keeps_every_regions_calls(caller, sample):
    out, docker = caller([{}, {}], bed=_TWO_CHUNKS)
    assert docker.started == 2
    assert out == str(sample["dir"] / "clair3")
    assert _ids(sample["dir"] / "clair3" / "snpsift_output.vcf") == {"NEW_SNV_1", "NEW_SNV_2"}
    assert _ids(sample["dir"] / "clair3" / "snpsift_indel_output.vcf") == {"NEW_INDEL_1", "NEW_INDEL_2"}


def test_a_failed_region_fails_the_run(caller):
    out, docker = caller([{}, {"exit": 1, "snv": False, "indel": False}], bed=_TWO_CHUNKS)
    assert docker.started == 2
    assert out == "", "a region that failed was left out of calls shown as complete"


def test_old_indel_calls_are_not_kept_when_the_caller_writes_none(caller, sample):
    clair = sample["dir"] / "clair3"
    _vcf_gz(str(clair / "region_1" / "indel.vcf.gz"), ["OLD_INDEL"], old=True)
    _vcf_gz(str(clair / "output_indel_done.vcf.gz"), ["OLD_INDEL"], old=True)
    out, _ = caller([{"indel": False}])
    assert out == "", "an INDEL file from an earlier run was shown as this run's"


def test_an_old_snv_file_is_not_taken_as_new_after_a_failed_single_pass(caller, sample):
    _vcf_gz(str(sample["dir"] / "clair3" / "snv.vcf.gz"), ["OLD_SNV"], old=True)
    _vcf_gz(str(sample["dir"] / "clair3" / "indel.vcf.gz"), ["OLD_INDEL"], old=True)
    out, _ = caller([{"exit": 1, "snv": False, "indel": False}], split=False)
    assert out == "", "an SNV file from an earlier run was shown as this run's"


class _UnreadableSnv(_FakeDocker):
    """Exits 0 but writes an SNV file that is not gzip: the VCF merge fails."""

    def start(self, container):
        super().start(container)
        with open(os.path.join(self._out, "snv.vcf.gz"), "wb") as f:
            f.write(b"not gzip")


def test_a_failed_vcf_merge_fails_the_run(caller, sample):
    _vcf_gz(str(sample["dir"] / "clair3" / "merged_snv.vcf.gz"), ["OLD_SNV"], old=True)
    out, _ = caller([{}], docker_cls=_UnreadableSnv)
    assert out == "", "an earlier run's merged SNV file was shown as this run's"


# 4. Calls are current only for target.bam as it is now ----------------------------------


def test_calls_are_not_reused_after_target_bam_changes_without_a_new_batch(sample, monkeypatch):
    import pysam

    _current_calls(sample, monkeypatch)
    tb = str(sample["dir"] / "target.bam")
    tmp = tb + ".half.bam"
    with pysam.AlignmentFile(tb, "rb") as src, pysam.AlignmentFile(tmp, "wb", header=src.header) as dst:
        for k, rec in enumerate(src):
            if k % 2 == 0:
                dst.write(rec)
    os.replace(tmp, tb)
    pysam.index(tb)
    assert ta_mod._merged_batch_names(tb), "the batch list must be unchanged for this case"
    assert not ta_mod.snp_calls_are_current(str(sample["dir"])), (
        "calls from the old target.bam count as current"
    )


# Round 3 (second external review) ------------------------------------------------------


def test_a_quick_retry_does_not_take_the_previous_runs_indel_file(caller):
    first, _ = caller([{}])
    assert first, "the first run should succeed"
    out, _ = caller([{"indel": False}])
    assert out == "", "the previous run's INDEL file was shown as this run's"


def test_old_run_folders_are_removed(caller, sample):
    caller([{}])
    caller([{}])
    assert len(list((sample["dir"] / "clair3").glob("run_*"))) == 1


def _drop_a_part(sample):
    """Remove a part of the newest staged set; its coverage file (written last) stays."""
    staging = sample["dir"] / "_staging"
    cov = sorted(staging.glob("coverage_*.parquet"))[-1]
    n = cov.name[len("coverage_"):-len(".parquet")]
    (staging / f"source_bam_{n}.txt").unlink()
    return cov.name


def test_finalize_reports_a_staged_set_missing_a_part(sample):
    _stage(sample, 1)
    name = _drop_a_part(sample)
    result = finalize_accumulation_for_sample(SAMPLE, str(sample["work"]), PANEL_NAME)
    assert result.get("incomplete_staged") == [name]
    assert (sample["dir"] / "_staging" / name).exists(), "the broken set was deleted"


@pytest.mark.parametrize("with_complete_set", [False, True], ids=["alone", "beside_a_complete_set"])
def test_an_snp_run_fails_on_a_staged_set_missing_a_part(sample, monkeypatch, with_complete_set):
    _stage(sample, 1)
    name = _drop_a_part(sample)
    if with_complete_set:
        _stage(sample, 2)
    seen = _record_snp_runs(monkeypatch)
    job = _run_handler(sample)
    assert seen == [], "variants were called without the broken set's reads"
    assert name in job.context.errors.get("snp_analysis", "")
