"""target.bam and the files made from it stay in step with the coverage (final review of
batch 6: Important 1 and 2, Minor 4). Synthetic BAMs only.
"""

import os
import shutil
import time

import pandas as pd
import pysam
import pytest

from target_fixtures import SAMPLE, bam_records, samtools_reference, write_bams, write_panel
from robin.analysis import target_analysis as ta_mod
from robin.analysis.target_analysis import (
    FileLock, TargetAnalysis, ensure_sorted_igv_bam, finalize_accumulation_for_sample,
)

PANEL_NAME = "TESTPANEL"


@pytest.fixture
def panel(tmp_path, monkeypatch):
    path = write_panel(str(tmp_path / "panel.bed"))
    monkeypatch.setattr(TargetAnalysis, "_find_target_bed", lambda self, p: path)
    return path


def _analysis(work):
    return TargetAnalysis(work_dir=str(work), target_panel=PANEL_NAME, batch_size=1, use_staging=True)


def _add(ta, bam):
    ta.process_file_with_staging(bam, {"sample_id": SAMPLE})
    return ta.accumulate_staged_files(SAMPLE, force=True)


def test_a_batch_whose_commit_fails_reaches_neither_the_coverage_nor_target_bam(tmp_path, panel, monkeypatch):
    """Important 1: a batch BAM is visible to finalize only once its coverage is committed,
    so a commit that fails (a lock timeout, an unreadable totals file) leaves the coverage
    and target.bam agreeing."""
    bams, _ = write_bams(str(tmp_path), 2)
    work = tmp_path / "work"
    ta = _analysis(work)
    assert _add(ta, bams[0])["status"] == "success"

    def failing_commit(self, *args, **kwargs):
        raise TimeoutError("lock held")

    monkeypatch.setattr(TargetAnalysis, "_commit_batch", failing_commit)
    assert _add(ta, bams[1])["status"] == "lock_timeout"
    monkeypatch.undo()
    finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)

    sd = work / SAMPLE
    only_first = samtools_reference(bams[0], panel)
    assert bam_records(str(sd / "target.bam")) == only_first["target_bam"]
    cov = pd.read_csv(sd / "coverage_main.csv").set_index("#rname")
    for name, (numreads, _, _) in only_first["contigs"].items():
        assert int(cov.loc[name, "numreads"]) == numreads


def test_finalize_does_not_hold_the_accumulation_lock_while_merging(tmp_path, panel, monkeypatch):
    """Important 1: a long merge must not make a worker's commit time out."""
    bams, _ = write_bams(str(tmp_path), 1)
    work = tmp_path / "work"
    ta = _analysis(work)
    _add(ta, bams[0])
    lock = ta._get_lock_file(SAMPLE, "accumulation")
    real_merge = ta_mod._merge_into_target_bam
    seen = []

    def merge_and_try_the_lock(*args, **kwargs):
        try:
            with FileLock(lock, timeout=0.5):
                seen.append("free")
        except TimeoutError:
            seen.append("held")
        return real_merge(*args, **kwargs)

    monkeypatch.setattr(ta_mod, "_merge_into_target_bam", merge_and_try_the_lock)
    finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
    assert seen == ["free"]


def _two_finalizes(tmp_path, panel):
    """Finalize after the first BAM, build the IGV BAM, then add the second BAM and finalize."""
    bams, merged = write_bams(str(tmp_path), 2)
    work = tmp_path / "work"
    sd = work / SAMPLE
    ta = _analysis(work)
    _add(ta, bams[0])
    finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
    igv = ensure_sorted_igv_bam(str(sd))
    assert igv and bam_records(igv) == bam_records(str(sd / "target.bam"))
    old = tmp_path / "igv_after_first"
    old.mkdir()
    shutil.copy2(igv, old / "igv_ready.bam")
    shutil.copy2(f"{igv}.bai", old / "igv_ready.bam.bai")
    _add(ta, bams[1])
    finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
    return sd, merged, old


def test_finalize_keeps_the_igv_bam_in_step_with_target_bam(tmp_path, panel):
    """Important 2: IGV and the per-gene variant check read igv_ready.bam; after a finalize
    it holds the same records as target.bam."""
    sd, merged, _ = _two_finalizes(tmp_path, panel)
    want = samtools_reference(merged, panel)["target_bam"]
    assert bam_records(str(sd / "target.bam")) == want
    assert bam_records(str(sd / "igv" / "igv_ready.bam")) == want


def test_a_stale_igv_bam_is_rebuilt_and_a_fresh_one_is_kept(tmp_path, panel):
    """Important 2: ensure_sorted_igv_bam rebuilds an IGV BAM (or a legacy Clair3 BAM) made
    from an older target.bam, and leaves one made from the current target.bam alone."""
    sd, merged, old = _two_finalizes(tmp_path, panel)
    want = bam_records(str(sd / "target.bam"))
    igv = sd / "igv" / "igv_ready.bam"

    # An IGV BAM from the first target.bam (a refresh that failed, or an older folder).
    shutil.copy2(old / "igv_ready.bam", igv)
    shutil.copy2(old / "igv_ready.bam.bai", f"{igv}.bai")
    assert bam_records(str(igv)) != want
    assert bam_records(ensure_sorted_igv_bam(str(sd))) == want

    # A fresh one is kept as it is.
    before = os.stat(igv).st_mtime_ns
    assert ensure_sorted_igv_bam(str(sd)) == str(igv)
    assert os.stat(igv).st_mtime_ns == before

    # No IGV BAM, and a legacy Clair3 BAM sorted from the first target.bam: not reused.
    os.remove(igv)
    os.remove(f"{igv}.bai")
    clair = sd / "clair3"
    clair.mkdir()
    shutil.copy2(old / "igv_ready.bam", clair / "sorted_targets_exceeding.bam")
    shutil.copy2(old / "igv_ready.bam.bai", clair / "sorted_targets_exceeding.bam.bai")
    assert bam_records(ensure_sorted_igv_bam(str(sd))) == want


def test_an_igv_bam_sorted_from_an_older_target_bam_is_rebuilt_even_if_newer(tmp_path, panel):
    """Important 2, the race: a sort that read the older target.bam finishes after the
    newer one was written, so its file is newer but its records are old. The merged-batch
    list in its header tells."""
    sd, _, old = _two_finalizes(tmp_path, panel)
    want = bam_records(str(sd / "target.bam"))
    igv = sd / "igv" / "igv_ready.bam"
    shutil.copy(old / "igv_ready.bam", igv)  # copy, not copy2: the file is new
    shutil.copy(old / "igv_ready.bam.bai", f"{igv}.bai")
    assert os.stat(igv).st_mtime_ns > os.stat(sd / "target.bam").st_mtime_ns
    assert bam_records(ensure_sorted_igv_bam(str(sd))) == want


def test_an_old_folder_igv_bam_follows_target_bam_by_time(tmp_path, panel):
    """A folder from before the merged-batch list: neither BAM has one, so an IGV BAM older
    than target.bam is rebuilt and a newer one is kept."""
    bams, merged = write_bams(str(tmp_path), 2)
    sd = tmp_path / "work" / SAMPLE
    (sd / "igv").mkdir(parents=True)
    igv = sd / "igv" / "igv_ready.bam"
    shutil.copy2(bams[0], igv)
    shutil.copy2(f"{bams[0]}.bai", f"{igv}.bai")
    shutil.copy(merged, sd / "target.bam")
    pysam.index(str(sd / "target.bam"))
    now = time.time_ns()
    os.utime(igv, ns=(now - 100 * 10**9,) * 2)  # the IGV BAM, then target.bam, in the past
    os.utime(sd / "target.bam", ns=(now - 50 * 10**9,) * 2)
    assert bam_records(ensure_sorted_igv_bam(str(sd))) == bam_records(merged)
    before = os.stat(igv).st_mtime_ns
    ensure_sorted_igv_bam(str(sd))
    assert os.stat(igv).st_mtime_ns == before


def test_finalize_rebuilds_a_stale_target_bam_index(tmp_path, panel):
    """Minor 4: if replacing the index failed after target.bam was replaced, the old index
    stays and region fetches go wrong; the next finalize rebuilds it."""
    bams, _ = write_bams(str(tmp_path), 2)
    work = tmp_path / "work"
    sd = work / SAMPLE
    ta = _analysis(work)
    _add(ta, bams[0])
    finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
    old_bai = tmp_path / "old.bai"
    shutil.copy2(sd / "target.bam.bai", old_bai)
    _add(ta, bams[1])
    finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)

    shutil.copy2(old_bai, sd / "target.bam.bai")  # copy2 keeps the old index's older mtime
    assert os.stat(sd / "target.bam.bai").st_mtime_ns < os.stat(sd / "target.bam").st_mtime_ns
    finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)  # nothing new to merge

    fresh = tmp_path / "fresh.bai"
    pysam.index(str(sd / "target.bam"), str(fresh))
    assert (sd / "target.bam.bai").read_bytes() == fresh.read_bytes()
