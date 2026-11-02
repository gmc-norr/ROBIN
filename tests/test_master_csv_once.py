"""A sample's master.csv counts each BAM once (batch 10d-1: JOB-2, the preprocessing step).

Synthetic BAMs (``target_fixtures``) through the real preprocessing handler.

- A BAM preprocessed again (a restart, the scan and the watcher both seeing it) is not
  counted again, and neither is a copy of it at another path: a BAM is known by its bytes.
- Different BAMs are each counted (control).
- master.csv that can't be read is kept and the failure is reported, never reset to zero.
- A failed master.csv update is reported, not only logged.
"""

import csv
import os
import shutil
from itertools import count

import pytest

pytest.importorskip("pysam")

from target_fixtures import fixture_reads, write_bam  # noqa: E402
from robin import workflow_ray as wr  # noqa: E402
from robin.analysis import bam_preprocessor  # noqa: E402
from robin.analysis.master_csv_manager import MasterCSVManager  # noqa: E402

_ids = count(1)
SAMPLE = "S"  # target_fixtures' read group LB
COUNTED = [  # every column a BAM adds to
    "counter_bam_passed", "counter_bam_failed",
    "counter_mapped_count", "counter_unmapped_count", "counter_bases_count",
    "counter_mapped_reads_num", "counter_unmapped_reads_num",
    "counter_mapped_bases", "counter_unmapped_bases",
    "counter_pass_mapped_count", "counter_pass_unmapped_count", "counter_pass_bases_count",
    "counter_pass_mapped_reads_num", "counter_pass_unmapped_reads_num",
    "counter_pass_mapped_bases", "counter_pass_unmapped_bases",
    "bam_tracking_counter", "bam_tracking_total_files",
]


def _bam(directory, name, part):
    """A pass BAM with every other fixture read (part 0 or 1), or all of them (None)."""
    reads = fixture_reads()
    os.makedirs(directory, exist_ok=True)
    return write_bam(os.path.join(str(directory), name), reads if part is None else reads[part::2])


def _run(bam, work):
    ctx = wr.WorkflowContext(filepath=str(bam), metadata={"work_dir": str(work), "target_panel": "rCNS2"})
    q = wr.job_queue_of("preprocessing")
    job = wr.Job(next(_ids), "preprocessing", q, [f"{q}:preprocessing"], 0, ctx)
    bam_preprocessor.bam_preprocessing_handler(job)
    return job


def _master(work):
    with open(os.path.join(str(work), SAMPLE, "master.csv"), newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 1
    return rows[0]


def _counts(work):
    row = _master(work)
    return {k: int(float(row[k])) for k in COUNTED}


def _sum(a, b):
    return {k: a[k] + b[k] for k in a}


def test_a_bam_preprocessed_again_is_counted_once(tmp_path):
    bam = _bam(tmp_path / "run" / "pass", "S_pass_0.bam", 0)
    _run(bam, tmp_path / "work")
    once = _counts(tmp_path / "work")
    assert once["counter_bam_passed"] == 1 and once["counter_mapped_count"] > 0
    _run(bam, tmp_path / "work")
    assert _counts(tmp_path / "work") == once, "a BAM preprocessed again was counted twice"


def test_a_copy_of_a_bam_at_another_path_is_counted_once(tmp_path):
    bam = _bam(tmp_path / "run" / "pass", "S_pass_0.bam", 0)
    copy = tmp_path / "moved" / "pass" / "S_pass_0.bam"
    copy.parent.mkdir(parents=True)
    shutil.copyfile(bam, copy)
    _run(bam, tmp_path / "work")
    once = _counts(tmp_path / "work")
    _run(copy, tmp_path / "work")
    assert _counts(tmp_path / "work") == once, "the same BAM at another path was counted again"


def test_different_bams_are_each_counted(tmp_path):
    # Control: each BAM alone, then both into one sample.
    a = _bam(tmp_path / "run" / "pass", "S_pass_0.bam", 0)
    b = _bam(tmp_path / "run" / "pass", "S_pass_1.bam", 1)
    _run(a, tmp_path / "wa")
    _run(b, tmp_path / "wb")
    _run(a, tmp_path / "both")
    _run(b, tmp_path / "both")
    assert _counts(tmp_path / "both") == _sum(_counts(tmp_path / "wa"), _counts(tmp_path / "wb"))


def test_other_writers_keep_the_record_of_counted_bams(tmp_path):
    # The GUI's overview and the panel update rewrite master.csv too.
    bam = _bam(tmp_path / "run" / "pass", "S_pass_0.bam", 0)
    work = tmp_path / "work"
    _run(bam, work)
    once = _counts(work)
    manager = MasterCSVManager(str(work))
    manager.update_sample_overview(SAMPLE, {"active_jobs": 1, "job_types": "cnv"})
    manager.update_analysis_panel(SAMPLE, "AML")
    _run(bam, work)
    assert _counts(work) == once, "a writer dropped the record of counted BAMs"


def test_a_master_csv_from_before_this_change_is_added_to(tmp_path):
    # Control: a row with no record of counted BAMs keeps its counts and gains the next BAM.
    a = _bam(tmp_path / "run" / "pass", "S_pass_0.bam", 0)
    b = _bam(tmp_path / "run" / "pass", "S_pass_1.bam", 1)
    _run(a, tmp_path / "wa")
    _run(b, tmp_path / "wb")
    work = tmp_path / "work"
    _run(a, work)
    path = os.path.join(str(work), SAMPLE, "master.csv")
    row = _master(work)
    row.pop("bams_added", None)
    row.pop("last_bam_added", None)  # round 2: the record moved to its own list file
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)
    listed = os.path.join(str(work), SAMPLE, "master_csv_bams.txt")
    if os.path.exists(listed):
        os.remove(listed)
    _run(b, work)
    assert _counts(work) == _sum(_counts(tmp_path / "wa"), _counts(tmp_path / "wb"))


def test_a_master_csv_that_cannot_be_read_is_kept_and_reported(tmp_path):
    a = _bam(tmp_path / "run" / "pass", "S_pass_0.bam", 0)
    b = _bam(tmp_path / "run" / "pass", "S_pass_1.bam", 1)
    work = tmp_path / "work"
    _run(a, work)
    path = os.path.join(str(work), SAMPLE, "master.csv")
    with open(path, "wb"):
        pass  # zero bytes
    job = _run(b, work)
    assert os.path.getsize(path) == 0, "the unreadable master.csv was replaced"
    failures = job.context.metadata.get("file_failures") or []
    assert [f["file"] for f in failures] == [str(b)], "the failed update was not reported"


def test_a_failed_master_csv_update_is_reported(tmp_path, monkeypatch):
    bam = _bam(tmp_path / "run" / "pass", "S_pass_0.bam", 0)

    def disk_full(self, data, csv_path):
        raise OSError("disk full")

    monkeypatch.setattr(MasterCSVManager, "_write_csv_internal", disk_full)
    job = _run(bam, tmp_path / "work")
    failures = job.context.metadata.get("file_failures") or []
    assert [f["file"] for f in failures] == [str(bam)], "the failed update was not reported"
    assert "disk full" in failures[0]["error"]


def test_a_bam_is_known_by_its_bytes(tmp_path):
    from robin.bam_identity import bam_key

    a = _bam(tmp_path / "x" / "pass", "S_pass_0.bam", 0)
    b = _bam(tmp_path / "y" / "pass", "S_pass_0.bam", 1)  # same header, other reads
    copy = tmp_path / "z" / "S_pass_0.bam"
    copy.parent.mkdir()
    shutil.copyfile(a, copy)
    assert bam_key(a) == bam_key(str(copy))
    assert bam_key(a) != bam_key(b)


# Round 2: the external review's finding 3 ----------------------------------------------------


def _add(manager, key, mapped):
    return manager.update_master_csv(SAMPLE, {"mapped_reads": mapped}, {"state": "pass"}, key)


def test_a_sample_with_4000_bams_keeps_a_master_csv_every_reader_can_read(tmp_path, monkeypatch):
    # The list of counted BAMs grew in one CSV cell: past 3,640 BAMs, Python's csv module (at its
    # default field limit) can't read the row, and the panel lookup falls back to rCNS2.
    assert csv.field_size_limit() <= 131072, "something raised the csv field limit: this test can't tell"
    monkeypatch.setattr(os, "fsync", lambda fd: None)  # speed only
    manager = MasterCSVManager(str(tmp_path))
    for i in range(4000):
        assert manager.update_master_csv(SAMPLE, {}, {"state": "pass"}, f"{i:032x}")
    manager.update_analysis_panel(SAMPLE, "AML")
    row = _master(tmp_path)  # csv.DictReader, as the panel lookup reads it
    assert (row["analysis_panel"], int(float(row["counter_bam_passed"]))) == ("AML", 4000)


def test_a_stop_between_master_csv_and_its_list_counts_the_bam_once(tmp_path, monkeypatch):
    manager = MasterCSVManager(str(tmp_path))
    real = MasterCSVManager._append_counted_key

    def stopped(self, path, key):
        raise OSError("synthetic: ROBIN stopped after writing master.csv")

    monkeypatch.setattr(MasterCSVManager, "_append_counted_key", stopped)
    with pytest.raises(OSError):
        _add(manager, "a" * 32, 5)
    monkeypatch.setattr(MasterCSVManager, "_append_counted_key", real)
    assert _add(manager, "a" * 32, 5) is False, "a BAM master.csv already counts was counted again"
    assert _add(manager, "b" * 32, 7) is True
    assert int(float(_master(tmp_path)["counter_mapped_count"])) == 12


def test_a_list_cut_off_mid_key_is_repaired(tmp_path):
    manager = MasterCSVManager(str(tmp_path))
    _add(manager, "a" * 32, 5)
    listed = os.path.join(str(tmp_path), SAMPLE, "master_csv_bams.txt")
    with open(listed, "r+b") as fh:
        fh.truncate(10)  # a stop part-way through writing the key
    assert _add(manager, "a" * 32, 5) is False, "a BAM master.csv already counts was counted again"
    assert _add(manager, "b" * 32, 7) is True
    assert int(float(_master(tmp_path)["counter_mapped_count"])) == 12
    with open(listed) as fh:
        assert fh.read().split() == ["a" * 32, "b" * 32]
