"""Batch 10e part 2: a coverage commit adds a batch to every file, or to none (TGT-3).

Adding a batch writes about ten files in the sample folder, then moves the batch BAM to
where finalize looks, then deletes the batch's claim folder. A commit that raised part-way,
or a process that died part-way, left some files with the batch and some without, for good.
After any stop, once the next accumulation or finalize has run, every file the commit writes
and target.bam must describe the same BAMs, and a BAM left out must be listed.

A crash is simulated by raising a BaseException (``_Killed``) from ``os.replace`` or
``shutil.rmtree``: no ``except Exception`` catches it, and the locks are let go, as when a
process dies. That measures a process dying, not a power loss: what must reach the disk
before what is checked by the order of the sync calls (C7). Synthetic BAMs only; samtools
(through pysam) is the reference.
"""

import glob
import json
import os
import shutil

import numpy as np
import pandas as pd
import pytest

import test_target_vs_samtools as tvs
from target_fixtures import (
    MERGED_PANEL, SAMPLE, fixture_reads, samtools_merged_depth, samtools_reference,
    write_bam, write_bams, write_panel,
)
from test_target_failed_files import _ledger_with_ended_target_steps, _rows_at_next_start
from robin.analysis import target_analysis as ta_mod
from robin.analysis.target_analysis import TargetAnalysis, finalize_accumulation_for_sample
from robin.analysis.target_coverage import covered_bases

PANEL_NAME = "TESTPANEL"
K = 3  # the fixture's reads split into 3 BAMs
# Measured with samtools: target B is over it after BAM 0; A and B after BAMs 0 and 1 (or 0
# and 2); A, B and C after all three. At the default 0.1x every target is over it after one
# BAM, so the threshold files would never change and a finish that skipped them would pass.
CALL_THRESHOLD = 1.5
LIST = os.path.join("_commit", "install.json")
# Every file a commit writes in the sample folder (target_coverage_time.csv is appended).
COMMIT_FILES = [
    "_covered_stretches.parquet",
    "coverage_main.csv",
    "bed_coverage_main.csv",
    "target_union_coverage.csv",
    "target_coverage.csv",
    "coverage_time_chart.npy",
    "_target_coverage_latest_reads.parquet",
    "target_coverage_time.csv",
    "targets_exceeding_threshold_count.txt",
    "targets_exceeding_threshold.bed",
]
ATOMIC_FILES = [f for f in COMMIT_FILES if f != "target_coverage_time.csv"]


class _Killed(BaseException):
    """A process dying: nothing catches it, and its locks are let go."""


@pytest.fixture
def panel(tmp_path, monkeypatch):
    path = write_panel(str(tmp_path / "panel.bed"))
    monkeypatch.setattr(TargetAnalysis, "_find_target_bed", lambda self, p: path)
    return path


@pytest.fixture(autouse=True)
def _threshold_crossed_at_each_commit(monkeypatch):
    real = TargetAnalysis._load_config
    monkeypatch.setattr(TargetAnalysis, "_load_config", lambda self: {**real(self), "call_threshold": CALL_THRESHOLD})


@pytest.fixture
def bams(tmp_path):
    (tmp_path / "bams").mkdir()
    return write_bams(str(tmp_path / "bams"), K)[0]


_REFS = {}


def _reference(tmp_path, panel, which):
    """samtools on the reads of the BAMs numbered in ``which``, merged into one BAM, with its
    depth over the hand-merged panel."""
    key = (str(tmp_path), panel, tuple(sorted(which)))
    if key not in _REFS:
        reads = fixture_reads()
        parts = [r for j in sorted(which) for r in reads[j::K]]
        name = "_".join(map(str, sorted(which)))
        merged = write_bam(str(tmp_path / f"ref_{name}.bam"), parts)
        ref = samtools_reference(merged, panel)
        folder = tmp_path / f"ref_{name}"
        folder.mkdir(exist_ok=True)
        ref["panel_depth"] = samtools_merged_depth(merged, MERGED_PANEL, str(folder))
        _REFS[key] = ref
    return _REFS[key]


def _analysis(work):
    return TargetAnalysis(work_dir=str(work), target_panel=PANEL_NAME, batch_size=1, use_staging=True)


def _stage(work, bam):
    meta, _ = _analysis(work).process_file_with_staging(bam, {"sample_id": SAMPLE})
    assert not meta.error_message, meta.error_message


def _add(work, bam):
    _stage(work, bam)
    return _analysis(work).accumulate_staged_files(SAMPLE, force=True)


def _snapshot(sample_dir):
    out = {}
    for name in COMMIT_FILES:
        path = os.path.join(sample_dir, name)
        out[name] = open(path, "rb").read() if os.path.exists(path) else None
    return out


def _time_csv_problems(sample_dir, ref, batches):
    """The time series: one block per counted batch, the last block's reads as samtools'."""
    path = os.path.join(sample_dir, "target_coverage_time.csv")
    if not os.path.exists(path):
        return ["target_coverage_time.csv missing"]
    try:
        df = pd.read_csv(path)
    except Exception as e:
        return [f"target_coverage_time.csv can't be read: {e}"]
    problems = []
    if df["timestamp"].nunique() != batches:
        problems.append(f"target_coverage_time.csv has {df['timestamp'].nunique()} blocks, {batches} batches counted")
    if len(df) != batches * len(ref["targets"]):
        problems.append(f"target_coverage_time.csv has {len(df)} rows, {batches * len(ref['targets'])} expected")
    last = df[df["timestamp"] == df["timestamp"].max()]
    for key, (_, _, reads) in ref["targets"].items():
        row = last[(last["chrom"] == key[0]) & (last["startpos"] == key[1]) & (last["name"] == key[3])]
        if len(row) != 1 or int(row["reads"].iloc[0]) != reads:
            problems.append(f"target_coverage_time.csv last block: {key[3]} reads {list(row['reads'])} vs {reads}")
    return problems


def _other_file_problems(work, sample_dir, robin, ref):
    """The files the totals check doesn't read: the mirrored targets, the threshold files,
    the covered stretches (the next commit's covered bases), and the panel's summary depth
    (through the union file)."""
    problems = []
    if not tvs._close(robin["panel_depth"], ref["panel_depth"]):
        problems.append(f"panel depth {robin['panel_depth']} vs samtools {ref['panel_depth']:.6f}")
    bed = pd.read_csv(os.path.join(sample_dir, "bed_coverage_main.csv"))
    try:
        mirror = pd.read_csv(os.path.join(sample_dir, "target_coverage.csv"))
        if not mirror.equals(bed):
            problems.append("target_coverage.csv differs from bed_coverage_main.csv")
    except Exception as e:
        problems.append(f"target_coverage.csv can't be read: {e}")
    threshold = _analysis(work).callthreshold
    want = sorted(
        (c, s, e) for (c, s, e, _), (_, depth, _) in ref["targets"].items() if depth >= threshold
    )
    try:
        count = int(open(os.path.join(sample_dir, "targets_exceeding_threshold_count.txt")).read())
        bed_path = os.path.join(sample_dir, "targets_exceeding_threshold.bed")
        got = []
        if os.path.getsize(bed_path):
            rows = pd.read_csv(bed_path, sep="\t", header=None)
            got = sorted((str(r[0]), int(r[1]), int(r[2])) for r in rows.itertuples(index=False))
    except Exception as e:
        return problems + [f"threshold files can't be read: {e}"]
    if count != len(want) or got != want:
        problems.append(f"threshold files: count {count}, rows {got} vs {want}")
    stretches = covered_bases(pd.read_parquet(os.path.join(sample_dir, "_covered_stretches.parquet")))
    for name, (_, covbases, _) in ref["contigs"].items():
        if stretches.get(name, 0) != covbases:
            problems.append(f"covered stretches: {name} {stretches.get(name, 0)} vs samtools {covbases}")
    return problems


def _agreement(tmp_path, panel, work, added, maybe, finalize_result, must_include=False):
    """Every file a commit writes and target.bam describe the same BAMs: ``added``, with or
    without ``maybe`` (a BAM whose commit stopped; with it only, when ``must_include``).
    Without it, finalize must list its claim folder (or its incomplete staged set)."""
    sample_dir = os.path.join(str(work), SAMPLE)
    robin = tvs.read_back(str(work))
    tried = {}
    options = [added | {maybe}] if must_include else [added | {maybe}, added]
    for which in options:
        ref = _reference(tmp_path, panel, which)
        problems = (
            tvs.total_problems(robin, ref)
            + tvs.target_bam_problems(robin, ref)
            + _time_csv_problems(sample_dir, ref, len(which))
            + _other_file_problems(work, sample_dir, robin, ref)
        )
        tried[tuple(sorted(which))] = problems
        if not problems:
            left = (finalize_result.get("claims_left") or []) + (finalize_result.get("incomplete_staged") or [])
            if maybe in which:
                return [] if not left else [f"BAM {maybe} counted and merged, but its claim folder is listed: {left}"]
            return [] if left else [f"BAM {maybe} left out, and no claim folder lists it"]
    return [f"no set of BAMs fits: {tried}"]


@pytest.mark.parametrize("stale", ATOMIC_FILES)
def test_the_checker_catches_a_file_from_an_earlier_commit(tmp_path, panel, bams, stale):
    """The checker is not blind: any one file put back as it was after the first BAM fails it
    (the external review built such a case against an earlier checker)."""
    work = tmp_path / "work"
    sample_dir = os.path.join(str(work), SAMPLE)
    assert _add(work, bams[0])["status"] == "success"
    first = _snapshot(sample_dir)
    assert _add(work, bams[1])["status"] == "success"
    result = finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
    assert _agreement(tmp_path, panel, work, {0}, 1, result) == []
    put_back = first[stale]
    assert put_back is not None
    assert put_back != _snapshot(sample_dir)[stale], f"{stale} is the same after one BAM as after two"
    with open(os.path.join(sample_dir, stale), "wb") as f:
        f.write(put_back)
    assert _agreement(tmp_path, panel, work, {0}, 1, result) != [], f"a stale {stale} passed"


# ------------------------------------------------------------------------------ a commit that raised


@pytest.mark.parametrize("failing", ATOMIC_FILES)
def test_a_commit_that_raises_part_way_changes_no_file(tmp_path, panel, bams, monkeypatch, failing):
    work = tmp_path / "work"
    sample_dir = os.path.join(str(work), SAMPLE)
    assert _add(work, bams[0])["status"] == "success"
    before = _snapshot(sample_dir)
    real = ta_mod.atomic_open

    def atomic_open(path, *args, **kwargs):
        if os.path.basename(str(path)).startswith(failing):
            raise OSError("synthetic disk full")
        return real(path, *args, **kwargs)

    monkeypatch.setattr(ta_mod, "atomic_open", atomic_open)
    assert _add(work, bams[1])["status"] == "error"
    monkeypatch.setattr(ta_mod, "atomic_open", real)
    after = _snapshot(sample_dir)
    changed = [name for name in COMMIT_FILES if before[name] != after[name]]
    assert changed == [], f"a failed commit changed {changed}"
    assert _add(work, bams[2])["status"] == "success"
    result = finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
    assert _agreement(tmp_path, panel, work, {0, 2}, 1, result) == []


# --------------------------------------------------------------------------------- a commit killed


def _in(path, folder):
    return f"{os.sep}{folder}" in str(path)


def _is_install_step(src, dst):
    """A rename that moves something out of _commit/ or out of a claim folder."""
    return (_in(src, "_commit") and not _in(dst, "_commit")) or (_in(src, "_batch_") and not _in(dst, "_batch_"))


def _kill_at(monkeypatch, n, only=None):
    """Die at the n-th os.replace or claim-folder removal from now (``only``: count only the
    install steps, only the batch BAM's moves, or only the claim folder's removal). Returns a
    function that undoes it and gives the number of calls counted."""
    real_replace, real_rmtree = os.replace, shutil.rmtree
    calls = {"n": 0}

    def count():
        calls["n"] += 1
        if calls["n"] == n:
            raise _Killed(f"killed at step {n}")

    def replace(src, dst, *args, **kwargs):
        if only is None or (only == "install" and _is_install_step(src, dst)) or (
            only == "batch_bam" and _in(src, "_batch_") and not _in(dst, "_batch_")
        ):
            count()
        return real_replace(src, dst, *args, **kwargs)

    def rmtree(path, *args, **kwargs):
        if os.path.basename(str(path).rstrip("/")).startswith("_batch_") and only in (None, "install", "drop_claim"):
            count()
        return real_rmtree(path, *args, **kwargs)

    monkeypatch.setattr(os, "replace", replace)
    monkeypatch.setattr(shutil, "rmtree", rmtree)

    def stop():
        monkeypatch.setattr(os, "replace", real_replace)
        monkeypatch.setattr(shutil, "rmtree", real_rmtree)
        return calls["n"]

    return stop


def _killed(monkeypatch, n, only, fn):
    """Run fn with a kill at its n-th counted step: (killed, steps counted)."""
    stop = _kill_at(monkeypatch, n, only)
    try:
        fn()
        killed = False
    except _Killed:
        killed = True
    return killed, stop()


@pytest.mark.parametrize("recovery", ["next_accumulation", "finalize"])
def test_a_commit_killed_at_any_step_leaves_the_files_and_target_bam_agreeing(tmp_path, panel, bams, monkeypatch, recovery):
    broken, n = [], 1
    while True:
        work = tmp_path / f"work_{n}"
        assert _add(work, bams[0])["status"] == "success"
        _stage(work, bams[1])  # staging a BAM is not part of the commit (the ledger covers it)
        killed, calls = _killed(monkeypatch, n, None, lambda: _analysis(work).accumulate_staged_files(SAMPLE, force=True))
        if not killed:
            assert calls < n
            break  # every step of the commit has been a kill point
        added = {0}
        if recovery == "next_accumulation":
            assert _add(work, bams[2])["status"] == "success"
            added = {0, 2}
        result = finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
        problems = _agreement(tmp_path, panel, work, added, 1, result)
        if problems:
            broken.append((n, problems[0][:300]))
        n += 1
    assert n > 10, f"only {n - 1} kill points: the commit was not reached"
    assert broken == [], f"{len(broken)} of {n - 1} kill points broke: {broken}"


def test_a_finish_killed_at_any_step_is_finished_by_the_next_one(tmp_path, panel, bams, monkeypatch):
    """A second stop, during the recovery: a commit killed during its install, then the
    finalize that finishes it killed too, at each of its install steps; the next finalize
    must finish it. Once its list was written, the commit counts: BAM 1 must be in."""
    def first_kill(work, n1):
        assert _add(work, bams[0])["status"] == "success"
        _stage(work, bams[1])
        return _killed(monkeypatch, n1, "install", lambda: _analysis(work).accumulate_staged_files(SAMPLE, force=True))

    _, total = first_kill(tmp_path / "count", 10_000)  # the install steps of one commit
    assert total >= 3, f"a commit has only {total} install steps"
    broken, with_a_list = [], 0
    for n1 in sorted({1, (total + 1) // 2, total}):
        m = 1
        while True:
            work = tmp_path / f"work_{n1}_{m}"
            killed, _ = first_kill(work, n1)
            assert killed
            if not os.path.exists(os.path.join(str(work), SAMPLE, LIST)):
                break  # no install list: nothing for a finish to do
            with_a_list += 1
            ran = []
            killed, _ = _killed(
                monkeypatch, m, "install",
                lambda: ran.append(finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)),
            )
            # Killed: the next finalize must finish it. Not killed: the finish ran through,
            # and its own result is checked (found by mutant Q13: a finish with no steps
            # was never checked).
            result = finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME) if killed else ran[0]
            problems = _agreement(tmp_path, panel, work, {0}, 1, result, must_include=True)
            if problems:
                broken.append((n1, m, problems[0][:300]))
            if not killed:
                break
            m += 1
    assert with_a_list > 0, "no kill point left an install list"
    assert broken == [], f"{len(broken)} broke: {broken}"


def test_a_time_series_append_cut_off_part_way_is_put_right(tmp_path, panel, bams, monkeypatch):
    """A crash part-way through the append: the list is still there, and the file holds part
    of the new rows. The finish cuts it back and appends them whole."""
    work = tmp_path / "work"
    sample_dir = os.path.join(str(work), SAMPLE)
    assert _add(work, bams[0])["status"] == "success"
    time_csv = os.path.join(sample_dir, "target_coverage_time.csv")
    size_before = os.path.getsize(time_csv)
    _stage(work, bams[1])
    killed, _ = _killed(monkeypatch, 1, "batch_bam", lambda: _analysis(work).accumulate_staged_files(SAMPLE, force=True))
    assert killed
    assert os.path.exists(os.path.join(sample_dir, LIST)), "no install list after the append"
    size_after = os.path.getsize(time_csv)
    assert size_after > size_before
    with open(time_csv, "r+b") as f:
        f.truncate(size_before + (size_after - size_before) // 2)  # only part of the rows landed
    result = finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
    assert _agreement(tmp_path, panel, work, {0}, 1, result, must_include=True) == []


def test_an_error_during_the_install_leaves_the_list_for_the_next_accumulation(tmp_path, panel, bams, monkeypatch):
    """C3: the list is written, then a rename fails (an error, not a crash). The commit counts
    as done: BAM 1 is not reported "not added", and the next accumulation finishes it."""
    work = tmp_path / "work"
    sample_dir = os.path.join(str(work), SAMPLE)
    assert _add(work, bams[0])["status"] == "success"
    _stage(work, bams[1])
    real_replace = os.replace
    failed = {"n": 0}

    def replace(src, dst, *args, **kwargs):
        if _is_install_step(src, dst) and not failed["n"]:
            failed["n"] += 1
            raise OSError("synthetic: a rename fails")
        return real_replace(src, dst, *args, **kwargs)

    monkeypatch.setattr(os, "replace", replace)
    result = _analysis(work).accumulate_staged_files(SAMPLE, force=True)
    monkeypatch.setattr(os, "replace", real_replace)
    assert failed["n"] == 1
    assert result["status"] == "error"
    assert not result.get("file_failures"), "a decided commit's BAM was reported not added"
    assert os.path.exists(os.path.join(sample_dir, LIST)), "the install list is gone"
    assert _add(work, bams[2])["status"] == "success"
    result = finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
    assert _agreement(tmp_path, panel, work, {0, 2}, 1, result, must_include=True) == []


def test_an_install_list_that_cannot_be_read_stops_the_accumulation(tmp_path, panel, bams):
    """C4: a list damaged from outside is not guessed at; the accumulation stops, names it,
    and changes no file."""
    work = tmp_path / "work"
    sample_dir = os.path.join(str(work), SAMPLE)
    assert _add(work, bams[0])["status"] == "success"
    before = _snapshot(sample_dir)
    os.makedirs(os.path.join(sample_dir, "_commit"), exist_ok=True)
    with open(os.path.join(sample_dir, LIST), "w") as f:
        f.write('{"steps": [')
    result = _add(work, bams[1])
    assert result["status"] == "error"
    assert "install.json" in result.get("error", "")
    after = _snapshot(sample_dir)
    assert [n for n in COMMIT_FILES if before[n] != after[n]] == []


@pytest.mark.parametrize("case", ["merged", "lost"])
def test_a_finish_whose_batch_bam_is_gone_looks_for_it_in_target_bam(tmp_path, panel, bams, monkeypatch, case):
    """C8: the commit moved its batch BAM into place and died before deleting its list. The
    batch BAM is then gone from the sample folder too:
    - merged: a finalize merged it and deleted it (it merges under its own lock). It is in
      target.bam: nothing is reported, and every file and target.bam agree, with BAM 1 in;
    - lost (from outside): BAM 1's reads are not in target.bam. It is reported as a failed
      file and in the ledger, and its claim folder is kept and listed."""
    work = tmp_path / "work"
    sample_dir = os.path.join(str(work), SAMPLE)
    _ledger_with_ended_target_steps(work, bams)  # each BAM's target job ended well
    assert _add(work, bams[0])["status"] == "success"
    _stage(work, bams[1])
    killed, _ = _killed(
        monkeypatch, 1, "drop_claim", lambda: _analysis(work).accumulate_staged_files(SAMPLE, force=True)
    )
    assert killed
    with open(os.path.join(sample_dir, LIST)) as f:
        step = next(s for s in json.load(f)["steps"] if s["do"] == "batch_bam")
    moved = os.path.join(sample_dir, step["to"])
    assert os.path.exists(moved), "the batch BAM was not in place before the claim folder's removal"
    if case == "merged":
        gone = sorted(glob.glob(os.path.join(sample_dir, "batch_*.bam")))
        ta_mod._merge_into_target_bam(
            gone, os.path.join(sample_dir, "target.bam"), [os.path.basename(b) for b in gone]
        )
    else:
        gone = [moved]
    for b in gone:
        os.remove(b)
        os.remove(f"{b}.bai")
    result = finalize_accumulation_for_sample(SAMPLE, str(work), PANEL_NAME)
    assert not os.path.exists(os.path.join(sample_dir, LIST)), "the finish did not end"
    if case == "merged":
        assert not result.get("file_failures"), result.get("file_failures")
        assert _agreement(tmp_path, panel, work, {0}, 1, result, must_include=True) == []
        assert _rows_at_next_start(work, bams[1]) == []
    else:
        assert [f["file"] for f in result.get("file_failures") or []] == [bams[1]]
        assert result.get("claims_left"), "the kept claim folder is not listed"
        assert _rows_at_next_start(work, bams[1]) == [("target", "failed_before")]
    assert _rows_at_next_start(work, bams[0]) == []


# ----------------------------------------------------------------------------- what reaches the disk


def _under(path, folder):
    return path == folder or path.startswith(folder + os.sep)


def _sync_order_problems(events, claim_kept):
    """Each sync in its place (steps 2, 3 and 5 of the plan). ``events``, in order, with paths
    relative to the sample folder: ("file"|"dir", p) syncs, ("replace", src, dst),
    ("remove", p), ("mkdir", p) and ("rmtree"|"rmdir", p)."""
    def find(pred, start=0, end=None):
        stop = len(events) if end is None else end
        return next((i for i in range(start, stop) if pred(events[i])), None)

    def last(pred, end):
        found = [i for i in range(end) if pred(events[i])]
        return found[-1] if found else None

    def synced(kind, path, start, end):
        return start is not None and end is not None and find(lambda e: e == (kind, path), start, end) is not None

    made = find(lambda e: e == ("mkdir", "_commit"))
    written = find(lambda e: e[0] == "replace" and e[2] == LIST)
    if made is None or written is None or made > written:
        return [f"_commit/ made at {made}, the list written at {written}"]
    deleted = find(lambda e: e == ("remove", LIST), written)
    if deleted is None:
        return ["the install list was not deleted"]
    steps = [
        i for i in range(written, deleted)
        if events[i][0] == "replace" and (
            (_under(events[i][1], "_commit") and not _under(events[i][2], "_commit"))
            or ("_batch_" in events[i][1] and "_batch_" not in events[i][2])
        )
    ]
    moved_out = [events[i][1] for i in steps if "_batch_" in events[i][1]]
    if not steps or len(moved_out) != 2:
        return [f"install renames {[events[i] for i in steps]}: expected the index and the batch BAM among them"]
    claim = os.path.dirname(moved_out[0])
    problems = []
    if not synced("dir", ".", made + 1, written):
        problems.append("the sample folder not synced after _commit/ was made, before the list")
    prepared = last(lambda e: e[0] == "replace" and _under(e[2], "_commit") and e[2] != LIST, written)
    if not synced("dir", "_commit", (prepared or 0) + 1, written):
        problems.append("_commit/ not synced after its prepared files, before the list")
    file_syncs = [find(lambda e, p=p: e == ("file", p), 0, written) for p in moved_out]
    if None in file_syncs:
        problems.append(f"{moved_out} not each synced before the list")
    else:
        if not synced("dir", claim, max(file_syncs) + 1, written):
            problems.append("the claim folder not synced after the batch BAM, before the list")
        if not synced("dir", "_staging", max(file_syncs) + 1, written):
            problems.append("_staging/ not synced after the batch BAM, before the list")
    if not synced("dir", "_commit", written + 1, steps[0]):
        problems.append("_commit/ not synced after the list, before the first install step")
    dropped = find(lambda e: e == ("rmtree", claim), written, deleted)
    if claim_kept:
        if dropped is not None:
            return problems + ["the claim folder was removed, but should be kept"]
        end = steps[-1]
        needs = [("file", "target_coverage_time.csv"), ("dir", "."), ("dir", claim)]
    else:
        if dropped is None:
            return problems + ["the claim folder was not removed during the install"]
        end = max(dropped, steps[-1])
        needs = [("file", "target_coverage_time.csv"), ("dir", "."), ("dir", "_staging")]
    for kind, path in needs:
        if not synced(kind, path, end + 1, deleted):
            problems.append(f"{(kind, path)} not synced after the install, before the list was deleted")
    gone = find(lambda e: e[0] in ("rmtree", "rmdir") and e[1] == "_commit", deleted)
    if not synced("dir", "_commit", deleted + 1, gone):
        problems.append("_commit/ not synced after the list was deleted, before _commit/ was")
    if gone is None or not synced("dir", ".", gone + 1, len(events)):
        problems.append("the sample folder not synced after _commit/ was deleted")
    return problems


@pytest.mark.parametrize("claim", ["dropped", "kept"])
def test_each_sync_comes_in_its_place(tmp_path, panel, bams, monkeypatch, claim):
    """A power loss can't be made here; the order of the syncs is checked instead (C7). The
    claim folder is kept when one of its sets can't be read (BAM 2's, here)."""
    work = tmp_path / "work"
    sample_dir = os.path.join(str(work), SAMPLE)
    assert _add(work, bams[0])["status"] == "success"
    _stage(work, bams[1])
    if claim == "kept":
        staged = glob.glob(os.path.join(sample_dir, "_staging", "coverage_*.parquet"))
        _stage(work, bams[2])
        (broken,) = set(glob.glob(os.path.join(sample_dir, "_staging", "coverage_*.parquet"))) - set(staged)
        with open(broken, "wb") as f:
            f.write(b"not a parquet file")
    events = []

    def rel(p):
        return os.path.relpath(os.path.abspath(str(p)), sample_dir)

    def spy(kind, real):
        def call(p, *args, **kwargs):
            events.append((kind, rel(p)))
            return real(p, *args, **kwargs)
        return call

    monkeypatch.setattr(ta_mod, "_sync_file", lambda p: events.append(("file", rel(p))), raising=False)
    monkeypatch.setattr(ta_mod, "_sync_dir", lambda p: events.append(("dir", rel(p))), raising=False)
    real_replace = os.replace

    def replace(src, dst, *args, **kwargs):
        events.append(("replace", rel(src), rel(dst)))
        return real_replace(src, dst, *args, **kwargs)

    monkeypatch.setattr(os, "replace", replace)
    monkeypatch.setattr(os, "remove", spy("remove", os.remove))
    monkeypatch.setattr(os, "mkdir", spy("mkdir", os.mkdir))  # os.makedirs calls it too
    monkeypatch.setattr(os, "rmdir", spy("rmdir", os.rmdir))
    monkeypatch.setattr(shutil, "rmtree", spy("rmtree", shutil.rmtree))
    result = _analysis(work).accumulate_staged_files(SAMPLE, force=True)
    monkeypatch.undo()
    assert result["status"] == "success", result
    assert result["files_processed"] == 1
    assert _sync_order_problems(events, claim == "kept") == []
