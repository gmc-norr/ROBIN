"""Step 6's locked rule (plan 7278fb9): ROBIN's target coverage, read counts and
target.bam equal samtools run once on all the sample's reads in one BAM.

The same synthetic reads, split into 1, 2 and 4 BAMs, go through the real entry
points (process_file_with_staging, accumulate_staged_files,
finalize_accumulation_for_sample, process_bam_reads) in five ways: in order,
reversed, all at once, two processes at the same time, and finalized twice. Every
run must give samtools' numbers, and the same numbers as every other run.
Problems are listed as short lines (a long list diff makes pytest look hung).
"""

import multiprocessing
import os
import time
from collections import Counter
from typing import Dict

import numpy as np
import pandas as pd
import pytest

from target_fixtures import (
    CONTIGS, SAMPLE, bam_records, samtools_reference, write_bams, write_panel,
)
from robin.analysis import target_analysis as ta_mod
from robin.analysis.bam_preprocessor import process_bam_reads
from robin.analysis.target_analysis import TargetAnalysis, finalize_accumulation_for_sample

PANEL_NAME = "TESTPANEL"
WAYS = ["in_order", "reversed", "all_at_once", "two_processes", "finalize_twice"]
_RUNS: Dict = {}


def _analysis(work):
    return TargetAnalysis(work_dir=work, target_panel=PANEL_NAME, batch_size=1, use_staging=True)


def _stage_and_accumulate(work, bams):
    ta = _analysis(work)
    for bam in bams:
        meta, _ = ta.process_file_with_staging(bam, {"sample_id": SAMPLE})
        assert not meta.error_message, meta.error_message
        ta.accumulate_staged_files(SAMPLE, force=True)


def _slow(fn):
    """Widen the gap between reading and writing the totals (outside the lock)."""
    def slow(*args, **kwargs):
        time.sleep(0.3)
        return fn(*args, **kwargs)
    return slow


def _worker(work, bams, panel):
    """One of two processes accumulating at the same time (a fresh interpreter: spawn)."""
    TargetAnalysis._find_target_bed = lambda self, p: panel
    ta_mod.get_read_counts_per_target = _slow(ta_mod.get_read_counts_per_target)
    if hasattr(ta_mod, "covered_bases"):
        # Also slow a step between reading and writing the totals inside the commit,
        # so a commit that isn't locked loses an update (a survivor of the first
        # mutation run; added after seeing it, rules unchanged).
        ta_mod.covered_bases = _slow(ta_mod.covered_bases)
    _stage_and_accumulate(work, bams)


def run_robin(work, bams, way, panel):
    os.makedirs(work, exist_ok=True)
    if way == "in_order":
        _stage_and_accumulate(work, bams)
    elif way == "reversed":
        _stage_and_accumulate(work, bams[::-1])
    elif way == "all_at_once":
        ta = _analysis(work)
        for bam in bams:
            ta.process_file_with_staging(bam, {"sample_id": SAMPLE})
        ta.accumulate_staged_files(SAMPLE, force=True)
    elif way == "two_processes":
        ctx = multiprocessing.get_context("spawn")
        procs = [ctx.Process(target=_worker, args=(work, half, panel)) for half in (bams[0::2], bams[1::2])]
        for p in procs:
            p.start()
        for p in procs:
            p.join(120)
        assert [p.exitcode for p in procs] == [0, 0]
    elif way == "finalize_twice":
        half = (len(bams) + 1) // 2
        _stage_and_accumulate(work, bams[:half])
        finalize_accumulation_for_sample(SAMPLE, work, PANEL_NAME)
        _stage_and_accumulate(work, bams[half:])
    finalize_accumulation_for_sample(SAMPLE, work, PANEL_NAME)
    return read_back(work)


def read_back(work) -> Dict:
    sd = os.path.join(work, SAMPLE)
    out = {"targets": {}, "contigs": {}, "global": None, "time_last": None, "target_bam": Counter()}
    cov_path = os.path.join(sd, "coverage_main.csv")
    if os.path.exists(cov_path):
        cov = pd.read_csv(cov_path)
        length = cov["endpos"] - cov["startpos"] + 1
        out["global"] = float((cov["meandepth"] * length).sum() / length.sum())
        for r in cov.itertuples(index=False):
            out["contigs"][r[0]] = (int(r.numreads), r.covbases, float(r.meandepth))
    reads = {}
    reads_path = os.path.join(sd, "_target_coverage_latest_reads.parquet")
    if os.path.exists(reads_path):
        for r in pd.read_parquet(reads_path).itertuples(index=False):
            reads[(r.chrom, int(r.startpos), int(r.endpos), r.name)] = int(r.reads)
    bed_path = os.path.join(sd, "bed_coverage_main.csv")
    if os.path.exists(bed_path):
        for r in pd.read_csv(bed_path).itertuples(index=False):
            key = (r.chrom, int(r.startpos), int(r.endpos), r.name)
            out["targets"][key] = (int(r.bases), float(r.coverage), reads.get(key))
    chart_path = os.path.join(sd, "coverage_time_chart.npy")
    if os.path.exists(chart_path):
        out["time_last"] = float(np.load(chart_path)[-1][1])
    target_bam = os.path.join(sd, "target.bam")
    if os.path.exists(target_bam):
        out["target_bam"] = bam_records(target_bam)
    return out


def _close(a, b) -> bool:
    return a is not None and abs(float(a) - b) <= 1e-9 * max(1.0, abs(b))


def total_problems(robin, ref):
    """Locked rules 1 and 2 (and the per-target read counts of rule 1)."""
    problems = []
    for key, (bases, depth, reads) in ref["targets"].items():
        got = robin["targets"].get(key)
        if got is None:
            problems.append(f"rule 1: {key[3]} missing")
        elif got[0] != bases or not _close(got[1], depth) or got[2] != reads:
            problems.append(f"rule 1: {key[3]} robin (bases, depth, reads) {got} vs samtools {(bases, round(depth, 6), reads)}")
    for name, (numreads, covbases, depth) in ref["contigs"].items():
        got = robin["contigs"].get(name)
        if got is None:
            problems.append(f"rule 1: contig {name} missing")
        elif got[0] != numreads or got[1] != covbases or not _close(got[2], depth):
            problems.append(f"rule 1: {name} robin (numreads, covbases, depth) {got} vs samtools {(numreads, covbases, round(depth, 6))}")
    if not _close(robin["global"], ref["global"]):
        problems.append(f"rule 2: global depth {robin['global']} vs samtools {ref['global']:.6f}")
    if not _close(robin["time_last"], ref["global"]):
        problems.append(f"rule 2: time chart's last value {robin['time_last']} vs samtools {ref['global']:.6f}")
    return problems


def target_bam_problems(robin, ref):
    got, want = robin["target_bam"], ref["target_bam"]
    if got == want:
        return []
    extra, missing = got - want, want - got
    return [f"rule 3: target.bam has {sum(got.values())} records, samtools view -L {sum(want.values())}; "
            f"extra {sorted(extra)[:4]}, missing {sorted(missing)[:4]}"]


def _run(tmp_path_factory, k, way):
    if (k, way) not in _RUNS:
        d = tmp_path_factory.mktemp(f"run_{k}_{way}")
        with pytest.MonkeyPatch.context() as mp:
            panel = write_panel(str(d / "panel.bed"))
            mp.setattr(TargetAnalysis, "_find_target_bed", lambda self, p: panel)
            bams, merged = write_bams(str(d), k)
            ref = samtools_reference(merged, panel)
            robin = run_robin(str(d / "work"), bams, way, panel)
        _RUNS[(k, way)] = (robin, ref)
    return _RUNS[(k, way)]


@pytest.mark.parametrize("way", WAYS)
@pytest.mark.parametrize("k", [1, 2, 4])
def test_totals_match_samtools(tmp_path_factory, k, way):
    robin, ref = _run(tmp_path_factory, k, way)
    problems = total_problems(robin, ref)
    assert not problems, "\n".join(problems)


@pytest.mark.parametrize("way", WAYS)
@pytest.mark.parametrize("k", [1, 2, 4])
def test_target_bam_matches_samtools_view_L(tmp_path_factory, k, way):
    robin, ref = _run(tmp_path_factory, k, way)
    problems = target_bam_problems(robin, ref)
    assert not problems, "\n".join(problems)


@pytest.mark.parametrize("k", [1, 2, 4])
def test_read_counts_match_samtools(tmp_path, k):
    bams, merged = write_bams(str(tmp_path), k)
    want = samtools_reference(merged, write_panel(str(tmp_path / "panel.bed")))["counts"]
    got = Counter()
    for bam in bams:
        counts = process_bam_reads(bam)
        for key in want:
            got[key] += counts[key]
    assert dict(got) == want


def test_every_run_gives_the_same_result(tmp_path_factory):
    """Locked rule 4: every k and every way writes the same numbers and target.bam."""
    first = None
    problems = []
    for k in (1, 2, 4):
        for way in WAYS:
            robin, _ = _run(tmp_path_factory, k, way)
            summary = (robin["targets"], robin["contigs"], robin["global"], robin["target_bam"])
            if first is None:
                first = (k, way, summary)
            elif summary != first[2]:
                problems.append(f"rule 4: k={k} {way} differs from k={first[0]} {first[1]}")
    assert not problems, "\n".join(problems)
