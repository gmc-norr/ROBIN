"""Step 5's locked rule (plan a95935d): ROBIN's methylation counts equal modkit's.

The same synthetic reads go through run_matkit + merge_modkit_files in 1, 2
and 4 batches, and modkit pileup (--cpg --combine-mods --combine-strands
--filter-threshold 0.73) runs on all of them at once, minus the reads whose
primary has qs < 12 (ROBIN's own rule). Per panel CpG: same sites, n_mod and
n_canonical exact, valid_cov their sum, percent within 0.005 of modkit's.

The reads hold every trap from the code review: both strands (METH-4), unsure
calls (METH-3), a 5mC/5hmC split, implicit MM modes (METH-13), low-QS reads
with supplementary alignments before and after the primary (METH-14), and a
panel with a duplicate row, a CGCG pair and a negative start (METH-5, METH-6).
"""

import os
from typing import Dict, List, Tuple

import pyarrow.parquet as pq
import pytest

from meth_fixtures import (
    CANON, CHROM, FAIL, MOD, ReadSpec, filter_bam, modkit_available, modkit_pileup,
    write_meth_bam, write_panel, write_reference,
)

pytestmark = pytest.mark.skipif(not modkit_available(), reason="modkit is not installed")

CPGS = [500, 900, 1300, 1700, 2100, 2500, 2502]
PANEL = [500, 900, 1700, 2100, 2500, 2502]  # 1300 is not in the panel
LOW_QS = {"X", "Y"}
S, E = 400, 2650


def fixture_reads() -> List[ReadSpec]:
    reads = []
    for i in range(6):
        reads.append(ReadSpec(f"pa{i}", "+", S, E, m={
            500: MOD, 900: MOD if i < 3 else CANON, 1300: MOD, 1700: CANON,
            2100: FAIL if i < 2 else MOD, 2500: MOD, 2502: CANON}))
    for i in range(4):
        reads.append(ReadSpec(f"ma{i}", "-", S, E, m={
            500: CANON, 900: MOD, 1300: CANON, 1700: MOD if i < 2 else FAIL,
            2100: MOD, 2500: CANON, 2502: MOD}))
    reads += [
        ReadSpec("mh", "+", S, E, m={900: 115}, h={900: 115}),  # 5mC 0.45 + 5hmC 0.45
        ReadSpec("imp", "+", S, E, m={500: MOD}, mode="."),  # skipped Cs are canonical
        ReadSpec("nomark", "-", S, E, m={1700: MOD}, mode=""),  # no marker means implicit
        ReadSpec("empty", "+", S, E, m={}, mode="."),  # every C canonical
        # qs 5: supplementary first in coordinate order, then primary first.
        ReadSpec("X", "+", 1600, E, m={1700: MOD, 2100: MOD, 2500: MOD, 2502: MOD}, qs=5),
        ReadSpec("X", "+", S, 1000, m={500: MOD, 900: MOD}, qs=None, supplementary=True),
        ReadSpec("Y", "+", S, 1000, m={500: MOD, 900: MOD}, qs=5),
        ReadSpec("Y", "+", 1600, E, m={1700: MOD, 2100: MOD, 2500: MOD, 2502: MOD}, qs=None, supplementary=True),
        # A good read split in two: each part counts at its own sites.
        ReadSpec("Z", "-", S, 1500, m={500: MOD, 900: CANON, 1300: MOD}),
        ReadSpec("Z", "-", 1600, E, m={1700: CANON, 2100: MOD, 2500: MOD, 2502: CANON}, qs=None, supplementary=True),
    ]
    return reads


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    d = str(tmp_path_factory.mktemp("meth"))
    ref = write_reference(d, CPGS)
    union = write_meth_bam(os.path.join(d, "union.bam"), ref, fixture_reads())
    kept = filter_bam(union, os.path.join(d, "kept.bam"), LOW_QS)
    combined = modkit_pileup(kept, ref, d, combine_strands=True)
    stranded = modkit_pileup(kept, ref, d, combine_strands=False)
    panel_rows = [(p, p + 2, f"cg{p}") for p in PANEL] + [(900, 902, "cg900dup"), (-1, 0, "cgneg")]
    panel = write_panel(os.path.join(d, "parquet_filter.txt"), panel_rows)
    return {
        "dir": d, "ref": ref, "union": union, "panel": panel,
        "combined": {k: v for k, v in combined.items() if k[0] in PANEL},
        "stranded": stranded,
    }


def _batches(k: int) -> List[List[ReadSpec]]:
    """Split the reads into k groups, keeping each read's records together."""
    names = list(dict.fromkeys(r.name for r in fixture_reads()))
    groups = [[] for _ in range(k)]
    for r in fixture_reads():
        groups[names.index(r.name) % k].append(r)
    return [g for g in groups if g]


def _robin_cumulative(world, k: int, tmp) -> Dict[Tuple[int, str], Tuple[int, int, int, float]]:
    from robin.analysis.temp_utilities import merge_modkit_files
    from robin.analysis.utilities.matkit import run_matkit

    work = str(tmp / f"work{k}")
    os.makedirs(os.path.join(work, "S"))
    cumulative = os.path.join(work, "S", "S.parquet")
    for i, group in enumerate(_batches(k)):
        bam = write_meth_bam(str(tmp / f"b{k}_{i}.bam"), world["ref"], group)
        per_bam = str(tmp / f"b{k}_{i}.parquet")
        run_matkit(bam, per_bam)
        merge_modkit_files([per_bam], cumulative, cumulative, world["panel"], "S", work, None, 1)
    table = pq.read_table(cumulative).to_pandas()
    sites = {}
    for row in table.itertuples():
        assert str(row.chrom) == CHROM
        sites[(int(row.chromStart), str(row.strand))] = (
            int(row.n_mod), int(row.n_canonical), int(row.valid_cov), float(row.percent_modified))
    return sites


def _assert_matches(robin, reference):
    assert sorted(robin) == sorted(reference)
    for key, (n_mod, n_can, valid, pct) in reference.items():
        r_mod, r_can, r_valid, r_pct = robin[key]
        assert (r_mod, r_can) == (n_mod, n_can), key
        assert r_valid == r_mod + r_can, key
        assert abs(r_pct - pct) <= 0.005, key
        assert abs(r_pct - 100 * r_mod / r_valid) <= 1e-4, key


@pytest.mark.parametrize("k", [1, 2, 4])
def test_the_cumulative_parquet_matches_modkit(world, k, tmp_path):
    _assert_matches(_robin_cumulative(world, k, tmp_path), world["combined"])


def test_splitting_into_batches_changes_nothing(world, tmp_path):
    one = _robin_cumulative(world, 1, tmp_path)
    assert _robin_cumulative(world, 2, tmp_path) == one
    assert _robin_cumulative(world, 4, tmp_path) == one


def test_the_mgmt_bed_matches_modkit_per_strand(world, tmp_path):
    """MGMT reads matkit's stranded 18-column BED; per strand-site it must match modkit."""
    from robin.analysis.utilities.matkit import run_matkit

    out = str(tmp_path / "mgmt.bed")
    run_matkit(world["union"], out)
    cpg_sites = {(p, "+") for p in CPGS} | {(p + 1, "-") for p in CPGS}
    robin = {}
    with open(out) as f:
        for line in f:
            col = line.split()
            key = (int(col[1]), col[5])
            if key in cpg_sites and int(col[9]) > 0:
                robin[key] = (int(col[11]), int(col[12]), int(col[9]), float(col[10]))
    _assert_matches(robin, world["stranded"])
