"""matkit's per-read calls, checked without modkit (so they run anywhere).

From the final review of step 5:
- An implicit MM group (``C+m.`` or ``C+m``) means the Cs it skips are
  unmethylated, but only CpG Cs matter: modkit pileup --cpg reports CpGs only,
  and a basecaller in ``?`` mode calls only there. Adding every C made a site
  of every C on every read (gigabytes on a real BAM) and gave non-CpG panel
  probes a 0% row.
- A record with MM but no ML has no probabilities; modkit drops it.
Synthetic reads only.
"""

import os

import pyarrow.parquet as pq
import pytest

from meth_fixtures import (
    CANON, MOD, ReadSpec, read_reference, write_meth_bam, write_panel, write_reference,
)
from robin.analysis.temp_utilities import merge_modkit_files
from robin.analysis.utilities.matkit import run_matkit

CPGS = [500, 900, 1300]
S, E = 400, 1400


def _ref(tmp_path) -> str:
    path = os.path.join(str(tmp_path), "ref.fa")
    return path if os.path.exists(path) else write_reference(str(tmp_path), CPGS)


def _per_bam(tmp_path, reads, name):
    bam = write_meth_bam(str(tmp_path / f"{name}.bam"), _ref(tmp_path), reads)
    out = str(tmp_path / f"{name}.parquet")
    run_matkit(bam, out)
    table = pq.read_table(out).to_pandas()
    return {
        (int(r.chromStart), r.strand.decode() if isinstance(r.strand, bytes) else str(r.strand)):
            (int(r.n_mod), int(r.n_canonical))
        for r in table.itertuples()
    }


@pytest.mark.parametrize("mode", [".", ""])
def test_an_implicit_tag_gives_the_calls_an_explicit_one_does(tmp_path, mode):
    implicit = [
        ReadSpec("f", "+", S, E, m={500: MOD}, mode=mode),
        ReadSpec("r", "-", S, E, m={900: MOD}, mode=mode),
    ]
    explicit = [
        ReadSpec("f", "+", S, E, m={500: MOD, 900: CANON, 1300: CANON}),
        ReadSpec("r", "-", S, E, m={500: CANON, 900: MOD, 1300: CANON}),
    ]
    want = {
        (500, "+"): (1, 0), (900, "+"): (0, 1), (1300, "+"): (0, 1),
        (501, "-"): (0, 1), (901, "-"): (1, 0), (1301, "-"): (0, 1),
    }
    assert _per_bam(tmp_path, explicit, "explicit") == want
    got = _per_bam(tmp_path, implicit, "implicit")
    assert len(got) == len(want), f"{len(got)} sites, not {len(want)}"
    assert got == want


def test_a_record_with_mm_but_no_ml_gives_no_calls(tmp_path):
    reads = [
        ReadSpec("q", "+", S, E, m={500: MOD, 900: MOD}, drop_ml=True),
        ReadSpec("d", "-", S, E, m={900: MOD}, mode=".", drop_ml=True),
    ]
    assert _per_bam(tmp_path, reads, "noml") == {}


def test_a_non_cpg_panel_row_gets_no_row(tmp_path):
    ref = read_reference(_ref(tmp_path))
    ch = next(i for i in range(S + 50, E - 50) if ref[i] == "C" and ref[i + 1] != "G")
    reads = [ReadSpec("f", "+", S, E, m={500: MOD}, mode=".")]
    bam = write_meth_bam(str(tmp_path / "f.bam"), _ref(tmp_path), reads)
    per_bam = str(tmp_path / "f.parquet")
    run_matkit(bam, per_bam)
    panel = write_panel(str(tmp_path / "parquet_filter.txt"), [(500, 502, "cg500"), (ch, ch + 1, "ch1")])
    work = tmp_path / "work"
    (work / "S").mkdir(parents=True)
    out = str(work / "S" / "S.parquet")
    merge_modkit_files([per_bam], out, out, panel, "S", str(work), None, 1)
    assert set(pq.read_table(out).to_pandas()["chromStart"]) == {500}
