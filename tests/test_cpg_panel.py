"""The CpG panel (parquet_filter.txt) is 0-based BED (METH-6), and its cache
can't outlive a code or panel change.

`parquet_filter.txt` lists each CpG as [C, C+2) in 0-based coordinates: three
probes read "CG" in hg38 at their listed starts. The merge subtracted 1 from
every start, so each CpG took in the base before it. Its PyRanges cache was
keyed by the file's basename only, so an old widened cache, or an old panel,
kept being used. Exact duplicate rows doubled the counts (METH-5).
Synthetic rows only.
"""

import os
import pickle

import pyarrow.parquet as pq
import pyranges as pr
import pandas as pd

from meth_fixtures import write_per_bam_parquet as _per_bam
from robin.analysis.temp_utilities import merge_modkit_files


def _panel(path, rows):
    with open(path, "w") as f:
        f.write("chr start end IlmnID\n")
        for start, end in rows:
            f.write(f"chr1\t{start}\t{end}\tcg{start}\n")
    return str(path)


def _merge(tmp_path, per_bam, panel):
    work = tmp_path / "work"
    (work / "S").mkdir(parents=True, exist_ok=True)
    out = str(work / "S" / "S.parquet")
    merge_modkit_files([per_bam], out, out, panel, "S", str(work), None, 1)
    return pq.read_table(out).to_pandas() if os.path.exists(out) else pd.DataFrame()


def test_the_base_before_a_cpg_is_not_kept(tmp_path):
    # 899 is the G of a CpG at 898, which is not in the panel.
    per_bam = _per_bam(tmp_path / "b.parquet", [(899, "-", 5, 0), (900, "+", 1, 1), (901, "-", 1, 1)])
    out = _merge(tmp_path, per_bam, _panel(tmp_path / "parquet_filter.txt", [(900, 902)]))
    assert set(out["chromStart"]) <= {900, 901}


def test_an_old_panel_cache_is_not_used(tmp_path):
    work_sample = tmp_path / "work" / "S"
    work_sample.mkdir(parents=True)
    old = pr.PyRanges(pd.DataFrame({"Chromosome": ["chr1"], "Start": [1299], "End": [1302]}))
    with open(work_sample / "parquet_filter.txt_1based.pgr_cache", "wb") as f:
        pickle.dump(old, f)
    per_bam = _per_bam(tmp_path / "b.parquet", [(900, "+", 2, 0), (1300, "+", 2, 0)])
    out = _merge(tmp_path, per_bam, _panel(tmp_path / "parquet_filter.txt", [(900, 902)]))
    assert set(out["chromStart"]) == {900}


def test_a_changed_panel_is_read_again(tmp_path):
    panel = tmp_path / "parquet_filter.txt"
    _merge(tmp_path, _per_bam(tmp_path / "a.parquet", [(900, "+", 1, 0)]), _panel(panel, [(900, 902)]))
    _panel(panel, [(900, 902), (1300, 1302)])
    os.utime(panel, (os.path.getmtime(panel) + 5,) * 2)
    out = _merge(tmp_path, _per_bam(tmp_path / "b.parquet", [(1300, "+", 1, 0)]), str(panel))
    assert 1300 in set(out["chromStart"])


def test_a_duplicated_panel_row_does_not_double_the_counts(tmp_path):
    per_bam = _per_bam(tmp_path / "b.parquet", [(900, "+", 3, 1)])
    out = _merge(tmp_path, per_bam, _panel(tmp_path / "parquet_filter.txt", [(900, 902), (900, 902)]))
    assert list(out["n_mod"]) == [3]
    assert list(out["n_canonical"]) == [1]


def test_a_negative_panel_row_is_harmless(tmp_path):
    per_bam = _per_bam(tmp_path / "b.parquet", [(900, "+", 1, 1)])
    out = _merge(tmp_path, per_bam, _panel(tmp_path / "parquet_filter.txt", [(-1, 0), (900, 902)]))
    assert set(out["chromStart"]) == {900}
