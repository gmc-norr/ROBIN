"""The cumulative methylation parquet holds counts, summed per CpG (METH-1,
METH-2, METH-4).

The merge averaged each batch's percent with the running percent (so the
newest batch weighed at least half), wrote the first batch without
aggregating it, and kept a CpG's + and - strand rows apart. Now every write
sums n_mod and n_canonical per CpG (the - row moved onto the C), and percent
is 100 * n_mod / (n_mod + n_canonical). A file written the old way is
migrated on its next merge. Synthetic rows only.
"""

import json
import os

import polars as pl
import pytest

from meth_fixtures import write_per_bam_parquet
from robin.analysis.temp_utilities import merge_modkit_files


@pytest.fixture
def setup(tmp_path):
    panel = tmp_path / "parquet_filter.txt"
    panel.write_text("chr start end IlmnID\nchr1\t900\t902\tcg900\nchr1\t1300\t1302\tcg1300\n")
    work = tmp_path / "work"
    (work / "S").mkdir(parents=True)
    out = str(work / "S" / "S.parquet")

    def merge(files):
        merge_modkit_files(files, out, out, str(panel), "S", str(work), None, len(files))
        return pl.read_parquet(out)

    return tmp_path, out, merge


def _row(df, pos):
    rows = df.filter(pl.col("chromStart") == pos)
    assert rows.height == 1, rows
    return rows.row(0, named=True)


def test_a_late_batch_does_not_outweigh_the_run(setup):
    # The code review's case: 50 reads at 80%, then 1 read at 0%.
    tmp, out, merge = setup
    merge([write_per_bam_parquet(tmp / "a.parquet", [(900, "+", 40, 10)])])
    row = _row(merge([write_per_bam_parquet(tmp / "b.parquet", [(900, "+", 0, 1)])]), 900)
    assert (row["n_mod"], row["n_canonical"], row["valid_cov"]) == (40, 11, 51)
    assert row["percent_modified"] == pytest.approx(100 * 40 / 51, abs=1e-4)


def test_the_first_write_is_summed_too(setup):
    tmp, out, merge = setup
    files = [write_per_bam_parquet(tmp / f"{i}.parquet", [(900, "+", 3, 1)]) for i in range(2)]
    row = _row(merge(files), 900)
    assert (row["n_mod"], row["n_canonical"], row["valid_cov"]) == (6, 2, 8)


def test_both_strands_of_a_cpg_are_one_row(setup):
    tmp, out, merge = setup
    df = merge([write_per_bam_parquet(tmp / "a.parquet", [(900, "+", 3, 0), (901, "-", 0, 1)])])
    assert df.height == 1
    row = _row(df, 900)
    assert (row["strand"], row["n_mod"], row["n_canonical"]) == (".", 3, 1)
    assert row["percent_modified"] == pytest.approx(75.0)


def test_an_old_format_file_is_migrated(setup):
    """A running sample's parquet from before this fix: stranded rows, failed
    calls in valid_cov, mean percents. Its counts are kept and re-summed."""
    tmp, out, merge = setup
    pl.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1"],
            "chromStart": [899, 900, 901],  # 899: the widened panel's extra base
            "mod_code": ["m", "m", "m"],
            "strand": ["-", "+", "-"],
            "valid_cov": [9, 12, 6],
            "percent_modified": [11.0, 41.0, 50.0],
            "n_mod": [1, 4, 3],
            "n_canonical": [8, 6, 3],
        },
        schema_overrides={"chromStart": pl.Int64, "valid_cov": pl.UInt32, "n_mod": pl.UInt32,
                          "n_canonical": pl.UInt32, "percent_modified": pl.Float32},
    ).write_parquet(out)
    with open(out.replace(".parquet", "_metadata.json"), "w") as f:
        json.dump({"bam_file_count": 3, "column_format": "optimized"}, f)

    df = merge([write_per_bam_parquet(tmp / "new.parquet", [(900, "+", 1, 0)])])
    assert df.height == 1
    row = _row(df, 900)
    assert (row["n_mod"], row["n_canonical"], row["valid_cov"]) == (8, 9, 17)
    assert row["percent_modified"] == pytest.approx(100 * 8 / 17, abs=1e-4)
    with open(out.replace(".parquet", "_metadata.json")) as f:
        meta = json.load(f)
    assert meta["bam_file_count"] == 4
    assert meta["column_format"] != "optimized"
