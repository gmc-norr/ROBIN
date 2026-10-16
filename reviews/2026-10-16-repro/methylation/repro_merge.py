"""Reproduce merge_modkit_files data-handling behaviour with synthetic per-BAM parquets.

Scenarios:
 A. Mean-of-means / exponential weighting across batches.
 B. Duplicate rows from the inner join when a per-BAM file has >1 row at (chrom, chromStart).
 C. Overlapping panel intervals -> intersect duplicates -> double counting.
 D. First-batch path writes un-grouped rows; later path groups -> inconsistent schema/semantics.
"""
import importlib.util
import logging
import os
import shutil
import sys
import tempfile

import pyarrow as pa
import pyarrow.parquet as pq
import polars as pl
import pandas as pd

logging.basicConfig(level=logging.WARNING)

SRC = "./src/robin/analysis/temp_utilities.py"
spec = importlib.util.spec_from_file_location("temp_utilities", SRC)
tu = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tu)

SCHEMA = pa.schema([
    ("chrom", pa.binary()),
    ("chromStart", pa.int64()),
    ("mod_code", pa.binary()),
    ("strand", pa.binary()),
    ("valid_cov", pa.uint32()),
    ("percent_modified", pa.float32()),
    ("n_mod", pa.uint32()),
    ("n_canonical", pa.uint32()),
])


def write_bam_parquet(path, rows):
    cols = {k: [] for k in SCHEMA.names}
    for chrom, pos, strand, cov, nmod in rows:
        cols["chrom"].append(chrom.encode())
        cols["chromStart"].append(pos)
        cols["mod_code"].append(b"m")
        cols["strand"].append(strand.encode())
        cols["valid_cov"].append(cov)
        cols["percent_modified"].append(round(100.0 * nmod / cov, 2))
        cols["n_mod"].append(nmod)
        cols["n_canonical"].append(cov - nmod)
    pq.write_table(pa.table(cols, schema=SCHEMA), path)


def run_batches(panel_path, batches, label):
    work = tempfile.mkdtemp(prefix=f"merge_{label}_")
    sample = "S1"
    os.makedirs(os.path.join(work, sample), exist_ok=True)
    out = os.path.join(work, sample, f"{sample}.parquet")
    for i, rows in enumerate(batches):
        p = os.path.join(work, sample, f"modkit_{i}.parquet")
        write_bam_parquet(p, rows)
        tu.merge_modkit_files([p], out, out, panel_path, sample, work, None, 1)
    df = pl.read_parquet(out).sort(["chrom", "chromStart"])
    shutil.rmtree(work)
    return df


def write_panel_txt(path, rows):
    with open(path, "w") as f:
        f.write("chr start end IlmnID\n")
        for r in rows:
            f.write("\t".join(map(str, r)) + "\n")


tmpd = tempfile.mkdtemp(prefix="panel_")
panel = os.path.join(tmpd, "panel_filter.txt")
# One CpG probe: C at 0-based 1000, G at 1001 (file stores 0-based [1000,1002) like parquet_filter.txt)
write_panel_txt(panel, [("chr1", 1000, 1002, "cgA")])

print("=== A. mean-of-means weighting across batches ===")
# Batch 1..9: 10 reads each, fully methylated (100%). Batch 10: 1 read, unmethylated (0%).
batches = [[("chr1", 1000, "+", 10, 10)] for _ in range(9)] + [[("chr1", 1000, "+", 1, 0)]]
df = run_batches(panel, batches, "A")
print(df)
r = df.row(0, named=True)
print(f"stored percent_modified={r['percent_modified']:.2f}  "
      f"coverage-weighted truth={100*r['n_mod']/r['valid_cov']:.2f}  (n_mod={r['n_mod']}, valid_cov={r['valid_cov']})")

print("\n=== B. duplicate (chrom,chromStart) rows in one per-BAM file -> join explosion ===")
# Same position, two strands (possible because matkit does no reference check)
batches = [[("chr1", 1000, "+", 5, 5), ("chr1", 1000, "-", 3, 0)]]
df = run_batches(panel, batches, "B1")
print("after FIRST batch (no group_by path):")
print(df)
batches2 = batches + [[("chr1", 1000, "+", 1, 1)]]
df = run_batches(panel, batches2, "B2")
print("after SECOND batch (group_by path): true valid_cov=9")
print(df)

print("\n=== C. overlapping panel intervals (adjacent probes / duplicate rows) -> double counting ===")
panel2 = os.path.join(tmpd, "panel_overlap.txt")
# CGCG: probe A C at 1000, probe B C at 1002 (0-based). After Start-1 widening: A=[999,1002) B=[1001,1004)
write_panel_txt(panel2, [("chr1", 1000, 1002, "cgA"), ("chr1", 1002, 1004, "cgB")])
batches = [[("chr1", 1001, "-", 4, 4)], [("chr1", 1001, "-", 4, 4)]]
df = run_batches(panel2, batches, "C")
print("two batches of 4 reads each at chr1:1001(-); true valid_cov=8")
print(df)

panel3 = os.path.join(tmpd, "panel_dup.txt")
write_panel_txt(panel3, [("chr1", 1000, 1002, "cgA"), ("chr1", 1000, 1002, "cgA_dup")])
batches = [[("chr1", 1000, "+", 4, 4)], [("chr1", 1000, "+", 4, 4)]]
df = run_batches(panel3, batches, "C2")
print("exact duplicate panel rows; two batches of 4 reads; true valid_cov=8")
print(df)

print("\n=== E. what the .txt Start-1 conversion does to a 0-based panel row ===")
cache = pr_df = None
import pyranges as pr
bed_df = pd.read_csv(panel, sep=r"\s+", header=0, dtype=str).rename(columns={"chr": "Chromosome", "start": "Start", "end": "End"})
bed_df["Start"] = bed_df["Start"].astype(int) - 1
bed_df["End"] = bed_df["End"].astype(int)
print(bed_df)
batches = [[("chr1", 999, "-", 4, 4), ("chr1", 1000, "+", 4, 0), ("chr1", 1001, "-", 4, 0)]]
df = run_batches(panel, batches, "E")
print("sites at 999(-, neighbouring CpG's minus strand), 1000(+), 1001(-): kept rows ->")
print(df)

shutil.rmtree(tmpd)
