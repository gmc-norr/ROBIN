"""Simulate realistic accumulation: 500 panel CpGs, true methylation 80%, one read per site per batch.

Runs the real temp_utilities.merge_modkit_files for each batch and compares the stored
percent_modified with the coverage-weighted value (n_mod / valid_cov) that the same file holds.
Reports the fraction of sites NanoDX (>=60 -> methylated) and Sturgeon (>0.7 methylated) would call
methylated under each.
"""
import importlib.util
import logging
import os
import shutil
import tempfile

import numpy as np
import polars as pl
import pyarrow as pa
import pyarrow.parquet as pq

logging.disable(logging.CRITICAL)
spec = importlib.util.spec_from_file_location(
    "temp_utilities", "./src/robin/analysis/temp_utilities.py")
tu = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tu)

SCHEMA = pa.schema([
    ("chrom", pa.binary()), ("chromStart", pa.int64()), ("mod_code", pa.binary()),
    ("strand", pa.binary()), ("valid_cov", pa.uint32()), ("percent_modified", pa.float32()),
    ("n_mod", pa.uint32()), ("n_canonical", pa.uint32()),
])
rng = np.random.default_rng(1)
N_SITES, N_BATCH, P_TRUE = 500, 20, 0.8
pos = np.arange(N_SITES) * 1000 + 10_000

work = tempfile.mkdtemp(prefix="sim_")
panel = os.path.join(work, "panel_filter.txt")
with open(panel, "w") as f:
    f.write("chr start end IlmnID\n")
    for i, p in enumerate(pos):
        f.write(f"chr1\t{p}\t{p+2}\tcg{i}\n")
os.makedirs(os.path.join(work, "S"), exist_ok=True)
out = os.path.join(work, "S", "S.parquet")

for b in range(N_BATCH):
    nmod = rng.binomial(1, P_TRUE, N_SITES)
    t = pa.table({
        "chrom": [b"chr1"] * N_SITES, "chromStart": pos.astype(np.int64),
        "mod_code": [b"m"] * N_SITES, "strand": [b"+"] * N_SITES,
        "valid_cov": np.ones(N_SITES, dtype=np.uint32),
        "percent_modified": (100.0 * nmod).astype(np.float32),
        "n_mod": nmod.astype(np.uint32), "n_canonical": (1 - nmod).astype(np.uint32),
    }, schema=SCHEMA)
    p = os.path.join(work, "S", f"modkit_{b}.parquet")
    pq.write_table(t, p)
    tu.merge_modkit_files([p], out, out, panel, "S", work, None, 1)
    os.remove(p)

df = pl.read_parquet(out).with_columns((100 * pl.col("n_mod") / pl.col("valid_cov")).alias("weighted"))
print(df.select(pl.col("valid_cov").min().alias("cov_min"), pl.col("valid_cov").max().alias("cov_max")))
for col in ["percent_modified", "weighted"]:
    v = df[col].to_numpy()
    print(f"{col:17s} mean={v.mean():6.2f}  sd={v.std():6.2f}  "
          f"NanoDX methylated (>=60): {np.mean(v >= 60)*100:5.1f}%  "
          f"Sturgeon methylated (>70): {np.mean(v > 70)*100:5.1f}%  "
          f"Sturgeon no-call (30..70): {np.mean((v >= 30) & (v <= 70))*100:5.1f}%")
shutil.rmtree(work)
