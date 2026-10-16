"""Reproduce how the cumulative parquet feeds NanoDX / Random Forest / Sturgeon inputs.

Uses robin's own functions (merge_bedmethyl.py) with a stub `sturgeon` module so the
module imports without the real dependency (only read_probes_file is imported there).
"""
import importlib.util
import os
import sys
import tempfile
import types

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

# --- stub sturgeon.utils.read_probes_file (unused by the functions under test)
st = types.ModuleType("sturgeon")
stu = types.ModuleType("sturgeon.utils")
stu.read_probes_file = lambda *a, **k: None
st.utils = stu
sys.modules["sturgeon"] = st
sys.modules["sturgeon.utils"] = stu

SRC = "./src/robin/analysis/utilities/merge_bedmethyl.py"
spec = importlib.util.spec_from_file_location("merge_bedmethyl", SRC)
mb = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mb)

tmp = tempfile.mkdtemp(prefix="clf_")

# Cumulative parquet as produced by matkit + merge for ONE CpG (cg13869341, C at 0-based 15864):
#   forward-strand reads -> row at 15864 '+', reverse-strand reads -> row at 15865 '-'
cum = pd.DataFrame({
    "chrom": pd.Categorical(["chr1", "chr1"]),
    "chromStart": np.array([15864, 15865], dtype=np.int64),
    "mod_code": pd.Categorical(["m", "m"]),
    "strand": pd.Categorical(["+", "-"]),
    "percent_modified": np.array([100.0, 0.0], dtype=np.float32),
    "valid_cov": np.array([3, 20], dtype=np.uint32),
    "n_mod": np.array([3, 0], dtype=np.uint32),
    "n_canonical": np.array([0, 20], dtype=np.uint32),
})
pqf = os.path.join(tmp, "S1.parquet")
cum.to_parquet(pqf)

df = mb.load_minimal_modkit_data(pqf)
print("load_minimal_modkit_data columns:", list(df.columns), "(valid_cov/n_mod dropped)")

print("\n=== NanoDX path (nanodx_analysis.process_parquet_file steps 3-5) ===")
nanodx_df = mb.collapse_minimal_bedmethyl(df.copy())
cpgs = pd.read_csv("./src/robin/resources/hglft_genome_260e9_91a970_clean.bed", sep="\t", header=None)
test_df = pd.merge(nanodx_df, cpgs, left_on=["chrom", "start_pos"], right_on=[0, 1])
test_df = test_df.rename(columns={3: "probe_id", "fraction": "methylation_call"})
print(nanodx_df[["chrom", "start_pos", "strand", "fraction"]])
print("rows reaching NN after merge with hglft bed:")
print(test_df[["chrom", "start_pos", "probe_id", "methylation_call"]])
print("-> 20x reverse-strand coverage (0% methylated) discarded; probe called from 3 forward reads only")

print("\n=== Random Forest path (random_forest_analysis.process_parquet_file steps 3-7) ===")
full = mb.reconstruct_full_bedmethyl_data(df.copy())
forest = mb.collapse_bedmethyl(full)
print(forest[["chrom", "start_pos", "strand", "Nvalid", "score"]])
print("-> cov column written to R is the constant 10, not real coverage (3 / 20)")

# First-batch duplicates (un-grouped write path) -> RF score summed
dup = pd.concat([cum.iloc[[0]], cum.iloc[[0]]], ignore_index=True)
pqf2 = os.path.join(tmp, "S2.parquet")
dup.to_parquet(pqf2)
d2 = mb.load_minimal_modkit_data(pqf2)
forest2 = mb.collapse_bedmethyl(mb.reconstruct_full_bedmethyl_data(d2.copy()))
print("duplicated row in first-batch parquet -> RF methylation_percent (score):")
print(forest2[["chrom", "start_pos", "Nvalid", "score"]])

print("\n=== Sturgeon path: modkit_pileup_file_to_bed with probe at 15864, margin 25 ===")
probes = os.path.join(tmp, "probes.bed")
with open(probes, "w") as f:
    f.write("chr start end ID_REF\nchr1 15864 15866 cg13869341\n")
res = mb.modkit_pileup_file_to_bed(df.copy(), os.path.join(tmp, "out.bed"), probes)
print(res)

print("\n=== percent_modified sensitivity: Sturgeon 0.3/0.7, NanoDX 60 thresholds ===")
for pm in [98.9, 50.0]:
    one = pd.DataFrame({"chrom": ["chr1"], "chromStart": [15864], "percent_modified": [pm], "mod_code": ["m"], "strand": ["+"]})
    r = mb.modkit_pileup_file_to_bed(one, os.path.join(tmp, "o2.bed"), probes)
    print(f"percent_modified={pm}: sturgeon rows={len(r)} methylation_call={list(r['methylation_call']) if len(r) else None}; "
          f"nanodx call={'1' if pm >= 60 else '-1'}")
