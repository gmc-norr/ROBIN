"""Time the per-BAM filter step of temp_utilities.merge_modkit_files on a realistic per-BAM parquet
(50k reads x ~8 kb -> ~4M CpG-strand rows genome-wide) against the real 900k-probe panel."""
import importlib.util, logging, os, tempfile, shutil, time
import numpy as np, pyarrow as pa, pyarrow.parquet as pq, polars as pl
logging.disable(logging.CRITICAL)
spec = importlib.util.spec_from_file_location("tu", "./src/robin/analysis/temp_utilities.py")
tu = importlib.util.module_from_spec(spec); spec.loader.exec_module(tu)
S = pa.schema([("chrom", pa.binary()), ("chromStart", pa.int64()), ("mod_code", pa.binary()), ("strand", pa.binary()),
               ("valid_cov", pa.uint32()), ("percent_modified", pa.float32()), ("n_mod", pa.uint32()), ("n_canonical", pa.uint32())])
rng = np.random.default_rng(0)
N = 4_000_000
chroms = np.array([f"chr{i}".encode() for i in range(1, 23)], dtype=object)
c = chroms[rng.integers(0, 22, N)]
pos = rng.integers(0, 150_000_000, N)
nmod = rng.integers(0, 2, N).astype(np.uint32)
w = tempfile.mkdtemp(); os.makedirs(f"{w}/S")
p = f"{w}/S/modkit_0.parquet"
pq.write_table(pa.table({"chrom": list(c), "chromStart": pos, "mod_code": [b"m"] * N, "strand": [b"+"] * N,
    "valid_cov": np.ones(N, np.uint32), "percent_modified": (nmod * 100).astype(np.float32), "n_mod": nmod,
    "n_canonical": (1 - nmod).astype(np.uint32)}, schema=S), p)
panel = "./src/robin/resources/parquet_filter.txt"
# warm the PyRanges cache so we time only the per-file work
tu.merge_modkit_files([], f"{w}/S/S.parquet", f"{w}/S/S.parquet", panel, "S", w, None, 1)
t = time.time()
df = pl.read_parquet(p, columns=["chrom", "mod_code", "strand"])
for col in ["chrom", "mod_code", "strand"]:
    df.with_columns(pl.col(col).map_elements(lambda x: x.decode("utf-8", errors="replace"), return_dtype=pl.Utf8))
t_decode = time.time() - t
t = time.time()
tu.merge_modkit_files([p], f"{w}/S/S.parquet", f"{w}/S/S.parquet", panel, "S", w, None, 1)
t_merge = time.time() - t
t = time.time(); df.with_columns([pl.col(x).cast(pl.Utf8) for x in ["chrom", "mod_code", "strand"]]); t_cast = time.time() - t
print(f"per-BAM merge_modkit_files: {t_merge:.1f}s  (of which map_elements decode ~{t_decode:.1f}s; vectorised cast: {t_cast:.2f}s)")
print("rows kept:", pl.read_parquet(f"{w}/S/S.parquet").height)
shutil.rmtree(w)
