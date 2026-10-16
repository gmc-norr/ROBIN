import importlib.util, os, sys, tempfile, warnings
warnings.simplefilter("default")
import polars as pl
spec = importlib.util.spec_from_file_location("tu", "./src/robin/analysis/temp_utilities.py")
tu = importlib.util.module_from_spec(spec); spec.loader.exec_module(tu)
d = tempfile.mkdtemp()
flt = os.path.join(d, "f.bed")
open(flt, "w").write("chr1\t100\t101\tcg1\n")
def batch(name, valid, nmod):
    p = os.path.join(d, name)
    pl.DataFrame({"chrom":["chr1"],"chromStart":pl.Series([100],dtype=pl.UInt32),"mod_code":["m"],"strand":["+"],
                  "valid_cov":pl.Series([valid],dtype=pl.UInt32),"percent_modified":pl.Series([100.0*nmod/valid],dtype=pl.Float32),
                  "n_mod":pl.Series([nmod],dtype=pl.UInt32),"n_canonical":pl.Series([valid-nmod],dtype=pl.UInt32)}).write_parquet(p)
    return p
out = os.path.join(d, "S", "S.parquet")
os.makedirs(os.path.dirname(out), exist_ok=True)
# batch1: 50 reads, 40 methylated (80%); batch2: 1 read, 0 methylated (0%)
tu.merge_modkit_files([batch("b1.parquet",50,40), batch("b3.parquet",10,0)], out, out, flt, "S", d, None, 2); print("after first call:", pl.read_parquet(out).shape)
tu.merge_modkit_files([batch("b2.parquet",1,0)], out, out, flt, "S", d, None, 1)
r = pl.read_parquet(out)
print(r)
print("percent_modified stored:", r["percent_modified"][0], " correct (100*n_mod/valid_cov):", 100*r["n_mod"][0]/r["valid_cov"][0])
