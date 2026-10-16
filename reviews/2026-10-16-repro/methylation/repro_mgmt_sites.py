"""extract_mgmt_site_rows_from_bed: percent-vs-fraction auto-detection on a low-methylation, deep MGMT bed."""
import importlib.util, sys, types, tempfile, os
for n in ["pysam"]:
    sys.modules[n] = types.ModuleType(n)
lc = types.ModuleType("robin.logging_config"); lc.get_job_logger = lambda *a, **k: None
r = types.ModuleType("robin"); r.__path__ = []; r.__file__ = "/tmp/x/robin/__init__.py"
sys.modules["robin"] = r; sys.modules["robin.logging_config"] = lc
spec = importlib.util.spec_from_file_location("mgmt", "./src/robin/analysis/mgmt_analysis.py")
mg = importlib.util.module_from_spec(spec); spec.loader.exec_module(mg)

def line(start, strand, cov, pct):
    # matkit write_bedmethyl_improved row (percent 0-100), then convert_to_mixed_delim merges col10/11
    f = ["chr10", str(start), str(start + 1), "m", str(cov), strand, str(start), str(start + 1), "255,0,0",
         f"{cov} {pct:.2f}", f"{pct:.2f}", "1", str(cov - 1), "0", "0", "0", "0", "0"]
    return "\t".join(f) + "\n"

d = tempfile.mkdtemp(); p = os.path.join(d, "1_mgmt.bed")
with open(p, "w") as fh:
    # unmethylated promoter, deep coverage: one methylated read out of 120 at site 1 (+), everything else 0%
    fh.write(line(129467254, "+", 120, 100 / 120))   # site 1 forward C (1-based 129467255)
    fh.write(line(129467255, "-", 110, 0.0))         # site 1 reverse C
    fh.write(line(129467257, "+", 118, 0.0))
    fh.write(line(129467258, "-", 105, 0.0))
for row in mg.extract_mgmt_site_rows_from_bed(p):
    print(row["site"], "meth%=", row["meth"], "cov=", row["cov_total"])
print("true site-1 methylation = 1/230 reads = 0.43%")
