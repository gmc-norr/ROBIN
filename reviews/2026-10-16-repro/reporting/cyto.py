import ast, numpy as np, pandas as pd, logging
src=open("./src/robin/reporting/sections/cnv.py").read()
tree=ast.parse(src)
fn=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="analyze_cytoband_cnv"][0]
ns={"np":np,"pd":pd,"logger":logging.getLogger("x")}
exec(ast.get_source_segment(src,fn), ns)
f=ns["analyze_cytoband_cnv"]
rng=np.random.default_rng(0)
bw=1_000_000
cnv={}
for c in range(1,23):
    cnv[f"chr{c}"]=rng.normal(0,0.05,100)
cnv["chr1"][40:50]+=0.12   # small +0.12 shift in one band, far below the 0.4 gain threshold
cyto=pd.DataFrame([{"chrom":"chr1","start_pos":i*10_000_000,"end_pos":(i+1)*10_000_000-1,"name":f"b{i}","stain":"x"} for i in range(10)])
cent=pd.DataFrame([{"chrom":"chr1","start_pos":0,"end_pos":0,"name":"c"}])
genes=pd.DataFrame([{"chrom":"chr1","start_pos":41_000_000,"end_pos":42_000_000,"gene":"GENEX"}])
out=f(cnv,"chr1",{"bin_width":bw},cyto,cent,genes,"XX")
print(out[["name","mean_cnv","cnv_state","genes"]].to_string())
