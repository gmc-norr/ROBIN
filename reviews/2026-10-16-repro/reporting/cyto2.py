import ast, numpy as np, pandas as pd, logging
src=open("./src/robin/reporting/sections/cnv.py").read()
tree=ast.parse(src)
fn=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=="analyze_cytoband_cnv"][0]
ns={"np":np,"pd":pd,"logger":logging.getLogger("x")}
exec(ast.get_source_segment(src,fn), ns)
f=ns["analyze_cytoband_cnv"]
cyto=pd.read_csv("./src/robin/resources/cytoBand.txt",sep="\t",names=["chrom","start_pos","end_pos","name","stain"])
cent=pd.read_csv("./src/robin/resources/cenSatRegions.bed",usecols=[0,1,2,3],names=["chrom","start_pos","end_pos","name"],header=None,sep=r"\s+")
genes=pd.read_csv("./src/robin/resources/unique_genes.bed",sep="\t",names=["chrom","start_pos","end_pos","gene"])
rng=np.random.default_rng(1)
bw=1_000_000
lens=cyto.groupby("chrom")["end_pos"].max()
# perfectly normal diploid genome: log-ratio-like noise around 0, smoothed with a 3-bin moving average (as the report does)
def uniform_filter1d(a,n,mode=None):
    p=np.pad(a,(n//2,n-1-n//2),mode="edge"); return np.convolve(p,np.ones(n)/n,mode="valid")
cnv={c:uniform_filter1d(rng.normal(0,0.15,int(lens[c]//bw)+1),3,mode="nearest") for c in [f"chr{i}" for i in range(1,23)]+["chrX"] if c in lens}
g=set(); l=set(); nchr=0
for c in cnv:
    out=f(cnv,c,{"bin_width":bw},cyto,cent,genes,"XX")
    if not out.empty:
        nchr+=1
        for _,r in out[out.cnv_state=="GAIN"].iterrows(): g.update(r["genes"])
        for _,r in out[out.cnv_state=="LOSS"].iterrows(): l.update(r["genes"])
print("chromosomes with GAIN/LOSS calls (=> per-chr plots shown):", nchr, "of", len(cnv))
print("Genes in Gained Regions:", len(g), " Genes in Lost Regions:", len(l))
print("max |bin| value:", max(abs(v).max() for v in cnv.values()), "(fixed GUI/central threshold is 0.4)")
