import ast, os, tempfile, types, numpy as np, pandas as pd
src=open("./src/robin/reporting/sections/coverage.py").read()
tree=ast.parse(src)
cls=[n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=="CoverageSection"][0]
fn=[n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=="_initialize_data"][0]
body="\n".join((l[4:] if l.startswith("    ") else l) for l in ast.get_source_segment(src,fn).splitlines())
ns={"os":os,"pd":pd,"np":np}; exec(body,ns); init=ns["_initialize_data"]
def run(files):
    d=tempfile.mkdtemp()
    for name,content in files.items(): open(os.path.join(d,name),"w").write(content)
    self=types.SimpleNamespace(report=types.SimpleNamespace(output=d))
    try: init(self); print("OK", getattr(self,"coverage_data",None))
    except Exception as e: print(type(e).__name__, str(e)[:90])
cov="#rname,startpos,endpos,numreads,covbases,coverage,meandepth,meanbaseq,meanmapq\nchr1,1,1000,5,500,50,0.5,20,60\n"
bed="chrom,startpos,endpos,name,bases,length\nchr1,1,100,EGFR,1000,100\n"
print("1) coverage_main.csv written, others not yet:"); run({"coverage_main.csv":cov})
print("2) target_coverage.csv header only:"); run({"coverage_main.csv":cov,"bed_coverage_main.csv":bed,"target_coverage.csv":"chrom,startpos,endpos,name,coverage\n"})
print("3) target_coverage.csv truncated to 0 bytes mid-write:"); run({"coverage_main.csv":cov,"bed_coverage_main.csv":bed,"target_coverage.csv":""})
