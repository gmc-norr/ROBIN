import ast, textwrap
src=open("./src/robin/reporting/sections/variants.py").read()
tree=ast.parse(src)
fn=[n for c in tree.body if isinstance(c,ast.ClassDef) and c.name=="VariantAnalysis" for n in c.body if isinstance(n,ast.FunctionDef) and n.name=="_is_pathogenic"][0]
code="import logging\nlogger=logging.getLogger('x')\n"+ast.get_source_segment(src,fn)
code=textwrap.dedent(code.replace("    def _is_pathogenic","def _is_pathogenic",1))
# dedent function body
lines=ast.get_source_segment(src,fn).splitlines()
body="\n".join(l[4:] if l.startswith("    ") else l for l in lines)
ns={}
exec("import logging\nlogger=logging.getLogger('x')\n"+body, ns)
f=ns["_is_pathogenic"]
cases={
 "conflicting, mostly benign":"CLNSIG=Conflicting_classifications_of_pathogenicity;CLNSIGCONF=Benign(12)|Likely_benign(3)|Uncertain_significance(1)",
 "benign":"CLNSIG=Benign",
 "pathogenic":"CLNSIG=Pathogenic",
}
for k,v in cases.items(): print(k, "->", f(None, v))
