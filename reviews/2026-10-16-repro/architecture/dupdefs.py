import ast, sys, os, collections, hashlib
root = "./src/robin"
defs = collections.defaultdict(list)
for dp, dn, fn in os.walk(root):
    if "submodules" in dp: continue
    for f in fn:
        if not f.endswith(".py"): continue
        p = os.path.join(dp, f)
        try: tree = ast.parse(open(p).read())
        except SyntaxError as e: print("SYNTAX", p, e); continue
        for node in tree.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
                src = ast.dump(node)
                defs[node.name].append((os.path.relpath(p, root), node.lineno, len(ast.unparse(node).splitlines()), hashlib.md5(src.encode()).hexdigest()[:8]))
for name, locs in sorted(defs.items()):
    files = {l[0] for l in locs}
    if len(files) > 1 and not name.startswith("__"):
        ident = len({l[3] for l in locs}) == 1
        print(f"{name}: {'IDENTICAL' if ident else 'differs'} " + "; ".join(f"{l[0]}:{l[1]}({l[2]}L)" for l in locs))
