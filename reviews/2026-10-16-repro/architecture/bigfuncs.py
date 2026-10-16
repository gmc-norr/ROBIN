import ast, os
root = "./src/robin"
out = []
tot_exc = 0; bare = 0; files = 0; total = 0
for dp, dn, fn in os.walk(root):
    for f in fn:
        if not f.endswith(".py"): continue
        p = os.path.join(dp, f); src = open(p).read(); total += src.count("\n")
        tree = ast.parse(src); files += 1
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                out.append((node.end_lineno - node.lineno + 1, os.path.relpath(p, root), node.name, node.lineno))
            if isinstance(node, ast.ExceptHandler):
                if node.type is None or (isinstance(node.type, ast.Name) and node.type.id in ("Exception","BaseException")):
                    tot_exc += 1
                    if len(node.body) == 1 and isinstance(node.body[0], (ast.Pass, ast.Continue)):
                        bare += 1
out.sort(reverse=True)
for r in out[:12]: print(*r)
print("files", files, "lines", total, "broad except handlers", tot_exc, "of which body is just pass/continue", bare)
