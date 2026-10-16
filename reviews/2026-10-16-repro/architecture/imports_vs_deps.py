import ast, os, sys, re, tomllib, collections
root = "./src/robin"
mods = collections.defaultdict(set)
for dp, dn, fn in os.walk(root):
    for f in fn:
        if f.endswith(".py"):
            p = os.path.join(dp, f)
            for node in ast.walk(ast.parse(open(p).read())):
                if isinstance(node, ast.Import):
                    for a in node.names: mods[a.name.split(".")[0]].add(os.path.relpath(p, root))
                elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                    mods[node.module.split(".")[0]].add(os.path.relpath(p, root))
std = set(sys.stdlib_module_names)
third = {m: v for m, v in mods.items() if m not in std and m != "robin"}
deps = tomllib.load(open("./pyproject.toml","rb"))["project"]["dependencies"]
names = [re.split(r"[<>=@!\[ ;]", d)[0].strip().lower().replace("_","-") for d in deps]
dup = [n for n, c in collections.Counter(names).items() if c > 1]
print("duplicate declared:", dup)
alias = {"sklearn":"scikit-learn","dateutil":"python-dateutil","cv2":"opencv-python","webview":"pywebview","zxingcpp":"zxing-cpp","tomli_w":"tomli-w","PyPDF2":"pypdf2","yaml":"pyyaml"}
imported = {alias.get(m, m).lower().replace("_","-"): sorted(v)[:3] for m, v in third.items()}
print("imported but not declared:", {k: v for k, v in imported.items() if k not in names})
print("declared but never imported:", sorted(set(names) - set(imported)))
