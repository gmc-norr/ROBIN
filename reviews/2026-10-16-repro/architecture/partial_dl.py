import sys, tempfile, json, importlib.util, types
from pathlib import Path
SRC = "./src/robin/utils/"
def load(name, fn):
    spec = importlib.util.spec_from_file_location(name, SRC + fn); m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m; spec.loader.exec_module(m); return m
# provide robin.utils.model_checker for the lazy import
pkg = types.ModuleType("robin"); pkg.__path__ = []; sys.modules["robin"] = pkg
up = types.ModuleType("robin.utils"); up.__path__ = []; sys.modules["robin.utils"] = up
mc = load("robin.utils.model_checker", "model_checker.py")
mu = load("robin.utils.model_updater", "model_updater.py")
d = Path(tempfile.mkdtemp())
manifest = json.load(open("./assets.json"))
mp = d / "assets.json"; mp.write_text(json.dumps(manifest))
calls = {"n": 0}
def fake_download(url, target_path, token, *, label=""):
    calls["n"] += 1
    with open(target_path, "wb") as f:
        f.write(b"\x80\x04truncated-pickle")   # partial content written ...
    raise ConnectionResetError("network dropped mid-transfer")
mu._download = fake_download
ok, msgs = mu.update_models(models_dir=d, manifest_path=str(mp))
print("run1 ok=", ok, msgs[-1][:80])
print("partial file left:", [(p.name, p.stat().st_size) for p in d.iterdir() if p.suffix in (".zip", ".pkl")])
# simulate 'retry' succeeding for the other files only
def good_download(url, target_path, token, *, label=""):
    raise AssertionError("should not be called for files that already exist")
for key, fn in mc.get_required_models():
    if not (d / fn).exists():
        (d / fn).write_bytes(b"x")  # pretend other files exist
mu._download = good_download
ok2, msgs2 = mu.update_models(models_dir=d, manifest_path=str(mp))
print("run2 ok=", ok2, [m for m in msgs2 if "Skipping" in m][:1])
mc.get_models_directory = lambda project_root=None: d
print("startup check_model_files:", mc.check_model_files())
