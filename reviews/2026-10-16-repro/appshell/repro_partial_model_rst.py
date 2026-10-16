"""Reproduce: interrupted model download leaves a truncated file that later passes check_model_files
and is skipped (never checksum-verified) by update_models on the next run."""
import importlib.util
import json
import socket
import sys
import tempfile
import threading
import types
from pathlib import Path

ROBIN = Path("./src/robin")


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# stub package hierarchy so model_updater's lazy "from robin.utils.model_checker import ..." works
sys.modules["robin"] = types.ModuleType("robin")
sys.modules["robin.utils"] = types.ModuleType("robin.utils")
checker = load("robin.utils.model_checker", ROBIN / "utils" / "model_checker.py")
updater = load("robin.utils.model_updater", ROBIN / "utils" / "model_updater.py")

# Tiny HTTP server that promises 1 MiB but drops the connection after 100 KiB
srv = socket.socket()
srv.bind(("127.0.0.1", 0))
srv.listen(8)
port = srv.getsockname()[1]


def serve():
    while True:
        c, _ = srv.accept()
        c.recv(65536)
        c.sendall(b"HTTP/1.1 200 OK\r\nContent-Length: 1048576\r\nConnection: close\r\n\r\n")
        c.sendall(b"x" * 100_000)
        import struct, time
        time.sleep(0.2)
        c.setsockopt(socket.SOL_SOCKET, socket.SO_LINGER, struct.pack("ii", 1, 0))
        c.close()


threading.Thread(target=serve, daemon=True).start()

work = Path(tempfile.mkdtemp(dir=sys.argv[1]))
models_dir = work / "models"
manifest = work / "assets.json"
manifest.write_text(json.dumps({"assets": {
    k: {"url": f"http://127.0.0.1:{port}/{fn}", "sha256": "0" * 64}
    for k, fn in checker.get_required_models()
}}))

ok, msgs = updater.update_models(models_dir=models_dir, manifest_path=str(manifest))
print("run 1 ok:", ok)
for m in msgs:
    print("   ", m)
for _, fn in checker.get_required_models():
    p = models_dir / fn
    print("   on disk:", fn, p.exists(), p.stat().st_size if p.exists() else None)

# simulate a second start: what check_model_files / update_models conclude
present_check = [((models_dir / fn).exists() and (models_dir / fn).stat().st_size > 0) for _, fn in checker.get_required_models()]
print("check_model_files-style presence (size>0):", present_check)
ok2, msgs2 = updater.update_models(models_dir=models_dir, manifest_path=str(manifest))
print("run 2 ok:", ok2)
for m in msgs2:
    print("   ", m)
