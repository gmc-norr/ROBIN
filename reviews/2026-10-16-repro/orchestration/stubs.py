"""Stub heavy deps (ray, watchdog, tqdm, GUI) so robin.workflow_ray can be imported for logic tests."""
import sys, types

SRC = "./src"
if SRC not in sys.path:
    sys.path.insert(0, SRC)

# ---- ray stub ----
ray = types.ModuleType("ray")
class _Remote:
    def __init__(self, obj): self._obj = obj
def _remote(obj=None, **kw):
    if obj is None:
        return lambda o: o
    return obj  # return the plain function/class so we can call it directly
ray.remote = _remote
def _get_actor(name): raise ValueError("no actor")
ray.get_actor = _get_actor
def _rt(): raise RuntimeError("no runtime")
ray.get_runtime_context = _rt
ray.get = lambda x, **k: x
ray.cancel = lambda *a, **k: None
ray.kill = lambda *a, **k: None
ray.wait = lambda refs, num_returns=1, timeout=None: ([], refs)
ray.cluster_resources = lambda: {}
ray.is_initialized = lambda: False
sys.modules["ray"] = ray

# ---- tqdm / watchdog stubs ----
tq = types.ModuleType("tqdm"); tq.tqdm = object; sys.modules["tqdm"] = tq
wd = types.ModuleType("watchdog"); wde = types.ModuleType("watchdog.events"); wdo = types.ModuleType("watchdog.observers")
class FileSystemEventHandler: pass
class Observer: pass
wde.FileSystemEventHandler = FileSystemEventHandler; wdo.Observer = Observer
sys.modules["watchdog"] = wd; sys.modules["watchdog.events"] = wde; sys.modules["watchdog.observers"] = wdo

# ---- GUI stub ----
import enum
gl = types.ModuleType("robin.gui_launcher")
class UpdateType(enum.Enum):
    LOG_MESSAGE = "log_message"; ERROR_UPDATE = "error_update"
GUI_SENT = []
gl.send_gui_update = lambda *a, **k: GUI_SENT.append((a, k))
gl.UpdateType = UpdateType
gl.launch_gui = lambda *a, **k: None
sys.modules["robin.gui_launcher"] = gl

# ---- avoid heavy robin.analysis import ----
an = types.ModuleType("robin.analysis")  # not a package -> handler imports fail -> None
sys.modules["robin.analysis"] = an

import robin.workflow_ray as wr  # noqa: E402
