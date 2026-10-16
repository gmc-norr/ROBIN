import sys, types, tempfile, pathlib, time, os
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent / "stubs2"))
sys.path.insert(0, "./src")
import logging; logging.basicConfig(level=logging.CRITICAL)
import nicegui
theme = types.ModuleType("robin.gui.theme")
theme.register_theme_sync_callback = lambda cb, **k: (lambda: None)
theme.styled_table = lambda **k: (None, nicegui.El())
sys.modules["robin.gui.theme"] = theme
try:
    import matplotlib; matplotlib.use("Agg")
except ImportError:
    # placeholder-plot path needs matplotlib; provide a minimal fake
    m = types.ModuleType("matplotlib"); p = types.ModuleType("matplotlib.pyplot")
    class Ax:
        def __getattr__(s, n): return lambda *a, **k: None
        transAxes = None
    class Fig:
        def get_axes(s): return []
        patch = Ax(); texts = []
    p.subplots = lambda *a, **k: (Fig(), Ax()); p.close = lambda *a: None
    m.pyplot = p; m.patches = types.ModuleType("matplotlib.patches")
    sys.modules.update({"matplotlib": m, "matplotlib.pyplot": p, "matplotlib.patches": m.patches})
from robin.gui.components import mgmt
d = pathlib.Path(tempfile.mkdtemp())
(d/"1_mgmt.csv").write_text("average,pred,status\n3.1,12.0,Unmethylated\n")
class Launcher: _mgmt_state = {}
L = Launcher()
def build():
    nicegui.TIMERS.clear()
    mgmt.add_mgmt_section(L, d)
    t = [t for t in nicegui.TIMERS if t.interval == 30.0][0]
    status = [c.cell_contents for c in t.cb.__closure__ if isinstance(c.cell_contents, nicegui.El) and "MGMT status" in c.cell_contents.text]
    return t.cb, status
# _refresh_mgmt closes over _refresh_mgmt_sync which closes over mgmt_status; walk closures
def find_status(fn):
    for c in fn.__closure__ or ():
        v = c.cell_contents
        if isinstance(v, nicegui.El) and v.text.startswith("MGMT status"): return v
        if callable(v) and getattr(v, "__closure__", None):
            r = find_status(v)
            if r: return r
cbA, _ = build(); stA = find_status(cbA)
cbB, _ = build(); stB = find_status(cbB)
cbA(); L._mgmt_state.clear()  # B connects -> _refresh_sample_plots clears state
cbB()
print("initial  A:", stA.text, "| B:", stB.text)
time.sleep(0.05)
(d/"2_mgmt.csv").write_text("average,pred,status\n41.7,88.0,Methylated\n")   # new batch result
cbA(); cbB()
print("after new batch, A ticks first -> A:", stA.text, "| B:", stB.text)
cbA(); cbB()
print("next tick -> A:", stA.text, "| B:", stB.text)
