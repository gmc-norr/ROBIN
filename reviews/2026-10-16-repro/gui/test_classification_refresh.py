import sys, types, tempfile, pathlib, time, os
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent / "stubs2"))
sys.path.insert(0, "./src")
import logging; logging.basicConfig(level=logging.ERROR, format="%(levelname)s %(message)s")
theme = types.ModuleType("robin.gui.theme")
theme.register_theme_sync_callback = lambda cb, **k: (lambda: None)
theme.styled_table = lambda **k: (None, None)
sys.modules["robin.gui.theme"] = theme
import nicegui
from robin.gui.components import classification as c
d = pathlib.Path(tempfile.mkdtemp())
def write(rows):
    with open(d/"sturgeon_scores.csv","w") as f:
        f.write("timestamp,Glioblastoma,Medulloblastoma,number_probes\n")
        for r in rows: f.write(",".join(map(str,r))+"\n")
write([(1000, 0.40, 0.60, 5000)])
(d/"NanoDX_scores.csv").write_text("timestamp,GBM,MB,number_probes\n1000,0.3,0.7,4000\n")
class L: workflow_steps = ["sturgeon", "nanodx", "cnv"]   # pannanodx / random_forest not run
c.add_classification_section(d, L())
refresh = [t for t in nicegui.TIMERS if t.interval == 30.0][0].cb
# find the Sturgeon summary label via closure
charts = [x for x in (cc.cell_contents for cc in refresh.__closure__) if isinstance(x, dict) and "Sturgeon" in x][0]
lbl = charts["Sturgeon"]["summary"]["class"]
refresh(); print("after tick 1:", lbl.text)
time.sleep(0.05)
write([(1000, 0.40, 0.60, 5000), (2000, 0.97, 0.03, 9000)])  # new batch: Glioblastoma now 97%
os.utime(d/"sturgeon_scores.csv", None)
refresh(); print("after tick 2 (new data written):", lbl.text)
refresh(); print("after tick 3:", lbl.text)
