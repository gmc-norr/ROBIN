"""Time merge_bedmethyl.map_methyl_calls_to_probes_chr (Sturgeon probe mapping) at panel scale."""
import importlib.util, sys, types, time
import numpy as np, pandas as pd
st = types.ModuleType("sturgeon"); stu = types.ModuleType("sturgeon.utils"); stu.read_probes_file = None
sys.modules["sturgeon"] = st; sys.modules["sturgeon.utils"] = stu
spec = importlib.util.spec_from_file_location("mb", "./src/robin/analysis/utilities/merge_bedmethyl.py")
mb = importlib.util.module_from_spec(spec); spec.loader.exec_module(mb)
rng = np.random.default_rng(0)
n_probes = 40_000   # roughly one chromosome's share of ~430k Sturgeon probes
probes = pd.DataFrame({"start": np.sort(rng.choice(250_000_000, n_probes, replace=False))})
pos = np.sort(np.concatenate([probes.start.values, probes.start.values + 1]))
calls = pd.DataFrame({"reference_pos": pos, "score": rng.random(len(pos))})
t = time.time()
mb.map_methyl_calls_to_probes_chr(probes.copy(), calls, 25, 0.3, 0.7)
dt = time.time() - t
print(f"{n_probes} probes: {dt:.2f}s -> extrapolated ~{dt*430_000/n_probes:.1f}s per Sturgeon run for full panel")
