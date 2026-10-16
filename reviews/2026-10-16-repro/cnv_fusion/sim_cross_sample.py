"""Cross-sample effect of the mutated cached reference (same Ray worker):
sample A (chr1 gain, chr2 loss) fully processed, then a NORMAL sample B."""
import importlib.util
spec = importlib.util.spec_from_file_location("s", __file__.replace("sim_cross_sample.py", "sim_ref_contamination.py"))
src = open(__file__.replace("sim_cross_sample.py", "sim_ref_contamination.py")).read().split("control_states =")[0]
ns = {}; exec(src, ns)
np = ns["np"]; draw = ns["draw_reads"]; add = ns["add"]; calc = ns["calc_cnv"]; diff = ns["diff"]; chroms = ns["chroms"]
normal = {c: 2 for c in chroms}
control = draw(400_000, normal)
ref_cached = {c: v.copy() for c, v in control.items()}
# sample A: 16 batches x 25k reads with chr1 gain / chr2 loss -> added to cached ref in pass 2
for b in range(16):
    add(ref_cached, draw(25_000, ns["copy_state"]))
# sample B (normal karyotype)
sB = {}
for b in range(1, 9):
    batch = draw(25_000, normal); add(sB, batch)
    s_cnv, bw = calc(sB)
    add(ref_cached, batch)
    r_cnv, _ = calc(ref_cached, bw)
    d = diff(s_cnv, r_cnv)
    if b in (1, 4, 8):
        print(f"normal sample B batch {b}: chr1 {np.mean(d['chr1']):+.3f}  chr2 {np.mean(d['chr2']):+.3f}  chr3 {np.mean(d['chr3']):+.3f}")
