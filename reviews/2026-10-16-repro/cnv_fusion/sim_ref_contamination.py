"""Simulate ROBIN direct-mode CNV pipeline to measure effect of the reference
dict being mutated in place by cnv_from_bam.iterate_bam_file (counts are
summed into the passed dict, per upstream lib.rs merge_and_reset/merge_vecs).

Re-implements cnv_from_bam.calculate_cnv: 1 kb base bins, bin_width =
genome_len / ceil(reads/100) rounded to 1000, chunk-sum / median * ploidy 2.
Then ROBIN: moving_average(n=3) of sample and reference, difference.
"""
import numpy as np

def uniform_filter1d(a, size, mode="nearest"):
    a = np.asarray(a, float)
    h = size // 2
    p = np.concatenate([np.repeat(a[:1], h), a, np.repeat(a[-1:], h)])
    return np.convolve(p, np.ones(size) / size, mode="valid")

rng = np.random.default_rng(0)
BIN = 1000
chroms = {f"chr{i}": 40_000_000 for i in range(1, 11)}
genome_len = sum(chroms.values())
copy_state = {c: 2 for c in chroms}
copy_state["chr1"] = 3  # gain
copy_state["chr2"] = 1  # loss


def draw_reads(n, states):
    w = np.array([chroms[c] * states[c] for c in chroms], float)
    per = rng.multinomial(n, w / w.sum())
    out = {}
    for (c, L), k in zip(chroms.items(), per):
        nb = L // BIN
        out[c] = np.bincount(rng.integers(0, nb, k), minlength=nb).astype(np.int64)
    return out


def add(acc, new):
    for c in new:
        acc[c] = acc.get(c, np.zeros_like(new[c])) + new[c]


def calc_cnv(counts, bin_width=None):
    nreads = sum(int(v.sum()) for v in counts.values())
    if bin_width is None:
        bw = genome_len / np.ceil(nreads / 100)
        bin_width = int((bw + BIN // 2) // BIN * BIN)
    k = max(1, bin_width // BIN)
    sums = {}
    for c, v in counts.items():
        n = int(np.ceil(len(v) / k))
        pad = np.zeros(n * k, np.int64)
        pad[: len(v)] = v
        sums[c] = pad.reshape(n, k).sum(1).astype(float)
    med = np.median(np.concatenate(list(sums.values())))
    med = 1.0 if med == 0 else round(med)
    return {c: s / med * 2 for c, s in sums.items()}, bin_width


def diff(sample_cnv, ref_cnv):
    return {
        c: uniform_filter1d(sample_cnv[c], 3, mode="nearest")
        - uniform_filter1d(ref_cnv[c], 3, mode="nearest")
        for c in sample_cnv
    }


control_states = {c: 2 for c in chroms}
control = draw_reads(400_000, control_states)  # fixed control counts

sample_acc = {}
ref_clean = {c: v.copy() for c, v in control.items()}  # never mutated (intended)
ref_cached = {c: v.copy() for c, v in control.items()}  # ROBIN direct-mode global cache

print("batch  sample_reads  chr1(gain) clean/contam   chr2(loss) clean/contam   chr3(neutral) clean/contam")
for b in range(1, 41):
    batch = draw_reads(25_000, copy_state)
    add(sample_acc, batch)
    s_cnv, bw = calc_cnv(sample_acc)
    # intended reference: control only at sample bin width
    r_clean, _ = calc_cnv(ref_clean, bw)
    # ROBIN direct mode: pass2 adds this batch into the cached reference in place
    add(ref_cached, batch)
    r_contam, _ = calc_cnv(ref_cached, bw)
    d_clean = diff(s_cnv, r_clean)
    d_contam = diff(s_cnv, r_contam)
    if b in (1, 2, 5, 10, 16, 20, 30, 40):
        nreads = sum(int(v.sum()) for v in sample_acc.values())
        print(
            f"{b:5d}  {nreads:12,d}  {np.mean(d_clean['chr1']):+.3f} / {np.mean(d_contam['chr1']):+.3f}"
            f"          {np.mean(d_clean['chr2']):+.3f} / {np.mean(d_contam['chr2']):+.3f}"
            f"          {np.mean(d_clean['chr3']):+.3f} / {np.mean(d_contam['chr3']):+.3f}"
        )
