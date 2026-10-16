"""DBSCAN on (primary_mid, supp_mid) as in _cluster_breakpoint_pairs_dbscan:
eps=7500 (manhattan), min_samples=3, then require >=3 distinct reads.
Same read model as sim_bp_keys.py. Naive DBSCAN (no sklearn available)."""
import numpy as np
rng = np.random.default_rng(1)
J1, J2 = 10_000_000, 50_000_000

def dbscan(X, eps, min_samples):
    n = len(X); labels = -np.ones(n, int); c = 0
    D = np.abs(X[:, None, :] - X[None, :, :]).sum(-1)
    nb = [np.where(D[i] <= eps)[0] for i in range(n)]
    core = [len(x) >= min_samples for x in nb]
    for i in range(n):
        if labels[i] != -1 or not core[i]:
            continue
        labels[i] = c; stack = list(nb[i])
        while stack:
            j = stack.pop()
            if labels[j] == -1:
                labels[j] = c
                if core[j]:
                    stack.extend(nb[j])
        c += 1
    return labels

for n_reads in (3, 5, 10, 20):
    trials = 200; found = 0; frac = []
    for t in range(trials):
        feats = []
        for r in range(n_reads):
            L = int(rng.lognormal(np.log(10_000), 0.5))
            a = int(L * rng.uniform(0.2, 0.8)); b = L - a
            left = (J1 - a, J1); right = (J2, J2 + b)
            if a >= b:
                prim, supp = left, right
                feats.append(((prim[0] + prim[1]) / 2, (supp[0] + 1 + supp[0] + 1 + L) / 2))
            else:
                prim, supp = right, left
                feats.append(((supp[0] + 1 + supp[0] + 1 + L) / 2, (prim[0] + prim[1]) / 2))  # grouped by chrom pair later; keep orientation chr1,chr5
        X = np.array(feats)
        lab = dbscan(X, 7500, 3)
        sizes = [int((lab == k).sum()) for k in set(lab) if k != -1]
        if sizes and max(sizes) >= 3:
            found += 1
        frac.append(max(sizes) / n_reads if sizes else 0)
    print(f"reads={n_reads:3d}: SV reported in {found}/{trials} trials; mean fraction of reads in best cluster={np.mean(frac):.2f}")
