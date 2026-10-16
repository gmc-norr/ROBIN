"""Simulate fusion_work master-BED breakpoint keying.

process_bam_single_pass builds rows: actual alignment (col4=master_bed_region)
and SA-derived rows with start=SA pos, end=SA pos + read.query_length
(fusion_work.py:1814-1821). _extract_master_bed_breakpoints_incremental keys each
(primary, supp) pair by the MIDPOINT of each interval binned to 500 bp
(fusion_work.py:4085-4099) and needs >=3 distinct reads per key.

Here: N reads all spanning ONE junction chr1:10,000,000 | chr5:50,000,000.
"""
import numpy as np
from collections import defaultdict

rng = np.random.default_rng(1)
J1, J2 = 10_000_000, 50_000_000
BIN = 500


def key(chrom, s, e, strand):
    mid = (s + e) // 2
    return (chrom, (mid // BIN) * BIN, strand)


for n_reads in (5, 10, 20, 50):
    support = defaultdict(set)
    for r in range(n_reads):
        L = int(rng.lognormal(np.log(10_000), 0.5))
        a = int(L * rng.uniform(0.2, 0.8))  # bases on chr1 (left of junction)
        b = L - a
        seg_left = ("chr1", J1 - a, J1)
        seg_right = ("chr5", J2, J2 + b)
        strand = "+"
        segs = [seg_left, seg_right]
        prim = max(segs, key=lambda s: s[2] - s[1])
        supp = seg_right if prim is seg_left else seg_left
        # rows emitted when iterating the primary record
        sa_supp = (supp[0], supp[1] + 1, supp[1] + 1 + L)  # SA pos is 1-based; end = pos + query_length
        k1 = key(*prim, strand)
        k2 = key(*sa_supp, strand)
        support[tuple(sorted((k1, k2)))].add(r)
    best = max(len(v) for v in support.values())
    n_ok = sum(1 for v in support.values() if len(v) >= 3)
    print(f"reads spanning the same junction={n_reads:3d}  distinct keys={len(support):3d}  "
          f"max reads/key={best}  keys passing min_read_support=3: {n_ok}")
