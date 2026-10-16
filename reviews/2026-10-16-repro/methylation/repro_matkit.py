"""Drive matkit.process_bam_counts_improved + write_counts_to_parquet with a fake pysam.

Checks:
 1. Fail calls (0 < max prob < 0.73) are counted in valid_cov -> percent_modified = n_mod/(mod+canon+fail).
 2. combine_mods uses max(P_m, P_h) rather than P_m + P_h (modkit --combine-mods sums them).
 3. Strand placement: forward read -> '+' at C; reverse read -> '-' at the G position (C+1).
"""
import importlib.util
import os
import sys
import tempfile
import types

import pyarrow.parquet as pq

# ---------------- stubs ----------------
class FakeRead:
    def __init__(self, name, ref_start, length, is_reverse, mods, qs=20):
        self.query_name = name
        self.is_unmapped = False
        self.is_secondary = False
        self.is_supplementary = False
        self.is_reverse = is_reverse
        self.reference_id = 0
        self._ref_start = ref_start
        self._len = length
        self.modified_bases = mods
        self._qs = qs

    def has_tag(self, t):
        return t == "qs"

    def get_tag(self, t):
        return self._qs

    def get_aligned_pairs(self, matches_only=True):
        return [(i, self._ref_start + i) for i in range(self._len)]


READS = []


class FakeAlignmentFile:
    def __init__(self, path, mode):
        self.references = ("chr1",)

    def fetch(self, contig=None, until_eof=False):
        return iter(READS)

    def close(self):
        pass


pysam = types.ModuleType("pysam")
pysam.AlignmentFile = FakeAlignmentFile
pysam.FastaFile = object
sys.modules["pysam"] = pysam
ap = types.ModuleType("alive_progress")
ap.alive_bar = lambda *a, **k: None
sys.modules["alive_progress"] = ap
for name in ["robin", "robin.analysis", "robin.analysis.utilities"]:
    m = types.ModuleType(name)
    m.__path__ = []
    sys.modules[name] = m
rb = types.ModuleType("robin.analysis.utilities.ReadBam")
rb.ReadBam = object
sys.modules["robin.analysis.utilities.ReadBam"] = rb
mf = types.ModuleType("robin.analysis.utilities.mnp_flex")
mf.APIClient = object
sys.modules["robin.analysis.utilities.mnp_flex"] = mf
res = types.ModuleType("robin.resources")
sys.modules["robin.resources"] = res
sys.modules["robin"].resources = res

spec = importlib.util.spec_from_file_location(
    "matkit", "./src/robin/analysis/utilities/matkit.py")
mk = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mk)


def run(reads, label):
    READS[:] = reads
    counts, sites, _ = mk.process_bam_counts_improved("x.bam", threshold=0.73, combine_mods=True)
    out = os.path.join(tempfile.mkdtemp(), "o.parquet")
    mk.write_counts_to_parquet(counts, sites, out)
    t = pq.read_table(out).to_pandas()
    for c in ("chrom", "mod_code", "strand"):
        t[c] = t[c].map(lambda b: b.decode())
    print(f"--- {label}")
    for k, c in sorted(counts.items()):
        print(f"   counts{k}: mod={c[0]} canon={c[1]} fail={c[4]} total={c[7]}")
    print(t.to_string(index=False))


# CpG with C at ref 1000 (0-based), G at 1001.
# 1. 8 forward reads: 6 confidently methylated (m=250), 2 confidently canonical (m=5), 2 'fail' (m=128)
reads = []
for i in range(6):
    reads.append(FakeRead(f"m{i}", 1000, 2, False, {("C", 0, "m"): [(0, 250)]}))
for i in range(2):
    reads.append(FakeRead(f"c{i}", 1000, 2, False, {("C", 0, "m"): [(0, 5)]}))
for i in range(2):
    reads.append(FakeRead(f"f{i}", 1000, 2, False, {("C", 0, "m"): [(0, 128)]}))
run(reads, "1. fail calls in denominator (modkit: 6/8 = 75.0%)")

# 2. read with P(5mC)=0.45, P(5hmC)=0.45 -> modkit --combine-mods: P(mod)=0.90 -> modified
reads = [FakeRead("mh", 1000, 2, False, {("C", 0, "m"): [(0, 115)], ("C", 0, "h"): [(0, 115)]})]
run(reads, "2. split m/h probabilities (modkit combine-mods: modified)")

# 3. reverse-strand read: pysam reports the call at the stored-SEQ G position -> ref 1001, strand flag 1
reads = [FakeRead("fwd", 1000, 2, False, {("C", 0, "m"): [(0, 250)]}),
         FakeRead("rev", 1000, 2, True, {("C", 1, "m"): [(1, 250)]})]
run(reads, "3. strand placement")
