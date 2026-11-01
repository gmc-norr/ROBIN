"""Synthetic BAMs for step 6 (target coverage, read counts, target.bam), and samtools
(through pysam) as the reference.

No real reads. Two contigs, a four-target panel (two targets overlap, one gets no
reads), and reads for every case the code review names: deletions, insertions and
soft clips, split reads with the supplementary on or off the panel, a secondary
alignment, unmapped reads, overlapping reads, and a read spanning two targets. A
read's records always go to the same BAM, as a basecaller writes them.
"""

import os
import random
from collections import Counter
from typing import Dict, List, Tuple

import pysam

CONTIGS = [("chr1", 3000), ("chr2", 2000)]
PANEL = [
    ("chr1", 100, 400, "GENE_A"),
    ("chr1", 300, 700, "GENE_B"),  # overlaps GENE_A
    ("chr1", 1500, 1800, "GENE_C"),
    ("chr2", 1700, 1900, "GENE_D"),  # no read reaches it
]
# Batch 8: the panel with its overlapping rows merged, written by hand (GENE_A + GENE_B).
MERGED_PANEL = [("chr1", 100, 700), ("chr1", 1500, 1800), ("chr2", 1700, 1900)]
# Batch 8: every way two rows can overlap, and the same panel merged by hand.
PANEL_OVERLAPS = PANEL + [
    ("chr1", 1550, 1650, "INSIDE_C"),  # a row inside another
    ("chr1", 1500, 1800, "GENE_C_ALIAS"),  # the same interval under another name
    ("chr1", 650, 1520, "BRIDGE_B_C"),  # a row bridging GENE_B and GENE_C
]
MERGED_PANEL_OVERLAPS = [("chr1", 100, 1800), ("chr2", 1700, 1900)]
SAMPLE = "S"
HEADER = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": n, "LN": ln} for n, ln in CONTIGS],
    "RG": [{
        "ID": "rg1", "DT": "2026-11-01T10:00:00Z", "LB": SAMPLE, "PL": "ONT",
        "DS": "runid=run1 basecall_model=dna_r10.4.1_e8.2_400bps_hac@v5.0.0",
    }],
}

Rec = Dict  # name, contig, start, cigar, flag (or qlen for unmapped)


def fixture_reads(seed: int = 11) -> List[List[Rec]]:
    """Each item is one read: its records (primary first)."""
    rng = random.Random(seed)
    reads: List[List[Rec]] = []

    def one(name, contig, start, cigar, flag=0):
        return {"name": name, "contig": contig, "start": start, "cigar": cigar, "flag": flag}

    for i in range(10):  # plain, either strand, anywhere on chr1
        ln = rng.randint(100, 250)
        reads.append([one(f"plain{i}", "chr1", rng.randint(0, 2700), f"{ln}M", 16 * (i % 2))])
    for i in range(4):  # deletions inside the panel
        reads.append([one(f"del{i}", "chr1", 120 + 140 * i, "60M7D80M")])
    for i in range(4):  # soft clip and insertion
        reads.append([one(f"ins{i}", "chr1", 1450 + 60 * i, "15S50M4I60M")])
    for i in range(3):  # primary in the panel, supplementary off it
        reads.append([one(f"spl_in{i}", "chr1", 150 + 40 * i, "80M60S"),
                      one(f"spl_in{i}", "chr2", 200 + 50 * i, "80H60M", 2048)])
    for i in range(3):  # primary off the panel, supplementary in it
        reads.append([one(f"spl_out{i}", "chr2", 300 + 70 * i, "70M50S"),
                      one(f"spl_out{i}", "chr1", 1550 + 30 * i, "70H50M", 2048)])
    for i in range(2):  # primary off the panel, secondary in it
        reads.append([one(f"sec{i}", "chr1", 2000 + 100 * i, "120M"),
                      one(f"sec{i}", "chr1", 320 + 30 * i, "120M", 256)])
    for i in range(4):  # unmapped
        reads.append([{"name": f"unm{i}", "flag": 4, "qlen": 150 + 10 * i}])
    for i in range(4):  # overlapping reads across GENE_A and GENE_B
        reads.append([one(f"span{i}", "chr1", 250 + 10 * i, "300M")])
    for i in range(3):  # chr2, before GENE_D
        reads.append([one(f"c2_{i}", "chr2", 100 + 400 * i, "200M", 16 * (i % 2))])
    return reads


def _segment(rec: Rec, rng: random.Random) -> pysam.AlignedSegment:
    seg = pysam.AlignedSegment()
    seg.query_name = rec["name"]
    seg.flag = rec["flag"]
    if rec["flag"] & 4:
        seg.reference_id, seg.reference_start, seg.mapping_quality = -1, -1, 0
        n = rec["qlen"]
    else:
        seg.reference_id = [c for c, _ in CONTIGS].index(rec["contig"])
        seg.reference_start = rec["start"]
        seg.mapping_quality = 60
        seg.cigarstring = rec["cigar"]
        n = seg.infer_query_length()
    seg.query_sequence = "".join(rng.choice("ACGT") for _ in range(n))
    seg.query_qualities = pysam.qualitystring_to_array("I" * n)
    return seg


def write_bam(path: str, reads: List[List[Rec]]) -> str:
    rng = random.Random(path)
    unsorted = path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=HEADER) as out:
        for read in reads:
            for rec in read:
                out.write(_segment(rec, rng))
    pysam.sort("-o", path, unsorted)
    os.remove(unsorted)
    pysam.index(path)
    return path


def write_bams(directory: str, k: int) -> Tuple[List[str], str]:
    """The fixture reads split into k BAMs (read i goes to BAM i % k), and the merged BAM."""
    reads = fixture_reads()
    bams = [write_bam(os.path.join(directory, f"{SAMPLE}_pass_{k}_{j}.bam"), reads[j::k]) for j in range(k)]
    merged = write_bam(os.path.join(directory, f"merged_{k}.bam"), reads)
    return bams, merged


def write_panel(path: str, rows=PANEL) -> str:
    with open(path, "w") as f:
        for chrom, start, end, name in rows:
            f.write(f"{chrom}\t{start}\t{end}\t{name}\n")
    return path


def samtools_merged_depth(merged_bam: str, merged_rows, directory: str) -> float:
    """samtools bedcov -j over a hand-merged panel, divided by its length: each base once."""
    rows = [(c, s, e, f"m{i}") for i, (c, s, e) in enumerate(merged_rows)]
    bed = write_panel(os.path.join(directory, "merged_panel.bed"), rows)
    bases = sum(int(line.split("\t")[-1]) for line in pysam.bedcov("-j", bed, merged_bam).splitlines())
    return bases / sum(e - s for _, s, e in merged_rows)


def record_key(line_or_seg) -> Tuple:
    """(name, flag, contig, 1-based pos, CIGAR) of a SAM line or an AlignedSegment."""
    if isinstance(line_or_seg, str):
        c = line_or_seg.split("\t")
        return (c[0], int(c[1]), c[2], int(c[3]), c[5])
    s = line_or_seg
    return (s.query_name, s.flag, s.reference_name or "*", s.reference_start + 1, s.cigarstring or "*")


def bam_records(path: str) -> Counter:
    with pysam.AlignmentFile(path, "rb") as f:
        return Counter(record_key(s) for s in f.fetch(until_eof=True))


def samtools_reference(merged: str, panel: str) -> Dict:
    """Locked rules 1-5, from samtools run once on the merged BAM."""
    targets = {}
    bedcov = pysam.bedcov("-j", panel, merged).splitlines()
    for (chrom, start, end, name), line in zip(PANEL, bedcov):
        bases = int(line.split("\t")[-1])
        reads = int(pysam.view("-c", "-F", "0x904", merged, f"{chrom}:{start + 1}-{end}").strip())
        targets[(chrom, start, end, name)] = (bases, bases / (end - start), reads)
    cov = {}
    for line in pysam.coverage(merged).splitlines()[1:]:
        c = line.split("\t")
        cov[c[0]] = (int(c[3]), int(c[4]))
    contigs, total_bases, total_len = {}, 0, 0
    for name, ln in CONTIGS:
        whole = os.path.join(os.path.dirname(merged), f"whole_{name}.bed")
        with open(whole, "w") as f:
            f.write(f"{name}\t0\t{ln}\t{name}\n")
        bases = int(pysam.bedcov("-j", whole, merged).split("\t")[-1])
        contigs[name] = (cov[name][0], cov[name][1], bases / ln)
        total_bases += bases
        total_len += ln
    stats = {}
    for line in pysam.stats(merged).splitlines():
        if line.startswith("SN\t"):
            key, value = line.split("\t")[1:3]
            stats[key.rstrip(":")] = value
    counts = {
        "mapped_reads": int(pysam.view("-c", "-F", "0x904", merged).strip()),
        "unmapped_reads": int(pysam.view("-c", "-f", "4", merged).strip()),
        "mapped_bases": int(stats["bases mapped"]),
        "yield_tracking": int(stats["total length"]),
    }
    target_bam = Counter(record_key(line) for line in pysam.view("-L", panel, merged).splitlines())
    return {
        "targets": targets, "contigs": contigs, "global": total_bases / total_len,
        "counts": counts, "target_bam": target_bam,
    }
