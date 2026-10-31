"""Synthetic methylation data: a small reference with planted CpGs, BAMs whose
reads carry MM/ML tags, and modkit as the reference implementation.

No real reads. A read is described by which CpGs it calls and how sure it is;
``write_meth_bam`` turns that into SEQ, CIGAR, MM and ML the way a basecaller
would (MM counts C bases on the read as sequenced, so a reverse read counts on
the reverse complement).
"""

import array
import os
import random
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

import pysam

CHROM = "chr1"
MOD, CANON, FAIL = 250, 5, 128  # ML bytes: sure 5mC, sure C, unsure
COMP = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]


def write_reference(directory: str, cpgs: Iterable[int], length: int = 3000, seed: int = 7) -> str:
    """Write ``ref.fa`` (+ .fai): random sequence whose only CG dinucleotides start at ``cpgs`` (0-based)."""
    rng = random.Random(seed)
    seq: List[str] = []
    for _ in range(length):
        base = rng.choice("ACGT")
        if seq and seq[-1] == "C" and base == "G":
            base = "A"
        seq.append(base)
    for p in sorted(cpgs):
        seq[p], seq[p + 1] = "C", "G"
        if seq[p - 1] == "C" and p - 1 not in cpgs:
            seq[p - 1] = "A"
    ref = "".join(seq)
    found = [i for i in range(len(ref) - 1) if ref[i:i + 2] == "CG"]
    assert found == sorted(cpgs), found
    path = os.path.join(directory, "ref.fa")
    with open(path, "w") as f:
        f.write(f">{CHROM}\n")
        f.write("\n".join(ref[i:i + 60] for i in range(0, len(ref), 60)) + "\n")
    pysam.faidx(path)
    return path


def read_reference(path: str) -> str:
    with pysam.FastaFile(path) as fa:
        return fa.fetch(CHROM)


@dataclass
class ReadSpec:
    """One alignment record. ``m``/``h`` map a CpG's C position (0-based, + strand)
    to the ML byte of that call; a ``-`` read calls the CpG's C on the minus strand."""

    name: str
    strand: str  # "+" or "-"
    start: int
    end: int
    m: Dict[int, int] = field(default_factory=dict)
    h: Optional[Dict[int, int]] = None
    mode: str = "?"  # "?" explicit, "." implicit, "" implicit (no marker)
    qs: Optional[int] = 20
    supplementary: bool = False


def _mm_ml(spec: ReadSpec, seg: str) -> Tuple[str, List[int]]:
    orig = seg if spec.strand == "+" else revcomp(seg)
    length = len(seg)
    called: Dict[int, int] = {}
    for p in spec.m:
        idx = p - spec.start if spec.strand == "+" else length - 1 - (p + 1 - spec.start)
        assert orig[idx] == "C", (spec.name, p)
        called[idx] = p
    c_idx = [i for i, b in enumerate(orig) if b == "C"]
    skips, order, last = [], [], -1
    for rank, i in enumerate(c_idx):
        if i in called:
            skips.append(rank - last - 1)
            order.append(called[i])
            last = rank
    body = "".join("," + str(s) for s in skips)
    groups, ml = [], []
    if spec.h is not None:
        groups.append(f"C+h{spec.mode}{body};")
        ml += [spec.h[p] for p in order]
    groups.append(f"C+m{spec.mode}{body};")
    ml += [spec.m[p] for p in order]
    return "".join(groups), ml


def write_meth_bam(path: str, ref_path: str, reads: List[ReadSpec]) -> str:
    """Write a sorted, indexed BAM of ``reads`` against the reference at ``ref_path``."""
    ref = read_reference(ref_path)
    header = {"HD": {"VN": "1.6", "SO": "unsorted"}, "SQ": [{"SN": CHROM, "LN": len(ref)}]}
    unsorted = path + ".unsorted.bam"
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for spec in reads:
            seg = ref[spec.start:spec.end]
            mm, ml = _mm_ml(spec, seg)
            rec = pysam.AlignedSegment()
            rec.query_name = spec.name
            rec.query_sequence = seg
            rec.flag = (16 if spec.strand == "-" else 0) | (2048 if spec.supplementary else 0)
            rec.reference_id = 0
            rec.reference_start = spec.start
            rec.mapping_quality = 60
            rec.cigarstring = f"{len(seg)}M"
            rec.query_qualities = pysam.qualitystring_to_array("I" * len(seg))
            rec.set_tag("MM", mm, "Z")
            rec.set_tag("ML", array.array("B", ml))
            if spec.qs is not None:
                rec.set_tag("qs", spec.qs, "i")
            out.write(rec)
    pysam.sort("-o", path, unsorted)
    os.remove(unsorted)
    pysam.index(path)
    return path


def write_panel(path: str, rows: List[Tuple[int, int, str]]) -> str:
    """Write a panel in parquet_filter.txt's format (header, space separated)."""
    with open(path, "w") as f:
        f.write("chr start end IlmnID\n")
        for start, end, name in rows:
            f.write(f"{CHROM}\t{start}\t{end}\t{name}\n")
    return path


def modkit_available() -> bool:
    return shutil.which("modkit") is not None


def modkit_pileup(bam: str, ref: str, out_dir: str, combine_strands: bool) -> Dict[Tuple[int, str], Tuple[int, int, int, float]]:
    """modkit's per-site counts: (pos, strand) -> (n_mod, n_canonical, valid_cov, percent)."""
    out = os.path.join(out_dir, "modkit_cs.bed" if combine_strands else "modkit.bed")
    cmd = [
        "modkit", "pileup", bam, out, "--ref", ref, "--cpg", "--combine-mods",
        "--modified-bases", "C", "--filter-threshold", "0.73",
        "--log-filepath", out + ".log",
    ]
    if combine_strands:
        cmd.append("--combine-strands")
    subprocess.run(cmd, check=True, capture_output=True)
    sites = {}
    with open(out) as f:
        for line in f:
            col = line.split()
            sites[(int(col[1]), col[5])] = (int(col[11]), int(col[12]), int(col[9]), float(col[10]))
    return sites


def filter_bam(src: str, dst: str, drop_names: Iterable[str]) -> str:
    """Copy ``src`` to ``dst`` without the records of the named reads, sorted and indexed."""
    drop = set(drop_names)
    with pysam.AlignmentFile(src, "rb") as inp, pysam.AlignmentFile(dst, "wb", template=inp) as out:
        for rec in inp:
            if rec.query_name not in drop:
                out.write(rec)
    pysam.index(dst)
    return dst
