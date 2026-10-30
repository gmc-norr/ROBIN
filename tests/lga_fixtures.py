"""Tiny synthetic inputs for the lightweight gene analysis (LGA) tests. No patient data."""

from pathlib import Path

import pysam

CONTIG_LENGTH = 1000


def write_bam(path, reads):
    """Write an indexed BAM on chr1. ``reads`` is a list of ``(start, sequence)``."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chr1", "LN": CONTIG_LENGTH}],
    }
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for i, (start, sequence) in enumerate(sorted(reads, key=lambda r: r[0])):
            read = pysam.AlignedSegment(bam.header)
            read.query_name = f"read{i}"
            read.reference_id = 0
            read.reference_start = start
            read.query_sequence = sequence
            read.cigarstring = f"{len(sequence)}M"
            read.query_qualities = pysam.qualitystring_to_array("?" * len(sequence))
            read.mapping_quality = 60
            bam.write(read)
    pysam.index(str(path))
    return path


def write_clinvar(path, sites):
    """Write a ClinVar-style VCF. ``sites`` is a list of ``(pos, ref, alt)`` on contig 1, gene TEST."""
    lines = ["##fileformat=VCFv4.2", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
    for i, (pos, ref, alt) in enumerate(sites, 1):
        lines.append(f"1\t{pos}\t{i}\t{ref}\t{alt}\t.\t.\tCLNSIG=Pathogenic;GENEINFO=TEST:1;CLNDN=test")
    Path(path).write_text("\n".join(lines) + "\n")
    return path


def write_bed(path, rows):
    """Write a BED on chr1. ``rows`` is a list of ``(start, end, name)``."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"chr1\t{start}\t{end}\t{name}\n" for start, end, name in rows))
    return path
