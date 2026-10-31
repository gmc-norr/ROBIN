"""Synthetic BAMs and controls for the CNV tests. No real reads."""

import os
from types import SimpleNamespace
from typing import Dict, List

import pysam

READ_LENGTH = 100


def write_bam(path: str, contigs: Dict[str, int], reads: Dict[str, List[int]]) -> str:
    """Write a sorted, indexed BAM of 100 M, MAPQ 60 reads at the given 0-based starts.

    ``contigs`` maps contig name to length (the header, in this order);
    ``reads`` maps contig name to read start positions. Returns ``path``.
    """
    header = {
        "HD": {"VN": "1.6", "SO": "unsorted"},
        "SQ": [{"SN": name, "LN": length} for name, length in contigs.items()],
    }
    names = list(contigs)
    unsorted = path + ".unsorted.bam"
    n = 0
    with pysam.AlignmentFile(unsorted, "wb", header=header) as out:
        for contig, starts in reads.items():
            for start in starts:
                read = pysam.AlignedSegment()
                read.query_name = f"read{n}"
                n += 1
                read.query_sequence = "A" * READ_LENGTH
                read.flag = 0
                read.reference_id = names.index(contig)
                read.reference_start = start
                read.mapping_quality = 60
                read.cigarstring = f"{READ_LENGTH}M"
                read.query_qualities = pysam.qualitystring_to_array("I" * READ_LENGTH)
                out.write(read)
    pysam.sort("-o", path, unsorted)
    os.remove(unsorted)
    pysam.index(path)
    return path


def every_kb(length: int, per_kb: int = 1) -> List[int]:
    """Read starts: ``per_kb`` reads at the start of each whole 1 kb bin."""
    return [start for start in range(0, length - 1000, 1000) for _ in range(per_kb)]


def flat_control(contigs: Dict[str, int], per_bin: int) -> Dict[str, List[int]]:
    """A control with ``per_bin`` reads in every 1 kb bin, sized as cnv_from_bam sizes them."""
    return {name: [per_bin] * (length // 1000 + 1) for name, length in contigs.items()}


def bam_metadata(sample_id: str, bam_path: str) -> SimpleNamespace:
    """The parts of a BamMetadata the CNV functions read."""
    return SimpleNamespace(
        extracted_data={"sample_id": sample_id},
        file_path=bam_path,
        file_size=os.path.getsize(bam_path),
    )
