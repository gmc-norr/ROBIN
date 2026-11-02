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


def cnv_worker(conn, contigs: Dict[str, int], control_per_bin: int, pause_seconds: float = 0.0) -> None:
    """A stand-in for a long-lived Ray worker process (run it in its own process).

    It runs the CNV jobs sent over ``conn`` as ``(kind, bams, sample_id, work_dir)``
    tuples, ``kind`` being "single" or "multiple", and keeps its module state between
    them, as a Ray worker does. It answers each with the job's error_message (None when
    it worked). ``None`` stops it. A pause after each BAM's analysis, before anything is
    saved, makes two workers overlap.
    """
    import logging
    import time

    from robin.analysis import cnv_analysis

    control = flat_control(contigs, control_per_bin)
    cnv_analysis.get_cached_ref_cnv_dict = lambda path, logger: control
    if pause_seconds:
        analyse = cnv_analysis.run_cnv_analysis_direct

        def slow(*args, **kwargs):
            result = analyse(*args, **kwargs)
            time.sleep(pause_seconds)
            return result

        cnv_analysis.run_cnv_analysis_direct = slow
    log = logging.getLogger("cnv_worker")
    while True:
        message = conn.recv()
        if message is None:
            return
        kind, bams, sample_id, work_dir = message
        if kind == "single":
            result = cnv_analysis.process_single_bam(
                bams[0], bam_metadata(sample_id, bams[0]), work_dir, log, threads=1
            )
        else:
            result = cnv_analysis.process_multiple_bams(
                bams, [bam_metadata(sample_id, b) for b in bams], work_dir, log, threads=1
            )
        conn.send(result.get("error_message"))
