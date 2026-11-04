"""A BAM that can't be read to its end fails preprocessing by name; it never becomes a sample.

Before this fix, a read error part-way through a BAM was swallowed: preprocessing reported
success with zero reads, and the sample id fell back to ``unknown_sample_<file>`` (pysam's
close error after a damaged block also broke the header read). A fake sample folder with an
all-zero master.csv appeared, and the later steps ran on it.

Synthetic BAMs only: the fixture header (read group LB "S"), enough reads for many BGZF
blocks, then one byte flipped in a middle block, or the file cut short.
"""

import os
import random
import struct
from itertools import count

import pytest

pysam = pytest.importorskip("pysam")

from target_fixtures import HEADER, SAMPLE  # noqa: E402
from robin import workflow_ray as wr  # noqa: E402
from robin.analysis import bam_preprocessor  # noqa: E402

_ids = count(1)
READS = 3000


def _write(path, n=READS):
    rng = random.Random(7)
    with pysam.AlignmentFile(path, "wb", header=HEADER) as out:
        for i, pos in enumerate(sorted(rng.randrange(0, 2700) for _ in range(n))):
            s = pysam.AlignedSegment()
            s.query_name = f"read{i:05d}"
            s.reference_id, s.reference_start, s.mapping_quality = 0, pos, 60
            s.cigarstring = "200M"
            s.query_sequence = "".join(rng.choice("ACGT") for _ in range(200))
            s.query_qualities = pysam.qualitystring_to_array("".join(rng.choice("5?I") for _ in range(200)))
            s.set_tag("RG", "rg1")
            out.write(s)
    return path


def _blocks(data):
    """(offset, size) of each BGZF block."""
    out, off = [], 0
    while off < len(data):
        size = struct.unpack_from("<H", data, off + 16)[0] + 1
        out.append((off, size))
        off += size
    return out


def _damaged(path):
    """One byte flipped inside a middle block's compressed data; the end block is kept."""
    data = bytearray(open(_write(path), "rb").read())
    blocks = _blocks(data)
    assert len(blocks) > 4, "the test BAM needs several blocks"
    off, size = blocks[len(blocks) // 2]
    data[off + size // 2] ^= 0xFF
    with open(path, "wb") as f:
        f.write(data)
    return path


def _cut(path):
    """The file cut at 60%: no end block."""
    data = open(_write(path), "rb").read()
    with open(path, "wb") as f:
        f.write(data[: int(len(data) * 0.6)])
    return path


def _run(bam, work):
    ctx = wr.WorkflowContext(filepath=str(bam), metadata={"work_dir": str(work), "target_panel": "rCNS2"})
    q = wr.job_queue_of("preprocessing")
    job = wr.Job(next(_ids), "preprocessing", q, [f"{q}:preprocessing"], 0, ctx)
    bam_preprocessor.bam_preprocessing_handler(job)
    return job


def _folders(work):
    return sorted(d for d in os.listdir(work) if os.path.isdir(os.path.join(work, d)))


@pytest.mark.parametrize("make", [_damaged, _cut], ids=["damaged_in_the_middle", "cut_off"])
def test_an_unreadable_bam_fails_preprocessing_by_name_and_makes_no_sample(tmp_path, make):
    bam = make(str(tmp_path / "FAB123_pass_0.bam"))
    work = tmp_path / "work"
    work.mkdir()

    job = _run(bam, work)

    errors = [e for e in job.context.errors if e["job_type"] == "preprocessing"]
    assert errors, "the job must fail, not report success with zero reads"
    assert "FAB123_pass_0.bam" in errors[0]["error"]
    assert "can't be read to its end" in errors[0]["error"]
    assert _folders(work) == [], "no sample folder (fake or real) for a BAM that wasn't counted"
    assert "preprocessing" not in (job.context.results or {}), "no success result"


def test_control_an_intact_bam_is_preprocessed_as_before(tmp_path):
    bam = _write(str(tmp_path / "FAB123_pass_0.bam"))
    work = tmp_path / "work"
    work.mkdir()

    job = _run(bam, work)

    assert not job.context.errors
    assert _folders(work) == [SAMPLE]
    result = job.context.results["preprocessing"]
    assert result["status"] == "success"
    assert result["mapped_reads"] == READS
