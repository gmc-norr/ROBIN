"""Coverage numbers that add up across a sample's BAMs.

Depth adds up across BAMs; "covered bases" (bases with at least one read) does
not. So each BAM gives exact integer totals (aligned bases per contig and per
target) and the stretches its reads cover, and the sample's covered bases come
from the union of those stretches. The rules are samtools': deletions and
reference skips are not aligned bases; unmapped, secondary, QC-fail and
duplicate records don't count; supplementary records do. The result equals
samtools run once on all of a sample's reads (step 6 of the code review).
"""

import os
import tempfile
from io import StringIO
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
import pysam

# samtools coverage's default --excl-flags: UNMAP, SECONDARY, QCFAIL, DUP
_SKIP_FLAGS = 0x4 | 0x100 | 0x200 | 0x400

CONTIG_COLUMNS = ["#rname", "startpos", "endpos", "numreads", "bases"]
TARGET_COLUMNS = ["chrom", "startpos", "endpos", "name", "bases"]
STRETCH_COLUMNS = ["chrom", "start", "end"]


def read_bed_targets(bed_path: str) -> List[Tuple[str, int, int, str]]:
    """The distinct (chrom, start, end, name) rows of a BED file, in file order."""
    rows, seen = [], set()
    with open(bed_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or not parts[0] or parts[0].startswith(("#", "track", "browser")):
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            name = parts[3] if len(parts) > 3 and parts[3].strip() else f"{chrom}:{start}-{end}"
            row = (chrom, start, end, name)
            if row not in seen:
                seen.add(row)
                rows.append(row)
    return rows


def _bedcov_j(bam_path: str, rows: List[Tuple[str, int, int, str]]) -> List[int]:
    """samtools bedcov -j (aligned bases, no deletions) for each row; 0 for a contig the BAM lacks."""
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        known = set(bam.references)
    present = [r for r in rows if r[0] in known]
    counts = {}
    if present:
        with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as f:
            for chrom, start, end, _ in present:
                f.write(f"{chrom}\t{start}\t{end}\n")
            bed = f.name
        try:
            lines = pysam.bedcov("-j", bed, bam_path).splitlines()
        finally:
            os.remove(bed)
        if len(lines) != len(present):
            raise RuntimeError(f"bedcov gave {len(lines)} rows for {len(present)} regions in {bam_path}")
        counts = {row: int(line.split("\t")[-1]) for row, line in zip(present, lines)}
    return [counts.get(r, 0) for r in rows]


def aligned_blocks(bam_path: str) -> Dict[str, np.ndarray]:
    """Per contig, the aligned blocks (0-based start, end) of every record samtools counts:
    deletions and reference skips split blocks; insertions and clips add none."""
    blocks: Dict[str, List[Tuple[int, int]]] = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            if read.flag & _SKIP_FLAGS or read.reference_id < 0:
                continue
            blocks.setdefault(read.reference_name, []).extend(read.get_blocks())
    return {chrom: np.asarray(pairs, dtype=np.int64).reshape(-1, 2) for chrom, pairs in blocks.items()}


def contig_totals(bam_path: str, blocks: Optional[Dict[str, np.ndarray]] = None) -> pd.DataFrame:
    """Per contig: samtools coverage's #rname, startpos, endpos (1-based) and numreads, plus
    ``bases``, the aligned bases (so mean depth = bases / contig length, exactly). ``bases``
    is the blocks' total length, which equals samtools bedcov -j over the whole contig (same
    records, no deletions; measured equal at 12,000x, so bedcov doesn't cap depth)."""
    if blocks is None:
        blocks = aligned_blocks(bam_path)
    cov = pd.read_csv(StringIO(pysam.coverage(bam_path)), sep="\t")
    cov = cov[["#rname", "startpos", "endpos", "numreads"]].copy()
    cov["bases"] = [
        int((blocks[c][:, 1] - blocks[c][:, 0]).sum()) if c in blocks else 0
        for c in cov["#rname"].astype(str)
    ]
    return cov[CONTIG_COLUMNS]


def target_bases(bam_path: str, bed_path: str) -> pd.DataFrame:
    """Aligned bases per distinct BED row (samtools bedcov -j)."""
    rows = read_bed_targets(bed_path)
    df = pd.DataFrame(rows, columns=["chrom", "startpos", "endpos", "name"])
    df["bases"] = _bedcov_j(bam_path, rows) if rows else []
    return df[TARGET_COLUMNS]


def _merge_intervals(starts: np.ndarray, ends: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    if len(starts) == 0:
        return starts, ends
    order = np.argsort(starts, kind="stable")
    s, e = starts[order], ends[order]
    run_end = np.maximum.accumulate(e)
    new = np.ones(len(s), dtype=bool)
    new[1:] = s[1:] > run_end[:-1]
    idx = np.flatnonzero(new)
    return s[idx], np.maximum.reduceat(e, idx)


def _merged_frame(per_contig: Dict[str, Tuple[np.ndarray, np.ndarray]]) -> pd.DataFrame:
    frames = []
    for chrom, (starts, ends) in per_contig.items():
        s, e = _merge_intervals(starts, ends)
        frames.append(pd.DataFrame({"chrom": chrom, "start": s, "end": e}))
    if not frames:
        return pd.DataFrame({"chrom": pd.Series(dtype=str), "start": pd.Series(dtype=np.int64),
                             "end": pd.Series(dtype=np.int64)})
    return pd.concat(frames, ignore_index=True)[STRETCH_COLUMNS]


def covered_stretches(bam_path: str, blocks: Optional[Dict[str, np.ndarray]] = None) -> pd.DataFrame:
    """The merged reference stretches that have at least one aligned base in this BAM."""
    if blocks is None:
        blocks = aligned_blocks(bam_path)
    return _merged_frame({chrom: (arr[:, 0], arr[:, 1]) for chrom, arr in blocks.items()})


def union_stretches(frames: Iterable[pd.DataFrame]) -> pd.DataFrame:
    """The union of several stretch tables."""
    frames = [f for f in frames if f is not None and len(f)]
    if not frames:
        return _merged_frame({})
    allf = pd.concat(frames, ignore_index=True)
    per_contig = {
        str(chrom): (g["start"].to_numpy(np.int64), g["end"].to_numpy(np.int64))
        for chrom, g in allf.groupby("chrom", sort=False)
    }
    return _merged_frame(per_contig)


def covered_bases(stretches: pd.DataFrame) -> Dict[str, int]:
    """Covered bases per contig: the total length of its stretches."""
    if stretches is None or not len(stretches):
        return {}
    lengths = (stretches["end"] - stretches["start"]).groupby(stretches["chrom"]).sum()
    return {str(k): int(v) for k, v in lengths.items()}


# --- Reading the sample's files (report and GUI) ---------------------------------


def target_length(df: pd.DataFrame) -> pd.Series:
    """Length of each 0-based, half-open BED row: end - start."""
    return (df["endpos"] - df["startpos"]).astype(float)


def with_target_depth(df: pd.DataFrame) -> pd.DataFrame:
    """A copy with ``length`` and ``coverage`` (mean depth) recomputed from ``bases``,
    so a file written with the old end - start + 1 length reads right too."""
    out = df.copy()
    out["length"] = target_length(out)
    if "bases" in out.columns:
        with np.errstate(divide="ignore", invalid="ignore"):
            out["coverage"] = out["bases"] / out["length"].replace(0, np.nan)
    return out


def target_depth(bed_df: pd.DataFrame):
    """Mean depth over all targets (aligned bases / target bases), or None."""
    if bed_df is None or "bases" not in bed_df.columns or not len(bed_df):
        return None
    length = target_length(bed_df).sum()
    return float(bed_df["bases"].sum() / length) if length > 0 else None


def contig_length(cov_df: pd.DataFrame) -> pd.Series:
    """Contig length from samtools coverage's 1-based startpos/endpos."""
    if "startpos" in cov_df.columns:
        return (cov_df["endpos"] - cov_df["startpos"] + 1).astype(float)
    return cov_df["endpos"].astype(float)


def global_depth(cov_df: pd.DataFrame):
    """Mean depth over the genome (sum of mean depth x length / genome length), or None.
    Exact for files written since step 6; samtools' printed mean depth before that."""
    if cov_df is None or not {"meandepth", "endpos"} <= set(cov_df.columns) or not len(cov_df):
        return None
    length = contig_length(cov_df)
    total = length.sum()
    return float((cov_df["meandepth"] * length).sum() / total) if total > 0 else None


def covered_bases_known(cov_df: pd.DataFrame) -> bool:
    """Covered bases are exact only in files written since step 6 (they carry ``bases``);
    before that they were per-BAM counts summed together."""
    return (
        cov_df is not None
        and "bases" in cov_df.columns
        and "covbases" in cov_df.columns
        and bool(cov_df["covbases"].notna().all())
    )
