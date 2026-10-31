#!/usr/bin/env python3
"""
Temporary utilities for robin BAM to parquet conversion.

This module provides utilities for merging modkit files and creating parquet files.
"""

import warnings
from typing import List, Optional
import gc
import os
import logging
import pickle
from datetime import datetime
import pandas as pd
import polars as pl
from contextlib import contextmanager
import tempfile
import json

from robin.utils.atomic_files import atomic_open

# Suppress pkg_resources deprecation warnings from sorted_nearest
warnings.filterwarnings(
    "ignore", message="pkg_resources is deprecated", category=UserWarning
)
# Suppress matplotlib tight_layout warnings
warnings.filterwarnings(
    "ignore", message="The figure layout has changed to tight", category=UserWarning
)

try:
    import pyranges as pr
    from robin.analysis.utilities.mnp_flex import APIClient as MnpFlexClient #ToDo: Maintain to future integration.
except ImportError as e:
    logging.warning(f"Some dependencies not available: {e}")



# Simple cross-process file lock using POSIX flock when available (no-op on unsupported platforms)
try:
    import fcntl  # type: ignore
except Exception:
    fcntl = None  # type: ignore


@contextmanager
def _exclusive_file_lock(lock_path: str):
    """Provide an exclusive lock on a file path. Best-effort no-op when flock is unavailable."""
    lock_dir = os.path.dirname(lock_path) or "."
    os.makedirs(lock_dir, exist_ok=True)
    lock_file = open(lock_path, "a+")
    if fcntl is not None:
        try:
            fcntl.flock(lock_file.fileno(), fcntl.LOCK_EX)
        except Exception:
            # Continue without locking if flock fails
            pass
    try:
        yield
    finally:
        if fcntl is not None:
            try:
                fcntl.flock(lock_file.fileno(), fcntl.LOCK_UN)
            except Exception:
                pass
        lock_file.close()


def _merge_failure_diagnostics(pl_df=None, gr=None, inter=None) -> str:
    """Describe merge state for a failed file: versions, shapes and contig names only (no positions or values)."""
    parts = [
        f"module={__file__}",
        f"pandas={pd.__version__}",
        f"polars={pl.__version__}",
        f"pyranges={getattr(globals().get('pr'), '__version__', 'not imported')}",
    ]
    checks = [
        ("rows", lambda: pl_df.height),
        ("schema", lambda: dict(pl_df.schema)),
        ("contigs", lambda: sorted(map(str, pl_df["chrom"].unique().to_list()))[:40]),
        ("ranges", lambda: len(gr)),
        ("intersect", lambda: f"{type(inter).__name__} rows={len(inter)} columns={list(inter.columns)}"),
    ]
    for name, get in checks:
        try:
            parts.append(f"{name}={get()}")
        except Exception:
            pass
    return "; ".join(parts)


# Bump when the panel's meaning changes, so old caches in sample folders are ignored.
# v2: 0-based BED as listed (no Start-1), negative starts and duplicates dropped.
_PANEL_CACHE_VERSION = "v2"


def _read_cpg_panel(filter_bed_file: str) -> pd.DataFrame:
    """Read the CpG panel as 0-based, half-open BED: Chromosome, Start, End.

    ``parquet_filter.txt`` (header "chr start end IlmnID") is already 0-based:
    each CpG is [C, C+2), and its starts read "CG" in hg38 (METH-6). Rows with a
    negative start (the file's first row is ``chr1 -1 0``) and exact duplicates
    are dropped (METH-5).
    """
    if filter_bed_file.endswith(".txt"):
        bed_df = pd.read_csv(filter_bed_file, sep=r"\s+", header=0, dtype=str)
        rename = {}
        for c in bed_df.columns:
            if c.lower() == "chr":
                rename[c] = "Chromosome"
            elif c.lower() == "start":
                rename[c] = "Start"
            elif c.lower() == "end":
                rename[c] = "End"
        bed_df = bed_df.rename(columns=rename)
    else:
        bed_df = pd.read_csv(
            filter_bed_file,
            sep="\t",
            header=None,
            names=["Chromosome", "Start", "End", "cg_label"],
            compression="gzip" if filter_bed_file.endswith(".gz") else None,
            dtype={"Chromosome": str},
        )
    bed_df = bed_df[["Chromosome", "Start", "End"]].astype({"Start": int, "End": int})
    bed_df = bed_df[bed_df["Start"] >= 0]
    return bed_df.drop_duplicates().reset_index(drop=True)


def load_cpg_panel(filter_bed_file: str, sample_output_dir: str) -> "pr.PyRanges":
    """Return the panel as PyRanges, cached in the sample folder.

    The cache name carries a format version and the panel file's size and
    modification time, so a code fix or an edited panel is never hidden by an
    old cache (METH-12). It is written whole or not at all.
    """
    stat = os.stat(filter_bed_file)
    cache_path = os.path.join(
        sample_output_dir,
        f"{os.path.basename(filter_bed_file)}.{_PANEL_CACHE_VERSION}"
        f".{stat.st_size}.{int(stat.st_mtime)}.pgr_cache",
    )
    if os.path.exists(cache_path):
        with open(cache_path, "rb") as f:
            return pickle.load(f)
    filter_ranges = pr.PyRanges(_read_cpg_panel(filter_bed_file))
    with atomic_open(cache_path, "wb") as f:
        pickle.dump(filter_ranges, f)
    return filter_ranges


# The sidecar JSON's column_format for a cumulative parquet of per-CpG counts
# (strands combined, percent from the counts). Older files are re-summed on
# their next merge by the same steps as new data.
CPG_COUNTS_FORMAT = "cpg_counts_v2"


def _panel_starts(filter_ranges: "pr.PyRanges") -> pl.DataFrame:
    """Each panel interval's start: the C of its CpG (0-based)."""
    df = filter_ranges.df[["Chromosome", "Start"]].drop_duplicates()
    return pl.DataFrame(
        {
            "chrom": df["Chromosome"].astype(str).tolist(),
            "chromStart": df["Start"].astype("int64").tolist(),
        },
        schema={"chrom": pl.Utf8, "chromStart": pl.Int64},
    )


def _sum_counts(df: pl.DataFrame) -> pl.DataFrame:
    return df.group_by(["chrom", "chromStart"]).agg(pl.sum("n_mod"), pl.sum("n_canonical"))


def _cpg_counts(df: pl.DataFrame, panel_starts: pl.DataFrame) -> pl.DataFrame:
    """Per-CpG counts: chrom, chromStart (the C), n_mod, n_canonical.

    A minus-strand row sits on the CpG's G, so it moves one base left onto the
    C, as ``modkit --combine-strands`` does (METH-4). Only rows on a panel
    CpG's C are kept, each once (a semi-join; METH-5), and counts are summed
    per CpG. Rows already combined (strand ".") pass through unchanged, which
    is what lets an old-format cumulative file be migrated by this same step.
    """
    shift = pl.when(pl.col("strand").cast(pl.Utf8) == "-").then(1).otherwise(0)
    df = df.select(
        pl.col("chrom").cast(pl.Utf8),
        (pl.col("chromStart").cast(pl.Int64) - shift).alias("chromStart"),
        pl.col("n_mod").cast(pl.Int64),
        pl.col("n_canonical").cast(pl.Int64),
    )
    return _sum_counts(df.join(panel_starts, on=["chrom", "chromStart"], how="semi"))


def _cpg_table(counts: pl.DataFrame) -> pl.DataFrame:
    """The cumulative parquet: one row per CpG, strand ".", valid_cov =
    n_mod + n_canonical and percent recomputed from them (METH-1, METH-2)."""
    valid = pl.col("n_mod") + pl.col("n_canonical")
    return (
        counts.filter(valid > 0)
        .sort(["chrom", "chromStart"])
        .select(
            pl.col("chrom"),
            pl.col("chromStart"),
            pl.lit("m").alias("mod_code"),
            pl.lit(".").alias("strand"),
            valid.cast(pl.UInt32).alias("valid_cov"),
            (100.0 * pl.col("n_mod") / valid).cast(pl.Float32).alias("percent_modified"),
            pl.col("n_mod").cast(pl.UInt32),
            pl.col("n_canonical").cast(pl.UInt32),
        )
    )


def merge_modkit_files(
    new_files: List[str],
    existing_file: str,
    output_file: str,
    filter_bed_file: str,
    sample_id: str,
    output_dir: str,
    mnpflex_config: Optional[dict],
    num_bam_files_seen: int,
) -> None:
    """
    Merge modkit files with optimized column set and improved caching.

    This function uses only essential columns to reduce memory usage and processing time:
    - chrom, chromStart: Required for all classifiers
    - mod_code: Required for Sturgeon filtering
    - strand: Required for proper aggregation
    - valid_cov: Required for coverage calculation
    - percent_modified: Primary methylation data (required)
    - n_mod, n_canonical: Required for modification counts

    Args:
        new_files (List[str]): List of new modkit files to merge
        existing_file (str): Path to existing parquet file
        output_file (str): Path to output parquet file
        filter_bed_file (str): Path to BED file for filtering
        sample_id (str): Sample ID for organizing output files
        output_dir (str): Base output directory
        mnpflex_config (Optional[dict]): Configuration for MNP-FLEX integration
        num_bam_files_seen (int): Number of BAM files being processed
    """
    # Create sample-specific output directory
    sample_output_dir = os.path.join(output_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)

    # Define optimized schema with only essential columns
    essential_cols = [
        "chrom",
        "chromStart",
        "mod_code",
        "strand",
        "valid_cov",
        "percent_modified",
        "n_mod",
        "n_canonical",
    ]

    categorical_cols = ["chrom", "mod_code", "strand"]
    unsigned_int_cols = ["chromStart", "valid_cov", "n_mod", "n_canonical"]
    float_cols = ["percent_modified"]

    try:
        # Track cumulative BAM file count
        cumulative_bam_file_count = 0

        # Check if existing file has metadata about BAM file count
        if os.path.exists(existing_file):
            try:
                # Read existing metadata
                metadata_file = existing_file.replace(".parquet", "_metadata.json")
                if os.path.exists(metadata_file):
                    with open(metadata_file, "r") as f:
                        metadata = json.load(f)
                        cumulative_bam_file_count = metadata.get("bam_file_count", 0)
                        logging.info(
                            f"Found existing metadata with {cumulative_bam_file_count} BAM files"
                        )
            except Exception as e:
                logging.warning(f"Could not read existing metadata: {str(e)}")
                cumulative_bam_file_count = 0

        # Add the number of new BAM files being processed
        cumulative_bam_file_count += num_bam_files_seen

        logging.info(
            f"Total cumulative BAM files contributing to parquet: {cumulative_bam_file_count} (added {num_bam_files_seen} new files)"
        )

        
        filter_ranges = load_cpg_panel(filter_bed_file, sample_output_dir)
        panel_starts = _panel_starts(filter_ranges)

        # BEDMethyl 18-column names for modkit output (tab-separated)
        full_cols = [
            "chrom",
            "chromStart",
            "chromEnd",
            "mod_code",
            "score_bed",
            "strand",
            "thickStart",
            "thickEnd",
            "color",
            "valid_cov",
            "percent_modified",
            "n_mod",
            "n_canonical",
            "n_othermod",
            "n_delete",
            "n_fail",
            "n_diff",
            "n_nocall",
        ]

        # Process new files: read with Polars, filter to CPG set, then combine
        new_frames: List[pl.DataFrame] = []
        for bed in new_files:
            pl_df = gr = inter = None
            try:
                if bed.endswith(".parquet"):
                    # Matkit wrote 8-column parquet; read directly (no CSV parse, types already correct)
                    pl_df = pl.read_parquet(bed, columns=essential_cols)
                    # Matkit may write chrom/mod_code/strand as binary; decode to Utf8 (Polars .str.decode is hex/base64 only)
                    def _binary_to_utf8(s: pl.Series) -> pl.Series:
                        return s.map_elements(
                            lambda x: x.decode("utf-8", errors="replace") if x is not None else None,
                            return_dtype=pl.Utf8,
                        )
                    for c in categorical_cols:
                        if pl_df.schema.get(c) == pl.Binary:
                            pl_df = pl_df.with_columns(_binary_to_utf8(pl.col(c)).alias(c))
                else:
                    # Legacy BEDMethyl text: read CSV, select essential columns, cast
                    pl_df = pl.read_csv(
                        bed,
                        separator="\t",
                        has_header=False,
                        new_columns=full_cols,
                        infer_schema_length=0,
                    ).select(essential_cols)
                    pl_df = pl_df.with_columns(
                        [pl.col(c).cast(pl.UInt32, strict=False) for c in unsigned_int_cols]
                        + [pl.col(c).cast(pl.Float32, strict=False) for c in float_cols]
                    )

                # Skip empty files (e.g. parquet from a batch where no reads passed the QS filter)
                if pl_df.is_empty():
                    continue

                counts = _cpg_counts(pl_df, panel_starts)
                # No sites on the panel
                if counts.is_empty():
                    continue
                new_frames.append(counts)
            except Exception as e:
                logging.error(
                    f"Error processing file {bed}: {str(e)} "
                    f"[{_merge_failure_diagnostics(pl_df, gr, inter)}]",
                    exc_info=True,
                )
                continue

        if not new_frames:
            logging.warning("No valid data to merge after filtering")
            return

        new_counts = _sum_counts(pl.concat(new_frames))

        # Critical section: read/merge/write parquet guarded by an exclusive lock
        lock_path = f"{output_file}.lock"
        with _exclusive_file_lock(lock_path):
            # The same steps on the first write and on every merge (METH-2). An
            # existing file is re-read through _cpg_counts, which leaves a
            # per-CpG file as it is and migrates an older stranded one.
            if os.path.exists(existing_file):
                existing = pl.read_parquet(
                    existing_file,
                    columns=["chrom", "chromStart", "strand", "n_mod", "n_canonical"],
                )
                counts = _sum_counts(pl.concat([_cpg_counts(existing, panel_starts), new_counts]))
            else:
                counts = new_counts
            table = _cpg_table(counts)

            # Atomic write: write to temp then replace
            fd, tmp_out = tempfile.mkstemp(
                prefix="lj_parquet_",
                suffix=".parquet",
                dir=os.path.dirname(output_file) or ".",
            )
            os.close(fd)
            try:
                table.write_parquet(tmp_out)
                os.replace(tmp_out, output_file)
            finally:
                try:
                    if os.path.exists(tmp_out):
                        os.remove(tmp_out)
                except Exception:
                    pass

            # Save metadata with updated cumulative BAM file count (atomic)
            metadata = {
                "bam_file_count": cumulative_bam_file_count,
                "last_updated": datetime.now().isoformat(),
                "sample_id": sample_id,
                "files_added_in_this_update": num_bam_files_seen,
                "column_format": CPG_COUNTS_FORMAT,
            }
            metadata_file = output_file.replace(".parquet", "_metadata.json")
            fdm, tmp_meta = tempfile.mkstemp(
                prefix="lj_meta_",
                suffix=".json",
                dir=os.path.dirname(metadata_file) or ".",
            )
            os.close(fdm)
            try:
                with open(tmp_meta, "w") as f:
                    json.dump(metadata, f, indent=2)
                os.replace(tmp_meta, metadata_file)
            finally:
                try:
                    if os.path.exists(tmp_meta):
                        os.remove(tmp_meta)
                except Exception:
                    pass

            logging.info(
                f"Updated methylation parquet: {table.height} CpGs from "
                f"{cumulative_bam_file_count} cumulative BAM files (added {num_bam_files_seen})"
            )

    except Exception as e:
        logging.error(f"Error in merge_modkit_files: {str(e)}")
        raise
    finally:
        gc.collect()
