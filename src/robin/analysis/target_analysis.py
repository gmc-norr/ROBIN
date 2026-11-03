#!/usr/bin/env python3
"""
Target Analysis Module for robin.

Requires Python 3.12+. Automated target analysis; integrates with robin's
workflow system and processes files for target-specific analysis.

Classes: TargetMetadata, TargetAnalysis.
"""

from __future__ import annotations

import sys
if sys.version_info < (3, 12):
    raise RuntimeError("robin target_analysis requires Python 3.12 or newer")

import os
import tempfile
import logging
import time
import gc
import json
import subprocess
import shutil
import glob
import fcntl
import uuid
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
from dataclasses import dataclass
from io import StringIO
import numpy as np
import pandas as pd
import pysam
from robin.logging_config import get_job_logger
from robin.bam_identity import bam_key
from robin.file_ledger import clear_bam_failures, record_bam_failure
from robin.analysis.snp_processing import build_snp_display_data
from robin.analysis.target_coverage import (
    UNION_FILE,
    aligned_blocks,
    contig_totals,
    covered_bases,
    covered_stretches,
    panel_union_bases,
    read_bed_targets,
    target_bases,
    union_stretches,
)
from robin.utils.atomic_files import atomic_open

# Optional import for Docker functionality
try:
    import docker
except ImportError:
    docker = None


def is_docker_available_for_snp_analysis() -> tuple[bool, str]:
    """
    Check if Docker is available for SNP analysis (Clair3 pipeline).
    Returns (True, "") if Docker is ready, (False, "error message") otherwise.
    """
    if docker is None:
        return False, "Docker Python package is not installed. Install it with: pip install docker"
    try:
        client = docker.from_env()
        client.ping()
        return True, ""
    except Exception as e:
        return False, f"Docker daemon is not available: {e}"


def json_serializable(obj):
    """Convert numpy/pandas types to JSON-serializable (3.12: type patterns)."""
    match obj:
        case int():
            return int(obj)
        case float():
            return float(obj)
        case _ if isinstance(obj, np.ndarray):
            return obj.tolist()
        case _ if isinstance(obj, pd.Series):
            return obj.tolist()
        case _ if isinstance(obj, pd.DataFrame):
            return obj.to_dict("records")
        case _:
            return obj


def _is_non_empty_file(path: Path) -> bool:
    """Return True when a path exists, is a file, and is non-empty."""
    return path.exists() and path.is_file() and path.stat().st_size > 0


def _resolve_clinvar_db_for_snpsift(logger: logging.Logger) -> Optional[str]:
    """
    Resolve ClinVar DB for SnpSift using bgzip + tabix index.

    Preferred and required format is `clinvar.vcf.gz` with a sibling `.tbi` index.
    If the index is missing, try to generate it in place.
    """
    try:
        # Use an alias to avoid shadowing module-level `resources`.
        from robin import resources as robin_resources

        resources_dir = Path(os.path.dirname(os.path.abspath(robin_resources.__file__)))
        clinvar_gz = resources_dir / "clinvar.vcf.gz"
        clinvar_tbi = resources_dir / "clinvar.vcf.gz.tbi"

        if not _is_non_empty_file(clinvar_gz):
            logger.warning(
                "ClinVar bgzipped VCF not found or empty at %s", clinvar_gz
            )
            return None

        if not _is_non_empty_file(clinvar_tbi):
            logger.info(
                "ClinVar tabix index missing at %s; attempting to generate it.",
                clinvar_tbi,
            )
            try:
                pysam.tabix_index(
                    str(clinvar_gz),
                    preset="vcf",
                    force=True,
                    keep_original=True,
                )
            except Exception as index_exc:
                logger.warning(
                    "Failed to create ClinVar tabix index for %s: %s",
                    clinvar_gz,
                    index_exc,
                )
                return None

        if not _is_non_empty_file(clinvar_tbi):
            logger.warning(
                "ClinVar tabix index still missing/empty after generation attempt: %s",
                clinvar_tbi,
            )
            return None

        logger.info(
            "Using ClinVar DB for SnpSift: %s (index: %s)",
            clinvar_gz,
            clinvar_tbi,
        )
        return str(clinvar_gz)
    except Exception as exc:
        logger.warning("Could not resolve ClinVar DB for SnpSift: %s", exc)
        return None


def _write_with_clinvar_header(src: str, dst: str, header_line: str) -> None:
    """Copy VCF ``src`` to ``dst`` with the ClinVar status line right after ``##fileformat=``,
    wherever that line is in the meta-information; first when there is none."""
    from robin.analysis.clinvar_status import CLINVAR_HEADER_KEY

    with open(src, "r") as fin, open(dst, "w") as fout:
        # The ## meta lines are few: read them all, then place the status line.
        meta = []
        first_other = None
        for line in fin:
            if line.startswith(CLINVAR_HEADER_KEY):
                continue  # never carry a stale status line over
            if not line.startswith("##"):
                first_other = line
                break
            meta.append(line)
        at = next(
            (i + 1 for i, line in enumerate(meta) if line.startswith("##fileformat=")), 0
        )
        meta.insert(at, header_line + "\n")
        fout.writelines(meta)
        if first_other is not None:
            fout.write(first_other)
        for line in fin:
            if line.startswith(CLINVAR_HEADER_KEY):
                continue
            fout.write(line)


def annotate_with_clinvar(
    snpeff_vcf: str,
    final_vcf: str,
    clinvar_db_path: Optional[str],
    logger: logging.Logger,
    *,
    verbose: bool = False,
    env: Optional[Dict[str, str]] = None,
    label: str = "SNP",
) -> str:
    """Make ``final_vcf`` from ``snpeff_vcf`` plus ClinVar, and record whether that worked.

    ``final_vcf`` always gets a ``##ROBIN_ClinVar=`` header line. If SnpSift ran
    and gave output it is ``annotated``; if the database is missing, SnpSift is
    not found, exits non-zero or gives nothing, the snpEff records are written
    unannotated and the line reads ``unavailable;reason=...``. The final file is
    only ever replaced by renaming a finished temp file from the same folder, so
    a reader never sees a half-written file or the previous one truncated.

    Returns ``"annotated"`` or ``"unavailable"``.
    """
    from robin.analysis.clinvar_status import (
        ANNOTATED,
        UNAVAILABLE,
        clean_reason,
        clinvar_header_line,
    )

    folder = os.path.dirname(os.path.abspath(final_vcf))
    reason = ""
    raw_path = None
    out_path = None
    try:
        fd, raw_path = tempfile.mkstemp(prefix=".clinvar_raw_", suffix=".vcf", dir=folder)
        os.close(fd)
        fd, out_path = tempfile.mkstemp(prefix=".clinvar_out_", suffix=".vcf", dir=folder)
        os.close(fd)

        if not clinvar_db_path:
            reason = "ClinVar database not found"
        else:
            cmd = ["SnpSift", "annotate"]
            if verbose:
                cmd.append("-v")
            cmd.extend([clinvar_db_path, snpeff_vcf])
            logger.info(f"Running SnpSift {label} command: {' '.join(cmd)}")
            try:
                with open(raw_path, "w") as fout:
                    result = subprocess.run(
                        cmd,
                        stdout=fout,
                        stderr=subprocess.PIPE,
                        text=True,
                        env=env,
                    )
            except FileNotFoundError:
                reason = "SnpSift not found"
                logger.error("SnpSift command not found in PATH")
            except Exception as exc:
                reason = f"SnpSift could not be run: {type(exc).__name__}: {exc}"
                logger.warning(f"SnpSift {label} annotation failed with exception: {exc}")
            else:
                if verbose and result.stderr:
                    logger.info(f"SnpSift {label} stderr output:\n{result.stderr}")
                if result.returncode != 0:
                    detail = (result.stderr or "").strip().splitlines()
                    reason = f"SnpSift exited with code {result.returncode}"
                    if detail:
                        reason += f": {detail[-1]}"
                    logger.warning(f"SnpSift ({label}) failed: {reason}")
                elif os.path.getsize(raw_path) == 0:
                    reason = "SnpSift produced no output"
                    logger.warning(f"SnpSift ({label}) produced no output")

        if reason:
            status = UNAVAILABLE
            source = snpeff_vcf
            header = clinvar_header_line(UNAVAILABLE, clean_reason(reason))
            logger.warning(
                f"ClinVar annotation unavailable for {label}: {reason}; "
                "writing the snpEff records marked as not ClinVar-annotated"
            )
        else:
            status = ANNOTATED
            source = raw_path
            header = clinvar_header_line(ANNOTATED)
            logger.info(f"SnpSift {label} completed successfully")

        _write_with_clinvar_header(source, out_path, header)
        os.replace(out_path, final_vcf)
        out_path = None
        return status
    finally:
        for leftover in (raw_path, out_path):
            if leftover and os.path.exists(leftover):
                try:
                    os.unlink(leftover)
                except OSError:
                    pass


class FileLock:
    """
    Simple file-based lock for coordinating access across processes/threads.
    Uses fcntl for POSIX systems.
    """
    
    def __init__(self, lock_file: str, timeout: float = 30.0):
        self.lock_file = lock_file
        self.timeout = timeout
        self.fd = None
    
    def __enter__(self):
        self.acquire()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.release()
    
    def acquire(self):
        """Acquire the lock with timeout"""
        os.makedirs(os.path.dirname(self.lock_file), exist_ok=True)
        self.fd = open(self.lock_file, 'w')
        
        start_time = time.time()
        while True:
            try:
                fcntl.flock(self.fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                return
            except IOError:
                if time.time() - start_time > self.timeout:
                    raise TimeoutError(f"Could not acquire lock {self.lock_file} within {self.timeout}s")
                time.sleep(0.1)
    
    def release(self):
        """Release the lock"""
        if self.fd:
            fcntl.flock(self.fd, fcntl.LOCK_UN)
            self.fd.close()
            self.fd = None


# Import robin resources for BED files
try:
    from robin import resources
except ImportError:
    resources = None


def _load_bed_regions(bedfile: str) -> List[Tuple[str, int, int]]:
    """Load BED regions into a list of (chrom, start, end) tuples."""
    regions: List[Tuple[str, int, int]] = []
    with open(bedfile, 'r') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                regions.append((chrom, start, end))
    return regions


def run_bedtools(bamfile, bedfile, tempbamfile, regions: Optional[List[Tuple[str, int, int]]] = None):
    """
    Extract target regions from BAM file, keeping all mappings (primary, secondary, supplementary)
    for reads that overlap the target regions.
    
    This function uses a two-step process:
    1. Extract read names from primary alignments that overlap target regions
    2. Extract ALL alignments (primary, secondary, supplementary) for those read names
    
    Parameters
    ----------
    bamfile : str
        Path to the input BAM file
    bedfile : str
        Path to the BED file defining regions
    tempbamfile : str
        Path where the output BAM file should be written
    regions : list of tuples, optional
        Pre-loaded BED regions to avoid repeated parsing
    """
    logger = logging.getLogger("robin.target")

    try:
        # Step 1: Extract read names from primary alignments overlapping target regions
        # Use bedtools to find primary alignments (filtering out supplementary and secondary)
        # Flag 2304 = 0x800 (supplementary) | 0x100 (secondary)
        # We use -F 2304 to exclude supplementary and secondary, keeping only primary alignments
        logger.debug(f"Step 1: Extracting read names from primary alignments overlapping {bedfile}")
        
        # Read BED file to get regions (unless preloaded)
        if regions is None:
            regions = _load_bed_regions(bedfile)
        
        if not regions:
            logger.warning(f"No valid regions found in BED file: {bedfile}")
            # Create empty BAM file
            with pysam.AlignmentFile(bamfile, "rb") as in_bam:
                header = in_bam.header.copy()
                with pysam.AlignmentFile(tempbamfile, "wb", header=header) as out_bam:
                    pass
            pysam.index(tempbamfile)
            return
        
        # Open input BAM and collect read names from primary alignments overlapping regions
        read_names = set()
        with pysam.AlignmentFile(bamfile, "rb") as in_bam:
            # Require BAM index for performance and correctness
            index_file = f"{bamfile}.bai"
            if not os.path.exists(index_file):
                raise FileNotFoundError(f"BAM index (.bai) not found for {bamfile}")
            
            # Use indexed access (faster)
            for chrom, start, end in regions:
                try:
                    # Fetch reads overlapping this region
                    for read in in_bam.fetch(chrom, start, end):
                        # Only consider primary alignments (not supplementary or secondary)
                        # Flag checks: not supplementary (0x800) and not secondary (0x100)
                        if not (read.flag & 0x800) and not (read.flag & 0x100):
                            read_names.add(read.query_name)
                except ValueError:
                    # Chromosome not found in BAM, skip
                    logger.debug(f"Chromosome {chrom} not found in BAM file, skipping")
                    continue
        
        logger.debug(f"Found {len(read_names)} unique read names overlapping target regions")
        
        if not read_names:
            logger.warning("No reads found overlapping target regions")
            # Create empty BAM file with same header
            with pysam.AlignmentFile(bamfile, "rb") as in_bam:
                header = in_bam.header.copy()
                with pysam.AlignmentFile(tempbamfile, "wb", header=header) as out_bam:
                    pass
            pysam.index(tempbamfile)
            return
        
        # Step 2: Extract ALL alignments (primary, secondary, supplementary) for those read names
        logger.debug(f"Step 2: Extracting all alignments for {len(read_names)} read names")
        
        reads_written = 0
        with pysam.AlignmentFile(bamfile, "rb") as in_bam:
            header = in_bam.header.copy()
            with pysam.AlignmentFile(tempbamfile, "wb", header=header) as out_bam:
                # Iterate through all reads in the BAM file
                for read in in_bam.fetch(until_eof=True):
                    if read.query_name in read_names:
                        out_bam.write(read)
                        reads_written += 1
                
                # Ensure all data is written to disk before closing
                out_bam.flush()
                logger.debug(f"Wrote {reads_written} alignments (including secondary/supplementary) to output BAM")
        
        # Verify the file was written successfully before indexing
        if not os.path.exists(tempbamfile):
            raise RuntimeError(f"Output BAM file was not created: {tempbamfile}")
        
        file_size = os.path.getsize(tempbamfile)
        if file_size == 0:
            logger.warning(f"Output BAM file is empty: {tempbamfile}")
        else:
            logger.debug(f"Output BAM file size: {file_size} bytes")
        
        # Index the output BAM (only if file has content)
        if file_size > 0:
            try:
                pysam.index(tempbamfile)
                logger.info(f"Successfully extracted target regions to {tempbamfile} ({reads_written} alignments)")
            except Exception as e:
                logger.error(f"Failed to index BAM file {tempbamfile}: {e}")
                # Try to verify if the BAM file is valid
                try:
                    with pysam.AlignmentFile(tempbamfile, "rb") as test_bam:
                        test_count = test_bam.count(until_eof=True)
                        logger.info(f"BAM file is readable, contains {test_count} reads")
                except Exception as verify_error:
                    logger.error(f"BAM file appears corrupted: {verify_error}")
                    raise
        else:
            logger.warning(f"Skipping indexing for empty BAM file: {tempbamfile}")
        
    except Exception as e:
        logger.error(f"Error in run_bedtools: {e}")
        import traceback
        logger.error(traceback.format_exc())


def get_covdfs(bamfile, bedfile=None):
    """
    Extract coverage information from a BAM file.

    This function calculates coverage statistics for both the entire genome
    and specific target regions defined in a BED file.

    Parameters
    ----------
    bamfile : str
        Path to the input BAM file.
    bedfile : str, optional
        Path to the target BED file. If not provided, will use unique_genes.bed.

    Returns
    -------
    tuple
        A tuple containing two pandas DataFrames:
        - Coverage statistics for the entire genome
        - Coverage statistics for target regions

    Notes
    -----
    Uses pysam for efficient BAM file processing.
    """
    logger = logging.getLogger("robin.target")

    try:
        # Get genome-wide coverage using pysam.coverage
        coverage_output = pysam.coverage(f"{bamfile}")

        newcovdf = pd.read_csv(StringIO(coverage_output), sep="\t")

        logger.debug(f"Raw pysam.coverage columns: {list(newcovdf.columns)}")
        logger.debug(f"Sample raw coverage data: {newcovdf.head(2).to_dict('records')}")

        newcovdf.drop(
            columns=["coverage", "meanbaseq", "meanmapq"],
            inplace=True,
        )

        logger.debug(f"After dropping columns: {list(newcovdf.columns)}")

        # Find target BED file - use provided bedfile or fallback to unique_genes.bed
        target_bed = bedfile
        if target_bed is None:
            # Fallback to unique_genes.bed for backward compatibility
            if resources is not None:
                try:
                    target_bed = os.path.join(
                        os.path.dirname(os.path.abspath(resources.__file__)),
                        "unique_genes.bed",
                    )
                    if not os.path.exists(target_bed):
                        target_bed = None
                except Exception:
                    pass

            # Fallback paths for unique_genes.bed
            if target_bed is None:
                possible_paths = [
                    "unique_genes.bed",
                    "data/unique_genes.bed",
                    "/usr/local/share/unique_genes.bed",
                ]
                for path in possible_paths:
                    if os.path.exists(path):
                        target_bed = path
                        break

        if target_bed is None:
            logger.warning("Target BED file not found, skipping bedcov analysis")
            bedcovdf = pd.DataFrame(
                columns=["chrom", "startpos", "endpos", "name", "bases"]
            )
        else:
            # Get target region coverage using pysam.bedcov
            bedcovdf = pd.read_csv(
                StringIO(
                    pysam.bedcov(
                        target_bed,
                        f"{bamfile}",
                    )
                ),
                names=["chrom", "startpos", "endpos", "name", "bases"],
                sep="\t",
            )

        logger.info(f"Successfully extracted coverage data from {bamfile}")
        logger.info(f"Genome coverage: {len(newcovdf)} regions")
        logger.info(f"Target coverage: {len(bedcovdf)} regions")

        return newcovdf, bedcovdf

    except Exception as e:
        logger.error(f"Error in get_covdfs: {str(e)}")
        return None, None


def get_read_counts_per_target(bamfile, bedfile):
    """
    Count reads overlapping each target region in a BED file.
    
    Parameters
    ----------
    bamfile : str
        Path to the input BAM file
    bedfile : str
        Path to the BED file defining target regions
    
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: chrom, startpos, endpos, name, reads
        Returns empty DataFrame if extraction fails
    """
    logger = logging.getLogger("robin.target")
    
    try:
        # Read BED file to get regions
        bed_regions = []
        with open(bedfile, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 4:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    name = parts[3] if parts[3].strip() else f"{chrom}:{start}-{end}"
                    bed_regions.append((chrom, start, end, name))
        
        if not bed_regions:
            logger.warning(f"No valid regions found in BED file: {bedfile}")
            return pd.DataFrame(columns=["chrom", "startpos", "endpos", "name", "reads"])
        
        # Count reads per region
        read_counts = []
        with pysam.AlignmentFile(bamfile, "rb") as bam:
            index_file = f"{bamfile}.bai"
            use_indexed = os.path.exists(index_file)
            
            for chrom, start, end, name in bed_regions:
                try:
                    # Count primary alignments only (not supplementary or secondary)
                    read_count = 0
                    if use_indexed:
                        # Use indexed access (faster)
                        for read in bam.fetch(chrom, start, end):
                            # Only count primary alignments
                            if not (read.flag & 0x800) and not (read.flag & 0x100):
                                read_count += 1
                    else:
                        # BAM is not indexed, count manually
                        for read in bam.fetch(chrom, start, end):
                            # Only count primary alignments
                            if not (read.flag & 0x800) and not (read.flag & 0x100):
                                read_count += 1
                    
                    read_counts.append({
                        'chrom': chrom,
                        'startpos': start,
                        'endpos': end,
                        'name': name,
                        'reads': read_count
                    })
                except ValueError:
                    logger.debug(f"Chromosome {chrom} not found in BAM file, skipping")
                    read_counts.append({
                        'chrom': chrom,
                        'startpos': start,
                        'endpos': end,
                        'name': name,
                        'reads': 0
                    })
                    continue
        
        df = pd.DataFrame(read_counts)
        logger.debug(f"Extracted read counts for {len(df)} target regions from {bamfile}")
        return df
        
    except Exception as e:
        logger.error(f"Error extracting read counts per target: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return pd.DataFrame(columns=["chrom", "startpos", "endpos", "name", "reads"])


def run_bedmerge(newcovdf, cov_df_main, bedcovdf, bedcov_df_main):
    """
    Merge coverage dataframes for incremental processing.

    Parameters
    ----------
    newcovdf : pd.DataFrame
        New genome coverage data
    cov_df_main : pd.DataFrame
        Existing genome coverage data
    bedcovdf : pd.DataFrame
        New target coverage data
    bedcov_df_main : pd.DataFrame
        Existing target coverage data

    Returns
    -------
    tuple
        Merged coverage dataframes
    """
    logger = logging.getLogger("robin.target")

    logger.debug("Starting coverage merge:")
    logger.debug(
        f"  New genome coverage: {len(newcovdf) if newcovdf is not None else 0} regions"
    )
    logger.debug(
        f"  Existing genome coverage: {len(cov_df_main) if cov_df_main is not None else 0} regions"
    )
    logger.debug(
        f"  New target coverage: {len(bedcovdf) if bedcovdf is not None else 0} targets"
    )
    logger.debug(
        f"  Existing target coverage: {len(bedcov_df_main) if bedcov_df_main is not None else 0} targets"
    )

    # Merge genome coverage data
    if cov_df_main is not None and not cov_df_main.empty:
        merged_df = pd.merge(
            newcovdf,
            cov_df_main,
            on=["#rname", "startpos", "endpos"],
            suffixes=("_df1", "_df2"),
        )
        merged_df["numreads"] = merged_df["numreads_df1"] + merged_df["numreads_df2"]
        merged_df["covbases"] = merged_df["covbases_df1"] + merged_df["covbases_df2"]
        merged_df["meandepth"] = merged_df["meandepth_df1"] + merged_df["meandepth_df2"]

        merged_df.drop(
            columns=[
                "numreads_df1",
                "numreads_df2",
                "meandepth_df1",
                "meandepth_df2",
                "covbases_df1",
                "covbases_df2",
            ],
            inplace=True,
        )

        logger.debug(
            f"Merged genome coverage: {len(cov_df_main)} + {len(newcovdf)} -> {len(merged_df)} regions"
        )
        logger.debug(
            f"Sample merged genome data: {merged_df.head(3).to_dict('records')}"
        )
    else:
        merged_df = newcovdf
        logger.debug(
            f"No existing genome coverage, using new data: {len(merged_df)} regions"
        )

    # Merge target coverage data
    if bedcov_df_main is not None and not bedcov_df_main.empty:
        merged_bed_df = pd.merge(
            bedcovdf,
            bedcov_df_main,
            on=["chrom", "startpos", "endpos", "name"],
            suffixes=("_df1", "_df2"),
        )
        merged_bed_df["bases"] = merged_bed_df["bases_df1"] + merged_bed_df["bases_df2"]
        merged_bed_df.drop(columns=["bases_df1", "bases_df2"], inplace=True)

        logger.debug(
            f"Merged target coverage: {len(bedcov_df_main)} + {len(bedcovdf)} -> {len(merged_bed_df)} targets"
        )
        logger.debug(
            f"Sample merged target data: {merged_bed_df.head(3).to_dict('records')}"
        )
    else:
        merged_bed_df = bedcovdf
        logger.debug(
            f"No existing target coverage, using new data: {len(merged_bed_df)} targets"
        )

    return merged_df, merged_bed_df


@dataclass(slots=True)
class TargetMetadata:
    """Container for target analysis metadata and results."""

    sample_id: str
    file_path: str
    analysis_timestamp: float
    target_data_path: Optional[str] = None
    target_plot_path: Optional[str] = None
    analysis_results: Optional[Dict] = None
    processing_steps: List[str] = None
    error_message: Optional[str] = None
    coverage_data: Optional[Dict] = None
    target_bam_path: Optional[str] = None
    coverage_over_time: Optional[np.ndarray] = None

    def __post_init__(self):
        if self.processing_steps is None:
            self.processing_steps = []
        if self.analysis_results is None:
            self.analysis_results = {}
        if self.coverage_data is None:
            self.coverage_data = {}


# A staged BAM is five files; the coverage file is written last (by rename).
_STAGED_KINDS = (
    ("coverage", "parquet"),
    ("bedcov", "parquet"),
    ("stretches", "parquet"),
    ("timestamp", "txt"),
    ("source_bam", "txt"),
)


def _regions(union: pd.DataFrame) -> frozenset:
    """The merged panel a union table covers: its (chrom, start, end) set."""
    return frozenset(
        (str(c), int(s), int(e)) for c, s, e in union[["chrom", "start", "end"]].itertuples(index=False)
    )


def _complete_staged_counters(staging_dir: str) -> List[str]:
    """Counters whose five staged files all exist, oldest first."""
    counters = []
    for path in sorted(glob.glob(os.path.join(staging_dir, "coverage_*.parquet"))):
        n = os.path.basename(path)[len("coverage_"):-len(".parquet")]
        if all(os.path.exists(os.path.join(staging_dir, f"{kind}_{n}.{ext}")) for kind, ext in _STAGED_KINDS):
            counters.append(n)
    return counters


class TargetAnalysis:
    """Target analysis worker"""

    def __init__(self, work_dir=None, config_path=None, threads=4, target_panel=None, 
                 batch_size=10, use_staging=True):
        logger = logging.getLogger("robin.target")

        self.work_dir = work_dir or os.getcwd()
        self.config_path = config_path
        self.threads = threads
        self.target_panel = target_panel
        self.batch_size = batch_size  # Trigger accumulation after N files
        self.use_staging = use_staging  # Enable/disable staging optimization

        # Initialize counter for incremental file naming
        self.file_counter = 1

        # Load configuration if provided
        self.config = self._load_config()

        # Note: State is now persisted to files instead of in-memory dictionaries

        # Find BED file for target extraction
        self.bedfile = self._find_target_bed(target_panel)

        # Configuration parameters
        self.callthreshold = self.config.get("call_threshold", 0.1)
        self.simtime = self.config.get("simtime", False)
        self.reference = self.config.get("reference", None)
        self.snp_calling = self.config.get("snp_calling", False)

        # Check environment variable for reference genome
        if not self.reference:
            self.reference = os.environ.get("robin_REFERENCE")

        # Log reference genome status
        if self.reference:
            logger.info(f"Reference genome configured: {self.reference}")
        else:
            logger.info(
                "No reference genome configured - SNP calling will not be available"
            )

        logger.info(f"Target Analysis initialized (staging={'enabled' if use_staging else 'disabled'}, batch_size={batch_size})")

    def _get_master_bed_path(self, sample_id: str) -> Optional[str]:
        """
        Get the path to the master BED file for a sample if it exists.
        Master BED includes target panel + CNV breakpoints + fusion breakpoints + master BED breakpoints.
        
        Args:
            sample_id: Sample ID
            
        Returns:
            Path to master BED file, or None if not found
        """
        logger = logging.getLogger("robin.target")
        try:
            import glob
            from robin.analysis.fusion_work import _load_analysis_counter
            
            sample_dir = os.path.join(self.work_dir, sample_id)
            bed_dir = os.path.join(sample_dir, "bed_files")
            
            if not os.path.exists(bed_dir):
                return None
            
            # Try to get analysis counter
            analysis_counter = _load_analysis_counter(sample_id, self.work_dir)
            master_bed_path = os.path.join(bed_dir, f"master_{analysis_counter:03d}.bed")
            
            if os.path.exists(master_bed_path):
                return master_bed_path
            
            # Try to find the latest master BED file if counter-based doesn't exist
            master_bed_files = glob.glob(os.path.join(bed_dir, "master_*.bed"))
            if master_bed_files:
                # Sort by modification time and return the latest
                latest = max(master_bed_files, key=os.path.getmtime)
                logger.debug(f"Using latest master BED file: {latest}")
                return latest
        except Exception as e:
            logger.debug(f"Error finding master BED file: {e}")
        
        return None

    def _find_target_bed(self, target_panel: str) -> str:
        """Find the target BED file from robin resources based on panel type"""
        logger = logging.getLogger("robin.target")
        
        # Determine the correct BED file name based on panel type
        bed_filename = None
        if target_panel == "rCNS2":
            bed_filename = "rCNS2_panel_name_uniq.bed"
        elif target_panel == "AML":
            bed_filename = "AML_panel_name_uniq.bed"
        else:
            # Check for custom panel
            bed_filename = f"{target_panel}_panel_name_uniq.bed"
            logger.info(f"Using custom panel: {target_panel}")
        
        if resources is not None:
            try:
                bed_path = os.path.join(
                    os.path.dirname(os.path.abspath(resources.__file__)),
                    bed_filename,
                )
                if os.path.exists(bed_path):
                    logger.info(f"Found {target_panel} BED file: {bed_path}")
                    return bed_path
            except Exception:
                pass

        # Fallback paths
        possible_paths = [
            bed_filename,
            f"data/{bed_filename}",
            f"/usr/local/share/{bed_filename}",
        ]

        for path in possible_paths:
            if os.path.exists(path):
                logger.info(f"Found {target_panel} BED file: {path}")
                return path

        # If not found, create a placeholder (this will cause an error later)
        logger.warning(f"Target BED file '{bed_filename}' not found for panel '{target_panel}', will use placeholder")
        return bed_filename

    def _load_config(self) -> Dict[str, Any]:
        """Load configuration from file or use defaults"""
        if self.config_path and os.path.exists(self.config_path):
            try:
                with open(self.config_path, "r") as f:
                    return json.load(f)
            except Exception as e:
                logging.warning(f"Error loading config from {self.config_path}: {e}")

        # Default configuration
        return {
            "analysis_type": "standard",
            "output_format": "json",
            "include_plots": True,
            "threshold": 0.5,
            "call_threshold": 0.1,
            "simtime": False,
            "snp_calling": False,
        }

    def _get_next_file_number(self) -> int:
        """Get the next file number for incremental processing"""
        return self.file_counter

    def _check_and_create_folder(self, base_dir: str, sample_id: str) -> str:
        """Create sample-specific directory and return path"""
        sample_dir = os.path.join(base_dir, sample_id)
        os.makedirs(sample_dir, exist_ok=True)
        return sample_dir
    
    def _get_staging_dir(self, sample_id: str) -> str:
        """Get staging directory for temporary per-file results"""
        staging_dir = os.path.join(self.work_dir, sample_id, "_staging")
        os.makedirs(staging_dir, exist_ok=True)
        return staging_dir
    
    def _get_lock_file(self, sample_id: str, lock_type: str = "counter") -> str:
        """Get lock file path for coordinating concurrent access"""
        lock_dir = os.path.join(self.work_dir, sample_id, "_locks")
        os.makedirs(lock_dir, exist_ok=True)
        return os.path.join(lock_dir, f"{lock_type}.lock")
    
    def _get_pending_count(self, sample_id: str) -> int:
        """Get count of files pending accumulation (thread-safe)"""
        staging_dir = self._get_staging_dir(sample_id)
        staging_files = glob.glob(os.path.join(staging_dir, "coverage_*.parquet"))
        return len(staging_files)
    
    def _atomic_counter_increment(self, sample_id: str) -> int:
        """
        Atomically increment and return the file counter for a sample.
        Uses file locking to prevent race conditions.
        
        Returns:
            The counter value to use for this file
        """
        lock_file = self._get_lock_file(sample_id, "counter")
        counter_file = os.path.join(self.work_dir, sample_id, "target_analysis_counter.txt")
        
        with FileLock(lock_file, timeout=30.0):
            # Read current counter
            if os.path.exists(counter_file):
                try:
                    with open(counter_file, "r") as f:
                        counter = int(f.read().strip())
                except (ValueError, IOError):
                    counter = 0
            else:
                counter = 0
            
            # Write incremented counter
            os.makedirs(os.path.dirname(counter_file), exist_ok=True)
            with open(counter_file, "w") as f:
                f.write(str(counter + 1))
            
            return counter
    
    def process_file_with_staging(
        self,
        file_path: str,
        metadata: Dict[str, Any],
        timestamp: Optional[float] = None,
    ) -> Tuple[TargetMetadata, bool]:
        """
        Fast per-file processing that saves results to staging area.
        Does NOT merge with accumulated data - much faster for large datasets.
        
        Args:
            file_path: Path to the input file
            metadata: File metadata from preprocessing
            timestamp: Optional timestamp for coverage tracking
        
        Returns:
            Tuple of (TargetMetadata, should_accumulate)
            - TargetMetadata: Results from this file
            - should_accumulate: True if batch accumulation should run now
        """
        logger = logging.getLogger("robin.target")
        
        logger.info(f"Processing file with staging: {file_path}")
        start_time = time.time()
        
        # Extract sample ID from metadata
        sample_id = metadata.get("sample_id", "unknown")
        logger.debug(f"Extracted sample_id: {sample_id}")
        
        target_result = TargetMetadata(
            sample_id=sample_id, file_path=file_path, analysis_timestamp=start_time
        )
        
        try:
            # Step 1: Validate input file
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"Input file not found: {file_path}")
            
            target_result.processing_steps.append("file_validation")
            
            # Step 2: Create sample-specific output directory
            sample_output_dir = self._check_and_create_folder(self.work_dir, sample_id)
            logger.info(f"Sample output directory: {sample_output_dir}")
            
            # Step 3: Get atomic counter (thread-safe)
            analysis_counter = self._atomic_counter_increment(sample_id)
            logger.info(f"Assigned file counter: {analysis_counter}")
            
            target_result.processing_steps.append("counter_assigned")
            
            # Step 4: This BAM's totals, which add up across BAMs (no accumulated data loaded)
            logger.info("Extracting coverage data...")
            blocks = aligned_blocks(file_path)  # one pass gives contig bases and stretches
            newcovdf = contig_totals(file_path, blocks)
            bedcovdf = target_bases(file_path, self.bedfile)
            union = panel_union_bases(file_path, self.bedfile)  # each panel base once
            stretches = covered_stretches(file_path, blocks)
            del blocks

            target_result.processing_steps.append("coverage_extracted")
            logger.info(
                f"Coverage data extracted: genome={newcovdf.shape}, targets={bedcovdf.shape}"
            )

            # Step 5: Save to staging. The coverage file is written last, by rename:
            # a claim takes only counters whose coverage file exists, so it never
            # meets a half-written set.
            staging_dir = self._get_staging_dir(sample_id)
            n = f"{analysis_counter:06d}"
            coverage_staging = os.path.join(staging_dir, f"coverage_{n}.parquet")
            bedcovdf.to_parquet(os.path.join(staging_dir, f"bedcov_{n}.parquet"), index=False)
            union.to_parquet(os.path.join(staging_dir, f"union_{n}.parquet"), index=False)
            stretches.to_parquet(os.path.join(staging_dir, f"stretches_{n}.parquet"), index=False)
            current_timestamp = timestamp * 1000 if (self.simtime and timestamp) else time.time() * 1000
            with open(os.path.join(staging_dir, f"timestamp_{n}.txt"), "w") as f:
                f.write(str(current_timestamp))
            with open(os.path.join(staging_dir, f"source_bam_{n}.txt"), "w") as f:
                f.write(file_path)
            # The BAM's key, read now while the BAM is here: a failure found later is then
            # recorded for this BAM even if it was renamed by then (batch 10e)
            with open(os.path.join(staging_dir, f"bam_key_{n}.txt"), "w") as f:
                f.write(bam_key(file_path))
            with atomic_open(coverage_staging, "wb") as f:
                newcovdf.to_parquet(f, index=False)

            target_result.processing_steps.append("saved_to_staging")
            logger.info(f"Saved to staging: {coverage_staging}")
            
            # Step 6: Store minimal coverage data in metadata (for logging)
            target_result.coverage_data = {
                "genome_coverage_shape": newcovdf.shape,
                "target_coverage_shape": bedcovdf.shape,
                "staging_file": coverage_staging,
            }
            
            # Step 7: Check if accumulation should run
            pending_count = self._get_pending_count(sample_id)
            should_accumulate = pending_count >= self.batch_size
            
            logger.info(
                f"File staged successfully. Pending files: {pending_count}/{self.batch_size}"
            )
            
            if should_accumulate:
                logger.info(
                    f"Accumulation threshold reached ({pending_count} >= {self.batch_size})"
                )
            
            # Step 8: Force garbage collection
            gc.collect()
            
            target_result.processing_steps.append("staging_complete")
            elapsed = time.time() - start_time
            logger.info(
                f"Staging complete for {sample_id} in {elapsed:.2f}s (vs ~{elapsed*10:.1f}s without staging)"
            )
            
            return target_result, should_accumulate
            
        except Exception as e:
            error_details = f"Error in staging for {sample_id}: {str(e)}"
            logger.error(error_details)
            target_result.error_message = error_details
            target_result.processing_steps.append("staging_failed")
            return target_result, False

    def process_file(
        self,
        file_path: str,
        metadata: Dict[str, Any],
        timestamp: Optional[float] = None,
    ) -> TargetMetadata:
        """
        Process a file for target analysis.

        Args:
            file_path: Path to the input file
            metadata: File metadata from preprocessing
            timestamp: Optional timestamp for coverage tracking

        Returns:
            TargetMetadata object with analysis results
        """
        logger = logging.getLogger("robin.target")

        logger.info(f"Processing file: {file_path} for target analysis")
        start_time = time.time()

        # Extract sample ID from metadata
        sample_id = metadata.get("sample_id", "unknown")
        logger.debug(f"Extracted sample_id: {sample_id}")

        target_result = TargetMetadata(
            sample_id=sample_id, file_path=file_path, analysis_timestamp=start_time
        )
        logger.debug(f"Created TargetMetadata object: {target_result}")

        logger.info(f"Starting target analysis for sample: {sample_id}")
        logger.debug(f"File: {os.path.basename(file_path)}")
        logger.debug(f"Work directory: {self.work_dir}")

        try:
            # Step 1: Validate input file
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"Input file not found: {file_path}")

            target_result.processing_steps.append("file_validation")

            # Step 2: Create sample-specific output directory
            sample_output_dir = self._check_and_create_folder(self.work_dir, sample_id)
            logger.info(f"Created sample output directory: {sample_output_dir}")

            target_result.processing_steps.append("directory_created")

            # Step 3: Load analysis counter from disk
            analysis_counter = self._load_analysis_counter(
                sample_id, self.work_dir, logger
            )
            logger.info(f"Loaded analysis counter for {sample_id}: {analysis_counter}")

            # Step 4: Load existing accumulated data from files
            existing_covdf = self._load_existing_coverage_data(
                sample_output_dir, "coverage_main.csv", logger
            )
            existing_bedcovdf = self._load_existing_coverage_data(
                sample_output_dir, "bed_coverage_main.csv", logger
            )
            existing_coverage_over_time = self._load_existing_coverage_over_time(
                sample_output_dir, logger
            )
            existing_target_bam_path = os.path.join(sample_output_dir, "target.bam")

            # Step 5: Create temporary directory for processing
            with tempfile.TemporaryDirectory(dir=sample_output_dir) as temp_dir:
                logger.info(f"Created temporary directory: {temp_dir}")

                # Step 6: Extract coverage data using get_covdfs
                logger.info("Extracting coverage data...")
                newcovdf, bedcovdf = get_covdfs(file_path, self.bedfile)

                if newcovdf is None or bedcovdf is None:
                    raise RuntimeError("Failed to extract coverage data from BAM file")

                # Store coverage data in metadata
                target_result.coverage_data = {
                    "genome_coverage_shape": newcovdf.shape,
                    "target_coverage_shape": bedcovdf.shape,
                }

                target_result.processing_steps.append("coverage_extracted")
                logger.info(
                    f"Coverage data extracted: genome={newcovdf.shape}, targets={bedcovdf.shape}"
                )

                # Step 7: Extract target regions using bedtools
                logger.info("Extracting target regions using bedtools...")
                tempbamfile = tempfile.NamedTemporaryFile(
                    dir=sample_output_dir, suffix=".bam"
                )

                # Prefer master BED file if available, otherwise use original target panel BED
                # Master BED includes target panel + CNV + fusion + master BED breakpoints
                sample_id = metadata.get("sample_id", "unknown")
                targets_bed = self._get_master_bed_path(sample_id) or self.bedfile
                
                # Run bedtools intersection
                run_bedtools(file_path, targets_bed, tempbamfile.name)

                # Check if target BAM has reads
                if (
                    pysam.AlignmentFile(tempbamfile.name, "rb").count(until_eof=True)
                    > 0
                ):
                    logger.info(f"Target regions found in {sample_id}")

                    # Handle target BAM file accumulation
                    if os.path.exists(existing_target_bam_path):
                        # Merge with existing target BAM
                        tempbamholder = tempfile.NamedTemporaryFile(
                            dir=sample_output_dir, suffix=".bam"
                        )
                        pysam.cat(
                            "-o",
                            tempbamholder.name,
                            existing_target_bam_path,
                            tempbamfile.name,
                        )
                        shutil.copy2(tempbamholder.name, existing_target_bam_path)
                        target_result.target_bam_path = existing_target_bam_path
                        logger.info(
                            f"Updated target BAM file: {existing_target_bam_path}"
                        )

                        # Clean up temporary files
                        try:
                            os.remove(f"{tempbamholder.name}.bai")
                        except FileNotFoundError:
                            pass
                    else:
                        # Create new target BAM file
                        shutil.copy2(tempbamfile.name, existing_target_bam_path)
                        target_result.target_bam_path = existing_target_bam_path
                        logger.info(
                            f"Created target BAM file: {existing_target_bam_path}"
                        )
                else:
                    logger.info(f"No target regions found in {sample_id}")

                # Clean up temporary BAM index
                try:
                    os.remove(f"{tempbamfile.name}.bai")
                except FileNotFoundError:
                    pass

                target_result.processing_steps.append("target_regions_extracted")

                # Step 8: Update coverage data
                if existing_covdf is not None:
                    updated_covdf, updated_bedcovdf = run_bedmerge(
                        newcovdf, existing_covdf, bedcovdf, existing_bedcovdf
                    )
                else:
                    updated_covdf, updated_bedcovdf = newcovdf, bedcovdf

                target_result.processing_steps.append("coverage_updated")

                # Step 9: Calculate coverage statistics
                bases = updated_covdf["covbases"].sum()
                genome = updated_covdf["endpos"].sum()
                coverage = bases / genome if genome > 0 else 0.0

                logger.info(
                    f"Coverage calculated: {coverage:.4f} ({bases} bases / {genome} genome)"
                )
                logger.debug(f"Coverage data shape: {updated_covdf.shape}")
                logger.debug(f"Coverage columns: {list(updated_covdf.columns)}")
                logger.debug("Sample coverage data:")
                logger.debug(
                    f"  First few rows: {updated_covdf.head(3).to_dict('records')}"
                )
                logger.debug(f"  covbases sum: {bases}")
                logger.debug(f"  endpos sum: {genome}")

                # Step 10: Track coverage over time
                if self.simtime and timestamp:
                    currenttime = timestamp * 1000
                else:
                    currenttime = time.time() * 1000

                if existing_coverage_over_time is not None:
                    updated_coverage_over_time = np.vstack(
                        [existing_coverage_over_time, [(currenttime, coverage)]]
                    )
                else:
                    updated_coverage_over_time = np.array([[currenttime, coverage]])

                target_result.coverage_over_time = updated_coverage_over_time
                target_result.processing_steps.append("coverage_tracked")

                # Step 11: Save updated coverage data to files
                np.save(
                    os.path.join(sample_output_dir, "coverage_time_chart.npy"),
                    updated_coverage_over_time,
                )

                updated_covdf.to_csv(
                    os.path.join(sample_output_dir, "coverage_main.csv"),
                    index=False,
                )

                # Calculate length and coverage for bed_coverage_main.csv (matching reference format)
                bed_coverage_main_df = updated_bedcovdf.copy()
                bed_coverage_main_df["length"] = (
                    bed_coverage_main_df["endpos"]
                    - bed_coverage_main_df["startpos"]
                    + 1
                )

                bed_coverage_main_df["coverage"] = (
                    bed_coverage_main_df["bases"] / bed_coverage_main_df["length"]
                )

                # Reorder columns to match reference format: chrom,startpos,endpos,name,length,coverage,bases
                bed_coverage_main_df = bed_coverage_main_df[
                    [
                        "chrom",
                        "startpos",
                        "endpos",
                        "name",
                        "length",
                        "coverage",
                        "bases",
                    ]
                ]

                bed_coverage_main_df.to_csv(
                    os.path.join(sample_output_dir, "bed_coverage_main.csv"),
                    index=False,
                )

                target_result.processing_steps.append("coverage_saved")

                # Step 12: Calculate target coverage
                target_coverage_df = updated_bedcovdf.copy()
                target_coverage_df["length"] = (
                    target_coverage_df["endpos"] - target_coverage_df["startpos"] + 1
                )

                target_coverage_df["coverage"] = (
                    target_coverage_df["bases"] / target_coverage_df["length"]
                )

                # Reorder columns to match reference format: chrom,startpos,endpos,name,length,coverage,bases
                target_coverage_df = target_coverage_df[
                    [
                        "chrom",
                        "startpos",
                        "endpos",
                        "name",
                        "length",
                        "coverage",
                        "bases",
                    ]
                ]

                target_coverage_df.to_csv(
                    os.path.join(sample_output_dir, "target_coverage.csv"),
                    index=False,
                )

                target_result.processing_steps.append("target_coverage_calculated")

                # Step 13b: Create target coverage with timestamp and read counts (optimized)
                logger.info("Calculating read counts per target...")
                new_read_counts_df = get_read_counts_per_target(file_path, self.bedfile)
                
                if not new_read_counts_df.empty:
                    # Use optimized approach: store latest cumulative reads in separate Parquet file for fast access
                    time_coverage_file = os.path.join(sample_output_dir, "target_coverage_time.csv")
                    latest_reads_cache = os.path.join(sample_output_dir, "_target_coverage_latest_reads.parquet")
                    
                    # Load previous cumulative reads from cache (much faster than reading entire CSV)
                    previous_cumulative_reads = None
                    if os.path.exists(latest_reads_cache):
                        try:
                            previous_cumulative_reads = pd.read_parquet(latest_reads_cache)
                            previous_cumulative_reads.rename(columns={'reads': 'previous_reads'}, inplace=True)
                        except Exception as e:
                            logger.debug(f"Could not load cached latest reads, will try CSV: {e}")
                            # Fallback to CSV if cache doesn't exist
                            if os.path.exists(time_coverage_file):
                                try:
                                    # Only read last chunk for efficiency
                                    existing_time_df = pd.read_csv(time_coverage_file)
                                    if not existing_time_df.empty:
                                        latest_timestamp = existing_time_df['timestamp'].max()
                                        previous_cumulative_reads = existing_time_df[
                                            existing_time_df['timestamp'] == latest_timestamp
                                        ][['chrom', 'startpos', 'endpos', 'name', 'reads']].copy()
                                        previous_cumulative_reads.rename(columns={'reads': 'previous_reads'}, inplace=True)
                                except Exception as e2:
                                    logger.warning(f"Error loading existing target_coverage_time.csv: {e2}")
                    
                    # Merge new read counts with target coverage data
                    target_coverage_with_reads = target_coverage_df.merge(
                        new_read_counts_df[['chrom', 'startpos', 'endpos', 'name', 'reads']],
                        on=['chrom', 'startpos', 'endpos', 'name'],
                        how='left'
                    )
                    # Fill missing reads with 0
                    target_coverage_with_reads['reads'] = target_coverage_with_reads['reads'].fillna(0).astype(int)
                    
                    # Accumulate with previous cumulative reads if available
                    if previous_cumulative_reads is not None:
                        target_coverage_with_reads = target_coverage_with_reads.merge(
                            previous_cumulative_reads,
                            on=['chrom', 'startpos', 'endpos', 'name'],
                            how='left'
                        )
                        target_coverage_with_reads['previous_reads'] = target_coverage_with_reads['previous_reads'].fillna(0).astype(int)
                        # Add new reads to previous cumulative reads
                        target_coverage_with_reads['reads'] = (
                            target_coverage_with_reads['reads'] + target_coverage_with_reads['previous_reads']
                        )
                        target_coverage_with_reads.drop(columns=['previous_reads'], inplace=True)
                    
                    # Calculate normalized reads (reads per length)
                    target_coverage_with_reads['reads_per_length'] = (
                        target_coverage_with_reads['reads'] / target_coverage_with_reads['length']
                    )
                    
                    # Add timestamp
                    if self.simtime and timestamp:
                        current_timestamp = timestamp * 1000
                    else:
                        current_timestamp = time.time() * 1000
                    target_coverage_with_reads['timestamp'] = current_timestamp
                    
                    # Reorder columns: chrom, startpos, endpos, name, length, coverage, bases, timestamp, reads, reads_per_length
                    target_coverage_with_reads = target_coverage_with_reads[
                        ['chrom', 'startpos', 'endpos', 'name', 'length', 'coverage', 'bases',
                         'timestamp', 'reads', 'reads_per_length']
                    ]
                    
                    # Save latest cumulative reads to cache for next time (fast access)
                    try:
                        target_coverage_with_reads[['chrom', 'startpos', 'endpos', 'name', 'reads']].to_parquet(
                            latest_reads_cache, index=False
                        )
                    except Exception as e:
                        logger.debug(f"Could not save latest reads cache: {e}")
                    
                    # Append to CSV using append mode (much faster than reading entire file)
                    try:
                        # Check if file exists to determine if we need header
                        file_exists = os.path.exists(time_coverage_file)
                        target_coverage_with_reads.to_csv(
                            time_coverage_file, 
                            mode='a', 
                            header=not file_exists,
                            index=False
                        )
                    except Exception as e:
                        logger.warning(f"Error appending to target_coverage_time.csv: {e}")
                        # Fallback: read and concat (slower but works)
                        if os.path.exists(time_coverage_file):
                            try:
                                existing_time_df = pd.read_csv(time_coverage_file)
                                target_coverage_with_reads = pd.concat(
                                    [existing_time_df, target_coverage_with_reads],
                                    ignore_index=True
                                )
                                target_coverage_with_reads.to_csv(time_coverage_file, index=False)
                            except Exception as e2:
                                logger.error(f"Error in fallback CSV write: {e2}")
                    
                    logger.info(f"Saved target coverage with timestamp and cumulative read counts: {time_coverage_file}")
                    target_result.processing_steps.append("target_coverage_time_saved")
                else:
                    logger.warning("No read counts extracted, skipping target_coverage_time.csv")

                # Step 13: Identify targets exceeding threshold
                run_list = target_coverage_df[
                    target_coverage_df["coverage"].ge(self.callthreshold)
                ]

                # Load existing targets exceeding threshold count
                targets_exceeding_file = os.path.join(
                    sample_output_dir, "targets_exceeding_threshold_count.txt"
                )
                existing_targets_count = 0
                if os.path.exists(targets_exceeding_file):
                    try:
                        with open(targets_exceeding_file, "r") as f:
                            existing_targets_count = int(f.read().strip())
                    except (ValueError, IOError):
                        existing_targets_count = 0

                # Always generate targets_exceeding_threshold.bed for potential SNP analysis
                if len(run_list) > 0:
                    logger.info(f"Found {len(run_list)} regions exceeding threshold")

                    # Save updated count
                    with open(targets_exceeding_file, "w") as f:
                        f.write(str(len(run_list)))

                    # Generate BED file for regions exceeding threshold
                    run_list[["chrom", "startpos", "endpos"]].to_csv(
                        os.path.join(
                            sample_output_dir, "targets_exceeding_threshold.bed"
                        ),
                        sep="\t",
                        header=None,
                        index=None,
                    )

                    logger.info(
                        f"Targets exceeding threshold saved: {len(run_list)} regions"
                    )

                    # Log the results for debugging (no storage to avoid memory leaks)
                    logger.info(
                        f"Generated BED file with {len(run_list)} regions exceeding threshold"
                    )
                else:
                    logger.info("No regions exceed coverage threshold")

                    # Create empty BED file to avoid errors
                    empty_bed_file = os.path.join(
                        sample_output_dir, "targets_exceeding_threshold.bed"
                    )
                    with open(empty_bed_file, "w") as f:
                        # Write minimal BED header
                        pass

                    logger.info("Created empty targets_exceeding_threshold.bed file")

                    # Log the empty file creation (no storage to avoid memory leaks)
                    logger.info("Created empty BED file - no regions exceed threshold")

                # Log reference genome status for SNP calling
                if self.reference:
                    logger.info(
                        f"Reference genome found: {self.reference} - SNP calling will be available"
                    )
                else:
                    logger.info(
                        "No reference genome provided - SNP calling will not be available"
                    )

                target_result.processing_steps.append("threshold_analysis_complete")

                # Step 13: Perform target analysis
                logger.info("Performing target analysis...")
                analysis_results = self._perform_target_analysis(
                    file_path, temp_dir, metadata, newcovdf, bedcovdf
                )
                target_result.analysis_results = analysis_results
                target_result.processing_steps.append("analysis_completed")

                # Step 14: Generate output files
                logger.info("Generating output files...")
                output_files = self._generate_output_files(
                    sample_output_dir,
                    analysis_counter,
                    analysis_results,
                    newcovdf,
                    bedcovdf,
                    logger,
                    metadata,
                )
                # No JSON files created - data is accumulated in CSV files
                target_result.target_data_path = None
                target_result.target_plot_path = None
                target_result.processing_steps.append("output_generated")

                # Step 15: Update analysis counter
                analysis_counter += 1
                self._save_analysis_counter(
                    sample_id, analysis_counter, self.work_dir, logger
                )
                target_result.processing_steps.append("counter_updated")

                # Step 16: Force garbage collection
                gc.collect()

                target_result.processing_steps.append("target_analysis_complete")
                logger.info(f"Target analysis complete for {sample_id}")
                logger.info(f"   Results: {len(analysis_results)} items")
                logger.info(f"   Output files: {len(output_files)} files")
                logger.info(
                    f"   Targets exceeding threshold: {len(run_list) if 'run_list' in locals() else 0}"
                )

                return target_result

        except Exception as e:
            error_details = f"Error in target analysis for {sample_id}: {str(e)}"
            logger.error(error_details)
            target_result.error_message = error_details
            target_result.processing_steps.append("analysis_failed")
            return target_result

    def _perform_target_analysis(
        self,
        file_path: str,
        temp_dir: str,
        metadata: Dict[str, Any],
        newcovdf: pd.DataFrame,
        bedcovdf: pd.DataFrame,
    ) -> Dict[str, Any]:
        """
        Perform the actual target analysis.

        Args:
            file_path: Path to input file
            temp_dir: Temporary directory for processing
            metadata: File metadata
            newcovdf: Genome coverage DataFrame
            bedcovdf: Target coverage DataFrame

        Returns:
            Dictionary containing analysis results
        """
        logger = logging.getLogger("robin.target")

        # This is where the actual target analysis logic would go
        # For now, we'll create a placeholder implementation that uses the coverage data

        results = {
            "analysis_type": self.config.get("analysis_type", "standard"),
            "input_file": os.path.basename(file_path),
            "sample_id": metadata.get("sample_id", "unknown"),
            "timestamp": time.time(),
            "targets_found": 0,
            "targets_analyzed": [],
            "coverage_summary": {
                "genome_regions": len(newcovdf) if newcovdf is not None else 0,
                "target_regions": len(bedcovdf) if bedcovdf is not None else 0,
                "total_bases_covered": 0,
                "average_coverage": 0.0,
            },
            "summary_stats": {
                "total_targets": 0,
                "significant_targets": 0,
                "average_score": 0.0,
            },
        }

        # Calculate coverage statistics
        if newcovdf is not None and not newcovdf.empty:
            results["coverage_summary"]["total_bases_covered"] = newcovdf.get(
                "numreads", pd.Series([0])
            ).sum()
            results["coverage_summary"]["average_coverage"] = newcovdf.get(
                "numreads", pd.Series([0])
            ).mean()

        # Analyze target regions
        if bedcovdf is not None and not bedcovdf.empty:
            results["targets_found"] = len(bedcovdf)
            results["summary_stats"]["total_targets"] = len(bedcovdf)

            # Calculate significance based on coverage
            significant_threshold = self.config.get("threshold", 0.5)
            significant_targets = bedcovdf[
                bedcovdf["bases"] > significant_threshold * bedcovdf["bases"].max()
            ]
            results["summary_stats"]["significant_targets"] = len(significant_targets)

            # Add target analysis results (list comp inlined in 3.12)
            max_bases = bedcovdf["bases"].max()
            results["targets_analyzed"] = [
                {
                    "target_id": row.get("name", f"target_{idx+1:03d}"),
                    "position": f"{row['chrom']}:{row['startpos']}-{row['endpos']}",
                    "coverage": row["bases"],
                    "significance": (
                        "high" if row["bases"] > significant_threshold * max_bases else "medium"
                    ),
                }
                for idx, row in bedcovdf.iterrows()
            ]

            # Calculate average score
            if results["targets_analyzed"]:
                scores = [t["coverage"] for t in results["targets_analyzed"]]
                results["summary_stats"]["average_score"] = sum(scores) / len(scores)

        logger.info(f"Analysis completed: {results['targets_found']} targets found")
        logger.info(f"Coverage summary: {results['coverage_summary']}")
        return results

    def _generate_output_files(
        self,
        sample_dir: str,
        analysis_counter: int,
        analysis_results: Dict[str, Any],
        newcovdf: pd.DataFrame,
        bedcovdf: pd.DataFrame,
        logger,
        metadata: Dict[str, Any],
    ) -> Dict[str, str]:
        """
        Generate minimal output files for the analysis.

        Args:
            sample_dir: Sample output directory
            analysis_counter: Current analysis counter
            analysis_results: Analysis results
            newcovdf: Genome coverage DataFrame
            bedcovdf: Target coverage DataFrame
            logger: Logger instance
            metadata: File metadata

        Returns:
            Dictionary with paths to generated files
        """
        output_files = {}

        try:
            # No JSON files needed - accumulated data is already saved in CSV files
            # The coverage data is persisted in:
            # - coverage_main.csv (genome coverage)
            # - bed_coverage_main.csv (target coverage)
            # - coverage_time_chart.npy (coverage over time)
            # - target_coverage.csv (calculated target coverage)

            logger.info(
                "No additional output files needed - data accumulated in CSV files"
            )
            return output_files

        except Exception as e:
            logger.error(f"Error in output file generation: {e}")
            return {}

    def _load_analysis_counter(self, sample_id: str, work_dir: str, logger) -> int:
        """Load the analysis counter for a sample from disk"""
        counter_file = os.path.join(work_dir, sample_id, "target_analysis_counter.txt")
        if os.path.exists(counter_file):
            try:
                with open(counter_file, "r") as f:
                    return int(f.read().strip())
            except (ValueError, IOError) as e:
                logger.warning(f"Error loading counter for {sample_id}: {e}")
        return 0

    def _save_analysis_counter(
        self, sample_id: str, counter: int, work_dir: str, logger
    ) -> None:
        """Save the analysis counter for a sample to disk"""
        counter_file = os.path.join(work_dir, sample_id, "target_analysis_counter.txt")
        try:
            # Create sample directory if it doesn't exist
            os.makedirs(os.path.dirname(counter_file), exist_ok=True)
            with open(counter_file, "w") as f:
                f.write(str(counter))
        except IOError as e:
            logger.error(f"Error saving counter for {sample_id}: {e}")

    def _load_existing_coverage_data(
        self, sample_dir: str, filename: str, logger
    ) -> Optional[pd.DataFrame]:
        """Load existing coverage data from CSV file if it exists"""
        file_path = os.path.join(sample_dir, filename)
        if os.path.exists(file_path):
            try:
                df = pd.read_csv(file_path)
                logger.debug(
                    f"Loaded existing coverage data from {filename}: {df.shape}"
                )

                # For bed_coverage_main.csv, extract only the raw data columns needed for merging
                if (
                    filename == "bed_coverage_main.csv"
                    and "length" in df.columns
                    and "coverage" in df.columns
                ):
                    # Extract only the raw data columns: chrom,startpos,endpos,name,bases
                    df = df[["chrom", "startpos", "endpos", "name", "bases"]]
                    logger.debug(
                        f"Extracted raw data columns from {filename}: {df.shape}"
                    )

                return df
            except Exception as e:
                logger.warning(
                    f"Error loading existing coverage data from {filename}: {e}"
                )
                return None
        return None

    def _load_existing_coverage_over_time(
        self, sample_dir: str, logger
    ) -> Optional[np.ndarray]:
        """Load existing coverage over time data from numpy file if it exists"""
        file_path = os.path.join(sample_dir, "coverage_time_chart.npy")
        if os.path.exists(file_path):
            try:
                data = np.load(file_path)
                logger.debug(f"Loaded existing coverage over time data: {data.shape}")
                return data
            except Exception as e:
                logger.warning(f"Error loading existing coverage over time data: {e}")
                return None
        return None
    
    def _extraction_regions(self, sample_id: str, out_path: str) -> str:
        """Write the regions whose reads go to target.bam: the panel, plus the master BED's
        regions when there is one. The panel is always in (TGT-2), even when the master
        BED was written without it."""
        rows = [(c, s, e) for c, s, e, _ in read_bed_targets(self.bedfile)]
        master = self._get_master_bed_path(sample_id)
        if master and os.path.exists(master) and os.path.abspath(master) != os.path.abspath(self.bedfile):
            rows += [(c, s, e) for c, s, e, _ in read_bed_targets(master)]
        with open(out_path, "w") as f:
            for chrom, start, end in rows:
                f.write(f"{chrom}\t{start}\t{end}\n")
        return out_path

    def _write_batch_bam(self, sample_id: str, batch_dir: str, source_bams: List[str]) -> Optional[str]:
        """Every source BAM's records in the extraction regions (samtools view -L), merged
        into one indexed BAM in the batch folder (TGT-1: every BAM, from its arrival). The
        commit moves it into the sample folder, where finalize looks, only once the
        batch's coverage is written; so target.bam never gets reads the coverage misses."""
        if not source_bams:
            return None
        regions = self._extraction_regions(sample_id, os.path.join(batch_dir, "regions.bed"))
        parts = []
        for i, bam in enumerate(source_bams):
            part = os.path.join(batch_dir, f"part_{i:04d}.bam")
            pysam.view("-b", "-L", regions, "-o", part, bam, catch_stdout=False)
            parts.append(part)
        tmp = os.path.join(batch_dir, "batch.bam")
        if len(parts) > 1:
            pysam.merge("-f", "-o", tmp, *parts)
        else:
            os.replace(parts[0], tmp)
        pysam.index(tmp)
        return tmp

    def _commit_batch(
        self,
        sample_dir: str,
        contigs: pd.DataFrame,
        targets: pd.DataFrame,
        stretches: pd.DataFrame,
        timestamps: List[float],
        batch_reads: Optional[pd.DataFrame],
        logger,
        union: Optional[pd.DataFrame] = None,
    ) -> Dict[str, Any]:
        """Add one batch to the sample's totals and prepare every file it changes in
        ``<sample>/_commit/`` (the caller makes it), writing nothing in the sample folder.
        Returns the coverage numbers, and in ``steps`` the install list's steps that put the
        prepared files in place (TGT-3). The caller holds the accumulation lock for the whole
        call (TGT-4). ``union`` is the batch's aligned bases over the merged panel, or None
        when a BAM in it was staged before batch 8."""
        prep = os.path.join(sample_dir, _COMMIT_DIR)
        cov_path = os.path.join(sample_dir, "coverage_main.csv")
        bed_path = os.path.join(sample_dir, "bed_coverage_main.csv")
        stretch_path = os.path.join(sample_dir, "_covered_stretches.parquet")
        existing_cov = self._load_existing_coverage_data(sample_dir, "coverage_main.csv", logger)
        existing_bed = self._load_existing_coverage_data(sample_dir, "bed_coverage_main.csv", logger)
        for path, df in ((cov_path, existing_cov), (bed_path, existing_bed)):
            if os.path.exists(path) and df is None:
                raise RuntimeError(f"Can't read {path}; not overwriting the sample's totals")

        # Contigs: numreads and aligned bases add up.
        keys = ["#rname", "startpos", "endpos"]
        frames = [contigs[keys + ["numreads", "bases"]]]
        if existing_cov is not None and len(existing_cov):
            old = existing_cov.copy()
            if "bases" not in old.columns:
                # Written before step 6: depth from samtools' printed mean depth.
                old["bases"] = (old["meandepth"] * (old["endpos"] - old["startpos"] + 1)).round()
            frames.insert(0, old[keys + ["numreads", "bases"]])
        totals = pd.concat(frames, ignore_index=True).groupby(keys, as_index=False, sort=False)[
            ["numreads", "bases"]
        ].sum()
        totals["numreads"] = totals["numreads"].astype("int64")
        totals["bases"] = totals["bases"].round().astype("int64")
        length = totals["endpos"] - totals["startpos"] + 1
        totals["meandepth"] = totals["bases"] / length

        # Covered bases: from the union of covered stretches. A folder written before
        # step 6 holds per-BAM counts summed together, which can't be undone, so its
        # covered bases stay unknown.
        breadth_known = existing_cov is None or os.path.exists(stretch_path)
        all_stretches = None
        if breadth_known:
            previous = pd.read_parquet(stretch_path) if os.path.exists(stretch_path) else None
            all_stretches = union_stretches([previous, stretches])
            per_contig = covered_bases(all_stretches)
            totals["covbases"] = [per_contig.get(str(c), 0) for c in totals["#rname"]]
        else:
            logger.warning(
                f"{cov_path} was written before covered bases could be summed; they stay unknown"
            )
            totals["covbases"] = np.nan
        totals = totals[["#rname", "startpos", "endpos", "numreads", "covbases", "meandepth", "bases"]]
        global_depth = float(totals["bases"].sum() / length.sum()) if length.sum() > 0 else 0.0

        # Targets: aligned bases add up; depth over a 0-based, half-open BED row.
        tkeys = ["chrom", "startpos", "endpos", "name"]
        tframes = [targets[tkeys + ["bases"]]]
        if existing_bed is not None and len(existing_bed):
            tframes.insert(0, existing_bed[tkeys + ["bases"]])
        target_df = pd.concat(tframes, ignore_index=True).groupby(tkeys, as_index=False, sort=False)[
            "bases"
        ].sum()
        target_df["bases"] = target_df["bases"].astype("int64")
        target_df["length"] = target_df["endpos"] - target_df["startpos"]
        target_df["coverage"] = target_df["bases"] / target_df["length"]
        target_df = target_df[["chrom", "startpos", "endpos", "name", "length", "coverage", "bases"]]
        run_list = target_df[target_df["coverage"].ge(self.callthreshold)]

        # The merged panel: each panel base once (batch 8). Known only while every batch of
        # the sample carried its merged-panel counts; a folder whose totals exist without the
        # file, or a batch without them, leaves it unknown and the readers use the row sum.
        union_path = os.path.join(sample_dir, UNION_FILE)
        union_total = None
        if union is not None and (existing_bed is None or os.path.exists(union_path)):
            ukeys = ["chrom", "start", "end"]
            batch_union = union[ukeys + ["bases"]].astype({"chrom": str})
            previous_union = None
            if os.path.exists(union_path):
                previous_union = pd.read_csv(union_path, dtype={"chrom": str})
                # The file is the previous generation only if it matches the row totals it was
                # written beside (a commit that failed after the rows were written leaves it
                # behind), and only if it covers the same merged panel (a panel changed mid-run).
                in_step = (
                    "rows_bases" in previous_union.columns
                    and existing_bed is not None
                    and "bases" in existing_bed.columns
                    and len(previous_union)
                    and int(previous_union["rows_bases"].iloc[0]) == int(existing_bed["bases"].sum())
                )
                if not in_step or _regions(previous_union) != _regions(batch_union):
                    logger.warning(
                        f"{union_path} doesn't match the sample's row totals or panel; the panel "
                        "depth falls back to the row sum"
                    )
                    previous_union = None
                    batch_union = None
            if batch_union is not None:
                uframes = [batch_union] if previous_union is None else [previous_union[ukeys + ["bases"]], batch_union]
                union_total = pd.concat(uframes, ignore_index=True).groupby(
                    ukeys, as_index=False, sort=False
                )["bases"].sum()
                union_total["bases"] = union_total["bases"].astype("int64")
                union_total["rows_bases"] = int(target_df["bases"].sum())

        chart_old = self._load_existing_coverage_over_time(sample_dir, logger)
        current_timestamp = max(timestamps) if timestamps else time.time() * 1000
        row = np.array([[current_timestamp, global_depth]])
        chart = np.vstack([chart_old, row]) if chart_old is not None else row

        # Each file is prepared whole in _commit/, under its own name; the install list's
        # steps put them in place: replaces, then a removal, then the time series' append.
        replaces, removes, appends = [], [], []

        def prepared(name: str, mode: str):
            replaces.append({"do": "replace", "file": name})
            return atomic_open(os.path.join(prep, name), mode)

        def write_csv(df: pd.DataFrame, name: str) -> None:
            with prepared(name, "w") as f:
                df.to_csv(f, index=False)

        if all_stretches is not None:
            with prepared(os.path.basename(stretch_path), "wb") as f:
                all_stretches.to_parquet(f, index=False)
        write_csv(totals, "coverage_main.csv")
        write_csv(target_df, "bed_coverage_main.csv")
        if union_total is not None:
            write_csv(union_total, UNION_FILE)
        elif os.path.exists(union_path):
            removes.append({"do": "remove", "file": UNION_FILE})  # no merged-panel counts: unknown from now on
        write_csv(target_df, "target_coverage.csv")
        with prepared("coverage_time_chart.npy", "wb") as f:
            np.save(f, chart)

        # Reads per target, cumulative, and the time series the GUI plots.
        time_coverage_file = os.path.join(sample_dir, "target_coverage_time.csv")
        latest_reads_cache = os.path.join(sample_dir, "_target_coverage_latest_reads.parquet")
        if batch_reads is not None:
            previous = None
            if os.path.exists(latest_reads_cache):
                try:
                    previous = pd.read_parquet(latest_reads_cache)
                except Exception:
                    previous = None
            if previous is None and os.path.exists(time_coverage_file):
                try:
                    time_df = pd.read_csv(time_coverage_file)
                    if not time_df.empty:
                        latest = time_df[time_df["timestamp"] == time_df["timestamp"].max()]
                        previous = latest[tkeys + ["reads"]].copy()
                except Exception:
                    previous = None
            with_reads = target_df.merge(batch_reads[tkeys + ["reads"]], on=tkeys, how="left")
            with_reads["reads"] = with_reads["reads"].fillna(0).astype(int)
            if previous is not None:
                previous = previous.rename(columns={"reads": "previous_reads"})
                with_reads = with_reads.merge(previous[tkeys + ["previous_reads"]], on=tkeys, how="left")
                with_reads["reads"] = with_reads["reads"] + with_reads["previous_reads"].fillna(0).astype(int)
                with_reads = with_reads.drop(columns=["previous_reads"])
            with_reads["reads_per_length"] = with_reads["reads"] / with_reads["length"]
            with_reads["timestamp"] = current_timestamp
            with_reads = with_reads[
                ["chrom", "startpos", "endpos", "name", "length", "coverage", "bases",
                 "timestamp", "reads", "reads_per_length"]
            ]
            with prepared(os.path.basename(latest_reads_cache), "wb") as f:
                with_reads[tkeys + ["reads"]].to_parquet(f, index=False)
            # Only the new rows: the install cuts the file back to its size now, then adds
            # them, so an append a stop cut off part-way is put right by doing it again.
            file_exists = os.path.exists(time_coverage_file)
            with atomic_open(os.path.join(prep, _TIME_SERIES), "w") as f:
                with_reads.to_csv(f, header=not file_exists, index=False)
            appends.append({
                "do": "append", "file": _TIME_SERIES,
                "size": os.path.getsize(time_coverage_file) if file_exists else 0,
            })

        with prepared("targets_exceeding_threshold_count.txt", "w") as f:
            f.write(str(len(run_list)))
        with prepared("targets_exceeding_threshold.bed", "w") as f:
            run_list[["chrom", "startpos", "endpos"]].to_csv(f, sep="\t", header=False, index=False)
        return {
            "coverage": global_depth,
            "targets_exceeding_threshold": len(run_list),
            "steps": replaces + removes + appends,
        }

    def _finish_commit(self, sample_dir: str, logger) -> List[Dict[str, str]]:
        """Finish a commit a stop left part-way (TGT-3). The caller holds the accumulation
        lock. With an install list, do its steps (each is safe to do twice), then its syncs,
        then delete it and _commit/; without one, delete _commit/: nothing was installed.
        A list that can't be read raises: nothing is guessed. Returns failed-file entries
        for BAMs whose reads did not reach target.bam. A claim folder it keeps holds
        incomplete.txt, which every finalize lists (_kept_claims)."""
        prep = os.path.join(sample_dir, _COMMIT_DIR)
        list_path = os.path.join(prep, _INSTALL_LIST)
        if os.path.exists(list_path):
            try:
                with open(list_path) as f:
                    steps = json.load(f)["steps"]
                if not isinstance(steps, list) or not all(isinstance(s, dict) and "do" in s for s in steps):
                    raise ValueError("no list of steps")
            except Exception as e:
                raise RuntimeError(
                    f"Can't read the install list {list_path} ({type(e).__name__}: {e}). The sample's "
                    "coverage files may be part-way through adding a batch, so nothing more is added "
                    "until it is dealt with"
                ) from e
            logger.warning(f"Finishing a coverage commit that a stop left part-way in {sample_dir}")
            entries, problem = _install(sample_dir, steps, self.work_dir, logger)
            _put_on_disk(sample_dir, steps)
            if problem is None:  # after the syncs: the clear holds only once the move is on disk
                _clear_left_out(self.work_dir, steps)
            _delete_list(sample_dir)  # after the clear: a clear that can't be written keeps the list
            return entries
        if os.path.isdir(prep):
            shutil.rmtree(prep)  # never listed: nothing was installed
        return []

    def accumulate_staged_files(self, sample_id: str, force: bool = False) -> Dict[str, Any]:
        """Add staged BAMs to the sample's totals.

        Claiming staged sets and committing a batch each hold the accumulation lock;
        the commit reads the totals, adds the batch and writes every file in one locked
        step, so two workers can't overwrite each other (TGT-4). Loading the batch,
        counting reads and writing its batch BAM run outside the lock.
        """
        logger = logging.getLogger("robin.target")
        lock_file = self._get_lock_file(sample_id, "accumulation")
        staging_dir = self._get_staging_dir(sample_id)
        sample_dir = os.path.join(self.work_dir, sample_id)
        start_time = time.time()
        claim = None
        claimed, batch_dir, committed = [], None, False
        taken = []  # the claimed sets with a file moved into the claim folder
        # A BAM whose coverage was not added, or whose reads are not in target.bam: one
        # failed-file entry each, which the handler shows (batch 10e)
        file_failures: List[Dict[str, str]] = []
        listed = set()
        # Each claimed set's BAM and key, read from _staging/ before any of its files move: a
        # set moved part-way keeps them in either place
        who: Dict[str, Tuple[str, Optional[str]]] = {}

        def failed(n: str, text: str) -> None:
            listed.add(n)
            error = f"{text} Its files are kept in _staging/{os.path.basename(batch_dir)}."
            bam, key = who.get(n) or (_staged_set_bam(batch_dir, n), _staged_set_key(batch_dir, n))
            file_failures.append({"file": bam, "error": error, "step": _COVERAGE_STEP})
            # In the ledger too, so every start reports it, whichever job found it
            _record_left_out(self.work_dir, key, bam, error)

        def not_added(result: Dict[str, Any]) -> Dict[str, Any]:
            # After the claim and before the commit, no claimed set's coverage was added; the
            # claim folder stays, and finalize lists it. Before the claim, nothing is lost:
            # the files are still staged for the next accumulation.
            if not committed:
                for n in taken:
                    if n not in listed:
                        failed(n, f"Its coverage was not added ({result['status']}: {result['error']}).")
            if batch_dir and not taken:
                try:  # nothing moved into it: no claim folder for finalize to list
                    os.rmdir(batch_dir)
                except OSError:
                    pass
            if file_failures:
                result["file_failures"] = file_failures
            return result

        def with_failures(result: Dict[str, Any]) -> Dict[str, Any]:
            if file_failures:
                result["file_failures"] = file_failures
            return result

        try:
            # --- Locked: finish a commit a stop left part-way, then claim complete staged sets ---
            with FileLock(lock_file, timeout=60.0):
                file_failures.extend(self._finish_commit(sample_dir, logger))
                counters = _complete_staged_counters(staging_dir)
                if not counters:
                    logger.info(f"No staged files to accumulate for {sample_id}")
                    return with_failures({"status": "no_files", "files_processed": 0})
                if not force and len(counters) < self.batch_size:
                    logger.info(
                        f"Skipping accumulation - only {len(counters)} files staged "
                        f"(threshold: {self.batch_size}, force={force})"
                    )
                    return with_failures({"status": "below_threshold", "files_pending": len(counters)})
                claimed = counters if force else counters[: self.batch_size]
                batch_dir = os.path.join(staging_dir, f"_batch_{uuid.uuid4()}")
                # Held until this claim ends, so finalize can wait for it (JOB-6).
                claim = FileLock(_claim_lock_path(batch_dir), timeout=0)
                claim.acquire()
                os.makedirs(batch_dir)
                for n in claimed:
                    who[n] = (_staged_set_bam(staging_dir, n), _staged_set_key(staging_dir, n))
                    for i, (kind, ext) in enumerate(_STAGED_KINDS):
                        name = f"{kind}_{n}.{ext}"
                        shutil.move(os.path.join(staging_dir, name), os.path.join(batch_dir, name))
                        if i == 0:
                            taken.append(n)  # its first file moved: no longer complete in staging
                    # Optional parts: a BAM staged before batch 8 has no merged-panel part, and
                    # one staged before batch 10e no key.
                    for name in (f"union_{n}.parquet", f"bam_key_{n}.txt"):
                        if os.path.exists(os.path.join(staging_dir, name)):
                            shutil.move(os.path.join(staging_dir, name), os.path.join(batch_dir, name))

            # --- Not locked: load the batch, count reads, write its batch BAM ---
            contigs, targets, stretch_frames, timestamps, source_bams = [], [], [], [], []
            union_frames = []
            used = []
            bam_refs = []  # each counted BAM and its key, for the install list's batch_bam step
            for n in claimed:
                def staged(kind: str, ext: str) -> str:
                    return os.path.join(batch_dir, f"{kind}_{n}.{ext}")
                try:
                    c = pd.read_parquet(staged("coverage", "parquet"))
                    t = pd.read_parquet(staged("bedcov", "parquet"))
                    s = pd.read_parquet(staged("stretches", "parquet"))
                    with open(staged("timestamp", "txt")) as f:
                        ts = float(f.read().strip())
                    with open(staged("source_bam", "txt")) as f:
                        src = f.read().strip()
                except Exception as e:
                    logger.warning(f"Error loading staged set {n} for {sample_id}: {e}")
                    failed(
                        n,
                        f"Its staged coverage could not be read ({type(e).__name__}: {e}), "
                        "so its coverage is not counted.",
                    )
                    continue
                # The merged-panel part is optional: unreadable counts as missing (the panel
                # depth falls back to the row sum); the BAM's other totals are still used.
                u = None
                if os.path.exists(staged("union", "parquet")):
                    try:
                        u = pd.read_parquet(staged("union", "parquet"))
                    except Exception as e:
                        logger.warning(f"Could not read the merged-panel part of staged set {n}: {e}")
                union_frames.append(u)
                contigs.append(c)
                targets.append(t)
                stretch_frames.append(s)
                timestamps.append(ts)
                source_bams.append(src)
                used.append(n)
                bam_refs.append({"file": src, "key": who[n][1]})
            if not contigs:
                logger.error("No staging files could be loaded successfully")
                _keep_incomplete_claim(batch_dir, [], claimed, None)
                return {"status": "load_failed", "files_attempted": len(claimed), "file_failures": file_failures}

            read_frames = []
            for bam_path in source_bams:
                counts = get_read_counts_per_target(bam_path, self.bedfile)
                if not counts.empty:
                    read_frames.append(counts)
            batch_reads = None
            if read_frames:
                batch_reads = pd.concat(read_frames, ignore_index=True).groupby(
                    ["chrom", "startpos", "endpos", "name"], as_index=False
                ).agg({"reads": "sum"})

            batch_bam, batch_bam_error, bam_problem = None, None, None
            try:
                batch_bam = self._write_batch_bam(sample_id, batch_dir, list(dict.fromkeys(source_bams)))
            except Exception as e:
                batch_bam_error = str(e) or type(e).__name__  # an error with no text is still one
                bam_problem = f"the batch BAM could not be written ({batch_bam_error})"
                logger.error(f"Could not write the batch BAM for {sample_id}; target.bam will miss these reads: {e}")

            # --- Locked: read the totals, add this batch, and prepare every file in _commit/;
            # write the install list, then do its steps: the files go in place, then the batch
            # BAM, then the claim folder is deleted. Until the list is written nothing has
            # changed, and the claim folder keeps the staged files: the coverage and
            # target.bam both miss them. Once it is written, the commit counts: a stop is
            # finished by the next accumulation or finalize (TGT-3).
            with FileLock(lock_file, timeout=120.0):
                file_failures.extend(self._finish_commit(sample_dir, logger))
                prep = os.path.join(sample_dir, _COMMIT_DIR)
                os.makedirs(prep)
                _sync_dir(sample_dir)  # _commit/'s name
                try:
                    summary = self._commit_batch(
                        sample_dir,
                        pd.concat(contigs, ignore_index=True),
                        pd.concat(targets, ignore_index=True),
                        union_stretches(stretch_frames),
                        timestamps,
                        batch_reads,
                        logger,
                        # None when any BAM in the batch has no merged-panel part, or the BAMs'
                        # merged panels differ (the panel changed while they were staged)
                        union=(
                            None
                            if any(u is None for u in union_frames)
                            or len({_regions(u) for u in union_frames}) > 1
                            else pd.concat(union_frames, ignore_index=True)
                        ),
                    )
                    steps = summary.pop("steps")
                    _sync_dir(prep)  # the prepared files' names
                    not_read = [n for n in claimed if n not in used]
                    if batch_bam_error or not_read:
                        # These reads are not in target.bam: keep them, and say so.
                        _keep_incomplete_claim(batch_dir, used, not_read, batch_bam_error)
                    if batch_bam:
                        for path in (batch_bam, f"{batch_bam}.bai"):
                            _sync_file(path)
                    # Before every list: the names in the claim folder (the batch BAM's,
                    # incomplete.txt's), then the claim folder's name and the sets moved into it
                    _sync_dir(batch_dir)
                    _sync_dir(staging_dir)
                    if batch_bam:
                        steps.append({
                            "do": "batch_bam",
                            "from": os.path.relpath(batch_bam, sample_dir),
                            "to": f"batch_{used[0]}_{uuid.uuid4().hex[:8]}.bam",
                            "bams": bam_refs,
                        })
                    elif bam_problem:
                        # The batch BAM could not be written: its BAMs' reads miss target.bam. In
                        # the list, so a finish after a stop reports them too.
                        steps.append({
                            "do": "keep_claim",
                            "dir": os.path.relpath(batch_dir, sample_dir),
                            "problem": bam_problem,
                            "bams": bam_refs,
                        })
                    if not (batch_bam_error or not_read):
                        steps.append({"do": "drop_claim", "dir": os.path.relpath(batch_dir, sample_dir)})
                    with atomic_open(os.path.join(prep, _INSTALL_LIST), "w") as f:
                        json.dump({"v": 1, "steps": steps}, f)
                    _sync_dir(prep)  # the list's name
                except Exception:
                    shutil.rmtree(prep, ignore_errors=True)  # never listed: nothing was installed
                    raise
                committed = True
                installed, moved_problem = _install(sample_dir, steps, self.work_dir, logger)
                file_failures.extend(installed)
                if moved_problem and not batch_bam_error:
                    batch_bam_error = moved_problem
                _retire(sample_dir, steps)

            gc.collect()
            elapsed = time.time() - start_time
            logger.info(f"Batch accumulation complete for {sample_id}: {len(used)} files in {elapsed:.2f}s")
            result = {
                "status": "success",
                "files_processed": len(used),
                "coverage": summary["coverage"],
                "targets_exceeding_threshold": summary["targets_exceeding_threshold"],
                "elapsed_time": elapsed,
            }
            if batch_bam_error:
                result["batch_bam_error"] = batch_bam_error
            return with_failures(result)

        except TimeoutError as e:
            logger.error(f"Could not acquire accumulation lock for {sample_id}: {e}")
            return not_added({"status": "lock_timeout", "error": str(e)})
        except Exception as e:
            logger.error(f"Error during batch accumulation for {sample_id}: {e}")
            import traceback
            logger.error(traceback.format_exc())
            return not_added({"status": "error", "error": str(e)})
        finally:
            if claim is not None:
                _end_claim(claim)


def process_single_file(
    file_path: str,
    metadata: Dict[str, Any],
    work_dir: str,
    logger,
    reference: Optional[str] = None,
    target_panel: str = "rCNS2",
) -> Dict[str, Any]:
    """
    Process a single file for target analysis using the complete pipeline.

    Args:
        file_path: Path to input file
        metadata: File metadata
        work_dir: Working directory
        logger: Logger instance
        reference: Optional reference genome path
        target_panel: Target panel type (rCNS2, AML, or custom panel)

    Returns:
        Dictionary with target analysis results
    """
    sample_id = metadata.get("sample_id", "unknown")
    logger.info(f"Starting target analysis for sample: {sample_id}")

    # Create sample-specific output directory
    sample_output_dir = os.path.join(work_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)
    logger.info(f"Created sample output directory: {sample_output_dir}")

    analysis_result = {
        "sample_id": sample_id,
        "file_path": file_path,
        "analysis_timestamp": time.time(),
        "target_data_path": None,
        "target_plot_path": None,
        "analysis_results": {},
        "processing_steps": [],
        "error_message": None,
        "coverage_data": {},
        "target_bam_path": None,
        "coverage_over_time": None,
    }

    try:
        # Initialize target analysis
        target_analysis = TargetAnalysis(work_dir=work_dir, target_panel=target_panel)

        # Set reference genome if provided
        if reference:
            target_analysis.reference = reference
            logger.info(f"Using reference genome: {reference}")

        # Process the file
        target_metadata = target_analysis.process_file(file_path, metadata)

        # Update analysis result with target metadata
        analysis_result.update(
            {
                "target_data_path": target_metadata.target_data_path,
                "target_plot_path": target_metadata.target_plot_path,
                "analysis_results": target_metadata.analysis_results,
                "processing_steps": target_metadata.processing_steps,
                "error_message": target_metadata.error_message,
                "coverage_data": target_metadata.coverage_data,
                "target_bam_path": target_metadata.target_bam_path,
                "coverage_over_time": target_metadata.coverage_over_time,
            }
        )

        if target_metadata.error_message:
            logger.error(f"Target analysis failed: {target_metadata.error_message}")
        else:
            logger.info(f"Target analysis complete for {sample_id}")
            logger.info(f"   Results: {len(target_metadata.analysis_results)} items")
            logger.info(
                f"   Processing steps: {', '.join(target_metadata.processing_steps)}"
            )
            logger.info(
                f"   Coverage data: genome={target_metadata.coverage_data.get('genome_coverage_shape', (0,0))}, targets={target_metadata.coverage_data.get('target_coverage_shape', (0,0))}"
            )
            if target_metadata.target_bam_path:
                logger.info(f"   Target BAM: {target_metadata.target_bam_path}")

        return analysis_result

    except Exception as e:
        error_details = f"Error in target analysis for {sample_id}: {str(e)}"
        logger.error(error_details)
        analysis_result["error_message"] = error_details
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def process_multiple_files(bam_paths, metadata_list, work_dir, logger, reference=None, target_panel=None):
    """
    Process multiple BAM files for target analysis using staged processing.
    
    This function processes multiple BAM files for the same sample using the
    existing staging infrastructure. Each file is processed individually and
    staged, then all staged files are accumulated in a single batch operation.

    Args:
        bam_paths: List of paths to BAM files
        metadata_list: List of metadata dictionaries (one per BAM file)
        work_dir: Working directory
        logger: Logger instance
        reference: Optional path to reference genome for SNP calling
        target_panel: Target panel type (rCNS2, AML, or custom panel)

    Returns:
        Dictionary with aggregated target analysis results
    """
    if not bam_paths or not metadata_list:
        raise ValueError("bam_paths and metadata_list must not be empty")
    
    if len(bam_paths) != len(metadata_list):
        raise ValueError("bam_paths and metadata_list must have the same length")
    
    # Get sample ID from first metadata (assuming all BAMs are from same sample)
    sample_id = metadata_list[0].get("sample_id", "unknown")
    
    logger.info(f"🎯 Starting multi-file target analysis for sample: {sample_id}")
    logger.info(f"Processing {len(bam_paths)} BAM files for sample {sample_id}")
    
    # Log essential metadata only
    for i, (bam_path, metadata) in enumerate(zip(bam_paths, metadata_list)):
        logger.debug(f"BAM file {i+1}: {os.path.basename(bam_path)}")

    analysis_result = {
        "sample_id": sample_id,
        "bam_paths": bam_paths,
        "analysis_timestamp": time.time(),
        "processing_steps": [],
        "error_message": None,
        "files_processed": 0,
        "total_files": len(bam_paths),
        "coverage_data": {},
        "target_bam_path": None,
        "coverage_over_time": None,
    }

    try:
        # Initialize target analysis with staging enabled
        target_analysis = TargetAnalysis(
            work_dir=work_dir, 
            target_panel=target_panel,
            batch_size=1,  # Force accumulation after each batch
            use_staging=True
        )
        
        # Set reference genome if provided
        if reference:
            target_analysis.reference = reference
            logger.info(f"Using reference genome: {reference}")

        # Process each BAM file individually using staging
        logger.info("Processing files with staging (fast path)")
        processed_files = 0
        # A BAM that could not be staged is not counted: a failed file (batch 10e)
        not_staged = []

        def not_counted(bam_path: str, why: str) -> None:
            not_staged.append({
                "file": bam_path,
                "error": f"Its coverage could not be prepared, so it is not counted: {why}",
                "step": _COVERAGE_STEP,
            })

        for i, (bam_path, metadata) in enumerate(zip(bam_paths, metadata_list)):
            logger.info(f"Processing BAM file {i+1}/{len(bam_paths)}: {os.path.basename(bam_path)}")

            try:
                # Process file with staging
                target_metadata, should_accumulate = target_analysis.process_file_with_staging(
                    bam_path, metadata
                )

                if target_metadata.error_message:
                    logger.warning(f"Error processing {os.path.basename(bam_path)}: {target_metadata.error_message}")
                    not_counted(bam_path, target_metadata.error_message)
                    continue

                processed_files += 1
                logger.debug(f"Successfully staged file {i+1}: {os.path.basename(bam_path)}")

            except Exception as e:
                logger.warning(f"Error processing {os.path.basename(bam_path)}: {e}")
                not_counted(bam_path, f"{type(e).__name__}: {e}")
                continue

        if processed_files == 0:
            analysis_result["error_message"] = "No files could be processed successfully"
            analysis_result["processing_steps"].append("no_files_processed")
            return analysis_result

        analysis_result["files_processed"] = processed_files
        analysis_result["processing_steps"].append("files_staged")

        # Force accumulation of all staged files
        logger.info(f"Accumulating {processed_files} staged files for sample {sample_id}")
        accumulation_result = target_analysis.accumulate_staged_files(
            sample_id, force=True
        )
        
        if accumulation_result.get("status") != "success":
            analysis_result["error_message"] = f"Accumulation failed: {accumulation_result.get('error', 'Unknown error')}"
            analysis_result["processing_steps"].append("accumulation_failed")
            return analysis_result

        analysis_result["processing_steps"].append("accumulation_complete")
        logger.info(f"Accumulation completed: {accumulation_result}")
        analysis_result["file_failures"] = not_staged + list(accumulation_result.get("file_failures") or [])

        # Load final accumulated data for result metadata
        sample_output_dir = os.path.join(work_dir, sample_id)
        
        # Load final coverage data
        final_covdf = target_analysis._load_existing_coverage_data(
            sample_output_dir, "coverage_main.csv", logger
        )
        final_bedcovdf = target_analysis._load_existing_coverage_data(
            sample_output_dir, "bed_coverage_main.csv", logger
        )
        final_coverage_over_time = target_analysis._load_existing_coverage_over_time(
            sample_output_dir, logger
        )
        
        # Store final results
        analysis_result["coverage_data"] = {
            "genome_coverage_shape": final_covdf.shape if final_covdf is not None else (0, 0),
            "target_coverage_shape": final_bedcovdf.shape if final_bedcovdf is not None else (0, 0),
            "coverage": accumulation_result.get("coverage", 0.0),
            "targets_exceeding_threshold": accumulation_result.get("targets_exceeding_threshold", 0),
        }
        
        analysis_result["target_bam_path"] = os.path.join(sample_output_dir, "target.bam")
        analysis_result["coverage_over_time"] = final_coverage_over_time
        
        analysis_result["processing_steps"].append("analysis_complete")
        logger.info(f"Multi-file target analysis completed for {sample_id}")
        logger.info(f"Files successfully processed: {analysis_result['files_processed']}/{analysis_result['total_files']}")
        logger.info(f"Final coverage: {analysis_result['coverage_data']['coverage']:.4f}")
        logger.info(f"Targets exceeding threshold: {analysis_result['coverage_data']['targets_exceeding_threshold']}")
        
        return analysis_result

    except Exception as e:
        logger.error(f"Error in multi-file target analysis for {sample_id}: {e}")
        analysis_result["error_message"] = str(e)
        analysis_result["processing_steps"].append("analysis_failed")
        return analysis_result


def target_handler(job, work_dir=None, reference=None, target_panel=None):
    """
    Handler function for target analysis jobs.
    This function processes files for target-specific analysis.

    Args:
        job: The workflow job containing file and metadata
        work_dir: Optional base directory for output (defaults to file directory)
        reference: Optional path to reference genome for SNP calling
        target_panel: Target panel type (required)
    """
    # Validate required parameters
    if not target_panel:
        raise ValueError("target_panel is required for target analysis")
    
    # Get job-specific logger
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
    
    # Check if this is a batched job
    batched_job = job.context.metadata.get("_batched_job")
    if batched_job:
        batch_size = batched_job.get_file_count()
        sample_id = batched_job.get_sample_id()
        batch_id = batched_job.batch_id
        logger.info(f"Processing target analysis batch: {batch_size} files for sample '{sample_id}' (batch_id: {batch_id})")
        
        # Get all filepaths in the batch
        filepaths = batched_job.get_filepaths()
        
        # Log individual files in the batch
        for i, filepath in enumerate(filepaths):
            logger.info(f"  Batch file {i+1}/{batch_size}: {os.path.basename(filepath)}")
        
        # Prepare metadata list for all BAM files in the batch (list comp inlined in 3.12)
        def _batch_metadata(i: int) -> dict:
            ctx = batched_job.contexts[i]
            sid = ctx.get_sample_id()
            return {**ctx.metadata.get("bam_metadata", {}), "sample_id": sid if sid != "unknown" else sample_id}
        metadata_list = [_batch_metadata(i) for i in range(len(filepaths))]

        # Determine work directory for the batch
        if work_dir is None:
            # Default to first BAM file directory
            batch_work_dir = os.path.dirname(filepaths[0])
        else:
            # Use specified work directory, create if it doesn't exist
            os.makedirs(work_dir, exist_ok=True)
            batch_work_dir = work_dir
            logger.debug(f"Using specified work directory: {batch_work_dir}")
        
        # Log and validate target panel
        job_panel = batched_job.contexts[0].metadata.get("target_panel", target_panel)
        if job_panel != target_panel:
            logger.warning(f"Panel mismatch: job metadata has '{job_panel}' but handler received '{target_panel}'. Using '{job_panel}' from metadata.")
            target_panel = job_panel
        logger.info(f"Using target panel: {target_panel}")
        
        # Debug: Log reference genome status
        if reference:
            logger.info(f"Reference genome provided to target_handler: {reference}")
        else:
            logger.info("No reference genome provided to target_handler")

        # Also check job metadata for reference
        job_reference = batched_job.contexts[0].metadata.get("reference")
        if job_reference:
            logger.info(f"Reference genome found in job metadata: {job_reference}")
            # Use job metadata reference if not provided directly
            if not reference:
                reference = job_reference
                logger.info(f"Using reference from job metadata: {reference}")
        else:
            logger.info("No reference genome found in job metadata")
        
        # Process all BAM files in the batch using the new aggregated function
        logger.info(f"Processing {batch_size} BAM files as aggregated batch for sample '{sample_id}'")
        batch_result = process_multiple_files(
            bam_paths=filepaths,
            metadata_list=metadata_list,
            work_dir=batch_work_dir,
            logger=logger,
            reference=reference,
            target_panel=target_panel
        )
        
        # Store batch results in job context (maintain compatibility with existing structure)
        job.context.add_metadata("target_analysis", {
            "batch_result": batch_result,  # Single aggregated result
            "batch_size": batch_size,
            "sample_id": sample_id,
            "batch_id": batch_id,
            "files_processed": batch_result.get("files_processed", batch_size),
            "total_files": batch_result.get("total_files", batch_size)
        })
        
        logger.info(f"Completed target analysis batch processing: {batch_size} files for sample '{sample_id}'")
        logger.info(f"Files successfully processed: {batch_result.get('files_processed', batch_size)}/{batch_result.get('total_files', batch_size)}")
        
        if batch_result.get("error_message"):
            logger.error(f"Batch processing completed with errors: {batch_result['error_message']}")
            job.context.add_error("target_analysis", batch_result["error_message"])
        else:
            logger.info("Batch processing completed successfully with aggregated target analysis")
            # BAMs whose coverage was not added: failed files, the job goes on (batch 10e)
            _add_file_failures(job, batch_result.get("file_failures"))
            job.context.add_result(
                "target_analysis",
                {
                    "status": "success",
                    "sample_id": sample_id,
                    "analysis_time": batch_result.get("analysis_timestamp", 0),
                    "targets_found": batch_result.get("coverage_data", {}).get("targets_exceeding_threshold", 0),
                    "processing_steps": batch_result.get("processing_steps", []),
                    "target_data_path": batch_result.get("target_data_path", ""),
                    "target_plot_path": batch_result.get("target_plot_path", ""),
                    "target_bam_path": batch_result.get("target_bam_path", ""),
                    "coverage_summary": batch_result.get("coverage_data", {}),
                    "files_processed": batch_result.get("files_processed", batch_size),
                    "total_files": batch_result.get("total_files", batch_size),
                },
            )
        
        # Update master.csv with panel information
        if sample_id and not batch_result.get("error_message"):
            try:
                from robin.analysis.master_csv_manager import MasterCSVManager
                
                # Update master.csv with panel information
                csv_manager = MasterCSVManager(batch_work_dir)
                csv_manager.update_analysis_panel(sample_id, target_panel)
                logger.info(f"Updated master.csv with panel '{target_panel}' for sample {sample_id}")
                
            except Exception as e:
                logger.warning(f"Could not update master.csv with panel info for {sample_id}: {e}")
        
        return
        
    else:
        # Single file processing (backward compatibility)
        try:
            file_path = job.context.filepath

            logger.info(f"Starting target analysis for: {os.path.basename(file_path)}")
            
            # Log and validate target panel
            job_panel = job.context.metadata.get("target_panel", target_panel)
            if job_panel != target_panel:
                logger.warning(f"Panel mismatch: job metadata has '{job_panel}' but handler received '{target_panel}'. Using '{job_panel}' from metadata.")
                target_panel = job_panel
            logger.info(f"Using target panel: {target_panel}")

            # Debug: Log reference genome status
            if reference:
                logger.info(f"Reference genome provided to target_handler: {reference}")
            else:
                logger.info("No reference genome provided to target_handler")

            # Also check job metadata for reference
            job_reference = job.context.metadata.get("reference")
            if job_reference:
                logger.info(f"Reference genome found in job metadata: {job_reference}")
                # Use job metadata reference if not provided directly
                if not reference:
                    reference = job_reference
                    logger.info(f"Using reference from job metadata: {reference}")
            else:
                logger.info("No reference genome found in job metadata")

            # Get metadata from preprocessing
            file_metadata = job.context.metadata.get("bam_metadata", {})

            # Determine work directory
            if work_dir is None:
                # Default to file directory
                work_dir = os.path.dirname(file_path)
            else:
                # Use specified work directory, create if it doesn't exist
                os.makedirs(work_dir, exist_ok=True)

            # Initialize target analysis with staging enabled
            target_analysis = TargetAnalysis(
                work_dir=work_dir, 
                target_panel=target_panel,
                batch_size=10,  # Accumulate every 10 files
                use_staging=True
            )
            
            # Set reference genome if provided
            if reference:
                target_analysis.reference = reference
            
            # Use fast staging-based processing
            logger.info("Using staging-based processing (fast path)")
            target_metadata, should_accumulate = target_analysis.process_file_with_staging(
                file_path, file_metadata
            )
            
            # Convert TargetMetadata to dict for storage
            result = {
                "sample_id": target_metadata.sample_id,
                "file_path": target_metadata.file_path,
                "analysis_timestamp": target_metadata.analysis_timestamp,
                "processing_steps": target_metadata.processing_steps,
                "error_message": target_metadata.error_message,
                "coverage_data": target_metadata.coverage_data,
            }
            
            # Store results in job context
            job.context.add_metadata("target_analysis", result)
            
            # Trigger accumulation if threshold reached
            if should_accumulate:
                logger.info("Accumulation threshold reached - running batch accumulation")
                accumulation_result = target_analysis.accumulate_staged_files(
                    target_metadata.sample_id, force=False
                )
                logger.info(f"Accumulation result: {accumulation_result}")
                job.context.add_metadata("accumulation_result", accumulation_result)
                # BAMs whose coverage was not added: failed files, the job goes on (batch 10e)
                _add_file_failures(job, accumulation_result.get("file_failures"))
            
            # Store flag for potential end-of-queue accumulation
            job.context.add_metadata("needs_final_accumulation", True)

            # Update master.csv with panel information
            if result.get("sample_id") and not result.get("error_message"):
                try:
                    from robin.analysis.master_csv_manager import MasterCSVManager
                    
                    # Determine work directory
                    if work_dir is None:
                        work_dir = os.path.dirname(file_path)
                    
                    # Update master.csv with panel information
                    csv_manager = MasterCSVManager(work_dir)
                    csv_manager.update_analysis_panel(result["sample_id"], target_panel)
                    logger.info(f"Updated master.csv with panel '{target_panel}' for sample {result['sample_id']}")
                    
                except Exception as e:
                    logger.warning(f"Could not update master.csv with panel info for {result.get('sample_id', 'unknown')}: {e}")

            if result.get("error_message"):
                job.context.add_error("target_analysis", result["error_message"])
                logger.error(f"Target analysis failed: {result['error_message']}")
            else:
                job.context.add_result(
                    "target_analysis",
                    {
                        "status": "success",
                        "sample_id": result.get("sample_id", "unknown"),
                        "analysis_time": result.get("analysis_timestamp", 0),
                        "targets_found": result.get("analysis_results", {}).get(
                            "targets_found", 0
                        ),
                        "processing_steps": result.get("processing_steps", []),
                        "target_data_path": result.get("target_data_path", ""),
                        "target_plot_path": result.get("target_plot_path", ""),
                        "target_bam_path": result.get("target_bam_path", ""),
                        "coverage_summary": result.get("analysis_results", {}).get(
                            "coverage_summary", {}
                        ),
                    },
                )
                logger.info(f"Target analysis complete for {os.path.basename(file_path)}")
                logger.info(f"Sample ID: {result.get('sample_id', 'unknown')}")
                logger.info(
                    f"Targets found: {result.get('analysis_results', {}).get('targets_found', 0)}"
                )

        except Exception as e:
            error_details = f"Error in target analysis for {job.context.filepath}: {str(e)}"
            job.context.add_error("target_analysis", error_details)
            logger.error(error_details)


_MERGED_TAG = "robin-merged-batch: "


def _merged_batch_names(target_bam_path: str) -> set:
    """The batch BAMs already merged into target.bam, from its header's @CO lines."""
    if not os.path.exists(target_bam_path):
        return set()
    with pysam.AlignmentFile(target_bam_path, "rb") as bam:
        comments = bam.header.to_dict().get("CO", [])
    return {c[len(_MERGED_TAG):] for c in comments if c.startswith(_MERGED_TAG)}


# A worker moves the staged files it claims into _staging/_batch_<id>/ and adds them
# outside the accumulation lock, holding this claim's own lock until it ends.
_CLAIM_WAIT_SECONDS = 300.0


def _claim_lock_path(claim_dir: str) -> str:
    sample_dir = os.path.dirname(os.path.dirname(claim_dir))
    return os.path.join(sample_dir, "_locks", f"claim{os.path.basename(claim_dir)}.lock")


def _end_claim(claim: "FileLock") -> None:
    try:
        os.remove(claim.lock_file)
    except FileNotFoundError:
        pass
    claim.release()


# The failed-file step name for a BAM whose coverage the target step could not add.
_COVERAGE_STEP = "target_coverage"


def _staged_set_bam(set_dir: str, n: str) -> str:
    """The BAM a staged set came from, or the set's own coverage file when that can't be read."""
    try:
        with open(os.path.join(set_dir, f"source_bam_{n}.txt")) as f:
            source = f.read().strip()
        if source:
            return source
    except OSError:
        pass
    return os.path.join(set_dir, f"coverage_{n}.parquet")


def _staged_set_key(set_dir: str, n: str) -> Optional[str]:
    """The key of the BAM a staged set came from: the one saved at staging, or for a set
    staged before batch 10e, read from the BAM's path while that still works. None if neither."""
    try:
        with open(os.path.join(set_dir, f"bam_key_{n}.txt")) as f:
            key = f.read().strip()
        if key:
            return key
    except OSError:
        pass
    try:
        with open(os.path.join(set_dir, f"source_bam_{n}.txt")) as f:
            return bam_key(f.read().strip())
    except OSError:
        return None


def _record_left_out(work_dir: str, key: Optional[str], bam: str, error: str) -> None:
    """Record in the work folder's ledger that the target step left this BAM out, so every
    start reports it, whichever job found it (batch 10e)."""
    if key is None:
        logging.getLogger("robin.target").warning(f"{bam} can't be identified: its failure is only logged")
        return
    record_bam_failure(work_dir, key, "target", error)


# A coverage commit prepares every file it changes in <sample>/_commit/, writes the list of
# steps that put them in place (install.json), then does them. A list a stop left is finished
# by the next accumulation or finalize, so a batch reaches every file, or none (TGT-3).
_COMMIT_DIR = "_commit"
_INSTALL_LIST = "install.json"
_TIME_SERIES = "target_coverage_time.csv"


def _sync_file(path: str) -> None:
    """Put a file's contents on disk."""
    fd = os.open(path, os.O_RDONLY)
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _sync_dir(path: str) -> None:
    """Put a folder's entries (the names of what is in it) on disk."""
    fd = os.open(path, os.O_RDONLY | getattr(os, "O_DIRECTORY", 0))
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _append_rows(path: str, size: int, rows_path: str) -> None:
    """Cut ``path`` back to ``size`` bytes, then add the rows in ``rows_path``. Safe to do
    twice, and it puts right an append a stop cut off part-way."""
    with open(rows_path, "rb") as f:
        rows = f.read()
    fd = os.open(path, os.O_WRONLY | os.O_CREAT, 0o666)
    with os.fdopen(fd, "wb") as f:
        end = f.seek(0, os.SEEK_END)
        if end < size:
            raise RuntimeError(
                f"{path} is shorter ({end} bytes) than when this batch's rows were prepared "
                f"({size} bytes); not guessing which rows it lost"
            )
        f.truncate(size)
        f.seek(size)
        f.write(rows)


def _claim_of(steps: List[dict]) -> Optional[str]:
    """The claim folder an install list moves the batch BAM from, keeps or deletes, relative
    to the sample folder; None when it names none."""
    for step in steps:
        if step["do"] == "batch_bam":
            return os.path.dirname(step["from"])
        if step["do"] in ("keep_claim", "drop_claim"):
            return step["dir"]
    return None


def _batch_bam_in_place(sample_dir: str, step: dict) -> bool:
    """The batch_bam step's BAM is in the sample folder, or target.bam lists it as merged (a
    finalize merges and deletes batch BAMs under its own lock, not the accumulation's, so it
    can do so before the list is deleted)."""
    return os.path.exists(os.path.join(sample_dir, step["to"])) or step["to"] in _merged_batch_names(
        os.path.join(sample_dir, "target.bam")
    )


def _install_batch_bam(sample_dir: str, step: dict) -> Optional[str]:
    """The install list's batch_bam step: move the index, then the BAM, from the claim folder
    to where finalize looks (finalize globs *.bam). Returns what went wrong, or None.

    Whether the BAM is in place is judged by where it is (_batch_bam_in_place), not by
    whether its move reported an error: over NFS a rename can take effect and still report
    one (rename(2), BUGS). A BAM gone from the claim folder and not in place is lost."""
    src, dst = os.path.join(sample_dir, step["from"]), os.path.join(sample_dir, step["to"])
    if os.path.exists(src):
        try:
            if os.path.exists(f"{src}.bai"):
                os.replace(f"{src}.bai", f"{dst}.bai")
            os.replace(src, dst)
        except Exception as e:
            if _batch_bam_in_place(sample_dir, step):
                logging.getLogger("robin.target").warning(
                    f"Moving {step['from']} into place reported an error, but it is in place: {e}"
                )
                return None
            if os.path.exists(f"{dst}.bai") and not os.path.exists(f"{src}.bai"):
                try:  # the claim folder keeps both
                    os.replace(f"{dst}.bai", f"{src}.bai")
                except OSError:
                    pass
            return f"the batch BAM could not be moved into place ({e})"
        return None
    if _batch_bam_in_place(sample_dir, step):
        return None
    return "the batch BAM is gone: it is neither in the sample folder nor merged into target.bam"


def _left_out_of_target_bam(work_dir: str, claim: str, bams, problem: str) -> List[Dict[str, str]]:
    """Failed-file entries, and ledger lines, for the BAMs of a batch whose batch BAM is not
    in target.bam: their coverage is counted, their reads are not in target.bam. Each names
    the claim folder (as the install list does), so that a later finish that puts the batch
    BAM in place clears it (_clear_left_out)."""
    error = (
        f"Its coverage is counted, but its reads are not in target.bam: {problem}. "
        f"Its files are kept in _staging/{os.path.basename(claim)}."
    )
    entries = []
    for bam in bams or []:
        entries.append({"file": bam["file"], "error": error, "step": _COVERAGE_STEP, "claim": claim})
        if bam.get("key"):
            record_bam_failure(work_dir, bam["key"], "target", error, claim=claim)
        else:
            logging.getLogger("robin.target").warning(
                f"{bam['file']} can't be identified: its failure is only logged"
            )
    return entries


def _install(sample_dir: str, steps: List[dict], work_dir: str, logger) -> Tuple[List[Dict[str, str]], Optional[str]]:
    """Do an install list's steps, in order; each is safe to do twice. Returns failed-file
    entries for BAMs whose reads did not reach target.bam, and why (None when they did)."""
    prep = os.path.join(sample_dir, _COMMIT_DIR)
    entries: List[Dict[str, str]] = []
    problem = None
    for step in steps:
        do = step["do"]
        if do == "replace":
            src = os.path.join(prep, step["file"])
            if os.path.exists(src):  # gone: renamed into place already
                os.replace(src, os.path.join(sample_dir, step["file"]))
        elif do == "remove":
            try:
                os.remove(os.path.join(sample_dir, step["file"]))
            except FileNotFoundError:
                pass
        elif do == "append":
            _append_rows(os.path.join(sample_dir, step["file"]), int(step["size"]), os.path.join(prep, step["file"]))
        elif do == "batch_bam":
            problem = _install_batch_bam(sample_dir, step)
            if problem:
                logger.error(f"Batch BAM {step['from']} in {sample_dir}: {problem}; target.bam misses its reads")
                entries += _left_out_of_target_bam(work_dir, os.path.dirname(step["from"]), step.get("bams"), problem)
                claim = os.path.join(sample_dir, os.path.dirname(step["from"]))
                os.makedirs(claim, exist_ok=True)
                counted = [os.path.basename(b["file"]) for b in step.get("bams") or []]
                _keep_incomplete_claim(claim, counted, [], problem)
        elif do == "keep_claim":  # its batch BAM could not be written; the claim folder stays
            problem = step["problem"]
            logger.error(f"Batch {step['dir']} in {sample_dir}: {problem}; target.bam misses its reads")
            entries += _left_out_of_target_bam(work_dir, step["dir"], step.get("bams"), problem)
            claim = os.path.join(sample_dir, step["dir"])
            if not os.path.exists(os.path.join(claim, "incomplete.txt")):  # its name did not reach the disk
                os.makedirs(claim, exist_ok=True)
                counted = [os.path.basename(b["file"]) for b in step.get("bams") or []]
                _keep_incomplete_claim(claim, counted, [], problem)
        elif do == "drop_claim":
            claim = os.path.join(sample_dir, step["dir"])
            if problem is None and os.path.isdir(claim):  # kept when its batch BAM is not in place
                shutil.rmtree(claim)
        else:
            raise ValueError(f"unknown install step {do!r}")
    return entries, problem


def _retire(sample_dir: str, steps: List[dict]) -> None:
    """Put what the install did on disk, then delete the list (and sync that), then _commit/.
    A sync that fails raises first, so the list stays and the next finish does it again."""
    _put_on_disk(sample_dir, steps)
    _delete_list(sample_dir)


def _clear_left_out(work_dir: str, steps: List[dict]) -> None:
    """A finish put the list's batch BAM in place, and on disk: the failures its commit
    recorded for the batch's BAMs (_left_out_of_target_bam) no longer hold. Raises when that
    can't be recorded, so the list stays and the next finish tries again."""
    for step in steps:
        if step["do"] == "batch_bam":
            keys = [b["key"] for b in step.get("bams") or [] if b.get("key")]
            if keys:
                clear_bam_failures(work_dir, keys, "target", os.path.dirname(step["from"]))


def _put_on_disk(sample_dir: str, steps: List[dict]) -> None:
    """Sync what the install did: the appended files, the sample folder's names (the renames,
    the removal, the moved batch BAM), then the claim folder and its parent."""
    for step in steps:
        if step["do"] == "append":
            _sync_file(os.path.join(sample_dir, step["file"]))
    _sync_dir(sample_dir)  # the renames, the removal, the moved batch BAM
    claim = _claim_of(steps)
    if claim:
        if os.path.isdir(os.path.join(sample_dir, claim)):
            _sync_dir(os.path.join(sample_dir, claim))  # kept: the batch BAM moved out of it
        _sync_dir(os.path.join(sample_dir, os.path.dirname(claim)))  # _staging/: the claim's removal


def _delete_list(sample_dir: str) -> None:
    """Delete the install list (and sync that), then _commit/."""
    prep = os.path.join(sample_dir, _COMMIT_DIR)
    os.remove(os.path.join(prep, _INSTALL_LIST))
    _sync_dir(prep)  # a power loss can't bring the list back, to be done again after a later commit
    shutil.rmtree(prep)
    _sync_dir(sample_dir)


def _add_file_failures(job, entries) -> None:
    """Add failed-file entries to the job's list, keeping those already there."""
    if entries:
        failures = list(job.context.metadata.get("file_failures") or [])
        job.context.add_metadata("file_failures", failures + list(entries))


def _keep_incomplete_claim(claim_dir: str, counted, not_read, batch_bam_error) -> None:
    """Keep a claim folder whose reads did not all reach target.bam, and say why: they
    count as reads not in target.bam until someone deals with them."""
    if batch_bam_error:
        bam = f"could not be written: {batch_bam_error}"
    else:
        bam = "written" if counted else "not written"
    with open(os.path.join(claim_dir, "incomplete.txt"), "w") as f:
        f.write(
            "The reads of these staged files are not all in target.bam.\n"
            f"Coverage counted for: {' '.join(counted) or 'none'}\n"
            f"Could not be read (coverage not counted): {' '.join(not_read) or 'none'}\n"
            f"Batch BAM: {bam}\n"
        )


def _claims_left(sample_dir: str) -> List[str]:
    """Wait for workers still adding claimed staged files, up to _CLAIM_WAIT_SECONDS in
    all. Returns the claim folders still there: a worker still adding them, or one that
    stopped or failed. Their reads are not in target.bam."""
    deadline = time.time() + _CLAIM_WAIT_SECONDS
    left = []
    for claim_dir in sorted(glob.glob(os.path.join(sample_dir, "_staging", "_batch_*"))):
        lock_path = _claim_lock_path(claim_dir)
        lock = FileLock(lock_path, timeout=max(0.0, deadline - time.time()))
        try:
            lock.acquire()
        except TimeoutError:
            lock.release()
            left.append(os.path.basename(claim_dir))  # still being added
            continue
        lock.release()
        try:
            os.remove(lock_path)
        except FileNotFoundError:
            pass
        if os.path.exists(claim_dir):
            left.append(os.path.basename(claim_dir))
    return left


def _kept_claims(sample_dir: str) -> List[str]:
    """The claim folders kept because their reads are not all in target.bam: each holds
    incomplete.txt, whichever job kept it. Call it under the accumulation lock."""
    return sorted(
        os.path.basename(os.path.dirname(path))
        for path in glob.glob(os.path.join(sample_dir, "_staging", "_batch_*", "incomplete.txt"))
    )


def _incomplete_staged_sets(staging_dir: str) -> List[str]:
    """Staged coverage files whose set misses a part. Staging writes the coverage file
    last, so such a set is broken, not still being written; no claim takes it."""
    complete = set(_complete_staged_counters(staging_dir))
    return [
        os.path.basename(path)
        for path in sorted(glob.glob(os.path.join(staging_dir, "coverage_*.parquet")))
        if os.path.basename(path)[len("coverage_"):-len(".parquet")] not in complete
    ]


def _merge_missed_reads(result: Dict[str, Any]) -> Optional[str]:
    """What kept reads out of target.bam in a finalize result; None when nothing did."""
    status = result.get("status")
    if status not in ("success", "no_files", "no_pending_files"):
        detail = result.get("error")
        return f"the accumulation ended with '{status}'" + (f": {detail}" if detail else "")
    if result.get("batch_bam_error"):
        return f"a batch BAM could not be written: {result['batch_bam_error']}"
    if result.get("claims_left"):
        folders = ", ".join(f"_staging/{name}" for name in result["claims_left"])
        return f"staged reads in {folders} were not added (a worker is still adding them, or it stopped)"
    if result.get("incomplete_staged"):
        return (
            "staged sets in _staging/ miss a part, so their reads were not added: "
            + ", ".join(result["incomplete_staged"])
        )
    return None


# Which target.bam the SNP calls in clair3/ came from, when the call ran.
_SNP_SOURCE = "target_bam_source.json"


def _snp_source_path(sample_dir: str) -> str:
    return os.path.join(sample_dir, "clair3", _SNP_SOURCE)


def _target_bam_identity(target_bam: str) -> Dict[str, Any]:
    """Which target.bam this is: its merged batches, size and modification time. Every
    way ROBIN writes target.bam changes it. A replacement from outside ROBIN that keeps
    the old size and modification time would not be noticed."""
    st = os.stat(target_bam)
    return {
        "batches": sorted(_merged_batch_names(target_bam)),
        "size": st.st_size,
        "mtime_ns": st.st_mtime_ns,
    }


def _reads_wait_to_merge(sample_dir: str) -> bool:
    """True when a batch BAM, a staged file or a claim folder holds reads not yet in
    target.bam."""
    staging = os.path.join(sample_dir, "_staging")
    return bool(
        glob.glob(os.path.join(sample_dir, "batch_*.bam"))
        or glob.glob(os.path.join(staging, "coverage_*.parquet"))
        or glob.glob(os.path.join(staging, "_batch_*"))
    )


def snp_calls_are_current(sample_dir: str) -> bool:
    """True when the SNP calls in clair3/ came from target.bam as it is now, and no read
    waits to be merged into it. Calls made before this record existed count as old."""
    clair = os.path.join(sample_dir, "clair3")
    target_bam = os.path.join(sample_dir, "target.bam")
    vcfs = [os.path.join(clair, n) for n in ("snpsift_output.vcf", "snpsift_indel_output.vcf")]
    if not os.path.exists(target_bam) or not all(os.path.exists(v) for v in vcfs):
        return False
    try:
        with open(_snp_source_path(sample_dir)) as f:
            recorded = json.load(f)
        current = _target_bam_identity(target_bam)
    except Exception:
        return False
    return recorded == current and not _reads_wait_to_merge(sample_dir)


def _write_snp_source(sample_dir: str, identity: Dict[str, Any]) -> None:
    path = _snp_source_path(sample_dir)
    tmp = f"{path}.tmp"
    with open(tmp, "w") as f:
        json.dump(identity, f)
    os.replace(tmp, path)


def _merge_into_target_bam(inputs: List[str], target_bam_path: str, merged_names: List[str]) -> None:
    """Merge ``inputs`` into target.bam, listing ``merged_names`` in its header, and
    replace target.bam (and its index) only once the new file is complete."""
    folder = os.path.dirname(target_bam_path)
    tag = uuid.uuid4().hex[:8]
    tmp_merged = os.path.join(folder, f".target_merge_{tag}.bam")
    tmp_header = os.path.join(folder, f".target_header_{tag}.sam")
    tmp_final = os.path.join(folder, f".target_final_{tag}.bam")
    try:
        pysam.merge("-f", "-o", tmp_merged, *inputs)
        with pysam.AlignmentFile(tmp_merged, "rb") as src:
            header = src.header.to_dict()
        header["CO"] = [c for c in header.get("CO", []) if not c.startswith(_MERGED_TAG)] + [
            f"{_MERGED_TAG}{name}" for name in merged_names
        ]
        with open(tmp_header, "w") as f:
            f.write(str(pysam.AlignmentHeader.from_dict(header)))
        # reheader copies the compressed records as they are: no second decode and encode.
        pysam.reheader("-P", tmp_header, tmp_merged, save_stdout=tmp_final)
        pysam.index(tmp_final)
        os.replace(tmp_final, target_bam_path)
        os.replace(f"{tmp_final}.bai", f"{target_bam_path}.bai")
    finally:
        for path in (tmp_merged, tmp_header, tmp_final, f"{tmp_final}.bai"):
            try:
                os.remove(path)
            except FileNotFoundError:
                pass


def _reindex_if_stale(bam_path: str, logger) -> None:
    """Index ``bam_path`` again when its index is missing or older than it: a finalize
    that replaced target.bam but failed to replace the index leaves the old one, and
    region fetches through it return the wrong records."""
    if not os.path.exists(bam_path):
        return
    bai = f"{bam_path}.bai"
    if os.path.exists(bai) and os.stat(bai).st_mtime_ns >= os.stat(bam_path).st_mtime_ns:
        return
    logger.warning(f"{os.path.basename(bai)} is missing or older than its BAM; indexing again")
    tmp = os.path.join(os.path.dirname(bam_path), f".reindex_{uuid.uuid4().hex[:8]}.bai")
    try:
        pysam.index(bam_path, tmp)
        os.replace(tmp, bai)
    finally:
        try:
            os.remove(tmp)
        except FileNotFoundError:
            pass


def _follows_target_bam(path: str, target_bam_path: str) -> bool:
    """True when ``path``, a BAM sorted from target.bam, holds target.bam's current records:
    its header lists the same merged batches (a sort keeps the @CO lines), and it isn't
    older than target.bam (for folders from before the list existed)."""
    if not os.path.exists(target_bam_path):
        return True
    try:
        return (
            _merged_batch_names(path) == _merged_batch_names(target_bam_path)
            and os.stat(path).st_mtime_ns >= os.stat(target_bam_path).st_mtime_ns
        )
    except Exception:
        return False


def finalize_accumulation_for_sample(
    sample_id: str, work_dir: str, target_panel: str, reference: Optional[str] = None
) -> Dict[str, Any]:
    """
    Force final accumulation of any remaining staged files for a sample
    and merge all batch BAMs into final target.bam.
    
    This should be called when:
    - All files for a sample have been processed
    - The workflow is completing
    - There are staged files that haven't been accumulated yet
    
    Args:
        sample_id: Sample identifier
        work_dir: Working directory containing sample data
        target_panel: Target panel type
    
    Returns:
        Dictionary with accumulation results
    """
    logger = logging.getLogger("robin.target")
    
    try:
        logger.info(f"Finalizing accumulation for sample {sample_id}")
        
        sample_output_dir = os.path.join(work_dir, sample_id)
        
        # Initialize target analysis
        target_analysis = TargetAnalysis(
            work_dir=work_dir,
            target_panel=target_panel,
            batch_size=1,  # Force accumulation regardless of count
            use_staging=True
        )
        if reference:
            target_analysis.reference = reference
            logger.info(f"Finalize path using reference genome: {reference}")
        
        # Check if there are pending files
        pending_count = target_analysis._get_pending_count(sample_id)
        
        if pending_count > 0:
            logger.info(f"Found {pending_count} pending files for {sample_id} - forcing accumulation")
            # Force accumulation of remaining files
            result = target_analysis.accumulate_staged_files(sample_id, force=True)
            logger.info(f"Final accumulation complete for {sample_id}: {result}")
        else:
            logger.info(f"No pending files for {sample_id} - skipping accumulation")
            result = {"status": "no_pending_files", "sample_id": sample_id}

        # Staged files another worker claimed are not in a batch BAM yet: wait for it,
        # so its reads are merged too (JOB-6).
        left = _claims_left(sample_output_dir)
        # A staged set with a missing part is never claimed: list it. Under the
        # accumulation lock, so a claim moving a set's files is never seen half done.
        # First finish a commit a stopped worker left part-way (TGT-3): after the wait, so
        # a worker that died while waited for is finished too, and its claim folder dropped.
        with FileLock(target_analysis._get_lock_file(sample_id, "accumulation"), timeout=60.0):
            finished = target_analysis._finish_commit(sample_output_dir, logger)
            # Every kept claim folder, whichever job kept it: this finish, or another job's
            # after the wait above (it holds incomplete.txt)
            kept = _kept_claims(sample_output_dir)
            broken = _incomplete_staged_sets(target_analysis._get_staging_dir(sample_id))
        if finished:
            result["file_failures"] = list(result.get("file_failures") or []) + finished
        # A claim folder seen in the wait and still there is listed, and so is every kept
        # one: their reads are not in target.bam.
        left = [name for name in left if os.path.isdir(os.path.join(sample_output_dir, "_staging", name))]
        left += [name for name in kept if name not in left]
        if left:
            logger.warning(f"Staged reads in {left} for {sample_id} are not in target.bam")
            result["claims_left"] = left
        if broken:
            logger.warning(f"Staged sets missing a part for {sample_id}: {broken}")
            result["incomplete_staged"] = broken

        # Merge every new batch BAM into target.bam, keeping the reads already in it
        # (TGT-1). target.bam's header lists the batch BAMs merged into it, so one left
        # behind by a crash after the merge is not merged twice. The merge has its own
        # lock, not the accumulation's: a batch BAM appears here only after its commit, by
        # rename, so the glob below is a safe snapshot, and a long merge never makes a
        # worker's commit time out.
        target_bam_path = os.path.join(sample_output_dir, "target.bam")
        lock_file = target_analysis._get_lock_file(sample_id, "finalize")
        with FileLock(lock_file, timeout=300.0):
            _reindex_if_stale(target_bam_path, logger)
            merged_before = _merged_batch_names(target_bam_path)
            batch_bams = sorted(glob.glob(os.path.join(sample_output_dir, "batch_*.bam")))
            new_bams = [b for b in batch_bams if os.path.basename(b) not in merged_before]
            if new_bams:
                logger.info(f"Merging {len(new_bams)} batch BAM files into target.bam for {sample_id}")
                inputs = ([target_bam_path] if os.path.exists(target_bam_path) else []) + new_bams
                names = sorted(merged_before | {os.path.basename(b) for b in new_bams})
                _merge_into_target_bam(inputs, target_bam_path, names)
                result["final_merge"] = "success"
                result["batch_files_merged"] = len(new_bams)
            else:
                logger.info(f"No new batch BAM files for {sample_id}")
                result["final_merge"] = "no_batch_files"
            for bam in batch_bams:  # merged now, or merged before a crash
                for path in (bam, f"{bam}.bai"):
                    try:
                        os.remove(path)
                    except FileNotFoundError:
                        pass
                    except OSError as e:
                        logger.warning(f"Could not remove {os.path.basename(path)}: {e}")
            # IGV and the per-gene variant check read igv_ready.bam: keep it in step.
            if os.path.exists(os.path.join(sample_output_dir, "igv", "igv_ready.bam")):
                if not ensure_sorted_igv_bam(sample_output_dir):
                    logger.warning(f"Could not refresh igv_ready.bam for {sample_id}")
        return result
        
    except Exception as e:
        logger.error(f"Error during final accumulation for {sample_id}: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return {"status": "error", "error": str(e), "sample_id": sample_id}


def target_bam_finalize_handler(job, work_dir: Optional[str] = None) -> None:
    """
    Workflow handler to finalize target BAMs for a sample.
    Runs accumulation and merges batch BAMs into target.bam.

    Required metadata:
    - sample_id: Sample identifier (optional; defaults to context sample ID)
    - work_dir: Parent directory containing sample folder (optional)
    - target_panel: Target panel type (optional; resolved from master.csv if missing)
    """
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)

    try:
        sid = (
            job.context.get_sample_id()
            if hasattr(job.context, "get_sample_id")
            else "unknown"
        )
        metadata = job.context.metadata or {}
        sample_id = metadata.get("sample_id") or sid
        sample_dir = metadata.get("sample_dir")
        base = work_dir or metadata.get("work_dir")

        if not sample_dir:
            if base and sample_id and sample_id != "unknown":
                sample_dir = os.path.join(base, sample_id)
            elif os.path.isdir(job.context.filepath):
                sample_dir = job.context.filepath
            else:
                sample_dir = os.path.dirname(job.context.filepath)

        if not base and sample_dir:
            base = os.path.dirname(sample_dir)

        if not base or not sample_id:
            raise RuntimeError(
                f"Missing work_dir/sample_id for target BAM finalization (work_dir={base}, sample_id={sample_id})"
            )

        target_panel = metadata.get("target_panel")
        reference = metadata.get("reference")
        if not target_panel and sample_dir and os.path.isdir(sample_dir):
            try:
                import csv
                master_csv = os.path.join(sample_dir, "master.csv")
                if os.path.exists(master_csv):
                    with open(master_csv, "r", newline="") as fh:
                        reader = csv.DictReader(fh)
                        first_row = next(reader, None)
                        if first_row:
                            panel = first_row.get("analysis_panel", "").strip()
                            if panel:
                                target_panel = panel
            except Exception:
                pass

        if not target_panel:
            target_panel = "rCNS2"

        logger.info(
            f"Finalizing target BAM: sample_id={sample_id}, work_dir={base}, target_panel={target_panel}, reference={reference}"
        )

        result = finalize_accumulation_for_sample(
            sample_id=sample_id, work_dir=base, target_panel=target_panel, reference=reference
        )

        if result.get("status") == "error":
            job.context.add_error("target_bam_finalize", result.get("error", "Unknown error"))
        else:
            # BAMs whose coverage its accumulation could not add: failed files (batch 10e)
            _add_file_failures(job, result.get("file_failures"))
            job.context.add_result("target_bam_finalize", result)

    except Exception as e:
        job.context.add_error("target_bam_finalize", str(e))
        logger.error(f"Error in target BAM finalization handler: {e}")


def ensure_sorted_igv_bam(
    sample_dir: str, threads: int = 4, force_regenerate: bool = False
) -> str:
    """
    Ensure an IGV-ready coordinate-sorted and indexed BAM exists for a sample.

    Preferred output: ``sample_dir/igv/igv_ready.bam`` (+ .bai).
    Backward compat: if an existing Clair3-sorted BAM is present, reuse it.
    If neither exists, sort and index ``sample_dir/target.bam`` into
    ``sample_dir/igv/igv_ready.bam``. An existing IGV or Clair3 BAM is used only while
    it holds target.bam's current records; after a finalize adds reads it is rebuilt.

    Args:
        sample_dir: Path to the sample output directory.
        threads: Number of threads to pass to pysam.sort.
        force_regenerate: If True, force recreation even if output exists.

    Returns:
        The path to the sorted BAM file (empty string on failure).
    """
    logger = logging.getLogger("robin.target")
    try:
        sample_dir = os.path.abspath(sample_dir)
        igv_dir = os.path.join(sample_dir, "igv")
        clair_dir = os.path.join(sample_dir, "clair3")
        os.makedirs(igv_dir, exist_ok=True)

        # Preferred IGV-ready path
        out_bam = os.path.join(igv_dir, "igv_ready.bam")
        out_bai = f"{out_bam}.bai"
        in_bam = os.path.join(sample_dir, "target.bam")

        def install(bam: str, bai: str) -> None:
            # Readers (IGV, the per-gene check) never see half a file.
            os.replace(bai, out_bai)
            os.replace(bam, out_bam)

        # One builder at a time per sample, so a sort of an older target.bam can't
        # replace one of the current target.bam after finalize checked it.
        with FileLock(os.path.join(sample_dir, "_locks", "igv.lock"), timeout=600.0):
            if (
                not force_regenerate
                and os.path.exists(out_bam)
                and os.path.exists(out_bai)
                and _follows_target_bam(out_bam, in_bam)
            ):
                logger.info(f"IGV BAM already present: {out_bam}")
                return out_bam

            tag = uuid.uuid4().hex[:8]
            tmp_bam = os.path.join(igv_dir, f".igv_ready_{tag}.bam")
            tmp_bai = f"{tmp_bam}.bai"
            try:
                # Backward compatibility: reuse Clair3-sorted BAM if present
                legacy_bams = [
                    os.path.join(clair_dir, "sorted_targets_exceeding_rerun.bam"),
                    os.path.join(clair_dir, "sorted_targets_exceeding.bam"),
                ]
                for legacy in legacy_bams:
                    if (
                        os.path.exists(legacy)
                        and os.path.exists(f"{legacy}.bai")
                        and _follows_target_bam(legacy, in_bam)
                    ):
                        logger.info(f"Reusing legacy Clair3 BAM for IGV: {legacy}")
                        # Optionally copy to IGV location to standardize
                        try:
                            shutil.copy2(legacy, tmp_bam)
                            shutil.copy2(f"{legacy}.bai", tmp_bai)
                            install(tmp_bam, tmp_bai)
                            return out_bam
                        except Exception:
                            # Fallback: just return legacy path
                            return legacy

                if not os.path.exists(in_bam):
                    logger.warning(
                        f"target.bam not found in {sample_dir}; cannot build IGV BAM"
                    )
                    return ""

                logger.info(f"Sorting BAM for IGV: {in_bam} -> {out_bam}")
                pysam.sort(f"-@{threads}", "-o", tmp_bam, in_bam)
                pysam.index(tmp_bam, tmp_bai)
                install(tmp_bam, tmp_bai)
                logger.info(f"IGV BAM ready: {out_bam}")
                return out_bam
            finally:
                for path in (tmp_bam, tmp_bai):
                    try:
                        os.remove(path)
                    except FileNotFoundError:
                        pass
    except Exception as e:
        logger.error(f"Failed to create IGV BAM in {sample_dir}: {e}")
        return ""


def igv_bam_handler(job, work_dir: Optional[str] = None) -> None:
    """
    Workflow handler to build an IGV-ready BAM for a sample.

    Expects a per-sample directory containing `target.bam`. Outputs to
    `sample_dir/igv/igv_ready.bam` (and .bai). On success, records the output
    path in the job context under results['igv_bam'].
    """
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)
    try:
        # Determine sample dir: prefer provided work_dir + sample_id
        sid = (
            job.context.get_sample_id()
            if hasattr(job.context, "get_sample_id")
            else "unknown"
        )
        base = work_dir or job.context.metadata.get("work_dir")
        sample_dir = None

        logger.info(f"IGV BAM job: sample_id={sid}, work_dir={work_dir}, base={base}")

        if base and sid and sid != "unknown":
            sample_dir = os.path.join(base, sid)
            logger.info(f"Using work_dir + sample_id path: {sample_dir}")
        else:
            # Fallback: directory of the current file
            sample_dir = os.path.dirname(job.context.filepath)
            logger.info(f"Using fallback path: {sample_dir}")

        # Additional fallback: if the filepath is already a sample directory, use it directly
        if os.path.isdir(job.context.filepath):
            sample_dir = job.context.filepath
            logger.info(f"Using filepath directly as sample directory: {sample_dir}")

        if not os.path.isdir(sample_dir):
            raise RuntimeError(f"Sample directory not found: {sample_dir}")

        logger.info(f"Final sample directory: {sample_dir}")
        logger.info(
            f"Looking for target.bam in: {os.path.join(sample_dir, 'target.bam')}"
        )

        # Check if force_regenerate is requested
        force_regenerate = job.context.metadata.get("force_regenerate", False)
        if force_regenerate:
            logger.info(
                "Force regenerate requested - will recreate IGV BAM even if it exists"
            )

        out = ensure_sorted_igv_bam(sample_dir, force_regenerate=force_regenerate)
        if not out:
            raise RuntimeError("Failed to create IGV-ready BAM")

        # Record result
        try:
            job.context.add_result("igv_bam", {"bam": out, "bai": f"{out}.bai"})
        except Exception:
            pass
        logger.info(f"IGV-ready BAM created: {out}")
    except Exception as e:
        job.context.add_error("igv_bam", str(e))
        logger.error(f"igv_bam handler failed: {e}")


def run_snp_analysis(
    sample_dir: str,
    threads: int = 4,
    force_regenerate: bool = False,
    reference: Optional[str] = None,
    annotation_only: bool = False,
    annotation_verbose: bool = False,
    target_panel: Optional[str] = None,
) -> str:
    """
    Run SNP analysis using Clair3 for a sample.

    This function runs the complete SNP calling pipeline including:
    - BAM sorting and preparation
    - Clair3 variant calling
    - snpEff annotation
    - SnpSift annotation with ClinVar data
    - Optional annotation-only reruns with verbose logging

    Args:
        sample_dir: Path to the sample output directory
        threads: Number of threads to use for processing
        force_regenerate: If True, force recreation even if output exists
        annotation_only: If True, run only the snpEff/SnpSift annotation steps
        annotation_verbose: If True, run annotation steps with verbose logging enabled
        target_panel: Target panel name used to resolve full target BED file

    Returns:
        The path to the output directory containing SNP results (empty string on failure)
    """
    logger = logging.getLogger("robin.target")

    # CRITICAL: Immediate logging to see if we even get here
    logger.info("ENTERING run_snp_analysis function")
    logger.info(
        "Parameters: sample_dir=%s, threads=%s, force_regenerate=%s, "
        "reference=%s, annotation_only=%s, annotation_verbose=%s, target_panel=%s",
        sample_dir,
        threads,
        force_regenerate,
        reference,
        annotation_only,
        annotation_verbose,
        target_panel,
    )

    # Also log to stderr to make sure we see it

    try:
        logger.info("STEP 1: Setting up sample directory")
        sample_dir = os.path.abspath(sample_dir)
        logger.info(f"Absolute sample_dir: {sample_dir}")

        clair_dir = os.path.join(sample_dir, "clair3")
        logger.info(f"Clair3 output directory: {clair_dir}")

        os.makedirs(clair_dir, exist_ok=True)
        logger.info(f"Clair3 directory created/verified: {clair_dir}")

        logger.info("STEP 2: Checking for existing SNP analysis")

        # Check if SNP analysis already exists and force_regenerate is not set
        snp_output_files = [
            os.path.join(clair_dir, "snpsift_output.vcf"),
            os.path.join(clair_dir, "snpsift_indel_output.vcf"),
        ]

        if annotation_only:
            logger.info("Annotation-only mode enabled; skipping existing output shortcut.")
        # Earlier calls are reused only when they came from target.bam as it is now, with
        # nothing waiting to be merged (JOB-6); otherwise they would hide newer reads.
        elif not force_regenerate and snp_calls_are_current(sample_dir):
            logger.info(f"SNP analysis already present in {clair_dir}")
            logger.info("Rebuilding SNP display JSON from existing snpsift output.")
            try:
                snp_display_vcf = Path(clair_dir) / "snpsift_output.vcf"
                snp_display_path = Path(clair_dir) / "snpsift_output_display.json"
                if snp_display_vcf.exists():
                    snp_display = build_snp_display_data(snp_display_vcf)
                    if snp_display is not None:
                        with snp_display_path.open("w", encoding="utf-8") as f_out:
                            json.dump(snp_display, f_out)
                        logger.info(
                            f"SNP display data refreshed at {snp_display_path}"
                        )
                    else:
                        logger.warning(
                            "Could not regenerate SNP display data from existing VCF."
                        )
                else:
                    logger.warning(
                        "snpsift_output.vcf not found while refreshing display JSON."
                    )
            except Exception as display_exc:
                logger.warning(f"Failed to refresh SNP display JSON: {display_exc}")
            return clair_dir

        if not annotation_only:
            # The calls are about to be redone: drop their record first, so a run that
            # stops part-way never leaves one that matches target.bam.
            try:
                os.remove(_snp_source_path(sample_dir))
            except FileNotFoundError:
                pass

        logger.info("STEP 3: Checking for required input files")
        target_bam = os.path.join(sample_dir, "target.bam")
        threshold_targets_bed = os.path.join(sample_dir, "targets_exceeding_threshold.bed")
        targets_bed = threshold_targets_bed
        output_snv_vcf = os.path.join(clair_dir, "output_done.vcf.gz")
        output_indel_vcf = os.path.join(clair_dir, "output_indel_done.vcf.gz")

        def resolve_full_targets_bed() -> str:
            # 1) Prefer target panel BED from resources when panel is known.
            if target_panel:
                if target_panel == "rCNS2":
                    bed_filename = "rCNS2_panel_name_uniq.bed"
                elif target_panel == "AML":
                    bed_filename = "AML_panel_name_uniq.bed"
                else:
                    bed_filename = f"{target_panel}_panel_name_uniq.bed"

                panel_paths = []
                if resources is not None:
                    try:
                        panel_paths.append(
                            os.path.join(
                                os.path.dirname(os.path.abspath(resources.__file__)),
                                bed_filename,
                            )
                        )
                    except Exception:
                        pass
                panel_paths.extend(
                    [bed_filename, f"data/{bed_filename}", f"/usr/local/share/{bed_filename}"]
                )
                for path in panel_paths:
                    if os.path.exists(path):
                        return path

            # 2) Fall back to sample-specific master BED if available.
            try:
                import glob

                bed_dir = os.path.join(sample_dir, "bed_files")
                if os.path.isdir(bed_dir):
                    master_bed_files = glob.glob(os.path.join(bed_dir, "master_*.bed"))
                    if master_bed_files:
                        return max(master_bed_files, key=os.path.getmtime)
            except Exception as e:
                logger.debug(f"Error resolving master BED from sample bed_files/: {e}")

            # 3) Final fallback: legacy threshold BED.
            return threshold_targets_bed

        logger.info("CHECKING INPUT FILES FOR CLAIR3")
        logger.info(f"  Sample directory: {sample_dir}")
        logger.info(f"  Target BAM path: {target_bam}")
        targets_bed = resolve_full_targets_bed()
        logger.info(f"  Threshold BED path: {threshold_targets_bed}")
        logger.info(f"  Targets BED path selected for ClairS: {targets_bed}")
        logger.info(f"  Reference path: {reference}")

        if annotation_only:
            missing_annotation_inputs = [
                path for path in [output_snv_vcf, output_indel_vcf] if not os.path.exists(path)
            ]
            if missing_annotation_inputs:
                for missing in missing_annotation_inputs:
                    logger.error(f"Annotation input not found: {missing}")
                logger.error(
                    "Annotation-only mode requires existing Clair3 outputs. "
                    "Run the full SNP pipeline first."
                )
                return ""
        else:
            if not os.path.exists(target_bam):
                logger.error(f"target.bam not found: {target_bam}")
                return ""

            if not os.path.exists(targets_bed):
                logger.error(f"Target BED file not found: {targets_bed}")
                return ""

        # Check for reference genome
        if not reference:
            # Try to find reference from common locations
            possible_refs = [
                os.path.join(sample_dir, "reference.fasta"),
                os.path.join(sample_dir, "reference.fa"),
                os.path.join(os.path.dirname(sample_dir), "reference.fasta"),
                os.path.join(os.path.dirname(sample_dir), "reference.fa"),
                "/usr/local/share/reference.fasta",
                "/usr/local/share/reference.fa",
            ]

            for ref_path in possible_refs:
                if os.path.exists(ref_path):
                    reference = ref_path
                    break

        if not reference:
            if annotation_only:
                logger.warning(
                    "Reference genome not found, continuing because annotation-only mode is enabled."
                )
            else:
                logger.error(
                    "Reference genome not found. SNP calling requires a reference genome."
                )
                return ""

        logger.info(f"Starting SNP analysis for sample: {os.path.basename(sample_dir)}")
        logger.info(f"Reference: {reference}")
        logger.info(f"Target BAM: {target_bam}")
        logger.info(f"Targets BED: {targets_bed}")

        if annotation_only:
            logger.info("Annotation-only mode enabled; skipping Clair3 variant calling.")
            sorted_bam = os.path.join(clair_dir, "sorted_targets_exceeding.bam")
        else:
            logger.info("STEP 4: Sorting target BAM for Clair3")
            sorted_bam = os.path.join(clair_dir, "sorted_targets_exceeding.bam")
            logger.info(f"Sorting BAM: {target_bam} -> {sorted_bam}")

            try:
                # Check if sorted BAM already exists and is valid
                if os.path.exists(sorted_bam) and os.path.exists(sorted_bam + ".bai"):
                    logger.info(f"Sorted BAM already exists: {sorted_bam}")
                    # Verify the file is readable
                    try:
                        test_bam = pysam.AlignmentFile(sorted_bam, "rb")
                        test_bam.close()
                        logger.info("Sorted BAM file is valid and readable")
                    except Exception as e:
                        logger.warning(
                            f"Existing sorted BAM appears corrupted, regenerating: {e}"
                        )
                        os.remove(sorted_bam)
                        if os.path.exists(sorted_bam + ".bai"):
                            os.remove(sorted_bam + ".bai")
                        raise Exception("Corrupted BAM file")
                else:
                    logger.info("Creating new sorted BAM file...")

                # Sort the BAM file
                logger.info(f"Sorting BAM with {threads} threads...")
                pysam.sort(f"-@{threads}", "-o", sorted_bam, target_bam)

                # Verify the sorted file was created
                if not os.path.exists(sorted_bam):
                    raise Exception(f"Sorted BAM file was not created: {sorted_bam}")

                # Create index
                logger.info("Creating BAM index...")
                pysam.index(sorted_bam)

                # Verify index was created
                if not os.path.exists(sorted_bam + ".bai"):
                    raise Exception(f"BAM index was not created: {sorted_bam}.bai")

                # Final verification
                file_size = os.path.getsize(sorted_bam)
                logger.info(
                    f"BAM sorting completed successfully: {sorted_bam} ({file_size / (1024**2):.1f} MB)"
                )

            except Exception as e:
                logger.error(f"BAM sorting failed: {e}")
                # Clean up any partial files
                for file_path in [sorted_bam, sorted_bam + ".bai"]:
                    if os.path.exists(file_path):
                        os.remove(file_path)
                        logger.info(f"Cleaned up partial file: {file_path}")
                return ""

            logger.info("STEP 5: Running Clair3 pipeline")
            logger.info("Running Clair3 pipeline...")

            # Run Clair3 pipeline using Docker
            logger.info("CHECKING DOCKER AVAILABILITY")
            if docker is None:
                logger.error("Docker module not available. Cannot run Clair3 pipeline.")
                return ""
            else:
                logger.info("Docker module is available")

            client = docker.from_env()

            # Test Docker connectivity
            try:
                client.ping()
                logger.info("Docker daemon is accessible")
            except Exception as docker_error:
                logger.error(f"Docker daemon not accessible: {docker_error}")
                return ""

            # Define container-specific paths using the original working pattern
            container_bamfile = f"/data/output/{os.path.basename(sorted_bam)}"  # BAM goes to output directory
            container_bedfile = f"/data/bed/{os.path.basename(targets_bed)}"  # BED goes to bed directory
            container_reference = f"/data/reference/{os.path.basename(reference)}"  # Reference goes to reference directory
            container_output = "/data/output"

            # Check if Clair3 image exists
            try:
                client.images.get("hkubal/clairs-to:latest")
                logger.info("Clair3 image found: hkubal/clairs-to:latest")
            except Exception:
                logger.warning(
                    "Clair3 image not found, attempting to pull: hkubal/clairs-to:latest"
                )
                try:
                    client.images.pull("hkubal/clairs-to:latest")
                    logger.info("Clair3 image pulled successfully")
                except Exception as pull_error:
                    logger.error(f"Failed to pull Clair3 image: {pull_error}")
                    return ""

            # Clean up old debugging - no longer needed with the fixed approach

            # Use the original working pattern: bind each directory to its own container path
            volume_bindings = {
                os.path.abspath(os.path.dirname(sorted_bam)): {
                    "bind": "/data/bam",
                    "mode": "ro",
                },  # BAM directory (clair3)
                os.path.abspath(os.path.dirname(targets_bed)): {
                    "bind": "/data/bed",
                    "mode": "ro",
                },  # BED directory (Sample_103)
                os.path.abspath(os.path.dirname(reference)): {
                    "bind": "/data/reference",
                    "mode": "ro",
                },  # Reference directory
                os.path.abspath(clair_dir): {
                    "bind": "/data/output",
                    "mode": "rw",
                },  # Output directory (clair3)
            }

            logger.info("Verifying volume bindings...")
            for source_path, binding_info in volume_bindings.items():
                if not os.path.exists(source_path):
                    logger.error(f"Volume source path does not exist: {source_path}")
                    return ""
                if not os.access(source_path, os.R_OK):
                    logger.error(f"Volume source path not readable: {source_path}")
                    return ""
                logger.info(
                    f"Volume binding verified: {source_path} -> {binding_info['bind']}"
                )

            # Verify the specific files exist in their directories
            logger.info("Verifying input files in volume directories...")
            
            if not os.path.exists(sorted_bam):
                logger.error(f"Sorted BAM file not found in directory: {sorted_bam}")
                return ""
            if not os.path.exists(targets_bed):
                logger.error(
                    f"Targets BED file not found in directory: {targets_bed}"
                )
                return ""
            if not os.path.exists(reference):
                logger.error(f"Reference file not found in directory: {reference}")
                return ""

            logger.info("All input files verified in their directories")

            host_config = client.api.create_host_config(binds=volume_bindings)

            # Function to split BED file into manageable chunks.
            # Chunks can include multiple chromosomes, constrained by covered genomic span.
            def split_bed_into_chunks(bed_file, max_chunk_size=150000000):
                """Split BED entries into chunks up to max_chunk_size covered genomic span."""
                chunks = []
                try:
                    # Read all BED entries and sort them
                    bed_entries = []
                    with open(bed_file, "r") as f:
                        for line_num, line in enumerate(f, 1):
                            raw_line = line.rstrip("\n")
                            line = raw_line.strip()
                            if not line or line.startswith("#"):
                                continue

                            parts = line.split("\t")
                            if len(parts) < 3:
                                logger.warning(
                                    f"Skipping invalid BED line {line_num}: {line}"
                                )
                                continue

                            chrom = parts[0]
                            start = int(parts[1])
                            end = int(parts[2])
                            if end <= start:
                                continue
                            entry_len = end - start
                            bed_entries.append((chrom, start, end, raw_line, entry_len))

                    # Sort by chromosome and start position
                    bed_entries.sort(key=lambda x: (x[0], x[1]))

                    def covered_span(chrom_bounds):
                        # Sum per-chromosome spans so sparse targets cannot inflate a chunk indefinitely.
                        return sum((end - start) for start, end in chrom_bounds.values())

                    current_chunk = []
                    current_bounds = {}

                    for chrom, start, end, raw_line, entry_len in bed_entries:
                        projected_bounds = dict(current_bounds)
                        if chrom in projected_bounds:
                            prev_start, prev_end = projected_bounds[chrom]
                            projected_bounds[chrom] = (min(prev_start, start), max(prev_end, end))
                        else:
                            projected_bounds[chrom] = (start, end)

                        if current_chunk and covered_span(projected_bounds) > max_chunk_size:
                            chunks.append(current_chunk)
                            current_chunk = []
                            current_bounds = {}

                        if chrom in current_bounds:
                            prev_start, prev_end = current_bounds[chrom]
                            current_bounds[chrom] = (min(prev_start, start), max(prev_end, end))
                        else:
                            current_bounds[chrom] = (start, end)
                        current_chunk.append((chrom, start, end, raw_line, entry_len))

                    if current_chunk:
                        chunks.append(current_chunk)

                    logger.info(
                        f"Combined BED entries into {len(chunks)} chunks (max chunk size: {max_chunk_size:,} bases)"
                    )
                    return chunks
                except Exception as e:
                    logger.error(f"Error splitting BED file into chunks: {e}")
                    return []

            # Default to region-split mode to keep INDEL calling memory usage bounded.
            # Set ROBIN_CLAIRS_SPLIT_REGIONS=0 to force single-pass mode.
            split_regions_env = os.environ.get("ROBIN_CLAIRS_SPLIT_REGIONS")
            if split_regions_env is None:
                use_split_regions = True
            else:
                use_split_regions = split_regions_env.lower() in ("1", "true", "yes", "on")

            if use_split_regions:
                logger.info(
                    "Using region-split ClairS mode (recommended for INDEL memory stability)."
                )
                chunked_entries = split_bed_into_chunks(targets_bed, max_chunk_size=250000000)
                if not chunked_entries:
                    logger.error("Failed to split BED file into chunks")
                    return ""
                chunk_bed_dir = os.path.join(clair_dir, "chunk_beds")
                os.makedirs(chunk_bed_dir, exist_ok=True)

                regions = []
                for i, chunk in enumerate(chunked_entries, 1):
                    chunk_host_bed = os.path.join(chunk_bed_dir, f"region_{i}.bed")
                    with open(chunk_host_bed, "w") as f_out:
                        for _, _, _, raw_line, _ in chunk:
                            f_out.write(raw_line + "\n")

                    first_chrom, first_start = chunk[0][0], chunk[0][1]
                    last_chrom, last_end = chunk[-1][0], chunk[-1][2]
                    chrom_bounds = {}
                    for chrom, start, end, _, _ in chunk:
                        if chrom in chrom_bounds:
                            prev_start, prev_end = chrom_bounds[chrom]
                            chrom_bounds[chrom] = (min(prev_start, start), max(prev_end, end))
                        else:
                            chrom_bounds[chrom] = (start, end)
                    span_bases = sum((end - start) for start, end in chrom_bounds.values())
                    label = f"{first_chrom}:{first_start+1}-{last_end}"
                    if first_chrom != last_chrom:
                        label = f"{first_chrom}:{first_start+1}..{last_chrom}:{last_end}"

                    regions.append(
                        {
                            "label": label,
                            "bed_container": f"{container_output}/chunk_beds/region_{i}.bed",
                            "is_split": True,
                            "span_bases": span_bases,
                        }
                    )
                    logger.info(
                        f"Chunk {i}: {label} (entries: {len(chunk)}, span: {span_bases:,} bases)"
                    )
            else:
                logger.info(
                    "Running ClairS in single-pass mode over full BED."
                )
                logger.info(
                    "Set ROBIN_CLAIRS_SPLIT_REGIONS=1 to re-enable region-splitting mode."
                )
                regions = [
                    {
                        "label": "full_target_bed",
                        "bed_container": container_bedfile,
                        "is_split": False,
                        "span_bases": None,
                    }
                ]

            # Process each region separately (or single pass)
            all_snv_vcfs = []
            all_indel_vcfs = []

            logger.info(f"Container input BAM: {container_bamfile}")
            logger.info(f"Container input BED: {container_bedfile}")
            logger.info(f"Container output dir: {container_output}")
            logger.info(f"Container reference: {container_reference}")
            logger.info(f"Volume bindings: {volume_bindings}")

            # Function to merge VCF files
            def merge_vcf_files(vcf_files, output_file, variant_type):
                """Merge multiple VCF files into a single output file"""
                try:
                    if not vcf_files:
                        logger.error(f"No {variant_type} VCF files to merge")
                        return False

                    logger.info(
                        f"Merging {len(vcf_files)} {variant_type} VCF files into {output_file}"
                    )

                    # Read headers from first file
                    import gzip

                    headers = []
                    first_vcf = vcf_files[0]

                    with gzip.open(first_vcf, "rt") as f:
                        for line in f:
                            if line.startswith("#"):
                                headers.append(line)
                            else:
                                break

                    # Write merged output
                    with gzip.open(output_file, "wt") as out_f:
                        # Write headers
                        for header in headers:
                            out_f.write(header)

                        # Write variants from each file
                        for vcf_file in vcf_files:
                            if os.path.exists(vcf_file):
                                with gzip.open(vcf_file, "rt") as in_f:
                                    for line in in_f:
                                        if not line.startswith("#"):
                                            out_f.write(line)

                    logger.info(
                        f"Successfully merged {variant_type} VCF files into {output_file}"
                    )
                    return True

                except Exception as e:
                    logger.error(f"Error merging {variant_type} VCF files: {e}")
                    return False

            # Log system resources before starting
            try:
                import psutil

                memory = psutil.virtual_memory()
                logger.info(
                    f"System memory: {memory.available / (1024**3):.1f} GB available out of {memory.total / (1024**3):.1f} GB"
                )
                logger.info(f"System CPU cores: {psutil.cpu_count()}")
            except ImportError:
                logger.info("psutil not available, skipping resource logging")

            # Each run writes Clair3's output into its own new folder, so a VCF there can
            # only be this run's: files from an earlier run are never taken as new. Earlier
            # run folders are deleted when possible (the container writes as root).
            for old_run in glob.glob(os.path.join(clair_dir, "run_*")):
                shutil.rmtree(old_run, ignore_errors=True)
            run_name = f"run_{uuid.uuid4().hex[:8]}"
            run_dir = os.path.join(clair_dir, run_name)
            os.makedirs(run_dir)

            # Process each region
            for i, region_spec in enumerate(regions):
                region_label = region_spec["label"]
                print(
                    f"Processing region {i+1}/{len(regions)}: {region_label} - PRINT STATEMENT"
                )
                logger.info(f"Processing region {i+1}/{len(regions)}: {region_label}")

                # Create region-specific output directory
                if region_spec["is_split"]:
                    region_output = f"{container_output}/{run_name}/region_{i+1}"
                else:
                    region_output = f"{container_output}/{run_name}"

                # Build Clair3 command for this region
                command = (
                    f"/opt/bin/run_clairs_to "
                    f"--tumor_bam_fn {container_bamfile} "
                    f"--ref_fn {container_reference} "
                    f"--threads {threads} "
                    f"--remove_intermediate_dir "
                    f"--platform ont_r10_guppy_hac_5khz "
                    f"--output_dir {region_output} "
                    f"-b {region_spec['bed_container']} "
                    f"--chunk_size 5000000 "
                    f"--disable_intermediate_phasing "
                )

                logger.info(f"Running Clair3 command for region {i+1}: {command}")

                # Create and start container for this region
                logger.info(f"Creating Clair3 container for region {i+1}...")
                container = client.api.create_container(
                    image="hkubal/clairs-to:latest",
                    command=command,
                    volumes=list(volume_bindings.keys()),
                    host_config=host_config,
                )
                logger.info(f"Container created with ID: {container.get('Id')}")

                # Start container
                client.api.start(container=container.get("Id"))

                # Stream container logs
                for line in client.api.logs(
                    container=container.get("Id"), stream=True, follow=True
                ):
                    logger.info(f"Clair3 Region {i+1}: {line.decode().strip()}")

                # Wait for completion
                result = client.api.wait(container=container.get("Id"))
                if result["StatusCode"] != 0:
                    status_code = result.get("StatusCode")
                    container_id = container.get("Id")

                    # Classify abnormal termination (OOM/SIGKILL/SIGTERM) for clearer diagnostics.
                    oom_killed = False
                    kill_reason = None
                    try:
                        inspect_data = client.api.inspect_container(container=container_id)
                        state = inspect_data.get("State", {}) if inspect_data else {}
                        oom_killed = bool(state.get("OOMKilled", False))
                    except Exception as inspect_exc:
                        logger.debug(
                            f"Could not inspect Clair3 container state for region {i+1}: {inspect_exc}"
                        )

                    if oom_killed:
                        kill_reason = "oom_killed"
                    elif status_code == 137:
                        kill_reason = "sigkill_or_oom"
                    elif status_code == 143:
                        kill_reason = "sigterm"

                    # Get detailed error logs before removing container
                    error_logs = client.api.logs(
                        container=container_id, stderr=True, stdout=False
                    )
                    error_details = error_logs.decode().strip()
                    if error_details:
                        logger.error(f"Clair3 region {i+1} error logs: {error_details}")

                    # Also get stdout logs for context
                    stdout_logs = client.api.logs(
                        container=container_id, stdout=True, stderr=False
                    )
                    stdout_details = stdout_logs.decode().strip()
                    if stdout_details:
                        logger.info(
                            f"Clair3 region {i+1} stdout logs: {stdout_details}"
                        )

                    if kill_reason is not None:
                        logger.error(
                            "Clair3 region %s terminated abnormally (%s, exit=%s). "
                            "This is often memory pressure (especially for oom_killed/sigkill_or_oom). "
                            "Consider lowering threads or reducing chunk span.",
                            i + 1,
                            kill_reason,
                            status_code,
                        )

                    # Calls without this region would look complete: the run fails.
                    logger.error(
                        f"Clair3 region {i+1} ({region_label}) failed with status code "
                        f"{status_code}; SNP calling failed"
                    )
                    return ""

                # Clean up container
                client.api.remove_container(container=container.get("Id"))

                # Check if output files were created for this region/pass
                if region_spec["is_split"]:
                    region_snv = f"{run_dir}/region_{i+1}/snv.vcf.gz"
                    region_indel = f"{run_dir}/region_{i+1}/indel.vcf.gz"
                else:
                    region_snv = f"{run_dir}/snv.vcf.gz"
                    region_indel = f"{run_dir}/indel.vcf.gz"

                for kind, path, found in (
                    ("SNV", region_snv, all_snv_vcfs),
                    ("INDEL", region_indel, all_indel_vcfs),
                ):
                    if not os.path.exists(path):
                        logger.error(
                            f"Region {i+1} ({region_label}) wrote no {kind} output: {path}; "
                            "SNP calling failed"
                        )
                        return ""
                    found.append(path)
                    logger.info(f"Region {i+1} {kind} output created: {path}")

            if use_split_regions:
                # Merge all split-region outputs
                logger.info(
                    f"Merging {len(all_snv_vcfs)} SNV VCF files and {len(all_indel_vcfs)} INDEL VCF files..."
                )

                if all_snv_vcfs:
                    merged_snv = f"{run_dir}/merged_snv.vcf.gz"
                    if not merge_vcf_files(all_snv_vcfs, merged_snv, "SNV"):
                        return ""
                    shutil.copy2(merged_snv, f"{clair_dir}/output_done.vcf.gz")
                    logger.info(f"Merged SNV VCFs into: {clair_dir}/output_done.vcf.gz")
                else:
                    logger.error("No SNV VCF files were successfully created")
                    return ""

                if all_indel_vcfs:
                    merged_indel = f"{run_dir}/merged_indel.vcf.gz"
                    if not merge_vcf_files(all_indel_vcfs, merged_indel, "INDEL"):
                        return ""
                    shutil.copy2(merged_indel, f"{clair_dir}/output_indel_done.vcf.gz")
                    logger.info(
                        f"Merged INDEL VCFs into: {clair_dir}/output_indel_done.vcf.gz"
                    )
                else:
                    logger.warning("No INDEL VCF files were successfully created")
                logger.info("Clair3 pipeline completed successfully for all regions")
            else:
                # Single-pass output normalization (no merge step needed).
                single_snv = f"{run_dir}/snv.vcf.gz"
                single_indel = f"{run_dir}/indel.vcf.gz"
                if not os.path.exists(single_snv):
                    logger.error(f"Single-pass SNV output not found: {single_snv}")
                    return ""
                shutil.copy2(single_snv, f"{clair_dir}/output_done.vcf.gz")
                logger.info(f"Single-pass SNV output copied to: {clair_dir}/output_done.vcf.gz")

                if os.path.exists(single_indel):
                    shutil.copy2(single_indel, f"{clair_dir}/output_indel_done.vcf.gz")
                    logger.info(
                        f"Single-pass INDEL output copied to: {clair_dir}/output_indel_done.vcf.gz"
                    )
                else:
                    logger.warning(f"Single-pass INDEL output not found: {single_indel}")
                logger.info("Clair3 pipeline completed successfully in single-pass mode")


        if annotation_only:
            logger.info(
                "Proceeding directly to annotation pipeline using existing Clair3 outputs."
            )

        logger.info("STEP 6: Running annotation pipeline")
        logger.info("Running annotation pipeline...")

        try:
            # SNP annotation with snpEff
            logger.info("Starting snpEff annotation for SNPs...")
            logger.info(f"Input file: {output_snv_vcf}")
            logger.info(f"Output file: {clair_dir}/snpeff_output.vcf")

            # Define output file paths
            snpeff_out = f"{clair_dir}/snpeff_output.vcf"
            snpsift_out = f"{clair_dir}/snpsift_output.vcf"
            snpeff_indel_out = f"{clair_dir}/snpeff_indel_output.vcf"
            snpsift_indel_out = f"{clair_dir}/snpsift_indel_output.vcf"

            # Check if input file exists and has content
            if os.path.exists(output_snv_vcf):
                file_size = os.path.getsize(output_snv_vcf)
                logger.info(f"Input VCF file exists, size: {file_size} bytes")

                # Check first few lines to verify it's a valid VCF
                try:
                    import gzip

                    with gzip.open(output_snv_vcf, "rt") as f:
                        first_lines = [next(f) for _ in range(5)]
                    logger.info("First 5 lines of input VCF:")
                    for i, line in enumerate(first_lines):
                        logger.info(f"  Line {i+1}: {line.strip()}")
                except Exception as e:
                    logger.warning(f"Could not read input VCF file: {e}")
            else:
                logger.error(
                    f"Input VCF file does not exist: {output_snv_vcf}"
                )

            snpeff_cmd = ["snpEff"]
            snpeff_cmd.append("-v" if annotation_verbose else "-q")
            snpeff_cmd.extend(["hg38", output_snv_vcf])
            logger.info(f"Running snpEff command: {' '.join(snpeff_cmd)}")

            annotation_env = os.environ.copy()
            annotation_env["_JAVA_OPTIONS"] = "-Xmx16g"

            # Check if snpEff is available
            try:
                snpeff_version = subprocess.run(
                    ["snpEff", "-version"], capture_output=True, text=True
                )
                if snpeff_version.returncode == 0:
                    logger.info(f"snpEff version: {snpeff_version.stdout.strip()}")
                else:
                    logger.warning(
                        f"Could not get snpEff version: {snpeff_version.stderr}"
                    )
            except FileNotFoundError:
                logger.error("snpEff command not found in PATH")

            with open(snpeff_out, "w") as fout:
                logger.info("Starting snpEff execution...")
                result = subprocess.run(
                    snpeff_cmd,
                    stdout=fout,
                    stderr=subprocess.PIPE,
                    text=True,
                    env=annotation_env,
                )
                logger.info(
                    f"snpEff execution completed with return code: {result.returncode}"
                )
                if annotation_verbose and result.stderr:
                    logger.info(f"snpEff stderr output:\n{result.stderr}")

            if result.returncode != 0:
                logger.warning(f"snpEff failed with return code: {result.returncode}")
                logger.warning(f"snpEff stderr output: {result.stderr}")
                logger.info("Using raw Clair3 output for SNP annotation")
                # Decompress the gzipped file when copying to .vcf extension
                import gzip

                with gzip.open(output_snv_vcf, "rt") as f_in:
                    with open(snpeff_out, "w") as f_out:
                        f_out.write(f_in.read())
                logger.info("Created uncompressed VCF from gzipped input")
            else:
                logger.info("snpEff completed successfully!")
                # Check the output file
                if os.path.exists(snpeff_out):
                    output_size = os.path.getsize(snpeff_out)
                    logger.info(
                        f"snpEff output file created, size: {output_size} bytes"
                    )

                    # Check first few lines of output to verify annotation
                    try:
                        with open(snpeff_out, "r") as f:
                            first_lines = [next(f) for _ in range(5)]
                        logger.info("First 5 lines of snpEff output:")
                        for i, line in enumerate(first_lines):
                            logger.info(f"  Line {i+1}: {line.strip()}")
                    except Exception as e:
                        logger.warning(f"Could not read snpEff output: {e}")
                else:
                    logger.error("snpEff output file was not created!")

            # SnpSift annotation
            logger.info("Starting SnpSift annotation for SNPs...")
            logger.info(f"Input file: {snpeff_out}")
            logger.info(f"Output file: {snpsift_out}")

            # Check if snpEff output exists and has content
            if os.path.exists(snpeff_out):
                snpeff_size = os.path.getsize(snpeff_out)
                logger.info(f"snpEff output file exists, size: {snpeff_size} bytes")
            else:
                logger.error(f"snpEff output file does not exist: {snpeff_out}")

            clinvar_db_path = _resolve_clinvar_db_for_snpsift(logger)
            snv_clinvar_status = annotate_with_clinvar(
                snpeff_out,
                snpsift_out,
                clinvar_db_path,
                logger,
                verbose=annotation_verbose,
                env=annotation_env.copy(),
                label="SNP",
            )
            logger.info(f"ClinVar status for SNPs: {snv_clinvar_status}")

            # INDEL annotation
            logger.info("Starting snpEff annotation for INDELs...")
            logger.info(f"Input file: {output_indel_vcf}")
            logger.info(f"Output file: {snpeff_indel_out}")

            # Check if input file exists and has content
            if os.path.exists(output_indel_vcf):
                file_size = os.path.getsize(output_indel_vcf)
                logger.info(f"INDEL input VCF file exists, size: {file_size} bytes")

                # Check first few lines to verify it's a valid VCF
                try:
                    import gzip

                    with gzip.open(output_indel_vcf, "rt") as f:
                        first_lines = [next(f) for _ in range(5)]
                    logger.info("First 5 lines of INDEL input VCF:")
                    for i, line in enumerate(first_lines):
                        logger.info(f"  Line {i+1}: {line.strip()}")
                except Exception as e:
                    logger.warning(f"Could not read INDEL input VCF file: {e}")
            else:
                logger.error(
                    f"INDEL input VCF file does not exist: {output_indel_vcf}"
                )

            snpeff_indel_cmd = ["snpEff"]
            snpeff_indel_cmd.append("-v" if annotation_verbose else "-q")
            snpeff_indel_cmd.extend(["hg38", output_indel_vcf])
            logger.info(f"Running snpEff INDEL command: {' '.join(snpeff_indel_cmd)}")

            with open(snpeff_indel_out, "w") as fout:
                logger.info("Starting snpEff INDEL execution...")
                result = subprocess.run(
                    snpeff_indel_cmd,
                    stdout=fout,
                    stderr=subprocess.PIPE,
                    text=True,
                    env=annotation_env,
                )
                logger.info(
                    f"snpEff INDEL execution completed with return code: {result.returncode}"
                )
                if annotation_verbose and result.stderr:
                    logger.info(f"snpEff INDEL stderr output:\n{result.stderr}")

            if result.returncode != 0:
                logger.warning(
                    f"snpEff (INDEL) failed with return code: {result.returncode}"
                )
                logger.warning(f"snpEff (INDEL) stderr output: {result.stderr}")
                logger.info("Using raw Clair3 INDEL output for annotation")
                # Decompress the gzipped file when copying to .vcf extension
                import gzip

                with gzip.open(output_indel_vcf, "rt") as f_in:
                    with open(snpeff_indel_out, "w") as f_out:
                        f_out.write(f_in.read())
                logger.info("Created uncompressed INDEL VCF from gzipped input")
            else:
                logger.info("snpEff INDEL completed successfully!")
                # Check the output file
                if os.path.exists(snpeff_indel_out):
                    output_size = os.path.getsize(snpeff_indel_out)
                    logger.info(
                        f"snpEff INDEL output file created, size: {output_size} bytes"
                    )

                    # Check first few lines of output to verify annotation
                    try:
                        with open(snpeff_indel_out, "r") as f:
                            first_lines = [next(f) for _ in range(5)]
                        logger.info("First 5 lines of snpEff INDEL output:")
                        for i, line in enumerate(first_lines):
                            logger.info(f"  Line {i+1}: {line.strip()}")
                    except Exception as e:
                        logger.warning(f"Could not read snpEff INDEL output: {e}")
                else:
                    logger.error("snpEff INDEL output file was not created!")

            # SnpSift annotation for INDELs
            logger.info("Starting SnpSift annotation for INDELs...")
            logger.info(f"Input file: {snpeff_indel_out}")
            logger.info(f"Output file: {snpsift_indel_out}")

            # Check if snpEff INDEL output exists and has content
            if os.path.exists(snpeff_indel_out):
                snpeff_indel_size = os.path.getsize(snpeff_indel_out)
                logger.info(
                    f"snpEff INDEL output file exists, size: {snpeff_indel_size} bytes"
                )
            else:
                logger.error(
                    f"snpEff INDEL output file does not exist: {snpeff_indel_out}"
                )

            indel_clinvar_status = annotate_with_clinvar(
                snpeff_indel_out,
                snpsift_indel_out,
                clinvar_db_path,
                logger,
                verbose=annotation_verbose,
                env=annotation_env.copy(),
                label="INDEL",
            )
            logger.info(f"ClinVar status for INDELs: {indel_clinvar_status}")

            logger.info("Annotation pipeline completed successfully")

        except Exception as e:
            logger.error(f"Error in annotation pipeline: {e}")
            return ""

        # Parse VCF files to CSV format
        logger.info("Parsing VCF files...")
        try:
            # Import parse_vcf function or implement basic parsing
            def basic_vcf_to_csv(vcf_file, csv_file):
                """Basic VCF to CSV conversion"""
                try:
                    with open(vcf_file, "r") as f:
                        lines = [line.strip() for line in f if not line.startswith("#")]

                    if not lines:
                        return

                    # Parse header from first non-comment line
                    header = lines[0].split("\t")

                    # Parse data
                    data = []
                    for line in lines[1:]:
                        fields = line.split("\t")
                        if len(fields) >= len(header):
                            data.append(fields[: len(header)])

                    # Create DataFrame and save
                    df = pd.DataFrame(data, columns=header)
                    df.to_csv(csv_file, index=False)
                    logger.info(f"VCF parsed to CSV: {csv_file}")

                except Exception as e:
                    logger.warning(f"Basic VCF parsing failed for {vcf_file}: {e}")

            # Parse both SNP and INDEL VCFs
            basic_vcf_to_csv(
                f"{clair_dir}/snpsift_output.vcf", f"{clair_dir}/snpsift_output.vcf.csv"
            )
            basic_vcf_to_csv(
                f"{clair_dir}/snpsift_indel_output.vcf",
                f"{clair_dir}/snpsift_indel_output.vcf.csv",
            )

            # Build pre-formatted SNP display data for the GUI
            try:
                snp_display_path = Path(clair_dir) / "snpsift_output_display.json"
                snp_display = build_snp_display_data(Path(clair_dir) / "snpsift_output.vcf")
                if snp_display is not None:
                    with snp_display_path.open("w", encoding="utf-8") as f_out:
                        json.dump(snp_display, f_out)
                    logger.info(f"SNP display data written to {snp_display_path}")
                else:
                    logger.warning("Could not generate SNP display data from VCF.")
            except Exception as display_exc:
                logger.warning(f"Failed to build SNP display data: {display_exc}")

            # Note: VCF files are now properly uncompressed when they have .vcf extensions
            # Compressed .vcf.gz files are only used as input to processing tools

        except Exception as e:
            logger.warning(f"VCF parsing failed: {e}")

        logger.info("FINAL SUCCESS: SNP analysis completed successfully!")
        logger.info(
            f"SNP analysis completed successfully for {os.path.basename(sample_dir)}"
        )
        logger.info(f"Output directory: {clair_dir}")
        return clair_dir

    except Exception as e:
        print(f"EXCEPTION in run_snp_analysis: {e} - PRINT STATEMENT")
        logger.error(f"Failed to run SNP analysis in {sample_dir}: {e}")
        import traceback

        print(f"TRACEBACK: {traceback.format_exc()} - PRINT STATEMENT")
        return ""


def snp_analysis_handler(job, work_dir: Optional[str] = None) -> None:
    """
    Workflow handler to run SNP analysis using Clair3 for a sample.

    Expects a per-sample directory containing `target.bam` and `targets_exceeding_threshold.bed`.
    Outputs to `sample_dir/clair3/` with annotated VCF files and CSV summaries.
    On success, records the output directory in the job context under results['snp_analysis'].

    Required metadata:
    - reference: Path to reference genome (optional, will auto-detect)
    - threads: Number of threads to use (default: 4)
    - force_regenerate: Whether to force regeneration of existing results (default: False)
    - annotation_only: Run only the annotation steps (default: False)
    - annotation_verbose: Enable verbose logging for snpEff/SnpSift (default: False)
    """
    logger = get_job_logger(str(job.job_id), job.job_type, job.context.filepath)

    try:
        # Determine sample dir: prefer provided work_dir + sample_id
        sid = (
            job.context.get_sample_id()
            if hasattr(job.context, "get_sample_id")
            else "unknown"
        )
        base = work_dir or job.context.metadata.get("work_dir")
        sample_dir = None

        logger.info(
            f"SNP analysis job: sample_id={sid}, work_dir={work_dir}, base={base}"
        )

        if base and sid and sid != "unknown":
            sample_dir = os.path.join(base, sid)
            logger.info(f"Using work_dir + sample_id path: {sample_dir}")
        else:
            # Fallback: directory of the current file
            sample_dir = os.path.dirname(job.context.filepath)
            logger.info(f"Using fallback path: {sample_dir}")

        # Additional fallback: if the filepath is already a sample directory, use it directly
        if os.path.isdir(job.context.filepath):
            sample_dir = job.context.filepath
            logger.info(f"Using filepath directly as sample directory: {sample_dir}")

        if not os.path.isdir(sample_dir):
            raise RuntimeError(f"Sample directory not found: {sample_dir}")

        logger.info(f"Final sample directory: {sample_dir}")

        # Let's see what's actually in this directory
        try:
            if os.path.exists(sample_dir):
                dir_contents = os.listdir(sample_dir)
                logger.info(f"Directory contents of {sample_dir}:")
                for item in dir_contents:
                    item_path = os.path.join(sample_dir, item)
                    if os.path.isfile(item_path):
                        size = os.path.getsize(item_path)
                        logger.info(f"  {item} ({size / (1024**2):.1f} MB)")
                    else:
                        logger.info(f"  {item}/")
            else:
                logger.error(f"Sample directory does not exist: {sample_dir}")
        except Exception as e:
            logger.warning(f"Could not list directory contents: {e}")

        target_bam = os.path.join(sample_dir, "target.bam")
        targets_bed = os.path.join(sample_dir, "targets_exceeding_threshold.bed")

        # Get job parameters
        threads = job.context.metadata.get("threads", 4)
        force_regenerate = job.context.metadata.get("force_regenerate", False)
        reference = job.context.metadata.get("reference")
        annotation_only = job.context.metadata.get("annotation_only", False)
        annotation_verbose = job.context.metadata.get("annotation_verbose", False)
        target_panel = job.context.metadata.get("target_panel")
        if not target_panel:
            try:
                master_csv = os.path.join(sample_dir, "master.csv")
                if os.path.exists(master_csv):
                    import csv

                    with open(master_csv, "r", newline="") as fh:
                        reader = csv.DictReader(fh)
                        first_row = next(reader, None)
                        if first_row:
                            panel = first_row.get("analysis_panel", "").strip()
                            if panel:
                                target_panel = panel
            except Exception:
                pass

        # Merge what is not yet in target.bam first: new batch BAMs and staged files
        # (JOB-6). Calls from a target.bam missing reads would look complete.
        sample_path = os.path.normpath(sample_dir)
        merge = finalize_accumulation_for_sample(
            sample_id=os.path.basename(sample_path),
            work_dir=os.path.dirname(sample_path),
            target_panel=target_panel or "rCNS2",  # as the finalize job does
            reference=reference,
        )
        missed = _merge_missed_reads(merge)
        if missed:
            raise RuntimeError(f"Not every read is in target.bam, so no variants were called: {missed}")

        # Check for required input files
        if not os.path.exists(target_bam):
            raise RuntimeError(f"target.bam not found: {target_bam}")

        if not os.path.exists(targets_bed):
            raise RuntimeError(
                f"targets_exceeding_threshold.bed not found: {targets_bed}"
            )
        read_from = _target_bam_identity(target_bam)

        # Debug logging for reference genome
        logger.info("=== SNP ANALYSIS HANDLER DEBUGGING ===")
        logger.info(f"Job metadata keys: {list(job.context.metadata.keys())}")
        logger.info(f"Job metadata: {job.context.metadata}")
        logger.info(f"Extracted reference: {reference}")
        logger.info(f"Annotation-only: {annotation_only}")
        logger.info(f"Annotation-verbose: {annotation_verbose}")
        logger.info(f"Target panel: {target_panel}")
        logger.info("=== END SNP ANALYSIS HANDLER DEBUGGING ===")

        if force_regenerate:
            logger.info(
                "Force regenerate requested - will recreate SNP analysis even if it exists"
            )

        if reference:
            logger.info(f"Using reference genome from job metadata: {reference}")
        else:
            logger.warning(
                "No reference genome in job metadata - will try to auto-detect"
            )

        # Run SNP analysis
        logger.info("ABOUT TO CALL run_snp_analysis function")
        logger.info(
            "Calling run_snp_analysis with: sample_dir=%s, threads=%s, "
            "force_regenerate=%s, reference=%s, annotation_only=%s, annotation_verbose=%s, target_panel=%s"
            % (
                sample_dir,
                threads,
                force_regenerate,
                reference,
                annotation_only,
                annotation_verbose,
                target_panel,
            )
        )

        output_dir = run_snp_analysis(
            sample_dir,
            threads=threads,
            force_regenerate=force_regenerate,
            reference=reference,
            annotation_only=annotation_only,
            annotation_verbose=annotation_verbose,
            target_panel=target_panel,
        )

        logger.info(f"run_snp_analysis returned: {output_dir}")

        if not output_dir:
            raise RuntimeError("Failed to complete SNP analysis")
        if not annotation_only:
            # Which reads these calls came from; earlier calls are reused only on a match.
            _write_snp_source(sample_dir, read_from)

        # Record result
        try:
            job.context.add_result(
                "snp_analysis",
                {
                    "output_dir": output_dir,
                    "snps_vcf": os.path.join(output_dir, "snpsift_output.vcf"),
                    "indels_vcf": os.path.join(output_dir, "snpsift_indel_output.vcf"),
                    "snps_csv": os.path.join(output_dir, "snpsift_output.vcf.csv"),
                    "indels_csv": os.path.join(
                        output_dir, "snpsift_indel_output.vcf.csv"
                    ),
                    "sample_id": sid,
                    "analysis_timestamp": time.time(),
                },
            )
        except Exception:
            pass

        logger.info(f"SNP analysis completed successfully for {sid}")
        logger.info(f"Output directory: {output_dir}")

    except Exception as e:
        job.context.add_error("snp_analysis", str(e))
        logger.error(f"snp_analysis handler failed: {e}")
