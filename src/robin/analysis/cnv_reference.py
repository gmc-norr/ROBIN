"""The CNV control profile, built from the stored control counts alone.

Never pass the control to ``cnv_from_bam.iterate_bam_file(copy_numbers=...)``:
that call adds the BAM's reads into the dict it is given, and it sums each bin
chunk as a 16-bit number, which wraps when a chunk holds more than 65,535
reads (a large control at the wide bins of an early run). This module applies
the library's ``calculate_cnv`` formula to the control in 64-bit instead.
"""

import math
from typing import Dict, List, Mapping, Sequence

import natsort
import numpy as np


def reference_cnv_profile(
    control_counts: Mapping[str, Sequence[int]],
    bin_width: int,
    bin_size: int = 1000,
    ploidy: int = 2,
) -> Dict[str, List[float]]:
    """Return the control's CNV profile at ``bin_width``, without changing it.

    ``control_counts`` holds read starts per ``bin_size`` bin for each contig.
    Like ``cnv_from_bam.calculate_cnv``: the bin width is rounded to a whole
    number of bins, each contig's bins are summed in chunks of that many (the
    last chunk may be partial), and each sum is divided by the median of all
    chunk sums, rounded half away from zero (0 becomes 1), times ``ploidy``.

    One difference: the library takes its median over chunks of all contigs
    laid end to end, in hash-map order, so a chunk can straddle two contigs.
    Here every chunk stays in its contig, as the returned values do.
    """
    rounded = (bin_width + bin_size // 2) // bin_size * bin_size
    chunk = max(1, rounded // bin_size)

    sums: Dict[str, np.ndarray] = {}
    for name in natsort.natsorted(control_counts):
        counts = np.asarray(control_counts[name], dtype=np.int64)
        n_chunks = -(-len(counts) // chunk)
        padded = np.zeros(n_chunks * chunk, dtype=np.int64)
        padded[: len(counts)] = counts
        sums[name] = padded.reshape(n_chunks, chunk).sum(axis=1)

    if not sums:
        return {}
    median = math.floor(float(np.median(np.concatenate(list(sums.values())))) + 0.5)
    median = median or 1
    return {name: (s / median * ploidy).tolist() for name, s in sums.items()}
