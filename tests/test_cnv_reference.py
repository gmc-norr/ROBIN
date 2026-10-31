"""CNV-1: the control profile is built from the control's counts alone.

``cnv_from_bam.iterate_bam_file(bam, copy_numbers=control)`` adds the BAM's
reads into ``control`` and sums each bin chunk as a 16-bit number, so ROBIN's
"control" grew with every patient and wrapped early in a run.
``reference_cnv_profile`` re-bins the stored control counts with the
library's formula (chunk sums / rounded median * ploidy), in 64-bit, and never
changes its input. All counts here are synthetic.
"""

import copy
import logging

import numpy as np
import pytest

from robin.analysis.cnv_reference import reference_cnv_profile


def _library_profile(control, bin_width):
    cnv_from_bam = pytest.importorskip("cnv_from_bam")
    result = cnv_from_bam.iterate_bam_file(
        None,
        copy_numbers={k: list(v) for k, v in control.items()},
        bin_width=bin_width,
        log_level=int(logging.ERROR),
    )
    return result.cnv


def test_matches_the_library_on_one_contig():
    rng = np.random.default_rng(1)
    control = {"chr1": rng.integers(0, 10, 2001).tolist()}
    for bin_width in (1_000, 7_000, 100_000):
        ours = reference_cnv_profile(control, bin_width)
        assert ours == _library_profile(control, bin_width)
        assert isinstance(ours["chr1"], list)


def test_matches_the_library_on_aligned_contigs():
    rng = np.random.default_rng(2)
    # Lengths are multiples of the 7-bin chunk, so no library chunk straddles
    # two contigs and both medians see the same chunk sums.
    control = {
        "chr10": rng.integers(0, 10, 700).tolist(),
        "chr2": rng.integers(0, 10, 1400).tolist(),
        "chr1": rng.integers(0, 10, 2100).tolist(),
    }
    ours = reference_cnv_profile(control, 7_000)
    library = _library_profile(control, 7_000)
    assert ours == library
    assert list(ours) == ["chr1", "chr2", "chr10"] == list(library)


def test_wide_bins_keep_exact_sums():
    # 2 Mb bins: chr1's first chunk holds 80,000 reads, past 16 bits.
    control = {"chr1": [40] * 2001, "chr2": [20] * 2001}
    profile = reference_cnv_profile(control, 2_000_000)
    median = 20_020  # chunk sums 80,000, 40, 40,000, 20
    assert profile["chr1"][0] == 80_000 / median * 2
    assert profile["chr1"][0] / profile["chr2"][0] == pytest.approx(2.0)


def test_a_half_median_rounds_up():
    # Median 2.5 rounds to 3 (Rust's f64::round), not to 2 (Python's round).
    assert reference_cnv_profile({"chr1": [2, 3]}, 1_000) == {"chr1": [4 / 3, 2.0]}


def test_an_all_zero_control_gives_zeros():
    assert reference_cnv_profile({"chr1": [0, 0, 0]}, 1_000) == {"chr1": [0.0, 0.0, 0.0]}


def test_no_contigs_give_an_empty_profile():
    assert reference_cnv_profile({}, 1_000) == {}


def test_a_bin_width_under_one_kb_uses_single_bins():
    assert len(reference_cnv_profile({"chr1": [1, 2, 3]}, 499)["chr1"]) == 3


def test_a_partial_last_chunk_is_kept():
    profile = reference_cnv_profile({"chr1": [1, 1, 1, 1, 1]}, 2_000)
    assert profile == {"chr1": [2.0, 2.0, 1.0]}


def test_the_control_is_not_changed():
    control = {"chr1": [5, 6, 7, 8], "chr2": [1, 2, 3, 4]}
    before = copy.deepcopy(control)
    from_lists = reference_cnv_profile(control, 2_000)
    assert control == before

    frozen = {}
    for name, counts in control.items():
        arr = np.asarray(counts, dtype=np.int64)
        arr.setflags(write=False)
        frozen[name] = arr
    assert reference_cnv_profile(frozen, 2_000) == from_lists
    assert all(not arr.flags.writeable for arr in frozen.values())
