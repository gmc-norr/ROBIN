"""CNV-5: every chromosome got a fake breakpoint at its end.

``ruptures``' ``predict()`` always returns the length of the signal as its final
element. That value is a terminator marking the end of the last segment, not a
detected change point, and ``run_ruptures`` turned every returned value into a
``(start, end)`` range. With about 24 chromosomes analysed per run, the master
BED gained about 24 spurious regions past the end of each contig.

These tests call ``run_ruptures`` directly on a synthetic signal, so they pin
the list comprehension in ``run_ruptures`` rather than the consumer's
``if start >= 0 and end > start`` guard, which cannot see the problem: the
terminator's range is positive and well-ordered, just fictional.
"""

from robin.analysis.cnv_analysis import run_ruptures


def test_the_series_end_does_not_become_a_region():
    """CNV-5: ruptures always returns len(signal) last; it is not a breakpoint."""
    import numpy as np

    signal = np.concatenate([np.zeros(100), np.ones(100) * 5.0])
    ranges = run_ruptures(signal, penalty_value=10, bin_width=1000)
    assert all(end <= len(signal) * 1000 - 1000 for _, end in ranges)
    assert len(ranges) == 1


def test_the_detected_shift_is_still_reported():
    import numpy as np

    signal = np.concatenate([np.zeros(100), np.ones(100) * 5.0])
    ((start, end),) = run_ruptures(signal, penalty_value=10, bin_width=1000)
    assert start == 100 * 1000 - 1000
    assert end == 100 * 1000 + 1000


def test_a_flat_signal_yields_no_regions():
    """Review Focus 5: ruptures returns only [len(signal)], so result[:-1] is empty."""
    import numpy as np

    assert run_ruptures(np.zeros(200), penalty_value=10, bin_width=1000) == []
