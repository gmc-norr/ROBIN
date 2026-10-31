"""CNV-1: no sample's reads ever reach the CNV control.

ROBIN's second ``iterate_bam_file`` pass handed the cached control to the
library as ``copy_numbers``, which adds the BAM's reads into it. So the
"control" grew with every BAM of every sample in a long-lived worker: real
gains and losses faded, and the next sample showed mirror-image false events.
That pass also summed chunks as 16-bit numbers, which wrapped at the wide bins
of an early run. All BAMs and controls here are synthetic.
"""

import copy
import logging
import os
import pickle
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("cnv_from_bam")

from cnv_fixtures import bam_metadata, every_kb, flat_control, write_bam  # noqa: E402
from robin.analysis import cnv_analysis  # noqa: E402
from robin.analysis.cnv_reference import reference_cnv_profile  # noqa: E402

LOG = logging.getLogger("test_cnv_control_reference")
CONTIGS = {"chr1": 2_000_000, "chr2": 2_000_000}


def _changed_contigs(after, before):
    """Contigs whose counts differ (a short failure message: pytest's diff of
    two long lists takes minutes)."""
    return sorted(c for c in set(after) | set(before) if list(after.get(c, [])) != list(before.get(c, [])))


def _direct(bam, control, copy_numbers=None):
    result = cnv_analysis.run_cnv_analysis_direct(
        bam, {} if copy_numbers is None else copy_numbers, control, LOG, threads=1
    )
    assert result["success"], result
    return result


def test_a_run_leaves_the_control_unchanged(tmp_path):
    bam = write_bam(str(tmp_path / "s.bam"), CONTIGS, {c: every_kb(n) for c, n in CONTIGS.items()})
    control = flat_control(CONTIGS, 3)
    before = copy.deepcopy(control)
    _direct(bam, control)
    assert _changed_contigs(control, before) == []


def test_the_reference_profile_is_the_control_alone(tmp_path):
    reads = {"chr1": every_kb(CONTIGS["chr1"], 2), "chr2": every_kb(CONTIGS["chr2"], 1)}
    bam = write_bam(str(tmp_path / "s.bam"), CONTIGS, reads)
    control = flat_control(CONTIGS, 3)
    result = _direct(bam, control)
    assert result["r2_cnv"] == reference_cnv_profile(flat_control(CONTIGS, 3), result["r_bin"])
    assert result["r2_cnv"]["chr1"] == result["r2_cnv"]["chr2"]


def test_an_early_run_reference_does_not_wrap(tmp_path):
    # 150 reads over 4 Mb give 2 Mb bins; the control's chr1 chunk then
    # holds 80,000 reads, past 16 bits.
    reads = {c: list(range(0, 1_500_000, 20_000)) for c in CONTIGS}
    bam = write_bam(str(tmp_path / "s.bam"), CONTIGS, reads)
    control = {"chr1": [40] * 2001, "chr2": [20] * 2001}
    result = _direct(bam, control)
    assert result["r_bin"] == 2_000_000
    r2 = result["r2_cnv"]
    assert r2["chr1"][0] / r2["chr2"][0] == pytest.approx(2.0)


@pytest.fixture
def cached_control(monkeypatch):
    """Patch the cached control to a synthetic one; return a setter."""
    holder = {}
    monkeypatch.setattr(cnv_analysis, "get_cached_ref_cnv_dict", lambda path, logger: holder["control"])
    cnv_analysis.clear_sample_cache()
    yield holder
    cnv_analysis.clear_sample_cache()


def _single(bam, sample_id, work_dir):
    result = cnv_analysis.process_single_bam(bam, bam_metadata(sample_id, bam), str(work_dir), LOG, threads=1)
    assert result["error_message"] is None, result
    return np.load(os.path.join(work_dir, sample_id, "CNV3.npy"), allow_pickle=True).item()


def test_a_normal_sample_after_an_abnormal_one_is_unchanged(tmp_path, cached_control):
    gain = {"chr1": every_kb(CONTIGS["chr1"], 3), "chr2": every_kb(CONTIGS["chr2"], 1)}
    normal = {c: every_kb(n) for c, n in CONTIGS.items()}
    bam_a = write_bam(str(tmp_path / "a.bam"), CONTIGS, gain)
    bam_b = write_bam(str(tmp_path / "b.bam"), CONTIGS, normal)

    cached_control["control"] = flat_control(CONTIGS, 1)
    _single(bam_a, "A", tmp_path / "run")
    after_a = _single(bam_b, "B", tmp_path / "run")

    cnv_analysis.clear_sample_cache()
    cached_control["control"] = flat_control(CONTIGS, 1)
    fresh = _single(bam_b, "B", tmp_path / "fresh")

    assert list(after_a) == list(fresh)
    for name in fresh:
        np.testing.assert_array_equal(after_a[name], fresh[name])


@pytest.fixture
def control_file(tmp_path):
    """A pickled synthetic control, with the module's control cache cleared around the test."""
    path = tmp_path / "control.pkl"
    with open(path, "wb") as f:
        pickle.dump(flat_control(CONTIGS, 3), f)
    cnv_analysis.clear_ref_cnv_dict_cache()
    yield str(path)
    cnv_analysis.clear_ref_cnv_dict_cache()


def test_the_cached_control_cannot_be_changed(control_file):
    cached = cnv_analysis.get_cached_ref_cnv_dict(control_file, LOG)
    assert _changed_contigs(cached, flat_control(CONTIGS, 3)) == []
    with pytest.raises(TypeError):
        cached["chr1"] = [0]
    with pytest.raises(ValueError):
        cached["chr1"][0] = 9


def test_the_cached_control_cannot_be_passed_to_the_library(tmp_path, control_file):
    import cnv_from_bam

    bam = write_bam(str(tmp_path / "s.bam"), CONTIGS, {c: every_kb(n) for c, n in CONTIGS.items()})
    cached = cnv_analysis.get_cached_ref_cnv_dict(control_file, LOG)
    with pytest.raises(TypeError):
        cnv_from_bam.iterate_bam_file(bam, copy_numbers=cached, log_level=int(logging.ERROR))
    assert _changed_contigs(cached, flat_control(CONTIGS, 3)) == []


def test_the_subprocess_runner_still_writes_the_control(tmp_path, control_file, monkeypatch):
    monkeypatch.setattr(cnv_analysis.subprocess, "run", lambda *a, **k: SimpleNamespace(returncode=1))
    out = tmp_path / "out"
    assert cnv_analysis.run_cnv_analysis_subprocess(str(tmp_path / "s.bam"), {}, control_file, str(out), LOG) is None
    with open(out / "ref_cnv_dict.pkl", "rb") as f:
        written = pickle.load(f)
    assert _changed_contigs(written, flat_control(CONTIGS, 3)) == []


@pytest.mark.parametrize("as_arrays", [False, True], ids=["lists", "arrays"])
def test_the_subprocess_script_reference_is_the_control_alone(tmp_path, as_arrays):
    from robin.analysis import cnv_subprocess

    reads = {"chr1": every_kb(CONTIGS["chr1"], 2), "chr2": every_kb(CONTIGS["chr2"], 1)}
    bam = write_bam(str(tmp_path / "s.bam"), CONTIGS, reads)
    control = flat_control(CONTIGS, 3)
    # The shipped control holds lists; run_cnv_analysis_subprocess writes
    # the cached control's arrays.
    stored = dict(cnv_analysis._read_only_control(control)) if as_arrays else control
    path = tmp_path / "control.pkl"
    with open(path, "wb") as f:
        pickle.dump(stored, f)
    result = cnv_subprocess.run_cnv_analysis(str(bam), ref_cnv_dict_path=str(path), output_dir=str(tmp_path / "out"))
    assert result["success"], result
    assert result["r2_cnv"] == reference_cnv_profile(control, result["r_bin"])


def test_a_batch_counts_each_bam_once(tmp_path, cached_control):
    # CNV-3: after the loop had added every BAM, a "final" call added the
    # last one again, and the saved counts kept it twice.
    first = write_bam(str(tmp_path / "1.bam"), CONTIGS, {"chr1": list(range(0, 300_000, 1000))})
    second = write_bam(str(tmp_path / "2.bam"), CONTIGS, {"chr2": list(range(0, 200_000, 1000))})
    cached_control["control"] = flat_control(CONTIGS, 1)
    bams = [first, second]
    result = cnv_analysis.process_multiple_bams(
        bams, [bam_metadata("S", b) for b in bams], str(tmp_path / "run"), LOG, threads=1
    )
    assert result["error_message"] is None, result
    with open(tmp_path / "run" / "S" / "S_copy_numbers.pkl", "rb") as f:
        counts = pickle.load(f)
    assert {name: sum(bins) for name, bins in counts.items()} == {"chr1": 300, "chr2": 200}


def test_a_contig_only_the_sample_has_is_skipped(tmp_path, cached_control):
    contigs = dict(CONTIGS, chrUn=100_000)
    bam = write_bam(str(tmp_path / "s.bam"), contigs, {c: every_kb(n) for c, n in contigs.items()})
    cached_control["control"] = flat_control(CONTIGS, 1)
    normalized = _single(bam, "S", tmp_path / "run")
    assert set(normalized) == {"chr1", "chr2"}
