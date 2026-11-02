"""A sample's CNV read counts are shared by every worker that adds to them (batch 10c:
CNV-2). Synthetic BAMs and control; real separate processes stand in for Ray workers.

- No worker keeps its own copy of a sample's counts: each run reads them from disk, so
  a batch another worker added is never lost.
- Two runs for one sample never overlap: a per-sample lock is held from reading the
  counts to saving them.
- Counts that can't be read fail the run and are kept, never reset to none.
- Counts are written whole or not at all.
"""

import logging
import multiprocessing
import os
import pickle

import pytest

pytest.importorskip("cnv_from_bam")

from cnv_fixtures import bam_metadata, cnv_worker, flat_control, write_bam  # noqa: E402
from robin.analysis import cnv_analysis  # noqa: E402

LOG = logging.getLogger("test_cnv_counts_shared")
CONTIGS = {"chr1": 2_000_000, "chr2": 2_000_000}
KINDS = ["single", "multiple"]


def _bam(tmp_path, name, contig, n):
    """A BAM of ``n`` reads on ``contig``, 1 kb apart."""
    return write_bam(str(tmp_path / name), CONTIGS, {contig: list(range(0, n * 1000, 1000))})


def _counts(work_dir, sample="S"):
    with open(work_dir / sample / f"{sample}_copy_numbers.pkl", "rb") as f:
        counts = pickle.load(f)
    return {name: int(sum(bins)) for name, bins in counts.items() if sum(bins)}


# Workers in their own processes ---------------------------------------------------------


@pytest.fixture
def workers():
    """Start stand-in Ray workers, each in its own process; stop them after the test."""
    ctx = multiprocessing.get_context("spawn")
    started = []

    def start(pause_seconds=0.0):
        parent, child = ctx.Pipe()
        proc = ctx.Process(target=cnv_worker, args=(child, CONTIGS, 1, pause_seconds), daemon=True)
        proc.start()
        started.append((proc, parent))
        return parent

    yield start
    for proc, conn in started:
        try:
            conn.send(None)
        except Exception:
            pass
        proc.join(10)
        if proc.is_alive():
            proc.kill()


def _send(conn, kind, bams, work_dir):
    conn.send((kind, bams, "S", str(work_dir)))


def _answer(conn):
    assert conn.poll(180), "the worker did not answer"
    return conn.recv()


@pytest.mark.parametrize("kind", KINDS)
def test_a_worker_keeps_the_batch_another_worker_added(tmp_path, workers, kind):
    # CNV-2's example: batch 1 on worker A, batch 2 on worker B, batch 3 on worker A.
    work = tmp_path / "run"
    first = _bam(tmp_path, "1.bam", "chr1", 100)
    second = _bam(tmp_path, "2.bam", "chr2", 200)
    third = _bam(tmp_path, "3.bam", "chr1", 50)
    a, b = workers(), workers()
    for conn, bam in ((a, first), (b, second), (a, third)):
        _send(conn, kind, [bam], work)
        assert _answer(conn) is None
    assert _counts(work) == {"chr1": 150, "chr2": 200}, "a batch another worker added was lost"
    with open(work / "S" / "cnv_analysis_counter.txt") as f:
        assert f.read().strip() == "3"


@pytest.mark.parametrize("kind", KINDS)
def test_two_runs_for_one_sample_at_once_keep_both_batches(tmp_path, workers, kind):
    work = tmp_path / "run"
    first = _bam(tmp_path, "1.bam", "chr1", 100)
    second = _bam(tmp_path, "2.bam", "chr2", 200)
    a, b = workers(pause_seconds=1.5), workers(pause_seconds=1.5)
    _send(a, kind, [first], work)
    _send(b, kind, [second], work)
    assert _answer(a) is None
    assert _answer(b) is None
    assert _counts(work) == {"chr1": 100, "chr2": 200}, "one run overwrote the other's counts"


# In this process ------------------------------------------------------------------------


@pytest.fixture
def control(monkeypatch):
    monkeypatch.setattr(cnv_analysis, "get_cached_ref_cnv_dict", lambda path, logger: flat_control(CONTIGS, 1))
    _forget_cached_counts()
    yield
    _forget_cached_counts()


def _forget_cached_counts():
    """What a new worker process starts with (before this batch: an empty cache)."""
    clear = getattr(cnv_analysis, "clear_sample_cache", None)
    if clear is not None:
        clear()


def _run(kind, bam, work):
    if kind == "single":
        return cnv_analysis.process_single_bam(bam, bam_metadata("S", bam), str(work), LOG, threads=1)
    return cnv_analysis.process_multiple_bams([bam], [bam_metadata("S", bam)], str(work), LOG, threads=1)


@pytest.mark.parametrize("kind", KINDS)
def test_counts_that_cannot_be_read_fail_the_run_and_are_kept(tmp_path, control, kind):
    work = tmp_path / "run"
    assert _run(kind, _bam(tmp_path, "1.bam", "chr1", 100), work)["error_message"] is None
    path = work / "S" / "S_copy_numbers.pkl"
    path.write_bytes(b"not a pickle")
    _forget_cached_counts()
    result = _run(kind, _bam(tmp_path, "2.bam", "chr2", 200), work)
    assert result["error_message"], "the run went on with no counts"
    assert path.read_bytes() == b"not a pickle", "the unreadable counts were replaced"


@pytest.mark.parametrize("kind", KINDS)
def test_a_failed_counts_write_keeps_the_last_good_counts(tmp_path, control, monkeypatch, kind):
    work = tmp_path / "run"
    assert _run(kind, _bam(tmp_path, "1.bam", "chr1", 100), work)["error_message"] is None
    _forget_cached_counts()

    def fail_part_way(obj, f, *args, **kwargs):
        f.write(b"partial")
        raise OSError("disk full")

    monkeypatch.setattr(cnv_analysis.pickle, "dump", fail_part_way)
    result = _run(kind, _bam(tmp_path, "2.bam", "chr2", 200), work)
    monkeypatch.undo()
    assert result["error_message"]
    assert _counts(work) == {"chr1": 100}, "a write that failed part-way destroyed the counts"


# The legacy multi-sample file (round 2: the external review) ---------------------------


def _legacy_file(tmp_path, work, holds):
    """Write ``<work>/S/update_cnv_dict.pkl``, the old multi-sample counts file: each
    sample in ``holds`` gets 100 reads on chr1, made by a real run in another folder."""
    seed = tmp_path / "seed"
    assert _run("single", _bam(tmp_path, "0.bam", "chr1", 100), seed)["error_message"] is None
    with open(seed / "S" / "S_copy_numbers.pkl", "rb") as f:
        counts = pickle.load(f)
    folder = work / "S"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / "update_cnv_dict.pkl"
    path.write_bytes(pickle.dumps({name: counts for name in holds}))
    return path


@pytest.mark.parametrize("kind", KINDS)
def test_a_failed_move_out_of_the_legacy_file_fails_the_run(tmp_path, control, monkeypatch, kind):
    work = tmp_path / "run"
    legacy = _legacy_file(tmp_path, work, ["S", "T"])
    before = legacy.read_bytes()
    real_save, calls = cnv_analysis._save_copy_numbers, []

    def first_write_fails(path, counts):  # the move's write; the run's own save works
        calls.append(path)
        if len(calls) == 1:
            raise OSError("rename failed")
        return real_save(path, counts)

    monkeypatch.setattr(cnv_analysis, "_save_copy_numbers", first_write_fails)
    result = _run(kind, _bam(tmp_path, "2.bam", "chr2", 200), work)
    assert result["error_message"], "the run went on without the sample's earlier counts"
    assert not (work / "S" / "S_copy_numbers.pkl").exists(), "counts were saved without the legacy ones"
    assert legacy.read_bytes() == before, "the legacy file was changed"

    # The next run moves the counts and adds to them; the failed run's BAM is not in them.
    monkeypatch.setattr(cnv_analysis, "_save_copy_numbers", real_save)
    assert _run(kind, _bam(tmp_path, "3.bam", "chr2", 50), work)["error_message"] is None
    assert _counts(work) == {"chr1": 100, "chr2": 50}


@pytest.mark.parametrize("kind", KINDS)
def test_a_legacy_file_that_cannot_be_read_fails_the_run_and_is_kept(tmp_path, control, kind):
    work = tmp_path / "run"
    (work / "S").mkdir(parents=True)
    legacy = work / "S" / "update_cnv_dict.pkl"
    legacy.write_bytes(b"not a pickle")
    result = _run(kind, _bam(tmp_path, "2.bam", "chr2", 200), work)
    assert result["error_message"], "the run went on with no counts"
    assert not (work / "S" / "S_copy_numbers.pkl").exists(), "counts were saved without the legacy ones"
    assert legacy.read_bytes() == b"not a pickle", "the legacy file was changed"


@pytest.mark.parametrize("kind", KINDS)
def test_a_legacy_file_is_moved_and_added_to(tmp_path, control, kind):
    # Control: the move works as before.
    work = tmp_path / "run"
    legacy = _legacy_file(tmp_path, work, ["S", "T"])
    assert _run(kind, _bam(tmp_path, "2.bam", "chr2", 200), work)["error_message"] is None
    assert _counts(work) == {"chr1": 100, "chr2": 200}
    assert legacy.exists()


@pytest.mark.parametrize("kind", KINDS)
def test_a_legacy_file_without_the_sample_starts_from_none(tmp_path, control, kind):
    # Control: a legacy file that doesn't hold this sample is not a failure.
    work = tmp_path / "run"
    _legacy_file(tmp_path, work, ["T"])
    assert _run(kind, _bam(tmp_path, "2.bam", "chr2", 200), work)["error_message"] is None
    assert _counts(work) == {"chr2": 200}


def test_the_workflow_still_imports_the_cnv_handler():
    from robin import workflow_ray

    assert workflow_ray._cnv_handler is not None, workflow_ray._HANDLER_IMPORT_ERRORS.get("cnv")
