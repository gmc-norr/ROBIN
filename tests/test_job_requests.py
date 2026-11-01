"""JOB-1 and JOB-5: a request for a sample's classifier or SNP run is never thrown away.

The rule (user's choice "A", 2026-11-01): for each sample and job type, at most one job runs and
at most one waits. A new request joins the waiting one if there is one, so it is never dropped,
and a merged run starts after the newest request arrived. Each classifier type still runs one
job at a time across all samples (as before), and waiting samples take turns, oldest first.

Every case goes through a real local Ray ``Coordinator`` with the synthetic handlers of
``failure_fixtures.py``. ``sleep`` keeps a job busy so the next requests arrive while it runs;
each test first checks that they did (a fresh Ray worker takes about 2 s to start, so a short
busy job can end before the next request lands). A classifier run reads the sample's whole
merged methylation file, so the check is: each sample's last classifier run started after its
last bed_conversion ended.

The tests that need a request to land while a job runs use the ``standard`` scheduler (the
default) only. On ``high`` the Coordinator answers nothing while a job runs (measured: a
``stats`` call took the job's full 4 s), so no request can land mid-run there. The tests without
that need run on both schedulers, with batching off and on.

On the old code every fault test runs and fails by assertion; the guards pass there.
"""

import os
import time

import pytest
import ray

from failure_fixtures import SCHEDULERS, Harness, init_ray, shutdown_ray

SLOW = 6.0  # seconds a busy job takes
WAIT = 60.0  # longest a test waits for the runs it expects


@pytest.fixture(scope="module")
def ray_cluster():
    init_ray()
    yield
    shutdown_ray()


def _harness(request, tmp_path, scheduler, enable_batching=False):
    harness = Harness(scheduler, str(tmp_path / "work"), enable_batching=enable_batching)
    request.addfinalizer(harness.close)
    return harness


@pytest.fixture
def harness(request, ray_cluster, tmp_path):
    """The standard scheduler (see the module doc for why not high)."""
    return _harness(request, tmp_path, "standard")


def _wait(predicate, timeout=WAIT):
    t0 = time.time()
    while time.time() - t0 < timeout:
        if predicate():
            return True
        time.sleep(0.2)
    return False


def _beds(harness, samples):
    """New methylation data: one bed_conversion per sample (sent together), each triggering
    Sturgeon. Returns once all of them have ended."""
    ids = [
        harness.submit("bed_conversion", s, downstream=("sturgeon",), sleep={"sturgeon": SLOW})
        for s in samples
    ]
    assert _wait(lambda: {r["job_id"] for r in harness.ended("bed_conversion")} >= set(ids))
    return max(r["end"] for r in harness.ended("bed_conversion") if r["job_id"] in ids)


def _last_bed_end(harness, sample):
    return max(r["end"] for r in harness.ended("bed_conversion", sample))


def _first_run_end(harness, job_type, sample):
    assert _wait(lambda: harness.ended(job_type, sample))
    return harness.ended(job_type, sample)[0]["end"]


# ------------------------------------------------------------------ JOB-1: classifiers


def test_every_sample_gets_a_classifier_run_after_its_data(harness):
    # Old code: A runs, B takes the one waiting slot, C's request is thrown away.
    _beds(harness, ["A"])
    assert _wait(lambda: harness.started("sturgeon", "A"))
    sent = _beds(harness, ["B", "C"])
    assert sent < _first_run_end(harness, "sturgeon", "A"), "B and C must arrive during A's run"
    _wait(lambda: {r["sample"] for r in harness.ended("sturgeon")} >= {"A", "B", "C"})
    ran_after = {
        s for s in ("A", "B", "C")
        if any(r["start"] > _last_bed_end(harness, s) for r in harness.ended("sturgeon", s))
    }
    assert ran_after == {"A", "B", "C"}


def test_new_data_while_another_sample_waits_gets_its_own_run(harness):
    # Old code: A runs, B waits, then A's new data finds the slot taken and is dropped, so A's
    # last classification is older than its data.
    _beds(harness, ["A"])
    assert _wait(lambda: harness.started("sturgeon", "A"))
    _beds(harness, ["B"])
    sent = _beds(harness, ["A"])
    assert sent < _first_run_end(harness, "sturgeon", "A"), "B, then A, must arrive during A's run"
    _wait(lambda: len(harness.ended("sturgeon", "A")) >= 2 and harness.ended("sturgeon", "B"))
    for sample in ("A", "B"):
        runs = harness.ended("sturgeon", sample)
        assert runs and runs[-1]["start"] > _last_bed_end(harness, sample), sample
    # Samples take turns: B waited first, so B runs before A's second run.
    assert harness.ended("sturgeon", "B")[0]["start"] < harness.ended("sturgeon", "A")[-1]["start"]


def test_guard_requests_for_one_sample_merge_into_one_run(harness):
    # Passes before and after: three new batches for A during its run give one more run.
    _beds(harness, ["A"])
    assert _wait(lambda: harness.started("sturgeon", "A"))
    sent = _beds(harness, ["A", "A", "A"])
    assert sent < _first_run_end(harness, "sturgeon", "A"), "the batches must arrive during A's run"
    _wait(lambda: len(harness.ended("sturgeon", "A")) >= 2)
    time.sleep(SLOW + 2)  # room for a wrongly queued third run to show
    runs = harness.ended("sturgeon", "A")
    assert len(runs) == 2
    assert runs[-1]["start"] > _last_bed_end(harness, "A")


# ------------------------------------------------------------------ JOB-5: SNP requests


@pytest.fixture(params=[(s, b) for s in SCHEDULERS for b in (False, True)],
                ids=[f"{s}-{'batched' if b else 'unbatched'}" for s in SCHEDULERS for b in (False, True)])
def snp_harness(request, ray_cluster, tmp_path):
    scheduler, batching = request.param
    return _harness(request, tmp_path, scheduler, enable_batching=batching)


@pytest.fixture(params=[False, True], ids=["standard-unbatched", "standard-batched"])
def snp_harness_standard(request, ray_cluster, tmp_path):
    return _harness(request, tmp_path, "standard", enable_batching=request.param)


def _snp_request(harness, sample):
    """What the GUI does: submit_snp_analysis_job on the sample's folder."""
    sample_dir = os.path.join(harness.work_dir, sample)
    os.makedirs(sample_dir, exist_ok=True)
    return ray.get(harness.coord.submit_snp_analysis_job.remote(sample_dir, sample))


def test_a_second_snp_request_runs(snp_harness):
    # Old code: the first request is never cleared, so every later one is skipped while
    # submit_snp_analysis_job still returns True ("queued").
    h = snp_harness
    assert _snp_request(h, "S") is True
    assert _wait(lambda: len(h.ended("snp_analysis", "S")) == 1)
    assert _snp_request(h, "S") is True
    _wait(lambda: len(h.ended("snp_analysis", "S")) == 2)
    assert len(h.ended("snp_analysis", "S")) == 2


def test_snp_requests_while_one_runs_give_one_more_run(snp_harness_standard):
    # Old code: both later requests are skipped. After: they merge into one run after the first,
    # and that run is the newest request (it carries the newest settings, e.g. force_regenerate).
    h = snp_harness_standard
    h.submit("snp_analysis", "S", sleep={"snp_analysis": SLOW})
    assert _wait(lambda: h.started("snp_analysis", "S"))
    sent = time.time()
    h.submit("snp_analysis", "S", sleep={"snp_analysis": SLOW})
    newest = h.submit("snp_analysis", "S", sleep={"snp_analysis": SLOW})
    sent_done = time.time()
    assert sent_done < _first_run_end(h, "snp_analysis", "S"), "requests must arrive during the run"
    _wait(lambda: len(h.ended("snp_analysis", "S")) >= 2)
    time.sleep(SLOW + 2)  # room for a wrongly queued third run to show
    runs = h.ended("snp_analysis", "S")
    assert len(runs) == 2
    assert runs[-1]["start"] > sent
    # With batching on too: a per-sample request skips the batcher, which ran a batch with
    # its first request's settings (batch 10a's external review)
    assert runs[-1]["job_id"] == newest


def test_guard_snp_requests_for_two_samples_both_run(snp_harness):
    # Passes before and after: different samples never held each other back.
    h = snp_harness
    assert _snp_request(h, "S") is True
    assert _snp_request(h, "T") is True
    _wait(lambda: h.ended("snp_analysis", "S") and h.ended("snp_analysis", "T"))
    assert len(h.ended("snp_analysis", "S")) == 1
    assert len(h.ended("snp_analysis", "T")) == 1
