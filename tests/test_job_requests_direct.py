"""Batch 10a's external review: an accepted request never waits for room inside a finish
step, and each SNP request reaches the one-running / one-waiting rule as it was sent.

These tests call the Coordinator's own methods directly: the class behind the Ray actor,
built without ``setup()``, so no Ray cluster, pools or background loops start. Processors
only record the jobs they are sent. So the backlog, the limits and the batcher's timeout are
set exactly, not by timing.

Why a finish step must never wait for room: a pool keeps its slot busy until the
Coordinator's ``_on_finish`` for that job returns (``Pool.enqueue``'s ``_wait``), and the
``high`` scheduler handles completions one at a time in ``_drain_loop``. So a finish step
that waits for room another job must make can wait for good. The tests use small limits; the
defaults are large (20,480 waiting jobs in all), but the cycle is the same.
"""

import asyncio
import itertools
import time
from types import SimpleNamespace

import pytest

from robin import workflow_ray as wr

CoordinatorClass = wr.Coordinator.__ray_metadata__.modified_class
STUCK = 2.0  # a finish step that takes longer than this is waiting for room
_ids = itertools.count(50_000_000)


def _coordinator(tmp_path, batching=False):
    c = CoordinatorClass(target_panel="rCNS2", preset="standard", enable_batching=batching)
    c.using_pools = True
    c._spool_enabled = False
    c.set_work_dir(str(tmp_path))
    sent = []
    proc = SimpleNamespace(enqueue=SimpleNamespace(remote=sent.append))
    c.processors = {t: proc for t in ("snp_analysis", "igv_bam", "mgmt")}
    return c, sent


def _job(job_type, sample="S", **metadata):
    q = wr.job_queue_of(job_type)
    ctx = wr.WorkflowContext(
        filepath=f"/data/{sample}/{job_type}.bam",
        metadata={
            "bam_metadata": {"sample_id": sample},
            "original_workflow": [f"{q}:{job_type}"],
            **metadata,
        },
    )
    return wr.Job(next(_ids), job_type, q, [f"{q}:{job_type}"], 0, ctx)


async def _finish(c, job):
    try:
        await asyncio.wait_for(c._on_finish(job, True, job.context, None), STUCK)
    except asyncio.TimeoutError:
        pytest.fail(f"finishing {job.job_type} waited {STUCK} s for room")


async def _flush_batcher(c):
    """What _batch_timeout_loop does once every batch's timeout has passed."""
    for by_type in c.sample_batcher.last_job_time.values():
        for job_type in by_type:
            by_type[job_type] = time.time() - 3600
    batches = await c.sample_batcher.check_timeouts()
    await c._submit_jobs_internal(c._convert_batched_jobs(batches))


def _tags(sent, job_type="snp_analysis"):
    return [j.context.metadata.get("tag") for j in sent if j.job_type == job_type]


# ------------------------------------------------------------------ never wait for room


def test_a_waiting_snp_request_starts_without_waiting_for_room(tmp_path):
    # Old code: S's waiting request was dispatched as a new submission, which first waits
    # until the backlog has room. Here the backlog is an igv_bam job that leaves it only when
    # the running igv_bam job ends, and in a real pool that job can't start before S's finish
    # step returns.
    async def run():
        c, sent = _coordinator(tmp_path)
        c.max_inflight_per_type = 1
        c.max_total_waiting = 1
        first = _job("snp_analysis", tag="first")
        waiting = _job("snp_analysis", tag="waiting")
        igv_running, igv_waiting = _job("igv_bam", "T"), _job("igv_bam", "T")
        for j in (first, waiting, igv_running, igv_waiting):
            await c.submit_jobs([j])
        assert sent == [first, igv_running]
        assert c.dedup_waiting[("S", "snp_analysis")] is waiting
        assert c.waiting_by_type_global["igv_bam"] == [igv_waiting]
        assert c._total_waiting() == c.max_total_waiting
        await _finish(c, first)
        assert sent[-1] is waiting
        assert c.dedup_running[("S", "snp_analysis")] == waiting.job_id

    asyncio.run(run())


def test_a_job_from_a_waiting_list_does_not_wait_for_queue_room(tmp_path):
    # Old code: a job taken off the global waiting list, whose queue was at its running
    # limit, waited for room in that queue's waiting list, inside another job's finish step.
    # The room comes only when the queue's running job ends; on `high` that finish is handled
    # after this one, so it never comes.
    async def run():
        c, sent = _coordinator(tmp_path)
        slow = wr.job_queue_of("igv_bam")
        assert wr.job_queue_of("mgmt") != slow
        c.max_total_inflight = 2
        c.max_inflight_per_queue[slow] = 1
        c.max_waiting_per_queue[slow] = 1
        snp = _job("snp_analysis")
        igv_a, other, igv_b = _job("igv_bam", "T"), _job("mgmt", "T"), _job("igv_bam", "T")
        for j in (snp, igv_a, other, igv_b):
            await c.submit_jobs([j])
        assert sent == [snp, other]
        assert c.waiting_by_queue[slow] == [igv_a]
        assert c.waiting_global == [igv_b]
        await _finish(c, other)
        # igv_b left the global list for its queue's list, where it waits its turn
        assert c.waiting_global == []
        assert c.waiting_by_queue[slow] == [igv_a, igv_b]

    asyncio.run(run())


# ------------------------------------------------------------------ batching on (the default)


def test_batched_snp_requests_keep_the_newest_request(tmp_path):
    # Old code: the two later requests became one batch that ran with its first request's
    # context, so the follow-up run had the older settings (no forced regeneration).
    async def run():
        c, sent = _coordinator(tmp_path, batching=True)
        await c.submit_jobs([_job("snp_analysis", tag="first")])
        await _flush_batcher(c)
        assert _tags(sent) == ["first"]
        running = sent[-1]
        await c.submit_jobs([_job("snp_analysis", tag="older", force_regenerate=False)])
        await c.submit_jobs([_job("snp_analysis", tag="newer", force_regenerate=True)])
        await _flush_batcher(c)
        await _finish(c, running)
        assert _tags(sent) == ["first", "newer"]
        assert sent[-1].context.metadata["force_regenerate"] is True

    asyncio.run(run())


def test_batched_snp_requests_during_a_run_give_one_more_run(tmp_path):
    # Old code: request 2 had left the batcher for the waiting slot, request 3 was still in
    # the batcher when run 1 ended, so 3 waited behind 2: three runs for requests that all
    # came during run 1.
    async def run():
        c, sent = _coordinator(tmp_path, batching=True)
        await c.submit_jobs([_job("snp_analysis", tag="r1")])
        await _flush_batcher(c)
        running = sent[-1]
        await c.submit_jobs([_job("snp_analysis", tag="r2")])
        await _flush_batcher(c)
        await c.submit_jobs([_job("snp_analysis", tag="r3")])
        await _finish(c, running)
        second = sent[-1]
        await _flush_batcher(c)
        await _finish(c, second)
        await _flush_batcher(c)
        assert _tags(sent) == ["r1", "r3"]

    asyncio.run(run())
