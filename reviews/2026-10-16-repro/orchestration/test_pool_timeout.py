import sys, os, asyncio, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stubs
wr = stubs.wr

async def main():
    loop = asyncio.get_running_loop()
    futures = {}
    class RF:
        def options(self, **kw): return self
        def remote(self, job):
            f = loop.create_future(); futures[job.job_id] = f; return f
    finished = []
    class Coord:
        class _on_finish:
            @staticmethod
            async def remote(job, ok, ctx, err, kind=None): finished.append((job.job_id, ok, kind))
    pool = wr.Pool("cnv", max_parallel=1, job_timeout_seconds=1)
    pool.handlers["cnv"] = (RF(), {})
    pool._coordinator = Coord
    mk = lambda i: wr.Job(i, "cnv", "cnv", ["cnv:cnv"], 0, wr.WorkflowContext(f"/r/{i}.bam"))
    await pool.enqueue(mk(1))
    print("after job1 start: running_count", pool._running_count)
    # simulate one sweep of Pool._timeout_loop body (job 1 exceeded timeout)
    for ref, (job, st) in list(pool._inflight_data.items()):
        pool._timed_out_jobs.add(job.job_id); pool._inflight_data.pop(ref, None)
        if pool._running_count > 0: pool._running_count -= 1
        await Coord._on_finish.remote(job, False, job.context, "timed out", "timed_out")
    await pool.enqueue(mk(2))          # slot looks free -> job 2 starts
    print("job2 started: running_count", pool._running_count)
    futures[1].set_result(wr.WorkflowContext("/r/1.bam"))   # job 1 (not actually killed) finishes later
    await asyncio.sleep(0.01)
    print("job1 late finish: running_count", pool._running_count, "(job2 still running)")
    t = asyncio.create_task(pool.enqueue(mk(3)))
    await asyncio.sleep(0.05)
    print("job3 started while job2 running:", 3 in futures, "| concurrently running:",
          [j for j, f in futures.items() if not f.done()], "| max_parallel:", pool.max_parallel)
asyncio.run(main())
