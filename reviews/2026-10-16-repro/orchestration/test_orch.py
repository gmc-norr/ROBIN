import sys, os, asyncio, pickle, time, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import stubs
wr = stubs.wr

WD = tempfile.mkdtemp(prefix="wd_")
PLAN = ["preprocessing:preprocessing", "bed_conversion:bed_conversion", "mgmt:mgmt",
        "cnv:cnv", "target:target", "fusion:fusion", "classification:sturgeon"]

def mkctx(fp, sid):
    return wr.WorkflowContext(fp, metadata={"work_dir": WD, "target_panel": "rCNS2",
            "bam_metadata": {"sample_id": sid}, "original_workflow": PLAN})

class _RemoteMethod:
    def __init__(self, sink): self.sink = sink
    def remote(self, job):
        # emulate Ray: arguments are serialised at call time
        self.sink.append(pickle.loads(pickle.dumps(job)))
class FakeProc:
    def __init__(self):
        self.jobs = []
        self.enqueue = _RemoteMethod(self.jobs)
        self.process = _RemoteMethod(self.jobs)

def new_coord(batching=False):
    c = wr.Coordinator(target_panel="rCNS2", enable_batching=batching)
    c.using_pools = True
    procs = {}
    for jt in ["preprocessing", "bed_conversion", "mgmt", "cnv", "target", "fusion",
               "sturgeon", "nanodx", "pannanodx", "random_forest", "igv_bam",
               "snp_analysis", "target_bam_finalize"]:
        procs[jt] = FakeProc()
    c.processors = procs
    c._spool_enabled = False
    return c, procs

# ---------------------------------------------------------------- T1
def t1_handler_rerun_on_exception():
    calls = []
    def h(job, work_dir=None):
        calls.append(work_dir)
        raise RuntimeError("real failure inside handler")
    impl = wr._wrap_real_handler(h, "mgmt")
    job = wr.Job(1, "mgmt", "mgmt", ["mgmt:mgmt"], 0, mkctx("/x/a.bam", "S1"))
    err = None
    try:
        impl(job)
    except Exception as e:
        err = e
    print("T1 handler invocations:", calls, "| final error:", repr(err))

# ---------------------------------------------------------------- T2
def t2_returned_error_counted_as_success():
    def h(job, work_dir=None):
        # mimic mgmt_handler: swallow and record under "<type>_analysis"
        job.context.add_error("mgmt_analysis", "methylartist crashed")
    impl = wr._wrap_real_handler(h, "mgmt")
    job = wr.Job(2, "mgmt", "mgmt", ["mgmt:mgmt"], 0, mkctx("/x/b.bam", "S1"))
    ctx = impl(job)
    # Pool._wait detection logic (workflow_ray.py:3516-3521)
    job_errors = [e for e in ctx.errors if isinstance(e, dict) and e.get("job_type") == job.job_type]
    print("T2 results[mgmt] =", ctx.results.get("mgmt"), "| Pool sees error:", bool(job_errors),
          "| errors recorded:", [e["job_type"] for e in ctx.errors])

# ---------------------------------------------------------------- T3
async def t3_classifier_trigger_dropped_for_other_sample():
    c, procs = new_coord()
    def bed_done(fp, sid, jid):
        ctx = mkctx(fp, sid)
        ctx.add_result("bed_conversion", {"status": "success", "parquet_path": f"{WD}/{sid}.parquet"})
        return wr.Job(jid, "bed_conversion", "bed_conversion", PLAN, 1, ctx), ctx
    j, ctx = bed_done("/r/A_1.bam", "A", 10); await c._on_finish(j, True, ctx, None)   # sturgeon(A) runs
    j, ctx = bed_done("/r/A_2.bam", "A", 11); await c._on_finish(j, True, ctx, None)   # waiting slot <- A
    j, ctx = bed_done("/r/B_1.bam", "B", 12); await c._on_finish(j, True, ctx, None)   # B: slot busy -> dropped
    # sturgeon(A) finishes
    first = procs["sturgeon"].jobs[0]
    await c._on_finish(first, True, first.context, None)
    # second sturgeon (A, from slot) finishes
    second = procs["sturgeon"].jobs[1]
    await c._on_finish(second, True, second.context, None)
    print("T3 sturgeon dispatched for samples:", [x.context.get_sample_id() for x in procs["sturgeon"].jobs],
          "| waiting slot:", c.classif_waiting_by_type.get("sturgeon"))

# ---------------------------------------------------------------- T4
async def t4_snp_dedup_pending_never_released():
    c, procs = new_coord(batching=False)
    ok1 = await c.submit_snp_analysis_job(f"{WD}/S1", sample_id="S1", reference="/ref.fa")
    j = procs["snp_analysis"].jobs[0]
    await c._on_finish(j, True, j.context, None)          # first run completes
    ok2 = await c.submit_snp_analysis_job(f"{WD}/S1", sample_id="S1", reference="/ref.fa", force_regenerate=True)
    print("T4 submit returned:", ok1, ok2, "| snp jobs dispatched:", len(procs["snp_analysis"].jobs),
          "| pending map:", c.pending, "| total_skipped:", c.total_skipped)

# ---------------------------------------------------------------- T5
async def t5_batched_job_overwrite_on_shared_ctx():
    c, procs = new_coord(batching=True)
    ctx = mkctx("/r/A_1.bam", "A")
    ctx.add_result("preprocessing", {"status": "success", "sample_id": "A"})
    prep = wr.Job(20, "preprocessing", "preprocessing", PLAN, 0, ctx)
    await c._on_finish(prep, True, ctx, None)
    # age the batcher so all types time out in the same sweep (as _batch_timeout_loop does)
    for sid in c.sample_batcher.last_job_time:
        for jt in c.sample_batcher.last_job_time[sid]:
            c.sample_batcher.last_job_time[sid][jt] -= 10
    batches = await c.sample_batcher.check_timeouts()
    regular = c._convert_batched_jobs(batches)
    await c._submit_jobs_internal(regular)
    out = []
    for jt in ["bed_conversion", "mgmt", "target", "fusion", "cnv"]:
        for j in procs[jt].jobs:
            out.append((j.job_type, j.context.metadata["_batched_job"].job_type))
    print("T5 (dispatched job_type, embedded _batched_job.job_type):", out)

# ---------------------------------------------------------------- T6
def t6_watcher_holds_files_until_next_event():
    sent = []
    class C:
        class submit_jobs:
            @staticmethod
            def remote(batch): sent.append([j.context.filepath for j in batch]); return None
    w = wr.RayFileWatcher(C, PLAN, "rCNS2", patterns=["*.bam"], work_dir=WD)
    w._handle("/run/pass/A_0.bam")          # immediately after start
    w._handle("/run/pass/A_0.bam")          # duplicate modified event (ignored, no flush)
    w._handle("/run/fail/A_0.bam")
    time.sleep(0.6)
    w._handle("/run/pass/A_0.bam")          # more modified events for known file: still no flush
    print("T6 submitted after events:", sent, "| still pending:", [j.context.filepath for j in w._pending_jobs],
          "| reference in watcher job metadata:", "reference" in w._pending_jobs[0].context.metadata if w._pending_jobs else None)

t1_handler_rerun_on_exception()
t2_returned_error_counted_as_success()
asyncio.run(t3_classifier_trigger_dropped_for_other_sample())
asyncio.run(t4_snp_dedup_pending_never_released())
asyncio.run(t5_batched_job_overwrite_on_shared_ctx())
t6_watcher_holds_files_until_next_event()
