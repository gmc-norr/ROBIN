"""Fault harness for the failure-visibility tests (batch 7).

It drives synthetic handlers through a REAL local Ray ``Coordinator`` on either scheduler:
``standard`` (grouped pool actors) and ``high`` (one TypeProcessor per job type). Nothing about
the scheduler is mocked; only the job handlers are synthetic, and each is wrapped exactly like
the real ones (``_wrap_real_handler``).

How it works, and why:
- A handler's fault switch travels in the job's metadata (``ctx.metadata["fault"]``): Ray workers
  do not see the driver's argv or any module global changed after the worker started.
- Every job carries a temp ``work_dir`` (otherwise the Coordinator writes ``job_durations.csv``
  into the current directory) and a ``target_panel`` (otherwise the wrapper takes its ValueError
  fallback path).
- ``original_workflow`` is set to just the job's own type, so a successful job does not trigger
  downstream jobs (which would run the real handlers).
- ``failed_jobs.csv`` is only written by ``Coordinator.shutdown()``, so ``Harness.close()`` returns
  its rows.
- Only one Coordinator can exist at a time (pools find it by the name ``robin_coordinator``).
"""

import csv
import itertools
import json
import os
import time
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import ray

from robin import workflow_ray as wr

SCHEDULERS = ("standard", "high")

# The 13 job types, in the order the real pipeline lists them.
JOB_TYPES = (
    "preprocessing", "bed_conversion", "mgmt", "cnv", "target", "fusion", "sturgeon",
    "nanodx", "pannanodx", "random_forest", "igv_bam", "snp_analysis", "target_bam_finalize",
)

# The name each real handler files its errors under (read on 92c7c00). cnv files none today (it
# only logs); Task 2 makes the real handler file "cnv_analysis", so the harness uses that.
ERROR_NAME = {
    "preprocessing": "preprocessing",
    "bed_conversion": "bed_conversion",
    "mgmt": "mgmt_analysis",
    "cnv": "cnv_analysis",
    "target": "target_analysis",
    "fusion": "fusion_analysis",
    "sturgeon": "sturgeon_analysis",
    "nanodx": "nanodx_analysis",
    "pannanodx": "pannanodx_analysis",
    "random_forest": "random_forest_analysis",
    "igv_bam": "igv_bam",
    "snp_analysis": "snp_analysis",
    "target_bam_finalize": "target_bam_finalize",
}

# Which pool (standard preset) hosts each job type.
POOL_OF = {
    "preprocessing": "prep", "bed_conversion": "bed_conversion", "mgmt": "analysis",
    "cnv": "cnv", "target": "analysis", "fusion": "analysis", "sturgeon": "classif",
    "nanodx": "classif", "pannanodx": "classif", "random_forest": "rf", "igv_bam": "slow",
    "snp_analysis": "slow", "target_bam_finalize": "slow",
}

# Types the wrapper treats as "expected failure" on a fail-only BAM submission.
FAIL_ONLY_TYPES = tuple(sorted(wr.CLASSIFICATION_TYPES | {"bed_conversion"}))

FAULTS = ("ok", "raise", "file_error", "skip")


def _synthetic_handler(job):
    """Runs inside a Ray worker. The fault comes from the job's metadata."""
    fault = job.context.metadata.get("fault", "ok")
    name = ERROR_NAME[job.job_type]
    if fault == "raise":
        raise RuntimeError(f"synthetic {job.job_type} crash")
    if fault == "file_error":
        job.context.add_error(name, f"synthetic {job.job_type} failure")
    elif fault == "skip":
        job.context.add_metadata("skip_downstream", True)


def _make_remote(job_type: str):
    inner = wr._wrap_real_handler(_synthetic_handler, job_type)

    def outer(job):
        # Record what the wrapped handler left on the context (results and errors), whether it
        # returned or raised, so the test can see exactly what the real wrapper produced.
        try:
            return inner(job)
        finally:
            out_dir = os.path.join(job.context.metadata["work_dir"], "outcomes")
            os.makedirs(out_dir, exist_ok=True)
            with open(os.path.join(out_dir, f"{job.job_id}.json"), "w") as fh:
                json.dump(
                    {"results": job.context.results, "errors": job.context.errors}, fh, default=str
                )

    return ray.remote(outer)


@dataclass
class Outcome:
    """What one submitted job did, as seen from outside the Coordinator."""

    job_type: str
    fault: str
    scheduler: str
    job_id: int
    d_completed: int = 0
    d_failed: int = 0
    d_skipped: int = 0
    d_failed_by_type: Dict[str, int] = field(default_factory=dict)
    d_completed_by_type: Dict[str, int] = field(default_factory=dict)
    result: Dict[str, Any] = field(default_factory=dict)  # ctx.results after the wrapper
    errors: List[Dict[str, Any]] = field(default_factory=list)  # ctx.errors after the wrapper
    sample_id: str = ""


class Harness:
    """One real Coordinator on one scheduler, with synthetic handlers for all 13 job types."""

    _ids = itertools.count(1_000_000)

    def __init__(self, scheduler: str, work_dir: str):
        assert scheduler in SCHEDULERS
        self.scheduler = scheduler
        self.work_dir = work_dir
        os.makedirs(work_dir, exist_ok=True)
        self.coord = wr.Coordinator.options(name="robin_coordinator").remote(
            target_panel="rCNS2", analysis_workers=1, preset=scheduler, enable_batching=False
        )
        ray.get(self.coord.setup.remote())
        ray.get(self.coord.set_work_dir.remote(work_dir))
        for jt in JOB_TYPES:
            remote_fn = _make_remote(jt)
            if scheduler == "standard":
                pool = ray.get_actor(f"pool_{POOL_OF[jt]}")
                ray.get(pool.register_handler.remote(jt, remote_fn, {}))
            else:
                ray.get(self.coord.register_handler.remote(jt, remote_fn))
        self._closed = False

    def _counters(self) -> Dict[str, Any]:
        s = ray.get(self.coord.stats.remote())
        return {
            "completed": int(s.get("completed", 0)),
            "failed": int(s.get("failed", 0)),
            "skipped": int(s.get("total_skipped", 0)),
            "failed_by_type": dict(s.get("failed_by_type", {})),
            "completed_by_type": dict(s.get("completed_by_type", {})),
        }

    @staticmethod
    def _delta(after: Dict[str, int], before: Dict[str, int]) -> Dict[str, int]:
        keys = set(after) | set(before)
        return {k: after.get(k, 0) - before.get(k, 0) for k in keys if after.get(k, 0) != before.get(k, 0)}

    def run_one(
        self,
        job_type: str,
        fault: str,
        tag: str = "",
        fail_only: bool = False,
        timeout: float = 60.0,
        pre_error: bool = False,
        downstream: bool = False,
    ) -> Outcome:
        """Submit one job of ``job_type`` whose handler does ``fault``; wait for it to finish.

        ``fault``: ok | raise | file_error | skip. ``fail_only`` makes the job a fail-only BAM
        submission (metadata flag plus a file name with "fail" and not "pass").
        ``pre_error`` puts an error from an "earlier job" in the context before submission.
        ``downstream`` lets the job's success trigger its real downstream jobs (wr.TRIGGERS), by
        listing them in original_workflow; they run on the synthetic handlers too.
        Each call uses its own sample id and file name, so dedup and the one-in-flight
        classifier rules never hold a job back. Raises (with the stats printed) on timeout.
        """
        assert fault in FAULTS, fault
        job_id = next(Harness._ids)
        sample = f"{self.scheduler}_{job_type}_{fault}{'_fo' if fail_only else ''}_{tag}_{job_id}"
        kind = "fail" if fail_only else "pass"
        steps = [f"{wr.job_queue_of(job_type)}:{job_type}"]
        if downstream:
            steps += [f"{wr.job_queue_of(t)}:{t}" for t in wr.TRIGGERS.get(job_type, [])]
        ctx = wr.WorkflowContext(
            filepath=os.path.join(self.work_dir, f"{sample}_{kind}_0.bam"),
            metadata={
                "bam_metadata": {"sample_id": sample},
                "work_dir": self.work_dir,
                "target_panel": "rCNS2",
                "original_workflow": steps,
                "fault": fault,
                "fail_only_bam_submission": bool(fail_only),
            },
        )
        if pre_error:
            ctx.add_error("earlier_job", "left over from an earlier job in this context")
        job = wr.Job(job_id, job_type, wr.job_queue_of(job_type), [f"{wr.job_queue_of(job_type)}:{job_type}"], 0, ctx)
        before = self._counters()
        ray.get(self.coord.submit_jobs.remote([job]))
        t0 = time.time()
        while True:
            now = self._counters()
            if (now["completed"] + now["failed"] + now["skipped"]) > (
                before["completed"] + before["failed"] + before["skipped"]
            ):
                break
            if time.time() - t0 > timeout:
                raise TimeoutError(
                    f"{self.scheduler}/{job_type}/{fault}: job {job_id} did not finish in "
                    f"{timeout}s; stats now {now}, before {before}"
                )
            time.sleep(0.1)
        rec = {"results": {}, "errors": []}
        rec_path = os.path.join(self.work_dir, "outcomes", f"{job_id}.json")
        t1 = time.time()
        while not os.path.exists(rec_path) and time.time() - t1 < 10:
            time.sleep(0.05)
        if os.path.exists(rec_path):
            with open(rec_path) as fh:
                rec = json.load(fh)
        return Outcome(
            job_type=job_type,
            fault=fault,
            scheduler=self.scheduler,
            job_id=job_id,
            d_completed=now["completed"] - before["completed"],
            d_failed=now["failed"] - before["failed"],
            d_skipped=now["skipped"] - before["skipped"],
            d_failed_by_type=self._delta(now["failed_by_type"], before["failed_by_type"]),
            d_completed_by_type=self._delta(now["completed_by_type"], before["completed_by_type"]),
            result=rec["results"].get(job_type, {}),
            errors=rec["errors"],
            sample_id=sample,
        )

    def finished_count(self, job_type: str, timeout: float = 0.0) -> int:
        """Completed + failed jobs of this type so far; with ``timeout``, waits for the first."""
        t0 = time.time()
        while True:
            c = self._counters()
            n = c["completed_by_type"].get(job_type, 0) + c["failed_by_type"].get(job_type, 0)
            if n or time.time() - t0 >= timeout:
                return n
            time.sleep(0.2)

    def close(self) -> List[Dict[str, str]]:
        """Shut the Coordinator down (this writes failed_jobs.csv) and return its rows."""
        if not self._closed:
            self._closed = True
            try:
                ray.get(self.coord.shutdown.remote(), timeout=120)
            finally:
                for name in ["robin_coordinator"] + [f"pool_{p}" for p in set(POOL_OF.values())]:
                    try:
                        ray.kill(ray.get_actor(name))
                    except Exception:  # noqa: BLE001 - already gone
                        pass
        path = os.path.join(self.work_dir, "failed_jobs.csv")
        if not os.path.exists(path):
            return []
        with open(path, newline="", encoding="utf-8") as fh:
            return list(csv.DictReader(fh))


def init_ray() -> None:
    if not ray.is_initialized():
        ray.init(num_cpus=4, include_dashboard=False, log_to_driver=False)
    # Ship the synthetic handler by value: the Ray workers need not import this tests module.
    import sys

    from ray import cloudpickle

    cloudpickle.register_pickle_by_value(sys.modules[__name__])


def shutdown_ray() -> None:
    if ray.is_initialized():
        ray.shutdown()
