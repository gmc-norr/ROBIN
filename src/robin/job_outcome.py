"""One rule for a job's outcome, shared by the Ray wrapper, both Ray judges and the
non-Ray WorkflowManager.

A job fails when its handler raises, or when the handler filed an error into the
context while it ran, whatever name the error was filed under. Only errors filed
*during* the job count: the caller notes ``len(context.errors)`` before the handler
runs and passes it as ``errors_before``, so an error an earlier job left in the same
context cannot fail a later job.

Fail-only BAM submissions (the downstream tools may not support them) turn a failure of
a classifier or ``bed_conversion`` job into an expected, not-failed outcome; the rule
is the same for an exception and for a filed error.

The outcome is stored in ``context.metadata["job_outcomes"]["<job_id>:<job_type>"]``.
``Job.next_job()`` reuses the job id and the context for the next step, so the type is
part of the key: a context carried through several jobs never hands job A's outcome to
job B. It travels back with the context from the Ray worker to the judges.
"""

from __future__ import annotations

import os
from typing import Any, Dict, Optional, Set

# Classification job types (single global pipeline per type)
CLASSIFICATION_TYPES: Set[str] = {"sturgeon", "nanodx", "pannanodx", "random_forest"}

# Job types whose failure on a fail-only BAM submission is expected.
FAIL_ONLY_TYPES: Set[str] = CLASSIFICATION_TYPES | {"bed_conversion"}

OUTCOMES_KEY = "job_outcomes"


class JobReportedFailure(Exception):
    """Raised by a runner after job_outcome() judged the job failed by a filed error."""


def bam_filename_kind(filepath: str) -> str:
    """Classify BAM filenames by substring: 'pass' | 'fail' | 'unknown'."""
    try:
        n = os.path.basename(filepath).lower()
    except Exception:
        return "unknown"
    if "pass" in n:
        return "pass"
    if "fail" in n:
        return "fail"
    return "unknown"


def is_expected_fail_only_failure(job_type: str, context: Any) -> bool:
    """True for a fail-only BAM submission's classifier / bed_conversion job.

    If *any* BAM in the submission contains 'pass', errors must stay visible.
    """
    try:
        if job_type in FAIL_ONLY_TYPES and bool(
            context.metadata.get("fail_only_bam_submission", False)
        ):
            return bam_filename_kind(getattr(context, "filepath", "") or "") == "fail"
    except Exception:
        pass
    return False


def job_outcome(
    job_type: str,
    context: Any,
    errors_before: int,
    exc: Optional[BaseException] = None,
) -> Dict[str, Any]:
    """Decide a job's outcome: {"status": success|failed|expected_failure, "ok", ...}.

    ``errors_before`` is ``len(context.errors)`` taken before the handler ran.
    """
    new_errors = list(getattr(context, "errors", None) or [])[errors_before:]
    if exc is None and not new_errors:
        return {"status": "success", "ok": True}
    if exc is not None:
        message = str(exc)
    else:
        last = new_errors[-1]
        message = (
            last.get("error") if isinstance(last, dict) else str(last)
        ) or "handler reported error"
    if is_expected_fail_only_failure(job_type, context):
        return {
            "status": "expected_failure",
            "ok": True,
            "reason": "fail_only_bam_submission",
            "message": message,
        }
    return {"status": "failed", "ok": False, "error": message}


def _key(job_id: Any, job_type: str) -> str:
    return f"{job_id}:{job_type}"


def record_outcome(context: Any, job_id: Any, job_type: str, outcome: Dict[str, Any]) -> None:
    """Store the outcome in the context, keyed by job id and type."""
    context.metadata.setdefault(OUTCOMES_KEY, {})[_key(job_id, job_type)] = outcome


def recorded_failure(context: Any, job_id: Any, job_type: str) -> Optional[str]:
    """The failure text if the stored outcome for this job is 'failed', else None."""
    try:
        outcome = (context.metadata.get(OUTCOMES_KEY) or {}).get(_key(job_id, job_type))
    except Exception:
        return None
    if isinstance(outcome, dict) and outcome.get("status") == "failed":
        return outcome.get("error") or "handler reported error"
    return None
