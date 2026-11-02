"""Which BAMs ROBIN has taken in, and which of their steps started and ended (JOB-2).

One append-only file per work folder, ``robin_bam_ledger.jsonl``, one JSON object per line:

- ``{"e": "take", "key": K, "file": path}``: the intake took the BAM in.
- ``{"e": "start", "key": K, "step": S, "sample": id, "file": path}``: a worker is about to
  run step S for the BAM (written before the handler runs).
- ``{"e": "end", "key": K, "step": S, "status": st, "outcome": o}``: step S ended for the BAM
  (after the handler). ``o`` is the job's outcome (success, failed or expected_failure), with
  ``"error"`` when it failed; ``"file_failed"`` when this BAM failed inside a job that went
  on; ``"downstream"`` for preprocessing: whether its steps follow.

K is the BAM's key (``robin.bam_identity.bam_key``): a BAM is known by its bytes. Only the
per-BAM steps are recorded (preprocessing and the five it starts); the per-sample ones
(classifiers, IGV BAM, SNP, finalize) are not.

After a restart, for each BAM: a step that ended is skipped (success or not), a step that
never started runs, and a step that started but never ended (it was running when ROBIN
stopped) is not run again, because it may have added its data. Preprocessing is the
exception: it counts a BAM once (batch 10d-1), so a BAM whose preprocessing never ended goes
through it again, and a restarted BAM is sent through its preprocessing told which steps to
start (``only_steps``). At every start, each step that failed or may be incomplete for the
BAM is reported as a failed file (none is run again).
"""

import fcntl
import json
import logging
import os
import threading
import time
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

from robin.bam_identity import bam_key

LEDGER_NAME = "robin_bam_ledger.jsonl"
PER_BAM_STEPS = ("preprocessing", "bed_conversion", "mgmt", "cnv", "target", "fusion")
AFTER_PREPROCESSING = PER_BAM_STEPS[1:]
# The 28-byte BGZF end-of-file block htslib writes when it closes a BAM
_BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")

_log = logging.getLogger("robin.file_ledger")


class LedgerError(RuntimeError):
    """The ledger can't be read or written: ROBIN can't tell which BAMs it has done."""


def append_whole_lines(path: str, data: bytes) -> None:
    """Append ``data`` (whole lines) to ``path`` under an exclusive flock, and fsync.

    A torn last line, left by a stop part-way through an earlier append, is cut off first:
    otherwise the new lines would join it, and the file could no longer be read.
    """
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "a+b") as f:
        fcntl.flock(f.fileno(), fcntl.LOCK_EX)
        try:
            end = f.seek(0, os.SEEK_END)
            keep = _whole_lines_end(f, end)
            if keep < end:
                _log.warning(
                    f"Cutting off an unfinished last line of {path} ({end - keep} bytes, left by a stop)"
                )
                f.truncate(keep)
            f.write(data)
            f.flush()
            os.fsync(f.fileno())
        finally:
            fcntl.flock(f.fileno(), fcntl.LOCK_UN)


def _whole_lines_end(f, end: int) -> int:
    """Where the file's whole lines end: just after its last newline (0 if it has none)."""
    if end == 0:
        return 0
    f.seek(end - 1)
    if f.read(1) == b"\n":
        return end
    pos = end
    while pos > 0:
        size = min(65536, pos)
        f.seek(pos - size)
        i = f.read(size).rfind(b"\n")
        if i >= 0:
            return pos - size + i + 1
        pos -= size
    return 0


def bam_is_complete(path: str) -> bool:
    """True when the BAM ends with the BGZF end-of-file block: its writer has closed it."""
    try:
        with open(path, "rb") as f:
            f.seek(0, os.SEEK_END)
            size = f.tell()
            if size < len(_BGZF_EOF):
                return False
            f.seek(size - len(_BGZF_EOF))
            return f.read() == _BGZF_EOF
    except OSError:
        return False


@dataclass
class StepState:
    started: bool = False
    ended: bool = False
    status: Optional[str] = None
    downstream: Optional[bool] = None
    outcome: Optional[str] = None  # the job's outcome: success, failed or expected_failure
    error: str = ""
    file_failed: Optional[str] = None  # this BAM failed inside a job that went on


@dataclass
class BamState:
    file: str = ""
    sample: str = ""
    steps: Dict[str, StepState] = field(default_factory=dict)


class FileLedger:
    """The ledger file of one work folder."""

    def __init__(self, work_dir: str):
        self.path = os.path.join(work_dir, LEDGER_NAME)

    def load(self) -> Dict[str, BamState]:
        """Each BAM's state, by key. Text after the last newline is a record a stop cut off and
        is ignored; any other line that can't be read raises LedgerError."""
        state: Dict[str, BamState] = {}
        if not os.path.exists(self.path):
            return state
        try:
            with open(self.path, "rb") as f:
                data = f.read()
        except OSError as e:
            raise LedgerError(f"Could not read the BAM ledger {self.path}: {e}") from e
        lines = data.split(b"\n")
        torn = lines.pop()
        for n, raw in enumerate(lines, 1):
            if not raw.strip():
                continue
            try:
                rec = json.loads(raw.decode("utf-8"))
                self._apply(state, rec)
            except (ValueError, KeyError, TypeError, AttributeError) as e:
                raise LedgerError(
                    f"Could not read line {n} of the BAM ledger {self.path} ({e}). ROBIN can't tell "
                    "which BAMs it has done; move the file away to start this work folder's "
                    "record again (BAMs done before may then be counted again)"
                ) from e
        if torn.strip():
            _log.warning(f"Ignoring an unfinished last line of {self.path} (cut off by a stop)")
        return state

    @staticmethod
    def _apply(state: Dict[str, BamState], rec: dict) -> None:
        kind, key = rec["e"], rec["key"]
        bam = state.setdefault(key, BamState())
        if kind == "take":
            bam.file = rec.get("file") or bam.file
        elif kind == "start":
            bam.steps[rec["step"]] = StepState(started=True)
            bam.sample = rec.get("sample") or bam.sample
            bam.file = bam.file or rec.get("file") or ""
        elif kind == "end":
            step = bam.steps.setdefault(rec["step"], StepState())
            step.ended, step.status = True, rec.get("status")
            step.downstream = rec.get("downstream")
            step.outcome = rec.get("outcome")
            step.error = str(rec.get("error") or "")
            step.file_failed = rec.get("file_failed")
        else:
            raise ValueError(f"unknown record kind {kind!r}")

    def append(self, records: List[dict]) -> None:
        """Append the records in one write, under an exclusive flock, and fsync (a torn last
        line is cut off first)."""
        stamp = time.time()
        data = "".join(
            json.dumps({"v": 1, "t": stamp, **r}, sort_keys=True) + "\n" for r in records
        ).encode("utf-8")
        append_whole_lines(self.path, data)


@dataclass
class Decision:
    """What the intake decided for one BAM."""

    run: bool
    key: Optional[str] = None
    only_steps: Optional[List[str]] = None  # None: every step of the plan
    cut_short: List[Dict[str, str]] = field(default_factory=list)  # rows to report as failed files
    reason: str = ""


class BamIntake:
    """Takes each complete BAM in once (JOB-2). One per session, shared by the start-up scan,
    the watcher and the GUI's "add folder"; ``decide`` is thread-safe."""

    def __init__(self, work_dir: Optional[str], plan: Iterable[str]):
        types = [str(s).split(":", 1)[-1] for s in plan]
        self.after = [t for t in AFTER_PREPROCESSING if t in types]
        self.ledger = FileLedger(work_dir) if work_dir else None
        self.state: Dict[str, BamState] = self.ledger.load() if self.ledger else {}
        self._seen: set = set()
        # Half-written BAMs: path -> [(size, mtime), when that last changed, reported]
        self._waiting: Dict[str, list] = {}
        self._lock = threading.Lock()

    def decide(self, path: str) -> Decision:
        """Whether to send the BAM at ``path``, and which of its steps to start.

        Not yet: a BAM still being written (no end-of-file block). Not again: a BAM this
        session already decided. A new BAM is taken in (all steps). A BAM the ledger knows
        gets only its unstarted steps; its cut-short steps are recorded as ended and returned
        to report.
        """
        if not bam_is_complete(path):
            sig, now = _file_sig(path), time.monotonic()
            with self._lock:
                held = self._waiting.get(path)
                if held is None:
                    self._waiting[path] = [sig, now, False]
                    _log.info(f"Waiting for {path} to be complete (no BAM end-of-file block yet)")
                elif held[0] != sig:
                    held[0], held[1] = sig, now  # still being written
            return Decision(False, reason="incomplete")
        key = bam_key(path)
        with self._lock:
            self._waiting.pop(path, None)
            if key in self._seen:
                return Decision(False, key, reason="seen")
            bam = self.state.get(key)
            if bam is None:
                self._append([{"e": "take", "key": key, "file": path}])
                self.state[key] = BamState(file=path)
                self._seen.add(key)
                return Decision(True, key, None, reason="new")
            decision = self._resume(key, bam, path)
            self._seen.add(key)
            return decision

    def _resume(self, key: str, bam: BamState, path: str) -> Decision:
        pre = bam.steps.get("preprocessing")
        to_run = [s for s in self.after if s not in bam.steps]
        cut = [s for s in AFTER_PREPROCESSING if s in bam.steps and bam.steps[s].started and not bam.steps[s].ended]
        if pre is not None and pre.ended and not (pre.status == "success" and pre.downstream):
            to_run = []  # preprocessing failed or skipped the BAM: no steps follow
        if cut:
            self._append([{"e": "end", "key": key, "step": s, "status": "cut_short"} for s in cut])
            for s in cut:
                bam.steps[s].ended, bam.steps[s].status = True, "cut_short"
        # Every start shows each step that failed or may be incomplete (none is run again)
        rows = []
        for s in PER_BAM_STEPS:
            step = bam.steps.get(s)
            found = _earlier_failure(s, step) if step is not None and step.ended else None
            if found is None:
                continue
            kind, text = found
            rows.append({
                "sample_id": bam.sample or "unknown",
                "filepath": path,
                "job_type": s,
                "error": text,
                "failure_kind": kind,
            })
            _log.error(f"{os.path.basename(path)}: {text}")
        run = bool(to_run) or pre is None or not pre.ended
        return Decision(run, key, list(to_run) if run else None, rows, reason="resume")

    def waiting_paths(self) -> List[str]:
        """The half-written BAMs seen so far, to check again (JOB-4)."""
        with self._lock:
            return list(self._waiting)

    def incomplete_too_long(self, seconds: float) -> List[str]:
        """Half-written BAMs whose size and time have not changed for ``seconds``, each listed
        once: they are marked reported."""
        now = time.monotonic()
        found = []
        with self._lock:
            for path, held in self._waiting.items():
                if not held[2] and now - held[1] >= seconds:
                    held[2] = True
                    found.append(path)
        return found

    def _append(self, records: List[dict]) -> None:
        if self.ledger is not None:
            self.ledger.append(records)


def _file_sig(path: str):
    """A file's size and modification time (None if it can't be read): a change means it is
    still being written."""
    try:
        st = os.stat(path)
        return st.st_size, st.st_mtime_ns
    except OSError:
        return None


def _earlier_failure(step: str, state: StepState) -> Optional[Tuple[str, str]]:
    """(failure_kind, text) for an ended step that was cut short or failed for this BAM; None
    otherwise (a fail-only BAM's expected failure is not a failure)."""
    if state.status == "cut_short":
        return "cut_short", (
            f"{step} was running for this BAM when ROBIN stopped: it may be incomplete for this "
            "BAM, and it is not run again (it may have added its data)"
        )
    if state.outcome == "failed":
        return "failed_before", (
            f"{step} failed for this BAM in an earlier session ({state.error or 'no error text'}); "
            "it is not run again"
        )
    if state.file_failed:
        return "failed_before", (
            f"{step} failed for this BAM in an earlier session, while its batch went on "
            f"({state.file_failed}); it is not run again"
        )
    return None


def _job_bams(job) -> List[Tuple[str, object]]:
    """(key, context) for each BAM the job works on: a batched job's contexts, or its own."""
    metadata = getattr(job.context, "metadata", None) or {}
    batched = metadata.get("_batched_job")
    contexts = list(getattr(batched, "contexts", None) or [job.context])
    found = []
    for c in contexts:
        key = (getattr(c, "metadata", None) or {}).get("bam_key")
        if key:
            found.append((key, c))
    return found


def _sample(context) -> str:
    bm = (getattr(context, "metadata", None) or {}).get("bam_metadata") or {}
    return str(bm.get("sample_id") or "") if isinstance(bm, dict) else ""


def _ledger_of(job) -> Optional[FileLedger]:
    work_dir = (getattr(job.context, "metadata", None) or {}).get("work_dir")
    return FileLedger(work_dir) if work_dir else None


def record_step_start(job, job_type: str) -> None:
    """Record that ``job_type`` is starting for each of the job's BAMs. Raises LedgerError
    when it can't: the step must not run unrecorded (a restart would run it again)."""
    if job_type not in PER_BAM_STEPS:
        return
    bams, ledger = _job_bams(job), _ledger_of(job)
    if not bams or ledger is None:
        return
    records = [
        {"e": "start", "key": k, "step": job_type, "sample": _sample(c), "file": getattr(c, "filepath", "")}
        for k, c in bams
    ]
    try:
        ledger.append(records)
    except Exception as e:
        raise LedgerError(
            f"Could not record the start of {job_type} in the BAM ledger ({e}); the step was not run"
        ) from e


def file_failure_count(job) -> int:
    """How many files the job's context lists as failed so far (taken before the handler, so
    that only the files this job adds count: sibling jobs can share a context)."""
    return len((getattr(job.context, "metadata", None) or {}).get("file_failures") or [])


def record_step_end(
    job, job_type: str, outcome: Optional[dict] = None, failures_before: int = 0
) -> None:
    """Record that ``job_type`` ended for each of the job's BAMs: the result's status, the
    job's outcome (``job_outcome``: success, failed or expected_failure) with its error, and
    for each BAM among the files the job reported as failed while it went on, that failure.
    A failure to write is logged: after a restart the step then shows as cut short."""
    if job_type not in PER_BAM_STEPS:
        return
    bams, ledger = _job_bams(job), _ledger_of(job)
    if not bams or ledger is None:
        return
    outcome = outcome or {"status": "failed"}
    result = (getattr(job.context, "results", None) or {}).get(job_type)
    status = (
        result.get("status") if isinstance(result, dict) and result.get("status")
        else outcome.get("status") or "failed"
    )
    extra: Dict[str, object] = {"outcome": outcome.get("status") or "failed"}
    if extra["outcome"] == "failed":
        extra["error"] = str(outcome.get("error") or "")[:500]
    if job_type == "preprocessing":
        skipped = bool((job.context.metadata or {}).get("skip_downstream"))
        extra["downstream"] = status == "success" and not skipped
    failed_files: Dict[str, str] = {}
    added = ((job.context.metadata or {}).get("file_failures") or [])[failures_before:]
    for f in added:
        if isinstance(f, dict) and f.get("file"):
            failed_files[os.path.realpath(str(f["file"]))] = str(f.get("error") or "")[:500]
    records = []
    for k, c in bams:
        rec = {"e": "end", "key": k, "step": job_type, "status": status, **extra}
        text = failed_files.get(os.path.realpath(str(getattr(c, "filepath", "") or "")))
        if text is not None:
            rec["file_failed"] = text or "no error text"
        records.append(rec)
    try:
        ledger.append(records)
    except Exception as e:
        _log.error(f"Could not record the end of {job_type} in the BAM ledger: {e}")
