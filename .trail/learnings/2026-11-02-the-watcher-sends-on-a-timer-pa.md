---
title: "The watcher sends on a timer, puts back a failed send, and stops before the Coordinator; half-written BAMs are checked again"
author: Pär
date: 2026-11-02
branch: fix/batch-10d-file-ledger
cost: { minutes: 17, tokens: 106832 }
id: l-18297593
scope: repo
status: active
src: agent
touches:
  - src/robin/file_ledger.py
  - src/robin/workflow_ray.py
  - tests/test_watcher_sends.py
  - reviews/plans/2026-11-02-batch-10d3-watcher-sends.md
  - CHANGELOG.md
  - reviews/2026-10-16-code-review.md
  - reviews/2026-10-16-code-review.pdf
  - reviews/plans/2026-11-02-batch-10d2-bam-ledger.md
---

## Decided
- The Ray watcher has a timer thread (`start`/`stop`/`tick`, every 0.5 s). It sends what waits, and checks half-written BAMs again. Before, sending happened only inside a file event, so the last BAMs of a run could wait for good, and a BAM half-written at the start-up scan was never looked at again.
- A failed send puts the batch back and is tried again. Only each BAM's first failure is logged.
- At shutdown, `_stop_watching` (observer, then watcher) runs before `coord.shutdown`, on both paths. Once told to stop, the Coordinator's `submit_jobs` returns silently. What still can't be sent is logged by name. The BAM ledger runs it at the next start: it was taken in, but no step started.
- A half-written BAM unchanged (same size and mtime) for 10 minutes is shown once as a failed file, `failure_kind` "incomplete". It is still taken in if it completes later.

## Watch out
- File events and the timer both call `_handle` and touch the buffer: one `RLock`. A tick can block in `ray.get` (the `high` preset's Coordinator answers nothing while a job runs), so `stop()` joins the timer with a 60 s limit.
- `flush_remaining()` makes one pass. A loop until the buffer is empty would spin for good now that a failed batch goes back.
- To check a call's order in `run()` from its source, use its full dotted name (`ast.unparse`). `coord.shutdown.remote()` is a call named `remote`.

In the code: `RayFileWatcher` (`start`, `tick`, `stop`, `_send_batch`, `flush_remaining`) and `_stop_watching` in `src/robin/workflow_ray.py`; `BamIntake.waiting_paths` and `incomplete_too_long` in `src/robin/file_ledger.py`. Tests: `tests/test_watcher_sends.py`.
