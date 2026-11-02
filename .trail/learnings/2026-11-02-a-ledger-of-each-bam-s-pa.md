---
title: "A ledger of each BAM's step starts and ends: after a restart skip finished, run unstarted, report cut-short"
author: Pär
date: 2026-11-02
branch: fix/batch-10d-file-ledger
cost: { minutes: 42, tokens: 267220 }
id: l-294c2aaa
scope: repo
status: active
src: agent
touches:
  - src/robin/file_ledger.py
  - src/robin/workflow_ray.py
  - src/robin/workflow_simple.py
  - tests/test_bam_ledger.py
---

## Decided
- The BAM ledger (`<work_dir>/robin_bam_ledger.jsonl`, append-only JSON lines under flock) records "take", and each per-BAM step's "start" and "end". The worker writes start right before the handler and end right after (in `_impl` for Ray, `WorkflowManager.worker` for non-Ray), so the record sits next to the work; the Coordinator's notifications are async and could lag a crash.
- After a restart: ended steps are skipped, unstarted ones run (the BAM goes back through preprocessing with `only_steps`; `master.csv` counts a BAM once, batch 10d-1), and a step started but not ended is not run again (it may have added its data) and is reported once as a failed file, `failure_kind` "cut_short" (its end is then written as "cut_short").
- One `BamIntake` per session, shared by the scan, the watcher and `add_watch_path` (Ray: built in `run()`, kept in `_GLOBAL_WATCH_CONTEXT["intake"]`). Two intakes in one process would each think a just-taken BAM was never started and send it twice.
- A BAM is taken in only once it ends with the 28-byte BGZF end-of-file block; before that it is not marked processed, so its next watch event decides it again. The content key only means something on a complete file.

## Watch out
- If the start record can't be written the step must not run (it raises): a step that ran unrecorded would run again after a restart. A missing end only makes it show as cut short.
- Tests that feed the watcher placeholder `.bam` files must make them complete and different (name + the EOF block), or the intake takes in none, or only one.
- A torn last ledger line (a stop mid-append) is ignored; that is safe because each record is written before the thing it announces happens. Any other bad line stops the intake.

In the code: `src/robin/file_ledger.py`; `workflow_ray.py` (`_impl`, `_on_finish` `only_steps`, `record_cut_short_steps`, `_jobs_for_file`, `submit_existing_paths`, `RayFileWatcher._handle`, `run`); `workflow_simple.py` (`worker`, `can_trigger_jobs`, `FileWatcher.handle_file`, `start`). Tests: `tests/test_bam_ledger.py`.
