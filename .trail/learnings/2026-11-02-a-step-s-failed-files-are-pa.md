---
title: "A step's failed files are found by an owner tag, not a count; failure reports wait until the Coordinator has them"
author: Pär
date: 2026-11-02
branch: fix/batch-10d-file-ledger
cost: { minutes: 24, tokens: 135114 }
id: l-ec95893e
scope: repo
status: active
src: agent
touches:
  - src/robin/analysis/master_csv_manager.py
  - src/robin/file_ledger.py
  - src/robin/workflow_ray.py
  - src/robin/workflow_simple.py
  - tests/test_bam_ledger.py
  - tests/test_master_csv_once.py
  - tests/test_watcher_sends.py
  - reviews/plans/2026-11-02-batch-10d-round3.md
  - CHANGELOG.md
  - reviews/2026-10-16-code-review.md
---

## Decided
- A step's failed files are found by an owner tag, not by a count. Both runners set the running step (`enter_step`/`leave_step` in `robin.file_ledger`) around the handler. `WorkflowContext.add_metadata("file_failures", …)` tags each new entry with it, and the ledger records only the job's own entries. A count taken before the handler broke: the BED step replaces the list, and in the non-Ray runner sibling steps share it at once.
- Failed-file rows the watcher sends wait in a queue until the Coordinator takes them. The Coordinator ignores a row whose kind, file and step it already holds, so sending one again is safe.
- `update_master_csv` stops (raises) when `master_csv_bams.txt` lists BAMs but `master.csv` is missing or has an empty `last_bam_added`: the listed BAMs' counts are gone, and counting from zero would hide it.

## Tried and ruled out
- Importing the ContextVar itself into `workflow_ray`: Ray pickles the Coordinator class with the module names its code uses, and a ContextVar can't be pickled. Every Ray actor test failed with "Could not serialize the actor class". Keep the ContextVar inside `file_ledger`, behind functions.

## Watch out
- Mark something "reported" only after the report got there (case 15, broken again in part 3: case 16).
- For a shared list or file, check every writer and how it writes: append, replace or delete.

In the code: `CURRENT_STEP`, `enter_step`, `leave_step`, `tag_file_failures` and `record_step_end` in `src/robin/file_ledger.py`; `add_metadata` in both `WorkflowContext`s; `RayFileWatcher._report` and `_deliver_rows`, and `Coordinator.record_cut_short_steps`, in `workflow_ray.py`; `update_master_csv` in `master_csv_manager.py`.
