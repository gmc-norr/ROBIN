---
title: A failure another job finds goes in the BAM ledger under the key saved at staging
author: Pär
date: 2026-11-02
branch: fix/batch-10e-recorded-items
cost: { minutes: 1, tokens: 1638 }
id: l-5ee67801
scope: repo
status: active
src: agent
touches:
  - src/robin/file_ledger.py
  - src/robin/analysis/target_analysis.py
  - reviews/plans/2026-11-02-batch-10e1-round2.md
---

## Decided
- A BAM that another job's add-up leaves out is written into the work folder's BAM ledger by the add-up itself (a `file_failed` line), not by the job that ran: a target job, a finalize job, an SNP run's merge, the end-of-run finalize and the GUI's finalize all add up, and only jobs run through the workers reach `record_step_end`.
- The line is keyed by the BAM's key saved at staging (`bam_key_<n>.txt`, an optional part of the set, like `union_<n>.parquet`). Reading the key from the staged path when the failure is found fails once the BAM is renamed, and the failure was then only logged (the plan check's case). The key is read at staging from the BAM's bytes; passing the job's key down was not chosen, because no test could tell the two apart.
- A later `end` keeps an earlier `file_failed`; only a `start` resets the step. A step with only a `file_failed` is reported at the next start even though it never ended.
- A staged set counts as taken once its first file has moved (part 1 marked it taken just before the move). A claim that took nothing removes its empty folder.

## Watch out
- The ledger keeps a line for any key (`_apply` uses `setdefault`). A failure line for a BAM the ledger never took in would make that BAM known, and the next start would skip its target step, which never ran. So the line is written only when the ledger has taken the key in.
- That key check also stops a line in a folder with no ledger, so removing the separate no-ledger check alone changes nothing a test can see (mutant P19 needed both edits).
- Running the whole suite from a copy of the tests in another folder gives false failures: some tests open source files by path, and some error texts get cut short by the longer path. Run on the branch instead, leaving out test files that belong to a later commit.

In the code: `record_bam_failure`, `FileLedger._apply` and `BamIntake._resume` in `src/robin/file_ledger.py`; `process_file_with_staging`, `accumulate_staged_files`, `_staged_set_key` and `_record_left_out` in `src/robin/analysis/target_analysis.py`. Tests: `tests/test_target_failed_files.py` (Y1, Y3, Y6, Y7), `tests/test_bam_ledger.py` (Y4, Y5). Plan: `reviews/plans/2026-11-02-batch-10e1-round2.md`.
