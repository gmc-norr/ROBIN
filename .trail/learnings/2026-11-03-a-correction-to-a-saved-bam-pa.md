---
title: A correction to a saved BAM failure is written after the install is on disk and before its list is deleted
author: Pär
date: 2026-11-03
branch: fix/batch-10e-recorded-items
cost: { minutes: 31, tokens: 607427 }
id: l-c65e7566
scope: repo
status: active
src: agent
touches:
  - src/robin/analysis/target_analysis.py
---

## Decided
- A failure the coverage commit records in the BAM ledger ("its reads are not in target.bam") names the commit's claim folder. When a later finish of the same install list puts the batch BAM in place, it writes a `cleared` line for that claim, and the reader drops every failure naming it, including ones that arrive later (the BAM's own end line can come after the finish).
- The clear goes in one exact place: after the install's syncs, before the list is deleted. Before the syncs, a power loss can keep the clear and lose the move, so a failure that holds again is hidden. After the list is deleted, a stop loses the clear for good, because no list is left to finish.
- A clear that can't be written raises, so the list stays and the next add-up or finalize tries again. No ledger, or a BAM the ledger never took in, is a quiet skip, not an error.
- A `cleared` line never adds a step or a BAM to the reader's state. The intake runs every step a BAM has no record of, so an added empty step would silently stop that step from running.

## Tried and ruled out
- Writing the clear inside the install step (the first draft): Codex's plan check found it could reach the disk before the batch BAM's move.
- A clear that only logs when it fails: the list was then deleted, and the stale failure stayed for good.
- One failure slot per step: an older unrelated failure, overwritten by the claim's failure, was lost when the claim was cleared.
- Matching the failure by its text: texts are cut to 500 characters.

## Watch out
- Finalize swaps in the new target.bam and deletes the batch BAMs without syncing either. That is older than this fix, puts reads at risk on a power loss, and is recorded for later. The clear rests on the same evidence as the list's deletion, so it is no weaker than that.
- Ledgers written before this fix have failures with no claim. They are never cleared.

In the code: `src/robin/file_ledger.py` (`_add_failure`, `FileLedger._apply`, `record_step_end`, `record_bam_failure(claim=)`, `clear_bam_failures`); `src/robin/analysis/target_analysis.py` (`_finish_commit`: `_install`, `_put_on_disk`, `_clear_left_out`, `_delete_list`). Tests: `tests/test_coverage_commit.py` C16–C18, `tests/test_bam_ledger.py` L1–L7. Plan: `reviews/plans/2026-11-03-batch-10e-round6.md`.
