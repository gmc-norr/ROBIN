---
title: "A failure a coverage commit finds goes in its install list, so a finish after a stop still reports it"
author: Pär
date: 2026-11-02
branch: fix/batch-10e-recorded-items
cost: { minutes: 25, tokens: 202182 }
id: l-2fd0c4a0
scope: repo
status: active
src: agent
touches:
  - src/robin/analysis/target_analysis.py
---

## Decided
- When a commit's batch BAM can't be written, the install list gets a `keep_claim` step: the problem text, and each of the batch's BAMs with its ledger key. Doing the step reports those BAMs as failed files and ledger lines. So they are reported whoever finishes the list: the worker itself, or a later add-up or finalize after a stop.
- Finalize adds the claim folder its own finish kept (a lost batch BAM) to `claims_left`. The SNP check already refuses on `claims_left`, so no SNP calls are made on a `target.bam` without those reads.
- The claim reads each set's BAM path and key from `_staging/` before it moves any of the set's files. A move that fails part-way still names the BAM and writes its ledger line.

## Tried and ruled out
- Reporting the failed write after the install, from the worker's memory: a stop between writing the list and that report lost it for good (found by an external code review).
- Recording the failure before the list is written: a commit that then failed before its list would say "its coverage is counted", which is false.
- Making the SNP check read `file_failures`: a finish's failed-file entries come once; a kept claim folder is listed by every later finalize until someone deals with it.
- Listing every claim folder after the finish: that could list a new worker's claim taken after finalize's wait.

## Watch out
- Anything a commit must report has to be a step in the list, or follow from one. The worker's memory dies with a stop.
- A new step kind needs `_claim_of` too, or its kept claim folder and `_staging/` are not synced before the list is deleted.

In the code: `src/robin/analysis/target_analysis.py`: `_install` (the `keep_claim` step), `_claim_of`, `_finish_commit` (returns `(entries, kept)`), `finalize_accumulation_for_sample` (`left += kept`), and the claim loop's `who`. Tests C9 and C10 in `tests/test_coverage_commit.py`, Y8 in `tests/test_target_failed_files.py`. Plan: `reviews/plans/2026-11-03-batch-10e-round3.md`.
