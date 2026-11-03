---
title: "Finalize finds a kept claim folder by its incomplete.txt, not from the job that finished the list"
author: Pär
date: 2026-11-03
branch: fix/batch-10e-recorded-items
cost: { minutes: 34, tokens: 165698 }
id: l-b078116f
scope: repo
status: active
src: agent
touches:
  - src/robin/analysis/target_analysis.py
---

This replaces round 3's way of listing a kept claim folder (learning l-2fd0c4a0, its second "Decided" point).

## Decided
- Any job can finish a coverage commit's install list: the worker, the next add-up, or any finalize. So a fact every later finalize needs is kept on disk where every finalize looks. A kept claim folder holds `incomplete.txt`. Finalize lists every `_staging/_batch_*` that holds it, read under the accumulation lock after its own finish, plus the claim folders it saw in its wait that are still there.
- A `keep_claim` step makes its claim folder again, with `incomplete.txt`, when that file is not there. The lost-batch-BAM step does the same.
- Before every install list, the claim folder and `_staging/` are synced, with or without a batch BAM.
- An error with no text is named by its type (`str(e) or type(e).__name__`).

## Tried and ruled out
- Round 3's `kept`, returned by the finish to its caller: when another finalize or an add-up finished the list after finalize's wait, finalize got nothing back, and the SNP calls ran on missing reads (found by an external code review).
- A new marker file in the sample folder: new saved state, and `incomplete.txt` already is that marker.
- Listing every claim folder after the finish: it would list a new worker's claim still in progress, which has no `incomplete.txt`.

## Watch out
- `str(e)` can be empty (`RuntimeError()`). Don't test an error's text for truth.
- A rule claimed "by reading" about every error or every caller needs a test of the odd one.

In the code: `src/robin/analysis/target_analysis.py`: `_kept_claims`, `finalize_accumulation_for_sample`, `_install` (the `keep_claim` step), and the second locked part of `accumulate_staged_files`. Tests C11-C14 in `tests/test_coverage_commit.py`. Plan: `reviews/plans/2026-11-03-batch-10e-round4.md`. Case 17 in `.claude/judgment-gate-cases.md`.
