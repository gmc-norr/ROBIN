---
title: A target coverage commit is an install list finished by the next add-up or finalize
author: Pär
date: 2026-11-02
branch: fix/batch-10e-recorded-items
cost: { minutes: 24, tokens: 126887 }
id: l-41e35d4f
scope: repo
status: active
src: agent
touches:
  - src/robin/analysis/target_analysis.py
  - tests/test_coverage_commit.py
  - reviews/plans/2026-11-02-batch-10e2-coverage-commit.md
---

## Decided
- A target coverage commit is an install list, not a rewrite in place: every file is prepared in `<sample>/_commit/`, then `install.json` lists the steps (replace, remove, append, batch_bam, drop_claim), then they run. Once the list is written the batch counts (part 1's `committed`). The next accumulation, or finalize, finishes a list it finds; each step is safe to do twice. Not chosen: a generation folder with a pointer (every reader would follow it), backup and roll back (copies every file and must put claimed sets back), rewriting the time series whole (it grows with the run).
- The time series is appended by cutting the file back to the size recorded in the list, then adding the rows: a half-written append is put right by doing it again.
- A batch BAM gone from its claim folder counts as installed only if it is in the sample folder or `target.bam`'s header lists it as merged: finalize merges and deletes batch BAMs under its own lock, so it can do so between a commit's batch_bam step and the list's deletion. Otherwise its BAMs are reported, and the claim folder is kept.
- Finalize finishes a list after it waits for claims (`_claims_left`), not before: a worker that dies during the wait leaves a list, and its counted claim folder would otherwise be listed.

## Watch out
- Syncing a file's contents does not put its name on disk. Before the list is written: the sample folder after `_commit/` is made, `_commit/` after its files, the batch BAM and index, then the claim folder and `_staging/`. After the list is deleted, `_commit/` again: a list a power loss brought back would be done after a later commit, and its append would cut off that commit's rows.
- A crash test that only checks after a kill misses a finish that runs through without one: mutant Q13 (a finish that deletes `_commit/` without doing its steps) passed the locked C5 until C5 also checked the unkilled finish.
- At the default call threshold (0.1x) every target in the test fixture is over it after one BAM, so the threshold files never change between commits and a test of them is blind; `tests/test_coverage_commit.py` uses 1.5x.

In the code: `_commit_batch`, `TargetAnalysis._finish_commit`, `accumulate_staged_files`, `finalize_accumulation_for_sample`, and `_install`, `_install_batch_bam`, `_append_rows`, `_retire`, `_sync_file`, `_sync_dir` in `src/robin/analysis/target_analysis.py`. Tests: `tests/test_coverage_commit.py` (C0-C8), `tests/test_panel_depth.py` (R1, its guard, R1 legacy). Plan: `reviews/plans/2026-11-02-batch-10e2-coverage-commit.md`.
