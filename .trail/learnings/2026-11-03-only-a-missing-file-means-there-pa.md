---
title: Only a missing file means there is no BAM ledger; a lookup error is not that
author: Pär
date: 2026-11-03
branch: fix/batch-10e-recorded-items
cost: { minutes: 12, tokens: 298198 }
id: l-066ecfe6
scope: repo
status: active
src: agent
touches:
  - src/robin/file_ledger.py
---

## Decided
- In the BAM ledger, only `FileNotFoundError` means "there is no ledger". `FileLedger.load()` opens the file: missing gives an empty state, and any other `OSError` raises `LedgerError`. `clear_bam_failures` relies on that, so a ledger that can't be looked up keeps the install list for the next finish.

## Watch out
- `os.path.exists` returns False for any error while looking the file up (an I/O error, a folder ROBIN may not search), not only for a missing file. Don't use it where "missing" and "can't tell" lead to different actions.
- A start whose ledger can't be looked up now stops, as one that can't be read does. It used to start as if there were no ledger and take every BAM in again.

In the code: `src/robin/file_ledger.py` (`FileLedger.load`, `record_bam_failure`, `clear_bam_failures`). Tests: `tests/test_coverage_commit.py` C19, `tests/test_bam_ledger.py` L8, L9. Plan: `reviews/plans/2026-11-03-batch-10e-round7.md`.
