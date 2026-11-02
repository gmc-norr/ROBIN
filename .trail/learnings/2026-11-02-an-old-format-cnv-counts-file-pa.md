---
title: "An old-format CNV counts file that can't be moved fails the run, never starts from none"
author: Pär
date: 2026-11-02
branch: fix/batch-10c-cnv-counts
cost: { minutes: 1, tokens: 1504 }
id: l-8d851cd5
scope: repo
status: active
src: agent
touches:
  - src/robin/analysis/cnv_analysis.py
  - tests/test_cnv_counts_shared.py
---

## Decided
- A sample's CNV counts have two ways in, and both fail the run on any error: the per-sample counts file, and the one-time move out of the old multi-sample `update_cnv_dict.pkl` in the sample's folder. A move that can't read the old file, or can't write the per-sample one, raises; the run fails before counting any BAM, the old file is kept, and the next run tries the move again.
- Going on after a failed move is not harmless: with no per-sample file the run starts from none, saves its own BAM, and from then on the per-sample file exists, so the old file is never read again and its counts are lost for good.

## Watch out
- The original code's move caught every exception with a warning. When a rule like "a failure says failed" is applied to a loader, check every `except` on every path into the data, not only the lines being changed: batch 10c's first round gated the per-sample read and missed this one.
- `atomic_open` adds a failure point (the rename) to any write it wraps, so a quiet `except` around a write that used to succeed can start firing.

In the code: `src/robin/analysis/cnv_analysis.py` (`_load_copy_numbers`). Tests: `tests/test_cnv_counts_shared.py` (`test_a_failed_move_out_of_the_legacy_file_fails_the_run`, `test_a_legacy_file_that_cannot_be_read_fails_the_run_and_is_kept`, and two controls).
