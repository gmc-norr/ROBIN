---
title: "Each BAM is added once only if the batcher, the runner and Ray's own retries are all covered; failures show at every start"
author: Pär
date: 2026-11-02
branch: fix/batch-10d-file-ledger
cost: { minutes: 1, tokens: 8444 }
id: l-2e772cc2
scope: repo
status: active
src: agent
touches:
  - src/robin/file_ledger.py
  - tests/test_bam_ledger.py
  - tests/test_master_csv_once.py
---

## Decided
- A batch's job gets its own copy of its first BAM's job details: its own metadata, with that BAM's results and errors shared. A BAM's steps share one context, so writing the batch list (`_batched_job`) into it let a later batch change the list of a batch still waiting to run. One BAM was then added twice and another left out.
- The six per-BAM Ray tasks (preprocessing, bed_conversion, mgmt, cnv, target, fusion) have `max_retries=0`. Ray re-runs a task whose worker died (default 3), and a step that had added a BAM's data added it again. A dead worker's step now fails, and after a restart shows as cut short. The per-sample steps keep Ray's default.
- `master.csv` holds only `last_bam_added`; the list of counted BAMs is `<sample>/master_csv_bams.txt`, one key per line. Write order: `master.csv` first, then the list. If the last key is missing from the list, the next update appends it.
- Failed and cut-short steps are reported at every start, and never run again (the user's choice). They come from the ledger's end records (`outcome`, `error`, `file_failed`). A fail-only BAM's expected failure is not a failure.

## Tried and ruled out
- Adding a newline before appending after a torn line: the torn fragment stays, as a bad line in the middle. The append cuts the torn tail off instead (Codex's point).
- Raising Python's csv field limit: other readers of `master.csv` still break. Excel cuts cells at 32,767 characters.
- A warning instead of the automatic resume (the user asked): it was kept. Without it, a restart means re-analysing the whole run, or numbers counted twice.

## Watch out
- pandas reads a hex key like `000…01` as a number: read `last_bam_added` as text (`dtype=str`).
- In the non-Ray runner, sibling jobs share a context and nothing removes `file_failures` from it. Count only the files a job adds: take a count before the handler.
- To check "each BAM is added once", list every layer that can run a step: the batcher, the runner and Ray itself. Size a growing record at a real run's scale.

In the code: `src/robin/file_ledger.py` (`append_whole_lines`, `_earlier_failure`, `record_step_end`, `file_failure_count`); `_convert_batched_jobs` in `workflow_ray.py` and `workflow_simple.py`; `_NO_RETRY` in `workflow_ray.py`; `MasterCSVManager.update_master_csv`, `_read_counted_keys`, `_append_counted_key`. Tests: the "Round 2" section of `tests/test_bam_ledger.py`, and the end of `tests/test_master_csv_once.py`.
