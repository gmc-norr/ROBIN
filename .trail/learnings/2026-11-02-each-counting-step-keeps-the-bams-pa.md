---
title: "master.csv counts each BAM once, known by its bytes; a ledger of step starts and ends does the rest of JOB-2"
author: Pär
date: 2026-11-02
branch: fix/batch-10d-file-ledger
cost: { minutes: 26, tokens: 248314 }
id: l-bfd9bb7a
scope: repo
status: active
src: agent
touches:
  - src/robin/bam_identity.py
  - src/robin/analysis/master_csv_manager.py
  - src/robin/analysis/bam_preprocessor.py
  - tests/test_master_csv_once.py
---

## Decided
- JOB-2 (a restart counts every BAM again): a ledger records when each of a BAM's steps starts and ends. After a restart, finished steps are skipped, unstarted ones run, and a step that was running at the stop is not run again but listed as possibly incomplete (it may have added its counts). The user's choice, 2026-11-02.
- A BAM is known by its bytes, `robin.bam_identity.bam_key`: blake2b-128 of a version tag, the size, the first 64 KiB and the last 64 KiB. A copy or a moved folder is the same BAM.
- `master.csv` counts each BAM once anyway (a `bams_added` JSON-list column, checked and written with the counts under `_locks/master_csv.lock`), because re-running a BAM's unstarted steps starts from its preprocessing. A failed update is a file failure (step "master_csv"), so the BAM's other steps still run.

## Tried and ruled out
- Making every counting step keep the BAMs it has added (exactly once even for a step cut short): mapped, it was six steps and seven batches, two large. Target coverage's commit writes about ten files one by one, and fusion appends part files with new ids per run, so neither can record "this BAM is in" in one write without a redesign.

## Watch out
- `MasterCSVManager._load_existing_data` used to return the defaults when `master.csv` couldn't be read, so every writer wrote over the sample's counts. It now raises.
- The non-Ray runner has no `file_failures` channel: per-file failures there are only logged.

In the code: `src/robin/bam_identity.py`, `src/robin/analysis/master_csv_manager.py` (`update_master_csv`, `_load_existing_data`), `src/robin/analysis/bam_preprocessor.py` (step 4). Tests: `tests/test_master_csv_once.py`.
