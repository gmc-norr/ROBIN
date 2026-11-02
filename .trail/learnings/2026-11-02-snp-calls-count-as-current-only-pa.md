---
title: "SNP calls count as current only if every read reached target.bam, including reads a worker is still adding"
author: Pär
date: 2026-11-02
branch: fix/batch-10b-finalize
cost: { minutes: 33, tokens: 466911 }
id: l-ce525a11
scope: repo
status: active
src: agent
touches:
  - reviews/plans/2026-11-02-batch-10b-finalize-and-snp.md
  - src/robin/analysis/target_analysis.py
  - src/robin/gui_launcher.py
  - tests/test_snp_every_read.py
  - tests/test_snp_reads_current.py
  - CHANGELOG.md
  - reviews/2026-10-16-code-review.md
  - reviews/2026-10-16-code-review.pdf
---

## Decided
- A claim folder (`_staging/_batch_<id>/`) counts as reads not in target.bam. Finalize waits for its worker (up to 300 s), and the SNP job fails while one is left. A worker holds a per-claim lock (`_locks/claim_batch_<id>.lock`) from claim to end, so a waiter can tell a live claim (lock held) from one left by a crash or a failure (lock free, folder still there).
- When a batch BAM can't be written, or a staged file can't be read, the claim folder is kept with an `incomplete.txt`, not deleted. Their coverage may already be counted, so they must never go back to the top of `_staging/`, or they would be counted twice.
- The SNP job goes ahead only on finalize results `success` (with no `batch_bam_error`), `no_files` or `no_pending_files`, and with no `claims_left`. It is an allow-list, so a new failure status fails safe.
- The SNP record (`clair3/target_bam_source.json`) is target.bam's batch list plus its size and mtime_ns. It replaced round 1's batch-names-only `target_bam_batches.txt` (never released). Any rewrite changes size or time.
- A Clair3 region counts only if it exits 0 and writes both `snv.vcf.gz` and `indel.vcf.gz` with an mtime at or after the run's start (less 2 s).

## Tried and ruled out
- Deleting old Clair3 outputs before a run: the container runs as root (no `user=`), so its region folders may not be ours to delete. Hence the mtime check.
- One shared reader-writer lock over all accumulations: a pending exclusive flock doesn't stop new shared holders, so finalize could starve during a busy live run. Per-claim locks wait only on the claims that exist when finalize looks.

## Watch out
- Not measured: whether ClairS-TO always writes `indel.vcf.gz` for a region with no candidates. If it doesn't, split runs now fail with a message naming the region.
- A kept claim folder fails that sample's SNP runs until someone deals with it. There is no automatic recovery yet.

In the code: `src/robin/analysis/target_analysis.py`: `_claims_left`, `_claim_lock_path`, `_keep_incomplete_claim`, `_merge_missed_reads`, `_target_bam_identity`, `accumulate_staged_files`, `finalize_accumulation_for_sample`, `run_snp_analysis` (its region loop). Tests: `tests/test_snp_every_read.py`.
