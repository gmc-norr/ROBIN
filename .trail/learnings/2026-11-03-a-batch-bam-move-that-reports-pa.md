---
title: A batch BAM move that reports an error is judged by where the BAM is
author: Pär
date: 2026-11-03
branch: fix/batch-10e-recorded-items
id: l-4cafa66e
scope: repo
status: active
src: agent
touches:
  - src/robin/analysis/target_analysis.py
---

## Decided
- When moving a batch BAM into the sample folder raises an error, the step is done if the BAM is in the sample folder or target.bam lists it as merged. That is the same question a finish asks about a BAM gone from its claim folder; one helper answers both. Over NFS a rename can take effect and still report an error (Linux rename(2), BUGS). Taking every error as "not moved" kept a claim folder with incomplete.txt whose reads finalize had merged, so every later finalize refused the SNP calls, and nothing could clear it (the install list was already gone).
- The "merged" half is needed, not only "in the sample folder": finalize merges and deletes batch BAMs under its own lock, not the accumulation lock, so it can do so between the move and its late error.

## Tried and ruled out
- Asking only whether the BAM left the claim folder: a BAM gone from both places is lost and must still be reported.
- Retrying the move: the state on disk after the error already answers it.

## Watch out
- An error from a file operation is not proof that nothing happened. Where an error becomes a lasting record (a failed file, a kept folder, a ledger line), check what is on disk first. The other install steps raise and leave the list, and the next finish does them again by what is on disk, so only the batch BAM step needed this.

In the code:
- `src/robin/analysis/target_analysis.py`: `_batch_bam_in_place`, `_install_batch_bam`.
- `tests/test_coverage_commit.py`: C15, `test_a_batch_bam_move_that_reports_an_error_is_judged_by_where_the_bam_is` (took_effect, merged_meanwhile, did_not). It makes the real rename, then raises ENOENT; no NFS fault is made.
- Plan: `reviews/plans/2026-11-03-batch-10e-round5.md`.
