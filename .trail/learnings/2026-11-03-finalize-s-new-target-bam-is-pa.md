---
title: "Finalize's new target.bam is synced by samtools before it is renamed into place"
author: Pär
date: 2026-11-03
branch: fix/batch-10e-recorded-items
cost: { minutes: 12, tokens: 298198 }
id: l-bca818aa
scope: repo
status: active
src: agent
touches:
  - src/robin/analysis/target_analysis.py
---

## Decided
- This corrects the "Watch out" of the round-6 learning on clearing a saved BAM failure (l-c65e7566). Finalize does sync target.bam's contents: samtools (pysam 0.24.1, samtools 1.24) `fdatasync`s the new target.bam and its index before `_merge_into_target_bam` renames them into place. An `fdatasync` error makes the merge raise before target.bam is replaced. Codex measured it with a syscall trace.

## Watch out
- What is left, not measured: finalize does not sync the sample folder between renaming target.bam into place and deleting the batch BAMs.
- Judge what reaches the disk by a syscall trace, not by the Python code: the native library below it does its own syncs.

In the code: `src/robin/analysis/target_analysis.py` (`_merge_into_target_bam`, and finalize's batch BAM deletion loop after it).
