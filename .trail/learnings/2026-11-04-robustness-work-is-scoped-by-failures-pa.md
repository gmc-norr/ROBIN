---
title: "Robustness work is scoped by failures seen in real runs, not by every fault a review can imagine"
author: Pär
date: 2026-11-04
branch: fix/unreadable-bam
cost: { minutes: 6, tokens: 65997 }
id: l-ec2b9a44
scope: repo
status: active
src: agent
touches:
  - reviews/plans/2026-11-04-rescue-at-door-gate-a.md
  - tests/test_unreadable_bam.py
  - src/robin/analysis/bam_preprocessor.py
---

## Decided
- The owner (2026-11-04): robustness work on ROBIN is scoped by failures seen, or clearly likely, in real runs. A fault a review can imagine is not by itself a reason to build a mechanism. Name it as a known limit instead.
- Fix a measured bug the smallest way, through channels ROBIN already has. A BAM that can't be read to its end now fails preprocessing with an error naming it (the existing failed-job path), so nothing runs after it. It was being counted as a zero-read sample called unknown_sample_<file>.

## Tried and ruled out
- A "never stop a sample" redesign of the target step (snapshot trial phase 2). Three plan revisions grew it from about 300 to about 580 lines and over 60 fault cases. Every fault came from reviews, none from a real run. Each revision added saved state the next review picked apart, and some fixes were only needed because of other fixes. It is parked on branch exp/snapshot-trial, not refuted.
- Rescuing the readable reads of a damaged or cut-off BAM into a copy for all analyses. How often MinKNOW leaves such a file is unknown; ONT writes to a temporary file and then moves it. Most of the design served the unlikely cut-off case (paused writer, file completing later, moved folders). Parked until a real run shows a damaged BAM.

## Watch out
- When a plan grows across review rounds, stop and check its scope with the owner before revising it again.
- process_bam_reads still returns None for a BAM with no read group in its header. That old fallback (a sample named unknown_sample_<file>, with zero counts) is unchanged and is a separate question.

In the code: bam_preprocessor.py `UnreadableBam` and process_bam_reads; tests/test_unreadable_bam.py.
