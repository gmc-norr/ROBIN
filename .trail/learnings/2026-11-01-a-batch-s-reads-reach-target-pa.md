---
title: "A batch's reads reach target.bam only after its coverage is committed, finalize merges under its own lock, and the IGV copy follows target.bam"
author: Pär
date: 2026-11-01
branch: fix/batch-6-target
cost: { minutes: 28, tokens: 134466 }
id: l-cd40c8f5
scope: repo
status: active
src: agent
touches:
  - src/robin/analysis/target_analysis.py
  - src/robin/analysis/target_coverage.py
---

## Decided
- A batch's reads become visible to finalize only after its coverage is committed: the batch BAM is built outside any lock in the batch folder, and renamed into the sample folder inside the accumulation lock, after the commit. A commit that fails (lock timeout, unreadable totals) then leaves the coverage and target.bam agreeing: both miss the batch.
- Finalize merges under its own "finalize" lock, not the accumulation lock. Batch BAMs appear only by rename after their commit, so its glob is a safe snapshot. Holding the accumulation lock made a long merge (about 10 s per 24k records) time out workers' 120 s commits and drop their coverage.
- igv/igv_ready.bam (read by IGV and the per-gene variant check) must follow target.bam: it counts as current only when its header's merged-batch @CO list equals target.bam's and it isn't older. Finalize refreshes an existing one; it is written by rename, one builder at a time.
- The merge sets the @CO list with samtools reheader (copies the compressed records), not a Python pass over every record.
- Contig bases come from the same aligned-block pass as the covered stretches. Equal to bedcov -j: bedcov doesn't cap depth (equal at 12,000x).

## Tried and ruled out
- Freshness of the IGV BAM by file time alone: a sort that read the older target.bam can finish after the newer one is written, so its file is newer but its reads old. The merged-batch list catches that; time is kept for folders from before the list.

## Watch out
- The earlier target-coverage learning says an old folder's covered bases "are shown as unknown". Nothing shows covered bases: in coverage_main.csv they are blank cells.
- In an old folder the coverage-over-time chart changes meaning at the upgrade: covered bases over genome length before, mean depth after (fixture: 0.4364, then 1.1414).
- The AML panel's "270 overlapping rows" counts rows that overlap an earlier row in file order; 473 rows overlap some row, and 3.59% of the summed length is counted twice in the single target-depth number.

In the code: `src/robin/analysis/target_analysis.py` (`accumulate_staged_files`, `_write_batch_bam`, `_merge_into_target_bam`, `_reindex_if_stale`, `_follows_target_bam`, `finalize_accumulation_for_sample`, `ensure_sorted_igv_bam`), `target_coverage.py` (`aligned_blocks`). Tests: `tests/test_target_finalize.py`. Plan addendum: `reviews/plans/2026-11-01-batch-6-target-coverage.md`.
