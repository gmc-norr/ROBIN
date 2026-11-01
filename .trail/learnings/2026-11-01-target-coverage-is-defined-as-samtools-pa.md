---
title: "Target coverage is defined as samtools on all reads, and covered bases come from a union of stretches"
author: Pär
date: 2026-11-01
branch: fix/batch-6-target
cost: { minutes: 66, tokens: 483095 }
id: l-21545c04
scope: repo
status: active
src: agent
touches:
  - src/robin/analysis/bam_preprocessor.py
---

## Decided
- Target coverage, read counts and target.bam mean what samtools gives when run once on all of a sample's reads (agreed with the user). Depth adds up across BAMs; covered bases don't, so each BAM stages its covered stretches and the sample keeps their union.
- Aligned bases follow `samtools coverage`: `bedcov -j` (no deletions), supplementary records count, secondary don't. Read counts and yield count primary records only.
- target.bam is `samtools view -L` over the panel plus the master BED's regions, for every BAM as it arrives; the call threshold only gates SNP calling. The final merge includes the existing target.bam and records merged batch names in its header (@CO), so a crash can't merge one twice.
- One locked step reads, adds and writes the totals; a staged BAM is five files and its coverage file is written last, by rename.
- Per-base panel depth (D4 via mosdepth) was offered and set aside as a separate feature; the user chose the stretch list for covered bases.

## Tried and ruled out
- Summing per-BAM `covbases`: breadth doesn't add up (4,048 vs 1,669 on the fixture).
- Keeping "all records of reads whose primary overlaps" for target.bam: not what samtools gives, and it carries off-panel pieces nobody reads.
- A two-process test that slows only work outside the lock: an unlocked commit survived it. Slow a step inside the commit too.

## Watch out
- A folder written before this change has summed covered bases; they are shown as unknown, not migrated.
- `_get_master_bed_path` imports fusion_work and swallows errors: without networkx or scikit-learn installed it silently returns None, and the master BED is ignored.
- `samtools coverage`'s startpos/endpos are 1-based (length = endpos − startpos + 1); a BED row's length is end − start.

In the code: `src/robin/analysis/target_coverage.py`, `target_analysis.py` (`process_file_with_staging`, `accumulate_staged_files`, `_commit_batch`, `finalize_accumulation_for_sample`), `bam_preprocessor.py`, `cnv_analysis.py` (`cnv_handler`, `save_cnv_files`). Tests: `tests/test_target_vs_samtools.py`, `tests/test_target_coverage.py`, `tests/test_master_bed_panel.py`, `tests/test_coverage_consumers.py`. Plan: `reviews/plans/2026-11-01-batch-6-target-coverage.md`.
