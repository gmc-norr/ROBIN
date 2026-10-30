---
title: "Batch 2: lightweight variant check (LGA) sample identity, interval and evidence rules"
author: Pär
date: 2026-10-30
branch: fix/batch-2-sample-identity
cost: { minutes: 35, tokens: 304131 }
id: l-baecec8a
scope: repo
status: active
src: agent
touches:
  - CHANGELOG.md
  - reviews/2026-10-16-code-review.md
  - reviews/2026-10-16-code-review.pdf
  - reviews/2026-10-16-repro/independent/README.md
  - reviews/2026-10-16-repro/independent/reproduce.py
  - src/robin/analysis/lightweight_gene_analysis.py
  - src/robin/gui/components/coverage.py
  - tests/test_lga_variant_support.py
  - tests/test_lga_intervals.py
  - tests/test_lga_sample_inputs.py
---

## Decided
- LGA (src/robin/analysis/lightweight_gene_analysis.py) reads only the sample's own folder: resolve_sample_folder(base_dir, sample_id), then <folder>/igv/igv_ready.bam and the given or <folder>/targets_exceeding_threshold.bed. Missing → FileNotFoundError naming the path. Never os.walk the monitored root.
- Results JSON records metadata.input_bam/input_bed/sample_id. The GUI shows a results file only if lga_results_belong_to_sample() says its input_bam is inside the sample's folder (old files are hidden).
- A ClinVar site counts as found only with 3+ alternate reads (is_supported: low/medium/high_frequency). needs_review: 1-2 reads, indels, same-length multi-base changes (mnv_requires_review), reads showing neither allele (reference_mismatch). Status "insufficient_evidence" when nothing is supported but something needs review.
- Same-named BED intervals: merge overlaps (_merge_intervals), pile up each, match each interval's sites to its own pileup.

## Tried and ruled out
- Concatenating per-interval pileups without merging: overlapping intervals double-count coverage.

## Watch out
- len(ref) != len(alt) is not the only uncallable case: CC>TT is same-length and a one-base pileup can't call it.
- The pileup's "reference" field is always "T" (max over (base, count) tuples compares letters). Don't trust it.
- reviews/2026-10-16-repro/independent/reproduce.py copies in only the functions it names; when LGA gains a helper, add it to the names lists or the tests NameError instead of showing the fix.
- tests/lga_fixtures.py writes tiny indexed BAMs, ClinVar VCFs and BEDs; keep ClinVar POS off the BED start.
