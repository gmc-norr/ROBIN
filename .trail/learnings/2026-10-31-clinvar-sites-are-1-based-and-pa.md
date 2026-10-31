---
title: "ClinVar sites are 1-based and BED starts are 0-based, so a site at the BED start is outside the region"
author: Pär
date: 2026-10-31
branch: fix/batch-3c-review-findings
cost: { minutes: 12, tokens: 103397 }
id: l-945288e5
scope: repo
status: active
src: agent
touches:
  - "src/robin/analysis/lightweight_gene_analysis.py#_find_intersecting_variants"
  - tests/test_lga_variant_support.py
---

## Decided
- The lightweight variant check now picks a ClinVar site for a BED region when start < POS <= end. POS is 1-based (VCF); the region is a 0-based half-open BED interval, i.e. 1-based start+1..end, which is exactly what the pysam pileup (truncate=True) covers. The old start <= POS picked the base before the region.
- A picked site with no pileup column now means "no read covers it" and is kept as evidence "no_reads" (needs review), so "None Detected" can't stand on unread sites. That is only safe because of the boundary fix: otherwise the base before every region would read "no_reads" while reads cover it.

## Tried and ruled out
- Treating zero usable bases as "needs review" for indels too: a column where every read is deleted supports a ClinVar deletion there.

## Watch out
- tests/lga_fixtures.py keeps ClinVar POS off the BED start for this reason; test_a_site_just_before_the_bed_start_is_not_checked pins the boundary.
- analyze_genes / _analyze_single_gene has no callers and is broken on chr-prefixed BAMs (it looks for contig "1"); don't trust a test that goes through it to exercise anything.
