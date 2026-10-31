---
title: "The per-gene variant path fails because its pileup starts one base late, not because of chromosome names"
author: Pär
date: 2026-10-31
branch: fix/batch-3c-review-findings
cost: { minutes: 14, tokens: 96895 }
id: l-033533a0
scope: repo
status: active
src: agent
touches:
  - "src/robin/analysis/lightweight_gene_analysis.py#_get_gene_coordinates"
  - tests/test_lga_variant_support.py
---

## Decided
- Correction to the 2026-10-31 learning on 1-based ClinVar sites: it says analyze_genes / _analyze_single_gene "is broken on chr-prefixed BAMs (it looks for contig 1)". That was wrong. _perform_targeted_pileup already maps "1" to "chr1" and back. The real fault was _get_gene_coordinates returning the 1-based min(position) as the 0-based pileup start, so the first site was never piled up and a one-site gene piled up an empty range (11..11). It now returns min - 1.

## Tried and ruled out
- Reading "Chromosome ... not found" into an empty pileup without checking: probe the pileup range first. A report whose chromosome field says "1" is only the label from ClinVar.

## Watch out
- Every coordinate hand-off in this module mixes 1-based VCF positions with 0-based pysam/BED ranges. Pin each with a one-site test.
