---
title: "The lightweight variant check's pileup reference is now the most common base, not always T"
author: Pär
date: 2026-10-31
branch: fix/batch-3-quick-wins
cost: { minutes: 1, tokens: 1901 }
id: l-d84f5206
scope: repo
status: active
src: agent
touches:
  - "src/robin/analysis/lightweight_gene_analysis.py#_perform_targeted_pileup"
  - tests/test_lga_variant_support.py
  - "tests/lga_fixtures.py#write_bam"
---

## Decided
- The batch 2 learning on the lightweight variant check says the pileup's "reference" field is always "T". That was true then and is no longer true. The fix: max() now compares the counts, not the (base, count) tuples, so the reference is the most common base at that position.
- That learning stays as it was. This note replaces its one "always T" bullet. Learnings are a record of what a session knew, so they are never edited; a new one says what changed.

## Tried and ruled out
- Editing the old bullet in place. The unattended batch 3 run did this; it was undone at merge for the rule above.

## Watch out
- "Most common base" is still a guess, not the genome reference. A site where most reads carry the alternate allele gets that allele as its "reference". The ClinVar REF is the real reference; use it when you need one.
- Before the fix, insertion sites read "reference_mismatch" instead of "insertion_requires_review", and the real reference base was listed as a variant. Results files written before batch 3 may carry both.

In the code: _perform_targeted_pileup in src/robin/analysis/lightweight_gene_analysis.py (ref_base = max(..., key=lambda kv: kv[1])); tests in tests/test_lga_variant_support.py; write_bam in tests/lga_fixtures.py takes an optional CIGAR as a third tuple element.
