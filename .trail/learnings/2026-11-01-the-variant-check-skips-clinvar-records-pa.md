---
title: "The variant check skips ClinVar records that have no alternate allele, because they used to wipe out their whole gene"
author: Pär
date: 2026-11-01
branch: fix/batch-9-clinvar
cost: { minutes: 13, tokens: 105710 }
id: l-43dd56d5
scope: repo
status: active
src: agent
touches:
  - src/robin/analysis/lightweight_gene_analysis.py
  - tests/test_clinvar_significance.py
  - reviews/plans/2026-11-01-batch-9-clinvar.md
  - CHANGELOG.md
  - reviews/2026-10-16-code-review.md
  - reviews/2026-10-16-code-review.pdf
---

## Decided
- The variant check skips a pathogenic ClinVar record whose ALT is not plain bases. ClinVar writes 13 records with ALT "." (HGVS `g.N=`): the record is the reference base itself, so no read can support or refute it. The results count them as `clinvar_records_not_checkable`, and the log gets one warning.
- Skip, not "needs review": as needs review, every PTEN and BRCA2 result with nothing found would read "insufficient evidence" for good, for a reason unrelated to the sample.

## Watch out
- An exception anywhere inside a region still drops that gene from the results with only a log line (`analyze_bed_regions` catches per region). Before this fix a lone ALT "." record did that to PTEN (AML panel) and BRCA2 (both panels) on every run.
- When measuring ClinVar input, tally every shape of each field the code reads (ALT by alphabet, not just "has a comma").

In the code: `_is_plain_allele` and the skip in `_parse_clinvar_vcf` (`src/robin/analysis/lightweight_gene_analysis.py`); tests `test_a_record_without_an_alternate_allele_does_not_drop_its_region` and `test_a_lone_record_without_an_alternate_allele_does_not_drop_its_region` in `tests/test_clinvar_significance.py`.
