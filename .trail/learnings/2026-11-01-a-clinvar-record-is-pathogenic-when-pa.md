---
title: "A ClinVar record is pathogenic when its germline, oncogenicity or somatic label says so, never on conflicting alone"
author: Pär
date: 2026-11-01
branch: fix/batch-9-clinvar
cost: { minutes: 35, tokens: 327762 }
id: l-09371d30
scope: repo
status: active
src: agent
touches:
  - reviews/plans/2026-11-01-batch-9-clinvar.md
  - src/robin/analysis/clinvar_significance.py
  - src/robin/analysis/lightweight_gene_analysis.py
  - src/robin/analysis/snp_processing.py
  - src/robin/gui/components/snp.py
  - src/robin/reporting/sections/variants.py
  - tests/test_clinvar_significance.py
  - CHANGELOG.md
  - reviews/2026-10-16-code-review.md
  - reviews/2026-10-16-code-review.pdf
---

## Decided
- A ClinVar record counts as pathogenic when ClinVar's own aggregate label says so in any of its three classifications. CLNSIG (germline): Pathogenic or Likely_pathogenic. ONC (oncogenicity): Oncogenic or Likely_oncogenic. SCI (somatic clinical impact): Tier_I_-_Strong or Tier_II_-_Potential. Labels are matched whole, after splitting the value on `|`, `/` and `,`. The user chose this on 2026-11-01, over a vote on CLNSIGCONF.
- One helper decides it for the variant check, the GUI's SNP table and the report. They used to disagree, and all three matched the word "pathogenic" anywhere.
- The variant check keys ClinVar records by (chrom, pos, ref, alt). In the 2026-10-25 ClinVar file no key repeats, while 53,674 positions hold more than one record the old rule called pathogenic.

## Tried and ruled out
- Matching exact CLNSIG labels only (the code review's suggested fix for VAR-2). BRAF V600E (7:140753336 A>T) and TERT C228T (5:1295113 G>A) have the germline label "Conflicting_classifications_of_pathogenicity". The old substring match counted them by accident; only their ONC and SCI labels call them. The strict fix would have dropped both.
- A vote over CLNSIGCONF (what the report's unreached code meant to do). That is ROBIN's verdict, not ClinVar's, and it keeps about 30,000 germline-disputed records.

## Watch out
- For a brain-tumour panel the somatic labels are the ones that carry the hotspots. Some sites have only SCI (TERT 5:1295113 G>T), or CLNSIG "other" with ONC (an H3-3A record at 1:226064454).
- CLNSIG can hold `Pathogenic,_low_penetrance`, so split on commas too.
- Saved SNP tables (`clair3/snpsift_output_display.json`) and saved variant-check results keep the flags they were made with until that step runs again. The report reads the VCF each time.
- Measure a rule change on the real file first: ClinVar's GRCh38 VCF (about 200 MB, the URL in `utils/clinvar_manager.py`). Read its labels and run all three checks over every record.

In the code: `src/robin/analysis/clinvar_significance.py` (`is_clinvar_pathogenic`, `clinvar_label_text`); used by `lightweight_gene_analysis.py`, `snp_processing.py`, `reporting/sections/variants.py`; tests in `tests/test_clinvar_significance.py`; plan `reviews/plans/2026-11-01-batch-9-clinvar.md`.
