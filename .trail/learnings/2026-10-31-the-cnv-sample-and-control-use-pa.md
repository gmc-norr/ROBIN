---
title: "The CNV sample and control use different median rules, up to 0.035 copies apart"
author: Pär
date: 2026-10-31
branch: fix/batch-4-cnv-reference
cost: { minutes: 5, tokens: 47184 }
id: l-85c8f78a
scope: repo
status: active
src: agent
touches:
  - tests/test_cnv_control_reference.py
---

## Decided
- The CNV control's median and the sample's median use different chunking rules: the control's chunks stay within each contig (our numpy profile), the sample's are the library's, laid end to end across contigs. Measured on the public control where nothing wraps, the two profiles differ by up to 0.035 copies (at 25,000 sample reads); 0.015 at 50,000, 0.0005 at 100,000, 0.024 at 200,000. It does not shrink steadily with depth. The earlier learning's "at most 0.024" came from five sampled read counts only.
- That gap is a baseline offset in the normalized profile, about a tenth of the ±0.4 call threshold, so it was left as is.

## Watch out
- To remove the offset, profile the sample with the same numpy function from its accumulated counts instead of the library's result, so both sides share one median rule and 64-bit sums.
- A bound measured at a few sample sizes is not a bound: scan the whole range before writing "at most".

In the code: `src/robin/analysis/cnv_reference.py` (`reference_cnv_profile`), `src/robin/analysis/cnv_analysis.py` (`run_cnv_analysis_direct`: `r_cnv` from the library, `r2_cnv` from `reference_cnv_profile`), `src/robin/analysis/cnv_subprocess.py`.
