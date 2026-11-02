---
title: "Each Clair3 run writes into its own new folder, and a staged set missing a part fails the SNP run"
author: Pär
date: 2026-11-02
branch: fix/batch-10b-finalize
cost: { minutes: 14, tokens: 93009 }
id: l-06cbae83
scope: repo
status: active
src: agent
touches:
  - src/robin/analysis/target_analysis.py
  - tests/test_snp_every_read.py
  - reviews/plans/2026-11-02-batch-10b-finalize-and-snp.md
  - CHANGELOG.md
  - reviews/2026-10-16-code-review.md
  - reviews/2026-10-16-code-review.pdf
---

## Decided
- `run_snp_analysis` makes a new folder, `clair3/run_<8 hex>/`, for each calling run, and reads only VCFs there. A time window ("written since the run started, less 2 s") let a file from a run that ended just before pass; a new folder makes that impossible. Earlier `run_*` folders are deleted when possible. The container writes as root, so some may stay.
- Finalize lists staged sets whose coverage file exists but a part is missing (`incomplete_staged`), and the SNP job fails on them. Staging writes coverage last, by rename, so such a set is broken, not in progress. The list is read under the accumulation lock, because a claim moves a set's files one by one.

## Watch out
- This corrects the earlier learning on SNP calls being current: not "any rewrite" of target.bam changes its record. Every way ROBIN writes it does. A replacement from outside ROBIN that keeps the old size and time is not noticed.
