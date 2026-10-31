---
title: "Methylation counts are defined as modkit's, and the parquet stores counts per CpG"
author: Pär
date: 2026-10-31
branch: fix/batch-5-methylation
cost: { minutes: 52, tokens: 382953 }
id: l-9ec27fea
scope: repo
status: active
src: agent
touches:
  - src/robin/analysis/utilities/merge_bedmethyl.py
---

## Decided
- ROBIN's methylation numbers mean what `modkit pileup --cpg --combine-mods --combine-strands --filter-threshold 0.73` gives for the same reads (agreed with the user before any code). A test runs an installed modkit against ROBIN on synthetic BAMs; it skips where modkit is missing.
- The cumulative parquet stores counts per CpG (minus rows moved onto the C); percent is always recomputed from them. Every write re-sums, so an old-format file migrates through the same step.
- Only primary alignments count, as modkit does by default; that also settles the low-QS read whose supplementary came first.
- Every panel interval's start is a CpG's C, 1 bp rows included.

## Tried and ruled out
- Expecting `--combine-mods` to add 5mC and 5hmC probabilities: modkit takes the most likely class, so 0.45 + 0.45 is a failed call.
- Comparing modkit pileup's Nfail column: 0.6.4 prints 0 for calls that `modkit extract calls` marks failed. Compare n_mod and n_canonical.
- Keeping rows by interval overlap plus an inner join: one row per overlapping panel interval.

## Watch out
- modkit 0.6 needs `--modified-bases` (pass `C` with `--combine-mods`).
- MM counts bases on the read as sequenced: a reverse read's C positions are counted on the reverse complement.
- `merge_bedmethyl` imports `sturgeon` at the top; without Sturgeon installed, tests need a stand-in module.
- Change `_PANEL_CACHE_VERSION` whenever the panel's meaning changes, or old caches in sample folders hide the change.

In the code: `src/robin/analysis/utilities/matkit.py` (`_with_implicit_canonical_calls`, `run_matkit`, both writers), `src/robin/analysis/temp_utilities.py` (`load_cpg_panel`, `_cpg_counts`, `_cpg_table`, `merge_modkit_files`), tests `tests/test_meth_vs_modkit.py`, `tests/meth_fixtures.py`, `tests/test_meth_consumers.py`. Plan with the locked rule: `reviews/plans/2026-11-01-batch-5-methylation.md`.
