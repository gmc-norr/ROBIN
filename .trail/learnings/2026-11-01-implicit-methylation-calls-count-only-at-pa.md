---
title: "Implicit methylation calls count only at CpGs read on the read, and records without probabilities are dropped"
author: Pär
date: 2026-11-01
branch: fix/batch-5-methylation
cost: { minutes: 42, tokens: 253178 }
id: l-3bb31140
scope: repo
status: active
src: agent
touches:
  - tests/test_matkit_reads.py
---

## Decided
- An implicit MM group (`C+m.` or `C+m`) adds an unmethylated call only at Cs the read reads as CpG (C then G on SEQ; on a reverse read, the G of a "CG"). Adding every C made a site of every C on every read (hundreds of MB per BAM) and gave non-CpG panel probes a 0% row.
- A record with MM but no ML is skipped, as modkit drops it. pysam reports its calls with probability -1, which matkit's arithmetic counted as canonical.
- Not every panel start is a CpG's C, correcting the earlier methylation learning: the 3,091 1 bp "ch" rows are non-CpG Cs, and some 1 bp "cg" rows sit one base off a CpG in hg38 (2 of 4 checked; 12 others sit on one). The merge still treats every start as a C; CpG calls give "ch" rows no data.

## Tried and ruled out
- Passing the reference to `run_matkit` so implicit calls follow the reference's CpGs, as modkit does: exact, but bed_conversion gets the reference only when one is configured, and implicit-mode BAMs are rare (dorado writes `?`). Measured cost of reading CpG from the read instead: at a read with an error at a CpG's other base, modkit counts the skipped C and ROBIN doesn't (same as `?` mode data).
- Trusting a mutation helper's "caught": pytest's exit 5 (no tests ran) and 4 (bad usage) are not failures. Count exit 1 only, and add a no-op control mutant that must survive.

## Watch out
- modkit `--cpg` decides CpG from the reference, not the read: in `.` mode a read reading CA at a reference CpG still counts there.
- `drop_duplicates()` in `_read_cpg_panel` is redundant (`_panel_starts` dedupes and a semi-join can't double rows); no test can fail if it goes.

In the code: `src/robin/analysis/utilities/matkit.py` (`_with_implicit_canonical_calls`, the ML check in `process_bam_counts_improved`); tests `tests/test_matkit_reads.py` (no modkit needed), `tests/test_meth_vs_modkit.py`. Measurements: plan `reviews/plans/2026-11-01-batch-5-methylation.md`, addendum.
