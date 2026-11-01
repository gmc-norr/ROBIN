# Batch 8: the panel's summary depth counts each base once

**Locked before code, 2026-11-01.** Base `c81cc37`. Branch `fix/batch-8-panel-depth`. Done by the controller (no agents, to save the user's tokens). An external review (Sol 6.1 / Astra) comes before the merge.

## Goal and rule

"Targets Estimated Coverage" is samtools `bedcov -j` over the panel with its overlapping rows merged, divided by the merged length, on all of a sample's reads, so each base counts once. The enrichment ratio and the GUI's quality label are made from it.

Each row's own depth, and "% of targets above 30×/20×/10×", don't change.

**Spec:** the user's "please fix it" (2026-11-01), and batch 6's known limit (TGT-3 Status in `reviews/2026-10-16-code-review.md`).

## Measured before this lock

- **Shipped panels:**
  - AML: 1,180 rows → 910 merged regions. 3,173,388 of 88,388,915 bases (3.59%) are counted more than once. One interval sits under two names (AK6/TAF9).
  - rCNS2: 243 rows, no overlaps, so nothing changes for it.
- **Error formula:** row-sum depth = merged depth + (D/L)·(overlap depth − merged depth), with D the length counted more than once and L the summed length.
- **Batch-6 fixture (GENE_A/GENE_B overlap at chr1:300-400):** the row sum gives 2.808333; samtools over the hand-merged panel gives 2.566364.
- **Where the number is made:** `target_depth(bed_df)` in `target_coverage.py`, the only place. It is called by the GUI coverage panel (`coverage.py` ~6663), the GUI summary card (`summary.py` ~1176) and the report (`reporting/sections/coverage.py` ~78).
- **Where it is written:** `process_file_with_staging` stages each BAM's totals, and `_commit_batch` adds them under the accumulation lock (TGT-4).

## Design

1. **Staging:** each BAM's staging step also writes `union_{n}.parquet`: aligned bases (`_bedcov_j`) over the panel's merged regions (`_merge_intervals`). It is written before the coverage file, which is still written last by rename.
2. **Claiming:** the union part is optional. A staged BAM from before this change has none and is still claimed.
3. **Commit:** `_commit_batch` adds the batch's union parts to `target_union_coverage.csv` (`chrom,start,end,bases`, written atomically) under the lock. The file is written only while every batch of the sample has carried its union part.
   - A batch without one removes the file.
   - A sample whose totals exist without the file stays without it.
   - With no file, the number is unknown, and the readers fall back to the row sum, the old number.
4. **Reading:** one helper, `panel_depth(sample_dir, bed_df)`, gives the merged-panel depth from the file, or the row sum without it (logging a warning if the file exists but can't be read). The three readers use it.

## Tests (measured on `c81cc37` before this lock: 19 failed, 4 passed as planned)

| # | Test | On `c81cc37` |
|---|---|---|
| U1 | `test_panel_depth_matches_samtools_over_the_merged_panel`: the batch-6 runs (1, 2 and 4 BAMs × in order, reversed, all at once, two processes, finalized twice). The helper's depth equals samtools over `MERGED_PANEL` (rel 1e-9). | 15 failed (no helper) |
| U1b | `test_every_overlap_shape_counts_each_base_once`: a panel adding a row inside another, the same interval under two names, and a bridging row, on 2 BAMs. The report, the GUI summary and the helper equal samtools over `MERGED_PANEL_OVERLAPS`. The fixture tells this apart from the row sum (> 0.01). U2: each row's bases equal `bedcov -j` for that row. | failed (no helper) |
| U1c | `test_the_gui_coverage_panel_uses_the_shared_helper` (source check; the panel's code sits in nested closures) | failed |
| C | `test_coverage_consumers.py`: the report and the GUI summary show samtools' target depth. This **existing test's reference changes** from the row sum to the merged panel (stricter); stated here before code. | 2 failed (2.808 vs 2.566) |
| U3 | guard: a BAM staged before this change (no union part), then a new one. No union file; report and summary equal the row sum. | passed |
| U4 | guard: a folder written before this change (no union file), then a new BAM. No union file; the report equals the row sum. | passed |
| K1 | control: a panel without overlaps. Samtools over the merged panel equals the row sum, and the report shows it. | passed |
| K2 | control: an old folder with no new BAMs. Report and summary equal the row sum. | passed |

The batch-6 tests (totals, target.bam, read counts, every run the same) pass on `c81cc37` and must still pass.

## Pass, all of

- U1, U1b, U1c and C pass after; U3, U4, K1 and K2 pass on both.
- One mutant per design part (1, 3 and 4, plus the fallback rule in 3) is caught, and the no-op control survives.
- The full suite passes.
- `robin.yml` and `pyproject.toml` are unchanged.

**Fail:** any one; then fix and re-run, or report it.

## Out of scope

- Per-base panel depth (D4, the later feature).
- A panel changed mid-run: merged regions are keyed by position, as rows are today.
- "% of targets above".

## Result (2026-11-01, on `c584ad3`): supported

- **On `c81cc37`:** the 19 fix tests failed and the 4 guards and controls passed (the table above). After the fix, all pass, along with batch 6's samtools tests.
- **One test added after the lock (guard):** `test_a_bam_staged_before_this_change_after_a_new_one_removes_the_union_file`. The rule "a batch without the parts removes the file" had no test, so its mutant would have survived. It passes on `c81cc37` and after.
- **Mutation check:** 5 of 5 caught (staging, commit sum, old-folder rule, removal rule, readers); the no-op control survived (76 passed).
- **Full suite:** 656 passed. `robin.yml` and `pyproject.toml` are unchanged.
- **Shipped panels, through the code's own merging:** AML 1,180 rows → 910 regions (3.59% counted more than once), as measured before the lock; rCNS2 243 → 243.
- **Not yet done:** the external review (Sol 6.1 / Astra), before the merge.
