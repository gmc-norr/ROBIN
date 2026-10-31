# Batch 5: methylation counts match modkit (step 5) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** every classifier gets, per CpG, methylated ÷ (methylated + unmethylated) calls over all reads seen so far, each read once, both strands together, unsure calls left out — what `modkit pileup --combine-strands --combine-mods` gives for the same reads.

**Architecture:** fix the per-read layer in `matkit` (fails out of coverage, implicit MM mode, QS decided per read), then the merge in `temp_utilities.merge_modkit_files` (0-based panel, semi-join, strands combined, counts summed and percent recomputed on every write, old files migrated). Consumers should then need no maths of their own; their inputs are tested at the classifier boundary.

**Tech Stack:** Python 3.12, polars, pyranges, pysam (synthetic BAMs with MM/ML tags), `modkit` 0.6.4 (bioconda, via pixi; tests that need it skip without it).

**Spec:** `reviews/2026-10-16-code-review.md` section 1 (METH-1..6, METH-13, METH-14), and the user's agreement on 2026-10-31 (Gate A): the principle above; modkit installed as the reference; METH-7 left for later; the case file built. Data path mapped on `2f29cf7` (Explore report, summarised in the Spec notes below).

**Spec notes (measured on `2f29cf7`):**
- Gate B (synthetic BAM, `?` mode, both strands): matkit and modkit agree exactly on n_mod and n_canonical at all 10 strand-sites. modkit 0.6.4's pileup prints Nfail = 0 even where `modkit extract calls` marks the same calls `fail=true`; both leave fails out of coverage, and modkit drops sites with no valid call. So Nfail is not compared.
- `parquet_filter.txt` is 0-based: cg13869341, cg02984631, cg14817997 read `CG` in hg38 at their listed starts (UCSC API). Its first row is `chr1 -1 0`.
- The live merge is `temp_utilities.merge_modkit_files`; the panel cache is `<sample>/parquet_filter.txt_1based.pgr_cache`, keyed by basename, so an old widened cache would survive a code fix.
- MGMT reads matkit's 18-column BED (`write_bedmethyl_improved`), stranded, not the parquet: per-read fixes reach it, merge fixes do not.
- Sturgeon, NanoDX's `NN_classifier`, and the RF/MGMT R scripts are not installed here (empty submodules); consumer tests capture the input at the call into them.

## Locked metric and decision rule (Gate C — committed before any code)

**Fixtures** (synthetic reference + BAMs, `tests/meth_fixtures.py`), covering every trap: calls on both strands; unsure calls (ML 128); a read whose 5mC and 5hmC each carry 0.45; MM in explicit `?` mode, implicit `.` mode, no-marker mode, and an implicit group with an empty list (`C+m.;`); a QS 5 read whose supplementary alignment comes before its primary, and one whose primary comes first; a panel with overlapping and exactly duplicated CpG rows and a negative-start row.

**Reference:** `modkit pileup <union BAM> out.bed --ref ref.fa --cpg --combine-mods --modified-bases C --filter-threshold 0.73 --combine-strands`, on the union of the fixture reads minus reads whose primary has `qs < 12` (ROBIN's own rule, which modkit doesn't have), restricted to the panel's CpGs.

**ROBIN:** the cumulative `<sample>.parquet` after the same reads go through `run_matkit` and `merge_modkit_files`, split into 1, 2 and 4 batches.

**Pass, all of:**
1. The same set of CpGs as modkit.
2. Per CpG, `n_mod` and `n_canonical` equal modkit's exactly.
3. `valid_cov == n_mod + n_canonical` exactly.
4. `percent_modified` within 0.005 of modkit's (it prints 2 decimals) and within 1e-4 of `100 * n_mod / valid_cov`.
5. The 1-, 2- and 4-batch results are identical.
6. At each classifier boundary (Sturgeon, NanoDX/PanNanoDX, random forest, MNP-Flex): one row per panel CpG, and the value equals the parquet's percent (on that consumer's scale) within 1e-4; random forest never above 100.
7. MGMT's per-BAM BED (stranded): per strand-site, `Nmod` and `Ncanonical` equal modkit's (strands not combined), and `Nvalid_cov` = their sum.

**Fail:** any one of these. Then the claim is refuted: fix and re-run, or report it.

**Commit order:** this plan (`plan: …`) → code (`code: …`, one commit per task) → `result: … — supported / refuted`, with the comparison's output in the message.

## Global Constraints

- Base `2f29cf7`; branch `fix/batch-5-methylation`; worktree `/home/parlar_ai/dev/robin/.claude/worktrees/batch-5`.
- Tests `/home/parlar_ai/robin-batch3-run/t5.sh` (wrap in `kill-after.sh`); Python `py5.sh`. Baseline 206 passed.
- `git diff 2f29cf7..HEAD -- robin.yml pyproject.toml` stays empty; no existing test weakened.
- Synthetic data only. Public reference data may be fetched (UCSC); never patient data.
- METH-7 (coverage to classifiers), METH-8..12 except what METH-6 needs from METH-12, and the dead `matkit.merge_modkit_files` are out of scope.
- Review document: Status lines, fork rows, Updated line, "Fix first"/step 5 markers, "Already done". Existing learnings are never edited.

## Review Focus

1. **A running sample's existing parquet** (old format: stranded rows, mean percents, fails in coverage): the next merge must migrate it, not mix formats. (Task 4)
2. **The old panel cache** in a sample folder must not be reused after the fix. (Task 3)
3. **A read with MM/ML on some records only**, or `Mm`/`Ml` tags: no crash, no invented calls. (Task 2)
4. **A batch whose BAMs have no panel CpG**: the cumulative file is unchanged and no error is logged (the existing test). (Task 4)
5. **MGMT**: the per-read fixes change its BED; the R script reads `Nvalid_cov` and percent from it. (Task 2)

---

### Task 1: Fixtures and the comparison harness
**Files:** create `tests/meth_fixtures.py`; test `tests/test_meth_vs_modkit.py`.
- `write_reference(dir, cpgs) -> str` (CpG-free random sequence with planted CpGs, faidx'd), `write_meth_bam(path, ref, reads) -> str` (reads given as strand, start/end, per-CpG ML values, MM mode, optional h values, qs, supplementary/primary, sorted + indexed), `modkit_pileup(bam, ref, combine_strands) -> Dict[(chrom,pos,strand), (n_mod, n_can, valid, pct)]` (skips when `modkit` is missing).
- The fixture set from the locked rule, and the comparison helper. Commit with the plan's tests marked xfail-free: they are expected to FAIL on the current code (RED for Tasks 2-5).

### Task 2: matkit per-read layer (METH-3, METH-13, METH-14)
**Files:** `src/robin/analysis/utilities/matkit.py` (`process_bam_counts_improved`, `write_counts_to_parquet`, `write_bedmethyl_improved`).
- `valid_cov = n_mod + n_canonical` and percent from that, in both writers; `Nfail` kept in the BED.
- Implicit MM groups (`.` or no marker): every candidate base the group skips is a canonical call; `?` skips stay no-calls. Parsed from the MM tag text; `Mm` accepted.
- QS: a first pass collects reads whose primary has `qs < 12`; the counting pass skips all their records.
- The 5mC/5hmC rule follows modkit on the m/h-split fixture.
- Tests: per strand-site vs modkit (rule 7); QS order independence; implicit/empty-list cases.

### Task 3: The CpG panel (METH-6, METH-5, the cache)
**Files:** `src/robin/analysis/temp_utilities.py` (panel load and cache); `docs/releases/v0.5.md` (the wrong note).
- Read `parquet_filter.txt` as 0-based BED (no `Start - 1`), drop rows with a negative start, drop duplicate intervals.
- Cache file name carries a format version and the panel's size and mtime; written atomically.
- Tests: no 3 bp widening; the negative row dropped; an old `_1based` cache ignored.

### Task 4: The merge (METH-1, METH-2, METH-4, METH-5 join)
**Files:** `src/robin/analysis/temp_utilities.py` (`merge_modkit_files`).
- Shift `-` rows to `pos - 1`, keep rows inside a panel CpG by semi-join, group by (chrom, pos) summing `n_mod`, `n_canonical`, `valid_cov`; strand `.`; percent from the sums; the same on first write and on every merge.
- An old-format existing file (sidecar `column_format` not the new value) is migrated the same way before merging.
- Tests: locked rules 1-5 against modkit for 1/2/4 batches; the review's 50 reads at 80% + 1 read at 0% → 78.43%; migration of an old-format file.

### Task 5: Classifier boundaries
**Files:** as the tests require (`merge_bedmethyl.py`, `nanodx_analysis.py`, `random_forest_analysis.py`, `sturgeon_analysis.py`, `matkit.reconstruct_full_bedmethyl_for_mnpflex`).
- Tests capture each classifier's input (patch the predictor / R call / upload) from a cumulative parquet built in Task 4, and check locked rule 6. Fix only what fails.

### Task 6: Result and records
- Run the full comparison; commit `result: … — supported / refuted` with its output.
- CHANGELOG; review Status lines for METH-1..6, 13, 14 and the fork table; v0.5 note; a learning; plan Outcome; `trail check`.
