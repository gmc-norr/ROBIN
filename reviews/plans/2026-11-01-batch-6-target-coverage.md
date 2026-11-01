# Batch 6: target coverage, read counts and SNP input match samtools (step 6) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** every per-sample coverage number, read count and the SNP-input BAM (`target.bam`) equal what samtools gives when run once on all the sample's reads in one BAM, however the reads were split into BAMs, in whatever order they arrived, and whether one or two workers ran.

**Architecture:** each BAM yields exact integer totals (aligned bases per contig and per target, read counts) and its covered stretches; the accumulation adds them to the sample's totals in one locked read-add-write step and unions the stretches; every BAM's panel reads go to a batch BAM when it arrives, and the final merge keeps everything already in `target.bam`. Consumers read depth from shared helpers.

**Tech Stack:** Python 3.12, pysam 0.24.1 (bundles samtools 1.24: `coverage`, `bedcov -j`, `view -L`, `stats`, `merge`), pandas, numpy.

**Spec:** `reviews/2026-10-16-code-review.md` section 3, "Target coverage, read counts and SNP input" (TGT-1 to TGT-5), and the user's answers on 2026-11-01 (Gate A): the principle above ("Same as one samtools run"); exact covered bases from a list of covered stretches (not D4; D4 per-base panel depth goes on the later list as a separate feature).

**Spec notes (measured on `db7269c`, synthetic BAMs, scratchpad `b6/`):**
- Depth adds up across BAMs: per-target `bedcov` sums over 3 BAMs equal the merged BAM's, with and without `-j`. `bedcov` counts deletions by default (2724 vs 2708 with `-j`); `samtools coverage` does not (`50M10D50M` → 100 covered bases, meandepth 0.1 on 1 kb). Contig bases from `bedcov -j` over the whole contig, divided by its length, equal `coverage`'s meandepth (5.185000 vs 5.185; 5.269333 vs 5.26933, which `coverage` prints to 6 significant digits).
- Covered bases don't add up: summed per BAM 4,048 vs 1,669 merged (c1). The union of each record's aligned blocks (`get_blocks()`, records not unmapped/secondary/QC-fail/duplicate) over 3 BAMs equals merged `covbases` (926 = 926). `coverage` counts supplementary records, not secondary, not deletions or `N`.
- A list of covered stretches stays small: simulated hg38-sized genome, lognormal reads (mean 8.3 kb): 43,722 stretches at 50k reads, 131,863 at 500k (~1.3×, ~2 MB), 9,485 at 2M, 10 at 5M.
- Read counts: primary mapped `view -c -F 0x904` and unmapped `-f 4` add up across BAMs (100, 20). Today's rule (`not is_secondary`) counts 120 mapped (supplementary included). Primary records' `query_length` sums equal `samtools stats` "bases mapped" (18,774) and "total length" (22,211).
- `samtools view -L` writes each record once even with overlapping and duplicated BED rows, with or without an index; it keeps supplementary and secondary records that overlap the regions. The rCNS2 panel BED has 243 rows, 4 columns, no exact duplicates.
- Today (from the code): genome numbers come from `pysam.coverage`, target numbers from `pysam.bedcov` (no `-j`); totals are read under the accumulation lock (`target_analysis.py:1798-1814`), the lock is released, and they are written under a second lock (`1909-2009`); filtered batch BAMs are made only once some target passes `callthreshold` (`1855`), for that batch's BAMs only; `finalize_accumulation_for_sample` deletes new batch BAMs without merging them when `target.bam` exists (`2603-2618`); target length is `end - start + 1` on a 0-based BED (`1830`, `1918`, and GUI `coverage.py:2514, 6192, 6274, 6440`); global coverage is Σcovbases/Σendpos (`1817`, report `coverage.py:70`, GUI `summary.py:1162`, `coverage.py:6663`); `cnv_handler` doesn't pass its panel on, so `save_cnv_files` takes it from fusion metadata or None (`cnv_analysis.py:1086`); `bam_preprocessor.py:352` counts supplementary records as reads.
- Learnings that bear on this: an empty result must still be written when a reader picks the newest file (master BED); LGA reads `igv/igv_ready.bam`, built from `target.bam`.

## Locked metric and decision rule (Gate C — committed before any code)

**Fixture** (`tests/target_fixtures.py`, synthetic): a reference with `chr1` (3,000 bp) and `chr2` (2,000 bp); a panel BED with three targets on `chr1`, two of them overlapping, and one on `chr2` that no read reaches; about 40 reads: plain, with deletions, with insertions and soft clips, split reads with the primary in the panel and the supplementary off it, split reads with the primary off the panel and the supplementary in it, a primary with a secondary in the panel, unmapped reads, reads overlapping each other, a read spanning two targets. A read's records stay in one BAM.

**Runs:** the reads split into 1, 2 and 4 BAMs, each processed through the real entry points (`process_file_with_staging` + `accumulate_staged_files`, `bam_preprocessor`'s read counting, `finalize_accumulation_for_sample`) in five ways: (a) one accumulation per BAM in order; (b) the same in reverse order; (c) all staged, then one forced accumulation; (d) two processes accumulating at the same time; (e) finalize after half the BAMs, then the rest, then finalize again.

**Reference:** samtools (pysam) run once on the merged BAM:
1. Per target: `bases` = `bedcov -j`; depth = bases / (end − start); reads = `view -c -F 0x904 merged.bam <target>`.
2. Per contig: `numreads` and `covbases` from `samtools coverage`; mean depth = `bedcov -j` over the whole contig / contig length.
3. Global depth = Σ contig bases / Σ contig length.
4. Sample counts: mapped reads = `view -c -F 0x904`; unmapped = `view -c -f 4`; mapped bases = `stats` "bases mapped"; yield = `stats` "total length".
5. `target.bam`'s records = `view -L <panel> merged.bam`, compared as a multiset of (name, flag, contig, position, CIGAR).

**Pass, all of:**
1. Rules 1, 2 and 4: every integer equal; depths within 1e-9 × max(1, value) of the reference ratio.
2. Rule 3: `coverage_main.csv`'s global depth and the last row of `coverage_time_chart.npy` within 1e-9 × max(1, value).
3. Rule 5: equal multisets, after the final finalize of every run.
4. Every run (k = 1, 2, 4; ways a–e) writes the same totals and the same `target.bam` records.
5. At the consumers (the PDF report's coverage section; the GUI summary's and coverage card's global and target depth): the numbers equal rules 1–3.
6. With no fusion metadata and a master BED holding only off-panel CNV breakpoints, rule 5 still holds; and `cnv_handler`'s master BED request carries the configured panel.

**Fail:** any one of these. Then the claim is refuted: fix and re-run, or report it.

**Commit order:** this plan (`plan: …`) → code (`code: …`, one commit per task) → `result: … — supported / refuted`, with the comparison's output in the message.

## Global Constraints

- Base `db7269c`; branch `fix/batch-6-target`; worktree `/home/parlar_ai/dev/robin/.claude/worktrees/batch-6`.
- Tests `/home/parlar_ai/robin-batch3-run/t6.sh` (wrap in `kill-after.sh`); Python `py6.sh`. Baseline 229 passed.
- `git diff db7269c..HEAD -- robin.yml pyproject.toml` stays empty; no new dependency; no existing test weakened.
- Synthetic data only; public reference data may be fetched; never patient data.
- Out of scope: D4 per-base depth (later list); recovering orphaned `_batch_*` folders after a crash; ClairS-TO itself (docker; tested at its input, `target.bam`); the dead `process_file` / `process_single_file` path and `ReadBam.py` (no callers).
- Review document: Status lines, fork rows, Updated line, step 6 marker, "Already done", D4 on the later list. Existing learnings are never edited.

## Review Focus

1. **A sample folder written by today's code** (no `bases` column in `coverage_main.csv`), opened in the GUI or report, or continued by a new BAM: depth still shown (from meandepth), covered bases shown as not available, no crash. (Tasks 4, 7)
2. **A staging set half-written** when another worker claims files: no coverage file paired with another BAM's source path. (Task 3)
3. **A BAM with no panel reads**, and a panel contig missing from the BAM's header: empty batch BAM or none, zeros in the totals, no exception. (Tasks 3, 5)
4. **A crash between replacing `target.bam` and deleting the merged batch BAMs**: the next finalize must not duplicate reads. (Task 5)
5. **A panel with overlapping targets**: each record once in `target.bam`; per-target totals per row. (Tasks 1, 4)

---

### Task 1: Fixtures and the samtools comparison (RED)
**Files:** create `tests/target_fixtures.py`, `tests/test_target_vs_samtools.py`.
- `write_reference(dir) -> str`, `fixture_reads() -> List[List[AlignedSegment-spec]]`, `write_bams(dir, k) -> List[str]` (sorted, indexed, a read's records together), `write_panel(dir) -> str`, `samtools_reference(merged_bam, panel, contigs) -> dict` (rules 1–5).
- `run_robin(work, bams, way) -> dict` drives the real entry points for ways a–e and reads back `coverage_main.csv`, `bed_coverage_main.csv`, `target_coverage.csv`, `coverage_time_chart.npy`, the read-count cache, the summed `bam_preprocessor` counters and `target.bam`.
- Tests: `test_totals_match_samtools[k-way]` (rules 1–3), `test_read_counts_match_samtools[k]` (rule 4 part), `test_target_bam_matches_samtools_view_L[k-way]` (rule 3), `test_every_run_gives_the_same_result` (rule 4). Expected to FAIL on today's code.

### Task 2: Read counts count primary records (TGT-5)
**Files:** `src/robin/analysis/bam_preprocessor.py:352`.
- Mapped/unmapped reads, bases, pass/fail splits and `yield_tracking` count primary records only (`not is_secondary and not is_supplementary`); the `supplementary_reads` counter stays. The 50,000-read gate then counts primary records too.
- Test: `test_read_counts_match_samtools[k]` turns green; a one-read + one-supplementary BAM counts 1 mapped read.

### Task 3: Per-BAM totals and covered stretches (TGT-3)
**Files:** create `src/robin/analysis/target_coverage.py`; modify `target_analysis.py` (`process_file_with_staging`).
- `contig_totals(bam) -> DataFrame[#rname, startpos, endpos, numreads, bases]` (`coverage` for numreads, `bedcov -j` over each contig for bases); `target_bases(bam, bed) -> DataFrame[chrom, startpos, endpos, name, bases]` (`bedcov -j`, one row per BED row); `covered_stretches(bam) -> DataFrame[chrom, start, end]`; `union_stretches(frames) -> DataFrame`; `covered_bases(stretches) -> Dict[str, int]`.
- Staging writes `bedcov_N`, `stretches_N`, `timestamp_N`, `source_bam_N`, then `coverage_N.parquet` last by atomic rename; a claim takes only counters whose five files exist.
- Tests: the helpers against samtools on the fixture; a half-written set is not claimed.

### Task 4: One locked accumulation (TGT-4, TGT-3, TGT-1 filtering)
**Files:** `target_analysis.py` (`accumulate_staged_files`, `run_bedmerge` callers).
- Outside the lock: load and sum the claimed sets; per source BAM, per-target read counts and a filtered BAM (`view -L` on the panel ∪ the master BED, every BAM, no threshold) merged into a uniquely named `batch_*.bam` written by atomic rename.
- One locked section: read the totals and the stretches, add the batch, recompute covbases from the stretches, meandepth = bases / length, target depth = bases / (end − start), write every file atomically, append the time rows, write `targets_exceeding_threshold.bed`, delete the batch folder.
- An old-format `coverage_main.csv` (no `bases`): bases = round(meandepth × length) to continue; covered bases become not available for that sample (no stretches to start from).
- Tests: locked rules 1–4 for ways a–d; two processes in parallel give the review's 200 bases, not 100; migration of an old folder.

### Task 5: `target.bam` keeps every batch (TGT-1)
**Files:** `target_analysis.py` (`finalize_accumulation_for_sample`).
- Under the accumulation lock: merge the existing `target.bam` and every `batch_*.bam` into a temp file, replace `target.bam` atomically, then delete the merged batch BAMs. The block that deletes batch BAMs when `target.bam` exists goes. A crash before the deletes must not duplicate reads on the next run (Review Focus 4): record the merged batch names in `target.bam`'s header (`@CO`) and skip any already listed.
- Tests: way e; a re-run after a simulated crash keeps each record once.

### Task 6: The panel always reaches the extraction regions (TGT-2)
**Files:** `cnv_analysis.py` (`cnv_handler`, `process_multiple_bams`, `save_cnv_files`); `target_analysis.py` (regions = panel ∪ master BED, Task 4).
- `cnv_handler` passes its `target_panel` down; `save_cnv_files` uses it, falling back to fusion metadata. The master BED cache key is left as is: the union in target analysis makes an old panel-less master BED harmless (ruling recorded).
- Tests: locked rule 6.

### Task 7: Consumers
**Files:** `reporting/sections/coverage.py`, `gui/components/summary.py`, `gui/components/coverage.py`, `reporting/plotting.py`; helpers in `target_coverage.py`: `target_length(df)`, `target_depth(df)`, `global_depth(cov_df)`, `covered_bases_known(cov_df)`.
- Every place that computed global coverage from covbases, or target length as end − start + 1, uses the helpers. Covered bases unknown (old folder) show "not available".
- Tests: locked rule 5 (report section values; the GUI's computing functions); a scan test that no `+ 1` target length or covbases-based global coverage remains in `src/robin/gui` and `src/robin/reporting`.

### Task 8: Result and records
- Run the full comparison; commit `result: … — supported / refuted` with its output.
- CHANGELOG; review Status lines for TGT-1..5 and the fork table; D4 per-base panel depth on the later list; a learning; plan Outcome; `trail check`.
