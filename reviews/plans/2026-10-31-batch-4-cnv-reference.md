# Batch 4: a clean CNV control reference (CNV-1, CNV-3) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The CNV control profile is built from the stored control counts only, exactly, and every BAM in a batch is counted once.

**Architecture:** A new numpy function computes the control's profile at the sample's bin width with 64-bit sums, using the same formula as `cnv_from_bam.calculate_cnv`. Both runners (direct and subprocess) use it in place of their second `iterate_bam_file` pass, so the control is never handed to the library again; the cached control becomes read-only so it cannot be. The multi-BAM loop keeps its last result instead of reading the last BAM again.

**Tech Stack:** Python 3.12, numpy, pytest, pysam (synthetic BAMs), `cnv_from_bam` 0.4.4 (installed in the test env from PyPI's wheel; tests that need it use `pytest.importorskip`).

**Outcome:** built as planned in `b682da0` (Task 1), `6ab6928` (Task 2), `a06c4a2` (Task 3), `2681f8b` (Task 4) and `391e159` (Task 5); suite 186 → 206. Differences: Task 1 gained a test for empty input (9, not 8); Task 2's contig-only-in-the-sample test failed on the old code (chrUn was compared with its own reads), so it is a behaviour change, not a guard; Task 4's test runs on a control stored as lists and as arrays (2, not 1). A failing `==` on two 2,001-item lists made pytest spend minutes in difflib, which looked like a hang; the tests compare the names of contigs that differ instead. The wrap was measured to last up to about 24,000 sample reads, not 20,000. Found and left open: in subprocess mode (off by default) a multi-BAM batch keeps only its last BAM. The final review (model fable: ready with fixes, 0 Critical, 1 Important, 6 Minor) found the median-rule gap understated; a full scan measured up to 0.035 copies (at 25,000 sample reads), not 0.024, so the records were corrected, a CHANGELOG sentence fixed, and a new learning added; no code changed, suite still 206.

**Spec:** `reviews/2026-10-16-code-review.md`, CNV-1 and CNV-3 (section 2), and these measurements on `b63d66a`:
- `cnv_from_bam` 0.4.4 `src/lib.rs`: `iterate_bam_file(bam, copy_numbers=d)` merges the BAM's counts into `d` and replaces its values (`merge_and_reset`); `calculate_cnv` sums each contig's 1 kb bins in chunks of `bin_width // 1000` as **`u16`** (release build, no overflow checks, so sums wrap), takes the median of all chunk sums (`f64::round`, half away from zero; 0 becomes 1) and returns `sum / median * 2`.
- Probe (synthetic BAM, current code): after one `run_cnv_analysis_direct`, a control of 3 reads per bin holds 4 — the BAM's reads were added.
- Probe (public control `HG01280_control_new.pkl`, 12.7 M reads, no-BAM mode): at the bin widths of a sample with 4,000 / 20,000 reads, 46/54 and 136/213 control chunk sums exceed 65,535 and the library's profile is off by up to 4.2 and 17.2 copies. At 50,000 reads and more nothing wraps (differences 0.015 at 50,000 and 0.024 at 200,000, from the median rule below; a full scan after the final review found up to 0.035, at 25,000). Today's code wraps the same way, since it sums control + BAM.
- Cost: the 64-bit re-binning takes ~6 ms per call on cached int64 arrays (~100 ms on lists); today's second pass re-reads the whole BAM.

## Global Constraints

- Base `b63d66a`; branch `fix/batch-4-cnv-reference`; worktree `/home/parlar_ai/dev/robin/.claude/worktrees/batch-4`.
- Tests: `/home/parlar_ai/robin-batch3-run/t4.sh`; Python: `/home/parlar_ai/robin-batch3-run/py4.sh <file>`. Baseline `186 passed, 1 warning`.
- `git diff b63d66a..HEAD -- robin.yml pyproject.toml` stays empty. No existing test is weakened or removed.
- Synthetic data only (fake BAMs and fake controls). The public control may be read for measurements, never patient data.
- Review document: only Status lines, fork-table rows, the "Updated" line and the "Already done" paragraph change. Existing learnings are never edited.
- Nothing is merged or pushed.

## Review Focus

1. **A contig only the sample has** (or only the control has): the run succeeds; the normalized profile skips that contig instead of comparing it with itself. (Task 2)
2. **A median that ends in .5** rounds up, as Rust's `f64::round` does, not to even as Python's `round` does. (Task 1)
3. **An all-zero control** gives zeros, not a division by zero. (Task 1)
4. **A bin width under 1 kb** uses chunks of one bin, as the library does. (Task 1)
5. **A control contig whose length is not a multiple of the chunk** ends in a partial chunk, as the library's does. (Task 1)

**Known difference (ruling):** the library chunks the flattened genome for its median, so a chunk can straddle two contigs; this function chunks per contig, as the profile values themselves are. The two medians differ slightly (measured on the public control where nothing wraps: up to 0.035 copies, at 25,000 sample reads; corrected after the final review from "≤ 0.024", which was only the largest of five sampled read counts). Mirroring the library is not possible: its contig order is Rust hash-map order.

---

### Task 1: `reference_cnv_profile`

**Files:** Create `src/robin/analysis/cnv_reference.py`; test `tests/test_cnv_reference.py`.

**Interfaces:**
- Produces: `reference_cnv_profile(control_counts: Mapping[str, Sequence[int]], bin_width: int, bin_size: int = 1000, ploidy: int = 2) -> Dict[str, List[float]]` — keys in natural order (`natsort.natsorted`), values plain lists of floats (the library's output types). Never mutates its input. `chunk = max(1, ((bin_width + bin_size // 2) // bin_size * bin_size) // bin_size)`; per-contig int64 chunk sums (last chunk partial); median of all chunk sums, `floor(m + 0.5)`, 0 → 1; value `sum / median * ploidy`. No contigs → `{}`.

- [ ] **Step 1: tests:**
  - `test_matches_the_library_on_one_contig` — importorskip; one contig of 2,001 bins with varied small counts, bin widths 1,000 / 7,000 / 100,000: equal to `iterate_bam_file(None, copy_numbers=dict(c), bin_width=bw).cnv` (exact).
  - `test_matches_the_library_on_aligned_contigs` — importorskip; three contigs whose lengths are multiples of the chunk: exact.
  - `test_wide_bins_keep_exact_sums` — chr1 40/bin, chr2 20/bin, 2,001 bins, bin width 2 Mb: chr1[0] / chr2[0] == 2.0 and chr1[0] == 80,000 / median * 2.
  - `test_a_half_median_rounds_up` — Review Focus 2: chunk sums [2, 3] → median 2.5 → 3; values [4/3, 2.0].
  - `test_an_all_zero_control_gives_zeros` — Review Focus 3.
  - `test_a_bin_width_under_one_kb_uses_single_bins` — Review Focus 4: bin width 499 → one value per bin.
  - `test_a_partial_last_chunk_is_kept` — Review Focus 5: 5 bins, chunk 2 → 3 values, last from one bin.
  - `test_the_control_is_not_changed` — lists and read-only arrays in, deep-equal after; read-only arrays accepted.
- [ ] **Step 2:** run; all fail (module missing).
- [ ] **Step 3:** implement.
- [ ] **Step 4:** `t4.sh tests/test_cnv_reference.py` → 8 passed; `t4.sh` → 194 passed. Commit.

### Task 2: The direct runner never adds a BAM to the control

**Files:** `src/robin/analysis/cnv_analysis.py` (`run_cnv_analysis_direct`); create `tests/cnv_fixtures.py` (synthetic BAM writer); test `tests/test_cnv_control_reference.py`.

**Interfaces:**
- Consumes: `reference_cnv_profile` (Task 1).
- Produces: `tests/cnv_fixtures.py`: `write_bam(path, contigs: Dict[str, int], reads: Dict[str, List[int]]) -> str` (sorted, indexed, MAPQ 60, 100 M reads at the given 0-based starts) and `flat_control(contigs, per_bin) -> Dict[str, List[int]]` (`length // 1000 + 1` bins, as the library sizes them).
- `run_cnv_analysis_direct` keeps its signature and result keys; `r2_cnv` becomes `reference_cnv_profile(ref_cnv_dict, result.bin_width)`; `timing["pass2_time"]` times that step.

- [ ] **Step 1: tests** (module-level importorskip of `cnv_from_bam`):
  - `test_a_run_leaves_the_control_unchanged` — RED (3 → 4 per bin today).
  - `test_the_reference_profile_is_the_control_alone` — flat control; sample with chr1 at twice chr2's read density: `r2_cnv == reference_cnv_profile(control, r_bin)`, chr1 bins equal chr2 bins. RED.
  - `test_an_early_run_reference_does_not_wrap` — 150 sample reads (2 Mb bins), control chr1 40/bin, chr2 20/bin: `r2["chr1"][0] / r2["chr2"][0] == 2.0`. RED.
  - `test_a_normal_sample_after_an_abnormal_one_is_unchanged` — through `process_single_bam` with the cached control patched to a synthetic one: sample A (chr1 gain), then sample B (normal); B's `CNV3.npy` equals B's on a fresh cache. RED.
  - `test_a_contig_only_the_sample_has_is_skipped` — Review Focus 1: the BAM has `chrUn` too; `process_single_bam` succeeds and `CNV3.npy` has no `chrUn`. Measure on old code first; a guard if it passes.
- [ ] **Step 2:** run on the old code; the first four fail.
- [ ] **Step 3:** implement; the docstring says why the control never goes to `iterate_bam_file`.
- [ ] **Step 4:** `t4.sh` → 199 passed. Commit.

### Task 3: The cached control is read-only

**Files:** `src/robin/analysis/cnv_analysis.py` (`get_cached_ref_cnv_dict`, `run_cnv_analysis_subprocess`); test `tests/test_cnv_control_reference.py`.

**Interfaces:**
- `get_cached_ref_cnv_dict` returns a `types.MappingProxyType` over read-only int64 arrays (`setflags(write=False)`), built once per path. `run_cnv_analysis_subprocess` pickles `{k: v for k, v in control.items()}` (a proxy cannot be pickled).

- [ ] **Step 1: tests:**
  - `test_the_cached_control_cannot_be_changed` — item assignment raises `TypeError`; writing into an array raises `ValueError`. RED.
  - `test_the_cached_control_cannot_be_passed_to_the_library` — importorskip; `iterate_bam_file(bam, copy_numbers=cached)` raises `TypeError`. RED.
  - `test_the_subprocess_runner_still_writes_the_control` — `subprocess.run` patched to return code 1; `ref_cnv_dict.pkl` in the temp dir loads equal to the control. Guard (passes before; must keep passing).
- [ ] **Step 2:** run; two fail.
- [ ] **Step 3:** implement.
- [ ] **Step 4:** `t4.sh` → 202 passed. Commit.

### Task 4: The subprocess script uses the same reference

**Files:** `src/robin/analysis/cnv_subprocess.py` (`run_cnv_analysis`); test `tests/test_cnv_control_reference.py`.

- [ ] **Step 1: test** `test_the_subprocess_script_reference_is_the_control_alone` — call `cnv_subprocess.run_cnv_analysis(bam, ref_cnv_dict_path=<pickled synthetic control>, output_dir=tmp)` in-process: `r2_cnv == reference_cnv_profile(control, r_bin)`. RED.
- [ ] **Step 2:** run; fails.
- [ ] **Step 3:** import `reference_cnv_profile` from `robin.analysis.cnv_reference` and replace pass 2.
- [ ] **Step 4:** `t4.sh` → 203 passed. Commit.

### Task 5: A multi-BAM batch counts each BAM once (CNV-3)

**Files:** `src/robin/analysis/cnv_analysis.py` (`process_multiple_bams`); test `tests/test_cnv_control_reference.py`.

- [ ] **Step 1: test** `test_a_batch_counts_each_bam_once` — two synthetic BAMs (300 and 200 reads) through `process_multiple_bams` (cached control patched): the saved `<sample>_copy_numbers.pkl` sums to 500. RED (700 today). Measure first that the call runs to the save in a temp dir; patch only what it cannot reach (the master BED's reference genome).
- [ ] **Step 2:** run; fails.
- [ ] **Step 3:** the loop keeps the last result as `final_result`; the second "final" call is removed (both modes).
- [ ] **Step 4:** `t4.sh` → 204 passed. Commit.

### Task 6: Records

- CHANGELOG `[Unreleased]` Fixed: the control reference (CNV-1 plus the 16-bit wrap) and the batch double count (CNV-3).
- Review: CNV-1 and CNV-3 Status lines, two fork-table rows, the "Updated" line, step 4 marked done in "Suggested order of work" and "Already done". PDF rebuilt with `/home/parlar_ai/robin-batch3-run/pdf/build.sh <worktree>`.
- A new trail learning (the library mutates `copy_numbers`, wraps at 16 bits, and its no-BAM mode).
- This plan's Outcome line. `trail check --base b63d66a`. Commit.

### After Task 6

Final whole-branch review (fable), then the fix pass, as in `superpowers:executing-plans`.
