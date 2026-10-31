# Batch 3c: the five findings of the 2026-10-31 fix-series review — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fix the five findings the user's review of `493f53f..7ae25ad` raised, each with a test that fails first.

**Architecture:** Five independent tasks, P1s first, then one recording task. Each changes one function (two in Task 1) and adds tests next to the existing ones for that area.

**Tech Stack:** Python 3.12, pytest, numpy, pandas, pysam (synthetic BAMs via `tests/lga_fixtures.py`), NiceGUI 3.18's real `Timer` loop, `/bin/sh`.

**Outcome:** built as planned in `4032f8b` (Task 1), `f3f29cd` (Task 2), `6acc03d` (Task 3), `b26c410` (Task 4) and `ab5151b` (Task 5); suite 139 → 161, as predicted. Three differences: Task 3's first opened-view test re-opened the view after every refresh and passed on the old code, so it was fixed to open once; Task 4's insertion guard measured `reference_mismatch`, not `insertion_requires_review`, so it pins "needs review" and the gene status instead of the label, and indels stay out of the zero-coverage rule because an all-deleted column is evidence for a ClinVar deletion; Task 5's two-level nested test passed on the old code (an odd count of `"` before `{file}` lands right by luck), so it uses three levels. After the final review (model fable: ready, 0 Critical, 2 Important, 6 Minor), the fix pass added `2ca41c7` (a failed CNV detection keeps the last good breakpoints) and `0964752` (no "(legacy location)" caption), and the user chose Task 7 (`65ab945`: sites no read covers need review), which also fixed the site-selection boundary; suite 161 → 170.

**Spec:** `/home/parlar_ai/dev/robin/reviews/2026-10-31-fix-series-review.md` (the user's review; untracked in the main checkout) and its probes in `/home/parlar_ai/dev/robin/reviews/2026-10-31-fix-series-repro.py`. All five reproduced on `7ae25ad` before this plan was written.

## Global Constraints

- Base `7ae25ad`; branch `fix/batch-3c-review-findings`; worktree `/home/parlar_ai/dev/robin/.claude/worktrees/batch-3c`.
- Tests: `/home/parlar_ai/robin-batch3-run/t3c.sh`; Python: `/home/parlar_ai/robin-batch3-run/py3c.sh <file>`. Baseline `139 passed, 1 warning`.
- `git diff 7ae25ad..HEAD -- robin.yml pyproject.toml` stays empty. No existing test is weakened or removed.
- Synthetic data only. The user's review files stay untouched in the main checkout.
- Review document: only Status lines, fork-table rows, the "Updated" line and the "Already done" paragraph change.

## Review Focus

1. **CNV: a run with breakpoints after a flat run** (empty → nonempty) must show the new breakpoints, not an empty leftover. (Task 1)
2. **LGA: a trusted file with a missing or unreadable `generated_at`** must not beat one that has a date, and must not crash the picker. (Task 2)
3. **LGA timer: the analysis still running** — the timer must keep polling until completion or error. (Task 3)
4. **LGA: zero usable coverage at an indel or MNV site** keeps its existing review level; only SNVs change. (Task 4)
5. **Templates: the quoting that works today** (unquoted, double, single, escaped quotes, `$((…))`, a subshell) must not change. (Task 5)

---

### Task 1: An empty CNV breakpoint result replaces the last one (finding 1)

**Files:** `src/robin/analysis/cnv_analysis.py` (the BED writer at ~L911, `save_cnv_files` at ~L973), `src/robin/analysis/master_bed_generator.py` (`_load_bed_file`); test `tests/test_cnv_breakpoint_files.py`.

- [ ] **Step 1: tests**, each saving a 200-bin `chr1` signal with one shift at bin 100 as counter 1 (via `detect_breakpoints_from_cnv` and `save_cnv_files`, as the review's `cnv_probe` does), then a flat 200-bin signal as counter 2:
  - `test_a_flat_run_clears_the_saved_breakpoint_array` — `cnv_data_array.npy` loads as length 0 with field names `("name", "start", "end")`.
  - `test_a_flat_run_writes_an_empty_breakpoint_bed_for_its_counter` — `bed_files/breakpoints_002.bed` exists and is empty.
  - `test_the_master_bed_for_a_flat_run_has_no_old_breakpoints` — `_build_master_bed_data("A", root, 2, None, missing_fa, log)` holds no `(99000, 101000)` region (measure what it returns when there are no regions and assert that).
  - `test_breakpoints_after_a_flat_run_are_shown` — Review Focus 1: flat as 1, shift as 2 → the array has the one breakpoint.
  - `test_an_empty_bed_loads_as_empty_without_a_warning` — `_load_bed_file(empty)` is empty and logs no "Error loading BED file" (caplog).
- [ ] **Step 2:** run on the old code; the first three and the last fail.
- [ ] **Step 3:** always write `breakpoints_NNN.bed` (empty when there are none) and always save the structured array (empty, same dtype). `_load_bed_file` returns an empty frame for a zero-byte file before calling pandas.
- [ ] **Step 4:** `t3c.sh` → 144 passed. Commit.

### Task 2: The newest trusted LGA results file wins (finding 2)

**Files:** `src/robin/analysis/lightweight_gene_analysis.py` (`pick_lga_results`); test `tests/test_lga_sample_inputs.py`.

**Decision:** among trusted candidates, the one whose `metadata.generated_at` is latest wins; a missing or unparseable date counts as oldest; ties keep candidate order (so `test_pick_lga_results_keeps_the_first_of_two_trusted_candidates` still holds unchanged). One canonical output location was the alternative; it would still need old files from the other layout to lose, which this does.

- [ ] **Step 1: tests:**
  - `test_pick_lga_results_prefers_the_newer_of_two_trusted_candidates` — nested `2026-10-30T10:00:00`, flat `2026-10-31T10:00:00` → flat. RED.
  - `test_pick_lga_results_keeps_a_newer_first_candidate` — reversed dates → nested. Guard.
  - `test_a_trusted_file_without_a_date_loses_to_one_with_a_date` — Review Focus 2: nested `generated_at` is `"not a date"`, flat has a real date → flat, no exception. RED.
  - `test_a_fresh_run_in_the_other_layout_is_the_one_shown` — the review's two real handler runs (20 ALT reads with `work_dir=<sample>`, then 20 REF reads with `work_dir=<root>`); `pick_lga_results` returns the flat file, whose status is `none_detected`. RED.
- [ ] **Step 2:** run on the old code; three fail.
- [ ] **Step 3:** implement; docstring says why.
- [ ] **Step 4:** `t3c.sh` → 148 passed. Commit.

### Task 3: The LGA status timer stops itself (finding 3)

**Files:** `src/robin/gui/components/coverage.py` (`check_lga_status` and its `ui.timer` at ~L6080-6107); test `tests/test_lga_status_timer.py`.

- [ ] **Step 1: tests.** Compile `check_lga_status` out of `coverage.py` with `ast` (as the review's `timer_probe` does) into a namespace holding stand-in widgets, a `_check_lga_results` that counts calls, and `status_timer` = a `nicegui.timer.Timer` subclass with `_skip_registration` → `True` and a bound of 5 ticks in `_should_stop`, run with `asyncio.run(timer._run_in_loop())`:
  - `test_completion_refreshes_once_and_stops_the_timer` — status `completed`: 1 tick, 1 refresh, timer cancelled. RED (5 ticks, 5 refreshes).
  - `test_an_opened_results_view_survives_after_completion` — the view opened after the refresh is still open when the loop ends. RED.
  - `test_an_error_stops_the_timer` — status `error`: 1 tick, 0 refreshes. RED.
  - `test_a_running_analysis_keeps_polling` — Review Focus 3: status `running` → 5 ticks (the bound), 0 refreshes. Guard.
  These lean on NiceGUI 3.18's private `_skip_registration`/`_should_stop`/`_run_in_loop` (pinned); say so in the docstring.
- [ ] **Step 2:** run on the old code; three fail.
- [ ] **Step 3:** `status_timer = ui.timer(0.5, check_lga_status, active=True)`; both end branches call `status_timer.cancel()` instead of `return False`; drop the return values.
- [ ] **Step 4:** import check of `robin.gui.components.coverage`; `t3c.sh` → 152 passed. Commit.

### Task 4: An SNV site with no usable bases needs review (finding 4)

**Files:** `src/robin/analysis/lightweight_gene_analysis.py` (`REVIEW_EVIDENCE`, the SNV branch of `_analyze_variant_evidence`); test `tests/test_lga_variant_support.py`.

**Scope:** only sites that have a pileup column with zero usable bases. Sites with no pileup column at all still drop out before export; whether to report them is a product decision left to the user (counting them would turn a gene "insufficient evidence" whenever any ClinVar site in it had no reads).

- [ ] **Step 1: tests** with 20 synthetic `10M1D39M` reads and ClinVar C>T at the deleted position 11 (the review's `zero_coverage_probe`):
  - `test_an_snv_site_with_only_deletions_needs_review` — evidence `no_usable_coverage`, `needs_review` true, `supported` false. RED.
  - `test_a_gene_whose_only_site_has_no_usable_bases_is_not_none_detected` — status `insufficient_evidence`. RED.
  - `test_an_indel_site_with_only_deletions_keeps_its_review_level` — Review Focus 4: a C>CA insertion at position 11 → `insertion_requires_review` (or the level it has today; measure). Guard.
- [ ] **Step 2:** run on the old code; two fail.
- [ ] **Step 3:** add `"no_usable_coverage"` to `REVIEW_EVIDENCE`; the SNV branch checks `total_coverage == 0` first.
- [ ] **Step 4:** `t3c.sh` → 155 passed. Commit.

### Task 5: `{file}` inside `$(…)` is quoted for its own context (finding 5)

**Files:** `src/robin/command_template.py` (`shell_command_for_file`, module docstring); test `tests/test_command_template.py`.

- [ ] **Step 1: tests**, exact arguments through `printf '%s\n'`, with `tmp_path / "a b" / "s.bam"`:
  - `test_a_quoted_placeholder_inside_command_substitution_stays_one_word` — `"$(dirname "{file}")/output.bai"` → `[<tmp>/a b/output.bai]`. RED.
  - `test_nested_command_substitutions_keep_the_path_whole` — `"$(basename "$(dirname "{file}")")"` → `["a b"]`. RED.
  - `test_a_hostile_name_inside_command_substitution_comes_back_exact` — `"$(basename "{file}")"` with the `HOSTILE` name → `[HOSTILE]`, no marker files. RED.
  - Review Focus 5 guards, added to the existing parametrized list: `"$((1+2))" "{file}"` → `["3", PLAIN]`; `( printf '%s\n' "{file}" )` → `[PLAIN]`; ``d=`dirname "{file}"`; printf '%s\n' "$d"`` → `["/d"]`.
- [ ] **Step 2:** run on the old code; the three named tests fail, the guards pass.
- [ ] **Step 3:** replace the single `quote` variable with a stack of contexts: unquoted, `'`, `"`, `$(`, `(` and backtick. `$(` and backtick open a new unquoted context from inside `"` or unquoted; `)` closes `$(` or `(`; a backtick closes a backtick context. `{file}` is emitted for the innermost context: `"${ROBIN_FILE}"` in any unquoted one, `${ROBIN_FILE}` in `"`, `'"${ROBIN_FILE}"'` in `'`. Backslash escapes as today. The file name still never becomes shell text. Docstring lists the handled contexts and the limits (`case` patterns inside `$(…)`, here-documents).
- [ ] **Step 4:** `t3c.sh` → 161 passed. Commit.

### Task 7: A ClinVar site with no reads needs review (added after the final review; the user chose this)

**Files:** `src/robin/analysis/lightweight_gene_analysis.py` (site selection ~L617, `_find_intersecting_variants`, `REVIEW_EVIDENCE`, the export summary, `_analyze_single_gene`); test `tests/test_lga_variant_support.py`.

**Decision (user, 2026-10-31):** a ClinVar site inside a checked region that no read covers is kept, with evidence `no_reads` (a review level), so the gene reads "insufficient evidence"; each gene summary gets `sites_without_reads`. **Boundary:** sites are picked with `start <= pos <= end` (1-based POS against the 0-based BED start), so a site at `pos == start` is the base before the region and is never piled up; it must become `start < pos <= end`, or it would be flagged `no_reads` falsely. The unused `_analyze_single_gene` path counts only supported sites, like the region path.

- [ ] **Step 1: tests** (BED `(0, 100)`, 20 reads covering 1-based 1–50, ClinVar SNVs at 11 and/or 80):
  - `test_a_site_with_no_reads_needs_review` — only site 80: evidence `no_reads`, `needs_review`, gene `insufficient_evidence`, `sites_without_reads == 1`. RED (gene absent).
  - `test_a_reference_only_site_and_an_unread_site_is_not_none_detected` — REF reads, sites 11 and 80 → `insufficient_evidence`. RED (`none_detected`).
  - `test_a_supported_site_still_wins_over_an_unread_one` — ALT reads, sites 11 and 80 → `present_with_good_coverage`, `sites_without_reads == 1`. RED (no count).
  - `test_a_site_just_before_the_bed_start_is_not_checked` — BED `(10, 60)`, site at POS 10 → not in the gene's variants. Guard (passes today; fails if the boundary is left).
  - `test_analyze_genes_counts_only_supported_sites` — the unused per-gene path, REF reads at site 11 → `pathogenic_snps_found == 0`. RED (counts 1).
- [ ] **Step 2:** run on the old code. **Step 3:** implement. **Step 4:** `t3c.sh` → 170 passed. Commit.

### Task 6: Record

- [ ] `CHANGELOG.md` `[Unreleased]` → `### Fixed`: one line per finding, saying what users notice.
- [ ] Review document: Status lines of CNV-5, VAR-9, VAR-8 and SEC-9 gain a follow-up sentence with the commit; fork-table rows; an Updated line; "Already done".
- [ ] `.gitignore:178` trailing space (the review's note).
- [ ] Rebuild the PDF; the plan's Outcome paragraph; `trail check`; a `trail note` if something took digging; commit.
- [ ] Final whole-branch review on model `fable`; then stop for the user.
