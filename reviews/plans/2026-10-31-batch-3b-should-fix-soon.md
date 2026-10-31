# Batch 3b: the six "should fix soon" follow-ups — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the six items the batch 3 final review marked "should fix soon".

**Architecture:** Four small, independent tasks plus one recording task. Two harden existing tests (the GUI timer test, and the NiceGUI stub the refresh test borrows from `reviews/`); one rewrites how `robin add-panel` reads its BED; one makes the LGA panel's "Hide" re-run the results check instead of rebuilding the view with a hard-coded path.

**Tech Stack:** Python 3.12, pytest, click's `CliRunner`, `ast`.

**Spec:** `/home/parlar_ai/robin-batch3-run/final-review.md` — the "should fix soon" rows of the deferred-minor triage (#8, #20, #21, #23), the "Known residues" rows for the 8-column BED and `_hide_lga_results`; and `reviews/2026-10-16-code-review.md` VAR-7 **Status** line (the 8-column problem, "still open").

## Global Constraints

- Base `6e0a842` on `main`; branch `fix/batch-3b-should-fix-soon`; worktree `/home/parlar_ai/dev/robin/.claude/worktrees/batch-3b`.
- Tests run with `/home/parlar_ai/robin-batch3-run/t3b.sh` (the batch 3 venv, `PYTHONPATH=src`). Baseline: `131 passed, 1 warning`.
- `git diff 6e0a842..HEAD -- robin.yml pyproject.toml` stays empty.
- No existing test is weakened or removed: `git diff 6e0a842..HEAD -- tests/ | grep '^-[^-]'` prints only lines a task below names and replaces with a stricter one.
- Synthetic data only.
- In the review document change only Status lines, fork-table rows, Fix-first / order-of-work marks and the "Already done" paragraph.
- Each new test is seen failing on the old code (scratch copy via `git archive 6e0a842 src tests | tar -x`), or is named a guard and why.

## Review Focus

1. **A standard BED12 panel** (UCSC export) — today `add-panel` exits 0 and writes chromosome `100`, gene `1`. It must write the panel's real 4 columns. (Task 3)
2. **A file mixing 4-column and 12-column rows** — the validator accepts it; today processing fails with a pandas tokenizer error. (Task 3)
3. **Mixed 4/6-column rows, CRLF line endings, a comment line between rows** — all work today (measured); the rewrite must not break them. (Task 3, guard tests)
4. **The NiceGUI stub missing** — the refresh test's subprocess must fail loudly, not fall back to the real NiceGUI. (Task 2)
5. **"Hide" after the results file stopped being trusted** — it must not offer "View" for that file. (Task 4: re-running `_check_lga_results` already handles it.)

---

### Task 1: Harden the GUI timer test (#20, #21)

**Files:** Modify `tests/test_gui_timer_cleanup.py`.

- [ ] **Step 1:** `test_no_refresh_timer_is_torn_down_on_disconnect` scans every `*.py` under `src/robin/` (sorted `rglob`), not the `GUI_FILES` list, and first asserts every `GUI_FILES` entry is among the scanned files. `offenders == []` as before.
- [ ] **Step 2:** `test_every_file_that_had_one_now_uses_on_delete` asserts `len(GUI_FILES) == 9`, and per file `"client.on_delete(" in text` (was `"on_delete" in text`, which `on_deleted` also matches).
- [ ] **Step 3: RED by mutation**, in a scratch copy: (a) add `ui.context.client.on_disconnect(lambda: None)` to `src/robin/gui/components/snp.py` (not in `GUI_FILES`) — old test passes, new test fails; (b) replace `summary.py`'s `ui.context.client.on_delete(lambda: refresh_timer.deactivate())` with `pass  # on_deleted` — old test passes, new test fails.
- [ ] **Step 4:** `t3b.sh` → 131 passed (no new tests). Commit.

### Task 2: The NiceGUI stub lives under `tests/` (#8)

**Files:** Create `tests/_stubs/nicegui/__init__.py` (verbatim copy of `reviews/2026-10-16-repro/gui/stubs2/nicegui/__init__.py`); modify `tests/test_classification_refresh.py`, `reviews/2026-10-16-repro/README.md:13`.

`stubs2/` stays: `reviews/2026-10-16-repro/gui/test_mgmt_shared_state.py` still uses it. Do not name the folder `tests/nicegui/`: `tests/` is on `sys.path` under pytest and that would shadow the real package.

- [ ] **Step 1:** `_STUBS = _REPO_ROOT / "tests" / "_stubs"`. `_stubbed_env()` asserts `(_STUBS / "nicegui" / "__init__.py").is_file()` so a missing stub fails loudly (Review Focus 4).
- [ ] **Step 2: RED:** in a scratch tree holding only `src/` and `tests/` (no `reviews/`), the old test fails; the new one passes. Record the old failure's actual message.
- [ ] **Step 3:** README line 13: the sentence about `stubs2/` and `tests/test_classification_refresh.py` now says the test has its own copy under `tests/_stubs/`.
- [ ] **Step 4:** `t3b.sh` → 131 passed. Commit.

### Task 3: `add-panel` reads every BED line the way the validator does (#23 and the 8-column defect)

**Decision:** columns after the 4th are allowed and ignored. The validator already accepts them, the docstring says "at least 4 columns", and BED12 is a standard format.

**Files:** Modify `src/robin/cli.py` (`_generate_unique_gene_bed`, `add_panel` docstring); test `tests/test_add_panel_regions.py`.

- [ ] **Step 1: tests.** A helper `_run_add_panel(tmp_path, monkeypatch, bed_text, name="P")` writes the BED as bytes, points `robin.resources.__file__` at `tmp_path / "res" / "__init__.py"` (so nothing is written into the package), sets `robin.cli._get_user_acknowledgment` to return `True`, runs `CliRunner().invoke(main, ["add-panel", str(bed), name])`, and returns `(result, lines of res/<name>_panel_name_uniq.bed or None)`. Expected lines are always `["chr17\t300\t400\tTP53", "chr7\t100\t200\tEGFR"]` (sort is by chromosome as text). Tests, each asserting `exit_code == 0` and those lines:
  - `test_the_command_writes_a_4_column_panel` — guard (passes today; the command-level test #23 asked for).
  - `test_an_8_column_bed_is_read_by_its_first_4_columns` — RED today: exit 1, `int('EGFR')`.
  - `test_a_bed12_panel_is_not_silently_garbled` — RED today: exit 0, one line `100\t200\t0\t1`.
  - `test_rows_with_4_and_12_columns_mix` — RED today: exit 1, tokenizer error.
  - `test_mixed_4_and_6_columns_crlf_and_comments_still_work` — guard (works today).
  - `test_a_failed_panel_leaves_no_files` — a row whose 4th column is empty: exit 1, and `res/` holds no `P_panel_name_uniq.bed` and no `P_panel_source.bed`. Guard.
- [ ] **Step 2:** run them on the old code; the three RED tests fail for the reasons above.
- [ ] **Step 3:** in `_generate_unique_gene_bed`, replace both `read_csv` calls and the `except ValueError` fallback with a line loop: `strip()`, skip blank and `#` lines, `split('\t')`, take `parts[0]`, `int(parts[1])`, `int(parts[2])`, `parts[3]`. Keep the comma split, dedupe, sort and 4-column write as they are. Update the docstring's "Supported formats" to say extra columns are ignored.
- [ ] **Step 4:** `t3b.sh` → 137 passed (131 + 6). Commit.

### Task 4: "Hide" re-runs the results check (A·2 Minor #3)

**Files:** Modify `src/robin/gui/components/coverage.py` (`_hide_lga_results`); create `tests/test_lga_hide_results.py`.

- [ ] **Step 1: test** `test_hide_re_runs_the_results_check`: `ast`-parse `coverage.py`, find the `FunctionDef` named `_hide_lga_results`, assert it contains a call to `_check_lga_results` and no string constant `"lightweight_gene_analysis_results.json"`. RED today (hard-coded path, no call). Lint-shaped: GUI code has no harness; say so in the commit.
- [ ] **Step 2:** body becomes `_check_lga_results()` with a docstring saying why. `_check_lga_results` picks the trusted file via `pick_lga_results` (tested in `tests/test_lga_sample_inputs.py`) and renders the same summary view through `_display_lga_results`.
- [ ] **Step 3:** `python -c "import robin.gui.components.coverage"` loads; `t3b.sh` → 138 passed. Commit.

### Task 5: Record

- [ ] `CHANGELOG.md` `[Unreleased]` → `### Fixed`: the BED12 / extra-columns fix (what users notice: a BED12 panel used to be saved garbled with a success message).
- [ ] Review VAR-7 Status line: the second problem is fixed, with the commit, and that BED12 was silently garbled, not only an error. README row 40 likewise.
- [ ] Rebuild the PDF: `bash /home/parlar_ai/robin-batch3-run/pdf/build.sh /home/parlar_ai/dev/robin/.claude/worktrees/batch-3b`.
- [ ] This plan gets an "Outcome" paragraph.
- [ ] `trail check`; a `trail note` if something took digging; commit.
- [ ] Final whole-branch review by a fresh reviewer (model `fable`), then stop for the user's merge decision.
