# Batch 3: Quick Wins Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Outcome (2026-10-31):** done on `fix/batch-3-quick-wins` (`50da100`..`849e078`), not yet merged. All nine tasks were built; the suite went 80 → 131. Every task review came back spec-compliant; four tasks needed one fix round each.

**This plan was reviewed before execution and four of its defects would have produced a false pass.** All four were corrected in the text above before Task 1 started:

- **Task 7's tabix test could not detect the bug.** Positions 100 and 200 share a 16 kb tabix bin, so a *stale* index still finds the new record — the test would have passed before and after. It now moves the record to position 5,000,000 and asserts the index file changed.
- **Task 4 had an untested crash.** `int(None)` at `summary.py:650` raises today and an import check never reaches it, so widening the defaults to `None` could have shipped a broken summary page with a green suite. A stub-`ui` test now reddens for exactly that.
- **Task 6's own verification step was impossible.** `robin.gui.components.mnpflex` cannot be imported in the test environment (`alive_progress` is missing and installing is forbidden), so it is byte-compiled instead and the other eight are imported.
- **Task 1's tests could not fail for GUI-1's reason.** The `KeyError` comes from the `is_fresh_visit` comprehension, not the first loop, so a half fix would have gone green *and* deleted the reproduction script. An end-to-end subprocess test was added, and retiring the script was made conditional on it.

What else was built differently:

- **Tasks 1 and 2 extract a helper rather than patching in place**, because the inline fix is two untestable lines inside a NiceGUI closure and a `try/except Exception: continue`.
- **REP-1's "raise on an unknown classifier" is scoped to a new resolver only.** Making the threshold lookups strict would break the summary page: `summary.py` passes section anchors (`coverage`, `cnv`, `mgmt`, `fusion`) and the literal `"default"`, all of which need the fallback.
- **GUI-2 was fixed at nine call sites, not the eight the finding lists** — `gui_launcher.py` is a ninth with the identical defect.
- **GUI-4 shows "Not available" rather than a fabricated `0`.** The finding gives no fix line; stopping at "remove the counter from the list" would leave the card asserting "Gained: 0" for every sample, which is the same class of defect as "Gained: 37".
- **JOB-6 is only partly fixed, and says so.** The constant and its configurability are done; the finding's other two consequences overlap JOB-7 and are outside the review's own quick-win scope. Its Fix-first mark reads *(partly fixed in our fork)*.
- **Task 3's brief demanded an unsatisfiable grep invariant** (a match count that the new function definition itself increments) and **Task 7's demanded a flaky assertion** (a `.tbi` mtime comparison that fails about a third of the time, because the filesystem stores mtime to ~1 ms). Both were reported with measurements rather than worked around.
- **A new defect was found and deliberately not fixed:** an 8-column BED makes `_generate_unique_gene_bed` return `False`, so `robin add-panel` rejects any BED with more than 6 columns *after* calling it valid. No review finding covers it; it is recorded in VAR-7's Status line.
- **The final whole-branch review found one Important defect in this work** — a forced tabix rebuild that *failed* left the stale index on disk for the next start-up to accept as good, reviving VAR-3 behind an error from an earlier run. Fixed in `849e078`.

**Goal:** Fix the eight findings in step 3 ("Quick wins") of the review's *Suggested order of work* — GUI-1, REP-1, JOB-6, GUI-4, CNV-5, GUI-2, VAR-3, VAR-7 — each as its own change with its own regression test.

**Architecture:** Eight independent one-to-few-line fixes in eight different subsystems. Nothing passes between them. Where the defective code sits inside a NiceGUI closure (GUI-1) or a GUI render path (GUI-4), the decision is lifted into a module-level pure function that a test can call, and the GUI keeps only the wiring — this repo has no NiceGUI test harness, so logic left inside a closure can only ever be covered by an import check.

**Tech Stack:** Python 3.12, pytest, pandas 2.1.4, numpy 1.26.4, pysam, ruptures, NiceGUI 3.18 (wiring only), click.

**Spec:** `reviews/2026-10-16-code-review.md` in the worktree. Each task quotes its finding in full; read the finding itself with `grep -n '^\*\*GUI-1 ' reviews/2026-10-16-code-review.md` and the following lines. Reproduction scripts for five of the eight live in `reviews/2026-10-16-repro/` — read the script for your task; it shows the bug on synthetic data.

## Global Constraints

- **Worktree:** `/home/parlar_ai/dev/robin/.claude/worktrees/batch-3`, branch `fix/batch-3-quick-wins`, base `eccbf18`. Every git command is `git -C /home/parlar_ai/dev/robin/.claude/worktrees/batch-3 …`; a bare `git` runs in the MAIN checkout and would commit to the wrong tree. Run `git -C … branch --show-current` before any commit; it must print `fix/batch-3-quick-wins`.
- **The main checkout at `/home/parlar_ai/dev/robin` stays untouched.** Write files only inside the worktree.
- **Tests:** `/home/parlar_ai/robin-batch3-run/t3.sh [pytest args]` (no args = all of `tests`). Suite at the start of this plan: **80 passed, 1 warning**. The one warning is a pre-existing polars `StringCache` deprecation; a new warning is a defect.
- **Each task's report states the new passed count and which tests account for the difference.** Never adjust an expected number to make the arithmetic close; find the term that changed and name it.
- **A verification that reddens nothing is a finding, not a pass.** Every new test must be seen failing on the old code, for the reason it targets, with the output captured verbatim. A test that passes before the fix is a *guard*: say so plainly and either strengthen it or record why it stays.
- **False-pass signature of this suite:** a test file whose import fails prints `N error(s) during collection` and `Interrupted`, with **no** `passed` line and no `FAILED`. Judge on the `passed` count, never on the absence of `FAILED`.
- **Do not create a venv and do not install anything.** The environment at `/home/parlar_ai/robin-batch3-run/venv` already has pandas 2.1.4, numpy 1.26.4, pysam 0.24.1, polars, nicegui 3.18.0 and ruptures 1.1.10. `bcftools` is on PATH; a standalone `tabix` binary is **not**.
- **`robin.yml` and `pyproject.toml` must not change.** `git diff eccbf18..HEAD -- robin.yml pyproject.toml` stays empty.
- **No existing test weakened, loosened or removed.** All 80 stay green.
- **Synthetic data only.** No patient data, and never go looking for any.
- **Scripted multi-line edits** (python/sed) must use anchors occurring exactly once, with the match count asserted, and `wc -l` plus `git diff --stat` checked afterwards. An earlier batch duplicated ~830 lines this way. Prefer ordinary single-hunk edits.
- **Never** `git checkout --`, `git reset --hard`, `git clean`, or `git stash` (the stash is shared with other sessions). To set work aside, make a WIP commit.
- **One commit per task**, message naming the finding. Every commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **No merge, push, PR, tag, or remote change. Ever.**
- **Reproduction scripts:** when every finding a script shows is fixed, its pytest regression test replaces it — delete the script and update the table in `reviews/2026-10-16-repro/README.md`. If the script also shows findings that are still open, keep it and note the fixed one in the table. Which is which is stated in each task. Note that `reviews/2026-10-16-repro/independent/reproduce.py` is inverted: there `OK` means the bug still reproduces and `FAIL` means it is fixed; it currently fails ir01, ir02, ir08 and ir09.

## Review Focus

1. **A classifier name REP-1's new resolver has never seen** must fail loudly rather than silently scoring against the 80/50 default — but `get_confidence_level` and `get_confidence_ui_tier` are deliberately called with non-classifier anchors (`"coverage"`, `"cnv"`, `"mgmt"`, `"fusion"`, and the literal default `"default"`), which must keep falling back to `DEFAULT_CONFIDENCE_THRESHOLDS`. A fix that makes those raise takes the whole summary page down. Pinned in Task 2.
2. **`ROBIN_COMPLETION_TIMEOUT_SECONDS` set to garbage, empty or negative.** A sample is marked Complete and finalised once this many seconds pass with no new BAM, and finalisation deletes unmerged batch BAMs. A value that silently becomes `0` would finalise every sample on the first tick. Pinned in Task 3.
3. **A ClinVar `.gz` replaced while a valid-looking old `.tbi` is still on disk** — including on the download-if-missing path, not only the update path. This is the whole of VAR-3: `_file_ok(tbi)` is true, so the old index survives and region queries return nothing with exit code 0. Pinned in Task 7.
4. **A panel listing the same gene on two chromosomes, or the same region twice.** VAR-7's fix must keep every *distinct* region while not reintroducing exact duplicate rows. Pinned in Task 8.
5. **A sample with no classifier enabled at all** (`charts` empty) and **a chromosome whose only change point is the series end** (`result[:-1]` empty). Both are the degenerate end of the collection the fix now iterates; neither may crash or emit a region. Pinned in Tasks 1 and 5.

---

### Task 1: GUI-1 — the classification panel stops updating when any classifier is disabled

> **Finding (`reviews/2026-10-16-code-review.md`, `^\*\*GUI-1 `):** critical, Confirmed (spot-checked). "The classification panel stops updating when any classifier is disabled. **Where:** `gui/components/classification.py:1003-1021`. **What's wrong:** The loop goes over all four tools in `tool_to_file`, but `charts` only has entries for the enabled ones. So `charts['PanNanoDX']` raises KeyError on every tick after the first. The panel keeps showing the first result, while the summary card moves on. **Fix:** in both loops of `_refresh_classification` (lines 1003-1008 and 1019-1021), loop over `charts` instead of `tool_to_file`, or skip tools that aren't in `charts`. The third loop (1041-1047) is already guarded."

**Repro script:** `reviews/2026-10-16-repro/gui/test_classification_refresh.py` (run with `python`, not pytest; it stubs NiceGUI via `gui/stubs2`). It shows **only GUI-1**, so when this task is done, delete it and update the `reviews/2026-10-16-repro/README.md` table row.

**Files:**
- Modify: `src/robin/gui/components/classification.py` — add a module-level helper; rewrite the body of `_refresh_classification` (def at `:989`) between `:997` and `:1027` to call it, and change the third loop at `:1041` to iterate `charts`.
- Test: `tests/test_classification_refresh.py` (create).

**Interfaces:**
- Consumes: nothing from other tasks.
- Produces: `classification_files_to_refresh(charts: Dict[str, Dict[str, Any]], sample_dir: Optional[Path]) -> Tuple[Dict[str, float], Dict[str, bool], bool]`, module-level in `src/robin/gui/components/classification.py`, returning `(mtimes, changed, is_fresh_visit)` keyed by the tool names **present in `charts`** and nothing else. `mtimes[tool]` is the file's `st_mtime`, or `0.0` when the file is missing or `sample_dir` is falsy. `changed[tool]` is True when the file exists and its mtime is greater than `charts[tool]["last_mtime"]`, treating a `last_mtime` of `None` as "changed". `is_fresh_visit` is True when any tool in `charts` still has `last_mtime is None`. Each entry of `charts` already carries its own `"file"` (set at `:587-595`), so the helper does not take `tool_to_file`.

- [ ] **Step 1: Write the failing tests** in `tests/test_classification_refresh.py`. Build `charts` as a plain dict of the shape the production code uses — only the keys the helper reads, i.e. `{"file": <name>, "last_mtime": <float or None>}` — and write real files under `tmp_path`.

```python
def test_only_charts_tools_are_considered(tmp_path):
    """A tool absent from charts must not be looked up at all (GUI-1: KeyError)."""
    (tmp_path / "sturgeon_scores.csv").write_text("x\n")
    charts = {"Sturgeon": {"file": "sturgeon_scores.csv", "last_mtime": None}}
    mtimes, changed, fresh = classification_files_to_refresh(charts, tmp_path)
    assert set(mtimes) == {"Sturgeon"}
    assert set(changed) == {"Sturgeon"}
    assert changed["Sturgeon"] is True
    assert fresh is True

def test_unchanged_file_is_not_flagged(tmp_path):
    p = tmp_path / "sturgeon_scores.csv"; p.write_text("x\n")
    charts = {"Sturgeon": {"file": "sturgeon_scores.csv", "last_mtime": p.stat().st_mtime}}
    mtimes, changed, fresh = classification_files_to_refresh(charts, tmp_path)
    assert changed["Sturgeon"] is False
    assert fresh is False

def test_missing_file_reports_zero_and_no_change(tmp_path):
    charts = {"NanoDX": {"file": "NanoDX_scores.csv", "last_mtime": 1.0}}
    mtimes, changed, fresh = classification_files_to_refresh(charts, tmp_path)
    assert mtimes["NanoDX"] == 0.0
    assert changed["NanoDX"] is False

def test_no_enabled_classifiers_is_not_a_fresh_visit(tmp_path):
    """Review Focus 5: charts empty. Must not crash and must not claim a fresh visit."""
    assert classification_files_to_refresh({}, tmp_path) == ({}, {}, False)

def test_missing_sample_dir_is_handled(tmp_path):
    charts = {"Sturgeon": {"file": "sturgeon_scores.csv", "last_mtime": None}}
    mtimes, changed, fresh = classification_files_to_refresh(charts, None)
    assert mtimes == {"Sturgeon": 0.0}
    assert changed == {"Sturgeon": False}
    assert fresh is True
```

- [ ] **Step 1b: Write the end-to-end test that actually reddens for GUI-1's reason.** The five tests above can only fail by `ImportError` — by construction they exercise a helper that does not exist yet, so **none of them can distinguish a half fix**. That matters here: the `KeyError` does **not** come from the first loop. Measured by running the repro script, it comes from the `is_fresh_visit` comprehension at `classification.py:1019-1020`:

```
File ".../classification.py", line 1020, in <genexpr>
    charts[tool_name].get("last_mtime") is None
KeyError: 'PanNanoDX'
```

The first loop at `:1003-1016` only touches `charts[tool_name]` inside `if file_path.exists()`, and a disabled classifier has written no score file — so it never raises. An implementer who guards only the first loop gets all five helper tests green, a green suite, and a still-frozen panel.

So add one test that drives `add_classification_section` with a disabled classifier and asserts the panel advances. **Run it in a subprocess** (`subprocess.run([sys.executable, ...], capture_output=True, text=True)` on a small script this test writes into `tmp_path`): the flow needs NiceGUI stubbed, and stubbing NiceGUI in-process would poison the real-NiceGUI imports that Tasks 4 and 6 rely on in the same pytest session. Model the script on `reviews/2026-10-16-repro/gui/test_classification_refresh.py`, which already does exactly this and whose `stubs2/` directory stays in the repo.

```python
def test_panel_keeps_refreshing_with_a_classifier_disabled(tmp_path):
    """GUI-1 end to end: workflow_steps omits pannanodx/random_forest, so charts
    lacks them; before the fix the second tick raised KeyError and the label froze."""
    script = tmp_path / "drive.py"
    script.write_text(_DRIVER_SCRIPT)          # module-level string in this test file
    r = subprocess.run([sys.executable, str(script), str(tmp_path)],
                       capture_output=True, text=True, env=_stubbed_env())
    assert r.returncode == 0, r.stderr
    assert "KeyError" not in r.stderr
    ticks = [l.split(": ", 1)[1] for l in r.stdout.splitlines() if l.startswith("tick")]
    assert ticks[0] != ticks[1], f"panel froze: {ticks}"
    assert "Glioblastoma" in ticks[1]
```

The driver writes `sturgeon_scores.csv` with one row, calls `add_classification_section`, prints the Sturgeon summary label after tick 1, appends a second row whose Glioblastoma score is 0.97, and prints the label after ticks 2 and 3. `_stubbed_env()` puts `reviews/2026-10-16-repro/gui/stubs2` and `src` on `PYTHONPATH`.

- [ ] **Step 2: Run and confirm RED.** Run `/home/parlar_ai/robin-batch3-run/t3.sh tests/test_classification_refresh.py -v`. Expected: the five helper tests fail by collection error / `ImportError` (the helper does not exist yet) — **say explicitly in your report that this is the only reason they redden**. The Step 1b test must fail for GUI-1's own reason: `KeyError` in stderr and both ticks showing the same frozen label. Capture both outputs verbatim. Also run the repro script directly (`PYTHONDONTWRITEBYTECODE=1 /home/parlar_ai/robin-batch3-run/venv/bin/python reviews/2026-10-16-repro/gui/test_classification_refresh.py` from the worktree root) and record its three printed lines — the label is frozen at "Medulloblastoma" across all three ticks.

- [ ] **Step 3: Implement `classification_files_to_refresh`** at module level in `src/robin/gui/components/classification.py`, above `add_classification_section`. Then rewrite `_refresh_classification`'s body: replace the `file_mtimes` / `files_changed` / `any_changes` loop (`:997-1016`) **and the `is_fresh_visit` comprehension at `:1019-1022` — that comprehension is where the `KeyError` is actually raised, so a fix that leaves it alone fixes nothing** — with one call to the helper plus `any_changes = any(files_changed.values())`. Then change the update loop at `:1041` to iterate `charts.items()` instead of `tool_to_file.items()` (each `charts` entry carries `"file"`, so `cfg["file"]` still reads correctly). The logging block between them indexes `tool_to_file[tool_name]` at `:1036`, so it must index `charts` too, or read `cfg["file"]`. `tool_to_file` keeps its role as the construction-time source at `:338-343`.

After editing, grep the function for every remaining `charts[` subscript and confirm each one is reached only for a key the helper returned: `grep -n 'charts\[' src/robin/gui/components/classification.py`. Report that list.

- [ ] **Step 4: Run and confirm GREEN.** `/home/parlar_ai/robin-batch3-run/t3.sh tests/test_classification_refresh.py -v` → 6 passed. Then the full suite → **86 passed** (80 + 6). Then `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src /home/parlar_ai/robin-batch3-run/venv/bin/python -c "import robin.gui.components.classification"` → exit 0. Then re-run the repro script and show the label now advancing to Glioblastoma. State in the report which lines of `_refresh_classification` are covered **only** by the import check.

- [ ] **Step 5: Retire the repro script — only if Step 1b's test exists and reddens.** The rule is that a pytest regression test *replaces* the script; the five helper tests do not, because they cannot fail for GUI-1's reason. If Step 1b landed and was seen RED, `git -C <worktree> rm reviews/2026-10-16-repro/gui/test_classification_refresh.py` and update its row in `reviews/2026-10-16-repro/README.md` ("Which script shows which finding" table). If Step 1b could not be made to work, **keep the script** and say so — a deleted script with no replacement is a net loss of coverage. Either way `reviews/2026-10-16-repro/gui/stubs2/` stays: Step 1b's driver depends on it. Report whether any other script also uses `stubs2/`.

- [ ] **Step 6: Commit.** Message names GUI-1 and says the panel now keeps refreshing with a classifier disabled.

---

### Task 2: REP-1 — Random Forest confidence uses the default cut-offs

> **Finding (`^\*\*REP-1 `):** high, Confirmed (spot-checked). "Random Forest confidence uses the default cut-offs. **Where:** `reporting/sections/classification.py:399`. **What's wrong:** It looks up `"random forest"`, but the key is `"random_forest"`, so the default 80/50 is used instead of 85/65. At 82%, the PDF says 'High' and the GUI says 'Medium'. **Fix:** map display names to config keys, and raise an error on an unknown classifier."

**Repro script:** `reviews/2026-10-16-repro/reporting/rf_conf.py`. It shows **only REP-1** → delete it in this task and update the README table row.

**Files:**
- Modify: `src/robin/classification_config.py` — add the display-name map and the resolver next to `CLASSIFIER_CONFIDENCE_THRESHOLDS` (`:11-32`).
- Modify: `src/robin/reporting/sections/classification.py:399` — resolve the display name before calling `get_confidence_status`; import the resolver at `:25`.
- Test: `tests/test_classifier_config_keys.py` (create).

**Interfaces:**
- Produces: `config_key_for_classifier(display_name: str) -> str` in `src/robin/classification_config.py`. Maps a human-facing classifier name to its `CLASSIFIER_CONFIDENCE_THRESHOLDS` key, case-insensitively and tolerating surrounding whitespace, and raises `ValueError` naming the unknown name and listing the known ones for anything else. Backed by a module-level `CLASSIFIER_DISPLAY_NAMES: Dict[str, str]` mapping `"Sturgeon"→"sturgeon"`, `"NanoDX"→"nanodx"`, `"PanNanoDX"→"pannanodx"`, `"Random Forest"→"random_forest"` — the same four display names `src/robin/gui/config.py:37-42` (`CLASSIFICATION_STEPS`) already uses, and the same four `reporting/sections/classification.py:340-343` reads score files for.
- **Do NOT make `get_confidence_level`, `get_confidence_ui_tier` or `get_confidence_status` raise on an unknown key.** Measured call sites that pass deliberate non-classifier keys and must keep falling back to `DEFAULT_CONFIDENCE_THRESHOLDS`: `gui/components/summary.py:423` passes `anchor_key`, which is `"coverage"`, `"cnv"`, `"mgmt"` or `"fusion"` (`:274`, `:285`, `:296`, `:310`); `gui/components/summary.py:1768` defaults `classifier="default"`. The raising happens in the new resolver only, at the one call site that has a display name.

- [ ] **Step 1: Write the failing tests** in `tests/test_classifier_config_keys.py`.

```python
def test_random_forest_display_name_maps_to_the_config_key():
    assert config_key_for_classifier("Random Forest") == "random_forest"

def test_every_display_name_resolves_to_a_real_threshold_entry():
    for display in CLASSIFIER_DISPLAY_NAMES:
        assert config_key_for_classifier(display) in CLASSIFIER_CONFIDENCE_THRESHOLDS

def test_resolution_is_case_and_whitespace_insensitive():
    assert config_key_for_classifier("  random forest ") == "random_forest"
    assert config_key_for_classifier("STURGEON") == "sturgeon"

def test_unknown_classifier_raises_and_names_itself():
    with pytest.raises(ValueError, match="Dinosaur"):
        config_key_for_classifier("Dinosaur")

def test_random_forest_at_82_percent_is_medium_not_high():
    """REP-1: the whole point. 82% is above the 80 default but below the RF 85."""
    status, _ = get_confidence_status(config_key_for_classifier("Random Forest"), 0.82)
    assert status == "Medium"

def test_the_report_resolves_the_display_name_before_looking_up_thresholds():
    """Covers reporting/sections/classification.py:399 itself — without this, nothing
    proves the call site was changed, only that the resolver exists."""
    src = (REPO_ROOT / "src/robin/reporting/sections/classification.py").read_text()
    assert "get_confidence_status(name.lower()" not in src
    assert "config_key_for_classifier" in src

def test_non_classifier_anchors_still_fall_back_to_the_defaults():
    """Review Focus 1: these must not raise — the summary cards depend on them."""
    for anchor in ("coverage", "cnv", "mgmt", "fusion", "default"):
        assert get_confidence_ui_tier(anchor, 82.0) == "high"   # 82 >= default high of 80
        assert get_confidence_level(anchor, 82.0) == "High confidence"
```

Derive `REPO_ROOT` as `Path(__file__).resolve().parents[1]`, not from the cwd.

Note: do **not** also assert that the old wrong key `"random forest"` still scores 82% as "High".
That would enshrine the fallback for a key nobody should pass, and it passes before and after the
fix, so it proves nothing.

- [ ] **Step 2: Run and confirm RED.** `/home/parlar_ai/robin-batch3-run/t3.sh tests/test_classifier_config_keys.py -v` — most tests fail by `ImportError` for the new names, which is the weak kind of red. The two that redden for REP-1's own reason are `test_the_report_resolves_the_display_name_before_looking_up_thresholds` (the call site still reads `get_confidence_status(name.lower()`) and `test_non_classifier_anchors_still_fall_back_to_the_defaults`, which is a **guard** that passes before and after — say so. Also run the repro script `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src /home/parlar_ai/robin-batch3-run/venv/bin/python reviews/2026-10-16-repro/reporting/rf_conf.py` and record what it prints; that is the independent RED for REP-1.

- [ ] **Step 3: Implement** `CLASSIFIER_DISPLAY_NAMES` and `config_key_for_classifier` in `src/robin/classification_config.py`, then change `reporting/sections/classification.py:399` from `get_confidence_status(name.lower(), confidence_value)` to resolve `name` through the new function first. Note that `:399` sits inside a `try:` whose `except Exception` at `:411` logs and `continue`s — so a raise there skips that classifier's row rather than aborting the report. Say in the report whether you consider that acceptable for an unknown classifier, and do not restructure the `try`.

- [ ] **Step 4: Run and confirm GREEN.** Per-file → 7 passed. Full suite → **93 passed** (86 + 7). Re-run the repro script and show the cut-offs now correct.

- [ ] **Step 5: Retire the repro script** (`git rm reviews/2026-10-16-repro/reporting/rf_conf.py`) and update its README table row.

- [ ] **Step 6: Commit.** Message names REP-1.

---

### Task 3: JOB-6 — the 60-second idle timeout is a test value

> **Finding (`^\*\*JOB-6 `):** high, Constant confirmed (spot-checked); impact Plausible. "The 60-second idle timeout is a test value. **Where:** `gui_launcher.py:420, 927-934, 7406, 7424-7426`, `target_analysis.py:2604-2618`. **What's wrong:** The code still says `completion_timeout_seconds = 60  # Set to 3600 for production`. A gap of 60 s or more between BAMs marks the sample Complete and finalises it. Later SNP runs then skip the merge and use the old `target.bam`. After a restart, unmerged batch BAMs are deleted."

**No repro script.**

**Files:**
- Modify: `src/robin/gui_launcher.py` — add a module-level helper near the top-level functions, and change the assignment at `:418-420`.
- Test: `tests/test_completion_timeout.py` (create).

**Interfaces:**
- Produces: `completion_timeout_seconds() -> int`, module-level in `src/robin/gui_launcher.py`. Returns the idle gap in seconds after which a sample with no active or pending jobs is marked Complete and finalised. Reads `ROBIN_COMPLETION_TIMEOUT_SECONDS`; returns `3600` when the variable is unset, empty, whitespace, not an integer, or less than `1`. Mirror the established idiom at `src/robin/workflow_ray.py:956-962` (`_job_timeout_seconds`) — same `os.environ.get(...).strip()` / `int(...)` / `except (ValueError, TypeError)` shape — so the codebase keeps one pattern for env-var timeouts.
- **One deliberate difference from that idiom, and the docstring must say it:** `_job_timeout_seconds` treats `0` as *disabled*; here `0` (and any value `< 1`) means *fall back to 3600*. This timeout triggers finalisation, and finalisation deletes unmerged batch BAMs, so a `0` that was honoured would destroy data on the first tick. Do not silently inherit the other function's meaning for `0`.
- Every existing reader already goes through `self.completion_timeout_seconds`, so the single assignment at `:420` is the whole wiring change. **Measured: 15 occurrences in the file — one write at `:420` and fourteen reads** (`:905, 927, 936, 941, 3447, 3465, 3487, 3825, 7869, 7872, 7891, 8176, 8181, 8412`). Verify with `grep -n completion_timeout_seconds src/robin/gui_launcher.py | wc -l` before and after and report both counts; the count must stay 15 and only `:420` may change.

- [ ] **Step 1: Write the failing tests** in `tests/test_completion_timeout.py`. Use pytest's `monkeypatch.setenv` / `delenv`; the helper must read the environment on each call so a test can set it.

```python
def test_default_is_one_hour(monkeypatch):
    monkeypatch.delenv("ROBIN_COMPLETION_TIMEOUT_SECONDS", raising=False)
    assert completion_timeout_seconds() == 3600

def test_env_var_overrides_the_default(monkeypatch):
    monkeypatch.setenv("ROBIN_COMPLETION_TIMEOUT_SECONDS", "120")
    assert completion_timeout_seconds() == 120

@pytest.mark.parametrize("raw", ["", "   ", "abc", "12.5", "-1", "0"])
def test_unusable_values_fall_back_to_the_default(monkeypatch, raw):
    """Review Focus 2: a value that became 0 would finalise every sample on the first tick."""
    monkeypatch.setenv("ROBIN_COMPLETION_TIMEOUT_SECONDS", raw)
    assert completion_timeout_seconds() == 3600

def test_surrounding_whitespace_is_tolerated(monkeypatch):
    monkeypatch.setenv("ROBIN_COMPLETION_TIMEOUT_SECONDS", " 900 ")
    assert completion_timeout_seconds() == 900
```

- [ ] **Step 2: Run and confirm RED.** Expected: `ImportError` — the helper does not exist. Record it.

- [ ] **Step 3: Implement** `completion_timeout_seconds()` and change `gui_launcher.py:420` to `self.completion_timeout_seconds = completion_timeout_seconds()`. Replace the two stale comment lines at `:418-419` ("Change to 60 for testing (1 minute), 3600 for production" / "Set to 3600 for production") with one line naming the env var and what the timeout actually triggers — a sample being marked Complete and its `target.bam` finalised. Do not change any of the ten read sites, and do not add a CLI option: the env var is the whole configuration surface for this task.

- [ ] **Step 4: Run and confirm GREEN.** Per-file → 9 passed (4 named tests, one of them parametrised 6 ways). Full suite → **102 passed** (93 + 9). Also run `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src /home/parlar_ai/robin-batch3-run/venv/bin/python -c "import robin.gui_launcher"` → exit 0, since `gui_launcher.py` is ~8000 lines of GUI code with no harness.

- [ ] **Step 5: Commit.** Message names JOB-6 and states both halves: the default is now an hour, and `ROBIN_COMPLETION_TIMEOUT_SECONDS` overrides it.

---

### Task 4: GUI-4 — the CNV summary card shows the batch counter as "Gained"

> **Finding (`^\*\*GUI-4 `):** high, Confirmed (spot-checked). "The CNV summary card shows the batch counter as 'Gained'. **Where:** `summary.py:1245-1292`. **What's wrong:** It reads `cnv_analysis_counter.txt` as the number of gains. After 37 batches the card says 'Gained: 37, Lost: 0' for every sample. The other files it looks for are never written."

**Repro script:** `reviews/2026-10-16-repro/gui/test_summary_cnv.py` shows GUI-4 **and GUI-7 and GUI-9, which stay open** → **keep the script** and only note in the README table that GUI-4 is fixed in our fork.

**Measured background:** `cnv_analysis_counter.txt` really is written, as a batch counter — by `save_analysis_counter` at `analysis/cnv_analysis.py:804-811`; the references at `:789` and `analysis/fusion_work.py:3084` *read* it. Of the other three names in the list, `cnv_results.csv` and `cnv_summary.txt` are written by nothing in `src/` (grep both names across `src/`: the only hits are this candidate list and `gui_launcher.py:64-65`, which also only reads). `cnv_analysis_results.pkl` is written by `analysis/cnv_subprocess.py:143` (`analysis/cnv_analysis.py:492` reads it), but the loop has no branch for `.pkl` — it matches, reads nothing, and `break`s. So after this fix the counts are always absent, which is the honest answer.

**Files:**
- Modify: `src/robin/gui/components/summary.py` — `_extract_cnv_data` (`:1204-1297`): drop `"cnv_analysis_counter.txt"` from the candidate list at `:1250-1255`, remove the now-dead `.txt` branch, and default `"gained"`/`"lost"` to `None` instead of `0` (`:1210-1211`). Then `_create_cnv_dashboard_card_with_data` (`:641-689`) must render an absent count as `"Not available"`, matching the idiom already used for `genetic_sex`, `bin_width` and `variance`, **and must not call `int()` on `None` at `:650`**. Call site at `:283-284` passes `cnv_data.get("gained", 0)` — change both defaults to `None`.
- Test: `tests/test_summary_cnv_counts.py` (create).

**Interfaces:**
- Consumes: nothing.
- Produces: `_extract_cnv_data(sample_dir: Path) -> Dict[str, Any]` keeps its name and signature; `"gained"` and `"lost"` become `Optional[int]`, `None` when no file supplied a real count. `_create_cnv_dashboard_card_with_data(genetic_sex: str, bin_width: str, variance: str, gained: Optional[int], lost: Optional[int], anchor_key: str = "cnv") -> None`. Its `total = int(gained) + int(lost)` at `:650` must treat `None` as contributing nothing, and `total_denom` stays `"—"` when no count is known, which is already what it prints for `total == 0`.
- **Explicit non-goal:** do not start reading `cnv_analysis_results.pkl` to compute real gain/loss counts. That is a feature, not this fix.

- [ ] **Step 1: Write the failing tests** in `tests/test_summary_cnv_counts.py`. `robin.gui.components.summary` imports cleanly under the test env (measured), and both functions are module-level, so call them directly. Only `_extract_cnv_data` needs asserting; for the card, assert the formatting helper's behaviour through `_extract_cnv_data`'s output shape rather than rendering NiceGUI.

```python
def test_batch_counter_is_not_reported_as_gains(tmp_path):
    """GUI-4: after 37 batches the card said 'Gained: 37'."""
    (tmp_path / "cnv_analysis_counter.txt").write_text("37\n")
    data = _extract_cnv_data(tmp_path)
    assert data["gained"] is None
    assert data["lost"] is None

def test_no_cnv_files_leaves_the_counts_unknown(tmp_path):
    data = _extract_cnv_data(tmp_path)
    assert data["gained"] is None and data["lost"] is None
    assert data["genetic_sex"] == "Not available"

def test_a_real_results_csv_still_supplies_the_counts(tmp_path):
    (tmp_path / "cnv_results.csv").write_text("gained,lost,genetic_sex\n3,2,XY\n")
    data = _extract_cnv_data(tmp_path)
    assert data["gained"] == 3
    assert data["lost"] == 2
    assert data["genetic_sex"] == "XY"
```

**And one test for the card itself — this one is not optional.** `_create_cnv_dashboard_card_with_data` does `total = int(gained) + int(lost)` at `summary.py:650`. Passing `None` raises `TypeError` **today**, and the import check never executes that line, so a forgotten guard would break every summary page while the suite stayed green. Drive the card with a stub `ui` object (monkeypatch `robin.gui.components.summary.ui` with a small fake whose `element`/`column`/`row`/`label`/`icon` return context managers and record the label texts):

```python
def test_the_card_renders_unknown_counts_without_crashing(monkeypatch):
    """GUI-4: `int(None)` at summary.py:650 would take down the whole summary page."""
    texts = _install_stub_ui(monkeypatch)        # helper in this test file
    _create_cnv_dashboard_card_with_data("XY", "1,000", "0.1", None, None)
    joined = " | ".join(texts)
    assert "Gained: 0" not in joined
    assert "Lost: 0" not in joined
    assert "Not available" in joined
```

- [ ] **Step 2: Run and confirm RED.** The first test must fail asserting `37 is None` — that is GUI-4 itself. The card test must fail with `TypeError: int() argument must be ... not 'NoneType'` from `summary.py:650` — that is the crash the fix must not ship. Record which of the four redden and which are guards; the third (`cnv_results.csv`) pins the surviving `.csv` branch and may already pass.

- [ ] **Step 3: Implement** the three changes named under Files. Keep the `except Exception: pass` style of the surrounding block rather than tightening it.

- [ ] **Step 4: Run and confirm GREEN.** Per-file → 4 passed. Full suite → **106 passed** (102 + 4). Then `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src /home/parlar_ai/robin-batch3-run/venv/bin/python -c "import robin.gui.components.summary"` → exit 0, and state which card-rendering lines are covered **only** by that import check (the stub-`ui` test now covers `:650` and the label calls, so the list should be shorter than you expect).

- [ ] **Step 5: Update the README table row** for `gui/test_summary_cnv.py` to say GUI-4 is fixed in our fork while GUI-7 and GUI-9 remain. **Keep the script.**

- [ ] **Step 6: Commit.** Message names GUI-4 and says the card now shows the counts as unavailable rather than inventing them.

---

### Task 5: CNV-5 — every chromosome gets a fake breakpoint at its end

> **Finding (`^\*\*CNV-5 `):** low, Confirmed (by reading). "Every chromosome gets a fake breakpoint at its end. **Where:** `cnv_analysis.py:574-580, 714-730`. **What's wrong:** ruptures always returns the series length as its last change point. That adds about 24 spurious regions to the master BED. **Fix:** use `result[:-1]`."

**No repro script.**

**Files:**
- Modify: `src/robin/analysis/cnv_analysis.py` — `run_ruptures` (`:574-580`), the list comprehension that turns change points into `(start, end)` pairs.
- Test: `tests/test_cnv_changepoints.py` (create).

**Interfaces:**
- `run_ruptures(r_cnv, penalty_value, bin_width)` keeps its signature and its `List[Tuple[int, int]]` return. The only change is that the final element of `algo_c.predict(...)`, which `ruptures` always sets to the series length, no longer produces a range. The consumer at `:714-730` keeps its `if start >= 0 and end > start` guard — do not remove it.
- `ruptures` 1.1.10 is installed in the test env, so `robin.analysis.cnv_analysis` imports and `run_ruptures` can be called directly on a synthetic signal.

- [ ] **Step 1: Write the failing tests** in `tests/test_cnv_changepoints.py`. Build a signal with one obvious mean shift so `ruptures` returns exactly `[shift, len(signal)]`, then assert on what `run_ruptures` emits.

```python
def test_the_series_end_does_not_become_a_region():
    """CNV-5: ruptures always returns len(signal) last; it is not a breakpoint."""
    import numpy as np
    signal = np.concatenate([np.zeros(100), np.ones(100) * 5.0])
    ranges = run_ruptures(signal, penalty_value=10, bin_width=1000)
    assert all(end <= len(signal) * 1000 - 1000 for _, end in ranges)
    assert len(ranges) == 1

def test_the_detected_shift_is_still_reported():
    import numpy as np
    signal = np.concatenate([np.zeros(100), np.ones(100) * 5.0])
    (start, end), = run_ruptures(signal, penalty_value=10, bin_width=1000)
    assert start == 100 * 1000 - 1000
    assert end == 100 * 1000 + 1000

def test_a_flat_signal_yields_no_regions():
    """Review Focus 5: ruptures returns only [len(signal)], so result[:-1] is empty."""
    import numpy as np
    assert run_ruptures(np.zeros(200), penalty_value=10, bin_width=1000) == []
```

- [ ] **Step 2: Run and confirm RED.** Expected on the old code: `len(ranges) == 2` in the first test and a non-empty list in the third, because the series-end change point still produces a range.

**These numbers were measured, not predicted** — on `ruptures` 1.1.10 in this env, `KernelCPD(kernel="linear").predict(pen=10)` returns `[100, 200]` for the step signal and `[200]` for the flat one, so `run_ruptures(..., bin_width=1000)` currently yields `[(99000, 101000), (199000, 201000)]` and `[(199000, 201000)]`, and with the `[:-1]` slice yields `[(99000, 101000)]` and `[]`. Re-run the measurement yourself before trusting it; if your numbers differ, use yours and say so in the report.

- [ ] **Step 3: Implement** the `[:-1]` slice in `run_ruptures`, with a comment saying why the last element is not a breakpoint.

- [ ] **Step 4: Run and confirm GREEN.** Per-file → 3 passed. Full suite → **109 passed** (106 + 3).

- [ ] **Step 5: Commit.** Message names CNV-5 and says roughly 24 spurious master-BED regions per run go away.

---

### Task 6: GUI-2 — a brief websocket drop permanently stops every refresh timer

> **Finding (`^\*\*GUI-2 `):** critical, Confirmed (against NiceGUI 3.6 source). "A brief websocket drop permanently stops every refresh timer. **Where:** `summary.py:352`, `classification.py:1068`, `coverage.py:6917`, `mgmt.py:721`, `cnv.py:2467`, `fusion.py:3421`, `mnpflex.py:860`, `bed_coverage.py:456`. **What's wrong:** Every component calls `timer.deactivate()` in `on_disconnect`, which in NiceGUI 3 also fires on reconnect, and nothing reactivates the timers. **Effect:** after a Wi-Fi blip, laptop sleep or proxy timeout, the page shows hours-old results. **Fix:** use `client.on_delete`."

**No repro script.** This is a batched same-shape change: nine one-line call-site edits and nothing else.

**Measured:** `nicegui.Client.on_delete(handler)` exists in the installed NiceGUI 3.18.0 with the same single-callable signature as `on_disconnect`. **There are nine sites, not the eight the finding lists** — `src/robin/gui_launcher.py:5130` is a ninth `on_disconnect(lambda: refresh_timer.deactivate())`, and `coverage.py` has moved from the finding's `:6917` to `:6938` because of earlier work on this branch. Confirm all nine for yourself with `grep -rn on_disconnect src/robin/gui/ src/robin/gui_launcher.py` and report the count you found.

All nine callbacks were read: every one is a timer teardown, so there is no legitimate non-timer `on_disconnect` in these files and the lint test below can assert zero occurrences. Four of them also unregister a theme-sync callback (`classification.py`, `coverage.py`, `mgmt.py`, `cnv.py`), `fusion.py` unregisters one across several lines, and `bed_coverage.py:456-460` deactivates **two** timers (`result_timer` and `refresh_timer`) plus a theme-sync. Preserve every callback body exactly; only the method name changes.

**Files:** modify only the `on_disconnect` → `on_delete` call at each of:
`src/robin/gui/components/summary.py:352`, `classification.py:1068`, `coverage.py:6938`, `mgmt.py:721`, `cnv.py:2467`, `fusion.py:3421`, `mnpflex.py:860`, `bed_coverage.py:456`, and `src/robin/gui_launcher.py:5130`.
- Test: `tests/test_gui_timer_cleanup.py` (create).

**Interfaces:** nothing is produced or consumed. Each callback body stays exactly as it is — `classification.py:1068` also calls `unregister_classification_theme_sync()`, which must be preserved.

- [ ] **Step 1: Write the failing test** in `tests/test_gui_timer_cleanup.py`. A NiceGUI page cannot be driven here, so pin the thing that is actually checkable: that no GUI source registers a timer teardown on `on_disconnect`. Read the files as text — this is a lint-shaped regression test, and it is honest about being one.

```python
GUI_FILES = [  # the nine sites GUI-2 names, plus gui_launcher
    "src/robin/gui/components/summary.py", "src/robin/gui/components/classification.py",
    "src/robin/gui/components/coverage.py", "src/robin/gui/components/mgmt.py",
    "src/robin/gui/components/cnv.py", "src/robin/gui/components/fusion.py",
    "src/robin/gui/components/mnpflex.py", "src/robin/gui/components/bed_coverage.py",
    "src/robin/gui_launcher.py",
]

def test_no_refresh_timer_is_torn_down_on_disconnect():
    """GUI-2: on_disconnect also fires on reconnect in NiceGUI 3, and nothing reactivates."""
    offenders = []
    for rel in GUI_FILES:
        text = (REPO_ROOT / rel).read_text()
        for n, line in enumerate(text.splitlines(), 1):
            if "on_disconnect" in line:
                offenders.append(f"{rel}:{n}")
    assert offenders == []

def test_every_file_that_had_one_now_uses_on_delete():
    for rel in GUI_FILES:
        assert "on_delete" in (REPO_ROOT / rel).read_text(), rel

def test_on_delete_exists_in_the_installed_nicegui():
    """If NiceGUI ever drops on_delete, this fix is silently a no-op."""
    from nicegui import Client
    assert callable(getattr(Client, "on_delete", None))
```

Derive `REPO_ROOT` from `__file__` (`Path(__file__).resolve().parents[1]`), not from the cwd.

- [ ] **Step 2: Run and confirm RED.** Expected: the first test fails listing all nine `file:line` offenders, and the second fails on the files that have no `on_delete` yet. Record the full offender list — it is the measurement of how many sites there are.

- [ ] **Step 3: Implement** the nine edits. One hunk per file, changing only the method name. **Do not use a repo-wide `sed`**: `on_disconnect` may appear in a comment or an unrelated handler, and the Global Constraints forbid unanchored scripted edits. Edit each site individually and check `git diff --stat` shows exactly nine files with one changed line each (plus whatever the comment above a site needs).

- [ ] **Step 4: Run and confirm GREEN.** Per-file → 3 passed. Full suite → **112 passed** (109 + 3).

Then check each of the nine modules still loads. **`robin.gui.components.mnpflex` cannot be imported in this environment** — measured: it raises `ModuleNotFoundError: No module named 'alive_progress'` via `matkit.py:27`, and installing is forbidden. So import the other eight and byte-compile mnpflex:

```
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src /home/parlar_ai/robin-batch3-run/venv/bin/python -c "import robin.gui.components.summary, robin.gui.components.classification, robin.gui.components.coverage, robin.gui.components.mgmt, robin.gui.components.cnv, robin.gui.components.fusion, robin.gui.components.bed_coverage, robin.gui_launcher"
PYTHONDONTWRITEBYTECODE=1 /home/parlar_ai/robin-batch3-run/venv/bin/python -m py_compile src/robin/gui/components/mnpflex.py && echo "mnpflex compiles"
```

Both must succeed. Say in the report that mnpflex got a syntax check only, and why.

State plainly that **no test proves a timer survives a reconnect** — that needs a live browser — and that the suite's coverage here is the lint test plus the import/compile check.

- [ ] **Step 5: Commit.** Message names GUI-2, says nine call sites, and notes that `gui_launcher.py:5130` was a ninth site the review did not list.

---

### Task 7: VAR-3 — a ClinVar refresh keeps the old tabix index

> **Finding (`^\*\*VAR-3 `):** high, Confirmed (reproduced with bcftools). "A ClinVar refresh keeps the old tabix index. **Where:** `utils/clinvar_manager.py:130-131, 221, 300`. **What's wrong:** After a new `clinvar.vcf.gz` is swapped in, `_ensure_tabix_index` returns early because the old `.tbi` exists. Region queries then return nothing with exit code 0, so pathogenic annotations silently disappear. **Fix:** always rebuild the index after replacing the file."

**Repro script:** `reviews/2026-10-16-repro/appshell/stale_tbi.sh` (needs `bcftools`, which is on PATH). It shows **only VAR-3** → delete it in this task and update the README table row.

**Files:**
- Modify: `src/robin/utils/clinvar_manager.py` — `_ensure_tabix_index` (`:126-131`) gains the force parameter; `ensure_clinvar_files` (`:220-221`) and `update_clinvar_if_newer` (`:298-300`) pass it where they have just written the `.gz`.
- Test: `tests/test_clinvar_index.py` (create).

**Interfaces:**
- `_ensure_tabix_index(gz_path: Path, tbi_path: Path, *, force: bool = False) -> None`. With `force=True` the `if _file_ok(tbi_path): return` short-circuit at `:130-131` is skipped and the index is rebuilt; the existing `pysam.tabix_index(..., force=True, keep_original=True)` call and its `tabix`-binary fallback are unchanged. Default `False` keeps today's behaviour for any caller that has not just replaced the file.
- Call sites to change: `update_clinvar_if_newer` at `:300`, immediately after `tmp_path.replace(gz_path)` — always `force=True`. And `ensure_clinvar_files` at `:221` — `force=True` **only when it just ran `_download_url_to_file` at `:220`**, because a fresh `.gz` can arrive beside a stale `.tbi`; when the `.gz` was already present and fine, leave `force=False` so an existing good index is not rebuilt on every start-up.
- There is a **third** place a `.gz` is written, at `update_clinvar_if_newer:273-277`, which calls `_ensure_tabix_index` not at all. It is practically unreachable on the current flow, so do not restructure it — but name it in your report so the next reader knows it was seen and left.
- `pysam` 0.24.1 is installed; there is no standalone `tabix` binary, so the `tabix`-binary fallback at `:146-168` cannot be exercised here. Say so rather than claiming it works.

- [ ] **Step 1: Write the failing tests** in `tests/test_clinvar_index.py`. Build two tiny single-record VCFs with different positions, bgzip and index them with `pysam`, then swap the `.gz` under a surviving `.tbi` and assert what a region query returns. Keep both VCFs synthetic and hand-written.

**CRITICAL — the obvious version of this test does not work.** A tabix index addresses records by
16 kb bin. Swapping a record at position 100 for one at position 200 leaves both in the *same* bin,
so the **stale** index still finds the new record. Measured in this env:

```
same-bin (100->200), stale index, fetch 150-250:                 ['200']
far-bin  (100->5e6), stale index, fetch 4999000-5001000:         []
far-bin  (100->5e6), REBUILT index, same fetch:                  ['5000000']
```

A same-bin test would pass before *and* after the fix and prove nothing. **Put the new record in a
different bin** — position 100 replaced by position 5,000,000, queried at 4,999,000-5,001,000 — and
additionally assert the `.tbi` file actually changed.

```python
def test_a_replaced_vcf_gets_a_rebuilt_index(tmp_path):
    """VAR-3: the old .tbi survives the swap, so the new records are invisible."""
    gz, tbi = _write_indexed_vcf(tmp_path, pos=100)            # helper in this file
    _overwrite_gz_only(gz, pos=5_000_000)                      # .tbi now stale
    _ensure_tabix_index(gz, tbi, force=True)
    with pysam.TabixFile(str(gz)) as tf:
        assert [r.split("\t")[1] for r in tf.fetch("chr1", 4_999_000, 5_001_000)] == ["5000000"]

def test_force_actually_rewrites_the_index_file(tmp_path):
    """Guards against a `force` parameter that is accepted and then ignored."""
    gz, tbi = _write_indexed_vcf(tmp_path, pos=100)
    _overwrite_gz_only(gz, pos=5_000_000)
    before = tbi.stat().st_mtime_ns
    _ensure_tabix_index(gz, tbi, force=True)
    assert tbi.stat().st_mtime_ns != before

def test_without_force_a_good_index_is_not_rebuilt(tmp_path):
    gz, tbi = _write_indexed_vcf(tmp_path, pos=100)
    before = tbi.stat().st_mtime_ns
    _ensure_tabix_index(gz, tbi)
    assert tbi.stat().st_mtime_ns == before

def test_the_stale_index_really_hides_the_new_record(tmp_path):
    """Pins the bug itself: this is what users saw — zero results, exit code 0."""
    gz, tbi = _write_indexed_vcf(tmp_path, pos=100)
    _overwrite_gz_only(gz, pos=5_000_000)
    with pysam.TabixFile(str(gz)) as tf:
        assert list(tf.fetch("chr1", 4_999_000, 5_001_000)) == []
```

Write `_write_indexed_vcf` with `pysam.tabix_compress` then `pysam.tabix_index(..., preset="vcf",
force=True, keep_original=True)`, and `_overwrite_gz_only` by compressing a second hand-written VCF
to a temp name and copying it over the `.gz` so the `.tbi` is left behind. Give the `##contig` header
a length large enough for position 5,000,000 (`length=250000000`).

- [ ] **Step 2: Run and confirm RED.** Tests 1 and 2 must fail with a `TypeError` for the unexpected `force` keyword. Then run the repro script (`bash reviews/2026-10-16-repro/appshell/stale_tbi.sh`) and record its output — that is the independent RED for VAR-3. Tests 3 and 4 are **guards**: test 3 pins that a good index is left alone, test 4 pins the old broken behaviour, and both pass before and after. Say so. Test 4 is what makes test 1's assertion meaningful, so keep it.

- [ ] **Step 3: Implement** the `force` parameter and the two call-site changes. The docstring must say which caller is expected to pass `force=True` and why a `.tbi` that passes `_file_ok` can still be wrong.

- [ ] **Step 4: Run and confirm GREEN.** Per-file → 4 passed. Full suite → **116 passed** (112 + 4). Re-run the repro script and show the query now returning the new record.

- [ ] **Step 5: Retire the repro script** (`git rm reviews/2026-10-16-repro/appshell/stale_tbi.sh`) and update its README table row. The README also mentions this script by name in its "How to run" section (`appshell/stale_tbi.sh also needs bcftools`) — remove that line too.

- [ ] **Step 6: Commit.** Message names VAR-3.

---

### Task 8: VAR-7 — `robin add-panel` keeps only the first region of each gene

> **Finding (`^\*\*VAR-7 `):** high, Confirmed. "`robin add-panel` keeps only the first region of each gene. **Where:** `cli.py:957`. **What's wrong:** `drop_duplicates(subset=['gene'])` drops every later exon or tile of a gene. In the test, a panel with 3 EGFR rows and 2 TP53 rows became 2 rows. Coverage, SNP calling and fusion detection then ignore the dropped regions."

**Repro script:** `reviews/2026-10-16-repro/appshell/repro_unique_bed.py <tmpdir>`. It shows **only VAR-7** → delete it in this task and update the README table row.

**Files:**
- Modify: `src/robin/cli.py:958` inside `_generate_unique_gene_bed` (def at `:916`).
- Test: `tests/test_add_panel_regions.py` (create).

**Interfaces:**
- `_generate_unique_gene_bed(input_bed_path: Path, output_bed_path: Path) -> bool` keeps its name, signature and return. It still writes a headerless 4-column BED of `chrom start end gene`, still expands a comma-separated `gene` cell into one row per gene (`:943-952`), and still sorts by `['chrom', 'start', 'end']` (`:961`). A prototype of this fix was run against the four test inputs below and produced 5 / 1 / 2 / 2 lines with no pandas warnings. The only change is the dedupe key: **drop duplicates on the whole region** (`subset=['chrom', 'start', 'end', 'gene']`) instead of on `['gene']` alone, so every distinct region survives while an exactly repeated row still collapses.
- `robin.cli` imports cleanly under the test env and the function is module-level, so call it directly; do not invoke the click command.

- [ ] **Step 1: Write the failing tests** in `tests/test_add_panel_regions.py`.

```python
def test_every_exon_of_a_gene_survives(tmp_path):
    """VAR-7: 3 EGFR + 2 TP53 rows became 2 rows."""
    src = tmp_path / "panel.bed"
    src.write_text(
        "chr7\t100\t200\tEGFR\n" "chr7\t300\t400\tEGFR\n" "chr7\t500\t600\tEGFR\n"
        "chr17\t100\t200\tTP53\n" "chr17\t300\t400\tTP53\n"
    )
    out = tmp_path / "out.bed"
    assert _generate_unique_gene_bed(src, out) is True
    lines = [l for l in out.read_text().splitlines() if l.strip()]
    assert len(lines) == 5
    assert sum(1 for l in lines if l.endswith("EGFR")) == 3

def test_exactly_repeated_rows_still_collapse(tmp_path):
    """Review Focus 4: keep every distinct region, but do not reintroduce duplicates."""
    src = tmp_path / "panel.bed"
    src.write_text("chr7\t100\t200\tEGFR\n" "chr7\t100\t200\tEGFR\n")
    out = tmp_path / "out.bed"
    _generate_unique_gene_bed(src, out)
    assert len([l for l in out.read_text().splitlines() if l.strip()]) == 1

def test_the_same_gene_on_two_chromosomes_keeps_both(tmp_path):
    """Review Focus 4: a pseudogene or a mis-annotated panel must not lose a locus."""
    src = tmp_path / "panel.bed"
    src.write_text("chr1\t100\t200\tFOO\n" "chr2\t100\t200\tFOO\n")
    out = tmp_path / "out.bed"
    _generate_unique_gene_bed(src, out)
    assert len([l for l in out.read_text().splitlines() if l.strip()]) == 2

def test_a_comma_separated_gene_cell_still_expands(tmp_path):
    src = tmp_path / "panel.bed"
    src.write_text("chr7\t100\t200\tEGFR,ERBB1\n")
    out = tmp_path / "out.bed"
    _generate_unique_gene_bed(src, out)
    genes = sorted(l.split("\t")[3] for l in out.read_text().splitlines() if l.strip())
    assert genes == ["EGFR", "ERBB1"]
```

- [ ] **Step 2: Run and confirm RED.** Expected: test 1 fails with 2 lines instead of 5, and test 3 fails with 1 line instead of 2. Tests 2 and 4 are guards on behaviour the fix must preserve and will pass before the fix — say so.

- [ ] **Step 3: Implement** the `subset` change at `:958`, with a comment saying a gene legitimately has many regions and only an identical region is a duplicate. Also fix the now-wrong comment at `:957` ("Remove duplicates based on gene name (keep first occurrence)").

- [ ] **Step 4: Run and confirm GREEN.** Per-file → 4 passed. Full suite → **120 passed** (116 + 4). Then run the repro script with a temp dir and show it now keeping all five regions.

- [ ] **Step 5: Retire the repro script** (`git rm reviews/2026-10-16-repro/appshell/repro_unique_bed.py`) and update its README table row.

- [ ] **Step 6: Commit.** Message names VAR-7.

---

### Task 9: Record the work on the branch

Batch 1 and batch 2 each recorded their own fixes in their own commit, and the review document is
the project's running record of what is fixed in this fork versus upstream. Leaving that out would
mean eight fixes that nobody reading the review can see. This task is bookkeeping only — **no
production code changes**, so no new tests and the suite must stay at **120 passed**.

**Files:**
- Modify: `CHANGELOG.md` — `[Unreleased]`.
- Modify: `reviews/2026-10-16-code-review.md` — only the places listed below.
- Modify: `reviews/2026-10-16-repro/README.md` — the table rows the earlier tasks changed, if any were missed.
- Regenerate: `reviews/2026-10-16-code-review.pdf`.

**What may change in the review document, and nothing else.** The findings keep their wording.
Permitted edits, exactly:
1. **Status lines** — add a `- **Status:**` line to each of GUI-1, REP-1, JOB-6, GUI-4, CNV-5, GUI-2, VAR-3, VAR-7 saying it is fixed in our fork, with the commit, and what remains. Follow the wording pattern of the existing Status line on VAR-8 (`grep -n '^- \*\*Status:\*\*' reviews/2026-10-16-code-review.md` to see them all). **JOB-6's Status line must say that only the constant and its configurability were addressed** — the finding's second and third bullets (a later SNP run skipping the merge, and unmerged batch BAMs being deleted after a restart) are consequences of the timeout firing early and are not fixed by changing it; they belong with JOB-7.
2. **Rows of the "This fork and the original ROBIN" table** — one row per finding, with the commit and what a user notices.
3. **The `*(fixed in our fork)*` marks in the Fix-first table.**
4. **The `*(done in our fork)*` mark on step 3 of "Suggested order of work".**
5. **The "Already done" paragraph.**
6. **The `**Updated:**` header line.**

- [ ] **Step 1: Add the `CHANGELOG.md` `[Unreleased]` entries** — one per finding, in the file's existing voice: what a user notices, not what the code now does.

- [ ] **Step 2: Edit the review document**, the six permitted kinds of change and nothing else. Then prove you changed nothing else: `git -C <worktree> diff -- reviews/2026-10-16-code-review.md` and read every hunk. Any hunk touching finding prose is a mistake to undo. Paste the `--stat` and the list of hunk headers into your report.

- [ ] **Step 3: Check the repro README** is consistent with what the earlier tasks did: `gui/test_classification_refresh.py`, `reporting/rf_conf.py`, `appshell/stale_tbi.sh` and `appshell/repro_unique_bed.py` should be gone from the table (unless Task 1 kept its script, in which case say so), `gui/test_summary_cnv.py` should remain with GUI-4 noted as fixed, and the "How to run" mention of `appshell/stale_tbi.sh also needs bcftools` should be gone. `ls reviews/2026-10-16-repro/*/` and compare against the table.

- [ ] **Step 4: Rebuild the PDF.** `bash /home/parlar_ai/robin-batch3-run/pdf/build.sh /home/parlar_ai/dev/robin/.claude/worktrees/batch-3`. It prints the page count — record it.

- [ ] **Step 5: Confirm nothing regressed.** Full suite → **120 passed, 1 warning**. `git -C <worktree> diff eccbf18..HEAD -- robin.yml pyproject.toml` → empty.

- [ ] **Step 6: Commit.** One commit, message saying the batch 3 fixes are recorded in the changelog and the review.

---

## Expected suite arithmetic

The chain starts at **80 passed** (71 at `eccbf18`, plus 9 tests added by the four batch 2b
follow-up commits already on this branch). Each task states its own count; the predicted totals are:

| After task | Tests added | Suite |
|---|---|---|
| 1 — GUI-1 | 6 (5 helper + 1 end-to-end subprocess) | 86 |
| 2 — REP-1 | 7 | 93 |
| 3 — JOB-6 | 9 (4 tests, one parametrised 6×) | 102 |
| 4 — GUI-4 | 4 (3 extract + 1 card render) | 106 |
| 5 — CNV-5 | 3 | 109 |
| 6 — GUI-2 | 3 | 112 |
| 7 — VAR-3 | 4 | 116 |
| 8 — VAR-7 | 4 | 120 |
| 9 — record the work | 0 | 120 |

These are predictions, not targets. If a count differs, find the term that changed and name it —
never adjust an expected number to make the arithmetic close.

## How this plan was reviewed before execution

A fresh reviewer on the most capable model checked every claim here against the spec, the repro
scripts and the real code before any task started, and found four defects that would have produced
a false pass. All four are fixed above, and each fix was independently re-measured:

1. **Task 7's tabix test was broken.** Positions 100 and 200 share a 16 kb tabix bin, so the *stale*
   index still returned the new record — the test would have passed before and after. Now uses
   position 5,000,000 and also asserts the `.tbi` mtime changed.
2. **Task 4 had an untested crash.** `int(None)` at `summary.py:650` raises today and the import
   check never reaches it, so the fix could have shipped a broken summary page with a green suite.
   Now has a stub-`ui` test that reddens for exactly that.
3. **Task 6's import check could not pass.** `robin.gui.components.mnpflex` cannot be imported in
   this env (`alive_progress` missing, installing forbidden). Now imports eight and byte-compiles
   the ninth.
4. **Task 1's tests could not fail for GUI-1's reason.** The `KeyError` comes from the
   `is_fresh_visit` comprehension at `:1020`, not the first loop, so a half fix would have gone
   green. Now has an end-to-end subprocess test, and the repro script is only retired if that test
   exists and was seen red.

The reviewer also confirmed, by running them: the 80-passed baseline, all five repro script outputs,
that the independent suite still fails ir01/ir02/ir08/ir09, that all nine `on_disconnect` callbacks
are timer teardowns, that NiceGUI 3.18's `Client.on_delete` exists and that `on_disconnect`'s own
docstring says it fires on reconnect, and that `ruptures` returns `[100, 200]` and `[200]` for the
two signals Task 5 uses.
