# Batch 2b: Lightweight Variant Check Follow-ups Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Outcome (2026-10-31):** done on `fix/batch-3-quick-wins` (`87984a2`, `80fda58`, `abfe903`, `a3eecaf`), not yet merged. All three tasks were built as planned; the suite went 71 → 80.

What was built differently, or found to be different:

- **Task 2's open question was already answered.** The plan asked whether `_check_lga_results` is in scope for `check_lga_status`, and said to find "the smallest way to make the refresh happen (e.g. pass a callback)" if not. It *is* in scope: both are locals of `add_coverage_section`, verified with Python's `ast`, so the change is a direct call and no callback, registry or parameter was added.
- **Task 2 needed a fix round for a defect the plan created.** Adding the refresh call made `_check_lga_results` reachable *after* a render for the first time, which exposed that its two non-display branches never cleared `lga_results_container`. Without that, a run whose results fail the trust check would set the label "Old variant results can't be confirmed as this sample's" while leaving a live "View Analysis Results" button pointing at that very file. Both branches now clear first (`abfe903`).
- **Task 3 needed the shared BAM fixture extended.** `tests/lga_fixtures.py::write_bam` could only emit all-match CIGARs, so it could not build the all-deletion pileup column the plan asks for. `write_bam` now accepts an optional third tuple element, a CIGAR string; existing two-tuple callers produce a byte-identical BAM (verified by md5).
- **Task 3's all-deletion test is a guard, not a regression test** — it passed before the fix, because the old code's `("T", 0)` satisfied its assertions too. It is kept because its `reference in ("A","T","G","C","N")` assertion is the only thing in the suite that would catch a wrong fix dropping the `!= "indel"` filter.
- **The plan's line number for `_perform_targeted_pileup` was stale** (~928; the function starts at 822). Noted, not treated as a defect.

Three deferred minors went to the final whole-branch review; none blocks merge.

**Goal:** Fix the three minors the batch 2 final review deferred: a wrong folder when the monitored root has the sample's name, the GUI showing a stale "run again" message (or an untrusted file over a trusted one), and the pileup's `reference` always being `"T"`.

**Architecture:** All in `src/robin/analysis/lightweight_gene_analysis.py` (LGA) plus the GUI's `_check_lga_results` in `src/robin/gui/components/coverage.py`. The GUI choice of results file moves into a tested LGA helper.

**Tech Stack:** Python 3.12, pysam, pytest, NiceGUI (wiring only).

**Spec:** the batch 2 final review's Minor findings 3, 4 and 5, quoted under each task. Background: `reviews/2026-10-16-code-review.md` VAR-8/VAR-9 and `reviews/plans/2026-10-30-batch-2-sample-identity.md`.

## Global Constraints

- Tests: `/home/parlar_ai/robin-batch3-run/t3.sh <pytest args>` (default `tests`). Baseline at `eccbf18`: **71 passed**.
- Synthetic data only. Reuse `tests/lga_fixtures.py` (`write_bam`, `write_clinvar`, `write_bed`).
- Keep every JSON key and every behaviour batch 2 pinned; all 71 existing tests stay green.
- One commit per task; the message names the follow-up it fixes.

## Review Focus

1. **Layout (b) — `work_dir=<root>/A`, sample `A`, no `<root>/A/A/igv`:** must still resolve to `<root>/A`. Pinned by the existing `test_sample_folder_passed_directly_still_works`.
2. **The handler creates `<root>/A/A/` as its output folder** in layout (b). A resolver that tests `isdir(join(base, id))` would flip after one run. Pinned in Task 1 (`test_output_folder_does_not_change_the_resolved_folder`).
3. **Two results files, old one untrusted, new one trusted:** the trusted one must be shown. Pinned in Task 2.
4. **A pileup column with zero counts everywhere** (all reads deleted at that position): choosing the reference must not crash. Pinned in Task 3.
5. **Insertion ClinVar site on reference reads:** must label `insertion_requires_review`, not `reference_mismatch`. Pinned in Task 3.

---

### Task 1: Root folder named like the sample

> Review: "`resolve_sample_folder` regresses the case where the monitored root is named like the sample (`/x/A` holding sample `/x/A/A`): the basename rule returns the root, and the handler fails with 'No igv_ready.bam … at /x/A/igv/igv_ready.bam'. Fix: resolve to `join(base, id)` when `join(base, id, "igv")` exists, else fall back to `base` on basename match. Don't test `isdir(join(base, id))` alone: the handler itself creates that folder as the output dir."

**Files:** LGA `resolve_sample_folder`; test `tests/test_lga_sample_inputs.py`.

**Interfaces:** `resolve_sample_folder(base_dir, sample_id) -> str` keeps its signature and its `ValueError` rules. New order: no ID/`"unknown"` → `base_dir`; `join(base, id, "igv")` is a directory → `join(base, id)`; basename matches → `base_dir`; else `join(base, id)`.

- [ ] **Step 1: Failing tests.**
```python
def test_root_named_like_the_sample_resolves_to_the_sample_folder(tmp_path):
    root = tmp_path / "A"                      # monitored root happens to be called "A"
    (root / "A" / "igv").mkdir(parents=True)   # sample A lives at /x/A/A
    assert resolve_sample_folder(str(root), "A") == str(root / "A")

def test_output_folder_does_not_change_the_resolved_folder(tmp_path):
    sample = tmp_path / "root" / "A"
    (sample / "igv").mkdir(parents=True)
    (sample / "A").mkdir()                     # the handler's output folder in layout (b)
    assert resolve_sample_folder(str(sample), "A") == str(sample)
```
- [ ] **Step 2: Run.** Expected: the first fails (returns `/x/A`); the second passes and guards the change.
- [ ] **Step 3: Implement** the new order. Docstring says why `igv` is the marker.
- [ ] **Step 4: All tests pass.** **Step 5: Commit.**

### Task 2: The GUI shows the first trusted results file, and refreshes after a run

> Review: "`_check_lga_results` prefers the first existing file over the first trusted one. An old untrusted file at `<sample>/<id>/…json` permanently hides a fresh trusted one at `<sample>/…json` … pick the first trusted candidate and show the old-results message only when none is trusted. Related (pre-existing): nothing re-runs `_check_lga_results()` after a run completes (only call at coverage.py ~5896), so a user who follows the new message's advice sees it unchanged until they reload; call it from `check_lga_status` on completion."

**Files:** LGA (new helper); `coverage.py` `_check_lga_results` (~4804) and `check_lga_status` (~6109); test `tests/test_lga_sample_inputs.py`.

**Interfaces:** `pick_lga_results(candidates: list[str], sample_folder: str) -> tuple[Optional[str], bool]` in LGA: returns `(first candidate that exists and lga_results_belong_to_sample says True, found_untrusted)` where `found_untrusted` is True when some candidate exists but none is trusted. `(None, False)` when none exist.

- [ ] **Step 1: Failing tests** (use the existing `_results_file` helper in that test file):
  - trusted file at the second candidate, untrusted at the first → returns the second, `False`.
  - only an untrusted file → `(None, True)`.
  - nothing on disk → `(None, False)`.
  - both trusted → the first.
- [ ] **Step 2: Run.** Expected: `ImportError`.
- [ ] **Step 3: Implement** the helper. In `_check_lga_results`, use it: trusted → display (keep the "(legacy location)" wording when it is the fallback path); untrusted only → the existing "Old variant results can't be confirmed…" message; nothing → "No lightweight variant analysis results found". In `check_lga_status`, when the status becomes `completed`, call `_check_lga_results()` if it is in scope; if it is not, find the smallest way to make the refresh happen (e.g. pass a callback) and record the choice in the ledger.
- [ ] **Step 4: All tests pass, plus** `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src /home/parlar_ai/robin-batch3-run/venv/bin/python -c "import robin.gui.components.coverage"` from the worktree (GUI wiring has no test harness; say so in the commit).
- [ ] **Step 5: Commit.**

### Task 3: The pileup's `reference` is the most common base

> Review: "Pileup `reference` is always `"T"` (`max((k, v) …)` compares base letters, so `("T", 0)` beats `("C", 20)`). Pre-existing. It makes every insertion site report `reference_mismatch` instead of `insertion_requires_review` … and `has_variants`/`variant_alleles` are meaningless. Fix: `max(…, key=lambda kv: kv[1])`."

**Files:** LGA `_perform_targeted_pileup` (~line 928); test `tests/test_lga_variant_support.py`.

- [ ] **Step 1: Failing tests.**
  - `_perform_targeted_pileup` on 20 reads of `"C"*50` (chr1 0-50): the entry at position 11 has `reference == "C"`, `variant_alleles == []`, `has_variants is False`.
  - ClinVar insertion `(11, "C", "CA")` on 20 reference reads → the variant's `evidence_level == "insertion_requires_review"` and the gene status stays `insufficient_evidence`.
  - A column where every read has a deletion: no crash (build reads with a `D` in the CIGAR at that position, or assert via a direct call with a crafted `base_counts` if the pileup can't produce it — record which).
- [ ] **Step 2: Run.** Expected: reference `"T"`; insertion labelled `reference_mismatch`.
- [ ] **Step 3: Implement** `max(..., key=lambda kv: kv[1])` over the non-indel bases. Comment why.
- [ ] **Step 4: All tests pass.** **Step 5: Commit.**
