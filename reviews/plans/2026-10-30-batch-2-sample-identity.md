# Batch 2: Sample Identity and Truthful Variant Summaries Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The lightweight variant check (LGA) uses only the sample's own files, counts a pathogenic variant only when reads support it, and analyses every BED interval of a gene.

**Architecture:** All three bugs live in `src/robin/analysis/lightweight_gene_analysis.py` (LGA below). VAR-9 adds one function that resolves the sample's folder and replaces both `os.walk` fallbacks with a lookup inside that folder; the results file records which BAM and BED were used, and the GUI refuses result files that don't prove they came from this sample. VAR-10 runs a pileup per interval inside `_analyze_single_region`. VAR-8 adds one rule for "supported" and bases every count and status on it.

**Tech Stack:** Python 3.12, pysam, pandas, pytest, NiceGUI (GUI wiring only).

**Spec:** [`reviews/2026-10-16-code-review.md`](../2026-10-16-code-review.md) VAR-8, VAR-9, VAR-10 (IR-01, IR-02, IR-09). The independent suite `reviews/2026-10-16-repro/independent/reproduce.py` reproduces all three.

**Outcome (2026-10-30):** done and merged to `main` (`b9ffef8`..`765f004`). The user chose "3 or more alternate reads". Changes from this plan:
- Task 3 merges overlapping intervals before the pileup, because concatenating pileups double-counted the overlap.
- Task 4 makes the GUI variant tables list sites that are found or need review.

The final review added:
- Same-length multi-base changes (`CC>TT`) and reads showing neither allele now give "insufficient evidence".
- A per-variant `needs_review` flag.
- The independent suite now loads the new helpers.

## Global Constraints

- Branch `fix/batch-2-sample-identity` in a git worktree cut from `main` (`01728bf`).
- Test env: Python 3.12 venv with numpy 1.26.4, pandas 2.1.4, pysam, pytest, pytest-cov. Run `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=src python -m pytest <path> -p no:cacheprovider --no-cov -q` from the repo root.
- Synthetic data only: tiny BAMs written with pysam and indexed with `pysam.index`; a hand-written ClinVar VCF passed as metadata `clinvar_path`. No patient data.
- Keep every JSON key the GUI reads today. New keys are additions.
- Each task is its own commit, message names its finding ID.
- **Decision — what counts as "supported":** 3 or more alternate reads, which is where LGA's own evidence levels stop calling it `low_support` (evidence levels `low_frequency`, `medium_frequency`, `high_frequency`). Sites with 1–2 alternate reads, and indel sites (always `*_requires_review` or `reference_mismatch`), are "insufficient evidence". Sites with 0 alternate reads are "not detected".

## Review Focus

1. **A sample ID with path parts** (`../B`, `A/B`, absolute path): must never reach another sample's folder. Pinned in Task 1 (`test_sample_id_cannot_leave_the_base_folder`).
2. **Result files from before this fix:** they may hold another sample's evidence and must not be shown as this sample's. Pinned in Task 2 (`test_old_results_without_inputs_are_not_trusted`).
3. **Overlapping intervals with the same name:** a ClinVar site in both must be counted once. Pinned in Task 3 (`test_site_in_overlapping_intervals_counts_once`).
4. **Indel ClinVar sites:** LGA can't call them, so they must read "insufficient evidence", never "present" or "none detected". Pinned in Task 4 (`test_indel_site_is_insufficient_evidence`).
5. **A sample folder passed with a trailing slash** (`/root/A/`): must resolve to the same folder. Pinned in Task 1 (`test_trailing_slash_resolves_to_the_same_folder`).

## Shared test fixture

`tests/lga_fixtures.py` (a helper module, not a test file):
- `write_bam(path, reads)`: `reads` is a list of `(start, sequence)` on contig `chr1` (length 1000), all base qualities 30, MAPQ 60; sorts and indexes with pysam.
- `write_clinvar(path, sites)`: `sites` is a list of `(pos, ref, alt)`; each line `1\t<pos>\t<id>\t<ref>\t<alt>\t.\t.\tCLNSIG=Pathogenic;GENEINFO=TEST:1;CLNDN=test`, with a `##fileformat=VCFv4.1` header.
- `write_bed(path, rows)`: `rows` are `(start, end, name)` on `chr1`.

Keep ClinVar positions off the BED start (LGA's intersection is `start <= POS <= end`, its pileup covers `start+1..end`).

---

### Task 1: The LGA handler uses only the sample's own folder (VAR-9)

**Files:**
- Modify: LGA `lightweight_gene_analysis_handler` (1503-1768): BAM search 1562-1624, BED search 1650-1680, JSON export call ~1730
- Modify: LGA `export_json(gene_reports, output_path, inputs=None)` (1317)
- Create: `tests/lga_fixtures.py`, `tests/test_lga_sample_inputs.py`

**Interfaces:**
- Produces: `resolve_sample_folder(base_dir: str, sample_id: Optional[str]) -> str` in LGA. No ID or `"unknown"` → `base_dir`. `basename(normpath(base_dir)) == sample_id` → `base_dir`. Otherwise `join(base_dir, sample_id)`. Raises `ValueError` if the ID is absolute or has a path separator or `..`.
- Produces: JSON `metadata.input_bam`, `metadata.input_bed`, `metadata.sample_id` (absolute paths), via `export_json(..., inputs={...})`; the same two paths in `add_result`.

- [ ] **Step 1: Write the failing tests** in `tests/test_lga_sample_inputs.py`. Layout under `tmp_path`: `root/A/targets_exceeding_threshold.bed`, `root/A/target.bam`, `root/B/igv/igv_ready.bam` (20 alt reads), `root/B/targets_exceeding_threshold.bed`. The job is a `SimpleNamespace` with `context.filepath`, `context.metadata={"clinvar_path": ...}`, `context.get_sample_id=lambda: "A"`, and recording `add_result`/`add_error`.

```python
def test_missing_bam_does_not_use_another_samples_bam(...):
    with pytest.raises(FileNotFoundError, match=r"A/igv/igv_ready\.bam"):
        lightweight_gene_analysis_handler(job, work_dir=str(root))
    assert not list(root.rglob("lightweight_gene_analysis_results.json"))

def test_missing_bed_does_not_use_another_samples_bed(...):  # A has igv bam, no BED, no bed_path
    with pytest.raises(FileNotFoundError, match=r"A/targets_exceeding_threshold\.bed"): ...

def test_results_record_the_inputs_used(...):  # A has both
    meta = json.loads((root / "A" / "lightweight_gene_analysis_results.json").read_text())["metadata"]
    assert meta["input_bam"] == str(root / "A" / "igv" / "igv_ready.bam")
    assert meta["input_bed"] == str(root / "A" / "targets_exceeding_threshold.bed")

def test_sample_folder_passed_directly_still_works(...):  # work_dir=root/A, ID "A"
    # uses root/A/igv/igv_ready.bam; output stays where the GUI reads it first: root/A/A/...json

def test_unknown_sample_id_looks_only_in_the_given_folder(...):  # get_sample_id -> "unknown", work_dir=root
    with pytest.raises(FileNotFoundError): ...  # B's BAM is never used

@pytest.mark.parametrize("sample_id", ["../B", "A/../B", "/abs/B"])
def test_sample_id_cannot_leave_the_base_folder(tmp_path, sample_id):
    with pytest.raises(ValueError):
        resolve_sample_folder(str(tmp_path / "root"), sample_id)

def test_trailing_slash_resolves_to_the_same_folder(tmp_path):
    assert resolve_sample_folder(str(tmp_path / "A") + "/", "A") == str(tmp_path / "A") + "/"  # basename matches
```

- [ ] **Step 2: Run them.** Expected: the two "missing" tests fail (B's files are used and a JSON is written under A); `test_results_record_the_inputs_used` fails with `KeyError: 'input_bam'`; the `resolve_sample_folder` tests fail with `ImportError`.

- [ ] **Step 3: Implement.**
  - Add `resolve_sample_folder`.
  - In the handler, resolve the folder once. BAM is `<folder>/igv/igv_ready.bam` only. BED is metadata `bed_path` if it exists, else `<folder>/targets_exceeding_threshold.bed`. Missing either: raise `FileNotFoundError` naming the expected path (the existing `except` already calls `add_error` and re-raises). Delete both `os.walk` loops and the dead second BAM check.
  - Leave the output-folder logic as it is (the GUI reads both locations).
  - Pass `inputs={"input_bam": ..., "input_bed": ..., "sample_id": ...}` to `export_json`, which merges it into `metadata`. Add the two paths to `add_result`.

- [ ] **Step 4: Run all tests.** Expected: pass.

- [ ] **Step 5: Commit** ("VAR-9: LGA uses only the sample's own BAM and BED, and records them").

### Task 2: The GUI doesn't show results that may belong to another sample (VAR-9)

**Why:** result files written before Task 1 may hold another sample's evidence. They have no `input_bam`, so they can't prove where they came from.

**Files:**
- Modify: LGA (new function), `src/robin/gui/components/coverage.py` `_check_lga_results` (4804-4839)
- Test: `tests/test_lga_sample_inputs.py`

**Interfaces:**
- Consumes: `metadata.input_bam` from Task 1.
- Produces: `lga_results_belong_to_sample(json_path: str, sample_folder: str) -> bool` in LGA: True only when `metadata.input_bam` exists and its real path is inside the real path of `sample_folder`. Unreadable JSON → False.

- [ ] **Step 1: Write the failing tests.**

```python
def test_old_results_without_inputs_are_not_trusted(tmp_path):  # metadata has no input_bam
    assert lga_results_belong_to_sample(str(json_path), str(sample_a)) is False

def test_results_from_another_samples_bam_are_not_trusted(tmp_path):  # input_bam under root/B
    assert ... is False

def test_results_from_this_samples_bam_are_trusted(tmp_path):  # input_bam under root/A
    assert ... is True
```

- [ ] **Step 2: Run them.** Expected: `ImportError`.
- [ ] **Step 3: Implement** the function. In `_check_lga_results`, call it on whichever file it found; when False, set the status to "Old variant results can't be confirmed as this sample's. Run the variant analysis again." and don't display them.
- [ ] **Step 4: Run all tests.** Expected: pass. Also run `python -c "import robin.gui.components.coverage"` in the test env to catch syntax errors in the GUI file.
- [ ] **Step 5: Commit** ("VAR-9: don't show LGA results that can't be tied to this sample").

### Task 3: Analyse every interval of a gene (VAR-10)

**Files:**
- Modify: LGA `_analyze_single_region` (676-735), `_determine_analysis_quality` (1068)
- Test: `tests/test_lga_intervals.py`

**Interfaces:**
- Keeps: `analyze_bed_regions(bam, bed) -> Dict[str, GeneReport]` keyed by BED name; `_find_bed_clinvar_intersections` still assigns each site to the first interval it falls in.

- [ ] **Step 1: Write the failing tests.** BED `chr1 0 30 TEST` and `chr1 100 130 TEST`; ClinVar C>T at 11 and 111; reads `"T"*30` ×10 at 0 and ×30 at 100.

```python
def test_every_interval_with_the_same_name_is_analysed(...):
    r = LightweightGeneAnalysis(str(tmp_path / "w"), str(vcf)).analyze_bed_regions(str(bam), str(bed))["TEST"]
    assert sorted(v["position"] for v in r.pathogenic_snps_details) == [11, 111]

def test_coverage_covers_every_interval(...):
    assert r.mean_coverage == pytest.approx(20.0)   # 30 positions at 10x, 30 at 30x

def test_site_in_overlapping_intervals_counts_once(...):  # BED 0-30 TEST and 5-30 TEST, site at 11
    assert [v["position"] for v in r.pathogenic_snps_details] == [11]
```

- [ ] **Step 2: Run them.** Expected: the first two fail (only position 11; mean 10.0). The overlap test may already pass; it guards the change.
- [ ] **Step 3: Implement.** Collect the distinct `intersecting_region` intervals of the group, in order. Pileup each and concatenate. `_find_intersecting_variants` and `_calculate_coverage_statistics` work on the combined pileup. Report `chromosome`, `start=min`, `end=max` (same as today for one interval). Give `_determine_analysis_quality` the summed interval length instead of `end - start`.
- [ ] **Step 4: Run all tests.** Expected: pass.
- [ ] **Step 5: Commit** ("VAR-10: analyse every BED interval of a gene").

### Task 4: Count a pathogenic variant only when reads support it (VAR-8)

**Files:**
- Modify: LGA `_analyze_single_region` (`pathogenic_snps_found`), `export_json` (variant blocks 1366-1394 and 1419-1447, summary 1448-1467, metadata 1474-1487)
- Modify: `coverage.py` summary cards (4926-4960)
- Test: `tests/test_lga_variant_support.py`

**Interfaces:**
- Produces: `is_supported(evidence: dict) -> bool` and `needs_review(evidence: dict) -> bool` in LGA, per the Global Constraints decision.
- JSON per variant: new `supported: bool`. Per gene `summary`: `total_variants`, `high_coverage_variants`, `low_coverage_variants` count **supported** sites only; new `sites_checked`, `insufficient_evidence_sites`; `pathogenic_status` is `present_with_good_coverage` (a supported site at ≥10×), `present_but_low_coverage`, `insufficient_evidence`, or `none_detected`. `metadata`: `total_variants_found`, `genes_with_pathogenic_variants`, `genes_with_good_coverage_variants` count supported only; new `total_sites_checked`, `genes_with_good_coverage` (genes with any checked site at ≥10×).
- The per-gene `variants` list keeps every checked site, with its `evidence_analysis`.

- [ ] **Step 1: Write the failing tests.** BED `chr1 0 50 TEST`; ClinVar C>T at 11. Helper runs `analyze_bed_regions` then `export_json` and returns `(report, json)`.

```python
def test_reference_only_site_is_not_reported_as_present(...):  # 20 reads of "C"*50
    assert report.pathogenic_snps_found == 0
    gene = data["genes"]["TEST"]
    assert gene["summary"]["pathogenic_status"] == "none_detected"
    assert gene["summary"]["total_variants"] == 0 and gene["summary"]["sites_checked"] == 1
    assert gene["variants"][0]["supported"] is False
    assert gene["variants"][0]["evidence_analysis"]["evidence_level"] == "no_support"
    assert data["metadata"]["genes_with_pathogenic_variants"] == 0
    assert data["metadata"]["genes_with_good_coverage"] == 1

def test_supported_site_is_reported(...):  # 20 reads with T at 11
    assert gene["summary"]["pathogenic_status"] == "present_with_good_coverage"
    assert data["metadata"]["total_variants_found"] == 1

def test_supported_site_at_low_coverage(...):  # 5 alt reads
    assert gene["summary"]["pathogenic_status"] == "present_but_low_coverage"

def test_two_alternate_reads_are_insufficient_evidence(...):  # 18 ref + 2 alt
    assert gene["summary"]["pathogenic_status"] == "insufficient_evidence"
    assert gene["summary"]["total_variants"] == 0

def test_indel_site_is_insufficient_evidence(...):  # ClinVar CA>C at 11, 20 ref reads
    assert gene["summary"]["pathogenic_status"] == "insufficient_evidence"
```

- [ ] **Step 2: Run them.** Expected: the reference-only, 2-read and indel tests fail (status `present_with_good_coverage`, counts 1); the supported tests may pass already and guard the change.
- [ ] **Step 3: Implement** the two rules and the counts above, once, used by both variant blocks in `export_json` and by `_analyze_single_region`. In `coverage.py`, the "Genes with Good Coverage (≥10x)" card and the "Low Coverage" card read `genes_with_good_coverage` instead of `genes_with_good_coverage_variants`. The gene table needs no change: it title-cases the status, so `insufficient_evidence` shows as "Insufficient Evidence".
- [ ] **Step 4: Run all tests**, plus the GUI import check. Expected: pass.
- [ ] **Step 5: Commit** ("VAR-8: count pathogenic variants only when reads support them").

### Task 5: Record the fixes

**Files:** `CHANGELOG.md`, `reviews/2026-10-16-code-review.md`, `reviews/2026-10-16-repro/independent/README.md`

- [ ] **Step 1: Run** `reviews/2026-10-16-repro/independent/reproduce.py`. Expected: `test_ir01`, `test_ir02` and `test_ir09` now **fail** (the bugs are gone), together with `test_ir08` from batch 1; the other 7 still reproduce.
- [ ] **Step 2: CHANGELOG `[Unreleased]` → Fixed:** VAR-9, VAR-10, VAR-8. Under **Changed:** the new statuses and JSON keys; old LGA result files are not shown until the analysis is run again.
- [ ] **Step 3: Review:** a **Status** line on VAR-8, VAR-9 and VAR-10; rows in "This fork and the original ROBIN"; the Fix-first table rows; step 2 of the order of work marked done. Rebuild the PDF.
- [ ] **Step 4: Independent README:** add the three tests to the "Fixed since" note.
- [ ] **Step 5: Commit** ("Record batch 2 fixes in the changelog and review").

---

**Not in this batch:** the dead `analyze_genes`/`_analyze_single_gene` path (no callers) and the unused HTML report (`generate_report`, no callers) keep the old counting; both belong to ARCH-4 (delete dead code). Building `igv_ready.bam` from `target.bam` when it's missing is a possible follow-up; this batch only makes the missing file an explicit error.

**After the batch:** a fresh reviewer reviews the whole branch, then merge to `main` and push to `gmc-norr`.
