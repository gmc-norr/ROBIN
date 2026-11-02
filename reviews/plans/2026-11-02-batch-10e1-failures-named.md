# Batch 10e part 1: the target step and the variant check name what they could not do

**Locked before code, 2026-11-02.** Base `fe10f4b`. Branch `fix/batch-10e-recorded-items`. Done by the controller (no agents). Small fixes that use existing channels (failed files, job errors), so no plan review before code (the user's rule, 2026-11-02).

Batch 10e takes the items batches 8 and 9 recorded for later. The user chose (2026-11-02, "go"): part 1 = the target step and the variant check name what they could not do; part 2 = the coverage commit made all-or-nothing, with a plan review first. Two more silent cases of the same kind were found while mapping part 1, and are shown to the user before the code: a BAM that can't be staged in a batch, and a BED line that isn't a BED row.

## Gate A

1. **Principle:** when a step could not do part of its work, ROBIN says which part. It never reports success, or "nothing found", for work it didn't do. (Batch 7's rule, for two steps it missed.)
2. **What would kill it:** any of these, after the fix.
   - A BAM whose coverage was not added, or whose reads are not in `target.bam`, with no failed file.
   - A failed file where nothing failed.
   - A region that hit an error, missing from the results file, or the job reporting success.
   - The GUI saying "completed successfully" for results with a failed region.
3. **Refuted before here?**
   - Batch 8 recorded: "an unreadable staged coverage file drops its BAM and still reports success".
   - Batch 10b keeps such a set in its claim folder (`_keep_incomplete_claim`). Finalize lists it, and SNP runs refuse. But nothing is shown when it happens.
   - Batch 9 recorded: "any other error inside a region still drops the region with only a log line".
   - No mechanism here was tried and dropped before.
4. **Simplest thing:**
   - The accumulation returns the BAMs it could not add, as failed-file entries. Both handler paths pass them on, through `file_failures`, as `bam_preprocessor` does for `master.csv`.
   - The variant check keeps a list of failed regions. The results file holds it, the handler turns it into a job error, and one helper gives the GUI's status text.
5. **Inputs** (read in code at `fe10f4b`):
   - **One-BAM path** (`target_handler`, `target_analysis.py:2663-2667`): the accumulation's result is only logged and stored as metadata. `load_failed`, `lock_timeout`, `error` and `batch_bam_error` never reach the job.
   - **Accumulation** (`accumulate_staged_files`):
     - An unreadable staged set is skipped with a warning (`:2110-2112`), and the result still says "success" (`:2186`).
     - A batch BAM that can't be written gives "success" with `batch_bam_error` (`:2144-2149`). The batched path checks only `status`, so it is silent too (`:2404`).
     - A failure after the claim keeps the claim folder (finalize lists it) and returns `lock_timeout` or `error`.
     - A lock timeout before the claim leaves the files staged; the next accumulation takes them.
   - **Batched path** (`process_multiple_files`, `:2369-2388`): a BAM whose staging fails is skipped with a warning; the job reports success if any other BAM was staged.
   - **Failed files:** `{"file", "error", "step"}` entries in `context.metadata["file_failures"]`. The Coordinator shows one row per entry. The ledger marks a BAM's step `file_failed` only for the job's own BAMs (`record_step_end`, owner-tagged).
   - **Variant check:**
     - `analyze_bed_regions` catches each region's error, logs it and drops the region (`lightweight_gene_analysis.py:409-412`).
     - `_load_bed_regions` returns no regions on any error, a bad line included (`:562-565`).
     - The GUI sets "completed successfully!" for any trusted results file (`coverage.py:4836-4839`).
     - ROBIN's own `targets_exceeding_threshold.bed` has no header line.

## Gate B: measured on `fe10f4b` before this lock

The 16 new tests: **13 failed, each on its own assertion; the 3 guards passed** (`before_10e1.txt`, in this session's scratchpad folder).

## Fixes

1. **Target step: a BAM whose coverage was not added is a failed file** (step `target_coverage`).
   - `accumulate_staged_files` returns `file_failures`, one entry per BAM. The file is the BAM the staged set came from (`source_bam_<n>.txt`), or the staged coverage file when that can't be read. An entry is added:
     - for each staged set that can't be read: its coverage is not counted;
     - for each claimed set, when the accumulation fails after the claim (a lock timeout, an error, every set unreadable): its coverage was not added;
     - for each counted set, when its batch BAM can't be written or moved into place: its coverage is counted, but its reads are not in `target.bam`.
   - Each text names the claim folder its files are kept in.
   - Moving the batch BAM into place gets its own `try`: a failure there is treated like a batch BAM that couldn't be written (the claim is kept), not as a failed commit.
   - `process_multiple_files` adds an entry for each BAM it could not stage, plus the accumulation's entries, when the job goes on.
   - Both paths of `target_handler` add the entries to `file_failures`, keeping any already there.
2. **Variant check: a region it could not check is named.**
   - `analyze_bed_regions` keeps `regions_failed`. Each entry has the region, its ClinVar genes, its intervals and the error.
   - `export_json` writes `metadata.regions_failed`, always: an empty list when none failed.
   - The handler writes the results file, then adds one job error naming the failed regions. The job counts as failed, and the other genes' results stay.
   - `_load_bed_regions` raises on a line with three or more tab-separated fields whose start or end is not a whole number, naming the line. Other errors raise too. The handler's existing `except` turns that into a job error.
   - New `lga_status_text(path)`: "could not be checked" with the region names when any failed, else the old "completed successfully" text. A results file written before this batch has no list and reads as before.
   - The GUI's `_check_lga_results` sets its status from `lga_status_text`.

## Rulings

- **A failed file, not a failed job, on the one-BAM path.** The job's own BAM was staged; it is the totals that missed a BAM. On the batched path, an accumulation that fails outright still fails the job (unchanged). Its failed files are added only when the job goes on.
- **A lock timeout before the claim is not a failure:** the files stay staged, and the next accumulation or finalize takes them.
- **A batch BAM that can't be moved into place** now leaves the job going on, with failed files, instead of failing the batched job. The coverage was counted, so it is the same case as a batch BAM that couldn't be written.
- **A region failure fails the job, and the results file is still written.** The variant check is per sample, so the BAM ledger is not involved.
- **BED lines:** lines with fewer than three tab-separated fields, and `#` lines, are still skipped (unchanged).
- **`analyze_genes`** has no caller in `src/`, and is left as it is.

## Tests (new)

`tests/test_target_failed_files.py`:

| # | Test | Before (`fe10f4b`) |
|---|------|------|
| F1 | one-BAM path: an unreadable staged set is a failed file | fail: none |
| F2 | one-BAM path: a failed commit is a failed file | fail: none |
| F3 | every staged set unreadable: the result names each BAM | fail: none |
| F4 | batched path: a batch BAM that can't be written is a failed file, for each BAM | fail: none |
| F5 | one-BAM path: a batch BAM that can't be moved into place is a failed file, and its claim folder stays | fail: none |
| F6 | batched path: a BAM that can't be staged is a failed file | fail: none |
| G1 | clean runs, both paths: no failed file, no error | pass |
| G2 | a lock timeout before the claim: no failed file, the set stays staged | pass |

`tests/test_lga_regions_failed.py`:

| # | Test | Before (`fe10f4b`) |
|---|------|------|
| L1 | a failed region is listed; the other genes are kept | fail: no list |
| L2 | a clean check lists no failed region | fail: no list |
| L3 | a BED line that isn't a BED row stops the check | fail: did not raise |
| L4 | the handler fails the job and still writes the other genes | fail: no error |
| L5 | the GUI status names the regions not checked | fail: no helper |
| L6 | the GUI status for a clean, and for an older, results file | fail: no helper |
| L7 | the GUI's results check uses the status text | fail: not used |
| G3 | the handler reports no error for a clean check | pass |

## Pass rules (locked)

- **S1:** on `fe10f4b`, the 13 fail and the 3 guards pass (measured).
- **S2:** after the fixes, both new files pass, and so do `test_target_finalize.py`, `test_snp_every_read.py`, `test_panel_depth.py`, `test_target_coverage.py`, `test_clinvar_significance.py`, `test_lga_variant_support.py`, `test_lga_sample_inputs.py`, `test_file_failures.py` and `test_bam_ledger.py`.
- **S3:** the full suite passes: 875 + 16 = 891.
- **S4:** each mutant below, applied alone, makes its test fail (pytest exit 1). The no-op control N0 leaves them passing.
  - P1: an unreadable staged set gets no entry → F1.
  - P2: a failure after the claim gets no entries → F2.
  - P3: a batch BAM that can't be written gets no entries → F4.
  - P4: a failed move into place is a failed commit again (no try of its own) → F5.
  - P5: the batched path drops its staging failures → F6.
  - P6: the one-BAM handler ignores the accumulation's entries → F1.
  - P7: a failed region is not listed → L1.
  - P8: `regions_failed` is written only when not empty → L2.
  - P9: the handler adds no job error → L4.
  - P10: `_load_bed_regions` swallows errors again → L3.
  - P11: `lga_status_text` ignores the list → L5.
  - P12: the GUI sets its fixed text again → L7.
  - N0: a no-op change → nothing fails.
- **S5:** `robin.yml` and `pyproject.toml` unchanged.

## Result (2026-11-02, on `3ab033c`): supported

Every pass rule holds. Three tests were added after the lock, one for each change found while coding; S3's count grows by them. The output files are in this session's scratchpad folder.

- **S1, before** (`fe10f4b`, `before_10e1_full.txt`):
  - The 16 planned tests: 13 failed, each on its own assertion, and the 3 guards passed, as locked.
  - The 3 added after the lock: L8 and F7 failed, and the guard G4 passed.
- **Found after the lock**, each with its own test and mutant:
  - **The GUI's own run** of the check (`run_lga_analysis_direct`) threw its job errors away (`add_error` was a no-op), then said "completed successfully!" in green. It now keeps them, and shows the last one in red. Test L8 (`test_the_gui_run_button_keeps_the_jobs_errors`), mutant P13.
  - **Finalize's own accumulation** only kept its result, so its BAMs that weren't added showed nowhere until an SNP run refused. It now passes the entries on. Test F7 (`test_finalize_shows_a_bam_its_accumulation_could_not_add`), mutant P14.
  - **A staged set counts as taken only once its files begin to move**, so a claim that fails before that reports nothing; its files stay staged. Guard G4 (`test_a_claim_that_fails_before_any_file_moves_reports_nothing`), mutant P15.
- **S2, after:** the 11 files: 259 passed (`after_10e1.txt`).
- **S3, full suite:** 894 passed: 875 + 19 (`suite_10e1.txt`).
  - The first full run had 2 failures in `test_lga_status_timer.py`. That test compiles `check_lga_status` on its own, and the new red-or-green choice read `lga_errors` from the function around it.
  - Fixed in the code, before the commit: the flag rides in `lga_status`. No test changed.
- **S4, mutants** (`mutate_10e1.py`; only pytest exit 1 counts as caught; `mutants_10e1.txt`): P1–P15 each caught by its own test. N0 survived (10 and 9 passed).
- **S5, pins:** `git diff fe10f4b -- robin.yml pyproject.toml` is empty.
- **Recorded, not fixed:** the GUI's SNP button also builds its job with a no-op `add_error` (`coverage.py`, the SNP section).
