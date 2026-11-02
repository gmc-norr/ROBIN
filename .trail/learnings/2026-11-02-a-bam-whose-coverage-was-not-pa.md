---
title: "A BAM whose coverage was not added is a failed file, and the variant check names regions it could not check"
author: Pär
date: 2026-11-02
branch: fix/batch-10e-recorded-items
cost: { minutes: 49, tokens: 437610 }
id: l-5e9df47f
scope: repo
status: active
src: agent
touches:
  - reviews/plans/2026-11-02-batch-10e1-failures-named.md
  - src/robin/analysis/lightweight_gene_analysis.py
  - src/robin/analysis/target_analysis.py
  - src/robin/gui/components/coverage.py
  - tests/test_lga_regions_failed.py
  - tests/test_target_failed_files.py
  - CHANGELOG.md
  - reviews/2026-10-16-code-review.md
  - reviews/2026-10-16-code-review.pdf
  - reviews/plans/2026-11-02-batch-10d-round3.md
---

## Decided
- A BAM whose coverage the target step could not add is a failed file (step `target_coverage`), not a failed job: the job's own BAM was staged, it is the totals that missed a BAM. The accumulation returns the entries; the one-BAM path, the batched path and finalize all pass them on through `file_failures`. On the batched path an accumulation that fails outright still fails the job, and then no entries are added.
- A staged set counts as taken only once its files begin to move into the claim folder. A failure before that leaves it staged for the next accumulation, so it is not reported.
- A region the variant check could not check is listed in the results (`metadata.regions_failed`) and fails the job; the file is still written with the other genes. A BED line that isn't a BED row stops the check.

## Watch out
- A step's "success" can be said in more than one place. The variant check's GUI has two: the results check (`_check_lga_results`) and the GUI's own run (`run_lga_analysis_direct`), which also threw its job errors away with a no-op `add_error`. Search for the success text and for every caller of the changed function before locking a plan.
- The GUI's SNP button builds its job with a no-op `add_error` too (`coverage.py`, the SNP section): not fixed here.

In the code: `accumulate_staged_files`, `_staged_set_bam`, `_install_batch_bam`, `_add_file_failures`, `process_multiple_files`, `target_handler`, `target_bam_finalize_handler` in `src/robin/analysis/target_analysis.py`; `analyze_bed_regions`, `_load_bed_regions`, `export_json`, `lga_status_text` and the handler in `src/robin/analysis/lightweight_gene_analysis.py`; `_check_lga_results` and `run_lga_analysis_direct` in `src/robin/gui/components/coverage.py`. Tests: `tests/test_target_failed_files.py`, `tests/test_lga_regions_failed.py`.
