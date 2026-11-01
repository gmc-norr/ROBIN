---
title: "A failed fusion build is kept in a .failed file, and a quiet build must rebuild while it is there"
author: Pär
date: 2026-11-01
branch: fix/batch-7-failures-visible
cost: { minutes: 162, tokens: 795153 }
id: l-501cc8f6
scope: repo
status: active
src: agent
touches:
  - reviews/plans/2026-11-01-batch-7-failures-visible.md
  - src/robin/analysis/fusion_work.py
  - src/robin/reporting/sections/fusion.py
  - src/robin/workflow_ray.py
  - tests/test_fusion_report_states.py
  - tests/test_gui_errors.py
  - src/robin/fusion_status.py
  - src/robin/gui/components/fusion.py
  - src/robin/bam_submission.py
  - src/robin/workflow_simple.py
---

## Decided
A fusion build that fails leaves a `.failed` file next to that type's processed pickle (target panel = master, genome-wide = all). The report and the GUI read it through one shared helper: a type with the file and no candidates reads "could not be read"; with earlier candidates they show with "The latest … fusion build failed; these results are from an earlier run." The file is removed only by a build that returns without raising.

## Watch out
A build with no new valid reads does not rebuild the pickle by itself. It must rebuild from the candidate CSV whenever the `.failed` file is there, or it clears the file over a stale pickle: the CSV was appended before the preprocessing failed. That bug shipped in one round. Its test passed because its failing build raised before the CSV append, so the pickle happened to match. Test such paths by comparing the pickle's gene pairs with the CSV's, not just whether the note is gone.

A fusion job that fails before its builds (staging, the parquet merge) writes no `.failed` file, so an earlier pickle still reads as current then. That is known and not done.

## In the code
- `src/robin/fusion_status.py`: `failed_build_marker`, `candidate_type_state`, `earlier_result_note`, `fusion_states`.
- `src/robin/analysis/fusion_work.py`: `_record_build_outcome`, and `_rebuild_missing` inside `_build_filtered_candidates`.
- Tests: `tests/test_fusion_report_states.py`, `test_a_quiet_build_*`.
