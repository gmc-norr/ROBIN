---
title: "Empty results must still be written, and NiceGUI timers ignore what the callback returns"
author: Pär
date: 2026-10-31
branch: fix/batch-3c-review-findings
cost: { minutes: 26, tokens: 227769 }
id: l-45e0ed45
scope: repo
status: active
src: agent
touches:
  - "src/robin/analysis/cnv_analysis.py#save_cnv_files"
  - "src/robin/analysis/master_bed_generator.py#_load_bed_file"
  - "src/robin/analysis/lightweight_gene_analysis.py#pick_lga_results"
  - src/robin/gui/components/coverage.py
  - "src/robin/command_template.py#shell_command_for_file"
---

## Decided
- When a reader picks "the newest file" (master BED's _get_latest_bed_file, the GUI's cnv_data_array.npy), a run with an empty result must still write its file, empty. Skipping the write on an empty result silently leaves the previous run's data as the newest. CNV breakpoints now always write breakpoints_NNN.bed and the structured array; _load_bed_file treats a zero-byte BED as no regions.
- LGA results can live in two current places (<sample>/<id>/ when the handler is given the sample folder, <sample>/ when given the root). pick_lga_results takes the trusted file with the latest metadata.generated_at, not the first candidate. Content dates, not mtimes, so copying files doesn't reorder them.
- A NiceGUI ui.timer stops only via .cancel() or .deactivate() on its handle. NiceGUI 3.18 uses the callback's return value only to decide whether to await it; "return False # stop" does nothing.

## Tried and ruled out
- Extending the zero-usable-coverage rule to indel sites: a column where every read is deleted is evidence FOR a ClinVar deletion there.
- A two-level nested "$(… "$(… "{file}")")" test for command templates: it passed on the broken one-level quote tracker by luck (odd number of '"' before {file}). Three levels exposes it.

## Watch out
- Any other "return False to stop polling" pattern in GUI code is the same bug; grep for it.
- Tests for timer behaviour: compile the nested callback out of the source with ast and run it in nicegui.timer.Timer's real loop with _skip_registration -> True and a tick bound in _should_stop (tests/test_lga_status_timer.py). Make the simulated user act once, or the test can't see a reset.

In the code: tests/test_cnv_breakpoint_files.py, tests/test_lga_sample_inputs.py (newest-trusted tests), tests/test_lga_status_timer.py, tests/test_lga_variant_support.py (no_usable_coverage), tests/test_command_template.py (context stack).
