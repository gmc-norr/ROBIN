---
title: "A sample's CNV counts are read from disk under a per-sample lock, never cached in a worker"
author: Pär
date: 2026-11-02
branch: fix/batch-10c-cnv-counts
cost: { minutes: 15, tokens: 158773 }
id: l-924ecef3
scope: repo
status: active
src: agent
touches:
  - src/robin/analysis/cnv_analysis.py
  - src/robin/workflow_ray.py
  - tests/cnv_fixtures.py
  - tests/test_cnv_control_reference.py
  - tests/test_cnv_counts_shared.py
  - reviews/plans/2026-11-02-batch-10c-cnv-counts.md
  - CHANGELOG.md
  - reviews/2026-10-16-code-review.md
  - reviews/2026-10-16-code-review.pdf
---

## Decided
- No per-worker cache of a sample's CNV counts or counter. Each Ray task can land on a different worker process, so a cache keyed on the path alone goes stale, and a batch another worker added is lost when the stale copy is saved. Every run reads `<sample>_copy_numbers.pkl` from disk.
- A per-sample file lock (`<sample>/_locks/cnv.lock`, the `FileLock` from master_bed_generator) is held from reading the counts to saving them. It wraps `process_single_bam` and `process_multiple_bams`, which keep their old bodies (`_process_*`). The wait is 3 h, not measured; past it, the run fails with an error_message.
- Unreadable counts raise, so the run fails and the file is kept. They used to reset to `{}`, which saved over every earlier BAM. Writes go through `atomic_open`.

## Tried and ruled out
- Serialising CNV in the scheduler with `SERIALIZE_BY_TYPE_PER_SAMPLE`: the set is empty, so that code has never run, and it would not cover a CNV run outside Ray.

## Watch out
- `workflow_ray.py` imports `cnv_handler` in the same try-block as any other names from cnv_analysis. A name removed from cnv_analysis but left in that import makes the whole import fail, and CNV then has no handler. `test_the_workflow_still_imports_the_cnv_handler` guards it.
- To test worker staleness, use real separate processes (multiprocessing "spawn", a long-lived worker loop in `tests/cnv_fixtures.py: cnv_worker`). Threads in one process share the module state, so they can't show it.

In the code: `src/robin/analysis/cnv_analysis.py` (`_load_copy_numbers`, `_save_copy_numbers`, `_locked_cnv_run`). Tests: `tests/test_cnv_counts_shared.py`.
