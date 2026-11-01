---
title: "Job failure tests drive a real Ray Coordinator, with each fault switch in the job's metadata"
author: Pär
date: 2026-11-01
branch: fix/batch-7-failures-visible
cost: { minutes: 126, tokens: 770800 }
id: l-0b875f98
scope: repo
status: active
src: agent
touches:
  - tests/test_fusion_report_states.py
  - src/robin/analysis/fusion_work.py
  - src/robin/analysis/fusion_analysis.py
  - src/robin/reporting/sections/fusion.py
  - src/robin/analysis/clinvar_status.py
  - src/robin/analysis/snp_processing.py
  - src/robin/analysis/target_analysis.py
  - src/robin/gui/components/snp.py
  - src/robin/reporting/sections/variants.py
  - tests/test_clinvar_status.py
---

## Decided
- Tests of whether a job counts as failed run a real local Ray Coordinator (enable_batching=False) on both schedulers, with synthetic handlers. Each scheduler judges a job in its own place (the pool actors for `standard`, the type processors for `high`), so a mocked judge would test neither.
- A job's outcome comes from the errors filed while it ran, whatever name they were filed under. Under Ray that is by position (errors past `errors_before`; each task has its own copy of the context). The non-Ray runner tags each error with the job that filed it, because sibling jobs share one context object on different threads and a positional rule let one sibling's error fail another.

## Watch out
- Ray workers don't see the driver's `sys.argv` or globals a test sets. A fault switch read from there does nothing, and the job "succeeds". Put the switch in the job's metadata. A control (a raising handler must count as failed) is what caught this in the first probe.
- Give every job a `work_dir` in its metadata, or the Coordinator writes `job_durations.csv` into the current folder, which is the repo.
- `failed_jobs.csv` is written only when the Coordinator shuts down. Read it after shutdown.
- A Ray worker has no GUI. `send_gui_update` there only fills a small outbox in that process; the handler wrapper moves it into the job's metadata, and the Coordinator forwards it. A message arrives when its job ends, not before.
- These tests start Ray, so a mutation check over them takes about 2 minutes a mutant. Run Ray test runs one at a time: a slow run that times out looks like a caught mutant.
- Before locking a control as "passes on the old code", run it there. Parts of three controls in batch 7 were new behaviour (the expected_failure label, the ClinVar status line, the empty fusion pickle).

In the code: tests/failure_fixtures.py (the driver: `Harness`, fault switch in `ctx.metadata["fault"]`), src/robin/job_outcome.py (`job_outcome`, `record_outcome`, `recorded_failure`), src/robin/workflow_ray.py (`_wrap_real_handler`, `Coordinator._on_finish`, `publish_failures`), src/robin/workflow_simple.py (`WorkflowManager.worker`, the owner tag).
