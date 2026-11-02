---
title: Scheduling tests must prove each request landed while the busy job ran
author: Pär
date: 2026-11-01
branch: fix/batch-10a-dropped-requests
cost: { minutes: 1, tokens: 1485 }
id: l-d2b16fcd
scope: repo
status: active
src: agent
touches:
  - tests/test_job_requests.py
  - tests/failure_fixtures.py
  - "src/robin/workflow_ray.py#Coordinator._on_finish"
  - "src/robin/workflow_ray.py#Coordinator._submit_jobs_internal"
---

## Decided
- Tests of what happens to a request sent while a job runs (JOB-1, JOB-5) make the busy job take 6 s, and each test first asserts that its requests were sent before that job ended. A timing slip then fails loudly on the test's own check, instead of quietly testing a request that came after the run.
- The synthetic handler tasks take no CPU (`num_cpus=0`). With 4 test CPUs, part held by the Coordinator's actors, jobs of different types queued for a CPU one after another, which a real run doesn't do.
- The tests that need a request to land mid-run use the `standard` scheduler only. The others run on both schedulers, with batching off and on.

## Tried and ruled out
- A 2 s busy job: a fresh Ray worker takes about 2.2 s to start, so every "request during the run" landed after it. On the old code a guard failed (3 and 4 runs, not 2) and two fault tests passed, so the tests said the opposite of the truth.
- Running the mid-run tests on `high`: there the Coordinator answers nothing while a job runs. A `stats` call took the job's full 4 s (0.01 s on `standard`), so no request can land mid-run. This is a finding in its own right (noted under step 8 in the code review), not fixed yet.

## Watch out
- Classifier runs (all four, random forest too) are started from `_on_finish` when bed_conversion ends, not through `_submit_jobs_internal`. So the per-sample request rule (`dedup_running` / `dedup_waiting`) only matters for `snp_analysis` in real runs.
- Don't exempt the batching-on case from a check. The sample batcher runs a batch with its first request's context, under a new id. Per-sample requests now skip it (an external review found the lost settings), so with batching on, the newest request's own id still runs.
- For exact limits and timing, call the Coordinator's methods directly: `Coordinator.__ray_metadata__.modified_class`, built without `setup()`, with processors that only record what they are sent (`tests/test_job_requests_direct.py`). No Ray is needed, and those tests take seconds.

In the code: tests/test_job_requests.py (SLOW = 6.0, the "must arrive during" asserts); tests/failure_fixtures.py (Harness: enable_batching, submit, started, ended; per-type sleep in metadata["sleep"]).
