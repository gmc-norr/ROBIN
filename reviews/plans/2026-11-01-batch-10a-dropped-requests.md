# Batch 10a: a request for a classifier or SNP run is never thrown away (JOB-1, JOB-5)

**Locked before code, 2026-11-01.** Base `0545e2f`. Branch `fix/batch-10a-dropped-requests`. Done by the controller (no agents). Step 8 of the review is split into five batches (10a–10e). This is the first: the dropped requests. An external review comes before the merge.

## Rule (agreed by the user, option "A", 2026-11-01)

- **One running, one waiting:** for each sample and job type, at most one job runs and at most one waits.
- **Nothing is dropped:** a new request joins the waiting one if there is one. A merged run starts after the newest request arrived.
- **Classifiers stay one at a time:** each classifier type (Sturgeon, NanoDX, PanNanoDX, Random Forest) still runs one job at a time across all samples, as before.
- **Waiting samples take turns,** oldest first.

**Spec:** JOB-1 and JOB-5 in `reviews/2026-10-16-code-review.md`, and the user's choice.

## Gate A

1. **Principle:** every request ROBIN accepts leads to a run that starts after it; only requests that would repeat a waiting run are merged.
2. **What would kill it:** any of these.
   - A sample with new methylation data whose last classifier run started before that data's bed_conversion ended.
   - A second SNP request that never runs.
   - Merged requests that each run separately: three new batches during a run giving three more runs, not one.
3. **Refuted before here?** No case-file entry is about scheduling. Batch 7 (ERR-1, JOB-3) touched `_on_finish` and the job wrapper, not these queues. Learning `2026-11-01-job-failure-tests-drive-a-real-pa`: test through a real Coordinator, with fault switches in the job's metadata.
4. **Simplest thing:**
   - JOB-1: the one global waiting slot per classifier type becomes one waiting job per sample, held in arrival order.
   - JOB-5: the per-sample request counter that is never decremented becomes "the request running now" (its job id) plus "the one waiting request".
5. **Inputs (read in code at `0545e2f`, measured below):**
   - **Classifier runs are triggered** when a bed_conversion finishes. `_on_finish` starts them straight on the type's processor, with one global "pending" per type and one global waiting `Job` per type (`classif_waiting_by_type`). A request finding that slot taken is dropped (no else branch).
   - **Classifier input:** a classifier reads the sample's whole merged file. Its context's `bed_conversion.parquet_path` is `<work_dir>/<sample>/<sample>.parquet` (`bed_conversion.py:411`). So a merged request loses nothing, as long as the run starts after the newest data.
   - **SNP requests** (`submit_snp_analysis_job`, and the generic `submit_sample_job`) go through `submit_jobs` (via the sample batcher when batching is on, as in real runs), then `_submit_jobs_internal`. For `DEDUP_TYPES` that adds one to `pending[(sample, type)]`.
   - Only `mark_running` takes it off again, and it has no callers. So `pend >= 1` stays true for good, and every later request is skipped (`total_skipped += 1`), while `submit_snp_analysis_job` returns True and logs "queued".
   - `_on_finish` takes one off `running[(sample, type)]`, which nothing ever adds to. It does this for every finished job of these types, triggered classifier runs too.
   - `_on_finish` has no early return, so code at its end runs for every finished job, failed ones too.

## Gate B: measured on `0545e2f` before this lock

The tests below were run on `0545e2f` (results in the table). The harness (`failure_fixtures.py`) gained four things:
- a per-type `sleep` in the job's metadata;
- a "started" record and start/end times in each job's outcome record;
- an `enable_batching` option (default off, as before);
- `num_cpus=0` on the synthetic handler tasks. With 4 test CPUs, part held by the Coordinator's actors, jobs of different types queued for a CPU one at a time, which a real run doesn't do.

The existing harness tests are unchanged.

**Three test-design faults, found by running on `0545e2f`, fixed before this lock:**
1. **The busy job was too short.** With 2 s, every "request during the run" landed after it, because a fresh Ray worker takes about 2.2 s to start (`b10/probe_block.py`). So a guard failed on the old code (3 and 4 runs, not 2), and two fault tests passed there. Busy jobs now take 6 s. Each test first asserts that its requests arrived before the busy run ended, so a timing slip fails loudly instead of testing the wrong thing.
2. **Jobs queued for a CPU,** fixed by the `num_cpus=0` change above.
3. **The high scheduler can't hold a request mid-run.**
   - **Measured** (`probe_block.py`, high): while a 4 s Sturgeon job ran, a `stats` call to the Coordinator took 4.00 s. On standard it took 0.01 s.
   - So on high the Coordinator handles nothing while a job runs, and the mid-run tests' own checks failed there.
   - Those tests (J1, J2, G1, S2) run on standard, the default, only. S1 and G2 run on both schedulers.
   - The freeze is a new finding, recorded below for a later batch, and is not fixed here.

## Fixes

1. **Classifiers (JOB-1):** `classif_waiting_by_type` becomes type → {sample: Job}, in arrival order.
   - A request finding its type busy stores its job under its sample. That replaces an older waiting job for the same sample, which keeps its place in line.
   - When a classifier job finishes, the first waiting sample's job starts. The start code is the same as today's.
2. **Per-sample requests (JOB-5):** `pending` / `running` become `dedup_running[(sample, type)] = job id` and `dedup_waiting[(sample, type)] = Job`.
   - A request whose key is running waits, and replaces an older waiting request.
   - When the job with that id finishes, the waiting request is dispatched, after the finished job's backpressure release. A triggered classifier run never touches these maps, because its id is never the stored one.
   - `mark_running` (no callers, and it kept the old counters) is removed.

## Tests: `tests/test_job_requests.py`

| # | Test | Params | On `0545e2f` |
|---|---|---|---|
| J1 | `test_every_sample_gets_a_classifier_run_after_its_data`: A runs, then new data for B and C arrives during it. Each sample's last Sturgeon run starts after its last bed_conversion ends. | standard | failed: `{'A', 'B'} == {'A', 'B', 'C'}` (C dropped) |
| J2 | `test_new_data_while_another_sample_waits_gets_its_own_run`: A runs, B waits, then A gets new data. Both get a run after their data, and B's run comes before A's second. | standard | failed: A's last run is older than its data |
| G1 | guard `test_guard_requests_for_one_sample_merge_into_one_run`: three new batches for A during its run give exactly one more run, after the last batch. | standard | passed |
| S1 | `test_a_second_snp_request_runs`: two `submit_snp_analysis_job` calls, the second after the first finished. Both run, and both return True. | 2 schedulers × batching off/on | failed ×4: `1 == 2` |
| S2 | `test_snp_requests_while_one_runs_give_one_more_run`: two requests during a run give exactly one more run. It starts after they were sent, and, without batching, it is the newest request (the batcher re-wraps jobs under new ids). | standard × batching off/on | failed ×2: `1 == 2` |
| G2 | guard `test_guard_snp_requests_for_two_samples_both_run`. | 2 × 2 | passed ×4 |

**On `0545e2f`:** 8 failed, all by assertion on the fault itself (not on a test's own precondition), and 5 passed (the guards). That run used these tests with high as well. The high cases of J1, J2, G1 and S2 were dropped afterwards, for reason 3 above; the remaining cases' results are unchanged.

The existing scheduling and failure tests (`test_job_failures.py`, `test_file_failures.py`, `test_classifier_failures.py`) must still pass.

## Mutants (one per fix part), each must be caught, and the no-op control must survive

- **M1:** a busy classifier type keeps only the first waiting sample: the old drop.
- **M2:** the finished classifier starts the newest waiting sample, not the oldest (turn-taking broken).
- **M3:** a finished per-sample request doesn't start its waiting one.
- **M4:** a request while one runs is dropped, not stored.
- **M5:** a waiting per-sample request is not replaced by a newer one.
- **M6:** the running per-sample request is never released: the old bug.
- **Control:** a no-op change.

Only exit 1 counts as caught.

A waiting classifier job is replaced too, but that can't be seen: every request for a sample reads the same merged file. So it has no mutant.

## Pass, all of

- J1, J2, S1 and S2 pass after. G1 and G2 pass before and after.
- M1–M6 are each caught.
- The full suite passes.
- `robin.yml` and `pyproject.toml` are unchanged.

## Out of scope (later batches of step 8)

- 10b: JOB-7 and the rest of JOB-6.
- 10c: CNV-2.
- 10d: JOB-2 and JOB-4.
- 10e: the items recorded in batches 8 and 9.
- **New finding: on the `high` preset the Coordinator stops answering while a job runs** (measured above). Submissions, and the GUI's status calls, wait until the job ends. `high` is opt-in (the default is `standard`). The cause is not found yet: `_drain_loop` waits only 0.1 s, and `TypeProcessor.process` is async. For a later batch.
- **Manual classifier runs** (`submit_sample_job` with a classifier type) get the per-sample rule but not the one-at-a-time rule. Such a job has no bed_conversion result and fails in the handler ("No parquet paths found"), so nothing in ROBIN sends one.
