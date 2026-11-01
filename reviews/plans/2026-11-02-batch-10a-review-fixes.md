# Batch 10a review fixes: an accepted request never waits for room, and each SNP request is kept as sent

**Locked before code, 2026-11-02.** Base `f77bc5d` (batch 10a's code). Branch `fix/batch-10a-dropped-requests`. Done by the controller (no agents). The findings are from the external review of batch 10a (`/home/parlar_ai/robin-batch3-run/review-requests/batch10a-report-codex.md`). Each was checked in the code before this plan, and each test below was run on `f77bc5d` before the lock.

## Rule (unchanged: the user's option "A" from batch 10a)

For each sample and job type, one runs and one waits; nothing is dropped; a merged run starts after the newest request. This plan makes the code keep that rule in two places where the review showed it doesn't.

## Findings and verdicts

1. **Critical, accepted:** a waiting SNP request is started from inside the finished job's `_on_finish` as a new submission, so it first waits for room in the backlog (`_wait_for_global_capacity`). A pool keeps its slot busy until `_on_finish` returns (`Pool.enqueue`'s `_wait`, `workflow_ray.py:3784-3790`), and the `high` scheduler handles completions one at a time (`_drain_loop`). So if the backlog can only shrink through a job that waits for that slot, nothing moves again. The same wait can happen one step later, for room in a full queue's waiting list (`_wait_for_queue_capacity`), and `_drain_waiting_jobs` (also called from `_on_finish`) can hit that one already, without this batch.
2. **Important, accepted:** with the sample batcher on (the default), SNP requests sent within 2 s become one batch, which runs with its **first** request's context (`_convert_batched_jobs`: `contexts[0]`). The SNP handler reads only that context, so a newer request's settings (forced regeneration, threads, reference) are lost. And a request still in the batcher when the running job ends joins later, behind the request just started: three runs where the rule gives two.
3. **Minor, no change:** if a classifier type's processor were missing, its oldest waiting sample would be dropped. Processors are set only in `setup()` (`workflow_ray.py:1507, 1525, 1544`) and never removed, and the start path itself skips a type with no processor. The reviewer found no way to reach it either. Recorded as a known limit.

Also noted by the reviewer, older than this batch, not fixed here: a classifier run started from `_on_finish` never adds to `inflight_by_queue` / `inflight_by_type`, but its end subtracts from them, so the `slow` queue's running limit is loosened. Recorded for a later batch of step 8.

## Gate A

1. **Principle:** a finish step never waits for something only a later finish can free; and the request waiting for a sample is the newest one, as it was sent.
2. **What would kill it:** a finish step that still waits (more than 2 s) with a full backlog; a follow-up SNP run with an older request's settings or id, with batching on; three runs from requests all sent during one run.
3. **Refuted before here?** Batch 10a's plan said a waiting request is the newest, with the newest settings, but its test S2 left out the id check with batching on ("the batcher re-wraps jobs under new ids") without measuring what the batcher does to the settings. That gap is finding 2 (case-file entry to add with the records). No entry is about waiting inside `_on_finish`.
4. **Simplest thing:**
   - Finding 1: the promotion passes `from_waiting=True` (already accepted: no global-room wait), and `_dispatch_ready_job` with `from_waiting=True` no longer waits for queue room either: a full queue's waiting list takes the job anyway.
   - Finding 2: requests of the per-sample types (`DEDUP_TYPES`) skip the sample batcher in `submit_jobs`, as `preprocessing` already does. Each request then meets the one-running / one-waiting rule as sent. Choosing `contexts[-1]` instead would keep the newest settings within one batch but not fix the three runs, so it is a mutant below, not the fix.
5. **Inputs (read in code at `f77bc5d`):**
   - `_dispatch_ready_job` waits in two places: for global room unless `from_waiting` (1603-1604), and for room in a full queue's waiting list (1643-1644), whatever `from_waiting` says. The other limits append to a waiting list without waiting.
   - Its callers inside `_on_finish`: `_drain_waiting_jobs` (3049, `from_waiting=True`), the promotion (3064, a new submission), and `submit_jobs` for triggered jobs (2813, 2860, 2865, new submissions, out of scope: they are new work).
   - `SampleJobBatcher`: `snp_analysis` batches up to 20 requests, or after 2 s without one (`BATCH_CONFIG`). The batch keeps its own id, which `dedup_running` holds (the reviewer checked this too).
   - The SNP handler (`_wrap_real_handler(_snp_analysis_handler, …)`) reads `job.context`, the primary context only. Nothing in `workflow_ray.py` reads `_batched_job` back.
   - The four classifier types are in `DEDUP_TYPES` too, but their runs come from `_on_finish`'s bed_conversion path, not `submit_jobs`, so skipping the batcher changes nothing for them.
   - The CNV batching inside `_on_finish` (2806-2857) is for `cnv` only and is untouched.

## Gate B: the tests, run on `f77bc5d` before this lock

**New file `tests/test_job_requests_direct.py`.** These tests call the Coordinator's own methods directly: the class behind the Ray actor, built without `setup()`, so no cluster, pools or background loops start. Processors only record what they are sent. The backlog, the limits and the batcher's timeout are set exactly, not by timing. The limits are small (one waiting job in all); the defaults are large, but the cycle is the same.

| # | Test | On `f77bc5d` |
|---|---|---|
| U1 | `test_a_waiting_snp_request_starts_without_waiting_for_room`: SNP S runs, a second waits; an igv_bam job runs and one more fills the backlog (limit 1). S's finish must return within 2 s and start the waiting request. | failed: "finishing snp_analysis waited 2.0 s for room" |
| U2 | `test_a_job_from_a_waiting_list_does_not_wait_for_queue_room`: a job taken off the global waiting list finds its queue at its running limit and the queue's waiting list full. The finish step must return, and the job must wait in that queue's list. | failed: "finishing mgmt waited 2.0 s for room" |
| U3 | `test_batched_snp_requests_keep_the_newest_request`: batching on; two requests during S's run, the newer with `force_regenerate=True`. The follow-up run is the newer one. | failed: `['first', 'older'] == ['first', 'newer']` |
| U4 | `test_batched_snp_requests_during_a_run_give_one_more_run`: batching on; request 2 has left the batcher for the waiting slot, request 3 is still in it when run 1 ends. Two runs in all, the second request 3. | failed: `['r1', 'r2', 'r3'] == ['r1', 'r3']` |
| S2 | `test_snp_requests_while_one_runs_give_one_more_run` (Ray, batch 10a's test): the follow-up run's id is the newest request's **with batching on too** (the exemption is removed: a stronger test, not a weaker one). | batched failed: `1001 == 1000005` (the batch's id); unbatched passed |

On `f77bc5d` the four direct tests failed in 7 s, each after its own preconditions held. The rest of batch 10a's tests, and the existing scheduling and failure tests, must still pass.

## Mutants (each must be caught; the no-op control must survive; only exit 1 counts)

- **N1:** the promotion is a new submission again (no `from_waiting=True`).
- **N2:** a job with `from_waiting=True` waits for queue room again.
- **N3:** per-sample request types go through the batcher again.
- **N4:** they go through the batcher, but a batch runs with its newest context (`contexts[-1]`): the half fix.
- **Control:** a no-op change.

## Pass, all of

- U1–U4 fail on `f77bc5d`, each by its own assertion (U1 and U2 by "waited 2 s for room"), and pass after. S2 with batching on fails on `f77bc5d` and passes after.
- N1–N4 are each caught.
- The full suite passes.
- `robin.yml` and `pyproject.toml` are unchanged.

## Out of scope

- Finding 3 (no change, above), and the loosened `slow` queue limit (later in step 8).
- Triggered jobs submitted from `_on_finish` (`submit_jobs` at 2813, 2860, 2865) still wait for global room, as new work. On `high`, where completions are handled one at a time, any such wait holds up every other completion, so these may block the same way. Older than this batch and not measured here; for a later batch of step 8.
- The reviewer's note on timing: the Ray tests use a 6 s busy job and check that their requests arrived during it, so a slow machine fails them loudly rather than passing them wrongly. Gates in place of sleeps would make them exact; not done here.
- Later batches of step 8 (10b–10e), and the `high` freeze.
