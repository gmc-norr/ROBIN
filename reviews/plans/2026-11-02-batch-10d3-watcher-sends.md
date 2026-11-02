# Batch 10d-3: the Ray watcher sends every BAM it takes in, and says so when it can't (JOB-4)

**Locked before code, 2026-11-02.** Base `b477842`. Branch `fix/batch-10d-file-ledger`. Part 3 of batch 10d. Done by the controller (no agents). An external review covers round 2 of part 2 and this part before the merge.

## Rule (agreed by the user, 2026-11-02, "go")

1. The watcher sends waiting jobs on a timer (every 0.5 s), not only when the next file arrives.
2. A send that fails is logged with the BAMs' names and tried again by the timer. Nothing is dropped.
3. At shutdown, ROBIN stops watching and sends what waits, and only then stops the Coordinator. What still can't be sent is logged by name; the BAM ledger runs those BAMs at the next start.
4. The timer checks half-written BAMs again. One that completes is taken in, even with no new file event. One that stays incomplete and unchanged for 10 minutes shows once as a failed file ("incomplete, not analysed"); if it completes later, it is still analysed.

**Spec:** JOB-4's first three bullets in `reviews/2026-10-16-code-review.md`, and "a BAM that never completes is only logged" from part 2's Not done.

## Gate A

1. **Principle:** every BAM the watcher takes in reaches the Coordinator, or the user sees which ones didn't.
2. **What would kill it:** any of these.
   - A taken-in BAM whose jobs wait with no further file event and are never sent.
   - A failed send whose jobs are dropped, or not logged by name.
   - Waiting jobs sent after the Coordinator was told to stop (it ignores them).
   - A half-written BAM that completes with no new file event and is never taken in.
   - A BAM that never completes and is never shown, or is shown while it is still being written, or shown more than once.
3. **Refuted before here?** Case 15 (round 2): list every layer between taking a BAM in and its step running. Here that is the watcher's buffer and the shutdown order. Learning `l-294c2aaa`: a BAM taken in but never started runs at the next start (the ledger); this part relies on that for what can't be sent at the end.
4. **Simplest thing:** one daemon thread per watcher, ticking every 0.5 s: send what waits, and check the intake's waiting BAMs again. `stop()` ends it and sends the rest. One helper in `run()` stops watching first, on both shutdown paths.
5. **Inputs** (read in code at `b477842`):
   - **The buffer:** `RayFileWatcher._flush_if_needed` (`workflow_ray.py`) is called only from `_handle`. It sends when 128 jobs wait or 0.5 s passed since the last send. So jobs added within 0.5 s of a send wait for the next file event.
   - **Failed sends:** it removes a batch from the buffer, then sends it. An exception is swallowed (`except Exception: pass`), so the batch is gone, and its path is already marked processed.
   - **Shutdown order** (`run()`, `finally`):
     - It calls `coord.shutdown` first, which sets `_shutdown_requested`.
     - `submit_jobs` returns at once when that is set, so the watcher's flush that follows sends into nothing.
     - The observer is stopped last, so new events still come in meanwhile.
     - The Ctrl-C path kills the Coordinator before the `finally` runs.
   - **Half-written BAMs:** the intake keeps them in `_waiting`, and they are decided again only on a file event. A BAM the start-up scan found half-written, or whose writer stopped, gets no further event.
   - **The non-Ray runner** has no buffer (it hands jobs to its own queue), so it is not changed here.

## Gate B: measured on `b477842` before this lock

The tests below: **9 failed** (`before_10d3.txt`, in this session's scratchpad folder).
- 3 failed on their own assertion: W2 (the failed send is not logged), W3 (the unsent BAM is not named), W4 (the Coordinator is stopped first).
- 6 failed because the timer doesn't exist yet (`RayFileWatcher.start`): W1, W5, W6, W7, and the two controls C7 and C8. The controls need the timer to run at all; they must pass after, and each is the target of a mutant.

## Fixes

1. **The timer.**
   - `RayFileWatcher.start()` starts one daemon thread. Every `_tick_s` (0.5 s) it calls `tick()`: check the waiting BAMs again, then send everything that waits.
   - `stop()` ends the thread, then sends what waits (one pass), and logs what is left.
   - A lock (`RLock`) guards the buffer, the processed set and `_handle`, which now run on two threads.
2. **Sending.**
   - A batch that fails to send goes back to the front of the buffer, and is tried again on the next tick.
   - Each BAM's first failure is logged as an error with its name and the reason.
   - `flush_remaining()` makes one pass. It logs any BAM still unsent by name, with what happens next: "the BAM ledger runs it at the next start", or, with no work folder, "it is not analysed". It returns the unsent jobs.
3. **Shutdown order.** A helper `_stop_watching(observer, watcher)` stops the observer, then the watcher (sending what waits). `run()` calls it first on the Ctrl-C path and in `finally`, before `coord.shutdown`. It does nothing the second time. `run()` starts the watcher's timer right after the observer.
4. **Half-written BAMs.** The intake keeps, for each waiting path, its size and modification time and when they last changed.
   - `waiting_paths()` lists them; the timer passes each to `_handle` again.
   - `incomplete_too_long(seconds)` returns those unchanged for that long and not yet reported, and marks them reported.
   - The watcher reports each once through `record_cut_short_steps`: `failure_kind` "incomplete", job type "preprocessing", text "…has no BAM end-of-file block and has not changed for 10 minutes: not analysed (it will be if it completes)".
   - `_incomplete_after_s` is 600.

## Rulings

- **A BAM shown as incomplete that completes later** is taken in and analysed. Its earlier failed-file row stays in the list.
- **The incomplete row's sample** is "unknown": the BAM was not read.
- **The non-Ray runner** is unchanged: it has no buffer, and its half-written BAMs are only logged, as before.
- **Retries of a failed send** happen every tick for as long as the session runs. Only each BAM's first failure is logged.

## Tests (new file `tests/test_watcher_sends.py`)

| # | Test | Before (`b477842`) |
|---|------|------|
| W1 | the last BAM is sent without another file event | fail: no timer |
| W2 | a failed send is logged and tried again | fail: not logged (and dropped) |
| W3 | what can't be sent at the end is logged, and runs at the next start | fail: not named |
| W4 | the Ray run stops watching before it stops the Coordinator (both shutdown paths, from `run()`'s source) | fail: Coordinator first |
| W5 | stop sends what waits and ends the timer | fail: no timer |
| W6 | a BAM the scan found half-written is taken in once it completes, with no file event | fail: no timer |
| W7 | a BAM still incomplete and unchanged after the set time shows once as a failed file | fail: no timer |
| C7 | a BAM still growing is not shown as incomplete (control) | fail: no timer |
| C8 | a BAM taken in is sent once however many ticks pass (control) | fail: no timer |

## Pass rules (locked)

- **P1:** on `b477842`, the 9 fail as shown (measured).
- **P2:** after the fixes, all 9 pass, and so do `test_bam_ledger.py`, `test_job_failures.py` and `test_file_failures.py`.
- **P3:** the full suite passes: 859 + 9 = 868.
- **P4:** each mutant below, applied alone, makes its test fail (pytest exit 1). The control N0 leaves them passing.
  - Y1: the tick doesn't send → W1.
  - Y2: a failed batch is not put back → W2.
  - Y3: unsent BAMs are not logged → W3.
  - Y4: `run()` stops the Coordinator before `_stop_watching` → W4.
  - Y5: `stop()` doesn't send what waits → W5.
  - Y6: the tick doesn't check waiting BAMs again → W6.
  - Y7: no incomplete report → W7.
  - Y8: the unchanged-since time is never reset → C7.
  - Y9: the send keeps sent jobs in the buffer → C8.
  - N0: a no-op change → nothing fails.
- **P5:** `robin.yml` and `pyproject.toml` unchanged.

## Result (2026-11-02, on `39039b7`): supported

Every pass rule holds. The output files are in this session's scratchpad folder.

- **P1, before** (`b477842`, `before_10d3.txt`): the 9 failed as shown.
- **P2, after:** all 9 pass. With `test_bam_ledger.py`, `test_job_failures.py` and `test_file_failures.py`: 250 passed (`after_10d3.txt`). The 9 passed in each of 5 more runs (the timer tests depend on timing).
- **P3, full suite:** 868 passed (`suite_10d3.txt`).
- **P4, mutants** (`mutate_10d3.py`; only pytest exit 1 counts as caught; `mutants_10d3.txt`): Y1–Y9 were each caught by the test the plan names. N0 survived (9 passed).
- **P5, pins:** `git diff b477842 -- robin.yml pyproject.toml` is empty.

### Ruling made after the lock: W4 checks each call's full name

W4 named each call by its last part. `coord.shutdown.remote()` is then a call named `remote`, so W4 would have missed the Ctrl-C path's Coordinator stop.
- **The change:** W4 now names each call by its full dotted name (`ast.unparse`), and checks that a call to `_stop_watching` comes before any call with "shutdown" in it, in each `except` and `finally` block of `run()` that has one.
- **Before the code:** the corrected W4 was run again on `b477842`'s code and still failed for its reason (`before_10d3_w4.txt`).
- **What stays:** its rule (watching stops before the Coordinator, on both paths) is unchanged.
