# Batch 10d-2: each BAM is taken in once, and after a restart only its unstarted steps run (JOB-2, JOB-8)

**Locked before code, 2026-11-02.** Base `98dc6ef` (batch 10d-1 on top of `01a5482`). Branch `fix/batch-10d-file-ledger`. Done by the controller.

## Rule (agreed by the user, 2026-11-02, "1")

A ledger in the work folder records each BAM ROBIN takes in, and when each of its steps starts and ends.

**After a restart**, for each BAM:
- A step that ended is skipped, whether it succeeded or failed (a failure was shown).
- A step that never started runs.
- A step that was running when ROBIN stopped is not run again, because it may have added its data. It is reported once, as a failed file that "may be incomplete".

**Taking a BAM in:** the start-up scan, the watcher and the GUI's "add folder" all use one intake. So does the non-Ray runner. The Ray runner starts watching before it scans.

**How a BAM is known:** by its bytes (`bam_key`, batch 10d-1), so it is taken in only once it is complete.

**JOB-8:** the watcher's BAMs get the reference genome, as the scan's do.

## Gate A

1. **Principle:** a BAM's per-BAM steps (preprocessing, methylation, MGMT, CNV, target, fusion) each run for it once, across restarts, scans and watchers. A step that may have run part-way is shown, never silently run again or left out.
2. **What would kill it:** any of these.
   - A BAM taken in again after its steps ended: by the next session, by the watcher after the scan, or as a copy.
   - A step that never started not run after a restart.
   - A step that was running at a stop run again, or not reported.
   - A half-written BAM taken in.
   - A watched BAM without the reference.
3. **Refuted before here?** The stricter plan (each counting step keeps its own list) was mapped and set aside: six steps, seven batches (batch 10d-1's result). `SQLiteStateTracker` (`state_tracker.py`) has no callers, and its key includes the size and time, so a growing file gets a new key each time. It is not used.
4. **Simplest thing:**
   - **The ledger:** one append-only JSON-lines file, `<work_dir>/robin_bam_ledger.jsonl`. Each append is under `flock`, then fsynced. One small module, `robin/file_ledger.py`.
   - **The records:** the driver writes "take". The worker writes "start" just before a handler runs and "end" just after, so the record sits next to the work.
   - **A restarted BAM** goes back through its preprocessing, which counts it once (10d-1), with `only_steps` set to the steps to start.
5. **Inputs (read in code at `98dc6ef`):**
   - **Ray: every handler** runs through `_wrap_real_handler`'s `_impl` (`workflow_ray.py:799-977`). Every registered remote is wrapped (`:1017-1033`, `:1378-1390`). The handler call is at `:881`. The fallback result is set in its `finally` (`:946-958`).
   - **Non-Ray: every handler** is called in `WorkflowManager.worker` at `workflow_simple.py:917`.
   - **Batched jobs:** a batched job carries its BAMs' contexts in `metadata["_batched_job"].contexts` (`workflow_ray.py:3168-3184`; non-Ray `workflow_simple.py:1311`). Downstream steps reuse the BAM's preprocessing context (`workflow_ray.py:2764`), so metadata set at take-in reaches them.
   - **Ray triggers:** preprocessing starts `TRIGGERS["preprocessing"]` (bed_conversion, mgmt, cnv, target, fusion; `workflow_ray.py:379-386`). It does so only when its result status is `success` and the BAM wasn't skipped (`:2752-2767`), filtered by the plan (`:2663-2664`, `:2749`).
   - **Non-Ray triggers:** the same five, through `job_dependencies` and `can_trigger_jobs` (`workflow_simple.py:185-245`).
   - **Steps the ledger ignores:** classifiers, IGV BAM, SNP and finalize are per sample, triggered or requested per sample. They are not in the ledger.
   - **Taking in (Ray):**
     - The scan (`submit_existing_paths`, `:3833-3946`) builds jobs with the reference and panel.
     - The watcher (`RayFileWatcher._handle`, `:4267-4286`) builds them without the reference (JOB-8). It marks a path done on its first event, before the file is complete.
     - `run()` scans (`:5056-5065`) before it starts the observer (`:5067-5084`).
     - `add_watch_path` schedules the watch first, then scans in another thread (`:4611-4652`).
     - The watcher and the scan keep separate in-memory sets.
   - **Taking in (non-Ray):** `FileWatcher.handle_file` (`workflow_simple.py:1549-1579`) has one set for the scan and the watcher. It scans before it starts the observer (`:1600-1606`). Its jobs get `work_dir` from the CLI's builder (`cli.py:2026-2035`).
   - **Failures:** the Coordinator records failed files as `_failed_jobs_log` rows, `file_failed` events and `file_failed_by_type` (`workflow_ray.py:2497-2537`). The non-Ray runner has no such channel.
   - **Completeness:** a BAM written by htslib/pysam ends with the 28-byte BGZF end-of-file block, written when the file is closed. A copy in progress doesn't end with it. This holds for the synthetic BAMs in the tests; for MinKNOW's own BAMs it is not measured here.
   - **Not supported:** two ROBIN sessions on one work folder at once. The intake's in-memory set is per process.

## Gate B: measured on `98dc6ef` before this lock

The tests below: **17 failed, 2 passed** (`before_10d2.txt`, in this session's scratchpad folder).
- 10 failed on their own assertion:
  - L1, L2, L7 and L13 failed with the BAM taken in again.
  - L4 failed: `only_steps` was None.
  - L8 failed: the half-written BAM was taken in.
  - L9 failed: no reference.
  - L12 and L18 failed: every step started.
  - L14 failed: the scan runs first.
- 7 failed because the code doesn't exist yet: no `robin.file_ledger` (L5, L6, L10, L11, L15), no `RayFileWatcher.intake` (L3), no `record_cut_short_steps` (L16).
- The 2 controls passed (C1, C2).

## Fixes

1. **`robin/file_ledger.py`:**
   - `bam_is_complete(path)`.
   - `FileLedger`:
     - `load()` returns per-BAM state. A torn last line (a stop mid-append) is ignored; any other unreadable line raises `LedgerError`.
     - `append(records)` appends under `flock` and fsyncs.
   - `BamIntake(work_dir, plan)`:
     - `decide(path)` returns whether to run the BAM and which steps (`only_steps`), plus any cut-short rows. It writes the "take" and the cut-short "end" records.
     - It takes nothing incomplete.
     - It keeps one in-memory set per session, under a thread lock.
   - `record_step_start(job, job_type)` and `record_step_end(job, job_type)`, for the six per-BAM steps and every BAM of a batch.
     - A start that can't be recorded raises: the step doesn't run.
     - An end that can't be recorded is logged. The step then shows as cut short after a restart.
2. **Ray:**
   - **The worker wrapper (`_impl`):** it records the start before the handler call and the end in `finally`, with the result's status. For preprocessing, the end also says whether its steps follow.
   - **`_on_finish`:** a preprocessing with `only_steps` starts only those.
   - **`Coordinator.record_cut_short_steps(rows)`:** a failed row (`failure_kind` "cut_short"), a `file_failed` event and the count.
   - **The scan and the watcher** build jobs with one helper (reference, panel, `work_dir`, fail-only, `bam_key`, `only_steps`), through the intake.
     - The watcher marks a path done only once it is decided. It fetches the reference from the Coordinator once (JOB-8).
   - **`run()`** makes one intake, starts the observer, then scans with that intake. `add_watch_path` uses the same intake.
3. **Non-Ray:**
   - **The worker** records the start and end around the handler call.
   - **`can_trigger_jobs`** honours `only_steps` for preprocessing.
   - **`FileWatcher.handle_file`** goes through an intake for the job's `work_dir` and plan, read from the job it builds. Cut-short steps are logged (there is no failure channel).
   - **`start()`** starts the observer before it scans.

## Rulings

- **Steps per sample** (classifiers, IGV BAM, SNP, finalize) are not in the ledger. A stop can leave a sample's last run of one undone until its next data or a manual run.
- **With `--no-process-existing`** nothing is scanned, so earlier BAMs' unstarted steps run only if the watcher sees those files again.
- **A BAM in the ledger that this session's folders don't hold** is not reported.
- **No work folder:** without `--work-dir` there is no ledger, only the in-memory set (as before), plus the completeness check.
- **A BAM that never completes** is never taken in. It is logged once as waiting. Showing it in the GUI is part 3 (JOB-4).
- **Restarted steps:** a restarted preprocessing writes its own start and end again (the latest one counts).

## Tests (new file `tests/test_bam_ledger.py`)

| # | Test | Before (`98dc6ef`) |
|---|------|------|
| L1 | a BAM scanned again after a restart is not taken in again | fail: taken again |
| C1 | different BAMs in one scan are each taken in, with the reference (control) | pass |
| L2 | a copy of a finished BAM is not taken in | fail: taken |
| L3 | a BAM seen by the scan and then the watcher is taken in once | fail: no shared intake |
| L4 | after a restart only the steps a BAM never started run (`only_steps`) | fail: None |
| L5 | a step cut short by a stop is not run again, is reported once, and not again the session after | fail: no module |
| L6 | a BAM whose preprocessing was cut short runs again, with every step | fail: no module |
| L7 | a BAM its preprocessing skipped runs nothing more | fail: taken again |
| L8 | a BAM still being written is taken in once it is complete | fail: taken half-written |
| L9 | the watcher gives its BAMs the reference (JOB-8) | fail: none |
| L10 | the worker records each BAM of a batch, and a failed step as failed | fail: no module |
| L11 | a step whose start can't be recorded doesn't run | fail: no module |
| L12 | a preprocessing told its steps starts only those (Coordinator) | fail: all started |
| L16 | the Coordinator reports a cut-short step as a failed file | fail: no method |
| L14 | the Ray run starts watching before it scans (order of the calls in `run()`) | fail: scan first |
| L15 | a ledger that can't be read stops the intake; a torn last line is ignored | fail: no module |
| L13 | the non-Ray watcher doesn't take a finished BAM in again | fail: taken again |
| C2 | a non-Ray preprocessing with no list starts every step of the plan (control) | pass |
| L18 | a non-Ray preprocessing told its steps starts only those | fail: all started |

## Pass rules (locked)

- **R1:** on `98dc6ef`, the 17 fail as shown, and C1 and C2 pass (measured).
- **R2:** after the fixes, all 19 pass, and so do `test_master_csv_once.py` and the batch-7 failure tests (`test_job_failures.py`, `test_file_failures.py`).
- **R3:** the full suite passes: 824 + 19 = 843.
- **R4:** each mutant below, applied alone to the fixed code, makes the file fail (pytest exit 1). The control N0 leaves it passing.
  - M1: the intake doesn't check the ledger (takes every BAM) → L1.
  - M2: the intake keys BAMs by path → L2.
  - M3: no completeness check → L8.
  - M4: a cut-short step is run again instead of reported → L5.
  - M5: `_on_finish` ignores `only_steps` → L12.
  - M6: the worker records no end → L1.
  - M7: a failed start record is logged and the step runs → L11.
  - M8: `run()` scans before it watches → L14.
  - M9: the watcher builds jobs without the reference → L9.
  - M10: the non-Ray `can_trigger_jobs` ignores `only_steps` → L18.
  - N0: a no-op change → nothing fails.
- **R5:** `robin.yml` and `pyproject.toml` unchanged.
