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

## Result (2026-11-02, on `4c60efb`): supported

Every pass rule holds.

- **R1, before** (`98dc6ef`, measured before the lock, `before_10d2.txt`): 17 failed as shown, and C1 and C2 passed.
- **R2, after:** all 19 pass. With `test_master_csv_once.py`, `test_job_failures.py` and `test_file_failures.py`: 236 passed (`after_10d2.txt`).
- **R3, full suite:** 843 passed (`suite_10d2.txt`).
- **R4, mutants** (`mutate_10d2.py`; only pytest exit 1 counts as caught; `mutants_10d2.txt`). Each was caught by the test the plan names:
  - M1 by L1, M2 by L2, M3 by L8, M4 by L5, M5 by L12.
  - M6 by L1, M7 by L11, M8 by L14, M9 by L9, M10 by L18.
  - N0 survived (19 passed).
- **R5, pins:** `git diff 98dc6ef -- robin.yml pyproject.toml` is empty.

The output files are in this session's scratchpad folder.

### Ruling made after the lock: the batch-7 stand-in BAMs

Four tests in `tests/test_job_failures.py` failed on the fixed code. Their stand-in BAMs were empty files, or paths under `/nonexistent`. A BAM is now taken in only once it is complete, and it is known by its bytes. So the intake took in none of the empty files, or only one, since they were all the same.
- **The change:** each stand-in now holds its own name, then the 28-byte BGZF end-of-file block (`_stand_in_bam`). The watch-session test writes its files under `tmp_path`.
- **What stays:** every assertion in those tests is unchanged. What they check (which BAMs go to fail-only submission, and that each failure is reported) is the same.

## Round 2: the external review of parts 1–2 (locked before code, 2026-11-02)

Codex reviewed `01a5482..4c60efb` (`batch10d-report-codex.md`, in the run notes folder) and found five problems. Claude found a sixth while writing the brief (B1). Each was checked by reading the code; Codex's probe scripts were not run. The tests below measure each one.

### The user's choices (2026-11-02)

- **Keep the automatic resume** ("1"). The user asked whether a warning alone would be better, and kept the resume: without it, a restart means re-analysing the run in a new work folder, or numbers counted twice.
- **A failed or cut-short step is reported at every start**, as long as its BAM is in the folders scanned. It is still never run again. Before, it was reported once (L5 changes with this; see the rulings).
- **The list of BAMs `master.csv` counts moves out of `master.csv`** into its own file ("1"). The other option was raising Python's CSV size limit in every ROBIN process, which leaves other readers of the file broken.

### Gate A

1. **Principle:** a BAM's data is added to its sample's totals at most once, and a step that failed or may be incomplete says so at every start.
2. **What would kill it:** any of these.
   - A step's batch list that holds another step's BAMs.
   - A per-BAM step run a second time after its worker died.
   - A `master.csv` that Python's `csv` module can't read at its default limit.
   - A failed or cut-short step with no failed-file row at a later start.
   - A fail-only BAM's expected failure reported as a failure.
   - A torn ledger line that makes the ledger unreadable a restart later.
3. **Refuted before here?**
   - The 10d-2 learning (`l-294c2aaa`) says a torn last line is safe to ignore. That holds for reading, but the next append was glued onto it (B1).
   - Case 14 (batch 10c): a second way into the same data. Here, Ray's own retry is a second way into a step.
   - Batch 7 decides what counts as a failure: the job's outcome (`job_outcome`). A fail-only BAM's expected failure is not one.
4. **Simplest thing:** for each finding, the fix listed below. None adds a dependency.
5. **Inputs** (read in code at `2905c04`):
   - **Batch lists:**
     - `_convert_batched_jobs` (Ray `workflow_ray.py:3206`, non-Ray `workflow_simple.py:1302`) uses the batch's first context as the job's context, and writes `_batched_job` into it.
     - A BAM's per-BAM jobs share one context: `_on_finish` builds them all with `ctx` (`:2802`).
     - With batching on, every one of them goes through the batcher (`submit_jobs`, `:1827`).
     - Handlers read `_batched_job` without checking its step (`mgmt_analysis.py:1477`, `cnv_analysis.py:1656`, `target_analysis.py:2471`), and so does the ledger's `_job_bams`.
     - A job is pickled for its worker when it is dispatched, which can be later than its conversion.
   - **Retries:**
     - The six per-BAM remotes (`workflow_ray.py:1022–1031`) set no `max_retries`, so Ray's default applies: 3 retries after the worker dies.
     - Both call sites (`TypeProcessor.process`, `:1088`; `:3735`) add only `RESOURCE_HINTS` (CPUs, memory).
     - No ROBIN code re-runs a failed job.
   - **`master.csv`:**
     - `bams_added` is a JSON list in one cell, 36 characters per key.
     - Python's `csv` field limit is 131,072, so the 3,641st key makes the row unreadable to `csv`.
     - The target finalize's panel lookup (`target_analysis.py:3095–3108`) reads it with `csv.DictReader`, swallows the error and uses `rCNS2`.
   - **Failures after a restart:**
     - `_resume` writes a cut-short step's end record when it decides. The scan reports the rows only after it has submitted its jobs (`workflow_ray.py:4016`).
     - An ended step makes no row, whatever its status.
     - The end record keeps the results status, not the job's outcome.
     - Files that failed inside a batch that went on are in `metadata["file_failures"]` (`{file, error, step}`). In the non-Ray runner, sibling jobs share a context and nothing removes that list.
   - **The ledger:** `FileLedger.append` opens the file to append, and writes after whatever is there.

### Gate B: measured on `2905c04` before this lock

The new tests and the ruled L5, on the unfixed code: **13 failed, 30 passed** (`before_10d2r2.txt`, in this session's scratchpad folder). Each failed for its own reason (table below). The 4 new controls and the ruled `master.csv` test passed.

### Fixes

1. **B1, a torn ledger line.** The ledger, and the `master.csv` list, append through one helper. Under the exclusive lock, it first cuts a torn tail back to the last newline (and logs it), then writes and fsyncs.
2. **Finding 1, batch lists.** Both `_convert_batched_jobs` give the job its own context: `dataclasses.replace(first, metadata=dict(first.metadata))`, with `_batched_job` set on that copy only. Its results, history and errors stay the first BAM's own (the same objects), so what a handler writes lands where it did.
3. **Finding 2, retries.** The six per-BAM remotes get `max_retries=0`.
   - A step whose worker dies then fails, and the Coordinator shows it as a failed job.
   - Its start has no end, so after a restart it is reported as cut short.
   - The per-sample steps keep Ray's default.
4. **Finding 3, `master.csv`.** It drops `bams_added` and gains `last_bam_added` (one key). The list is `<sample>/master_csv_bams.txt`, one key per line. Under the per-sample lock:
   - Read the list. A torn last line is ignored; any other line that isn't a key raises (the file is kept and reported).
   - If `last_bam_added` is not on the list, append it first: the last update stopped between its two writes.
   - A key on the list changes nothing.
   - Otherwise: add the counts, set `last_bam_added`, write `master.csv`, then append the key.
5. **Findings 4 and 5, failures after a restart.**
   - The worker's end record also holds:
     - the job's `outcome` (success, failed or expected_failure);
     - its `error` (shortened);
     - `file_failed`, for each BAM named in the files the job reported as failed. Only files added during this job count: a count is taken before the handler, in both runners.
   - At every start, `_resume` returns a row for each per-BAM step of the BAM, preprocessing too, that:
     - was cut short: it started with no end. Its end is then written as cut_short, as before.
     - ended cut_short in an earlier session;
     - ended with the outcome "failed";
     - ended with a file failure.
   - Rows have `failure_kind` "cut_short" or "failed_before", and the error says which, with the earlier text. An expected failure makes no row. No step is run again.
   - `Coordinator.record_cut_short_steps` takes each row's `failure_kind` (default "cut_short").
   - The non-Ray runner logs the rows, as before.

### Rulings

- **L5** (`tests/test_bam_ledger.py`): its last check changes from "not reported the session after" to "reported again", with the user's rule change. Its name changes with it (`..._is_reported_at_each_start`).
- **"A `master.csv` from before this change is added to"** (`tests/test_master_csv_once.py`): its setup also removes `last_bam_added` and the list file, since the record moved. Its check is unchanged.
- **End records from before this round** have no `outcome`, and are read as not failed (only cut_short shows). Such records exist only on this unmerged branch.
- **A failure to append the key after `master.csv` was written** raises, so the BAM's preprocessing shows a `master_csv` failed file. The next update repairs the list.

### Tests (new)

| # | Test | Before (`2905c04`) |
|---|------|------|
| Q1 | a torn last line is cut off before the next record | fail: `LedgerError` |
| Q2 | two batches that share a BAM keep their own lists (Ray; and the ledger records the CNV batch's BAMs) | fail: CNV list became the MGMT batch's |
| Q3 | the same, non-Ray | fail: same |
| C3 | a converted batch sees its batch, and its results reach the first BAM (control) | pass |
| Q4 | a per-BAM step is not run again when its worker dies (`max_retries` 0) | fail: None |
| C4 | a classifier keeps Ray's retries (control) | pass |
| Q5 | a step that failed is reported at each start and not run | fail: no row |
| Q6 | a preprocessing that failed is reported at the next start and not run | fail: no row |
| Q7 | a BAM that failed inside a batch that went on is reported at the next start (only that BAM) | fail: no row |
| C5 | a fail-only BAM's expected failure is not reported at the next start (control) | pass |
| C6 | a BAM whose steps all ended well reports nothing at the next start (control) | pass |
| Q8 | a cut-short step whose report was lost to a second stop is reported at the next start | fail: no row |
| Q9 | the Coordinator keeps the kind of an earlier failure | fail: "cut_short" |
| Q10 | a sample with 4,000 BAMs keeps a `master.csv` every reader can read (with its panel) | fail: `csv.Error`, field limit |
| Q11 | a stop between `master.csv` and its list counts the BAM once | fail: no list helper |
| Q12 | a list cut off mid-key is repaired | fail: no list file |
| L5 (ruled) | a step cut short is reported at each start | fail: not reported the session after |

### Pass rules (locked)

- **R6:** on `2905c04`, the tests fail as shown, and the controls pass (measured).
- **R7:** after the fixes, `test_bam_ledger.py` and `test_master_csv_once.py` all pass, and so do `test_job_failures.py` and `test_file_failures.py`.
- **R8:** the full suite passes: 843 + 16 = 859.
- **R9:** each mutant below, applied alone, makes its test fail (pytest exit 1). The control N0 leaves them passing.
  - X1: the append doesn't cut a torn tail → Q1.
  - X2: Ray's `_convert_batched_jobs` writes into the shared first context → Q2.
  - X3: the non-Ray one does → Q3.
  - X4: the CNV remote keeps Ray's retries → Q4.
  - X5: no row for a step whose outcome failed → Q5.
  - X6: a row for any outcome but success → C5.
  - X7: the worker ignores file failures → Q7.
  - X8: no row for an earlier cut_short → Q8.
  - X9: the Coordinator ignores the row's `failure_kind` → Q9.
  - X10: `master.csv` also keeps every key in a column → Q10.
  - X11: no repair when `last_bam_added` is not on the list → Q11.
  - X12: no row for a failed preprocessing → Q6.
  - X13: the converted batch's context is a deep copy → C3.
  - N0: a no-op change → nothing fails.
- **R10:** `robin.yml` and `pyproject.toml` unchanged.

### Round 2 result (2026-11-02, on `80d457f`): supported

Every pass rule holds. The output files are in this session's scratchpad folder.

- **R6, before** (`2905c04`, measured before the lock, `before_10d2r2.txt`): 13 failed as shown, and the 4 controls and the ruled `master.csv` test passed.
- **R7, after:** with `test_job_failures.py` and `test_file_failures.py`, 252 passed (`after_10d2r2.txt`).
  - The first run had 1 failure, in the code, not the test: `_append_counted_key` was a static method, and Q11 (written first) treats it as a method. Made it one; the test is unchanged.
- **R8, full suite:** 859 passed (`suite_10d2r2.txt`).
- **R9, mutants** (`mutate_10d2r2.py`; only pytest exit 1 counts as caught; `mutants_10d2r2_ledger.txt`, `mutants_10d2r2_master.txt`):
  - Each was caught by the test the plan names: X1 by Q1, X2 by Q2, X3 by Q3, X4 by Q4, X5 by Q5, X6 by C5, X7 by Q7, X9 by Q9, X10 by Q10, X11 by Q11, X12 by Q6, X13 by C3.
  - X8 was caught first by L5 (the runner stops at the first failure). Run against Q8 alone, it is caught by Q8 too.
  - N0 survived both files (32 and 11 passed).
- **R10, pins:** `git diff 2905c04 -- robin.yml pyproject.toml` is empty.
