# Batch 10c: a sample's CNV read counts are shared by every worker (CNV-2)

**Locked before code, 2026-11-02.** Base `0ed84c4`. Branch `fix/batch-10c-cnv-counts`. Done by the controller (no agents). Third of step 8's five batches (10a–10e).

## Rule

- **No worker keeps its own copy of a sample's counts.** Each CNV run reads the sample's counts from disk.
- **Two CNV runs for one sample never overlap.** One holds a per-sample lock from reading the counts to saving them; the other waits.
- **Counts that can't be read fail the run and are kept.** They are never reset to none (the rule since batch 7: a failure says "failed").
- **Counts are written whole or not at all.**

**Spec:** CNV-2 in `reviews/2026-10-16-code-review.md` ("always reload from disk under a per-sample lock, and write with tmp + `os.replace`"), and step 8's line "CNV-2: no worker caches, and atomic writes".

## Gate A

1. **Principle:** a sample's CNV counts hold every BAM that a run reported as added, exactly once, whichever worker added it. A run that can't add its BAMs says so.
2. **What would kill it:** any of these.
   - A batch that another worker added is missing from the counts after a later run.
   - Two runs for one sample at once, and the saved counts hold only one run's BAMs.
   - Unreadable counts and a run that reports success.
   - A write that fails part-way, and the last good counts are gone.
3. **Refuted before here?**
   - Case 4's lesson (batch 4, CNV-1): the library adds a BAM's reads into the dict it is given. So a cached dict is changed in place, and it must never be shared.
   - `SERIALIZE_BY_TYPE_PER_SAMPLE` (`workflow_ray.py:410`) is empty, and its code paths (`:1916`, `:2965`) have never run. Turning it on would serialise CNV in the scheduler. It is not used here: it is untested, and it would not cover a CNV run outside Ray.
4. **Simplest thing:**
   - Delete the per-sample cache.
   - Two small helpers read and write the counts.
   - A thin wrapper takes the lock around each of the two CNV entry points, so their bodies are not re-indented.
5. **Inputs (read in code at `0ed84c4`):**
   - **The cache:** `_sample_cache` (`cnv_analysis.py:116`) is a module-level dict, one per Ray worker process.
     - `process_single_bam` (`:1195`) and `process_multiple_bams` (`:1510`) use its counts whenever its path matches, without reading the disk.
     - `load_analysis_counter` (`:813`) uses its counter in the same way.
   - **Reading:** a load error sets the counts to `{}` (`:1229`, `:1541`). The run then goes on and saves them, so every earlier BAM is lost.
   - **Writing:** the counts are written with `open(path, "wb")` (`:1301`, `:1645`, and the legacy migration at `:1214`, `:1530`).
   - **Concurrency:**
     - The CNV pool runs `analysis_workers` jobs at once (`workflow_ray.py:1453`) as Ray tasks, and nothing serialises a sample's CNV jobs.
     - `cnv_handler` always calls `process_multiple_bams` (`:1842`). `process_single_bam` has no caller in `src/`.
   - **The Coordinator's cleanup:** `_cleanup_cnv_cache` (`workflow_ray.py:2214`) runs in the Coordinator's process, whose cache is not the workers'. So it clears nothing.
     - The import block that brings in the cache functions (`:160-171`) also brings in `cnv_handler`. Removing a name there without editing the block would leave CNV with no handler.
   - **The result:** the handler reads its result only with `.get` (`:1852-1867`), and files `error_message` as a failure (`:1865-1868`).
   - **The lock:** `master_bed_generator.FileLock` (fcntl `flock`, polling, closes its file on a timeout) is already imported by `cnv_analysis` (`:1076`). It works across processes and across threads.

## Gate B: measured on `0ed84c4` before this lock

The tests below: **8 failed, 1 passed** (`before_10c.txt`, in this session's scratchpad folder), each failure with the message it was written for. The import guard passed.

The process tests run stand-in Ray workers in their own processes (`multiprocessing` with "spawn"), each keeping its module state between jobs, with a synthetic control.

## Fixes

1. **No cache.** Deleted:
   - In `cnv_analysis.py`: `_sample_cache`, `_current_sample_id`, `get_sample_cache`, `update_sample_cache`, `clear_sample_cache`, `set_current_sample`, `cleanup_sample_cache_on_completion` and `get_sample_cache_stats`.
   - In `workflow_ray.py`: the Coordinator's `_cleanup_cnv_cache`, its scheduling and counters, and `get_cnv_cache_stats`. The import block keeps `cnv_handler` only.
   - The counter is read from disk each time.
2. **Lock.** `process_single_bam` and `process_multiple_bams` each take `<sample>/_locks/cnv.lock` and run their old bodies, unchanged, under it.
   - The wait is up to `_CNV_LOCK_WAIT_SECONDS` (3 h). That number is not measured; it only has to exceed a CNV run, whose own timeout starts at 1 h.
   - On a timeout, the run returns an `error_message`, and its BAMs are not added.
3. **Reading.** `_load_copy_numbers` reads the counts, or migrates the legacy file as before.
   - No file means no counts yet (`{}`).
   - An unreadable file raises, so the run fails with the path in its message and the file is left as it is.
4. **Writing.** `_save_copy_numbers` writes through `atomic_open`, and so does the legacy migration.

## Rulings

- `tests/test_cnv_control_reference.py` calls `clear_sample_cache()` three times (its `cached_control` fixture, and once mid-test). Those calls are removed, because there is no cache to clear. No assertion changes.
- **Not done here (recorded for later):**
  - `cnv_analysis_counter.txt` is still written in place.
  - In subprocess mode (off by default), a batch still keeps only its last BAM (CNV-3's note).
  - `SERIALIZE_BY_TYPE_PER_SAMPLE` stays empty and untested.

## Tests (new file `tests/test_cnv_counts_shared.py`)

| # | Test | Before (`0ed84c4`) |
|---|------|------|
| T1 | a worker keeps the batch another worker added (single; multiple) | fail ×2: the batch lost |
| T2 | two runs for one sample at once keep both batches (single; multiple) | fail ×2: one overwrote the other |
| T3 | counts that cannot be read fail the run and are kept (single; multiple) | fail ×2: the run went on |
| T4 | a failed counts write keeps the last good counts (single; multiple) | fail ×2: the counts unreadable |
| T5 | the workflow still imports the CNV handler | pass: guard |

- T1 is CNV-2's example. Worker A runs batch 1, worker B batch 2, worker A batch 3. It also checks that the counter is 3.
- T2 pauses each worker 1.5 s after its analysis and before it saves, so the two overlap.

## Pass rules (locked)

- **P1:** on `0ed84c4`, the 8 fail for the reason shown, and T5 passes (measured).
- **P2:** after the fixes, all 9 pass.
- **P3:** the full suite passes: 799 + 9 = 808.
- **P4:** each mutant below, applied alone to the fixed code, makes the new file fail (pytest exit 1). The control N0 leaves it passing (exit 0).
  - M1: a per-process cache of the counts returns (`_load_copy_numbers` reuses a module-level copy) → T1.
  - M2: `process_multiple_bams` runs without the lock → T2[multiple].
  - M3: `process_single_bam` runs without the lock → T2[single].
  - M4: unreadable counts read as none → T3.
  - M5: the counts are written in place → T4.
  - N0: a no-op change → nothing fails.
- **P5:** `robin.yml` and `pyproject.toml` unchanged.

## Result (2026-11-02): supported

Every pass rule holds.

- **P1, before** (`0ed84c4`, measured before the lock, `before_10c.txt`): 8 failed, 1 passed. Each failure gave the message it was written for:
  - T1: "a batch another worker added was lost".
  - T2: "one run overwrote the other's counts".
  - T3: "the run went on with no counts".
  - T4: the counts file was left truncated (`UnpicklingError: pickle data was truncated`).
  - T5, the import guard, passed.
- **P2, after** (`1f2c5b1`): the 9 new tests pass, and so do the 20 in `test_cnv_control_reference.py` and `test_cnv_reference.py` (29 passed, `after_10c.txt`). Re-run on `1f2c5b1` while writing this result: all 7 CNV test files, 49 passed.
- **P3, full suite:** 808 passed (799 before, plus the 9), `suite_10c.txt`. It ran on the `src/` and `tests/` committed as `1f2c5b1`; neither has changed since.
- **P4, mutants** (`mutate_10c.py`, each alone on the fixed code; only pytest exit 1 counts as caught; `mutants_10c.txt`):
  - M1, a per-process cache returns: caught by T1[single].
  - M2, `process_multiple_bams` without the lock: caught by T2[multiple].
  - M3, `process_single_bam` without the lock: caught by T2[single].
  - M4, unreadable counts read as none: caught by T3[single].
  - M5, the counts written in place: caught by T4[single].
  - N0, the no-op control, survived (9 passed).
- **P5, pins:** `git diff 0ed84c4 -- robin.yml pyproject.toml` is empty.

The output files are in this session's scratchpad folder, as in Gate B.
