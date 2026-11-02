# Batch 10d-1: a sample's master.csv counts each BAM once (JOB-2, preprocessing)

**Locked before code, 2026-11-02.** Base `01a5482`. Branch `fix/batch-10d-file-ledger`. Done by the controller (no agents for code). First of batch 10d's eight parts.

## Rule (agreed by the user, 2026-11-02)

The user chose "re-run only its unfinished steps" for a BAM that was part-way through when ROBIN stopped: exactly once, not "never run it again".

For that, every step that adds a BAM's data to a sample's totals keeps the list of BAMs it has added in the same file as the totals, written in one go. A step run again for a BAM on its list changes nothing. A ledger (part 7) then says which steps each BAM has not finished, and only those run after a restart.

**Order** (agreed, "go"):
1. Preprocessing (`master.csv`): this part.
2. CNV.
3. Methylation.
4. MGMT.
5. Fusion.
6. Target coverage.
7. The ledger, JOB-8.
8. The watcher (JOB-4).

After parts 1–6, a restart can't double anything; it only redoes work. After part 7 it redoes only what's unfinished.

**A BAM is known by its bytes:** its size and its first and last 64 KiB, hashed (`robin.bam_identity.bam_key`). So a moved folder, or a copy of a BAM, is the same BAM.

## Gate A

1. **Principle:** a sample's `master.csv` counts each BAM once, whatever runs it again: a restart, the scan and the watcher both seeing it, or a copy of it at another path. A BAM it can't count is reported.
2. **What would kill it:** any of these.
   - A BAM preprocessed twice, from the same path or a copy, and its counts added twice.
   - Two different BAMs, and one not counted.
   - An unreadable `master.csv` replaced by a fresh one.
   - A failed update that is only logged.
3. **Refuted before here?** No. Case 14 applies: every path into the counts must fail loudly. Here that's `_load_existing_data`, which today starts from zero when the file can't be read.
4. **Simplest thing:**
   - A `bams_added` column in `master.csv`'s one row: a JSON list of BAM keys. It is checked and extended under the existing lock, and written in the same atomic write as the counts.
   - The key lives in a new small module, `robin/bam_identity.py`.
   - The handler reports a failed update as a file failure (the batch 7 channel), so the job and its other steps go on.
5. **Inputs (read in code at `01a5482`):**
   - **The file:** `master.csv` is one row, written whole (temp file + `os.replace`, `master_csv_manager.py:34-54`), under `<sample>/_locks/master_csv.lock` (10 s, `:64-66`).
   - **The writers:** every one goes through `MasterCSVManager`, and loading keeps every column the file has (`:104-106`), so a new column survives the others.
     - `update_master_csv`: `bam_preprocessor.py:908`.
     - `update_analysis_panel`: `bam_preprocessor.py:915`, `target_analysis.py:2579`, `:2683`.
     - `update_sample_overview`: `gui_launcher.py:963`.
   - **The readers** pick named columns:
     - `csv.DictReader` at `workflow_ray.py:2040`, `:2128`, `workflow_simple.py:2053`, `:2127`, `gui/components/mnpflex.py:47`, `target_analysis.py:3100`, `:4465`, `gui_launcher.py:853`, `:7426`.
     - pandas for the GUI summary (`gui/components/summary.py:99-104`) and the report (`reporting/sections/run_data.py:116-176`).
     - So an extra column changes nothing they show.
   - **Unreadable file:** an unreadable `master.csv` prints a warning and returns the defaults (`:118-119`). The update then writes this BAM alone over it. The other two writers do the same with their own fields.
   - **Errors:** `update_master_csv` logs and swallows its own errors (`:80-82`). The handler's own catch (`bam_preprocessor.py:919-920`) never sees them.
   - **One BAM per job:** preprocessing runs one BAM per job (`workflow_ray.py:243`, `max_batch_size: 1`; `:1837`, submitted directly). The BAM is `job.context.filepath` (`bam_preprocessor.py:714`).
   - **Pass or fail:** decided by "pass" in the path (`bam_preprocessor.py:236`).
   - **File failures:** the Ray runner records `file_failures` as failed rows and events (`workflow_ray.py:2497-2537`). The non-Ray runner has no such channel (no `file_failures` in `workflow_simple.py`).
   - **The key's assumptions:**
     - **Not verified here:** a BAM is not changed after it is complete.
     - **Argued, not measured on real BAMs:** two different BAMs differ in size or in their first or last 64 KiB. The last 64 KiB holds the last reads and the end-of-file block, and read names differ between BAMs. The test below checks two synthetic BAMs with the same header and different reads.

## Gate B: measured on `01a5482` before this lock

The tests below: **6 failed, 2 passed** (`before_10d1.txt`, in this session's scratchpad folder).
- T1, T2 and T4 failed with `counter_bam_passed` 2 instead of 1.
- T6 failed: the zero-byte `master.csv` was replaced (1,172 bytes).
- T7 failed: no file failure was reported.
- T8 failed: there is no `robin.bam_identity` yet.
- The controls T3 and T5 passed.

## Fixes

1. **`robin/bam_identity.py`:** `bam_key(path)` returns blake2b-128 of a version tag, the size, the first 64 KiB and the last 64 KiB. A BAM of 128 KiB or less is hashed whole.
2. **`MasterCSVManager.update_master_csv(sample_id, bam_data, bam_info, bam_key)`:**
   - Under the lock it loads the row. A key already in `bams_added` changes nothing (logged). Otherwise it adds the counts, appends the key and writes everything in one write.
   - Any failure raises.
3. **`_load_existing_data`:** a `master.csv` that exists but doesn't give exactly one row, or whose `bams_added` isn't a JSON list, raises. It no longer returns the defaults.
   - The two other writers already catch and log, so they now leave such a file as it is instead of resetting it.
4. **The handler:** it computes the key and calls the update.
   - A failure is added to the job's `file_failures` (`{"file": <bam>, "error": …, "step": "master_csv"}`), and the panel update is skipped.
   - The job and its other steps go on.

## Rulings

- **Before this change:** a `master.csv` written earlier has no `bams_added`. Its counts are kept and the next BAMs are added (T5). But the BAMs it already counted are not on any list. So a restart before part 7 can count each of them once more. The CHANGELOG says so: for exact numbers, use a new work folder.
- **The other writers:** `update_analysis_panel` and `update_sample_overview` still log and go on when they fail. Unchanged, except that they no longer reset an unreadable file.
- **The non-Ray runner** has no file-failure channel, so there the failure is only logged (recorded for later). The count is still never added twice, and the file is still never reset.
- **Not done here:** the other steps still add a BAM again (parts 2–6).

## Tests (new file `tests/test_master_csv_once.py`)

All go through the real `bam_preprocessing_handler`, with synthetic BAMs from `target_fixtures`.

| # | Test | Before (`01a5482`) |
|---|------|------|
| T1 | a BAM preprocessed again is counted once | fail: counted twice |
| T2 | a copy of a BAM at another path is counted once | fail: counted twice |
| T3 | different BAMs are each counted (control: A alone + B alone = A then B) | pass |
| T4 | other writers (overview, panel) keep the record of counted BAMs | fail: counted twice |
| T5 | a `master.csv` from before this change (no record) keeps its counts and is added to (control) | pass |
| T6 | a `master.csv` that can't be read (zero bytes) is kept, and the failure is reported for that BAM | fail: replaced |
| T7 | a failed update (the write raises) is reported for that BAM, with its error | fail: not reported |
| T8 | a BAM is known by its bytes: a copy has the same key, and the same header with other reads has another | fail: no module |

## Pass rules (locked)

- **P1:** on `01a5482`, T1, T2, T4, T6, T7 and T8 fail for the reasons shown, and T3 and T5 pass (measured).
- **P2:** after the fixes, all 8 pass.
- **P3:** the full suite passes: 816 + 8 = 824.
- **P4:** each mutant below, applied alone to the fixed code, makes the new file fail (pytest exit 1). The control N0 leaves it passing.
  - M1: the record isn't checked (a counted BAM is added again) → T1.
  - M2: the key is the path, not the bytes → T2.
  - M3: an unreadable `master.csv` starts from zero again → T6.
  - M4: `update_master_csv` swallows its errors again → T7.
  - M5: the handler logs a failed update instead of reporting it → T7.
  - N0: a no-op change → nothing fails.
- **P5:** `robin.yml` and `pyproject.toml` unchanged.

## Result (2026-11-02, on `52b5169`): supported

Every pass rule holds.

- **P1, before** (`01a5482`, measured before the lock, `before_10d1.txt`): T1, T2, T4, T6, T7 and T8 failed, each for its own reason, and T3 and T5 passed.
- **P2, after:** all 8 pass.
- **P3, full suite:** 824 passed (`suite_10d1.txt`).
- **P4, mutants** (`mutate_10d1.py`; only pytest exit 1 counts as caught; `mutants_10d1.txt`):
  - M1: caught by T1.
  - M2: caught by T2.
  - M3: caught by T6.
  - M4: caught by T7.
  - M5: caught first by T6 (the runner stops at the first failure). Run against T7 alone, it is caught by T7 too.
  - N0 survived (8 passed).
- **P5, pins:** `git diff 01a5482 -- robin.yml pyproject.toml` is empty.

The output files are in this session's scratchpad folder.

## The plan changed after this part (the user's choice, 2026-11-02)

The 8-part order above is replaced, and this part stays.
- **Why:** mapping parts 2–6 showed six counting steps and seven batches, two of them large (target coverage commits about ten files one by one, and fusion appends part files with new ids on each run). The user asked for a simpler way and chose it.
- **The new plan:** a ledger records when each of a BAM's steps starts and ends. After a restart, finished steps are skipped and unstarted ones run. A step that was running at the stop is not run again (it may have added its counts) and is listed as possibly incomplete.
- **The remaining parts:**
  - Part 2: the ledger, with JOB-8.
  - Part 3: JOB-4.
- **Why this part stays:** re-running a BAM's unstarted steps starts from its preprocessing, which must not count the BAM in `master.csv` again.
