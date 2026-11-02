# Batch 10e part 2: a coverage commit adds a batch to every file, or to none (TGT-3)

**Locked (committed) before code, 2026-11-02, after its plan review and a short check of the revision.** Base `01d9edc` (part 1 done). Branch `fix/batch-10e-recorded-items`. Done by the controller (no agents). It adds saved state (an install list, and a folder of prepared files), so by the user's rule its plan is reviewed before code. The user chose (2026-11-02, "1"): one review for part 1's code and this plan, then a short check of both revised plans.

Batch 8 recorded, under TGT-3: "a commit is not a transaction: one that fails part-way leaves its files at different generations, and its batch folder stranded." Batch 6 recorded the same for a crash between a batch's coverage commit and its BAM's move.

**What the plan review changed** (`batch10e-report-codex.md`, findings 3 and 4):
- **The crash test's checker** missed five of the files a commit writes, and the panel depth; Codex built a mixed state it passed. It now checks every file. A test shows it catches a stale copy of each one.
- **New tests:** a second stop during the recovery; an append cut off part-way; an older folder's union file a generation behind.
- **The syncs** are now spelled out, in order, with what they guard against. What is measured is a process dying; a power loss can't be made here, so the order of the syncs is tested instead.

**What the short check changed** (`batch10e-r2-report-codex.md`, findings 2 and 3; the user said "go"):
- **Names on disk before the list:** syncing a file's contents does not put its name on disk. So `_commit/`, the claim folder and `_staging/` are now synced before the list is written. After a power loss, a list on disk could otherwise point at a batch BAM whose name was lost.
- **A batch BAM gone from its claim folder** is no longer taken as installed just because it is gone. The finish looks for it in the sample folder and in `target.bam`'s list of merged batch BAMs. If it is in neither, its BAMs are reported (C8).
- **The sync test (C7) was too loose:** it passed a wrong order, and syncs of the wrong files. It now checks each sync's place and the exact files moved. A second case keeps the claim folder.
- **The test's threshold files never changed:** at the default threshold every target is over it after one BAM. The test file now uses 1.5×, so the threshold files change at each commit, and C0 checks real stale copies (no made-up contents).
- **After the list is deleted, `_commit/` is synced** (it was in the plan, now with its reason and a test). A list a power loss brought back would be done again after a later commit, and its append would cut off that commit's rows.

## Gate A

1. **Principle:** a batch is added to all of a sample's coverage files, and its batch BAM handed to finalize, together or not at all. After a process stops, once the next accumulation or finalize has run, every file the commit writes and `target.bam` describe the same BAMs, and a BAM left out is listed. After a power loss the same holds, provided the disk keeps what was synced; that part is designed in (the syncs below), not measured.
2. **What would kill it:** any of these.
   - After a crash at any step of a commit (each `os.replace`, and the claim folder's removal), the next accumulation or finalize leaves any of the commit's files, or `target.bam`, describing different BAMs. The files are: the contig and target totals, the mirrored target file, the threshold count and BED, the covered stretches, the time chart, the latest read counts, the time series, and the union file through the panel's summary depth.
   - The same after a second crash, during the recovery.
   - A BAM counted and merged whose claim folder is still listed (an SNP run then refuses for good).
   - A commit that raised changes any file.
   - A time-series block missing, or counted twice, after a recovery, including one cut off part-way.
   - A sync the recovery needs that comes after the list is deleted, or a list written before what it moves, and those files' names, are synced.
   - A batch BAM gone from its claim folder counted as installed when it is neither in the sample folder nor in `target.bam`; or one merged into `target.bam` reported as lost.
3. **Refuted before here?**
   - Batch 8's `rows_bases` detects a union file a generation behind its rows. It stays, as a second guard for folders written before this change, with its own test (R1 legacy).
   - Batches 6 and 10b recorded the stranded claim folder. Part 1 (`3ab033c`) reports the BAMs of a failed commit as "not added"; for a commit that failed part-way that text was wrong, because some files had them.
   - Cases 15 and 16: mark something done only after it is. Here the list is deleted only after its steps, and after the syncs.
   - This plan's first draft claimed recovery "after any stop" and checked part of the files: its review found both.
4. **Simplest thing:** prepare every new file in a folder beside the real ones, write a list of the renames, then do them. A later accumulation or finalize finishes a list it finds. Considered and not chosen:
   - **A generation folder with a pointer file:** every reader of these files (the GUI, the report, the SNP and variant steps) would have to follow the pointer.
   - **Back up and roll back:** copies every file at each commit, and must put the claimed sets back in staging. More I/O and more steps.
   - **Rewriting the time series whole at each commit:** it grows with the run (targets × batches rows). An append that first cuts the file back to a recorded size is safe to repeat, and costs only the new rows.
   - **Ordering the writes so that one last rename publishes the commit:** there are ten files and a BAM, and readers read each file by its own name, so no single rename covers them.
5. **Inputs** (read in code at `01d9edc`):
   - **One commit at a time per sample:** the accumulation lock is `fcntl.flock` (`target_analysis.py:301-314`); the kernel lets it go when a process dies.
   - **Files a commit writes** (`_commit_batch`, `:1880-2050`): `_covered_stretches.parquet`, `coverage_main.csv`, `bed_coverage_main.csv`, `target_union_coverage.csv` (written, or removed), `target_coverage.csv`, `coverage_time_chart.npy`, `_target_coverage_latest_reads.parquet`, `target_coverage_time.csv` (appended in place, `mode="a"`), `targets_exceeding_threshold_count.txt`, `targets_exceeding_threshold.bed`. Each but the time series is written by `atomic_open`: a temp file, its contents synced, renamed. The folder is not synced.
   - **The batch BAM and its index** are written by `pysam` into the claim folder as `batch.bam` and `batch.bam.bai` (`_write_batch_bam`, `:1859-1878`), and not synced.
   - **After the commit, still under the lock** (`accumulate_staged_files`): the batch BAM's index, then the BAM, move to `<sample>/batch_<n>_<id>.bam` (`_install_batch_bam`, part 1). Then the claim folder is deleted, or kept with `incomplete.txt` when a set was unreadable or the batch BAM failed.
   - **Other writers:** `process_file`, the old non-staging path, writes the same files. Its only caller, `process_single_file`, has no caller in `src/` (dead code, ARCH-4).
   - **Readers** don't lock (unchanged), and none searches a sample folder recursively for these names (checked: every `rglob` and `os.walk` in `src/`).
   - **Finalize** takes the accumulation lock already, to list staged sets missing a part (`:3101-3102`), after it waits for claims (`_claims_left`, `:3095`).
   - **Finalize merges under its own lock**, not the accumulation lock (`:3113-3137`). So it can merge a batch BAM, and delete it, after a commit moved it into place and before that commit's list is deleted. `target.bam` is replaced by rename, and its header lists every batch BAM merged into it (`_merged_batch_names`, `:2786`).
   - **The call threshold** is 0.1× by default (`:842`, from the config). With the tests' reads every target is over it after one BAM, so the threshold files are the same after one, two and three BAMs (measured, as the plan check found). At 1.5×: target B is over it after BAM 0; A and B after BAMs 0 and 1, or 0 and 2; A, B and C after all three (measured with samtools).

## Gate B: measured on `01d9edc` before this plan

`tests/test_coverage_commit.py` (new) and `tests/test_panel_depth.py`'s ruled and new tests, on the code as it is (`before_10e_locked_all.txt`, in this session's scratchpad folder, with part 1 round 2's tests: 26 failed, 68 passed; the draft's numbers were in `before_10e_r2_all.txt`):

- **C0, the checker catches a stale copy of each of the 9 files:** 9 pass. It checks the test, so it passes before and after.
- **C1, a commit that raises at one file (9 cases):** 8 fail, each because earlier files were changed. The case that raises at the first file passes.
- **C2, a crash at each step:** a commit has 13 crash points. 11 broke, with the next accumulation as recovery and with finalize alone. Most: no set of BAMs fits both the files and `target.bam`. With finalize alone, the threshold files are now among the mismatches.
- **C3–C8:** each fails (C7 and C8 in both cases), mostly because there is no install list.
- **The suite collects 931 tests** (`collect_10e_locked.txt`).
- **R1, ruled:** fails. Right after the failed commit the shown depth jumps from 0.57× to 1.34×: it counts a BAM that is not in `target.bam`, and that part 1 reports as "not added".
- **R1 guard, ruled:** fails: the failed first commit wrote `bed_coverage_main.csv`.
- **R1 legacy (new guard):** passes.

## The change

All under the accumulation lock, as now.

1. **Finish an earlier commit first** (`_finish_commit(sample_dir)`), at the start of each locked part of an accumulation, and in finalize before it looks for claim folders (whether or not anything is staged):
   - `<sample>/_commit/install.json` exists: do its steps (each is safe to do twice), then the syncs of step 5, then delete the list, then `_commit/`.
   - `_commit/` exists without a list: delete it. Nothing was installed.
   - A list that can't be read raises, naming it. Nothing is guessed (see the rulings).
   - A finish's failed-file entries join the result of the accumulation or finalize that ran it.
2. **Prepare:**
   - Make `<sample>/_commit/`, then sync the sample folder, so the new folder's name is on disk.
   - `_commit_batch` keeps its maths, and writes each new file into `_commit/` under its own name, its contents synced. `_commit/` is inside the sample folder, so each rename stays on one filesystem.
   - For the time series it writes only the new rows (with a header when the file is new), and records the real file's size before.
   - Sync `_commit/`, so the prepared files' names are on disk.
   - When there is a batch BAM: sync `batch.bam` and `batch.bam.bai`, then the claim folder (their names), then `_staging/` (the claim folder's name, and the sets moved out of `_staging/` into it).
3. **Write the install list,** `_commit/install.json`, whole (temp file, contents synced, renamed), then sync `_commit/`, so the list's name is on disk. Its layout (paths relative to the sample folder): `{"v": 1, "steps": [...]}`, the steps in order:
   - `{"do": "replace", "file": F}`: `_commit/F` over `F`, for each prepared file. Done once `_commit/F` is gone.
   - `{"do": "remove", "file": "target_union_coverage.csv"}`: when this batch makes the merged total unknown.
   - `{"do": "append", "file": "target_coverage_time.csv", "size": N}`: cut the file back to N bytes, then add `_commit/target_coverage_time.csv`'s rows.
   - `{"do": "batch_bam", "from": "_staging/_batch_<id>/batch.bam", "to": "batch_<n>_<id>.bam", "bams": [{"file": path, "key": K}, ...]}`: move the index, then the BAM. `bams` are the BAMs whose reads it holds, with their keys (part 1 round 2's key files), for a finish to report.
   - `{"do": "drop_claim", "dir": "_staging/_batch_<id>"}`: delete the claim folder. Listed only when every set was counted and the batch BAM was written (as now).
4. **Install:** do the steps in order. In `batch_bam`, the index and the BAM each move if still in the claim folder. When the BAM is gone from the claim folder:
   - it is in place when `to` exists, or `target.bam`'s header lists `to` as merged (a finalize merged it and deleted it, under its own lock);
   - otherwise it is lost: handled as a `batch_bam` step that fails (below), with that reason;
   - a `target.bam` that can't be read stops the finish with the error; the list stays.
5. **Sync, then retire the list:**
   - Sync the appended time series, the sample folder (the renames, the removal, the moved BAM), the claim folder when it was kept, and `_staging/` after the claim folder's removal.
   - Then delete the list, sync `_commit/` (so a power loss can't bring the list back), delete `_commit/`, and sync the sample folder.
   - A sync that fails raises before the list is deleted, so the list stays and the next finish does it all again.

**What a stop means now:**
- **Before the list is written** (an error, or a crash): `_commit/` is deleted, now or by the next finish. Nothing changed. The claim folder stays, and part 1 reports its BAMs as "not added" (now true).
- **After the list is written:** the commit counts as done (part 1's `committed`), so no "not added" entries. If the install raised, or the process died, the next accumulation or finalize finishes it, and a stop during that finish is finished by the one after.
- **A `batch_bam` step that fails** (part 1's "could not be moved into place", or lost): its index goes back when it moved, the claim folder is kept with its note, `drop_claim` is skipped, and the list's `bams` are reported as not in `target.bam`: failed-file entries, and ledger lines (part 1 round 2). The other steps go on. A finish does the same, and logs it.

## Rulings

- **R1 and its guard change** (`test_panel_depth.py`; batch 8's tests of a commit that fails on the union file). They checked that, after a commit wrote the rows and then failed, the row sum showed instead of a partial merged total. A failed commit now changes nothing. So:
  - R1 (`test_a_commit_that_fails_part_way_changes_nothing`) checks that the shown numbers right after the failed commit are those before it, and that after two more commits the merged total is samtools' over BAMs 1, 3 and 4. The failed BAM's set is kept and listed, not counted.
  - The guard (`test_a_first_commit_that_fails_part_way_writes_nothing`) checks that a failed first commit writes nothing, and that the next commit's merged total is known and right.
  - Both are stricter than before: exact numbers, where they checked for a fallback.
  - The fallback itself keeps a test, R1 legacy (`test_an_older_folders_union_file_a_generation_behind_falls_back_to_the_row_sum`). It builds such a folder by hand, as one written before part 2 could be.
- **An install list that can't be read** stops the sample's accumulations, with an error naming it, as an unreadable totals file does now. It is written whole and synced, so this needs damage from outside. Guessing could count a batch twice or not at all.
- **A step that fails every time** (a rename that keeps failing) makes each later accumulation of that sample fail at its finish, with the error. Skipping it would leave the files mixed, which is the failure this change removes.
- **A prepared file gone from `_commit/`** was renamed into place. A rename leaves a file under one of its two names, and `_commit/` is synced before the list, so nothing else removes it. (The real file can't be told apart from the prepared one by its name; the batch BAM can, so it is looked for.)
- **A lost batch BAM is reported, not a stop.** Its BAMs' coverage is counted, and their reads are not in `target.bam`: the same as a batch BAM that could not be moved, which part 1 reports and goes on from. Stopping every later accumulation would not bring the reads back.
- **What is measured:** a process dying, at every rename and the claim folder's removal (C2), and a second time during the recovery (C5); an append cut off part-way (C6); a batch BAM gone from both places, merged or not (C8). **A power loss is not measured.** The syncs are in the order steps 2, 3 and 5 give, which C7 checks.
- **Readers are unchanged:** during an install's renames (milliseconds), a reader can still see two generations, as today. After a crash, the files stay mixed until the next accumulation or finalize (each SNP run and the end of a run do one). Finishing at start-up was not chosen: nothing at start-up knows the sample folders before their first BAM.
- **`process_file`** (dead code) is left as it is, with ARCH-4.

## Tests

`tests/test_coverage_commit.py`:

| # | Test | Before (`01d9edc`) |
|---|------|------|
| C0 | the checker catches a stale copy of each of the 9 files; each differs from the current one (9 cases) | pass (it checks the test) |
| C1 | a commit that raises at any one of its 9 files changes no file; with the next BAM, every file and `target.bam` agree, and the failed BAM is listed (9 cases) | 8 fail, 1 passes |
| C2 | a commit killed at any of its steps leaves every file and `target.bam` agreeing, once the next accumulation, or finalize, has run (2 cases; every step tried in each) | fail: 11 of 13 broke, in each |
| C3 | an error during the install leaves the list; the BAM is not reported "not added", and the next accumulation finishes the commit | fail |
| C4 | an install list that can't be read stops the accumulation with an error naming it, and changes no file | fail |
| C5 | a commit killed during its install, then the finish killed at each of its install steps: the next finish completes it, with the BAM in | fail: no list |
| C6 | a time-series append cut off part-way is put right by the finish | fail: no list |
| C7 | each sync in its place (steps 2, 3 and 5), for the exact files the install moves; the claim folder dropped, or kept (2 cases) | fail: no list |
| C8 | a finish whose batch BAM is gone from the claim folder and the sample folder: merged into `target.bam`, nothing reported and every file agrees; or lost, its BAM reported (failed file, ledger line) and its claim folder kept (2 cases) | fail: no list |

`tests/test_panel_depth.py`: R1 and its guard, ruled above (both fail before); R1 legacy (new guard, passes before).

The test file sets the call threshold to 1.5× (through the config), so the threshold files change at each commit.

The kill in C2, C5 and C8 raises a `BaseException` from `os.replace` or the claim folder's `shutil.rmtree`: nothing catches it, and the locks are let go, as when a process dies. The checker compares every file the commit writes, the time series (one block per counted batch, its row count, the last block's reads), the panel's summary depth, and `target.bam` with samtools on one set of BAMs. It checks that finalize lists the left-out BAM's claim folder (or its incomplete staged set) only when that BAM was left out.

C7 records each sync, rename, removal, folder made and folder deleted, in order. It checks:
- after `_commit/` is made and before the list: the sample folder;
- after the last prepared file and before the list: `_commit/`;
- before the list: the batch BAM and index the install later moves (those exact files), then the claim folder, and `_staging/`;
- after the list and before the first install step: `_commit/`;
- after the install's last step (the claim folder's removal, or the last rename when it is kept) and before the list is deleted: the time series, the sample folder, and `_staging/` (dropped) or the claim folder (kept);
- after the list is deleted: `_commit/`; after `_commit/` is deleted: the sample folder.

## Pass rules (locked after the check)

- **S1:** on `01d9edc`, the tests fail or pass as measured above.
- **S2:** after the change, `test_coverage_commit.py`, `test_panel_depth.py`, `test_target_failed_files.py`, `test_target_finalize.py`, `test_snp_every_read.py`, `test_target_coverage.py`, `test_target_vs_samtools.py` and `test_bam_ledger.py` pass.
- **S3:** the full suite passes: 894, plus part 1 round 2's 8, plus this part's 29 (C0–C8: 28; R1 legacy: 1; R1 and its guard replace two) = 931 (counted with `--collect-only` before locking).
- **S4:** each mutant below, alone, makes its test fail (pytest exit 1); the no-op control N0 leaves them passing.
  - Q1: files written straight to the sample folder (no prepare) → C1.
  - Q2: no finish at the start of an accumulation → C2 (next accumulation).
  - Q3: finalize doesn't finish → C2 (finalize).
  - Q4: the list deleted before its steps → C2.
  - Q5: the time rows appended without cutting back first → C6.
  - Q6: the commit not counted as done once its list is written → C3.
  - Q7: the claim folder deleted after the list, not as a step in it → C2.
  - Q8: a failing `batch_bam` step stops the install → `test_one_bam_path_a_batch_bam_that_cannot_be_moved_into_place_is_a_failed_file` (part 1).
  - Q9: a list that can't be read is skipped → C4.
  - Q10: a finish skips the threshold files' steps (a commit's own install still does them) → C2 (finalize).
  - Q11: the install skips the stretches' step → C2.
  - Q12: the install skips the union file's step → C2.
  - Q13: a finish deletes `_commit/` without doing the list's steps → C5.
  - Q14: the list deleted before the syncs of step 5 → C7.
  - Q15: the batch BAM not synced before the list is written → C7.
  - Q16: the sample folder synced before `_commit/` is made, not after → C7.
  - Q17: `_staging/` synced before the claim folder is removed, not after → C7 (dropped).
  - Q18: the claim folder and `_staging/` not synced before the list is written → C7.
  - Q19: `_commit/` not synced between its prepared files and the list → C7.
  - Q20: `_commit/` not synced after the list is deleted → C7.
  - Q21: a finish takes a batch BAM gone from its claim folder as installed, without looking → C8 (lost).
  - Q22: a finish takes it as lost without looking at `target.bam`'s header → C8 (merged).
  - Q23: a kept claim folder not synced after the install → C7 (kept).
- **S5:** `robin.yml` and `pyproject.toml` unchanged.
