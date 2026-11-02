# Batch 10b: one finalize per request, and SNP calls from every read (JOB-7, JOB-6)

**Locked before code, 2026-11-02.** Base `79c8870`. Branch `fix/batch-10b-finalize`. Done by the controller (no agents). Second of step 8's five batches (10a–10e). An external review comes before the merge.

## Rule (agreed by the user, option "A", 2026-11-02)

- **One request, one merge:** when the finalize job is queued, the GUI does not also run the merge itself.
- **No session-long skip:** a sample finalized once in a session still gets its later reads merged.
- **SNP calls from every read:** an SNP run first merges what is not yet in `target.bam`.
- **Earlier SNP calls are reused only when they came from those same reads;** otherwise the variants are called again.

**Spec:** JOB-6 (its second and third bullets) and JOB-7 in `reviews/2026-10-16-code-review.md`, and the user's choice.

## Gate A

1. **Principle:** an SNP run calls on every read the sample has when it starts, and shows earlier calls only if they came from those reads; one finalize request merges once.
2. **What would kill it:** any of these.
   - A queued finalize that the GUI also runs itself.
   - A batch BAM or a staged file that never reaches `target.bam` because the sample was finalized before in the session.
   - An SNP run that reads a `target.bam` missing reads the sample already has.
   - Earlier calls reused after `target.bam` gained reads, or while reads wait to be merged.
3. **Refuted before here?** Case 10 (batch 10a): don't exempt the default path from a check. Here the default is `force_regenerate=False` (all three SNP buttons send it), so the reuse rule is tested on that path. Learning `2026-11-01-a-batch-s-reads-reach-target-pa`: a batch BAM appears only after its coverage commit, and finalize merges under its own lock; the IGV copy follows `target.bam` by the `@CO` batch list, not by file time. The SNP record below uses the same batch list, for the same reason.
4. **Simplest thing:** one condition in the GUI (merge itself only when the job was not queued); delete the in-memory "finalized" set (the on-disk check already does its job); one call to the existing finalize at the start of the SNP job; a small text file next to the SNP calls with the batch names they came from.
5. **Inputs (read in code at `79c8870`; the test results below measure each):**
   - **JOB-7:** `gui_launcher.py:7522` sets `finalization_succeeded = True` after a queued finalize job, and `:7554` then runs `finalize_accumulation_for_sample` in the GUI as well (without the reference).
   - **JOB-6, second bullet:** `_finalized_samples` is read at `:7454` (`already_finalized`: skip the merge) and `:953` (no finalize on a later Live → Complete). After a sample's first finalize, later reads stay out of `target.bam` for the rest of the session.
   - **`_is_target_bam_finalize_redundant`** (`:719`) calls a sample done when `target.bam` and its index exist and no `batch_*.bam` is left. It ignores staged files (`_staging/coverage_*.parquet`), which finalize would merge.
   - **The SNP job** (`snp_analysis_handler`, `target_analysis.py:4190`) never merges: it reads `target.bam` as it is.
   - **All three SNP entry points** (the finalize+SNP button `gui_launcher.py:7707`, the coverage panel `coverage.py:4626`, bulk SNP `:7103`) send `force_regenerate=False`, and `run_snp_analysis` (`:3191`) then reuses earlier calls whenever the two `snpsift` VCFs exist. Bulk SNP and its sample list skip any sample with those files (`_sample_has_snp_calling_outputs`, `:6724`); the "SNP complete" poller (`:6931`) uses the same check.
   - **Batch BAM names are unique** (`batch_<n>_<8 hex>.bam`, `:2168`), and `_merged_batch_names` reads the merged ones from `target.bam`'s `@CO` lines.
   - **JOB-6, third bullet** (unmerged batch BAMs deleted when `target.bam` exists) is already fixed by `24ebbc8` (batch 6, TGT-1): the old cleanup (`eccbf18`, `target_analysis.py:2602-2618`) is gone, and finalize merges every new batch BAM into the existing `target.bam`. Only the review's Status line is out of date.

## Gate B: measured on `79c8870` before this lock

**JOB-7's damage, two finalizes at once** (`scratchpad/b10b/gateb_two_finalizes.py`, in the batch-7 session folder). Two separate processes finalize one sample at the same moment; each trial checks `target.bam` against samtools on all reads, and that no batch BAM or staged file is left. 20 trials with no earlier finalize and 20 after one.
- First run: lock on 0/40 bad, but the control (the `finalize` lock switched off) was also 0/40 bad. The fixture's merge takes milliseconds, so the two never overlapped: the check could not see damage.
- With a random 0.1–0.6 s pause before each merge: **lock on 0 bad of 40; lock off 39 bad of 40.** Without the lock, one finalize fails ("samtools merge: fail to open batch_….bam": the other deleted it), JOB-7's symptom.
- So batch 6's lock and unique temp names already stop the damage. What is left of JOB-7 is the double run (wasted work, and the GUI's copy runs without the reference).

**The tests below, run on `79c8870`:** 13 failed, 5 passed (`scratchpad/before_10b.txt`, in this session's folder). Each failure is the one the test was written for (messages in the table).

## Fixes

1. **G1 (JOB-7):** the GUI runs the merge itself only when there is no runner, the runner has no finalize submit, or the submit fails.
2. **G2 (JOB-6):** `_finalized_samples` is deleted, with its seeding from disk: its two reads (`:7454`, `:953`) are the skips above, and the on-disk check already stops a finalize that has nothing to do. Every Live → Complete starts a finalize.
3. **G3 (JOB-6):** `_is_target_bam_finalize_redundant` also needs no staged file.
4. **M1 (JOB-6):** the SNP job first calls `finalize_accumulation_for_sample` on its sample (its panel, defaulting to `rCNS2` as the finalize job does, and its reference). If that returns an error, the job fails with it and no variants are called.
5. **M2 (JOB-6):** a finished SNP call writes `clair3/target_bam_batches.txt`: the batch names in `target.bam` right after the job's merge, one per line, written under a temporary name and renamed. Not after an annotation-only run (it calls nothing).
   - `snp_calls_are_current(sample_dir)`: both `snpsift` VCFs exist, the record exists and equals `target.bam`'s batch list, and no batch BAM or staged file waits.
   - `run_snp_analysis` reuses earlier calls only when that is true (and `force_regenerate` is off). When it does not reuse them and it is not annotation-only, it deletes the record before calling, so a run that stops part-way never leaves a record that matches.
   - The GUI's `_sample_has_snp_calling_outputs` returns `snp_calls_are_current`. So bulk SNP and its list include a sample whose calls are older than its reads, and the poller waits for this run's calls.
   - A folder from before this batch has no record, so its next SNP run calls again, once.

## Tests (new files `tests/test_finalize_once.py`, `tests/test_snp_reads_current.py`)

| # | Test | Before (`79c8870`) | Why it failed there |
|---|------|------|------|
| F1 | a queued finalize is not also run by the GUI | fail | GUI merged as well |
| F2 | a queued finalize then SNP sends one finalize and one SNP job | fail | GUI merged as well |
| F3 | the GUI merges itself when the job can't be queued (no runner; queue refused) | pass ×2 | control |
| F4 | a second finalize in one session merges new batches | fail | second finalize skipped |
| F5 | files still staged are merged after an earlier finalize | fail | called "done" with a staged file |
| F6 | a sample that completes twice is finalized twice | fail | second completion started none |
| S1 | an SNP run first merges reads not yet in target.bam | fail | read the old target.bam |
| S2 | an SNP run fails when its merge fails | fail | called variants anyway |
| S3 | a finished call records the batches it read | fail | no record |
| S4 | an annotation-only run records nothing | pass | guard |
| S5 | calls from the current target.bam are reused | pass | control |
| S6 | calls from an older target.bam are not reused | fail | reused |
| S7 | calls without a record are not reused | fail | reused |
| S8 | calls are not reused while reads wait (batch; staged) | fail ×2 | reused |
| S9 | bulk SNP lists a sample whose calls are older than target.bam | fail | left out |
| S10 | bulk SNP leaves out a sample whose calls are current | pass | control |

The GUI tests build a real `GUILauncher` without NiceGUI, and swap `threading.Thread` for one that runs its target at once. The SNP tests never run Clair3: the handler tests replace `run_snp_analysis` with one that notes `target.bam`'s reads, and the reuse tests remove the targets BED, so a run that does not reuse stops at its input checks and returns "".

## Pass rules (locked)

- **P1:** on `79c8870`, the 13 above fail for the reason shown and the 5 pass (measured).
- **P2:** after the fixes, all 18 pass.
- **P3:** the full suite passes: 757 + 18 = 775. No existing test is changed.
- **P4:** each mutant below, applied alone to the fixed code, makes the two new files fail (pytest exit 1); the control N0 leaves them passing (exit 0).
  - N1: G1 undone (merge in the GUI after a queued job) → F1, F2.
  - N2: the "already finalized" skip restored → F4.
  - N3: the scan starts a finalize only on a sample's first completion → F6.
  - N4: G3 undone (staged files ignored) → F5.
  - N5: the SNP job does not merge → S1.
  - N6: a failed merge is ignored → S2.
  - N7: reuse by the VCFs alone → S6, S7.
  - N8: reuse ignores reads waiting to merge → S8.
  - N9: the record is kept when the calls may be redone → S6.
  - N10: no record is written → S3.
  - N11: a record is written after annotation-only → S4.
  - N12: the GUI check uses the VCFs alone → S9.
  - N0: a no-op change → nothing fails.
- **P5:** `robin.yml` and `pyproject.toml` unchanged.

## Not done here (recorded for later)

- **A Complete without a finalize:** when the workflow's job types are known and their outputs exist, the background scan's first Complete check (`:926-943`) marks a timed-out sample Complete without starting a finalize, so the later check never sees a Live → Complete. Its batch BAMs then wait for the next SNP run or a manual finalize. Which of the GUI's scans sees the change first in a real run is not measured.
- **The poller while reads keep arriving:** with new reads after an SNP run, the "SNP complete" poller waits until its 24 h limit and then says "still running", though the run has ended (its calls are no longer current).

## Result (2026-11-02, on `c93a711`): supported

Measured in this session's scratchpad folder.

- **P1, before:** on `79c8870` src with the new tests: 13 failed, 5 passed (`before_10b.txt`). Each failure has the message the test was written for: "the GUI ran the merge itself" (F1, F2), "the second finalize skipped the merge" (F4), the "done" check true with a staged file (F5), "the second completion started no finalize" (F6), the old `target.bam`'s reads (S1), "variants were called on a target.bam missing reads" (S2), no record (S3), the calls reused (S6, S7, S8 ×2), the sample left out of bulk SNP (S9). F3 ×2, S4, S5 and S10 passed.
- **P2, after:** 18 passed (`after_10b.txt`).
- **P3, suite:** 775 passed (`suite_10b.txt`): 757 + 18. No existing test was changed. The target, completion-timeout and GUI-error tests ran on their own as well: 99 passed.
- **P4, mutants:** N1–N12 each caught (pytest exit 1); N0 survived, 18 passed (`mutants10b.txt`, runner `mutate10b.py`). The runner stops at the first failing test, so each line names one failure, not every test a mutant breaks.
- **P5, pins:** `git diff 79c8870 -- robin.yml pyproject.toml` is empty.

Not changed, noticed while fixing: the GUI's own merge (now only when the job can't be queued) still runs without the reference, as before. The coverage panel shows SNP calls without saying when they are older than `target.bam`; the poller and bulk SNP now know, the panel doesn't.

## Round 2: fixes from the external review (locked before code, 2026-11-02)

Codex reviewed `c93a711` (`/home/parlar_ai/robin-batch3-run/review-requests/batch10b-report-codex.md`) and found four gaps. Each was confirmed by reading the code at `54fe364`; Codex's probe scripts were read, not run. The user chose to fix all four on this branch ("1", 2026-11-02).

### Gate A

1. **Principle (unchanged):** an SNP run calls on every read the sample has when it starts, and shows calls as current only if they came from exactly those reads. When it can't, it fails and says why.
2. **What would kill it:** any of these.
   - A finalize or SNP run that goes ahead while a worker still holds claimed reads (R1).
   - A claim folder deleted while its reads are not in `target.bam` (R1).
   - An SNP run that goes ahead after an accumulation that timed out, could not read its files or could not write its batch BAM (R2).
   - A Clair3 run counted complete when a region failed, a VCF is missing, or a VCF is from an earlier run (R3).
   - Calls counted current after `target.bam` changed (R4).
3. **Refuted before here?** Round 1 counted only top-level staged files and batch BAMs as "waiting", and only `"error"` as a failed merge. Both were inferred, not checked against `accumulate_staged_files`.
4. **Simplest thing:**
   - R1: a lock per claim, held by its worker from claim to end. Finalize waits on the claims it finds. A claim folder still there once its lock is free is kept, and counts as reads not in `target.bam`.
   - R2: a list of the merge results that allow calling.
   - R3: a start time and a check per region. No deletions: the Clair3 container runs as root, so its region folders may not be ours to delete.
   - R4: `target.bam`'s size and modification time join the batch list in the record.
5. **Inputs (read in code at `54fe364`):**
   - **Claims:** accumulation moves claimed sets into `_staging/_batch_<uuid>/` under the lock (`target_analysis.py:2079-2089`) and works on them outside it.
     - The folder is deleted at the end (`:2172`), even when the batch BAM failed (`:2140-2145`) or a set could not be read (`:2106-2108`).
     - With no set read, it is deleted and `load_failed` is returned (`:2124-2127`).
     - A commit that raises leaves the folder.
   - **Waiting checks:** `_get_pending_count` (`:986`) and `_reads_wait_to_merge` (`:2731`) look only at top-level staged files. So does the GUI's `_is_target_bam_finalize_redundant` (`gui_launcher.py:731`).
   - **Results:**
     - Accumulation returns `no_files`, `below_threshold`, `load_failed`, `success` (with `batch_bam_error` when the batch BAM failed), `lock_timeout` or `error`.
     - Finalize adds `no_pending_files`. It returns `error` for its own exceptions, including its merge lock timeout.
     - The SNP job stops only on `error` (`:4342`).
   - **Clair3 region loop** (`:3742-3871`):
     - A non-zero exit is logged and skipped (`continue`).
     - A missing per-region VCF only warns.
     - `merge_vcf_files` returning False is ignored (`:3879`, `:3888`).
     - Single-pass copies `snv.vcf.gz` whenever it exists (`:3897`).
     - Files from earlier runs stay in `clair3/`.
   - **Container user:** the container runs as root (no `user=` in `create_host_config`, `:3532`).
   - **Not measured:** whether ClairS-TO always writes `indel.vcf.gz` when it succeeds.
     - A run with no INDEL file already fails today, at annotation (`gzip.open` on the missing file, `:4108`, then `return ""`). So requiring the file adds no new failure for a first run.
     - If ClairS-TO skips the file for a region with no candidates, split runs will now fail, with a message naming the region.
   - **Other writers:** `TargetAnalysis.process_file` (`:1254-1265`) rewrites `target.bam` without updating its `@CO` lines. The workflow does not call it (`target_handler` → `process_multiple_files` → `process_file_with_staging`), but it shows that `target.bam` can change without a new batch.

### Gate B: measured on `54fe364` before this lock

The 19 tests below: **17 failed, 2 passed** (`before_r2.txt`, in this session's scratchpad folder). Each failure has the message it was written for. The 2 controls passed.

### Fixes

1. **R1 (finding 1, and the lost reads in finding 2):**
   - `accumulate_staged_files` takes a claim lock (`_locks/claim<folder name>.lock`) inside the accumulation lock, before it makes the claim folder. It releases and removes the lock when it ends (`finally`).
   - The claim folder is kept, not deleted, when its batch BAM failed or a set could not be read. A file `incomplete.txt` in it says which sets' coverage was counted, and why the folder was kept.
   - Finalize, after its own accumulation, waits on each claim folder's lock, up to `_CLAIM_WAIT_SECONDS` (300 s) in all.
     - A folder still there after that is listed in its result as `claims_left`.
     - Its reads are neither merged nor deleted.
   - `_reads_wait_to_merge` also counts claim folders. The GUI's `_is_target_bam_finalize_redundant` now uses it.
2. **R2 (finding 2):** the SNP job goes ahead only when finalize returns `success` without `batch_bam_error`, `no_files` or `no_pending_files`, and no `claims_left`. Otherwise it fails with what went wrong, and calls nothing.
3. **R3 (finding 3):** a run that calls notes its start time.
   - Each region must exit 0, and must write `snv.vcf.gz` and `indel.vcf.gz` during this run: modified at or after the start, less 2 s for file-system time rounding.
   - A VCF merge that fails fails the run.
   - Otherwise the run returns "" (failed), so no record is written.
4. **R4 (finding 4):** the record, renamed `clair3/target_bam_source.json`, holds `target.bam`'s batch list, size and modification time (ns), taken after the job's merge. Calls are current only when all three match.

### Rulings

- **Round 1's tests** are on this branch and not merged. `test_snp_reads_current.py`'s helpers now write and read the record in its new form; the assertions keep their meaning. No test on main is changed.
- **Not done here (recorded for later):**
  - Recovering a kept claim folder. It needs knowing whether its coverage was counted; `incomplete.txt` records that.
  - The finalize job still fails only on `"error"`.
  - The coverage panel still shows calls from before a failed run without saying so.

### Tests (new file `tests/test_snp_every_read.py`)

| # | Test | Before (`54fe364`) |
|---|------|------|
| W1 | finalize waits for a worker still adding reads | fail: merged before the worker added its reads |
| W2 | calls are not current while a worker adds reads | fail: current |
| W3 | the GUI does not call finalize redundant while a worker adds reads | fail: redundant |
| W4 | an SNP run fails while a worker still adds reads after the wait | fail: variants called |
| W5 | finalize reports a claim folder a stopped worker left | fail: not reported |
| W6 | an SNP run fails on a claim folder a stopped worker left | fail: variants called |
| W7 | a finished claim leaves no folder and no lock | pass: control |
| A1 | a batch BAM that cannot be written keeps its staged files | fail: deleted |
| A2 | a staged file that cannot be read keeps its folder | fail: deleted |
| A3 | an SNP run fails again after a batch BAM could not be written | fail: variants called |
| A4 | an SNP run fails when its merge missed reads (lock timeout; load failed; batch BAM error) | fail ×3: variants called |
| C1 | a complete run keeps every region's calls | pass: control |
| C2 | a failed region fails the run | fail: shown complete |
| C3 | old INDEL calls are not kept when the caller writes none | fail: old file used |
| C4 | an old SNV file is not taken as new after a failed single pass | fail: old file used |
| C5 | a failed VCF merge fails the run | fail: old merged file used |
| I1 | calls are not reused after target.bam changes without a new batch | fail: current |

The C tests run the real `run_snp_analysis`, with a stand-in for the `docker` module (each container writes the VCFs a region would) and for snpEff and SnpSift (they pass the calls through unchanged).

### Pass rules (locked)

- **Q1:** on `54fe364`, the 17 fail for the reason shown and the 2 controls pass (measured).
- **Q2:** after the fixes, the 19 pass, and round 1's 18 pass.
- **Q3:** the full suite passes: 775 + 19 = 794.
- **Q4:** each mutant below, applied alone to the fixed code, makes the new file fail (pytest exit 1). The control N0 leaves it passing (exit 0).
  - K1: finalize does not wait for claims → W1.
  - K2: `_reads_wait_to_merge` ignores claim folders → W2, W3.
  - K3: the SNP job ignores `claims_left` → W6.
  - K4: the SNP job stops only on `"error"` → A4.
  - K5: the claim folder is deleted after a batch BAM error → A1, A3.
  - K6: the claim folder is deleted when a set could not be read → A2.
  - K7: no claim lock is taken → W1.
  - K8: the claim lock file is not removed → W7.
  - K9: a failed region is skipped again → C2.
  - K10: per-region VCFs are not checked against the start time → C3.
  - K11: a failed VCF merge is ignored → C5.
  - K12: the record has no size and time → I1.
  - N0: a no-op change → nothing fails.
- **Q5:** `robin.yml` and `pyproject.toml` unchanged.

### Round 2 result (2026-11-02, on `e139331`): supported

Measured in this session's scratchpad folder.

- **Q1, before:** on `54fe364` with the 19 new tests: 17 failed, 2 passed (`before_r2.txt`). Each failure had the message it was written for. W7 and C1 passed.
- **Q2, after:** the 19 new tests and round 1's 18 tests: 37 passed (`after_r2.txt`).
- **Q3, suite:** 794 passed (`suite_r2.txt`): 775 + 19.
- **Q4, mutants:** K1–K12 each caught (pytest exit 1). N0 survived, with 19 passed (`mutants_r2.txt`, runner `mutate_r2.py`). The runner stops at the first failing test:
  - K3 failed first at W4. W4 also needs `claims_left`, and it runs before W6.
  - Every other mutant failed first at the test the plan named.
- **Q5, pins:** `git diff 54fe364 -- robin.yml pyproject.toml` is empty.

Not changed, noticed while fixing: the coverage panel's own SNP fallback (`coverage.py:4691`) runs through `snp_analysis_handler`, so it gets the same merge, checks and record.

## Round 3: fixes from the second external review (locked before code, 2026-11-02)

Codex reviewed `e139331` (`/home/parlar_ai/robin-batch3-run/review-requests/batch10b-r2-report-codex.md`) and found two gaps, both rated P2. Both were confirmed by reading the code at `75c9e7f`. The user chose to fix both, then merge and push ("1", 2026-11-02). Codex found no fault in the claim locks or the kept folders.

### Gate A

1. **Principle (unchanged):** an SNP run calls on every read the sample has when it starts, and shows calls as current only if they came from exactly those reads. When it can't, it fails and says why.
2. **What would kill it:**
   - A run that takes a VCF this run did not write (D).
   - An SNP run that goes ahead while a staged set with a missing part sits in `_staging/` (E).
3. **Refuted before here?**
   - The 2 s allowance in round 2 was for "file-system time rounding" that was never measured. It opened the gap Codex found: a file from a run that ended less than 2 s earlier passes.
   - `no_files` means "no complete staged set", not "no staged reads" (`_complete_staged_counters`, `target_analysis.py:806-813`). That is the same mistake as case 12: an upstream value taken for what it seemed to mean.
4. **Simplest thing:**
   - D: each calling run writes Clair3's output into its own new folder `clair3/run_<8 hex>/`, so a file there can only be this run's, and the time check goes. Earlier `run_*` folders are deleted when possible (the container writes as root, so some files may stay).
   - E: finalize lists the staged sets whose coverage file exists but a part is missing (`incomplete_staged`), and the SNP job's check fails on them.
5. **Inputs (read in code at `75c9e7f`):**
   - **Staging order:** staging writes the coverage file last, by rename (`:1086-1101`), so a coverage file with a missing part is broken, not still being written.
   - **Moves:** a claim moves a set's files one at a time, under the accumulation lock (`:2082-2089`). Listing broken sets under the same lock never sees a set half moved.
   - **Readers:** nothing outside `run_snp_analysis` reads Clair3's raw `snv.vcf.gz`, `indel.vcf.gz` or `region_*` files (grep of `src/`).
   - **Output folders:** in split mode the container makes each `region_<i>` folder inside the output folder it is given. In single-pass mode it writes into a folder that exists. The new run folder is made by ROBIN, so both patterns stay the same.

### Gate B: measured on `75c9e7f` before this lock

The 5 tests below: **5 failed** (`before_r3.txt`, in this session's scratchpad folder), each with the message it was written for. The 19 round-2 tests passed.
- D1 depends on timing at `75c9e7f`: it fails there because the second run starts less than 2 s after the first. After the fix it does not depend on timing.

### Fixes

1. **D:** `run_snp_analysis` makes `clair3/run_<8 hex>/` before its region loop, after deleting earlier `run_*` folders (errors ignored).
   - Region outputs go to `run_<id>/region_<i>/` (split) or `run_<id>/` (single pass).
   - The merged VCFs go to `run_<id>/`.
   - Each region must exit 0 and leave both VCFs there. The time check is removed.
   - `output_done.vcf.gz` and `output_indel_done.vcf.gz` are still copied from this run's files, as before.
2. **E:**
   - `_incomplete_staged_sets(staging_dir)`: the coverage files whose set is not complete, read under the accumulation lock.
   - Finalize puts them in its result as `incomplete_staged`.
   - `_merge_missed_reads` fails on them, naming the files.
   - The files are left in place.
3. **Wording:** the claim that any rewrite changes `target.bam`'s identity is too strong. Every way ROBIN writes `target.bam` changes it; a replacement from outside that keeps the old size and modification time would not be noticed (Codex's note). The docstring says so. Round 2's learning, which says otherwise, is not edited; a new note corrects it.

### Tests (added to `tests/test_snp_every_read.py`)

| # | Test | Before (`75c9e7f`) |
|---|------|------|
| D1 | a quick retry does not take the previous run's INDEL file | fail: the previous run's file used |
| D2 | old run folders are removed | fail: no run folder |
| E1 | finalize reports a staged set missing a part | fail: not reported |
| E2 | an SNP run fails on a staged set missing a part (alone; beside a complete set) | fail ×2: variants called |

### Pass rules (locked)

- **R1:** on `75c9e7f`, the 5 fail for the reason shown (measured).
- **R2:** after the fixes, the whole file passes (24), and round 1's 18 pass.
- **R3:** the full suite passes: 794 + 5 = 799.
- **R4:** each mutant below, applied alone to the fixed code, makes the file fail (pytest exit 1). The control N0 leaves it passing (exit 0).
  - L1: region outputs go back to `clair3/region_<i>/`, with no time check → D1.
  - L2: earlier run folders are not deleted → D2.
  - L3: finalize does not list broken staged sets → E1.
  - L4: the SNP job's check ignores `incomplete_staged` → E2.
  - N0: a no-op change → nothing fails.
- **R5:** `robin.yml` and `pyproject.toml` unchanged.
