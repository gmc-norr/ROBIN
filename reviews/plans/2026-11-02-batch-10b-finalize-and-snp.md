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
