# Batch 7: failures are visible (step 7) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** when a step fails or never runs, what the user sees says "failed" or "not available"; "none found" appears only when the step ran and found nothing.

**Architecture:** the handler wrapper decides each job's outcome from the errors filed during that job, whatever name they were filed under; both schedulers (`standard`/`p2i` pools and the `high` preset's processors) and the non-Ray runner read that one outcome. Failures and worker-side GUI messages travel back with the job, the Coordinator keeps the recent ones, and the GUI process sends them on. Classifiers report their errors instead of logging them. The ClinVar step writes its own status line into the VCF it produces; the report and GUI read it. The fusion report tells "not run", "unreadable" and "none found" apart.

**Outcome:** (filled in after the result)

**Tech Stack:** Python 3.12, Ray 2.59.0 (local, in tests), pandas, ReportLab, NiceGUI 3.18 (GUI code: no harness).

**Spec:** `reviews/2026-10-16-code-review.md` step 7 — ERR-1 (layers 1–5), METH-8, VAR-4, REP-6 — and the user's answers on 2026-11-01 (Gate A): the principle above; all four findings in one batch, job counting first, worker messages last.

**Spec notes (measured or read on `92c7c00`; probes in the session scratchpad `b7/`):**
- **Measured, real local Ray, a real `Coordinator` (`enable_batching=False`), one MGMT job whose synthetic handler files its error under `mgmt_analysis` as the real one does:** `standard` (pools) → completed 1, failed 0; `high` (processors) → completed 1, failed 0. **Control:** the same handler raising → completed 0, failed 1 on both. Probe pitfall, found by the control: Ray workers don't see the driver's `sys.argv`; a test's fault switch must travel in the job's metadata.
- Error names (read): mgmt, target, fusion, sturgeon, nanodx, pannanodx and random_forest file under `<type>_analysis`; CNV only logs a batch error (`cnv_analysis.py:1864-1867`); preprocessing, bed_conversion, igv_bam, snp_analysis and target_bam_finalize file under the job type. The pool judge matches the job type name only (`workflow_ray.py:3519-3526`); the `high` judge looks at exceptions only (`2905-2915`); the wrapper's `finally` adds `{"status": "success"}` unless an error matches the job type name (`803-811`). The non-Ray runner counts exceptions only (`workflow_simple.py:948-958`).
- A failed job does not trigger its next jobs (`_on_finish`). Among the 8 mismatched types only `target` has one (`igv_bam`, `TRIGGERS`); SNP calling and finalize are scheduled separately, and finalize rebuilds an existing IGV BAM (batch 6).
- Fail-only BAM submissions: the batch paths of sturgeon, nanodx, pannanodx, random_forest and bed_conversion skip filing an expected error; the single-file paths don't (`sturgeon_analysis.py:493`, `nanodx_analysis.py:771`, `random_forest_analysis.py:847`); the wrapper applies the same rule to exceptions (`workflow_ray.py:770-797`).
- GUI (read): `_update_errors` sets the three error-count labels from the payload (`gui_launcher.py:1794-1812`); nothing sends `ERROR_UPDATE`; the "Recent errors" label is created as "No errors detected." and never changed (`6587`). In a Ray worker `send_gui_update` finds no launcher and drops the update silently (`gui/app.py:64-67`), e.g. the "BAM Alignment Warning" (`bam_preprocessor.py:635-660`) and "Skipped (too many reads)" (`workflow_ray.py:755-767`).
- ERR-1 layer 1 (read): the modkit merge logs a failing file and goes on (`temp_utilities.py:360-366`); the parallel BAM step does the same (`bed_conversion.py:252-262`); the job reports success. `failed_jobs.csv` is written only at Coordinator shutdown (`workflow_ray.py:3296-3300`).
- METH-8 (read): Sturgeon's `sturgeon_bam_background_work` catches every exception (`585-587`); NanoDX's `_direct_classification` returns empty arrays and writes `errordf.csv` (the sample's methylation calls) into the worker's current directory (`180-184`), after which the store gets a row with `number_probes` 0 and no scores; MGMT's `run_matkit` returns True after an R failure (`850-869`), and both callers only log a False (`1095-1102`, `1407-1414`).
- VAR-4 (read): `snpsift_output.vcf` is opened with `"w"` before SnpSift runs (`target_analysis.py:3830`); snpEff output is copied over it when the ClinVar DB is missing, SnpSift exits non-zero, or raises (`3845-3887`), and the same for INDELs. The report prints "No pathogenic variants were identified." whenever the files exist (`variants.py:326-335`); the GUI shows a pathogenic count only when it is above 0 (`snp.py:244`). Real SnpSift is not available here (an install was refused by the session's permission check), so the status must not rest on SnpSift's own header lines: ROBIN writes its own.
- REP-6 (read): the report prints "No fusion candidates were identified." when the processed pickles are missing, unreadable or empty (`fusion.py:130-190, 637-645`). Fusion writes no processed pickle when no read passes the filters (`fusion_work.py:2826-2833`); `sv_count.txt` is written whenever the step finishes (`2850-2853`). The methods text says MAPQ > 50, > 200 bp and 4 or more reads (`fusion.py:508`); the code keeps MAPQ > 49, span > 249 and ≥ 3 distinct reads (`classification_config.py:308-320`, `fusion_work.py:1150-1152, 6016-6017`).
- Learnings that bear on this: a new test counts only once seen failing on the old code; a tested helper is not a tested fix (test the call site); GUI code has no harness (say which lines are covered only by import); empty results must still be written.

## Locked metric and decision rule (Gate C — committed before any code)

**Faults.** Each is injected on purpose, synthetic data only. Job-level faults run through a real local Ray `Coordinator` on both schedulers (`standard` and `high`); the synthetic handler is registered on the job type's pool (`standard`) or with `register_handler` (`high`), and reads its fault switch from the job's metadata.

| # | Fault | Must show |
|---|---|---|
| F1 | For each of the 13 job types, the handler files an error under the name the real handler uses (table above) | failed count +1; `failed_by_type[type]` +1; a `failed_jobs.csv` row after shutdown with the error text; an `ERROR_UPDATE` payload naming the job type, the sample and the error |
| F2 | For each of the 13, the handler raises | as F1 (this works today; it must keep working) |
| F3 | CNV batch processing returns an `error_message` | the CNV handler files an error; judged failed |
| F4 | One of two files fails in the modkit merge, and one in the parallel BAM step | the job completes (its next jobs still run); one `failed_jobs.csv` row of kind `file_failed` per failed file; the `ERROR_UPDATE` payload lists it |
| F5 | A worker sends a GUI warning (the alignment warning, the too-many-reads skip) | the GUI process sends the same update type with the same message after the job ends |
| M1 | Sturgeon's classifier raises inside its background step | the Sturgeon handler files an error; judged failed |
| M2 | NanoDX's network raises | judged failed; the store CSV gains no row; no `errordf.csv` anywhere |
| M3 | MGMT's R step: non-zero exit; `Rscript` missing; script missing — in the single and the batch path | the MGMT handler files an error; judged failed; the MGMT BED output is still written |
| V1 | ClinVar DB missing | `snpsift_output.vcf` (and the INDEL file) carry `##ROBIN_ClinVar=unavailable` with the reason; the report says ClinVar annotation was unavailable and does not say "No pathogenic variants were identified."; the SNP display data's summary says `clinvar: unavailable` |
| V2 | SnpSift (a stand-in on `PATH`) exits non-zero | as V1 |
| V3 | SnpSift not found | as V1 |
| V4 | While SnpSift runs, the previous `snpsift_output.vcf` is read | it is the previous file, whole (written to a temp file, then renamed) |
| V5 | A `snpsift_output.vcf` with no `##ROBIN_ClinVar` line (a folder from before this batch) | the report says ClinVar status was not recorded, not "none found" |
| R1 | No fusion output at all | "Fusion analysis has not been run for this sample." |
| R2 | A processed fusion pickle that can't be read | says the fusion results could not be read |
| R3 | Fusion's preprocessing raises | the fusion handler files an error; judged failed |

**Controls (must pass on `92c7c00` and after):**

| # | Run | Must show |
|---|---|---|
| C1 | For each of the 13 job types on both schedulers, a handler that succeeds | completed +1; no failed count, no `failed_jobs.csv` row, no `ERROR_UPDATE` failure |
| C2 | A classifier or bed_conversion handler that files an error on a fail-only BAM submission (`fail_only_bam_submission`, file name with "fail" and not "pass") | not failed: `expected_failure`, as for an exception today |
| C3 | Preprocessing marks `skip_downstream` | counted as skipped, not failed |
| C4 | SnpSift stand-in succeeds; no pathogenic record | `##ROBIN_ClinVar=annotated`; the report says "No pathogenic variants were identified." |
| C5 | Fusion ran and no read passed the filters (new run: an empty processed pickle; old folder: `sv_count.txt`, no pickle, no candidate rows) | "No fusion candidates were identified." |
| C6 | Successful Sturgeon, NanoDX and MGMT runs | NanoDX's row has its scores; MGMT's CSV is written; no error filed |

**Also locked:** the fusion methods sentence states MAPQ ≥ 50, mapping span ≥ 250 bp and ≥ 3 supporting reads, built from `FUSION_DETECTION_THRESHOLDS`; a test checks it against the config and checks the filter at the boundaries (MAPQ 49 dropped, 50 kept; span 249 dropped, 250 kept). The GUI's error panel: `_update_errors` given an `ERROR_UPDATE` payload sets the three counts and replaces "No errors detected." with the newest failures (job type, sample, message), tested on a stand-in object holding the labels.

**Pass, all of:** every F, M, V and R row shows what the table says; every control holds; every new F/M/V/R test fails on `92c7c00`'s code and every control passes there; mutation check (one fix undone at a time, only pytest exit 1 counts as caught, a no-op control mutant must survive) catches every fix, re-run on the result commit.

**Fail:** any one of these. Then the claim is refuted: fix and re-run, or report it.

**Commit order:** this plan (`plan: …`) → code (`code: …`, one commit per task) → `result: … — supported / refuted`, with the test output in the message → records.

## Global Constraints

- Base `92c7c00`; branch `fix/batch-7-failures-visible`; worktree `/home/parlar_ai/dev/robin/.claude/worktrees/batch-7`.
- Tests `/home/parlar_ai/robin-batch3-run/t7.sh` (wrap in `kill-after.sh`); Python `py7.sh`. Baseline 286 passed.
- `git diff 92c7c00..HEAD -- robin.yml pyproject.toml` stays empty; no new dependency; no existing test weakened.
- Synthetic data only; never patient data.
- Out of scope: keeping a failed file's data after the merge (METH-9's deletion; only its visibility is here); marking the classification panel's last result as stale; ClairS-TO's failed regions (VAR-5); import errors (ERR-2); the broad `except` count (ERR-3); GUI forwarding in the non-Ray runner (it runs in the GUI process; only its job judging changes here).
- Review document: Status lines for ERR-1, METH-8, VAR-4, REP-6; fork rows; Updated line; step 7 marker; "Already done". Existing learnings are never edited; add a new note.

## Review Focus

1. **An error filed by an earlier job in the same context** must not fail a later job: only errors filed during the job count. (Task 2)
2. **Fail-only BAM submissions** raise no false alarm on any path, including the single-file paths. (Task 2)
3. **A `target` job that fails** no longer starts `igv_bam`; finalize still refreshes an existing IGV BAM. (Task 2)
4. **Worker messages from two jobs in one worker process**: each still reaches the GUI once. (Task 3)
5. **A sample folder from before this batch** (no `##ROBIN_ClinVar` line; no empty fusion pickle): the report says "not recorded" or reads `sv_count.txt`, never crashes. (Tasks 6, 7)

---

### Task 1: The fault harness and the tests (RED)
**Files:** create `tests/failure_fixtures.py`, `tests/test_job_failures.py`, `tests/test_classifier_failures.py`, `tests/test_clinvar_status.py`, `tests/test_fusion_report_states.py`, `tests/test_gui_errors.py`.
- A module-scoped local Ray (`num_cpus=4`, no dashboard); a Coordinator per scheduler; `run_one(job_type, fault, preset)` registers a synthetic handler (fault switch in the job's metadata), submits one job, waits for it, returns stats and, after shutdown, `failed_jobs.csv`. Every test job carries a temp `work_dir` in its metadata: without one the Coordinator writes `job_durations.csv` into the current directory (seen in the probe).
- The tests for F1–F5, M1–M3, V1–V5, R1–R3, C1–C6, the methods sentence and the error panel. Expected to fail on today's code except the controls.

### Task 2: One rule for a job's outcome (ERR-1 layers 2 and 5; CNV)
**Files:** `workflow_ray.py` (`_wrap_real_handler`, `Pool.enqueue._wait`, the `high` drain loop), `workflow_simple.py` (`WorkflowManager.worker`), `cnv_analysis.py` (`cnv_handler`).
- The wrapper notes how many errors the context held when the job started; afterwards, the errors filed since then decide the outcome, whatever their name; the fail-only rule applies to them as it does to exceptions; the outcome is stored in the context and the fallback result says `failed`. Both Ray judges and the non-Ray runner read it through one helper. CNV files its batch error.
- Tests: F1–F3, C1–C3 turn green.

### Task 3: Failures and worker messages reach the GUI (ERR-1 layers 3 and 4)
**Files:** `gui/app.py` (`send_gui_update`), `workflow_ray.py` (`_wrap_real_handler`, `Coordinator._on_finish`, `stats`, `_publish_gui`), `gui_launcher.py` (`_update_errors`).
- In a process with no GUI, `send_gui_update` keeps the update in a bounded outbox; the wrapper moves it into the job's context; `_on_finish` keeps the recent failures and messages with a running number; `stats` returns them with failure counts by category (preprocessing: preprocessing, bed_conversion; classification: sturgeon, nanodx, pannanodx, random_forest; analysis: the rest). A module-level `publish_failures(stats, state, send)` sends new ones as `ERROR_UPDATE` and as the original update types; `_publish_gui` calls it. `_update_errors` also fills "Recent errors".
- Tests: the payload part of F1, F5, the error panel; a check that `_publish_gui` calls `publish_failures`.

### Task 4: Failed files in the merge are reported (ERR-1 layer 1)
**Files:** `temp_utilities.py` (`merge_modkit_files`), `bed_conversion.py`.
- Both loops return the files that failed, with the error; the bed_conversion handler puts them in the context as file failures; `_on_finish` writes one `failed_jobs` row of kind `file_failed` each and passes them to the GUI; the job's outcome is unchanged.
- Tests: F4.

### Task 5: Classifiers report their errors (METH-8)
**Files:** `sturgeon_analysis.py`, `nanodx_analysis.py`, `mgmt_analysis.py`.
- Sturgeon's background step raises; NanoDX's classification raises (no `errordf.csv`, no empty row); `run_matkit` returns False on any R failure and both callers turn False into the result's error message.
- Tests: M1–M3, C6.

### Task 6: ClinVar status in the VCF (VAR-4)
**Files:** `target_analysis.py` (one `annotate_with_clinvar(...)` used for SNVs and INDELs; `build_snp_display_data`), `reporting/sections/variants.py`, `gui/components/snp.py`.
- Output goes to a temp file, then is renamed; the file's header carries `##ROBIN_ClinVar=annotated` or `=unavailable;reason=…`; the report and the display data read it, and a file without the line is "not recorded". The GUI shows a label when ClinVar is not annotated (import check only; said so in the commit).
- Tests: V1–V5, C4; `run_snp_analysis` calls `annotate_with_clinvar` for both files and no direct copy of snpEff output over the SnpSift path remains.

### Task 7: Fusion report states and its methods text (REP-6)
**Files:** `fusion_work.py` (write an empty processed pickle when nothing passes; a failed preprocessing is filed), `fusion_analysis.py` (handler), `reporting/sections/fusion.py`.
- The report tells not run / unreadable / none found apart, reading `sv_count.txt` and the candidate CSV for folders from before this batch; the methods sentence comes from the config.
- Tests: R1–R3, C5, the methods sentence and the boundaries.

### Task 8: Result and records
- Run the full suite and the mutation check on the result commit; commit `result: … — supported / refuted` with the output.
- CHANGELOG; review Status lines and the fork table; step 7 marker; a learning (Ray workers don't see the driver's `sys.argv`; how to drive a real Coordinator in a test); plan Outcome; `trail check`.
