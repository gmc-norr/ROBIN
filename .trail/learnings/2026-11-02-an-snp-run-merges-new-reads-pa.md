---
title: "An SNP run merges new reads first, and earlier SNP calls are reused only when they came from the same target.bam"
author: Pär
date: 2026-11-02
branch: fix/batch-10b-finalize
cost: { minutes: 40, tokens: 348146 }
id: l-491f45c6
scope: repo
status: active
src: agent
touches:
  - reviews/plans/2026-11-02-batch-10b-finalize-and-snp.md
  - src/robin/analysis/target_analysis.py
  - src/robin/gui_launcher.py
  - tests/test_finalize_once.py
  - tests/test_snp_reads_current.py
  - CHANGELOG.md
  - reviews/2026-10-16-code-review.md
  - reviews/2026-10-16-code-review.pdf
  - reviews/plans/2026-11-02-batch-10a-review-fixes.md
  - src/robin/workflow_ray.py
---

## Decided
- The SNP job itself merges new batch BAMs and staged files into target.bam before calling variants, so all three SNP buttons (finalize+SNP, the coverage panel, bulk SNP) get the same fresh input. Fixing only the GUI's finalize path would have left the other two reading a stale target.bam.
- Earlier calls are reused only when a record written by the finished call (the batch names its target.bam held) equals target.bam's current batch list and no batch BAM or staged file waits. This is the same @CO batch list batch 6 uses for the IGV copy; file times were ruled out there and here, because a merge can land while a call runs.
- The record is deleted before a real (non-annotation-only) call starts, so a run that stops part-way never leaves a matching record. Annotation-only runs write none: they call nothing.
- The GUI's in-memory set of finalized samples was deleted, not fixed: its only reads were the two skips that caused JOB-6. The on-disk "nothing to merge" check does that job, and now also counts staged files.

## Tried and ruled out
- Checking JOB-7's damage with the tiny fixture BAMs as they are: the merge takes milliseconds, so two finalizes never overlapped and the control (lock off) was clean too (0 of 40). A random 0.1-0.6 s pause before each merge made them overlap: 39 of 40 bad without the lock, 0 of 40 with it. A race check needs a control that shows the damage first.

## Watch out
- The SNP job now needs the sample's panel BED (the merge builds a TargetAnalysis); a panel that can't be found fails the SNP job with "Could not merge new reads into target.bam".
- Folders from before batch 10b have no record, so each sample's next SNP run calls again, once.
- GUI tests can build a real GUILauncher() without NiceGUI; swap threading.Thread for a class whose start() runs the target, after building it.

In the code: src/robin/analysis/target_analysis.py (snp_calls_are_current, _write_snp_source, clair3/target_bam_batches.txt; snp_analysis_handler's merge; run_snp_analysis's reuse check); src/robin/gui_launcher.py (_trigger_target_bam_finalization's merge_here and nothing_to_merge; _is_target_bam_finalize_redundant; _sample_has_snp_calling_outputs). Tests: tests/test_finalize_once.py, tests/test_snp_reads_current.py.
