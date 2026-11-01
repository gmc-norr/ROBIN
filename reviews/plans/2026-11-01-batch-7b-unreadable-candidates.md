# Batch 7b: a fusion candidate file that can't be read is a failure, not "none found"

**Locked before code, 2026-11-01.** Base `287172f` (batch 7 merged). Branch `fix/batch-7b-unreadable-candidates`.

## Why

Finding 2 of the external review of batch 7 (`/home/parlar_ai/robin-batch3-run/review-requests/batch7-round4-report-codex.md`), given its own batch by the user. The rule is batch 7's: when a step fails, the user sees "failed" or "not available"; "none found" only when it ran and found nothing. A new false alarm is a failure too.

## What was measured before this lock (Gate B, synthetic reads, the real writer, no Ray)

Probe: `scratchpad/b7b/probe_unreadable.py`.
- **A damaged staged target file** (next to good genome-wide staged files) gives:
  - status "success";
  - the target panel reads "No fusion networks identified within the target panel.";
  - the damaged file is deleted by the cleanup, so that BAM's target-panel reads are gone for good.
- **A damaged older part** of the accumulated dataset did no harm: its reads were already counted in the saved filter state.
- **The reviewer's measurement:** an `OSError` reading the target dataset in either pass of the build gives "success", an empty target pickle, no record, and "No fusion networks identified within the target panel.".

## Why the errors were skipped (Gate A, read)

- `_append_staged_files` and `_iter_fusion_candidates_parquet_batches` catch every read error, log it and go on.
- Staged files (`process_bam_with_staging`) and dataset parts (`_append_fusion_candidates_parquet`) are written straight to their final names.
- Accumulation takes no lock: any worker whose staging count reaches the threshold may run it ("Another worker may have already accumulated").
- So a reader can meet a half-written file. Skipping hides that, and hides real damage too. Raising without changing the writes would turn that race into a false alarm, so the writes change as well.
- Two accumulations at once are a wider problem (double appends, state races). It belongs to step 8 (JOB-1, JOB-2, CNV-2) and is not fixed here.

## Locked

- **F1, writes:** staged files and dataset parts are written under a temporary name (`<name>.parquet.tmp`, which no `*.parquet` glob matches) and then renamed into place, so a reader never sees a half-written file. A helper `_write_parquet_atomically` serves both writers.
- **F2, staged files:** before anything is appended, every staged file of all three kinds is opened. If one can't be read, the accumulation raises and names it. Nothing is appended, and the staged files stay for the next try.
- **F3, the dataset:** `_iter_fusion_candidates_parquet_batches` raises a read error, naming the file, instead of skipping it. That type's build fails, and its `.failed` record is written.
- **F4, the CSV:** the build's second pass appends to the candidate CSV only after every file has been read. A read failure part-way leaves the CSV as it was, so a retry can't append the same rows twice.
- **F5, an accumulation that fails before its builds** (staging or append) writes the `.failed` record for both types. The report then says the latest build failed beside earlier results, or "could not be read" with none. A failure inside one type's build still marks only that type.

## Tests (each fails on `287172f` unless marked control or guard)

1. A damaged staged target file, first build: status "error", the damaged file still staged, the target-panel summary not "No fusion networks identified within the target panel.".
2. The same, after an earlier good build: the summary carries "The latest target-panel fusion build failed; these results are from an earlier run." (F5).
3. A dataset read error (`pq.ParquetFile` raising for the target dataset), parametrized over the first and second pass, first build: status "error", the target record written, the target panel not "none" (F3).
4. F4: a good build 1, then build 2 with two new parts and a read error on the second part in the second pass, then build 3 with nothing failing. After build 3, the target CSV holds each read's rows once, not twice.
5. **Guard:** after a staged-file failure, repairing the file and building again gives both fusions, with each read's rows once in both CSVs (F2: nothing was appended by the failed build).
6. F1: `_append_fusion_candidates_parquet` and `process_bam_with_staging` write through `_write_parquet_atomically`, which renames a `.tmp` file into place. Its test checks the rename, and the writers are checked on their source.
7. **Control:** a failure inside the target panel's build only (the round-3 setup) leaves no record for the genome-wide type. The existing `test_earlier_candidates_after_a_failed_build_say_so` pins that.

## Pass, all of

- Each test fails on `287172f` and passes after; controls and guards pass on both.
- One mutant per fix (F1–F5) is caught, and the no-op control survives.
- The full suite passes.
- `robin.yml` and `pyproject.toml` are unchanged.

**Fail:** any one; then fix and re-run, or report it.
