# Judgment-gate case file: claims this repo got wrong

Each entry: what was claimed, what it rested on, and which gate would have caught it. Add one
whenever a claim is withdrawn or a finding turns out to be an artifact. Concrete, not advice.

## 1. "The CpG filter is 1-based; shift it to 0-based" (v0.5)
- **Claimed:** `b20c2d5` and `docs/releases/v0.5.md:19`: `parquet_filter.txt` was 1-based Illumina,
  so `Start - 1` was added and users were told to rebuild their parquets.
- **Rested on:** an assumption about the file's origin. The file was already 0-based:
  cg13869341 sits at `15864-15866` in it and in the 0-based NanoDX and MNP-Flex BEDs, and about
  486k probes match the NanoDX file exactly (code review METH-6). The shift widened every CpG to
  3 bp, and the first row became `chr1 -1 0`.
- **Gate:** A.5 — compare a few coordinates with an independent 0-based source before shifting.

## 2. "The per-gene variant path fails on chr naming" (batch 3c, 2026-10-31)
- **Claimed:** in conversation and in learning `2026-10-31-clinvar-sites-are-1-based-and-pa.md`,
  `analyze_genes` failed because `1` was not mapped to `chr1`.
- **Rested on:** reasoning from a log line. It does map the name; its pileup range was empty
  (`11..11`), because the gene's start was one base late. Corrected in a new learning.
- **Gate:** D.3 — probe the pileup range before naming a cause.

## 3. Tests that passed on the broken code (batch 3c)
- **Claimed:** a timer test and a two-level nested `$(…)` template test pinned their bugs.
- **Rested on:** never running them on the old code. The timer test re-opened the view every
  tick; the two-level template passed because an odd number of `"` landed right by luck.
- **Gate:** D.2 — run every new test on the old code and watch it fail.

## 4. A predicted label instead of a measured one (batch 3c)
- **Claimed:** an insertion site would read `insertion_requires_review`.
- **Rested on:** reading the code. Measured: `reference_mismatch`.
- **Gate:** D.3 — measure, then write the expected value.

## 5. Bounds from a few sampled points (batch 4, 2026-10-31)
- **Claimed:** the CNV median-rule gap is "at most 0.024 copies"; the u16 wrap ends "under
  about 20,000 reads".
- **Rested on:** five sampled read counts. A full scan found 0.035 copies (at 25,000 reads) and
  wraps up to 24,000 reads.
- **Gate:** D.3 — "at most" needs a scan of the whole range.

## 6. A mutation check that could not report a survivor (batch 5, 2026-11-01)
- **Claimed:** in the plan's Result (Gate D.2), each of nine fixes, undone alone, fails a test.
- **Rested on:** a helper that called any non-zero pytest exit "caught". Given several test
  files as one argument, pytest ran no tests (exit 5) and every mutant printed "caught". One
  claim ("no dedupe") was true at Task 3 and went stale at Task 4, when a semi-join made the
  dedupe redundant; the Result repeated it without a re-run. The final review found both.
- **Gate:** D.2 applied to the checker itself — add a no-op control mutant that must survive,
  count only "N failed" (exit 1) as caught, and re-run every mutant on the commit the claim is
  about.

## 7. Controls described from reading, not from running (batch 7, 2026-11-01)
- **Claimed:** in the locked plan (`31894f0`), controls C2, C4 and C5 "must pass on 92c7c00 and
  after". C2 said a fail-only BAM's filed error is `expected_failure`, "as for an exception today".
- **Rested on:** reading the code. Run on 92c7c00: C2's result is `success` for the classifiers
  and `failed` for bed_conversion (on `standard` also counted failed, though only with a synthetic
  handler: the real one never files there); C4's `##ROBIN_ClinVar=annotated` line and C5's empty
  pickle for a new run don't exist yet. Parts of three controls were new behaviour. Task 1's commit
  saw C2 fail and called it "on purpose", but no ruling squared that with the plan; the Result's
  old-code run found it.
- **Gate:** C, with D.3 — run every control on the old code before locking it. A control that
  can't pass there is a fault row.

## 8. "Cleared after the next good build", tested on a build that rebuilt nothing (batch 7, 2026-11-01)
- **Claimed:** round 3 (`97b34ee`): after a failed fusion build, the earlier-result note "goes away
  after the next good build"; its test passed, and the build loop cleared the `.failed` record on
  any build that returned.
- **Rested on:** a test whose failing build raised in the filter step, before the CSV was appended,
  so the old pickle still matched the CSV. When preprocessing fails after the append, a quiet
  next build rebuilds nothing, clears the record, and the stale pickle reads as current (a
  target-panel fusion then shows as "genome-wide"). Found by the second review's probe, which
  compared the pickle's pairs with the CSV's.
- **Gate:** D.1 and D.2 — when a check says "state X is cleared", test it against the source of
  truth (pickle vs CSV), and on the path where the state is actually stale, not one where it
  happens to be fine.

## 9. "Every ClinVar record is checked", with one input shape measured (batch 9, 2026-11-01)
- **Claimed:** in the locked plan (`4a7f106`, Gate A.5), the real ClinVar file's alleles were
  measured: "multi-ALT records: 0, repeated keys: 0, so the new key is unique". The result
  (`81cf075`) said each record is checked on its own allele.
- **Rested on:** counting only the shape that was feared (a comma in ALT). Among the 367,875
  selected records, 13 have ALT `.` (the reference base itself). The evidence code fails on them
  and drops the whole gene's results. A lone one did this before batch 9 (PTEN, BRCA2), and
  keeping every record at a position added BRCA1. Codex's review found it by scanning every
  selected record's ALT.
- **Gate:** A.5 — for each field the code reads, tally every value shape in the real input
  (here, ALT by alphabet), not the one shape you expect to break.

## 10. "The newest request carries the newest settings", with the default setting exempted (batch 10a, 2026-11-02)
- **Claimed:** in the locked plan (`182c509`), a request sent while one runs replaces the older
  waiting one, because the newest carries the newest settings. The result (`f4ef637`) said
  supported.
- **Rested on:** test S2 checked that the follow-up run was the newest request only with
  batching off. With batching on, the default, it skipped the check, with the comment "the
  batcher re-wraps jobs under new ids". That treated the batcher as a test nuisance, and never
  measured what it does to a request: a batch ran with its *first* request's settings, and a
  request still in the batcher missed the waiting slot, giving three runs. Codex's review found
  both by sending requests with different settings.
- **Gate:** A.5 and D.2 — when a test leaves out an assertion for one parameter, that is an
  input assumption: measure what that path does before exempting it, and above all for the
  default setting.

## 11. "Unmerged batch BAMs are still deleted", a day after they weren't (batch 6 → batch 10b, 2026-11-02)
- **Claimed:** JOB-6's Status line in the code review, and the CHANGELOG's JOB-6 entry, said the
  batch-BAM deletion "when `target.bam` already exists (`target_analysis.py:2604-2618`)" was not
  fixed. Both were still saying so when batch 10b started.
- **Rested on:** the line written in batch 3. Batch 6's `24ebbc8` (TGT-1) deleted that cleanup
  and made finalize merge every new batch BAM, but its records updated TGT-1 only, not the other
  finding that cited the same lines. Found by reading the old code at `eccbf18` beside today's.
- **Gate:** D.4, for records. When a fix removes or rewrites code, search the review for every
  finding that cites that file and those lines, and re-check each one's Status in the same
  records commit.

## 12. "An SNP run calls on every read", with two places reads wait and one way a merge fails (batch 10b, 2026-11-02)
- **Claimed:** in the locked plan (`e027717`), an SNP run first merges every read the sample
  has, and earlier calls count as current only when nothing waits to be merged. The result
  (`910a296`) said supported.
- **Rested on:**
  - "Waiting" meant top-level staged files and batch BAMs. A worker moves the files it claims
    into `_staging/_batch_<id>/` and adds them outside the lock, and nothing looked there.
  - "Failed" meant `status == "error"`. Accumulation also returns `lock_timeout`,
    `load_failed`, and `success` with `batch_bam_error`, and the last two had deleted the
    reads.
  - The tests injected only `"error"`, and staged files only at the top level. Below the
    merge, a failed Clair3 region was skipped, and old VCFs were used.
  - Codex's review found all of it by pausing a real worker mid-claim and failing a real
    batch BAM write.
- **Gate:** A.5. List every place the input can be (here, every folder reads pass through)
  and every value an upstream function returns (read each of its `return`s). Then test each
  one, not the one you expect.

## 13. A 2-second allowance for a rounding nobody measured (batch 10b round 2 → 3, 2026-11-02)
- **Claimed:** in round 2's plan (`8b07d81`), a Clair3 VCF counts only if this run wrote it:
  modified at or after the run's start, "less 2 s for file-system time rounding".
- **Rested on:** an allowance for rounding that was never measured, for a gap no measured
  case needed. It let in what the check was there to stop: a file from a run that ended less
  than 2 s earlier. Codex's second review found it with a quick retry. Round 3 dropped time
  checks and gave each run its own new folder.
- **Gate:** A.4 and D.3. A tolerance is a mechanism. Measure the need before adding it, and
  prefer making the wrong input impossible (a new folder) to making it unlikely (a time window).

## 14. "Counts that can't be read fail the run", with a second way into the counts (batch 10c, 2026-11-02)
- **Claimed:** in the locked plan (`3ff9e57`), a sample's CNV counts that can't be read fail
  the run and are kept, never reset to none. The result (`e108753`) said supported.
- **Rested on:**
  - "The counts" meant the per-sample file. The loader has a second way in: moving the
    sample's counts out of the old multi-sample `update_cnv_dict.pkl`.
  - That move's `try/except`, kept from the original code, logged a warning and went on. So
    a failed move gave no counts, the run succeeded, and it saved over the legacy counts for
    good.
  - The batch's own atomic write added a way to hit it: a failed rename. The tests took only
    the per-sample path.
  - Codex's review found it by failing that rename, through both entry points.
- **Gate:** A.5 again, as in case 12. List every way the input gets in (here, both files the
  loader reads), and check every `except` on those paths against the rule ("a failure says
  failed"), not just the lines being changed.

## 15. "Each BAM's data is added once", with a step's batch and Ray's retries as other ways in (batch 10d, 2026-11-02)
- **Claimed:** in the locked plan (`aad9b4e`) and its result (`41b5c64`): with the ledger, a
  BAM's data is added to its sample's totals once, and a step cut short is reported once.
- **Rested on:**
  - "Which BAMs a step works on" was read from the job's batch list. That list lives in the
    first BAM's job details, which all of that BAM's steps share, so a later batch could
    overwrite an earlier one's before it ran. Older code, but the ledger trusted it.
  - "A step runs once per start record" ignored Ray's own retry of a task whose worker died
    (default 3).
  - The `master.csv` key list was tested with a handful of BAMs, not with a run's worth
    (3,641 keys break Python's `csv`).
  - "Reported once" was enforced by writing the end record before the report was delivered,
    and an ended failure made no report at all.
  - Codex found all four, each through the real batcher, worker or Ray task.
  - Claude found a fifth while writing the brief: the claimed-safe torn-line case broke the
    file one restart later.
- **Gate:** A.5 (inputs verified, not assumed). List every way a step can run, and every
  layer that runs it: the batcher, the runner and Ray itself. Size any growing record at a
  real run's scale, not a test's. And a record that says "reported" must come after the
  report, not before.

## 16. Case 15's lesson, broken in the same day's next part (batch 10d part 3, 2026-11-02)
- **Claimed:** part 3's plan (`721507c`) and result (`f7742f0`): a half-written BAM that
  never completes "shows once as a failed file". Round 2's result (`6b4e991`): the ledger
  records a step's own failed files.
- **Rested on:**
  - The intake marked the BAM reported before the watcher delivered the row, and a failed
    delivery was only logged. That is case 15's last point, written a few hours earlier.
  - "A step's own failed files are the ones after the count taken before it" assumed an
    append-only list owned by one step. The BED step replaces the list, and in the non-Ray
    runner sibling steps share it at once. Errors already used an owner tag for the same
    reason (batch 7).
  - "Counted BAMs are on the list" assumed the list and `master.csv` go together. A
    deleted `master.csv` left the list saying BAMs were counted.
  - Codex's second review found all three, each with the real handler or manager.
- **Gate:** A.3. Re-read the case file's newest entries against the plan being locked, not
  only the old ones. And A.5: for any shared list or file, check who else writes it, and
  how (append, replace, delete).
