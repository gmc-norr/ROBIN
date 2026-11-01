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
