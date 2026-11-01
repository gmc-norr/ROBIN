# Batch 9: which ClinVar records count as pathogenic, and that each one is checked (VAR-1, VAR-2)

**Locked before code, 2026-11-01.** Base `cf4f84e`. Branch `fix/batch-9-clinvar`. Done by the controller (no agents, to save the user's tokens). An external review (Codex / Sol 6.1 / Astra) comes before the merge.

## Rule (agreed by the user, option "A", 2026-11-01)

A ClinVar record counts as pathogenic when ClinVar's own aggregate label says so in any of its three classifications:

- **CLNSIG** (germline): `Pathogenic` or `Likely_pathogenic`.
- **ONC** (oncogenicity): `Oncogenic` or `Likely_oncogenic`.
- **SCI** (somatic clinical impact): `Tier_I_-_Strong` or `Tier_II_-_Potential`.

Labels are matched whole, never as substrings. A value is split into single labels on `|`, `/` and `,`, so `Pathogenic/Likely_pathogenic` and `Pathogenic,_low_penetrance` count. `Conflicting_classifications_of_pathogenicity` never counts by itself. There is no vote over CLNSIGCONF: that would be ROBIN's verdict, not ClinVar's.

Every ClinVar record is checked on its own allele (VAR-1). A record is keyed by `(chrom, pos, ref, alt)`, not `(chrom, pos)`.

**Spec:** VAR-1 and VAR-2 in `reviews/2026-10-16-code-review.md`, and the user's choice above. The review's own suggested fix ("match exact CLNSIG tokens") is not used as written: measured below, it drops BRAF V600E and TERT C228T.

## Gate A

1. **Principle:** ROBIN calls a site pathogenic only when ClinVar does, in any of ClinVar's three classifications, and checks every record it loads.
2. **What would kill it:** any of these.
   - BRAF V600E or TERT C228T stops counting at any of the three checks.
   - A conflicting-only record (TP53, NF1 below) still counts.
   - An R132C-type sample still reads anything but found, with the IDH1 c.394 records in ClinVar's file order.
3. **Refuted before here?** No case-file entry is about the ClinVar rule. The learnings on 1-based sites and the boundary rule are about where sites are picked, which this batch doesn't change. The review's "exact CLNSIG tokens" is refuted by the measurement below.
4. **Simplest thing:**
   - One small helper with three label sets.
   - It is used by the variant check (LGA), the GUI's SNP table and the report.
   - One key change for VAR-1.
5. **Inputs (measured on ClinVar's GRCh38 VCF, fileDate 2026-10-25, 4,599,696 records):**
   - **CLNSIG labels**, after splitting on `| / ,`:
     - Uncertain_significance, Likely_benign, Benign
     - Pathogenic (235,021), Conflicting_classifications_of_pathogenicity (182,151), Likely_pathogenic (176,753)
     - not_provided, drug_response, other, no_classification_for_the_single_variant, risk_factor, association
     - Uncertain_risk_allele, `_low_penetrance`, Affects, no_classifications_from_unflagged_records, Likely_risk_allele
     - VUS-high/mid/low, protective, confers_sensitivity, Established_risk_allele, association_not_found
   - **ONC values:** Likely_oncogenic, no_classification_for_the_single_variant, Uncertain_significance, Oncogenic, Benign, Likely_benign, Conflicting_classifications_of_oncogenicity, Oncogenic/Likely_oncogenic.
   - **SCI values:** Tier_I_-_Strong, no_classification_for_the_single_variant, Tier_II_-_Potential, Tier_III_-_Unknown, Tier_IV_-_Benign/Likely_benign.
   - **Alleles:** multi-ALT records: 0. Repeated `(chrom, pos, ref, alt)` keys: 0, so the new key is unique.
   - **Not measured:** SnpSift copies ONC, SCI and CLNSIGCONF into `snpsift_output.vcf`. ROBIN runs `SnpSift annotate <db> <vcf>` with no `-info`, and SnpSift's documented default copies every INFO field. SnpSift is not installed here (and must not be).

Scripts and outputs: session scratchpad `b9/` (`measure_clinvar.py`, `hotspots.py`, `rules.py`, `examples.py`, `probe_var1.py`). The ClinVar file is in `b9/clinvar_dl/`.

## Gate B: measured on `cf4f84e` before this lock

- **Today's three checks** are the variant check's `_is_pathogenic`, `snp_processing._is_pathogenic` and the report's `VariantAnalysis._is_pathogenic`.
  - 182,153 records have "conflicting" in CLNSIG or ONC. The variant check calls all of them pathogenic; the other two call 182,152 of them pathogenic.
  - In all, the variant check calls 549,186 records pathogenic, and the other two call 548,359.
- **The rule above, compared with the variant check today:**
  - It drops 182,012 records (conflicting, with no oncogenic or Tier I/II label).
  - It adds 701 records: their oncogenic or Tier I/II label is one today's check ignores.
  - It keeps 367,174.
  - A CLNSIGCONF vote would keep about 30,000 more conflicting records (option B, not chosen).
- **Hotspots** (ClinVar labels, GRCh38):

  | Site | CLNSIG | ONC | SCI | Counts under the rule |
  |---|---|---|---|---|
  | BRAF V600E (7:140753336 A>T) | Conflicting | Oncogenic | Tier I | yes, by ONC/SCI |
  | TERT C228T (5:1295113 G>A) | Conflicting | Oncogenic | Tier I | yes, by ONC/SCI |
  | TERT 5:1295113 G>T | — | — | Tier I | yes (today: no) |
  | an H3-3A record (1:226064454 G>A) | other | Likely_oncogenic | Tier I | yes (today: no) |
  | IDH1 c.395 C>T / C>A, c.394 G>A / G>T | Pathogenic or P/LP | Oncogenic | Tier I/II | yes |

- **VAR-1:**
  - 53,674 positions hold more than one record that the variant check calls pathogenic today. 67,276 of its 549,186 records are overwritten.
  - At IDH1 c.394, today's code keeps only G>T. A synthetic sample with 20 reads all carrying A (R132C-type) reads "reference_mismatch" and the gene "insufficient_evidence".
  - At KRAS c.35, today's code keeps only CC>AA, which a one-base pileup can't call.

## Fixes

1. **New** `src/robin/analysis/clinvar_significance.py` (no pandas):
   - `is_clinvar_pathogenic(info)` takes the INFO string or a dict (values as strings, or lists of strings) and applies the rule above. Case is ignored.
   - `clinvar_label_text(info, default="unknown")` gives the labels for display: CLNSIG as written, then `oncogenicity: <ONC>` and `somatic impact: <SCI>` when present, joined by `; `. A CLNSIG-only record reads exactly as before.
2. **Variant check** (`lightweight_gene_analysis.py`):
   - `_is_pathogenic` uses the helper. This also drops the old "ONC only when there's no CLNSIG" branch.
   - The record key becomes `(chrom, pos, ref, alt)`, and `_find_bed_clinvar_intersections` reads the 4-tuple.
   - `clinical_significance` becomes `clinvar_label_text(...)`. The GUI and the HTML/JSON outputs show it as they do today.
3. **SNP table** (`snp_processing.py`):
   - `_is_pathogenic` uses the helper.
   - `parse_vcf` groups with `dropna=False`, so a record without ANN stays.
4. **Report** (`reporting/sections/variants.py`):
   - `_is_pathogenic` uses the helper. The CLNSIGCONF vote branch goes: it was never reached, and its `int("pathogenic(3")` would have raised.
   - The note under the table states the rule.
5. **GUI SNP table** (`gui/components/snp.py`): `ONC` and `SCI` follow `CLNSIG` in `preferred_display_fields`. A row flagged by ONC or SCI while CLNSIG says "Conflicting" then shows why.

## Tests: `tests/test_clinvar_significance.py`, written and run on `cf4f84e` before this lock

The INFO values are copied from the ClinVar file above (classification keys and GENEINFO only). The reads are synthetic.

On `cf4f84e`: **53 failed, 19 passed**, as planned.
- 32 failed by `ModuleNotFoundError` (the helper's table, which can't exist on the old code).
- 21 failed by assertion.
- The 19 that passed are the guards and controls below.

| # | Test | On `cf4f84e` |
|---|---|---|
| T1 | `test_the_rule_on_clinvar_labels`: 32 labels across CLNSIG, ONC and SCI, plus labels under CLNSIGCONF/CLNSIGINCL, which never count. String and dict input. | 32 failed (no module) |
| T2 | `test_every_check_follows_the_rule`: 10 records × the three checks. | 14 failed, 16 passed |
| T2 faults | TERT SCI-only and H3-3A (×3: today "not pathogenic"); KRAS ONC-only (SNP table and report); TP53 and NF1 conflicting-only (×3: today "pathogenic"). | failed |
| T2 guards | BRAF V600E and TERT C228T (×3: pathogenic today by the substring accident, and must stay so); plain Pathogenic, Benign and not-in-ClinVar (×3); KRAS ONC-only in the variant check. | passed |
| T3 | `test_the_report_lists_the_records_clinvar_calls_pathogenic`: 7 rows. Listed: BRAF, TERT ×2, H3-3A, plain. Not listed: TP53, NF1. Genes, and the note states the rule. | failed |
| T4 | `test_the_snp_table_flags_the_records_clinvar_calls_pathogenic`: the same 7 rows with ANN; the flags and counts. | failed |
| T5 | `test_a_record_without_a_gene_annotation_stays_in_the_snp_table`: 3 rows, one without ANN. `parse_vcf` and the table keep all 3. | failed (2 of 3) |
| T6 | `test_the_gui_snp_table_shows_the_oncogenicity_and_somatic_labels` (source check: the columns are picked inside a nested GUI function). | failed |
| T7 | `test_every_clinvar_allele_at_one_position_is_checked`: IDH1 c.394, three records in file order, reads all A. All three alleles are listed. G>A has 20 alternate reads, high_frequency, supported. Gene: present_with_good_coverage. | failed |
| T8 | `test_a_snv_next_to_a_multi_base_record_at_the_same_position_is_checked`: KRAS c.35, four records, reads all T. C>T supported; CC>AA mnv_requires_review. | failed |
| T9 | `test_the_variant_check_uses_the_rule_and_shows_why`: a TP53-type conflicting record is not checked. An H3-3A-type record is checked and supported, with its labels shown in `clinical_significance`. | failed |
| K1 | `test_control_report_with_plain_pathogenic_rows`: Pathogenic listed, Benign not. | passed |
| K2 | `test_control_every_record_annotated_keeps_its_rows`: with ANN on every row, nothing changes. | passed |
| K3 | `test_control_a_plain_pathogenic_record_keeps_its_label`: one record, CLNSIG only. Supported, and the label reads `Pathogenic` as before. | passed |

No existing test changes.

## Mutants (one per fix part), each must be caught, and the no-op control must survive

- **M1:** the variant-check key is `(chrom, pos)` again.
- **M2:** the helper matches substrings, not whole labels.
- **M3:** the helper ignores ONC.
- **M4:** the helper ignores SCI.
- **M5:** the helper doesn't split on `,`.
- **M6, M7, M8:** each of the three checks goes back to its old code, one at a time.
- **M9:** `dropna=False` is removed.
- **M10:** `clinical_significance` is CLNSIG only again.
- **M11:** the GUI's ONC/SCI fields are removed.
- **M12:** the report's note is back to the old text.
- **Control:** a no-op change.

Only exit 1 ("N failed") counts as caught (case 6).

## Pass, all of

- T1–T9 pass after. K1–K3 and the T2 guards pass on `cf4f84e` and after.
- M1–M12 are each caught, and the control survives.
- The full suite passes (661 on `cf4f84e`, plus these).
- `robin.yml` and `pyproject.toml` are unchanged.

**Fail:** any one; then fix and re-run, or report it.

## Out of scope (recorded)

- **Saved results from before this change keep their old flags.**
  - The GUI's saved SNP table (`clair3/snpsift_output_display.json`) changes when the sample's SNP step runs again.
  - Saved variant-check results (`lightweight_gene_analysis_results.json`) change when the check runs again.
  - The report reads the VCF each time, so it follows the new rule at once.
- **A ClinVar site inside two differently named BED regions** is checked under the first region only (`break` in `_find_bed_clinvar_intersections`).
- **The report's CSV export** fails (and logs) when no listed row has HGVS fields. This is pre-existing. Real ClinVar-annotated rows carry CLNHGVS and snpEff's ANN.
- **`reviews/2026-10-16-repro/independent/reproduce.py`:**
  - Its variant-check repros (IR-01, IR-02, IR-09) already fail on `cf4f84e`, because they reproduce bugs that are now fixed.
  - IR-01 fails there by `NameError: datetime`. Not touched.
- **Whether SnpSift matches ClinVar records by REF/ALT as well as position** (it decides which record's labels a sample row gets). Not checked here.

## Result: supported (code `502fd72`)

Every pass condition holds. The commands and outputs follow.

- **New tests on the base:** `tests/test_clinvar_significance.py` at `502fd72`, run on `cf4f84e`'s `src` (`git archive cf4f84e src`): `53 failed, 19 passed`. The failures are 32 `ModuleNotFoundError` and 21 by assertion, as locked.
- **New tests after:** `t9.sh tests/test_clinvar_significance.py`: `72 passed`.
- **Mutants** (`b9/scripts/mutate9.py`; only exit 1 counts):
  - M1–M12 are each `CAUGHT`, each with 1 failed under `-x`.
  - The no-op control `SURVIVED` (163 passed).
  - The mutated files: M1, M6, M10 the variant check; M2–M5 the helper; M7, M9 `snp_processing.py`; M8, M12 the report; M11 the GUI.
- **Full suite:** `733 passed` (661 on `cf4f84e` + 72).
- **Pins:** `git diff cf4f84e -- robin.yml pyproject.toml` is empty.
- **No existing test changed.**

**Real ClinVar file, re-measured on `502fd72`** (`measure_clinvar.py` with the worktree's `src`):
- All three checks flag **367,875** records. That is 367,174 kept + 701 added, as Gate B's `rules.py` predicted.
- Of the 182,153 records with "conflicting" in CLNSIG or ONC, 141 still count (by an oncogenic or Tier I/II label), and 182,012 don't.
- The three checks now give the same count. Before, the variant check's count was 827 higher than the other two.

**IDH1 probe on `502fd72`** (`probe_var1.py`):
- All 3 records are kept. G>A has 20 alternate reads and reads high_frequency, so the gene is present_with_good_coverage.
- Measured side effect: G>C and G>T read `reference_mismatch` (needs review) there, because every read carries A, which is neither allele.
  - In a sample whose reads show G as well, they read `no_support`. Measured at the records step (`probe_mixed.py`, 10 G + 10 A reads): G>A has 10 alternate reads and is supported; G>C and G>T read `no_support`, are not supported and need no review; the gene reads present_with_good_coverage.
  - The gene's status follows the supported site, so it is unchanged. It is recorded, not changed: a sibling allele's evidence is still judged on its own two alleles.

Gate D:
1. **Could the pipeline make this up?** The rule is checked on ClinVar's own labels (copied verbatim) and on the whole real file, not only on hand-made labels.
2. **Do the controls tell working from broken?** Yes. The no-op control survived, and every mutant of a fix part was caught.
3. **Measured or inferred?** Measured, except that SnpSift copies ONC/SCI into `snpsift_output.vcf` (its documented default; not run here).
4. **Fresh output:** above.

## External review (Codex, `review-requests/batch9-report-codex.md`): fix list locked before code

Two findings, both verified here with my own measurements. I did not run the reviewer's scripts.

### Finding 1 (Important): a record with no alternate allele drops its whole region

- **Measured on the real file** (`no_alt.py`): 13 of the 367,875 selected records have ALT `.`; every other selected ALT is plain bases.
  - ClinVar's HGVS for them is `g.N=`: the record is the reference base itself.
  - Examples: PTEN 10:87864144 (AML panel), BRCA2 13:32337606 (both panels), BRCA1 17:43067696.
  - Three share a position with other selected records: UGT1A 2:233767046, BRCA1 17:43067696 (with C>A, C>G, C>T) and COL1A1 17:50201450.
- **Why it drops the region:** such a record gets `alternate=None`. Both evidence paths call `len()` on it, with reads or without. `analyze_bed_regions` catches the error per region, logs it and drops the region, so the JSON has no entry for that gene.
- **Measured with the new tests** (synthetic reads, the real labels):
  - **On `cf4f84e`:** a lone PTEN-type G>. drops the PTEN region, covered or not. This predates batch 9. The BRCA1 quartet "works" there only because the position key hid C>. behind C>T, along with C>A and C>G.
  - **On `502fd72`:** the BRCA1 region is dropped too. Batch 9 made this reachable for the three shared positions.
- **Fix R1:** the loader skips a selected record whose ALT is not plain bases (A, C, G, T, N).
  - It counts them, logs one warning with the count, and records the count in the results metadata as `clinvar_records_not_checkable`.
  - **Ruling (mine):** skip, not "needs review". The record is the reference base itself, so no read can support or refute it. As "needs review", every PTEN and BRCA2 result with nothing found would read "insufficient evidence" for good, for a reason that has nothing to do with the sample.
- **Recorded, not fixed:** any other error inside a region still drops the region with only a log line. This predates batch 9 and belongs with the "failures visible" follow-ups. After R1, no known record shape triggers it.

### Finding 2 (Minor): the HTML export counts every checked record as found

- `generate_report` / `_generate_html_report` count `len(variants)` and coverage, not the evidence.
  - On `502fd72`, three IDH1 sibling records with reads all A read "PATHOGENIC SNPS PRESENT (3 with good coverage)".
  - Reference-only reads read the same.
- This predates batch 9 (VAR-8 fixed the JSON only); VAR-1 makes it larger. Nothing in ROBIN calls it (graft: 0 callers), but it is public, as `analyze_genes` was when batch 3c fixed it.
- **Fix R2:** one module function, `gene_pathogenic_status(high, total, review)`, gives the gene status for both exports.
  - The JSON summary calls it with the counts it makes today, so its output is unchanged.
  - The HTML counts found variants with `is_supported`, and sites to review with `needs_review`, and shows the status the JSON shows. The list's heading says "ClinVar sites checked".

### Tests (added to `tests/test_clinvar_significance.py`; run on `cf4f84e` and `502fd72` before this lock)

| # | Test | `cf4f84e` | `502fd72` |
|---|---|---|---|
| T10 ×2 | `test_a_record_without_an_alternate_allele_does_not_drop_its_region` [site covered / uncovered]: the BRCA1 quartet in file order, plus a supported C>T at 121. BRCA1 is present, the three plain siblings and the neighbour are listed, C>. is not, the status is present_with_good_coverage, and the metadata count is 1. | failed (only C>T and the neighbour listed) | failed (`'BRCA1' in {}`) |
| T11 ×2 | `test_a_lone_record_without_an_alternate_allele_does_not_drop_its_region` [covered / uncovered]: PTEN G>. (and the Benign G>C), plus a supported G>A at 121. PTEN is present, only the neighbour is listed, and the count is 1. | failed (`'PTEN' in {}`) | failed (`'PTEN' in {}`) |
| T12 ×3 | `test_the_html_export_counts_what_the_json_counts` [all alt / reference only / two alt reads] on the IDH1 trio. The JSON status and found count are as expected, and the HTML shows the same status phrase and "Total Pathogenic Variants: n". | failed | failed |

All 7 fail by assertion on both. No existing test changes. The existing variant-check tests guard the JSON summary refactor.

### Mutants (round 2)

- **N1:** the loader keeps records whose ALT is not plain bases.
- **N2:** the metadata count is not written.
- **N3:** the HTML counts every checked record as found.
- **N4:** the HTML ignores sites to review.
- **N5:** `gene_pathogenic_status` drops its "insufficient evidence" branch, which the JSON uses too.
- **Control:** a no-op change.

All of M1–M12 are re-run on the new code too.

### Pass, all of

- T10–T12 pass.
- N1–N5 and M1–M12 are each caught, and the control survives.
- The full suite passes.
- The pins are unchanged.

### Result: supported (code `b97a709`)

- **Before:** T10–T12 failed on `cf4f84e` and on `502fd72` (7 each, by assertion; measured before the lock).
- **After:** the variant-check tests (`test_clinvar_significance.py`, `test_lga_*`, `test_atomic_writes.py`): `147 passed`.
- **Mutants** (`b9/scripts/mutate9b.py`; only exit 1 counts): N1–N5 and M1–M12 are each `CAUGHT`, and the no-op control `SURVIVED` (177 passed).
- **Full suite:** `740 passed` (733 + 7).
- **Pins:** `git diff cf4f84e -- robin.yml pyproject.toml` is empty.
- **No existing test changed.**
- **The real ClinVar file, loaded with the variant check's own loader** (`load_real.py`): 367,862 records are kept and 13 skipped, which is 367,875 in all, as the rule selects. No ALT `.` record is loaded, and BRCA1 17:43067696 keeps A, G and T. The log has one warning: "13 pathogenic ClinVar records have no alternate allele (ALT '.') and can't be checked with a pileup; skipped".
- **Case file:** entry 9 added. The plan's Gate A.5 measured only the shape it feared in ALT (a comma), and missed `.`.
