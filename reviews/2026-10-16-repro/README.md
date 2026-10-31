# Reproduction scripts for the 2026-10-16 code review

These are the scripts behind the "Confirmed" findings in [`../2026-10-16-code-review.md`](../2026-10-16-code-review.md).

- **Synthetic data only.** No patient data. The VCFs in `cnv_fusion/vcf/` are hand-written ClinVar-style records.
- **Scripts, not tests.** Some files are named `test_*.py`, but they are plain scripts that print what they find. Run them with `python`, not `pytest`.
- **Version.** They were written against `493f53f`.
- **Re-checked.** The 12 scripts marked ✓ below were re-run on 2026-10-16 against `5be9e04` and still show their bug.
- **Next step.** When a finding is fixed, turn its script into a pytest regression test in `tests/` that fails before the fix and passes after it. Then delete the script.

## How to run

Run from the repository root, in a Python 3.12 environment with pandas 2.1.4, numpy < 2, polars, pyranges 0.1.4 and pyarrow 16.1. Ray, NiceGUI, pysam and sturgeon are not needed: the scripts stub them out (see `orchestration/stubs.py` and `gui/stubs*/`). `gui/stubs2/` stays because `gui/test_mgmt_shared_state.py` still uses it; `tests/test_classification_refresh.py` has its own copy under `tests/_stubs/`, so retiring these scripts does not break it.

```
python reviews/2026-10-16-repro/architecture/merge_mean.py
python reviews/2026-10-16-repro/cnv_fusion/test_snp.py reviews/2026-10-16-repro/cnv_fusion/vcf/test.vcf
```

The independent review's suite (IR-01 to IR-11) is in [`independent/`](independent/README.md). It is a single `unittest` file that also needs pysam.

## Which script shows which finding

| Script | Finding | Re-checked |
|---|---|---|
| `architecture/merge_mean.py` | METH-1: mean of means (stored 40.0%, true 78.4%) | ✓ |
| `architecture/merge_first.py` | METH-2: first write not aggregated | |
| `methylation/repro_merge.py` | METH-1, METH-2, METH-5 (join duplicates), METH-6 (Start−1 widening) | |
| `methylation/repro_mean_bias_sim.py` | METH-1: effect on NanoDX/Sturgeon calls over 20 batches | |
| `methylation/repro_matkit.py` | METH-3 (fails in denominator), METH-4 (minus-strand position) | ✓ |
| `methylation/repro_classifier_inputs.py` | METH-4, METH-7 (constant `cov=10`), METH-2 (RF score 200) | |
| `methylation/repro_mgmt_sites.py` | METH-10: 0.43% shown as 43% | |
| `methylation/perf_merge_filter.py` | PERF-8: per-BAM merge timing | |
| `methylation/perf_sturgeon_map.py` | PERF-8: Sturgeon probe-mapping timing | |
| `cnv_fusion/sim_ref_contamination.py` | CNV-1: control polluted within one patient | |
| `cnv_fusion/sim_cross_sample.py` | CNV-1: false events in the next patient | ✓ |
| `cnv_fusion/test_lw.py` + `vcf/clinvar_test.vcf` | VAR-1: ClinVar alleles overwrite each other | ✓ |
| `cnv_fusion/test_snp.py` + `vcf/test.vcf` | VAR-2: "conflicting" counted as pathogenic | ✓ |
| `cnv_fusion/sim_bp_keys.py`, `sim_dbscan.py` | FUS-1: breakpoints never reach support | |
| _(retired)_ `appshell/repro_unique_bed.py <tmpdir>` | VAR-7: `add-panel` keeps one region per gene — **fixed**, replaced by `tests/test_add_panel_regions.py`. Its 8-column case showed a second problem, under no finding ID and not addressed by VAR-7: `read_csv` is given 6 names, so an 8-column BED shifts every field and `_generate_unique_gene_bed` returns `False` on `int('EGFR')`. That is fixed too (`5af0a47`), with command-level tests in the same file, including BED12, which used to be saved garbled. | |
| _(retired)_ `appshell/stale_tbi.sh` | VAR-3: stale tabix index returns nothing — **fixed**, replaced by `tests/test_clinvar_index.py` | |
| `appshell/repro_partial_model_rst.py <tmpdir>` | SEC-3: partial model download accepted | |
| `architecture/partial_dl.py` | SEC-3 (same, with a fake downloader) | |
| `orchestration/test_orch.py` | T1 JOB-3, T2 ERR-1, T3 JOB-1, T4 JOB-5, T5 JOB-13, T6 JOB-4 and JOB-8 | ✓ |
| `orchestration/test_pool_timeout.py` | JOB-9: timeout frees a slot twice | |
| _(retired)_ `gui/test_classification_refresh.py` | GUI-1: panel freezes with a disabled classifier — **fixed**, replaced by `tests/test_classification_refresh.py` | ✓ |
| `gui/test_mgmt_shared_state.py` | GUI-3: second viewer misses updates | |
| `gui/test_summary_cnv.py` | GUI-4 (counter shown as gains) — **fixed** in our fork, now prints `'gained': None, 'lost': None`, regression test `tests/test_summary_cnv_counts.py`; GUI-7 and GUI-9 remain, so the script stays | ✓ |
| `gui/time_fusion_validate.py`, `time_tostring.py` | PERF-1: fusion refresh cost | |
| _(retired)_ `reporting/rf_conf.py` | REP-1: wrong Random Forest cut-offs — **fixed**, replaced by `tests/test_classifier_config_keys.py` | ✓ |
| `reporting/patho.py` | VAR-2 (report side) | ✓ |
| `reporting/cyto.py`, `cyto2.py` | REP-2: old ±1 SD rule | ✓ (`cyto.py`) |
| `reporting/fusion_logic.py` | REP-3 (weakest cluster), REP-4 (hyphenated genes) | |
| `reporting/cov_init.py` | REP-9: coverage file aborts the report | |
| `reporting/ylim.py` | REP-7: amplifications clipped | ✓ |
| `architecture/bigfuncs.py`, `dupdefs.py`, `imports_vs_deps.py` | ERR-3, ARCH-2, ARCH-4, PKG-6 (code counts) | |
