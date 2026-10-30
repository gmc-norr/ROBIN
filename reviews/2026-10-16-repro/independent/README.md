# Independent review test suite (IR-01 to IR-11)

These tests come from the independent review by GPT 6.1 Sol. Its findings are merged into [the code review](../../2026-10-16-code-review.md). Each `test_irNN` matches the finding marked *(IR-NN)* there.

| Test | Finding in the review |
|---|---|
| `test_ir01` | VAR-8: reference-only samples summarised as pathogenic |
| `test_ir02` | VAR-9: missing BAM falls back to another sample's BAM |
| `test_ir03` | TGT-3: global depth depends on how reads are split into BAMs |
| `test_ir04` | TGT-4: concurrent target accumulations overwrite coverage |
| `test_ir05` | METH-13: implicit canonical MM calls ignored |
| `test_ir06` | METH-14: QS filtering depends on alignment order |
| `test_ir07` | TGT-5: supplementary alignments counted as reads |
| `test_ir08` | SEC-8: MGMT runs the sample path through the shell |
| `test_ir09` | VAR-10: repeated BED names reduced to the first interval |
| `test_ir10` | TGT-2: CNV master BED can drop the selected panel |
| `test_ir11` | TGT-1: early BAMs never reach the SNP input |

## How to run

Use Python 3.12 on Linux. The concurrency test uses POSIX `fcntl` locks.

```sh
python3.12 -m venv /tmp/robin-independent-review
/tmp/robin-independent-review/bin/python -m pip install -r reviews/2026-10-16-repro/independent/requirements.txt
/tmp/robin-independent-review/bin/python reviews/2026-10-16-repro/independent/reproduce.py -v
```

## How it works

- The suite pulls the real function definitions out of the source with Python's `ast`, so Ray, NiceGUI, models, R and Docker are not needed.
- These parts are real, not stubbed:
  - BAMs and their indexes
  - coverage commands and pileups
  - parquet files and file locks
  - the shell command that the MGMT test triggers
- Every file it creates stays in a temporary folder, which is removed afterwards.

**The assertions confirm the existing bugs**, so `OK` means each bug was reproduced. When a finding is fixed, turn its test into a normal regression test in `tests/` that requires the correct behaviour.

**Verified:** 11 of 11 reproduced, once by the independent reviewer at `5be9e04` and again during the merge at `6593336`. Both runs used Python 3.12.3, NumPy 1.26.4, pandas 2.1.4, pysam 0.23.3 and PyArrow 16.1.0. No patient data, live Ray scheduling or GUI was used.

**Fixed since:** in our fork, these tests now fail. That is expected: the bugs they reproduce are gone.

| Test | Fixed by | Regression test |
|---|---|---|
| `test_ir08` (SEC-8) | `2d3577e` | `tests/test_mgmt_rscript.py` |
| `test_ir02` (VAR-9) | `b9ffef8` | `tests/test_lga_sample_inputs.py` |
| `test_ir09` (VAR-10) | `20ea771` | `tests/test_lga_intervals.py` |
| `test_ir01` (VAR-8) | `aeb5339` | `tests/test_lga_variant_support.py` |

Each of the four fails on its own assertion, which shows the fixed behaviour: `test_ir01` gets `none_detected`, `test_ir02` stops at A's missing `igv_ready.bam`, `test_ir08` finds no marker file, and `test_ir09` finds both sites, 11 and 111. The other 7 tests still reproduce.

This suite copies in only the functions it names. The fixed code calls new helpers, so `is_supported`, `needs_review` and `resolve_sample_folder` were added to its lists in `reproduce.py`. Without them the LGA tests stop with a `NameError` or `KeyError` instead. No assertion was changed.
