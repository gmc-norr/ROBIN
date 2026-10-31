---
title: "The CNV library changes the counts it is given and adds in 16 bits, so the control profile is built in numpy"
author: Pär
date: 2026-10-31
branch: fix/batch-4-cnv-reference
cost: { minutes: 1, tokens: 2659 }
id: l-089e9553
scope: repo
status: active
src: agent
touches:
  - src/robin/analysis/cnv_subprocess.py
---

## Decided
- The CNV control's profile is built in numpy from the control's stored counts (64-bit sums), not by cnv_from_bam. `iterate_bam_file(bam, copy_numbers=d)` merges the BAM's counts into `d` and replaces its values, so any dict you pass is changed. The cached control is a read-only mapping of read-only arrays: the library refuses it with a TypeError, which is the guard.
- The numpy formula mirrors the library's `calculate_cnv`: round the bin width to whole 1 kb bins, sum each contig's bins in chunks (last chunk partial), divide by the median of all chunk sums rounded half away from zero (Rust `f64::round`; Python's `round` goes to even), 0 becomes 1, times ploidy 2. On one contig, or contigs whose lengths are whole chunks, it matches the library exactly.
- The median chunks stay within each contig. The library chunks the genome end to end in Rust hash-map order, so its chunks can straddle contigs; that order can't be mirrored. Measured on the public control: at most 0.024 copies apart.

## Tried and ruled out
- The library's own no-BAM mode, `iterate_bam_file(None, copy_numbers=d, bin_width=w)`: it reads `d` without changing it and uses the library's exact maths, but it sums chunks as u16 (release build, no overflow checks). With the 12.7 M-read control and a sample of up to ~24,000 reads (wide bins), chunk sums pass 65,535 and wrap; at 4,000 and 20,000 reads the profile was off by up to 4 and 17 copies. ROBIN's old pass had the same wrap.
- A `copy.deepcopy` of the control per call: the copy still gets this BAM's reads added.

## Watch out
- cnv_from_bam installs a Ctrl-C/SIGTERM handler, so `timeout` (SIGTERM) does not stop a run stuck in it; use SIGKILL.
- A failing pytest `==` on two long lists (a 2,001-bin control) spends minutes in difflib building the diff and looks exactly like a hang. Compare something short (the names of contigs that differ) instead.

In the code: `src/robin/analysis/cnv_reference.py` (`reference_cnv_profile`), `src/robin/analysis/cnv_analysis.py` (`_read_only_control`, `get_cached_ref_cnv_dict`, `run_cnv_analysis_direct`, `process_multiple_bams`), `src/robin/analysis/cnv_subprocess.py`. Library source read: cnv_from_bam 0.4.4 `src/lib.rs` (`calculate_cnv`, `merge_and_reset`, `iterate_bam_file`). Tests: `tests/test_cnv_reference.py`, `tests/test_cnv_control_reference.py`, synthetic BAMs from `tests/cnv_fixtures.py`.
