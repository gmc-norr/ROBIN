---
title: "The ClinVar site boundary rule lives in the BED-intersection step, and writers must keep the umask"
author: Pär
date: 2026-10-31
branch: fix/batch-3c-review-findings
cost: { minutes: 11, tokens: 74796 }
id: l-75436380
scope: repo
status: active
src: agent
touches:
  - "src/robin/analysis/lightweight_gene_analysis.py#_find_bed_clinvar_intersections"
  - "src/robin/utils/atomic_files.py#atomic_open"
---

## Decided
- Correction to the 2026-10-31 learning on 1-based ClinVar sites: the start < POS <= end rule is in _find_bed_clinvar_intersections (where sites are picked for BED regions), not _find_intersecting_variants (which matches picked sites to pileup columns and adds the no_reads / pileup_failed records).
- Result files are written with robin.utils.atomic_files.atomic_open: a temp file made with os.open(O_CREAT|O_EXCL, 0o666) so the umask applies exactly as with open(), then os.replace.

## Tried and ruled out
- tempfile.mkstemp for atomic writes: it creates 0600 regardless of the umask, and os.replace carries that onto the result, so group members lose read access to results in shared folders.

## Watch out
- _perform_targeted_pileup returns None when it could not run (no index, contig missing, error) and [] when it ran and saw no reads; keep the two apart, they label sites differently.
