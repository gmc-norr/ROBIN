---
title: Reading a BED with pandas and fixed column names shifts every field when the file is wider
author: Pär
date: 2026-10-31
branch: fix/batch-3b-should-fix-soon
cost: { minutes: 18, tokens: 184333 }
id: l-a5b8a814
scope: repo
status: active
src: agent
touches:
  - "src/robin/cli.py#_generate_unique_gene_bed"
  - tests/test_add_panel_regions.py
---

## Decided
- robin add-panel reads the uploaded BED line by line, the same way its validator does (strip, skip blank and # lines, split on tabs), and keeps chrom, start, end and gene. Columns after the 4th are allowed and ignored. Parsing the same way as the validator means the two cannot disagree about a file.

## Tried and ruled out
- pd.read_csv(header=None, names=[six names]) on a wider file. When names is shorter than the row, pandas silently uses the extra leading columns as the index, so every field shifts left. Measured: 7 or 8 columns fail on int('EGFR') right after "validation passed"; BED12 succeeds, prints "added successfully" and saves thickStart, thickEnd, itemRgb and blockCount as chrom, start, end and gene. A file mixing 4- and 12-column rows fails on the tokenizer, because pandas sizes rows from the first line.

## Watch out
- Other read_csv(..., names=[...]) calls exist in src. master_bed_generator._load_bed_file counts columns before choosing names, so it is safe; the CNV GUI and report read the 4-column file add-panel writes. Check any new one against a wider file.
- Test add-panel through the command with click's CliRunner, pointing robin.resources.__file__ at tmp_path and replacing _get_user_acknowledgment; otherwise it writes the panel into the installed package.
- Custom panels added from a BED wider than 6 columns before this fix are wrong on disk and need re-adding.

In the code: _generate_unique_gene_bed and _validate_bed_file in src/robin/cli.py; _run_add_panel in tests/test_add_panel_regions.py.
