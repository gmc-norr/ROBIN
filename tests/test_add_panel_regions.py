"""VAR-7: ``robin add-panel`` must keep every region of a gene, not just the first.

A target panel legitimately lists many regions per gene -- one per exon, or per
tile. ``_generate_unique_gene_bed`` deduplicated on the gene name alone, so a
panel with 3 EGFR rows and 2 TP53 rows was written out as 2 rows and coverage,
SNP calling and fusion detection never saw the dropped regions. The dedupe key
is the whole region now, so only an identical region collapses.

The function is module level and ``robin.cli`` imports cleanly, so it is called
directly rather than through the click command. Every BED here is synthetic.

The tests from ``_run_add_panel`` down go through the command itself. Columns
after the 4th (score, strand, the rest of BED12) are allowed and ignored: the
validator already accepts them, and the processing step used to shift every
field when there were more than 6.
"""

from click.testing import CliRunner

import robin.cli as cli
import robin.resources
from robin.cli import _generate_unique_gene_bed


def test_every_exon_of_a_gene_survives(tmp_path):
    """VAR-7: 3 EGFR + 2 TP53 rows became 2 rows."""
    src = tmp_path / "panel.bed"
    src.write_text(
        "chr7\t100\t200\tEGFR\n" "chr7\t300\t400\tEGFR\n" "chr7\t500\t600\tEGFR\n"
        "chr17\t100\t200\tTP53\n" "chr17\t300\t400\tTP53\n"
    )
    out = tmp_path / "out.bed"
    assert _generate_unique_gene_bed(src, out) is True
    lines = [l for l in out.read_text().splitlines() if l.strip()]
    assert len(lines) == 5
    assert sum(1 for l in lines if l.endswith("EGFR")) == 3


def test_exactly_repeated_rows_still_collapse(tmp_path):
    """Review Focus 4: keep every distinct region, but do not reintroduce duplicates."""
    src = tmp_path / "panel.bed"
    src.write_text("chr7\t100\t200\tEGFR\n" "chr7\t100\t200\tEGFR\n")
    out = tmp_path / "out.bed"
    _generate_unique_gene_bed(src, out)
    assert len([l for l in out.read_text().splitlines() if l.strip()]) == 1


def test_the_same_gene_on_two_chromosomes_keeps_both(tmp_path):
    """Review Focus 4: a pseudogene or a mis-annotated panel must not lose a locus."""
    src = tmp_path / "panel.bed"
    src.write_text("chr1\t100\t200\tFOO\n" "chr2\t100\t200\tFOO\n")
    out = tmp_path / "out.bed"
    _generate_unique_gene_bed(src, out)
    assert len([l for l in out.read_text().splitlines() if l.strip()]) == 2


def test_a_comma_separated_gene_cell_still_expands(tmp_path):
    src = tmp_path / "panel.bed"
    src.write_text("chr7\t100\t200\tEGFR,ERBB1\n")
    out = tmp_path / "out.bed"
    _generate_unique_gene_bed(src, out)
    genes = sorted(l.split("\t")[3] for l in out.read_text().splitlines() if l.strip())
    assert genes == ["EGFR", "ERBB1"]


# --- through the command -----------------------------------------------------

EXPECTED = ["chr17\t300\t400\tTP53", "chr7\t100\t200\tEGFR"]  # sorted as text
EGFR = "chr7\t100\t200\tEGFR"
TP53 = "chr17\t300\t400\tTP53"
BED12_EGFR = EGFR + "\t0\t+\t100\t200\t0\t1\t100,\t0,"
BED12_TP53 = TP53 + "\t0\t-\t300\t400\t0\t1\t100,\t0,"


def _run_add_panel(tmp_path, monkeypatch, bed_text, name="P"):
    """Run ``robin add-panel``. Panels go to tmp_path/res, never into the package."""
    res = tmp_path / "res"
    res.mkdir()
    monkeypatch.setattr(robin.resources, "__file__", str(res / "__init__.py"))
    monkeypatch.setattr(cli, "_get_user_acknowledgment", lambda: True)
    bed = tmp_path / "upload.bed"
    bed.write_bytes(bed_text.encode())
    cli_result = CliRunner().invoke(cli.main, ["add-panel", str(bed), name])
    out = res / f"{name}_panel_name_uniq.bed"
    lines = out.read_text().splitlines() if out.exists() else None
    return cli_result, lines


def test_the_command_writes_a_4_column_panel(tmp_path, monkeypatch):
    result, lines = _run_add_panel(tmp_path, monkeypatch, f"{EGFR}\n{TP53}\n")
    assert result.exit_code == 0, result.output
    assert lines == EXPECTED


def test_an_8_column_bed_is_read_by_its_first_4_columns(tmp_path, monkeypatch):
    """Used to exit 1 on int('EGFR'): pandas took the 2 extra columns as the index."""
    extra = "\t0\t+\t100\t200"
    result, lines = _run_add_panel(tmp_path, monkeypatch, f"{EGFR}{extra}\n{TP53}{extra}\n")
    assert result.exit_code == 0, result.output
    assert lines == EXPECTED


def test_a_bed12_panel_is_not_silently_garbled(tmp_path, monkeypatch):
    """Used to exit 0, say "added successfully", and save columns 7-10 as the panel."""
    result, lines = _run_add_panel(tmp_path, monkeypatch, f"{BED12_EGFR}\n{BED12_TP53}\n")
    assert result.exit_code == 0, result.output
    assert lines == EXPECTED


def test_rows_with_4_and_12_columns_mix(tmp_path, monkeypatch):
    """The validator accepts this; processing used to stop on a tokenizer error."""
    result, lines = _run_add_panel(tmp_path, monkeypatch, f"{EGFR}\n{BED12_TP53}\n")
    assert result.exit_code == 0, result.output
    assert lines == EXPECTED


def test_mixed_4_and_6_columns_crlf_and_comments_still_work(tmp_path, monkeypatch):
    text = f"# panel v1\r\n{EGFR}\t0\t+\r\n\r\n# note\r\n{TP53}\r\n"
    result, lines = _run_add_panel(tmp_path, monkeypatch, text)
    assert result.exit_code == 0, result.output
    assert lines == EXPECTED


def test_a_failed_panel_leaves_no_files(tmp_path, monkeypatch):
    """A gene cell of only commas passes the validator but names no gene."""
    result, lines = _run_add_panel(tmp_path, monkeypatch, "chr7\t100\t200\t,\n")
    assert result.exit_code == 1, result.output
    assert lines is None
    assert sorted(p.name for p in (tmp_path / "res").iterdir()) == []
