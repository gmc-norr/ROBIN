"""VAR-7: ``robin add-panel`` must keep every region of a gene, not just the first.

A target panel legitimately lists many regions per gene -- one per exon, or per
tile. ``_generate_unique_gene_bed`` deduplicated on the gene name alone, so a
panel with 3 EGFR rows and 2 TP53 rows was written out as 2 rows and coverage,
SNP calling and fusion detection never saw the dropped regions. The dedupe key
is the whole region now, so only an identical region collapses.

The function is module level and ``robin.cli`` imports cleanly, so it is called
directly rather than through the click command. Every BED here is synthetic.
"""

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
