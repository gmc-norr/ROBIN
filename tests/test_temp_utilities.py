import logging

import polars as pl

from robin.analysis.temp_utilities import merge_modkit_files


def _write_modkit_parquet(path, chrom, chrom_start):
    pl.DataFrame(
        {
            "chrom": [chrom],
            "chromStart": [chrom_start],
            "mod_code": ["m"],
            "strand": ["+"],
            "valid_cov": [10],
            "percent_modified": [50.0],
            "n_mod": [5],
            "n_canonical": [5],
        },
        schema={
            "chrom": pl.Utf8,
            "chromStart": pl.UInt32,
            "mod_code": pl.Utf8,
            "strand": pl.Utf8,
            "valid_cov": pl.UInt32,
            "percent_modified": pl.Float32,
            "n_mod": pl.UInt32,
            "n_canonical": pl.UInt32,
        },
    ).write_parquet(path)


def test_file_with_no_panel_sites_is_skipped_without_error(tmp_path, caplog):
    filter_bed = tmp_path / "panel.bed"
    filter_bed.write_text("chr1\t100\t101\tcg0001\n")

    on_panel = tmp_path / "on_panel.parquet"
    off_panel = tmp_path / "off_panel.parquet"
    _write_modkit_parquet(on_panel, "chr1", 100)
    _write_modkit_parquet(off_panel, "chr2", 500)

    output_file = tmp_path / "out.parquet"
    with caplog.at_level(logging.ERROR):
        merge_modkit_files(
            new_files=[str(on_panel), str(off_panel)],
            existing_file=str(output_file),
            output_file=str(output_file),
            filter_bed_file=str(filter_bed),
            sample_id="sample",
            output_dir=str(tmp_path),
            mnpflex_config=None,
            num_bam_files_seen=2,
        )

    assert [r.getMessage() for r in caplog.records if r.levelno >= logging.ERROR] == []
    merged = pl.read_parquet(output_file)
    assert merged.select("chrom", "chromStart").rows() == [("chr1", 100)]
