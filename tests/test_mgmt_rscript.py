"""SEC-8: the MGMT R script runs without a shell."""

import logging
import os
import sys
import types
from pathlib import Path

import pytest

from robin.analysis import mgmt_analysis
from robin.analysis.mgmt_analysis import run_matkit

# One BEDMethyl line: convert_to_mixed_delim needs 18 columns.
BEDMETHYL_LINE = "chr10\t1\t2\tm\t0\t+\t1\t2\t0,0,0\t10\t50.0\t5\t5\t0\t0\t0\t0\t0\n"

# Records each argument on its own line, and writes the status CSV the real
# R script writes into --out_dir.
FAKE_RSCRIPT = """#!/bin/sh
: > "$ARGS_FILE"
for a in "$@"; do
  printf '%s\\n' "$a" >> "$ARGS_FILE"
  case "$a" in --out_dir=*) printf 'status\\n' > "${a#--out_dir=}/live_analysis_mgmt_status.csv";; esac
done
"""


@pytest.fixture
def mgmt_inputs(tmp_path, monkeypatch):
    """Fake the BAM sort and matkit, and lay out hv_path with the R script."""

    def fake_sort_and_index(input_bam, output_bam, **kwargs):
        Path(output_bam).touch()
        return True

    def fake_matkit(in_bam, out_bed):
        Path(out_bed).write_text(BEDMETHYL_LINE)

    monkeypatch.setattr(mgmt_analysis, "safely_sort_and_index_bam", fake_sort_and_index)
    matkit = types.ModuleType("robin.analysis.utilities.matkit")
    matkit.run_matkit = fake_matkit
    monkeypatch.setitem(sys.modules, "robin.analysis.utilities.matkit", matkit)

    hv = tmp_path / "hv"
    (hv / "bin").mkdir(parents=True)
    (hv / "bin" / "mgmt_pred_v0.3.R").write_text("# not run\n")
    # A shell injection would run in the current folder; keep that folder in tmp.
    monkeypatch.chdir(tmp_path)
    return hv


def _put_fake_rscript_on_path(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    rscript = bin_dir / "Rscript"
    rscript.write_text(FAKE_RSCRIPT)
    rscript.chmod(0o755)
    args_file = tmp_path / "rscript_args.txt"
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("ARGS_FILE", str(args_file))
    return args_file


def test_sample_dir_with_shell_syntax_is_one_literal_argument(tmp_path, monkeypatch, mgmt_inputs):
    args_file = _put_fake_rscript_on_path(tmp_path, monkeypatch)
    sample_dir = tmp_path / "s $(touch MARKER) x"
    sample_dir.mkdir()

    assert run_matkit(str(tmp_path), "in.bam", str(mgmt_inputs), str(sample_dir), 1) is True

    args = args_file.read_text().splitlines()
    assert f"--out_dir={sample_dir}" in args
    assert not (tmp_path / "MARKER").exists()
    # The R script's status file reached the sample folder and was renamed.
    assert (sample_dir / "1_mgmt.csv").exists()


def test_missing_rscript_logs_warning_and_continues(tmp_path, monkeypatch, mgmt_inputs, caplog):
    empty_bin = tmp_path / "empty_bin"
    empty_bin.mkdir()
    monkeypatch.setenv("PATH", str(empty_bin))
    sample_dir = tmp_path / "sample"
    sample_dir.mkdir()

    with caplog.at_level(logging.WARNING, logger="robin.mgmt"):
        assert run_matkit(str(tmp_path), "in.bam", str(mgmt_inputs), str(sample_dir), 1) is True
    assert "R script failed" in caplog.text
    # The BED file is still kept when R is missing.
    assert (sample_dir / "1_mgmt.bed").exists()
