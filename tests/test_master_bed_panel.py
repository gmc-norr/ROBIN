"""TGT-2 (locked rule 6): the configured panel always reaches target.bam's regions, and
the CNV step asks for a master BED that includes the panel, even before fusion has
written its metadata. Synthetic data only.
"""

import logging
import types
from collections import Counter

import numpy as np
import pysam

from target_fixtures import PANEL, SAMPLE, bam_records, record_key, write_bams, write_panel
from robin.analysis import cnv_analysis
from robin.analysis.target_analysis import TargetAnalysis, finalize_accumulation_for_sample

LOG = logging.getLogger("test_master_bed_panel")


def test_a_master_bed_without_the_panel_still_gives_the_panel_reads(tmp_path, monkeypatch):
    panel = write_panel(str(tmp_path / "panel.bed"))
    monkeypatch.setattr(TargetAnalysis, "_find_target_bed", lambda self, p: panel)
    bams, merged = write_bams(str(tmp_path), 2)
    work = tmp_path / "work"
    bed_dir = work / SAMPLE / "bed_files"
    bed_dir.mkdir(parents=True)
    # What the CNV step wrote without fusion metadata: breakpoints only, no panel.
    (bed_dir / "master_001.bed").write_text("chr2\t100\t300\tcnv_breakpoint\n")
    ta = TargetAnalysis(work_dir=str(work), target_panel="TESTPANEL", batch_size=1, use_staging=True)
    for bam in bams:
        ta.process_file_with_staging(bam, {"sample_id": SAMPLE})
        ta.accumulate_staged_files(SAMPLE, force=True)
    finalize_accumulation_for_sample(SAMPLE, str(work), "TESTPANEL")
    union = write_panel(str(tmp_path / "union.bed"), PANEL + [("chr2", 100, 300, "cnv_breakpoint")])
    want = Counter(record_key(line) for line in pysam.view("-L", union, merged).splitlines())
    assert bam_records(str(work / SAMPLE / "target.bam")) == want


def test_the_cnv_handler_passes_its_panel_on(tmp_path, monkeypatch):
    captured = {}

    def fake_process_multiple_bams(**kwargs):
        captured.update(kwargs)
        return {"files_processed": 1, "total_files": 1}

    monkeypatch.setattr(cnv_analysis, "process_multiple_bams", fake_process_multiple_bams)
    bam = str(tmp_path / "x.bam")
    context = types.SimpleNamespace(metadata={"bam_metadata": {}}, get_sample_id=lambda: SAMPLE)
    batched = types.SimpleNamespace(
        get_file_count=lambda: 1, get_sample_id=lambda: SAMPLE, batch_id="b1",
        get_filepaths=lambda: [bam], contexts=[context],
    )
    job = types.SimpleNamespace(
        job_id=1, job_type="cnv",
        context=types.SimpleNamespace(
            filepath=bam, metadata={"_batched_job": batched}, add_metadata=lambda key, value: None,
        ),
    )
    cnv_analysis.cnv_handler(job, work_dir=str(tmp_path / "work"), target_panel="rCNS2")
    assert captured.get("target_panel") == "rCNS2"


def test_the_cnv_step_asks_for_a_master_bed_with_the_panel(tmp_path, monkeypatch):
    from robin.analysis import master_bed_generator

    captured = {}
    monkeypatch.setattr(master_bed_generator, "generate_master_bed_async", lambda **kw: captured.update(kw))
    sample_dir = tmp_path / SAMPLE
    sample_dir.mkdir()
    signal = {"chr1": np.r_[np.zeros(100), np.ones(100) * 5]}
    cnv_analysis.save_cnv_files(
        str(sample_dir), 1, signal, signal, signal, 1000, 0.1, [], "Unknown", {}, LOG,
        generate_master_bed=True, target_panel="rCNS2",
    )
    assert captured.get("target_panel") == "rCNS2"
