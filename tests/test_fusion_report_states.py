"""REP-6: the fusion report tells "not run", "unreadable" and "none found" apart, and its
methods text matches the filters (batch 7, Task 7).

Rows from the locked plan:
  R1  No fusion output at all -> "Fusion analysis has not been run for this sample."
  R2  A processed fusion pickle that can't be read -> says the fusion results could not be read
  R3  Fusion's preprocessing raises -> the fusion handler files an error; judged failed
  C5  Fusion ran and no read passed the filters (new run: an empty processed pickle; old folder:
      ``sv_count.txt``, no pickle, no candidate rows) -> "No fusion candidates were identified."
  Also: the methods sentence states MAPQ >= 50, mapping span >= 250 bp and >= 3 supporting reads,
  built from ``FUSION_DETECTION_THRESHOLDS``; the filter is checked at the boundaries.

The report section is the real ``FusionSection`` on a temp sample folder; the writer is the real
``accumulate_fusion_candidates`` on synthetic staged rows; the R3 handler path is the real
``fusion_handler`` with only the preprocessing replaced by a function that raises. "Judged failed"
is the shared outcome helper applied to the job's context.

Which tests can fail on the old code by assertion: every report-text test (the old section says
"No fusion candidates were identified." for all three cases and ">50"/">200bp"/"4 or more" in the
methods), the empty-pickle tests (no pickle was written), and R3 (the preprocessing failure was
only a warning). ``test_methods_sentence_follows_the_config`` and the sentence helper tests can only
fail by ImportError there.
"""

import os
import pickle
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pysam
import pytest

from robin import classification_config as cc
from robin import workflow_ray as wr
from robin.analysis import fusion_analysis, fusion_work
from robin.job_outcome import job_outcome

SAMPLE = "S1"

NOT_RUN = "Fusion analysis has not been run for this sample."
UNREADABLE = "The fusion results could not be read."
NONE_FOUND = "No fusion candidates were identified."
NONE_FOUND_SUMMARY = "No fusion networks identified within the target panel."


# ------------------------------------------------------------------ synthetic data
def _rows(n_reads, gene_a="GENEA", gene_b="GENEB", prefix="r", mapq=60, span=400):
    """n_reads reads, each aligned to gene_a (primary) and gene_b (supplementary)."""
    rows = []
    for i in range(n_reads):
        common = dict(
            col2=1000, col3=5000, read_id=f"{prefix}{i}", mapping_quality=mapq, strand="+",
            is_secondary=False, mapping_span=span,
        )
        rows.append(dict(common, col1="chr1", col4=gene_a, reference_id="chr1",
                         reference_start=2000 + i, reference_end=2000 + i + span,
                         read_start=0, read_end=span, is_supplementary=False))
        rows.append(dict(common, col1="chr2", col4=gene_b, reference_id="chr2",
                         reference_start=3000 + i, reference_end=3000 + i + span,
                         read_start=span, read_end=2 * span, is_supplementary=True))
    return pd.DataFrame(rows)


def _sample_dir(tmp_path):
    d = tmp_path / SAMPLE
    d.mkdir()
    return d


def _report_text(sample_dir):
    from robin.reporting.sections.fusion import FusionSection
    from robin.reporting.styling.styles import ReportStyles

    report = SimpleNamespace(
        output=str(sample_dir), styles=ReportStyles(tempfile.mkdtemp())
    )
    section = FusionSection(report)
    section.add_content()
    body = "\n".join(
        el.getPlainText() for el in section.elements if hasattr(el, "getPlainText")
    )
    summary = "\n".join(
        el.getPlainText() for el in section.summary_elements if hasattr(el, "getPlainText")
    )
    return body, summary


def _write_pickles(sample_dir, df, names=("master", "all")):
    for n in names:
        fusion_work.preprocess_fusion_data_standalone(
            df, str(sample_dir / f"fusion_candidates_{n}_processed.pkl")
        )


def _stage(work_dir, df, kinds=("target", "genome"), n=1):
    staging = Path(fusion_work._get_staging_dir(str(work_dir), SAMPLE))
    for kind in kinds:
        df.to_parquet(staging / f"{kind}_{n}.parquet", index=False)


def _accumulate(work_dir):
    return fusion_work.accumulate_fusion_candidates(
        str(work_dir), SAMPLE, "rCNS2", force=True, batch_size=1
    )


def _load(path):
    with open(path, "rb") as f:
        return pickle.load(f)


# ------------------------------------------------------------------ R1, R2
def test_r1_no_output_says_not_run(tmp_path):
    body, summary = _report_text(_sample_dir(tmp_path))
    assert NOT_RUN in body
    assert NOT_RUN in summary
    assert NONE_FOUND not in body
    assert NONE_FOUND_SUMMARY not in summary


@pytest.mark.parametrize("content", [b"", b"not a pickle at all", b"\x80\x04\x95"])
def test_r2_unreadable_pickle_says_could_not_be_read(tmp_path, content):
    d = _sample_dir(tmp_path)
    (d / "fusion_candidates_master_processed.pkl").write_bytes(content)
    body, summary = _report_text(d)
    assert UNREADABLE in body
    assert UNREADABLE in summary
    assert NONE_FOUND not in body
    assert NOT_RUN not in body
    assert NONE_FOUND_SUMMARY not in summary


def test_r2_unreadable_beats_sv_count(tmp_path):
    d = _sample_dir(tmp_path)
    (d / "sv_count.txt").write_text("0")
    (d / "fusion_candidates_all_processed.pkl").write_bytes(b"garbage")
    body, _ = _report_text(d)
    assert UNREADABLE in body
    assert NONE_FOUND not in body


# ------------------------------------------------------------------ C5
def test_c5_new_run_where_nothing_passes_writes_empty_pickles(tmp_path):
    # Two reads for the gene pair: below the 3-read support.
    _stage(tmp_path, _rows(2))
    result = _accumulate(tmp_path)
    assert result["status"] == "success"
    d = tmp_path / SAMPLE
    for name in ("master", "all"):
        pkl = d / f"fusion_candidates_{name}_processed.pkl"
        assert pkl.exists(), f"no processed pickle for {name}"
        data = _load(pkl)
        assert data["annotated_data"].empty
        assert len(data["goodpairs"]) == 0
        assert data["gene_pairs"] == [] and data["gene_groups"] == []
        assert data["candidate_count"] == 0
    assert not list(d.glob("*.tmp")), "temp files left behind"
    body, summary = _report_text(d)
    assert NONE_FOUND in body
    assert NONE_FOUND_SUMMARY in summary
    assert NOT_RUN not in body and UNREADABLE not in body


def test_c5_empty_pickle_replaces_a_stale_populated_one(tmp_path):
    d = _sample_dir(tmp_path)
    _write_pickles(d, _rows(4, prefix="old"))
    assert not _load(d / "fusion_candidates_master_processed.pkl")["annotated_data"].empty
    _stage(tmp_path, _rows(2))
    _accumulate(tmp_path)
    for name in ("master", "all"):
        data = _load(d / f"fusion_candidates_{name}_processed.pkl")
        assert data["annotated_data"].empty
    body, _ = _report_text(d)
    assert NONE_FOUND in body


def test_c5_old_folder_with_sv_count_and_no_pickle_says_none_found(tmp_path):
    d = _sample_dir(tmp_path)
    (d / "sv_count.txt").write_text("0")
    body, summary = _report_text(d)
    assert NONE_FOUND in body
    assert NONE_FOUND_SUMMARY in summary
    assert NOT_RUN not in body


def test_c5_old_folder_with_header_only_csv_says_none_found(tmp_path):
    d = _sample_dir(tmp_path)
    (d / "sv_count.txt").write_text("0")
    _rows(1).iloc[0:0].to_csv(d / "fusion_candidates_master.csv", index=False)
    body, _ = _report_text(d)
    assert NONE_FOUND in body
    assert UNREADABLE not in body


def test_old_folder_with_candidate_rows_but_no_pickle_is_unreadable(tmp_path):
    d = _sample_dir(tmp_path)
    (d / "sv_count.txt").write_text("1")
    _rows(4).to_csv(d / "fusion_candidates_master.csv", index=False)
    body, summary = _report_text(d)
    assert UNREADABLE in body
    assert UNREADABLE in summary
    assert NONE_FOUND not in body


# ------------------------------------------------------------------ control: candidates
def test_candidates_are_reported_with_the_methods_sentence(tmp_path):
    d = _sample_dir(tmp_path)
    _write_pickles(d, _rows(4))
    body, summary = _report_text(d)
    assert "GENEA" in summary and "GENEB" in summary
    assert NONE_FOUND not in body and NOT_RUN not in body and UNREADABLE not in body
    assert fusion_report_sentence() in body
    assert "(>50)" not in body and "(>200bp)" not in body


# ------------------------------------------------------------------ methods sentence
def fusion_report_sentence():
    from robin.reporting.sections.fusion import fusion_methods_sentence

    return fusion_methods_sentence()


def test_methods_sentence_states_the_locked_values():
    s = fusion_report_sentence()
    assert "mapping quality of at least 50" in s
    assert "mapping span of at least 250 bp" in s
    assert "at least 3 supporting reads" in s


def test_methods_sentence_follows_the_config(monkeypatch):
    t = cc.FUSION_DETECTION_THRESHOLDS
    monkeypatch.setitem(t, "mapping_quality", dict(t["mapping_quality"], min_threshold=59))
    monkeypatch.setitem(t, "mapping_span", dict(t["mapping_span"], min_threshold=399))
    monkeypatch.setitem(t, "read_support", dict(t["read_support"], min_threshold=5))
    s = fusion_report_sentence()
    assert "mapping quality of at least 60" in s
    assert "mapping span of at least 400 bp" in s
    assert "at least 5 supporting reads" in s


# ------------------------------------------------------------------ the filters at the boundaries
GENES = {
    "chr1": [fusion_work.GeneRegion(0, 100000, "GENEA")],
    "chr2": [fusion_work.GeneRegion(0, 100000, "GENEB")],
}


def _bam(tmp_path, mapq, span):
    path = str(tmp_path / "t.bam")
    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 200000}, {"SN": "chr2", "LN": 200000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for flag, ref, pos, other_ref, other_pos in (
            (0, 0, 2000, "chr2", 3000),
            (2048, 1, 3000, "chr1", 2000),
        ):
            a = pysam.AlignedSegment(out.header)
            a.query_name = "read1"
            a.query_sequence = "A" * span
            a.flag = flag
            a.reference_id = ref
            a.reference_start = pos
            a.mapping_quality = mapq
            a.cigarstring = f"{span}M"
            a.set_tag("SA", f"{other_ref},{other_pos + 1},+,{span}M,{mapq},0;")
            out.write(a)
    return path


@pytest.mark.parametrize(
    "mapq,span,kept",
    [(49, 400, False), (50, 400, True), (60, 249, False), (60, 250, True)],
)
def test_filter_boundaries_match_the_stated_values(tmp_path, mapq, span, kept):
    df = fusion_work._process_reads_for_fusions(_bam(tmp_path, mapq, span), {"read1"}, GENES)
    assert (df is not None and not df.empty) is kept


def test_master_bed_filter_boundaries_match_the_stated_values(tmp_path):
    for mapq, span, kept in [(49, 400, False), (50, 400, True), (60, 249, False), (60, 250, True)]:
        df = fusion_work._process_reads_for_master_bed_fusions(
            _bam(tmp_path, mapq, span), {"read1"}
        )
        assert (df is not None and not df.empty) is kept, (mapq, span)


@pytest.mark.parametrize("n_reads,good", [(2, False), (3, True)])
def test_read_support_boundary(n_reads, good):
    _, goodpairs = fusion_work._annotate_results(_rows(n_reads))
    assert bool(goodpairs.any()) is good


# ------------------------------------------------------------------ R3 through the real handler
def _batched_job(tmp_path):
    bam = str(tmp_path / "in_0.bam")
    ctx = wr.WorkflowContext(
        filepath=bam,
        metadata={
            "bam_metadata": {"sample_id": SAMPLE, "has_supplementary_reads": True},
            "target_panel": "rCNS2",
        },
    )
    ctx.metadata["_batched_job"] = wr.BatchedJob(
        job_id=1, job_type="fusion", origin="fusion", workflow=["fusion:fusion"], step=0,
        contexts=[ctx], batch_id="b1", sample_id=SAMPLE,
    )
    return wr.Job(1, "fusion", "fusion", ["fusion:fusion"], 0, ctx)


def _run_handler(tmp_path, monkeypatch, preprocess=None):
    def staged(bam_path, temp_dir, metadata, fusion_metadata, target_panel, **kw):
        _stage(kw["work_dir"], _rows(4), kinds=("target", "genome"))
        return {}, False

    monkeypatch.setattr(fusion_analysis, "process_bam_with_staging", staged)
    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    if preprocess is not None:
        monkeypatch.setattr(fusion_work, "preprocess_fusion_data_standalone", preprocess)
    job = _batched_job(tmp_path)
    before = len(job.context.errors)
    fusion_analysis.fusion_handler(job, work_dir=str(tmp_path / "work"), target_panel="rCNS2")
    return job, job.context.errors[before:], job_outcome("fusion", job.context, before)


def test_r3_a_raising_preprocessing_is_filed_and_the_job_is_judged_failed(tmp_path, monkeypatch):
    def boom(df, output_file):
        raise RuntimeError("synthetic preprocessing crash")

    _, filed, outcome = _run_handler(tmp_path, monkeypatch, preprocess=boom)
    assert [e["job_type"] for e in filed] == ["fusion_analysis"]
    assert "synthetic preprocessing crash" in filed[0]["error"]
    assert outcome["status"] == "failed"


def test_r3_control_a_working_preprocessing_files_nothing(tmp_path, monkeypatch):
    job, filed, outcome = _run_handler(tmp_path, monkeypatch)
    assert filed == []
    assert outcome["status"] == "success"
    assert (tmp_path / "work" / SAMPLE / "fusion_candidates_master_processed.pkl").exists()


# ------------------------------------------------------------------ a failed build (Ruling J)
def _boom(df, output_file):
    raise RuntimeError("synthetic preprocessing crash")


def test_a_raising_preprocessing_still_finishes_the_steps_other_outputs(tmp_path, monkeypatch):
    _, filed, outcome = _run_handler(tmp_path, monkeypatch, preprocess=_boom)
    assert outcome["status"] == "failed" and filed
    d = tmp_path / "work" / SAMPLE
    assert (d / "sv_count.txt").exists()
    assert list((d / "bed_files").glob("fusion_breakpoints_*.bed")), "breakpoint BED not written"
    assert not list((d / "_fusion_staging").glob("*.parquet")), "staging not cleaned"
    assert not (d / "fusion_candidates_master_processed.pkl").exists()
    # The report does not call it "none found": results exist but were not processed.
    body, _ = _report_text(d)
    assert UNREADABLE in body


def test_a_failed_build_is_rebuilt_on_the_next_run_without_double_counting(tmp_path, monkeypatch):
    work = tmp_path / "work"
    d = work / SAMPLE
    monkeypatch.setattr(fusion_work, "_ensure_gene_regions_loaded", lambda panel: None)
    _stage(work, _rows(4))
    real = fusion_work.preprocess_fusion_data_standalone
    monkeypatch.setattr(fusion_work, "preprocess_fusion_data_standalone", _boom)
    assert _accumulate(work)["status"] == "error"
    assert not (d / "fusion_candidates_master_processed.pkl").exists()
    csv_rows = len(pd.read_csv(d / "fusion_candidates_master.csv"))
    assert csv_rows == 8

    # Second run: preprocessing works again; the only new rows are one read that does not pass.
    monkeypatch.setattr(fusion_work, "preprocess_fusion_data_standalone", real)
    _stage(work, _rows(1, "GENEC", "GENED", prefix="x"), n=2)
    assert _accumulate(work)["status"] == "success"
    for name in ("master", "all"):
        data = _load(d / f"fusion_candidates_{name}_processed.pkl")
        assert not data["annotated_data"].empty
        assert data["annotated_data"]["read_id"].nunique() == 4
    assert len(pd.read_csv(d / "fusion_candidates_master.csv")) == csv_rows
    body, summary = _report_text(d)
    assert "GENEA" in summary


def test_a_step_inside_the_real_preprocessing_that_raises_fails_the_job(tmp_path, monkeypatch):
    # Pins preprocess_fusion_data_standalone's own re-raise: the heavy step it calls raises,
    # the function itself is real.
    def boom(df):
        raise RuntimeError("synthetic annotate crash")

    monkeypatch.setattr(fusion_work, "_annotate_results", boom)
    _, filed, outcome = _run_handler(tmp_path, monkeypatch)
    assert [e["job_type"] for e in filed] == ["fusion_analysis"]
    assert "synthetic annotate crash" in filed[0]["error"]
    assert outcome["status"] == "failed"
    body, _ = _report_text(tmp_path / "work" / SAMPLE)
    assert UNREADABLE in body
