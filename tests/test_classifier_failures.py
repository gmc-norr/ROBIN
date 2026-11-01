"""Sturgeon, NanoDX and MGMT report their failures (batch 7, METH-8).

Rows from the locked plan:
  M1  Sturgeon's classifier raises inside its background step -> the Sturgeon handler files an
      error; the job is judged failed.
  M2  NanoDX's network raises -> judged failed; the store CSV gains no row; no ``errordf.csv``
      anywhere.
  M3  MGMT's R step: non-zero exit, ``Rscript`` missing, script missing, on the single and the
      batch path -> the MGMT handler files an error; judged failed; the MGMT BED is still written.
  C6  Successful runs -> NanoDX's row has its scores; MGMT's CSV is written; no error filed.
  Control: "No MGMT spanning reads found" stays an expected condition (not failed).

The real handlers run on synthetic inputs. Only the heavy external step is replaced: the Sturgeon
model, the NanoDX network, and ``subprocess.run`` for bedtools / Rscript. "Judged failed" is the
shared outcome helper (``robin.job_outcome.job_outcome``) applied to the job's context.
"""

import importlib
import subprocess
import sys
import types
from pathlib import Path

import pandas as pd
import pysam
import pytest

from robin import workflow_ray as wr
from robin.analysis import mgmt_analysis
from robin.job_outcome import job_outcome

SAMPLE = "S1"
BEDMETHYL_LINE = "chr10\t1\t2\tm\t0\t+\t1\t2\t0,0,0\t10\t50.0\t5\t5\t0\t0\t0\t0\t0\n"


# ---------------------------------------------------------------- shared helpers
def _job(job_type, contexts, batched):
    """A Job whose context is the first of ``contexts``; batched wraps them in a BatchedJob."""
    ctx = contexts[0]
    if batched:
        ctx.metadata["_batched_job"] = wr.BatchedJob(
            job_id=1, job_type=job_type, origin=job_type, workflow=[f"{job_type}:{job_type}"],
            step=0, contexts=contexts, batch_id="b1", sample_id=SAMPLE,
        )
    return wr.Job(1, job_type, job_type, [f"{job_type}:{job_type}"], 0, ctx)


def _run(handler, job, **kwargs):
    """Run the real handler; return (errors it filed, its outcome as the job runner judges it)."""
    before = len(job.context.errors)
    handler(job, **kwargs)
    filed = job.context.errors[before:]
    return filed, job_outcome(job.job_type, job.context, before)


def _parquet_contexts(tmp_path, n):
    contexts = []
    for i in range(n):
        pq = tmp_path / f"in_{i}.parquet"
        pq.write_bytes(b"x")
        contexts.append(
            wr.WorkflowContext(
                filepath=str(tmp_path / f"in_{i}.bam"),
                metadata={"bam_metadata": {"sample_id": SAMPLE}},
                results={"bed_conversion": {"parquet_path": str(pq)}},
            )
        )
    return contexts


# ---------------------------------------------------------------- Sturgeon
@pytest.fixture
def sturgeon_mod(monkeypatch):
    """robin.analysis.sturgeon_analysis imported against stub ``sturgeon`` / ``onnxruntime``
    packages (neither is installed here); only the Sturgeon model call is then patched."""
    stubs = {
        "sturgeon": types.ModuleType("sturgeon"),
        "sturgeon.prediction": types.ModuleType("sturgeon.prediction"),
        "sturgeon.utils": types.ModuleType("sturgeon.utils"),
        "sturgeon.constants": types.ModuleType("sturgeon.constants"),
        "onnxruntime": types.ModuleType("onnxruntime"),
    }
    stubs["sturgeon"].__file__ = "/nonexistent/sturgeon/__init__.py"
    stubs["sturgeon.prediction"].bed_to_numpy = None
    stubs["sturgeon.utils"].load_bed_file = None
    stubs["sturgeon.utils"].softmax = None
    stubs["sturgeon.utils"].read_probes_file = None
    stubs["sturgeon.constants"].NOMEASURE_VALUE = -1
    stubs["onnxruntime"].InferenceSession = object
    for name, mod in stubs.items():
        monkeypatch.setitem(sys.modules, name, mod)
    names = ("robin.analysis.utilities.merge_bedmethyl", "robin.analysis.sturgeon_analysis")
    for name in names:
        monkeypatch.delitem(sys.modules, name, raising=False)
    sa = importlib.import_module("robin.analysis.sturgeon_analysis")
    monkeypatch.setattr(sa, "load_modkit_data", lambda path: pd.DataFrame({"x": [1]}))
    monkeypatch.setattr(sa, "modkit_pileup_file_to_bed", lambda df, out, probes: df)
    yield sa
    for name in names:  # the stub-bound modules must not leak into other tests
        sys.modules.pop(name, None)


@pytest.mark.parametrize("batched", [False, True], ids=["single", "batch"])
def test_m1_sturgeon_classifier_that_raises_is_filed_and_judged_failed(
    sturgeon_mod, tmp_path, monkeypatch, batched
):
    def boom(df):
        raise RuntimeError("sturgeon model blew up")

    monkeypatch.setattr(sturgeon_mod, "predict_sample_from_dataframe", boom)
    job = _job("sturgeon", _parquet_contexts(tmp_path, 2 if batched else 1), batched)

    filed, outcome = _run(sturgeon_mod.sturgeon_handler, job, work_dir=str(tmp_path / "w"))

    assert [e["job_type"] for e in filed] == ["sturgeon_analysis"]
    assert outcome["status"] == "failed"
    if not batched:  # the single path carries the model's own message
        assert "sturgeon model blew up" in filed[0]["error"]
    assert not (tmp_path / "w" / SAMPLE / "sturgeon_scores.csv").exists()


@pytest.mark.parametrize("batched", [False, True], ids=["single", "batch"])
def test_c6_successful_sturgeon_writes_scores_and_files_no_error(
    sturgeon_mod, tmp_path, monkeypatch, batched
):
    monkeypatch.setattr(
        sturgeon_mod, "predict_sample_from_dataframe",
        lambda df: pd.DataFrame({"glioma": [0.9], "other": [0.1]}),
    )
    job = _job("sturgeon", _parquet_contexts(tmp_path, 2 if batched else 1), batched)

    filed, outcome = _run(sturgeon_mod.sturgeon_handler, job, work_dir=str(tmp_path / "w"))

    assert filed == []
    assert outcome["status"] == "success"
    scores = pd.read_csv(tmp_path / "w" / SAMPLE / "sturgeon_scores.csv", index_col=0)
    assert list(scores.columns) == ["glioma", "other"]
    assert len(scores) == (2 if batched else 1)


# ---------------------------------------------------------------- NanoDX / PanNanoDX
class _FakeNet:
    """Stands in for the neural network. ``mode``: ok | raise."""

    mode = "ok"

    def __init__(self, modelfile):
        pass

    def predict(self, df):
        if _FakeNet.mode == "raise":
            raise RuntimeError("network blew up")
        return [0.8, 0.2], ["glioma", "other"], len(df)


@pytest.fixture
def nanodx_mod(monkeypatch, tmp_path):
    from robin.analysis import nanodx_analysis as nd

    calls = pd.DataFrame({"chrom": ["chr1", "chr1", "chr2"], "start_pos": [1, 2, 3],
                          "fraction": [90.0, 10.0, 80.0]})
    monkeypatch.setattr(nd, "NN_classifier", _FakeNet, raising=False)
    monkeypatch.setattr(nd, "load_modkit_data", lambda path: calls, raising=False)
    monkeypatch.setattr(nd, "collapse_minimal_bedmethyl", lambda df: df, raising=False)
    cpgs = pd.DataFrame({0: ["chr1", "chr1", "chr2"], 1: [1, 2, 3], 2: [2, 3, 4],
                         3: ["cg1", "cg2", "cg3"]})
    monkeypatch.setattr(nd.NanodxAnalysis, "_load_cpgs_data", lambda self: cpgs)
    _FakeNet.mode = "ok"
    # A write to the current directory (the old errordf.csv) must be visible, and stay in tmp.
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    yield nd
    _FakeNet.mode = "ok"


NANODX_CASES = [("nanodx", "nanodx_handler", "NanoDX_scores.csv"),
                ("pannanodx", "pannanodx_handler", "PanNanoDX_scores.csv")]


@pytest.mark.parametrize("batched", [False, True], ids=["single", "batch"])
@pytest.mark.parametrize("job_type,handler_name,store", NANODX_CASES)
def test_m2_nanodx_network_that_raises_is_filed_with_no_row_and_no_errordf(
    nanodx_mod, tmp_path, job_type, handler_name, store, batched
):
    _FakeNet.mode = "raise"
    work = tmp_path / "w"
    store_path = work / SAMPLE / store
    store_path.parent.mkdir(parents=True)
    # An earlier good run left one row; a failed run must not add another.
    pd.DataFrame({"glioma": [0.7], "other": [0.3], "number_probes": [3]},
                 index=pd.Index([111.0], name="timestamp")).to_csv(store_path)
    job = _job(job_type, _parquet_contexts(tmp_path, 2 if batched else 1), batched)

    filed, outcome = _run(getattr(nanodx_mod, handler_name), job, work_dir=str(work))

    assert [e["job_type"] for e in filed] == [f"{job_type}_analysis"]
    assert outcome["status"] == "failed"
    if not batched:
        assert "network blew up" in filed[0]["error"]
    assert len(pd.read_csv(store_path, index_col=0)) == 1
    assert list(tmp_path.rglob("errordf.csv")) == []


def test_m2_nanodx_without_its_network_is_filed_not_an_empty_row(nanodx_mod, tmp_path, monkeypatch):
    monkeypatch.setattr(nanodx_mod, "NN_classifier", None)
    job = _job("nanodx", _parquet_contexts(tmp_path, 1), False)

    filed, outcome = _run(nanodx_mod.nanodx_handler, job, work_dir=str(tmp_path / "w"))

    assert [e["job_type"] for e in filed] == ["nanodx_analysis"]
    assert outcome["status"] == "failed"
    assert not (tmp_path / "w" / SAMPLE / "NanoDX_scores.csv").exists()
    assert list(tmp_path.rglob("errordf.csv")) == []


@pytest.mark.parametrize("batched", [False, True], ids=["single", "batch"])
@pytest.mark.parametrize("job_type,handler_name,store", NANODX_CASES)
def test_c6_successful_nanodx_row_has_its_scores(
    nanodx_mod, tmp_path, job_type, handler_name, store, batched
):
    job = _job(job_type, _parquet_contexts(tmp_path, 2 if batched else 1), batched)

    filed, outcome = _run(getattr(nanodx_mod, handler_name), job, work_dir=str(tmp_path / "w"))

    assert filed == []
    assert outcome["status"] == "success"
    rows = pd.read_csv(tmp_path / "w" / SAMPLE / store, index_col=0)
    assert len(rows) == (2 if batched else 1)
    assert list(rows.columns) == ["glioma", "other", "number_probes"]
    assert (rows["glioma"] == 0.8).all() and (rows["number_probes"] == 3).all()
    assert list(tmp_path.rglob("errordf.csv")) == []


# ---------------------------------------------------------------- MGMT
def _mgmt_bam(path):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr10", "LN": 135_000_000}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        a = pysam.AlignedSegment()
        a.query_name = "r1"
        a.query_sequence = "ACGTACGTAC"
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 129_467_000
        a.mapping_quality = 60
        a.cigar = [(0, 10)]
        a.query_qualities = pysam.qualitystring_to_array("IIIIIIIIII")
        out.write(a)


@pytest.fixture
def mgmt_env(tmp_path, monkeypatch):
    """MGMT with only the external steps replaced: bedtools and Rscript (``subprocess.run``),
    the BAM sort, matkit, the plot, and where the R script lives. ``env.r`` picks what R does:
    ok | exit | missing."""
    env = types.SimpleNamespace(r="ok", rscript_calls=0)
    src = tmp_path / "in.bam"
    _mgmt_bam(src)
    env.bam = str(src)

    def fake_run(cmd, *args, **kwargs):
        if cmd[0] == "bedtools":
            kwargs["stdout"].write(Path(cmd[3]).read_bytes())
            return subprocess.CompletedProcess(cmd, 0, b"", b"")
        assert cmd[0] == "Rscript", cmd
        env.rscript_calls += 1
        if env.r == "missing":
            raise FileNotFoundError(2, "No such file or directory: 'Rscript'")
        if env.r == "exit":
            return subprocess.CompletedProcess(cmd, 1, "", "R: model could not be fitted")
        out_dir = next(a.split("=", 1)[1] for a in cmd if a.startswith("--out_dir="))
        Path(out_dir, "live_analysis_mgmt_status.csv").write_text("status\nunmethylated\n")
        return subprocess.CompletedProcess(cmd, 0, "", "")

    def fake_sort_and_index(input_bam, output_bam, **kwargs):
        Path(output_bam).touch()
        return True

    def fake_matkit(in_bam, out_bed):
        Path(out_bed).write_text(BEDMETHYL_LINE)

    hv = tmp_path / "hv"
    (hv / "bin").mkdir(parents=True)
    (hv / "bin" / "mgmt_pred_v0.3.R").write_text("# not run\n")
    env.hv = hv

    monkeypatch.setattr(mgmt_analysis.subprocess, "run", fake_run)
    monkeypatch.setattr(mgmt_analysis, "safely_sort_and_index_bam", fake_sort_and_index)
    monkeypatch.setattr(mgmt_analysis, "_find_hv_path", lambda work_dir: str(hv))
    monkeypatch.setattr(
        mgmt_analysis, "run_methylartist_safely",
        lambda **kw: {"success": False, "fallback_used": False, "error_message": "no plot here"},
    )
    matkit = types.ModuleType("robin.analysis.utilities.matkit")
    matkit.run_matkit = fake_matkit
    monkeypatch.setitem(sys.modules, "robin.analysis.utilities.matkit", matkit)
    monkeypatch.chdir(tmp_path)
    return env


def _mgmt_job(env, batched, has_reads=True):
    meta = {"sample_id": SAMPLE, "has_mgmt_reads": has_reads, "mgmt_read_count": 1}
    contexts = [wr.WorkflowContext(filepath=env.bam, metadata={"bam_metadata": dict(meta)})]
    return _job("mgmt", contexts, batched)


@pytest.mark.parametrize("batched", [False, True], ids=["single", "batch"])
@pytest.mark.parametrize("failure", ["exit", "missing", "no_script"])
def test_m3_mgmt_r_failure_is_filed_judged_failed_and_keeps_the_bed(
    mgmt_env, tmp_path, failure, batched
):
    if failure == "no_script":
        (mgmt_env.hv / "bin" / "mgmt_pred_v0.3.R").unlink()
    else:
        mgmt_env.r = failure
    job = _mgmt_job(mgmt_env, batched)

    filed, outcome = _run(mgmt_analysis.mgmt_handler, job, work_dir=str(tmp_path / "w"))

    assert [e["job_type"] for e in filed] == ["mgmt_analysis"]
    assert outcome["status"] == "failed"
    sample_dir = tmp_path / "w" / SAMPLE
    bed = sample_dir / ("final_mgmt.bed" if batched else "1_mgmt.bed")
    assert bed.exists(), sorted(p.name for p in sample_dir.iterdir())
    assert not list(sample_dir.glob("*_mgmt.csv"))
    assert mgmt_env.rscript_calls == (0 if failure == "no_script" else 1)


@pytest.mark.parametrize("batched", [False, True], ids=["single", "batch"])
def test_c6_successful_mgmt_writes_its_csv_and_files_no_error(mgmt_env, tmp_path, batched):
    job = _mgmt_job(mgmt_env, batched)

    filed, outcome = _run(mgmt_analysis.mgmt_handler, job, work_dir=str(tmp_path / "w"))

    assert filed == []
    assert outcome["status"] == "success"
    sample_dir = tmp_path / "w" / SAMPLE
    prefix = "final" if batched else "1"
    assert (sample_dir / f"{prefix}_mgmt.csv").exists()
    assert (sample_dir / f"{prefix}_mgmt.bed").exists()


@pytest.mark.parametrize("batched", [False, True], ids=["single", "batch"])
def test_no_mgmt_reads_is_still_an_expected_condition_not_a_failure(mgmt_env, tmp_path, batched):
    job = _mgmt_job(mgmt_env, batched, has_reads=False)

    filed, outcome = _run(mgmt_analysis.mgmt_handler, job, work_dir=str(tmp_path / "w"))

    assert filed == []
    assert outcome["status"] == "success"
    assert job.context.results["mgmt_analysis"]["status"] == "no_mgmt_reads"
    assert mgmt_env.rscript_calls == 0
