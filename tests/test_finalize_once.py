"""The GUI's target.bam finalize (batch 10b: JOB-7, JOB-6). Synthetic BAMs only.

JOB-7: when the finalize job is queued, the GUI must not run the same merge itself too.
JOB-6: a sample finalized once in a session must still have later reads merged, both
batch BAMs and files still staged.

A real GUILauncher is built without starting NiceGUI. The finalize runs in a thread, so
``threading.Thread`` is swapped for one that runs its target at once.
"""

import os
import threading
import time

import pytest

from target_fixtures import SAMPLE, bam_records, samtools_reference, write_bams, write_panel
from robin.analysis import target_analysis as ta_mod
from robin.analysis.target_analysis import TargetAnalysis
from robin.gui_launcher import GUILauncher

PANEL_NAME = "TESTPANEL"


class _RunsAtOnce:
    """Stands in for threading.Thread: start() runs the target in the caller's thread."""

    def __init__(self, target=None, args=(), kwargs=None, daemon=None, name=None):
        self._target, self._args, self._kwargs = target, args, kwargs or {}
        self.name = name or "runs-at-once"

    def start(self):
        self._target(*self._args, **self._kwargs)


class _Runner:
    """A workflow runner that records what the GUI sends it."""

    def __init__(self, queues_finalize):
        self.queues_finalize = queues_finalize
        self.finalize, self.snp = [], []

    def submit_target_bam_finalize_job(self, **kwargs):
        self.finalize.append(kwargs)
        return self.queues_finalize

    def submit_snp_analysis_job(self, **kwargs):
        self.snp.append(kwargs)
        return True


@pytest.fixture
def panel(tmp_path, monkeypatch):
    path = write_panel(str(tmp_path / "panel.bed"))
    monkeypatch.setattr(TargetAnalysis, "_find_target_bed", lambda self, p: path)
    return path


def _master_csv(sample_dir, last_seen):
    with open(sample_dir / "master.csv", "w") as f:
        f.write(
            "analysis_panel,samples_overview_last_seen,samples_overview_active_jobs,"
            "samples_overview_pending_jobs\n"
        )
        f.write(f"{PANEL_NAME},{last_seen},0,0\n")


@pytest.fixture
def gui(tmp_path, monkeypatch):
    g = GUILauncher()  # before the swap: building it must not run anything
    work = tmp_path / "work"
    sample_dir = work / SAMPLE
    sample_dir.mkdir(parents=True)
    _master_csv(sample_dir, time.time())
    g.monitored_directory = str(work)
    g.workflow_runner = None
    monkeypatch.setattr(g, "send_update", lambda *a, **k: None)
    monkeypatch.setattr(threading, "Thread", _RunsAtOnce)
    return g


def _count_finalizes(monkeypatch, real=True):
    """Count the GUI's own merges; run the real one, or just answer success."""
    calls = []
    real_finalize = ta_mod.finalize_accumulation_for_sample

    def counted(*args, **kwargs):
        calls.append(kwargs.get("sample_id", args[0] if args else None))
        if real:
            return real_finalize(*args, **kwargs)
        return {"status": "success", "sample_id": SAMPLE}

    monkeypatch.setattr(ta_mod, "finalize_accumulation_for_sample", counted)
    return calls


def _analysis(gui):
    return TargetAnalysis(work_dir=gui.monitored_directory, target_panel=PANEL_NAME,
                          batch_size=1, use_staging=True)


def _stage(ta, bam):
    ta.process_file_with_staging(bam, {"sample_id": SAMPLE})


def _batch(ta, bam):
    _stage(ta, bam)
    ta.accumulate_staged_files(SAMPLE, force=True)


def _target_bam(gui):
    return os.path.join(gui.monitored_directory, SAMPLE, "target.bam")


# JOB-7 ------------------------------------------------------------------------------


def test_a_queued_finalize_is_not_also_run_by_the_gui(gui, monkeypatch):
    calls = _count_finalizes(monkeypatch, real=False)
    gui.workflow_runner = _Runner(queues_finalize=True)
    gui._trigger_target_bam_finalization(SAMPLE)
    assert len(gui.workflow_runner.finalize) == 1
    assert calls == [], "the GUI ran the merge itself as well as queueing it"


def test_a_queued_finalize_then_snp_sends_one_finalize_and_one_snp_job(gui, monkeypatch):
    calls = _count_finalizes(monkeypatch, real=False)
    gui.workflow_runner = _Runner(queues_finalize=True)
    open(_target_bam(gui), "w").close()  # an older target.bam; no index, so not "done"
    monkeypatch.setattr(gui, "_wait_for_target_bam_or_timeout", lambda *a, **k: True)
    monkeypatch.setattr(gui, "_wait_for_snp_outputs_or_timeout", lambda *a, **k: True)
    monkeypatch.setattr(gui, "_locate_reference_for_sample", lambda d: "/ref.fa")
    monkeypatch.setattr(ta_mod, "is_docker_available_for_snp_analysis", lambda: (True, ""))
    gui._trigger_target_bam_finalization(SAMPLE, trigger_snp=True)
    assert len(gui.workflow_runner.finalize) == 1
    assert len(gui.workflow_runner.snp) == 1
    assert calls == [], "the GUI ran the merge itself as well as queueing it"


@pytest.mark.parametrize("runner", [None, _Runner(queues_finalize=False)], ids=["no runner", "queue refused"])
def test_the_gui_merges_itself_when_the_job_cannot_be_queued(gui, monkeypatch, runner):
    calls = _count_finalizes(monkeypatch, real=False)
    gui.workflow_runner = runner
    gui._trigger_target_bam_finalization(SAMPLE)
    assert calls == [SAMPLE]


# JOB-6 ------------------------------------------------------------------------------


def test_a_second_finalize_in_one_session_merges_new_batches(gui, panel, tmp_path, monkeypatch):
    calls = _count_finalizes(monkeypatch)
    bams, merged = write_bams(str(tmp_path), 2)
    ta = _analysis(gui)
    _batch(ta, bams[0])
    gui._trigger_target_bam_finalization(SAMPLE)
    assert bam_records(_target_bam(gui)) == samtools_reference(bams[0], panel)["target_bam"]
    _batch(ta, bams[1])
    gui._trigger_target_bam_finalization(SAMPLE)
    assert len(calls) == 2, "the second finalize skipped the merge"
    assert bam_records(_target_bam(gui)) == samtools_reference(merged, panel)["target_bam"]


def test_files_still_staged_are_merged_after_an_earlier_finalize(gui, panel, tmp_path, monkeypatch):
    _count_finalizes(monkeypatch)
    bams, merged = write_bams(str(tmp_path), 2)
    ta = _analysis(gui)
    _batch(ta, bams[0])
    gui._trigger_target_bam_finalization(SAMPLE)
    _stage(ta, bams[1])  # staged, not yet made into a batch BAM
    assert not gui._is_target_bam_finalize_redundant(SAMPLE)
    gui._trigger_target_bam_finalization(SAMPLE)
    assert bam_records(_target_bam(gui)) == samtools_reference(merged, panel)["target_bam"]


def test_a_sample_that_completes_twice_is_finalized_twice(gui, monkeypatch):
    """The background scan starts a finalize on each Live -> Complete, not only the first."""
    triggered = []
    monkeypatch.setattr(gui, "_trigger_target_bam_finalization",
                        lambda sid, **k: triggered.append(sid))
    sample_dir = os.path.join(gui.monitored_directory, SAMPLE)
    old = time.time() - 2 * gui.completion_timeout_seconds

    from pathlib import Path
    _master_csv(Path(sample_dir), old)
    gui._background_scan_samples_sync()
    assert gui._samples_master_record[SAMPLE].origin == "Complete"
    assert triggered == [SAMPLE]

    _master_csv(Path(sample_dir), time.time())  # new reads: the sample is Live again
    gui._background_scan_samples_sync()
    assert gui._samples_master_record[SAMPLE].origin == "Live"

    _master_csv(Path(sample_dir), old)
    gui._background_scan_samples_sync()
    assert gui._samples_master_record[SAMPLE].origin == "Complete"
    assert triggered == [SAMPLE, SAMPLE], "the second completion started no finalize"
