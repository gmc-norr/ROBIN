"""A step that fails or never runs must be counted as failed, on both schedulers (batch 7, ERR-1).

Every case goes through a real local Ray ``Coordinator`` (see ``failure_fixtures.py``), once on
``standard`` (pool actors) and once on ``high`` (one processor per job type). Only the job
handlers are synthetic; each is wrapped exactly like the real ones.

Rows from the locked plan:
  F1  the handler FILES an error under the name the real handler uses  -> counted failed
  F2  the handler RAISES                                               -> counted failed
  C1  the handler succeeds                                             -> counted completed
  C2  a classifier / bed_conversion handler files an error (or raises) on a fail-only BAM
      submission                                                       -> not failed; result
                                                                          is expected_failure
  C3  preprocessing marks skip_downstream                              -> counted skipped

Each scheduler's whole matrix runs once, in a module-scoped fixture, then the Coordinator is shut
down (which is what writes ``failed_jobs.csv``) and the tests assert on what was recorded.
"""

import pytest

from failure_fixtures import (
    FAIL_ONLY_TYPES, JOB_TYPES, SCHEDULERS, Harness, init_ray, shutdown_ray,
)


@pytest.fixture(scope="module")
def ray_cluster():
    init_ray()
    yield
    shutdown_ray()


class Campaign:
    def __init__(self):
        self.outcomes = {}  # (case, job_type) -> Outcome
        self.csv_rows = []  # failed_jobs.csv rows after shutdown

    def rows_for(self, outcome):
        return [r for r in self.csv_rows if r["sample_id"] == outcome.sample_id]


@pytest.fixture(scope="module", params=SCHEDULERS)
def campaign(request, ray_cluster, tmp_path_factory):
    scheduler = request.param
    work = str(tmp_path_factory.mktemp(f"failures_{scheduler}"))
    harness = Harness(scheduler, work)
    c = Campaign()
    c.scheduler = scheduler
    try:
        for jt in JOB_TYPES:
            c.outcomes[("F1", jt)] = harness.run_one(jt, "file_error", tag="F1")
            c.outcomes[("F2", jt)] = harness.run_one(jt, "raise", tag="F2")
            c.outcomes[("C1", jt)] = harness.run_one(jt, "ok", tag="C1")
        for jt in FAIL_ONLY_TYPES:
            c.outcomes[("C2", jt)] = harness.run_one(jt, "file_error", tag="C2", fail_only=True)
            c.outcomes[("C2raise", jt)] = harness.run_one(jt, "raise", tag="C2", fail_only=True)
        c.outcomes[("C3", "preprocessing")] = harness.run_one("preprocessing", "skip", tag="C3")
    finally:
        c.csv_rows = harness.close()
    return c


def _counted_failed(c, o, text):
    assert o.d_failed == 1, f"{c.scheduler}/{o.job_type}: failed count moved by {o.d_failed}, expected 1"
    assert o.d_completed == 0, f"{c.scheduler}/{o.job_type}: also counted completed"
    assert o.d_failed_by_type == {o.job_type: 1}, (
        f"{c.scheduler}/{o.job_type}: failed_by_type moved by {o.d_failed_by_type}"
    )
    rows = c.rows_for(o)
    assert len(rows) == 1, f"{c.scheduler}/{o.job_type}: {len(rows)} failed_jobs.csv rows, expected 1"
    assert rows[0]["job_type"] == o.job_type
    assert text in rows[0]["error"], f"csv error text {rows[0]['error']!r} lacks {text!r}"


@pytest.mark.parametrize("job_type", JOB_TYPES)
def test_F1_filed_error_counts_as_failed(campaign, job_type):
    o = campaign.outcomes[("F1", job_type)]
    _counted_failed(campaign, o, f"synthetic {job_type} failure")


@pytest.mark.parametrize("job_type", JOB_TYPES)
def test_F2_raised_exception_counts_as_failed(campaign, job_type):
    o = campaign.outcomes[("F2", job_type)]
    _counted_failed(campaign, o, f"synthetic {job_type} crash")


@pytest.mark.parametrize("job_type", JOB_TYPES)
def test_C1_success_counts_as_completed(campaign, job_type):
    o = campaign.outcomes[("C1", job_type)]
    assert o.d_completed == 1, f"{campaign.scheduler}/{job_type}: completed moved by {o.d_completed}"
    assert o.d_failed == 0 and o.d_failed_by_type == {}
    assert o.d_completed_by_type == {job_type: 1}
    assert campaign.rows_for(o) == []


@pytest.mark.parametrize("job_type", FAIL_ONLY_TYPES)
def test_C2_filed_error_on_fail_only_bam_is_not_failed(campaign, job_type):
    o = campaign.outcomes[("C2", job_type)]
    assert o.d_failed == 0 and o.d_failed_by_type == {}, (
        f"{campaign.scheduler}/{job_type}: a fail-only BAM's filed error counted as failed"
    )
    assert o.d_completed == 1
    assert campaign.rows_for(o) == []


@pytest.mark.parametrize("job_type", FAIL_ONLY_TYPES)
def test_C2_filed_error_on_fail_only_bam_result_is_expected_failure(campaign, job_type):
    o = campaign.outcomes[("C2", job_type)]
    assert o.result.get("status") == "expected_failure", (
        f"{campaign.scheduler}/{job_type}: result is {o.result!r}"
    )
    assert o.result.get("reason") == "fail_only_bam_submission"


@pytest.mark.parametrize("job_type", FAIL_ONLY_TYPES)
def test_C2_raised_exception_on_fail_only_bam_is_expected_failure(campaign, job_type):
    """The behaviour the filed-error case must match: it works today for an exception."""
    o = campaign.outcomes[("C2raise", job_type)]
    assert o.d_failed == 0 and o.d_completed == 1
    assert o.result.get("status") == "expected_failure"
    assert campaign.rows_for(o) == []


def test_C3_deliberate_skip_is_skipped_not_failed(campaign):
    o = campaign.outcomes[("C3", "preprocessing")]
    assert o.d_skipped == 1, f"{campaign.scheduler}: skipped moved by {o.d_skipped}"
    assert o.d_failed == 0 and o.d_completed == 0 and o.d_failed_by_type == {}
    assert campaign.rows_for(o) == []
