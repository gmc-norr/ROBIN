"""JOB-6: the sample-completion idle timeout must default to an hour, not 60 s.

The timeout decides when a quiet sample is marked Complete and *finalised* —
its batch BAMs merged into ``target.bam`` and then deleted. A 60-second gap
between BAMs is routine on a real run, so the old test value finalised samples
mid-run. See ``completion_timeout_seconds`` in ``robin.gui_launcher``.
"""

import pytest

from robin.gui_launcher import completion_timeout_seconds


def test_default_is_one_hour(monkeypatch):
    monkeypatch.delenv("ROBIN_COMPLETION_TIMEOUT_SECONDS", raising=False)
    assert completion_timeout_seconds() == 3600


def test_env_var_overrides_the_default(monkeypatch):
    monkeypatch.setenv("ROBIN_COMPLETION_TIMEOUT_SECONDS", "120")
    assert completion_timeout_seconds() == 120


@pytest.mark.parametrize("raw", ["", "   ", "abc", "12.5", "-1", "0"])
def test_unusable_values_fall_back_to_the_default(monkeypatch, raw):
    """Review Focus 2: a value that became 0 would finalise every sample on the first tick."""
    monkeypatch.setenv("ROBIN_COMPLETION_TIMEOUT_SECONDS", raw)
    assert completion_timeout_seconds() == 3600


def test_surrounding_whitespace_is_tolerated(monkeypatch):
    monkeypatch.setenv("ROBIN_COMPLETION_TIMEOUT_SECONDS", " 900 ")
    assert completion_timeout_seconds() == 900
