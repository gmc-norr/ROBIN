"""GUI-4: the CNV summary card showed the batch counter as the number of gains.

``_extract_cnv_data`` looked for ``cnv_analysis_counter.txt`` first and read it as
the gained-region count, so after 37 batches every sample's card said
"Gained: 37, Lost: 0". The counter is a batch counter, written by
``save_analysis_counter`` in ``robin/analysis/cnv_analysis.py``.

Nothing in ``src/`` writes a per-sample gain/loss count, so the honest answer is
that the counts are unknown — ``None`` — and the card must say "Not available"
instead of inventing a number.

``test_the_card_renders_unknown_counts_without_crashing`` is the test that pins
the wiring: ``_create_cnv_dashboard_card_with_data`` did ``int(gained) + int(lost)``,
which raises ``TypeError`` on ``None``. Widening the defaults without guarding that
line would break every sample's summary page at render time while an import check
stayed happy, so the card is driven here with a stub ``ui``.
"""

from typing import Any, List, Optional

from robin.gui.components import summary
from robin.gui.components.summary import (
    _create_cnv_dashboard_card_with_data,
    _extract_cnv_data,
)


# --------------------------------------------------------------------------- #
# A stub ``ui`` just rich enough for the CNV card: every builder returns a node
# that is chainable (``classes``/``props``/``on``) and usable as a context
# manager, and every label records its text.
# --------------------------------------------------------------------------- #


class _StubNode:
    def __init__(self, texts: List[str], text: Optional[str] = None) -> None:
        self._texts = texts
        if text is not None:
            texts.append(text)

    def classes(self, *args: Any, **kwargs: Any) -> "_StubNode":
        return self

    def props(self, *args: Any, **kwargs: Any) -> "_StubNode":
        return self

    def on(self, *args: Any, **kwargs: Any) -> "_StubNode":
        return self

    def __enter__(self) -> "_StubNode":
        return self

    def __exit__(self, *exc: Any) -> bool:
        return False


class _StubUi:
    def __init__(self) -> None:
        self.texts: List[str] = []

    def element(self, *args: Any, **kwargs: Any) -> _StubNode:
        return _StubNode(self.texts)

    def column(self, *args: Any, **kwargs: Any) -> _StubNode:
        return _StubNode(self.texts)

    def row(self, *args: Any, **kwargs: Any) -> _StubNode:
        return _StubNode(self.texts)

    def icon(self, *args: Any, **kwargs: Any) -> _StubNode:
        return _StubNode(self.texts)

    def label(self, text: Any = "", *args: Any, **kwargs: Any) -> _StubNode:
        return _StubNode(self.texts, str(text))


def _install_stub_ui(monkeypatch) -> List[str]:
    stub = _StubUi()
    monkeypatch.setattr(summary, "ui", stub)
    return stub.texts


# --------------------------------------------------------------------------- #
# _extract_cnv_data
# --------------------------------------------------------------------------- #


def test_batch_counter_is_not_reported_as_gains(tmp_path):
    """GUI-4: after 37 batches the card said 'Gained: 37'."""
    (tmp_path / "cnv_analysis_counter.txt").write_text("37\n")
    data = _extract_cnv_data(tmp_path)
    assert data["gained"] is None
    assert data["lost"] is None


def test_no_cnv_files_leaves_the_counts_unknown(tmp_path):
    data = _extract_cnv_data(tmp_path)
    assert data["gained"] is None and data["lost"] is None
    assert data["genetic_sex"] == "Not available"


def test_a_real_results_csv_still_supplies_the_counts(tmp_path):
    (tmp_path / "cnv_results.csv").write_text("gained,lost,genetic_sex\n3,2,XY\n")
    data = _extract_cnv_data(tmp_path)
    assert data["gained"] == 3
    assert data["lost"] == 2
    assert data["genetic_sex"] == "XY"


def test_a_counter_file_does_not_shadow_the_results_csv(tmp_path):
    """The counter name must leave the candidate list, not just lose its branch.

    Keeping the name and only deleting its ``.txt`` branch would still ``break``
    the loop on the counter and never reach the real counts in the ``.csv``.
    """
    (tmp_path / "cnv_analysis_counter.txt").write_text("37\n")
    (tmp_path / "cnv_results.csv").write_text("gained,lost,genetic_sex\n3,2,XY\n")
    data = _extract_cnv_data(tmp_path)
    assert data["gained"] == 3
    assert data["lost"] == 2


# --------------------------------------------------------------------------- #
# _create_cnv_dashboard_card_with_data
# --------------------------------------------------------------------------- #


def test_the_card_renders_unknown_counts_without_crashing(monkeypatch):
    """GUI-4: `int(None)` at summary.py:650 would take down the whole summary page."""
    texts = _install_stub_ui(monkeypatch)
    _create_cnv_dashboard_card_with_data("XY", "1,000", "0.1", None, None)
    joined = " | ".join(texts)
    assert "Gained: 0" not in joined
    assert "Lost: 0" not in joined
    assert "Not available" in joined
    # "Gained vs total: Not available / —" would be worse than no line at all.
    assert "Gained vs total: —" in joined


def test_the_card_still_renders_real_counts(monkeypatch):
    """Guard: a sample that does have counts must still see them (passes pre-fix)."""
    texts = _install_stub_ui(monkeypatch)
    _create_cnv_dashboard_card_with_data("XY", "1,000", "0.1", 3, 2)
    joined = " | ".join(texts)
    assert "Gained vs total: 3 / 5" in joined
    assert "Gained: 3" in joined
    assert "Lost: 2" in joined
    assert "Not available" not in joined
