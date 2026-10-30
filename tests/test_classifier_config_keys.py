"""REP-1: the report must resolve a classifier's display name to its config key.

Each classifier has its own confidence band. Random Forest's is 85/65, not the
80/50 default. The PDF report used to look the band up by lowercasing the
display name, which gives ``"random forest"`` with a space — not the
``"random_forest"`` config key — so Random Forest was silently scored against
the default band and 82% read "High" in the PDF while the GUI said "Medium".
"""

from pathlib import Path

import pytest

from robin.classification_config import (
    CLASSIFIER_CONFIDENCE_THRESHOLDS,
    CLASSIFIER_DISPLAY_NAMES,
    config_key_for_classifier,
    get_confidence_level,
    get_confidence_status,
    get_confidence_ui_tier,
)

REPO_ROOT = Path(__file__).resolve().parents[1]


def test_random_forest_display_name_maps_to_the_config_key():
    assert config_key_for_classifier("Random Forest") == "random_forest"


def test_every_display_name_resolves_to_a_real_threshold_entry():
    for display in CLASSIFIER_DISPLAY_NAMES:
        assert config_key_for_classifier(display) in CLASSIFIER_CONFIDENCE_THRESHOLDS


def test_resolution_is_case_and_whitespace_insensitive():
    assert config_key_for_classifier("  random forest ") == "random_forest"
    assert config_key_for_classifier("STURGEON") == "sturgeon"


def test_unknown_classifier_raises_and_names_itself():
    with pytest.raises(ValueError, match="Dinosaur"):
        config_key_for_classifier("Dinosaur")


def test_random_forest_at_82_percent_is_medium_not_high():
    """REP-1: the whole point. 82% is above the 80 default but below the RF 85."""
    status, _ = get_confidence_status(config_key_for_classifier("Random Forest"), 0.82)
    assert status == "Medium"


def test_the_report_resolves_the_display_name_before_looking_up_thresholds():
    """Covers reporting/sections/classification.py:399 itself — without this, nothing
    proves the call site was changed, only that the resolver exists."""
    src = (REPO_ROOT / "src/robin/reporting/sections/classification.py").read_text()
    assert "get_confidence_status(name.lower()" not in src
    assert "config_key_for_classifier" in src


def test_non_classifier_anchors_still_fall_back_to_the_defaults():
    """Review Focus 1: these must not raise — the summary cards depend on them."""
    for anchor in ("coverage", "cnv", "mgmt", "fusion", "default"):
        assert get_confidence_ui_tier(anchor, 82.0) == "high"   # 82 >= default high of 80
        assert get_confidence_level(anchor, 82.0) == "High confidence"
