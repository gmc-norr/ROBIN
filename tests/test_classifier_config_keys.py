"""REP-1: the report must resolve a classifier's display name to its config key.

Each classifier has its own confidence band. Random Forest's is 85/65, not the
80/50 default. The PDF report used to look the band up by lowercasing the
display name, which gives ``"random forest"`` with a space — not the
``"random_forest"`` config key — so Random Forest was silently scored against
the default band and 82% read "High" in the PDF while the GUI said "Medium".

``test_the_report_bands_random_forest_at_82_percent_as_medium`` is the test that
reddens for REP-1's own reason: it drives ``ClassificationSection.add_content``
over a synthetic Random Forest score file and records the key the section
actually passes to ``get_confidence_status`` along with the band that comes
back. Against ``50da100`` it records ``('random forest', 0.82, 'High')``.

That test runs in a subprocess because the section's module imports ``seaborn``,
which is not installed here, so it needs a stub — and installing a fake module
into ``sys.modules`` would leak into every other test in the pytest session.

The other tests cover the resolver directly and so can only fail by
``ImportError`` on the old code. ``test_non_classifier_anchors_...`` is a guard:
it passes before and after, and pins that the three threshold lookups stay
permissive for the summary page's section anchors.
"""

import os
import subprocess
import sys
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
from robin.gui.config import CLASSIFICATION_STEPS

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


def test_the_display_names_match_the_guis_classification_steps():
    """The resolver's names must stay the exact inverse of the GUI's, or a renamed
    classifier silently loses its band again — today that drift shows up only at
    runtime, as a missing report row plus a log line."""
    assert CLASSIFIER_DISPLAY_NAMES == {v: k for k, v in CLASSIFICATION_STEPS.items()}


def test_the_report_resolves_the_display_name_before_looking_up_thresholds():
    """Covers reporting/sections/classification.py:399 itself. Cheap, and it catches a
    failure the behavioural test below does not: someone dropping the import. It does
    NOT prove the fix on its own — renaming the ``name`` variable would satisfy it."""
    src = (REPO_ROOT / "src/robin/reporting/sections/classification.py").read_text()
    assert "get_confidence_status(name.lower()" not in src
    assert "config_key_for_classifier" in src


_DRIVER_SCRIPT = '''\
"""Drive ClassificationSection.add_content over one Random Forest score file.

Records every (classifier_key, confidence, status) the section passes to
get_confidence_status and prints them. argv[1] is the report output dir.
"""
import logging
import pathlib
import sys
import types

logging.disable(logging.CRITICAL)

# The section module imports seaborn, which is not installed here. Only
# _create_time_plot touches it, and that call is already wrapped in its own
# try/except, so an attribute-less stub is enough to reach the lookup at :399.
sys.modules.setdefault("seaborn", types.ModuleType("seaborn"))

from robin.reporting.sections import classification as cs
from robin.reporting.styling.styles import ReportStyles

out = pathlib.Path(sys.argv[1])
# Random Forest scores are 0-100, so 82.0 becomes a confidence of 0.82: above
# the 80 default high cut-off, below Random Forest's own 85.
(out / "random_forest_scores.csv").write_text(
    "timestamp,Glioblastoma,Medulloblastoma,number_probes\\n1000,82.0,18.0,5000\\n"
)


class Rep:
    """The two attributes ReportSection.__init__ and the loop actually read."""

    output = str(out)
    styles = ReportStyles(str(out / "no-fonts"))  # falls back to Helvetica


recorded = []
real_get_confidence_status = cs.get_confidence_status


def recorder(classifier, confidence):
    status, color = real_get_confidence_status(classifier, confidence)
    recorded.append((classifier, confidence, status))
    return status, color


cs.get_confidence_status = recorder

try:
    cs.ClassificationSection(Rep()).add_content()
except AttributeError as exc:
    # add_content runs past the summary loop into the detailed section, which
    # reads self.report.doc.width; this stub report has no doc. The lookup under
    # test has already happened by then, so the recorder is the assertion target.
    if "doc" not in str(exc):
        raise

print("RECORDED %r" % (recorded,))
'''


def test_the_report_bands_random_forest_at_82_percent_as_medium(tmp_path):
    """REP-1 end to end. The section must hand get_confidence_status the config key
    'random_forest', not 'random forest'. Against 50da100 this records
    ('random forest', 0.82, 'High') - the PDF saying High where the GUI said Medium."""
    script = tmp_path / "drive.py"
    script.write_text(_DRIVER_SCRIPT)
    out_dir = tmp_path / "out"
    out_dir.mkdir()

    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(
        [str(REPO_ROOT / "src"), env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    env["PYTHONDONTWRITEBYTECODE"] = "1"

    r = subprocess.run(
        [sys.executable, str(script), str(out_dir)],
        capture_output=True,
        text=True,
        env=env,
    )
    assert r.returncode == 0, r.stderr
    recorded = [
        line.split(" ", 1)[1] for line in r.stdout.splitlines() if line.startswith("RECORDED ")
    ]
    assert len(recorded) == 1, r.stdout + r.stderr
    assert recorded[0] == "[('random_forest', 0.82, 'Medium')]", recorded[0]


def test_non_classifier_anchors_still_fall_back_to_the_defaults():
    """Review Focus 1: these must not raise — the summary cards depend on them."""
    for anchor in ("coverage", "cnv", "mgmt", "fusion", "default"):
        assert get_confidence_ui_tier(anchor, 82.0) == "high"   # 82 >= default high of 80
        assert get_confidence_level(anchor, 82.0) == "High confidence"
