"""PKG-3/PKG-6: the pins that keep a fresh or updated environment working.

Each pin here stopped a real break on the clinical server: pip under Python
3.14, conda moving scikit-learn past the models' version, and conda moving
numpy to 2.x under pandas 2.1.4 ("numpy.dtype size changed").
"""

import tomllib
from pathlib import Path

from packaging.requirements import Requirement
from packaging.specifiers import SpecifierSet
from packaging.version import Version

ROOT = Path(__file__).resolve().parents[1]

# Both NN model pickles (Capper_et_al_NN.pkl, pancan_devel_v5i_NN.pkl) record
# _sklearn_version 1.8.0.
MODEL_SKLEARN_VERSION = "1.8.0"


def _pyproject():
    return tomllib.loads((ROOT / "pyproject.toml").read_text())


def _pip_requirements(name):
    deps = [Requirement(d) for d in _pyproject()["project"]["dependencies"]]
    return [r for r in deps if r.name.lower() == name]


def _conda_pin(name):
    """Return the exact version robin.yml pins for a package, or None if it isn't pinned."""
    for line in (ROOT / "robin.yml").read_text().splitlines():
        entry = line.split("#")[0].strip()
        if entry.startswith("- "):
            spec = entry[2:].strip()
            pkg, _, version = spec.partition("=")
            if pkg == name:
                return version if version and not version.startswith(("=", "<", ">")) else None
    return None


def test_only_python_312_is_accepted():
    spec = SpecifierSet(_pyproject()["project"]["requires-python"])
    assert "3.12.3" in spec
    assert "3.11.9" not in spec
    assert "3.13.0" not in spec
    assert "3.14.0" not in spec


def test_scikit_learn_matches_the_models_in_pip_and_conda():
    reqs = _pip_requirements("scikit-learn")
    assert len(reqs) == 1
    assert reqs[0].specifier == SpecifierSet(f"=={MODEL_SKLEARN_VERSION}")
    assert _conda_pin("scikit-learn") == MODEL_SKLEARN_VERSION


def test_conda_numpy_and_pandas_are_pinned_inside_the_pip_ranges():
    numpy, pandas = _conda_pin("numpy"), _conda_pin("pandas")
    assert numpy is not None and pandas is not None
    # pandas 2.1.4 is built for numpy 1.x
    assert Version(numpy).major == 1
    assert all(numpy in r.specifier for r in _pip_requirements("numpy"))
    assert all(pandas in r.specifier for r in _pip_requirements("pandas"))
