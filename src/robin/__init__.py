"""ROBIN and robin - A Python CLI tool with file watching capabilities."""

import os

# Keep NiceGUI's storage (logins, disclaimer state) out of the folder Robin is
# started from. NiceGUI reads this once, when first imported, so it must be set
# before anything imports NiceGUI. An existing value (including the one Ray
# gives its workers) wins.
os.environ.setdefault(
    "NICEGUI_STORAGE_PATH", os.path.join(os.path.expanduser("~"), ".robin", "nicegui")
)

# Suppress pkg_resources deprecation warnings from sorted_nearest
import warnings
warnings.filterwarnings(
    "ignore", message="pkg_resources is deprecated", category=UserWarning
)
# Suppress matplotlib tight_layout warnings
warnings.filterwarnings(
    "ignore", message="The figure layout has changed to tight", category=UserWarning
)

__version__ = "0.5.1"
__author__ = "Matt Loose"
__email__ = "matt.loose@nottingham.ac.uk"

#
# Import analysis modules via consolidated analysis package.
#
# Some analysis dependencies are optional and may not be installed in minimal
# environments (e.g. when only utility modules are imported).
try:
    from . import analysis  # noqa: F401
except ModuleNotFoundError:
    analysis = None  # type: ignore[assignment]

# Export key modules
__all__ = [
    "analysis",
]
