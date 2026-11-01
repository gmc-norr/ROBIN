from __future__ import annotations

import threading
from collections import deque
from typing import Any, Optional, Dict, List, Tuple
from pathlib import Path

try:
    from nicegui import ui
except ImportError:  # pragma: no cover
    ui = None

# We intentionally import the heavy launcher class from the legacy module
# to keep code changes minimal while we migrate features into this package.
from ..gui_launcher import GUILauncher  # type: ignore

_current_gui_launcher: Optional[GUILauncher] = None


def launch_gui(
    host: str = "0.0.0.0",
    port: int = 8081,
    show: bool = False,
    workflow_runner: Any = None,
    workflow_steps: Optional[list] = None,
    monitored_directory: str = "",
    reload: bool = False,
    center: str = None,
) -> GUILauncher:
    """Central entrypoint for launching the robin GUI (refactored).

    This function mirrors the legacy signature and behavior while living under
    the package `robin.gui`. It returns the launcher instance and sets a
    global reference for update dispatching.
    """
    if ui is None:
        raise ImportError(
            "NiceGUI is not available. Please install it with: pip install nicegui"
        )

    launcher = GUILauncher(host, port, reload)

    # Normalize monitored directory to an absolute path to avoid CWD issues
    abs_mon_dir = (
        str(Path(monitored_directory).resolve()) if monitored_directory else ""
    )
    success = launcher.launch_gui(workflow_runner, workflow_steps, abs_mon_dir, reload, center)
    if not success:
        raise RuntimeError("Failed to launch GUI")

    global _current_gui_launcher
    _current_gui_launcher = launcher

    if show:
        import webbrowser  # local to avoid unnecessary import at module import time

        webbrowser.open(launcher.get_gui_url())

    return launcher


def get_gui_launcher() -> Optional[GUILauncher]:
    return _current_gui_launcher


# In a process with no GUI (a Ray worker) updates are kept here instead of being dropped; the job
# wrapper moves them into the job's context so they travel back with the finished job. Bounded:
# only the newest OUTBOX_MAX are kept.
OUTBOX_MAX = 100
_outbox: "deque[Tuple[Any, Dict[str, Any], int]]" = deque(maxlen=OUTBOX_MAX)
_outbox_lock = threading.Lock()


def send_gui_update(update_type: Any, data: Dict[str, Any], priority: int = 0) -> None:
    launcher = get_gui_launcher()
    if launcher:
        launcher.send_update(update_type, data, priority)
        return
    # No GUI here: keep the update (as plain values, safe to pickle) for the job to carry back.
    with _outbox_lock:
        _outbox.append((getattr(update_type, "value", update_type), data, priority))


def drain_gui_outbox() -> List[Tuple[Any, Dict[str, Any], int]]:
    """Return and clear the updates kept in this process, oldest first."""
    with _outbox_lock:
        items = list(_outbox)
        _outbox.clear()
    return items
