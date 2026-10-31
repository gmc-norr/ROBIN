"""GUI-2: refresh timers must be torn down on client delete, not on disconnect.

Every ROBIN GUI page starts a periodic `ui.timer` and deactivates it from a
client callback. In NiceGUI 3 `Client.on_disconnect` also fires when a client
*reconnects* ("Updated in version 3.0.0: The handler is also called when a
client reconnects"), and nothing in ROBIN reactivates the timers afterwards.
A Wi-Fi blip, a laptop sleep or a proxy timeout therefore left an open page
showing hours-old results with no sign that they were stale.

`Client.on_delete` is the correct hook: it runs only once the client is really
discarded, after NiceGUI's `reconnect_timeout` (60 s, set in gui_launcher.py).

These are lint-shaped tests: they read the GUI sources as text. Driving a real
reconnect needs a live browser and a real websocket, neither of which exists
here, and a stubbed `nicegui` in `sys.modules` would poison the real-NiceGUI
imports other tests in the same session depend on. So this pins the one thing
that is checkable without a browser -- that no GUI source hangs a timer
teardown off `on_disconnect`. `tests/test_dependency_pins.py` reads files as
text for the same reason.

NOTE: no test here proves that a timer survives a reconnect.
"""

from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

GUI_FILES = [  # the eight sites GUI-2 names, plus gui_launcher as a ninth
    "src/robin/gui/components/summary.py",
    "src/robin/gui/components/classification.py",
    "src/robin/gui/components/coverage.py",
    "src/robin/gui/components/mgmt.py",
    "src/robin/gui/components/cnv.py",
    "src/robin/gui/components/fusion.py",
    "src/robin/gui/components/mnpflex.py",
    "src/robin/gui/components/bed_coverage.py",
    "src/robin/gui_launcher.py",
]


def test_no_refresh_timer_is_torn_down_on_disconnect():
    """GUI-2: on_disconnect also fires on reconnect in NiceGUI 3, and nothing reactivates.

    Scans every source file, not only GUI_FILES, so a new page that hangs its
    timer off on_disconnect is caught too, and shortening the list cannot make
    this pass."""
    scanned = sorted((REPO_ROOT / "src" / "robin").rglob("*.py"))
    scanned_rel = {p.relative_to(REPO_ROOT).as_posix() for p in scanned}
    assert set(GUI_FILES) <= scanned_rel
    offenders = []
    for path in scanned:
        for n, line in enumerate(path.read_text().splitlines(), 1):
            if "on_disconnect" in line:
                offenders.append(f"{path.relative_to(REPO_ROOT)}:{n}")
    assert offenders == []


def test_every_file_that_had_one_now_uses_on_delete():
    # Nine sites: a shorter list would quietly test less.
    assert len(GUI_FILES) == 9
    for rel in GUI_FILES:
        # The call itself: a bare "on_delete" would also match "on_deleted".
        assert "client.on_delete(" in (REPO_ROOT / rel).read_text(), rel


def test_on_delete_exists_in_the_installed_nicegui():
    """If NiceGUI ever drops on_delete, this fix is silently a no-op."""
    from nicegui import Client

    assert callable(getattr(Client, "on_delete", None))
