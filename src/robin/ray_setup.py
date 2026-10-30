"""Shared settings for starting Ray, used by every ``ray.init`` call in Robin."""

import os
import tempfile

DEFAULT_DASHBOARD_HOST = "127.0.0.1"


def ray_dashboard_host() -> str:
    """Return the address the Ray dashboard listens on.

    The dashboard's Jobs API runs any command it is sent and has no login, so it
    listens on localhost only. An operator who wants it elsewhere sets
    ``RAY_DASHBOARD_HOST``.
    """
    return os.environ.get("RAY_DASHBOARD_HOST") or DEFAULT_DASHBOARD_HOST


def ray_worker_runtime_env() -> dict:
    """Return the Ray runtime environment for Robin's workers.

    Every Ray worker imports NiceGUI, and NiceGUI rewrites its storage file at
    import through one fixed ``.tmp`` name. Workers sharing the GUI's storage
    folder race on that file. Each run gives its workers a fresh, empty folder,
    where NiceGUI finds nothing to write.
    """
    storage = tempfile.mkdtemp(prefix="robin-ray-nicegui-")
    return {"env_vars": {"NICEGUI_STORAGE_PATH": storage}}


def ray_init_options(include_dashboard: bool) -> dict:
    """Return the dashboard and worker settings to merge into ``ray.init`` kwargs.

    ``include_dashboard`` is always passed: Ray's own default starts the
    dashboard whenever its dependencies are installed.
    """
    options = {
        "include_dashboard": include_dashboard,
        "runtime_env": ray_worker_runtime_env(),
    }
    if include_dashboard:
        options["dashboard_host"] = ray_dashboard_host()
    return options
