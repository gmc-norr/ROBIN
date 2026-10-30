"""SEC-1 and the NiceGUI storage race: how Robin starts Ray."""

import os
import re
import subprocess
import sys
import types
from pathlib import Path

import pytest

from robin.ray_setup import ray_dashboard_host, ray_init_options, ray_worker_runtime_env

SRC = Path(__file__).resolve().parents[1] / "src" / "robin"


def test_dashboard_host_defaults_to_localhost(monkeypatch):
    monkeypatch.delenv("RAY_DASHBOARD_HOST", raising=False)
    assert ray_dashboard_host() == "127.0.0.1"


def test_dashboard_host_honours_env(monkeypatch):
    monkeypatch.setenv("RAY_DASHBOARD_HOST", "0.0.0.0")
    assert ray_dashboard_host() == "0.0.0.0"


def test_worker_runtime_env_points_nicegui_at_empty_dir():
    path = ray_worker_runtime_env()["env_vars"]["NICEGUI_STORAGE_PATH"]
    assert os.path.isdir(path) and os.listdir(path) == []


def test_dashboard_is_off_unless_asked_for():
    # Ray starts the dashboard by default (include_dashboard=None), so "off" has to be said.
    options = ray_init_options(include_dashboard=False)
    assert options["include_dashboard"] is False
    assert "dashboard_host" not in options


def test_importing_nicegui_with_empty_storage_dir_writes_nothing(tmp_path):
    # NiceGUI 3.18 rewrites storage-general.json at import. The worker fix relies on
    # it writing nothing when its storage folder is empty.
    pytest.importorskip("nicegui")
    gui_dir = tmp_path / "gui"
    gui_dir.mkdir()
    gui_file = gui_dir / "storage-general.json"
    gui_file.write_text('{"k": 1}')
    before = gui_file.stat().st_mtime_ns
    worker_dir = tmp_path / "worker"
    worker_dir.mkdir()
    env = {**os.environ, "NICEGUI_STORAGE_PATH": str(worker_dir)}
    subprocess.run([sys.executable, "-c", "import nicegui"], env=env, check=True, cwd=gui_dir)
    assert list(worker_dir.iterdir()) == []
    assert gui_file.stat().st_mtime_ns == before


class _FakeRay(types.ModuleType):
    def __init__(self):
        super().__init__("ray")
        self.init_kwargs = None

    def is_initialized(self):
        return False

    def init(self, **kwargs):
        self.init_kwargs = kwargs


@pytest.fixture
def fake_ray(monkeypatch):
    ray = _FakeRay()
    monkeypatch.setitem(sys.modules, "ray", ray)
    monkeypatch.delenv("RAY_DASHBOARD_HOST", raising=False)
    return ray


def test_cli_starts_ray_without_dashboard_by_default(fake_ray):
    from robin.cli import _initialize_ray

    _initialize_ray(2)
    assert fake_ray.init_kwargs["include_dashboard"] is False
    assert "dashboard_host" not in fake_ray.init_kwargs
    storage = fake_ray.init_kwargs["runtime_env"]["env_vars"]["NICEGUI_STORAGE_PATH"]
    assert os.listdir(storage) == []


def test_cli_dashboard_when_asked_listens_on_localhost(fake_ray):
    from robin.cli import _initialize_ray

    _initialize_ray(2, include_dashboard=True)
    assert fake_ray.init_kwargs["include_dashboard"] is True
    assert fake_ray.init_kwargs["dashboard_host"] == "127.0.0.1"


def test_cli_workflow_dashboard_flag_defaults_off():
    from robin.cli import workflow

    defaults = {p.name: p.default for p in workflow.params}
    assert defaults["ray_dashboard"] is False


def test_workflow_ray_dashboard_flag_defaults_off(monkeypatch):
    from robin.workflow_ray import parse_args

    monkeypatch.setattr(sys, "argv", ["workflow_ray", "--center", "X"])
    assert parse_args().ray_dashboard is False
    monkeypatch.setattr(sys, "argv", ["workflow_ray", "--center", "X", "--ray-dashboard"])
    assert parse_args().ray_dashboard is True


def test_no_ray_dashboard_binds_all_interfaces_by_default():
    # Guard for the two ray.init sites a unit test can't reach (cli `workflow` and
    # workflow_ray's __main__): every dashboard bind goes through ray_dashboard_host().
    pattern = re.compile(r"""RAY_DASHBOARD_HOST["']\s*,\s*["']0\.0\.0\.0""")
    offenders = [p.name for p in SRC.rglob("*.py") if pattern.search(p.read_text())]
    assert offenders == []
