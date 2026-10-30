"""The GUI's NiceGUI storage lives under ~/.robin, not in the folder Robin starts from."""

import os
import subprocess
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"


def _import_robin_and_print_storage_path(extra_env, home):
    env = {k: v for k, v in os.environ.items() if k != "NICEGUI_STORAGE_PATH"}
    env.update({"PYTHONPATH": str(SRC), "HOME": str(home)}, **extra_env)
    result = subprocess.run(
        [sys.executable, "-c", "import os, robin; print(os.environ['NICEGUI_STORAGE_PATH'])"],
        env=env,
        cwd=home,
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip()


def test_default_storage_path_is_under_home(tmp_path):
    out = _import_robin_and_print_storage_path({}, home=tmp_path)
    assert out == str(tmp_path / ".robin" / "nicegui")


def test_existing_storage_path_is_kept(tmp_path):
    out = _import_robin_and_print_storage_path(
        {"NICEGUI_STORAGE_PATH": "/custom/place"}, home=tmp_path
    )
    assert out == "/custom/place"
