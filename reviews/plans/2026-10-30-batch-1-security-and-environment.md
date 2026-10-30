# Batch 1: Security and Environment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close the two shell/network code-execution holes (SEC-1, SEC-8, plus SEC-9 found while planning), stop the NiceGUI storage-file race between Ray workers and the GUI, and pin the versions that have broken the clinical server's environment three times.

**Architecture:** Two tiny new modules carry the shared logic, so no fix is copied across call sites: `robin/ray_setup.py` (dashboard host and Ray worker runtime environment, used by all three `ray.init` sites) and `robin/command_template.py` (safe `{file}` substitution, used by both command runners). The MGMT fix swaps one shell string for an argument list in place. Pins live in `robin.yml` and `pyproject.toml`.

**Tech Stack:** Python 3.12, Ray, NiceGUI 3.18, click/argparse, pytest.

**Spec:** [`reviews/2026-10-16-code-review.md`](../2026-10-16-code-review.md) sections SEC-1, SEC-8, PKG-3, PKG-6, ARCH-5. SEC-9 and the NiceGUI race were found after the review; their evidence is in Tasks 4 and 1.

**Outcome (2026-10-30):** done and merged to `main` (`3623e93`..`75485ce`). Task 4 was built differently from this plan: `{file}` becomes a quoted `$ROBIN_FILE` reference and the path travels in the environment (`robin.command_template.run_file_command`), because pasting `shlex.quote(path)` into the template leaves `"{file}.bai"` injectable. The final review added three fixes: the default storage folder is created, an Rscript that can't start warns and continues, and the SEC-9 docs say where the protection stops.

**Batch order:** this is batch 1 of the review's "Suggested order of work". Batch 2 is sample identity (VAR-9, VAR-8, VAR-10), and gets its own plan.

## Global Constraints

- Python `>=3.12,<3.13`; `numpy=1.26.4`, `pandas=2.1.4`, `scikit-learn=1.8.0` (both NN model pickles record `_sklearn_version` 1.8.0).
- Run tests in a Python 3.12 env that has the repo's pins: `PYTHONPATH=src python -m pytest <path> -p no:cacheprovider --no-cov -q` from the repo root.
- No new runtime dependencies.
- Don't change behaviour beyond each task's stated fix. Each task is its own commit, with a message that names its finding ID.
- Work on branch `fix/batch-1-security-environment` in a git worktree, so the uncommitted merged review on `main` is untouched.

## Review Focus

1. **Paths with spaces or shell characters** (`$(…)`, `;`, quotes) as a sample folder or watched file: they must reach the program as one literal argument, and nothing in them may run. Pinned in Task 3 (`test_sample_dir_with_shell_syntax_is_one_literal_argument`) and Task 4 (`test_filename_with_shell_syntax_does_not_run`).
2. **User command templates that already quote `{file}`** (`cp "{file}" /x`, `cp '{file}' /x`) must keep working. Pinned in Task 4 (`test_quoted_placeholders_are_replaced_whole`).
3. **An operator who sets `RAY_DASHBOARD_HOST=0.0.0.0` on purpose** must still get that. Pinned in Task 1 (`test_dashboard_host_honours_env`).
4. **`Rscript` missing from PATH:** the old shell call returned 127 and MGMT logged a warning and carried on. An argument-list call raises `FileNotFoundError`, which would turn that into a failed MGMT step. Keep the warning-and-continue behaviour. Pinned in Task 3 (`test_missing_rscript_logs_warning_and_continues`).
5. **A NICEGUI_STORAGE_PATH someone already set** (and the value Ray gives workers) must not be overridden by Robin's new default. Pinned in Task 2 (`test_existing_storage_path_is_kept`).

---

### Task 1: Ray dashboard opt-in on localhost, and an empty NiceGUI storage folder for Ray workers (SEC-1, NiceGUI race)

**Why:** The dashboard's Jobs API runs arbitrary commands and binds `0.0.0.0` with no auth (SEC-1). Separately, NiceGUI 3.18 calls `app.storage.general.initialize_sync()` at import (`nicegui/nicegui.py:51`). That rewrites `.nicegui/storage-general.json` through a fixed `storage-general.json.tmp`, and every Ray worker imports NiceGUI, so they collide (`FileNotFoundError … storage-general.json.tmp`). `FilePersistentDict.backup` returns without writing when the file is absent and the dict is empty, so pointing workers at an empty folder means they never write.

**Files:**
- Create: `src/robin/ray_setup.py`
- Modify: `src/robin/cli.py:1674-1747` (`_initialize_ray`), `src/robin/cli.py:2302-2306` (option), `src/robin/cli.py:2429-2454` (`workflow`)
- Modify: `src/robin/workflow_ray.py:5124-5128` (argparse), `src/robin/workflow_ray.py:5152-5169` (`ray.init`)
- Test: `tests/test_ray_setup.py`

**Interfaces:**
- Produces: `ray_dashboard_host() -> str` and `ray_worker_runtime_env() -> dict` in `robin.ray_setup`; constant `DEFAULT_DASHBOARD_HOST = "127.0.0.1"`.

- [ ] **Step 1: Install NiceGUI 3.18.0 into the test env** (`pip install nicegui==3.18.0`), so the NiceGUI test below runs instead of skipping.

- [ ] **Step 2: Write the failing tests** in `tests/test_ray_setup.py`:

```python
def test_dashboard_host_defaults_to_localhost(monkeypatch):
    monkeypatch.delenv("RAY_DASHBOARD_HOST", raising=False)
    assert ray_dashboard_host() == "127.0.0.1"

def test_dashboard_host_honours_env(monkeypatch):
    monkeypatch.setenv("RAY_DASHBOARD_HOST", "0.0.0.0")
    assert ray_dashboard_host() == "0.0.0.0"

def test_worker_runtime_env_points_nicegui_at_empty_dir():
    path = ray_worker_runtime_env()["env_vars"]["NICEGUI_STORAGE_PATH"]
    assert os.path.isdir(path) and os.listdir(path) == []

def test_importing_nicegui_with_empty_storage_dir_writes_nothing(tmp_path):
    pytest.importorskip("nicegui")
    gui_dir = tmp_path / "gui"; gui_dir.mkdir()
    gui_file = gui_dir / "storage-general.json"; gui_file.write_text('{"k": 1}')
    before = gui_file.stat().st_mtime_ns
    worker_dir = tmp_path / "worker"; worker_dir.mkdir()
    env = {**os.environ, "NICEGUI_STORAGE_PATH": str(worker_dir)}
    subprocess.run([sys.executable, "-c", "import nicegui"], env=env, check=True, cwd=gui_dir.parent)
    assert list(worker_dir.iterdir()) == []
    assert gui_file.stat().st_mtime_ns == before

def test_no_ray_dashboard_binds_all_interfaces_by_default():
    # guard: every dashboard bind goes through ray_dashboard_host()
    pattern = re.compile(r"""RAY_DASHBOARD_HOST["']\s*,\s*["']0\.0\.0\.0""")
    offenders = [p for p in Path("src/robin").rglob("*.py") if pattern.search(p.read_text())]
    assert offenders == []
```

- [ ] **Step 3: Run them and confirm they fail.** `ImportError: robin.ray_setup` for the first three. The last test lists `cli.py` and `workflow_ray.py`. The NiceGUI test passes already, because it documents NiceGUI's behaviour that the fix relies on.

- [ ] **Step 4: Create `src/robin/ray_setup.py`** with `DEFAULT_DASHBOARD_HOST`, `ray_dashboard_host()` (returns `RAY_DASHBOARD_HOST` if set, else the default) and `ray_worker_runtime_env()` (makes a fresh `tempfile.mkdtemp(prefix="robin-ray-nicegui-")` and returns `{"env_vars": {"NICEGUI_STORAGE_PATH": <it>}}`). Each docstring says why.

- [ ] **Step 5: Use it at all three `ray.init` sites.**
  - Replace each `os.environ.get("RAY_DASHBOARD_HOST", "0.0.0.0")` with `ray_dashboard_host()`.
  - Add `init_kwargs["runtime_env"] = ray_worker_runtime_env()`. In `workflow_ray.py`, also pass it in the `TypeError` fallback `ray.init(address=…)`.
  - Fix the comments that say "0.0.0.0" or "all interfaces".

- [ ] **Step 6: Make the dashboard opt-in.**
  - `cli.py`: change the `--ray-dashboard/--no-ray-dashboard` default to `False`, with help text "Enable the Ray dashboard on localhost (default: off). Set RAY_DASHBOARD_HOST to bind elsewhere."
  - `cli.py`: change the `_initialize_ray(..., include_dashboard: bool = False)` default.
  - `workflow_ray.py`: replace `--no-ray-dashboard` with `p.add_argument("--ray-dashboard", action=argparse.BooleanOptionalAction, default=False, help=…)`, and read `args.ray_dashboard`.

- [ ] **Step 7: Run** `tests/test_ray_setup.py` and the whole `tests/` folder. Expected: all pass. Run `trail grep "0.0.0.0"` and confirm that none of the remaining hits is a Ray dashboard bind.

- [ ] **Step 8: Commit** with a message like "SEC-1: Ray dashboard off by default and localhost-only; give Ray workers an empty NiceGUI storage folder".

### Task 2: Keep the GUI's NiceGUI storage out of the working directory

**Why:** NiceGUI defaults its storage to `./.nicegui`, so starting Robin from the checkout writes login-session data into the repo (and `.nicegui/` isn't gitignored). NiceGUI reads `NICEGUI_STORAGE_PATH` once, when it's first imported, and `robin/__init__.py` runs before any NiceGUI import.

**Files:**
- Modify: `src/robin/__init__.py` (top, before `from . import analysis`)
- Modify: `.gitignore`
- Test: `tests/test_nicegui_storage_path.py`

**Interfaces:**
- Produces: env default `NICEGUI_STORAGE_PATH = ~/.robin/nicegui`, set with `os.environ.setdefault`.

- [ ] **Step 1: Write the failing tests.** The helper `_import_robin_and_print_storage_path(extra_env, home)` runs `[sys.executable, "-c", "import os, robin; print(os.environ['NICEGUI_STORAGE_PATH'])"]` with `PYTHONPATH=src`, `HOME=home`, `NICEGUI_STORAGE_PATH` removed from the inherited env, and `extra_env` added. It returns stripped stdout.

```python
def test_default_storage_path_is_under_home(tmp_path):
    out = _import_robin_and_print_storage_path({}, home=tmp_path)
    assert out == str(tmp_path / ".robin" / "nicegui")

def test_existing_storage_path_is_kept(tmp_path):
    out = _import_robin_and_print_storage_path({"NICEGUI_STORAGE_PATH": "/custom/place"}, home=tmp_path)
    assert out == "/custom/place"
```

- [ ] **Step 2: Run them.** Expected: the first fails with `KeyError: 'NICEGUI_STORAGE_PATH'`, and the second passes.

- [ ] **Step 3: Add the `os.environ.setdefault(...)` line to `src/robin/__init__.py`,** with a comment that it must come before anything imports NiceGUI. Add `.nicegui/` to `.gitignore`.

- [ ] **Step 4: Run all tests.** Expected: pass.

- [ ] **Step 5: Commit.** The message notes the one-time side effect: users see the disclaimer and log in again once, because the storage moves.

### Task 3: MGMT runs Rscript without a shell (SEC-8)

**Why:** `run_matkit` builds `f"Rscript {r_script_path} --input=… --out_dir={sample_dir} …"` and runs it with `shell=True`, so a sample folder name with `$(…)` executes (reproduced by `reviews/2026-10-16-repro/independent/reproduce.py::test_ir08`).

**Files:**
- Modify: `src/robin/analysis/mgmt_analysis.py:834-852`
- Test: `tests/test_mgmt_rscript.py`

**Interfaces:**
- Consumes: `mgmt_analysis.run_matkit(temp_dir, mgmt_bamfile, hv_path, sample_dir, file_number) -> bool`, unchanged.

- [ ] **Step 1: Write the failing tests.**
  - Stub `mgmt_analysis.safely_sort_and_index_bam` so it touches `output_bam` and returns True.
  - Put a stub `robin.analysis.utilities.matkit` module in `sys.modules` (monkeypatch) whose `run_matkit(in, out)` writes one 18-column line.
  - Create `hv/bin/mgmt_pred_v0.3.R`.
  - Put a fake `Rscript` on PATH (monkeypatch `PATH` to a tmp bin folder). It's a shell script that writes its arguments, one per line, to `$ARGS_FILE`.

```python
def test_sample_dir_with_shell_syntax_is_one_literal_argument(tmp_path, ...):
    sample_dir = tmp_path / "s $(touch MARKER) x"; sample_dir.mkdir()
    assert run_matkit(str(tmp_path), "in.bam", str(hv), str(sample_dir), 1) is True
    args = args_file.read_text().splitlines()
    assert f"--out_dir={sample_dir}" in args
    assert not (tmp_path / "MARKER").exists() and not Path("MARKER").exists()

def test_missing_rscript_logs_warning_and_continues(tmp_path, caplog, ...):
    # PATH has no Rscript
    assert run_matkit(...) is True
    assert "R script failed" in caplog.text
```

- [ ] **Step 2: Run them.** Expected: the first fails, because the args are split at the space and the shell runs `touch MARKER`. The second fails if the list call raises.

- [ ] **Step 3: Replace the string and `shell=True`** with the list `["Rscript", r_script_path, f"--input={output_mgmt_bed}", f"--out_dir={sample_dir}", f"--probes={os.path.join(hv_path, 'bin', 'mgmt_probes.Rdata')}", f"--model={os.path.join(hv_path, 'bin', 'mgmt_137sites_mean_model.Rdata')}", "--sample=live_analysis"]`. Wrap the call in `try/except FileNotFoundError` and log `R script failed: Rscript not found on PATH`, so it keeps the old warning-and-continue path.

- [ ] **Step 4: Run all tests, plus `reviews/2026-10-16-repro/independent/reproduce.py`.** Expected: `test_ir08` now **fails**, because the bug is gone. That's correct. Note it in the commit message.

- [ ] **Step 5: Commit** ("SEC-8: run the MGMT R script without a shell").

### Task 4: Quote watched file paths in user command templates (SEC-9)

**Why:** `watcher.FileChangeHandler._run_command` (`watcher.py:116`) and `workflow_simple.command_handler` (`workflow_simple.py:1756`) do `template.replace("{file}", path)` and run the result with `shell=True`. The template is the operator's own shell command, so a shell is intended, but the file name comes from the watched folder: a BAM named `a;rm -rf ~.bam` runs `rm`.

**Files:**
- Create: `src/robin/command_template.py`
- Modify: `src/robin/watcher.py:116`, `src/robin/workflow_simple.py:1756`
- Test: `tests/test_command_template.py`

**Interfaces:**
- Produces: `fill_file_placeholder(template: str, file_path: str) -> str`. It replaces every `{file}`, `"{file}"` and `'{file}'` with `shlex.quote(file_path)`.

- [ ] **Step 1: Write the failing tests.**

```python
def test_bare_placeholder_is_quoted():
    assert fill_file_placeholder("echo {file}", "/d/a b.bam") == "echo '/d/a b.bam'"

def test_quoted_placeholders_are_replaced_whole():
    assert fill_file_placeholder('cp "{file}" /x', "/d/a b.bam") == "cp '/d/a b.bam' /x"
    assert fill_file_placeholder("cp '{file}' /x", "/d/a b.bam") == "cp '/d/a b.bam' /x"

def test_filename_with_shell_syntax_does_not_run(tmp_path):
    path = str(tmp_path / "x$(touch MARKER);y.bam")
    out = subprocess.run(fill_file_placeholder("printf %s {file}", path), shell=True,
                         capture_output=True, text=True, cwd=tmp_path).stdout
    assert out == path and not (tmp_path / "MARKER").exists()
```

- [ ] **Step 2: Run them.** Expected: `ImportError`.

- [ ] **Step 3: Implement `fill_file_placeholder`** with one `re.sub(r"""(["']?)\{file\}\1""", lambda m: shlex.quote(file_path), template)`. Use it in both call sites, which keep `shell=True` on purpose (the template is the operator's shell command). Add a one-line comment saying so.

- [ ] **Step 4: Run all tests.** Expected: pass.

- [ ] **Step 5: Commit** ("SEC-9: quote watched file paths in command templates"). Add SEC-9 to the review under section 8 in the same commit.

### Task 5: Pin the versions the server needs

**Why:** The clinical server's env broke three ways: pip ran under Python 3.14; conda moved scikit-learn to 1.9.1 while the model pickles were saved with 1.8.0; and conda moved numpy to 2.5.4 under pip's pandas 2.1.4 ("numpy.dtype size changed"). conda-forge has `numpy 1.26.4`, `pandas 2.1.4` (requires `numpy <2`) and `scikit-learn 1.8.0` for linux-64 py312 (checked 2026-10-20).

**Files:**
- Modify: `robin.yml:28,32` (+ add `pandas`), `pyproject.toml:21,49,66`, `CHANGELOG.md` (`[Unreleased]`)
- Test: `tests/test_dependency_pins.py`

- [ ] **Step 1: Write the failing tests.**
  - Parse `pyproject.toml` with `tomllib`: `requires-python == ">=3.12,<3.13"`; `scikit-learn==1.8.0` appears exactly once, and there's no bare `scikit-learn`.
  - Read `robin.yml` line by line: it contains `- numpy=1.26.4`, `- pandas=2.1.4` and `- scikit-learn=1.8.0`.
  - Add a comment in the test pointing to the pickles' recorded version.

- [ ] **Step 2: Run them.** Expected: they fail on each missing pin.

- [ ] **Step 3: Apply the pins.** Remove the duplicate `scikit-learn` line from `pyproject.toml`. Next to the pins in `robin.yml`, add a comment saying why (model pickles, and pandas built for numpy 1.x).

- [ ] **Step 4: Add `[Unreleased]` entries to `CHANGELOG.md`** for all five tasks. Under **Security:** SEC-1, SEC-8, SEC-9. Under **Fixed:** the NiceGUI storage race, and the storage location move with its one-time re-login. Under **Changed:** the pins and the Python upper bound.

- [ ] **Step 5: Run all tests.** Expected: pass. Then commit ("Pin numpy, pandas, scikit-learn and Python 3.12 (PKG-3, PKG-6)").

---

**After the batch:** a fresh reviewer reviews the whole branch. Then push the branch to `gmc-norr` and open a PR to gmc-norr `main`. SEC-1, SEC-8 and SEC-9 also need a private note to LooseLab.
