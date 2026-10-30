---
title: "Batch 1: Ray dashboard, NiceGUI storage, shell-safe commands, version pins"
author: Pär
date: 2026-10-30
branch: fix/batch-1-security-environment
cost: { minutes: 29, tokens: 282556 }
id: l-05340372
scope: repo
status: active
src: agent
touches:
  - CHANGELOG.md
  - src/robin/cli.py
  - src/robin/command_template.py
  - src/robin/analysis/mgmt_analysis.py
  - tests/test_mgmt_rscript.py
  - src/robin/__init__.py
  - tests/test_nicegui_storage_path.py
  - docs/getting-started/installation.md
  - pyproject.toml
  - robin.yml
---

## Decided
- All ray.init calls take robin.ray_setup.ray_init_options(include_dashboard). Always pass include_dashboard: Ray's default (None) starts the dashboard whenever ray[default] is installed, so leaving it out never turned it off.
- Ray workers get NICEGUI_STORAGE_PATH = a fresh empty temp folder via runtime_env. NiceGUI 3.18 rewrites storage-general.json at import through one fixed .tmp name; with no file and an empty dict it writes nothing.
- robin/__init__.py sets NICEGUI_STORAGE_PATH to ~/.robin/nicegui (only if unset) and creates it. NiceGUI's mkdir has no parents=True, and it logs instead of raising when a write fails.
- Command templates ({file}) go through robin.command_template.run_file_command: {file} becomes a quoted ${ROBIN_FILE} reference; the path travels in the environment.
- Pins: numpy=1.26.4, pandas=2.1.4, scikit-learn=1.8.0 (the NN model pickles record 1.8.0), Python >=3.12,<3.13.

## Tried and ruled out
- shlex.quote pasted into the template string: "{file}.bai" or "pre {file} post" inside double quotes stays injectable.

## Watch out
- {file} is only safe in the shell Robin starts. sh -c / eval / ssh templates re-parse the expanded name; use "$ROBIN_FILE" in the inner command.
- Installing scipy/matplotlib into a test env can pull numpy 2 and break pandas 2.1.4 ("numpy.dtype size changed"). Re-pin numpy==1.26.4 in the same install.
- NiceGUI storage tests must check the file on disk: NiceGUI swallows write errors.
