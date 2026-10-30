---
title: "Four ways a test in this repo passes without testing anything, and how to tell"
author: Pär
date: 2026-10-30
branch: fix/batch-3-quick-wins
cost: { minutes: 190, tokens: 2748 }
id: l-ce6fa244
scope: repo
status: active
src: agent
touches:
  - tests/test_classification_refresh.py
  - tests/test_clinvar_index.py
  - tests/test_completion_timeout.py
  - tests/test_classifier_config_keys.py
---

## Decided

A new test here is not believed until it has been seen failing against the old code, and the old
code is obtained with `git archive <sha> src | tar -x` into a scratchpad, never with `git stash`,
`git checkout --`, `git reset --hard` or `git clean` (the stash is shared with other sessions and
the rest destroy work). Several of this batch's tests went further and mutated the one edited line
in a scratchpad copy of `src` to confirm which test reddens; that is cheap and it is what caught
two tests that were pinning the wrong thing.

Judge a run on its `passed` count, never on the absence of `FAILED`. A test file whose import
fails prints `N error(s) during collection` and `Interrupted` with no `passed` line at all, which
reads like success at a glance.

## Tried and ruled out

**A stale-tabix test with records at positions 100 and 200 does not work.** A tabix index
addresses records by 16 kb bin, so both positions land in the same bin, the record's virtual
offset is unchanged, and a deliberately *stale* index still returns the new record. Measured:
same-bin stale fetch returns `['200']`, identical to a rebuilt index. Move the replacement record
across a bin boundary — position 100 replaced by 5,000,000, queried at 4,999,000-5,001,000 —
and the stale index returns `[]` while a rebuilt one returns `['5000000']`. The shell reproduction
script avoided this only by accident, by writing 200,000 records so the offsets drifted.

**Asserting that a rebuilt index's mtime changed is flaky, not strict.** Both filesystems here
store mtime to about 1 ms and a rebuild often lands inside the same millisecond: 122 failures in
300 runs, independently re-measured at 123/300. Backdate the surviving file by an hour with
`os.utime` first — which is also what the real situation looks like, an index written before the
file it indexes. 0/300 after.

**Stubbing NiceGUI in `sys.modules` inside the pytest session poisons other tests.** The real
NiceGUI 3.18 is installed and several test files import GUI modules for real. A GUI flow that
needs stubs has to run in a subprocess; `reviews/2026-10-16-repro/gui/stubs2/` is the stub set,
and `tests/test_classification_refresh.py` drives it that way.

**A tested helper is not a tested fix.** Twice in this batch a new module-level function had
thorough unit tests while the one line that *calls* it had none, so reverting that line left the
whole suite green. Both times the unit tests could only ever redden by `ImportError`. If you
extract a helper to make something testable, add the test that fails when the call site is
reverted — not only the ones that fail when the helper is missing.

## Watch out

GUI code has no harness, so an import check proves a file parses and nothing more. Say in the
commit message which lines are covered only that way. Two modules cannot even be imported in the
sanctioned test environment — `robin.gui.components.mnpflex` needs `alive_progress`, and
`robin.reporting.sections.classification` needs `seaborn`; neither is installed and installing is
not allowed. Byte-compile the first; stub the second inside a subprocess.

In the code:

- `/home/parlar_ai/robin-batch3-run/t3.sh` is the only runner that uses the right environment.
- `tests/test_clinvar_index.py` — the far-bin swap and the `os.utime` backdating.
- `tests/test_classification_refresh.py` — the subprocess driver and the `stubs2` path.
- `tests/test_completion_timeout.py` — helper tests plus the two that construct a `GUILauncher`,
  which is the wiring the helper tests miss.
- `tests/test_classifier_config_keys.py` — a subprocess driver that records what the real call
  site passes, alongside the weaker source-text assertion it replaced.
