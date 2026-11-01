"""The state of a sample's fusion results, per candidate type (REP-6, M8).

Shared by the PDF report and the GUI's fusion panel, and by the writer for the failed-build
record. Nothing here imports ReportLab, NiceGUI or the analysis modules.

A candidate type is the target panel (``"master"``) or genome-wide (``"all"``). Its state:

- ``found``: its processed pickle can be read and holds candidates;
- ``unreadable``: the pickle can't be read, or its last build failed (the writer leaves a
  ``<pickle>.failed`` record next to it, removed by the next good build), or the pickle is
  missing while the type's candidate CSV has rows (the results exist but were not processed);
- ``none``: the pickle can be read and is empty, or it is missing and the step ran
  (``sv_count.txt``) with no candidate rows;
- ``not_run``: no trace of the step for this type.
"""

from __future__ import annotations

import logging
import os
import pickle
from typing import Any, Dict, Optional, Tuple

import pandas as pd

logger = logging.getLogger(__name__)

STATE_NOT_RUN = "not_run"
STATE_UNREADABLE = "unreadable"
STATE_NONE = "none"
STATE_FOUND = "found"

KINDS = ("master", "all")
PICKLE_NAMES = {
    "master": "fusion_candidates_master_processed.pkl",
    "all": "fusion_candidates_all_processed.pkl",
}
CSV_NAMES = {
    "master": "fusion_candidates_master.csv",
    "all": "fusion_candidates_all.csv",
}
FAILED_BUILD_SUFFIX = ".failed"

# What the GUI's fusion panel says for a type with no candidates to show.
NONE_FOUND_MESSAGE = "No fusion candidates were identified"
_MESSAGES = {
    ("master", STATE_NOT_RUN): "Target panel fusion analysis has not been run yet",
    ("master", STATE_UNREADABLE): "The target panel fusion results could not be read",
    ("all", STATE_NOT_RUN): "Genome-wide fusion analysis has not been run yet",
    ("all", STATE_UNREADABLE): "The genome-wide fusion results could not be read",
}


# Shown beside a type's candidates when they come from an earlier run: the pickle still holds
# an earlier good build's candidates and the latest build of that type failed.
_EARLIER_RESULT_NOTES = {
    "master": "The latest target-panel fusion build failed; these results are from an earlier run.",
    "all": "The latest genome-wide fusion build failed; these results are from an earlier run.",
}


def failed_build_marker(pickle_path: Any) -> str:
    """The record the writer leaves next to a pickle whose last build failed."""
    return os.fspath(pickle_path) + FAILED_BUILD_SUFFIX


def earlier_result_note(sample_dir: Any, kind: str, state: Optional[str]) -> Optional[str]:
    """The line to show beside a type's candidates when they are from an earlier run, else None.

    That is a type whose state is found while its latest build failed (its ``.failed`` record
    is there). Without candidates the record already makes the type "could not be read".
    """
    if state != STATE_FOUND:
        return None
    path = os.path.join(os.fspath(sample_dir), PICKLE_NAMES[kind])
    if os.path.exists(failed_build_marker(path)):
        return _EARLIER_RESULT_NOTES[kind]
    return None


def read_processed_pickle(path: Any) -> Tuple[pd.DataFrame, bool]:
    """Return (candidates, readable) for one processed pickle.

    candidates is the good-pair frame (empty when the pickle holds none). A pickle that
    cannot be read, or holds something else than the processed structure, is not readable.
    """
    try:
        with open(path, "rb") as f:
            processed_data = pickle.load(f)
        annotated_data = processed_data.get("annotated_data", pd.DataFrame())
        goodpairs = processed_data.get("goodpairs", pd.Series(dtype=bool))
        if not isinstance(annotated_data, pd.DataFrame):
            raise TypeError("annotated_data is not a DataFrame")
        # Only keep the good pairs (same as GUI does)
        if not annotated_data.empty and not goodpairs.empty:
            candidates = annotated_data[goodpairs]
        else:
            candidates = annotated_data
        logger.debug(
            f"Loaded fusion candidates from {path} "
            f"({len(candidates)} good pairs from {len(annotated_data)} total candidates)"
        )
        return candidates, True
    except Exception as e:
        logger.warning(f"Could not read fusion pickle {path}: {e}")
        return pd.DataFrame(), False


def _missing_pickle_state(sample_dir: str, kind: str) -> str:
    """State of a type whose processed pickle is missing (and whose last build did not fail).

    CSV rows but no pickle: the results exist but were not processed (a folder from before
    the writer always made one, or this type's preprocessing failed). No rows: with
    ``sv_count.txt`` the step ran and nothing passed; without it, never run.
    """
    path = os.path.join(sample_dir, CSV_NAMES[kind])
    if os.path.exists(path):
        try:
            if not pd.read_csv(path, nrows=1).empty:
                return STATE_UNREADABLE
        except pd.errors.EmptyDataError:
            pass
        except Exception as e:
            logger.warning(f"Could not read fusion candidates {path}: {e}")
            return STATE_UNREADABLE
    if os.path.exists(os.path.join(sample_dir, "sv_count.txt")):
        return STATE_NONE
    return STATE_NOT_RUN


def candidate_type_state(
    sample_dir: Any, kind: str, pickle_read: Optional[Tuple[bool, bool]] = None
) -> str:
    """State of one candidate type (``"master"`` or ``"all"``) in ``sample_dir``.

    ``pickle_read``: ``(readable, has_candidates)`` when the caller has already read the
    processed pickle; otherwise it is read here.
    """
    sample_dir = os.fspath(sample_dir)
    path = os.path.join(sample_dir, PICKLE_NAMES[kind])
    failed = os.path.exists(failed_build_marker(path))
    if os.path.exists(path):
        if pickle_read is None:
            candidates, readable = read_processed_pickle(path)
            pickle_read = (readable, readable and not candidates.empty)
        readable, has_candidates = pickle_read
        if readable and has_candidates:
            return STATE_FOUND  # (from an earlier good build if the last one failed)
        if failed or not readable:
            return STATE_UNREADABLE
        return STATE_NONE
    if failed:
        return STATE_UNREADABLE
    return _missing_pickle_state(sample_dir, kind)


def sample_state(master_state: str, all_state: str) -> str:
    """The sample's state: found if either type is, else unreadable, else none, else not run."""
    states = {master_state, all_state}
    for state in (STATE_FOUND, STATE_UNREADABLE, STATE_NONE):
        if state in states:
            return state
    return STATE_NOT_RUN


def fusion_states(sample_dir: Any) -> Dict[str, Any]:
    """Read both processed pickles and decide each type's state and the sample's.

    Keys: ``master_candidates``/``all_candidates`` (the good-pair frame, None when the pickle
    could not be read or is missing), ``master_path``/``all_path``, ``master_state``/
    ``all_state``, ``master_note``/``all_note`` (``earlier_result_note``), ``state`` and
    ``unreadable`` (pickles that exist but can't be read).
    """
    sample_dir = os.fspath(sample_dir)
    out: Dict[str, Any] = {"unreadable": []}
    for kind in KINDS:
        path = os.path.join(sample_dir, PICKLE_NAMES[kind])
        out[f"{kind}_path"] = path
        out[f"{kind}_candidates"] = None
        pickle_read = None
        if os.path.exists(path):
            candidates, readable = read_processed_pickle(path)
            if readable:
                out[f"{kind}_candidates"] = candidates
            else:
                out["unreadable"].append(path)
            pickle_read = (readable, readable and not candidates.empty)
        out[f"{kind}_state"] = candidate_type_state(sample_dir, kind, pickle_read)
        out[f"{kind}_note"] = earlier_result_note(sample_dir, kind, out[f"{kind}_state"])
    out["state"] = sample_state(out["master_state"], out["all_state"])
    return out


def status_message(kind: str, state: Optional[str]) -> str:
    """What the GUI's fusion panel says for a type with nothing to show."""
    if state in (STATE_FOUND, STATE_NONE):
        return NONE_FOUND_MESSAGE
    return _MESSAGES.get((kind, state or STATE_UNREADABLE), _MESSAGES[(kind, STATE_UNREADABLE)])
