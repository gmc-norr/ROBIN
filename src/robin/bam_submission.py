"""Which files a submission holds, and whether it is a fail-only BAM submission.

Shared by the Ray runner (``workflow_ray``) and the non-Ray runner (``workflow_simple``), so
both mark fail-only submissions by the same two rules (the review's M6):

- a batch of existing files is fail-only when it holds BAMs, all of them "fail" and none
  "pass" (``detect_fail_only_bam_submission``);
- in a watch session, a "fail" BAM that arrives before any "pass" BAM is fail-only
  (``FailOnlyWatchSession``).

A fail-only submission's classifier and bed_conversion failures are expected (job_outcome).
Nothing here imports Ray.
"""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional

from robin.job_outcome import bam_filename_kind


def matches_any_pattern(p: Path, patterns: Optional[List[str]]) -> bool:
    if not patterns:
        return True
    try:
        return any(p.match(pat) for pat in patterns)
    except Exception:
        return True


def matches_no_ignores(p: Path, ignore_patterns: Optional[List[str]]) -> bool:
    # Always ignore hidden dotfiles/directories by default.
    # This keeps `robin workflow` from getting noisy when watch directories contain
    # files like `.DS_Store` or editor artifacts.
    if p.name.startswith("."):
        return False
    if not ignore_patterns:
        return True
    try:
        return not any(p.match(pat) for pat in ignore_patterns)
    except Exception:
        return True


def is_bam_path(p: Path) -> bool:
    try:
        return p.suffix.lower() == ".bam"
    except Exception:
        return False


def detect_fail_only_bam_submission(
    paths: List[str],
    patterns: Optional[List[str]],
    ignore_patterns: Optional[List[str]],
    recursive: bool,
) -> bool:
    """
    True iff this submission contains BAMs and *all* BAM filenames contain 'fail'
    and none contain 'pass'. Used to suppress expected downstream failures on
    fail-only read sets.
    """
    bam_total = 0
    bam_fail = 0
    bam_pass = 0
    for p in paths:
        pth = Path(p)
        if pth.is_file():
            files = [pth]
        elif pth.is_dir():
            try:
                files = list(pth.rglob("*") if recursive else pth.glob("*"))
            except Exception:
                files = []
        else:
            files = []
        for f in files:
            try:
                if not f.is_file():
                    continue
                if not matches_any_pattern(f, patterns) or not matches_no_ignores(
                    f, ignore_patterns
                ):
                    continue
                if not is_bam_path(f):
                    continue
                bam_total += 1
                k = bam_filename_kind(str(f))
                if k == "fail":
                    bam_fail += 1
                elif k == "pass":
                    bam_pass += 1
            except Exception:
                continue
    return bool(bam_total and bam_pass == 0 and bam_fail == bam_total)


class FailOnlyWatchSession:
    """The watch-session rule: a "fail" BAM seen before any "pass" BAM is fail-only.

    Once a "pass" BAM has been seen, fail-BAM errors stay visible for the rest of the session.
    """

    def __init__(self) -> None:
        self.seen_bam = False
        self.seen_pass_bam = False

    def is_fail_only(self, filepath: str) -> bool:
        """Note this file, and say whether its jobs belong to a fail-only submission."""
        try:
            if not is_bam_path(Path(filepath)):
                return False
            self.seen_bam = True
            kind = bam_filename_kind(filepath)
            if kind == "pass":
                self.seen_pass_bam = True
            return kind == "fail" and not self.seen_pass_bam
        except Exception:
            return False
