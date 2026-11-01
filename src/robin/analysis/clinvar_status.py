"""Record, and read back, whether the ClinVar annotation step ran (VAR-4).

The ClinVar step writes one header line into the annotated VCF, right after
``##fileformat=``::

    ##ROBIN_ClinVar=annotated
    ##ROBIN_ClinVar=unavailable;reason=<short reason, no newline>

A file with no such line (made before this was recorded) reads as
``not_recorded``. Nothing here imports pandas or the heavy analysis modules,
so the report and the GUI can use it freely.
"""

from __future__ import annotations

import os
from typing import Tuple

CLINVAR_HEADER_KEY = "##ROBIN_ClinVar="
ANNOTATED = "annotated"
UNAVAILABLE = "unavailable"
NOT_RECORDED = "not_recorded"

_MAX_REASON = 200


def clean_reason(reason: str) -> str:
    """One line, no tabs, bounded length: safe to put in a VCF header."""
    text = " ".join(str(reason).split())
    return text[:_MAX_REASON]


def clinvar_header_line(status: str, reason: str = "") -> str:
    """The header line (without the newline) for ``status``."""
    if status == ANNOTATED:
        return f"{CLINVAR_HEADER_KEY}{ANNOTATED}"
    return f"{CLINVAR_HEADER_KEY}{UNAVAILABLE};reason={clean_reason(reason)}"


def clinvar_status(vcf_path) -> Tuple[str, str]:
    """Return ``(status, reason)`` for an annotated VCF.

    ``status`` is ``"annotated"``, ``"unavailable"`` (``reason`` says why) or
    ``"not_recorded"`` (no readable ``##ROBIN_ClinVar`` line, a missing file, or
    a value this code does not know). ``reason`` is ``""`` otherwise.
    """
    try:
        with open(os.fspath(vcf_path), "r", errors="replace") as handle:
            for line in handle:
                if not line.startswith("##"):
                    break
                if not line.startswith(CLINVAR_HEADER_KEY):
                    continue
                value = line[len(CLINVAR_HEADER_KEY):].rstrip("\r\n")
                if value == ANNOTATED:
                    return ANNOTATED, ""
                if value == UNAVAILABLE or value.startswith(UNAVAILABLE + ";"):
                    reason = ""
                    if ";reason=" in value:
                        reason = value.split(";reason=", 1)[1]
                    return UNAVAILABLE, reason
                return NOT_RECORDED, ""
    except OSError:
        pass
    return NOT_RECORDED, ""


def describe_not_annotated(kind: str, status: str, reason: str = "") -> str:
    """Sentence for the report/GUI when ``kind`` ('SNP'/'INDEL') was not annotated."""
    if status == UNAVAILABLE:
        why = f" ({reason})" if reason else ""
        return (
            f"ClinVar annotation was unavailable for the {kind} file{why}; "
            f"pathogenic {kind} variants could not be assessed."
        )
    return (
        f"ClinVar annotation status was not recorded for the {kind} file "
        f"(it predates status recording); pathogenic {kind} variants could not be assessed."
    )
