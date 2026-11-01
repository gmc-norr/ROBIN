"""Whether ClinVar calls a record pathogenic (VAR-2), from ClinVar's own labels.

A record counts when ClinVar's aggregate label says so in any of its three classifications:

- ``CLNSIG`` (germline): ``Pathogenic`` or ``Likely_pathogenic``
- ``ONC`` (oncogenicity): ``Oncogenic`` or ``Likely_oncogenic``
- ``SCI`` (somatic clinical impact): ``Tier_I_-_Strong`` or ``Tier_II_-_Potential``

Labels are matched whole. ``Conflicting_classifications_of_pathogenicity`` contains the word
"pathogenic" but is not a pathogenic call, so it never counts by itself. A value is split into
single labels on ``|`` (several labels), ``/`` (``Pathogenic/Likely_pathogenic``) and ``,``
(``Pathogenic,_low_penetrance``).

The variant check, the GUI's SNP table and the report all use this. It imports nothing heavy.
"""

from __future__ import annotations

import re
from typing import Any, Dict, Mapping, Union

_COUNTING_LABELS = {
    "CLNSIG": frozenset({"pathogenic", "likely_pathogenic"}),
    "ONC": frozenset({"oncogenic", "likely_oncogenic"}),
    "SCI": frozenset({"tier_i_-_strong", "tier_ii_-_potential"}),
}
_DISPLAY_PREFIX = {"CLNSIG": "", "ONC": "oncogenicity: ", "SCI": "somatic impact: "}
_LABEL_SEPARATORS = re.compile(r"[|/,]")

Info = Union[str, Mapping[str, Any]]


def parse_info(info: str) -> Dict[str, str]:
    """``KEY=VALUE`` pairs of a VCF INFO string. Flags (no ``=``) are left out."""
    fields: Dict[str, str] = {}
    for part in info.split(";"):
        key, sep, value = part.partition("=")
        if sep:
            fields[key] = value
    return fields


def _value(info: Mapping[str, Any], key: str) -> str:
    value = info.get(key)
    if isinstance(value, (list, tuple)):
        return "|".join(str(item) for item in value)
    return value if isinstance(value, str) else ""


def _as_mapping(info: Info) -> Mapping[str, Any]:
    return parse_info(info) if isinstance(info, str) else info


def is_clinvar_pathogenic(info: Info) -> bool:
    """True when CLNSIG, ONC or SCI holds one of the labels that count (see the module doc)."""
    fields = _as_mapping(info)
    for key, counting in _COUNTING_LABELS.items():
        labels = {label.casefold() for label in _LABEL_SEPARATORS.split(_value(fields, key))}
        if labels & counting:
            return True
    return False


def clinvar_label_text(info: Info, default: str = "unknown") -> str:
    """ClinVar's labels for display.

    CLNSIG as written, then ``oncogenicity: <ONC>`` and ``somatic impact: <SCI>`` when present,
    joined by ``; ``. A record with only CLNSIG reads exactly as its CLNSIG.
    """
    fields = _as_mapping(info)
    parts = [
        prefix + _value(fields, key)
        for key, prefix in _DISPLAY_PREFIX.items()
        if _value(fields, key)
    ]
    return "; ".join(parts) or default
