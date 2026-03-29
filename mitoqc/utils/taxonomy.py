"""
taxonomy.py – Taxon group classification from NCBI lineage strings.

Maps the free-text lineage field from GenBank/RefSeq records to one of
the 12 canonical MitoQC taxon groups used for threshold selection.

Classification is keyword-based (order matters — more specific groups
are checked before broader ones).  The fallback is "Other_Metazoa".
"""

from __future__ import annotations

# Ordered list of (taxon_group, keywords_that_must_appear_in_lineage)
# Keywords are matched case-insensitively against the full lineage string.
_RULES: list[tuple[str, list[str]]] = [
    ("Mammalia",       ["Mammalia"]),
    ("Aves",           ["Aves"]),
    ("Reptilia",       ["Lepidosauria", "Testudines", "Crocodilia", "Squamata",
                        "Rhynchocephalia"]),
    ("Amphibia",       ["Amphibia"]),
    ("Actinopterygii", ["Actinopterygii"]),
    ("Insecta",        ["Insecta"]),
    ("Arachnida",      ["Arachnida"]),
    ("Crustacea",      ["Malacostraca", "Branchiopoda", "Maxillopoda",
                        "Ostracoda", "Crustacea"]),
    ("Mollusca",       ["Mollusca"]),
    ("Nematoda",       ["Nematoda"]),
    ("Echinodermata",  ["Echinodermata"]),
]


def classify_taxon(lineage: str) -> str:
    """Return the MitoQC taxon group for a given NCBI lineage string.

    Parameters
    ----------
    lineage:
        Semicolon-separated lineage string as returned by Biopython's
        ``record.annotations["taxonomy"]`` joined with ``"; "``.

    Returns
    -------
    str
        One of the 12 canonical taxon group names, or ``"Other_Metazoa"``
        if no rule matches.

    Examples
    --------
    >>> classify_taxon("Eukaryota; Metazoa; Chordata; Mammalia; Primates")
    'Mammalia'
    >>> classify_taxon("Eukaryota; Metazoa; Arthropoda; Insecta; Diptera")
    'Insecta'
    >>> classify_taxon("Eukaryota; Metazoa; Platyhelminthes; Trematoda")
    'Other_Metazoa'
    """
    lin_lower = lineage.lower()
    for group, keywords in _RULES:
        if any(kw.lower() in lin_lower for kw in keywords):
            return group
    return "Other_Metazoa"
