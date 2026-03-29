"""
M4 – COI Extractability and Species Name Validation

Performs two checks on the deposited organism identity:

1. **COI extractability** – attempts to extract the cytochrome c oxidase
   subunit I (COI/cox1) sequence from the annotated CDS features.  COI
   is the universal DNA barcoding locus (Hebert et al. 2003 Proc R Soc B
   270:313); its absence from a complete mitogenome annotation is a
   significant quality concern.

2. **Species name validation** – checks that the organism name field
   contains a valid binomial (genus + species epithet) and is not an
   ambiguous placeholder (e.g. "sp.", "cf.", "aff.", numeric strain IDs).

THRESHOLDS
----------
- COI minimum length: 600 bp (half-barcode; full barcode is 658 bp)
- Ambiguous name patterns: "sp.", "cf.", "aff.", "nr.", "gen. nov.",
  "sp. nov.", and organism names consisting of only one word.

NOTES
-----
This module (v1) performs local checks only.  MitoQC v0.2 introduces
``m4_species_v2.py`` with an optional NCBI BLAST verification layer
(disabled by default; requires network access and an NCBI API key for
high-throughput use).
"""

from __future__ import annotations

import re

from mitoqc.utils.constants import COI_CANONICAL_NAMES, GENE_NORM

_COI_MIN_LEN = 600   # bp

# Patterns that indicate an ambiguous/unresolved species name
_AMBIG_PATTERNS = re.compile(
    r"\bsp\.\b|\bcf\.\b|\baff\.\b|\bnr\.\b|"
    r"\bgen\.\s*nov\.\b|\bsp\.\s*nov\.\b",
    re.IGNORECASE,
)
_NUMERIC_STRAIN = re.compile(r"\b\d{4,}\b")   # 4+ digit number in name


def _normalise(name: str) -> str:
    return GENE_NORM.get(name.lower().strip(), name.lower().strip())


def _is_coi(gene_name: str) -> bool:
    return _normalise(gene_name) in COI_CANONICAL_NAMES


def _extract_coi(rec: dict) -> tuple[str | None, int]:
    """Return (coi_sequence, length) or (None, 0) if not found."""
    seq = rec["seq"]
    for feat in rec["features"]:
        if feat["type"] != "CDS":
            continue
        if not _is_coi(feat["gene"]):
            continue
        start, end, strand = feat["start"], feat["end"], feat["strand"]
        if strand == -1:
            from Bio.Seq import Seq
            coi_seq = str(Seq(seq[start:end]).reverse_complement())
        else:
            coi_seq = seq[start:end]
        return coi_seq, len(coi_seq)
    return None, 0


def run(rec: dict, config: dict | None = None) -> dict:
    """Run M4 COI extractability and species name checks.

    Parameters
    ----------
    rec:
        Parsed record dict.  Must include ``features``, ``seq``,
        and ``organism``.
    config:
        Optional overrides:
        - ``coi_min_len`` (int): minimum COI length in bp (default 600).

    Returns
    -------
    dict
        Keys: ``m4_coi_found``, ``m4_coi_len``, ``m4_organism``,
        ``m4_flags``, ``m4_pass``.
    """
    cfg         = config or {}
    coi_min_len = cfg.get("coi_min_len", _COI_MIN_LEN)
    organism    = rec.get("organism", "")

    coi_seq, coi_len = _extract_coi(rec)
    coi_found = coi_seq is not None and coi_len >= coi_min_len

    flags: list[str] = []

    # COI check
    if coi_seq is None:
        flags.append("COI_NOT_EXTRACTED")
    elif coi_len < coi_min_len:
        flags.append(f"COI_TOO_SHORT({coi_len}<{coi_min_len})")

    # Species name check
    words = organism.strip().split()
    if len(words) < 2:
        flags.append(f"NON_BINOMIAL_NAME('{organism}')")
    elif _AMBIG_PATTERNS.search(organism):
        flags.append(f"AMBIGUOUS_TAXON('{organism}')")
    elif _NUMERIC_STRAIN.search(organism):
        flags.append(f"NUMERIC_STRAIN_ID('{organism}')")

    return {
        "m4_coi_found": coi_seq is not None,
        "m4_coi_len":   coi_len,
        "m4_organism":  organism,
        "m4_flags":     "|".join(flags) if flags else "PASS",
        "m4_pass":      len(flags) == 0,
    }
