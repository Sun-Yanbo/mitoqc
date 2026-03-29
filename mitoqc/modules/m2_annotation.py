"""
M2 – Annotation Completeness Proxy

Checks whether the GenBank feature table contains the expected complement
of mitochondrial genes for a complete metazoan mitogenome:

  - 13 protein-coding genes (PCGs): nad1-6, nad4l, cox1-3, atp6, atp8, cytb
    (atp8 is exempt for Nematoda — biologically absent in this clade)
  - 2 ribosomal RNA genes: 12S rRNA (rrnS) and 16S rRNA (rrnL)
  - ≥18 tRNA genes (conservative lower bound; standard is 22)

This is a *proxy* check: it tests annotation completeness, not sequence
completeness.  A record can fail M2 because (a) the gene is genuinely
absent (assembly gap, NUMT), or (b) the gene is present but unannotated
or annotated under a non-standard alias.  Gene name normalisation via
``GENE_NORM`` mitigates (b) for common aliases.

NEMATODA EXCEPTION
------------------
Nematode mitogenomes lack atp8 (Lavrov & Brown 2001 Genetics 157:621).
MitoQC exempts atp8 from the expected PCG set for Nematoda records.

EXTRA PCG FLAG
--------------
Records with more than 13 annotated PCGs (after normalisation) are
flagged with EXTRA_PCGS.  This can indicate duplicated annotations,
NUMTs included in the feature table, or annotation errors.
"""

from __future__ import annotations

from mitoqc.utils.constants import (
    EXPECTED_PCGS, GROUPS_LACKING_ATP8, MIN_TRNA_COUNT, GENE_NORM
)

_RRNA_CANONICAL = frozenset({"rrnS", "rrnL"})


def _normalise(name: str) -> str:
    return GENE_NORM.get(name.lower().strip(), name.lower().strip())


def run(rec: dict, config: dict | None = None) -> dict:
    """Run M2 annotation completeness checks on a single record.

    Parameters
    ----------
    rec:
        Parsed record dict.  Must include ``features`` list and
        ``taxon_group`` string.
    config:
        Optional overrides:
        - ``min_trna`` (int): minimum tRNA count (default 18).

    Returns
    -------
    dict
        Keys: ``m2_n_pcg``, ``m2_n_rrna``, ``m2_n_trna``,
        ``m2_missing_pcgs``, ``m2_extra_pcgs``, ``m2_flags``, ``m2_pass``.
    """
    cfg      = config or {}
    group    = rec["taxon_group"]
    min_trna = cfg.get("min_trna", MIN_TRNA_COUNT)

    # Build expected PCG set (exempt atp8 for Nematoda)
    expected = set(EXPECTED_PCGS)
    if group in GROUPS_LACKING_ATP8:
        expected.discard("atp8")

    # Collect annotated genes by type
    pcg_found:  set[str] = set()
    rrna_found: set[str] = set()
    n_trna = 0

    for feat in rec["features"]:
        norm = _normalise(feat["gene"])
        if feat["type"] == "CDS":
            if norm in expected:
                pcg_found.add(norm)
        elif feat["type"] == "rRNA":
            if norm in _RRNA_CANONICAL:
                rrna_found.add(norm)
        elif feat["type"] == "tRNA":
            n_trna += 1

    missing_pcgs = sorted(expected - pcg_found)
    extra_pcgs   = sorted(pcg_found - expected)

    flags: list[str] = []
    if missing_pcgs:
        flags.append(f"MISSING_PCG({','.join(missing_pcgs)})")
    if extra_pcgs:
        flags.append(f"EXTRA_PCG({','.join(extra_pcgs)})")
    if len(rrna_found) < 2:
        missing_rrna = sorted(_RRNA_CANONICAL - rrna_found)
        flags.append(f"MISSING_RRNA({','.join(missing_rrna)})")
    if n_trna < min_trna:
        flags.append(f"LOW_TRNA({n_trna}<{min_trna})")

    return {
        "m2_n_pcg":       len(pcg_found),
        "m2_n_rrna":      len(rrna_found),
        "m2_n_trna":      n_trna,
        "m2_missing_pcgs": ",".join(missing_pcgs) if missing_pcgs else "",
        "m2_extra_pcgs":   ",".join(extra_pcgs)   if extra_pcgs   else "",
        "m2_flags":        "|".join(flags) if flags else "PASS",
        "m2_pass":         len(flags) == 0,
    }
