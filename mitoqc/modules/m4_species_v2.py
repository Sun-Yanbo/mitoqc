"""
M4 v2 – COI Extractability, Species Name Validation, and BLAST Verification

Extends M4 v1 with an optional NCBI BLAST verification layer.

LAYER 1 (inherited from M4 v1)
-------------------------------
- COI extractability: extract cox1 CDS from feature table
- COI length capture: record actual extracted length
- Species name validation: check for binomial name, ambiguous placeholders

LAYER 2 (new in v2, disabled by default)
-----------------------------------------
Submit the extracted COI sequence to NCBI BLAST (blastn, nt database) and
verify:
  1. Top hit identity ≥ 97% (Hebert et al. 2003 species threshold)
  2. Top hit organism genus matches the deposited organism genus

Enabled via ``config={"run_blast": True}``.  Requires network access.
With an NCBI API key (~10 req/s), processing 5,000 records takes ~8 min.
Without a key (~3 req/s), ~28 min.

RATE LIMITS
-----------
NCBI allows 3 requests/second without an API key, 10/second with one.
Set ``config={"api_key": "YOUR_KEY"}`` to use an API key.
Estimated runtime for 5,000 records: 28–83 min depending on key usage.

REFERENCES
----------
Hebert et al. 2003 Proc R Soc B 270:313 — COI barcoding threshold (97%)
"""

from __future__ import annotations

import re
import time

from mitoqc.utils.constants import COI_CANONICAL_NAMES, GENE_NORM

_COI_MIN_LEN = 600

_AMBIG_PATTERNS = re.compile(
    r"\bsp\.\b|\bcf\.\b|\baff\.\b|\bnr\.\b|"
    r"\bgen\.\s*nov\.\b|\bsp\.\s*nov\.\b",
    re.IGNORECASE,
)
_NUMERIC_STRAIN = re.compile(r"\b\d{4,}\b")


def _normalise(name: str) -> str:
    return GENE_NORM.get(name.lower().strip(), name.lower().strip())


def _is_coi(gene_name: str) -> bool:
    return _normalise(gene_name) in COI_CANONICAL_NAMES


def _extract_coi(rec: dict) -> tuple[str | None, int, str]:
    """Return (coi_sequence, length, gene_name) or (None, 0, '') if not found."""
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
        return coi_seq, len(coi_seq), feat["gene"]
    return None, 0, ""


def _genus_match(organism1: str, organism2: str) -> bool:
    """Return True if both names share the same genus (first word)."""
    w1 = organism1.strip().split()
    w2 = organism2.strip().split()
    if not w1 or not w2:
        return False
    return w1[0].lower() == w2[0].lower()


def _blast_verify(coi_seq: str, organism: str,
                  api_key: str | None = None,
                  identity_threshold: float = 97.0) -> dict:
    """Submit COI to NCBI BLAST and return verification result."""
    from Bio.Blast import NCBIWWW, NCBIXML

    kwargs = {
        "program":   "blastn",
        "database":  "nt",
        "sequence":  coi_seq,
        "hitlist_size": 5,
        "entrez_query": "Metazoa[Organism]",
    }
    if api_key:
        kwargs["api_key"] = api_key

    try:
        handle  = NCBIWWW.qblast(**kwargs)
        records = list(NCBIXML.parse(handle))
        handle.close()
    except Exception as exc:
        return {
            "m4_blast_run":      True,
            "m4_blast_identity": None,
            "m4_blast_hit":      None,
            "m4_blast_genus_ok": None,
            "m4_blast_error":    str(exc),
        }

    if not records or not records[0].alignments:
        return {
            "m4_blast_run":      True,
            "m4_blast_identity": None,
            "m4_blast_hit":      "NO_HIT",
            "m4_blast_genus_ok": False,
            "m4_blast_error":    None,
        }

    top_aln = records[0].alignments[0]
    top_hsp = top_aln.hsps[0]
    identity_pct = top_hsp.identities / top_hsp.align_length * 100
    hit_title    = top_aln.title
    genus_ok     = _genus_match(organism, hit_title)

    return {
        "m4_blast_run":      True,
        "m4_blast_identity": round(identity_pct, 2),
        "m4_blast_hit":      hit_title[:80],
        "m4_blast_genus_ok": genus_ok,
        "m4_blast_error":    None,
    }


def run(rec: dict, config: dict | None = None) -> dict:
    """Run M4 v2 COI + species name + optional BLAST checks.

    Parameters
    ----------
    rec:
        Parsed record dict.
    config:
        Optional overrides:
        - ``coi_min_len`` (int): minimum COI length (default 600).
        - ``run_blast`` (bool): enable BLAST layer (default False).
        - ``api_key`` (str): NCBI API key for higher rate limits.
        - ``identity_threshold`` (float): BLAST identity threshold (default 97.0).

    Returns
    -------
    dict
        All M4 v1 keys plus BLAST fields if run_blast=True.
    """
    cfg                = config or {}
    coi_min_len        = cfg.get("coi_min_len", _COI_MIN_LEN)
    run_blast          = cfg.get("run_blast", False)
    api_key            = cfg.get("api_key", None)
    identity_threshold = cfg.get("identity_threshold", 97.0)
    organism           = rec.get("organism", "")

    coi_seq, coi_len, coi_name = _extract_coi(rec)
    coi_found = coi_seq is not None and coi_len >= coi_min_len

    flags: list[str] = []

    # COI check
    if coi_seq is None:
        flags.append("COI_NOT_FOUND")
    elif coi_len < coi_min_len:
        flags.append(f"COI_TOO_SHORT({coi_len}<{coi_min_len})")

    # Species name check
    words = organism.strip().split()
    if len(words) < 2:
        flags.append(f"NON_BINOMIAL_NAME('{organism}')")
    elif _AMBIG_PATTERNS.search(organism):
        flags.append("AMBIGUOUS_TAXON")
    elif _NUMERIC_STRAIN.search(organism):
        flags.append("NUMERIC_STRAIN_ID")

    result = {
        "m4_coi_found":    coi_seq is not None,
        "m4_coi_name":     coi_name,
        "m4_coi_len":      coi_len,
        "m4_organism_ok":  len(flags) == 0,
        "m4_blast_run":    False,
        "m4_blast_identity": None,
        "m4_blast_hit":    None,
        "m4_blast_genus_ok": None,
        "m4_blast_error":  None,
        "m4_flags":        "|".join(flags) if flags else "PASS",
        "m4_pass":         len(flags) == 0,
    }

    # Optional BLAST layer
    if run_blast and coi_found:
        blast_result = _blast_verify(
            coi_seq, organism,
            api_key=api_key,
            identity_threshold=identity_threshold,
        )
        result.update(blast_result)
        # Add BLAST-based flags
        if blast_result.get("m4_blast_identity") is not None:
            if blast_result["m4_blast_identity"] < identity_threshold:
                result["m4_flags"] += f"|LOW_BLAST_IDENTITY({blast_result['m4_blast_identity']:.1f}%)"
                result["m4_pass"] = False
            if blast_result.get("m4_blast_genus_ok") is False:
                result["m4_flags"] += "|BLAST_GENUS_MISMATCH"
                result["m4_pass"] = False

    return result
