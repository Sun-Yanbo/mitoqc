"""
M2 v2 – Annotation Completeness with ORF Scan Cross-Check

Extends M2 v1 with an independent ORF scan layer that validates annotated
CDS boundaries against *de novo* open reading frames found in the sequence.

TWO-LAYER APPROACH
------------------
Layer 1 (inherited from M2 v1):
  Check presence of expected PCGs, rRNAs, and tRNAs in the feature table.

Layer 2 (new in v2):
  For each annotated CDS, scan all 6 reading frames in the surrounding
  region for ORFs ≥ ``min_orf_aa`` amino acids (default 100 aa = 300 bp).
  If no valid ORF overlaps the annotated CDS coordinates, flag as
  CDS_NO_ORF — indicating a likely boundary error or NUMT.

CONFIGURATION
-------------
ORF scan is enabled by default for sequences ≥ 5,000 bp.  It can be
disabled via ``config={"run_orf_scan": False}`` for speed.

MITOS2 INTEGRATION (PLANNED v0.3)
----------------------------------
Full re-annotation via MITOS2 (Bernt et al. 2013 Mol Phylogenet Evol
69:313) is planned for v0.3.  MITOS2 provides independent gene boundary
prediction and would replace the ORF scan heuristic.  Integration requires
either a local Docker installation or an accessible REST endpoint.
"""

from __future__ import annotations

from Bio.Seq import Seq

from mitoqc.utils.constants import (
    EXPECTED_PCGS, GROUPS_LACKING_ATP8, MIN_TRNA_COUNT, GENE_NORM,
    GENETIC_CODE_TABLE,
)

_RRNA_CANONICAL = frozenset({"rrnS", "rrnL"})
_MIN_ORF_AA     = 100   # minimum ORF length in amino acids
_MIN_SEQ_LEN    = 5_000  # skip ORF scan for very short sequences


def _normalise(name: str) -> str:
    return GENE_NORM.get(name.lower().strip(), name.lower().strip())


def _find_orfs(seq: str, table: int, min_aa: int = _MIN_ORF_AA) -> list[tuple[int, int]]:
    """Return list of (start, end) for all ORFs ≥ min_aa in all 6 frames."""
    orfs: list[tuple[int, int]] = []
    seq_len = len(seq)
    rc_seq  = str(Seq(seq).reverse_complement())

    for strand_seq, is_rc in [(seq, False), (rc_seq, True)]:
        for frame in range(3):
            i = frame
            while i + 3 <= len(strand_seq):
                codon = strand_seq[i:i+3]
                # Find start codon
                if codon == "ATG":
                    j = i + 3
                    while j + 3 <= len(strand_seq):
                        stop = strand_seq[j:j+3]
                        if stop in ("TAA", "TAG", "TGA"):
                            orf_len_aa = (j - i) // 3
                            if orf_len_aa >= min_aa:
                                if is_rc:
                                    # Convert RC coordinates back to forward
                                    fwd_end   = seq_len - i
                                    fwd_start = seq_len - j - 3
                                    orfs.append((fwd_start, fwd_end))
                                else:
                                    orfs.append((i, j + 3))
                            break
                        j += 3
                i += 3
    return orfs


def _overlaps(feat_start: int, feat_end: int,
              orfs: list[tuple[int, int]]) -> bool:
    """Return True if any ORF overlaps the feature by ≥ 50%."""
    feat_len = feat_end - feat_start
    if feat_len <= 0:
        return False
    for orf_start, orf_end in orfs:
        overlap = min(feat_end, orf_end) - max(feat_start, orf_start)
        if overlap / feat_len >= 0.5:
            return True
    return False


def run(rec: dict, config: dict | None = None) -> dict:
    """Run M2 v2 annotation completeness + ORF scan checks.

    Parameters
    ----------
    rec:
        Parsed record dict.
    config:
        Optional overrides:
        - ``min_trna`` (int): minimum tRNA count (default 18).
        - ``run_orf_scan`` (bool): enable ORF scan layer (default True).
        - ``min_orf_aa`` (int): minimum ORF length in aa (default 100).

    Returns
    -------
    dict
        All M2 v1 keys plus: ``m2_orf_scan_run``, ``m2_cds_no_orf``.
    """
    cfg          = config or {}
    group        = rec["taxon_group"]
    min_trna     = cfg.get("min_trna", MIN_TRNA_COUNT)
    run_orf_scan = cfg.get("run_orf_scan", True)
    min_orf_aa   = cfg.get("min_orf_aa", _MIN_ORF_AA)
    seq          = rec["seq"]
    table        = GENETIC_CODE_TABLE.get(group, 5)

    # ── Layer 1: annotation proxy (same as M2 v1) ────────────────────────────
    expected = set(EXPECTED_PCGS)
    if group in GROUPS_LACKING_ATP8:
        expected.discard("atp8")

    pcg_found:  set[str] = set()
    rrna_found: set[str] = set()
    n_trna = 0
    cds_features: list[dict] = []

    for feat in rec["features"]:
        norm = _normalise(feat["gene"])
        if feat["type"] == "CDS":
            if norm in expected:
                pcg_found.add(norm)
            cds_features.append(feat)
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

    # ── Layer 2: ORF scan cross-check ────────────────────────────────────────
    orf_scan_run = False
    cds_no_orf   = 0

    if run_orf_scan and len(seq) >= _MIN_SEQ_LEN:
        orf_scan_run = True
        orfs = _find_orfs(seq, table, min_aa=min_orf_aa)
        for feat in cds_features:
            if not _overlaps(feat["start"], feat["end"], orfs):
                cds_no_orf += 1
        if cds_no_orf > 0:
            flags.append(f"CDS_NO_ORF({cds_no_orf})")

    return {
        "m2_n_pcg":        len(pcg_found),
        "m2_n_rrna":       len(rrna_found),
        "m2_n_trna":       n_trna,
        "m2_missing_pcgs": ",".join(missing_pcgs) if missing_pcgs else "",
        "m2_extra_pcgs":   ",".join(extra_pcgs)   if extra_pcgs   else "",
        "m2_orf_scan_run": orf_scan_run,
        "m2_cds_no_orf":   cds_no_orf,
        "m2_flags":        "|".join(flags) if flags else "PASS",
        "m2_pass":         len(flags) == 0,
    }
