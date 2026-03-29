"""
M5 – Protein-Coding Gene (PCG) Integrity

Translates each annotated CDS feature and checks for premature internal
stop codons, which are the most reliable sequence-level indicator of:
  - Assembly errors (frameshifts, mis-joins)
  - Annotation boundary errors (wrong start/stop position)
  - NUMTs (nuclear mitochondrial DNA segments) mis-annotated as mtDNA

GENETIC CODE HANDLING
---------------------
Metazoan mitogenomes use non-standard genetic codes.  The correct table
is selected per taxon group:

  Table  2  – Vertebrate Mitochondrial   (Mammalia, Aves, Reptilia,
               Amphibia, Actinopterygii)
  Table  5  – Invertebrate Mitochondrial  (Insecta, Arachnida, Crustacea,
               Mollusca, Nematoda, Echinodermata, Other_Metazoa)

INCOMPLETE STOP CODONS
-----------------------
Many metazoan PCGs end with a single T or TA that is completed to TAA by
post-transcriptional polyadenylation.  Biopython's ``translate(to_stop=False)``
will render these as '*' at the final position.  We therefore only flag
*internal* stop codons (position < last codon), not terminal ones.

AVES nad3 FRAMESHIFT (TAXON-AWARE EXCEPTION)
---------------------------------------------
Avian nad3 contains a well-documented +1 frameshift at a conserved position
that is corrected by post-transcriptional RNA editing (Mindell et al. 1998
Mol Biol Evol 15:1568).  GenBank annotations of avian nad3 therefore
routinely include an internal stop codon when translated naively.
MitoQC v0.2 empirically confirmed this: 47% of RefSeq Aves and 10.6% of
GenBank Aves records fail M5 exclusively on nad3.

To avoid systematic false positives, nad3 is excluded from M5 stop-codon
checking when taxon_group == "Aves".  This is stored in
``M5_TAXON_GENE_EXCEPTIONS`` in constants.py and applied at runtime.

All other Aves PCGs are still checked normally.
"""

from __future__ import annotations

from Bio.Seq import Seq

from mitoqc.utils.constants import (
    GENETIC_CODE_TABLE, GENE_NORM, M5_TAXON_GENE_EXCEPTIONS
)


def _normalise(name: str) -> str:
    return GENE_NORM.get(name.lower().strip(), name.lower().strip())


def _translate_cds(seq_str: str, table: int) -> tuple[str, bool]:
    """Translate a CDS string; return (protein, has_internal_stop).

    Handles sequences whose length is not a multiple of 3 by trimming the
    trailing 1–2 bases (incomplete terminal codon).
    """
    trim = len(seq_str) - (len(seq_str) % 3)
    cds  = Seq(seq_str[:trim]) if trim > 0 else Seq(seq_str)

    try:
        protein = str(cds.translate(table=table, to_stop=False))
    except Exception:
        return "", False

    # Internal stop = '*' anywhere except the last position
    internal = "*" in protein[:-1]
    return protein, internal


def run(rec: dict, config: dict | None = None) -> dict:
    """Run M5 PCG integrity checks on a single record.

    Parameters
    ----------
    rec:
        Parsed record dict.  Must include ``seq`` and ``features``.
    config:
        Optional overrides:
        - ``genetic_code`` (int): force a specific NCBI translation table.

    Returns
    -------
    dict
        Keys: ``m5_n_cds_checked``, ``m5_n_internal_stop``,
        ``m5_affected_genes``, ``m5_flags``, ``m5_pass``.
    """
    cfg   = config or {}
    seq   = rec["seq"]
    group = rec["taxon_group"]
    table = cfg.get("genetic_code") or GENETIC_CODE_TABLE.get(group, 5)

    # Genes to skip for this taxon (e.g. avian nad3 RNA-editing frameshift)
    skip_genes: frozenset[str] = M5_TAXON_GENE_EXCEPTIONS.get(group, frozenset())

    n_checked      = 0
    n_internal     = 0
    affected_genes: list[str] = []

    for feat in rec["features"]:
        if feat["type"] != "CDS":
            continue

        start  = feat["start"]
        end    = feat["end"]
        strand = feat["strand"]
        gene   = _normalise(feat["gene"])

        # Skip taxon-specific known RNA-editing frameshifts
        if gene in skip_genes:
            continue

        if strand == -1:
            cds_seq = str(Seq(seq[start:end]).reverse_complement())
        else:
            cds_seq = seq[start:end]

        if len(cds_seq) < 3:
            continue

        n_checked += 1
        _, has_internal = _translate_cds(cds_seq, table)
        if has_internal:
            n_internal += 1
            affected_genes.append(gene)

    flags: list[str] = []
    if n_internal > 0:
        flags.append(
            f"INTERNAL_STOP({n_internal}/{n_checked},"
            f"genes={','.join(affected_genes)})"
        )

    return {
        "m5_n_cds_checked":   n_checked,
        "m5_n_internal_stop": n_internal,
        "m5_affected_genes":  ",".join(affected_genes) if affected_genes else "",
        "m5_flags":           "|".join(flags) if flags else "PASS",
        "m5_pass":            len(flags) == 0,
    }
