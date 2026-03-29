"""
constants.py – Biological constants used across all MitoQC modules.

All numeric thresholds are documented with their source.  Where values
are derived empirically from the RefSeq Release 234 pilot dataset, this
is stated explicitly so future users can recalibrate on larger corpora.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Expected mitogenome length ranges (bp) per taxon group
# Source: empirical mean ± 3 SD from NCBI RefSeq Release 234 (n=16,075)
#         cross-checked against:
#         - Boore 1999 Nucleic Acids Res 27:1767 (metazoan mtDNA overview)
#         - Cameron 2014 Curr Biol 24:R488 (insect mtDNA)
# ---------------------------------------------------------------------------
EXPECTED_LEN: dict[str, tuple[int, int]] = {
    "Mammalia":        (15_000, 18_000),
    "Aves":            (14_500, 18_500),
    "Reptilia":        (14_000, 20_000),
    "Amphibia":        (14_000, 20_000),
    "Actinopterygii":  (14_000, 20_000),
    "Insecta":         (13_000, 22_000),
    "Arachnida":       (12_000, 22_000),
    "Crustacea":       (13_000, 22_000),
    "Mollusca":        (13_000, 25_000),
    "Nematoda":        ( 9_000, 22_000),   # wide range; nematodes are variable
    "Echinodermata":   (14_000, 20_000),
    "Other_Metazoa":   (10_000, 30_000),
}

# ---------------------------------------------------------------------------
# Expected protein-coding genes (PCGs) – canonical 13-gene set
# Source: Anderson et al. 1981 Nature 290:457 (human mtDNA reference)
#         Boore 1999 Nucleic Acids Res 27:1767
# ---------------------------------------------------------------------------
EXPECTED_PCGS: frozenset[str] = frozenset({
    "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6",
    "cox1", "cox2", "cox3",
    "atp6", "atp8",
    "cytb",
})

# Groups where atp8 is known to be absent (biological, not annotation error)
# Source: Hu et al. 2003 Int J Parasitol 33:1171 (Strongyloides atp8 absence)
#         Lavrov & Brown 2001 Genetics 157:621 (nematode mtDNA)
GROUPS_LACKING_ATP8: frozenset[str] = frozenset({
    "Nematoda",
})

# Minimum tRNA count expected in a complete metazoan mitogenome
# Source: Wolstenholme 1992 Int Rev Cytol 141:173 (22 tRNAs standard)
#         Reduced to 18 as conservative lower bound for annotation proxy
MIN_TRNA_COUNT: int = 18

# ---------------------------------------------------------------------------
# NCBI genetic code tables per taxon group
# Source: NCBI Taxonomy genetic codes
#         https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
# ---------------------------------------------------------------------------
GENETIC_CODE_TABLE: dict[str, int] = {
    # Vertebrate Mitochondrial (Table 2)
    "Mammalia":       2,
    "Aves":           2,
    "Reptilia":       2,
    "Amphibia":       2,
    "Actinopterygii": 2,
    # Invertebrate Mitochondrial (Table 5)
    "Insecta":        5,
    "Arachnida":      5,
    "Crustacea":      5,
    "Mollusca":       5,
    "Nematoda":       5,
    "Echinodermata":  5,
    "Other_Metazoa":  5,
}

# ---------------------------------------------------------------------------
# Chimera detection: GC-SD thresholds per taxon group
# Source: empirically calibrated as mean + 3 SD of per-record GC-window SD
#         from NCBI RefSeq Release 234 pilot dataset (n=4,985)
#         Window size: 500 bp non-overlapping
# ---------------------------------------------------------------------------
CHIMERA_GC_SD_THRESHOLDS: dict[str, float] = {
    "Mammalia":        7.5,
    "Aves":            8.0,
    "Reptilia":        9.0,
    "Amphibia":        9.5,
    "Actinopterygii":  9.5,
    "Insecta":        14.0,
    "Arachnida":      13.0,
    "Crustacea":      17.5,
    "Mollusca":       11.0,
    "Nematoda":       12.0,
    "Echinodermata":  10.0,
    "Other_Metazoa":  11.0,
}

# ---------------------------------------------------------------------------
# M5 taxon-aware gene exceptions: genes to skip stop-codon checking per group
# Source: Mindell et al. 1998 Mol Biol Evol 15:1568 — avian nad3 +1 frameshift
#         corrected by post-transcriptional RNA editing; confirmed empirically
#         in MitoQC v0.2 (47% RefSeq Aves M5 failures exclusively on nad3,
#         n=4,985 + n=4,946 cross-dataset validation, 2025-03)
# ---------------------------------------------------------------------------
M5_TAXON_GENE_EXCEPTIONS: dict[str, frozenset[str]] = {
    "Aves": frozenset({"nad3"}),
}

# ---------------------------------------------------------------------------
# Empirical quality-grade and module-fail rates
# Source: MitoQC v0.2 analysis of n=9,931 metazoan mitogenomes
#         (RefSeq Release 234 pilot n=4,985 + GenBank non-RefSeq n=4,946)
#         Stratified sample across 12 taxon groups, processed 2025-03
# ---------------------------------------------------------------------------
EMPIRICAL_GRADE_RATES: dict[str, float] = {
    "Gold":   0.759,   # 75.9% — passes all 5 modules
    "Silver": 0.151,   # 15.1% — fails exactly 1 module
    "Bronze": 0.083,   #  8.3% — fails exactly 2 modules
    "Fail":   0.008,   #  0.8% — fails 3+ modules
}

# Per-module empirical fail rates (combined RefSeq + GenBank)
EMPIRICAL_MODULE_FAIL_RATES: dict[str, float] = {
    "M1_seq_integrity":  0.091,   # 9.1%  — length/N-content/GC outliers
    "M2_annotation":     0.107,   # 10.7% — missing PCG/rRNA/tRNA annotations
    "M3_chimera":        0.004,   # 0.4%  — GC-heterogeneity chimera proxy
    "M4_coi_species":    0.063,   # 6.3%  — COI not extractable / ambiguous taxon
    "M5_pcg_integrity":  0.046,   # 4.6%  — internal stop codons in PCGs
}

# Per-database Gold rates (for downstream ML training set selection guidance)
EMPIRICAL_GOLD_RATES_BY_DB: dict[str, float] = {
    "RefSeq":            0.758,
    "GenBank_nonRefSeq": 0.759,
}

# ---------------------------------------------------------------------------
# COI canonical names (post-normalisation)
# Covers all common aliases found in GenBank annotations
# ---------------------------------------------------------------------------
COI_CANONICAL_NAMES: frozenset[str] = frozenset({
    "cox1", "coi", "co1", "coxi",
    "cytochrome oxidase subunit 1",
    "cytochrome oxidase subunit i",
    "cytochrome c oxidase subunit 1",
    "cytochrome c oxidase subunit i",
})

# ---------------------------------------------------------------------------
# Gene name normalisation map
# Maps common aliases / typos → canonical lower-case name
# Sources:
#   - NCBI RefSeq annotation conventions
#   - Bernt et al. 2013 Mol Phylogenet Evol (MITOS gene names)
#   - Iwasaki et al. 2013 Mol Biol Evol (MitoFish conventions)
# ---------------------------------------------------------------------------
GENE_NORM: dict[str, str] = {
    # ND subunits
    "nd1": "nad1", "nd2": "nad2", "nd3": "nad3",
    "nd4": "nad4", "nd4l": "nad4l",
    "nd5": "nad5", "nd6": "nad6",
    "nadh1": "nad1", "nadh2": "nad2", "nadh3": "nad3",
    "nadh4": "nad4", "nadh4l": "nad4l",
    "nadh5": "nad5", "nadh6": "nad6",
    "nadh dehydrogenase subunit 1": "nad1",
    "nadh dehydrogenase subunit 2": "nad2",
    "nadh dehydrogenase subunit 3": "nad3",
    "nadh dehydrogenase subunit 4": "nad4",
    "nadh dehydrogenase subunit 4l": "nad4l",
    "nadh dehydrogenase subunit 5": "nad5",
    "nadh dehydrogenase subunit 6": "nad6",
    # COX subunits
    "cox1": "cox1", "coi": "cox1", "co1": "cox1", "coxi": "cox1",
    "cox2": "cox2", "coii": "cox2", "co2": "cox2",
    "cox3": "cox3", "coiii": "cox3", "co3": "cox3",
    "cytochrome oxidase subunit 1": "cox1",
    "cytochrome oxidase subunit i": "cox1",
    "cytochrome c oxidase subunit 1": "cox1",
    "cytochrome c oxidase subunit i": "cox1",
    "cytochrome oxidase subunit 2": "cox2",
    "cytochrome oxidase subunit ii": "cox2",
    "cytochrome c oxidase subunit 2": "cox2",
    "cytochrome c oxidase subunit ii": "cox2",
    "cytochrome oxidase subunit 3": "cox3",
    "cytochrome oxidase subunit iii": "cox3",
    "cytochrome c oxidase subunit 3": "cox3",
    "cytochrome c oxidase subunit iii": "cox3",
    # ATP synthase
    "atp6": "atp6", "atpase6": "atp6", "atpase 6": "atp6",
    "atp8": "atp8", "atpase8": "atp8", "atpase 8": "atp8",
    "atp synthase subunit 6": "atp6",
    "atp synthase subunit 8": "atp8",
    "atp synthase f0 subunit 6": "atp6",
    "atp synthase f0 subunit 8": "atp8",
    # Cytochrome b
    "cytb": "cytb", "cob": "cytb", "cytochrome b": "cytb",
    "cytochrome b apoenzyme": "cytb",
    # rRNA
    "12s rrna": "rrnS", "12s ribosomal rna": "rrnS",
    "small subunit ribosomal rna": "rrnS",
    "16s rrna": "rrnL", "16s ribosomal rna": "rrnL",
    "large subunit ribosomal rna": "rrnL",
    "s-rrna": "rrnS", "l-rrna": "rrnL",
}
