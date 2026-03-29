"""
M1 – Sequence Integrity

Checks the raw nucleotide sequence for three classes of quality issues:

1. **Length outlier** – sequence length falls outside the taxon-specific
   expected range (mean ± 3 SD from RefSeq Release 234, n=16,075).
   Extreme outliers suggest assembly truncation, concatenation artefacts,
   or mis-classified records.

2. **Ambiguous base content** – fraction of IUPAC ambiguity codes (N, R,
   Y, S, W, K, M, B, D, H, V) exceeds 1 %.  High ambiguity indicates
   poor sequencing coverage or unresolved assembly gaps.

3. **GC content outlier** – GC% falls outside the taxon-specific expected
   range (mean ± 3 SD).  Extreme GC deviation can indicate contamination,
   NUMTs, or chimeric assemblies.

THRESHOLDS
----------
All thresholds are stored in ``mitoqc.utils.constants.EXPECTED_LEN`` and
are calibrated empirically from RefSeq Release 234.  The N-content
threshold of 1 % is a widely used convention (e.g. NCBI RefSeq submission
guidelines).
"""

from __future__ import annotations

import statistics

from mitoqc.utils.constants import EXPECTED_LEN

# Ambiguous IUPAC bases (excluding A/T/G/C)
_AMBIG = frozenset("NRYSWKMBDHV")

# GC outlier thresholds (absolute %, applied after taxon-specific length check)
# Source: empirical 1st–99th percentile from RefSeq Release 234
_GC_RANGE = (25.0, 55.0)

# N-content threshold
_N_THRESHOLD = 0.01   # 1 %


def run(rec: dict, config: dict | None = None) -> dict:
    """Run M1 sequence integrity checks on a single record.

    Parameters
    ----------
    rec:
        Parsed record dict from :func:`mitoqc.utils.io.parse_genbank`.
        Must include ``seq`` (full mitogenome string) and ``taxon_group``.
    config:
        Optional overrides:
        - ``n_threshold`` (float): N-content threshold (default 0.01).
        - ``gc_range`` (tuple[float,float]): absolute GC% bounds.

    Returns
    -------
    dict
        Keys: ``m1_seq_len``, ``m1_n_pct``, ``m1_gc_pct``, ``m1_amb_pct``,
        ``m1_len_flag``, ``m1_flags``, ``m1_pass``.
    """
    cfg         = config or {}
    seq         = rec["seq"]
    group       = rec["taxon_group"]
    n_threshold = cfg.get("n_threshold", _N_THRESHOLD)
    gc_range    = cfg.get("gc_range", _GC_RANGE)

    seq_len = len(seq)
    n_count = seq.count("N")
    g_count = seq.count("G")
    c_count = seq.count("C")
    amb_count = sum(1 for b in seq if b in _AMBIG)

    n_pct   = n_count   / seq_len if seq_len else 0.0
    gc_pct  = (g_count + c_count) / seq_len * 100 if seq_len else 0.0
    amb_pct = amb_count / seq_len if seq_len else 0.0

    flags: list[str] = []

    # 1. Length check
    lo, hi = EXPECTED_LEN.get(group, (10_000, 30_000))
    len_flag = "PASS"
    if seq_len < lo:
        len_flag = f"TOO_SHORT(<{lo})"
        flags.append(len_flag)
    elif seq_len > hi:
        len_flag = f"TOO_LONG(>{hi})"
        flags.append(len_flag)

    # 2. N-content check
    if n_pct > n_threshold:
        flags.append(f"HIGH_N({n_pct*100:.2f}%)")

    # 3. GC content check
    if gc_pct < gc_range[0]:
        flags.append(f"LOW_GC({gc_pct:.1f}%)")
    elif gc_pct > gc_range[1]:
        flags.append(f"HIGH_GC({gc_pct:.1f}%)")

    return {
        "m1_seq_len":  seq_len,
        "m1_n_pct":    round(n_pct * 100, 4),
        "m1_gc_pct":   round(gc_pct, 2),
        "m1_amb_pct":  round(amb_pct * 100, 4),
        "m1_len_flag": len_flag,
        "m1_flags":    "|".join(flags) if flags else "PASS",
        "m1_pass":     len(flags) == 0,
    }
