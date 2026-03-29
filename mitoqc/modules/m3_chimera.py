"""
M3 – Chimera Detection Proxy

Detects potential chimeric assemblies by measuring GC-content heterogeneity
across non-overlapping windows of the mitogenome.  A chimeric sequence
(formed by joining segments from two different species) typically shows
an abrupt GC shift at the junction, producing an elevated standard
deviation of per-window GC content.

ALGORITHM
---------
1. Divide the sequence into non-overlapping windows of ``window_size`` bp
   (default 500 bp).  Discard the trailing partial window.
2. Compute GC% for each window.
3. Compute the standard deviation (SD) of all window GC values.
4. Flag the record if SD exceeds the taxon-specific threshold stored in
   ``CHIMERA_GC_SD_THRESHOLDS``.

VALIDATION (MitoQC v0.2)
------------------------
Validated on a synthetic chimera set (n=248: 39 chimeras + 209 negatives)
constructed by splicing real RefSeq sequences at 25/50/75% breakpoints
across 8 taxon-pair combinations.

  ROC-AUC = 0.984   Average Precision = 0.921
  Optimal threshold (best F1 = 0.875):
    GC-SD = 13.6%  →  Sensitivity=0.897, Specificity=0.971,
                       PPV=0.854, NPV=0.981

KNOWN LIMITATION
----------------
Within-group splices (e.g. Insecta × Insecta) produce GC-SD 16–20% due
to natural D-loop vs. coding-region GC heterogeneity, causing false
positives at low thresholds.  Taxon-specific thresholds are calibrated
to minimise this effect.

THRESHOLDS
----------
Stored in ``mitoqc.utils.constants.CHIMERA_GC_SD_THRESHOLDS``.
Calibrated as mean + 3 SD of per-record GC-window SD from RefSeq
Release 234 pilot (n=4,985).
"""

from __future__ import annotations

import statistics

from mitoqc.utils.constants import CHIMERA_GC_SD_THRESHOLDS

_DEFAULT_WINDOW = 500
_MIN_WINDOWS    = 5   # require at least 5 windows for a meaningful SD


def _gc_pct(window: str) -> float:
    n = len(window)
    if n == 0:
        return 0.0
    return (window.count("G") + window.count("C")) / n * 100


def run(rec: dict, config: dict | None = None) -> dict:
    """Run M3 chimera detection on a single record.

    Parameters
    ----------
    rec:
        Parsed record dict.  Must include ``seq`` and ``taxon_group``.
    config:
        Optional overrides:
        - ``window_size`` (int): window size in bp (default 500).
        - ``gc_sd_threshold`` (float): override taxon-specific threshold.

    Returns
    -------
    dict
        Keys: ``m3_gc_sd``, ``m3_gc_threshold``, ``m3_n_windows``,
        ``m3_flags``, ``m3_pass``.
    """
    cfg         = config or {}
    seq         = rec["seq"]
    group       = rec["taxon_group"]
    window_size = cfg.get("window_size", _DEFAULT_WINDOW)
    threshold   = cfg.get("gc_sd_threshold",
                          CHIMERA_GC_SD_THRESHOLDS.get(group, 11.0))

    # Compute per-window GC%
    windows = [
        seq[i:i + window_size]
        for i in range(0, len(seq) - window_size + 1, window_size)
    ]
    n_windows = len(windows)

    if n_windows < _MIN_WINDOWS:
        return {
            "m3_gc_sd":        None,
            "m3_gc_threshold": threshold,
            "m3_n_windows":    n_windows,
            "m3_flags":        "INSUFFICIENT_WINDOWS",
            "m3_pass":         True,   # cannot assess; do not penalise
        }

    gc_values = [_gc_pct(w) for w in windows]
    gc_sd     = statistics.stdev(gc_values)

    flags: list[str] = []
    if gc_sd > threshold:
        flags.append(f"HIGH_GC_SD({gc_sd:.2f}>{threshold})")

    return {
        "m3_gc_sd":        round(gc_sd, 4),
        "m3_gc_threshold": threshold,
        "m3_n_windows":    n_windows,
        "m3_flags":        "|".join(flags) if flags else "PASS",
        "m3_pass":         len(flags) == 0,
    }
