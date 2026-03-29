"""
scorer.py – Quality grade assignment for MitoQC.

Aggregates the pass/fail results from all five QC modules (M1–M5) into
a single quality grade and numeric score.

GRADING SCHEME
--------------
The scheme is designed to be informative for downstream use:

  Gold   – passes all 5 modules.  Suitable for high-confidence analyses
            (phylogenomics, ML training sets, population genomics).

  Silver – fails exactly 1 non-critical module.  Usable with caution;
            the specific flag should be reviewed.

  Bronze – fails exactly 1 critical module, OR fails exactly 2 modules
            (any combination).  Requires careful evaluation before use.

  Fail   – fails 3 or more modules.  Not recommended for downstream
            analyses without manual curation.

CRITICAL MODULES
----------------
M1 (sequence integrity) and M5 (PCG integrity) are designated critical
because their failures indicate fundamental sequence-level problems that
cannot be resolved by re-annotation.  A single critical failure → Bronze.

QUALITY SCORE
-------------
A continuous score in [0, 1] is also computed:
  score = (sum of module pass values) / 5
This allows ranking within grades and is useful for ML training set
selection (e.g. select top-N records by score within Gold grade).
"""

from __future__ import annotations

_CRITICAL_MODULES = frozenset({"m1_pass", "m5_pass"})
_ALL_MODULES      = ("m1_pass", "m2_pass", "m3_pass", "m4_pass", "m5_pass")


def grade(result: dict, config: dict | None = None) -> dict:
    """Assign a quality grade to a record based on module pass/fail results.

    Parameters
    ----------
    result:
        Dict containing boolean values for ``m1_pass`` through ``m5_pass``.
        Missing keys are treated as passed (not penalised).
    config:
        Optional overrides:
        - ``critical_modules`` (set[str]): override the critical module set.

    Returns
    -------
    dict
        Keys: ``quality_grade``, ``n_modules_failed``, ``failed_modules``,
        ``quality_score``.
    """
    cfg      = config or {}
    critical = cfg.get("critical_modules", _CRITICAL_MODULES)

    passed  = [result.get(m, True) for m in _ALL_MODULES]
    n_fail  = sum(1 for p in passed if not p)
    n_pass  = len(_ALL_MODULES) - n_fail

    failed_names = [
        m.replace("_pass", "").upper()
        for m, p in zip(_ALL_MODULES, passed) if not p
    ]

    # Determine grade
    if n_fail == 0:
        quality_grade = "Gold"
    elif n_fail == 1:
        failed_mod = [m for m, p in zip(_ALL_MODULES, passed) if not p][0]
        if failed_mod in critical:
            quality_grade = "Bronze"
        else:
            quality_grade = "Silver"
    elif n_fail == 2:
        quality_grade = "Bronze"
    else:
        quality_grade = "Fail"

    quality_score = round(n_pass / len(_ALL_MODULES), 4)

    return {
        "quality_grade":    quality_grade,
        "n_modules_failed": n_fail,
        "failed_modules":   ",".join(failed_names) if failed_names else None,
        "quality_score":    quality_score,
    }
