"""
pipeline.py – MitoQC main pipeline orchestrator.

Runs all five QC modules (M1–M5) on a single record dict and returns
a flat result dict suitable for appending to a pandas DataFrame.

Usage
-----
>>> from mitoqc.pipeline import run_pipeline
>>> result = run_pipeline(rec)
>>> print(result["quality_grade"])
'Gold'

For batch processing, use the CLI (``python -m mitoqc``) or call
``run_pipeline`` in a loop over records from
:func:`mitoqc.utils.io.iter_genbank_file` or
:func:`mitoqc.utils.io.fetch_genbank_records`.
"""

from __future__ import annotations

from mitoqc.modules import (
    m1_integrity,
    m2_annotation,
    m3_chimera,
    m4_species,
    m5_pcg,
)
from mitoqc.modules.scorer import grade

# Module registry — order matters for result dict column ordering
_MODULES = [
    ("m1", m1_integrity),
    ("m2", m2_annotation),
    ("m3", m3_chimera),
    ("m4", m4_species),
    ("m5", m5_pcg),
]


def run_pipeline(rec: dict, config: dict | None = None) -> dict:
    """Run all MitoQC modules on a single record.

    Parameters
    ----------
    rec:
        Standardised record dict from :func:`mitoqc.utils.io.parse_genbank`.
    config:
        Optional per-module configuration overrides.  Keys are module
        prefixes (``"m1"``, ``"m2"``, etc.); values are dicts passed to
        the corresponding module's ``run()`` function.

        Example::

            config = {
                "m1": {"n_threshold": 0.02},
                "m4": {"run_blast": True, "api_key": "YOURKEY"},
            }

    Returns
    -------
    dict
        Flat result dict with all module outputs plus quality grade fields.
        Suitable for direct use as a pandas DataFrame row.
    """
    cfg = config or {}

    # Metadata passthrough
    result: dict = {
        "accession":   rec.get("accession", ""),
        "organism":    rec.get("organism", ""),
        "taxon_group": rec.get("taxon_group", "Other_Metazoa"),
        "seq_len":     rec.get("seq_len", 0),
        "gc_pct":      round(
            (rec["seq"].count("G") + rec["seq"].count("C"))
            / max(rec.get("seq_len", 1), 1) * 100, 2
        ) if rec.get("seq") else None,
        "date":        rec.get("date", ""),
        "source":      rec.get("source", ""),
    }

    # Run each module
    for mod_name, mod in _MODULES:
        mod_cfg = cfg.get(mod_name)
        try:
            result.update(mod.run(rec, config=mod_cfg))
        except Exception as exc:
            result[f"{mod_name}_pass"]  = None
            result[f"{mod_name}_flags"] = f"ERROR:{exc}"

    # Assign quality grade
    result.update(grade(result))
    return result
