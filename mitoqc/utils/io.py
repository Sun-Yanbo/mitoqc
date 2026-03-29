"""
io.py – GenBank/RefSeq record parsing utilities for MitoQC.

Provides a single entry point ``parse_genbank`` that converts a
Biopython SeqRecord into the standardised dict format consumed by
all MitoQC QC modules.
"""

from __future__ import annotations

from typing import Iterator
from pathlib import Path
from io import StringIO

from Bio import SeqIO, Entrez
from Bio.Seq import UndefinedSequenceError

from mitoqc.utils.taxonomy import classify_taxon
from mitoqc.utils.constants import GENE_NORM


def _normalise_gene(name: str) -> str:
    return GENE_NORM.get(name.lower().strip(), name.lower().strip())


def _parse_features(record) -> list[dict]:
    """Extract relevant features from a Biopython SeqRecord."""
    feats = []
    for feat in record.features:
        if feat.type not in ("gene", "CDS", "rRNA", "tRNA",
                              "D-loop", "misc_feature"):
            continue
        gene_name = (
            feat.qualifiers.get("gene", [""])[0]
            or feat.qualifiers.get("product", [""])[0]
            or feat.type
        )
        feats.append({
            "type":   feat.type,
            "gene":   gene_name,
            "start":  int(feat.location.start),
            "end":    int(feat.location.end),
            "strand": feat.location.strand,
        })
    return feats


def parse_genbank(record, source: str = "GenBank") -> dict | None:
    """Convert a Biopython SeqRecord to a MitoQC record dict.

    Parameters
    ----------
    record:
        A ``Bio.SeqRecord.SeqRecord`` object parsed from GenBank format.
    source:
        Label for the data source (e.g. ``"RefSeq"``, ``"GenBank_nonRefSeq"``).

    Returns
    -------
    dict or None
        Standardised record dict, or ``None`` if the sequence is undefined.
    """
    try:
        seq_str = str(record.seq).upper()
    except UndefinedSequenceError:
        return None

    taxonomy = record.annotations.get("taxonomy", [])
    lineage  = "; ".join(taxonomy)

    return {
        "accession":   record.id,
        "organism":    record.annotations.get("organism", ""),
        "lineage":     lineage,
        "taxon_group": classify_taxon(lineage),
        "seq_len":     len(record.seq),
        "seq":         seq_str,
        "date":        record.annotations.get("date", ""),
        "platform":    "Unknown",
        "features":    _parse_features(record),
        "n_features":  len(record.features),
        "source":      source,
    }


def iter_genbank_file(path: str | Path, source: str = "GenBank") -> Iterator[dict]:
    """Yield MitoQC record dicts from a local GenBank flat file.

    Parameters
    ----------
    path:
        Path to a ``.gb`` or ``.gbk`` file (may be gzip-compressed).
    source:
        Source label attached to each record.

    Yields
    ------
    dict
        Parsed record dict (skips records with undefined sequences).
    """
    import gzip
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        for rec in SeqIO.parse(fh, "genbank"):
            parsed = parse_genbank(rec, source=source)
            if parsed is not None:
                yield parsed


def fetch_genbank_records(
    accessions: list[str],
    email: str = "mitoqc@example.com",
    batch_size: int = 50,
    sleep: float = 0.35,
) -> Iterator[dict]:
    """Fetch records from NCBI Entrez and yield MitoQC record dicts.

    Parameters
    ----------
    accessions:
        List of accession strings or NCBI UIDs.
    email:
        Email address for NCBI Entrez (required by NCBI policy).
    batch_size:
        Number of records to fetch per API call.
    sleep:
        Seconds to sleep between batches (respect NCBI rate limits).

    Yields
    ------
    dict
        Parsed record dict.
    """
    import time

    Entrez.email = email
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i:i + batch_size]
        try:
            handle  = Entrez.efetch(
                db="nucleotide", id=",".join(batch),
                rettype="gb", retmode="text"
            )
            text = handle.read()
            handle.close()
        except Exception as exc:
            import warnings
            warnings.warn(f"Entrez fetch error for batch {i//batch_size}: {exc}")
            time.sleep(5)
            continue

        for rec in SeqIO.parse(StringIO(text), "genbank"):
            parsed = parse_genbank(rec, source="GenBank")
            if parsed is not None:
                yield parsed

        time.sleep(sleep)
