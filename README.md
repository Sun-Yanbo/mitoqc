# MitoQC v0.2

**Systematic quality control for metazoan mitochondrial genomes**

MitoQC is a Python package that assigns quality grades to complete metazoan mitogenomes from NCBI GenBank/RefSeq. It runs five independent QC modules covering sequence integrity, annotation completeness, chimera detection, COI extractability, and protein-coding gene integrity.

## Installation

```bash
pip install -e .
# or
pip install biopython>=1.79 pandas>=1.5 numpy>=1.23
```

**Requirements**: Python ≥ 3.9, Biopython ≥ 1.79, pandas ≥ 1.5

## Quick Start

```bash
# Run on a local GenBank file
python -m mitoqc --input sequences.gb --output results.csv --verbose

# Run on NCBI accessions
python -m mitoqc --accessions NC_012920 NC_001643 --output results.csv

# Use v2 modules (ORF scan + BLAST interface)
python -m mitoqc --input sequences.gb --output results.csv --v2-modules
```

## Python API

```python
from mitoqc.utils.io import iter_genbank_file
from mitoqc.pipeline import run_pipeline
import pandas as pd

records = list(iter_genbank_file("sequences.gb"))
results = [run_pipeline(rec) for rec in records]
df = pd.DataFrame(results)
print(df["quality_grade"].value_counts())
```

## Quality Grades

| Grade  | Criteria                          | Recommended use                        |
|--------|-----------------------------------|----------------------------------------|
| Gold   | All 5 modules pass                | Phylogenomics, ML training, population genomics |
| Silver | Exactly 1 non-critical module fails | Use with caution; review specific flag |
| Bronze | 1 critical module fails OR 2 modules fail | Manual curation recommended    |
| Fail   | 3+ modules fail                   | Not recommended without curation       |

**Critical modules**: M1 (sequence integrity) and M5 (PCG integrity).

## QC Modules

| Module | Name                  | What it checks                                      |
|--------|-----------------------|-----------------------------------------------------|
| M1     | Sequence Integrity    | Length outliers, N-content >1%, GC outliers         |
| M2     | Annotation Completeness | Missing PCGs, rRNAs, tRNAs in feature table       |
| M3     | Chimera Detection     | GC-content heterogeneity across 500 bp windows      |
| M4     | COI/Species Validation | COI extractability, binomial name, ambiguous taxa  |
| M5     | PCG Integrity         | Internal stop codons in translated CDS features     |

## v0.2 Upgrades

- **M2 v2**: ORF scan cross-check validates annotated CDS boundaries against *de novo* ORFs
- **M4 v2**: Optional NCBI BLAST species verification (disabled by default)
- **M5**: Taxon-aware exception for avian *nad3* RNA-editing frameshift (Mindell et al. 1998)
- **Empirical baselines**: `constants.py` updated with rates from n=9,931 records

## Empirical Baselines (n=9,931)

| Grade  | Rate  | | Module | Fail rate |
|--------|-------|-|--------|-----------|
| Gold   | 75.9% | | M1     | 9.1%      |
| Silver | 15.1% | | M2     | 10.7%     |
| Bronze | 8.3%  | | M3     | 0.4%      |
| Fail   | 0.8%  | | M4     | 6.3%      |
|        |       | | M5     | 4.6%      |

## Output Columns

The output CSV contains one row per record with columns:
- `accession`, `organism`, `taxon_group`, `seq_len`, `gc_pct`, `date`, `source`
- Per-module flags and pass/fail: `m1_*`, `m2_*`, `m3_*`, `m4_*`, `m5_*`
- `quality_grade`, `n_modules_failed`, `failed_modules`, `quality_score`

## Running Tests

```bash
pip install pytest
pytest tests/ -v
# Expected: 46+ tests passing
```

## Citation

If you use MitoQC in your research, please cite:

> Liu J, Sun Y.B. (2026) MitoQC: A systematic quality control framework for
> metazoan mitochondrial genomes reveals taxon-specific error patterns across
> NCBI GenBank and RefSeq. *Unpublished*.

## License

MIT License. See LICENSE file.
