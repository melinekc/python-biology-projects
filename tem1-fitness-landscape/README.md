# TEM-1 beta-lactamase fitness landscape

Python notebook reproducing and extending the deep mutational scanning
analysis of TEM-1 beta-lactamase from Firnberg et al. (2014, *Molecular
Biology and Evolution*). The dataset covers 18,081 codon-level variants
across 13 ampicillin concentrations (0.25 to 1024 µg/mL), available on
MaveDB as `urn:mavedb:00000070-a-3`.

## What this notebook does

- Loads and validates the raw MaveDB count file
- Computes normalized fitness scores (`w`) relative to wild-type
- Applies a count-based quality filter, justified on this dataset's own
  noise-versus-count trade-off
- Validates internal controls (synonymous, nonsense variants)
- Validates the resulting fitness scores and per-position tolerance
  against the values published in the original paper's supplementary
  data (Spearman rho of 0.965 and 0.932)
- Builds a full position-by-amino-acid fitness heatmap and identifies
  the least and most tolerant positions, cross-referenced with known
  catalytic residues and UniProt secondary structure annotations
- Tests whether intolerant and beneficial positions cluster along the
  sequence, with a permutation test
- Extends the analysis with a codon-level comparison: does the specific
  codon used to encode a substitution affect its fitness, and does this
  depend on the position's tolerance or its location at the start of
  the gene
- Extends this further with an mRNA folding analysis (ViennaRNA) around
  the translation start site

The "Key results" section at the top of the notebook summarises the main
findings. Full results, statistical tests, and their interpretation are
detailed step by step throughout the notebook, with limitations discussed
both inline and in the dedicated "Limitations" section near the end.

## Requirements

Python 3.10 or later, with pandas, numpy, matplotlib, scipy,
scikit-posthocs, ViennaRNA, and openpyxl installed.

## Data

Two files must be downloaded and placed in the same directory as the
notebook.

1. The MaveDB count file, from
   https://www.mavedb.org/score-sets/urn:mavedb:00000070-a-3/
   Save it as `urn_mavedb_00000070-a-3_counts.csv`.

2. The original publication's supplementary data (Data S1-S4), from the
   article's page on Oxford Academic (search for Firnberg et al. 2014,
   MBE, doi:10.1093/molbev/msu081, Supplementary Data). Save it as
   `Data_S1-S4.xlsx`.

## Usage

Run all cells in order. Steps depend on variables defined in earlier
cells.

## References

Firnberg E, Labonte JW, Gray JJ, Ostermeier M. A comprehensive,
high-resolution map of a gene's fitness landscape. Mol Biol Evol.
2014;31(6):1581-1592. doi:10.1093/molbev/msu081

MaveDB. urn:mavedb:00000070-a-3.
https://www.mavedb.org/score-sets/urn:mavedb:00000070-a-3/

See the notebook's own References section for the remaining sources.
