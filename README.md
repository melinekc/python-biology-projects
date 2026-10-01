# python-biology-projects

Python tools for parsing and analysing biological sequence data.

This repository collects small, self-contained projects written while moving
from wet-lab process engineering towards computational biology. Each project
sits in its own folder with a short README describing the question it answers,
how to run it, and what it produces.

## Projects

| Folder | Description |
|---|---|
| [`tem1-fitness-landscape`](tem1-fitness-landscape/) | Rebuilds the fitness landscape of TEM-1 beta-lactamase from raw deep mutational scanning counts (18,081 codon-level variants, 13 ampicillin concentrations) and validates it against the published scores (Spearman 0.965). Extends the analysis to synonymous codons and mRNA folding near the start of the gene. |
| [`cai`](cai/) | Implements the codon adaptation index from first principles and asks what the number means. The choice of reference table can change it by more than 0.3, while the exact definition barely matters. Across 1,000 synonymous luciferase variants, CAI and GC content are tightly coupled (Spearman 0.96). |
| [`genbank2fasta`](genbank2fasta/) | Extracts every annotated gene from a GenBank file into FASTA, handling sense and antisense strands, with a check of the extracted genome length against the LOCUS line. Follows an exercise from a Python course for biology, with three additions of my own. |

## Requirements

Python 3.10 or later. Each project lists its own dependencies in its README:
`cai.py` and `genbank2fasta` use only the standard library, while the notebooks
rely on pandas, numpy, scipy and matplotlib, plus a few extra packages for
TEM-1 (scikit-posthocs, ViennaRNA, openpyxl).

## Author

Méline Kuric, molecular biologist and mRNA process engineer.
