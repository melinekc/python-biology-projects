# genbank2fasta

Extract every annotated gene from a GenBank file and write the sequences in
FASTA format.

## What it does

A GenBank file bundles annotations and genome sequence in a single text file.
This script separates the two: it locates every `gene` feature, retrieves the
corresponding stretch of the genome, reverse-complements it when the gene lies
on the antisense strand, and writes the result as FASTA.

Sense and antisense features are written differently in GenBank:

```
     gene            58..272
     gene            complement(55979..56935)
```

Both forms are recognised, as are the partial variants marked `<` and `>`.

## Usage

```
python genbank2fasta.py NC_001133.gbk
```

All genes are written to a single multi-FASTA file, `genes.fasta`, with headers
of the form:

```
>organism|gene_n|start|end|strand
```

Passing `single_file=None` to `extract_genes()` writes one file per gene
instead.

## Example

Run on `NC_001133.gbk`, chromosome I of *Saccharomyces cerevisiae* S288C
(230,218 bp, downloaded from NCBI):

```
$ python genbank2fasta.py NC_001133.gbk
...
```

The extracted genome length is checked against the size declared on the first
line of the GenBank file. A silent off-by-one in the sequence parsing would
otherwise shift every gene.

## Notes and limitations

- GenBank coordinates are 1-based and inclusive; the script converts them to
  Python slicing.
- Genes flagged `<` or `>` are **partial**: the start or stop codon falls
  outside the assembled region. Their sequences are extracted as annotated, so
  they are not complete reading frames and should not be translated as such.
  This information is not carried into the FASTA headers.
- Only `gene` features are extracted; CDS, tRNA and other feature types are
  ignored.
- The genome sequence is identified by line layout rather than by the
  `ORIGIN` / `//` delimiters. This works on standard GenBank files.

## Dependencies

Standard library only (`re`, `sys`, `pathlib`).

## Context

This project follows exercise 27.2.3 of *Cours de Python pour la biologie*
(Fuchs & Poulain, Université Paris Cité). The implementation is my own, with
three additions beyond the original exercise:

- the extracted sequence length is checked against the size declared on the
  LOCUS line, so that a parsing error stops the script instead of silently
  shifting every gene;
- genes can be written either to one file each or to a single multi-FASTA
  file, which is the format most downstream tools expect;
- sense and antisense genes are counted separately, as a quick sanity check on
  the annotation parsing.
