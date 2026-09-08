"""Convert a GenBank file into FASTA files, one entry per annotated gene.

Each gene is extracted from the genome sequence, reverse-complemented when it
lies on the antisense strand, and written in FASTA format. The extracted
sequence length is checked against the size declared on the LOCUS line, since a
silent off-by-one would shift every gene.

Usage:
    python genbank2fasta.py genome.gbk

Output:
    genes.fasta, containing one FASTA entry per gene, with headers of the form
    >organism|gene_n|start|end|strand
"""

import re
import sys
from pathlib import Path


BASE_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
FASTA_LINE_WIDTH = 80


def read_file(filename):
    """Read a text file and return its content as a list of stripped lines.

    Args:
        filename (str): path to the file to read.

    Returns:
        list of str: the lines of the file, without trailing whitespace.
    """
    with open(filename, "r") as f_in:
        return [line.strip() for line in f_in]


def extract_organism(lines):
    """Extract the organism name from the ORGANISM line of a GenBank file.

    Args:
        lines (list of str): content of the GenBank file.

    Returns:
        str: the organism name.

    Raises:
        ValueError: if no ORGANISM line is found.
    """
    organism_pattern = re.compile(r"ORGANISM\s+(.+)")
    for line in lines:
        match = organism_pattern.search(line)
        if match:
            return match.group(1).strip()
    raise ValueError("No ORGANISM line found in the file.")


def extract_declared_length(lines):
    """Read the sequence length declared on the LOCUS line.

    Args:
        lines (list of str): content of the GenBank file.

    Returns:
        int: the number of base pairs declared in the file header.

    Raises:
        ValueError: if no LOCUS line carrying a length is found.
    """
    locus_pattern = re.compile(r"LOCUS\s+\S+\s+(\d+)\s+bp")
    for line in lines:
        match = locus_pattern.search(line)
        if match:
            return int(match.group(1))
    raise ValueError("No LOCUS line with a sequence length found in the file.")


def find_genes(lines):
    """Collect the coordinates and strand of every annotated gene.

    Both sense genes (``gene  58..272``) and antisense genes
    (``gene  complement(55979..56935)``) are recognised. The ``<`` and ``>``
    symbols marking partial genes are ignored.

    Args:
        lines (list of str): content of the GenBank file.

    Returns:
        list of list: one entry per gene, as [start, end, strand], where
        start and end are integers and strand is "sense" or "antisense".
    """
    genes = []
    gene_pattern = re.compile(r"gene\s+(complement\()?<?(\d+)\.\.>?(\d+)")
    for line in lines:
        match = gene_pattern.search(line)
        if match:
            strand = "antisense" if match.group(1) else "sense"
            genes.append([int(match.group(2)), int(match.group(3)), strand])
    return genes


def count_strands(genes):
    """Count how many genes lie on each strand.

    Args:
        genes (list of list): one entry per gene, as [start, end, strand].

    Returns:
        tuple of int: number of sense genes, number of antisense genes.
    """
    sense = sum(1 for gene in genes if gene[2] == "sense")
    return sense, len(genes) - sense


def extract_sequence(lines):
    """Extract the genome sequence from a GenBank file.

    Sequence lines are identified by their layout: a position number followed
    by blocks of nucleotides. Spaces are removed and the sequence is returned
    in upper case.

    Args:
        lines (list of str): content of the GenBank file.

    Returns:
        str: the full genome sequence.
    """
    sequence_pattern = re.compile(r"^\d+\s+([gatc ]+)")
    sequence = ""
    for line in lines:
        match = sequence_pattern.search(line)
        if match:
            sequence += match.group(1)
    return sequence.replace(" ", "").upper()


def reverse_complement(dna):
    """Return the reverse complement of a DNA sequence.

    Args:
        dna (str): DNA sequence, in upper or lower case.

    Returns:
        str: the reverse complement, in upper case.
    """
    complement = "".join(BASE_COMPLEMENT[base.upper()] for base in dna)
    return complement[::-1]


def write_fasta(filename, header, sequence, mode="w"):
    """Write a sequence to a file in FASTA format.

    Args:
        filename (str): path to the output file.
        header (str): header line, written after the ">" character.
        sequence (str): the sequence to write.
        mode (str): "w" to overwrite the file, "a" to append to it.
    """
    with open(filename, mode) as f_out:
        f_out.write(f">{header}\n")
        for start in range(0, len(sequence), FASTA_LINE_WIDTH):
            f_out.write(sequence[start:start + FASTA_LINE_WIDTH] + "\n")


def extract_genes(genes, genome, organism, single_file=None):
    """Write every gene sequence in FASTA format.

    GenBank coordinates are 1-based and inclusive, so the start position is
    shifted by one to match Python indexing.

    Args:
        genes (list of list): one entry per gene, as [start, end, strand],
            with 1-based inclusive coordinates.
        genome (str): the full genome sequence.
        organism (str): organism name, used in the FASTA headers.
        single_file (str): if given, all genes are written to this one file;
            otherwise one file per gene is created.
    """
    if single_file:
        Path(single_file).unlink(missing_ok=True)

    for number, (start, end, strand) in enumerate(genes, start=1):
        if single_file:
            filename = single_file
            mode = "a"
        else:
            filename = f"gene_{number}.fasta"
            mode = "w"

        sequence = genome[start - 1:end]
        if strand == "antisense":
            sequence = reverse_complement(sequence)

        header = f"{organism}|gene_{number}|{start}|{end}|{strand}"
        write_fasta(filename, header, sequence, mode)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("Usage: python genbank2fasta.py genome.gbk")
    if not Path(sys.argv[1]).exists():
        sys.exit(f"Error: file {sys.argv[1]} not found.")

    content = read_file(sys.argv[1])
    print(f"{len(content)} lines read")

    organism = extract_organism(content)
    print(f"Organism: {organism}")

    genes = find_genes(content)
    sense, antisense = count_strands(genes)
    print(f"{len(genes)} genes found: {sense} sense, {antisense} antisense")

    genome = extract_sequence(content)
    declared_length = extract_declared_length(content)
    if len(genome) != declared_length:
        sys.exit(f"Error: extracted {len(genome)} bases, "
                 f"but {declared_length} are declared on the LOCUS line.")
    print(f"Genome length: {len(genome)} bases (matches the LOCUS line)")

    extract_genes(genes, genome, organism, single_file="genes.fasta")
    print(f"{len(genes)} genes written to genes.fasta")
