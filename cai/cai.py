"""Calculate the Codon Adaptation Index (CAI) of coding sequences.

The CAI measures how closely the codon usage of a coding sequence matches the
preferred codons of a reference organism (Sharp & Li, 1987). Values range from
0 to 1, where 1 means every codon is the most frequent one of its synonymous
family.

Following the original definition, single-codon amino acids (Met, Trp) and stop
codons are excluded by default, since their weight is fixed and carries no
information about codon choice. Both can be included, which is one of the
reasons why tools disagree on the CAI of a given sequence.

The module also provides the building blocks needed to explore the space of
synonymous variants of a protein: translation, a synonym table and a random
variant generator.

Usage:
    python cai.py codon_table.txt sequences.fasta

Output:
    One line per sequence: identifier, tab, CAI rounded to four decimals.
    Sequences that cannot be scored are reported and skipped.
"""

import math
import random
import re
import sys
from collections import defaultdict
from pathlib import Path


VALID_BASES = set("AUGC")
EXPECTED_CODONS = 64
CODON_LENGTH = 3


def read_kazusa_table(filename):
    """Read a codon usage table in Kazusa format.

    The expected layout is one entry per codon, as produced by the Kazusa
    database with the "Codon Usage Table with Amino Acids" option:
    ``UUU F 0.46 17.6 (714298)``. Only the codon, the amino acid and the
    fraction are used; header lines are ignored.

    Args:
        filename (str): path to the table file.

    Returns:
        dict: codon (str) mapped to a tuple (amino acid as a one-letter code,
        fraction of use within its synonymous family).

    Raises:
        ValueError: if the table does not contain exactly 64 codons.
    """
    codon_pattern = re.compile(r"([AUGC]{3}) ([A-Z*]) ([\d.]+)")
    with open(filename, "r", encoding="utf-8") as f_in:
        entries = codon_pattern.findall(f_in.read())

    table = {codon: (amino_acid, float(fraction))
             for codon, amino_acid, fraction in entries}

    if len(table) != EXPECTED_CODONS:
        raise ValueError(f"Expected {EXPECTED_CODONS} codons, "
                         f"found {len(table)}.")
    return table


def create_aa_table(codon_table):
    """Summarise a codon table by amino acid.

    Args:
        codon_table (dict): codon mapped to (amino acid, fraction).

    Returns:
        dict: amino acid (str) mapped to a tuple (highest fraction within the
        family, number of synonymous codons).
    """
    families = defaultdict(list)
    for amino_acid, fraction in codon_table.values():
        families[amino_acid].append(fraction)

    return {amino_acid: (max(fractions), len(fractions))
            for amino_acid, fractions in families.items()}


def create_weight_table(codon_table):
    """Compute the relative adaptiveness of every codon.

    The weight of a codon is its fraction of use divided by the fraction of the
    most frequent codon of the same amino acid, so the preferred codon of each
    family has a weight of 1.

    Args:
        codon_table (dict): codon mapped to (amino acid, fraction).

    Returns:
        dict: codon (str) mapped to a tuple (amino acid, weight, number of
        synonymous codons).
    """
    aa_table = create_aa_table(codon_table)
    weights = {}
    for codon, (amino_acid, fraction) in codon_table.items():
        max_fraction, family_size = aa_table[amino_acid]
        weights[codon] = (amino_acid, fraction / max_fraction, family_size)
    return weights


def create_synonym_table(codon_table):
    """List the synonymous codons of every amino acid.

    Args:
        codon_table (dict): codon mapped to (amino acid, fraction).

    Returns:
        dict: amino acid (str) mapped to the list of codons encoding it. Stop
        codons are grouped under "*".
    """
    families = defaultdict(list)
    for codon, (amino_acid, _) in codon_table.items():
        families[amino_acid].append(codon)
    return dict(families)


def read_fasta(filename):
    """Read one or several sequences from a FASTA file.

    Args:
        filename (str): path to the FASTA file.

    Returns:
        dict: header without the ">" character, mapped to its sequence.
    """
    sequences = {}
    seq_id = None
    seq = ""
    with open(filename, "r", encoding="utf-8") as f_in:
        for line in f_in:
            if line.startswith(">"):
                if seq_id:
                    sequences[seq_id] = seq
                seq_id = line[1:].strip()
                seq = ""
            else:
                seq += line.strip()
    if seq_id:
        sequences[seq_id] = seq
    return sequences


def check_sequence(sequence):
    """Normalise a coding sequence and check that it can be read as codons.

    Thymines are converted to uracils, so DNA sequences are accepted.

    Args:
        sequence (str): RNA or DNA coding sequence.

    Returns:
        str: the sequence in upper case, with T replaced by U.

    Raises:
        ValueError: if the length is not a multiple of three, or if the
            sequence contains characters other than A, U, G, C and T.
    """
    if len(sequence) % CODON_LENGTH != 0:
        raise ValueError(f"Sequence length ({len(sequence)}) "
                         f"is not a multiple of {CODON_LENGTH}.")

    sequence = sequence.upper().replace("T", "U")
    unknown = set(sequence) - VALID_BASES
    if unknown:
        raise ValueError(f"Invalid characters in sequence: {unknown}")
    return sequence


def translate(sequence, codon_table):
    """Translate a coding sequence into its protein sequence.

    Args:
        sequence (str): RNA or DNA coding sequence.
        codon_table (dict): codon mapped to (amino acid, fraction).

    Returns:
        str: the protein sequence in one-letter code, stop codons written "*".

    Raises:
        ValueError: if the sequence is not a valid coding sequence.
    """
    sequence = check_sequence(sequence)
    protein = []
    for start in range(0, len(sequence), CODON_LENGTH):
        codon = sequence[start:start + CODON_LENGTH]
        amino_acid, _ = codon_table[codon]
        protein.append(amino_acid)
    return "".join(protein)


def calculate_cai(sequence, weights_table, stop=False, unique=False):
    """Calculate the CAI of a coding sequence.

    The CAI is the geometric mean of the weights of the codons kept in the
    calculation. It is computed through logarithms, since multiplying hundreds
    of values below 1 would underflow.

    Args:
        sequence (str): RNA or DNA coding sequence.
        weights_table (dict): codon mapped to (amino acid, weight, number of
            synonymous codons), as returned by create_weight_table().
        stop (bool): if True, stop codons are included in the calculation.
        unique (bool): if True, amino acids encoded by a single codon are
            included in the calculation.

    Returns:
        float: the CAI, between 0 and 1.

    Raises:
        ValueError: if the sequence is not a valid coding sequence, or if no
            codon is left once the exclusions are applied.
    """
    sequence = check_sequence(sequence)

    sum_logs = 0
    scored_codons = 0
    for start in range(0, len(sequence), CODON_LENGTH):
        codon = sequence[start:start + CODON_LENGTH]
        amino_acid, weight, family_size = weights_table[codon]
        if family_size == 1 and not unique:
            continue
        if amino_acid == "*" and not stop:
            continue
        sum_logs += math.log(weight)
        scored_codons += 1

    if scored_codons == 0:
        raise ValueError("No codon left once the exclusions are applied.")
    return math.exp(sum_logs / scored_codons)


def random_variant(protein, synonym_table, weight_table, bias):
    """Build a random synonymous coding sequence for a protein.

    At every position, the codon is drawn among the synonyms of the amino acid,
    with a probability proportional to its weight raised to the power of
    bias. A bias of 0 gives a uniform draw over the synonyms. A positive bias
    favours the preferred codon of each family, pushing the CAI up; a negative
    bias favours the least frequent codons instead, pushing the CAI down.
    Stop codons are drawn like any other family.

    Args:
        protein (str): protein sequence in one-letter code.
        synonym_table (dict): amino acid mapped to its list of codons, as
            returned by create_synonym_table().
        weight_table (dict): codon mapped to (amino acid, weight, number of
            synonymous codons), as returned by create_weight_table().
        bias (float): exponent applied to codon weights before sampling.
            Zero gives a uniform draw; positive values favour the preferred
            codon, negative values favour the rarest ones.

    Returns:
        str: a coding sequence translating into the given protein.
    """
    sequence = []
    for amino_acid in protein:
        codons = synonym_table[amino_acid]
        weights = [weight_table[c][1] ** bias for c in codons]
        sequence.append(random.choices(codons, weights=weights)[0])
    return "".join(sequence)


def gc_content(sequence):
    """Compute the GC content of a nucleotide sequence.

    Args:
        sequence (str): nucleotide sequence, RNA or DNA.

    Returns:
        float: the fraction of G and C bases, between 0 and 1.

    Raises:
        ValueError: if the sequence is empty.
    """
    if not sequence:
        raise ValueError("Cannot compute GC content of an empty sequence.")
    sequence = sequence.upper()
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("Usage: python cai.py codon_table.txt sequences.fasta")
    for argument in sys.argv[1:]:
        if not Path(argument).exists():
            sys.exit(f"Error: file {argument} not found.")

    codon_table = read_kazusa_table(sys.argv[1])
    weights_table = create_weight_table(codon_table)
    sequences = read_fasta(sys.argv[2])

    for identifier, coding_sequence in sequences.items():
        try:
            cai = calculate_cai(coding_sequence, weights_table)
            print(f"{identifier}\t{cai:.4f}")
        except ValueError as error:
            print(f"{identifier}\tskipped: {error}")
