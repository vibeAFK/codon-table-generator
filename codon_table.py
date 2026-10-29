"""
Codon Table Generator
Translates a DNA (or RNA) sequence into amino acids using the standard genetic code.
"""

# The standard genetic code: every DNA codon mapped to its amino acid.
CODON_TABLE = {
    "TTT": "Phenylalanine",
    "TTC": "Phenylalanine",
    "TTA": "Leucine",
    "TTG": "Leucine",
    "TCT": "Serine",
    "TCC": "Serine",
    "TCA": "Serine",
    "TCG": "Serine",
    "TAT": "Tyrosine",
    "TAC": "Tyrosine",
    "TAA": "Stop",
    "TAG": "Stop",
    "TGT": "Cysteine",
    "TGC": "Cysteine",
    "TGA": "Stop",
    "TGG": "Tryptophan",
    "CTT": "Leucine",
    "CTC": "Leucine",
    "CTA": "Leucine",
    "CTG": "Leucine",
    "CCT": "Proline",
    "CCC": "Proline",
    "CCA": "Proline",
    "CCG": "Proline",
    "CAT": "Histidine",
    "CAC": "Histidine",
    "CAA": "Glutamine",
    "CAG": "Glutamine",
    "CGT": "Arginine",
    "CGC": "Arginine",
    "CGA": "Arginine",
    "CGG": "Arginine",
    "ATT": "Isoleucine",
    "ATC": "Isoleucine",
    "ATA": "Isoleucine",
    "ATG": "Methionine",
    "ACT": "Threonine",
    "ACC": "Threonine",
    "ACA": "Threonine",
    "ACG": "Threonine",
    "AAT": "Asparagine",
    "AAC": "Asparagine",
    "AAA": "Lysine",
    "AAG": "Lysine",
    "AGT": "Serine",
    "AGC": "Serine",
    "AGA": "Arginine",
    "AGG": "Arginine",
    "GTT": "Valine",
    "GTC": "Valine",
    "GTA": "Valine",
    "GTG": "Valine",
    "GCT": "Alanine",
    "GCC": "Alanine",
    "GCA": "Alanine",
    "GCG": "Alanine",
    "GAT": "Aspartic acid",
    "GAC": "Aspartic acid",
    "GAA": "Glutamic acid",
    "GAG": "Glutamic acid",
    "GGT": "Glycine",
    "GGC": "Glycine",
    "GGA": "Glycine",
    "GGG": "Glycine",
}

# Safety check: there must be exactly 64 codons (4 x 4 x 4).
assert len(CODON_TABLE) == 64, "The codon table must contain 64 codons"


VALID_DNA = set("ATGC")


def validate_sequence(raw):
    """Clean the input and check that it only contains A, T, G and C.

    Returns the cleaned sequence in upper case.
    Raises ValueError if the input is empty or has invalid characters.
    """
    seq = "".join(raw.split()).upper()  # remove all whitespace, make upper case
    if not seq:
        raise ValueError("The sequence is empty.")
    bad = sorted(set(seq) - VALID_DNA)
    if bad:
        raise ValueError("Invalid character(s): " + ", ".join(bad)
                         + ". Only A, T, G and C are allowed.")
    return seq


def split_codons(seq):
    """Split a sequence into complete codons (groups of three nucleotides)."""
    codons = []
    for i in range(0, len(seq) - 2, 3):
        codons.append(seq[i:i + 3])
    return codons


def translate(codons):
    """Translate each codon into its amino acid.

    Returns a list of (position, codon, amino_acid) tuples.
    """
    results = []
    for position, codon in enumerate(codons, start=1):
        results.append((position, codon, CODON_TABLE[codon]))
    return results


START_CODON = "ATG"
STOP_CODONS = ("TAA", "TAG", "TGA")


def mark_codons(results):
    """Add a note to each result: start codon, stop codon or nothing."""
    marked = []
    for position, codon, amino_acid in results:
        note = ""
        if codon == START_CODON:
            note = "Start codon"
        elif codon in STOP_CODONS:
            note = "Stop codon"
        marked.append((position, codon, amino_acid, note))
    return marked


def get_leftover(seq):
    """Return the 1 or 2 trailing nucleotides that do not form a full codon."""
    return seq[len(seq) - len(seq) % 3:]


# Three-letter abbreviations used in the codon table grid.
ABBREVIATIONS = {
    "Phenylalanine": "Phe", "Leucine": "Leu", "Isoleucine": "Ile",
    "Methionine": "Met", "Valine": "Val", "Serine": "Ser", "Proline": "Pro",
    "Threonine": "Thr", "Alanine": "Ala", "Tyrosine": "Tyr", "Histidine": "His",
    "Glutamine": "Gln", "Asparagine": "Asn", "Lysine": "Lys",
    "Aspartic acid": "Asp", "Glutamic acid": "Glu", "Cysteine": "Cys",
    "Tryptophan": "Trp", "Arginine": "Arg", "Glycine": "Gly", "Stop": "Stop",
}


def display_table():
    """Print all 64 codons as a grid (rows = 1st base, columns = 2nd base)."""
    bases = "TCAG"
    print("\nThe Standard Genetic Code (DNA codons)")
    print("=" * 65)
    print("1st    " + "      ".join("{:<8}".format("2nd: " + b) for b in bases) + "     3rd")
    print("-" * 65)
    for first in bases:
        for third in bases:
            cells = []
            for second in bases:
                codon = first + second + third
                cells.append("{} {:<4}".format(codon, ABBREVIATIONS[CODON_TABLE[codon]]))
            label = first if third == "T" else " "
            print(" {}     ".format(label) + "      ".join(cells) + "     " + third)
        print("-" * 65)
