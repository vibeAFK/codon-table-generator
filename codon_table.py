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
