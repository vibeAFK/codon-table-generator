"""Compare my codon table (codon_table.py) with the reference table (reference_table.py)."""

from codon_table import CODON_TABLE
from reference_table import STANDARD_CODE

# One-letter code -> the full name used in CODON_TABLE.
FULL_NAMES = {
    "F": "Phenylalanine", "L": "Leucine", "I": "Isoleucine", "M": "Methionine",
    "V": "Valine", "S": "Serine", "P": "Proline", "T": "Threonine",
    "A": "Alanine", "Y": "Tyrosine", "H": "Histidine", "Q": "Glutamine",
    "N": "Asparagine", "K": "Lysine", "D": "Aspartic acid", "E": "Glutamic acid",
    "C": "Cysteine", "W": "Tryptophan", "R": "Arginine", "G": "Glycine",
    "*": "Stop",
}


def flip_reference():
    """Turn the reference (amino acid -> codons) into codon -> amino acid name."""
    flipped = {}
    for letter, codons in STANDARD_CODE.items():
        for codon in codons:
            flipped[codon] = FULL_NAMES[letter]
    return flipped


def main():
    reference = flip_reference()

    print("Structure of my table:        codon -> amino acid name")
    print("Structure of reference table: amino acid letter -> list of codons")
    print()
    print("Codons in my table:        ", len(CODON_TABLE))
    print("Codons in reference table: ", len(reference))

    missing = sorted(set(reference) - set(CODON_TABLE))
    extra = sorted(set(CODON_TABLE) - set(reference))
    mismatches = [c for c in sorted(CODON_TABLE)
                  if c in reference and CODON_TABLE[c] != reference[c]]

    print("Missing from my table:     ", missing if missing else "none")
    print("Not in reference table:    ", extra if extra else "none")
    for codon in mismatches:
        print("Mismatch:", codon, "mine =", CODON_TABLE[codon], "| reference =", reference[codon])
    if not mismatches:
        print("Mismatching amino acids:    none")

    print()
    print("Tables match:", not (missing or extra or mismatches))


if __name__ == "__main__":
    main()
