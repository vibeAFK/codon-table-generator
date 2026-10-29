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
VALID_RNA = set("AUGC")


def detect_mode(seq):
    """Return "RNA" if the sequence contains U, otherwise "DNA"."""
    return "RNA" if "U" in seq else "DNA"


def validate_sequence(raw):
    """Clean the input and check it is a valid DNA or RNA sequence.

    Returns (sequence, mode) where mode is "DNA" or "RNA".
    Raises ValueError if the input is empty, mixes T and U,
    or contains characters other than A, T/U, G and C.
    """
    seq = "".join(raw.split()).upper()  # remove all whitespace, make upper case
    if not seq:
        raise ValueError("The sequence is empty.")
    if "T" in seq and "U" in seq:
        raise ValueError("The sequence contains both T and U. Use DNA (T) or RNA (U), not both.")
    mode = detect_mode(seq)
    allowed = VALID_RNA if mode == "RNA" else VALID_DNA
    bad = sorted(set(seq) - allowed)
    if bad:
        raise ValueError("Invalid character(s): " + ", ".join(bad)
                         + ". Only " + ", ".join(sorted(allowed)) + " are allowed for " + mode + ".")
    return seq, mode


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


def display_result(results, leftover):
    """Print the translation as a table, followed by a short summary."""
    if not results:
        print("\nNo complete codons found in the sequence.")
    else:
        print("\n{:<10}{:<8}{:<16}{}".format("Position", "Codon", "Amino acid", "Note"))
        print("-" * 46)
        for position, codon, amino_acid, note in results:
            print("{:<10}{:<8}{:<16}{}".format(position, codon, amino_acid, note))

        starts = [r[1] for r in results if r[3] == "Start codon"]
        stops = [r[1] for r in results if r[3] == "Stop codon"]
        print()
        print("Codons:", " - ".join(r[1] for r in results))
        print("Amino acids:", " - ".join(r[2] for r in results))
        print("Start codon:", starts[0] if starts else "None")
        print("Stop codon:", stops[0] if stops else "None")

    if leftover:
        print("Incomplete codon, not translated:", leftover)


def main():
    """Ask for a DNA or RNA sequence, translate it and show the results."""
    display_table()
    while True:
        raw = input("\nEnter a DNA or RNA sequence: ")
        try:
            seq, mode = validate_sequence(raw)
            break
        except ValueError as error:
            print("Error:", error)

    print("Detected sequence type:", mode)
    dna = seq.replace("U", "T")  # the codon table is written with T
    leftover = get_leftover(dna)
    results = mark_codons(translate(split_codons(dna)))

    if mode == "RNA":  # show the codons in their RNA form (U instead of T)
        results = [(p, c.replace("T", "U"), aa, note) for p, c, aa, note in results]
        leftover = leftover.replace("T", "U")

    display_result(results, leftover)


if __name__ == "__main__":
    main()
