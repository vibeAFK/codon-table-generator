# Codon Table Generator

A Python program that builds the standard genetic code (all 64 codons) and translates a DNA or RNA sequence into its amino acid sequence, codon by codon. It was written for Python Assessment 2 (MSc Bioinformatics and Biotechnology, Chanakya University).

**Author:** Vaibhav S. Vijayan

## Objective

Generate a codon table and identify the amino acid for each codon of a sequence entered by the user, marking start and stop codons and handling incomplete codons at the end of the sequence.

## Features

- Dictionary of all 64 codons and their amino acids (standard genetic code), with a built-in check that it has exactly 64 entries
- Input validation: only A, T, G, C (DNA) or A, U, G, C (RNA); empty input, invalid characters and mixed T/U are rejected with a clear message and the user is asked again
- Splits the sequence into codons of three nucleotides
- Translates every codon and keeps going after a stop codon, marking each stop
- Identifies the start codon (ATG) and the stop codons (TAA, TAG, TGA)
- Reports an incomplete codon (1-2 leftover nucleotides) instead of translating it
- Prints the full codon table as a grid and the result as a table
- Additional challenge: RNA input is detected automatically (U instead of T) and translated; codons are shown in their RNA form (AUG, GCU, ...)
- Compares the codon table against an independent reference table (`compare_tables.py`)
- Unit tests (`test_codon_table.py`)

## How the genetic code works

DNA is read in groups of three nucleotides called codons. Four bases in three positions give 4 x 4 x 4 = 64 codons. They code for 20 amino acids plus stop signals, so most amino acids are encoded by more than one codon. ATG codes for methionine and is the usual start codon; TAA, TAG and TGA tell the ribosome to stop. In RNA, thymine (T) is replaced by uracil (U), so ATG becomes AUG.

## Requirements

- Python 3.8 or newer
- No extra packages (standard library only)

## How to run

```
python codon_table.py
```

The program prints the codon table, then asks for a sequence. To compare the tables and run the tests:

```
python compare_tables.py
python -m unittest -v
```

## Sample run (DNA)

Input: `ATGGCTTTTTAA`

```
Position  Codon   Amino acid      Note
----------------------------------------------
1         ATG     Methionine      Start codon
2         GCT     Alanine         
3         TTT     Phenylalanine   
4         TAA     Stop            Stop codon

Codons: ATG - GCT - TTT - TAA
Amino acids: Methionine - Alanine - Phenylalanine - Stop
Start codon: ATG
Stop codon: TAA
```

## RNA mode (additional challenge)

Input: `AUGGCUUUUUAA`

```
Detected sequence type: RNA

Position  Codon   Amino acid      Note
----------------------------------------------
1         AUG     Methionine      Start codon
2         GCU     Alanine         
3         UUU     Phenylalanine   
4         UAA     Stop            Stop codon

Codons: AUG - GCU - UUU - UAA
Amino acids: Methionine - Alanine - Phenylalanine - Stop
Start codon: AUG
Stop codon: UAA
```

## Edge cases handled

Incomplete codon at the end (`ATGGCTTTTTAAGC`): the two leftover bases are reported and not translated.

```
Codons: ATG - GCT - TTT - TAA
Amino acids: Methionine - Alanine - Phenylalanine - Stop
Start codon: ATG
Stop codon: TAA
Incomplete codon, not translated: GC
```

Invalid input is rejected and the program asks again:

```
Error: Invalid character(s): X. Only A, C, G, T are allowed for DNA.
Error: The sequence contains both T and U. Use DNA (T) or RNA (U), not both.
Error: The sequence is empty.
```

Lower case and spaces are accepted (`atg gct` is read as `ATGGCT`). If a sequence has no stop codon, the stop codon line shows `None`.

## Tests

```
python -m unittest -v
```

Nine tests cover the sample input, lower case input, invalid characters, empty input, mixed T and U, incomplete codons, RNA input, the 64-codon check and the three stop codons. Result:

```
Ran 9 tests in 0.001s

OK
```

## Comparison with a reference codon table

`reference_table.py` holds an independent copy of the standard genetic code (NCBI translation table 1). It follows the structure used by Biopython's `CodonTable` module (`Bio/Data/CodonTable.py` in the biopython repository on GitHub): one-letter amino acid code mapped to a list of codons. It was written out by hand and is not downloaded or imported from Biopython.

My table is structured the other way round: codon mapped to the full amino acid name. `compare_tables.py` converts the reference into my structure and checks both tables codon by codon:

```
Structure of my table:        codon -> amino acid name
Structure of reference table: amino acid letter -> list of codons

Codons in my table:         64
Codons in reference table:  64
Missing from my table:      none
Not in reference table:     none
Mismatching amino acids:    none

Tables match: True
```

The two layouts hold the same information. The codon-to-amino-acid layout is faster for translating a sequence (one lookup per codon), while the amino-acid-to-codons layout is better for seeing which codons share an amino acid (synonymous codons).

## Project structure

```
codon-table-generator/
|-- codon_table.py        main program
|-- reference_table.py    reference table (standard genetic code)
|-- compare_tables.py     compares my table with the reference
|-- test_codon_table.py   unit tests
|-- README.md
|-- .gitignore
```

## Learning outcomes

Python dictionaries, loops, conditional statements, string manipulation, input validation, and basic biological sequence analysis.
