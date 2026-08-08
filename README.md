# Combinatorial Robustness of the Genetic Code

This repository contains the computational part of our study of the standard genetic code as a Hamming graph. A codon is represented by a vertex of `H(3,4)`, and two codons are adjacent when they differ in exactly one nucleotide. Once each codon is assigned to an amino acid or stop signal, an edge is called *silent* when both endpoints have the same output.

The main program builds this graph directly from the four-letter DNA alphabet and studies how the arrangement of synonymous codons affects the number of silent edges. No biological similarity scores or mutation weights are used; the calculation is purely combinatorial.

## What the program does

[`codon_robustness.py`](codon_robustness.py) carries out the following calculations:

1. Constructs all 64 codons and the 288 edges of `H(3,4)`.
2. Encodes the standard genetic code, including the stop signal, and checks its degeneracy profile.
3. Counts the silent edges in the biological assignment.
4. Examines every swap of two codons with different outputs. These swaps preserve all class sizes.
5. Applies the two improving Ser/Cys swaps and examines the full swap neighborhood of the resulting assignment.
6. Checks the explicit zero-silent-edge assignment given in the paper's appendix.
7. Calculates the random fixed-profile expectation and the classwise upper bound discussed in the paper.

The exhaustive search considers 1,926 nontrivial codon swaps. The graph and genetic-code data are defined in the same file, so there are no external datasets or hidden preprocessing steps.

## Running the analysis

Python 3.9 or newer is sufficient; the numerical analysis uses only the standard library.

```bash
python codon_robustness.py
```

The program prints the graph size, random expectation, biological robustness, improving swaps, and the result of the second local search. For use by another program, the same information can be written as JSON:

```bash
python codon_robustness.py --json
```

## Biological-code visualization

[`visualize_code.py`](visualize_code.py) draws the standard genetic code on the Hamming graph. Nodes are colored by output, silent edges are darkened, and nonsilent edges are shown faintly.

```bash
python -m pip install -r requirements.txt
python visualize_code.py --output biological_codon_graph.png --labels
```

The image is generated from the same codon assignment used in the analysis. It is not needed for the numerical calculations and is therefore not stored in the repository.

## Repository contents

- `codon_robustness.py` contains the genetic code, graph construction, robustness calculation, and exhaustive swap search.
- `visualize_code.py` produces the optional graph visualization.
- `requirements.txt` lists the packages used only for plotting.

Exploratory toy models, non-degeneracy-preserving reassignments, presentation scripts, virtual environments, and generated images are not part of this repository.

## Author

William Anklam
