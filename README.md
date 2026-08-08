# Combinatorial Robustness of the Genetic Code

This repository contains the Python code used to study the standard genetic code as a Hamming graph. Codons are the 64 vertices of `H(3,4)`, and two codons are connected when they differ in exactly one nucleotide. An edge is silent when its two codons have the same assigned output.

The calculation is deliberately unweighted: it does not use amino-acid similarity, mutation probabilities, codon usage, or substitution matrices.

## Computations

[`codon_robustness.py`](codon_robustness.py) implements the computational work described in the paper:

1. It constructs the 64 codons and lists every unordered pair at Hamming distance one, producing the 288 edges of `H(3,4)`.
2. It encodes the standard genetic code and counts its silent edges.
3. It calculates the exact expected robustness for a uniformly random assignment with the biological degeneracy profile.
4. It can generate and analyze random fixed-profile reassignments directly.
5. It examines all unordered codon pairs, excludes synonymous pairs, and reconstructs and recounts every nontrivial degeneracy-preserving swap.
6. It applies the two improving swaps found by the search and repeats the exhaustive search from the resulting assignment.

The genetic-code assignment is included in the program, so no external dataset or preprocessing step is required.

## Running the analysis

The main calculation uses only the Python standard library and requires Python 3.9 or newer.

```bash
python codon_robustness.py
```

Machine-readable output is available with:

```bash
python codon_robustness.py --json
```

To perform random reassignment trials, specify the number of trials and a seed:

```bash
python codon_robustness.py --random-trials 10000 --seed 2026
```

Each trial draws uniformly from assignments with the biological degeneracy profile. The seed makes the sample reproducible.

## Visualization

[`visualize_code.py`](visualize_code.py) draws the standard genetic code on the Hamming graph. Nodes are colored by output, silent edges are darkened, and nonsilent edges are shown faintly.

```bash
python -m pip install -r requirements.txt
python visualize_code.py --output biological_codon_graph.png --labels
```

The generated image is not required for the numerical calculations and is not stored in the repository.

## Repository contents

- `codon_robustness.py`: graph construction, robustness calculation, random reassignment, and exhaustive codon-swap search.
- `visualize_code.py`: optional visualization of the biological assignment.
- `requirements.txt`: packages used only for visualization.

## Authors

William Anklam, Naomi Hekman, and Michael Sill

## License

The code in this repository is available under the [MIT License](LICENSE).
