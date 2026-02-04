# Codon Graph Visualization

A visualization tool for exploring mappings from codons to amino acids (colors) using Hamming distance graphs.

## Overview

This project models the relationship between genetic codons and their encoded amino acids as a graph where:
- **Nodes** represent individual codons
- **Edges** connect codons that differ by exactly one nucleotide (Hamming distance = 1)
- **Node colors** represent the assigned amino acid

This structure is useful for studying:
- Robustness of genetic code assignments
- How mutations (single nucleotide changes) affect amino acid translation
- Properties of codon partitioning schemes

## Features

- **Configurable parameters**: Adjust alphabet size, codon length, and number of amino acids
- **Flexible assignment modes**: Random or lexicographic codon-to-amino acid assignment
- **Automatic partitioning**: Evenly distribute codons among amino acids
- **Biological scale option**: Pre-configured for real DNA codons (4 nucleotides, triplets, 21 amino acids)
- **Graph statistics**: Compute robustness metrics (intra-amino-acid edge fractions)
- **Interactive visualization**: NetworkX + Matplotlib visualization with legend

## Installation

1. Clone the repository
2. Create a virtual environment (optional but recommended):
   ```bash
   python -m venv myenv
   myenv\Scripts\activate  # Windows
   # or: source myenv/bin/activate  # Unix/macOS
   ```
3. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

Run the visualization:
```bash
python codon_graph.py
```

### Configuration

Edit the CONFIG section at the top of `codon_graph.py` to customize:

```python
ALPHABET_SIZE = 3           # Number of unique nucleotides
CODON_LENGTH = 2            # Length of each codon
NUM_AMINO_ACIDS = 4         # Number of amino acid classes
PARTITION = [3, 3, 2, 1]    # Distribution of codons per amino acid (optional)
ASSIGNMENT_MODE = 'random'  # 'random' or 'lexicographic'
SEED = 18                   # Random seed for reproducibility
SHOW_NODE_LABELS = True     # Display codon labels on graph
```

#### Biological Scale

Set `BIOLOGICAL_SCALE = True` to use standard DNA parameters:
- 4 nucleotides (A, T, G, C)
- Triplet codons (64 total)
- 21 amino acids (standard genetic code partition)

### Example Output

The script generates a colored network graph and prints robustness statistics:
```
Intra-amino-acid edges: 12/18 (0.667)
Average robustness over 1 runs: 0.667
```

Higher robustness indicates that most edge connections remain within the same amino acid class.

## Algorithm

1. Generate all possible codons from the alphabet
2. Partition codons among amino acids according to the specified distribution
3. Build a graph with Hamming distance = 1 edges
4. Assign colors based on amino acid partition
5. Visualize and compute robustness metrics

## Requirements

- Python 3.7+
- networkx
- matplotlib
- numpy

## Author

William Anklam
