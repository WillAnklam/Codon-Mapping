"""
A visualization tool for exploring mappings f: C -> A,
where C is a set of codons and A is a set of amino acids (colors).

Each codon is represented as a vertex, and edges connect codons
that differ by exactly one symbol (Hamming distance = 1).

You can:
 - Control alphabet size and codon length.
 - Choose number of amino acids (classes/colors).
 - Automatically partition codons among amino acids.
 - Visualize the resulting colored graph.
"""

import itertools
import random
import string
import networkx as nx
import matplotlib.pyplot as plt

# ===============================================================
# ------------------------- CONFIG -------------------------------
# ===============================================================

# Set to true to enable full sized calculations
BIOLOGICAL_SCALE = False

# Number of unique letters in the alphabet (e.g., 4 for DNA, 3 for ternary)
ALPHABET_SIZE = 3

# Length of each codon (e.g., 3 -> triplets)
CODON_LENGTH = 2

# total codons = ALPHABET_SIZE ^ CODON_LENGTH

# Number of amino acids / colors to map to (e.g., 4 different colors)
NUM_AMINO_ACIDS = 4

# Optional explicit partition (must sum to total codons).
# If None, a balanced partition will be generated automatically.
PARTITION = [3, 3, 2, 1]  # Example: for 9 codons and 4 amino acids: [3, 3, 2, 1]

# Whether codon assignment should be random or lexicographic
ASSIGNMENT_MODE = 'random'  # 'random' | 'lexicographic'

# Random seed for reproducibility
SEED = 18

# Whether to show codon labels (turn off for large graphs)
SHOW_NODE_LABELS = True

if BIOLOGICAL_SCALE:
    ALPHABET_SIZE = 4  # DNA
    CODON_LENGTH = 3   # Triplets
    NUM_AMINO_ACIDS = 21
    PARTITION = [6, 6, 6, 4, 4, 4, 4, 4, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1]  # Standard biological partition
    ASSIGNMENT_MODE = 'random'
    SHOW_NODE_LABELS = True


# ===============================================================
# --------------------- HELPER FUNCTIONS -------------------------
# ===============================================================

def generate_alphabet(size):
    """Return an alphabet N of the given size."""
    # For size <= 26, use single letters A–Z.
    # For larger alphabets, label them A1, A2, A3, ...
    if size <= 26:
        return list(string.ascii_uppercase[:size])
    return [f"A{i+1}" for i in range(size)]


def generate_codons(alphabet, length):
    """Generate all codons = all possible strings of given length over N."""
    return [''.join(p) for p in itertools.product(alphabet, repeat=length)]


def hamming_distance(x, y):
    """Compute the Hamming distance between two codons (count of differing symbols)."""
    return sum(a != b for a, b in zip(x, y))


def auto_partition(total, k):
    """
    Evenly divide 'total' codons among 'k' amino acids.
    If it doesn't divide evenly, some get one extra.
    """
    base = total // k
    r = total % k
    return [base + 1 if i < r else base for i in range(k)]


def make_labels(k):
    """Return a list of amino-acid labels like ['a1', 'a2', ..., 'ak']."""
    return [f"a{i+1}" for i in range(k)]


def assign_codons_to_amino_acids(codons, partition, aa_labels, mode='random', seed=None):
    """
    Given a set of codons and a partition vector, assign codons to amino acids.
    The partition defines how many codons belong to each amino acid.
    """
    if seed is not None:
        random.seed(seed)

    # Check that the partition sums to total codon count
    if sum(partition) != len(codons):
        raise ValueError(
            f"Partition sum {sum(partition)} != number of codons {len(codons)}"
        )

    # Check that there are enough amino-acid labels
    if len(aa_labels) < len(partition):
        raise ValueError("Not enough amino-acid labels for the partition.")

    # Copy codon list so we don’t modify the original
    codons_order = codons.copy()

    # Optionally randomize codon order before assigning
    if mode == 'random':
        random.shuffle(codons_order)

    # Build the mapping f: C -> A
    mapping = {}
    idx = 0
    for count, label in zip(partition, aa_labels):
        for _ in range(count):
            mapping[codons_order[idx]] = label
            idx += 1

    return mapping


def build_hamming_graph(codons):
    """
    Construct an undirected graph where:
    - Nodes are codons.
    - Edges connect codons with Hamming distance = 1.
    """
    G = nx.Graph()
    G.add_nodes_from(codons)

    # Compare each pair once (upper triangle)
    for i, c1 in enumerate(codons):
        for c2 in codons[i + 1:]:
            if hamming_distance(c1, c2) == 1:
                G.add_edge(c1, c2)
    return G


def plot_codon_graph(G, mapping, show_labels=True):
    """
    Draw the codon graph using NetworkX + Matplotlib.
    Each amino acid is shown as a color.
    """
    # Get unique amino acids
    amino_acids = sorted(set(mapping.values()))

    # Use a built-in Matplotlib colormap
    colors = plt.cm.tab20(range(len(amino_acids)))
    color_map = {aa: colors[i % len(colors)] for i, aa in enumerate(amino_acids)}

    # Map node colors
    node_colors = [color_map[mapping[c]] for c in G.nodes()]

    pos = nx.spring_layout(G, seed=0)

    plt.figure(figsize=(9, 9))
    # Draw faint edges for Hamming-1 connections
    nx.draw_networkx_edges(G, pos, alpha=0.25, width=0.8)
    # Draw nodes with black borders and colored fills
    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=node_colors,
        node_size=260,
        edgecolors='black',
        linewidths=0.6,
    )
    # Add codon labels if enabled
    if show_labels:
        nx.draw_networkx_labels(G, pos, font_size=7)

    # Draw color legend (amino acid -> color)
    for aa, color in color_map.items():
        plt.scatter([], [], color=color, label=aa)
    plt.legend(title="Amino acids", loc='upper left', frameon=False, fontsize=8)

    plt.axis('off')
    plt.title("Codon -> Amino Acid Assignment (Hamming-1 Graph)")
    plt.tight_layout()
    plt.show()


# ===============================================================
# -------------------------- MAIN -------------------------------
# ===============================================================

def main():
    """
    Main execution function that orchestrates the codon graph generation,
    assignment, and visualization pipeline.
    
    The function:
    1. Generates alphabet and codons based on configuration
    2. Creates a partition of codons to amino acids
    3. Assigns codons to amino acids
    4. Builds a Hamming graph
    5. Visualizes the colored graph
    6. Computes and reports robustness statistics
    
    When RUN_MULTIPLE is True, performs multiple runs with different random
    seeds and reports average statistics.
    """
    #Set true to run a bunch of times
    RUN_MULTIPLE = False

    total_rob_to_avg = 0
    total_edges_to_avg = 0

    for run in range(1 if not RUN_MULTIPLE else 1000):
        # Step 1: Generate alphabet and codons
        alphabet = generate_alphabet(ALPHABET_SIZE)
        codons = generate_codons(alphabet, CODON_LENGTH)
        total = len(codons)  # |N|^L

        # Step 2: Create or verify partition
        part = PARTITION if PARTITION is not None else auto_partition(total, NUM_AMINO_ACIDS)
        labels = make_labels(len(part))

        # print(f"Alphabet ({ALPHABET_SIZE}): {alphabet}")
        # print(f"Codon length L = {CODON_LENGTH} -> total codons = {total}")
        # print(f"Partition = {part} (sum={sum(part)})")
        # print(f"Labels    = {labels}\n")

        # Step 3: Create assignment f: C -> A
        mapping = assign_codons_to_amino_acids(
            codons, part, labels, mode=ASSIGNMENT_MODE, seed=run
        )

        # Step 4: Build Hamming graph
        G = build_hamming_graph(codons)

        # Step 5: Compute robustness statistic (fraction of intra-color edges)
        same, total_edges = 0, G.number_of_edges()
        for u, v in G.edges():
            if mapping[u] == mapping[v]:
                same += 1
        # if total_edges > 0:
            # print(f"Intra-amino-acid edges: {same}/{total_edges} "
                # f"({same / total_edges})\n")

        
        # Step 6: Visualize graph
        plot_codon_graph(G, mapping, show_labels=SHOW_NODE_LABELS)
        total_rob_to_avg += (same / total_edges) if total_edges > 0 else 0
        total_edges_to_avg += same

    # run a bunch of times

    print(f"Average robustness over {1000 if RUN_MULTIPLE else 1} runs: {total_rob_to_avg / (1000 if RUN_MULTIPLE else 1)}")
    print(f"Average intra-amino-acid edges over {1000 if RUN_MULTIPLE else 1} runs: {total_edges_to_avg / (1000 if RUN_MULTIPLE else 1)}")
# Run main if file is executed directly
if __name__ == "__main__":
    main()
