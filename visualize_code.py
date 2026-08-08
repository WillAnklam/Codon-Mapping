"""Generate a visualization of the standard biological codon assignment."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import networkx as nx

from codon_robustness import BIOLOGICAL_CLASSES, assignment_from_classes, hamming_edges


def draw_biological_code(output: Path, show_labels: bool = False, show: bool = False) -> None:
    assignment = assignment_from_classes(BIOLOGICAL_CLASSES)
    graph = nx.Graph()
    graph.add_nodes_from(sorted(assignment))
    graph.add_edges_from(hamming_edges(assignment))

    labels = sorted(set(assignment.values()))
    palette = plt.colormaps["tab20"].resampled(len(labels))
    colors = {label: palette(index) for index, label in enumerate(labels)}
    node_colors = [colors[assignment[codon]] for codon in graph.nodes]
    silent = [(u, v) for u, v in graph.edges if assignment[u] == assignment[v]]
    nonsilent = [(u, v) for u, v in graph.edges if assignment[u] != assignment[v]]
    positions = nx.spring_layout(graph, seed=7, k=0.55, iterations=300)

    figure, axis = plt.subplots(figsize=(14, 14))
    nx.draw_networkx_edges(graph, positions, edgelist=nonsilent, edge_color="lightgray", alpha=0.25, width=0.8, ax=axis)
    nx.draw_networkx_edges(graph, positions, edgelist=silent, edge_color="black", alpha=0.75, width=2.4, ax=axis)
    nx.draw_networkx_nodes(
        graph,
        positions,
        node_color=node_colors,
        node_size=430,
        edgecolors="black",
        linewidths=0.7,
        ax=axis,
    )
    if show_labels:
        nx.draw_networkx_labels(graph, positions, font_size=7, font_weight="bold", ax=axis)

    handles = [
        plt.Line2D([], [], marker="o", linestyle="", color=colors[label], label=label)
        for label in labels
    ]
    axis.legend(handles=handles, title="Output", loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False)
    axis.set_title("Standard genetic code in H(3,4)\nBlack edges are silent single-nucleotide substitutions")
    axis.set_axis_off()
    figure.tight_layout()
    figure.savefig(output, dpi=300, bbox_inches="tight")
    print(f"Saved {output} ({len(silent)}/{graph.number_of_edges()} silent edges)")
    if show:
        plt.show()
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("biological_codon_graph.png"))
    parser.add_argument("--labels", action="store_true", help="draw codon labels")
    parser.add_argument("--show", action="store_true", help="also open an interactive plot window")
    args = parser.parse_args()
    draw_biological_code(args.output, show_labels=args.labels, show=args.show)


if __name__ == "__main__":
    main()
