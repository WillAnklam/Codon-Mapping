"""Study the combinatorial robustness of the standard genetic code.

Codons are vertices of H(3,4), and two codons are joined when they differ in
one nucleotide.  The program counts edges whose endpoints encode the same
output and explores how that count changes when two codon assignments are
swapped.  It also checks the explicit zero-edge construction from the paper.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb
from typing import Iterable, Mapping


NUCLEOTIDES = ("A", "C", "G", "T")
CODON_LENGTH = 3

BIOLOGICAL_CLASSES = {
    "Phe": ("TTT", "TTC"),
    "Leu": ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
    "Ile": ("ATT", "ATC", "ATA"),
    "Met": ("ATG",),
    "Val": ("GTT", "GTC", "GTA", "GTG"),
    "Ser": ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
    "Pro": ("CCT", "CCC", "CCA", "CCG"),
    "Thr": ("ACT", "ACC", "ACA", "ACG"),
    "Ala": ("GCT", "GCC", "GCA", "GCG"),
    "Tyr": ("TAT", "TAC"),
    "His": ("CAT", "CAC"),
    "Gln": ("CAA", "CAG"),
    "Asn": ("AAT", "AAC"),
    "Lys": ("AAA", "AAG"),
    "Asp": ("GAT", "GAC"),
    "Glu": ("GAA", "GAG"),
    "Cys": ("TGT", "TGC"),
    "Trp": ("TGG",),
    "Arg": ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
    "Gly": ("GGT", "GGC", "GGA", "GGG"),
    "Stop": ("TAA", "TAG", "TGA"),
}

# The explicit zero-silent-edge assignment in the paper's appendix.
APPENDIX_ZERO_EDGE_CLASSES = (
    ("TTC", "TCT", "CTT", "CCA", "CAC", "ACC"),
    ("AAT", "ATA", "TAA", "GCA", "GAC", "AGC"),
    ("TTA", "TAT", "ATT", "CAA", "ACA", "AAC"),
    ("GAA", "AGA", "AAG", "ATC"),
    ("ACT", "CTA", "CAT", "TCA"),
    ("GGA", "GAG", "AGG", "TCG"),
    ("ACG", "CCT", "CTC", "TCC"),
    ("GAT", "GTA", "ATG", "AGT"),
    ("TGA", "TAG", "GGC"),
    ("TTG", "TGT", "GTT"),
    ("TGC", "GTC"),
    ("AAA", "GGG"),
    ("TAC", "GGT"),
    ("GTG", "TGG"),
    ("GCG", "CGG"),
    ("CAG", "CGA"),
    ("GCT", "CGT"),
    ("CGC", "CCG"),
    ("CTG", "GCC"),
    ("CCC",),
    ("TTT",),
)

EXPECTED_PROFILE = tuple(sorted((6,) * 3 + (4,) * 5 + (3,) * 2 + (2,) * 9 + (1,) * 2))
IMPROVING_SWAPS = (("AGC", "TGT"), ("AGT", "TGC"))


def all_codons() -> tuple[str, ...]:
    """Return the 64 vertices of H(3,4)."""
    return tuple("".join(symbols) for symbols in product(NUCLEOTIDES, repeat=CODON_LENGTH))


def hamming_distance(left: str, right: str) -> int:
    if len(left) != len(right):
        raise ValueError("Hamming distance requires strings of equal length")
    return sum(a != b for a, b in zip(left, right))


def hamming_edges(codons: Iterable[str]) -> tuple[tuple[str, str], ...]:
    """Construct every unordered Hamming-distance-one edge exactly once."""
    return tuple((u, v) for u, v in combinations(sorted(codons), 2) if hamming_distance(u, v) == 1)


def assignment_from_classes(classes: Mapping[str, Iterable[str]]) -> dict[str, str]:
    assignment: dict[str, str] = {}
    for label, codons in classes.items():
        for codon in codons:
            if codon in assignment:
                raise ValueError(f"Codon {codon} appears in more than one class")
            assignment[codon] = label
    return assignment


def appendix_zero_edge_assignment() -> dict[str, str]:
    classes = {f"A_{index}": codons for index, codons in enumerate(APPENDIX_ZERO_EDGE_CLASSES, 1)}
    return assignment_from_classes(classes)


def degeneracy_profile(assignment: Mapping[str, str]) -> tuple[int, ...]:
    return tuple(sorted(Counter(assignment.values()).values()))


def silent_edge_count(assignment: Mapping[str, str], edges: Iterable[tuple[str, str]]) -> int:
    return sum(assignment[u] == assignment[v] for u, v in edges)


def swap(assignment: Mapping[str, str], u: str, v: str) -> dict[str, str]:
    swapped = dict(assignment)
    swapped[u], swapped[v] = swapped[v], swapped[u]
    return swapped


def enumerate_nontrivial_swaps(
    assignment: Mapping[str, str], edges: tuple[tuple[str, str], ...]
) -> dict[str, object]:
    """Exhaustively evaluate swaps of codons having different outputs."""
    baseline = silent_edge_count(assignment, edges)
    improving: list[dict[str, object]] = []
    increasing = decreasing = unchanged = 0
    best_count = -1

    for u, v in combinations(sorted(assignment), 2):
        if assignment[u] == assignment[v]:
            continue

        candidate_count = silent_edge_count(swap(assignment, u, v), edges)
        delta = candidate_count - baseline
        if delta > 0:
            increasing += 1
            improving.append(
                {
                    "codons": (u, v),
                    "labels": (assignment[u], assignment[v]),
                    "delta": delta,
                    "silent_edges": candidate_count,
                }
            )
        elif delta < 0:
            decreasing += 1
        else:
            unchanged += 1
        best_count = max(best_count, candidate_count)

    return {
        "considered": increasing + decreasing + unchanged,
        "increasing": increasing,
        "decreasing": decreasing,
        "unchanged": unchanged,
        "best_silent_edges": best_count,
        "improving_swaps": improving,
    }


def analyze_standard_code() -> dict[str, object]:
    """Compute the graph, robustness, and complete codon-swap landscape."""
    codons = all_codons()
    edges = hamming_edges(codons)
    biological = assignment_from_classes(BIOLOGICAL_CLASSES)

    assert len(codons) == 64
    assert len(edges) == 288
    degrees = Counter(endpoint for edge in edges for endpoint in edge)
    assert set(degrees.values()) == {9}
    assert set(biological) == set(codons)
    assert degeneracy_profile(biological) == EXPECTED_PROFILE

    same_class_pairs = sum(comb(size, 2) for size in EXPECTED_PROFILE)
    random_expectation = Fraction(same_class_pairs, comb(len(codons), 2))
    biological_silent = silent_edge_count(biological, edges)
    biological_neighborhood = enumerate_nontrivial_swaps(biological, edges)

    assert same_class_pairs == 90
    assert biological_silent == 69
    assert biological_neighborhood["considered"] == 1926
    found_swaps = tuple(item["codons"] for item in biological_neighborhood["improving_swaps"])
    assert found_swaps == IMPROVING_SWAPS
    assert all(item["silent_edges"] == 70 for item in biological_neighborhood["improving_swaps"])

    improved = dict(biological)
    for u, v in IMPROVING_SWAPS:
        improved = swap(improved, u, v)
    improved_silent = silent_edge_count(improved, edges)
    improved_neighborhood = enumerate_nontrivial_swaps(improved, edges)

    assert degeneracy_profile(improved) == EXPECTED_PROFILE
    assert improved_silent == 71
    assert improved_neighborhood["considered"] == 1926
    assert improved_neighborhood["increasing"] == 0

    zero_edge = appendix_zero_edge_assignment()
    assert set(zero_edge) == set(codons)
    assert degeneracy_profile(zero_edge) == EXPECTED_PROFILE
    zero_edge_silent = silent_edge_count(zero_edge, edges)
    assert zero_edge_silent == 0

    fixed_profile_upper_bound = 3 * 9 + 5 * 6 + 2 * 3 + 9 * 1 + 2 * 0
    assert fixed_profile_upper_bound == 72

    return {
        "vertices": len(codons),
        "edges": len(edges),
        "degree": 9,
        "degeneracy_profile": list(reversed(EXPECTED_PROFILE)),
        "within_class_codon_pairs": same_class_pairs,
        "random_expected_robustness": {
            "fraction": f"{random_expectation.numerator}/{random_expectation.denominator}",
            "decimal": float(random_expectation),
        },
        "biological": {
            "silent_edges": biological_silent,
            "robustness": biological_silent / len(edges),
            "nontrivial_swaps": biological_neighborhood["considered"],
            "improving_swaps": biological_neighborhood["improving_swaps"],
        },
        "improved_assignment": {
            "silent_edges": improved_silent,
            "robustness": improved_silent / len(edges),
            "further_improving_swaps": improved_neighborhood["increasing"],
        },
        "appendix_zero_edge_assignment": {"silent_edges": zero_edge_silent},
        "fixed_profile_upper_bound": fixed_profile_upper_bound,
    }


def print_report(results: Mapping[str, object]) -> None:
    biological = results["biological"]
    improved = results["improved_assignment"]
    expectation = results["random_expected_robustness"]
    print("Combinatorial robustness of the standard genetic code")
    print(f"H(3,4): {results['vertices']} vertices, {results['edges']} edges, degree {results['degree']}")
    print(
        "Random fixed-profile expectation: "
        f"{expectation['fraction']} = {expectation['decimal']:.10f}"
    )
    print(
        "Biological code: "
        f"{biological['silent_edges']}/288 = {biological['robustness']:.10f}"
    )
    print(f"Nontrivial swaps checked: {biological['nontrivial_swaps']}")
    for item in biological["improving_swaps"]:
        u, v = item["codons"]
        left, right = item["labels"]
        print(f"  improving swap: {u} ({left}) <-> {v} ({right}); silent edges = {item['silent_edges']}")
    print(
        "Both swaps applied: "
        f"{improved['silent_edges']}/288 = {improved['robustness']:.10f}; "
        f"further improving swaps = {improved['further_improving_swaps']}"
    )
    print(f"Appendix construction: {results['appendix_zero_edge_assignment']['silent_edges']} silent edges")
    print(f"Fixed-profile upper bound from classwise bounds: {results['fixed_profile_upper_bound']} silent edges")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    args = parser.parse_args()
    results = analyze_standard_code()
    if args.json:
        print(json.dumps(results, indent=2))
    else:
        print_report(results)


if __name__ == "__main__":
    main()
