"""Computational analysis of Hamming robustness in the standard genetic code.

The program constructs the codon graph, counts silent mutation edges, samples
random assignments with the biological degeneracy profile, and exhaustively
searches the degeneracy-preserving codon-swap neighborhood.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction
from itertools import combinations, product
from math import comb
from random import Random
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

EXPECTED_PROFILE = tuple(sorted((6,) * 3 + (4,) * 5 + (3,) * 2 + (2,) * 9 + (1,) * 2))
EXPECTED_IMPROVING_SWAPS = (("AGC", "TGT"), ("AGT", "TGC"))


def all_codons() -> tuple[str, ...]:
    """Construct the 64 vertices of H(3,4)."""
    return tuple("".join(symbols) for symbols in product(NUCLEOTIDES, repeat=CODON_LENGTH))


def hamming_distance(left: str, right: str) -> int:
    if len(left) != len(right):
        raise ValueError("Hamming distance requires strings of equal length")
    return sum(a != b for a, b in zip(left, right))


def hamming_edges(codons: Iterable[str]) -> tuple[tuple[str, str], ...]:
    """List all unordered codon pairs at Hamming distance one."""
    return tuple(
        (left, right)
        for left, right in combinations(sorted(codons), 2)
        if hamming_distance(left, right) == 1
    )


def assignment_from_classes(classes: Mapping[str, Iterable[str]]) -> dict[str, str]:
    assignment: dict[str, str] = {}
    for label, codons in classes.items():
        for codon in codons:
            if codon in assignment:
                raise ValueError(f"Codon {codon} appears in more than one class")
            assignment[codon] = label
    return assignment


def degeneracy_profile(assignment: Mapping[str, str]) -> tuple[int, ...]:
    return tuple(sorted(Counter(assignment.values()).values()))


def silent_edge_count(
    assignment: Mapping[str, str], edges: Iterable[tuple[str, str]]
) -> int:
    """Count edges whose endpoints have the same assigned output."""
    return sum(assignment[left] == assignment[right] for left, right in edges)


def swap_assignment(
    assignment: Mapping[str, str], left: str, right: str
) -> dict[str, str]:
    swapped = dict(assignment)
    swapped[left], swapped[right] = swapped[right], swapped[left]
    return swapped


def random_fixed_profile_assignment(rng: Random) -> dict[str, str]:
    """Draw uniformly from assignments with the biological profile.

    The profile sizes are shuffled among the 21 output labels, the codons are
    independently shuffled, and consecutive blocks are assigned to the labels.
    This gives every labeled assignment with the required profile the same
    probability.
    """
    labels = sorted(BIOLOGICAL_CLASSES)
    sizes = list(EXPECTED_PROFILE)
    codons = list(all_codons())
    rng.shuffle(sizes)
    rng.shuffle(codons)

    assignment: dict[str, str] = {}
    start = 0
    for label, size in zip(labels, sizes):
        for codon in codons[start : start + size]:
            assignment[codon] = label
        start += size

    assert start == len(codons)
    assert degeneracy_profile(assignment) == EXPECTED_PROFILE
    return assignment


def sample_random_reassignments(
    trials: int, seed: int, edges: tuple[tuple[str, str], ...]
) -> dict[str, object]:
    """Estimate robustness by repeated uniform fixed-profile reassignment."""
    if trials <= 0:
        raise ValueError("trials must be positive")

    rng = Random(seed)
    silent_counts = [
        silent_edge_count(random_fixed_profile_assignment(rng), edges)
        for _ in range(trials)
    ]
    return {
        "trials": trials,
        "seed": seed,
        "mean_silent_edges": sum(silent_counts) / trials,
        "mean_robustness": sum(silent_counts) / (trials * len(edges)),
        "minimum_silent_edges": min(silent_counts),
        "maximum_silent_edges": max(silent_counts),
    }


def enumerate_nontrivial_swaps(
    assignment: Mapping[str, str], edges: tuple[tuple[str, str], ...]
) -> dict[str, object]:
    """Construct and recount every swap of codons with different outputs."""
    baseline = silent_edge_count(assignment, edges)
    improving: list[dict[str, object]] = []
    increasing = decreasing = unchanged = skipped_synonymous = 0

    for left, right in combinations(sorted(assignment), 2):
        if assignment[left] == assignment[right]:
            skipped_synonymous += 1
            continue

        candidate = swap_assignment(assignment, left, right)
        candidate_count = silent_edge_count(candidate, edges)
        delta = candidate_count - baseline
        if delta > 0:
            increasing += 1
            improving.append(
                {
                    "codons": (left, right),
                    "labels": (assignment[left], assignment[right]),
                    "delta": delta,
                    "silent_edges": candidate_count,
                }
            )
        elif delta < 0:
            decreasing += 1
        else:
            unchanged += 1

    return {
        "all_pairs": skipped_synonymous + increasing + decreasing + unchanged,
        "skipped_synonymous": skipped_synonymous,
        "considered": increasing + decreasing + unchanged,
        "increasing": increasing,
        "decreasing": decreasing,
        "unchanged": unchanged,
        "improving_swaps": improving,
    }


def analyze_standard_code(random_trials: int = 0, seed: int = 2026) -> dict[str, object]:
    """Run the graph, robustness, random-reassignment, and swap calculations."""
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
    all_codon_pairs = comb(len(codons), 2)
    random_expectation = Fraction(same_class_pairs, all_codon_pairs)
    biological_silent = silent_edge_count(biological, edges)
    biological_neighborhood = enumerate_nontrivial_swaps(biological, edges)

    assert all_codon_pairs == 2016
    assert same_class_pairs == 90
    assert biological_silent == 69
    assert biological_neighborhood["all_pairs"] == all_codon_pairs
    assert biological_neighborhood["skipped_synonymous"] == same_class_pairs
    assert biological_neighborhood["considered"] == 1926

    found_swaps = tuple(
        item["codons"] for item in biological_neighborhood["improving_swaps"]
    )
    assert found_swaps == EXPECTED_IMPROVING_SWAPS
    assert all(
        item["silent_edges"] == 70
        for item in biological_neighborhood["improving_swaps"]
    )

    improved = dict(biological)
    for left, right in found_swaps:
        improved = swap_assignment(improved, left, right)
    improved_silent = silent_edge_count(improved, edges)
    improved_neighborhood = enumerate_nontrivial_swaps(improved, edges)

    assert degeneracy_profile(improved) == EXPECTED_PROFILE
    assert improved_silent == 71
    assert improved_neighborhood["considered"] == 1926
    assert improved_neighborhood["increasing"] == 0

    results: dict[str, object] = {
        "vertices": len(codons),
        "edges": len(edges),
        "degree": 9,
        "degeneracy_profile": list(reversed(EXPECTED_PROFILE)),
        "random_expected_robustness": {
            "within_class_pairs": same_class_pairs,
            "all_codon_pairs": all_codon_pairs,
            "fraction": f"{random_expectation.numerator}/{random_expectation.denominator}",
            "decimal": float(random_expectation),
        },
        "biological": {
            "silent_edges": biological_silent,
            "robustness": biological_silent / len(edges),
            "all_codon_pairs": biological_neighborhood["all_pairs"],
            "synonymous_pairs_skipped": biological_neighborhood["skipped_synonymous"],
            "nontrivial_swaps": biological_neighborhood["considered"],
            "improving_swaps": biological_neighborhood["improving_swaps"],
        },
        "improved_assignment": {
            "silent_edges": improved_silent,
            "robustness": improved_silent / len(edges),
            "further_improving_swaps": improved_neighborhood["increasing"],
        },
    }
    if random_trials:
        results["random_reassignment_sample"] = sample_random_reassignments(
            random_trials, seed, edges
        )
    return results


def print_report(results: Mapping[str, object]) -> None:
    expectation = results["random_expected_robustness"]
    biological = results["biological"]
    improved = results["improved_assignment"]

    print("Combinatorial robustness of the standard genetic code")
    print(
        f"H(3,4): {results['vertices']} vertices, "
        f"{results['edges']} edges, degree {results['degree']}"
    )
    print(
        "Random fixed-profile expectation: "
        f"{expectation['within_class_pairs']}/{expectation['all_codon_pairs']} = "
        f"{expectation['fraction']} = {expectation['decimal']:.10f}"
    )
    print(
        "Biological code: "
        f"{biological['silent_edges']}/{results['edges']} = "
        f"{biological['robustness']:.10f}"
    )
    print(
        "Codon pairs: "
        f"{biological['all_codon_pairs']} total - "
        f"{biological['synonymous_pairs_skipped']} synonymous = "
        f"{biological['nontrivial_swaps']} nontrivial swaps checked"
    )
    for item in biological["improving_swaps"]:
        left, right = item["codons"]
        left_label, right_label = item["labels"]
        print(
            f"  improving swap: {left} ({left_label}) <-> "
            f"{right} ({right_label}); silent edges = {item['silent_edges']}"
        )
    print(
        "Both swaps applied: "
        f"{improved['silent_edges']}/{results['edges']} = "
        f"{improved['robustness']:.10f}; "
        f"further improving swaps = {improved['further_improving_swaps']}"
    )

    if "random_reassignment_sample" in results:
        sample = results["random_reassignment_sample"]
        print(
            f"Random reassignment sample ({sample['trials']} trials, "
            f"seed {sample['seed']}): mean robustness = "
            f"{sample['mean_robustness']:.10f}"
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", action="store_true", help="emit machine-readable JSON")
    parser.add_argument(
        "--random-trials",
        type=int,
        default=0,
        metavar="N",
        help="also sample N uniform random fixed-profile reassignments",
    )
    parser.add_argument("--seed", type=int, default=2026, help="seed for random reassignment")
    args = parser.parse_args()
    if args.random_trials < 0:
        parser.error("--random-trials must be nonnegative")

    results = analyze_standard_code(random_trials=args.random_trials, seed=args.seed)
    if args.json:
        print(json.dumps(results, indent=2))
    else:
        print_report(results)


if __name__ == "__main__":
    main()
