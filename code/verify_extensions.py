#!/usr/bin/env python3
"""Exact checks for cut vertices, kernel symmetries and convergence rates.

Run from the repository root:
    python3 code/verify_extensions.py --output extension_results.json
Uses only the Python standard library and the companion verification module.
"""
import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from math import comb
from pathlib import Path

from verify import KERNELS, QUASI, compositions, direct_density, rank, subdivide


# A loop belongs to the suppressed kernel only. Its replacement has length >= 2.
CUT_KERNELS = {
    'F': (1, [(0, 0), (0, 0), (0, 0)]),
    'H': (2, [(0, 0), (1, 1), (0, 1), (0, 1)]),
    'J': (2, [(0, 0), (0, 1), (0, 1), (0, 1)]),
    'M': (3, [(2, 2), (0, 2), (2, 1), (0, 1), (0, 1)]),
}


def loopless_lengths(edges, lengths):
    return all(u != v or length >= 2 for (u, v), length in zip(edges, lengths))


def cut_counts(name, length):
    """Pair of labelled and unlabelled counts in the four cut kernels."""
    t, residue = divmod(length, 3)
    if residue or t < 2:
        return 0, 0
    if name == 'F':
        return 1, 1
    if name == 'H':
        return t - 1, t // 2
    triples = (t - 1) * (t + 4) // 2
    if name == 'J':
        partitions = sum((t - 3 * k) // 2 + 1 for k in range(t // 3 + 1))
        return triples, partitions - 1
    assert name == 'M'
    labelled = triples * (2 * t - 3) // 3
    orbits = sum((j // 2) * (j // 2 + 1) for j in range(2, t + 1)) - t // 2
    return labelled, orbits


def kernel_constraints(vertices, edges):
    """Cycle counts on every edge subset, including loops and parallel copies."""
    result = []
    for size in range(1, len(edges) + 1):
        for indices in itertools.combinations(range(len(edges)), size):
            nullity = size - rank(vertices, [edges[i] for i in indices])
            if nullity:
                result.append((indices, nullity))
    return result


def admissible_by_subsets(lengths, constraints):
    total = sum(lengths)
    return all(3 * sum(lengths[i] for i in subset) >= total * nullity
               for subset, nullity in constraints)


def edge_actions(vertices, edges):
    """Distinct edge permutations induced by vertex maps and parallel-copy permutations."""
    endpoints = [tuple(sorted(edge)) for edge in edges]
    actions = set()
    permutations = list(itertools.permutations(range(len(edges))))
    for vertex_map in itertools.permutations(range(vertices)):
        targets = [tuple(sorted((vertex_map[u], vertex_map[v]))) for u, v in edges]
        for perm in permutations:
            if all(targets[i] == endpoints[perm[i]] for i in range(len(edges))):
                actions.add(perm)
    return sorted(actions)


def permutation_cycles(perm):
    unseen = set(range(len(perm)))
    cycles = []
    while unseen:
        i = min(unseen)
        cycle = []
        while i in unseen:
            cycle.append(i)
            unseen.remove(i)
            i = perm[i]
        cycles.append(cycle)
    return cycles


def fixed_count(total, cycles, edge_count, predicate):
    """Count fixed lengths directly on the disjoint cycles of one permutation."""
    vector = [0] * edge_count
    sizes = [len(cycle) for cycle in cycles]
    def count(index, budget):
        if index == len(cycles):
            return int(budget == 0 and predicate(vector))
        cap = (budget - sum(sizes[index + 1:])) // sizes[index]
        result = 0
        for value in range(1, cap + 1):
            for edge in cycles[index]:
                vector[edge] = value
            result += count(index + 1, budget - sizes[index] * value)
        return result
    return count(0, total)


def orbit_checks(max_total):
    """Compare Burnside's sum with canonical orbit representatives for all eight kernels."""
    result = {}
    expected_orders = dict(A=24, B=8, C=16, D=24, F=6, H=4, J=6, M=4)
    for name, (vertices, edges) in (KERNELS | CUT_KERNELS).items():
        actions = edge_actions(vertices, edges)
        assert len(actions) == expected_orders[name]
        cycles = [permutation_cycles(perm) for perm in actions]
        constraints = kernel_constraints(vertices, edges)
        def admissible(vector):
            return loopless_lengths(edges, vector) and admissible_by_subsets(vector, constraints)
        rows = []
        for total in range(1, max_total + 1):
            representatives = set()
            labelled = 0
            for length in compositions(total, len(edges)):
                if not admissible(length):
                    continue
                labelled += 1
                representatives.add(min(tuple(length[i] for i in perm) for perm in actions))
            fixed = sum(fixed_count(total, c, len(edges), admissible) for c in cycles)
            assert fixed == len(actions) * len(representatives)
            if name in CUT_KERNELS:
                assert (labelled, len(representatives)) == cut_counts(name, total)
            rows.append(dict(length=total, labelled=labelled, unlabelled=len(representatives)))
        result[name] = dict(group_order=len(actions), counts=rows)
    return result


def cut_checks(max_total, subset_total):
    """Test closed forms by kernel-subset inequalities and exhaustive subsets of small graphs."""
    result = {}
    for name, (vertices, edges) in CUT_KERNELS.items():
        constraints = kernel_constraints(vertices, edges)
        checked = small = subsets = 0
        for total in range(1, max_total + 1):
            actual = 0
            for lengths in compositions(total, len(edges)):
                if not loopless_lengths(edges, lengths):
                    continue
                expected = admissible_by_subsets(lengths, constraints)
                if total <= subset_total:
                    n, graph = subdivide(vertices, edges, lengths)
                    assert expected == direct_density(n, graph), (name, lengths)
                    small += 1
                    subsets += 2**total - 1
                actual += expected
                checked += 1
            assert actual == cut_counts(name, total)[0], (name, total)
        result[name] = dict(length_vectors=checked, expanded_graphs=small, edge_subsets=subsets)
    return result


def small_graph_orbits(orbit_report):
    """Enumerate small multigraphs directly, independently of the list of kernels."""
    report = []
    for vertices in range(2, 6):
        total = vertices + 2
        pairs = list(itertools.combinations(range(vertices), 2))
        maps = list(itertools.permutations(range(vertices)))
        classes = set()
        uniform = connected = 0
        for graph in itertools.combinations_with_replacement(pairs, total):
            if rank(vertices, graph) != vertices - 1:
                continue
            connected += 1
            if not direct_density(vertices, graph):
                continue
            uniform += 1
            key = min(tuple(sorted(tuple(sorted((p[u], p[v]))) for u, v in graph)) for p in maps)
            classes.add(key)
        predicted = sum(item['counts'][total - 1]['unlabelled'] for item in orbit_report.values())
        assert len(classes) == predicted, (vertices, len(classes), predicted)
        report.append(dict(vertices=vertices, total_edges=total, connected=connected,
                           uniform_vertex_labelled=uniform, uniform_nonisomorphic=len(classes)))
    return report


def theta_checks(max_total):
    for total in range(3, max_total + 1):
        actual = sum(2 * max(lengths) <= total for lengths in compositions(total, 3))
        t, residue = divmod(total, 2)
        expected = t * (t + 1) // 2 if residue else (t - 1) * (t + 4) // 2
        assert actual == expected
    return dict(min_total=3, max_total=max_total)


def rate_checks():
    expected = dict(A=['24', '12', '0'], B=['128/5', '48/5', '8'],
                    C=['300/11', '80/11', '100/11'], D=['135/7', '10', '65/7'])
    for name, rows in QUASI.items():
        for residue, (_, co) in enumerate(rows):
            h = len(co) - 1
            coefficient = 3 * Fraction(co[-2], co[-1]) - h * residue + comb(h + 1, 2)
            assert str(coefficient) == expected[name][residue]
    return expected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='extension_results.json')
    args = parser.parse_args()
    results = {'cut_max_total': 24, 'expanded_max_total': 12, 'orbit_max_total': 15}
    results['cut_checks'] = cut_checks(results['cut_max_total'], results['expanded_max_total'])
    results['orbit_checks'] = orbit_checks(results['orbit_max_total'])
    results['independent_graph_orbits'] = small_graph_orbits(results['orbit_checks'])
    results['theta_checks'] = theta_checks(100)
    results['first_order_coefficients'] = rate_checks()
    results['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(args.output).write_text(json.dumps(results, indent=2) + '\n')
    print('Cut counts, expanded graphs, orbit counts, theta counts and rates passed.')


if __name__ == '__main__':
    main()
