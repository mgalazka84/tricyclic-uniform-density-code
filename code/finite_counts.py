#!/usr/bin/env python3
"""Reproduce the finite proof table, simple examples and U(L) term list.

Run from the repository root:
    python3 code/finite_counts.py --output finite_counts.json

All sums are integer counts. The identity element in Burnside's average
uses the certified labelled counts. Other elements enumerate fixed lengths.
"""
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path

from verify import (KERNELS, compositions, explicit_criterion, formula,
                    quasipolynomial, recognize_tricyclic, subdivide)
from verify_extensions import (CUT_KERNELS, cut_counts, edge_actions,
                               orbit_checks, permutation_cycles)


def simple_subdivision(edges, lengths):
    # All kernels in this call are loopless. Parallel edges can survive only
    # if both paths have length one and share their pair of endpoints.
    single = [tuple(sorted(edge)) for edge, length in zip(edges, lengths) if length == 1]
    return len(single) == len(set(single))


def finite_table():
    triples, simple_examples = {}, {}
    needed = 0
    for name, (_, edges) in KERNELS.items():
        d = len(edges)
        rows, simple = [], {}
        for q in range(1, 7):
            row = []
            for residue in range(3):
                total = 3 * q + residue
                admissible = simple_count = 0
                for lengths in compositions(total, d):
                    if explicit_criterion(name, lengths):
                        admissible += 1
                        simple_count += simple_subdivision(edges, lengths)
                assert admissible == formula(name, total) == quasipolynomial(name, total)
                if q <= d:
                    needed += 1
                expected_simple = admissible if total >= 7 else int(name == "D" and total == 6)
                assert simple_count == expected_simple
                simple[total] = simple_count
                row.append(admissible)
            rows.append(row)
        triples[name] = rows
        simple_examples[name] = simple
    assert needed == 63
    return dict(determining_values=needed, triples=triples, simple=simple_examples)


def fixed_lengths_by_total(name, cycles, maximum):
    """One length per permutation cycle, counting every positive assignment once."""
    counts = [0] * (maximum + 1)
    sizes = [len(cycle) for cycle in cycles]
    vector = [0] * sum(sizes)
    suffix = [sum(sizes[i:]) for i in range(len(cycles) + 1)]
    def visit(i, used):
        if i == len(cycles):
            if explicit_criterion(name, vector):
                counts[used] += 1
            return
        # Every admissible coordinate is at most floor(L/3).
        upper = min(maximum // 3, (maximum - used - suffix[i + 1]) // sizes[i])
        for value in range(1, upper + 1):
            for edge in cycles[i]:
                vector[edge] = value
            visit(i + 1, used + value * sizes[i])
    visit(0, 0)
    return counts


def aggregate_sequence(maximum):
    by_kernel = {}
    for name, (vertices, edges) in KERNELS.items():
        actions = edge_actions(vertices, edges)
        fixed_sums = [quasipolynomial(name, total) for total in range(maximum + 1)]
        partitions = Counter()
        for permutation in actions:
            cycles = permutation_cycles(permutation)
            if len(cycles) < len(edges):
                partitions[tuple(sorted(tuple(sorted(c)) for c in cycles))] += 1
        for cycles, multiplicity in partitions.items():
            fixed = fixed_lengths_by_total(name, cycles, maximum)
            for total, count in enumerate(fixed):
                fixed_sums[total] += multiplicity * count
        assert all(n % len(actions) == 0 for n in fixed_sums)
        by_kernel[name] = [n // len(actions) for n in fixed_sums]
    for name in CUT_KERNELS:
        by_kernel[name] = [cut_counts(name, total)[1] for total in range(maximum + 1)]
    aggregate = [sum(values[total] for values in by_kernel.values())
                 for total in range(1, maximum + 1)]
    # Check the Burnside implementation against explicit canonical representatives.
    checked_to = min(20, maximum)
    direct = orbit_checks(checked_to)
    for name, item in direct.items():
        actual = [r['unlabelled'] for r in item['counts']]
        assert by_kernel[name][1:checked_to + 1] == actual
    if maximum >= 12:
        assert aggregate[11] == 75
        assert [by_kernel[k][12] for k in 'ABCDFHJM'] == [4, 20, 17, 20, 1, 2, 3, 8]
    return dict(first_total=1, last_total=maximum, checked_by_canonical_orbits_to=checked_to,
                counts=aggregate)


def density_example():
    n, edges = KERNELS['A']
    size, graph = subdivide(n, edges, [3, 3, 3, 7])
    result = recognize_tricyclic(size, graph)
    assert not result['uniform'] and tuple(result['density']) == (9, 7)
    assert result['girth'] == 6 and result['largest_series_class'] == 7
    assert set(result['densest_set']) == set(range(9))
    return dict(lengths=[3, 3, 3, 7], whole_density=[16, 13],
                cycle_density=[6, 5], maximum_density=[9, 7], certificate_edges=9)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--max-total', type=int, default=60)
    parser.add_argument('--output', default='finite_counts.json')
    args = parser.parse_args()
    if args.max_total < 1:
        parser.error('--max-total must be positive')
    report = {'finite_table': finite_table(), 'aggregate_U': aggregate_sequence(args.max_total),
              'density_example': density_example()}
    report['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(args.output).write_text(json.dumps(report, indent=2) + '\n')
    print('Finite counts, simple subdivisions, orbit totals and witness example verified.')


if __name__ == '__main__':
    main()
