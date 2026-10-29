#!/usr/bin/env python3
"""Exact, dependency-free checks for tricyclic uniform density and counts.

Run: python3 code/verify.py --output verification_results.json
All graph and counting comparisons use integers. Tests supplement proofs.
"""
from __future__ import annotations
import argparse
import hashlib
import itertools as it
import json
import math
import time
from pathlib import Path


def rank(n, edges):
    parent = list(range(n))
    def root(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v
    r = 0
    for u, v in edges:
        a, b = root(u), root(v)
        if a != b:
            parent[a] = b
            r += 1
    return r


def signatures(n, edges):
    """Coordinates of edge evaluations on a fundamental cycle basis."""
    adjacency = [[] for _ in range(n)]
    for i, (u, v) in enumerate(edges):
        adjacency[u].append((v, i))
        adjacency[v].append((u, i))
    tree = [[] for _ in range(n)]
    seen, tree_indices, stack = {0}, set(), [0]
    while stack:
        u = stack.pop()
        for v, i in adjacency[u]:
            if v not in seen:
                seen.add(v)
                stack.append(v)
                tree_indices.add(i)
                tree[u].append((v, i))
                tree[v].append((u, i))
    assert len(seen) == n, 'The input must be connected.'
    chords = [i for i in range(len(edges)) if i not in tree_indices]
    sig = [0] * len(edges)
    for j, idx in enumerate(chords):
        u, v = edges[idx]
        prev = {u: None}
        stack = [u]
        while stack and v not in prev:
            x = stack.pop()
            for y, e in tree[x]:
                if y not in prev:
                    prev[y] = (x, e)
                    stack.append(y)
        sig[idx] |= 1 << j
        x = v
        while x != u:
            x, e = prev[x]
            sig[e] |= 1 << j
    return sig, len(chords)


def recognize_tricyclic(n, edges):
    """Recognize a connected loopless tricyclic multigraph; return certificates.

    Edges are indexed by their positions, including distinct parallel copies.
    For a bridgeless input, density and densest_set give the exact maximum.
    For an input with a bridge, obstruction is a violating edge-index set.
    """
    sig, k = signatures(n, edges)
    assert k == 3
    if 0 in sig:
        bridge = sig.index(0)
        return dict(uniform=False, obstruction=[i for i in range(len(edges))
                                               if i != bridge],
                    density=None, densest_set=None, bridge=bridge)
    weights = [sig.count(i) for i in range(8)]
    largest = max(range(1, 8), key=lambda c: weights[c])
    cycle = min(range(1, 8), key=lambda z: sum(weights[c] for c in range(1, 8)
                                             if (c & z).bit_count() % 2))
    shortest_set = [i for i, c in enumerate(sig) if (c & cycle).bit_count() % 2]
    complement = [i for i, c in enumerate(sig) if c != largest]
    m, g, s = len(edges), len(shortest_set), weights[largest]
    candidates = [((m, m - 3), list(range(m))),
                  ((g, g - 1), shortest_set),
                  ((m - s, m - s - 2), complement)]
    best, witness = candidates[0]
    for value, subset in candidates[1:]:
        if value[0] * best[1] > best[0] * value[1]:
            best, witness = value, subset
    uniform = 3 * s <= m <= 3 * g
    return dict(uniform=uniform, obstruction=None if uniform else witness,
                density=best, densest_set=witness, girth=g, largest_series_class=s)


def signature_criterion(n, edges):
    return recognize_tricyclic(n, edges)['uniform']


def direct_density(n, edges, return_max=False):
    m = len(edges)
    # Iterate over all nonempty subsets, even after a violation was found.
    # This makes the reported number of comparisons exact and reproducible.
    okay = True
    maximum = (1, 1)
    for mask in range(1, 1 << m):
        chosen = [e for i, e in enumerate(edges) if mask >> i & 1]
        r = rank(n, chosen)
        if maximum[0] * r < len(chosen) * maximum[1]:
            maximum = (len(chosen), r)
        if m * r < (n - 1) * len(chosen):
            okay = False
    return (okay, maximum) if return_max else okay


def check_tricyclic_graph(n, edges):
    """Compare signatures and certificates with exhaustive density and deletions."""
    m = len(edges)
    actual, maximum = direct_density(n, edges, return_max=True)
    answer = recognize_tricyclic(n, edges)
    assert actual == answer['uniform'], (n, edges, actual, answer)
    if not actual:
        obstruction = [edges[i] for i in answer['obstruction']]
        assert m * rank(n, obstruction) < (n - 1) * len(obstruction)
    sig, _ = signatures(n, edges)
    for i in range(m):
        no_i = [e for j, e in enumerate(edges) if j != i]
        assert (sig[i] == 0) == (rank(n, no_i) < n - 1)
    if 0 not in sig:
        p, q = answer['density']
        assert p * maximum[1] == maximum[0] * q
        chosen = [edges[i] for i in answer['densest_set']]
        assert len(chosen) * maximum[1] == maximum[0] * rank(n, chosen)
        for i, j in it.combinations(range(m), 2):
            remainder = [e for h, e in enumerate(edges) if h not in (i, j)]
            assert (sig[i] == sig[j]) == (rank(n, remainder) < n - 1)
    return actual


def exhaustive_graphs():
    records = []
    comparisons = 0
    for n in range(4, 7):
        m = n + 2
        available = list(it.combinations(range(n), 2))
        connected = uniform = 0
        for edges in it.combinations(available, m):
            if rank(n, edges) != n - 1:
                continue
            connected += 1
            uniform += check_tricyclic_graph(n, edges)
            comparisons += (1 << m) - 1
        records.append(dict(vertices=n, connected_graphs=connected,
                            uniformly_dense=uniform))
    return dict(records=records, subset_comparisons=comparisons,
                total_graphs=sum(x['connected_graphs'] for x in records))


def exhaustive_multigraphs():
    """All connected loopless tricyclic edge-multiplicity vectors on 2--5 vertices."""
    records, comparisons = [], 0
    for n in range(2, 6):
        m = n + 2
        pairs = list(it.combinations(range(n), 2))
        connected = uniform = nonsimple = 0
        for edges in it.combinations_with_replacement(pairs, m):
            if rank(n, edges) != n - 1:
                continue
            connected += 1
            nonsimple += len(set(edges)) < m
            uniform += check_tricyclic_graph(n, edges)
            comparisons += (1 << m) - 1
        records.append(dict(vertices=n, connected_multigraphs=connected,
                            with_parallel_edges=nonsimple, uniformly_dense=uniform))
    return dict(records=records, subset_comparisons=comparisons,
                total_multigraphs=sum(x['connected_multigraphs'] for x in records),
                convention='labelled vertices; one input per edge-multiplicity vector')


# Label order matches the manuscript and kernel diagrams exactly.
KERNELS = {
    'A': (2, [(0, 1)] * 4),
    'B': (3, [(0, 1), (0, 1), (1, 2), (1, 2), (0, 2)]),
    'C': (4, [(0, 1), (0, 1), (1, 2), (2, 3), (2, 3), (3, 0)]),
    'D': (4, [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]),
}


def explicit_criterion(name, lengths):
    L = sum(lengths)
    if name == 'A':
        return 3 * max(lengths) <= L
    if name == 'B':
        a, b, c, d, e = lengths
        return (3 * max(lengths) <= L and
                3 * min(a + b, c + d) >= L)
    if name == 'C':
        a, b, c, d, e, f = lengths
        return (3 * max(a, b, d, e, c + f) <= L and
                3 * min(a + b, d + e) >= L)
    a, b, c, d, e, f = lengths
    return (3 * max(lengths) <= L and
            3 * min(a + b + d, a + c + e, b + c + f, d + e + f) >= L)


def kernel_checks(max_length=6):
    result = {}
    for name, (n, edges) in KERNELS.items():
        constraints = []
        for mask in range(1, 1 << len(edges)):
            indices = [i for i in range(len(edges)) if mask >> i & 1]
            beta = len(indices) - rank(n, [edges[i] for i in indices])
            if beta:
                constraints.append((indices, beta))
        uniform = total = 0
        for lengths in it.product(range(1, max_length + 1), repeat=len(edges)):
            L = sum(lengths)
            brute = all(3 * sum(lengths[i] for i in I) >= L * b
                        for I, b in constraints)
            explicit = explicit_criterion(name, lengths)
            assert brute == explicit, (name, lengths, brute, explicit)
            uniform += brute
            total += 1
        result[name] = dict(length_assignments=total, uniformly_dense=uniform)
    return result


def choose(n, k):
    return math.comb(n, k) if n >= k >= 0 else 0


def capped(d, L, q):
    return sum((-1) ** j * choose(d, j) * choose(L - j * q - 1, d - 1)
               for j in range(3))


def formula(name, L):
    q, r = divmod(L, 3)
    if q == 0:
        return 0
    u = q - 1 + int(r > 0)
    if name == 'A':
        return capped(4, L, q)
    if name == 'B':
        return capped(5, L, q) - 2 * sum(
            (R - 1) * choose(R - r + 2, 2) for R in range(2, u + 1))
    if name == 'C':
        A = (choose(L - 1, 5) - 5 * choose(L - q - 1, 5)
             + 10 * choose(L - 2 * q - 1, 5)
             - q * choose(L - q - 1, 4)
             + 4 * q * choose(L - 2 * q - 1, 4))
        B = sum((R - 1) * ((q - 1) * choose(R - r + 2, 2)
                          - choose(R - r + 2, 3)) for R in range(2, u + 1))
        return A - 2 * B
    return capped(6, L, q) - 4 * sum(
        choose(R - 1, 2) * choose(R - r + 2, 2) for R in range(3, u + 1))


# Ascending coefficients of the numerator, separately for L = 3q + epsilon.
# These explicit polynomials are checked against exhaustive compositions below.
QUASI = {
    'A': [(6, [-18, 11, 6, 1]), (6, [0, 2, 3, 1]), (6, [0, -1, 0, 1])],
    'B': [(24, [96, -122, -5, 26, 5]), (24, [0, -6, -5, 6, 5]),
          (24, [0, 2, 7, 10, 5])],
    'C': [(120, [-480, 834, -285, -125, 45, 11]),
          (120, [0, 24, 10, -35, -10, 11]), (120, [0, -6, -15, -5, 15, 11])],
    'D': [(60, [-360, 634, -150, -175, 30, 21]),
          (20, [0, 8, 0, -15, 0, 7]), (20, [0, -2, -10, -5, 10, 7])],
}


def quasipolynomial(name, L):
    if L < 3:
        return 0
    q, epsilon = divmod(L, 3)
    denominator, coefficients = QUASI[name][epsilon]
    numerator = sum(c * q**i for i, c in enumerate(coefficients))
    assert numerator % denominator == 0
    return numerator // denominator


def quasipolynomial_checks():
    from derive_quasipolynomials import derive_all
    from fractions import Fraction
    identities = derive_all()
    for name in KERNELS:
        for epsilon in range(3):
            denominator, coefficients = QUASI[name][epsilon]
            assert identities[name][epsilon] == [Fraction(c, denominator)
                                                 for c in coefficients]
    checked = 0
    for name in KERNELS:
        for q in range(1, 101):
            for epsilon in range(3):
                L = 3 * q + epsilon
                assert quasipolynomial(name, L) == formula(name, L), (name, L)
                checked += 1
    return dict(symbolic_polynomial_identities=12,
                max_q=100, residue_classes=[0, 1, 2], comparisons=checked,
                numerator_coefficients_ascending=QUASI)


def compositions(total, parts):
    if parts == 1:
        if total >= 1:
            yield (total,)
        return
    for first in range(1, total - parts + 2):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def counting_checks(max_total=33):
    records, assignments = [], {}
    for name, (_, edges) in KERNELS.items():
        checked = 0
        for L in range(1, max_total + 1):
            actual = 0
            for lengths in compositions(L, len(edges)):
                checked += 1
                actual += explicit_criterion(name, lengths)
            expected = formula(name, L)
            assert actual == expected, (name, L, actual, expected)
            assert actual == quasipolynomial(name, L), (name, L, actual)
            records.append(dict(kernel=name, total_length=L, count=actual))
        assignments[name] = checked
    return dict(max_total=max_total, assignments=assignments, counts=records)


def subdivide(n, edges, lengths):
    result = []
    for (u, v), length in zip(edges, lengths):
        start = u
        for _ in range(length - 1):
            result.append((start, n))
            start = n
            n += 1
        result.append((start, v))
    return n, result


def counterexample():
    edges = [(0, 1), (0, 1), (0, 2), (1, 2), (0, 3), (1, 3), (2, 3)]
    lengths = (2, 2, 1, 1, 3, 3, 3)
    n, E = subdivide(4, edges, lengths)
    sig, k = signatures(n, E)
    s = max(sig.count(c) for c in set(sig))
    g = min(sum((c & z).bit_count() % 2 for c in sig) for z in range(1, 1 << k))
    assert k == 4 and 0 not in sig and 4 * s <= len(E) <= 4 * g
    assert not direct_density(n, E)
    assert len(set(tuple(sorted(e)) for e in E)) == len(E)
    for v in range(n):
        vertices = [x for x in range(n) if x != v]
        mapping = {x: i for i, x in enumerate(vertices)}
        rest = [(mapping[a], mapping[b]) for a, b in E if v not in (a, b)]
        assert rank(n - 1, rest) == n - 2
    return dict(vertices=n, edges=len(E), cycle_rank=k, largest_series_class=s,
                girth=g, witness_edges=sum(lengths[:4]), witness_cycle_rank=2,
                two_connected=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', default='verification_results.json')
    args = parser.parse_args()
    started = time.time()
    result = {}
    for name, fun in [('graphs', exhaustive_graphs), ('multigraphs', exhaustive_multigraphs),
                      ('kernels', kernel_checks), ('counts', counting_checks),
                      ('quasipolynomials', quasipolynomial_checks),
                      ('counterexample', counterexample)]:
        result[name] = fun()
        print(name, 'passed', flush=True)
    result['elapsed_seconds'] = round(time.time() - started, 3)
    result['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result['source_sha256_by_file'] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in [Path(__file__), Path(__file__).with_name('derive_quasipolynomials.py')]}
    Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ('counts', 'quasipolynomials')}, indent=2))


if __name__ == '__main__':
    main()
