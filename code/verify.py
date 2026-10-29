#!/usr/bin/env python3
"""Exact, dependency-free checks for the tricyclic uniform-density draft.

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
    parent = list(range(n))
    tree = [[] for _ in range(n)]
    chords = []
    def root(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v
    for i, (u, v) in enumerate(edges):
        a, b = root(u), root(v)
        if a != b:
            parent[a] = b
            tree[u].append((v, i))
            tree[v].append((u, i))
        else:
            chords.append(i)
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


def signature_criterion(n, edges):
    sig, k = signatures(n, edges)
    assert k == 3
    if 0 in sig:
        return False
    weights = [sig.count(i) for i in range(8)]
    shortest = min(sum(weights[c] for c in range(1, 8)
                       if (c & z).bit_count() % 2) for z in range(1, 8))
    return 3 * max(weights) <= len(edges) <= 3 * shortest


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
            a, maximum = direct_density(n, edges, return_max=True)
            b = signature_criterion(n, edges)
            assert a == b, (n, edges, a, b)
            uniform += a
            comparisons += (1 << m) - 1
            sig, _ = signatures(n, edges)
            # Independently verify bridge and series-class interpretations.
            for i in range(m):
                no_i = [e for j, e in enumerate(edges) if j != i]
                assert (sig[i] == 0) == (rank(n, no_i) < n - 1)
            if 0 not in sig:
                s = max(sig.count(c) for c in set(sig))
                g = min(sum((c & z).bit_count() % 2 for c in sig)
                        for z in range(1, 8))
                candidates = [(m, m - 3), (g, g - 1), (m - s, m - s - 2)]
                assert all(p * maximum[1] <= maximum[0] * q for p, q in candidates)
                assert any(p * maximum[1] == maximum[0] * q for p, q in candidates)
                for i, j in it.combinations(range(m), 2):
                    remainder = [e for h, e in enumerate(edges) if h not in (i, j)]
                    assert (sig[i] == sig[j]) == (rank(n, remainder) < n - 1)
        records.append(dict(vertices=n, connected_graphs=connected,
                            uniformly_dense=uniform))
    return dict(records=records, subset_comparisons=comparisons,
                total_graphs=sum(x['connected_graphs'] for x in records))


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


def kernel_checks(max_length=5):
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


def compositions(total, parts):
    if parts == 1:
        if total >= 1:
            yield (total,)
        return
    for first in range(1, total - parts + 2):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def counting_checks(max_total=30):
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
    for name, fun in [('graphs', exhaustive_graphs), ('kernels', kernel_checks),
                      ('counts', counting_checks), ('counterexample', counterexample)]:
        result[name] = fun()
        print(name, 'passed', flush=True)
    result['elapsed_seconds'] = round(time.time() - started, 3)
    result['source_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(args.output).write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'counts'}, indent=2))


if __name__ == '__main__':
    main()
