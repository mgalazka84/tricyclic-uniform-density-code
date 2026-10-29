#!/usr/bin/env python3
"""Exact symbolic expansion of the counting decompositions (standard library).

Polynomials are ascending coefficient lists over Fraction. Finite differences
convert each polynomial summand to the binomial basis; the hockey-stick
identity then sums it symbolically. No interpolation of subdivision counts
or floating-point arithmetic is used.
"""
from fractions import Fraction as F
from math import comb
import json


def trim(p):
    p = list(map(F, p))
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(*polynomials):
    out = [F(0)] * max(map(len, polynomials))
    for p in polynomials:
        for i, a in enumerate(p):
            out[i] += a
    return trim(out)


def scale(p, a):
    return trim([a * x for x in p])


def mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return trim(out)


def value(p, x):
    return sum(a * x**i for i, a in enumerate(p))


def binomial(p, k):
    """The polynomial binom(p(q), k), with k a fixed nonnegative integer."""
    out = [F(1)]
    for j in range(k):
        out = scale(mul(out, add(p, [-j])), F(1, j + 1))
    return out


def symbolic_sum(summand, upper, lower):
    """Sum a polynomial in R from lower to upper(q), as a polynomial in q."""
    differences = [value(summand, i) for i in range(len(summand))]
    out, order = [F(0)], 0
    while differences:
        out = add(out, scale(binomial(add(upper, [1]), order + 1), differences[0]))
        differences = [b - a for a, b in zip(differences, differences[1:])]
        order += 1
    return add(out, [-sum(value(summand, i) for i in range(lower))])


def derive_all():
    out = {name: [] for name in 'ABCD'}
    q, R = [0, 1], [0, 1]
    for epsilon in range(3):
        upper = [-int(epsilon == 0), 1]
        L = [epsilon, 3]
        def capped(d):
            return add(*[scale(binomial(add(L, scale(q, -j), [-1]), d - 1),
                               (-1)**j * comb(d, j)) for j in range(3)])
        pair_summand = mul(add(R, [-1]), binomial(add(R, [2 - epsilon]), 2))
        pair_sum = symbolic_sum(pair_summand, upper, 2)
        connector_summand = mul(add(R, [-1]), binomial(add(R, [2 - epsilon]), 3))
        connector_sum = symbolic_sum(connector_summand, upper, 2)
        triangle_summand = mul(binomial(add(R, [-1]), 2),
                               binomial(add(R, [2 - epsilon]), 2))
        triangle_sum = symbolic_sum(triangle_summand, upper, 3)
        W = add(binomial(add(L, [-1]), 5),
                scale(binomial(add(L, scale(q, -1), [-1]), 5), -5),
                scale(binomial(add(L, scale(q, -2), [-1]), 5), 10),
                scale(mul(q, binomial(add(L, scale(q, -1), [-1]), 4)), -1),
                scale(mul(q, binomial(add(L, scale(q, -2), [-1]), 4)), 4))
        out['A'].append(capped(4))
        out['B'].append(add(capped(5), scale(pair_sum, -2)))
        out['C'].append(add(W, scale(mul(add(q, [-1]), pair_sum), -2),
                            scale(connector_sum, 2)))
        out['D'].append(add(capped(6), scale(triangle_sum, -4)))
    return out


if __name__ == '__main__':
    from verify import QUASI
    polynomials = derive_all()
    for name, rows in polynomials.items():
        for epsilon, row in enumerate(rows):
            denominator, coefficients = QUASI[name][epsilon]
            assert row == [F(c, denominator) for c in coefficients]
    print('All 12 symbolic polynomial identities passed.')
    print(json.dumps({name: [[str(x) for x in row] for row in rows]
                      for name, rows in polynomials.items()}, indent=2))
