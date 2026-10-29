# Uniformly dense tricyclic subdivisions: counts and certificates

Exact computational checks accompanying the manuscript by Marek Gałązka and Hanna Wdowicka.

The code implements a linear-time recognition algorithm with density and obstruction certificates, checks the corank-three density criterion, and derives twelve explicit polynomials counting positive labelled subdivision lengths for four tricyclic kernels. The criterion is a consequence of matroid duality; it is not claimed here as an independent new structural theorem. Finite verification supplements the proofs and does not establish mathematical priority.

## Run

Requires Python 3.10 or newer. No third-party packages are needed.

From the repository root:

```sh
python3 code/verify.py --output verification_results.json
```

A successful run prints `passed` for `graphs`, `multigraphs`, `kernels`, `counts`, `symbolic`, and `counterexample`, then writes the results as JSON. A failed comparison raises an assertion error. Run Python normally, without `-O`, which disables assertions.

The symbolic derivation can also be run separately:

```sh
python3 code/derive_quasipolynomials.py
```

This expands the exact bounded-composition sums over rational coefficients and uses the binomial-basis summation identity. It does not fit polynomials to subdivision counts. For each kernel, the three coefficient lists correspond to `L = 3q + epsilon`, with `epsilon = 0, 1, 2` and `q >= 1`; coefficients are in ascending powers of `q`. Counts at `L = 1, 2` are zero.

The committed `verification_results.json` is an executed reference run. Its mathematical records are deterministic; `elapsed_seconds` depends on the machine. `source_sha256_by_file` identifies both exact scripts, while `source_sha256` retains the main script checksum.

## Coverage

| Check | Scope |
| --- | --- |
| Simple graphs | All 6,286 connected labelled simple tricyclic graphs on 4–6 vertices; 1,587,378 nonempty edge subsets |
| Multigraphs | All 7,540 connected loopless tricyclic edge-multiplicity vectors on 2–5 labelled vertices; 933,916 nonempty edge subsets |
| Four kernel inequalities | All 102,384 assignments of path lengths in `{1,2,3,4,5,6}` |
| Subdivision counts | All 2,493,392 positive labelled length vectors of total at most 33, including every residue modulo three |
| Symbolic formulas | Exact coefficient comparison for all twelve polynomials; a further 1,200 sum/formula comparisons for `1 <= q <= 100` |
| Boundary at cycle rank four | A simple, two-connected counterexample, checked directly |

The multigraph group includes 121 simple inputs already present in the first group; these are not disjoint graph totals. Parallel copies are distinct when enumerating edge subsets. The multigraph convention enumerates each edge-multiplicity vector once, without separately labelling indistinguishable copies as input graphs.

For bridgeless inputs, direct exhaustive density checks validate the exact maximum and the edge set returned by `recognize_tricyclic`. For every nonuniform input, the returned obstruction is checked directly. Bridges and series classes are also tested independently by edge deletion and connectivity. Kernel criteria are compared with all kernel-subset density inequalities.

## Conventions and files

Graphs are finite, undirected, and loopless. Kernel edges are labelled and may be parallel. Counts concern positive lengths assigned to those labelled edges, not graph isomorphism classes. Kernel labels A–D and their edge order are defined in `KERNELS`; the fifth edge of kernels B and C is denoted `p` in the manuscript.

- `code/verify.py`: recognition algorithm, independent exhaustive checks, and explicit residue-class formulas.
- `code/derive_quasipolynomials.py`: exact symbolic derivation with Python's standard-library rational arithmetic.
- `verification_results.json`: executed reference results and script checksums.
- `CITATION.cff`: software citation metadata.
- `.zenodo.json`: metadata prepared for a Zenodo software archive.
- `LICENSE`: MIT license.

Exhaustive verification is intentionally limited to the stated ranges. It is separate from the linear-time recognition algorithm, whose input contract is a connected loopless tricyclic graph. All comparisons use exact integers or rational numbers.

The manuscript and private research notes are not part of this code repository. A Zenodo DOI has not yet been assigned; the metadata file alone is not an archived release.

## License

MIT. Copyright (c) 2026 Marek Gałązka and Hanna Wdowicka.
