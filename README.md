# Uniform density of tricyclic graphs: verification code

Exact computational checks accompanying the working manuscript *Uniform density of tricyclic graphs and their subdivisions* by Marek Gałązka and Hanna Wdowicka.

The code checks an explicit density criterion for graphs with three independent cycles, subdivision-length counting formulas for four kernels, and a counterexample at cycle rank four. All mathematical comparisons use integers. Finite verification supplements the proofs; it does not establish a theorem or priority of a result.

## Run

Requires Python 3.10 or newer. No third-party packages are needed.

From the repository root:

```sh
python3 code/verify.py --output verification_results.json
```

A successful run prints `passed` for `graphs`, `kernels`, `counts`, and `counterexample`, then writes the results as JSON. A failed comparison raises an assertion error. Run Python normally, without `-O`, which disables assertions.

The committed `verification_results.json` is a reference run. Its mathematical records are deterministic; `elapsed_seconds` depends on the machine. `source_sha256` identifies the exact verification script.

## Coverage

| Check | Scope |
| --- | --- |
| Graph criterion | All 6,286 connected labelled simple tricyclic graphs on 4–6 vertices |
| Direct density comparisons | All 1,587,378 nonempty edge subsets of those graphs |
| Four kernel inequalities | All 35,000 assignments of path lengths in `{1,2,3,4,5}` |
| Subdivision counts | All 1,357,461 positive labelled length vectors of total at most 30, including every residue modulo three |
| Boundary at cycle rank four | A simple, two-connected counterexample, checked directly |

For bridgeless inputs, the graph checks compare the exact maximum density with the three-term formula. Bridges and series classes are also tested independently by deleting edges and checking connectivity. Kernel criteria are compared with all kernel-subset density inequalities.

## Conventions

Graphs are finite, undirected, and loopless. Kernel edges are labelled and may be parallel. Subdivision counts concern positive lengths assigned to those labelled edges, not graph isomorphism classes. Kernel labels A–D and their edge order are defined in `KERNELS` in the script.

The script intentionally uses exhaustive checks for independent validation; it is not a benchmark implementation of the manuscript's linear-time recognition algorithm. Exhaustive graph checks stop at six vertices and do not assert verification of all graphs.

## Files

- `code/verify.py`: self-contained verification script.
- `verification_results.json`: executed reference results and script checksum.
- `LICENSE`: MIT license.

The manuscript and private research notes are not part of this code repository.

## License

MIT. Copyright (c) 2026 Marek Gałązka and Hanna Wdowicka.
