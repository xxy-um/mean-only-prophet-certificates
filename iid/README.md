# I.I.D. lower-bound certificate

This directory is the publication package for the numerical certificate in Section 3.3 and Appendix A.1 of **Mean-only Prophet Inequalities**.

The certificate constructs a feasible solution of the paper's original fixed-coefficient factor-revealing LP. Its exact objective is

```
815733368877 / 1099511627776
```

It uses 40 threshold slopes, 640 time intervals, 1,200 state subdivisions in each coordinate, and 131,093 rate intervals, including the last unbounded interval. No optimality claim is made.

## Start here

- Open `index.html` for the public file index.
- Read `docs/verification.pdf` for the mathematical and numerical correctness proof.
- Use `reproduce.py` to check the frozen inputs and reproduce the backward evaluation.
- Read `REPRODUCE.md` for exact commands, requirements, and the scope of each check.
- Inspect `validation/` for the independent release-validation records.

`validation/README.md` explains the scope and dates of the evidence;
`validation/RELEASE_CHECKS.json` records the current packaging checks.
The dated full-run and row-check records retain the historical executions of
the unchanged numerical evaluator and inputs.

The proof follows the paper's notation. In particular, `N=40` denotes the number of threshold pieces in the mathematical document. The C++ constant named `N` denotes `N_t=640`; the source-to-paper mapping is given in the proof.

## Reproduce

From this directory, with Python 3.9 or later and GNU C++17 with OpenMP and signed 128-bit integers:

```sh
python3 reproduce.py --check-only
python3 reproduce.py --full --threads 8
```

The first command checks integrity, the exact input relationships, and all frozen coefficient enclosures. The second additionally compiles and runs all 640 backward layers and compares the exact result with the recorded certificate. The full run writes its own checkpoints and logs; reserve approximately 1.5 GB of free disk space and allow tens of minutes, depending on hardware.

## Files and provenance

| Path | Purpose |
| --- | --- |
| `docs/verification.pdf` | Primary correctness document for this release |
| `certificate/` | Frozen evaluator, exact inputs, checkpoints, guarded launchers, and audit records |
| `reproduce.py` | Release preflight checks and reproducibility entry point |
| `REPRODUCE.md` | Requirements, commands, and interpretation of the checks |
| `validation/` | Evidence produced while preparing this release |
| `PAPER_REFERENCE.json` | Manuscript title and frozen-source identity |
| `checksums.sha256` | SHA-256 checksums for the complete release |

The primary proof is `docs/verification.pdf`, distributed as a compiled PDF only. It supplies the implementation justification and numerical error analysis. The evaluator, exact numerical inputs, and frozen checkpoints are unchanged. The launchers reject disabled Python assertions; documentation and historical path labels have been normalized for this release. The current manifests cover the distributed files.

The final rational inputs are supplied. The historical heuristic search that discovered the threshold is not needed to verify its guarantee and is not reconstructed by this package.

The root README and proof concern the i.i.d. lower-bound certificate only. The paper supplies the separate analytical transfer to every finite input length.
