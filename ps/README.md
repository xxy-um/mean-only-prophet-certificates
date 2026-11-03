# Prophet Secretary: a certified 160-segment threshold

This supplementary repository contains an independently checkable finite certificate for the prophet secretary lower bound in *Mean-only Prophet Inequalities*. It certifies

$$
\gamma\ge\frac{49839735470020217}{72057594037927936}
=0.69166527325054429076\ldots>0.6916.
$$

The threshold has 160 equal intervals and strictly increasing positive rational drop rates, with

$$
\theta(0)=\frac{567246577807}{800000000000}=0.70905822225875.
$$

The certificate covers the entire stopping-time domain using 5,120 intervals and uses 10,240 common probability-averaging cells. The verifier constructs one fixed affine lower bound on each whole stopping-time interval and checks all 10,241 suffix constraints, including the empty suffix.

## Mathematical verification

Read [Verification.pdf](Verification.pdf) for the complete supplementary verification, including the parameter interpretation, coefficient bounds, finite minimization, and arithmetic enclosures. Its editable source is the single self-contained file [Verification.tex](Verification.tex).

The proof specializes the paper's original reward criterion and Lemmas 4.7, 4.8, B.1, and B.2. In particular, the survival coefficient includes the constant deduction prescribed by the coefficient-error inequality in Lemma B.1. The choices of lower-bound planes remain fixed throughout each stopping-time interval. The certificate proves feasibility and the stated lower bound; it makes no optimality claim.

## Run the checks

Requirements: Python 3.9 or later and a C++17 compiler available as `g++` with support for `__int128`. Verification and proposal reconstruction use no numerical optimizer or external Python package. Run all commands from the `ps/` directory containing `replay.py`.

First check all inputs, integrity hashes, common parameters, plane mass budgets, and complete interval coverage:

```sh
python3 replay.py --audit-only
```

Then recompile the verifier and recompute every arithmetic record:

```sh
python3 replay.py --jobs 8
```

The final report must have `status: VERIFIED`, `all_5120_arithmetic_checks_replayed: true`, and `arithmetic_checks_replayed: 5120`. Acceptance uses the exact rational target `1729/2500`, not a rounded decimal or an optimization solver's status. The number of jobs may be reduced to suit the available machine.

A shorter diagnostic run is also available:

```sh
python3 replay.py --jobs 4 --critical-only
```

It audits every input and the whole interval manifest, then recomputes four selected intervals. Its report explicitly identifies the partial arithmetic coverage; use the full command above for complete verification.

The validation regression tests check rejection of altered targets, incorrect
checksums, invalid job counts, and inconsistent result metadata under ordinary
Python, `python3 -O`, and `PYTHONOPTIMIZE=1`:

```sh
python3 -m unittest discover -s tests -v
```

## Rebuild the proposal files

```sh
python3 generation/rebuild_proposals.py --output-dir regenerated_proposals --compare-dir proposals
```

This reconstructs all 5,120 final proposal files from the included frozen rational candidate inputs and selection schedule. It checks their recorded SHA-256 hashes and compares every output byte for byte with the supplied proposals. The output directory must differ from `proposals/`.

This is deterministic reconstruction from supplied inputs. The earlier optimization search that discovered the threshold and candidate planes is outside this reconstruction scope. Mathematical validity is checked directly from the final rational data and does not depend on that search. See [generation/README.md](generation/README.md) for details.

## Compile the verification document

With a LaTeX installation providing the packages declared in `Verification.tex`, run:

```sh
pdflatex -interaction=nonstopmode -halt-on-error Verification.tex
pdflatex -interaction=nonstopmode -halt-on-error Verification.tex
```

Repeat the compilation if LaTeX requests another run for cross-references. All mathematical content is in that one source file; no separate TeX chapters or bibliography file are required. TeX versions and PDF metadata can change the bytes of a recompiled PDF. Preserve the distributed files when checking their published integrity manifest, or compile a copy in another directory.

## Repository contents

| Path | Purpose |
| --- | --- |
| `Verification.pdf`, `Verification.tex` | Complete supplementary verification and its single-file source |
| `parameters.json` | Exact threshold and allocation parameters |
| `proposals/` | Fixed plane choices and masses for all 5,120 stopping-time intervals |
| `intervals.json` | Complete stopping-time interval assignments |
| `records.json` | Directed-arithmetic records to be reproduced |
| `sources/` | C++ verifier, interval arithmetic, and logarithmic-plane checks |
| `replay.py` | Input audit and fresh arithmetic replay |
| `tests/` | Validation regression tests |
| `RESULTS.json` | Exact minimum, target, and certificate dimensions |
| `certificate.json` | SHA-256 integrity manifest for the distributed files |
| `generation/` | Frozen candidate inputs and deterministic proposal reconstruction |
| `audit/` | Supplied verification and audit records |
| `RELEASE_NOTES.md` | Certified values, manuscript consistency notes, and reproducibility scope |

The certificate can be verified locally using the commands above.
