# Mean-only Prophet Inequalities: numerical certificates

Exact certificate inputs, verification programs, and supplementary correctness proofs.

| Directory | Certified result | Correctness proof |
| --- | --- | --- |
| `iid/` | `815733368877 / 1099511627776 > 0.7419` | [IID verification](iid/docs/verification.pdf) |
| `ps/` | `49839735470020217 / 72057594037927936 > 0.6916` | [PS verification](ps/Verification.pdf) |
| `c2/` | `U_hat <= 0.741966695100 < 0.741967` | [C2 verification](c2/docs/Verification%20of%20the%20Appendix%20C.2%20Certificate.pdf) |

## Requirements

Python 3.10 or later; no external Python packages. IID and PS also require
GNU C++17 with signed 128-bit integers (`g++`); IID additionally requires OpenMP.
IID requires binary64 round-to-nearest arithmetic with gradual underflow.
Run Python without `-O` or `PYTHONOPTIMIZE`. Use the supplied compiler options;
do not enable fast-math, floating-point contraction, or disabled C++ assertions.

## Verification

Run from the repository directory containing this README:

```sh
python3 iid/reproduce.py --full --threads 8
python3 ps/replay.py --jobs 8
python3 c2/verify_certificate.py --output-dir c2/runs/recheck
```

The IID command validates every frozen coefficient, reconstructs all 640 backward
layers, and checks the exact objective and final checkpoint digests. Reserve about
1.5 GB for generated outputs and allow tens of minutes. Its output directory is
`iid/runs/reproduction/`; resume an interrupted run with the same command plus
`--resume`, or select a fresh location with `--output`.

The PS command verifies all 5,120 complete stopping-time intervals and all
10,241 suffix constraints per interval. It computes fresh arithmetic records and
checks their exact minimum. Outputs go to `ps/runs/recheck/`; select a fresh
location with `--output-dir` for another run.

The C2 command computes all 1,025 weighted coefficient bounds with integer
interval arithmetic. It writes the bound, coefficient intervals, and log to its
specified output directory. The cache may require substantial memory.

All three commands must exit successfully. The exact accepted results are shown
in the table above. Runtime output is generated locally and is not an input.

For input checks only, use `python3 iid/reproduce.py --check-only` or
`python3 ps/replay.py --audit-only`. These shorter commands do not perform the
complete numerical verification.

## Mathematical scope

The IID certificate proves feasibility of the original FRLP in Section 3.3 and
Appendix A.1. Appendix A.2 of the paper supplies the transfer to every finite n.
The PS certificate verifies the original pointwise criterion in Section 4 and
Appendix B.2. The C2 certificate supplies the finite dual witness in Lemma C.11;
Lemma C.12 gives the finite-n bound `exp(1/n) U + eta`, for every positive `eta`,
under its size condition (`n >= 1,000,000,001` for this witness). Its asymptotic
upper bound concerns the specified precommitted threshold class.

The three proof PDFs are part of the verification argument. No optimality of a
numerical witness or reconstruction of the exploratory search is claimed.

## File integrity

From this directory, run `sha256sum -c SHA256SUMS`, or on macOS,
`shasum -a 256 -c SHA256SUMS`. The verifiers also check their required numerical
inputs. Generated output directories are not part of the distributed package.
