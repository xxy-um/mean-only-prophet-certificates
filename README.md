# Mean-only Prophet Inequalities

This repository contains the numerical certificates and supplementary verification proofs for **Mean-only Prophet Inequalities**. Each directory is self-contained and has its own reproduction instructions and integrity manifest.

| Directory | Certified result | Supplementary proof | Manuscript location |
| --- | --- | --- | --- |
| [iid/](iid/) | I.I.D. lower bound: `815733368877 / 1099511627776 > 0.7419` | [Verification PDF](iid/docs/verification.pdf) | Section 3.3; Appendix A.1 |
| [ps/](ps/) | Prophet secretary lower bound: `49839735470020217 / 72057594037927936 > 0.6916` | [Verification PDF](ps/Verification.pdf) | Section 4; Appendix B.2 |
| [c2/](c2/) | Finite dual witness: `U_hat <= 0.741966695100 < 0.741967` | [Verification PDF](c2/docs/Verification%20of%20the%20Appendix%20C.2%20Certificate.pdf) | Appendix C.2; Lemmas C.11 and C.12 |

## Verification

Run these commands from the repository root. Python 3.10 or later suffices for all three packages. The IID and PS packages also require GNU C++17 with signed 128-bit integers; IID additionally uses OpenMP. C2 requires only the Python standard library.

```sh
# IID: validate the frozen inputs, then reproduce all 640 backward layers.
python3 iid/reproduce.py --check-only
python3 iid/reproduce.py --full --threads 8

# PS: audit all inputs, then reproduce all 5,120 arithmetic records.
python3 ps/replay.py --audit-only
python3 ps/replay.py --jobs 8

# C2: recompute all 1,025 weighted coefficient bounds.
python3 c2/verify_certificate.py --output-dir c2/runs/recheck
```

The full IID run takes tens of minutes and needs approximately 1.5 GB of free disk space. See [IID reproduction instructions](iid/REPRODUCE.md), [PS instructions](ps/README.md), and [C2 instructions](c2/README.md) for requirements, output locations, additional checks, and the meaning of successful verification. Input checks and selected-row diagnostics do not replace a full arithmetic replay.

## Mathematical scope

The IID certificate proves feasibility for the paper's original fixed-coefficient factor-revealing LP. The analytical transfer to every finite input length is established separately in Appendix A.2 of the manuscript.

The PS certificate checks one fixed affine lower bound over each complete stopping-time interval and all probability suffix constraints, including the empty suffix. Its final rational proposals can also be reconstructed from the supplied frozen candidate inputs and schedule.

The C2 witness concerns the paper's precommitted deterministic, continuous, strictly decreasing threshold rules under independent uniform arrival times. The value `0.741967` is an asymptotic upper bound for this class. For finite input size, Lemma C.12 gives `exp(1/n) U + eta`, for every positive `eta`, subject to its size condition; the supplied witness requires `n >= 1,000,000,001`.

The supplementary PDFs prove the numerical implementation and error bounds used by the verifiers. They are part of the certification argument. None of the packages claims optimality of its numerical witness or reconstruction of the exploratory optimization search.

## Files and integrity

The three supplementary verification proofs are supplied as PDFs. No LaTeX installation is needed for numerical verification.

Exact inputs, numerical verifier sources, and recorded numerical evidence are included. Check the distributed files against [IID checksums](iid/checksums.sha256), [PS certificate manifest](ps/certificate.json), and [C2 checksums](c2/SHA256SUMS) before changing them.

Manuscript metadata identifies the work by its title, **Mean-only Prophet Inequalities**.
