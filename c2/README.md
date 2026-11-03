# Appendix C.2: finite dual certificate

This directory contains the numerical certificate used in Appendix C.2,
**“Precommitted Thresholds in the I.I.D. Setting,”** of *Mean-only Prophet
Inequalities*. It certifies the finite coefficient conditions
of Lemma C.11 with

$$
N=1024,\qquad |\mathcal Q|=1014,\qquad
\widehat U\le 0.741966695100 < U=0.741967.
$$

The certificate concerns the **upper bound for precommitted thresholds**: the
paper considers deterministic, continuous, strictly decreasing thresholds with
$\theta(1)=0$, fixed before the values and their independent uniform arrival
times are sampled. The threshold may depend on $n$. The mathematical reduction
from these coefficient inequalities to hard i.i.d. inputs is in the paper. For finite $n$, Lemma C.12 gives the bound $e^{1/n}U+\eta$, for every $\eta>0$, under its size condition;
$U$ is the resulting asymptotic upper bound.

## Verify the certificate

Download and extract the complete repository or release ZIP. In a terminal,
change to the `c2` directory, then run:

```sh
python3 verify_certificate.py --output-dir runs/recheck
```

Use **Python 3.10 or later**. The verifier requires only the Python standard
library; there is no package installation, solver, network access, compilation,
or random seed. The default input is `data/upper_1024.json.gz`, resolved relative
to the script. A different input path can be passed using `--data`.

A successful run exits with status zero only after checking the data and the
target inequalities. It writes `verification.json`, `coefficient_bounds.json`,
and `verification.log` in the specified output directory. A failed check exits
with nonzero status. Use a fresh output directory to distinguish separate runs.

The reproducible mathematical output is

```text
N = 1024
rate vectors = 1014
interval scale = 1461501637330902918203684832716283019655932542976
upper numerator = 1084385539732742468991572064431813102297740351518
upward-rounded upper bound = 0.741966695100
target U = 0.741967
```

The numerator divided by the interval scale is the exact rational
$\widehat U$ produced by the interval computation. Timing, timestamps, Python
version, and platform information may differ between machines. The integer
result does not depend on those fields. The archived run for this release is
in [`results/`](results/).

The verifier evaluates more than one million rate entries. Its cache can use
substantial memory; this is a local computation, not an LP optimization.

## Files

| File | Purpose |
| --- | --- |
| [`data/upper_1024.json.gz`](data/upper_1024.json.gz) | Exact rational certificate, preserved byte-for-byte from the supplied artifact bundle. |
| [`verify_certificate.py`](verify_certificate.py) | Standalone verifier with exact input validation and outward-rounded interval arithmetic. |
| [`docs/Verification of the Appendix C.2 Certificate.pdf`](docs/Verification%20of%20the%20Appendix%20C.2%20Certificate.pdf) | Correctness proof of the released verification procedure. |
| [`docs/C2_Certificate_Verification.tex`](docs/C2_Certificate_Verification.tex) | LaTeX source of the proof note. |
| [`results/verification.json`](results/verification.json) | Actual run record, exact bound, input/code hashes, and target checks. |
| [`results/coefficient_bounds.json`](results/coefficient_bounds.json) | Exact interval endpoints for all 1025 weighted coefficient averages. |
| [`results/verification.log`](results/verification.log) | Human-readable output of the archived run. |
| [`tests/test_verifier.py`](tests/test_verifier.py) | Focused tests of interval arithmetic, exponential bounds, and input validation. |
| [`PROVENANCE.json`](PROVENANCE.json) | Source artifact identifiers and hashes, manuscript references, and packaging details. |
| [`SHA256SUMS`](SHA256SUMS) | SHA-256 checksums of the released files, excluding the checksum file itself. |

## Exact data format and notation

The compressed file contains the complete JSON certificate. It can be read
directly with Python's standard `gzip` module or decompressed with `gzip -dc`;
a second uncompressed copy is not required for verification.

The JSON contains the following fields. Integer arrays are not floating-point
approximations; the indicated divisions define the exact rates and weights.

| JSON field | Meaning in the paper |
| --- | --- |
| `N` | Number of equal time subintervals; here 1024. |
| `rate_denominator` | Common rate denominator, $10^8$. |
| `profiles` | Each row stores the numerators of one $q=(q_0,\ldots,q_{N-1})$. |
| `weight_denominator` | Common weight denominator, $10^{12}$. |
| `weights` | Numerators of $\lambda_q$, in the same order as the rows of `profiles`. |
| `cutoffs` | Auxiliary metadata retained from the original artifact; it is not used by the verifier or by the certificate proof. |

The verifier checks all dimensions and integer types, the expected denominators,
nonnegative nondecreasing rate vectors, nonnegative weights, and a weight-
numerator sum exactly equal to $10^{12}$. It checks that all rates are at most
$10^9$ and reports the actual maximum. In this certificate that maximum is
$10^9$, so Lemma C.12 applies for $n\ge 1{,}000{,}000{,}001$.

The mathematical notation in the proof follows the paper:
$S_q$, $a_j$, $A_q$, $k_{q,j}$, $\lambda_q$, and $\widehat U$.
The arithmetic scale is $B=2^{160}$. Interval endpoints in the output are
integers divided by $B$; no interval endpoint is a decimal approximation.

## What is certified

The program computes enclosures for

$$
\sum_{q\in\mathcal Q}\lambda_q A_q(0)
\quad\text{and}\quad
\sum_{q\in\mathcal Q}\lambda_q(1-k_{q,j}),\qquad 0\le j<N.
$$

It includes the $A_q(0)$ coefficient as well as all 1024 coefficients involving
$k_{q,j}$. The exact upper endpoint used for the latter is one minus the
lower endpoint of the weighted $k_{q,j}$ sum, using the verified equality
$\sum_q\lambda_q=1$.

All arithmetic affecting the certificate uses integers and directed rounding.
The exponential enclosure is proved in the PDF; the verifier does not call a
floating-point exponential function. The final comparisons are exact integer
comparisons with $741966695100/10^{12}$ and $741967/10^6$.

The proof note establishes that passing these checks implies all hypotheses
of Lemma C.11. The stored run records the checks on the released data. Neither
a solver's feasibility status nor a pre-entered expected answer is used as
a substitute for computing the coefficients.

## Additional checks

Run the focused implementation tests with:

```sh
python3 -m unittest discover -s tests -v
```

On macOS, check the release file hashes from the `c2` directory with:

```sh
shasum -a 256 -c SHA256SUMS
```

On systems with GNU coreutils, the corresponding command is:

```sh
sha256sum -c SHA256SUMS
```

Running the verifier into `runs/recheck` does not alter the archived
`results/` files. Those archived files are included in the checksum list.

## Provenance and scope

The certificate was extracted from the supplied `Mean_Only_IID_Artifacts(2).zip`
artifact bundle. Its compressed bytes are unchanged. The released verifier
is a readable standalone implementation of the original integer interval
calculation, with explicit input errors, exact target checks, and full
coefficient output. Its evaluation order is retained to reproduce the original
exact numerator. `PROVENANCE.json` records an independent replay of the
original verifier on the same certificate and the matching exact bound; the
original verifier and source archive are historical references, not required
inputs to this package.

The certificate-generation/LP-search procedure is not needed to verify a fixed
finite dual witness and is not included in this focused C.2 release. No claim
of optimality of that witness is required or made.
