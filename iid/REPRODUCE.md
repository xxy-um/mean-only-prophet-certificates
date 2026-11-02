# Reproduction instructions

## Requirements

- Python 3.9 or later; no third-party Python packages.
- A GNU-compatible C++17 compiler with OpenMP and signed `__int128`.
- IEEE binary64 floating-point arithmetic with round-to-nearest, gradual underflow (no flush-to-zero or denormals-are-zero mode), and enabled assertions.
- About 1.5 GB of free disk space for the default full-run checkpoints and logs.

The launcher selects `-O3 -std=c++17 -fopenmp -ffp-contract=off`. Do not add fast-math, reassociation, contraction, or `-DNDEBUG` options. The frozen evaluator rejects fast-math and disabled assertions and checks the rounding mode on its OpenMP workers.

## Check the supplied files and exact input bounds

Run at the release root:

```sh
python3 reproduce.py --check-only
```

This checks the original SHA-256 manifest and additional fixed-instance assumptions, including equality of the actual text state grid and the validated JSON grid, the exact threshold data, the rate guards, the lower constant for `log(2)`, and the integer-range bound. It then runs the original exact-rational coefficient validator. This may take several minutes.

Passing these checks establishes the stated input enclosures. It does not by itself check every backward LP constraint; the correctness proof and the complete backward evaluation supply that step.

## Recompute the complete certificate

```sh
python3 reproduce.py --full --threads 8
```

The program checks the inputs, compiles the original evaluator, and evaluates all layers from `k=639` through `k=0`. It accepts the numerical result only if the exact numerator and denominator match the recorded certificate. Decimal output is for display.

The original launcher refuses a nonempty output directory unless resuming an existing run. To select an output directory or compiler, see:

```sh
python3 reproduce.py --help
```

For the original launcher's lower-level interface:

```sh
python3 certificate/run.py --help
```

The output includes the exact final result, source/input identities, compiler command, layer log, and periodic binary table checkpoints. Checkpoints may be retained locally for independent row audits.

## Independent checks of individual original rows

The original `certificate/audit/README.md` documents the long-double diagnostic scanner and the exact-rational interval row checker. The diagnostic scanner inspects every rate row at selected states. The rational checker rigorously encloses specified rows. Neither selected-state procedure replaces the complete backward proof and full evaluation.

The release's `validation/` directory records what was actually executed while preparing these files. Historical records inside `certificate/audit/` belong to the supplied archive; the new release records are identified separately.

## Build the proof PDF

With a standard LaTeX installation containing `fontenc`, `lmodern`, `amsmath`, `amssymb`, `amsthm`, `mathtools`, `booktabs`, `longtable`, `array`, `microtype`, `geometry`, `hyperref`, and `enumitem`:

```sh
cd docs
pdflatex -interaction=nonstopmode -halt-on-error verification.tex
pdflatex -interaction=nonstopmode -halt-on-error verification.tex
```

The supplied PDF is already compiled. No LaTeX installation is needed to run the certificate.
