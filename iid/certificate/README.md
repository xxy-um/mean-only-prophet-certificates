# I.I.D. lower bound: original Appendix A.1 FRLP

A completed backward run constructs a feasible solution of the original Appendix A.1 FRLP with

\[
\gamma=\frac{815733368877}{1099511627776}
      =0.74190517705301318\ldots>0.7419.
\]

The target comparison is exact:

```text
10000 * 815733368877 - 7419 * 1099511627776 = 56922299856 > 0.
```

The full run took about **31 minutes on eight threads**. This is a feasible lower bound; no LP optimality is claimed.

The LP is the one constructed in Appendix A.1: each rate interval uses its lower endpoint for the A.2 continuation state and its upper endpoint for the applicable A.3 running-cost bounds, followed by the original grid's triangular interpolation. The fixed threshold has 40 slopes. The numerical dimensions are

- 640 time intervals;
- 1,200 state subdivisions in each coordinate;
- 131,093 rate intervals, including the final unbounded interval.

The backward evaluator uses conservative arithmetic and pruning that bounds the **original leaf rows**. It does not substitute ancestor rows, and it does not raise table entries using same-layer concavity or monotonicity. A final monotonicity pass can only lower the current layer; because the next layer is already fixed, this preserves its LP constraints. `FRLP_PROOF.md` explains these statements.

## Files

- `verify.cpp`: the frozen evaluator source used for the completed run.
- `data/`: the exact threshold, state/rate grids, and frozen binary coefficients.
- `validate_inputs.py`: exact rational validation of those coefficient inputs.
- `run.py`: compilation and replay launcher.
- `audit/INPUT_VALIDATION.json`: completed input-validation record.
- `RESULTS.json`: exact final objective, parameters, and target comparison.
- `audit/FINAL_RUN_*`: the completed run result, metadata, and log.
- `audit/README.md`: portable independent row-check instructions.

## Replay

Requires Python 3.9 or later and a GNU-compatible C++17 compiler with OpenMP and signed `__int128`. The evaluator requires binary64 arithmetic and enabled assertions. The launcher supplies strict compiler options; do not add `-ffast-math`, `-DNDEBUG`, or floating contraction/reassociation options.

From this directory:

```sh
python3 run.py --threads 8
```

This validates the input coefficients, compiles the evaluator, and performs the full backward run. The exact final objective is reported as an integer numerator over `2^40`; the decimal is only a display value. The launcher compares the final fraction with `7419/10000` using integers.

A run writes its table checkpoints and logs under `runs/run`. To resume that same run:

```sh
python3 run.py --threads 8 --resume
```

The source, inputs, and numerical parameters must match the saved run metadata. Each table checkpoint uses signed 64-bit integer entries in triangular order `i*(i+1)/2+j`, with scale `2^40`.

This package concerns the asymptotic FRLP lower bound. The article's separate analytical transfer theorem supplies the all-finite-`n` statement; no historical finite-`n` loss estimate is used here.


## Numerical value reported in the manuscript

The theorem-level guarantee **0.7419** is reached with the same 40 slopes and
the same displayed dimensions. The finer objective associated with this
certificate is **0.74190517705301318...**. If Section 3.3 reports a finer
numerical value for this instance, it should use this value instead of the
old strengthened-table value `0.74191589...`.
