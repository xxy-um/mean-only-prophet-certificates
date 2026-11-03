# Validation records

`RELEASE_CHECKS.json` describes the currently distributed files and proof PDF.
`PREFLIGHT.json` records the current integrity and exact-input checks. These checks
validate the 41 entries in the current `certificate/SHA256SUMS.json`; they do not
claim that a new complete backward run was performed during packaging.

The remaining files are dated evidence from the historical validation of the
same numerical evaluator, exact input data, and final checkpoints:

- `FULL_REPRODUCTION.json` and `full_run.log`: all 640 backward layers, exact
  objective, and the three final checkpoint comparisons.
- `INPUT_VALIDATION.json` and `exact_assumptions.json`: coefficient enclosures
  and additional exact assumptions checked for that run.
- `INDEPENDENT_CHECKS.json` and `MIDPOINT_CHECKS.json`: explicitly scoped
  selected-state and exact-row checks.
- `ENVIRONMENT.json`: the compiler, interpreter, commands, and source identities
  used for that historical execution.

Historical statements such as `launcher_source_unchanged` and the count of 42
original manifest entries describe that execution. The current release adds
Python-optimization guards to the launcher, coefficient validator, and exact-row
checker, updates documentation, removes a superseded proof note, and normalizes
historical workspace path prefixes. The current certificate manifest therefore
has 41 entries and a new digest. The mathematical evaluator, numerical inputs,
recorded arithmetic results, timings, and checkpoints remain unchanged.

The complete feasibility argument is in `../docs/verification.pdf`; selected-row
checks are supplementary diagnostics and do not replace that argument.
