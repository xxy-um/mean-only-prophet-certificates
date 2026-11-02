# Rebuilding the final proposal files

This directory contains the deterministic final assembly used for the reported
0.6916 certificate. It takes 165 frozen candidate plane lists, a fixed rational
seed schedule, and the common rational threshold and allocation parameters.
It rebuilds all 5,120 final proposal files using only Python's standard library.

From the certificate's root directory, run:

```sh
python3 generation/rebuild_proposals.py --output-dir regenerated_proposals --compare-dir proposals
```

Every generated file is checked against its recorded SHA-256. The optional
`--compare-dir` also compares all files byte for byte with the shipped final
proposals. The generator writes a `REBUILD_REPORT.json` in its output directory.

For each complete stopping-time interval, the rule chooses two candidate lists
from `plane_choice_seeds.json`, evaluates their rational mixture weight at the
interval midpoint, refines their owner grids, and rounds integer mixture masses
down. That chosen mixture is then fixed throughout the interval. The verifier
separately establishes the bounds over the whole interval; the midpoint used to
choose a candidate is not a sampled-time correctness test.

`parameters.json` supplies the common parameters. `proposals/` contains the
165 candidate lists used as inputs. `SEED_SHA256SUMS.json` checks these inputs;
`FINAL_PROPOSAL_SHA256SUMS.json` records the expected 5,120 output checksums.

## Scope

This reproduces the final certificate generation from frozen candidate inputs.
It does not reproduce the earlier exploratory optimization that selected the
threshold and candidate plane lists, and makes no claim of optimality or of
restarting the entire historical search from no supplied parameters.
The verifier and the mathematical validity of the final certificate do not
depend on how those candidate inputs were discovered.
