#!/usr/bin/env python3
"""Check the frozen IID certificate, or reproduce its complete backward run.

  python3 reproduce.py --check-only
  python3 reproduce.py --full --threads 8

The default is --check-only.  The evaluation arithmetic is unchanged.
Check-only verifies release integrity, exact input assumptions, and
the exact-rational coefficient validator; it does not assert that
the complete backward computation has been rerun.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time


MANIFEST_SHA256 = "2405fa74f9309e6bb7e5803da99dbb59071c327221488acea2ba30dc1bef8432"
SOURCE_SHA256 = "cbc968020a858132d77278eb8fa97b4f95d481e8ccc5c681dcf6b5f5310de9cf"
F = Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def digest(path: Path) -> str:
    result = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def read_json(path: Path):
    return json.loads(path.read_text())


def integrity(bundle: Path) -> dict:
    manifest_path = bundle / "SHA256SUMS.json"
    require(digest(manifest_path) == MANIFEST_SHA256,
            "The frozen SHA256SUMS.json differs from this release.")
    manifest = read_json(manifest_path)
    require(isinstance(manifest, dict) and len(manifest) == 41,
            "Unexpected frozen manifest structure.")
    for relative, expected in manifest.items():
        path = (bundle / relative).resolve()
        require(path.is_relative_to(bundle), "Manifest path escapes bundle: " + relative)
        require(path.is_file(), "Missing frozen file: " + relative)
        require(digest(path) == expected, "Frozen file hash mismatch: " + relative)
    require(digest(bundle / "verify.cpp") == SOURCE_SHA256,
            "Unexpected evaluator source.")
    return {"status": "PASS", "manifest_sha256": MANIFEST_SHA256,
            "files_checked": len(manifest), "source_sha256": SOURCE_SHA256}


def exact_assumptions(bundle: Path) -> dict:
    data = bundle / "data"
    grid = read_json(data / "grid1200.json")
    runtime_grid = list(map(int, (data / "grids/grid1200.txt").read_text().split()))
    require(runtime_grid == [grid["m"]] + grid["numerators"],
            "Runtime text grid differs from the grid validated by validate_inputs.py.")
    require(grid["m"] == 1200 and grid["denominator"] == 2**30,
            "Unexpected state-grid dimensions.")
    knots = grid["numerators"]
    require(len(knots) == 1201 and knots[0] == 0 and knots[-1] == 2**30
            and all(a < b for a, b in zip(knots, knots[1:])),
            "Invalid state-grid knots.")

    threshold = read_json(data / "threshold40.json")
    require(threshold["steps"] == 40, "Unexpected threshold dimension.")
    raw = [F(n, threshold["a_denominator"]) for n in threshold["a_numerators"]]
    require(len(raw) == 40, "Expected 40 threshold parameters.")
    average = sum(raw) / 40
    slopes = [a / (1 + average) for a in raw]
    height = average / (1 + average)
    cap = max(slopes[j] + sum(slopes[j:]) / 40 - 1 for j in range(40))
    require(F(threshold["H"]) == height and F(threshold["cap"]) == cap,
            "Threshold height or cap does not match exact reconstruction.")
    direct_threshold = read_json(data / "threshold_slopes.json")
    direct_slopes = [F(row["numerator"], row["denominator"]) for row in direct_threshold["b"]]
    direct_height = direct_threshold["initial_height"]
    require(direct_threshold["steps"] == 40 and direct_slopes == slopes
            and F(direct_height["numerator"], direct_height["denominator"]) == height,
            "The direct rational slopes/height differ from the evaluator's threshold reconstruction.")
    require(all(F(1, 4) < b < 2 for b in slopes)
            and all(a <= b for a, b in zip(slopes, slopes[1:]))
            and F(85, 100) < cap < F(9, 10) and 0 < height <= 1,
            "Threshold assumptions in the arithmetic proof fail.")

    rate_data = read_json(data / "rates_refined_Q131093.json")
    rates = [F(n, rate_data["denominator"]) for n in rate_data["numerators"]]
    require(rate_data["Q"] == len(rates) == 131093 and rates[0] == 0
            and all(a < b for a, b in zip(rates, rates[1:]))
            and rates[-1] == 8192, "Unexpected rate-grid data.")
    require(F(64) in rates and all(r == 64 or abs(r - 64) > F(1, 10**6) for r in rates),
            "The rate-64 guard separation assumption fails.")

    source = (bundle / "verify.cpp").read_text()
    match = re.search(r"ln2_lower\s*=\s*(0x[0-9a-f.]+p[+-]?\d+)", source)
    require(match is not None, "Cannot find the frozen log(2) lower constant.")
    ln2_constant = F.from_float(float.fromhex(match.group(1)))
    w = F(1, 3)
    ln2_lower = 2 * sum(w**(2*j + 1) / (2*j + 1) for j in range(40))
    ln2_tail = 2 * w**81 / (81 * (1 - w*w))
    require(ln2_constant + F(1, 2**53) < ln2_lower,
            "The hexadecimal log(2) constant lacks the proof's 2^-53 rational safety margin.")

    # This checks the exact recurrence stated in the accompanying proof.
    # It does not replace that proof's derivation of the recurrence.
    eps, unit = F(1, 2**32), F(1, 2**40)
    multiplier = 1 + 2*eps + 2*unit
    additive = F(3, 640) + 4*eps + 8*unit + F(1, 2**28) + F(1, 2**16)
    denominator = 1 - 640 * (multiplier - 1)
    require(denominator > 0, "The table-magnitude upper-bound denominator is nonpositive.")
    magnitude = 640 * additive / denominator
    require(magnitude == F(1846693964, 613566573) and magnitude < F(301, 100),
            "The accompanying proof's conservative table-magnitude bound fails.")
    require(4 * 2**40 * 2**60 < 2**127, "The signed-128 barycentric bound fails.")
    return {
        "status": "PASS", "runtime_grid_matches_validated_json": True,
        "direct_rational_threshold_matches_reconstruction": True,
        "threshold_height": str(height), "cap": str(cap),
        "slope_count": 40, "state_subdivisions": 1200, "rate_intervals": len(rates),
        "last_finite_rate": 8192, "rate_64_guard_separated": True,
        "ln2_constant_below_rational_lower_enclosure": True,
        "ln2_constant_rational_safety_margin_exceeds_2_pow_minus_53": True,
        "ln2_lower_gap_display": float(ln2_lower - ln2_constant),
        "ln2_enclosure_width_display": float(ln2_tail),
        "table_recurrence_bound_below_3_01": True,
        "table_recurrence_bound_exact": str(magnitude),
        "table_recurrence_bound_display": float(magnitude),
        "signed128_barycentric_range_check": True,
        "scope": "Exact checks of stated input assumptions; not a replacement for the mathematical error/pruning proofs."
    }


def main() -> None:
    root = Path(__file__).resolve().parent
    default_bundle = root / "certificate" if (root / "certificate").is_dir() else root
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check-only", action="store_true", help="integrity and exact-input checks only (default)")
    mode.add_argument("--full", action="store_true", help="also construct all 640 backward layers")
    parser.add_argument("--bundle", type=Path, default=default_bundle)
    parser.add_argument("--threads", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--compiler", default="g++")
    parser.add_argument("--output", type=Path, default=root / "runs/reproduction")
    parser.add_argument("--resume", action="store_true", help="resume an interrupted --full run with matching inputs")
    args = parser.parse_args()
    require(args.threads > 0, "--threads must be positive.")
    require(not args.resume or args.full, "--resume requires --full.")
    require(sys.flags.optimize == 0 and not os.environ.get("PYTHONOPTIMIZE"),
            "Assertions must remain enabled: do not use python -O or PYTHONOPTIMIZE.")
    bundle, output = args.bundle.resolve(), args.output.resolve()
    require(bundle.is_dir(), "Certificate directory does not exist.")
    started = time.monotonic()
    report = {"started_utc": datetime.now(timezone.utc).isoformat(),
              "mode": "full" if args.full else "check-only"}
    report["integrity"] = integrity(bundle)
    print("PASS: all 41 certificate file hashes.", flush=True)
    report["exact_assumptions"] = exact_assumptions(bundle)
    print("PASS: exact grid, threshold, rate-guard, logarithm-constant, and range checks.", flush=True)
    result = subprocess.run([sys.executable, str(bundle / "validate_inputs.py")],
                            check=True, text=True, capture_output=True)
    report["exact_coefficient_validation"] = json.loads(result.stdout)
    require(report["exact_coefficient_validation"].get("validation", "").startswith("PASS"),
            "Exact coefficient validator did not report PASS.")
    print(result.stdout, end="", flush=True)
    report["preflight_seconds"] = time.monotonic() - started
    report["complete_backward_run_replayed"] = False
    output.parent.mkdir(parents=True, exist_ok=True)
    preflight_path = output.with_name(output.name + ".preflight.json")
    preflight_path.write_text(json.dumps(report, indent=2) + "\n")
    if not args.full:
        print("CHECK-ONLY PASS. No backward layers were replayed.", flush=True)
        print("Report: " + str(preflight_path), flush=True)
        return

    command = [sys.executable, str(bundle / "run.py"), "--threads", str(args.threads),
               "--compiler", args.compiler, "--output", str(output), "--skip-input-validation"]
    if args.resume:
        command.append("--resume")
    print("Starting the full evaluator. Progress is written to "
          + str(output / "replay.stderr"), flush=True)
    subprocess.run(command, check=True)
    final, expected = read_json(output / "RESULT.json"), read_json(bundle / "RESULTS.json")
    require(final.get("last_k") == 0, "The evaluator did not reach layer zero.")
    require((int(final["numerator"]), int(final["denominator"])) ==
            (int(expected["gamma_numerator"]), int(expected["gamma_denominator"])),
            "The completed objective differs from the frozen exact result.")
    require(final["m"] == 1200 and final["N"] == 640 and final["Q"] == 131093
            and final["slack_bits"] == 28 and final["potential_bits"] == 16,
            "The completed evaluator used unexpected numerical parameters.")
    require(10000 * int(final["numerator"]) >= 7419 * int(final["denominator"]),
            "The exact target comparison failed.")
    checkpoint_checks = {}
    for suffix in ("table", "next", "hints"):
        name = "layer_0." + suffix
        expected_digest = digest(bundle / "checkpoint_0" / name)
        actual_digest = digest(output / name)
        require(actual_digest == expected_digest,
                "The complete layer-zero checkpoint differs from the frozen witness: " + name)
        checkpoint_checks[name] = actual_digest
    report.update({"complete_backward_run_replayed": True,
                   "exact_objective_matches_frozen_result": True,
                   "layer_zero_checkpoint_matches_frozen": True,
                   "layer_zero_checkpoint_sha256": checkpoint_checks,
                   "result": final, "finished_utc": datetime.now(timezone.utc).isoformat(),
                   "total_seconds": time.monotonic() - started})
    (output / "REPRODUCTION_CHECK.json").write_text(json.dumps(report, indent=2) + "\n")
    print("FULL REPRODUCTION PASS: all 640 layers completed; exact objective matches RESULTS.json.")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print("REPRODUCTION FAILED: " + str(error), file=sys.stderr)
        raise SystemExit(1)
