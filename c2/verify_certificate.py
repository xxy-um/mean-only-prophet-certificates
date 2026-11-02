#!/usr/bin/env python3
"""Verify the finite dual certificate in Appendix C.2 using exact integers.

Only the Python standard library is required.  No optimization is performed.
The input's ``cutoffs`` field is auxiliary metadata retained from the source
artifact; it is not used in verification. The mathematical input consists only of N, the rate
and weight denominators, profiles, and weights.

An interval (lo, hi) means [lo/B, hi/B], where B is the scale in the proof
note.  Every operation affecting a certificate bound uses integers only.
Floating-point values occur only in elapsed-time measurements.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from fractions import Fraction
import gzip
import hashlib
import json
from math import factorial
from pathlib import Path
import platform
import sys
import time


B = 1 << 160
CLAIMED_CEILING_NUMERATOR = 741966695100
CLAIMED_CEILING_DENOMINATOR = 10**12
TARGET_NUMERATOR = 741967
TARGET_DENOMINATOR = 10**6
ROOT = Path(__file__).resolve().parent


class CertificateError(ValueError):
    """An input or a required certificate inequality failed verification."""


def ceil_div(numerator: int, denominator: int) -> int:
    """Mathematical ceiling, including for negative numerators."""
    if denominator <= 0:
        raise CertificateError("A division denominator must be positive.")
    return -((-numerator) // denominator)


def add(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
    return a[0] + b[0], a[1] + b[1]


def sub(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
    return a[0] - b[1], a[1] - b[0]


def mul(a: tuple[int, int], b: tuple[int, int]) -> tuple[int, int]:
    products = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return min(products) // B, ceil_div(max(products), B)


def scale(a: tuple[int, int], numerator: int,
          denominator: int = 1) -> tuple[int, int]:
    """Multiply an interval by the exact rational numerator/denominator."""
    if denominator <= 0:
        raise CertificateError("A scaling denominator must be positive.")
    lo, hi = a if numerator >= 0 else (a[1], a[0])
    return (lo * numerator // denominator,
            ceil_div(hi * numerator, denominator))


def expneg(numerator: int, denominator: int = 1) -> tuple[int, int]:
    """Enclose exp(-numerator/denominator) by the Appendix C.2 procedure."""
    if type(numerator) is not int or type(denominator) is not int:
        raise CertificateError("Exponential arguments must be exact integers.")
    if numerator < 0 or denominator <= 0:
        raise CertificateError("The exponential argument must be nonnegative.")
    if numerator == 0:
        return B, B
    if numerator >= 160 * denominator:
        return 0, 1  # e > 2 implies exp(-160) < 2**(-160).

    squarings = 0
    while numerator > denominator:
        denominator *= 2
        squarings += 1

    # Degree 50, hence 51 Taylor terms including the constant term.
    lo = hi = lower_sum = upper_sum = B
    for m in range(1, 51):
        lo = lo * numerator // (denominator * m)
        hi = ceil_div(hi * numerator, denominator * m)
        lower_sum += lo
        upper_sum += hi
    upper_sum += ceil_div(3 * B, factorial(51))

    # Reciprocal of the positive enclosure of exp(z).
    lo, hi = B * B // upper_sum, ceil_div(B * B, lower_sum)
    for _ in range(squarings):
        lo, hi = lo * lo // B, ceil_div(hi * hi, B)
    return lo, hi


def _object_without_duplicate_keys(pairs: list[tuple[str, object]]) -> dict:
    result = {}
    for key, value in pairs:
        if key in result:
            raise CertificateError(f"Duplicate JSON object key: {key}")
        result[key] = value
    return result


def load_certificate(path: Path) -> tuple[dict, str]:
    """Load exact JSON integers and hash the original (possibly gzip) bytes."""
    raw = path.read_bytes()
    payload = gzip.decompress(raw) if path.suffix == ".gz" else raw
    data = json.loads(payload.decode("utf-8"),
                      object_pairs_hook=_object_without_duplicate_keys)
    if not isinstance(data, dict):
        raise CertificateError("The certificate must be a JSON object.")
    return data, hashlib.sha256(raw).hexdigest()


def validate_certificate(data: dict, fixed_instance: bool = True) -> dict:
    """Check the exact input; fixed_instance also checks the released dimensions."""
    if not isinstance(data, dict):
        raise CertificateError("The certificate must be an object.")
    for key in ("N", "rate_denominator", "weight_denominator"):
        if type(data.get(key)) is not int or data[key] <= 0:
            raise CertificateError(f"{key} must be a positive integer.")
    N = data["N"]
    rate_denominator = data["rate_denominator"]
    weight_denominator = data["weight_denominator"]
    weights, profiles = data.get("weights"), data.get("profiles")
    if not isinstance(weights, list) or not isinstance(profiles, list):
        raise CertificateError("weights and profiles must be lists.")
    if not profiles or len(profiles) != len(weights):
        raise CertificateError("There must be one weight for every rate profile.")
    if any(type(w) is not int or w < 0 for w in weights):
        raise CertificateError("All weight numerators must be nonnegative integers.")
    if sum(weights) != weight_denominator:
        raise CertificateError("The weights do not sum exactly to one.")

    if fixed_instance:
        expected = {"N": 1024, "rate_denominator": 10**8,
                    "weight_denominator": 10**12}
        for key, value in expected.items():
            if data[key] != value:
                raise CertificateError(f"The released instance requires {key}={value}.")
        if len(profiles) != 1014:
            raise CertificateError("The released instance requires 1014 rate profiles.")

    maximum_rate_numerator = 0
    for row, rates in enumerate(profiles):
        if not isinstance(rates, list) or len(rates) != N:
            raise CertificateError(f"Profile {row} must have exactly {N} entries.")
        if any(type(r) is not int or r < 0 for r in rates):
            raise CertificateError(f"Profile {row} has a noninteger or negative rate.")
        if any(left > right for left, right in zip(rates, rates[1:])):
            raise CertificateError(f"Profile {row} is not nondecreasing.")
        maximum_rate_numerator = max(maximum_rate_numerator, rates[-1])
    # This also establishes the finite-n cutoff reported in the manuscript.
    if maximum_rate_numerator > 10**9 * rate_denominator:
        raise CertificateError("A rate exceeds the released certificate's cap of 10^9.")

    return {
        "N": N,
        "profiles": len(profiles),
        "rate_denominator": rate_denominator,
        "weight_denominator": weight_denominator,
        "weight_numerator_sum": sum(weights),
        "max_rate_numerator": maximum_rate_numerator,
        "max_rate": str(Fraction(maximum_rate_numerator, rate_denominator)),
        "all_rates_nonnegative_and_nondecreasing": True,
        "all_weights_nonnegative": True,
        "weights_sum_exactly_to_one": True,
        "cutoffs_used": False,
    }


def evaluate_coefficients(data: dict, progress=None) -> tuple[tuple[int, int], list]:
    """Enclose the weighted A_q(0) and k_{q,j} averages.

    validate_certificate must be called first. The operation order matches
    the original verifier. For q_j > 0 the cell integral is evaluated as
    S_q(j/N) * (1-exp(-q_j/N))/q_j, algebraically the paper's formula for a_j.
    """
    N = data["N"]
    rate_denominator = data["rate_denominator"]
    weight_denominator = data["weight_denominator"]
    weights, profiles = data["weights"], data["profiles"]
    mean_A = (0, 0)
    mean_k = [(0, 0) for _ in range(N)]
    cache = {}

    for row, (weight, rates) in enumerate(zip(weights, profiles)):
        survival = (B, B)  # S_q(0).
        cell_integrals = []  # a_j.
        cell_end_survival = []  # S_q((j+1)/N).
        exp_rate = []  # exp(-q_j).
        for rate in rates:
            if rate not in cache:
                exp_cell = expneg(rate, rate_denominator * N)
                exp_full = expneg(rate, rate_denominator)
                if rate:
                    integral_factor = scale(sub((B, B), exp_cell),
                                            rate_denominator, rate)
                else:
                    integral_factor = (B // N, ceil_div(B, N))
                cache[rate] = exp_cell, exp_full, integral_factor
            exp_cell, exp_full, integral_factor = cache[rate]
            cell_integrals.append(mul(survival, integral_factor))
            survival = mul(survival, exp_cell)
            cell_end_survival.append(survival)
            exp_rate.append(exp_full)

        A = (sum(a[0] for a in cell_integrals),
             sum(a[1] for a in cell_integrals))  # A_q(0).
        tail = (0, 0)  # Sum of a_l over l > j in the backward loop.
        mean_A = add(mean_A, scale(A, weight, weight_denominator))
        for j in range(N - 1, -1, -1):
            k = sub(sub(cell_end_survival[j],
                        scale(tail, rates[j], rate_denominator)),
                    mul(A, exp_rate[j]))
            mean_k[j] = add(mean_k[j], scale(k, weight, weight_denominator))
            tail = add(tail, cell_integrals[j])
        if progress is not None and (row % max(1, len(profiles) // 8) == 0
                                    or row + 1 == len(profiles)):
            progress(f"Evaluated {row + 1}/{len(profiles)} rate profiles.")
    return mean_A, mean_k


def certify_bounds(mean_A: tuple[int, int], mean_k: list) -> tuple[dict, dict]:
    """Check the claimed decimal ceiling and U using exact integer comparisons."""
    intervals = [mean_A] + [(B - hi, B - lo) for lo, hi in mean_k]
    if any(lo > hi for lo, hi in intervals):
        raise CertificateError("An invalid coefficient interval was computed.")
    upper_numerator = max(hi for _, hi in intervals)
    at_most_claimed = (upper_numerator * CLAIMED_CEILING_DENOMINATOR
                       <= CLAIMED_CEILING_NUMERATOR * B)
    claimed_below_target = (CLAIMED_CEILING_NUMERATOR * TARGET_DENOMINATOR
                            < TARGET_NUMERATOR * CLAIMED_CEILING_DENOMINATOR)
    below_target = upper_numerator * TARGET_DENOMINATOR < TARGET_NUMERATOR * B
    coefficient_checks = [hi * TARGET_DENOMINATOR <= TARGET_NUMERATOR * B
                          for _, hi in intervals]
    if not at_most_claimed:
        raise CertificateError("The computed bound exceeds 0.741966695100.")
    if not claimed_below_target or not below_target or not all(coefficient_checks):
        raise CertificateError("The coefficient inequalities for U=0.741967 failed.")

    rounded = ceil_div(upper_numerator * 10**12, B)
    quotient, remainder = divmod(rounded, 10**12)
    result = {
        "upper_numerator": str(upper_numerator),
        "upper_denominator": str(B),
        "upper_rounded_up_12_decimals": rounded,
        "upper_rounded_up_decimal": f"{quotient}.{remainder:012d}",
        "claimed_ceiling_numerator": CLAIMED_CEILING_NUMERATOR,
        "claimed_ceiling_denominator": CLAIMED_CEILING_DENOMINATOR,
        "target_U_numerator": TARGET_NUMERATOR,
        "target_U_denominator": TARGET_DENOMINATOR,
        "bound_at_most_claimed_ceiling": at_most_claimed,
        "claimed_ceiling_strictly_below_target": claimed_below_target,
        "computed_bound_strictly_below_target": below_target,
        "all_coefficient_inequalities_verified": all(coefficient_checks),
        "coefficient_count": len(intervals),
    }
    coefficients = []
    for index, (lo, hi) in enumerate(intervals):
        entry = {"coefficient": "A_q(0)" if index == 0 else "1-k_{q,j}",
                 "lower_numerator": str(lo), "upper_numerator": str(hi),
                 "at_most_U": coefficient_checks[index]}
        if index:
            entry["j"] = index - 1
        coefficients.append(entry)
    bounds = {"status": "PASS", "scale": str(B),
              "interpretation": "Weighted coefficient averages; divide endpoints by scale.",
              "coefficients": coefficients}
    return result, bounds


def _write_json(path: Path, content: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(content, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path,
                        default=ROOT / "data" / "upper_1024.json.gz",
                        help="Exact JSON or gzip-compressed JSON certificate.")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "runs" / "latest",
                        help="Directory for verification.json, coefficient_bounds.json and verification.log.")
    args = parser.parse_args(argv)
    output_dir = args.output_dir.resolve()
    data_path = args.data.resolve()
    source_path = Path(__file__).resolve()
    destinations = [output_dir / name for name in
                    ("verification.json", "coefficient_bounds.json", "verification.log")]
    if data_path in destinations or source_path in destinations:
        parser.error("The output files must not overwrite the input or verifier.")

    started_utc = datetime.now(timezone.utc).isoformat()
    start = time.perf_counter()
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        with (output_dir / "verification.log").open("w", encoding="utf-8") as log:
            def report(message):
                print(message, flush=True)
                log.write(message + "\n")
                log.flush()

            report("Appendix C.2 finite dual certificate verification")
            report(f"Started UTC: {started_utc}")
            report("All certificate arithmetic uses exact integers and outward intervals.")
            # An interrupted rerun must not leave previous PASS reports looking
            # like the output of the new run. Only these report files are replaced.
            pending = {"status": "RUNNING", "started_utc": started_utc}
            _write_json(output_dir / "verification.json", pending)
            _write_json(output_dir / "coefficient_bounds.json", pending)
            try:
                data, data_hash = load_certificate(data_path)
                metadata = validate_certificate(data)
                metadata.update({"file": data_path.name, "sha256": data_hash})
                report(f"Input checks passed: N={metadata['N']}, profiles={metadata['profiles']}.")
                report(f"Exact weight numerator sum: {metadata['weight_numerator_sum']}.")
                report(f"Computed maximum rate: {metadata['max_rate']}.")
                mean_A, mean_k = evaluate_coefficients(data, report)
                result, bounds = certify_bounds(mean_A, mean_k)
                elapsed = time.perf_counter() - start
                verification = {
                    "status": "PASS", "certificate": metadata,
                    "arithmetic": {"interval_scale": str(B), "precision_bits": 160,
                                   "exponential_taylor_degree": 50,
                                   "certificate_comparisons": "exact integer arithmetic"},
                    "result": result,
                    "execution": {
                        "started_utc": started_utc,
                        "completed_utc": datetime.now(timezone.utc).isoformat(),
                        "elapsed_seconds": elapsed, "python": sys.version,
                        "platform": platform.platform(),
                        "verifier_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
                        "command_arguments": sys.argv[1:] if argv is None else argv,
                    },
                }
                _write_json(output_dir / "coefficient_bounds.json", bounds)
                _write_json(output_dir / "verification.json", verification)
                report("Computed upper numerator: " + result["upper_numerator"])
                report("Interval scale / upper denominator: " + str(B))
                report("Upward-rounded bound: " + result["upper_rounded_up_decimal"])
                report("Exact check: U_hat <= 0.741966695100 < U = 0.741967.")
                report(f"All {result['coefficient_count']} coefficient inequalities verified.")
                report(f"Elapsed seconds (timing only): {elapsed:.6f}")
                report("PASS")
                return 0
            except (CertificateError, OSError, ValueError, KeyError, TypeError) as error:
                failure = {"status": "FAIL", "error": str(error),
                           "started_utc": started_utc,
                           "completed_utc": datetime.now(timezone.utc).isoformat()}
                # Replace previous reports for these filenames, so a failed rerun
                # cannot leave an old successful coefficient report in their place.
                _write_json(output_dir / "verification.json", failure)
                _write_json(output_dir / "coefficient_bounds.json", failure)
                report("FAIL: " + str(error))
                return 1
    except OSError as error:
        print("FAIL: could not write verification output: " + str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
