"""Small exact-arithmetic and rejection tests; no full certificate run."""

from copy import deepcopy
from fractions import Fraction
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SOURCE = Path(__file__).resolve().parents[1] / "verify_certificate.py"
SPEC = importlib.util.spec_from_file_location("verify_certificate", SOURCE)
v = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(v)


def reference_expneg(p, q):
    """Independent alternating-series enclosure, with exact rational terms."""
    z = Fraction(p, q)
    h = 0
    while z > 1:
        z /= 2
        h += 1
    term = total = Fraction(1)
    for m in range(1, 201):
        term *= -z / m
        total += term
    upper = total  # Even truncation degree.
    lower = total + term * (-z / 201)  # Odd truncation degree.
    for _ in range(h):
        lower *= lower
        upper *= upper
    return lower, upper


class ArithmeticTests(unittest.TestCase):
    def test_signed_ceil_div(self):
        for numerator in (-11, -10, -9, -1, 0, 1, 9, 10, 11):
            result = v.ceil_div(numerator, 5)
            exact = Fraction(numerator, 5)
            self.assertLess(result - 1, exact)
            self.assertLessEqual(exact, result)

    def test_signed_product_contains_corner_extrema(self):
        B = v.B
        cases = [((-3*B-7, 2*B+11), (-5*B+3, -B+1)),
                 ((B//3, B//2+1), (B//7, 3*B//5)),
                 ((-B//3, B//7), (-B//11, B//13))]
        for a, b in cases:
            lo, hi = v.mul(a, b)
            exact = [Fraction(x*y, B*B) for x in a for y in b]
            self.assertLessEqual(Fraction(lo, B), min(exact))
            self.assertLessEqual(max(exact), Fraction(hi, B))
            self.assertLess(min(exact)-Fraction(lo, B), Fraction(1, B))
            self.assertLess(Fraction(hi, B)-max(exact), Fraction(1, B))

    def test_subtraction_and_negative_scaling(self):
        B = v.B
        a, b = (-B-1, 2*B+3), (-3*B+2, B+4)
        self.assertEqual(v.sub(a, b), (a[0]-b[1], a[1]-b[0]))
        for p in (-7, 0, 7):
            lo, hi = v.scale(a, p, 3)
            exact = [Fraction(x*p, 3*B) for x in a]
            self.assertLessEqual(Fraction(lo, B), min(exact))
            self.assertLessEqual(max(exact), Fraction(hi, B))

    def test_exp_against_independent_alternating_series(self):
        for p, q in ((1, 1000), (1, 3), (1, 1), (3, 2), (8, 1)):
            lo, hi = v.expneg(p, q)
            lower, upper = reference_expneg(p, q)
            self.assertLessEqual(Fraction(lo, v.B), lower)
            self.assertLessEqual(upper, Fraction(hi, v.B))

    def test_exp_zero_and_cutoff(self):
        self.assertEqual(v.expneg(0), (v.B, v.B))
        self.assertEqual(v.expneg(160), (0, 1))
        self.assertEqual(v.expneg(160000, 1000), (0, 1))
        self.assertEqual(v.expneg(10**9), (0, 1))
        with self.assertRaises(v.CertificateError):
            v.expneg(-1)


class CertificateTests(unittest.TestCase):
    def setUp(self):
        self.data = {"N": 2, "rate_denominator": 10,
                     "weight_denominator": 10,
                     "weights": [4, 6], "profiles": [[0, 1], [2, 3]],
                     "cutoffs": [0, 0]}

    def test_structural_rejections(self):
        changes = [("N", True), ("rate_denominator", 0),
                   ("weights", [4, 5]), ("weights", [True, 9]),
                   ("weights", [-1, 11]), ("profiles", [[1, 0], [2, 3]]),
                   ("profiles", [[0, 1.0], [2, 3]]),
                   ("profiles", [[0], [2, 3]]),
                   ("profiles", [[0, 1], [2, 10**11]])]
        for key, value in changes:
            data = deepcopy(self.data)
            data[key] = value
            with self.subTest(key=key, value=value):
                with self.assertRaises(v.CertificateError):
                    v.validate_certificate(data, fixed_instance=False)
        with self.assertRaises(v.CertificateError):
            v.validate_certificate(self.data)

    def test_metadata_does_not_change_coefficients(self):
        v.validate_certificate(self.data, fixed_instance=False)
        expected = v.evaluate_coefficients(self.data)
        changed = deepcopy(self.data)
        changed["cutoffs"] = [123, 456]
        self.assertEqual(v.evaluate_coefficients(changed), expected)

    def test_zero_rate_branch(self):
        data = {"N": 2, "rate_denominator": 1, "weight_denominator": 1,
                "profiles": [[0, 0]], "weights": [1]}
        v.validate_certificate(data, fixed_instance=False)
        mean_A, mean_k = v.evaluate_coefficients(data)
        self.assertEqual(mean_A, (v.B, v.B))
        self.assertEqual(mean_k, [(0, 0), (0, 0)])
        with self.assertRaises(v.CertificateError):
            v.certify_bounds(mean_A, mean_k)

    def test_endpoint_direction_and_target_rejection(self):
        # The 1-k coefficient uses B minus the LOWER endpoint of the k average.
        low_k = v.B // 2
        result, bounds = v.certify_bounds((0, v.B//4), [(low_k, v.B)])
        self.assertEqual(int(result["upper_numerator"]), v.B-low_k)
        self.assertEqual(bounds["coefficients"][1]["upper_numerator"], str(v.B-low_k))
        with self.assertRaises(v.CertificateError):
            v.certify_bounds((0, v.B), [(0, v.B)])

    def test_duplicate_json_keys_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "bad.json"
            path.write_text('{"N": 2, "N": 3}', encoding="utf-8")
            with self.assertRaises(v.CertificateError):
                v.load_certificate(path)


if __name__ == "__main__":
    unittest.main()
