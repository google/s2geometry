"""Tests for S1Angle pybind11 bindings."""

import math
import unittest
import s2geometry_pybind as s2


class TestS1Angle(unittest.TestCase):
    """Test cases for S1Angle bindings."""

    # Constructors

    def test_default_constructor(self):
        a = s2.S1Angle()
        self.assertEqual(a.radians(), 0.0)

    def test_constructor_from_two_points(self):
        x = s2.S2Point(1.0, 0.0, 0.0)
        y = s2.S2Point(0.0, 1.0, 0.0)
        a = s2.S1Angle(x, y)
        self.assertAlmostEqual(a.radians(), math.pi / 2)

    def test_constructor_from_same_point(self):
        p = s2.S2Point(1.0, 0.0, 0.0)
        a = s2.S1Angle(p, p)
        self.assertEqual(a.radians(), 0.0)

    def test_constructor_from_antipodal_points(self):
        x = s2.S2Point(1.0, 0.0, 0.0)
        y = s2.S2Point(-1.0, 0.0, 0.0)
        a = s2.S1Angle(x, y)
        self.assertAlmostEqual(a.radians(), math.pi)

    # Factory methods

    def test_from_radians(self):
        a = s2.S1Angle.from_radians(math.pi / 4)
        self.assertAlmostEqual(a.radians(), math.pi / 4)

    def test_from_degrees(self):
        a = s2.S1Angle.from_degrees(45.0)
        self.assertAlmostEqual(a.degrees(), 45.0)

    def test_from_e5(self):
        a = s2.S1Angle.from_e5(4500000)
        self.assertAlmostEqual(a.degrees(), 45.0)

    def test_from_e6(self):
        a = s2.S1Angle.from_e6(45000000)
        self.assertAlmostEqual(a.degrees(), 45.0)

    def test_from_e7(self):
        a = s2.S1Angle.from_e7(450000000)
        self.assertAlmostEqual(a.degrees(), 45.0)

    def test_zero(self):
        a = s2.S1Angle.zero()
        self.assertEqual(a.radians(), 0.0)

    def test_infinity(self):
        a = s2.S1Angle.infinity()
        self.assertTrue(math.isinf(a.radians()))
        self.assertTrue(a.radians() > 0)

    # Geometric operations

    def test_radians(self):
        a = s2.S1Angle.from_radians(1.5)
        self.assertEqual(a.radians(), 1.5)

    def test_degrees(self):
        a = s2.S1Angle.from_degrees(180.0)
        self.assertAlmostEqual(a.radians(), math.pi)
        self.assertAlmostEqual(a.degrees(), 180.0)

    def test_degrees_radians_conversion(self):
        # Exact conversions per the C++ header documentation.
        self.assertEqual(s2.S1Angle.from_degrees(180.0).radians(), math.pi)
        for k in range(9):
            self.assertEqual(
                s2.S1Angle.from_degrees(45 * k).radians(),
                k * math.pi / 4)

    def test_e5(self):
        a = s2.S1Angle.from_degrees(45.0)
        self.assertEqual(a.e5(), 4500000)

    def test_e6(self):
        a = s2.S1Angle.from_degrees(45.0)
        self.assertEqual(a.e6(), 45000000)

    def test_e7(self):
        a = s2.S1Angle.from_degrees(45.0)
        self.assertEqual(a.e7(), 450000000)

    def test_e5_e6_e7_consistency(self):
        # Degrees(n) == E6(1000000 * n) == E7(10000000 * n) for integer n.
        for n in [0, 1, 45, 90, 180, -180, -90]:
            self.assertEqual(s2.S1Angle.from_degrees(n).radians(),
                             s2.S1Angle.from_e6(1000000 * n).radians())
            self.assertEqual(s2.S1Angle.from_degrees(n).radians(),
                             s2.S1Angle.from_e7(10000000 * n).radians())

    def test_e5_not_normalized_raises(self):
        a = s2.S1Angle.from_degrees(270.0)
        with self.assertRaises(ValueError) as cm:
            a.e5()
        self.assertEqual(str(cm.exception),
                         "Angle 270.000000 degrees is not in the "
                         "normalized range (-180, 180]")

    def test_e6_not_normalized_raises(self):
        a = s2.S1Angle.from_degrees(270.0)
        with self.assertRaises(ValueError) as cm:
            a.e6()
        self.assertEqual(str(cm.exception),
                         "Angle 270.000000 degrees is not in the "
                         "normalized range (-180, 180]")

    def test_e7_not_normalized_raises(self):
        a = s2.S1Angle.from_degrees(270.0)
        with self.assertRaises(ValueError) as cm:
            a.e7()
        self.assertEqual(str(cm.exception),
                         "Angle 270.000000 degrees is not in the "
                         "normalized range (-180, 180]")

    def test_e6_minus_180_not_normalized_raises(self):
        # -180 degrees is outside the normalized range (-180, 180].
        a = s2.S1Angle.from_degrees(-180.0)
        with self.assertRaises(ValueError) as cm:
            a.e6()
        self.assertEqual(str(cm.exception),
                         "Angle -180.000000 degrees is not in the "
                         "normalized range (-180, 180]")

    def test_e6_plus_180_is_valid(self):
        # +180 degrees is inside the normalized range (-180, 180].
        a = s2.S1Angle.from_degrees(180.0)
        self.assertEqual(a.e6(), 180000000)

    def test_abs(self):
        pos = s2.S1Angle.from_degrees(45.0)
        neg = s2.S1Angle.from_degrees(-45.0)
        self.assertAlmostEqual(abs(pos).degrees(), 45.0)
        self.assertAlmostEqual(abs(neg).degrees(), 45.0)

    def test_abs_zero(self):
        a = s2.S1Angle.zero()
        self.assertEqual(abs(a).radians(), 0.0)

    def test_normalized(self):
        a = s2.S1Angle.from_degrees(270.0)
        self.assertAlmostEqual(a.normalized().degrees(), -90.0)

    def test_normalized_negative(self):
        a = s2.S1Angle.from_degrees(-270.0)
        self.assertAlmostEqual(a.normalized().degrees(), 90.0)

    def test_normalized_minus_180_becomes_180(self):
        # Normalized range is (-180, 180]; -180 maps to 180.
        a = s2.S1Angle.from_degrees(-180.0)
        self.assertAlmostEqual(a.normalized().degrees(), 180.0)

    # Operators

    def test_equality(self):
        a = s2.S1Angle.from_degrees(45.0)
        b = s2.S1Angle.from_degrees(45.0)
        c = s2.S1Angle.from_degrees(90.0)
        self.assertTrue(a == b)
        self.assertTrue(a != c)
        self.assertFalse(a != b)
        self.assertFalse(a == c)

    def test_comparison(self):
        small = s2.S1Angle.from_degrees(10.0)
        large = s2.S1Angle.from_degrees(20.0)
        self.assertTrue(small < large)
        self.assertTrue(large > small)
        self.assertTrue(small <= large)
        self.assertTrue(large >= small)
        self.assertTrue(small <= small)
        self.assertTrue(small >= small)
        self.assertFalse(small > large)
        self.assertFalse(large < small)

    def test_negation(self):
        a = s2.S1Angle.from_degrees(45.0)
        neg = -a
        self.assertAlmostEqual(neg.degrees(), -45.0)

    def test_addition(self):
        a = s2.S1Angle.from_degrees(30.0)
        b = s2.S1Angle.from_degrees(60.0)
        result = a + b
        self.assertAlmostEqual(result.degrees(), 90.0)

    def test_subtraction(self):
        a = s2.S1Angle.from_degrees(90.0)
        b = s2.S1Angle.from_degrees(30.0)
        result = a - b
        self.assertAlmostEqual(result.degrees(), 60.0)

    def test_scalar_multiplication(self):
        a = s2.S1Angle.from_degrees(45.0)
        result1 = a * 2.0
        result2 = 2.0 * a
        self.assertAlmostEqual(result1.degrees(), 90.0)
        self.assertAlmostEqual(result2.degrees(), 90.0)

    def test_scalar_division(self):
        a = s2.S1Angle.from_degrees(90.0)
        result = a / 2.0
        self.assertAlmostEqual(result.degrees(), 45.0)

    def test_angle_division(self):
        a = s2.S1Angle.from_degrees(90.0)
        b = s2.S1Angle.from_degrees(45.0)
        ratio = a / b
        self.assertAlmostEqual(ratio, 2.0)

    def test_in_place_addition(self):
        a = s2.S1Angle.from_degrees(30.0)
        a += s2.S1Angle.from_degrees(60.0)
        self.assertAlmostEqual(a.degrees(), 90.0)

    def test_in_place_subtraction(self):
        a = s2.S1Angle.from_degrees(90.0)
        a -= s2.S1Angle.from_degrees(30.0)
        self.assertAlmostEqual(a.degrees(), 60.0)

    def test_in_place_scalar_multiplication(self):
        a = s2.S1Angle.from_degrees(45.0)
        a *= 2.0
        self.assertAlmostEqual(a.degrees(), 90.0)

    def test_in_place_scalar_division(self):
        a = s2.S1Angle.from_degrees(90.0)
        a /= 2.0
        self.assertAlmostEqual(a.degrees(), 45.0)

    # String representation

    def test_repr(self):
        a = s2.S1Angle.from_degrees(45.0)
        self.assertEqual(repr(a), "S1Angle(45.0000000)")

    def test_str(self):
        a = s2.S1Angle.from_degrees(45.0)
        self.assertEqual(str(a), "45.0000000")

    def test_repr_zero(self):
        a = s2.S1Angle.zero()
        self.assertEqual(repr(a), "S1Angle(0.0000000)")

    # Module-level functions

    def test_sin(self):
        a = s2.S1Angle.from_degrees(90.0)
        self.assertAlmostEqual(s2.sin(a), 1.0)

    def test_cos(self):
        a = s2.S1Angle.from_degrees(0.0)
        self.assertAlmostEqual(s2.cos(a), 1.0)

    def test_tan(self):
        a = s2.S1Angle.from_degrees(45.0)
        self.assertAlmostEqual(s2.tan(a), 1.0)



if __name__ == "__main__":
    unittest.main()
