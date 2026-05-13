"""Tests for S1ChordAngle pybind11 bindings."""

import math
import unittest
import s2geometry_pybind as s2


class TestS1ChordAngle(unittest.TestCase):
    """Test cases for S1ChordAngle bindings."""

    # Constructors

    def test_default_constructor(self):
        a = s2.S1ChordAngle()
        self.assertTrue(a.is_zero())

    def test_constructor_from_two_points(self):
        x = s2.S2Point(1.0, 0.0, 0.0)
        y = s2.S2Point(0.0, 1.0, 0.0)
        a = s2.S1ChordAngle(x, y)
        self.assertAlmostEqual(a.radians, math.pi / 2)

    def test_constructor_from_same_point(self):
        p = s2.S2Point(1.0, 0.0, 0.0)
        a = s2.S1ChordAngle(p, p)
        self.assertTrue(a.is_zero())

    def test_constructor_from_antipodal_points(self):
        x = s2.S2Point(1.0, 0.0, 0.0)
        y = s2.S2Point(-1.0, 0.0, 0.0)
        a = s2.S1ChordAngle(x, y)
        self.assertAlmostEqual(a.radians, math.pi)

    def test_constructor_non_unit_length_raises(self):
        bad = s2.S2Point(2.0, 0.0, 0.0)
        ok = s2.S2Point(1.0, 0.0, 0.0)
        with self.assertRaises(ValueError):
            s2.S1ChordAngle(bad, ok)
        with self.assertRaises(ValueError):
            s2.S1ChordAngle(ok, bad)

    def test_constructor_from_s1_angle(self):
        angle = s2.S1Angle.from_degrees(90.0)
        a = s2.S1ChordAngle(angle)
        self.assertAlmostEqual(a.radians, math.pi / 2)

    def test_constructor_from_s1_angle_negative(self):
        a = s2.S1ChordAngle(s2.S1Angle.from_radians(-1.0))
        self.assertTrue(a.is_negative())

    def test_constructor_from_s1_angle_infinity(self):
        a = s2.S1ChordAngle(s2.S1Angle.infinity())
        self.assertTrue(a.is_infinity())

    def test_constructor_from_s1_angle_larger_than_pi_clamps(self):
        a = s2.S1ChordAngle(s2.S1Angle.from_radians(math.pi + 1.0))
        self.assertEqual(a, s2.S1ChordAngle.straight())

    # Factory methods

    def test_from_radians(self):
        a = s2.S1ChordAngle.from_radians(math.pi / 2)
        self.assertAlmostEqual(a.radians, math.pi / 2)

    def test_from_degrees(self):
        a = s2.S1ChordAngle.from_degrees(90.0)
        self.assertAlmostEqual(a.degrees, 90.0)

    def test_zero(self):
        a = s2.S1ChordAngle.zero()
        self.assertTrue(a.is_zero())

    def test_right(self):
        a = s2.S1ChordAngle.right()
        self.assertAlmostEqual(a.radians, math.pi / 2)

    def test_straight(self):
        a = s2.S1ChordAngle.straight()
        self.assertAlmostEqual(a.radians, math.pi)

    def test_infinity(self):
        a = s2.S1ChordAngle.infinity()
        self.assertTrue(a.is_infinity())
        self.assertTrue(a.is_special())

    def test_negative(self):
        a = s2.S1ChordAngle.negative()
        self.assertTrue(a.is_negative())
        self.assertTrue(a.is_special())

    # Properties

    def test_radians(self):
        a = s2.S1ChordAngle.from_radians(1.0)
        self.assertAlmostEqual(a.radians, 1.0)

    def test_degrees(self):
        a = s2.S1ChordAngle.from_degrees(45.0)
        self.assertAlmostEqual(a.degrees, 45.0)

    def test_e5_e6_e7(self):
        a = s2.S1ChordAngle.from_degrees(45.0)
        self.assertEqual(a.e5, 4500000)
        self.assertEqual(a.e6, 45000000)
        self.assertEqual(a.e7, 450000000)

    # Predicates

    def test_is_zero(self):
        self.assertTrue(s2.S1ChordAngle.zero().is_zero())
        self.assertFalse(s2.S1ChordAngle.right().is_zero())

    def test_is_negative(self):
        self.assertTrue(s2.S1ChordAngle.negative().is_negative())
        self.assertFalse(s2.S1ChordAngle.zero().is_negative())

    def test_is_infinity(self):
        self.assertTrue(s2.S1ChordAngle.infinity().is_infinity())
        self.assertFalse(s2.S1ChordAngle.straight().is_infinity())

    def test_is_special(self):
        self.assertTrue(s2.S1ChordAngle.negative().is_special())
        self.assertTrue(s2.S1ChordAngle.infinity().is_special())
        self.assertFalse(s2.S1ChordAngle.zero().is_special())
        self.assertFalse(s2.S1ChordAngle.straight().is_special())

    def test_is_valid(self):
        self.assertTrue(s2.S1ChordAngle.zero().is_valid())
        self.assertTrue(s2.S1ChordAngle.straight().is_valid())
        self.assertTrue(s2.S1ChordAngle.negative().is_valid())
        self.assertTrue(s2.S1ChordAngle.infinity().is_valid())

    # Geometric operations

    def test_to_angle(self):
        a = s2.S1ChordAngle.from_degrees(90.0)
        angle = a.to_angle()
        self.assertAlmostEqual(angle.degrees, 90.0)

    def test_sin_cos_tan(self):
        a = s2.S1ChordAngle.from_degrees(30.0)
        self.assertAlmostEqual(a.sin(), 0.5, places=10)
        self.assertAlmostEqual(a.cos(), math.sqrt(3) / 2, places=10)
        self.assertAlmostEqual(a.tan(), 1.0 / math.sqrt(3), places=10)

    # Operators

    def test_equality(self):
        a = s2.S1ChordAngle.from_degrees(45.0)
        b = s2.S1ChordAngle.from_degrees(45.0)
        c = s2.S1ChordAngle.from_degrees(90.0)
        self.assertTrue(a == b)
        self.assertTrue(a != c)

    def test_comparison(self):
        a = s2.S1ChordAngle.from_degrees(30.0)
        b = s2.S1ChordAngle.from_degrees(60.0)
        self.assertTrue(a < b)
        self.assertTrue(b > a)
        self.assertTrue(a <= a)
        self.assertTrue(a >= a)

    def test_comparison_special_values(self):
        # Negative < Zero < finite < Infinity.
        neg = s2.S1ChordAngle.negative()
        zero = s2.S1ChordAngle.zero()
        finite = s2.S1ChordAngle.from_degrees(45.0)
        inf = s2.S1ChordAngle.infinity()
        self.assertLess(neg, zero)
        self.assertLess(zero, finite)
        self.assertLess(finite, inf)

    def test_add(self):
        a = s2.S1ChordAngle.from_degrees(30.0)
        b = s2.S1ChordAngle.from_degrees(60.0)
        total = a + b
        self.assertAlmostEqual(total.degrees, 90.0, places=5)

    def test_add_clamps_to_pi(self):
        a = s2.S1ChordAngle.from_degrees(120.0)
        b = s2.S1ChordAngle.from_degrees(120.0)
        total = a + b
        self.assertEqual(total, s2.S1ChordAngle.straight())

    def test_add_special_raises(self):
        a = s2.S1ChordAngle.from_degrees(30.0)
        with self.assertRaises(ValueError):
            _ = a + s2.S1ChordAngle.negative()
        with self.assertRaises(ValueError):
            _ = s2.S1ChordAngle.infinity() + a

    def test_sub(self):
        a = s2.S1ChordAngle.from_degrees(60.0)
        b = s2.S1ChordAngle.from_degrees(30.0)
        diff = a - b
        self.assertAlmostEqual(diff.degrees, 30.0, places=5)

    def test_sub_clamps_at_zero(self):
        a = s2.S1ChordAngle.from_degrees(30.0)
        b = s2.S1ChordAngle.from_degrees(60.0)
        diff = a - b
        self.assertTrue(diff.is_zero())

    def test_sub_special_raises(self):
        a = s2.S1ChordAngle.from_degrees(30.0)
        with self.assertRaises(ValueError):
            _ = a - s2.S1ChordAngle.negative()

    def test_iadd(self):
        a = s2.S1ChordAngle.from_degrees(30.0)
        a += s2.S1ChordAngle.from_degrees(60.0)
        self.assertAlmostEqual(a.degrees, 90.0, places=5)

    def test_isub(self):
        a = s2.S1ChordAngle.from_degrees(60.0)
        a -= s2.S1ChordAngle.from_degrees(30.0)
        self.assertAlmostEqual(a.degrees, 30.0, places=5)

    def test_hash(self):
        a = s2.S1ChordAngle.from_degrees(45.0)
        b = s2.S1ChordAngle.from_degrees(45.0)
        self.assertEqual(hash(a), hash(b))
        s = {a, b}
        self.assertEqual(len(s), 1)

    # String representation

    def test_repr(self):
        a = s2.S1ChordAngle.from_degrees(90.0)
        r = repr(a)
        self.assertTrue(r.startswith("S1ChordAngle("))
        self.assertTrue(r.endswith(")"))

    def test_str(self):
        a = s2.S1ChordAngle.from_degrees(90.0)
        s = str(a)
        self.assertIn("90", s)


if __name__ == "__main__":
    unittest.main()
