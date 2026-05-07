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
        self.assertEqual(a.length2, 0.0)

    def test_constructor_from_two_points(self):
        x = s2.S2Point(1.0, 0.0, 0.0)
        y = s2.S2Point(0.0, 1.0, 0.0)
        a = s2.S1ChordAngle(x, y)
        # Chord length between orthogonal unit vectors is sqrt(2);
        # length2 = 2.
        self.assertAlmostEqual(a.length2, 2.0)
        self.assertAlmostEqual(a.radians, math.pi / 2)

    def test_constructor_from_same_point(self):
        p = s2.S2Point(1.0, 0.0, 0.0)
        a = s2.S1ChordAngle(p, p)
        self.assertTrue(a.is_zero())

    def test_constructor_from_antipodal_points(self):
        x = s2.S2Point(1.0, 0.0, 0.0)
        y = s2.S2Point(-1.0, 0.0, 0.0)
        a = s2.S1ChordAngle(x, y)
        # Length2 for antipodal is 4.
        self.assertAlmostEqual(a.length2, 4.0)
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
        self.assertAlmostEqual(a.length2, 2.0)

    def test_constructor_from_s1_angle_negative(self):
        # Negative angles map to Negative().
        a = s2.S1ChordAngle(s2.S1Angle.from_radians(-1.0))
        self.assertTrue(a.is_negative())

    def test_constructor_from_s1_angle_infinity(self):
        a = s2.S1ChordAngle(s2.S1Angle.infinity())
        self.assertTrue(a.is_infinity())

    def test_constructor_from_s1_angle_larger_than_pi_clamps(self):
        a = s2.S1ChordAngle(s2.S1Angle.from_radians(math.pi + 1.0))
        # Clamped to Straight() (length2 == 4).
        self.assertAlmostEqual(a.length2, 4.0)

    # Factory methods

    def test_from_radians(self):
        a = s2.S1ChordAngle.from_radians(math.pi / 2)
        self.assertAlmostEqual(a.length2, 2.0)

    def test_from_degrees(self):
        a = s2.S1ChordAngle.from_degrees(90.0)
        self.assertAlmostEqual(a.length2, 2.0)

    def test_from_e5(self):
        a = s2.S1ChordAngle.from_e5(9000000)
        self.assertAlmostEqual(a.degrees, 90.0, places=3)

    def test_from_e6(self):
        a = s2.S1ChordAngle.from_e6(90000000)
        self.assertAlmostEqual(a.degrees, 90.0, places=5)

    def test_from_e7(self):
        a = s2.S1ChordAngle.from_e7(900000000)
        self.assertAlmostEqual(a.degrees, 90.0, places=7)

    def test_fast_upper_bound_from(self):
        angle = s2.S1Angle.from_radians(0.1)
        bound = s2.S1ChordAngle.fast_upper_bound_from(angle)
        exact = s2.S1ChordAngle(angle)
        # The fast bound is an upper bound on the exact value.
        self.assertGreaterEqual(bound.length2, exact.length2)

    def test_from_length2(self):
        a = s2.S1ChordAngle.from_length2(2.0)
        self.assertAlmostEqual(a.length2, 2.0)

    def test_from_length2_clamps_at_4(self):
        # Values above 4 are clamped (handles roundoff).
        a = s2.S1ChordAngle.from_length2(5.0)
        self.assertAlmostEqual(a.length2, 4.0)

    def test_from_length2_negative_raises(self):
        with self.assertRaises(ValueError):
            s2.S1ChordAngle.from_length2(-0.1)

    def test_zero(self):
        a = s2.S1ChordAngle.zero()
        self.assertTrue(a.is_zero())
        self.assertEqual(a.length2, 0.0)

    def test_right(self):
        a = s2.S1ChordAngle.right()
        self.assertAlmostEqual(a.length2, 2.0)
        self.assertAlmostEqual(a.radians, math.pi / 2)

    def test_straight(self):
        a = s2.S1ChordAngle.straight()
        self.assertAlmostEqual(a.length2, 4.0)
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

    def test_length2(self):
        a = s2.S1ChordAngle.from_length2(1.5)
        self.assertAlmostEqual(a.length2, 1.5)

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

    def test_sin2(self):
        a = s2.S1ChordAngle.from_degrees(30.0)
        self.assertAlmostEqual(a.sin2(), 0.25, places=10)

    def test_successor(self):
        # Negative().successor() == Zero().
        self.assertEqual(
            s2.S1ChordAngle.negative().successor(), s2.S1ChordAngle.zero())
        # Straight().successor() == Infinity().
        self.assertEqual(
            s2.S1ChordAngle.straight().successor(), s2.S1ChordAngle.infinity())
        # Infinity().successor() == Infinity() (fixed point).
        self.assertEqual(
            s2.S1ChordAngle.infinity().successor(), s2.S1ChordAngle.infinity())
        # Successor is strictly greater for non-special values.
        a = s2.S1ChordAngle.from_degrees(45.0)
        self.assertGreater(a.successor(), a)

    def test_predecessor(self):
        # Infinity().predecessor() == Straight().
        self.assertEqual(
            s2.S1ChordAngle.infinity().predecessor(),
            s2.S1ChordAngle.straight())
        # Zero().predecessor() == Negative().
        self.assertEqual(
            s2.S1ChordAngle.zero().predecessor(), s2.S1ChordAngle.negative())
        # Negative().predecessor() == Negative() (fixed point).
        self.assertEqual(
            s2.S1ChordAngle.negative().predecessor(),
            s2.S1ChordAngle.negative())
        a = s2.S1ChordAngle.from_degrees(45.0)
        self.assertLess(a.predecessor(), a)

    def test_plus_error(self):
        a = s2.S1ChordAngle.from_degrees(45.0)
        adjusted = a.plus_error(0.0)
        self.assertAlmostEqual(a.length2, adjusted.length2)
        adjusted2 = a.plus_error(1e-10)
        self.assertGreater(adjusted2.length2, a.length2)

    def test_constructor_max_error(self):
        a = s2.S1ChordAngle.from_degrees(45.0)
        self.assertGreater(a.s2_point_constructor_max_error(), 0.0)
        self.assertGreater(a.s1_angle_constructor_max_error(), 0.0)

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
        # Adding two large chord angles clamps the result to Straight().
        a = s2.S1ChordAngle.from_degrees(120.0)
        b = s2.S1ChordAngle.from_degrees(120.0)
        total = a + b
        self.assertAlmostEqual(total.length2, 4.0)

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
        # Subtracting a larger chord angle clamps to zero.
        a = s2.S1ChordAngle.from_degrees(30.0)
        b = s2.S1ChordAngle.from_degrees(60.0)
        diff = a - b
        self.assertAlmostEqual(diff.length2, 0.0)

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
        # The stream inserter prints the equivalent S1Angle.
        s = str(a)
        self.assertIn("90", s)


if __name__ == "__main__":
    unittest.main()
