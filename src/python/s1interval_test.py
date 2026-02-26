"""Tests for S1Interval pybind11 bindings."""

import unittest
import math
import s2geometry_pybind as s2


class TestS1Interval(unittest.TestCase):
    """Test cases for S1Interval bindings."""

    # Constructors

    def test_default_constructor(self):
        interval = s2.S1Interval()
        self.assertTrue(interval.is_empty())

    def test_constructor_with_bounds(self):
        interval = s2.S1Interval(0.0, math.pi / 2)
        self.assertAlmostEqual(interval.lo, 0.0)
        self.assertAlmostEqual(interval.hi, math.pi / 2)
        self.assertFalse(interval.is_empty())

    def test_constructor_converts_minus_pi_to_pi(self):
        # "The value -π is converted internally to π except for the Full() and Empty() intervals"
        # (from s1interval.h constructor documentation)
        interval = s2.S1Interval(-math.pi, 0.0)
        self.assertEqual(interval.lo, math.pi)  # -π converted to π
        self.assertEqual(interval.hi, 0.0)
        self.assertTrue(interval.is_inverted())

    # Static factory methods

    def test_empty_constructor(self):
        empty = s2.S1Interval.empty()
        self.assertTrue(empty.is_empty())

    def test_full_constructor(self):
        full = s2.S1Interval.full()
        self.assertTrue(full.is_full())

    def test_from_point(self):
        point_interval = s2.S1Interval.from_point(1.0)
        self.assertEqual(point_interval.lo, 1.0)
        self.assertEqual(point_interval.hi, 1.0)

    def test_from_point_pair(self):
        pair_interval = s2.S1Interval.from_point_pair(0.5, 2.0)
        self.assertTrue(pair_interval.contains(1.0))

    def test_from_point_pair_invalid_raises(self):
        with self.assertRaises(ValueError):
            s2.S1Interval.from_point_pair(0.5, 4.0)  # 4.0 > π
        with self.assertRaises(ValueError):
            s2.S1Interval.from_point_pair(-4.0, 0.5)  # -4.0 < -π

    # Properties

    def test_properties_lo_hi(self):
        interval = s2.S1Interval(0.0, math.pi / 2)
        self.assertAlmostEqual(interval.lo, 0.0)
        self.assertAlmostEqual(interval.hi, math.pi / 2)

    def test_lo_hi_are_readonly(self):
        interval = s2.S1Interval(0.0, math.pi / 2)
        with self.assertRaises(AttributeError):
            interval.lo = 0.5
        with self.assertRaises(AttributeError):
            interval.hi = 0.5

    def test_bounds(self):
        interval = s2.S1Interval(0.5, 2.0)
        bounds = interval.bounds()
        self.assertEqual(bounds, (0.5, 2.0))
        self.assertAlmostEqual(bounds[0], interval.lo)
        self.assertAlmostEqual(bounds[1], interval.hi)

    # Predicates

    def test_invalid_bounds_too_large(self):
        with self.assertRaises(ValueError):
            s2.S1Interval(0.0, 4.0)  # 4.0 > π

    def test_invalid_bounds_too_small(self):
        with self.assertRaises(ValueError):
            s2.S1Interval(-4.0, 0.0)  # -4.0 < -π

    def test_invalid_both_bounds_out_of_range(self):
        with self.assertRaises(ValueError):
            s2.S1Interval(-4.0, 4.0)

    def test_is_full_and_empty(self):
        empty = s2.S1Interval.empty()
        self.assertTrue(empty.is_empty())
        self.assertFalse(empty.is_full())

        full = s2.S1Interval.full()
        self.assertTrue(full.is_full())
        self.assertFalse(full.is_empty())

        normal = s2.S1Interval(0.0, math.pi)
        self.assertFalse(normal.is_empty())
        self.assertFalse(normal.is_full())

    def test_is_inverted(self):
        # Normal interval
        normal = s2.S1Interval(0.0, math.pi)
        self.assertFalse(normal.is_inverted())

        # Inverted interval (wraps around)
        inverted = s2.S1Interval(math.pi, 0.0)
        self.assertTrue(inverted.is_inverted())

    # Geometric operations

    def test_center(self):
        interval = s2.S1Interval(0.0, math.pi)
        self.assertAlmostEqual(interval.center(), math.pi / 2)

    def test_length(self):
        interval = s2.S1Interval(0.0, math.pi)
        self.assertAlmostEqual(interval.length(), math.pi)

    def test_complement_center(self):
        interval = s2.S1Interval(0.0, math.pi / 2)
        complement_center = interval.complement_center()
        self.assertIsNotNone(complement_center)

    def test_contains_point(self):
        interval = s2.S1Interval(0.0, math.pi)
        self.assertTrue(interval.contains(math.pi / 2))
        self.assertTrue(interval.contains(0.0))
        self.assertTrue(interval.contains(math.pi))
        self.assertFalse(interval.contains(-math.pi / 2))

    def test_contains_point_invalid_raises(self):
        interval = s2.S1Interval(0.0, math.pi)
        with self.assertRaises(ValueError):
            interval.contains(4.0)  # 4.0 > π

    def test_interior_contains_point(self):
        interval = s2.S1Interval(0.0, math.pi)
        self.assertTrue(interval.interior_contains(math.pi / 2))
        self.assertFalse(interval.interior_contains(0.0))
        self.assertFalse(interval.interior_contains(math.pi))

    def test_interior_contains_point_invalid_raises(self):
        interval = s2.S1Interval(0.0, math.pi)
        with self.assertRaises(ValueError):
            interval.interior_contains(4.0)  # 4.0 > π

    def test_contains_interval(self):
        interval1 = s2.S1Interval(0.0, math.pi)
        interval2 = s2.S1Interval(0.5, 2.0)
        interval3 = s2.S1Interval(-math.pi, -2.0)
        self.assertTrue(interval1.contains(interval2))
        self.assertFalse(interval1.contains(interval3))

    def test_interior_contains_interval(self):
        interval1 = s2.S1Interval(0.0, math.pi)
        interval2 = s2.S1Interval(0.5, 2.0)
        self.assertTrue(interval1.interior_contains(interval2))

    def test_intersects(self):
        interval1 = s2.S1Interval(0.0, 2.0)
        interval2 = s2.S1Interval(1.0, 3.0)
        interval3 = s2.S1Interval(-3.0, -2.5)
        self.assertTrue(interval1.intersects(interval2))
        self.assertFalse(interval1.intersects(interval3))

    def test_interior_intersects(self):
        interval1 = s2.S1Interval(0.0, 2.0)
        interval2 = s2.S1Interval(1.0, 3.0)
        interval3 = s2.S1Interval(-3.0, -2.5)
        self.assertTrue(interval1.interior_intersects(interval2))
        self.assertFalse(interval1.interior_intersects(interval3))

    def test_add_point(self):
        interval = s2.S1Interval(0.0, 1.0)
        interval.add_point(1.5)
        self.assertTrue(interval.contains(1.5))

    def test_add_point_invalid_raises(self):
        interval = s2.S1Interval(0.0, 1.0)
        with self.assertRaises(ValueError):
            interval.add_point(4.0)  # 4.0 > π

    def test_project(self):
        interval = s2.S1Interval(0.0, 1.5)
        self.assertAlmostEqual(interval.project(1.0), 1.0)  # Point inside interval

    def test_project_on_empty_raises(self):
        with self.assertRaises(ValueError):
            s2.S1Interval.empty().project(1.0)

    def test_project_invalid_point_raises(self):
        interval = s2.S1Interval(0.0, 1.5)
        with self.assertRaises(ValueError):
            interval.project(4.0)  # 4.0 > π

    def test_expanded(self):
        interval = s2.S1Interval(0.5, 1.5)
        expanded = interval.expanded(0.5)
        self.assertAlmostEqual(expanded.length(), 
                             interval.length() + 1.0)

    def test_union(self):
        interval1 = s2.S1Interval(0.0, 1.0)
        interval2 = s2.S1Interval(0.5, 1.5)
        union = interval1.union(interval2)
        self.assertTrue(union.contains(0.5))
        self.assertTrue(union.contains(1.2))

    def test_intersection(self):
        interval1 = s2.S1Interval(0.0, 2.0)
        interval2 = s2.S1Interval(1.0, 3.0)
        intersection = interval1.intersection(interval2)
        self.assertAlmostEqual(intersection.lo, 1.0)
        self.assertAlmostEqual(intersection.hi, 2.0)

    def test_complement(self):
        interval = s2.S1Interval(0.0, math.pi / 2)
        complement = interval.complement()
        self.assertTrue(complement.is_inverted())
        self.assertFalse(complement.contains(math.pi / 4))
        self.assertTrue(complement.contains(math.pi))

    def test_directed_hausdorff_distance(self):
        interval1 = s2.S1Interval(0.0, 1.0)
        interval2 = s2.S1Interval(0.5, 1.5)
        distance = interval1.directed_hausdorff_distance(interval2)
        self.assertIsNotNone(distance)

    def test_approx_equals(self):
        interval1 = s2.S1Interval(0.0, 2.0)
        interval2 = s2.S1Interval(0.0, 2.0)
        interval3 = s2.S1Interval(0.0, 2.0 + 1e-16)  # Very small difference
        interval4 = s2.S1Interval(0.0, 2.1)  # Larger difference
        
        # Test with default max_error (1e-15)
        self.assertTrue(interval1.approx_equals(interval2))
        self.assertTrue(interval1.approx_equals(interval3))  # Within default tolerance
        self.assertFalse(interval1.approx_equals(interval4))  # Outside default tolerance
        
        # Test with explicit max_error values
        self.assertTrue(interval1.approx_equals(interval3, 1e-15))
        self.assertFalse(interval1.approx_equals(interval4, 0.05))
        self.assertTrue(interval1.approx_equals(interval4, 0.2))

    # Operators

    def test_equality(self):
        interval1 = s2.S1Interval(0.0, 2.0)
        interval2 = s2.S1Interval(0.0, 2.0)
        interval3 = s2.S1Interval(0.0, 3.0)
        self.assertTrue(interval1 == interval2)
        self.assertTrue(interval1 != interval3)

    # String representation

    def test_string_representation(self):
        interval = s2.S1Interval(0.0, 2.0)
        self.assertEqual(repr(interval), "S1Interval([0, 2])")
        self.assertEqual(str(interval), "[0, 2]")

        empty = s2.S1Interval.empty()
        # Empty interval has lo=pi, hi=-pi
        self.assertEqual(repr(empty), "S1Interval([3.14159, -3.14159])")
        self.assertEqual(str(empty), "[3.14159, -3.14159]")
        
        full = s2.S1Interval.full()
        # Full interval spans from -pi to pi
        self.assertEqual(repr(full), "S1Interval([-3.14159, 3.14159])")
        self.assertEqual(str(full), "[-3.14159, 3.14159]")


if __name__ == "__main__":
    unittest.main()
