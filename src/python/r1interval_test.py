"""Tests for R1Interval pybind11 bindings."""

import unittest
import math
import s2geometry_pybind as s2


class TestR1Interval(unittest.TestCase):
    """Test cases for R1Interval bindings."""

    # Constructors

    def test_default_constructor(self):
        interval = s2.R1Interval()
        self.assertTrue(interval.is_empty())

    def test_constructor_with_bounds(self):
        interval = s2.R1Interval(1.0, 3.0)
        self.assertEqual(interval.lo, 1.0)
        self.assertEqual(interval.hi, 3.0)
        self.assertFalse(interval.is_empty())

    def test_constructor_lo_greater_than_hi_is_empty(self):
        interval = s2.R1Interval(3.0, 1.0)
        self.assertTrue(interval.is_empty())

    def test_constructor_negative_values(self):
        interval = s2.R1Interval(-5.0, -2.0)
        self.assertEqual(interval.lo, -5.0)
        self.assertEqual(interval.hi, -2.0)
        self.assertFalse(interval.is_empty())

    def test_constructor_large_values(self):
        interval = s2.R1Interval(-1e10, 1e10)
        self.assertEqual(interval.lo, -1e10)
        self.assertEqual(interval.hi, 1e10)

    # Static factory methods

    def test_empty(self):
        empty = s2.R1Interval.empty()
        self.assertTrue(empty.is_empty())

    def test_from_point(self):
        interval = s2.R1Interval.from_point(5.0)
        self.assertEqual(interval.lo, 5.0)
        self.assertEqual(interval.hi, 5.0)
        self.assertFalse(interval.is_empty())

    def test_from_point_pair(self):
        interval = s2.R1Interval.from_point_pair(5.0, 2.0)
        self.assertEqual(interval.lo, 2.0)
        self.assertEqual(interval.hi, 5.0)
        self.assertTrue(interval.contains(3.0))

    def test_from_point_pair_ordered(self):
        interval = s2.R1Interval.from_point_pair(2.0, 5.0)
        self.assertEqual(interval.lo, 2.0)
        self.assertEqual(interval.hi, 5.0)

    # Properties

    def test_properties_lo_hi(self):
        interval = s2.R1Interval(1.0, 3.0)
        self.assertEqual(interval.lo, 1.0)
        self.assertEqual(interval.hi, 3.0)

    def test_lo_hi_are_readonly(self):
        interval = s2.R1Interval(1.0, 3.0)
        with self.assertRaises(AttributeError):
            interval.lo = 0.5
        with self.assertRaises(AttributeError):
            interval.hi = 0.5

    def test_bounds(self):
        interval = s2.R1Interval(1.0, 3.0)
        bounds = interval.bounds()
        self.assertEqual(bounds, (1.0, 3.0))

    # Predicates

    def test_is_empty(self):
        self.assertTrue(s2.R1Interval().is_empty())
        self.assertTrue(s2.R1Interval.empty().is_empty())
        self.assertFalse(s2.R1Interval(1.0, 3.0).is_empty())
        self.assertFalse(s2.R1Interval(1.0, 1.0).is_empty())  # single point

    # Geometric operations

    def test_center(self):
        interval = s2.R1Interval(1.0, 3.0)
        self.assertAlmostEqual(interval.center(), 2.0)

    def test_length(self):
        interval = s2.R1Interval(1.0, 3.0)
        self.assertAlmostEqual(interval.length(), 2.0)

    def test_length_empty_is_negative(self):
        empty = s2.R1Interval.empty()
        self.assertLess(empty.length(), 0)

    def test_contains_point(self):
        interval = s2.R1Interval(1.0, 3.0)
        self.assertTrue(interval.contains(1.0))   # lo boundary
        self.assertTrue(interval.contains(2.0))   # interior
        self.assertTrue(interval.contains(3.0))   # hi boundary
        self.assertFalse(interval.contains(0.5))  # below
        self.assertFalse(interval.contains(3.5))  # above

    def test_interior_contains_point(self):
        interval = s2.R1Interval(1.0, 3.0)
        self.assertFalse(interval.interior_contains(1.0))  # lo boundary
        self.assertTrue(interval.interior_contains(2.0))   # interior
        self.assertFalse(interval.interior_contains(3.0))  # hi boundary

    def test_contains_interval(self):
        outer = s2.R1Interval(1.0, 5.0)
        inner = s2.R1Interval(2.0, 4.0)
        disjoint = s2.R1Interval(6.0, 8.0)
        self.assertTrue(outer.contains(inner))
        self.assertFalse(inner.contains(outer))
        self.assertFalse(outer.contains(disjoint))

    def test_contains_empty_interval(self):
        interval = s2.R1Interval(1.0, 3.0)
        empty = s2.R1Interval.empty()
        self.assertTrue(interval.contains(empty))

    def test_interior_contains_interval(self):
        outer = s2.R1Interval(1.0, 5.0)
        inner = s2.R1Interval(2.0, 4.0)
        boundary = s2.R1Interval(1.0, 5.0)
        self.assertTrue(outer.interior_contains(inner))
        self.assertFalse(outer.interior_contains(boundary))

    def test_intersects(self):
        i1 = s2.R1Interval(1.0, 3.0)
        i2 = s2.R1Interval(2.0, 4.0)
        i3 = s2.R1Interval(4.0, 6.0)
        self.assertTrue(i1.intersects(i2))
        self.assertFalse(i1.intersects(i3))

    def test_interior_intersects(self):
        i1 = s2.R1Interval(1.0, 3.0)
        i2 = s2.R1Interval(2.0, 4.0)
        touching = s2.R1Interval(3.0, 5.0)
        self.assertTrue(i1.interior_intersects(i2))
        self.assertFalse(i1.interior_intersects(touching))

    def test_add_point(self):
        interval = s2.R1Interval(1.0, 3.0)
        interval.add_point(5.0)
        self.assertEqual(interval.hi, 5.0)
        interval.add_point(-1.0)
        self.assertEqual(interval.lo, -1.0)

    def test_add_point_to_empty(self):
        interval = s2.R1Interval()
        interval.add_point(2.0)
        self.assertEqual(interval.lo, 2.0)
        self.assertEqual(interval.hi, 2.0)

    def test_add_interval(self):
        i1 = s2.R1Interval(1.0, 3.0)
        i2 = s2.R1Interval(2.0, 5.0)
        i1.add_interval(i2)
        self.assertEqual(i1.lo, 1.0)
        self.assertEqual(i1.hi, 5.0)

    def test_project(self):
        interval = s2.R1Interval(1.0, 3.0)
        self.assertAlmostEqual(interval.project(2.0), 2.0)  # inside
        self.assertAlmostEqual(interval.project(0.0), 1.0)  # below
        self.assertAlmostEqual(interval.project(5.0), 3.0)  # above

    def test_project_on_empty_raises(self):
        with self.assertRaises(ValueError):
            s2.R1Interval.empty().project(1.0)

    def test_expanded(self):
        interval = s2.R1Interval(1.0, 3.0)
        expanded = interval.expanded(0.5)
        self.assertAlmostEqual(expanded.lo, 0.5)
        self.assertAlmostEqual(expanded.hi, 3.5)

    def test_expanded_shrink(self):
        interval = s2.R1Interval(1.0, 3.0)
        shrunk = interval.expanded(-0.5)
        self.assertAlmostEqual(shrunk.lo, 1.5)
        self.assertAlmostEqual(shrunk.hi, 2.5)

    def test_expanded_empty_stays_empty(self):
        empty = s2.R1Interval.empty()
        expanded = empty.expanded(1.0)
        self.assertTrue(expanded.is_empty())

    def test_union(self):
        i1 = s2.R1Interval(1.0, 3.0)
        i2 = s2.R1Interval(5.0, 7.0)
        union = i1.union(i2)
        self.assertEqual(union.lo, 1.0)
        self.assertEqual(union.hi, 7.0)

    def test_intersection(self):
        i1 = s2.R1Interval(1.0, 5.0)
        i2 = s2.R1Interval(3.0, 7.0)
        intersection = i1.intersection(i2)
        self.assertEqual(intersection.lo, 3.0)
        self.assertEqual(intersection.hi, 5.0)

    def test_intersection_disjoint(self):
        i1 = s2.R1Interval(1.0, 3.0)
        i2 = s2.R1Interval(5.0, 7.0)
        intersection = i1.intersection(i2)
        self.assertTrue(intersection.is_empty())

    def test_directed_hausdorff_distance(self):
        i1 = s2.R1Interval(1.0, 3.0)
        i2 = s2.R1Interval(2.0, 5.0)
        dist = i1.directed_hausdorff_distance(i2)
        # max of (1-2, 3-5) clamped to 0 = max(0, 0) = 0
        # Actually: max(0, max(hi-y.hi, y.lo-lo)) = max(0, max(3-5, 2-1)) = max(0, 1) = 1
        self.assertAlmostEqual(dist, 1.0)

    def test_directed_hausdorff_distance_empty(self):
        empty = s2.R1Interval.empty()
        i = s2.R1Interval(1.0, 3.0)
        self.assertAlmostEqual(empty.directed_hausdorff_distance(i), 0.0)

    def test_approx_equals(self):
        i1 = s2.R1Interval(1.0, 3.0)
        i2 = s2.R1Interval(1.0, 3.0)
        i3 = s2.R1Interval(1.0, 3.0 + 1e-16)
        i4 = s2.R1Interval(1.0, 3.1)

        self.assertTrue(i1.approx_equals(i2))
        self.assertTrue(i1.approx_equals(i3))
        self.assertFalse(i1.approx_equals(i4))
        self.assertTrue(i1.approx_equals(i4, 0.2))

    # Operators

    def test_equality(self):
        i1 = s2.R1Interval(1.0, 3.0)
        i2 = s2.R1Interval(1.0, 3.0)
        i3 = s2.R1Interval(1.0, 4.0)
        self.assertTrue(i1 == i2)
        self.assertTrue(i1 != i3)

    def test_empty_intervals_are_equal(self):
        e1 = s2.R1Interval.empty()
        e2 = s2.R1Interval()
        self.assertTrue(e1 == e2)

    # String representation

    def test_string_representation(self):
        interval = s2.R1Interval(1.0, 3.0)
        self.assertEqual(repr(interval), "R1Interval([1, 3])")
        self.assertEqual(str(interval), "[1, 3]")

    def test_string_representation_empty(self):
        empty = s2.R1Interval.empty()
        self.assertEqual(repr(empty), "R1Interval([1, 0])")
        self.assertEqual(str(empty), "[1, 0]")


if __name__ == "__main__":
    unittest.main()
