"""Tests for R2Rect pybind11 bindings."""

import unittest
import s2geometry_pybind as s2


class TestR2Rect(unittest.TestCase):
    """Test cases for R2Rect bindings."""

    # Constructors

    def test_default_constructor(self):
        rect = s2.R2Rect()
        self.assertTrue(rect.is_empty())

    def test_constructor_from_points(self):
        lo = s2.R2Point(1.0, 2.0)
        hi = s2.R2Point(3.0, 4.0)
        rect = s2.R2Rect(lo, hi)
        self.assertEqual(rect.lo, lo)
        self.assertEqual(rect.hi, hi)
        self.assertFalse(rect.is_empty())

    def test_constructor_from_intervals(self):
        x = s2.R1Interval(1.0, 3.0)
        y = s2.R1Interval(2.0, 4.0)
        rect = s2.R2Rect(x, y)
        self.assertEqual(rect.x, x)
        self.assertEqual(rect.y, y)

    def test_constructor_from_intervals_both_empty(self):
        x = s2.R1Interval.empty()
        y = s2.R1Interval.empty()
        rect = s2.R2Rect(x, y)
        self.assertTrue(rect.is_empty())

    def test_constructor_from_intervals_mixed_empty_raises(self):
        x = s2.R1Interval(1.0, 3.0)
        y = s2.R1Interval.empty()
        with self.assertRaises(ValueError):
            s2.R2Rect(x, y)

    def test_constructor_from_points_mixed_empty_raises(self):
        # lo.x > hi.x but lo.y < hi.y -> one axis empty, other not
        lo = s2.R2Point(3.0, 1.0)
        hi = s2.R2Point(1.0, 4.0)
        with self.assertRaises(ValueError):
            s2.R2Rect(lo, hi)

    # Static factory methods

    def test_empty(self):
        rect = s2.R2Rect.empty()
        self.assertTrue(rect.is_empty())

    def test_from_point(self):
        p = s2.R2Point(1.0, 2.0)
        rect = s2.R2Rect.from_point(p)
        self.assertEqual(rect.lo, p)
        self.assertEqual(rect.hi, p)
        self.assertFalse(rect.is_empty())

    def test_from_point_pair(self):
        p1 = s2.R2Point(3.0, 4.0)
        p2 = s2.R2Point(1.0, 2.0)
        rect = s2.R2Rect.from_point_pair(p1, p2)
        self.assertAlmostEqual(rect.lo.x, 1.0)
        self.assertAlmostEqual(rect.lo.y, 2.0)
        self.assertAlmostEqual(rect.hi.x, 3.0)
        self.assertAlmostEqual(rect.hi.y, 4.0)

    def test_from_center_size(self):
        center = s2.R2Point(2.0, 3.0)
        size = s2.R2Point(4.0, 6.0)
        rect = s2.R2Rect.from_center_size(center, size)
        self.assertAlmostEqual(rect.lo.x, 0.0)
        self.assertAlmostEqual(rect.lo.y, 0.0)
        self.assertAlmostEqual(rect.hi.x, 4.0)
        self.assertAlmostEqual(rect.hi.y, 6.0)

    def test_from_center_size_negative_raises(self):
        center = s2.R2Point(0.0, 0.0)
        with self.assertRaises(ValueError):
            s2.R2Rect.from_center_size(center, s2.R2Point(-1.0, 2.0))
        with self.assertRaises(ValueError):
            s2.R2Rect.from_center_size(center, s2.R2Point(1.0, -2.0))

    # Properties

    def test_lo_hi(self):
        rect = s2.R2Rect(s2.R2Point(1.0, 2.0), s2.R2Point(3.0, 4.0))
        self.assertAlmostEqual(rect.lo.x, 1.0)
        self.assertAlmostEqual(rect.lo.y, 2.0)
        self.assertAlmostEqual(rect.hi.x, 3.0)
        self.assertAlmostEqual(rect.hi.y, 4.0)

    def test_x_y(self):
        rect = s2.R2Rect(s2.R2Point(1.0, 2.0), s2.R2Point(3.0, 4.0))
        self.assertEqual(rect.x, s2.R1Interval(1.0, 3.0))
        self.assertEqual(rect.y, s2.R1Interval(2.0, 4.0))

    def test_properties_are_readonly(self):
        rect = s2.R2Rect(s2.R2Point(1.0, 2.0), s2.R2Point(3.0, 4.0))
        with self.assertRaises(AttributeError):
            rect.lo = s2.R2Point(0.0, 0.0)
        with self.assertRaises(AttributeError):
            rect.hi = s2.R2Point(0.0, 0.0)
        with self.assertRaises(AttributeError):
            rect.x = s2.R1Interval(0.0, 1.0)
        with self.assertRaises(AttributeError):
            rect.y = s2.R1Interval(0.0, 1.0)

    # Predicates

    def test_is_valid(self):
        rect = s2.R2Rect(s2.R2Point(1.0, 2.0), s2.R2Point(3.0, 4.0))
        self.assertTrue(rect.is_valid())
        self.assertTrue(s2.R2Rect.empty().is_valid())

    def test_is_empty(self):
        self.assertTrue(s2.R2Rect().is_empty())
        self.assertTrue(s2.R2Rect.empty().is_empty())
        rect = s2.R2Rect(s2.R2Point(1.0, 2.0), s2.R2Point(3.0, 4.0))
        self.assertFalse(rect.is_empty())

    # Geometric operations

    def test_vertex(self):
        rect = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(1.0, 1.0))
        # CCW: lower-left, lower-right, upper-right, upper-left
        v0 = rect.vertex(0)
        v1 = rect.vertex(1)
        v2 = rect.vertex(2)
        v3 = rect.vertex(3)
        self.assertEqual(v0, s2.R2Point(0.0, 0.0))
        self.assertEqual(v1, s2.R2Point(1.0, 0.0))
        self.assertEqual(v2, s2.R2Point(1.0, 1.0))
        self.assertEqual(v3, s2.R2Point(0.0, 1.0))

    def test_vertex_ij(self):
        rect = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(1.0, 1.0))
        self.assertEqual(rect.vertex_ij(0, 0), s2.R2Point(0.0, 0.0))
        self.assertEqual(rect.vertex_ij(1, 0), s2.R2Point(1.0, 0.0))
        self.assertEqual(rect.vertex_ij(0, 1), s2.R2Point(0.0, 1.0))
        self.assertEqual(rect.vertex_ij(1, 1), s2.R2Point(1.0, 1.0))

    def test_center(self):
        rect = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(2.0, 4.0))
        center = rect.center()
        self.assertAlmostEqual(center.x, 1.0)
        self.assertAlmostEqual(center.y, 2.0)

    def test_size(self):
        rect = s2.R2Rect(s2.R2Point(1.0, 2.0), s2.R2Point(3.0, 6.0))
        size = rect.size()
        self.assertAlmostEqual(size.x, 2.0)
        self.assertAlmostEqual(size.y, 4.0)

    def test_contains_point(self):
        rect = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(2.0, 2.0))
        self.assertTrue(rect.contains(s2.R2Point(1.0, 1.0)))   # interior
        self.assertTrue(rect.contains(s2.R2Point(0.0, 0.0)))   # corner
        self.assertTrue(rect.contains(s2.R2Point(2.0, 2.0)))   # corner
        self.assertFalse(rect.contains(s2.R2Point(3.0, 1.0)))  # outside

    def test_interior_contains_point(self):
        rect = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(2.0, 2.0))
        self.assertTrue(rect.interior_contains(s2.R2Point(1.0, 1.0)))
        self.assertFalse(rect.interior_contains(s2.R2Point(0.0, 0.0)))
        self.assertFalse(rect.interior_contains(s2.R2Point(2.0, 1.0)))

    def test_contains_rect(self):
        outer = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(4.0, 4.0))
        inner = s2.R2Rect(s2.R2Point(1.0, 1.0), s2.R2Point(3.0, 3.0))
        disjoint = s2.R2Rect(s2.R2Point(5.0, 5.0), s2.R2Point(6.0, 6.0))
        self.assertTrue(outer.contains(inner))
        self.assertFalse(inner.contains(outer))
        self.assertFalse(outer.contains(disjoint))

    def test_interior_contains_rect(self):
        outer = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(4.0, 4.0))
        inner = s2.R2Rect(s2.R2Point(1.0, 1.0), s2.R2Point(3.0, 3.0))
        same = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(4.0, 4.0))
        self.assertTrue(outer.interior_contains(inner))
        self.assertFalse(outer.interior_contains(same))

    def test_intersects(self):
        r1 = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(2.0, 2.0))
        r2 = s2.R2Rect(s2.R2Point(1.0, 1.0), s2.R2Point(3.0, 3.0))
        r3 = s2.R2Rect(s2.R2Point(3.0, 3.0), s2.R2Point(4.0, 4.0))
        self.assertTrue(r1.intersects(r2))
        self.assertFalse(r1.intersects(r3))

    def test_interior_intersects(self):
        r1 = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(2.0, 2.0))
        r2 = s2.R2Rect(s2.R2Point(1.0, 1.0), s2.R2Point(3.0, 3.0))
        touching = s2.R2Rect(s2.R2Point(2.0, 2.0), s2.R2Point(3.0, 3.0))
        self.assertTrue(r1.interior_intersects(r2))
        self.assertFalse(r1.interior_intersects(touching))

    def test_add_point(self):
        rect = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(1.0, 1.0))
        rect.add_point(s2.R2Point(3.0, 3.0))
        self.assertAlmostEqual(rect.hi.x, 3.0)
        self.assertAlmostEqual(rect.hi.y, 3.0)

    def test_add_rect(self):
        r1 = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(1.0, 1.0))
        r2 = s2.R2Rect(s2.R2Point(2.0, 2.0), s2.R2Point(3.0, 3.0))
        r1.add_rect(r2)
        self.assertAlmostEqual(r1.lo.x, 0.0)
        self.assertAlmostEqual(r1.lo.y, 0.0)
        self.assertAlmostEqual(r1.hi.x, 3.0)
        self.assertAlmostEqual(r1.hi.y, 3.0)

    def test_project(self):
        rect = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(2.0, 2.0))
        # Point inside
        p = rect.project(s2.R2Point(1.0, 1.0))
        self.assertEqual(p, s2.R2Point(1.0, 1.0))
        # Point outside
        p = rect.project(s2.R2Point(5.0, 5.0))
        self.assertEqual(p, s2.R2Point(2.0, 2.0))
        # Point below
        p = rect.project(s2.R2Point(-1.0, -1.0))
        self.assertEqual(p, s2.R2Point(0.0, 0.0))

    def test_project_on_empty_raises(self):
        with self.assertRaises(ValueError):
            s2.R2Rect.empty().project(s2.R2Point(1.0, 1.0))

    def test_expanded(self):
        rect = s2.R2Rect(s2.R2Point(1.0, 2.0), s2.R2Point(3.0, 4.0))
        expanded = rect.expanded(s2.R2Point(0.5, 1.0))
        self.assertAlmostEqual(expanded.lo.x, 0.5)
        self.assertAlmostEqual(expanded.lo.y, 1.0)
        self.assertAlmostEqual(expanded.hi.x, 3.5)
        self.assertAlmostEqual(expanded.hi.y, 5.0)

    def test_union(self):
        r1 = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(1.0, 1.0))
        r2 = s2.R2Rect(s2.R2Point(2.0, 2.0), s2.R2Point(3.0, 3.0))
        union = r1.union(r2)
        self.assertAlmostEqual(union.lo.x, 0.0)
        self.assertAlmostEqual(union.lo.y, 0.0)
        self.assertAlmostEqual(union.hi.x, 3.0)
        self.assertAlmostEqual(union.hi.y, 3.0)

    def test_intersection(self):
        r1 = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(2.0, 2.0))
        r2 = s2.R2Rect(s2.R2Point(1.0, 1.0), s2.R2Point(3.0, 3.0))
        inter = r1.intersection(r2)
        self.assertAlmostEqual(inter.lo.x, 1.0)
        self.assertAlmostEqual(inter.lo.y, 1.0)
        self.assertAlmostEqual(inter.hi.x, 2.0)
        self.assertAlmostEqual(inter.hi.y, 2.0)

    def test_approx_equals(self):
        r1 = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(1.0, 1.0))
        r2 = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(1.0, 1.0))
        r3 = s2.R2Rect(s2.R2Point(0.0, 0.0),
                        s2.R2Point(1.0, 1.0 + 1e-16))
        r4 = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(1.0, 1.1))
        self.assertTrue(r1.approx_equals(r2))
        self.assertTrue(r1.approx_equals(r3))
        self.assertFalse(r1.approx_equals(r4))
        self.assertTrue(r1.approx_equals(r4, 0.2))

    # Operators

    def test_equality(self):
        r1 = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(1.0, 1.0))
        r2 = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(1.0, 1.0))
        r3 = s2.R2Rect(s2.R2Point(0.0, 0.0), s2.R2Point(2.0, 2.0))
        self.assertTrue(r1 == r2)
        self.assertTrue(r1 != r3)

    def test_empty_rects_are_equal(self):
        e1 = s2.R2Rect.empty()
        e2 = s2.R2Rect()
        self.assertTrue(e1 == e2)

    # String representation

    def test_string_representation(self):
        rect = s2.R2Rect(s2.R2Point(1.0, 2.0), s2.R2Point(3.0, 4.0))
        repr_str = repr(rect)
        str_str = str(rect)
        self.assertTrue(repr_str.startswith("R2Rect("))
        # str should not have the class name prefix
        self.assertFalse(str_str.startswith("R2Rect("))


if __name__ == "__main__":
    unittest.main()
