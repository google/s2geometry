"""Tests for S2Cap pybind11 bindings."""

import math
import unittest
import s2geometry_pybind as s2


class TestS2Cap(unittest.TestCase):
    """Test cases for S2Cap bindings."""

    # Constructors

    def test_default_constructor_is_empty(self):
        cap = s2.S2Cap()
        self.assertTrue(cap.is_empty())
        self.assertFalse(cap.is_full())

    def test_constructor_s1angle(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        radius = s2.S1Angle.from_radians(math.pi / 4)
        cap = s2.S2Cap(center, radius)
        self.assertFalse(cap.is_empty())
        self.assertFalse(cap.is_full())

    def test_constructor_s1chord_angle(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        radius = s2.S1ChordAngle(s2.S1Angle.from_radians(math.pi / 4))
        cap = s2.S2Cap(center, radius)
        self.assertFalse(cap.is_empty())

    def test_constructor_negative_radius_is_empty(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(-1.0))
        self.assertTrue(cap.is_empty())

    def test_constructor_pi_radius_is_full(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(math.pi))
        self.assertTrue(cap.is_full())

    def test_constructor_normalizes_center(self):
        p = s2.S2Point(2.0, 0.0, 0.0)  # not unit length
        cap = s2.S2Cap(p, s2.S1Angle.from_degrees(10.0))
        self.assertAlmostEqual(cap.center.norm(), 1.0, places=15)

    # Factory methods

    def test_from_point(self):
        p = s2.S2Point(0.0, 1.0, 0.0)
        cap = s2.S2Cap.from_point(p)
        self.assertTrue(cap.contains_point(p))

    def test_from_center_height(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap.from_center_height(center, 0.5)
        self.assertAlmostEqual(cap.height, 0.5)

    def test_from_center_height_negative_is_empty(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        self.assertTrue(s2.S2Cap.from_center_height(center, -1.0).is_empty())

    def test_from_center_height_two_is_full(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        self.assertTrue(s2.S2Cap.from_center_height(center, 2.0).is_full())

    def test_from_center_area(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        area = 2 * math.pi  # hemisphere
        cap = s2.S2Cap.from_center_area(center, area)
        self.assertAlmostEqual(cap.area(), area, places=10)

    def test_empty(self):
        cap = s2.S2Cap.empty()
        self.assertTrue(cap.is_empty())
        self.assertFalse(cap.is_full())

    def test_full(self):
        cap = s2.S2Cap.full()
        self.assertTrue(cap.is_full())
        self.assertFalse(cap.is_empty())

    def test_from_points_empty(self):
        cap = s2.S2Cap.from_points([])
        self.assertTrue(cap.is_empty())

    def test_from_points_single(self):
        p = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap.from_points([p])
        self.assertTrue(cap.contains_point(p))

    def test_from_points_multiple(self):
        p1 = s2.S2Point(1.0, 0.0, 0.0)
        p2 = s2.S2Point(0.0, 1.0, 0.0)
        p3 = s2.S2Point(0.0, 0.0, 1.0)
        cap = s2.S2Cap.from_points([p1, p2, p3])
        self.assertTrue(cap.contains_point(p1))
        self.assertTrue(cap.contains_point(p2))
        self.assertTrue(cap.contains_point(p3))

    def test_from_caps_empty(self):
        cap = s2.S2Cap.from_caps([])
        self.assertTrue(cap.is_empty())

    def test_from_caps_single(self):
        c = s2.S2Cap(s2.S2Point(1.0, 0.0, 0.0), s2.S1Angle.from_degrees(10.0))
        result = s2.S2Cap.from_caps([c])
        self.assertTrue(result.contains_point(c.center))

    def test_from_caps_multiple(self):
        c1 = s2.S2Cap(s2.S2Point(1.0, 0.0, 0.0), s2.S1Angle.from_degrees(10.0))
        c2 = s2.S2Cap(s2.S2Point(0.0, 1.0, 0.0), s2.S1Angle.from_degrees(10.0))
        result = s2.S2Cap.from_caps([c1, c2])
        self.assertTrue(result.contains_point(c1.center))
        self.assertTrue(result.contains_point(c2.center))

    # Properties

    def test_center(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap.from_point(center)
        self.assertEqual(cap.center, center)

    def test_radius(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        chord = s2.S1ChordAngle(s2.S1Angle.from_radians(1.0))
        cap = s2.S2Cap(center, chord)
        self.assertEqual(cap.radius, chord)

    def test_height(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap.from_center_height(center, 0.5)
        self.assertAlmostEqual(cap.height, 0.5)

    # Predicates

    def test_is_empty(self):
        self.assertTrue(s2.S2Cap.empty().is_empty())
        self.assertFalse(s2.S2Cap.full().is_empty())

    def test_is_full(self):
        self.assertTrue(s2.S2Cap.full().is_full())
        self.assertFalse(s2.S2Cap.empty().is_full())

    # Geometric operations

    def test_radius_angle(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        self.assertAlmostEqual(cap.radius_angle().radians, 1.0, places=14)

    def test_area(self):
        self.assertAlmostEqual(s2.S2Cap.full().area(), 4 * math.pi)
        self.assertAlmostEqual(s2.S2Cap.empty().area(), 0.0)

    def test_centroid(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap.from_center_height(center, 0.5)
        centroid = cap.centroid()
        # Centroid lies along the same axis as center.
        self.assertGreater(centroid.x, 0.0)
        self.assertAlmostEqual(centroid.y, 0.0)
        self.assertAlmostEqual(centroid.z, 0.0)

    def test_complement_of_empty_is_full(self):
        self.assertTrue(s2.S2Cap.empty().complement().is_full())

    def test_complement_of_full_is_empty(self):
        self.assertTrue(s2.S2Cap.full().complement().is_empty())

    def test_complement_roundtrip(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        comp = cap.complement()
        self.assertFalse(comp.is_empty())
        self.assertFalse(comp.is_full())

    def test_expanded(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap.from_point(center)
        delta = s2.S1Angle.from_radians(0.1)
        expanded = cap.expanded(delta)
        self.assertGreater(expanded.radius_angle().radians,
                           cap.radius_angle().radians)

    def test_expanded_empty_stays_empty(self):
        self.assertTrue(
            s2.S2Cap.empty().expanded(s2.S1Angle.from_radians(1.0)).is_empty())

    def test_union_contains_both_centers(self):
        p1 = s2.S2Point(1.0, 0.0, 0.0)
        p2 = s2.S2Point(math.cos(0.1), math.sin(0.1), 0.0)
        cap1 = s2.S2Cap(p1, s2.S1Angle.from_radians(0.01))
        cap2 = s2.S2Cap(p2, s2.S1Angle.from_radians(0.01))
        u = cap1.union(cap2)
        # Centers lie strictly inside the union radius, so contains_point is reliable.
        self.assertTrue(u.contains_point(p1))
        self.assertTrue(u.contains_point(p2))

    def test_full_contains_any_cap(self):
        cap = s2.S2Cap(s2.S2Point(1.0, 0.0, 0.0), s2.S1Angle.from_radians(1.0))
        self.assertTrue(s2.S2Cap.full().contains(cap))

    def test_empty_contained_by_any_cap(self):
        cap = s2.S2Cap(s2.S2Point(1.0, 0.0, 0.0), s2.S1Angle.from_radians(1.0))
        self.assertTrue(cap.contains(s2.S2Cap.empty()))

    def test_contains_point(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        self.assertTrue(cap.contains_point(center))
        self.assertFalse(cap.contains_point(s2.S2Point(-1.0, 0.0, 0.0)))

    def test_contains_cell(self):
        cell = s2.S2Cell(s2.S2CellId(s2.S2Point(1.0, 0.0, 0.0)))
        self.assertTrue(s2.S2Cap.full().contains_cell(cell))
        self.assertFalse(s2.S2Cap.empty().contains_cell(cell))

    def test_intersects(self):
        cap = s2.S2Cap(s2.S2Point(1.0, 0.0, 0.0), s2.S1Angle.from_radians(1.0))
        self.assertTrue(cap.intersects(cap))
        self.assertFalse(cap.intersects(s2.S2Cap.empty()))

    def test_interior_intersects(self):
        cap = s2.S2Cap(s2.S2Point(1.0, 0.0, 0.0), s2.S1Angle.from_radians(1.0))
        self.assertTrue(cap.interior_intersects(cap))

    def test_interior_contains_point(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        self.assertTrue(cap.interior_contains_point(center))

    def test_may_intersect(self):
        cell = s2.S2Cell(s2.S2CellId(s2.S2Point(1.0, 0.0, 0.0)))
        self.assertTrue(s2.S2Cap.full().may_intersect(cell))

    def test_cap_bound(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        bound = cap.cap_bound()
        self.assertTrue(bound.approx_equals(cap))

    def test_cell_union_bound(self):
        cell_ids = s2.S2Cap.full().cell_union_bound()
        self.assertGreater(len(cell_ids), 0)

    # Operators

    def test_equality(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap1 = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        cap2 = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        self.assertTrue(cap1 == cap2)
        self.assertFalse(cap1 != cap2)

    def test_inequality(self):
        cap1 = s2.S2Cap.from_point(s2.S2Point(1.0, 0.0, 0.0))
        cap2 = s2.S2Cap.from_point(s2.S2Point(0.0, 1.0, 0.0))
        self.assertTrue(cap1 != cap2)
        self.assertFalse(cap1 == cap2)

    def test_approx_equals(self):
        cap = s2.S2Cap(s2.S2Point(1.0, 0.0, 0.0), s2.S1Angle.from_radians(1.0))
        self.assertTrue(cap.approx_equals(cap))

    def test_approx_equals_with_max_error(self):
        cap = s2.S2Cap(s2.S2Point(1.0, 0.0, 0.0), s2.S1Angle.from_radians(1.0))
        self.assertTrue(cap.approx_equals(cap, s2.S1Angle.from_radians(1e-10)))

    def test_approx_equals_false(self):
        cap1 = s2.S2Cap(s2.S2Point(1.0, 0.0, 0.0), s2.S1Angle.from_radians(1.0))
        cap2 = s2.S2Cap(s2.S2Point(0.0, 1.0, 0.0), s2.S1Angle.from_radians(1.0))
        self.assertFalse(cap1.approx_equals(cap2))

    def test_hash(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap1 = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        cap2 = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        self.assertEqual(hash(cap1), hash(cap2))
        s = {cap1, cap2}
        self.assertEqual(len(s), 1)

    # String representation

    def test_repr(self):
        cap = s2.S2Cap.empty()
        r = repr(cap)
        self.assertTrue(r.startswith("S2Cap("))
        self.assertIn("Center=", r)
        self.assertIn("Radius=", r)

    def test_str(self):
        cap = s2.S2Cap.empty()
        s = str(cap)
        self.assertIn("Center=", s)
        self.assertIn("Radius=", s)


if __name__ == "__main__":
    unittest.main()
