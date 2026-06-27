"""Tests for S2Cap pybind11 bindings."""

import math
import unittest
import s2geometry_pybind as s2


class TestS2Cap(unittest.TestCase):

    # --- Constructors ---

    def test_default_constructor_is_empty(self):
        cap = s2.S2Cap()
        self.assertTrue(cap.is_empty())
        self.assertFalse(cap.is_full())

    def test_constructor_s1angle(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        radius = s2.S1Angle.from_radians(math.pi / 4)
        cap = s2.S2Cap(center, radius)
        self.assertTrue(cap.is_valid())
        self.assertFalse(cap.is_empty())
        self.assertFalse(cap.is_full())

    def test_constructor_s1chord_angle(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        radius = s2.S1ChordAngle(s2.S1Angle.from_radians(math.pi / 4))
        cap = s2.S2Cap(center, radius)
        self.assertTrue(cap.is_valid())

    def test_constructor_negative_radius_is_empty(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(-1.0))
        self.assertTrue(cap.is_empty())

    def test_constructor_pi_radius_is_full(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(math.pi))
        self.assertTrue(cap.is_full())

    # --- Static factories ---

    def test_from_point(self):
        p = s2.S2Point(0.0, 1.0, 0.0)
        cap = s2.S2Cap.from_point(p)
        self.assertTrue(cap.contains_point(p))
        self.assertTrue(cap.is_valid())

    def test_from_center_height(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap.from_center_height(center, 0.5)
        self.assertTrue(cap.is_valid())
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
        self.assertTrue(cap.is_valid())
        self.assertAlmostEqual(cap.get_area(), area, places=10)

    def test_empty(self):
        cap = s2.S2Cap.empty()
        self.assertTrue(cap.is_empty())
        self.assertFalse(cap.is_full())

    def test_full(self):
        cap = s2.S2Cap.full()
        self.assertTrue(cap.is_full())
        self.assertFalse(cap.is_empty())

    # --- Properties ---

    def test_center(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap.from_point(center)
        self.assertEqual(cap.center, center)

    def test_radius_type(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        chord = s2.S1ChordAngle(s2.S1Angle.from_radians(1.0))
        cap = s2.S2Cap(center, chord)
        self.assertIsInstance(cap.radius, s2.S1ChordAngle)

    def test_height(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap.from_center_height(center, 0.5)
        self.assertAlmostEqual(cap.height, 0.5)

    def test_get_radius(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        angle = s2.S1Angle.from_radians(1.0)
        cap = s2.S2Cap(center, angle)
        self.assertIsInstance(cap.get_radius(), s2.S1Angle)
        self.assertAlmostEqual(cap.get_radius().radians, 1.0, places=14)

    def test_get_area(self):
        # Full cap has area 4*pi
        self.assertAlmostEqual(s2.S2Cap.full().get_area(), 4 * math.pi)
        # Empty cap has area 0
        self.assertAlmostEqual(s2.S2Cap.empty().get_area(), 0.0)

    def test_get_centroid(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap.from_center_height(center, 0.5)
        centroid = cap.get_centroid()
        self.assertIsInstance(centroid, s2.S2Point)
        # Centroid lies along the same axis as center
        self.assertGreater(centroid.x, 0.0)
        self.assertAlmostEqual(centroid.y, 0.0)
        self.assertAlmostEqual(centroid.z, 0.0)

    # --- Predicates ---

    def test_is_valid(self):
        self.assertTrue(s2.S2Cap.empty().is_valid())
        self.assertTrue(s2.S2Cap.full().is_valid())
        center = s2.S2Point(1.0, 0.0, 0.0)
        self.assertTrue(s2.S2Cap(center, s2.S1Angle.from_radians(1.0)).is_valid())

    def test_is_empty(self):
        self.assertTrue(s2.S2Cap.empty().is_empty())
        self.assertFalse(s2.S2Cap.full().is_empty())

    def test_is_full(self):
        self.assertTrue(s2.S2Cap.full().is_full())
        self.assertFalse(s2.S2Cap.empty().is_full())

    # --- Geometric operations ---

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
        self.assertGreater(expanded.get_radius().radians,
                           cap.get_radius().radians)

    def test_expanded_empty_stays_empty(self):
        self.assertTrue(
            s2.S2Cap.empty().expanded(s2.S1Angle.from_radians(1.0)).is_empty())

    def test_union_contains_both_centers(self):
        p1 = s2.S2Point(1.0, 0.0, 0.0)
        p2 = s2.S2Point(math.cos(0.1), math.sin(0.1), 0.0)
        cap1 = s2.S2Cap(p1, s2.S1Angle.from_radians(0.01))
        cap2 = s2.S2Cap(p2, s2.S1Angle.from_radians(0.01))
        u = cap1.union(cap2)
        # The centers of both input caps are strictly inside the union radius,
        # so contains_point is reliable here.
        self.assertTrue(u.contains_point(p1))
        self.assertTrue(u.contains_point(p2))

    # --- Containment / intersection ---

    def test_full_contains_any_cap(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        self.assertTrue(s2.S2Cap.full().contains(cap))

    def test_empty_contained_by_any_cap(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        self.assertTrue(cap.contains(s2.S2Cap.empty()))

    def test_contains_point(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        self.assertTrue(cap.contains_point(center))
        antipode = s2.S2Point(-1.0, 0.0, 0.0)
        self.assertFalse(cap.contains_point(antipode))

    def test_contains_cell(self):
        cell = s2.S2Cell(s2.S2CellId(s2.S2Point(1.0, 0.0, 0.0)))
        self.assertTrue(s2.S2Cap.full().contains_cell(cell))
        self.assertFalse(s2.S2Cap.empty().contains_cell(cell))

    def test_intersects(self):
        p = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(p, s2.S1Angle.from_radians(1.0))
        self.assertTrue(cap.intersects(cap))
        self.assertFalse(cap.intersects(s2.S2Cap.empty()))

    def test_interior_intersects(self):
        p = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(p, s2.S1Angle.from_radians(1.0))
        self.assertTrue(cap.interior_intersects(cap))

    def test_interior_contains_point(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        self.assertTrue(cap.interior_contains_point(center))

    def test_may_intersect_cell(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cell = s2.S2Cell(s2.S2CellId(center))
        self.assertTrue(s2.S2Cap.full().may_intersect(cell))

    # --- Mutation ---

    def test_add_point(self):
        cap = s2.S2Cap.empty()
        p = s2.S2Point(1.0, 0.0, 0.0)
        cap.add_point(p)
        self.assertTrue(cap.contains_point(p))

    def test_add_cap(self):
        cap1 = s2.S2Cap.empty()
        p = s2.S2Point(1.0, 0.0, 0.0)
        cap2 = s2.S2Cap.from_point(p)
        cap1.add_cap(cap2)
        self.assertTrue(cap1.contains_point(p))

    # --- S2Region interface ---

    def test_get_cap_bound(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        bound = cap.get_cap_bound()
        self.assertIsInstance(bound, s2.S2Cap)
        self.assertTrue(bound.approx_equals(cap))

    def test_cell_union_bound(self):
        cap = s2.S2Cap.full()
        cell_ids = cap.cell_union_bound()
        self.assertIsInstance(cell_ids, list)
        self.assertGreater(len(cell_ids), 0)
        for cid in cell_ids:
            self.assertIsInstance(cid, s2.S2CellId)

    # --- Operators ---

    def test_equality(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap1 = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        cap2 = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        self.assertEqual(cap1, cap2)

    def test_inequality(self):
        cap1 = s2.S2Cap.from_point(s2.S2Point(1.0, 0.0, 0.0))
        cap2 = s2.S2Cap.from_point(s2.S2Point(0.0, 1.0, 0.0))
        self.assertNotEqual(cap1, cap2)

    def test_approx_equals(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        self.assertTrue(cap.approx_equals(cap))

    def test_approx_equals_with_max_error(self):
        center = s2.S2Point(1.0, 0.0, 0.0)
        cap = s2.S2Cap(center, s2.S1Angle.from_radians(1.0))
        self.assertTrue(cap.approx_equals(cap, s2.S1Angle.from_radians(1e-10)))

    # --- String representation ---

    def test_repr(self):
        cap = s2.S2Cap.empty()
        r = repr(cap)
        self.assertIn("Center=", r)
        self.assertIn("Radius=", r)

    def test_str(self):
        cap = s2.S2Cap.empty()
        self.assertIsInstance(str(cap), str)

    # --- S2Cell.get_cap_bound ---

    def test_s2cell_get_cap_bound(self):
        cell = s2.S2Cell(s2.S2CellId.from_face(0))
        cap = cell.get_cap_bound()
        self.assertIsInstance(cap, s2.S2Cap)
        self.assertTrue(cap.is_valid())
        # The cap must contain the cell's center.
        self.assertTrue(cap.contains_point(cell.center()))


if __name__ == "__main__":
    unittest.main()
