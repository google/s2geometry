"""Tests for S2LatLngRect pybind11 bindings."""

import math
import unittest
import s2geometry_pybind as s2


class TestS2LatLngRect(unittest.TestCase):

    # --- Constructors ---

    def test_default_constructor_is_empty(self):
        rect = s2.S2LatLngRect()
        self.assertTrue(rect.is_empty())
        self.assertFalse(rect.is_full())

    def test_constructor_from_lo_hi(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertTrue(rect.is_valid())
        self.assertFalse(rect.is_empty())

    def test_constructor_from_intervals(self):
        lat = s2.R1Interval(-0.5, 0.5)
        lng = s2.S1Interval(-1.0, 1.0)
        rect = s2.S2LatLngRect(lat, lng)
        self.assertTrue(rect.is_valid())

    # --- Static factories ---

    def test_from_center_size(self):
        center = s2.S2LatLng.from_degrees(10.0, 20.0)
        size = s2.S2LatLng.from_degrees(4.0, 6.0)
        rect = s2.S2LatLngRect.from_center_size(center, size)
        self.assertTrue(rect.is_valid())
        self.assertTrue(rect.contains_latlng(center))

    def test_from_point(self):
        p = s2.S2LatLng.from_degrees(45.0, 90.0)
        rect = s2.S2LatLngRect.from_point(p)
        self.assertTrue(rect.is_point())
        self.assertTrue(rect.contains_latlng(p))

    def test_from_point_pair(self):
        p1 = s2.S2LatLng.from_degrees(-10.0, -20.0)
        p2 = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect.from_point_pair(p1, p2)
        self.assertTrue(rect.contains_latlng(p1))
        self.assertTrue(rect.contains_latlng(p2))

    def test_empty(self):
        self.assertTrue(s2.S2LatLngRect.empty().is_empty())
        self.assertFalse(s2.S2LatLngRect.empty().is_full())

    def test_full(self):
        self.assertTrue(s2.S2LatLngRect.full().is_full())
        self.assertFalse(s2.S2LatLngRect.full().is_empty())

    def test_full_lat(self):
        lat = s2.S2LatLngRect.full_lat()
        self.assertIsInstance(lat, s2.R1Interval)
        self.assertAlmostEqual(lat.lo, -math.pi / 2)
        self.assertAlmostEqual(lat.hi, math.pi / 2)

    def test_full_lng(self):
        lng = s2.S2LatLngRect.full_lng()
        self.assertIsInstance(lng, s2.S1Interval)
        self.assertTrue(lng.is_full())

    # --- Properties ---

    def test_lat_property(self):
        lo = s2.S2LatLng.from_degrees(-10.0, 0.0)
        hi = s2.S2LatLng.from_degrees(10.0, 0.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertIsInstance(rect.lat, s2.R1Interval)
        self.assertAlmostEqual(rect.lat.lo, lo.lat.radians)
        self.assertAlmostEqual(rect.lat.hi, hi.lat.radians)

    def test_lng_property(self):
        lo = s2.S2LatLng.from_degrees(0.0, -30.0)
        hi = s2.S2LatLng.from_degrees(0.0, 30.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertIsInstance(rect.lng, s2.S1Interval)

    def test_lat_lo_lat_hi(self):
        lo = s2.S2LatLng.from_degrees(-10.0, 0.0)
        hi = s2.S2LatLng.from_degrees(20.0, 0.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertAlmostEqual(rect.lat_lo.degrees, -10.0)
        self.assertAlmostEqual(rect.lat_hi.degrees, 20.0)

    def test_lng_lo_lng_hi(self):
        lo = s2.S2LatLng.from_degrees(0.0, -30.0)
        hi = s2.S2LatLng.from_degrees(0.0, 60.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertAlmostEqual(rect.lng_lo.degrees, -30.0)
        self.assertAlmostEqual(rect.lng_hi.degrees, 60.0)

    def test_lo_hi(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertAlmostEqual(rect.lo.lat.degrees, -10.0)
        self.assertAlmostEqual(rect.hi.lat.degrees, 10.0)

    # --- Predicates ---

    def test_is_valid(self):
        self.assertTrue(s2.S2LatLngRect.empty().is_valid())
        self.assertTrue(s2.S2LatLngRect.full().is_valid())

    def test_is_empty(self):
        self.assertTrue(s2.S2LatLngRect.empty().is_empty())
        self.assertFalse(s2.S2LatLngRect.full().is_empty())

    def test_is_full(self):
        self.assertTrue(s2.S2LatLngRect.full().is_full())
        self.assertFalse(s2.S2LatLngRect.empty().is_full())

    def test_is_point(self):
        p = s2.S2LatLng.from_degrees(45.0, 90.0)
        self.assertTrue(s2.S2LatLngRect.from_point(p).is_point())
        self.assertFalse(s2.S2LatLngRect.full().is_point())

    def test_is_inverted(self):
        # A rectangle spanning the 180-degree line has lo.lng > hi.lng
        lo = s2.S2LatLng.from_degrees(0.0, 170.0)
        hi = s2.S2LatLng.from_degrees(10.0, -170.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertTrue(rect.is_inverted())

    def test_inverted_rect(self):
        # A rect where lng_lo > lng_hi spans the antimeridian (is "inverted").
        lo = s2.S2LatLng.from_degrees(-10, 160)
        hi = s2.S2LatLng.from_degrees(10, -160)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertTrue(rect.is_inverted())
        # A point at lng=175 is inside (between 160 and -160 via antimeridian).
        self.assertTrue(rect.contains_latlng(s2.S2LatLng.from_degrees(0, 175)))
        # A point at lng=0 is outside.
        self.assertFalse(rect.contains_latlng(s2.S2LatLng.from_degrees(0, 0)))

    # --- Geometric accessors ---

    def test_vertex(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        for k in range(4):
            v = rect.vertex(k)
            self.assertIsInstance(v, s2.S2LatLng)

    def test_center(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        c = rect.center()
        self.assertAlmostEqual(c.lat.degrees, 0.0)
        self.assertAlmostEqual(c.lng.degrees, 0.0)

    def test_size(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        sz = rect.size()
        self.assertAlmostEqual(sz.lat.degrees, 20.0)
        self.assertAlmostEqual(sz.lng.degrees, 40.0)

    def test_area(self):
        self.assertGreater(s2.S2LatLngRect.full().area(), 0.0)
        self.assertAlmostEqual(s2.S2LatLngRect.empty().area(), 0.0)

    def test_centroid(self):
        rect = s2.S2LatLngRect.full()
        centroid = rect.centroid()
        self.assertIsInstance(centroid, s2.S2Point)

    # --- Containment ---

    def test_contains_rect(self):
        outer = s2.S2LatLngRect.full()
        inner = s2.S2LatLngRect.from_point(s2.S2LatLng.from_degrees(0.0, 0.0))
        self.assertTrue(outer.contains(inner))
        self.assertFalse(inner.contains(outer))

    def test_contains_latlng(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertTrue(rect.contains_latlng(s2.S2LatLng.from_degrees(0.0, 0.0)))
        self.assertFalse(rect.contains_latlng(s2.S2LatLng.from_degrees(45.0, 0.0)))

    def test_contains_point(self):
        rect = s2.S2LatLngRect.full()
        self.assertTrue(rect.contains_point(s2.S2Point(1.0, 0.0, 0.0)))

    def test_contains_cell(self):
        rect = s2.S2LatLngRect.full()
        cell = s2.S2Cell(s2.S2CellId.from_face(0))
        self.assertTrue(rect.contains_cell(cell))

    def test_interior_contains_point(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertTrue(rect.interior_contains_point(s2.S2Point(1.0, 0.0, 0.0)))

    def test_interior_contains_latlng(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertTrue(
            rect.interior_contains_latlng(s2.S2LatLng.from_degrees(0.0, 0.0)))
        self.assertFalse(
            rect.interior_contains_latlng(s2.S2LatLng.from_degrees(10.0, 0.0)))

    def test_interior_contains_rect(self):
        outer = s2.S2LatLngRect.full()
        inner = s2.S2LatLngRect.from_point(s2.S2LatLng.from_degrees(0.0, 0.0))
        self.assertTrue(outer.interior_contains(inner))

    # --- Intersection ---

    def test_intersects_rect(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertTrue(rect.intersects(rect))
        self.assertFalse(rect.intersects(s2.S2LatLngRect.empty()))

    def test_intersects_cell(self):
        rect = s2.S2LatLngRect.full()
        cell = s2.S2Cell(s2.S2CellId.from_face(0))
        self.assertTrue(rect.intersects_cell(cell))

    def test_interior_intersects(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertTrue(rect.interior_intersects(rect))

    def test_boundary_intersects(self):
        rect = s2.S2LatLngRect.full()
        p1 = s2.S2Point(1.0, 0.0, 0.0)
        p2 = s2.S2Point(0.0, 1.0, 0.0)
        # Full rect has no boundary on the sphere, so this returns False.
        self.assertFalse(rect.boundary_intersects(p1, p2))

        rect2 = s2.S2LatLngRect(s2.S2LatLng.from_degrees(-10, -20),
                                 s2.S2LatLng.from_degrees(10, 20))
        # An edge crossing the boundary: one point inside, one outside.
        p_in = s2.S2LatLng.from_degrees(0, 0).to_point()
        p_out = s2.S2LatLng.from_degrees(0, 45).to_point()
        self.assertTrue(rect2.boundary_intersects(p_in, p_out))
        # An edge entirely inside: both points inside.
        p_a = s2.S2LatLng.from_degrees(-5, -10).to_point()
        p_b = s2.S2LatLng.from_degrees(5, 10).to_point()
        self.assertFalse(rect2.boundary_intersects(p_a, p_b))

    def test_may_intersect(self):
        rect = s2.S2LatLngRect.full()
        cell = s2.S2Cell(s2.S2CellId.from_face(0))
        self.assertTrue(rect.may_intersect(cell))

    # --- Mutation ---

    def test_add_point_s2point(self):
        rect = s2.S2LatLngRect.empty()
        rect.add_point(s2.S2Point(1.0, 0.0, 0.0))
        self.assertFalse(rect.is_empty())

    def test_add_point_latlng(self):
        rect = s2.S2LatLngRect.empty()
        ll = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect.add_point(ll)
        self.assertTrue(rect.contains_latlng(ll))

    # --- Set operations ---

    def test_expanded(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        margin = s2.S2LatLng.from_degrees(5.0, 5.0)
        expanded = rect.expanded(margin)
        self.assertTrue(expanded.contains(rect))

    def test_polar_closure_of_full_is_full(self):
        self.assertTrue(s2.S2LatLngRect.full().polar_closure().is_full())

    def test_union(self):
        p1 = s2.S2LatLng.from_degrees(-10.0, -10.0)
        p2 = s2.S2LatLng.from_degrees(10.0, 10.0)
        r1 = s2.S2LatLngRect.from_point(p1)
        r2 = s2.S2LatLngRect.from_point(p2)
        u = r1.union(r2)
        self.assertTrue(u.contains_latlng(p1))
        self.assertTrue(u.contains_latlng(p2))

    def test_intersection(self):
        lo1 = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi1 = s2.S2LatLng.from_degrees(10.0, 20.0)
        lo2 = s2.S2LatLng.from_degrees(0.0, 0.0)
        hi2 = s2.S2LatLng.from_degrees(20.0, 30.0)
        r1 = s2.S2LatLngRect(lo1, hi1)
        r2 = s2.S2LatLngRect(lo2, hi2)
        inter = r1.intersection(r2)
        self.assertTrue(r1.contains(inter))
        self.assertTrue(r2.contains(inter))

    def test_expanded_by_distance(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        distance = s2.S1Angle.from_radians(0.1)
        expanded = rect.expanded_by_distance(distance)
        self.assertTrue(expanded.contains(rect))

    # --- Distance ---

    def test_distance_to_rect(self):
        lo1 = s2.S2LatLng.from_degrees(10.0, 0.0)
        hi1 = s2.S2LatLng.from_degrees(20.0, 10.0)
        lo2 = s2.S2LatLng.from_degrees(30.0, 0.0)
        hi2 = s2.S2LatLng.from_degrees(40.0, 10.0)
        r1 = s2.S2LatLngRect(lo1, hi1)
        r2 = s2.S2LatLngRect(lo2, hi2)
        d = r1.distance_to_rect(r2)
        self.assertIsInstance(d, s2.S1Angle)
        self.assertGreater(d.radians, 0.0)

    def test_distance_to_latlng(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -10.0)
        hi = s2.S2LatLng.from_degrees(10.0, 10.0)
        rect = s2.S2LatLngRect(lo, hi)
        p = s2.S2LatLng.from_degrees(0.0, 0.0)
        d = rect.distance_to_latlng(p)
        self.assertAlmostEqual(d.radians, 0.0)

    def test_directed_hausdorff_distance(self):
        r1 = s2.S2LatLngRect.full()
        lo = s2.S2LatLng.from_degrees(-10.0, -10.0)
        hi = s2.S2LatLng.from_degrees(10.0, 10.0)
        r2 = s2.S2LatLngRect(lo, hi)
        d = r2.directed_hausdorff_distance(r1)
        self.assertIsInstance(d, s2.S1Angle)

    def test_hausdorff_distance(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -10.0)
        hi = s2.S2LatLng.from_degrees(10.0, 10.0)
        rect = s2.S2LatLngRect(lo, hi)
        d = rect.hausdorff_distance(rect)
        self.assertAlmostEqual(d.radians, 0.0)

    # --- S2Region interface ---

    def test_rect_bound(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        bound = rect.rect_bound()
        self.assertIsInstance(bound, s2.S2LatLngRect)
        self.assertTrue(bound.approx_equals(rect))

    def test_cell_union_bound(self):
        rect = s2.S2LatLngRect.full()
        cell_ids = rect.cell_union_bound()
        self.assertIsInstance(cell_ids, list)
        self.assertGreater(len(cell_ids), 0)
        for cid in cell_ids:
            self.assertIsInstance(cid, s2.S2CellId)

    # --- Operators ---

    def test_equality(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        r1 = s2.S2LatLngRect(lo, hi)
        r2 = s2.S2LatLngRect(lo, hi)
        self.assertEqual(r1, r2)

    def test_inequality(self):
        r1 = s2.S2LatLngRect.empty()
        r2 = s2.S2LatLngRect.full()
        self.assertNotEqual(r1, r2)

    def test_approx_equals(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertTrue(rect.approx_equals(rect))

    def test_approx_equals_latlng(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        tol = s2.S2LatLng.from_degrees(1e-10, 1e-10)
        self.assertTrue(rect.approx_equals_latlng(rect, tol))

    def test_hash(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        r1 = s2.S2LatLngRect(lo, hi)
        r2 = s2.S2LatLngRect(lo, hi)
        self.assertEqual(hash(r1), hash(r2))

    # --- String representation ---

    def test_repr(self):
        rect = s2.S2LatLngRect.empty()
        r = repr(rect)
        self.assertIn("Lo", r)
        self.assertIn("Hi", r)
        self.assertTrue(repr(rect).startswith("S2LatLngRect("))

    def test_str(self):
        self.assertIsInstance(str(s2.S2LatLngRect.empty()), str)


if __name__ == "__main__":
    unittest.main()
