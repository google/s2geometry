"""Tests for S2LatLngRect pybind11 bindings."""

import math
import unittest
import s2geometry_pybind as s2


class TestS2LatLngRect(unittest.TestCase):
    """Test cases for S2LatLngRect bindings."""

    # Constructors

    def test_default_constructor_is_empty(self):
        rect = s2.S2LatLngRect()
        self.assertTrue(rect.is_empty())
        self.assertFalse(rect.is_full())

    def test_constructor_from_lo_hi(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertEqual(rect.lo, lo)
        self.assertEqual(rect.hi, hi)

    def test_constructor_rejects_reversed_lat(self):
        lo = s2.S2LatLng.from_degrees(10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(-10.0, 20.0)
        with self.assertRaises(ValueError):
            s2.S2LatLngRect(lo, hi)

    def test_constructor_from_intervals(self):
        lat = s2.R1Interval(-0.5, 0.5)
        lng = s2.S1Interval(-1.0, 1.0)
        rect = s2.S2LatLngRect(lat, lng)
        self.assertEqual(rect.lat, lat)
        self.assertEqual(rect.lng, lng)

    def test_constructor_from_both_empty_intervals_is_canonical(self):
        # Non-canonical empty R1Interval must produce canonical Empty(), so
        # that == and __hash__ are consistent.
        rect = s2.S2LatLngRect(s2.R1Interval(0.5, -0.5), s2.S1Interval.empty())
        self.assertEqual(rect, s2.S2LatLngRect.empty())
        self.assertEqual(hash(rect), hash(s2.S2LatLngRect.empty()))

    def test_constructor_rejects_lo_lng_180_hi_lng_neg180(self):
        lo = s2.S2LatLng.from_degrees(0.0, 180.0)
        hi = s2.S2LatLng.from_degrees(10.0, -180.0)
        with self.assertRaises(ValueError):
            s2.S2LatLngRect(lo, hi)

    def test_constructor_rejects_lat_out_of_range(self):
        lat = s2.R1Interval(0.0, 2.0)  # hi > pi/2 (~1.5708)
        lng = s2.S1Interval(0.0, 1.0)
        with self.assertRaises(ValueError):
            s2.S2LatLngRect(lat, lng)

    def test_constructor_rejects_nan_lat(self):
        lat = s2.R1Interval(float('nan'), float('nan'))
        lng = s2.S1Interval(0.0, 1.0)
        with self.assertRaises(ValueError):
            s2.S2LatLngRect(lat, lng)

    def test_constructor_rejects_empty_mismatch(self):
        lat = s2.R1Interval.empty()
        lng = s2.S1Interval(0.0, 1.0)  # non-empty lng with empty lat
        with self.assertRaises(ValueError):
            s2.S2LatLngRect(lat, lng)

    # Factory methods

    def test_from_center_size(self):
        center = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect.from_center_size(
            center, s2.S1Angle.from_degrees(4.0), s2.S1Angle.from_degrees(6.0))
        self.assertTrue(rect.contains_latlng(center))
        self.assertAlmostEqual(rect.lat_lo.degrees, 8.0)
        self.assertAlmostEqual(rect.lat_hi.degrees, 12.0)
        self.assertAlmostEqual(rect.lng_lo.degrees, 17.0)
        self.assertAlmostEqual(rect.lng_hi.degrees, 23.0)

    def test_from_center_size_exceeds_latlng_validity_range(self):
        # Example from s2latlng_rect.h: sizes beyond 180°/360° clamp to full.
        center = s2.S2LatLng.from_degrees(10.0, 40.0)
        rect = s2.S2LatLngRect.from_center_size(
            center, s2.S1Angle.from_degrees(210.0),
            s2.S1Angle.from_degrees(400.0))
        self.assertTrue(rect.is_full())

    def test_from_center_size_round_trips_size(self):
        rect = s2.S2LatLngRect(s2.S2LatLng.from_degrees(-80.0, -170.0),
                               s2.S2LatLng.from_degrees(70.0, 150.0))
        rebuilt = s2.S2LatLngRect.from_center_size(rect.center(), *rect.size())
        self.assertTrue(rebuilt.approx_equals(rect, s2.S1Angle.from_radians(1e-14)))

    def test_from_center_size_rejects_nan(self):
        center = s2.S2LatLng.from_degrees(0.0, 0.0)
        nan = s2.S1Angle.from_radians(math.nan)
        one = s2.S1Angle.from_degrees(1.0)
        with self.assertRaises(ValueError):
            s2.S2LatLngRect.from_center_size(center, nan, one)
        with self.assertRaises(ValueError):
            s2.S2LatLngRect.from_center_size(center, one, nan)

    def test_from_point(self):
        p = s2.S2LatLng.from_degrees(45.0, 90.0)
        rect = s2.S2LatLngRect.from_point(p)
        self.assertTrue(rect.is_point())
        self.assertTrue(rect.contains_latlng(p))

    def test_from_point_pair(self):
        # Ordered points.
        p1 = s2.S2LatLng.from_degrees(-10.0, -20.0)
        p2 = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect.from_point_pair(p1, p2)
        self.assertTrue(rect.contains_latlng(p1))
        self.assertTrue(rect.contains_latlng(p2))
        # Reversed order must produce the same result.
        self.assertEqual(s2.S2LatLngRect.from_point_pair(p2, p1), rect)

    def test_from_latlngs_empty(self):
        rect = s2.S2LatLngRect.from_latlngs([])
        self.assertTrue(rect.is_empty())

    def test_from_latlngs_single(self):
        ll = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect.from_latlngs([ll])
        self.assertTrue(rect.contains_latlng(ll))
        self.assertTrue(rect.is_point())

    def test_from_latlngs_multiple(self):
        lls = [
            s2.S2LatLng.from_degrees(-10.0, -20.0),
            s2.S2LatLng.from_degrees(10.0, 20.0),
            s2.S2LatLng.from_degrees(5.0, 0.0),
        ]
        rect = s2.S2LatLngRect.from_latlngs(lls)
        for ll in lls:
            self.assertTrue(rect.contains_latlng(ll))

    def test_empty(self):
        self.assertTrue(s2.S2LatLngRect.empty().is_empty())
        self.assertFalse(s2.S2LatLngRect.empty().is_full())

    def test_full(self):
        self.assertTrue(s2.S2LatLngRect.full().is_full())
        self.assertFalse(s2.S2LatLngRect.full().is_empty())

    def test_full_lat(self):
        lat = s2.S2LatLngRect.full_lat()
        self.assertAlmostEqual(lat.lo, -math.pi / 2)
        self.assertAlmostEqual(lat.hi, math.pi / 2)

    def test_full_lng(self):
        lng = s2.S2LatLngRect.full_lng()
        self.assertTrue(lng.is_full())

    # Properties

    def test_lat_property(self):
        lo = s2.S2LatLng.from_degrees(-10.0, 0.0)
        hi = s2.S2LatLng.from_degrees(10.0, 0.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertAlmostEqual(rect.lat.lo, lo.lat.radians)
        self.assertAlmostEqual(rect.lat.hi, hi.lat.radians)

    def test_lng_property(self):
        lo = s2.S2LatLng.from_degrees(0.0, -30.0)
        hi = s2.S2LatLng.from_degrees(0.0, 30.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertAlmostEqual(rect.lng.lo, lo.lng.radians)
        self.assertAlmostEqual(rect.lng.hi, hi.lng.radians)

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
        self.assertEqual(rect.lo, lo)
        self.assertEqual(rect.hi, hi)

    # Predicates

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
        # A point at lng=175 is inside (between 160 and -160 via antimeridian).
        self.assertTrue(rect.contains_latlng(s2.S2LatLng.from_degrees(0, 175)))
        # A point at lng=0 is outside.
        self.assertFalse(rect.contains_latlng(s2.S2LatLng.from_degrees(0, 0)))

    # Geometric operations

    def test_vertex(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        # Vertices in CCW order: lower-left, lower-right, upper-right, upper-left.
        expected = [(-10.0, -20.0), (-10.0, 20.0), (10.0, 20.0), (10.0, -20.0)]
        for k, (lat, lng) in enumerate(expected):
            v = rect.vertex(k)
            self.assertAlmostEqual(v.lat.degrees, lat)
            self.assertAlmostEqual(v.lng.degrees, lng)

    def test_vertex_modulo_4(self):
        rect = s2.S2LatLngRect.full()
        self.assertEqual(rect.vertex(0), rect.vertex(4))
        self.assertEqual(rect.vertex(3), rect.vertex(-1))

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
        lat_ext, lng_ext = rect.size()
        self.assertAlmostEqual(lat_ext.degrees, 20.0)
        self.assertAlmostEqual(lng_ext.degrees, 40.0)

    def test_size_full_exceeds_latlng_validity_range(self):
        # full().size() returns (180°, 360°), which is outside S2LatLng validity
        # bounds — confirmed by checking the pair of S1Angles directly.
        lat_ext, lng_ext = s2.S2LatLngRect.full().size()
        self.assertAlmostEqual(lat_ext.degrees, 180.0)
        self.assertAlmostEqual(lng_ext.degrees, 360.0)

    def test_area(self):
        self.assertAlmostEqual(s2.S2LatLngRect.full().area(), 4 * math.pi)
        self.assertAlmostEqual(s2.S2LatLngRect.empty().area(), 0.0)

    def test_centroid(self):
        # Full sphere is symmetric; centroid (area-weighted) is the zero vector.
        centroid = s2.S2LatLngRect.full().centroid()
        self.assertAlmostEqual(centroid.x, 0.0)
        self.assertAlmostEqual(centroid.y, 0.0)
        self.assertAlmostEqual(centroid.z, 0.0)

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
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertTrue(rect.contains_point(s2.S2Point(1.0, 0.0, 0.0)))
        # Point on the equator at lng=45° is outside.
        outside = s2.S2LatLng.from_degrees(0.0, 45.0).to_point()
        self.assertFalse(rect.contains_point(outside))

    def test_contains_point_accepts_non_unit_length(self):
        rect = s2.S2LatLngRect(s2.S2LatLng.from_degrees(-10.0, -20.0),
                               s2.S2LatLng.from_degrees(10.0, 20.0))
        self.assertTrue(rect.contains_point(s2.S2Point(2.0, 0.0, 0.0)))
        self.assertFalse(rect.contains_point(s2.S2Point(-2.0, 0.0, 0.0)))

    def test_contains_cell(self):
        rect = s2.S2LatLngRect.full()
        cell = s2.S2Cell(s2.S2CellId.from_face(0))
        self.assertTrue(rect.contains_cell(cell))
        small = s2.S2LatLngRect(s2.S2LatLng.from_degrees(0.0, 0.0),
                                 s2.S2LatLng.from_degrees(1.0, 1.0))
        self.assertFalse(small.contains_cell(cell))

    def test_interior_contains_point(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertTrue(rect.interior_contains_point(s2.S2Point(1.0, 0.0, 0.0)))
        # A point on the boundary is not in the interior.
        boundary = s2.S2LatLng.from_degrees(10.0, 0.0).to_point()
        self.assertFalse(rect.interior_contains_point(boundary))

    def test_interior_contains_point_accepts_non_unit_length(self):
        rect = s2.S2LatLngRect(s2.S2LatLng.from_degrees(-10.0, -20.0),
                               s2.S2LatLng.from_degrees(10.0, 20.0))
        self.assertTrue(rect.interior_contains_point(s2.S2Point(2.0, 0.0, 0.0)))
        self.assertFalse(
            rect.interior_contains_point(s2.S2Point(-2.0, 0.0, 0.0)))

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
        # A rect does not interior-contain itself (boundary touches boundary).
        self.assertFalse(inner.interior_contains(inner))

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
        self.assertFalse(s2.S2LatLngRect.empty().intersects_cell(cell))

    def test_interior_intersects(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertTrue(rect.interior_intersects(rect))
        # A rect whose boundary only touches does not interior-intersect.
        touching = s2.S2LatLngRect(s2.S2LatLng.from_degrees(10.0, -20.0),
                                    s2.S2LatLng.from_degrees(20.0, 20.0))
        self.assertFalse(rect.interior_intersects(touching))

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

    def test_boundary_intersects_rejects_non_unit_length(self):
        rect = s2.S2LatLngRect(s2.S2LatLng.from_degrees(-10, -20),
                               s2.S2LatLng.from_degrees(10, 20))
        unit = s2.S2Point(1.0, 0.0, 0.0)
        non_unit = s2.S2Point(2.0, 0.0, 0.0)
        with self.assertRaises(ValueError):
            rect.boundary_intersects(non_unit, unit)
        with self.assertRaises(ValueError):
            rect.boundary_intersects(unit, non_unit)

    def test_may_intersect(self):
        rect = s2.S2LatLngRect.full()
        cell = s2.S2Cell(s2.S2CellId.from_face(0))
        self.assertTrue(rect.may_intersect(cell))
        self.assertFalse(s2.S2LatLngRect.empty().may_intersect(cell))

    def test_expanded(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        expanded = rect.expanded(s2.S1Angle.from_degrees(5.0),
                                 s2.S1Angle.from_degrees(7.0))
        self.assertTrue(expanded.contains(rect))
        self.assertAlmostEqual(expanded.lat_hi.degrees, 15.0)
        self.assertAlmostEqual(expanded.lng_hi.degrees, 27.0)

    def test_expanded_lat_margin_beyond_90(self):
        # A 180° lat margin grows a south-pole point to the full latitude range.
        south = s2.S2LatLngRect.from_point(s2.S2LatLng.from_degrees(-90.0, 0.0))
        expanded = south.expanded(s2.S1Angle.from_degrees(180.0),
                                  s2.S1Angle.from_degrees(0.0))
        self.assertEqual(expanded.lat, s2.S2LatLngRect.full_lat())

    def test_expanded_shrinks_full_to_empty(self):
        shrunk = s2.S2LatLngRect.full().expanded(s2.S1Angle.from_degrees(-91.0),
                                                 s2.S1Angle.from_degrees(0.0))
        self.assertTrue(shrunk.is_empty())

    def test_expanded_rejects_nan(self):
        rect = s2.S2LatLngRect.full()
        nan = s2.S1Angle.from_radians(math.nan)
        zero = s2.S1Angle.from_degrees(0.0)
        with self.assertRaises(ValueError):
            rect.expanded(nan, zero)
        with self.assertRaises(ValueError):
            rect.expanded(zero, nan)

    def test_polar_closure(self):
        self.assertTrue(s2.S2LatLngRect.full().polar_closure().is_full())
        # A rect touching the north pole with partial lng should expand to full lng.
        polar = s2.S2LatLngRect(s2.S2LatLng.from_degrees(80.0, 10.0),
                                 s2.S2LatLng.from_degrees(90.0, 20.0))
        self.assertTrue(polar.polar_closure().lng.is_full())
        # A rect not touching a pole should be unchanged.
        mid = s2.S2LatLngRect(s2.S2LatLng.from_degrees(0.0, 0.0),
                               s2.S2LatLng.from_degrees(10.0, 10.0))
        self.assertEqual(mid.polar_closure(), mid)

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

    def test_expanded_by_distance_infinity_is_full(self):
        rect = s2.S2LatLngRect.from_point(s2.S2LatLng.from_degrees(10.0, 20.0))
        self.assertTrue(
            rect.expanded_by_distance(s2.S1Angle.infinity()).is_full())

    def test_expanded_by_distance_rejects_nan(self):
        with self.assertRaises(ValueError):
            s2.S2LatLngRect.full().expanded_by_distance(
                s2.S1Angle.from_radians(math.nan))

    def test_distance(self):
        # Two rects separated by exactly 10° in latitude; closest points are
        # hi1.lat=20° and lo2.lat=30°, so distance is 10° in radians.
        lo1 = s2.S2LatLng.from_degrees(10.0, 0.0)
        hi1 = s2.S2LatLng.from_degrees(20.0, 10.0)
        lo2 = s2.S2LatLng.from_degrees(30.0, 0.0)
        hi2 = s2.S2LatLng.from_degrees(40.0, 10.0)
        r1 = s2.S2LatLngRect(lo1, hi1)
        r2 = s2.S2LatLngRect(lo2, hi2)
        d = r1.distance(r2)
        self.assertAlmostEqual(d.degrees, 10.0, places=10)

    def test_distance_rejects_empty(self):
        full = s2.S2LatLngRect.full()
        empty = s2.S2LatLngRect.empty()
        with self.assertRaises(ValueError):
            empty.distance(full)
        with self.assertRaises(ValueError):
            full.distance(empty)

    def test_distance_latlng(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -10.0)
        hi = s2.S2LatLng.from_degrees(10.0, 10.0)
        rect = s2.S2LatLngRect(lo, hi)
        # Interior point: distance is 0.
        self.assertAlmostEqual(
            rect.distance_latlng(s2.S2LatLng.from_degrees(0.0, 0.0)).radians, 0.0)
        # Exterior point due north: closest boundary is at lat=10°, distance=10°.
        d = rect.distance_latlng(s2.S2LatLng.from_degrees(20.0, 0.0))
        self.assertAlmostEqual(d.degrees, 10.0, places=10)

    def test_distance_latlng_rejects_empty(self):
        p = s2.S2LatLng.from_degrees(0.0, 0.0)
        with self.assertRaises(ValueError):
            s2.S2LatLngRect.empty().distance_latlng(p)

    def test_directed_hausdorff_distance(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -10.0)
        hi = s2.S2LatLng.from_degrees(10.0, 10.0)
        r2 = s2.S2LatLngRect(lo, hi)
        # Every point in r2 is in the full rect, so the directed distance is 0.
        d = r2.directed_hausdorff_distance(s2.S2LatLngRect.full())
        self.assertAlmostEqual(d.radians, 0.0)

    def test_hausdorff_distance(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -10.0)
        hi = s2.S2LatLng.from_degrees(10.0, 10.0)
        rect = s2.S2LatLngRect(lo, hi)
        d = rect.hausdorff_distance(rect)
        self.assertAlmostEqual(d.radians, 0.0)

    def test_cap_bound(self):
        bound = s2.S2LatLngRect.full().cap_bound()
        self.assertTrue(bound.is_full())

    def test_rect_bound(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertTrue(rect.rect_bound().approx_equals(rect))

    def test_cell_union_bound(self):
        cell_ids = s2.S2LatLngRect.full().cell_union_bound()
        self.assertEqual(cell_ids, [s2.S2CellId.from_face(i) for i in range(6)])

    # Operators

    def test_equality(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        r1 = s2.S2LatLngRect(lo, hi)
        r2 = s2.S2LatLngRect(lo, hi)
        self.assertTrue(r1 == r2)
        self.assertFalse(r1 != r2)

    def test_inequality(self):
        r1 = s2.S2LatLngRect.empty()
        r2 = s2.S2LatLngRect.full()
        self.assertTrue(r1 != r2)
        self.assertFalse(r1 == r2)

    def test_approx_equals(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        rect = s2.S2LatLngRect(lo, hi)
        self.assertTrue(rect.approx_equals(rect))
        self.assertFalse(rect.approx_equals(s2.S2LatLngRect.full()))

    def test_hash(self):
        lo = s2.S2LatLng.from_degrees(-10.0, -20.0)
        hi = s2.S2LatLng.from_degrees(10.0, 20.0)
        r1 = s2.S2LatLngRect(lo, hi)
        r2 = s2.S2LatLngRect(lo, hi)
        self.assertEqual(hash(r1), hash(r2))

    # String representation

    def test_repr(self):
        rect = s2.S2LatLngRect.empty()
        r = repr(rect)
        self.assertTrue(r.startswith("S2LatLngRect("))
        self.assertIn("Lo", r)
        self.assertIn("Hi", r)

    def test_str(self):
        rect = s2.S2LatLngRect.empty()
        s = str(rect)
        self.assertIn("Lo", s)
        self.assertIn("Hi", s)


if __name__ == "__main__":
    unittest.main()
