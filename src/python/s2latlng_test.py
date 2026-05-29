"""Tests for S2LatLng pybind11 bindings."""

import math
import unittest
import s2geometry_pybind as s2


class TestS2LatLng(unittest.TestCase):
    """Test cases for S2LatLng bindings."""

    # Constructors

    def test_default_constructor(self):
        ll = s2.S2LatLng()
        self.assertEqual(ll.lat.radians, 0.0)
        self.assertEqual(ll.lng.radians, 0.0)

    def test_constructor_from_angles(self):
        lat = s2.S1Angle.from_degrees(37.0)
        lng = s2.S1Angle.from_degrees(-122.0)
        ll = s2.S2LatLng(lat, lng)
        self.assertAlmostEqual(ll.lat.degrees, 37.0)
        self.assertAlmostEqual(ll.lng.degrees, -122.0)

    def test_constructor_from_angles_at_bounds(self):
        ll = s2.S2LatLng(
            s2.S1Angle.from_degrees(90.0), s2.S1Angle.from_degrees(180.0)
        )
        self.assertAlmostEqual(ll.lat.degrees, 90.0)
        self.assertAlmostEqual(ll.lng.degrees, 180.0)

    def test_constructor_from_angles_invalid_lat_raises(self):
        with self.assertRaises(ValueError) as cm:
            s2.S2LatLng(
                s2.S1Angle.from_degrees(91.0), s2.S1Angle.from_degrees(0.0)
            )
        self.assertEqual(
            str(cm.exception),
            "Invalid S2LatLng: (91, 0) "
            "(latitude must be in [-90, 90], longitude in [-180, 180])",
        )

    def test_constructor_from_angles_invalid_lng_raises(self):
        with self.assertRaises(ValueError):
            s2.S2LatLng(
                s2.S1Angle.from_degrees(0.0), s2.S1Angle.from_degrees(181.0)
            )

    def test_constructor_from_point(self):
        p = s2.S2Point(1.0, 0.0, 0.0)
        ll = s2.S2LatLng(p)
        self.assertAlmostEqual(ll.lat.degrees, 0.0)
        self.assertAlmostEqual(ll.lng.degrees, 0.0)

    def test_constructor_from_nan_point_raises(self):
        p = s2.S2Point(float("nan"), 0.0, 0.0)
        with self.assertRaises(ValueError):
            s2.S2LatLng(p)

    def test_constructor_from_point_north_pole(self):
        p = s2.S2Point(0.0, 0.0, 1.0)
        ll = s2.S2LatLng(p)
        self.assertAlmostEqual(ll.lat.degrees, 90.0)

    # Factory methods

    def test_from_radians(self):
        ll = s2.S2LatLng.from_radians(math.pi / 4, -math.pi / 2)
        self.assertAlmostEqual(ll.lat.degrees, 45.0)
        self.assertAlmostEqual(ll.lng.degrees, -90.0)

    def test_from_radians_invalid_raises(self):
        with self.assertRaises(ValueError):
            s2.S2LatLng.from_radians(2.0, 0.0)  # lat > pi/2

    def test_normalized_from_radians_clamps_latitude(self):
        ll = s2.S2LatLng.normalized_from_radians(2.0, 0.0)  # lat > pi/2
        self.assertAlmostEqual(ll.lat.radians, math.pi / 2)
        self.assertAlmostEqual(ll.lng.radians, 0.0)

    def test_normalized_from_radians_wraps_longitude(self):
        # 3*pi/2 wraps to -pi/2
        ll = s2.S2LatLng.normalized_from_radians(0.0, 3 * math.pi / 2)
        self.assertAlmostEqual(ll.lat.radians, 0.0)
        self.assertAlmostEqual(ll.lng.radians, -math.pi / 2)

    def test_normalized_from_radians_in_range_is_unchanged(self):
        ll = s2.S2LatLng.normalized_from_radians(0.5, -1.0)
        self.assertAlmostEqual(ll.lat.radians, 0.5)
        self.assertAlmostEqual(ll.lng.radians, -1.0)

    def test_normalized_from_radians_non_finite_raises(self):
        with self.assertRaises(ValueError):
            s2.S2LatLng.normalized_from_radians(float("nan"), 0.0)
        with self.assertRaises(ValueError):
            s2.S2LatLng.normalized_from_radians(0.0, float("inf"))

    def test_from_degrees(self):
        ll = s2.S2LatLng.from_degrees(45.0, -90.0)
        self.assertAlmostEqual(ll.lat.degrees, 45.0)
        self.assertAlmostEqual(ll.lng.degrees, -90.0)

    def test_from_degrees_invalid_raises(self):
        with self.assertRaises(ValueError):
            s2.S2LatLng.from_degrees(91.0, 0.0)
        with self.assertRaises(ValueError):
            s2.S2LatLng.from_degrees(0.0, 181.0)

    def test_from_e5(self):
        ll = s2.S2LatLng.from_e5(4500000, -9000000)
        self.assertAlmostEqual(ll.lat.degrees, 45.0)
        self.assertAlmostEqual(ll.lng.degrees, -90.0)

    def test_from_e5_invalid_raises(self):
        with self.assertRaises(ValueError):
            s2.S2LatLng.from_e5(9100000, 0)  # 91 degrees

    def test_from_e6(self):
        ll = s2.S2LatLng.from_e6(45000000, -90000000)
        self.assertAlmostEqual(ll.lat.degrees, 45.0)
        self.assertAlmostEqual(ll.lng.degrees, -90.0)

    def test_from_e6_invalid_raises(self):
        with self.assertRaises(ValueError):
            s2.S2LatLng.from_e6(91000000, 0)  # 91 degrees

    def test_from_e7(self):
        ll = s2.S2LatLng.from_e7(450000000, -900000000)
        self.assertAlmostEqual(ll.lat.degrees, 45.0)
        self.assertAlmostEqual(ll.lng.degrees, -90.0)

    def test_from_e7_invalid_raises(self):
        with self.assertRaises(ValueError):
            s2.S2LatLng.from_e7(910000000, 0)  # 91 degrees

    def test_latitude_static(self):
        p = s2.S2Point(0.0, 0.0, 1.0)
        lat = s2.S2LatLng.latitude(p)
        self.assertAlmostEqual(lat.degrees, 90.0)

    def test_longitude_static(self):
        p = s2.S2Point(0.0, 1.0, 0.0)
        lng = s2.S2LatLng.longitude(p)
        self.assertAlmostEqual(lng.degrees, 90.0)

    # Properties

    def test_lat_lng_are_s1angles(self):
        ll = s2.S2LatLng.from_degrees(37.0, -122.0)
        self.assertIsInstance(ll.lat, s2.S1Angle)
        self.assertIsInstance(ll.lng, s2.S1Angle)

    def test_lat_lng_readonly(self):
        ll = s2.S2LatLng.from_degrees(37.0, -122.0)
        with self.assertRaises(AttributeError):
            ll.lat = s2.S1Angle.from_degrees(0.0)
        with self.assertRaises(AttributeError):
            ll.lng = s2.S1Angle.from_degrees(0.0)

    def test_coords(self):
        ll = s2.S2LatLng.from_radians(0.5, -1.0)
        coords = ll.coords
        self.assertIsInstance(coords, s2.R2Point)
        self.assertAlmostEqual(coords.x, 0.5)
        self.assertAlmostEqual(coords.y, -1.0)

    # Geometric operations

    def test_to_point(self):
        ll = s2.S2LatLng.from_degrees(0.0, 0.0)
        p = ll.to_point()
        self.assertAlmostEqual(p.x, 1.0)
        self.assertAlmostEqual(p.y, 0.0)
        self.assertAlmostEqual(p.z, 0.0)

    def test_to_point_north_pole(self):
        ll = s2.S2LatLng.from_degrees(90.0, 0.0)
        p = ll.to_point()
        self.assertAlmostEqual(p.x, 0.0, places=15)
        self.assertAlmostEqual(p.y, 0.0, places=15)
        self.assertAlmostEqual(p.z, 1.0)

    def test_to_point_roundtrip(self):
        ll = s2.S2LatLng.from_degrees(37.7749, -122.4194)
        p = ll.to_point()
        ll2 = s2.S2LatLng(p)
        self.assertAlmostEqual(ll.lat.radians, ll2.lat.radians)
        self.assertAlmostEqual(ll.lng.radians, ll2.lng.radians)

    def test_get_distance(self):
        # Distance between two points on the equator, 90 degrees apart.
        a = s2.S2LatLng.from_degrees(0.0, 0.0)
        b = s2.S2LatLng.from_degrees(0.0, 90.0)
        self.assertAlmostEqual(a.get_distance(b).degrees, 90.0)

    def test_get_distance_same_point(self):
        ll = s2.S2LatLng.from_degrees(37.0, -122.0)
        self.assertAlmostEqual(ll.get_distance(ll).radians, 0.0)

    def test_to_string_in_degrees(self):
        ll = s2.S2LatLng.from_degrees(37.794, -122.395)
        self.assertEqual(ll.to_string_in_degrees(), "37.794000,-122.395000")

    def test_approx_equals(self):
        a = s2.S2LatLng.from_degrees(37.0, -122.0)
        b = s2.S2LatLng.from_degrees(37.0, -122.0)
        self.assertTrue(a.approx_equals(b))

    def test_approx_equals_within_tolerance(self):
        a = s2.S2LatLng.from_degrees(37.0, -122.0)
        b = s2.S2LatLng.from_degrees(37.0 + 1e-16, -122.0)
        self.assertTrue(a.approx_equals(b))

    def test_approx_equals_outside_default_tolerance(self):
        a = s2.S2LatLng.from_degrees(37.0, -122.0)
        b = s2.S2LatLng.from_degrees(37.1, -122.0)
        self.assertFalse(a.approx_equals(b))

    def test_approx_equals_custom_tolerance(self):
        a = s2.S2LatLng.from_degrees(37.0, -122.0)
        b = s2.S2LatLng.from_degrees(37.1, -122.0)
        self.assertTrue(a.approx_equals(b, s2.S1Angle.from_degrees(0.2)))
        self.assertFalse(a.approx_equals(b, s2.S1Angle.from_degrees(0.05)))

    # Operators

    def test_equality(self):
        a = s2.S2LatLng.from_degrees(37.0, -122.0)
        b = s2.S2LatLng.from_degrees(37.0, -122.0)
        c = s2.S2LatLng.from_degrees(38.0, -122.0)
        self.assertTrue(a == b)
        self.assertTrue(a != c)
        self.assertFalse(a != b)
        self.assertFalse(a == c)

    def test_comparison(self):
        a = s2.S2LatLng.from_degrees(10.0, 20.0)
        b = s2.S2LatLng.from_degrees(20.0, 10.0)
        self.assertTrue(a < b)
        self.assertTrue(b > a)
        self.assertTrue(a <= b)
        self.assertTrue(b >= a)
        self.assertTrue(a <= a)
        self.assertTrue(a >= a)

    def test_addition(self):
        a = s2.S2LatLng.from_degrees(10.0, 20.0)
        b = s2.S2LatLng.from_degrees(30.0, 40.0)
        result = a + b
        self.assertAlmostEqual(result.lat.degrees, 40.0)
        self.assertAlmostEqual(result.lng.degrees, 60.0)

    def test_addition_normalizes(self):
        # 80 + 20 = 100 degrees latitude, should clamp to 90.
        a = s2.S2LatLng.from_degrees(80.0, 170.0)
        b = s2.S2LatLng.from_degrees(20.0, 30.0)
        result = a + b
        self.assertAlmostEqual(result.lat.degrees, 90.0)
        # 170 + 30 = 200 degrees longitude, should wrap to -160.
        self.assertAlmostEqual(result.lng.degrees, -160.0)

    def test_subtraction(self):
        a = s2.S2LatLng.from_degrees(30.0, 40.0)
        b = s2.S2LatLng.from_degrees(10.0, 20.0)
        result = a - b
        self.assertAlmostEqual(result.lat.degrees, 20.0)
        self.assertAlmostEqual(result.lng.degrees, 20.0)

    def test_subtraction_normalizes(self):
        # -80 - 20 = -100 degrees latitude, should clamp to -90.
        a = s2.S2LatLng.from_degrees(-80.0, -170.0)
        b = s2.S2LatLng.from_degrees(20.0, 30.0)
        result = a - b
        self.assertAlmostEqual(result.lat.degrees, -90.0)
        # -170 - 30 = -200 degrees longitude, should wrap to 160.
        self.assertAlmostEqual(result.lng.degrees, 160.0)

    def test_scalar_multiplication(self):
        ll = s2.S2LatLng.from_degrees(10.0, 20.0)
        result1 = ll * 2.0
        result2 = 2.0 * ll
        for result in [result1, result2]:
            self.assertAlmostEqual(result.lat.degrees, 20.0)
            self.assertAlmostEqual(result.lng.degrees, 40.0)

    def test_scalar_multiplication_inf_raises(self):
        ll = s2.S2LatLng.from_degrees(10.0, 20.0)
        with self.assertRaises(ValueError):
            ll * float("inf")

    def test_scalar_multiplication_nan_raises(self):
        ll = s2.S2LatLng.from_degrees(10.0, 20.0)
        with self.assertRaises(ValueError):
            ll * float("nan")

    def test_rmul_inf_raises(self):
        ll = s2.S2LatLng.from_degrees(10.0, 20.0)
        with self.assertRaises(ValueError):
            float("inf") * ll

    def test_scalar_multiplication_normalizes(self):
        ll = s2.S2LatLng.from_degrees(60.0, 120.0)
        result = ll * 2.0
        # 120 degrees latitude -> clamped to 90
        self.assertAlmostEqual(result.lat.degrees, 90.0)
        # 240 degrees longitude -> wrapped to -120
        self.assertAlmostEqual(result.lng.degrees, -120.0)

    # String representation

    def test_repr(self):
        ll = s2.S2LatLng.from_degrees(45.0, -90.0)
        self.assertEqual(repr(ll), "S2LatLng([45.0000000, -90.0000000])")

    def test_str(self):
        ll = s2.S2LatLng.from_degrees(45.0, -90.0)
        self.assertEqual(str(ll), "[45.0000000, -90.0000000]")

    def test_repr_zero(self):
        ll = s2.S2LatLng()
        self.assertEqual(repr(ll), "S2LatLng([0.0000000, 0.0000000])")


if __name__ == "__main__":
    unittest.main()
