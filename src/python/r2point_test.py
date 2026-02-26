"""Tests for R2Point pybind11 bindings."""

import math
import unittest
import s2geometry_pybind as s2


class TestR2Point(unittest.TestCase):
    """Test cases for R2Point bindings."""

    # Constructors

    def test_default_constructor(self):
        p = s2.R2Point()
        self.assertEqual(p.x, 0.0)
        self.assertEqual(p.y, 0.0)

    def test_constructor_with_coordinates(self):
        p = s2.R2Point(1.0, 2.0)
        self.assertEqual(p.x, 1.0)
        self.assertEqual(p.y, 2.0)

    def test_constructor_from_tuple(self):
        p = s2.R2Point((3.0, 4.0))
        self.assertEqual(p.x, 3.0)
        self.assertEqual(p.y, 4.0)

    def test_constructor_from_tuple_wrong_length(self):
        with self.assertRaises(ValueError):
            s2.R2Point((1.0,))
        with self.assertRaises(ValueError):
            s2.R2Point((1.0, 2.0, 3.0))

    # Properties

    def test_properties_are_readonly(self):
        p = s2.R2Point(1.0, 2.0)
        with self.assertRaises(AttributeError):
            p.x = 5.0
        with self.assertRaises(AttributeError):
            p.y = 5.0

    def test_data(self):
        p = s2.R2Point(1.0, 2.0)
        coords = p.data()
        self.assertEqual(coords, (1.0, 2.0))

    # Vector operations

    def test_norm(self):
        p = s2.R2Point(3.0, 4.0)
        self.assertAlmostEqual(p.norm(), 5.0)
        self.assertAlmostEqual(p.norm2(), 25.0)

    def test_normalize(self):
        p = s2.R2Point(3.0, 4.0)
        normalized = p.normalize()
        self.assertAlmostEqual(normalized.norm(), 1.0)
        self.assertAlmostEqual(normalized.x, 0.6)
        self.assertAlmostEqual(normalized.y, 0.8)

    def test_dot_product(self):
        p1 = s2.R2Point(1.0, 0.0)
        p2 = s2.R2Point(0.0, 1.0)
        p3 = s2.R2Point(2.0, 3.0)
        self.assertAlmostEqual(p1.dot_prod(p2), 0.0)
        self.assertAlmostEqual(p1.dot_prod(p3), 2.0)

    def test_cross_product(self):
        p1 = s2.R2Point(1.0, 0.0)
        p2 = s2.R2Point(0.0, 1.0)
        # 2D cross product returns scalar: x1*y2 - y1*x2
        self.assertAlmostEqual(p1.cross_prod(p2), 1.0)

    def test_angle(self):
        p1 = s2.R2Point(1.0, 0.0)
        p2 = s2.R2Point(0.0, 1.0)
        self.assertAlmostEqual(p1.angle(p2), math.pi / 2)

    def test_ortho(self):
        p = s2.R2Point(1.0, 0.0)
        o = p.ortho()
        # Rotated 90 degrees CCW: (1, 0) -> (0, 1)
        self.assertAlmostEqual(o.x, 0.0)
        self.assertAlmostEqual(o.y, 1.0)

    def test_ortho_general(self):
        p = s2.R2Point(3.0, 4.0)
        o = p.ortho()
        # Should be perpendicular
        self.assertAlmostEqual(p.dot_prod(o), 0.0)
        # Same norm
        self.assertAlmostEqual(o.norm(), p.norm())

    def test_fabs(self):
        p = s2.R2Point(-3.0, -4.0)
        a = p.fabs()
        self.assertAlmostEqual(a.x, 3.0)
        self.assertAlmostEqual(a.y, 4.0)

    def test_fabs_mixed(self):
        p = s2.R2Point(-1.0, 2.0)
        a = p.fabs()
        self.assertAlmostEqual(a.x, 1.0)
        self.assertAlmostEqual(a.y, 2.0)

    # Operators

    def test_addition(self):
        p1 = s2.R2Point(1.0, 2.0)
        p2 = s2.R2Point(3.0, 4.0)
        result = p1 + p2
        self.assertAlmostEqual(result.x, 4.0)
        self.assertAlmostEqual(result.y, 6.0)

    def test_subtraction(self):
        p1 = s2.R2Point(3.0, 4.0)
        p2 = s2.R2Point(1.0, 2.0)
        result = p1 - p2
        self.assertAlmostEqual(result.x, 2.0)
        self.assertAlmostEqual(result.y, 2.0)

    def test_scalar_multiplication(self):
        p = s2.R2Point(1.0, 2.0)
        result1 = p * 3.0
        result2 = 3.0 * p
        for result in [result1, result2]:
            self.assertAlmostEqual(result.x, 3.0)
            self.assertAlmostEqual(result.y, 6.0)

    def test_scalar_division(self):
        p = s2.R2Point(4.0, 6.0)
        result = p / 2.0
        self.assertAlmostEqual(result.x, 2.0)
        self.assertAlmostEqual(result.y, 3.0)

    def test_negation(self):
        p = s2.R2Point(1.0, -2.0)
        result = -p
        self.assertAlmostEqual(result.x, -1.0)
        self.assertAlmostEqual(result.y, 2.0)

    def test_equality(self):
        p1 = s2.R2Point(1.0, 2.0)
        p2 = s2.R2Point(1.0, 2.0)
        p3 = s2.R2Point(1.0, 3.0)
        self.assertEqual(p1, p2)
        self.assertNotEqual(p1, p3)

    def test_in_place_addition(self):
        p = s2.R2Point(1.0, 2.0)
        p += s2.R2Point(3.0, 4.0)
        self.assertAlmostEqual(p.x, 4.0)
        self.assertAlmostEqual(p.y, 6.0)

    def test_in_place_subtraction(self):
        p = s2.R2Point(3.0, 4.0)
        p -= s2.R2Point(1.0, 2.0)
        self.assertAlmostEqual(p.x, 2.0)
        self.assertAlmostEqual(p.y, 2.0)

    def test_in_place_multiplication(self):
        p = s2.R2Point(1.0, 2.0)
        p *= 3.0
        self.assertAlmostEqual(p.x, 3.0)
        self.assertAlmostEqual(p.y, 6.0)

    def test_in_place_division(self):
        p = s2.R2Point(4.0, 6.0)
        p /= 2.0
        self.assertAlmostEqual(p.x, 2.0)
        self.assertAlmostEqual(p.y, 3.0)

    # String representation

    def test_string_representation(self):
        p = s2.R2Point(1.0, 2.0)
        self.assertEqual(repr(p), "R2Point([1, 2])")
        self.assertEqual(str(p), "[1, 2]")


if __name__ == "__main__":
    unittest.main()
