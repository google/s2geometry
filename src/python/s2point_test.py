"""Tests for S2Point pybind11 bindings."""

import unittest
import s2geometry_pybind as s2


class TestS2Point(unittest.TestCase):
    """Test cases for S2Point bindings."""

    def test_default_constructor(self):
        p = s2.S2Point()
        self.assertEqual(p.x, 0.0)
        self.assertEqual(p.y, 0.0)
        self.assertEqual(p.z, 0.0)

    def test_constructor_with_coordinates(self):
        p = s2.S2Point(1.0, 2.0, 3.0)
        self.assertEqual(p.x, 1.0)
        self.assertEqual(p.y, 2.0)
        self.assertEqual(p.z, 3.0)

    def test_norm(self):
        p = s2.S2Point(3.0, 4.0, 0.0)
        self.assertAlmostEqual(p.norm(), 5.0)
        self.assertAlmostEqual(p.norm2(), 25.0)

    def test_normalize(self):
        p = s2.S2Point(3.0, 4.0, 0.0)
        normalized = p.normalize()
        self.assertAlmostEqual(normalized.norm(), 1.0)
        self.assertAlmostEqual(normalized.x, 0.6)
        self.assertAlmostEqual(normalized.y, 0.8)
        self.assertAlmostEqual(normalized.z, 0.0)

    def test_dot_product(self):
        p1 = s2.S2Point(1.0, 0.0, 0.0)
        p2 = s2.S2Point(0.0, 1.0, 0.0)
        p3 = s2.S2Point(1.0, 0.0, 0.0)
        self.assertAlmostEqual(p1.dot_prod(p2), 0.0)
        self.assertAlmostEqual(p1.dot_prod(p3), 1.0)

    def test_cross_product(self):
        p1 = s2.S2Point(1.0, 0.0, 0.0)
        p2 = s2.S2Point(0.0, 1.0, 0.0)
        cross = p1.cross_prod(p2)
        self.assertAlmostEqual(cross.x, 0.0)
        self.assertAlmostEqual(cross.y, 0.0)
        self.assertAlmostEqual(cross.z, 1.0)

    def test_addition(self):
        p1 = s2.S2Point(1.0, 2.0, 3.0)
        p2 = s2.S2Point(4.0, 5.0, 6.0)
        result = p1 + p2
        self.assertAlmostEqual(result.x, 5.0)
        self.assertAlmostEqual(result.y, 7.0)
        self.assertAlmostEqual(result.z, 9.0)

    def test_subtraction(self):
        p1 = s2.S2Point(4.0, 5.0, 6.0)
        p2 = s2.S2Point(1.0, 2.0, 3.0)
        result = p1 - p2
        self.assertAlmostEqual(result.x, 3.0)
        self.assertAlmostEqual(result.y, 3.0)
        self.assertAlmostEqual(result.z, 3.0)

    def test_scalar_multiplication(self):
        p = s2.S2Point(1.0, 2.0, 3.0)
        result1 = p * 2.0
        result2 = 2.0 * p
        for result in [result1, result2]:
            self.assertAlmostEqual(result.x, 2.0)
            self.assertAlmostEqual(result.y, 4.0)
            self.assertAlmostEqual(result.z, 6.0)

    def test_scalar_division(self):
        p = s2.S2Point(2.0, 4.0, 6.0)
        result = p / 2.0
        self.assertAlmostEqual(result.x, 1.0)
        self.assertAlmostEqual(result.y, 2.0)
        self.assertAlmostEqual(result.z, 3.0)

    def test_negation(self):
        p = s2.S2Point(1.0, -2.0, 3.0)
        result = -p
        self.assertAlmostEqual(result.x, -1.0)
        self.assertAlmostEqual(result.y, 2.0)
        self.assertAlmostEqual(result.z, -3.0)

    def test_equality(self):
        p1 = s2.S2Point(1.0, 2.0, 3.0)
        p2 = s2.S2Point(1.0, 2.0, 3.0)
        p3 = s2.S2Point(1.0, 2.0, 4.0)

        self.assertEqual(p1, p2)
        self.assertNotEqual(p1, p3)

    def test_repr(self):
        p = s2.S2Point(1.0, 2.0, 3.0)
        self.assertEqual(repr(p), "S2Point(1.000000, 2.000000, 3.000000)")


if __name__ == "__main__":
    unittest.main()
