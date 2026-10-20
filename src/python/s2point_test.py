"""Tests for S2Point pybind11 bindings."""

import math
import unittest
import s2geometry_pybind as s2


class TestS2Point(unittest.TestCase):
    """Test cases for S2Point bindings."""

    # Constructors

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

    def test_unnormalized_point(self):
        # Unnormalized points are accepted by the constructor.
        p = s2.S2Point(3.0, 4.0, 0.0)
        self.assertAlmostEqual(p.x, 3.0)
        self.assertAlmostEqual(p.y, 4.0)
        self.assertAlmostEqual(p.z, 0.0)
        self.assertAlmostEqual(p.norm(), 5.0)  # Not on unit sphere

        # User may call normalize() to obtain a unit vector.
        normalized = p.normalize()
        self.assertAlmostEqual(normalized.norm(), 1.0)

    def test_constructor_from_tuple(self):
        p = s2.S2Point((1.0, 2.0, 3.0))
        self.assertEqual(p.x, 1.0)
        self.assertEqual(p.y, 2.0)
        self.assertEqual(p.z, 3.0)

    def test_constructor_from_tuple_wrong_length(self):
        with self.assertRaises(ValueError):
            s2.S2Point((1.0, 2.0))
        with self.assertRaises(ValueError):
            s2.S2Point((1.0, 2.0, 3.0, 4.0))

    def test_constructor_from_sequence(self):
        # List of floats
        p1 = s2.S2Point([1.0, 2.0, 3.0])
        self.assertEqual(p1.x, 1.0)
        self.assertEqual(p1.y, 2.0)
        self.assertEqual(p1.z, 3.0)

        # List of ints
        p2 = s2.S2Point([1, 2, 3])
        self.assertEqual(p2.x, 1.0)
        self.assertEqual(p2.y, 2.0)
        self.assertEqual(p2.z, 3.0)

    def test_constructor_from_sequence_wrong_length(self):
        with self.assertRaises(ValueError):
            s2.S2Point([])
        with self.assertRaises(ValueError):
            s2.S2Point([1.0, 2.0])
        with self.assertRaises(ValueError):
            s2.S2Point([1.0, 2.0, 3.0, 4.0])
        with self.assertRaises(ValueError):
            s2.S2Point(["a", "b"])

    def test_constructor_from_sequence_wrong_types(self):
        # Strings and bytes of any length consistently raise TypeError.
        with self.assertRaises(TypeError):
            s2.S2Point("")
        with self.assertRaises(TypeError):
            s2.S2Point("ab")
        with self.assertRaises(TypeError):
            s2.S2Point("abc")
        with self.assertRaises(TypeError):
            s2.S2Point(b"12")
        with self.assertRaises(TypeError):
            s2.S2Point(b"123")
        with self.assertRaises(TypeError):
            s2.S2Point([1.0, 2.0, "three"])

    # Sequence protocol

    def test_len(self):
        p = s2.S2Point(1.0, 2.0, 3.0)
        self.assertEqual(len(p), 3)

    def test_getitem(self):
        p = s2.S2Point(1.0, 2.0, 3.0)
        self.assertEqual(p[0], 1.0)
        self.assertEqual(p[1], 2.0)
        self.assertEqual(p[2], 3.0)
        self.assertEqual(p[-1], 3.0)
        self.assertEqual(p[-2], 2.0)
        self.assertEqual(p[-3], 1.0)

        with self.assertRaises(IndexError):
            _ = p[3]
        with self.assertRaises(IndexError):
            _ = p[-4]
        # Large indices raise IndexError rather than overflow TypeError.
        with self.assertRaises(IndexError):
            _ = p[10**18]
        with self.assertRaises(IndexError):
            _ = p[-10**18]
        # Slicing is not supported on S2Point and raises TypeError.
        with self.assertRaises(TypeError):
            _ = p[0:2]

    def test_unpacking(self):
        p = s2.S2Point(1.0, 2.0, 3.0)
        x, y, z = p
        self.assertEqual((x, y, z), (1.0, 2.0, 3.0))

        with self.assertRaises(ValueError):
            a, b = p
        with self.assertRaises(ValueError):
            a, b, c, d = p

    def test_iter(self):
        p = s2.S2Point(1.0, 2.0, 3.0)
        self.assertEqual(list(p), [1.0, 2.0, 3.0])
        self.assertEqual(tuple(p), (1.0, 2.0, 3.0))
        self.assertEqual([c for c in p], [1.0, 2.0, 3.0])

    # Properties

    def test_data(self):
        p = s2.S2Point(1.0, 2.0, 3.0)
        coords = p.data()
        self.assertEqual(coords, (1.0, 2.0, 3.0))

    # Vector operations

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

    def test_angle(self):
        p1 = s2.S2Point(1.0, 0.0, 0.0)
        p2 = s2.S2Point(0.0, 1.0, 0.0)
        self.assertAlmostEqual(p1.angle(p2), math.pi / 2)

    # Operators

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

    def test_hash(self):
        a = s2.S2Point(1.0, 0.0, 0.0)
        b = s2.S2Point(1.0, 0.0, 0.0)
        self.assertEqual(hash(a), hash(b))
        self.assertEqual(len({a, b}), 1)

    # String representation

    def test_string_representation(self):
        p = s2.S2Point(1.0, 2.0, 3.0)
        self.assertEqual(repr(p), "S2Point([1, 2, 3])")
        self.assertEqual(str(p), "[1, 2, 3]")


if __name__ == "__main__":
    unittest.main()
