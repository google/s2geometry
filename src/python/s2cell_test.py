"""Tests for S2Cell pybind11 bindings."""

import math
import unittest
import s2geometry_pybind as s2


class TestS2Cell(unittest.TestCase):
    """Test cases for S2Cell bindings."""

    # Constructors

    def test_constructor_from_cell_id(self):
        cell_id = s2.S2CellId.from_face(0)
        cell = s2.S2Cell(cell_id)
        self.assertEqual(cell.id, cell_id)
        self.assertEqual(cell.face, 0)
        self.assertEqual(cell.level, 0)

    def test_constructor_from_point(self):
        p = s2.S2Point(1.0, 0.0, 0.0)
        cell = s2.S2Cell(p)
        self.assertTrue(cell.is_leaf())

    def test_constructor_from_latlng(self):
        ll = s2.S2LatLng.from_degrees(0.0, 0.0)
        cell = s2.S2Cell(ll)
        self.assertTrue(cell.is_leaf())

    # Factory methods

    def test_from_face(self):
        for face in range(6):
            cell = s2.S2Cell.from_face(face)
            self.assertEqual(cell.face, face)
            self.assertEqual(cell.level, 0)

    def test_from_face_out_of_range_raises(self):
        with self.assertRaises(ValueError) as cm:
            s2.S2Cell.from_face(6)
        self.assertEqual(str(cm.exception), "Face 6 out of range [0, 5]")
        with self.assertRaises(ValueError):
            s2.S2Cell.from_face(-1)

    def test_from_face_pos_level(self):
        cell = s2.S2Cell.from_face_pos_level(0, 0, 0)
        self.assertEqual(cell.level, 0)
        self.assertEqual(cell.face, 0)

    def test_from_face_pos_level_out_of_range_raises(self):
        with self.assertRaises(ValueError):
            s2.S2Cell.from_face_pos_level(6, 0, 0)
        with self.assertRaises(ValueError):
            s2.S2Cell.from_face_pos_level(0, s2.S2CellId.MAX_POSITION + 1, 0)
        with self.assertRaises(ValueError):
            s2.S2Cell.from_face_pos_level(0, 0, 31)

    # Constants

    def test_boundary_constants(self):
        self.assertEqual(int(s2.S2Cell.BOTTOM_EDGE), 0)
        self.assertEqual(int(s2.S2Cell.RIGHT_EDGE), 1)
        self.assertEqual(int(s2.S2Cell.TOP_EDGE), 2)
        self.assertEqual(int(s2.S2Cell.LEFT_EDGE), 3)

    # Properties

    def test_id(self):
        cell_id = s2.S2CellId.from_face(3)
        cell = s2.S2Cell(cell_id)
        self.assertEqual(cell.id, cell_id)

    def test_face(self):
        for f in range(6):
            self.assertEqual(s2.S2Cell.from_face(f).face, f)

    def test_level(self):
        face = s2.S2Cell.from_face(0)
        self.assertEqual(face.level, 0)
        child = s2.S2Cell(face.id.child(0))
        self.assertEqual(child.level, 1)

    def test_orientation(self):
        # Face cells all have orientation 0 (Hilbert curve default).
        for f in range(6):
            self.assertEqual(s2.S2Cell.from_face(f).orientation, 0)
        # Orientation is a bitmask in [0, 3]; non-face cells may differ.
        child = s2.S2Cell(s2.S2CellId.from_face(0).child(0))
        self.assertIn(child.orientation, range(4))

    # Predicates

    def test_is_leaf(self):
        face = s2.S2Cell.from_face(0)
        self.assertFalse(face.is_leaf())
        leaf = s2.S2Cell(s2.S2Point(1.0, 0.0, 0.0))
        self.assertTrue(leaf.is_leaf())

    # Geometric operations

    def test_size_ij(self):
        # A face cell spans 2^kMaxLevel in (i,j).
        face = s2.S2Cell.from_face(0)
        self.assertEqual(face.size_ij(), 1 << s2.S2CellId.MAX_LEVEL)

    def test_size_st(self):
        # A face cell spans the full [0,1] range in (s,t).
        face = s2.S2Cell.from_face(0)
        self.assertAlmostEqual(face.size_st(), 1.0)

    def test_vertex_count_and_normalization(self):
        face = s2.S2Cell.from_face(0)
        for k in range(4):
            v = face.vertex(k)
            # Normalized vertices should have unit norm.
            self.assertAlmostEqual(v.norm(), 1.0)

    def test_vertex_mod_4(self):
        face = s2.S2Cell.from_face(0)
        # Argument is reduced modulo 4.
        v0 = face.vertex(0)
        v4 = face.vertex(4)
        self.assertAlmostEqual(v0.x, v4.x)
        self.assertAlmostEqual(v0.y, v4.y)
        self.assertAlmostEqual(v0.z, v4.z)

    def test_edge(self):
        face = s2.S2Cell.from_face(0)
        for k in range(4):
            e = face.edge(k)
            self.assertAlmostEqual(e.norm(), 1.0)

    def test_uv_coord_of_edge(self):
        face = s2.S2Cell.from_face(0)
        # Face cell UV bounds are [-1, 1] x [-1, 1]. Bottom/top edges are
        # constant in V; left/right constant in U.
        bottom = face.uv_coord_of_edge(s2.S2Cell.Boundary.BOTTOM_EDGE)
        top = face.uv_coord_of_edge(s2.S2Cell.Boundary.TOP_EDGE)
        self.assertAlmostEqual(bottom, -1.0)
        self.assertAlmostEqual(top, 1.0)
        right = face.uv_coord_of_edge(s2.S2Cell.Boundary.RIGHT_EDGE)
        left = face.uv_coord_of_edge(s2.S2Cell.Boundary.LEFT_EDGE)
        self.assertAlmostEqual(right, 1.0)
        self.assertAlmostEqual(left, -1.0)

    def test_ij_coord_of_edge(self):
        face = s2.S2Cell.from_face(0)
        # Face cell in (i,j): bottom/top edges are constant in J (0 and 2^30),
        # left/right edges are constant in I (0 and 2^30).
        self.assertEqual(face.ij_coord_of_edge(s2.S2Cell.Boundary.BOTTOM_EDGE), 0)
        self.assertEqual(face.ij_coord_of_edge(s2.S2Cell.Boundary.TOP_EDGE), 1 << 30)
        self.assertEqual(face.ij_coord_of_edge(s2.S2Cell.Boundary.LEFT_EDGE), 0)
        self.assertEqual(face.ij_coord_of_edge(s2.S2Cell.Boundary.RIGHT_EDGE), 1 << 30)

    def test_center(self):
        face = s2.S2Cell.from_face(0)
        c = face.center()
        # Face 0 center is (1, 0, 0).
        self.assertAlmostEqual(c.norm(), 1.0)
        self.assertAlmostEqual(c.x, 1.0, places=5)
        self.assertAlmostEqual(c.y, 0.0, places=5)
        self.assertAlmostEqual(c.z, 0.0, places=5)

    def test_average_area_for_level(self):
        # AverageArea halves (roughly) with each subdivision by 4: area at
        # level L is ~ (4*pi / 6) / 4^L.
        expected_0 = 4.0 * math.pi / 6.0
        self.assertAlmostEqual(
            s2.S2Cell.average_area_for_level(0), expected_0, places=5)
        # Sum over all cells at level 5 should still be 4*pi.
        total = s2.S2Cell.average_area_for_level(5) * 6 * (4 ** 5)
        self.assertAlmostEqual(total, 4.0 * math.pi, places=5)

    def test_average_area_for_level_out_of_range_raises(self):
        with self.assertRaises(ValueError):
            s2.S2Cell.average_area_for_level(31)
        with self.assertRaises(ValueError):
            s2.S2Cell.average_area_for_level(-1)

    def test_approx_area(self):
        face = s2.S2Cell.from_face(0)
        expected = 4.0 * math.pi / 6.0
        # approx_area is accurate to within 3%.
        self.assertAlmostEqual(face.approx_area(), expected, delta=expected * 0.03)

    def test_exact_area_sums_to_sphere(self):
        # Summing exact_area over all 6 face cells should give 4*pi.
        total = sum(s2.S2Cell.from_face(f).exact_area() for f in range(6))
        self.assertAlmostEqual(total, 4.0 * math.pi, places=10)

    def test_bound_uv_face(self):
        face = s2.S2Cell.from_face(0)
        uv = face.bound_uv()
        self.assertAlmostEqual(uv.lo.x, -1.0)
        self.assertAlmostEqual(uv.lo.y, -1.0)
        self.assertAlmostEqual(uv.hi.x, 1.0)
        self.assertAlmostEqual(uv.hi.y, 1.0)

    def test_distance_to_point_outside(self):
        face = s2.S2Cell.from_face(0)
        # A point on the opposite face should be at a positive distance.
        far = s2.S2Point(-1.0, 0.0, 0.0)
        d = face.distance(far)
        self.assertGreater(d.radians, 0.0)

    def test_distance_to_point_inside(self):
        face = s2.S2Cell.from_face(0)
        # Center of face 0 is inside; distance should be zero.
        center = face.center()
        self.assertAlmostEqual(face.distance(center).radians, 0.0)

    def test_boundary_distance(self):
        face = s2.S2Cell.from_face(0)
        center = face.center()
        # Interior point has positive boundary distance.
        self.assertGreater(face.boundary_distance(center).radians, 0.0)
        # A far exterior point also has positive boundary distance.
        far = s2.S2Point(-1.0, 0.0, 0.0)
        self.assertGreater(face.boundary_distance(far).radians, 0.0)

    def test_max_distance_to_point(self):
        face = s2.S2Cell.from_face(0)
        center = face.center()
        # Max distance from a face cell to its own center should be positive
        # (distance to the farthest corner).
        self.assertGreater(face.max_distance(center).radians, 0.0)

    def test_distance_to_edge(self):
        face = s2.S2Cell.from_face(0)
        # An edge entirely on the opposite face should be at positive distance.
        a = s2.S2Point(-1.0, 0.1, 0.0)
        b = s2.S2Point(-1.0, -0.1, 0.0)
        self.assertGreater(face.distance_to_edge(a, b).radians, 0.0)

    def test_max_distance_to_edge(self):
        face = s2.S2Cell.from_face(0)
        a = s2.S2Point(-1.0, 0.1, 0.0)
        b = s2.S2Point(-1.0, -0.1, 0.0)
        self.assertGreater(face.max_distance_to_edge(a, b).radians, 0.0)

    def test_distance_to_cell(self):
        face0 = s2.S2Cell.from_face(0)
        face1 = s2.S2Cell.from_face(1)
        # Distance from a cell to itself is zero.
        self.assertAlmostEqual(face0.distance_to_cell(face0).radians, 0.0)
        # Distance between disjoint faces is positive.
        self.assertGreater(face0.distance_to_cell(face1).radians, 0.0)

    def test_max_distance_to_cell(self):
        face0 = s2.S2Cell.from_face(0)
        face1 = s2.S2Cell.from_face(1)
        self.assertGreater(face0.max_distance_to_cell(face1).radians, 0.0)

    def test_cell_union_bound(self):
        cell = s2.S2Cell.from_face(0)
        bound = cell.cell_union_bound()
        self.assertEqual(len(bound), 1)
        self.assertEqual(bound[0], cell.id)

    def test_contains_cell(self):
        face = s2.S2Cell.from_face(0)
        child = s2.S2Cell(face.id.child(0))
        self.assertTrue(face.contains(child))
        self.assertFalse(child.contains(face))

    def test_contains_point(self):
        face = s2.S2Cell.from_face(0)
        # Face 0 center is (1, 0, 0), which is inside face 0.
        self.assertTrue(face.contains_point(s2.S2Point(1.0, 0.0, 0.0)))
        # A point on the opposite face is not in face 0.
        self.assertFalse(face.contains_point(s2.S2Point(-1.0, 0.0, 0.0)))

    def test_may_intersect(self):
        face0 = s2.S2Cell.from_face(0)
        face1 = s2.S2Cell.from_face(1)
        # A cell intersects itself.
        self.assertTrue(face0.may_intersect(face0))
        # Distinct faces (disjoint cell id ranges) do not intersect.
        self.assertFalse(face0.may_intersect(face1))
        # A cell contained within another may_intersect the outer cell.
        child = s2.S2Cell(face0.id.child(0))
        self.assertTrue(face0.may_intersect(child))
        self.assertTrue(child.may_intersect(face0))

    # Traversal

    def test_subdivide(self):
        face = s2.S2Cell.from_face(0)
        children = face.subdivide()
        self.assertEqual(len(children), 4)
        for i, child in enumerate(children):
            self.assertEqual(child.level, 1)
            self.assertEqual(child.face, 0)
            # Subdivide should produce the same children as constructing
            # from the S2CellId children.
            self.assertEqual(child.id, face.id.child(i))

    def test_subdivide_leaf_raises(self):
        leaf = s2.S2Cell(s2.S2Point(1.0, 0.0, 0.0))
        with self.assertRaises(ValueError) as cm:
            leaf.subdivide()
        self.assertEqual(str(cm.exception), "Leaf cell has no children")

    # Operators

    def test_equality(self):
        a = s2.S2Cell.from_face(0)
        b = s2.S2Cell.from_face(0)
        c = s2.S2Cell.from_face(1)
        self.assertTrue(a == b)
        self.assertTrue(a != c)

    def test_comparison(self):
        a = s2.S2Cell.from_face(0)
        b = s2.S2Cell.from_face(1)
        self.assertTrue(a < b)
        self.assertTrue(b > a)
        self.assertTrue(a <= a)
        self.assertTrue(a >= a)

    def test_hash(self):
        a = s2.S2Cell.from_face(0)
        b = s2.S2Cell.from_face(0)
        self.assertEqual(hash(a), hash(b))
        s = {a, b}
        self.assertEqual(len(s), 1)

    # String representation

    def test_repr(self):
        cell = s2.S2Cell.from_face(0)
        self.assertEqual(repr(cell), "S2Cell(0/)")

    def test_str(self):
        cell = s2.S2Cell.from_face(0)
        self.assertEqual(str(cell), "0/")

    def test_repr_child(self):
        cell = s2.S2Cell(s2.S2CellId.from_face(3).child(0).child(2))
        self.assertEqual(repr(cell), "S2Cell(3/02)")


if __name__ == "__main__":
    unittest.main()
