"""Tests for S2CellId pybind11 bindings."""

import unittest
import s2geometry_pybind as s2


class TestS2CellId(unittest.TestCase):
    """Test cases for S2CellId bindings."""

    # Constructors

    def test_constructor_from_point(self):
        p = s2.S2Point(1.0, 0.0, 0.0)
        cell = s2.S2CellId(p)
        self.assertTrue(cell.is_leaf())

    def test_constructor_from_latlng(self):
        ll = s2.S2LatLng.from_degrees(0.0, 0.0)
        cell = s2.S2CellId(ll)
        self.assertTrue(cell.is_leaf())

    def test_constructor_from_uint64(self):
        # Construct from a face cell, then use its id to reconstruct.
        face_cell = s2.S2CellId.from_face(0)
        cell = s2.S2CellId(face_cell.id)
        self.assertEqual(cell, face_cell)

    def test_constructor_from_invalid_uint64_raises(self):
        with self.assertRaises(ValueError):
            s2.S2CellId(0)

    # Constants

    def test_max_level(self):
        self.assertEqual(s2.S2CellId.kMaxLevel, 30)

    def test_num_faces(self):
        self.assertEqual(s2.S2CellId.kNumFaces, 6)

    # Factory methods

    def test_from_face(self):
        for face in range(6):
            cell = s2.S2CellId.from_face(face)
            self.assertEqual(cell.face(), face)
            self.assertTrue(cell.is_face())

    def test_from_face_out_of_range_raises(self):
        with self.assertRaises(ValueError) as cm:
            s2.S2CellId.from_face(6)
        self.assertEqual(str(cm.exception), "Face 6 out of range [0, 5]")
        with self.assertRaises(ValueError):
            s2.S2CellId.from_face(-1)

    def test_from_face_pos_level(self):
        cell = s2.S2CellId.from_face_pos_level(0, 0, 0)
        self.assertEqual(cell.level(), 0)
        self.assertEqual(cell.face(), 0)

    def test_from_token_roundtrip(self):
        cell = s2.S2CellId.from_face(3)
        token = cell.to_token()
        self.assertEqual(token, "7")
        cell2 = s2.S2CellId.from_token(token)
        self.assertEqual(cell, cell2)

    def test_from_token_invalid_raises(self):
        with self.assertRaises(ValueError) as cm:
            s2.S2CellId.from_token("invalid_token")
        self.assertEqual(str(cm.exception),
                         "Invalid S2CellId token: 'invalid_token'")

    def test_from_debug_string(self):
        cell = s2.S2CellId.from_debug_string("3/02")
        self.assertEqual(cell.face(), 3)
        self.assertEqual(cell.level(), 2)

    def test_from_debug_string_invalid_raises(self):
        with self.assertRaises(ValueError) as cm:
            s2.S2CellId.from_debug_string("bad")
        self.assertEqual(str(cm.exception),
                         "Invalid S2CellId debug string: 'bad'")

    def test_from_face_ij(self):
        cell = s2.S2CellId.from_face_ij(0, 0, 0)
        self.assertTrue(cell.is_leaf())

    def test_cells_level_out_of_range_raises(self):
        with self.assertRaises(ValueError) as cm:
            s2.S2CellId.cells(31)
        self.assertEqual(str(cm.exception), "Level 31 out of range [0, 30]")

    # Properties

    def test_id_property(self):
        cell = s2.S2CellId.from_face(0)
        self.assertIsInstance(cell.id, int)
        self.assertGreater(cell.id, 0)

    # Predicates

    def test_is_leaf(self):
        p = s2.S2Point(1.0, 0.0, 0.0)
        leaf = s2.S2CellId(p)
        self.assertTrue(leaf.is_leaf())
        face = s2.S2CellId.from_face(0)
        self.assertFalse(face.is_leaf())

    def test_is_face(self):
        face = s2.S2CellId.from_face(0)
        self.assertTrue(face.is_face())
        child = face.child(0)
        self.assertFalse(child.is_face())

    # Geometric operations

    def test_face(self):
        for f in range(6):
            self.assertEqual(s2.S2CellId.from_face(f).face(), f)

    def test_level(self):
        face = s2.S2CellId.from_face(0)
        self.assertEqual(face.level(), 0)
        child = face.child(0)
        self.assertEqual(child.level(), 1)

    def test_to_point(self):
        cell = s2.S2CellId.from_face(0)
        p = cell.to_point()
        # Face 0 center is approximately (1, 0, 0).
        self.assertAlmostEqual(p.x, 1.0, places=1)
        self.assertAlmostEqual(p.y, 0.0, places=1)
        self.assertAlmostEqual(p.z, 0.0, places=1)

    def test_to_lat_lng(self):
        ll = s2.S2LatLng.from_degrees(45.0, 90.0)
        cell = s2.S2CellId(ll)
        ll2 = cell.to_lat_lng()
        # The cell center may not exactly match the original point due to
        # quantization to the cell center.
        self.assertAlmostEqual(ll.lat.degrees(), ll2.lat.degrees(), places=5)
        self.assertAlmostEqual(ll.lng.degrees(), ll2.lng.degrees(), places=5)

    def test_get_center_st(self):
        cell = s2.S2CellId.from_face(0)
        st = cell.get_center_st()
        self.assertIsInstance(st, s2.R2Point)

    def test_get_bound_st(self):
        cell = s2.S2CellId.from_face(0)
        bound = cell.get_bound_st()
        self.assertIsInstance(bound, s2.R2Rect)

    def test_get_center_uv(self):
        cell = s2.S2CellId.from_face(0)
        uv = cell.get_center_uv()
        self.assertIsInstance(uv, s2.R2Point)

    def test_get_bound_uv(self):
        cell = s2.S2CellId.from_face(0)
        bound = cell.get_bound_uv()
        self.assertIsInstance(bound, s2.R2Rect)

    def test_get_size_ij(self):
        face = s2.S2CellId.from_face(0)
        self.assertEqual(face.get_size_ij(), 1 << 30)

    def test_get_size_ij_for_level(self):
        self.assertEqual(s2.S2CellId.get_size_ij_for_level(0), 1 << 30)
        self.assertEqual(s2.S2CellId.get_size_ij_for_level(30), 1)

    def test_get_size_st(self):
        face = s2.S2CellId.from_face(0)
        self.assertGreater(face.get_size_st(), 0.0)

    def test_get_center_si_ti(self):
        cell = s2.S2CellId.from_face(0)
        face, si, ti = cell.get_center_si_ti()
        self.assertEqual(face, 0)
        self.assertIsInstance(si, int)
        self.assertIsInstance(ti, int)

    def test_to_face_ij_orientation(self):
        cell = s2.S2CellId.from_face(0)
        face, i, j, orientation = cell.to_face_ij_orientation()
        self.assertEqual(face, 0)

    def test_child_position(self):
        face = s2.S2CellId.from_face(0)
        for pos in range(4):
            child = face.child(pos)
            self.assertEqual(child.child_position(), pos)

    def test_child_position_face_cell_raises(self):
        face = s2.S2CellId.from_face(0)
        with self.assertRaises(ValueError) as cm:
            face.child_position()
        self.assertEqual(str(cm.exception),
                         "Level 0 out of range [1, 30]")

    def test_child_position_at_level(self):
        cell = s2.S2CellId.from_debug_string("0/012")
        self.assertEqual(cell.child_position_at_level(1), 0)
        self.assertEqual(cell.child_position_at_level(2), 1)
        self.assertEqual(cell.child_position_at_level(3), 2)

    def test_to_token(self):
        cell = s2.S2CellId.from_face(0)
        self.assertEqual(cell.to_token(), "1")

    def test_to_string(self):
        cell = s2.S2CellId.from_debug_string("3/02")
        self.assertEqual(cell.to_string(), "3/02")

    # Traversal

    def test_range_min_max(self):
        cell = s2.S2CellId.from_face(0)
        self.assertTrue(cell.range_min() <= cell.range_max())

    def test_contains(self):
        face = s2.S2CellId.from_face(0)
        child = face.child(0)
        self.assertTrue(face.contains(child))
        self.assertFalse(child.contains(face))

    def test_intersects(self):
        face = s2.S2CellId.from_face(0)
        child = face.child(0)
        self.assertTrue(face.intersects(child))
        self.assertTrue(child.intersects(face))
        # Different faces don't intersect.
        face1 = s2.S2CellId.from_face(1)
        self.assertFalse(face.intersects(face1))

    def test_parent(self):
        face = s2.S2CellId.from_face(0)
        child = face.child(0)
        self.assertEqual(child.parent(), face)

    def test_parent_of_face_raises(self):
        face = s2.S2CellId.from_face(0)
        with self.assertRaises(ValueError) as cm:
            face.parent()
        self.assertEqual(str(cm.exception),
                         "Face cell has no parent")

    def test_parent_at_level(self):
        cell = s2.S2CellId.from_debug_string("0/012")
        self.assertEqual(cell.parent_at_level(0), s2.S2CellId.from_face(0))
        self.assertEqual(cell.parent_at_level(3), cell)

    def test_child(self):
        face = s2.S2CellId.from_face(0)
        for pos in range(4):
            child = face.child(pos)
            self.assertEqual(child.level(), 1)

    def test_child_of_leaf_raises(self):
        leaf = s2.S2CellId(s2.S2Point(1.0, 0.0, 0.0))
        with self.assertRaises(ValueError) as cm:
            leaf.child(0)
        self.assertEqual(str(cm.exception),
                         "Leaf cell has no children")

    def test_child_position_out_of_range_raises(self):
        face = s2.S2CellId.from_face(0)
        with self.assertRaises(ValueError) as cm:
            face.child(4)
        self.assertEqual(str(cm.exception),
                         "Child position 4 out of range [0, 3]")

    def test_children(self):
        face = s2.S2CellId.from_face(0)
        children = list(face.children())
        self.assertEqual(len(children), 4)
        for child in children:
            self.assertEqual(child.level(), 1)

    def test_children_positions(self):
        face = s2.S2CellId.from_face(0)
        children = list(face.children())
        for i, child in enumerate(children):
            self.assertEqual(child.child_position(), i)

    def test_children_at_level(self):
        face = s2.S2CellId.from_face(0)
        # Level 2 descendants = 4^2 = 16 cells.
        children = list(face.children(2))
        self.assertEqual(len(children), 16)
        for child in children:
            self.assertEqual(child.level(), 2)

    def test_children_of_leaf_raises(self):
        leaf = s2.S2CellId(s2.S2Point(1.0, 0.0, 0.0))
        with self.assertRaises(ValueError) as cm:
            list(leaf.children())
        self.assertEqual(str(cm.exception),
                         "Leaf cell has no children")

    def test_children_for_loop(self):
        face = s2.S2CellId.from_face(0)
        count = 0
        for child in face.children():
            count += 1
        self.assertEqual(count, 4)

    def test_cells_level_0(self):
        cells = list(s2.S2CellId.cells(0))
        self.assertEqual(len(cells), 6)
        for i, cell in enumerate(cells):
            self.assertEqual(cell.face(), i)

    def test_cells_level_1(self):
        cells = list(s2.S2CellId.cells(1))
        # 6 faces * 4 children = 24 cells.
        self.assertEqual(len(cells), 24)
        for cell in cells:
            self.assertEqual(cell.level(), 1)

    def test_iter_from_cell(self):
        # Iterating from face 3 should yield faces 3, 4, 5.
        face3 = s2.S2CellId.from_face(3)
        cells = list(face3)
        self.assertEqual(len(cells), 3)
        self.assertEqual(cells[0].face(), 3)
        self.assertEqual(cells[1].face(), 4)
        self.assertEqual(cells[2].face(), 5)

    def test_iter_from_last_face(self):
        face5 = s2.S2CellId.from_face(5)
        cells = list(face5)
        self.assertEqual(len(cells), 1)
        self.assertEqual(cells[0].face(), 5)

    def test_reversed_from_cell(self):
        # Reversed from face 3 should yield faces 3, 2, 1, 0.
        face3 = s2.S2CellId.from_face(3)
        cells = list(reversed(face3))
        self.assertEqual(len(cells), 4)
        self.assertEqual(cells[0].face(), 3)
        self.assertEqual(cells[1].face(), 2)
        self.assertEqual(cells[2].face(), 1)
        self.assertEqual(cells[3].face(), 0)

    def test_reversed_from_first_face(self):
        face0 = s2.S2CellId.from_face(0)
        cells = list(reversed(face0))
        self.assertEqual(len(cells), 1)
        self.assertEqual(cells[0].face(), 0)

    def test_children_len(self):
        face = s2.S2CellId.from_face(0)
        self.assertEqual(len(face.children()), 4)
        self.assertEqual(len(face.children(2)), 16)

    def test_children_getitem(self):
        face = s2.S2CellId.from_face(0)
        children = face.children()
        self.assertEqual(children[0].child_position(), 0)
        self.assertEqual(children[3].child_position(), 3)
        # Negative indexing.
        self.assertEqual(children[-1].child_position(), 3)

    def test_children_getitem_out_of_range_raises(self):
        face = s2.S2CellId.from_face(0)
        children = face.children()
        with self.assertRaises(IndexError):
            children[4]
        with self.assertRaises(IndexError):
            children[-5]

    def test_children_contains(self):
        face = s2.S2CellId.from_face(0)
        children = face.children()
        child0 = face.child(0)
        self.assertIn(child0, children)
        # A cell from a different face is not in the range.
        other = s2.S2CellId.from_face(1).child(0)
        self.assertNotIn(other, children)
        # A cell at the wrong level is not in the range.
        self.assertNotIn(face, children)

    def test_children_reversed(self):
        face = s2.S2CellId.from_face(0)
        children = list(face.children())
        rev_children = list(reversed(face.children()))
        self.assertEqual(rev_children, list(reversed(children)))

    def test_cells_len(self):
        self.assertEqual(len(s2.S2CellId.cells(0)), 6)
        self.assertEqual(len(s2.S2CellId.cells(1)), 24)

    def test_cells_getitem(self):
        cells = s2.S2CellId.cells(0)
        self.assertEqual(cells[0].face(), 0)
        self.assertEqual(cells[5].face(), 5)
        self.assertEqual(cells[-1].face(), 5)

    def test_cells_reversed(self):
        cells = list(s2.S2CellId.cells(0))
        rev_cells = list(reversed(s2.S2CellId.cells(0)))
        self.assertEqual(rev_cells, list(reversed(cells)))

    def test_get_common_ancestor_level(self):
        cell1 = s2.S2CellId.from_debug_string("0/012")
        cell2 = s2.S2CellId.from_debug_string("0/013")
        self.assertEqual(cell1.get_common_ancestor_level(cell2), 2)

    def test_get_common_ancestor_level_different_faces(self):
        cell1 = s2.S2CellId.from_face(0)
        cell2 = s2.S2CellId.from_face(1)
        self.assertEqual(cell1.get_common_ancestor_level(cell2), -1)

    def test_get_edge_neighbors(self):
        cell = s2.S2CellId.from_face(0)
        neighbors = cell.get_edge_neighbors()
        self.assertEqual(len(neighbors), 4)

    def test_get_vertex_neighbors(self):
        cell = s2.S2CellId.from_debug_string("0/012")
        neighbors = cell.get_vertex_neighbors(1)
        self.assertGreaterEqual(len(neighbors), 3)
        self.assertLessEqual(len(neighbors), 4)

    def test_get_vertex_neighbors_level_too_high_raises(self):
        cell = s2.S2CellId.from_face(0)
        with self.assertRaises(ValueError):
            cell.get_vertex_neighbors(0)  # level must be < cell level

    def test_get_all_neighbors(self):
        cell = s2.S2CellId.from_debug_string("0/0")
        neighbors = cell.get_all_neighbors(1)
        self.assertGreater(len(neighbors), 0)

    # Operators

    def test_equality(self):
        a = s2.S2CellId.from_face(0)
        b = s2.S2CellId.from_face(0)
        c = s2.S2CellId.from_face(1)
        self.assertTrue(a == b)
        self.assertTrue(a != c)

    def test_comparison(self):
        a = s2.S2CellId.from_face(0)
        b = s2.S2CellId.from_face(1)
        self.assertTrue(a < b)
        self.assertTrue(b > a)
        self.assertTrue(a <= a)
        self.assertTrue(a >= a)

    def test_hash(self):
        a = s2.S2CellId.from_face(0)
        b = s2.S2CellId.from_face(0)
        self.assertEqual(hash(a), hash(b))
        # Can be used in sets/dicts.
        s = {a, b}
        self.assertEqual(len(s), 1)

    # String representation

    def test_repr(self):
        cell = s2.S2CellId.from_face(0)
        self.assertEqual(repr(cell), "S2CellId(0/)")

    def test_str(self):
        cell = s2.S2CellId.from_face(0)
        self.assertEqual(str(cell), "0/")

    def test_repr_child(self):
        cell = s2.S2CellId.from_debug_string("3/02")
        self.assertEqual(repr(cell), "S2CellId(3/02)")


if __name__ == "__main__":
    unittest.main()
