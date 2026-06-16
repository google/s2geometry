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
        self.assertEqual(s2.S2CellId.MAX_LEVEL, 30)

    def test_num_faces(self):
        self.assertEqual(s2.S2CellId.NUM_FACES, 6)

    # Factory methods

    def test_from_face(self):
        for face in range(6):
            cell = s2.S2CellId.from_face(face)
            self.assertEqual(cell.face, face)
            self.assertTrue(cell.is_face())

    def test_from_face_out_of_range_raises(self):
        with self.assertRaises(ValueError) as cm:
            s2.S2CellId.from_face(6)
        self.assertEqual(str(cm.exception), "Face 6 out of range [0, 5]")
        with self.assertRaises(ValueError):
            s2.S2CellId.from_face(-1)

    def test_from_face_pos_level(self):
        cell = s2.S2CellId.from_face_pos_level(0, 0, 0)
        self.assertEqual(cell.level, 0)
        self.assertEqual(cell.face, 0)

    def test_from_face_pos_level_pos_out_of_range_raises(self):
        with self.assertRaises(ValueError):
            s2.S2CellId.from_face_pos_level(0, s2.S2CellId.MAX_POSITION + 1, 0)

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

    def test_from_face_ij(self):
        cell = s2.S2CellId.from_face_ij(0, 0, 0)
        self.assertTrue(cell.is_leaf())

    def test_from_face_ij_i_out_of_range_raises(self):
        with self.assertRaises(ValueError):
            s2.S2CellId.from_face_ij(0, -1, 0)
        with self.assertRaises(ValueError):
            s2.S2CellId.from_face_ij(0, 1 << s2.S2CellId.MAX_LEVEL, 0)

    def test_from_face_ij_j_out_of_range_raises(self):
        with self.assertRaises(ValueError):
            s2.S2CellId.from_face_ij(0, 0, -1)
        with self.assertRaises(ValueError):
            s2.S2CellId.from_face_ij(0, 0, 1 << s2.S2CellId.MAX_LEVEL)

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
            self.assertEqual(s2.S2CellId.from_face(f).face, f)

    def test_level(self):
        face = s2.S2CellId.from_face(0)
        self.assertEqual(face.level, 0)
        child = face.child(0)
        self.assertEqual(child.level, 1)

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
        self.assertAlmostEqual(ll.lat.degrees, ll2.lat.degrees, places=5)
        self.assertAlmostEqual(ll.lng.degrees, ll2.lng.degrees, places=5)

    def test_center_st(self):
        # (s,t)-space covers [0,1] x [0,1] per face; a face cell's center
        # is at (0.5, 0.5).
        st = s2.S2CellId.from_face(0).center_st()
        self.assertAlmostEqual(st.x, 0.5)
        self.assertAlmostEqual(st.y, 0.5)

    def test_bound_st(self):
        # A face cell in (s,t) spans the full [0,1] x [0,1] range.
        bound = s2.S2CellId.from_face(0).bound_st()
        self.assertAlmostEqual(bound.lo.x, 0.0)
        self.assertAlmostEqual(bound.lo.y, 0.0)
        self.assertAlmostEqual(bound.hi.x, 1.0)
        self.assertAlmostEqual(bound.hi.y, 1.0)

    def test_center_uv(self):
        # (u,v)-space covers [-1,1] x [-1,1] per face; a face cell's center
        # is at the origin.
        uv = s2.S2CellId.from_face(0).center_uv()
        self.assertAlmostEqual(uv.x, 0.0)
        self.assertAlmostEqual(uv.y, 0.0)

    def test_bound_uv(self):
        # A face cell in (u,v) spans the full [-1,1] x [-1,1] range.
        bound = s2.S2CellId.from_face(0).bound_uv()
        self.assertAlmostEqual(bound.lo.x, -1.0)
        self.assertAlmostEqual(bound.lo.y, -1.0)
        self.assertAlmostEqual(bound.hi.x, 1.0)
        self.assertAlmostEqual(bound.hi.y, 1.0)

    def test_size_ij(self):
        # In (i,j)-space a face spans 2^kMaxLevel = 2^30 units.
        face = s2.S2CellId.from_face(0)
        self.assertEqual(face.size_ij(), 1 << s2.S2CellId.MAX_LEVEL)

    def test_size_ij_for_level(self):
        # Level 0 cells span 2^kMaxLevel in (i,j); leaf cells span 1.
        self.assertEqual(s2.S2CellId.size_ij_for_level(0), 1 << s2.S2CellId.MAX_LEVEL)
        self.assertEqual(s2.S2CellId.size_ij_for_level(30), 1)

    def test_size_st(self):
        # In (s,t)-space a face spans the full unit interval, so size is 1.0.
        self.assertAlmostEqual(s2.S2CellId.from_face(0).size_st(), 1.0)

    def test_size_st_for_level(self):
        # Level-0 cells span the full unit interval; level-30 cells span 1/2^30.
        self.assertAlmostEqual(s2.S2CellId.size_st_for_level(0), 1.0)
        self.assertAlmostEqual(s2.S2CellId.size_st_for_level(30),
                               1.0 / (1 << 30))

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
        cell = s2.S2CellId.from_face(0).child(0).child(1).child(2)
        self.assertEqual(cell.child_position_at_level(1), 0)
        self.assertEqual(cell.child_position_at_level(2), 1)
        self.assertEqual(cell.child_position_at_level(3), 2)

    def test_to_token(self):
        cell = s2.S2CellId.from_face(0)
        self.assertEqual(cell.to_token(), "1")

    def test_range_min_max_leaf(self):
        # A leaf cell's range is just itself.
        leaf = s2.S2CellId(s2.S2Point(1.0, 0.0, 0.0))
        self.assertEqual(leaf.range_min(), leaf)
        self.assertEqual(leaf.range_max(), leaf)

    def test_range_min_max_covers_children(self):
        # range_min/range_max bound the leaf cells of this cell; every
        # direct child's range must sit within the parent's range.
        face = s2.S2CellId.from_face(0)
        self.assertLess(face.range_min(), face.range_max())
        children = [face.child(i) for i in range(4)]
        for child in children:
            self.assertGreaterEqual(child.range_min(), face.range_min())
            self.assertLessEqual(child.range_max(), face.range_max())
        # The first child's range_min equals the parent's, and likewise
        # the last child's range_max equals the parent's.
        self.assertEqual(children[0].range_min(), face.range_min())
        self.assertEqual(children[-1].range_max(), face.range_max())

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

    def test_common_ancestor_level(self):
        # Two cells sharing the first two Hilbert curve steps on face 0
        # diverge at level 2, so their common ancestor is at level 2.
        base = s2.S2CellId.from_face(0).child(0).child(1)
        cell1 = base.child(2)
        cell2 = base.child(3)
        self.assertEqual(cell1.common_ancestor_level(cell2), 2)

    def test_common_ancestor_level_different_faces(self):
        cell1 = s2.S2CellId.from_face(0)
        cell2 = s2.S2CellId.from_face(1)
        self.assertEqual(cell1.common_ancestor_level(cell2), -1)

    # Traversal

    def test_parent(self):
        face = s2.S2CellId.from_face(0)
        child = face.child(0)
        self.assertEqual(child.parent(), face)

    def test_parent_of_face_raises(self):
        face = s2.S2CellId.from_face(0)
        with self.assertRaises(ValueError) as cm:
            face.parent()
        self.assertEqual(str(cm.exception),
                         "Function invalid for face cells")

    def test_parent_at_level(self):
        cell = s2.S2CellId.from_face(0).child(0).child(1).child(2)
        self.assertEqual(cell.parent_at_level(0), s2.S2CellId.from_face(0))
        self.assertEqual(cell.parent_at_level(3), cell)

    def test_child(self):
        face = s2.S2CellId.from_face(0)
        for pos in range(4):
            child = face.child(pos)
            self.assertEqual(child.level, 1)

    def test_child_of_leaf_raises(self):
        leaf = s2.S2CellId(s2.S2Point(1.0, 0.0, 0.0))
        with self.assertRaises(ValueError) as cm:
            leaf.child(0)
        self.assertEqual(str(cm.exception),
                         "Function invalid for leaf cells")

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
            self.assertEqual(child.level, 1)

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
            self.assertEqual(child.level, 2)

    def test_children_of_leaf_raises(self):
        leaf = s2.S2CellId(s2.S2Point(1.0, 0.0, 0.0))
        with self.assertRaises(ValueError) as cm:
            list(leaf.children())
        self.assertEqual(str(cm.exception),
                         "Function invalid for leaf cells")

    def test_children_at_own_level_raises(self):
        cell = s2.S2CellId.from_face(0)
        with self.assertRaises(ValueError):
            cell.children(0)

    def test_cells_level_0(self):
        cells = list(s2.S2CellId.cells(0))
        self.assertEqual(len(cells), 6)
        for i, cell in enumerate(cells):
            self.assertEqual(cell.face, i)

    def test_cells_level_1(self):
        cells = list(s2.S2CellId.cells(1))
        # 6 faces * 4 children = 24 cells.
        self.assertEqual(len(cells), 24)
        for cell in cells:
            self.assertEqual(cell.level, 1)

    def test_cells_level_out_of_range_raises(self):
        with self.assertRaises(ValueError) as cm:
            s2.S2CellId.cells(31)
        self.assertEqual(str(cm.exception), "Level 31 out of range [0, 30]")

    def test_edge_neighbors(self):
        cell = s2.S2CellId.from_face(0)
        neighbors = cell.edge_neighbors()
        self.assertEqual(len(neighbors), 4)

    def test_vertex_neighbors(self):
        cell = s2.S2CellId.from_face(0).child(0).child(1).child(2)
        neighbors = cell.vertex_neighbors(1)
        self.assertGreaterEqual(len(neighbors), 3)
        self.assertLessEqual(len(neighbors), 4)

    def test_vertex_neighbors_face_cell_raises(self):
        face = s2.S2CellId.from_face(0)
        with self.assertRaises(ValueError):
            face.vertex_neighbors(0)

    def test_vertex_neighbors_level_out_of_range_raises(self):
        cell = s2.S2CellId.from_face(0).child(0).child(1).child(2)
        with self.assertRaises(ValueError):
            cell.vertex_neighbors(3)  # level must be < self.level() (3)

    def test_all_neighbors(self):
        cell = s2.S2CellId.from_face(0).child(0)
        neighbors = cell.all_neighbors(1)
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
        cell = s2.S2CellId.from_face(3).child(0).child(2)
        self.assertEqual(repr(cell), "S2CellId(3/02)")


class TestS2CellIdRange(unittest.TestCase):
    """Test cases for the S2CellIdRange class returned by children()/cells()."""

    def test_len(self):
        face = s2.S2CellId.from_face(0)
        self.assertEqual(len(face.children()), 4)
        self.assertEqual(len(face.children(2)), 16)

    def test_size_matches_len(self):
        face = s2.S2CellId.from_face(0)
        self.assertEqual(face.children().size(), 4)
        self.assertEqual(s2.S2CellId.cells(1).size(), 24)

    def test_size_large_range(self):
        # cells(20) = 6 * 4^20 ≈ 6.6e12; fits in int64 but iterating would
        # be prohibitive. size() returns it without materializing.
        self.assertEqual(s2.S2CellId.cells(20).size(), 6 * (4 ** 20))

    def test_getitem(self):
        face = s2.S2CellId.from_face(0)
        children = face.children()
        self.assertEqual(children[0].child_position(), 0)
        self.assertEqual(children[3].child_position(), 3)
        # Negative indexing.
        self.assertEqual(children[-1].child_position(), 3)

    def test_getitem_out_of_range_raises(self):
        children = s2.S2CellId.from_face(0).children()
        with self.assertRaises(IndexError):
            children[4]
        with self.assertRaises(IndexError):
            children[-5]

    def test_slice_basic(self):
        face = s2.S2CellId.from_face(0)
        children = face.children()
        sliced = children[1:3]
        self.assertEqual(len(sliced), 2)
        self.assertEqual(list(sliced),
                         [face.child(1), face.child(2)])

    def test_slice_open_ended(self):
        face = s2.S2CellId.from_face(0)
        children = face.children()
        self.assertEqual(list(children[:2]),
                         [face.child(0), face.child(1)])
        self.assertEqual(list(children[2:]),
                         [face.child(2), face.child(3)])
        self.assertEqual(list(children[:]), list(children))

    def test_slice_negative_indices(self):
        face = s2.S2CellId.from_face(0)
        children = face.children()
        self.assertEqual(list(children[-2:]),
                         [face.child(2), face.child(3)])

    def test_slice_clamps_out_of_range(self):
        children = s2.S2CellId.from_face(0).children()
        # Stop past the end is clamped.
        self.assertEqual(len(children[0:100]), 4)
        # Start past the end produces an empty range.
        self.assertEqual(len(children[100:]), 0)

    def test_slice_step_not_one_raises(self):
        children = s2.S2CellId.from_face(0).children()
        with self.assertRaises(ValueError) as cm:
            children[::2]
        self.assertIn("step 1", str(cm.exception))
        with self.assertRaises(ValueError):
            children[::-1]

    def test_slice_returns_range(self):
        # Slicing returns another S2CellIdRange, not a list.
        children = s2.S2CellId.from_face(0).children()
        sliced = children[1:3]
        self.assertIsInstance(sliced, s2.S2CellIdRange)

    def test_contains(self):
        face = s2.S2CellId.from_face(0)
        children = face.children()
        child0 = face.child(0)
        self.assertIn(child0, children)
        # A cell from a different face is not in the range.
        other = s2.S2CellId.from_face(1).child(0)
        self.assertNotIn(other, children)
        # A cell at the wrong level is not in the range.
        self.assertNotIn(face, children)

    def test_reversed(self):
        face = s2.S2CellId.from_face(0)
        children = list(face.children())
        rev_children = list(reversed(face.children()))
        self.assertEqual(rev_children, list(reversed(children)))


if __name__ == "__main__":
    unittest.main()
