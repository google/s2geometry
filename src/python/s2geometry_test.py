#
# Copyright 2006 Google Inc. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS-IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import math
import unittest
from collections import defaultdict

import s2geometry as s2


class PyWrapS2TestCase(unittest.TestCase):

  def testS2PointFromRawToNamedCorrectly(self):
      x = 1.0
      y = 2.0
      z = 3.0
      point = s2.S2Point_FromRaw(x, y, z)
      self.assertEqual(x, point.x())
      self.assertEqual(y, point.y())
      self.assertEqual(z, point.z())

  def testContainsIsWrappedCorrectly(self):
    london = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(51.3368602, 0.4931979),
                             s2.S2LatLng.FromDegrees(51.7323965, 0.1495211))
    e14lj = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(51.5213527, -0.0476026),
                            s2.S2LatLng.FromDegrees(51.5213527, -0.0476026))
    self.assertTrue(london.Contains(e14lj))

  def testS2CellIdEqualsIsWrappedCorrectly(self):
    london = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    cell = s2.S2CellId(london)
    same_cell = s2.S2CellId(london)
    self.assertEqual(cell, same_cell)

  def testS2CellIdComparsionIsWrappedCorrectly(self):
    london = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    cell = s2.S2CellId(london)
    self.assertLess(cell, cell.next())
    self.assertGreater(cell.next(), cell)

  def testS2CellIdFromToTokenIsWrappedCorrectly(self):
    cell = s2.S2CellId.FromToken("487604c489f841c3")
    self.assertEqual(cell.ToToken(), "487604c489f841c3")
    self.assertEqual(cell.id(), 0x487604c489f841c3)

    cell = s2.S2CellId.FromToken("487")
    self.assertEqual(cell.ToToken(), "487")
    self.assertEqual(cell.id(), 0x4870000000000000)

    cell = s2.S2CellId.FromToken("this is invalid")
    self.assertEqual(cell.ToToken(), "X")
    self.assertEqual(cell.id(), 0)

  def testS2CellIdGetEdgeNeighborsIsWrappedCorrectly(self):
    cell = s2.S2CellId(0x466d319000000000)
    expected_neighbors = [s2.S2CellId(0x466d31b000000000),
                          s2.S2CellId(0x466d317000000000),
                          s2.S2CellId(0x466d323000000000),
                          s2.S2CellId(0x466d31f000000000)]
    neighbors = cell.GetEdgeNeighbors()
    self.assertCountEqual(neighbors, expected_neighbors)

  def testS2CellIdGetVertexNeighborsIsWrappedCorrectly(self):
    cell = s2.S2CellId(0x466d319000000000)
    expected_neighbors = [s2.S2CellId(0x466d31c000000000),
                          s2.S2CellId(0x466d314000000000),
                          s2.S2CellId(0x466d324000000000),
                          s2.S2CellId(0x466d33c000000000)]
    self.assertEqual(cell.level(), 12)
    # Requires level < cell.level.
    neighbors = cell.GetVertexNeighbors(11)
    self.assertCountEqual(neighbors, expected_neighbors)

  def testS2CellIdGetAllNeighborsIsWrappedCorrectly(self):
    cell = s2.S2CellId(0x466d319000000000)
    expected_neighbors = [s2.S2CellId(0x466d31d000000000),
                          s2.S2CellId(0x466d311000000000),
                          s2.S2CellId(0x466d31b000000000),
                          s2.S2CellId(0x466d323000000000),
                          s2.S2CellId(0x466d31f000000000),
                          s2.S2CellId(0x466d317000000000),
                          s2.S2CellId(0x466d321000000000),
                          s2.S2CellId(0x466d33d000000000)]
    self.assertEqual(cell.level(), 12)
    # Requires level >= cell.level.
    neighbors = cell.GetAllNeighbors(12)
    self.assertCountEqual(neighbors, expected_neighbors)

  def testS2CellIdChild(self):
    valid = s2.S2CellId(0x89c259c000000000)
    invalid = s2.S2CellId(0)
    self.assertTrue(valid.is_valid())
    self.assertFalse(invalid.is_valid())

    self.assertEqual(valid.child(0).parent().id(), valid.id())

    with self.assertRaises(ValueError):
      valid.child(-1)

    with self.assertRaises(ValueError):
      valid.child(4)

    with self.assertRaises(ValueError):
      invalid.child(0)

    leaf = s2.S2CellId(s2.S2LatLng.FromDegrees(10.0, 20.0))
    with self.assertRaises(ValueError):
      leaf.child(0)

  def testS2CellIdChildLevelIsWrappedCorrectly(self):
    cell = s2.S2CellId(0x876bec2688e50000)
    self.assertEqual(cell.child_position(3), 2)
    with self.assertRaises(ValueError):
      cell.child_position(-1)
    with self.assertRaises(ValueError):
      cell.child_position(0)
    with self.assertRaises(ValueError):
      cell.child_position(40)

  def testS2CellIdContainsInvalidRaises(self):
    valid = s2.S2CellId(0x89c259c000000000)
    invalid = s2.S2CellId(0)
    self.assertTrue(valid.is_valid())
    self.assertFalse(invalid.is_valid())

    self.assertTrue(valid.contains(valid))

    with self.assertRaises(ValueError):
      valid.contains(invalid)

    with self.assertRaises(ValueError):
      invalid.contains(valid)

  def testS2CellIdGetAllNeighborsIsWrappedCorrectly(self):
    cell = s2.S2CellId(0x6aa7590000000000)
    expected_neighbors = (s2.S2CellId(0x2ab3530000000000),
                          s2.S2CellId(0x2ab34b0000000000),
                          s2.S2CellId(0x2ab34d0000000000),
                          s2.S2CellId(0x6aa75b0000000000),
                          s2.S2CellId(0x6aa7570000000000),
                          s2.S2CellId(0x6aa75f0000000000),
                          s2.S2CellId(0x6aa7510000000000),
                          s2.S2CellId(0x6aa75d0000000000))
    neighbors = cell.GetAllNeighbors(cell.level())
    self.assertEqual(neighbors, expected_neighbors)

  def testS2CellIdIntersectsIsTrueForOverlap(self):
    cell1 = s2.S2CellId(0x89c259c000000000)
    cell2 = s2.S2CellId(0x89c2590000000000)
    self.assertTrue(cell1.intersects(cell2))

  def testS2CellIdIntersectsIsFalseForNonOverlap(self):
    cell1 = s2.S2CellId(0x89c259c000000000)
    cell2 = s2.S2CellId(0x89e83d0000000000)
    self.assertFalse(cell1.intersects(cell2))

  def testS2CellIdIntersectsInvalidRaises(self):
    valid = s2.S2CellId(0x89c259c000000000)
    invalid = s2.S2CellId(0)
    self.assertTrue(valid.is_valid())
    self.assertFalse(invalid.is_valid())

    with self.assertRaises(ValueError):
      valid.intersects(invalid)

    with self.assertRaises(ValueError):
      invalid.intersects(valid)

  def testS2CellIdLevel(self):
    leaf = s2.S2CellId(s2.S2LatLng.FromDegrees(10.0, 20.0))
    self.assertEqual(leaf.level(), 30)

    with self.assertRaises(ValueError):
      s2.S2CellId(0).level()

  def testS2CellIdParent(self):
    leaf = s2.S2CellId(s2.S2LatLng.FromDegrees(10.0, 20.0))
    self.assertEqual(leaf.level(), 30)
    self.assertEqual(leaf.parent().level(), 29)

    level8 = leaf.parent(8)
    self.assertEqual(level8.level(), 8)

    self.assertEqual(level8.parent(0).level(), 0)

    # Error to have negative level.
    with self.assertRaises(ValueError):
      level8.parent(-1)

    # Same level is ok.
    self.assertEqual(level8.parent(8).level(), 8)

    # Error to ask for parent with lower level.
    with self.assertRaises(ValueError):
      level8.parent(9)

    # Parent of invalid is an error
    with self.assertRaises(ValueError):
      s2.S2CellId(0).parent()

  def testS2HashingIsWrappedCorrectly(self):
    london = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    cell = s2.S2CellId(london)
    same_cell = s2.S2CellId(london)
    self.assertEqual(hash(cell), hash(same_cell))

  def testCovererIsWrappedCorrectly(self):
    london = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(51.3368602, 0.4931979),
                             s2.S2LatLng.FromDegrees(51.7323965, 0.1495211))
    e14lj = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(51.5213527, -0.0476026),
                            s2.S2LatLng.FromDegrees(51.5213527, -0.0476026))
    coverer = s2.S2RegionCoverer()
    coverer.set_max_cells(6)
    self.assertEqual(6, coverer.max_cells())
    covering = coverer.GetCovering(e14lj)
    self.assertLessEqual(len(covering), 6)
    for cellid in covering:
      self.assertTrue(london.Contains(s2.S2Cell(cellid)))
    interior = coverer.GetInteriorCovering(e14lj)
    for cellid in interior:
      self.assertTrue(london.Contains(s2.S2Cell(cellid)))

  def testS2CellUnionIsWrappedCorrectly(self):
    cell_union = s2.S2CellUnion()
    cell_union.Init([0x466d319000000000, 0x466d31b000000000])
    self.assertEqual(cell_union.num_cells(), 2)
    trondheim = s2.S2LatLng.FromDegrees(63.431052, 10.395083)
    self.assertTrue(cell_union.Contains(s2.S2CellId(trondheim)))

    # Init() calls Normalized, so cell_ids() are normalized.
    cell_union2 = s2.S2CellUnion.FromNormalized(cell_union.cell_ids())
    # There is no S2CellUnion::Equals, and cell_ids is a non-iterable
    # SWIG object, so just perform the same checks again.
    self.assertEqual(cell_union2.num_cells(), 2)
    self.assertTrue(cell_union2.Contains(s2.S2CellId(trondheim)))

  def testS2PolygonIsWrappedCorrectly(self):
    london = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    polygon = s2.S2Polygon(s2.S2Cell(s2.S2CellId(london)))
    self.assertEqual(polygon.num_loops(), 1)
    point = london.ToPoint()
    self.assertTrue(polygon.Contains(point))

  def testS2LoopIsWrappedCorrectly(self):
    london = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    polygon = s2.S2Polygon(s2.S2Cell(s2.S2CellId(london)))
    loop = polygon.loop(0)
    self.assertTrue(loop.IsValid())
    self.assertEqual(0, loop.depth())
    self.assertFalse(loop.is_hole())
    self.assertEqual(4, loop.num_vertices())
    self.assertTrue(loop.IsNormalized())
    point = london.ToPoint()
    self.assertTrue(loop.Contains(point))

  def testS2LoopUsesValueEquality(self):
    self.assertEqual(s2.S2Loop(), s2.S2Loop())

  def testS2PolygonCopiesLoopInConstructorBecauseItTakesOwnership(self):
    london = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    loop = s2.S2Loop(s2.S2Cell(s2.S2CellId(london)))
    s2.S2Polygon(loop)

  def testS2LoopAreaIsWrappedCorrectly(self):
    london = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    loop = s2.S2Loop(s2.S2Cell(s2.S2CellId(london)))
    equivalent_polygon = s2.S2Polygon(loop)
    self.assertAlmostEqual(loop.GetArea(), equivalent_polygon.GetArea())

  def testS2PolygonInitNestedIsWrappedCorrectly(self):
    london = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    small_loop = s2.S2Loop(s2.S2Cell(s2.S2CellId(london)))
    big_loop = s2.S2Loop(s2.S2Cell(s2.S2CellId(london).parent(1)))
    polygon = s2.S2Polygon()
    polygon.InitNested([big_loop, small_loop])

  def testS2PolygonInitNestedWithIncorrectTypeIsWrappedCorrectly(self):
    london = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    loop = s2.S2Loop(s2.S2Cell(s2.S2CellId(london)))
    polygon = s2.S2Polygon()
    with self.assertRaises(TypeError):
      polygon.InitNested([loop, s2.S2CellId()])

  def testS2PolygonGetAreaIsWrappedCorrectly(self):
    # Cell at level 10 containing central London.

    london_level_10 = s2.S2CellId(
        s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)).parent(10)
    polygon = s2.S2Polygon(s2.S2Cell(london_level_10))
    # Because S2Cell.ExactArea() isn't swigged, compare S2Polygon.GetArea() with
    # S2CellUnion.ExactArea().
    cell_union = s2.S2CellUnion()
    cell_union.Init([london_level_10.id()])
    self.assertAlmostEqual(cell_union.ExactArea(), polygon.GetArea(), places=10)

  def testS2PolygonGetOverlapFractions(self):
    # Matches S2Polygon, OverlapFractions test from cs/s2polygon_test.cc
    a = s2.S2Polygon()
    b = s2.S2Polygon()
    r1, r2 = s2.S2Polygon.GetOverlapFractions(a, b)
    self.assertAlmostEqual(1.0, r1)
    self.assertAlmostEqual(1.0, r2)

    def verts2loop(vs):
      loop = s2.S2Loop()
      loop.Init([s2.S2LatLng.FromDegrees(*v).ToPoint() for v in vs])
      return loop

    loop1verts = [(-10, 10), (0, 10), (0, -10), (-10, -10), (-10, 0)]
    b = s2.S2Polygon(verts2loop(loop1verts))
    r1, r2 = s2.S2Polygon.GetOverlapFractions(a, b)
    self.assertAlmostEqual(1.0, r1)
    self.assertAlmostEqual(0.0, r2)

    loop2verts = [(-10, 0), (10, 0), (10, -10), (-10, -10)]
    a = s2.S2Polygon(verts2loop(loop2verts))
    r1, r2 = s2.S2Polygon.GetOverlapFractions(a, b)
    self.assertAlmostEqual(0.5, r1)
    self.assertAlmostEqual(0.5, r2)

  def testGetS2LatLngVertexIsWrappedCorrectly(self):
    london = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    polygon = s2.S2Polygon(s2.S2Cell(s2.S2CellId(london)))
    loop = polygon.loop(0)
    first_vertex = loop.GetS2LatLngVertex(0)
    self.assertIsInstance(first_vertex, s2.S2LatLng)
    self.assertEqual("51.500152,-0.126235", first_vertex.ToStringInDegrees())
    second_vertex = loop.GetS2LatLngVertex(1)
    self.assertIsInstance(second_vertex, s2.S2LatLng)
    self.assertEqual("51.500153,-0.126235", second_vertex.ToStringInDegrees())

  def testGetLastDescendant(self):
    def verts2loop(vs):
      loop = s2.S2Loop()
      loop.Init([s2.S2LatLng.FromDegrees(*v).ToPoint() for v in vs])
      return loop

    loop1 = verts2loop([(0, 0), (0, 10), (10, 10), (10, 0)])    # Shell
    loop2 = verts2loop([(2, 2), (2, 5), (5, 5), (5, 2)])        # Hole
    loop3 = verts2loop([(0, 20), (0, 30), (10, 30), (10, 20)])  # Another shell

    polygon = s2.S2Polygon()
    polygon.InitNested([loop1, loop2, loop3])

    self.assertEqual(1, polygon.GetLastDescendant(0))
    self.assertEqual(1, polygon.GetLastDescendant(1))
    self.assertEqual(2, polygon.GetLastDescendant(2))

  def testS2PolylineInitFromS2LatLngs(self):
    e7_10deg = 0x5f5e100
    list_ll = []
    for lat, lng in [(0, 0), (0, e7_10deg), (e7_10deg, e7_10deg)]:
      list_ll.append(s2.S2LatLng.FromE7(lat, lng))
    line = s2.S2Polyline()
    line.InitFromS2LatLngs(list_ll)
    self.assertAlmostEqual(20.0, line.GetLength().degrees())

  def testS2PolylineInitFromS2Points(self):
    e7_10deg = 0x5f5e100
    list_points = []
    for lat, lng in [(0, 0), (0, e7_10deg), (e7_10deg, e7_10deg)]:
      list_points.append(s2.S2LatLng.FromE7(lat, lng).ToPoint())
    line = s2.S2Polyline()
    line.InitFromS2Points(list_points)
    self.assertAlmostEqual(20.0, line.GetLength().degrees())

  def testS2PolylineUsesValueEquality(self):
    self.assertEqual(s2.S2Polyline(), s2.S2Polyline())

  def testS2PointsCanBeNormalized(self):
    line = s2.S2Polyline()
    line.InitFromS2LatLngs([s2.S2LatLng.FromDegrees(37.794484, -122.394871),
                            s2.S2LatLng.FromDegrees(37.762699, -122.435158)])
    self.assertNotAlmostEqual(line.GetCentroid().Norm(), 1.0)
    self.assertAlmostEqual(line.GetCentroid().Normalize().Norm(), 1.0)

  def testS1AngleComparsionIsWrappedCorrectly(self):
    ten_degrees = s2.S1Angle.Degrees(10)
    one_hundred_degrees = s2.S1Angle.Degrees(100)
    self.assertLess(ten_degrees, one_hundred_degrees)
    self.assertGreater(one_hundred_degrees, ten_degrees)

  def testS2PolygonIntersectsWithPolyline(self):
    london = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    polygon = s2.S2Polygon(s2.S2Cell(s2.S2CellId(london).parent(15)))
    line = s2.S2Polyline()
    line.InitFromS2LatLngs([s2.S2LatLng.FromDegrees(51.5, -0.128),
                            s2.S2LatLng.FromDegrees(51.5, -0.125)])
    intersections = polygon.IntersectWithPolyline(line)
    self.assertEqual(1, len(intersections))

  def testS2PolygonBoundaryNearIsSame(self):
    london_1 = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    polygon_1 = s2.S2Polygon(s2.S2Loop(s2.S2Cell(s2.S2CellId(london_1))))
    london_2 = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    polygon_2 = s2.S2Polygon(s2.S2Loop(s2.S2Cell(s2.S2CellId(london_2))))
    self.assertTrue(polygon_1.BoundaryNear(polygon_2))

  def testS2PolygonBoundaryNearIsTotallyDifferent(self):
    london = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    polygon_1 = s2.S2Polygon(s2.S2Loop(s2.S2Cell(s2.S2CellId(london))))
    seattle = s2.S2LatLng.FromDegrees(47.6062, -122.3321)
    polygon_2 = s2.S2Polygon(s2.S2Loop(s2.S2Cell(s2.S2CellId(seattle))))
    self.assertFalse(polygon_1.BoundaryNear(polygon_2))

  def testS2PolygonBoundaryNearIsNear(self):

    def verts2loop(vs):
      loop = s2.S2Loop()
      loop.Init([s2.S2LatLng.FromDegrees(*v).ToPoint() for v in vs])
      return loop

    vertices_1 = [(-10, 10), (0, 10), (0, -10), (-10, -10), (-10, 0)]
    polygon_1 = s2.S2Polygon(verts2loop(vertices_1))
    vertices_2 = [(-10, 10), (0, 10), (0, -10.1), (-10, -10), (-10, 0)]
    polygon_2 = s2.S2Polygon(verts2loop(vertices_2))
    self.assertTrue(polygon_1.BoundaryNear(polygon_2, s2.S1Angle.Degrees(1)))

  def testS2PolygonUsesValueEquality(self):
    self.assertEqual(s2.S2Polygon(), s2.S2Polygon())

  def testCrossingSign(self):
    a = s2.S2LatLng.FromDegrees(-1, 0).ToPoint()
    b = s2.S2LatLng.FromDegrees(1, 0).ToPoint()
    c = s2.S2LatLng.FromDegrees(0, -1).ToPoint()
    d = s2.S2LatLng.FromDegrees(0, 1).ToPoint()
    # SWIG flattens namespaces, so this is just s2.CrossingSign,
    # not s2.S2.CrossingSign.
    self.assertEqual(1, s2.CrossingSign(a, b, c, d))

  def testGetIntersection(self):
    a = s2.S2LatLng.FromDegrees(-1, 0).ToPoint()
    b = s2.S2LatLng.FromDegrees(1, 0).ToPoint()
    c = s2.S2LatLng.FromDegrees(0, -1).ToPoint()
    d = s2.S2LatLng.FromDegrees(0, 1).ToPoint()
    # SWIG namespace flattening as above.
    intersection = s2.GetIntersection(a, b, c, d)
    self.assertEqual(
        "0.000000,0.000000", s2.S2LatLng(intersection).ToStringInDegrees())

  def testS2CellDistance(self):
    # Level-0 cell (i.e. face) centered at (0, 0)
    cell = s2.S2Cell(s2.S2CellId(0x1000000000000000))

    p1 = s2.S2LatLng.FromDegrees(0, 0).ToPoint()
    self.assertTrue(cell.Contains(p1))
    d1 = cell.GetDistance(p1).ToAngle().degrees()
    # Inside, so distance is 0, but boundary distance is not.
    self.assertEqual(0.0, d1)
    bd1 = cell.GetBoundaryDistance(p1).ToAngle().degrees()
    self.assertEqual(45.0, bd1)

    p2 = s2.S2LatLng.FromDegrees(0, 90).ToPoint()
    self.assertFalse(cell.Contains(p2))
    d2 = cell.GetDistance(p2).ToAngle().degrees()
    self.assertAlmostEqual(45.0, d2)
    bd2 = cell.GetBoundaryDistance(p2).ToAngle().degrees()
    # Outside, so distance and boundary distance are the same.
    self.assertAlmostEqual(45.0, bd2)

  def testS2Rotate(self):
    mtv_a = s2.S2LatLng.FromDegrees(37.4402777, -121.9638888).ToPoint()
    mtv_b = s2.S2LatLng.FromDegrees(37.3613888, -121.9283333).ToPoint()
    angle = s2.S1Angle.Radians(0.039678)
    point = s2.Rotate(mtv_a, mtv_b, angle)
    self.assertEqual("37.439095,-121.967802",
                     s2.S2LatLng(point).ToStringInDegrees())

  def testS2TurnAngle(self):
    mtv_a = s2.S2LatLng.FromDegrees(37.4402777, -121.9638888).ToPoint()
    mtv_b = s2.S2LatLng.FromDegrees(37.3613888, -121.9283333).ToPoint()
    mtv_c = s2.S2LatLng.FromDegrees(37.3447222, -122.0308333).ToPoint()
    angle = s2.TurnAngle(mtv_a, mtv_b, mtv_c)
    self.assertAlmostEqual(-1.7132025, angle)

  def testEncodeDecode(self):
    london = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    polygon = s2.S2Polygon(s2.S2Cell(s2.S2CellId(london).parent(15)))
    self.assertEqual(polygon.num_loops(), 1)

    encoder = s2.Encoder()
    polygon.Encode(encoder)

    encoded = encoder.buffer()
    decoder = s2.Decoder(encoded)
    decoded_polygon = s2.S2Polygon()
    self.assertTrue(decoded_polygon.Decode(decoder))

    self.assertEqual(decoded_polygon.num_loops(), 1)
    self.assertTrue(decoded_polygon.Equals(polygon))

  def testS2CapRegion(self):
    center = s2.S2LatLng.FromDegrees(2.0, 3.0).ToPoint()
    cap = s2.S2Cap(center, s2.S1Angle.Degrees(1.0))

    inside = s2.S2LatLng.FromDegrees(2.1, 2.9).ToPoint()
    outside = s2.S2LatLng.FromDegrees(0.0, 0.0).ToPoint()
    self.assertTrue(cap.Contains(inside))
    self.assertFalse(cap.Contains(outside))
    self.assertTrue(cap.Contains(s2.S2Cell(inside)))
    self.assertFalse(cap.Contains(s2.S2Cell(outside)))
    self.assertTrue(cap.MayIntersect(s2.S2Cell(inside)))
    self.assertFalse(cap.MayIntersect(s2.S2Cell(outside)))

    self.assertTrue(cap.ApproxEquals(cap.GetCapBound()))

    rect_bound = cap.GetRectBound()
    self.assertTrue(rect_bound.Contains(inside))
    self.assertFalse(rect_bound.Contains(outside))

  def testS2LatLngRectRegion(self):
    rect = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(1.0, 2.0),
                           s2.S2LatLng.FromDegrees(3.0, 4.0))

    inside = s2.S2LatLng.FromDegrees(2.0, 3.0).ToPoint()
    outside = s2.S2LatLng.FromDegrees(0.0, 0.0).ToPoint()

    self.assertTrue(rect.Contains(inside))
    self.assertFalse(rect.Contains(outside))
    self.assertTrue(rect.Contains(s2.S2Cell(inside)))
    self.assertFalse(rect.Contains(s2.S2Cell(outside)))
    self.assertTrue(rect.MayIntersect(s2.S2Cell(inside)))
    self.assertFalse(rect.MayIntersect(s2.S2Cell(outside)))

    cap_bound = rect.GetCapBound()
    self.assertTrue(cap_bound.Contains(inside))
    self.assertFalse(cap_bound.Contains(outside))

    self.assertTrue(rect.ApproxEquals(rect.GetRectBound()))

  def testS2CellRegion(self):
    cell = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8))

    inside = s2.S2LatLng.FromDegrees(3.0, 4.0).ToPoint()
    outside = s2.S2LatLng.FromDegrees(30.0, 40.0).ToPoint()

    self.assertTrue(cell.Contains(inside))
    self.assertFalse(cell.Contains(outside))
    self.assertTrue(cell.Contains(s2.S2Cell(inside)))
    self.assertFalse(cell.Contains(s2.S2Cell(outside)))
    self.assertTrue(cell.MayIntersect(s2.S2Cell(inside)))
    self.assertFalse(cell.MayIntersect(s2.S2Cell(outside)))

    cap_bound = cell.GetCapBound()
    self.assertTrue(cap_bound.Contains(inside))
    self.assertFalse(cap_bound.Contains(outside))

    rect_bound = cell.GetRectBound()
    self.assertTrue(rect_bound.Contains(inside))
    self.assertFalse(rect_bound.Contains(outside))

  def testS2CellUnionRegion(self):
    cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8)
    cell_union = s2.S2CellUnion()
    cell_union.Init([cell_id.id()])

    inside = s2.S2LatLng.FromDegrees(3.0, 4.0).ToPoint()
    outside = s2.S2LatLng.FromDegrees(30.0, 40.0).ToPoint()

    self.assertTrue(cell_union.Contains(inside))
    self.assertFalse(cell_union.Contains(outside))
    self.assertTrue(cell_union.Contains(s2.S2Cell(inside)))
    self.assertFalse(cell_union.Contains(s2.S2Cell(outside)))
    self.assertTrue(cell_union.MayIntersect(s2.S2Cell(inside)))
    self.assertFalse(cell_union.MayIntersect(s2.S2Cell(outside)))

    cap_bound = cell_union.GetCapBound()
    self.assertTrue(cap_bound.Contains(inside))
    self.assertFalse(cap_bound.Contains(outside))

    rect_bound = cell_union.GetRectBound()
    self.assertTrue(rect_bound.Contains(inside))
    self.assertFalse(rect_bound.Contains(outside))

  def testS2CellUnionEmpty(self):
    empty_cell_union = s2.S2CellUnion()
    self.assertTrue(empty_cell_union.empty())

    cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8)
    cell_union = s2.S2CellUnion()
    cell_union.Init([cell_id.id()])
    self.assertFalse(cell_union.empty())

  def testS2CellUnionIntersectionWithS2CellUnion(self):
    cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0))
    cell_union = s2.S2CellUnion()
    cell_union.Init([cell_id.id()])

    # No intersection.
    outside_cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(5.0, 6.0))
    outside_cell_union = s2.S2CellUnion()
    outside_cell_union.Init([outside_cell_id.id()])
    empty_intersection = cell_union.Intersection(outside_cell_union)
    self.assertTrue(empty_intersection.empty())

    # Complete overlap.
    self_intersection = cell_union.Intersection(cell_union)
    self.assertTrue(self_intersection.Contains(cell_union))
    self.assertTrue(cell_union.Contains(self_intersection))

    # Some intersection.
    joint_cell_union = s2.S2CellUnion()
    joint_cell_union.Init([cell_id.id(), outside_cell_id.id()])
    outside_intersection = joint_cell_union.Intersection(outside_cell_union)
    self.assertTrue(outside_intersection.Contains(outside_cell_id))
    self.assertFalse(outside_intersection.Contains(cell_id))

  def testS2CellUnionIntersectionWithS2CellId(self):
    cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0))
    cell_union = s2.S2CellUnion()
    cell_union.Init([cell_id.id()])

    # No intersection.
    outside_cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(4.0, 5.0))
    empty_intersection = cell_union.Intersection(outside_cell_id)
    self.assertTrue(empty_intersection.empty())

    # Complete overlap.
    intersection = cell_union.Intersection(cell_id)
    self.assertTrue(intersection.Contains(cell_id))

    # Some intersection.
    joint_cell_union = s2.S2CellUnion()
    joint_cell_union.Init([cell_id.id(), outside_cell_id.id()])
    outside_intersection = joint_cell_union.Intersection(outside_cell_id)
    self.assertTrue(outside_intersection.Contains(outside_cell_id))
    self.assertFalse(outside_intersection.Contains(cell_id))

  def testS2CellUnionIsNormalized(self):
    empty_cell_union = s2.S2CellUnion()
    self.assertTrue(empty_cell_union.IsNormalized())

    london = s2.S2LatLng.FromDegrees(51.5001525, -0.1262355)
    london_cell_id = s2.S2CellId(london)
    normalized_union = s2.S2CellUnion()
    normalized_union.Init([london_cell_id.id()])
    self.assertTrue(normalized_union.IsNormalized())

  def testS2CellUnionNormalizeS2CellUnion(self):
    empty_cell_union = s2.S2CellUnion()
    empty_cell_union.NormalizeS2CellUnion()
    self.assertTrue(empty_cell_union.IsNormalized())

    cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8)
    cell_union = s2.S2CellUnion()
    cell_union.Init([cell_id.id()])
    cell_union.NormalizeS2CellUnion()
    self.assertTrue(cell_union.IsNormalized())

  def testS2LoopRegion(self):
    cell = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8))
    loop = s2.S2Loop(cell)

    inside = s2.S2LatLng.FromDegrees(3.0, 4.0).ToPoint()
    outside = s2.S2LatLng.FromDegrees(30.0, 40.0).ToPoint()

    self.assertTrue(loop.Contains(inside))
    self.assertFalse(loop.Contains(outside))
    self.assertTrue(loop.Contains(s2.S2Cell(inside)))
    self.assertFalse(loop.Contains(s2.S2Cell(outside)))
    self.assertTrue(loop.MayIntersect(s2.S2Cell(inside)))
    self.assertFalse(loop.MayIntersect(s2.S2Cell(outside)))

    cap_bound = loop.GetCapBound()
    self.assertTrue(cap_bound.Contains(inside))
    self.assertFalse(cap_bound.Contains(outside))

    rect_bound = loop.GetRectBound()
    self.assertTrue(rect_bound.Contains(inside))
    self.assertFalse(rect_bound.Contains(outside))

  def testS2PolygonRegion(self):
    cell = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8))
    polygon = s2.S2Polygon(cell)

    inside = s2.S2LatLng.FromDegrees(3.0, 4.0).ToPoint()
    outside = s2.S2LatLng.FromDegrees(30.0, 40.0).ToPoint()

    self.assertTrue(polygon.Contains(inside))
    self.assertFalse(polygon.Contains(outside))
    self.assertTrue(polygon.Contains(s2.S2Cell(inside)))
    self.assertFalse(polygon.Contains(s2.S2Cell(outside)))
    self.assertTrue(polygon.MayIntersect(s2.S2Cell(inside)))
    self.assertFalse(polygon.MayIntersect(s2.S2Cell(outside)))

    cap_bound = polygon.GetCapBound()
    self.assertTrue(cap_bound.Contains(inside))
    self.assertFalse(cap_bound.Contains(outside))

    rect_bound = polygon.GetRectBound()
    self.assertTrue(rect_bound.Contains(inside))
    self.assertFalse(rect_bound.Contains(outside))

  def testS2PolylineRegion(self):
    polyline = s2.S2Polyline()
    polyline.InitFromS2LatLngs([s2.S2LatLng.FromDegrees(0.0, 0.0),
                                s2.S2LatLng.FromDegrees(1.0, 1.0)])

    # Contains(S2Point) always return false.
    self.assertFalse(
        polyline.Contains(s2.S2LatLng.FromDegrees(0.0, 0.0).ToPoint()))
    self.assertFalse(
        polyline.Contains(s2.S2Cell(s2.S2LatLng.FromDegrees(0.0, 0.0))))

    self.assertTrue(
        polyline.MayIntersect(s2.S2Cell(s2.S2LatLng.FromDegrees(0.0, 0.0))))
    self.assertFalse(
        polyline.MayIntersect(s2.S2Cell(s2.S2LatLng.FromDegrees(3.0, 4.0))))

    cap_bound = polyline.GetCapBound()
    self.assertTrue(
        cap_bound.Contains(s2.S2LatLng.FromDegrees(0.0, 0.0).ToPoint()))
    self.assertFalse(
        cap_bound.Contains(s2.S2LatLng.FromDegrees(2.0, 2.0).ToPoint()))

    rect_bound = polyline.GetRectBound()
    self.assertTrue(
        rect_bound.Contains(s2.S2LatLng.FromDegrees(0.0, 0.0).ToPoint()))
    self.assertFalse(
        rect_bound.Contains(s2.S2LatLng.FromDegrees(2.0, 2.0).ToPoint()))

  def testS2CellIdCenterSiTi(self):
    cell = s2.S2CellId.FromFacePosLevel(3, 0x12345678, s2.S2CellId.kMaxLevel)

    # Check that the (si, ti) coordinates of the center end in a
    # 1 followed by (30 - level) 0s.

    # Leaf level, 30.
    face, si, ti = cell.GetCenterSiTi()
    self.assertEqual(3, face)
    self.assertEqual(1 << 0, si & 1)
    self.assertEqual(1 << 0, ti & 1)

    # Level 29.
    face, si, ti = cell.parent(s2.S2CellId.kMaxLevel - 1).GetCenterSiTi()
    self.assertEqual(3, face)
    self.assertEqual(1 << 1, si & 3)
    self.assertEqual(1 << 1, ti & 3)

    # Level 28.
    face, si, ti = cell.parent(s2.S2CellId.kMaxLevel - 2).GetCenterSiTi()
    self.assertEqual(3, face)
    self.assertEqual(1 << 2, si & 7)
    self.assertEqual(1 << 2, ti & 7)

    # Level 20.
    face, si, ti = cell.parent(s2.S2CellId.kMaxLevel - 10).GetCenterSiTi()
    self.assertEqual(3, face)
    self.assertEqual(1 << 10, si & ((1 << 11) - 1))
    self.assertEqual(1 << 10, ti & ((1 << 11) - 1))

    # Level 10.
    face, si, ti = cell.parent(s2.S2CellId.kMaxLevel - 20).GetCenterSiTi()
    self.assertEqual(3, face)
    self.assertEqual(1 << 20, si & ((1 << 21) - 1))
    self.assertEqual(1 << 20, ti & ((1 << 21) - 1))

    # Level 0.
    face, si, ti = cell.parent(0).GetCenterSiTi()
    self.assertEqual(3, face)
    self.assertEqual(1 << 30, si & ((1 << 31) - 1))
    self.assertEqual(1 << 30, ti & ((1 << 31) - 1))

  def testS2CellIdToFromFaceIJ(self):
    cell = s2.S2CellId.FromFaceIJ(3, 1234, 5678)
    face, i, j, _ = cell.ToFaceIJOrientation()
    self.assertEqual(3, face)
    self.assertEqual(1234, i)
    self.assertEqual(5678, j)

  def testS2EarthMetricRadians(self):
    radius_rad = s2.S2Earth.KmToRadians(12.34)
    self.assertAlmostEqual(radius_rad, 0.0019368985451286374)
    angle = s2.S1Angle.Radians(radius_rad)
    radius_m = s2.S2Earth.RadiansToMeters(angle.radians())
    self.assertEqual(radius_m, 12340.0)

  def testS2Point_ToFromRaw(self):
    p = s2.S2Point_FromRaw(1.0, 1.0, 1.0)
    p = p.Normalize()
    p_raw = s2.S2Point_ToRaw(p)
    self.assertAlmostEqual(p_raw[0], 0.57735027)
    self.assertAlmostEqual(p_raw[1], 0.57735027)
    self.assertAlmostEqual(p_raw[2], 0.57735027)

  def testS2Interpolate(self):
    p1 = s2.S2LatLng.FromDegrees(3.0, 4.0).ToPoint()
    p2 = s2.S2LatLng.FromDegrees(4.0, 5.0).ToPoint()

    p3 = s2.Interpolate(p1, p2, 0.5)
    self.assertEqual("3.500133,4.499733", s2.S2LatLng(p3).ToStringInDegrees())

  def testS2predOrderedCCW(self):
    pC = s2.S2LatLng.FromDegrees(3.0, 4.0).ToPoint()
    p1 = s2.S2LatLng.FromDegrees(3.0, 5.0).ToPoint()
    p2 = s2.S2LatLng.FromDegrees(2.7, 3.0).ToPoint()
    p3 = s2.S2LatLng.FromDegrees(3.3, 3.0).ToPoint()

    self.assertFalse(s2.OrderedCCW(p1, p2, p3, pC))
    self.assertTrue(s2.OrderedCCW(p2, p1, p3, pC))
    self.assertTrue(s2.OrderedCCW(p1, p3, p2, pC))
    self.assertTrue(s2.OrderedCCW(p3, p2, p1, pC))

  def testS2UpdateMinDistance(self):
    pC = s2.S2LatLng.FromDegrees(3.0, 4.0).ToPoint()
    p1 = s2.S2LatLng.FromDegrees(3.0, 5.0).ToPoint()
    p2 = s2.S2LatLng.FromDegrees(2.7, 3.0).ToPoint()

    md = s2.S1ChordAngle.Infinity()
    self.assertTrue(s2.UpdateMinDistance(pC, p1, p2, md))
    self.assertAlmostEqual(0.1478883, md.degrees())

    pC = s2.S2LatLng.FromDegrees(3.0, 40.0).ToPoint()
    self.assertFalse(s2.UpdateMinDistance(pC, p1, p2, md))
    self.assertAlmostEqual(0.1478883, md.degrees())

  def testS1AngleToString(self):
    a = s2.S1Angle.Degrees(123)
    self.assertEqual("123.0000000", str(a))

  def testS2CellIdToString(self):
    a = s2.S2CellId(12345)
    self.assertEqual("0/000000000000000000000001200130", str(a))

class S1AngleTest(unittest.TestCase):
  def testE5(self):
    a = s2.S1Angle.E5(1000000)
    self.assertAlmostEqual(10.0, a.degrees())

  def testE6(self):
    a = s2.S1Angle.E6(10000000)
    self.assertAlmostEqual(10.0, a.degrees())

  def testE7(self):
    a = s2.S1Angle.E7(100000000)
    self.assertAlmostEqual(10.0, a.degrees())

  def testAccessors_e6_e7(self):
    a = s2.S1Angle.Degrees(45.0)
    self.assertEqual(45000000, a.e6())
    self.assertEqual(450000000, a.e7())

  def testUnsignedE6(self):
    a = s2.S1Angle.UnsignedE6(10000000)
    self.assertAlmostEqual(10.0, a.degrees())

  def testNormalize(self):
    a = s2.S1Angle.Degrees(360 + 45)
    a.Normalize()
    self.assertAlmostEqual(45.0, a.degrees())

  def testNormalized(self):
    a = s2.S1Angle.Degrees(360 + 45)
    b = a.Normalized()
    self.assertAlmostEqual(45.0, b.degrees())

  def testAbs(self):
    a = s2.S1Angle.Degrees(-45.0)
    b = a.abs()
    self.assertAlmostEqual(45.0, b.degrees())

  def testRadiansAccessor(self):
    a = s2.S1Angle.Degrees(180.0)
    self.assertAlmostEqual(3.14159265, a.radians(), places=5)

class S2CapTest(unittest.TestCase):
  def testEmpty(self):
    cap = s2.S2Cap.Empty()
    self.assertTrue(cap.is_empty())

  def testFull(self):
    cap = s2.S2Cap.Full()
    self.assertTrue(cap.is_valid())
    self.assertFalse(cap.is_empty())

  def testFromPoint(self):
    p = s2.S2LatLng.FromDegrees(1.0, 2.0).ToPoint()
    cap = s2.S2Cap.FromPoint(p)
    self.assertTrue(cap.is_valid())
    self.assertTrue(cap.Contains(p))

  def testFromCenterHeight(self):
    p = s2.S2LatLng.FromDegrees(1.0, 2.0).ToPoint()
    cap = s2.S2Cap.FromCenterHeight(p, 0.5)
    self.assertTrue(cap.is_valid())
    self.assertAlmostEqual(0.5, cap.height())

  def testCenterAndHeight(self):
    center = s2.S2LatLng.FromDegrees(2.0, 3.0).ToPoint()
    cap = s2.S2Cap(center, s2.S1Angle.Degrees(1.0))
    c = cap.center()
    self.assertIsInstance(c, s2.S2Point)
    self.assertGreater(cap.height(), 0.0)

  def testGetCentroid(self):
    center = s2.S2LatLng.FromDegrees(2.0, 3.0).ToPoint()
    cap = s2.S2Cap(center, s2.S1Angle.Degrees(1.0))
    centroid = cap.GetCentroid()
    self.assertIsInstance(centroid, s2.S2Point)

  def testExpanded(self):
    center = s2.S2LatLng.FromDegrees(2.0, 3.0).ToPoint()
    cap = s2.S2Cap(center, s2.S1Angle.Degrees(1.0))
    expanded = cap.Expanded(s2.S1Angle.Degrees(1.0))
    self.assertGreater(expanded.height(), cap.height())

  def testUnion(self):
    p1 = s2.S2LatLng.FromDegrees(1.0, 1.0).ToPoint()
    p2 = s2.S2LatLng.FromDegrees(10.0, 10.0).ToPoint()
    cap1 = s2.S2Cap(p1, s2.S1Angle.Degrees(1.0))
    cap2 = s2.S2Cap(p2, s2.S1Angle.Degrees(1.0))
    union = cap1.Union(cap2)
    self.assertTrue(union.Contains(p1))
    self.assertTrue(union.Contains(p2))

  def testIntersects(self):
    p1 = s2.S2LatLng.FromDegrees(1.0, 1.0).ToPoint()
    cap1 = s2.S2Cap(p1, s2.S1Angle.Degrees(5.0))
    cap2 = s2.S2Cap(p1, s2.S1Angle.Degrees(1.0))
    self.assertTrue(cap1.Intersects(cap2))

  def testEncodeDecode(self):
    center = s2.S2LatLng.FromDegrees(2.0, 3.0).ToPoint()
    cap = s2.S2Cap(center, s2.S1Angle.Degrees(1.0))
    encoder = s2.Encoder()
    cap.Encode(encoder)
    decoder = s2.Decoder(encoder.buffer())
    cap2 = s2.S2Cap()
    self.assertTrue(cap2.Decode(decoder))
    # S2Cap.ApproxEquals is not exposed via SWIG, so compare heights.
    self.assertAlmostEqual(cap.height(), cap2.height())

  def testAddPoint(self):
    cap = s2.S2Cap.Empty()
    p = s2.S2LatLng.FromDegrees(1.0, 2.0).ToPoint()
    cap.AddPoint(p)
    self.assertTrue(cap.Contains(p))

class S2CellTest(unittest.TestCase):
  def testApproxArea(self):
    cell = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(0, 0)).parent(10))
    self.assertGreater(cell.ApproxArea(), 0.0)

  def testExactArea(self):
    cell = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(0, 0)).parent(10))
    self.assertGreater(cell.ExactArea(), 0.0)

  def testAverageArea(self):
    area = s2.S2Cell.AverageArea(10)
    self.assertGreater(area, 0.0)

  def testGetCenter(self):
    cell = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(0, 0)).parent(10))
    center = cell.GetCenter()
    self.assertIsInstance(center, s2.S2Point)

  def testGetVertex(self):
    cell = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(0, 0)).parent(10))
    v = cell.GetVertex(0)
    self.assertIsInstance(v, s2.S2Point)

  def testGetS2LatLngEdge(self):
    cell = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(0, 0)).parent(10))
    edge = cell.GetS2LatLngEdge(0)
    self.assertIsInstance(edge, s2.S2LatLng)

  def testFaceAndLevel(self):
    cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(0, 0)).parent(10)
    cell = s2.S2Cell(cell_id)
    self.assertEqual(10, cell.level())
    self.assertIn(cell.face(), range(6))

  def testIdRoundTrip(self):
    cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(0, 0)).parent(10)
    cell = s2.S2Cell(cell_id)
    self.assertEqual(cell_id, cell.id())

  def testEncodeDecode(self):
    cell = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(0, 0)).parent(10))
    encoder = s2.Encoder()
    cell.Encode(encoder)
    decoder = s2.Decoder(encoder.buffer())
    cell2 = s2.S2Cell()
    self.assertTrue(cell2.Decode(decoder))
    self.assertEqual(cell.id(), cell2.id())

class S2CellIdTest(unittest.TestCase):
  def testBeginEnd(self):
    begin = s2.S2CellId.Begin(0)
    end = s2.S2CellId.End(0)
    self.assertTrue(begin.is_valid())
    self.assertLess(begin, end)

  def testChildBeginEnd(self):
    cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(1.0, 2.0)).parent(10)
    begin = cell_id.child_begin()
    end = cell_id.child_end()
    self.assertLess(begin, end)

  def testIsLeaf(self):
    leaf = s2.S2CellId(s2.S2LatLng.FromDegrees(1.0, 2.0))
    self.assertTrue(leaf.is_leaf())
    non_leaf = leaf.parent(10)
    self.assertFalse(non_leaf.is_leaf())

  def testIsFace(self):
    face = s2.S2CellId(s2.S2LatLng.FromDegrees(1.0, 2.0)).parent(0)
    self.assertTrue(face.is_face())
    non_face = s2.S2CellId(s2.S2LatLng.FromDegrees(1.0, 2.0)).parent(10)
    self.assertFalse(non_face.is_face())

  def testPos(self):
    cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(1.0, 2.0))
    pos = cell_id.pos()
    self.assertIsInstance(pos, int)

  def testRangeMinMax(self):
    cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(1.0, 2.0)).parent(10)
    rmin = cell_id.range_min()
    rmax = cell_id.range_max()
    self.assertLess(rmin, rmax)

  def testNextPrev(self):
    cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(1.0, 2.0)).parent(10)
    nxt = cell_id.next()
    prv = cell_id.prev()
    self.assertNotEqual(cell_id, nxt)
    self.assertNotEqual(cell_id, prv)
    self.assertEqual(cell_id, nxt.prev())
    self.assertEqual(cell_id, prv.next())

  def testKMaxLevel(self):
    self.assertEqual(30, s2.S2CellId.kMaxLevel)

class S2CellUnionTest(unittest.TestCase):
  def testApproxArea(self):
    cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(1.0, 2.0)).parent(10)
    cu = s2.S2CellUnion()
    cu.Init([cell_id.id()])
    self.assertGreater(cu.ApproxArea(), 0.0)

  def testDenormalize(self):
    cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(1.0, 2.0)).parent(10)
    cu = s2.S2CellUnion()
    cu.Init([cell_id.id()])
    denormed = cu.Denormalize(0, 1)
    self.assertGreater(len(denormed), 0)

  def testCellId(self):
    cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(1.0, 2.0)).parent(10)
    cu = s2.S2CellUnion()
    cu.Init([cell_id.id()])
    self.assertEqual(cell_id, cu.cell_id(0))

  def testEncodeDecode(self):
    cell_id = s2.S2CellId(s2.S2LatLng.FromDegrees(1.0, 2.0)).parent(10)
    cu = s2.S2CellUnion()
    cu.Init([cell_id.id()])
    encoder = s2.Encoder()
    cu.Encode(encoder)
    decoder = s2.Decoder(encoder.buffer())
    cu2 = s2.S2CellUnion()
    self.assertTrue(cu2.Decode(decoder))
    self.assertEqual(cu.num_cells(), cu2.num_cells())

class S2LoopTest(unittest.TestCase):
  def _makeLoop(self):
    cell = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8))
    return s2.S2Loop(cell)

  def testClone(self):
    loop = self._makeLoop()
    clone = loop.Clone()
    self.assertTrue(loop.Equals(clone))

  def testEncodeDecode(self):
    loop = self._makeLoop()
    encoder = s2.Encoder()
    loop.Encode(encoder)
    decoder = s2.Decoder(encoder.buffer())
    loop2 = s2.S2Loop()
    self.assertTrue(loop2.Decode(decoder))
    self.assertTrue(loop.Equals(loop2))

  def testGetCentroid(self):
    loop = self._makeLoop()
    centroid = loop.GetCentroid()
    self.assertIsInstance(centroid, s2.S2Point)

  def testGetDistance(self):
    loop = self._makeLoop()
    p = s2.S2LatLng.FromDegrees(3.0, 4.0).ToPoint()
    d = loop.GetDistance(p)
    self.assertIsInstance(d, s2.S1Angle)

  def testIntersects(self):
    loop1 = self._makeLoop()
    loop2 = self._makeLoop()
    self.assertTrue(loop1.Intersects(loop2))

  def testProject(self):
    loop = self._makeLoop()
    p = s2.S2LatLng.FromDegrees(3.0, 4.0).ToPoint()
    projected = loop.Project(p)
    self.assertIsInstance(projected, s2.S2Point)

  def testIsEmpty(self):
    loop = self._makeLoop()
    self.assertFalse(loop.is_empty())

  def testSign(self):
    loop = self._makeLoop()
    sign = loop.sign()
    self.assertEqual(1, sign)

  def testVertex(self):
    loop = self._makeLoop()
    v = loop.vertex(0)
    self.assertIsInstance(v, s2.S2Point)

class S2PolylineTest(unittest.TestCase):
  def _makePolyline(self):
    line = s2.S2Polyline()
    line.InitFromS2LatLngs([s2.S2LatLng.FromDegrees(0.0, 0.0),
                            s2.S2LatLng.FromDegrees(1.0, 1.0)])
    return line

  def testClone(self):
    line = self._makePolyline()
    clone = line.Clone()
    self.assertTrue(line.ApproxEquals(clone))

  def testEncodeDecode(self):
    line = self._makePolyline()
    encoder = s2.Encoder()
    line.Encode(encoder)
    decoder = s2.Decoder(encoder.buffer())
    line2 = s2.S2Polyline()
    self.assertTrue(line2.Decode(decoder))
    self.assertTrue(line.ApproxEquals(line2))

  def testApproxEquals(self):
    line1 = self._makePolyline()
    line2 = self._makePolyline()
    self.assertTrue(line1.ApproxEquals(line2))

  def testGetSuffix(self):
    line = self._makePolyline()
    point, next_vertex = line.GetSuffix(0.5)
    self.assertIsInstance(point, s2.S2Point)
    self.assertIsInstance(next_vertex, int)

  def testInterpolateAndUnInterpolate(self):
    line = self._makePolyline()
    p = line.Interpolate(0.5)
    self.assertIsInstance(p, s2.S2Point)
    f = line.UnInterpolate(p, 1)
    self.assertAlmostEqual(0.5, f, places=5)

  def testIntersects(self):
    line1 = s2.S2Polyline()
    line1.InitFromS2LatLngs([s2.S2LatLng.FromDegrees(-1.0, 0.0),
                             s2.S2LatLng.FromDegrees(1.0, 0.0)])
    line2 = s2.S2Polyline()
    line2.InitFromS2LatLngs([s2.S2LatLng.FromDegrees(0.0, -1.0),
                             s2.S2LatLng.FromDegrees(0.0, 1.0)])
    self.assertTrue(line1.Intersects(line2))

  def testIsOnRight(self):
    line = self._makePolyline()
    p = s2.S2LatLng.FromDegrees(0.0, 1.0).ToPoint()
    result = line.IsOnRight(p)
    self.assertIsInstance(result, bool)

  def testIsValid(self):
    line = self._makePolyline()
    self.assertTrue(line.IsValid())

  def testProject(self):
    line = self._makePolyline()
    p = s2.S2LatLng.FromDegrees(0.5, 0.5).ToPoint()
    projected = line.Project(p)
    self.assertIsInstance(projected[0], s2.S2Point)

  def testReverse(self):
    line = self._makePolyline()
    v0_before = s2.S2LatLng(line.vertex(0)).ToStringInDegrees()
    line.Reverse()
    v_last = s2.S2LatLng(line.vertex(0)).ToStringInDegrees()
    self.assertNotEqual(v0_before, v_last)

  def testNumVertices(self):
    line = self._makePolyline()
    self.assertEqual(2, line.num_vertices())

  def testVertex(self):
    line = self._makePolyline()
    v = line.vertex(0)
    self.assertIsInstance(v, s2.S2Point)

class S2EarthTest(unittest.TestCase):
  def testGetDistanceLatLng(self):
    a = s2.S2LatLng.FromDegrees(0.0, 0.0)
    b = s2.S2LatLng.FromDegrees(1.0, 0.0)
    # GetDistance returns an opaque SWIG Meters* object.
    d = s2.S2Earth.GetDistance(a, b)
    self.assertTrue(d is not None)

  def testGetDistanceS2Point(self):
    a = s2.S2LatLng.FromDegrees(0.0, 0.0).ToPoint()
    b = s2.S2LatLng.FromDegrees(1.0, 0.0).ToPoint()
    # GetDistance returns an opaque SWIG Meters* object.
    d = s2.S2Earth.GetDistance(a, b)
    self.assertTrue(d is not None)

  def testGetDistanceKmLatLng(self):
    a = s2.S2LatLng.FromDegrees(0.0, 0.0)
    b = s2.S2LatLng.FromDegrees(1.0, 0.0)
    km = s2.S2Earth.GetDistanceKm(a, b)
    self.assertGreater(km, 100.0)

  def testGetDistanceKmS2Point(self):
    a = s2.S2LatLng.FromDegrees(0.0, 0.0).ToPoint()
    b = s2.S2LatLng.FromDegrees(1.0, 0.0).ToPoint()
    km = s2.S2Earth.GetDistanceKm(a, b)
    self.assertGreater(km, 100.0)

  def testGetDistanceMetersLatLng(self):
    a = s2.S2LatLng.FromDegrees(0.0, 0.0)
    b = s2.S2LatLng.FromDegrees(1.0, 0.0)
    m = s2.S2Earth.GetDistanceMeters(a, b)
    self.assertGreater(m, 100000.0)

  def testGetDistanceMetersS2Point(self):
    a = s2.S2LatLng.FromDegrees(0.0, 0.0).ToPoint()
    b = s2.S2LatLng.FromDegrees(1.0, 0.0).ToPoint()
    m = s2.S2Earth.GetDistanceMeters(a, b)
    self.assertGreater(m, 100000.0)

  def testGetInitialBearing(self):
    a = s2.S2LatLng.FromDegrees(0.0, 0.0)
    b = s2.S2LatLng.FromDegrees(1.0, 0.0)
    bearing = s2.S2Earth.GetInitialBearing(a, b)
    self.assertIsInstance(bearing, s2.S1Angle)

  def testRadiusMethods(self):
    km = s2.S2Earth.RadiusKm()
    self.assertAlmostEqual(6371.01, km, places=1)
    m = s2.S2Earth.RadiusMeters()
    self.assertAlmostEqual(6371010.0, m, places=-1)

  def testAltitudeMethods(self):
    hm = s2.S2Earth.HighestAltitudeMeters()
    self.assertGreater(hm, 0.0)
    lom = s2.S2Earth.LowestAltitudeMeters()
    self.assertLess(lom, 0.0)

  def testMetersToRadians(self):
    r = s2.S2Earth.MetersToRadians(1000.0)
    self.assertGreater(r, 0.0)

  def testRadiansToKm(self):
    km = s2.S2Earth.RadiansToKm(0.001)
    self.assertGreater(km, 0.0)

  def testSquareConversions(self):
    sr = s2.S2Earth.SquareKmToSteradians(1.0)
    self.assertGreater(sr, 0.0)
    sr2 = s2.S2Earth.SquareMetersToSteradians(1000000.0)
    self.assertGreater(sr2, 0.0)
    km2 = s2.S2Earth.SteradiansToSquareKm(sr)
    self.assertAlmostEqual(1.0, km2, places=5)
    m2 = s2.S2Earth.SteradiansToSquareMeters(sr2)
    self.assertAlmostEqual(1000000.0, m2, places=0)

  def testToDistance(self):
    a = s2.S1Angle.Degrees(1.0)
    # ToDistance returns an opaque SWIG Meters* object.
    d = s2.S2Earth.ToDistance(a)
    self.assertTrue(d is not None)

  def testToDistanceChordAngle(self):
    ca = s2.S1ChordAngle(s2.S1Angle.Degrees(1.0))
    # ToDistance returns an opaque SWIG Meters* object.
    d = s2.S2Earth.ToDistance(ca)
    self.assertTrue(d is not None)

  def testToKm(self):
    a = s2.S1Angle.Degrees(1.0)
    km = s2.S2Earth.ToKm(a)
    self.assertGreater(km, 100.0)

  def testToKmChordAngle(self):
    ca = s2.S1ChordAngle(s2.S1Angle.Degrees(1.0))
    km = s2.S2Earth.ToKm(ca)
    self.assertGreater(km, 100.0)

  def testToMeters(self):
    a = s2.S1Angle.Degrees(1.0)
    m = s2.S2Earth.ToMeters(a)
    self.assertGreater(m, 100000.0)

  def testToMetersChordAngle(self):
    ca = s2.S1ChordAngle(s2.S1Angle.Degrees(1.0))
    m = s2.S2Earth.ToMeters(ca)
    self.assertGreater(m, 100000.0)

class S2LatLngTest(unittest.TestCase):
  def testFromE6(self):
    ll = s2.S2LatLng.FromE6(10000000, 20000000)
    self.assertAlmostEqual(10.0, ll.lat().degrees())
    self.assertAlmostEqual(20.0, ll.lng().degrees())

  def testFromE7(self):
    ll = s2.S2LatLng.FromE7(100000000, 200000000)
    self.assertAlmostEqual(10.0, ll.lat().degrees())
    self.assertAlmostEqual(20.0, ll.lng().degrees())

  def testFromRadians(self):
    ll = s2.S2LatLng.FromRadians(0.0, math.pi / 2)
    self.assertAlmostEqual(0.0, ll.lat().degrees())
    self.assertAlmostEqual(90.0, ll.lng().degrees())

  def testFromUnsignedE6(self):
    ll = s2.S2LatLng.FromUnsignedE6(10000000, 20000000)
    self.assertAlmostEqual(10.0, ll.lat().degrees())
    self.assertAlmostEqual(20.0, ll.lng().degrees())

  def testFromUnsignedE7(self):
    ll = s2.S2LatLng.FromUnsignedE7(100000000, 200000000)
    self.assertAlmostEqual(10.0, ll.lat().degrees())
    self.assertAlmostEqual(20.0, ll.lng().degrees())

  def testCoords(self):
    ll = s2.S2LatLng.FromDegrees(10.0, 20.0)
    # coords() returns an opaque SWIG R2Point object.
    c = ll.coords()
    self.assertTrue(c is not None)

class S2LatLngRectTest(unittest.TestCase):
  def testFromPoint(self):
    r = s2.S2LatLngRect.FromPoint(s2.S2LatLng.FromDegrees(1.0, 2.0))
    self.assertTrue(r.is_point())

  def testFromPointPair(self):
    r = s2.S2LatLngRect.FromPointPair(
        s2.S2LatLng.FromDegrees(1.0, 2.0),
        s2.S2LatLng.FromDegrees(3.0, 4.0))
    self.assertFalse(r.is_empty())

  def testFromCenterSize(self):
    center = s2.S2LatLng.FromDegrees(1.0, 2.0)
    size = s2.S2LatLng.FromDegrees(2.0, 4.0)
    r = s2.S2LatLngRect.FromCenterSize(center, size)
    self.assertFalse(r.is_empty())

  def testEmpty(self):
    r = s2.S2LatLngRect.Empty()
    self.assertTrue(r.is_empty())

  def testFull(self):
    r = s2.S2LatLngRect.Full()
    self.assertFalse(r.is_empty())

  def testExpandedByDistance(self):
    r = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(1.0, 2.0),
                        s2.S2LatLng.FromDegrees(3.0, 4.0))
    expanded = r.ExpandedByDistance(s2.S1Angle.Degrees(1.0))
    self.assertGreater(expanded.Area(), r.Area())

  def testArea(self):
    r = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(0.0, 0.0),
                        s2.S2LatLng.FromDegrees(1.0, 1.0))
    self.assertGreater(r.Area(), 0.0)

  def testGetCenter(self):
    r = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(0.0, 0.0),
                        s2.S2LatLng.FromDegrees(2.0, 4.0))
    c = r.GetCenter()
    self.assertAlmostEqual(1.0, c.lat().degrees())
    self.assertAlmostEqual(2.0, c.lng().degrees())

  def testGetCentroid(self):
    r = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(0.0, 0.0),
                        s2.S2LatLng.FromDegrees(2.0, 4.0))
    centroid = r.GetCentroid()
    self.assertIsInstance(centroid, s2.S2Point)

  def testGetDistance(self):
    r = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(0.0, 0.0),
                        s2.S2LatLng.FromDegrees(1.0, 1.0))
    p = s2.S2LatLng.FromDegrees(5.0, 5.0)
    d = r.GetDistance(p)
    self.assertIsInstance(d, s2.S1Angle)

  def testGetSize(self):
    r = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(0.0, 0.0),
                        s2.S2LatLng.FromDegrees(2.0, 4.0))
    size = r.GetSize()
    self.assertAlmostEqual(2.0, size.lat().degrees())
    self.assertAlmostEqual(4.0, size.lng().degrees())

  def testGetVertex(self):
    r = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(0.0, 0.0),
                        s2.S2LatLng.FromDegrees(2.0, 4.0))
    for i in range(4):
      v = r.GetVertex(i)
      self.assertIsInstance(v, s2.S2LatLng)

  def testIntersectionAndUnion(self):
    r1 = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(0.0, 0.0),
                         s2.S2LatLng.FromDegrees(2.0, 2.0))
    r2 = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(1.0, 1.0),
                         s2.S2LatLng.FromDegrees(3.0, 3.0))
    inter = r1.Intersection(r2)
    self.assertFalse(inter.is_empty())
    union = r1.Union(r2)
    self.assertGreater(union.Area(), r1.Area())

  def testLatLngAccessors(self):
    r = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(1.0, 2.0),
                        s2.S2LatLng.FromDegrees(3.0, 4.0))
    self.assertAlmostEqual(1.0, r.lat_lo().degrees())
    self.assertAlmostEqual(3.0, r.lat_hi().degrees())
    self.assertAlmostEqual(2.0, r.lng_lo().degrees())
    self.assertAlmostEqual(4.0, r.lng_hi().degrees())
    # lat()/lng() return opaque SWIG R1Interval/S1Interval objects.
    self.assertTrue(r.lat() is not None)
    self.assertTrue(r.lng() is not None)
    self.assertIsInstance(r.lo(), s2.S2LatLng)
    self.assertIsInstance(r.hi(), s2.S2LatLng)

  def testEncodeDecode(self):
    r = s2.S2LatLngRect(s2.S2LatLng.FromDegrees(1.0, 2.0),
                        s2.S2LatLng.FromDegrees(3.0, 4.0))
    encoder = s2.Encoder()
    r.Encode(encoder)
    decoder = s2.Decoder(encoder.buffer())
    r2 = s2.S2LatLngRect()
    self.assertTrue(r2.Decode(decoder))
    self.assertTrue(r.ApproxEquals(r2))

  def testAddPoint(self):
    r = s2.S2LatLngRect.Empty()
    r.AddPoint(s2.S2LatLng.FromDegrees(1.0, 2.0))
    self.assertFalse(r.is_empty())

class EncoderDecoderTest(unittest.TestCase):
  def testPutGet8(self):
    enc = s2.Encoder()
    enc.Ensure(1)
    enc.put8(42)
    dec = s2.Decoder(enc.buffer())
    self.assertEqual(42, dec.get8())

  def testPutGet16(self):
    enc = s2.Encoder()
    enc.Ensure(2)
    enc.put16(1234)
    dec = s2.Decoder(enc.buffer())
    self.assertEqual(1234, dec.get16())

  def testPutGet32(self):
    enc = s2.Encoder()
    enc.Ensure(4)
    enc.put32(12345678)
    dec = s2.Decoder(enc.buffer())
    self.assertEqual(12345678, dec.get32())

  def testPutGet64(self):
    enc = s2.Encoder()
    enc.Ensure(8)
    enc.put64(123456789012345)
    dec = s2.Decoder(enc.buffer())
    self.assertEqual(123456789012345, dec.get64())

  def testPutGetFloat(self):
    enc = s2.Encoder()
    enc.Ensure(4)
    enc.putfloat(3.14)
    dec = s2.Decoder(enc.buffer())
    self.assertAlmostEqual(3.14, dec.getfloat(), places=5)

  def testPutGetDouble(self):
    enc = s2.Encoder()
    enc.Ensure(8)
    enc.putdouble(3.141592653589793)
    dec = s2.Decoder(enc.buffer())
    self.assertAlmostEqual(3.141592653589793, dec.getdouble())

  def testEncoderLengthAndAvail(self):
    enc = s2.Encoder()
    self.assertEqual(0, enc.length())
    enc.Ensure(1)
    enc.put8(1)
    self.assertEqual(1, enc.length())
    self.assertGreaterEqual(enc.avail(), 0)

  def testEncoderClear(self):
    enc = s2.Encoder()
    enc.Ensure(1)
    enc.put8(1)
    enc.clear()
    self.assertEqual(0, enc.length())

  def testEncoderEnsure(self):
    enc = s2.Encoder()
    enc.Ensure(100)
    self.assertGreaterEqual(enc.avail(), 100)

  def testDecoderPosAndAvail(self):
    enc = s2.Encoder()
    enc.Ensure(2)
    enc.put8(1)
    enc.put8(2)
    dec = s2.Decoder(enc.buffer())
    self.assertEqual(0, dec.pos())
    self.assertEqual(2, dec.avail())
    dec.get8()
    self.assertEqual(1, dec.pos())
    self.assertEqual(1, dec.avail())

class RegionTermIndexerTest(unittest.TestCase):
  def _randomCaps(self, query_type, **indexer_options):
    # This function creates an index consisting either of points (if
    # options.index_contains_points_only() is true) or S2Caps of random size.
    # It then executes queries consisting of points (if query_type == POINT)
    # or S2Caps of random size (if query_type == CAP).
    #
    # indexer_options are set on both the indexer & coverer (if relevant)
    # eg. _randomCaps('cap', min_level=0) calls indexer.set_min_level(0)
    ITERATIONS = 400

    indexer = s2.S2RegionTermIndexer()
    coverer = s2.S2RegionCoverer()

    # set indexer options
    for opt_key, opt_value in indexer_options.items():
      setter = "set_%s" % opt_key
      getattr(indexer, setter)(opt_value)
      if hasattr(coverer, setter):
        getattr(coverer, setter)(opt_value)

    caps = []
    coverings = []
    index = defaultdict(set)

    index_terms = 0
    query_terms = 0
    for i in range(ITERATIONS):
      # Choose the region to be indexed: either a single point or a cap
      # of random size (up to a full sphere).
      terms = []
      if indexer.index_contains_points_only():
        cap = s2.S2Cap.FromPoint(s2.S2Testing.RandomPoint())
        terms = indexer.GetIndexTerms(cap.center(), "")
      else:
        cap = s2.S2Testing.GetRandomCap(
          0.3 * s2.S2Cell.AverageArea(indexer.max_level()),
          4.0 * s2.S2Cell.AverageArea(indexer.min_level())
        )
        terms = indexer.GetIndexTerms(cap, "")

      caps.append(cap)
      coverings.append(s2.S2CellUnion(coverer.GetCovering(cap)))
      for term in terms:
        index[term].add(i)

      index_terms += len(terms)

    for i in range(ITERATIONS):
      # Choose the region to be queried: either a random point or a cap of
      # random size.
      terms = []

      if query_type == 'cap':
        cap = s2.S2Cap.FromPoint(s2.S2Testing.RandomPoint())
        terms = indexer.GetQueryTerms(cap.center(), "")
      else:
        cap = s2.S2Testing.GetRandomCap(
          0.3 * s2.S2Cell.AverageArea(indexer.max_level()),
          4.0 * s2.S2Cell.AverageArea(indexer.min_level())
        )
        terms = indexer.GetQueryTerms(cap, "")

      # Compute the expected results of the S2Cell query by brute force.
      covering = s2.S2CellUnion(coverer.GetCovering(cap))
      expected, actual = set(), set()
      for j in range(len(caps)):
        if covering.Intersects(coverings[j]):
          expected.add(j)

      for term in terms:
        actual |= index[term]

      self.assertEqual(expected, actual)
      query_terms += len(terms)

      print("Index terms/doc: %0.2f, Query terms/doc: %0.2f" % (
        float(index_terms) / ITERATIONS,
        float(query_terms) / ITERATIONS)
      )

  # We run one test case for each combination of space vs. time optimization,
  # and indexing regions vs. only points.

  def testIndexRegionsQueryRegionsOptimizeTime(self):
    self._randomCaps("cap",
      optimize_for_space=False,
      min_level=0,
      max_level=16,
      max_cells=20,
    )

  def testIndexRegionsQueryPointsOptimizeTime(self):
    self._randomCaps("point",
      optimize_for_space=False,
      min_level=0,
      max_level=16,
      max_cells=20,
    )

  def testIndexRegionsQueryRegionsOptimizeTimeWithLevelMod(self):
    self._randomCaps("cap",
      optimize_for_space=False,
      min_level=6,
      max_level=12,
      level_mod=3,
    )

  def testIndexRegionsQueryRegionsOptimizeSpace(self):
    self._randomCaps("cap",
      optimize_for_space=True,
      min_level=4,
      max_level=s2.S2CellId.kMaxLevel,
      max_cells=8,
    )

  def testIndexPointsQueryRegionsOptimizeTime(self):
    self._randomCaps("cap",
      optimize_for_space=False,
      min_level=0,
      max_level=s2.S2CellId.kMaxLevel,
      level_mod=2,
      max_cells=20,
      index_contains_points_only=True,
    )

  def testIndexPointsQueryRegionsOptimizeSpace(self):
    self._randomCaps("cap",
      optimize_for_space=True,
      index_contains_points_only=True,
    )

  def testMaxLevelSetLoosely(self):
    # Test that correct terms are generated even when (max_level - min_level)
    # is not a multiple of level_mod.
    indexer1 = s2.S2RegionTermIndexer()
    indexer1.set_min_level(1)
    indexer1.set_level_mod(2)
    indexer1.set_max_level(19)

    indexer2 = s2.S2RegionTermIndexer()
    indexer2.set_min_level(1)
    indexer2.set_level_mod(2)
    indexer2.set_max_level(19)
    indexer2.set_max_level(20)

    point = s2.S2Testing.RandomPoint()

    self.assertEqual(
      indexer1.GetIndexTerms(point, ""),
      indexer2.GetIndexTerms(point, "")
    )
    self.assertEqual(
      indexer1.GetQueryTerms(point, ""),
      indexer2.GetQueryTerms(point, "")
    )

    cap = s2.S2Testing.GetRandomCap(0.0, 1.0)
    self.assertEqual(
      indexer1.GetIndexTerms(cap, ""),
      indexer2.GetIndexTerms(cap, "")
    )
    self.assertEqual(
      indexer1.GetQueryTerms(cap, ""),
      indexer2.GetQueryTerms(cap, "")
    )

  def testS2CellIdFromDebugString(self):
    cell = s2.S2CellId.FromDebugString("5/31200")
    self.assertTrue(cell.is_valid())
    self.assertEqual("5/31200", cell.ToString())

  def testCovererSetFixedLevel(self):
    coverer = s2.S2RegionCoverer()
    coverer.set_fixed_level(10)
    self.assertEqual(10, coverer.min_level())
    self.assertEqual(10, coverer.max_level())

  def testCovererTrueMaxLevel(self):
    coverer = s2.S2RegionCoverer()
    coverer.set_max_level(20)
    self.assertEqual(20, coverer.true_max_level())

  def testIndexerMarkerCharacter(self):
    indexer = s2.S2RegionTermIndexer()
    default_marker = indexer.marker_character()
    self.assertIsInstance(default_marker, str)

  def testIndexerTrueMaxLevel(self):
    indexer = s2.S2RegionTermIndexer()
    indexer.set_max_level(20)
    self.assertEqual(20, indexer.true_max_level())

  def testIndexerSetFixedLevel(self):
    indexer = s2.S2RegionTermIndexer()
    indexer.set_fixed_level(10)
    self.assertEqual(10, indexer.min_level())
    self.assertEqual(10, indexer.max_level())

class S2PolygonTestCase(unittest.TestCase):

  def _makePolygon(self):
    cell = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8))
    return s2.S2Polygon(cell)
  def testInitToUnionSame(self):
    cell1 = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8))
    polygon1 = s2.S2Polygon(cell1)

    cell2 = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8))
    polygon2 = s2.S2Polygon(cell2)

    polygon3 = s2.S2Polygon()
    polygon3.InitToUnion(polygon1, polygon2)

    self.assertEqual(1, polygon3.num_loops())
    loop = polygon3.loop(0)
    self.assertEqual(4, loop.num_vertices())

  def testInitToUnionDistinct(self):
    cell1 = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8))
    polygon1 = s2.S2Polygon(cell1)

    cell2 = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(13.0, 4.0)).parent(8))
    polygon2 = s2.S2Polygon(cell2)

    polygon3 = s2.S2Polygon()
    polygon3.InitToUnion(polygon1, polygon2)

    self.assertEqual(2, polygon3.num_loops())
    loop = polygon3.loop(0)
    self.assertEqual(4, loop.num_vertices())
    loop = polygon3.loop(1)
    self.assertEqual(4, loop.num_vertices())

  def testClone(self):
    poly = self._makePolygon()
    clone = poly.Clone()
    self.assertTrue(poly.Equals(clone))

  def testCopy(self):
    poly = self._makePolygon()
    copy = s2.S2Polygon()
    copy.Copy(poly)
    self.assertTrue(poly.Equals(copy))

  def testGetCentroid(self):
    poly = self._makePolygon()
    centroid = poly.GetCentroid()
    self.assertIsInstance(centroid, s2.S2Point)

  def testGetDistance(self):
    poly = self._makePolygon()
    p = s2.S2LatLng.FromDegrees(3.0, 4.0).ToPoint()
    d = poly.GetDistance(p)
    self.assertIsInstance(d, s2.S1Angle)

  def testProject(self):
    poly = self._makePolygon()
    p = s2.S2LatLng.FromDegrees(3.0, 4.0).ToPoint()
    projected = poly.Project(p)
    self.assertIsInstance(projected, s2.S2Point)

  def testIsEmpty(self):
    empty = s2.S2Polygon()
    self.assertTrue(empty.is_empty())
    poly = self._makePolygon()
    self.assertFalse(poly.is_empty())

  def testNumVertices(self):
    poly = self._makePolygon()
    self.assertEqual(4, poly.num_vertices())

  def testIntersects(self):
    poly1 = self._makePolygon()
    poly2 = self._makePolygon()
    self.assertTrue(poly1.Intersects(poly2))

  def testIsValid(self):
    poly = self._makePolygon()
    self.assertTrue(poly.IsValid())

class S2ChordAngleTest(unittest.TestCase):
  def testBasic(self):
    ca = s2.S1ChordAngle(s2.S1Angle.Degrees(100))
    self.assertAlmostEqual(100, ca.degrees())

  def testArithmetic(self):
    ca1 = s2.S1ChordAngle(s2.S1Angle.Degrees(10))
    ca2 = s2.S1ChordAngle(s2.S1Angle.Degrees(20))
    ca3 = ca1 + ca2
    self.assertAlmostEqual(30, ca3.degrees())
    ca4 = ca2 - ca1
    self.assertAlmostEqual(10, ca4.degrees())

  def testComparison(self):
    ca1 = s2.S1ChordAngle(s2.S1Angle.Degrees(10))
    ca2 = s2.S1ChordAngle(s2.S1Angle.Degrees(20))
    self.assertTrue(ca1 < ca2)
    self.assertTrue(ca2 > ca1)
    self.assertFalse(ca1 > ca2)
    self.assertFalse(ca2 < ca1)

    ca3 = s2.S1ChordAngle(s2.S1Angle.Degrees(10))
    self.assertTrue(ca1 == ca3)
    self.assertFalse(ca1 == ca2)
    self.assertFalse(ca1 != ca3)
    self.assertTrue(ca1 != ca2)

  def testInfinity(self):
    ca1 = s2.S1ChordAngle(s2.S1Angle.Degrees(179))
    ca2 = s2.S1ChordAngle.Infinity()
    self.assertTrue(ca2 > ca1)

  def testCopy(self):
    ca1 = s2.S1ChordAngle(s2.S1Angle.Degrees(100))
    ca2 = s2.S1ChordAngle(ca1)
    self.assertAlmostEqual(100, ca2.degrees())

  def testToAngle(self):
    ca = s2.S1ChordAngle(s2.S1Angle.Degrees(90.0))
    a = ca.ToAngle()
    self.assertAlmostEqual(90.0, a.degrees())

class S2BufferOperationTest(unittest.TestCase):
  def setUp(self):
    self.opts = s2.S2BufferOperationOptions()
    self.result = s2.S2Polygon()
    self.layer = s2.S2PolygonLayer(self.result)

  def testDefaults(self):
    op = s2.S2BufferOperation(self.layer, self.opts)

    cell1 = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8))
    op.AddPolygon(s2.S2Polygon(cell1))
    op.Build()

    self.assertEqual(1, self.result.num_loops())
    loop = self.result.loop(0)
    self.assertEqual(4, loop.num_vertices())

  def testRadius(self):
    self.opts.set_buffer_radius(s2.S1Angle.Degrees(0.001))
    op = s2.S2BufferOperation(self.layer, self.opts)

    cell1 = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8))
    op.AddPolygon(s2.S2Polygon(cell1))
    op.Build()

    self.assertEqual(1, self.result.num_loops())
    loop = self.result.loop(0)
    self.assertEqual(20, loop.num_vertices())

  def testRadiusAndError(self):
    self.opts.set_buffer_radius(s2.S1Angle.Degrees(0.001))
    self.opts.set_error_fraction(0.1)
    op = s2.S2BufferOperation(self.layer, self.opts)

    cell1 = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8))
    op.AddPolygon(s2.S2Polygon(cell1))
    op.Build()

    self.assertEqual(1, self.result.num_loops())
    loop = self.result.loop(0)
    self.assertEqual(12, loop.num_vertices())

  def testPoint(self):
    self.opts.set_buffer_radius(s2.S1Angle.Degrees(0.001))
    op = s2.S2BufferOperation(self.layer, self.opts)

    op.AddPoint(s2.S2LatLng.FromDegrees(14.0, 15.0).ToPoint())
    op.Build()

    self.assertEqual(1, self.result.num_loops())
    loop = self.result.loop(0)
    self.assertEqual(16, loop.num_vertices())

class S2BooleanOperationTest(unittest.TestCase):
  def setUp(self):
    self.index1 = s2.MutableS2ShapeIndex()
    self.index2 = s2.MutableS2ShapeIndex()

  def operate(self):
    poly = s2.S2Polygon()
    layer = s2.S2PolygonLayer(poly)
    op = s2.S2BooleanOperation(s2.S2BooleanOperation.OpType_UNION,
                               layer)
    op.Build(self.index1, self.index2)
    return poly

  def testUnionSame(self):
    cell1 = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8))
    self.index1.Add(s2.S2Polygon(cell1))

    cell2 = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8))
    self.index2.Add(s2.S2Polygon(cell2))

    result = self.operate()
    self.assertEqual(1, result.num_loops())
    loop = result.loop(0)
    self.assertEqual(4, loop.num_vertices())

  def testUnionDistinct(self):
    cell1 = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(3.0, 4.0)).parent(8))
    self.index1.Add(s2.S2Polygon(cell1))

    cell2 = s2.S2Cell(s2.S2CellId(s2.S2LatLng.FromDegrees(13.0, 4.0)).parent(8))
    self.index2.Add(s2.S2Polygon(cell2))

    result = self.operate()
    self.assertEqual(2, result.num_loops())
    loop = result.loop(0)
    self.assertEqual(4, loop.num_vertices())
    loop = result.loop(1)
    self.assertEqual(4, loop.num_vertices())

# See https://github.com/google/s2geometry/issues/376
@unittest.skip("CHECK-fails in debug mode.")
class S2BuilderTest(unittest.TestCase):
  def setUp(self):
    self.p1 = s2.S2LatLng.FromDegrees(10.0, 10.0).ToPoint()
    self.p2 = s2.S2LatLng.FromDegrees(10.0, 11.0).ToPoint()
    self.p3 = s2.S2LatLng.FromDegrees(11.0, 10.0).ToPoint()
    self.p4 = s2.S2LatLng.FromDegrees(11.0, 11.0).ToPoint()

  def testSquareDirected(self):
    opts = s2.S2PolygonLayerOptions()
    opts.set_edge_type(s2.S2Builder.EdgeType_DIRECTED)

    result = s2.S2Polygon()
    b = s2.S2Builder()
    b.StartLayer(s2.S2PolygonLayer(result, opts))

    b.AddEdge(self.p1, self.p2)
    b.AddEdge(self.p2, self.p4)
    b.AddEdge(self.p4, self.p3)
    b.AddEdge(self.p3, self.p1)

    b.Build()

    self.assertEqual(1, result.num_loops())
    loop = result.loop(0)
    self.assertEqual(4, loop.num_vertices())

  def testSquareUndirected(self):
    opts = s2.S2PolygonLayerOptions()
    opts.set_edge_type(s2.S2Builder.EdgeType_UNDIRECTED)

    result = s2.S2Polygon()
    b = s2.S2Builder()
    b.StartLayer(s2.S2PolygonLayer(result, opts))

    b.AddEdge(self.p1, self.p2)
    b.AddEdge(self.p1, self.p3)
    b.AddEdge(self.p2, self.p4)
    b.AddEdge(self.p3, self.p4)

    b.Build()

    self.assertEqual(1, result.num_loops())
    loop = result.loop(0)
    self.assertEqual(4, loop.num_vertices())

  def testException(self):
    opts = s2.S2PolygonLayerOptions()
    opts.set_edge_type(s2.S2Builder.EdgeType_DIRECTED)

    result = s2.S2Polygon()
    b = s2.S2Builder()
    b.StartLayer(s2.S2PolygonLayer(result, opts))

    b.AddEdge(self.p1, self.p2)
    b.AddEdge(self.p1, self.p3)
    b.AddEdge(self.p2, self.p4)
    b.AddEdge(self.p3, self.p4)

    with self.assertRaises(ValueError):
      b.Build()

if __name__ == "__main__":
  unittest.main()
