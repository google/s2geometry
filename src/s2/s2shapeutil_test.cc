// Copyright 2013 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS-IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

// Author: ericv@google.com (Eric Veach)

#include "s2/s2shapeutil.h"

#include <memory>
#include <vector>

#include <gtest/gtest.h>
#include "s2/third_party/absl/memory/memory.h"
#include "s2/s2cap.h"
#include "s2/s2loop.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"
#include "s2/s2testing.h"
#include "s2/s2textformat.h"

using std::unique_ptr;
using std::vector;

namespace s2shapeutil {

TEST(LaxLoop, EmptyLoop) {
  // Test S2Loop constructor.
  LaxLoop shape;
  shape.Init(S2Loop(S2Loop::kEmpty()));
  EXPECT_EQ(0, shape.num_vertices());
  EXPECT_EQ(0, shape.num_edges());
  EXPECT_EQ(0, shape.num_chains());
  EXPECT_EQ(2, shape.dimension());
  EXPECT_FALSE(shape.contains_origin());
}

TEST(LaxLoop, NonEmptyLoop) {
  // Test vector<S2Point> constructor.
  vector<S2Point> vertices = s2textformat::ParsePoints("0:0, 0:1, 1:1, 1:0");
  LaxLoop shape(vertices);
  EXPECT_EQ(vertices.size(), shape.num_vertices());
  EXPECT_EQ(vertices.size(), shape.num_edges());
  EXPECT_EQ(1, shape.num_chains());
  EXPECT_EQ(0, shape.chain(0).start);
  EXPECT_EQ(vertices.size(), shape.chain(0).length);
  for (int i = 0; i < vertices.size(); ++i) {
    EXPECT_EQ(vertices[i], shape.vertex(i));
    auto edge = shape.edge(i);
    EXPECT_EQ(vertices[i], edge.v0);
    EXPECT_EQ(vertices[(i + 1) % vertices.size()], edge.v1);
  }
  EXPECT_EQ(2, shape.dimension());
  EXPECT_TRUE(shape.has_interior());
  EXPECT_FALSE(shape.contains_origin());
}

TEST(ClosedLaxPolyline, NoInterior) {
  vector<S2Point> vertices = s2textformat::ParsePoints("0:0, 0:1, 1:1, 1:0");
  ClosedLaxPolyline shape(vertices);
  EXPECT_EQ(1, shape.dimension());
  EXPECT_FALSE(shape.has_interior());
  EXPECT_FALSE(shape.contains_origin());
}

TEST(LaxPolygon, EmptyPolygon) {
  LaxPolygon shape((S2Polygon()));
  EXPECT_EQ(0, shape.num_loops());
  EXPECT_EQ(0, shape.num_vertices());
  EXPECT_EQ(0, shape.num_edges());
  EXPECT_EQ(0, shape.num_chains());
  EXPECT_EQ(2, shape.dimension());
  EXPECT_TRUE(shape.has_interior());
  EXPECT_FALSE(shape.contains_origin());
}

TEST(LaxPolygon, FullPolygon) {
  LaxPolygon shape(S2Polygon(s2textformat::MakeLoop("full")));
  EXPECT_EQ(1, shape.num_loops());
  EXPECT_EQ(0, shape.num_vertices());
  EXPECT_EQ(0, shape.num_edges());
  EXPECT_EQ(1, shape.num_chains());
  EXPECT_EQ(2, shape.dimension());
  EXPECT_TRUE(shape.has_interior());
  EXPECT_TRUE(shape.contains_origin());
}

TEST(LaxPolygon, SingleVertexPolygon) {
  // S2Polygon doesn't support single-vertex loops, so we need to construct
  // the LaxPolygon directly.
  vector<vector<S2Point>> loops;
  loops.push_back(s2textformat::ParsePoints("0:0"));
  LaxPolygon shape(loops);
  EXPECT_EQ(1, shape.num_loops());
  EXPECT_EQ(1, shape.num_vertices());
  EXPECT_EQ(1, shape.num_edges());
  EXPECT_EQ(1, shape.num_chains());
  EXPECT_EQ(0, shape.chain(0).start);
  EXPECT_EQ(1, shape.chain(0).length);
  auto edge = shape.edge(0);
  EXPECT_EQ(loops[0][0], edge.v0);
  EXPECT_EQ(loops[0][0], edge.v1);
  EXPECT_TRUE(edge == shape.chain_edge(0, 0));
  EXPECT_EQ(2, shape.dimension());
  EXPECT_TRUE(shape.has_interior());
  EXPECT_FALSE(shape.contains_origin());
}

TEST(LaxPolygon, SingleLoopPolygon) {
  // Test S2Polygon constructor.
  vector<S2Point> vertices = s2textformat::ParsePoints("0:0, 0:1, 1:1, 1:0");
  LaxPolygon shape(S2Polygon(absl::MakeUnique<S2Loop>(vertices)));
  EXPECT_EQ(1, shape.num_loops());
  EXPECT_EQ(vertices.size(), shape.num_vertices());
  EXPECT_EQ(vertices.size(), shape.num_loop_vertices(0));
  EXPECT_EQ(vertices.size(), shape.num_edges());
  EXPECT_EQ(1, shape.num_chains());
  EXPECT_EQ(0, shape.chain(0).start);
  EXPECT_EQ(vertices.size(), shape.chain(0).length);
  for (int i = 0; i < vertices.size(); ++i) {
    EXPECT_EQ(vertices[i], shape.loop_vertex(0, i));
    auto edge = shape.edge(i);
    EXPECT_EQ(vertices[i], edge.v0);
    EXPECT_EQ(vertices[(i + 1) % vertices.size()], edge.v1);
    EXPECT_EQ(edge.v0, shape.chain_edge(0, i).v0);
    EXPECT_EQ(edge.v1, shape.chain_edge(0, i).v1);
  }
  EXPECT_EQ(2, shape.dimension());
  EXPECT_TRUE(shape.has_interior());
  EXPECT_FALSE(shape.contains_origin());
}

TEST(LaxPolygon, MultiLoopPolygon) {
  // Test vector<vector<S2Point>> constructor.  Make sure that the loops are
  // oriented so that the interior of the polygon is always on the left.
  vector<LaxPolygon::Loop> loops = {
    s2textformat::ParsePoints("0:0, 0:3, 3:3"),  // CCW
    s2textformat::ParsePoints("1:1, 2:2, 1:2")   // CW
  };
  LaxPolygon shape(loops);

  EXPECT_EQ(loops.size(), shape.num_loops());
  int num_vertices = 0;
  EXPECT_EQ(loops.size(), shape.num_chains());
  for (int i = 0; i < loops.size(); ++i) {
    EXPECT_EQ(loops[i].size(), shape.num_loop_vertices(i));
    EXPECT_EQ(num_vertices, shape.chain(i).start);
    EXPECT_EQ(loops[i].size(), shape.chain(i).length);
    for (int j = 0; j < loops[i].size(); ++j) {
      EXPECT_EQ(loops[i][j], shape.loop_vertex(i, j));
      auto edge = shape.edge(num_vertices + j);
      EXPECT_EQ(loops[i][j], edge.v0);
      EXPECT_EQ(loops[i][(j + 1) % loops[i].size()], edge.v1);
    }
    num_vertices += loops[i].size();
  }
  EXPECT_EQ(num_vertices, shape.num_vertices());
  EXPECT_EQ(num_vertices, shape.num_edges());
  EXPECT_EQ(2, shape.dimension());
  EXPECT_TRUE(shape.has_interior());
  EXPECT_FALSE(shape.contains_origin());
}

TEST(LaxPolygon, DegenerateLoops) {
  vector<LaxPolygon::Loop> loops = {
    s2textformat::ParsePoints("1:1, 1:2, 2:2, 1:2, 1:3, 1:2, 1:1"),
    s2textformat::ParsePoints("0:0, 0:3, 0:6, 0:9, 0:6, 0:3, 0:0"),
    s2textformat::ParsePoints("5:5, 6:6")
  };
  LaxPolygon shape(loops);
  EXPECT_FALSE(shape.contains_origin());
}

TEST(LaxPolygon, InvertedLoops) {
  vector<LaxPolygon::Loop> loops = {
    s2textformat::ParsePoints("1:2, 1:1, 2:2"),
    s2textformat::ParsePoints("3:4, 3:3, 4:4")
  };
  LaxPolygon shape(loops);
  EXPECT_TRUE(shape.contains_origin());
}

TEST(LaxPolygon, PartiallyDegenerateLoops) {
  for (int iter = 0; iter < 100; ++iter) {
    // First we construct a long convoluted edge chain that follows the
    // S2CellId Hilbert curve.  At some random point along the curve, we
    // insert a small triangular loop.
    vector<LaxPolygon::Loop> loops(1);
    LaxPolygon::Loop* loop = &loops[0];
    int const num_vertices = 100;
    S2CellId start = S2Testing::GetRandomCellId(S2CellId::kMaxLevel - 1);
    S2CellId end = start.advance_wrap(num_vertices);
    S2CellId loop_cellid = start.advance_wrap(
        S2Testing::rnd.Uniform(num_vertices - 2) + 1);
    vector<S2Point> triangle;
    for (S2CellId cellid = start; cellid != end; cellid = cellid.next_wrap()) {
      if (cellid == loop_cellid) {
        // Insert a small triangular loop.  We save the loop so that we can
        // test whether it contains the origin later.
        triangle.push_back(cellid.child(0).ToPoint());
        triangle.push_back(cellid.child(1).ToPoint());
        triangle.push_back(cellid.child(2).ToPoint());
        loop->insert(loop->end(), triangle.begin(), triangle.end());
        loop->push_back(cellid.child(0).ToPoint());
      } else {
        loop->push_back(cellid.ToPoint());
      }
    }
    // Now we retrace our steps, except that we skip the three edges that form
    // the triangular loop above.
    for (S2CellId cellid = end; cellid != start; cellid = cellid.prev_wrap()) {
      if (cellid == loop_cellid) {
        loop->push_back(cellid.child(0).ToPoint());
      } else {
        loop->push_back(cellid.ToPoint());
      }
    }
    LaxPolygon shape(loops);
    S2Loop triangle_loop(triangle);
    EXPECT_EQ(triangle_loop.Contains(S2::Origin()), shape.contains_origin());
  }
}

void CompareS2LoopToShape(S2Loop const& loop, S2Shape* shape) {
  S2ShapeIndex index;
  index.Add(shape);
  S2Cap cap = loop.GetCapBound();
  for (int iter = 0; iter < 100; ++iter) {
    S2Point point = S2Testing::SamplePoint(cap);
    EXPECT_EQ(loop.Contains(point), index.ShapeContains(shape, point));
  }
}

TEST(LaxPolygon, CompareToS2Loop) {
  for (int iter = 0; iter < 100; ++iter) {
    S2Testing::Fractal fractal;
    fractal.set_max_level(S2Testing::rnd.Uniform(5));
    fractal.set_fractal_dimension(1 + S2Testing::rnd.RandDouble());
    S2Point center = S2Testing::RandomPoint();
    unique_ptr<S2Loop> loop(fractal.MakeLoop(
        S2Testing::GetRandomFrameAt(center), S1Angle::Degrees(5)));

    // Compare S2Loop to LaxLoop.
    CompareS2LoopToShape(*loop, new LaxLoop(*loop));

    // Compare S2Loop to LaxPolygon.
    vector<LaxPolygon::Loop> loops(
        1, vector<S2Point>(&loop->vertex(0),
                           &loop->vertex(0) + loop->num_vertices()));
    CompareS2LoopToShape(*loop, new LaxPolygon(loops));
  }
}

TEST(LaxPolyline, NoVertices) {
  vector<S2Point> vertices;
  LaxPolyline shape(vertices);
  EXPECT_EQ(0, shape.num_edges());
  EXPECT_EQ(0, shape.num_chains());
  EXPECT_EQ(1, shape.dimension());
}

TEST(LaxPolyline, OneVertex) {
  vector<S2Point> vertices = {S2Point(1, 0, 0)};
  LaxPolyline shape(vertices);
  EXPECT_EQ(0, shape.num_edges());
  EXPECT_EQ(0, shape.num_chains());
  EXPECT_EQ(1, shape.dimension());
}

TEST(LaxPolyline, EdgeAccess) {
  vector<S2Point> vertices = s2textformat::ParsePoints("0:0, 0:1, 1:1");
  LaxPolyline shape(vertices);
  EXPECT_EQ(2, shape.num_edges());
  EXPECT_EQ(1, shape.num_chains());
  EXPECT_EQ(0, shape.chain(0).start);
  EXPECT_EQ(2, shape.chain(0).length);
  EXPECT_EQ(1, shape.dimension());
  auto edge0 = shape.edge(0);
  EXPECT_EQ(vertices[0], edge0.v0);
  EXPECT_EQ(vertices[1], edge0.v1);
  auto edge1 = shape.edge(1);
  EXPECT_EQ(vertices[1], edge1.v0);
  EXPECT_EQ(vertices[2], edge1.v1);
}

TEST(EdgeVectorShape, EdgeAccess) {
  EdgeVectorShape shape;
  S2Testing::rnd.Reset(FLAGS_s2_random_seed);
  int const kNumEdges = 100;
  for (int i = 0; i < kNumEdges; ++i) {
    S2Point a = S2Testing::RandomPoint();  // Control the evaluation order
    shape.Add(a, S2Testing::RandomPoint());
  }
  EXPECT_EQ(kNumEdges, shape.num_edges());
  EXPECT_EQ(kNumEdges, shape.num_chains());
  EXPECT_EQ(1, shape.dimension());
  S2Testing::rnd.Reset(FLAGS_s2_random_seed);
  for (int i = 0; i < kNumEdges; ++i) {
    EXPECT_EQ(i, shape.chain(i).start);
    EXPECT_EQ(1, shape.chain(i).length);
    auto edge = shape.edge(i);
    EXPECT_EQ(S2Testing::RandomPoint(), edge.v0);
    EXPECT_EQ(S2Testing::RandomPoint(), edge.v1);
  }
}

TEST(EdgeVectorShape, SingletonConstructor) {
  S2Point a(1, 0, 0), b(0, 1, 0);
  EdgeVectorShape shape(a, b);
  EXPECT_EQ(1, shape.num_edges());
  EXPECT_EQ(1, shape.num_chains());
  auto edge = shape.edge(0);
  EXPECT_EQ(a, edge.v0);
  EXPECT_EQ(b, edge.v1);
}

TEST(PointVectorShape, ConstructionAndAccess) {
  std::vector<S2Point> points;
  S2Testing::rnd.Reset(FLAGS_s2_random_seed);
  int const kNumPoints = 100;
  for (int i = 0; i < kNumPoints; ++i) {
    points.push_back(S2Testing::RandomPoint());
  }
  PointVectorShape shape(&points);

  EXPECT_EQ(kNumPoints, shape.num_edges());
  EXPECT_EQ(kNumPoints, shape.num_chains());
  EXPECT_EQ(0, shape.dimension());
  S2Testing::rnd.Reset(FLAGS_s2_random_seed);
  for (int i = 0; i < kNumPoints; ++i) {
    EXPECT_EQ(i, shape.chain(i).start);
    EXPECT_EQ(1, shape.chain(i).length);
    auto edge = shape.edge(i);
    S2Point pt = S2Testing::RandomPoint();
    EXPECT_EQ(pt, edge.v0);
    EXPECT_EQ(pt, edge.v1);
  }
}

TEST(VertexIdLaxLoop, EmptyLoop) {
  VertexIdLaxLoop shape(vector<int32>(), nullptr);
  EXPECT_EQ(0, shape.num_edges());
  EXPECT_EQ(0, shape.num_vertices());
  EXPECT_EQ(1, shape.num_chains());
  EXPECT_EQ(2, shape.dimension());
  EXPECT_TRUE(shape.has_interior());
  EXPECT_FALSE(shape.contains_origin());
}

TEST(VertexIdLaxLoop, InvertedLoop) {
  vector<S2Point> vertex_array =
      s2textformat::ParsePoints("0:0, 0:1, 1:1, 1:0");
  vector<int32> vertex_ids { 0, 3, 2, 1 };  // Inverted.
  VertexIdLaxLoop shape(vertex_ids, &vertex_array[0]);
  EXPECT_EQ(4, shape.num_edges());
  EXPECT_EQ(4, shape.num_vertices());
  EXPECT_EQ(1, shape.num_chains());
  EXPECT_EQ(0, shape.chain(0).start);
  EXPECT_EQ(4, shape.chain(0).length);
  EXPECT_EQ(&vertex_array[0], &shape.vertex(0));
  EXPECT_EQ(&vertex_array[3], &shape.vertex(1));
  EXPECT_EQ(&vertex_array[2], &shape.vertex(2));
  EXPECT_EQ(&vertex_array[1], &shape.vertex(3));
  EXPECT_EQ(2, shape.dimension());
  EXPECT_TRUE(shape.has_interior());
  EXPECT_TRUE(shape.contains_origin());
}

TEST(S2LoopOwningShape, Ownership) {
  // Debug mode builds will catch any memory leak below.
  auto loop = absl::MakeUnique<S2Loop>(S2Loop::kEmpty());
  S2LoopOwningShape shape(std::move(loop));
}

TEST(S2PolygonOwningShape, Ownership) {
  // Debug mode builds will catch any memory leak below.
  vector<unique_ptr<S2Loop>> loops;
  auto polygon = absl::MakeUnique<S2Polygon>(std::move(loops));
  S2PolygonOwningShape shape(std::move(polygon));
}

TEST(S2PolylineOwningShape, Ownership) {
  // Debug mode builds will catch any memory leak below.
  vector<S2Point> vertices;
  auto polyline = absl::MakeUnique<S2Polyline>(vertices);
  S2PolylineOwningShape shape(std::move(polyline));
}

// A set of edge pairs within an S2ShapeIndex.
using EdgePairVector = std::vector<std::pair<ShapeEdgeId, ShapeEdgeId>>;

EdgePairVector GetCrossings(S2ShapeIndex const& index, CrossingType type) {
  EdgePairVector edge_pairs;
  VisitCrossings(index, type,
                 [&edge_pairs](ShapeEdge const& a, ShapeEdge const& b, bool) {
                   edge_pairs.push_back(std::make_pair(a.id(), b.id()));
                   return true;  // Continue visiting.
                 });
  if (edge_pairs.size() > 1) {
    std::sort(edge_pairs.begin(), edge_pairs.end());
    edge_pairs.erase(std::unique(edge_pairs.begin(), edge_pairs.end()),
                     edge_pairs.end());
  }
  return edge_pairs;
}

// An iterator that advances through all edges in an S2ShapeIndex.
// TODO(ericv): Consider moving this to S2ShapeIndex, or s2shapeutil?
class EdgeIterator {
 public:
  explicit EdgeIterator(S2ShapeIndex const& index)
      : index_(index), shape_id_(-1), num_edges_(0), edge_id_(-1) {
    Next();
  }
  int32 shape_id() const { return shape_id_; }
  int32 edge_id() const { return edge_id_; }
  ShapeEdgeId shape_edge_id() const {
    return ShapeEdgeId(shape_id_, edge_id_);
  }
  S2Shape::Edge edge() {
    return index_.shape(shape_id_)->edge(edge_id_);
  }
  bool Done() const {
    return shape_id() >= index_.num_shape_ids();
  }
  void Next() {
    while (++edge_id_ >= num_edges_) {
      if (++shape_id_ >= index_.num_shape_ids()) break;
      S2Shape* shape = index_.shape(shape_id_);
      num_edges_= (shape == nullptr) ? 0 : shape->num_edges();
      edge_id_ = -1;
    }
  }

 private:
  S2ShapeIndex const& index_;
  int32 shape_id_;
  int32 num_edges_;
  int32 edge_id_;
};

EdgePairVector GetCrossingEdgePairsBruteForce(S2ShapeIndex const& index,
                                              CrossingType type) {
  EdgePairVector result;
  int min_sign = (type == CrossingType::ALL) ? 0 : 1;
  for (EdgeIterator a_iter(index); !a_iter.Done(); a_iter.Next()) {
    auto a = a_iter.edge();
    EdgeIterator b_iter = a_iter;
    for (b_iter.Next(); !b_iter.Done(); b_iter.Next()) {
      auto b = b_iter.edge();
      if (S2EdgeUtil::CrossingSign(a.v0, a.v1, b.v0, b.v1) >= min_sign) {
        result.push_back(
            std::make_pair(a_iter.shape_edge_id(), b_iter.shape_edge_id()));
      }
    }
  }
  return result;
}

std::ostream& operator<<(std::ostream& os,
                         std::pair<ShapeEdgeId, ShapeEdgeId> const& pair) {
  return os << "(" << pair.first << "," << pair.second << ")";
}

void TestGetCrossingEdgePairs(S2ShapeIndex const& index, CrossingType type) {
  EdgePairVector expected = GetCrossingEdgePairsBruteForce(index, type);
  EdgePairVector actual = GetCrossings(index, type);
  if (actual != expected) {
    ADD_FAILURE() << "Unexpected edge pairs; see details below."
                  << "\nExpected number of edge pairs: " << expected.size()
                  << "\nActual number of edge pairs: " << actual.size();
    for (auto const& edge_pair : expected) {
      if (std::count(actual.begin(), actual.end(), edge_pair) != 1) {
        std::cout << "Missing value: " << edge_pair << std::endl;
      }
    }
    for (auto const& edge_pair : actual) {
      if (std::count(expected.begin(), expected.end(), edge_pair) != 1) {
        std::cout << "Extra value: " << edge_pair << std::endl;
      }
    }
  }
}

TEST(GetCrossingEdgePairs, NoIntersections) {
  S2ShapeIndex index;
  TestGetCrossingEdgePairs(index, CrossingType::ALL);
  TestGetCrossingEdgePairs(index, CrossingType::INTERIOR);
}

TEST(GetCrossingEdgePairs, EdgeGrid) {
  int const kGridSize = 10;  // (kGridSize + 1) * (kGridSize + 1) crossings
  S2ShapeIndex index;
  EdgeVectorShape* shape = new EdgeVectorShape;
  for (int i = 0; i <= kGridSize; ++i) {
    shape->Add(S2LatLng::FromDegrees(0, i).ToPoint(),
              S2LatLng::FromDegrees(kGridSize, i).ToPoint());
    shape->Add(S2LatLng::FromDegrees(i, 0).ToPoint(),
              S2LatLng::FromDegrees(i, kGridSize).ToPoint());
  }
  index.Add(shape);
  TestGetCrossingEdgePairs(index, CrossingType::ALL);
  TestGetCrossingEdgePairs(index, CrossingType::INTERIOR);
}

class TestLaxLoop : public LaxLoop {
 public:
  explicit TestLaxLoop(string const& vertex_str) {
    vector<S2Point> vertices = s2textformat::ParsePoints(vertex_str);
    Init(vertices);
  }
};

TEST(ResolveComponents, NoComponents) {
  vector<vector<S2Shape*>> faces, components;
  ResolveComponents(components, &faces);
  EXPECT_EQ(0, faces.size());
}

TEST(ResolveComponents, OneLoop) {
  TestLaxLoop a0("0:0, 1:0, 0:1");  // Outer face
  TestLaxLoop a1("0:0, 0:1, 1:0");
  vector<vector<S2Shape*>> faces, components = {{&a0, &a1}};
  ResolveComponents(components, &faces);
  EXPECT_EQ(2, faces.size());
}

TEST(ResolveComponents, TwoLoopsSameComponent) {
  TestLaxLoop a0("0:0, 1:0, 0:1");  // Outer face
  TestLaxLoop a1("0:0, 0:1, 1:0");
  TestLaxLoop a2("1:0, 0:1, 1:1");
  vector<vector<S2Shape*>> faces, components = {{&a0, &a1, &a2}};
  ResolveComponents(components, &faces);
  EXPECT_EQ(3, faces.size());
}

TEST(ResolveComponents, TwoNestedLoops) {
  TestLaxLoop a0("0:0, 3:0, 0:3");  // Outer face
  TestLaxLoop a1("0:0, 0:3, 3:0");
  TestLaxLoop b0("1:1, 2:0, 0:2");  // Outer face
  TestLaxLoop b1("1:1, 0:2, 2:0");
  vector<vector<S2Shape*>> faces, components = {{&a0, &a1}, {&b0, &b1}};
  ResolveComponents(components, &faces);
  EXPECT_EQ(3, faces.size());
  EXPECT_EQ((vector<S2Shape*>{&b0, &a1}), faces[0]);
}

TEST(ResolveComponents, TwoLoopsDifferentComponents) {
  TestLaxLoop a0("0:0, 1:0, 0:1");  // Outer face
  TestLaxLoop a1("0:0, 0:1, 1:0");
  TestLaxLoop b0("0:2, 1:2, 0:3");  // Outer face
  TestLaxLoop b1("0:2, 0:3, 1:2");
  vector<vector<S2Shape*>> faces, components = {{&a0, &a1}, {&b0, &b1}};
  ResolveComponents(components, &faces);
  EXPECT_EQ(3, faces.size());
  EXPECT_EQ((vector<S2Shape*>{&a0, &b0}), faces[2]);
}

TEST(ResolveComponents, OneDegenerateLoop) {
  TestLaxLoop a0("0:0, 1:0, 0:0");
  vector<vector<S2Shape*>> faces, components = {{&a0}};
  ResolveComponents(components, &faces);
  EXPECT_EQ(1, faces.size());
}

TEST(ResolveComponents, TwoDegenerateLoops) {
  TestLaxLoop a0("0:0, 1:0, 0:0");
  TestLaxLoop b0("2:0, 3:0, 2:0");
  vector<vector<S2Shape*>> faces, components = {{&a0}, {&b0}};
  ResolveComponents(components, &faces);
  EXPECT_EQ(1, faces.size());
  EXPECT_EQ(2, faces[0].size());
}

static void SortFaces(vector<vector<S2Shape*>>* faces) {
  for (auto& face : *faces) {
    std::sort(face.begin(), face.end());
  }
  std::sort(faces->begin(), faces->end());
}

TEST(ResolveComponents, ComplexTest1) {
  // Loops at index 0 are the outer (clockwise) loops.
  // Component "a" consists of 4 adjacent squares forming a larger square.
  TestLaxLoop a0("0:0, 25:0, 50:0, 50:25, 50:50, 25:50, 0:50, 0:50");
  TestLaxLoop a1("0:0, 0:25, 25:25, 25:0");
  TestLaxLoop a2("0:25, 0:50, 25:50, 25:25");
  TestLaxLoop a3("25:0, 25:25, 50:25, 50:0");
  TestLaxLoop a4("25:25, 25:50, 50:50, 50:25");
  // Component "b" consists of a degenerate loop to the left of "a".
  TestLaxLoop b0("0:-10, 10:-10");
  // Components "a1_a", "a1_b", and "a1_c" are located within "a1".
  TestLaxLoop a1_a0("5:5, 20:5, 20:10, 5:10");
  TestLaxLoop a1_a1("5:5, 5:10, 10:10, 10:5");
  TestLaxLoop a1_a2("10:5, 10:10, 15:10, 15:5");
  TestLaxLoop a1_a3("15:5, 15:10, 20:10, 20:5");
  TestLaxLoop a1_b0("5:15, 20:15, 20:20, 5:20");
  TestLaxLoop a1_b1("5:15, 5:20, 20:20, 20:15");
  TestLaxLoop a1_c0("2:5, 2:10, 2:5");
  // Two components located inside "a1_a2" and "a1_a3".
  TestLaxLoop a1_a2_a0("11:6, 14:6, 14:9, 11:9");
  TestLaxLoop a1_a2_a1("11:6, 11:9, 14:9, 14:6");
  TestLaxLoop a1_a3_a0("16:6, 19:9, 16:6");
  // Five component located inside "a3" and "a4".
  TestLaxLoop a3_a0("30:5, 45:5, 45:20, 30:20");
  TestLaxLoop a3_a1("30:5, 30:20, 45:20, 45:5");
  TestLaxLoop a4_a0("30:30, 40:30, 30:30");
  TestLaxLoop a4_b0("30:35, 40:35, 30:35");
  TestLaxLoop a4_c0("30:40, 40:40, 30:40");
  TestLaxLoop a4_d0("30:45, 40:45, 30:45");
  vector<vector<S2Shape*>> components = {
    {&a0, &a1, &a2, &a3, &a4},
    {&b0},
    {&a1_a0, &a1_a1, &a1_a2, &a1_a3},
    {&a1_b0, &a1_b1},
    {&a1_c0},
    {&a1_a2_a0, &a1_a2_a1},
    {&a1_a3_a0},
    {&a3_a0, &a3_a1},
    {&a4_a0},
    {&a4_b0},
    {&a4_c0},
    {&a4_d0}
  };
  vector<vector<S2Shape*>> expected_faces = {
    {&a0, &b0},
    {&a1, &a1_a0, &a1_b0, &a1_c0},
    {&a1_a1},
    {&a1_a2, &a1_a2_a0},
    {&a1_a2_a1},
    {&a1_a3, &a1_a3_a0},
    {&a1_b1},
    {&a2},
    {&a3, &a3_a0},
    {&a3_a1},
    {&a4, &a4_a0, &a4_b0, &a4_c0, &a4_d0}
  };
  vector<vector<S2Shape*>> faces;
  ResolveComponents(components, &faces);
  EXPECT_EQ(expected_faces.size(), faces.size());
  SortFaces(&expected_faces);
  SortFaces(&faces);
  EXPECT_EQ(expected_faces, faces);
}

}  // namespace s2shapeutil
