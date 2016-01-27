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

#include "s2shapeutil.h"

#include <memory>
#include <vector>
#include <gtest/gtest.h>
#include "s2cap.h"
#include "s2loop.h"
#include "s2polygon.h"
#include "s2polyline.h"
#include "s2testing.h"
#include "s2textformat.h"

using std::unique_ptr;
using std::unique_ptr;
using std::vector;

namespace s2shapeutil {

TEST(S2EdgeVectorShape, EdgeAccess) {
  S2EdgeVectorShape shape;
  S2Testing::rnd.Reset(FLAGS_s2_random_seed);
  int const kNumEdges = 100;
  for (int i = 0; i < kNumEdges; ++i) {
    S2Point a = S2Testing::RandomPoint();  // Control the evaluation order
    shape.Add(a, S2Testing::RandomPoint());
  }
  EXPECT_EQ(kNumEdges, shape.num_edges());
  S2Testing::rnd.Reset(FLAGS_s2_random_seed);
  for (int i = 0; i < kNumEdges; ++i) {
    S2Point const *a, *b;
    shape.GetEdge(i, &a, &b);
    EXPECT_EQ(S2Testing::RandomPoint(), *a);
    EXPECT_EQ(S2Testing::RandomPoint(), *b);
  }
}

TEST(S2EdgeVectorShape, SingletonConstructor) {
  S2Point a(1, 0, 0), b(0, 1, 0);
  S2EdgeVectorShape shape(a, b);
  EXPECT_EQ(1, shape.num_edges());
  S2Point const *pa, *pb;
  shape.GetEdge(0, &pa, &pb);
  EXPECT_EQ(a, *pa);
  EXPECT_EQ(b, *pb);
}

TEST(S2LoopOwningShape, Ownership) {
  // Debug mode builds will catch any memory leak below.
  S2Loop* loop = new S2Loop(S2Loop::kEmpty());
  S2LoopOwningShape shape(loop);
}

TEST(S2PolygonOwningShape, Ownership) {
  // Debug mode builds will catch any memory leak below.
  vector<S2Loop*> loops;
  S2Polygon* polygon = new S2Polygon(&loops);
  S2PolygonOwningShape shape(polygon);
}

TEST(S2PolylineOwningShape, Ownership) {
  // Debug mode builds will catch any memory leak below.
  vector<S2Point> vertices;
  S2Polyline* polyline = new S2Polyline(vertices);
  S2PolylineOwningShape shape(polyline);
}

TEST(FaceLoopShape, EmptyLoop) {
  // Test S2Loop constructor.
  FaceLoopShape shape;
  shape.Init(S2Loop(S2Loop::kEmpty()));
  EXPECT_EQ(0, shape.num_vertices());
  EXPECT_EQ(0, shape.num_edges());
  EXPECT_TRUE(shape.has_interior());
  EXPECT_FALSE(shape.contains_origin());
}

TEST(FaceLoopShape, NonEmptyLoop) {
  // Test vector<S2Point> constructor.
  vector<S2Point> vertices;
  s2textformat::ParsePoints("0:0, 0:1, 1:1, 1:0", &vertices);
  FaceLoopShape shape(vertices);
  EXPECT_EQ(vertices.size(), shape.num_vertices());
  for (int i = 0; i < vertices.size(); ++i) {
    EXPECT_EQ(vertices[i], shape.vertex(i));
    S2Point const *v0, *v1;
    shape.GetEdge(i, &v0, &v1);
    EXPECT_EQ(vertices[i], *v0);
    EXPECT_EQ(vertices[(i + 1) % vertices.size()], *v1);
  }
  EXPECT_EQ(vertices.size(), shape.num_edges());
  EXPECT_TRUE(shape.has_interior());
  EXPECT_FALSE(shape.contains_origin());
}

TEST(FaceShape, EmptyPolygon) {
  S2Polygon a;
  FaceShape shape(a);
  EXPECT_EQ(0, shape.num_loops());
  EXPECT_EQ(0, shape.num_vertices());
  EXPECT_EQ(0, shape.num_edges());
  EXPECT_TRUE(shape.has_interior());
  EXPECT_FALSE(shape.contains_origin());
}

TEST(FaceShape, SingleLoopPolygon) {
  // Test S2Polygon constructor.
  vector<S2Point> vertices;
  s2textformat::ParsePoints("0:0, 0:1, 1:1, 1:0", &vertices);
  FaceShape shape(S2Polygon(new S2Loop(vertices)));
  EXPECT_EQ(1, shape.num_loops());
  EXPECT_EQ(vertices.size(), shape.num_vertices());
  EXPECT_EQ(vertices.size(), shape.num_loop_vertices(0));
  for (int i = 0; i < vertices.size(); ++i) {
    EXPECT_EQ(vertices[i], shape.loop_vertex(0, i));
    S2Point const *v0, *v1;
    shape.GetEdge(i, &v0, &v1);
    EXPECT_EQ(vertices[i], *v0);
    EXPECT_EQ(vertices[(i + 1) % vertices.size()], *v1);
  }
  EXPECT_EQ(vertices.size(), shape.num_edges());
  EXPECT_TRUE(shape.has_interior());
  EXPECT_FALSE(shape.contains_origin());
}

TEST(FaceShape, MultiLoopPolygon) {
  // Test vector<vector<S2Point>> constructor.
  vector<FaceShape::Loop> loops(2);
  s2textformat::ParsePoints("1:1, 1:2, 2:2", &loops[0]);
  s2textformat::ParsePoints("0:0, 0:3, 3:3", &loops[1]);
  FaceShape shape(loops);

  EXPECT_EQ(loops.size(), shape.num_loops());
  int num_vertices = 0;
  for (int i = 0; i < loops.size(); ++i) {
    EXPECT_EQ(loops[i].size(), shape.num_loop_vertices(i));
    for (int j = 0; j < loops[i].size(); ++j) {
      EXPECT_EQ(loops[i][j], shape.loop_vertex(i, j));
      S2Point const *v0, *v1;
      shape.GetEdge(num_vertices + j, &v0, &v1);
      EXPECT_EQ(loops[i][j], *v0);
      EXPECT_EQ(loops[i][(j + 1) % loops[i].size()], *v1);
    }
    num_vertices += loops[i].size();
  }
  EXPECT_EQ(num_vertices, shape.num_vertices());
  EXPECT_EQ(num_vertices, shape.num_edges());
  EXPECT_TRUE(shape.has_interior());
  EXPECT_FALSE(shape.contains_origin());
}

TEST(FaceShape, DegenerateLoops) {
  vector<FaceShape::Loop> loops(3);
  s2textformat::ParsePoints("1:1, 1:2, 2:2, 1:2, 1:3, 1:2, 1:1", &loops[0]);
  s2textformat::ParsePoints("0:0, 0:3, 0:6, 0:9, 0:6, 0:3, 0:0", &loops[1]);
  s2textformat::ParsePoints("5:5, 6:6",  &loops[2]);
  FaceShape shape(loops);
  EXPECT_FALSE(shape.contains_origin());
}

TEST(FaceShape, InvertedLoops) {
  vector<FaceShape::Loop> loops(2);
  s2textformat::ParsePoints("1:2, 1:1, 2:2", &loops[0]);
  s2textformat::ParsePoints("3:4, 3:3, 4:4", &loops[1]);
  FaceShape shape(loops);
  EXPECT_TRUE(shape.contains_origin());
}

TEST(FaceShape, PartiallyDegenerateLoops) {
  for (int iter = 0; iter < 100; ++iter) {
    // First we construct a long convoluted edge chain that follows the
    // S2CellId Hilbert curve.  At some random point along the curve, we
    // insert a small triangular loop.
    vector<FaceShape::Loop> loops(1);
    FaceShape::Loop* loop = &loops[0];
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
    FaceShape shape(loops);
    S2Loop triangle_loop(triangle);
    EXPECT_EQ(triangle_loop.Contains(S2::Origin()), shape.contains_origin());
  }
}

void CompareS2LoopToShape(S2Loop const* loop, S2Shape* shape) {
  S2ShapeIndex index;
  index.Add(shape);
  S2Cap cap = loop->GetCapBound();
  for (int iter = 0; iter < 100; ++iter) {
    S2Point point = S2Testing::SamplePoint(cap);
    EXPECT_EQ(loop->Contains(point), index.ShapeContains(shape, point));
  }
}

TEST(FaceShapes, CompareToS2Loop) {
  for (int iter = 0; iter < 100; ++iter) {
    S2Testing::Fractal fractal;
    fractal.set_max_level(S2Testing::rnd.Uniform(5));
    fractal.set_fractal_dimension(1 + S2Testing::rnd.RandDouble());
    S2Point center = S2Testing::RandomPoint();
    unique_ptr<S2Loop> loop(fractal.MakeLoop(
        S2Testing::GetRandomFrameAt(center), S1Angle::Degrees(5)));

    // Compare S2Loop to FaceLoopShape.
    CompareS2LoopToShape(loop.get(), new FaceLoopShape(*loop));

    // Compare S2Loop to FaceShape.
    vector<FaceShape::Loop> loops(
        1, vector<S2Point>(&loop->vertex(0),
                           &loop->vertex(0) + loop->num_vertices()));
    CompareS2LoopToShape(loop.get(), new FaceShape(loops));
  }
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
  void GetEdge(S2Point const** v0, S2Point const** v1) {
    index_.shape(shape_id_)->GetEdge(edge_id_, v0, v1);
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

void GetCrossingEdgePairsBruteForce(S2ShapeIndex const& index,
                                    CrossingType type,
                                    EdgePairList* edge_pairs) {
  int min_crossing_value = (type == CrossingType::ALL) ? 0 : 1;
  for (EdgeIterator a_iter(index); !a_iter.Done(); a_iter.Next()) {
    S2Point const *a0, *a1;
    a_iter.GetEdge(&a0, &a1);
    EdgeIterator b_iter = a_iter;
    for (b_iter.Next(); !b_iter.Done(); b_iter.Next()) {
      S2Point const *b0, *b1;
      b_iter.GetEdge(&b0, &b1);
      if (S2EdgeUtil::CrossingSign(*a0, *a1, *b0, *b1) >= min_crossing_value) {
        edge_pairs->push_back(
            std::make_pair(a_iter.shape_edge_id(), b_iter.shape_edge_id()));
      }
    }
  }
}

std::ostream& operator<<(std::ostream& os,
                         std::pair<ShapeEdgeId, ShapeEdgeId> const& pair) {
  return os << "(" << pair.first << "," << pair.second << ")";
}

void TestGetCrossingEdgePairs(S2ShapeIndex const& index, CrossingType type) {
  EdgePairList expected, actual;
  GetCrossingEdgePairsBruteForce(index, type, &expected);
  GetCrossingEdgePairs(index, type, &actual);
  if (actual != expected) {
    ADD_FAILURE() << "Unexpected edge pairs; see details below."
                  << "\nExpected number of edge pairs: " << expected.size()
                  << "\nActual number of edge pairs: " << actual.size();
    for (auto const& edge_pair : expected) {
      if (std::count(actual.begin(), actual.end(), edge_pair) != 1) {
        std::cout << "Missing value: " << edge_pair << "\n";
      }
    }
    for (auto const& edge_pair : actual) {
      if (std::count(expected.begin(), expected.end(), edge_pair) != 1) {
        std::cout << "Extra value: " << edge_pair << "\n";
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
  S2EdgeVectorShape* shape = new S2EdgeVectorShape;
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

class TestLoopShape : public FaceLoopShape {
 public:
  TestLoopShape(string const& vertex_str) {
    vector<S2Point> vertices;
    s2textformat::ParsePoints(vertex_str, &vertices);
    Init(vertices);
  }
};

TEST(ResolveLoopContainment, NoComponents) {
  vector<vector<S2Shape*>> faces, components;
  ResolveLoopContainment(components, &faces);
  EXPECT_EQ(0, faces.size());
}

TEST(ResolveLoopContainment, OneLoop) {
  TestLoopShape a0("0:0, 1:0, 0:1");  // Outer face
  TestLoopShape a1("0:0, 0:1, 1:0");
  vector<vector<S2Shape*>> faces, components = {{&a0, &a1}};
  ResolveLoopContainment(components, &faces);
  EXPECT_EQ(2, faces.size());
}

TEST(ResolveLoopContainment, TwoLoopsSameComponent) {
  TestLoopShape a0("0:0, 1:0, 0:1");  // Outer face
  TestLoopShape a1("0:0, 0:1, 1:0");
  TestLoopShape a2("1:0, 0:1, 1:1");
  vector<vector<S2Shape*>> faces, components = {{&a0, &a1, &a2}};
  ResolveLoopContainment(components, &faces);
  EXPECT_EQ(3, faces.size());
}

TEST(ResolveLoopContainment, TwoNestedLoops) {
  TestLoopShape a0("0:0, 3:0, 0:3");  // Outer face
  TestLoopShape a1("0:0, 0:3, 3:0");
  TestLoopShape b0("1:1, 2:0, 0:2");  // Outer face
  TestLoopShape b1("1:1, 0:2, 2:0");
  vector<vector<S2Shape*>> faces, components = {{&a0, &a1}, {&b0, &b1}};
  ResolveLoopContainment(components, &faces);
  EXPECT_EQ(3, faces.size());
  EXPECT_EQ((vector<S2Shape*>{&b0, &a1}), faces[0]);
}

TEST(ResolveLoopContainment, TwoLoopsDifferentComponents) {
  TestLoopShape a0("0:0, 1:0, 0:1");  // Outer face
  TestLoopShape a1("0:0, 0:1, 1:0");
  TestLoopShape b0("0:2, 1:2, 0:3");  // Outer face
  TestLoopShape b1("0:2, 0:3, 1:2");
  vector<vector<S2Shape*>> faces, components = {{&a0, &a1}, {&b0, &b1}};
  ResolveLoopContainment(components, &faces);
  EXPECT_EQ(3, faces.size());
  EXPECT_EQ((vector<S2Shape*>{&a0, &b0}), faces[2]);
}

TEST(ResolveLoopContainment, OneDegenerateLoop) {
  TestLoopShape a0("0:0, 1:0, 0:0");
  vector<vector<S2Shape*>> faces, components = {{&a0}};
  ResolveLoopContainment(components, &faces);
  EXPECT_EQ(1, faces.size());
}

TEST(ResolveLoopContainment, TwoDegenerateLoops) {
  TestLoopShape a0("0:0, 1:0, 0:0");
  TestLoopShape b0("2:0, 3:0, 2:0");
  vector<vector<S2Shape*>> faces, components = {{&a0}, {&b0}};
  ResolveLoopContainment(components, &faces);
  EXPECT_EQ(1, faces.size());
  EXPECT_EQ(2, faces[0].size());
}

static void SortFaces(vector<vector<S2Shape*>>* faces) {
  for (auto& face : *faces) {
    std::sort(face.begin(), face.end());
  }
  std::sort(faces->begin(), faces->end());
}

TEST(ResolveLoopContainment, ComplexTest1) {
  // Loops at index 0 are the outer (clockwise) loops.
  // Component "a" consists of 4 adjacent squares forming a larger square.
  TestLoopShape a0("0:0, 25:0, 50:0, 50:25, 50:50, 25:50, 0:50, 0:50");
  TestLoopShape a1("0:0, 0:25, 25:25, 25:0");
  TestLoopShape a2("0:25, 0:50, 25:50, 25:25");
  TestLoopShape a3("25:0, 25:25, 50:25, 50:0");
  TestLoopShape a4("25:25, 25:50, 50:50, 50:25");
  // Component "b" consists of a degenerate loop to the left of "a".
  TestLoopShape b0("0:-10, 10:-10");
  // Components "a1_a", "a1_b", and "a1_c" are located within "a1".
  TestLoopShape a1_a0("5:5, 20:5, 20:10, 5:10");
  TestLoopShape a1_a1("5:5, 5:10, 10:10, 10:5");
  TestLoopShape a1_a2("10:5, 10:10, 15:10, 15:5");
  TestLoopShape a1_a3("15:5, 15:10, 20:10, 20:5");
  TestLoopShape a1_b0("5:15, 20:15, 20:20, 5:20");
  TestLoopShape a1_b1("5:15, 5:20, 20:20, 20:15");
  TestLoopShape a1_c0("2:5, 2:10, 2:5");
  // Two components located inside "a1_a2" and "a1_a3".
  TestLoopShape a1_a2_a0("11:6, 14:6, 14:9, 11:9");
  TestLoopShape a1_a2_a1("11:6, 11:9, 14:9, 14:6");
  TestLoopShape a1_a3_a0("16:6, 19:9, 16:6");
  // Five component located inside "a3" and "a4".
  TestLoopShape a3_a0("30:5, 45:5, 45:20, 30:20");
  TestLoopShape a3_a1("30:5, 30:20, 45:20, 45:5");
  TestLoopShape a4_a0("30:30, 40:30, 30:30");
  TestLoopShape a4_b0("30:35, 40:35, 30:35");
  TestLoopShape a4_c0("30:40, 40:40, 30:40");
  TestLoopShape a4_d0("30:45, 40:45, 30:45");
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
  ResolveLoopContainment(components, &faces);
  EXPECT_EQ(expected_faces.size(), faces.size());
  SortFaces(&expected_faces);
  SortFaces(&faces);
  EXPECT_EQ(expected_faces, faces);
}

}  // namespace s2shapeutil
