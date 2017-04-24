// Copyright 2017 Google Inc. All Rights Reserved.
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

#include "s2/s2boundary_operation.h"

#include <memory>
#include <gtest/gtest.h>
#include "s2/third_party/absl/memory/memory.h"
#include "s2/s2builder.h"
#include "s2/s2builder_graph.h"
#include "s2/s2builder_layer.h"
#include "s2/s2builderutil_snap_functions.h"
#include "s2/s2shapeutil.h"
#include "s2/s2textformat.h"

namespace {

using std::make_pair;
using std::pair;
using std::unique_ptr;
using std::vector;

using Graph = S2Builder::Graph;
using GraphOptions = S2Builder::GraphOptions;
using DegenerateEdges = GraphOptions::DegenerateEdges;
using DuplicateEdges = GraphOptions::DuplicateEdges;
using SiblingPairs = GraphOptions::SiblingPairs;

using OpType = S2BoundaryOperation::OpType;

// Returns an S2ShapeIndex containing the given loops, polylines, and points.
// Loops should be directed so that the region's interior is on the left.
// Loops can be degenerate (they do not need to meet S2Loop requirements).
unique_ptr<S2ShapeIndex> MakeIndex(
    vector<string> const& loops,
    vector<string> const& polylines = vector<string>(),
    string const& points = string()) {
  auto index = absl::MakeUnique<S2ShapeIndex>();
  if (!loops.empty()) {
    vector<vector<S2Point>> polygon;
    for (string const& str : loops) {
      polygon.push_back(s2textformat::ParsePoints(str));
    }
    index->Add(new s2shapeutil::LaxPolygon(polygon));
  }
  for (string const& str : polylines) {
    auto vertices = s2textformat::ParsePoints(str);
    index->Add(new s2shapeutil::LaxPolyline(vertices));
  }
  // TODO(ericv): Create an S2Shape that represents a set of points.
  for (S2Point const& point : s2textformat::ParsePoints(points)) {
    vector<S2Point> vertices = { point, point };
    index->Add(new s2shapeutil::LaxPolyline(vertices));
  }
  return index;
}

S2Error::Code INDEXES_DO_NOT_MATCH = S2Error::USER_DEFINED_START;

class IndexMatchingLayer : public S2Builder::Layer {
 public:
  explicit IndexMatchingLayer(S2ShapeIndex const& index) : index_(index) {
    // Don't do any preprocessing of the graph edges.
    graph_options_.set_edge_type(EdgeType::DIRECTED);
    graph_options_.set_degenerate_edges(DegenerateEdges::KEEP);
    graph_options_.set_duplicate_edges(DuplicateEdges::KEEP);
    graph_options_.set_sibling_pairs(SiblingPairs::KEEP);
  }

  GraphOptions graph_options() const override {
    return graph_options_;
  }

  void Build(Graph const& g, S2Error* error) override;

 private:
  using EdgeVector = vector<pair<S2Point, S2Point>>;
  static string ToString(EdgeVector const& edges);
  S2ShapeIndex const& index_;
  GraphOptions graph_options_;
};

string IndexMatchingLayer::ToString(EdgeVector const& edges) {
  string msg;
  for (auto const& edge : edges) {
    vector<S2Point> vertices = { edge.first, edge.second };
    msg += s2textformat::ToString(vertices);
    msg += "; ";
  }
  return msg;
}

void IndexMatchingLayer::Build(Graph const& g, S2Error* error) {
  vector<pair<S2Point, S2Point> > actual, expected;
  for (int e = 0; e < g.num_edges(); ++e) {
    Graph::Edge const& edge = g.edge(e);
    actual.push_back(make_pair(g.vertex(edge.first), g.vertex(edge.second)));
  }
  for (int s = 0; s < index_.num_shape_ids(); ++s) {
    S2Shape* shape = index_.shape(s);
    if (shape == nullptr) continue;
    for (int e = shape->num_edges(); --e >= 0; ) {
      S2Point const *v0, *v1;
      shape->GetEdge(e, &v0, &v1);
      expected.push_back(make_pair(*v0, *v1));
    }
  }
  std::sort(actual.begin(), actual.end());
  std::sort(expected.begin(), expected.end());

  // The edges are a multiset, so we can't use std::set_difference.
  vector<pair<S2Point, S2Point>> missing, extra;
  for (auto ai = actual.begin(), ei = expected.begin();
       ai != actual.end() || ei != expected.end(); ) {
    if (ei == expected.end() || (ai != actual.end() && *ai < *ei)) {
      extra.push_back(*ai++);
    } else if (ai == actual.end() || *ei < *ai) {
      missing.push_back(*ei++);
    } else {
      ++ai;
      ++ei;
    }
  }
  if (!missing.empty() || !extra.empty()) {
    error->Init(INDEXES_DO_NOT_MATCH,
                "Missing edges: %s\nExtra edges: %s\n",
                ToString(missing).c_str(), ToString(extra).c_str());
  }
}

void ExpectResult(S2BoundaryOperation::OpType op_type ,
                  S2BoundaryOperation::Options const& options,
                  S2ShapeIndex const& a, S2ShapeIndex const& b,
                  S2ShapeIndex const& expected) {
  S2BoundaryOperation op(
      op_type, absl::MakeUnique<IndexMatchingLayer>(expected), options);
  S2Error error;
  EXPECT_TRUE(op.Build(a, b, &error))
      << S2BoundaryOperation::OpTypeToString(op_type) << " failed:\n"
      << error.text();

  // Now try the same thing with boolean output.
  bool result_empty;
  S2BoundaryOperation op2(op_type, &result_empty, options);
  EXPECT_TRUE(op2.Build(a, b, &error)) << "Boolean "
      << S2BoundaryOperation::OpTypeToString(op_type) << " failed:\n"
      << error.text();
  EXPECT_EQ(expected.num_shape_ids() == 0, result_empty);
}

}  // namespace

// The intersections in the "expected" data below were computed in lat-lng
// space (i.e., the rectangular projection), while the actual intersections
// are computed using geodesics.  We can compensate for this by rounding the
// intersection points to a fixed precision in degrees (e.g., 2 decimals).
static S2BoundaryOperation::Options RoundToE(int exp) {
  S2BoundaryOperation::Options options;
  options.set_snap_function(s2builderutil::IntLatLngSnapFunction(exp));
  return options;
}

TEST(S2BoundaryOperation, TwoDisjointTriangles) {
  S2BoundaryOperation::Options options;
  auto a = MakeIndex({"4:2, 3:1, 3:3"});
  auto b = MakeIndex({"2:0, 0:0, 2:2"});
  ExpectResult(OpType::UNION, options, *a, *b,
               *MakeIndex({"4:2, 3:1, 3:3", "2:0, 0:0, 2:2"}));
  ExpectResult(OpType::INTERSECTION, options, *a, *b,
               *MakeIndex({}));
  ExpectResult(OpType::DIFFERENCE, options, *a, *b,
               *MakeIndex({"4:2, 3:1, 3:3"}));
  ExpectResult(OpType::SYMMETRIC_DIFFERENCE, options, *a, *b,
               *MakeIndex({"4:2, 3:1, 3:3", "2:0, 0:0, 2:2"}));
}

TEST(S2BoundaryOperation, TwoTrianglesSharingEdge) {
  S2BoundaryOperation::Options options;
  auto a = MakeIndex({"4:2, 3:1, 3:3"});
  auto b = MakeIndex({"3:1, 2:2, 3:3"});
  ExpectResult(OpType::UNION, options, *a, *b,
               *MakeIndex({"4:2, 3:1, 2:2, 3:3"}));
  ExpectResult(OpType::INTERSECTION, options, *a, *b,
               *MakeIndex({}));
  ExpectResult(OpType::DIFFERENCE, options, *a, *b,
               *MakeIndex({"4:2, 3:1, 3:3"}));

  // For symmetric difference, the raw graph output should have all 6 edges.
  // (Sibling edge pairs can be removed by an S2Builder layer if desired.)
  ExpectResult(OpType::SYMMETRIC_DIFFERENCE, options, *a, *b,
               *MakeIndex({"4:2, 3:1, 3:3", "3:1, 2:2, 3:3"}));
}

TEST(S2BoundaryOperation, ThreeOverlappingBars) {
  // Two vertical bars and a horizontal bar that overlaps both of the other
  // bars and connects them.

  // Round intersection points to E2 precision because the expected results
  // were computed in lat/lng space rather than using geodesics.
  S2BoundaryOperation::Options options = RoundToE(2);
  auto a = MakeIndex({"0:0, 0:2, 3:2, 3:0", "0:3, 0:5, 3:5, 3:3"});
  auto b = MakeIndex({"1:1, 1:4, 2:4, 2:1"});
  ExpectResult(OpType::UNION, options, *a, *b, *MakeIndex(
    {"0:0, 0:2, 1:2, 1:3, 0:3, 0:5, 3:5, 3:3, 2:3, 2:2, 3:2, 3:0"}));
  ExpectResult(OpType::INTERSECTION, options, *a, *b, *MakeIndex(
    {"1:1, 1:2, 2:2, 2:1", "1:3, 1:4, 2:4, 2:3"}));
  ExpectResult(OpType::DIFFERENCE, options, *a, *b, *MakeIndex(
    {"0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0",
     "0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3"}));
  ExpectResult(OpType::SYMMETRIC_DIFFERENCE, options, *a, *b, *MakeIndex(
    {"0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0",
     "0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3",
     "1:2, 1:3, 2:3, 2:2"}));
}

TEST(S2BoundaryOperation, FourOverlappingBars) {
  // Two vertical bars and two horizontal bars.

  // Round intersection points to E2 precision because the expected results
  // were computed in lat/lng space rather than using geodesics.
  S2BoundaryOperation::Options options = RoundToE(2);
  auto a = MakeIndex({"1:88, 1:93, 2:93, 2:88", "-1:88, -1:93, 0:93, 0:88"});
  auto b = MakeIndex({"-2:89, -2:90, 3:90, 3:89", "-2:91, -2:92, 3:92, 3:91"});
  ExpectResult(OpType::UNION, options, *a, *b, *MakeIndex(
    {"-1:88, -1:89, -2:89, -2:90, -1:90, -1:91, -2:91, -2:92, -1:92, "
     "-1:93, 0:93, 0:92, 1:92, 1:93, 2:93, 2:92, 3:92, 3:91, 2:91, "
     "2:90, 3:90, 3:89, 2:89, 2:88, 1:88, 1:89, 0:89, 0:88",
     "0:90, 1:90, 1:91, 0:91" /*CW*/ }));
  ExpectResult(OpType::INTERSECTION, options, *a, *b, *MakeIndex(
    {"1:89, 1:90, 2:90, 2:89", "1:91, 1:92, 2:92, 2:91",
     "-1:89, -1:90, 0:90, 0:89", "-1:91, -1:92, 0:92, 0:91"}));
  ExpectResult(OpType::DIFFERENCE, options, *a, *b, *MakeIndex(
    {"1:88, 1:89, 2:89, 2:88", "1:90, 1:91, 2:91, 2:90",
     "1:92, 1:93, 2:93, 2:92", "-1:88, -1:89, 0:89, 0:88",
     "-1:90, -1:91, 0:91, 0:90", "-1:92, -1:93, 0:93, 0:92"}));
  ExpectResult(OpType::SYMMETRIC_DIFFERENCE, options, *a, *b, *MakeIndex(
    {"1:88, 1:89, 2:89, 2:88", "-1:88, -1:89, 0:89, 0:88",
     "1:90, 1:91, 2:91, 2:90", "-1:90, -1:91, 0:91, 0:90",
     "1:92, 1:93, 2:93, 2:92", "-1:92, -1:93, 0:93, 0:92",
     "-2:89, -2:90, -1:90, -1:89", "-2:91, -2:92, -1:92, -1:91",
     "0:89, 0:90, 1:90, 1:89", "0:91, 0:92, 1:92, 1:91",
     "2:89, 2:90, 3:90, 3:89", "2:91, 2:92, 3:92, 3:91"}));
}

TEST(S2BoundaryOperation, OverlappingDoughnuts) {
  // Two overlapping square doughnuts whose holes do not overlap.
  // This means that the union polygon has only two holes rather than three.

  // Round intersection points to E2 precision because the expected results
  // were computed in lat/lng space rather than using geodesics.
  S2BoundaryOperation::Options options = RoundToE(1);
  auto a = MakeIndex({"-1:-93, -1:-89, 3:-89, 3:-93",
                      "0:-92, 2:-92, 2:-90, 0:-90" /*CW*/ });
  auto b = MakeIndex({"-3:-91, -3:-87, 1:-87, 1:-91",
                      "-2:-90, 0:-90, 0:-88, -2:-88" /*CW*/ });
  ExpectResult(OpType::UNION, options, *a, *b, *MakeIndex(
    {"-1:-93, -1:-91, -3:-91, -3:-87, 1:-87, 1:-89, 3:-89, 3:-93",
     "0:-92, 2:-92, 2:-90, 1:-90, 1:-91, 0:-91" /*CW */,
     "-2:-90, -1:-90, -1:-89, 0:-89, 0:-88, -2:-88" /* CW */ }));
  ExpectResult(OpType::INTERSECTION, options, *a, *b, *MakeIndex(
    {"-1:-91, -1:-90, 0:-90, 0:-91",
     "0:-90, 0:-89, 1:-89, 1:-90"}));
  ExpectResult(OpType::DIFFERENCE, options, *a, *b, *MakeIndex(
    {"-1:-93, -1:-91, 0:-91, 0:-92, 2:-92, 2:-90, 1:-90, 1:-89, 3:-89, 3:-93",
     "-1:-90, -1:-89, 0:-89, 0:-90"}));
  ExpectResult(OpType::SYMMETRIC_DIFFERENCE, options, *a, *b, *MakeIndex(
    {"-1:-93, -1:-91, 0:-91, 0:-92, 2:-92, 2:-90, 1:-90, 1:-89, 3:-89, 3:-93",
     "-3:-91, -3:-87, 1:-87, 1:-89, 0:-89, 0:-88,-2:-88,-2:-90,-1:-90,-1:-91",
     "-1:-90, -1:-89, 0:-89, 0:-90",
     "1:-91, 0:-91, 0:-90, 1:-90"}));
}

TEST(S2BoundaryOperation, PolylineOverlappingRectangle) {
  // A polyline that crosses from the outside to the inside of a rectangle at
  // one of its vertices.
  S2BoundaryOperation::Options options = RoundToE(1);
  auto a = MakeIndex({}, {"0:0, 2:2"});
  auto b = MakeIndex({"1:1, 1:3, 3:3, 3:1"});
  ExpectResult(OpType::UNION, options, *a, *b, *MakeIndex(
    {"1:1, 1:3, 3:3, 3:1"}, {"0:0, 1:1"}));
  ExpectResult(OpType::INTERSECTION, options, *a, *b, *MakeIndex(
    {}, {"1:1, 2:2"}));
  ExpectResult(OpType::DIFFERENCE, options, *a, *b, *MakeIndex(
    {}, {"0:0, 1:1"}));
  ExpectResult(OpType::SYMMETRIC_DIFFERENCE, options, *a, *b, *MakeIndex(
    {"1:1, 1:3, 3:3, 3:1"}, {"0:0, 1:1"}));
}

TEST(S2BoundaryOperation, PolylineFollowingRectangle) {
  // A polyline that starts outside a rectangle, follows one of its edges,
  // and continues outside the rectangle.
  S2BoundaryOperation::Options options = RoundToE(1);
  auto a = MakeIndex({}, {"0:0, 1:1, 1:3, 0:4"});
  auto b = MakeIndex({"1:1, 1:3, 3:3, 3:1"});
#if 0
  // Polyline/polygon union operations are not supported yet.
  ExpectResult(OpType::UNION, options, *a, *b, *MakeIndex(
    {"1:1, 1:3, 3:3, 3:1"}, {"0:0, 1:1", "1:3, 0:4"}));
#endif
  ExpectResult(OpType::INTERSECTION, options, *a, *b, *MakeIndex(
    {}, {"1:1, 1:3"}));
  ExpectResult(OpType::DIFFERENCE, options, *a, *b, *MakeIndex(
    {}, {"0:0, 1:1", "1:3, 0:4"}));
  ExpectResult(OpType::SYMMETRIC_DIFFERENCE, options, *a, *b, *MakeIndex(
    {"1:1, 1:3, 3:3, 3:1"}, {"0:0, 1:1", "1:3, 0:4"}));
}

TEST(S2BoundaryOperation, PolylineCrossingRectangleTwice) {
  // A polyline that crosses a rectangle in one direction, then moves to a
  // different side and crosses the rectangle in the other direction.  Note
  // that an extra vertex is added where the two polyline edges cross.
  S2BoundaryOperation::Options options = RoundToE(1);
  auto a = MakeIndex({}, {"0:-5, 0:5, 5:0, -5:0"});
  auto b = MakeIndex({"1:1, 1:-1, -1:-1, -1:1"});
  ExpectResult(OpType::UNION, options, *a, *b, *MakeIndex(
      {"1:1, 1:0, 1:-1, 0:-1, -1:-1, -1:0, -1:1, 0:1"},
      {"0:-5, 0:-1", "0:1, 0:5, 5:0, 1:0", "-1:0, -5:0"}));
  ExpectResult(OpType::INTERSECTION, options, *a, *b, *MakeIndex(
      {}, {"0:-1, 0:0, 0:1", "1:0, 0:0, -1:0"}));
  ExpectResult(OpType::DIFFERENCE, options, *a, *b, *MakeIndex(
      {}, {"0:-5, 0:-1", "0:1, 0:5, 5:0, 1:0", "-1:0, -5:0"}));
  ExpectResult(OpType::SYMMETRIC_DIFFERENCE, options, *a, *b, *MakeIndex(
      {"1:1, 1:0, 1:-1, 0:-1, -1:-1, -1:0, -1:1, 0:1"},
      {"0:-5, 0:-1", "0:1, 0:5, 5:0, 1:0", "-1:0, -5:0"}));
}

#if 0
// Polyline/polyline operations are not supported yet.
TEST(S2BoundaryOperation, IntersectingPolylines) {
  S2BoundaryOperation::Options options = RoundToE(1);
  auto a = MakeIndex({}, {"0:0, 2:2"});
  auto b = MakeIndex({}, {"2:0, 0:2"});
  ExpectResult(OpType::UNION, options, *a, *b, *MakeIndex(
    {}, {"0:0, 1:1, 2:2", "2:0, 1:1, 0:2"}));
  ExpectResult(OpType::INTERSECTION, options, *a, *b, *MakeIndex(
    {}, {}));
  ExpectResult(OpType::DIFFERENCE, options, *a, *b, *MakeIndex(
    {}, {"0:0, 2:2"}));
  ExpectResult(OpType::SYMMETRIC_DIFFERENCE, options, *a, *b, *MakeIndex(
    {}, {"0:0, 1:1, 2:2", "2:0, 1:1, 0:2"}));
}

TEST(S2BoundaryOperation, OverlappingPolylines) {
  S2BoundaryOperation::Options options = RoundToE(1);
  auto a = MakeIndex({}, {"0:0, 0:2, 2:2"});
  auto b = MakeIndex({}, {"0:2, 2:2, 2:4"});
  ExpectResult(OpType::UNION, options, *a, *b, *MakeIndex(
    {}, {"0:0, 0:2, 2:2, 2:4"}));
  ExpectResult(OpType::INTERSECTION, options, *a, *b, *MakeIndex(
    {}, {"0:2, 2:2"}));
  ExpectResult(OpType::DIFFERENCE, options, *a, *b, *MakeIndex(
    {}, {"0:0, 0:2"}));
  ExpectResult(OpType::SYMMETRIC_DIFFERENCE, options, *a, *b, *MakeIndex(
    {}, {"0:0, 0:2", "2:2, 2:4"}));
}
#endif
