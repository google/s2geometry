// Copyright 2016 Google Inc. All Rights Reserved.
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

#include "s2/s2builder.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <gflags/gflags.h>
#include <glog/log_severity.h>
#include "s2/base/timer.h"
#include "s2/strings/join.h"
#include <gtest/gtest.h>
#include "s2/third_party/absl/memory/memory.h"
#include "s2/s2builder_layer.h"
#include "s2/s2builderutil_layers.h"
#include "s2/s2builderutil_snap_functions.h"
#include "s2/s2cap.h"
#include "s2/s2cellid.h"
#include "s2/s2debug.h"
#include "s2/s2latlng.h"
#include "s2/s2loop.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"
#include "s2/s2predicates.h"
#include "s2/s2testing.h"
#include "s2/s2textformat.h"

using gtl::MakeUnique;
using std::cout;
using std::endl;
using std::min;
using std::unique_ptr;
using std::vector;
using s2builderutil::IdentitySnapFunction;
using s2builderutil::S2CellIdSnapFunction;
using s2builderutil::IntLatLngSnapFunction;
using s2builderutil::S2PolygonLayer;
using s2builderutil::S2PolylineLayer;
using s2builderutil::S2PolylineVectorLayer;
using s2textformat::MakePolygon;
using s2textformat::MakePolyline;
using EdgeType = S2Builder::EdgeType;
using GraphOptions = S2Builder::GraphOptions;;

DEFINE_int32(iteration_multiplier, 1,
             "Iteration multiplier for randomized tests");

namespace {

void ExpectPolygonsEqual(S2Polygon const& expected,
                         S2Polygon const& actual) {
  EXPECT_TRUE(expected.Equals(&actual))
      << "\nExpected:\n" << s2textformat::ToString(expected)
      << "\nActual:\n" << s2textformat::ToString(actual);
}

void ExpectPolygonsApproxEqual(S2Polygon const& expected,
                               S2Polygon const& actual,
                               S1Angle tolerance) {
  EXPECT_TRUE(expected.BoundaryApproxEquals(actual, tolerance))
      << "\nExpected:  " << s2textformat::ToString(expected)
      << "\nActual:    " << s2textformat::ToString(actual)
      << "\nTolerance: " << tolerance.degrees();
}

void ExpectPolylinesEqual(S2Polyline const& expected,
                          S2Polyline const& actual) {
  EXPECT_TRUE(expected.Equals(&actual))
      << "\nExpected:\n" << s2textformat::ToString(expected)
      << "\nActual:\n" << s2textformat::ToString(actual);
}

TEST(S2Builder, SimpleVertexMerging) {
  // When IdentitySnapFunction is used (i.e., no special requirements on
  // vertex locations), check that vertices closer together than the snap
  // radius are merged together.

  S1Angle snap_radius = S1Angle::Degrees(0.5);
  S2Builder builder((S2Builder::Options(IdentitySnapFunction(snap_radius))));
  S2Polygon output;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output));
  unique_ptr<S2Polygon> input(MakePolygon(
      "0:0, 0.2:0.2, 0.1:0.2, 0.1:0.9, 0:1, 0.1:1.1, 0.9:1, 1:1, 1:0.9"));
  builder.AddPolygon(*input);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  unique_ptr<S2Polygon> expected(MakePolygon("0:0, 0:1, 1:0.9"));
  ExpectPolygonsApproxEqual(*expected, output, snap_radius);
}

TEST(S2Builder, SimpleS2CellIdSnapping) {
  // When S2CellIdSnapFunction is used, check that all output vertices are the
  // centers of S2CellIds at the specified level level.

  int level = S2CellIdSnapFunction::LevelForMaxSnapRadius(S1Angle::Degrees(1));
  S2CellIdSnapFunction snap_function(level);
  S2Builder builder((S2Builder::Options(snap_function)));
  S2Polygon output;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output));
  unique_ptr<S2Polygon> input(MakePolygon(
      "2:2, 3:4, 2:6, 4:5, 6:6, 5:4, 6:2, 4:3"));
  builder.AddPolygon(*input);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  ASSERT_EQ(1, output.num_loops());
  S2Loop const* loop = output.loop(0);
  for (int i = 0; i < loop->num_vertices(); ++i) {
    EXPECT_EQ(S2CellId(loop->vertex(i)).parent(level).ToPoint(),
              loop->vertex(i));
  }
  ExpectPolygonsApproxEqual(*input, output, snap_function.snap_radius());
}

TEST(S2Builder, SimpleIntLatLngSnapping) {
  S2Builder builder(S2Builder::Options(IntLatLngSnapFunction(0)));  // E0 coords
  S2Polygon output;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output));
  unique_ptr<S2Polygon> input(MakePolygon(
      "2.01:2.09, 3.24:4.49, 1.78:6.25, 3.51:5.49, 6.11:6.11, "
      "5.22:3.88, 5.55:2.49, 4.49:2.51"));
  unique_ptr<S2Polygon> expected(MakePolygon(
      "2:2, 3:4, 2:6, 4:5, 6:6, 5:4, 6:2, 4:3"));
  builder.AddPolygon(*input);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  ASSERT_EQ(1, output.num_loops());
  ExpectPolygonsEqual(*expected, output);
}

TEST(S2Builder, VerticesMoveLessThanSnapRadius) {
  // Check that chains of closely spaced vertices do not collapse into a
  // single vertex.

  S1Angle snap_radius = S1Angle::Degrees(1);
  S2Builder builder((S2Builder::Options(IdentitySnapFunction(snap_radius))));
  S2Polygon output;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output));
  // The spacing between input vertices is about 2*pi*20/1000 = 0.125 degrees.
  // The output vertices are spaced between 1 and 2 degrees apart; the average
  // spacing is about 1.33 degrees.
  S2Polygon input(
      S2Loop::MakeRegularLoop(S2Point(1, 0, 0), S1Angle::Degrees(20), 1000));
  builder.AddPolygon(input);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  ASSERT_EQ(1, output.num_loops());
  EXPECT_GE(output.loop(0)->num_vertices(), 90);
  EXPECT_LE(output.loop(0)->num_vertices(), 100);
  EXPECT_TRUE(output.BoundaryNear(input, snap_radius));
}

TEST(S2Builder, MinEdgeVertexSeparation) {
  // Check that edges are separted from non-incident vertices by at least
  // min_edge_vertex_separation().  This requires adding new vertices (not
  // present in the input) in some cases.

  // The input is a skinny right triangle with two legs of length 10 and 1,
  // and whose diagonal is subdivided into 10 short edges.  Using a snap
  // radius of 0.5, about half of the long leg is snapped onto the diagonal
  // (which causes that part of the polygon to be removed).  But the real
  // problem is that the remaining part of the long leg gets too close to the
  // remaining vertices on the diagonal, i.e. it would violate the minimum
  // edge-vertex separation guarantee.  S2Builder handles this by creating at
  // least one vertex along the original long leg, to keep the snapped edge
  // far enough away from the diagonal.
  unique_ptr<S2Polygon> input(MakePolygon(
      "0:0, 0:1, 1:.9, 2:.8, 3:.7, 4:.6, 5:.5, 6:.4, 7:.3, 8:.2, 9:.1, 10:0"));
  unique_ptr<S2Polygon> expected(MakePolygon(
      "0:0, 0:1, 1:.9, 2:.8, 3:.7, 4:.6, 5:.5, 4.00021862252687:0"));
  S2Builder::Options options(IdentitySnapFunction(S1Angle::Degrees(0.5)));
  S2Builder builder(options);
  S2Polygon output;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output));
  builder.AddPolygon(*input);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  ExpectPolygonsApproxEqual(*expected, output, S1Angle::Radians(1e-15));
}

TEST(S2Builder, IdempotencySnapsUnsnappedVertices) {
  // When idempotency is requested, no snapping is done unless S2Builder finds
  // at least one vertex or edge that could not be the output of a previous
  // snapping operation.  This test checks that S2Builder detects vertices
  // that are not at a valid location returned by the given snap function.

  // In this example we snap two vertices to integer lat/lng coordinates.  The
  // two vertices are far enough apart (more than min_vertex_separation) so
  // that they might be the result of a previous snapping operation, but one
  // of the two vertices does not have integer lat/lng coordinates.  We use
  // internal knowledge of how snap sites are chosen (namely, that candidates
  // are considered in S2CellId order) to construct two different cases, one
  // where the snapped vertex is processed first and one where the unsnapped
  // vertex is processed first.  This exercises two different code paths.
  IntLatLngSnapFunction snap_function(0);
  EXPECT_GE(snap_function.snap_radius(), S1Angle::Degrees(0.7));
  EXPECT_LE(snap_function.min_vertex_separation(), S1Angle::Degrees(0.35));
  S2Builder builder((S2Builder::Options(snap_function)));

  // In this example, the snapped vertex (0, 0) is processed first and is
  // selected as a Voronoi site (i.e., output vertex).  The second vertex is
  // closer than snap_radius(), therefore it is snapped to the first vertex
  // and the polyline becomes degenerate.
  S2Point a = S2LatLng::FromDegrees(0, 0).ToPoint();
  S2Point b = S2LatLng::FromDegrees(0.01, 0.6).ToPoint();
  EXPECT_LT(S2CellId(a), S2CellId(b));
  S2Polyline input1(vector<S2Point>{a, b}), output1;
  builder.StartLayer(MakeUnique<S2PolylineLayer>(&output1));
  builder.AddPolyline(input1);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error));
  EXPECT_EQ("", s2textformat::ToString(output1));

  // In this example the unsnapped vertex is processed first and is snapped to
  // (0, 0).  The second vertex is further than snap_radius() away, so it is
  // also snapped (which does nothing) and is left at (0, 1).
  S2Point c = S2LatLng::FromDegrees(0.01, 0.4).ToPoint();
  S2Point d = S2LatLng::FromDegrees(0, 1).ToPoint();
  EXPECT_LT(S2CellId(c), S2CellId(d));
  S2Polyline input2(vector<S2Point>{c, d}), output2;
  builder.StartLayer(MakeUnique<S2PolylineLayer>(&output2));
  builder.AddPolyline(input2);
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  EXPECT_EQ("0:0, 0:1", s2textformat::ToString(output2));
}

TEST(S2Builder, IdempotencySnapsEdgesWithTinySnapRadius) {
  // When idempotency is requested, no snapping is done unless S2Builder finds
  // at least one vertex or edge that could not be the output of a previous
  // snapping operation.  This test checks that S2Builder detects edges that
  // are too close to vertices even when the snap radius is very small
  // (e.g., S2EdgeUtil::kIntersectionError).
  //
  // Previously S2Builder used a conservative approximation to decide whether
  // edges were too close to vertices; unfortunately this meant that when the
  // snap radius was very small then no snapping would be done at all, because
  // even an edge/vertex distance of zero was considered far enough apart.
  //
  // This tests that the current code (which uses exact predicates) handles
  // this situation correctly (i.e., that an edge separated from a
  // non-incident vertex by a distance of zero cannot be the output of a
  // previous snapping operation).
  S2Builder::Options options;
  options.set_snap_function(
      s2builderutil::IdentitySnapFunction(S2EdgeUtil::kIntersectionError));
  S2PolylineVectorLayer::Options layer_options;
  layer_options.set_duplicate_edges(
      S2PolylineVectorLayer::Options::DuplicateEdges::MERGE);
  S2Builder builder(options);
  vector<unique_ptr<S2Polyline>> output;
  builder.StartLayer(MakeUnique<S2PolylineVectorLayer>(&output, layer_options));
  builder.AddPolyline(*s2textformat::MakePolyline("0:0, 0:10"));
  builder.AddPolyline(*s2textformat::MakePolyline("0:5, 0:7"));
  S2Error error;
  ASSERT_TRUE(builder.Build(&error));
  ASSERT_EQ(1, output.size());
  EXPECT_EQ("0:0, 0:5, 0:7, 0:10", s2textformat::ToString(*output[0]));
}

TEST(S2Builder, IdempotencyDoesNotSnapAdequatelySeparatedEdges) {
  // When idempotency is requested, no snapping is done unless S2Builder finds
  // at least one vertex or edge that could not be the output of a previous
  // snapping operation.  This test checks that when an edge is further away
  // than min_edge_vertex_separation() then no snapping is done.
  S2Builder::Options options(IntLatLngSnapFunction(0));
  options.set_idempotent(true);  // Test fails if this is "false".
  S2Builder builder(options);
  S2Polygon output1, output2;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output1));
  builder.AddPolygon(*MakePolygon("1.49:0, 0:2, 0.49:3"));
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  const char* expected = "1:0, 0:2, 0:3";
  EXPECT_EQ(expected, s2textformat::ToString(output1));
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output2));
  builder.AddPolygon(output1);
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  EXPECT_EQ(expected, s2textformat::ToString(output2));
}

TEST(S2Builder, kMaxSnapRadiusCanSnapAtLevel0) {
  // Verify that kMaxSnapRadius will allow snapping at S2CellId level 0.
  EXPECT_LE(S2CellIdSnapFunction::MinSnapRadiusForLevel(0),
            S2Builder::SnapFunction::kMaxSnapRadius());
}

TEST(S2Builder, S2CellIdSnappingAtAllLevels) {
  unique_ptr<const S2Polygon> input(MakePolygon(
      "0:0, 0:2, 2:0; 0:0, 0:-2, -2:-2, -2:0"));
  for (int level = 0; level <= S2CellId::kMaxLevel; ++level) {
    S2CellIdSnapFunction snap_function(level);
    S2Builder builder((S2Builder::Options(snap_function)));
    S2Polygon output;
    builder.StartLayer(MakeUnique<S2PolygonLayer>(&output));
    builder.AddPolygon(*input);
    S2Error error;
    ASSERT_TRUE(builder.Build(&error)) << error.text();
    EXPECT_TRUE(output.IsValid());
    // The ApproxContains calls below are not guaranteed to succeed in general
    // because ApproxContains works by snapping both polygons together using
    // the given tolerance and then checking for containment.  Since
    // ApproxContains snaps to an arbitrary subset of the input vertices
    // rather than to S2CellId centers at the current level, this means that
    // corresponding vertices in "input" and "output" can snap to different
    // sites, which causes the containment test to fail.  Nevertheless, by
    // using a larger tolerance of 2 * snap_radius, all calls in this test
    // succeed (and would be likely to succeed in other similar tests).
    // (To guarantee correctness we would need to use S2CellIdSnapFunction
    // within the ApproxContains implementation.)
    S1Angle tolerance = min(2 * snap_function.snap_radius(),
                            snap_function.kMaxSnapRadius());
    EXPECT_TRUE(output.ApproxContains(input.get(), tolerance));
    EXPECT_TRUE(input->ApproxContains(&output, tolerance));
  }
}

TEST(S2Builder, SnappingDoesNotRotateVertices) {
  // This is already tested extensively elsewhere.
  unique_ptr<S2Polygon> input(MakePolygon(
      "49.9305505:-124.8345463, 49.9307448:-124.8299657, "
      "49.9332101:-124.8301996, 49.9331224:-124.8341368; "
      "49.9311087:-124.8327042, 49.9318176:-124.8312621, "
      "49.9318866:-124.8334451"));
  S2Builder::Options options((S2CellIdSnapFunction()));
  S2Builder builder(options);
  S2Polygon output1, output2;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output1));
  builder.AddPolygon(*input);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  // This checks that the vertices are in the same cyclic order, and that
  // vertices have not moved by more than "snap_radius".
  ExpectPolygonsApproxEqual(*input, output1,
                            options.snap_function().snap_radius());

  // Check that snapping twice doesn't rotate the vertices.  This also
  // verifies that S2Builder can be used again after Build() is called.
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output2));
  builder.AddPolygon(output1);
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  ExpectPolygonsEqual(output1, output2);
}

TEST(S2Builder, SelfIntersectingPolyline) {
  // Check that when two edges of a polyline cross, the intersection point is
  // added to both edges.

  S2Builder::Options options;
  IntLatLngSnapFunction snap_function(1);  // Snap to E1 coordinates
  options.set_snap_function(snap_function);
  options.set_split_crossing_edges(true);
  S2Builder builder(options);
  S2Polyline output;
  builder.StartLayer(MakeUnique<S2PolylineLayer>(&output));
  unique_ptr<S2Polyline> input(MakePolyline("3:1, 1:3, 1:1, 3:3"));
  unique_ptr<S2Polyline> expected(MakePolyline("3:1, 2:2, 1:3, 1:1, 2:2, 3:3"));
  builder.AddPolyline(*input);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  ExpectPolylinesEqual(*expected, output);
}

TEST(S2Builder, SelfIntersectingPolygon) {
  // Check that when two edge of a polygon cross, the intersection point is
  // added to both edges, and that the resulting (undirected) edges can be
  // assembled into a valid polygon.

  IntLatLngSnapFunction snap_function(1);  // Snap to E1 coordinates
  S2Builder::Options options;
  options.set_snap_function(snap_function);
  options.set_split_crossing_edges(true);
  S2Builder builder(options);
  S2Polygon output;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(
      &output, S2PolygonLayer::Options(EdgeType::UNDIRECTED)));
  unique_ptr<S2Polyline> input(MakePolyline("3:1, 1:3, 1:1, 3:3, 3:1"));
  unique_ptr<S2Polygon> expected(MakePolygon("1:1, 1:3, 2:2; 3:1, 2:2, 3:3"));
  builder.AddPolyline(*input);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  ExpectPolygonsEqual(*expected, output);
}

TEST(S2Builder, TieBreakingIsConsistent) {
  // Check that when an edge passes between two equally distant vertices, that
  // the choice of which one to snap to does not depend on the edge direction.

  S2Builder::Options options(IdentitySnapFunction(S1Angle::Degrees(2)));
  options.set_idempotent(false);
  S2Builder builder(options);
  builder.ForceVertex(S2LatLng::FromDegrees(1, 0).ToPoint());
  builder.ForceVertex(S2LatLng::FromDegrees(-1, 0).ToPoint());
  S2Polyline output1, output2;
  builder.StartLayer(gtl::MakeUnique<S2PolylineLayer>(&output1));
  builder.AddPolyline(*MakePolyline("0:-5, 0:5"));
  builder.StartLayer(gtl::MakeUnique<S2PolylineLayer>(&output2));
  builder.AddPolyline(*MakePolyline("0:5, 0:-5"));
  S2Error error;
  EXPECT_TRUE(builder.Build(&error)) << error.text();
  EXPECT_EQ(3, output1.num_vertices());
  EXPECT_EQ(3, output2.num_vertices());
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(output1.vertex(i), output2.vertex(2 - i));
  }
}

void TestPolylineLayers(
    vector<char const*> const& input_strs,
    vector<char const*> const& expected_strs,
    S2PolylineLayer::Options const& layer_options,
    S2Builder::Options const& builder_options = S2Builder::Options()) {
  SCOPED_TRACE(layer_options.edge_type() == EdgeType::DIRECTED ?
               "DIRECTED" : "UNDIRECTED");
  S2Builder builder(builder_options);
  vector<unique_ptr<S2Polyline>> output;
  for (auto input_str : input_strs) {
    output.emplace_back(new S2Polyline);
    builder.StartLayer(MakeUnique<S2PolylineLayer>(output.back().get(),
                                                   layer_options));
    builder.AddPolyline(*MakePolyline(input_str));
  }
  S2Error error;
  ASSERT_TRUE(builder.Build(&error));
  vector<string> output_strs;
  for (auto const& polyline : output) {
    output_strs.push_back(s2textformat::ToString(*polyline));
  }
  EXPECT_EQ(strings::Join(expected_strs, "; "),
            strings::Join(output_strs, "; "));
}

void TestPolylineVector(
    vector<char const*> const& input_strs,
    vector<char const*> const& expected_strs,
    S2PolylineVectorLayer::Options const& layer_options,
    S2Builder::Options const& builder_options = S2Builder::Options()) {
  S2Builder builder(builder_options);
  vector<unique_ptr<S2Polyline>> output;
  builder.StartLayer(MakeUnique<S2PolylineVectorLayer>(&output, layer_options));
  for (auto input_str : input_strs) {
    builder.AddPolyline(*MakePolyline(input_str));
  }
  S2Error error;
  ASSERT_TRUE(builder.Build(&error));
  vector<string> output_strs;
  for (auto const& polyline : output) {
    output_strs.push_back(s2textformat::ToString(*polyline));
  }
  EXPECT_EQ(strings::Join(expected_strs, "; "),
            strings::Join(output_strs, "; "));
}

void TestPolylineLayersBothEdgeTypes(
    vector<char const*> const& input_strs,
    vector<char const*> const& expected_strs,
    S2PolylineLayer::Options layer_options,  // by value
    S2Builder::Options const& builder_options = S2Builder::Options()) {
  layer_options.set_edge_type(EdgeType::DIRECTED);
  TestPolylineLayers(input_strs, expected_strs, layer_options, builder_options);
  layer_options.set_edge_type(EdgeType::UNDIRECTED);
  TestPolylineLayers(input_strs, expected_strs, layer_options, builder_options);
}

TEST(S2Builder, SimplifyOneEdge) {
  // Simplify a perturbed edge chain into a single edge.

  S2Builder::Options options(IdentitySnapFunction(S1Angle::Degrees(1)));
  options.set_simplify_edge_chains(true);
  TestPolylineLayersBothEdgeTypes({"0:0, 1:0.5, 2:-0.5, 3:0.5, 4:-0.5, 5:0"},
                                  {"0:0, 5:0"},
                                  S2PolylineLayer::Options(), options);
}

TEST(S2Builder, SimplifyTwoLayers) {
  // Construct two layers, each containing a polyline that could be simplified
  // to a single edge on its own.  However the two polylines actually cross,
  // so make sure that the output still contains the intersection vertex.

  S2Builder::Options options(IdentitySnapFunction(S1Angle::Degrees(0.5)));
  options.set_split_crossing_edges(true);
  options.set_simplify_edge_chains(true);
  TestPolylineLayersBothEdgeTypes(
      {"-2:-1, -1:0, 1:0, 2:1", "1:-2, 0:-1, 0:1, -1:2"},
      {"-2:-1, 0:0, 2:1", "1:-2, 0:0, -1:2"},
      S2PolylineLayer::Options(), options);
}

TEST(S2Builder, SimplifyOneLoop) {
  // Simplify a regular loop with 1000 vertices and a radius of 20 degrees.
  // Turning on edge chain simplification yields a dramatically smaller number
  // of vertices than snapping alone (10 vertices vs 95 vertices using a snap
  // radius of 1 degree).  This is because snapping alone yields vertices that
  // stay within 1 degree of the input *vertices*, while simplifying edge
  // chains yields edges that stay within 1 degree of the input *edges*.

  for (int i = 0; i < 2; ++i) {
    EdgeType edge_type = static_cast<EdgeType>(i);
    S1Angle snap_radius = S1Angle::Degrees(1);
    S2Builder::Options options((IdentitySnapFunction(snap_radius)));
    options.set_simplify_edge_chains(true);
    S2Builder builder(options);
    S2Polygon output;
    builder.StartLayer(MakeUnique<S2PolygonLayer>(
        &output, S2PolygonLayer::Options(edge_type)));
    // Spacing between vertices: approximately 2*pi*20/1000 = 0.125 degrees.
    S2Polygon input(
        S2Loop::MakeRegularLoop(S2Point(1, 0, 0), S1Angle::Degrees(20), 1000));
    builder.AddPolygon(input);
    S2Error error;
    ASSERT_TRUE(builder.Build(&error)) << error.text();
    ASSERT_EQ(1, output.num_loops());
    EXPECT_GE(output.loop(0)->num_vertices(), 10);
    EXPECT_LE(output.loop(0)->num_vertices(), 12);
    EXPECT_TRUE(output.BoundaryNear(input, snap_radius));
  }
}

TEST(S2Builder, SimplifyOppositeDirections) {
  // We build two layers with two polylines that follow the same circular arc
  // in opposite directions, and verify that they are snapped identically.
  // (The snap radius is adjusted so that the arc is simplified into a long
  // edge and a short edge, and therefore we would get a different result if
  // the two layers followed the edge chain in different directions.)

  S2Builder::Options options(IdentitySnapFunction(S1Angle::Degrees(0.5)));
  options.set_simplify_edge_chains(true);
  TestPolylineLayersBothEdgeTypes(
      {"-4:0.83, -3:0.46, -2:0.2, -1:0.05, 0:0, 1:0.5, 2:0.2, 3:0.46, 4:0.83",
            "4:.83, 3:.46, 2:.2, 1:.05, 0:0, -1:.5, -2:.2, -3:.46, -4:.83"},
      {"-4:0.83, -2:0.2, 4:0.83", "4:0.83, -2:0.2, -4:0.83"},
      S2PolylineLayer::Options(), options);
}

TEST(S2Builder, SimplifyKeepsEdgeVertexSeparation) {
  // We build two layers each containing a polyline, such that the polyline in
  // the first layer could be simplified to a straight line except that then
  // it would create an intersection with the second polyline.

  S2Builder::Options options(IdentitySnapFunction(S1Angle::Degrees(1.0)));
  options.set_simplify_edge_chains(true);
  TestPolylineLayersBothEdgeTypes(
      {"0:-10, 0.99:0, 0:10", "-5:-5, -0.2:0, -5:5"},
      {"0:-10, 0.99:0, 0:10", "-5:-5, -0.2:0, -5:5"},
      S2PolylineLayer::Options(), options);
}

TEST(S2Builder, SimplifyBacktrackingEdgeChain) {
  // Test simplifying an edge chain that backtracks on itself.
  S2Builder::Options options(IdentitySnapFunction(S1Angle::Degrees(0.5)));
  options.set_simplify_edge_chains(true);
  TestPolylineLayersBothEdgeTypes(
      {"0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 4:0, 3:0, "
            "2:0, 3:0, 4:0, 5:0, 6:0, 7:0"},
      {"0:0, 2:0, 5:0, 2:0, 5:0, 7:0"},
      S2PolylineLayer::Options(), options);
}

TEST(S2Builder, SimplifyLimitsEdgeDeviation) {
  // Make sure that simplification does not create long edges such that the
  // midpoint of the edge might be further than max_edge_deviation() from an
  // input edge.  In the example below, vertices are snapped to integer
  // lat/lng coordinates, and the snap radius is approximately 0.707 degrees.
  // Snapping moves the input vertices perpendicular to the input edge by just
  // slightly less than the snap radius (0.693 degrees).  Now the midpoint of
  // the snapped edge is about 0.98 degrees from the input edge, which causes
  // an extra site to be added at the midpoint of the original edge.
  //
  // When simplify_edge_chains() is enabled, then usually an extra site like
  // this would be simplified away (because the simplified edge would still be
  // within snap_radius() of all the input vertices) except that there is an
  // explicit check in S2Builder that prevents this.  (If the check is removed
  // then this test fails.)

  S2Builder::Options options(IntLatLngSnapFunction(0));  // E0 coordinates
  options.set_simplify_edge_chains(true);
  TestPolylineLayersBothEdgeTypes(
      {"-30.49:-29.51, 29.51:30.49"}, {"-30:-30, -1:1, 30:30"},
      S2PolylineLayer::Options(), options);
}

TEST(S2Builder, SimplifyPreservesTopology) {
  // Crate several nested concentric loops, and verify that the loops are
  // still nested after simplification.

  int const kNumLoops = 20;
  int const kNumVerticesPerLoop = 1000;
  S1Angle const kBaseRadius = S1Angle::Degrees(5);
  S1Angle const kSnapRadius = S1Angle::Degrees(0.1);
  S2Builder::Options options((IdentitySnapFunction(kSnapRadius)));
  options.set_simplify_edge_chains(true);
  S2Builder builder(options);
  vector<unique_ptr<S2Polygon>> input, output;
  for (int j = 0; j < kNumLoops; ++j) {
    // Spacing between vertices: approximately 2*pi*20/1000 = 0.125 degrees.
    S1Angle radius = kBaseRadius + 0.7 * j * j / kNumLoops * kSnapRadius;
    input.emplace_back(new S2Polygon(S2Loop::MakeRegularLoop(
        S2Point(1, 0, 0), radius, kNumVerticesPerLoop)));
    output.emplace_back(new S2Polygon);
    builder.StartLayer(MakeUnique<S2PolygonLayer>(output.back().get()));
    builder.AddPolygon(*input.back());
  }
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  for (int j = 0; j < kNumLoops; ++j) {
    EXPECT_TRUE(output[j]->BoundaryNear(*input[j], kSnapRadius));
    if (j > 0) EXPECT_TRUE(output[j]->Contains(output[j - 1].get()));
  }
}

TEST(S2Builder, SimplifyRemovesSiblingPairs) {
  S2Builder::Options options(IntLatLngSnapFunction(0));  // E0 coords
  S2PolylineVectorLayer::Options layer_options;
  layer_options.set_sibling_pairs(GraphOptions::SiblingPairs::DISCARD);

  // Check that there is no sibling pair without simplification.
  TestPolylineVector(
      {"0:0, 0:10", "0:10, 0.6:5, 0:0"},
      {"0:0, 0:10, 1:5, 0:0"}, layer_options, options);

  // Now check that (1) simplification produces a sibling pair,
  // and (2) the sibling pair is removed (since we requested it).
  options.set_simplify_edge_chains(true);
  TestPolylineVector(
  {"0:0, 0:10", "0:10, 0.6:5, 0:0"},
  {}, layer_options, options);
}

TEST(S2Builder, SimplifyMergesDuplicateEdges) {
  S2Builder::Options options(IntLatLngSnapFunction(0));  // E0 coords
  S2PolylineVectorLayer::Options layer_options;
  layer_options.set_duplicate_edges(GraphOptions::DuplicateEdges::MERGE);

  // Check that there are no duplicate edges without simplification.
  TestPolylineVector(
      {"0:0, 0:10", "0:0, 0.6:5, 0:10"},
      {"0:0, 0:10", "0:0, 1:5, 0:10"}, layer_options, options);

  // Now check that (1) simplification produces a duplicate edge pair,
  // and (2) the duplicate pair is merged (since we requested it).
  options.set_simplify_edge_chains(true);
  TestPolylineVector(
      {"0:0, 0:10", "0:0, 0.6:5, 0:10"},
      {"0:0, 0:10"}, layer_options, options);
}

TEST(S2Builder, HighPrecisionPredicates) {
  // To produce correct output in this example, the algorithm needs fall back
  // to high precision predicates when the output of the normal predicates is
  // uncertain.
  vector<S2Point> vertices = {
    {-0.1053119128423491, -0.80522217121852213, 0.58354661852470235},
    {-0.10531192039134209, -0.80522217309706012, 0.58354661457019508},
    {-0.10531192039116592, -0.80522217309701472, 0.58354661457028933},
  };
  S2Polyline input(vertices);
  S1Angle snap_radius = S2EdgeUtil::kIntersectionMergeRadius;
  S2Builder::Options options((IdentitySnapFunction(snap_radius)));
  options.set_idempotent(false);
  S2Builder builder(options);
  S2Polyline output;
  builder.StartLayer(gtl::MakeUnique<S2PolylineLayer>(&output));
  builder.ForceVertex(S2Point(
      -0.10531192039134191, -0.80522217309705857, 0.58354661457019719));
  builder.AddPolyline(input);
  S2Error error;
  EXPECT_TRUE(builder.Build(&error)) << error.text();
}

// Chooses a random S2Point that is often near the intersection of one of the
// coodinates planes or coordinate axes with the unit sphere.  (It is possible
// to represent very small perturbations near such points.)
S2Point ChoosePoint() {
  S2Point x = S2Testing::RandomPoint();
  for (int i = 0; i < 3; ++i) {
    if (S2Testing::rnd.OneIn(3)) {
      x[i] *= pow(1e-50, S2Testing::rnd.RandDouble());
    }
  }
  return x.Normalize();
}

TEST(S2Builder, HighPrecisionStressTest) {
  // This test constructs many small, random inputs such that the output is
  // likely to be inconsistent unless high-precision predicates are used.

  S1Angle snap_radius = S2EdgeUtil::kIntersectionMergeRadius;
  // Some S2Builder calculations use an upper bound that takes into account
  // S1ChordAngle errors.  We sometimes try perturbing points by very close to
  // that distance in an attempt to expose errors.
  S1ChordAngle ca(snap_radius);
  S1Angle snap_radius_with_error = ca.PlusError(
      ca.GetS1AngleConstructorMaxError() +
      S2EdgeUtil::GetUpdateMinDistanceMaxError(ca)).ToAngle();

  auto& rnd = S2Testing::rnd;
  int non_degenerate = 0;
  int const kIters = 8000 * FLAGS_iteration_multiplier;
  for (int iter = 0; iter < kIters; ++iter) {
    // TODO(ericv): This test fails with a random seed of 96.  Change this
    // back to "iter + 1" once all the exact predicates are implemented.
    rnd.Reset(iter + 1);  // Easier to reproduce a specific case.

    // We construct a nearly degenerate triangle where one of the edges is
    // sometimes very short.  Then we add a forced vertex somewhere near the
    // shortest edge.  Then after snapping, we check that (1) the edges still
    // form a loop, and (2) if the loop is non-degenerate, then it has the
    // same orientation as the original triangle.
    //
    // v1 is located randomly.  (v0,v1) is the longest of the three edges.
    // v2 is located along (v0,v1) but is perturbed by up to 2 * snap_radius.
    S2Point v1 = ChoosePoint(), v0_dir = ChoosePoint();
    double d0 = pow(1e-16, rnd.RandDouble());
    S2Point v0 = S2EdgeUtil::InterpolateAtDistance(S1Angle::Radians(d0),
                                                   v1, v0_dir);
    double d2 = 0.5 * d0 * pow(1e-16, pow(rnd.RandDouble(), 2));
    S2Point v2 = S2EdgeUtil::InterpolateAtDistance(S1Angle::Radians(d2),
                                                   v1, v0_dir);
    v2 = S2Testing::SamplePoint(S2Cap(v2, 2 * snap_radius));
    // Vary the edge directions by randomly swapping v0 and v2.
    if (rnd.OneIn(2)) std::swap(v0, v2);

    // The forced vertex (v3) is either located near the (v1, v2) edge.
    // We perturb it either in a random direction from v1 or v2, or
    // perpendicular to (v1, v2) starting from an interior edge point.
    S1Angle d3 = rnd.OneIn(2) ? snap_radius : snap_radius_with_error;
    if (rnd.OneIn(3)) d3 = 1.5 * rnd.RandDouble() * d3;
    S2Point v3;
    if (rnd.OneIn(5)) {
      v3 = rnd.OneIn(2) ? v1 : v2;
      v3 = S2EdgeUtil::InterpolateAtDistance(d3, v3, ChoosePoint());
    } else {
      v3 = S2EdgeUtil::Interpolate(pow(1e-16, rnd.RandDouble()), v1, v2);
      v3 = S2EdgeUtil::InterpolateAtDistance(d3, v3,
                                             v1.CrossProd(v2).Normalize());
    }
    S2Builder::Options options((IdentitySnapFunction(snap_radius)));
    options.set_idempotent(false);
    S2Builder builder(options);
    S2Polygon output;
    output.set_s2debug_override(S2Debug::DISABLE);
    builder.StartLayer(gtl::MakeUnique<S2PolygonLayer>(&output));
    builder.ForceVertex(v3);
    builder.AddEdge(v0, v1);
    builder.AddEdge(v1, v2);
    builder.AddEdge(v2, v0);
    S2Error error;
    if (!builder.Build(&error)) {
      LOG(ERROR) << "d0=" << d0 << ", d2=" << d2 << ", d3=" << d3;
    }
    if (error.ok() && !output.is_empty()) {
      EXPECT_EQ(1, output.num_loops());
      if (output.num_loops() == 1) {
        EXPECT_TRUE(output.IsValid());
        EXPECT_EQ(s2pred::Sign(v0, v1, v2) > 0, output.loop(0)->IsNormalized())
            << "d0=" << d0 << ", d2=" << d2 << ", d3=" << d3;
        ++non_degenerate;
      }
    }
  }
  LOG(INFO) << non_degenerate << " non-degenerate out of " << kIters;
  EXPECT_GE(non_degenerate, kIters / 10);
}

TEST(S2Builder, SelfIntersectionStressTest) {
  int const kIters = 50 * FLAGS_iteration_multiplier;
  for (int iter = 0; iter < kIters; ++iter) {
    S2Testing::rnd.Reset(iter + 1);  // Easier to reproduce a specific case.
    CycleTimer timer;
    timer.Start();

    // The minimum radius is about 36cm on the Earth's surface.  The
    // performance is reduced for radii much smaller than this because
    // S2ShapeIndex only indexes regions down to about 1cm across.
    S2Cap cap = S2Testing::GetRandomCap(1e-14, 1e-2);

    S2Builder::Options options;
    options.set_split_crossing_edges(true);
    if (S2Testing::rnd.OneIn(2)) {
      S1Angle radius = cap.GetRadius();
      int min_exp = IntLatLngSnapFunction::ExponentForMaxSnapRadius(radius);
      int exponent = min(IntLatLngSnapFunction::kMaxExponent,
                         min_exp + S2Testing::rnd.Uniform(5));
      options.set_snap_function(IntLatLngSnapFunction(exponent));
    }
    S2Builder builder(options);

    // Note that the number of intersections (and the running time) is
    // quadratic in the number of vertices.  With 200 input vertices, the
    // output consists of about 2300 loops and 9000 vertices.
    S2Polygon output;
    builder.StartLayer(MakeUnique<S2PolygonLayer>(
        &output, S2PolygonLayer::Options(EdgeType::UNDIRECTED)));
    vector<S2Point> vertices(google::DEBUG_MODE ? 50 : 200);
    for (S2Point& vertex : vertices) {
      vertex = S2Testing::SamplePoint(cap);
    }
    vertices.back() = vertices.front();
    S2Polyline input(vertices);
    builder.AddPolyline(input);
    S2Error error;
    EXPECT_TRUE(builder.Build(&error)) << error.text();
    EXPECT_FALSE(output.FindValidationError(&error)) << error.text();
    if (iter == -1) {
      cout << "S2Polyline: " << s2textformat::ToString(input) << endl;
      cout << "S2Polygon: " << s2textformat::ToString(output) << endl;
    }
    if (iter < 50) {
      printf("iter=%4d: ms=%4lld, radius=%8.3g, loops=%d, vertices=%d\n",
             iter, timer.GetInMs(), cap.GetRadius().radians(),
             output.num_loops(), output.num_vertices());
    }
  }
}

TEST(S2Builder, FractalStressTest) {
  int const kIters = (google::DEBUG_MODE ? 100 : 1000) * FLAGS_iteration_multiplier;
  for (int iter = 0; iter < kIters; ++iter) {
    S2Testing::rnd.Reset(iter + 1);  // Easier to reproduce a specific case.
    S2Testing::Fractal fractal;
    fractal.SetLevelForApproxMaxEdges(google::DEBUG_MODE ? 800 : 12800);
    fractal.SetLevelForApproxMinEdges(12);
    fractal.set_fractal_dimension(1.5 + 0.5 * S2Testing::rnd.RandDouble());
    S2Polygon input(fractal.MakeLoop(S2Testing::GetRandomFrame(),
                                     S1Angle::Degrees(20)));
    S2Builder::Options options;
    if (S2Testing::rnd.OneIn(3)) {
      int exponent = S2Testing::rnd.Uniform(11);
      options.set_snap_function(IntLatLngSnapFunction(exponent));
    } else if (S2Testing::rnd.OneIn(2)) {
      int level = S2Testing::rnd.Uniform(20);
      options.set_snap_function(S2CellIdSnapFunction(level));
    } else {
      options.set_snap_function(IdentitySnapFunction(
          S1Angle::Degrees(10 * pow(1e-4, S2Testing::rnd.RandDouble()))));
    }
    S2Builder builder(options);
    S2Polygon output;
    builder.StartLayer(MakeUnique<S2PolygonLayer>(&output));
    builder.AddPolygon(input);
    S2Error error;
    EXPECT_TRUE(builder.Build(&error)) << error.text();
    EXPECT_FALSE(output.FindValidationError(&error)) << error.text();
    if (iter == -1) {
      cout << "S2Polygon: " << s2textformat::ToString(input) << endl;
      cout << "S2Polygon: " << s2textformat::ToString(output) << endl;
    }
    if (iter < 50) {
      printf("iter=%4d: in_vertices=%d, out_vertices=%d\n",
             iter, input.num_vertices(), output.num_vertices());
    }
  }
}

void TestSnappingWithForcedVertices(char const* input_str,
                                    S1Angle snap_radius,
                                    char const* vertices_str,
                                    char const* expected_str) {
  S2Builder builder((S2Builder::Options(IdentitySnapFunction(snap_radius))));
  vector<S2Point> vertices = s2textformat::ParsePoints(vertices_str);
  for (auto const& vertex : vertices) {
    builder.ForceVertex(vertex);
  }
  S2Polyline output;
  builder.StartLayer(MakeUnique<S2PolylineLayer>(&output));
  builder.AddPolyline(*MakePolyline(input_str));
  S2Error error;
  EXPECT_TRUE(builder.Build(&error)) << error.text();
  EXPECT_EQ(expected_str, s2textformat::ToString(output));
}

TEST(S2Builder, AdjacentCoverageIntervalsSpanMoreThan90Degrees) {
  // The test for whether one Voronoi site excludes another along a given
  // input edge boils down to a test of whether two angle intervals "a" and
  // "b" overlap.  Let "ra" and "rb" be the semi-widths of the two intervals,
  // and let "d" be the angle between their centers.  Then "a" contains "b" if
  // (rb + d <= ra), and "b" contains "a" if (rb - d >= ra).  However the
  // actual code uses the sines of the angles, e.g.  sin(rb + d) <= sin(ra).
  // This works fine most of the time, but the first condition (rb + d <= ra)
  // also needs to check that rb + d < 90 degrees.  This test verifies that
  // case.

  // The following 3 tests have d < 90, d = 90, and d > 90 degrees, but in all
  // 3 cases rb + d > 90 degrees.
  TestSnappingWithForcedVertices("0:0, 0:80", S1Angle::Degrees(60),
                                 "0:0, 0:70", "0:0, 0:70");
  TestSnappingWithForcedVertices("0:0, 0:80", S1Angle::Degrees(60),
                                 "0:0, 0:90", "0:0, 0:90");
  TestSnappingWithForcedVertices("0:0, 0:80", S1Angle::Degrees(60),
                                 "0:0, 0:110", "0:0, 0:110");

  // This test has d = 180 degrees, i.e. the two sites project to points that
  // are 180 degrees apart along the input edge.  The snapped edge doesn't
  // stay within max_edge_deviation() of the input edge, so an extra site is
  // added and it is snapped again (yielding two edges).  The case we are
  // testing here is the first call to SnapEdge() before adding the site.
  TestSnappingWithForcedVertices("0:10, 0:170", S1Angle::Degrees(50),
                                 "47:0, 49:180", "47:0, 0:90, 49:180");

  // This test has d = 220 degrees, i.e. when the input edge is snapped it
  // goes the "wrong way" around the sphere.  Again, the snapped edge is too
  // far from the input edge so an extra site is added and it is resnapped.
  TestSnappingWithForcedVertices("0:10, 0:170", S1Angle::Degrees(70),
                                 "0:-20, 0:-160", "0:-20, 0:90, 0:-160");

  // Without using forced vertices, the maximum angle between the coverage
  // interval centers is d = 300 degrees.  This would use an edge 180 degrees
  // long, and then place two sites 60 degrees past either endpoint.  With
  // forced vertices we can increase the snap radius to 70 degrees and get an
  // angle of up to d = 320 degrees, but the sites are only 40 degrees apart
  // (which is why it requires forced vertices).  The test below is an
  // approximation of this situation with d = 319.6 degrees.
  TestSnappingWithForcedVertices("0:0.1, 0:179.9", S1Angle::Degrees(70),
                                 "0:-69.8, 0:-110.2",
                                 "0:-69.8, 0:90, 0:-110.2");
}

TEST(S2Builder, OldS2PolygonBuilderBug) {
  // This is a polygon that caused the obsolete S2PolygonBuilder class to
  // generate an invalid output polygon (duplicate edges).
  unique_ptr<S2Polygon> input(MakePolygon(
      "32.2983095:72.3416582, 32.2986281:72.3423059, "
      "32.2985238:72.3423743, 32.2987176:72.3427807, "
      "32.2988174:72.3427056, 32.2991269:72.3433480, "
      "32.2991881:72.3433077, 32.2990668:72.3430462, "
      "32.2991745:72.3429778, 32.2995078:72.3436725, "
      "32.2996075:72.3436269, 32.2985465:72.3413832, "
      "32.2984558:72.3414530, 32.2988015:72.3421839, "
      "32.2991552:72.3429416, 32.2990498:72.3430073, "
      "32.2983764:72.3416059"));
  ASSERT_TRUE(input->IsValid());

  S1Angle snap_radius = S2Testing::MetersToAngle(20/0.866);
  S2Builder builder((S2Builder::Options(IdentitySnapFunction(snap_radius))));
  S2Polygon output;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output));
  builder.AddPolygon(*input);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  EXPECT_TRUE(output.IsValid());
  unique_ptr<S2Polygon> expected(MakePolygon(
      "32.2991552:72.3429416, 32.2991881:72.3433077, 32.2996075:72.3436269; "
      "32.2988015:72.3421839, 32.2985465:72.3413832, 32.2983764:72.3416059, "
      "32.2985238:72.3423743, 32.2987176:72.3427807"));
  ExpectPolygonsEqual(*expected, output);
}

}  // namespace
