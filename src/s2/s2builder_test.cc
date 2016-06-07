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
#include <glog/log_severity.h>
#include "s2/base/timer.h"
#include <gtest/gtest.h>
#include "s2/s2builder_layer.h"
#include "s2/s2builderutil_layers.h"
#include "s2/s2builderutil_snap_functions.h"
#include "s2/s2cap.h"
#include "s2/s2cellid.h"
#include "s2/s2latlng.h"
#include "s2/s2loop.h"
#include "s2/s2polygon.h"
#include "s2/s2polygonbuilder.h"
#include "s2/s2polyline.h"
#include "s2/s2testing.h"
#include "s2/s2textformat.h"
#include "s2/util/gtl/ptr_util.h"

using gtl::MakeUnique;
using std::cout;
using std::min;
using std::unique_ptr;
using std::vector;
using s2builderutil::IdentitySnapFunction;
using s2builderutil::S2CellIdSnapFunction;
using s2builderutil::IntLatLngSnapFunction;
using s2builderutil::S2PolygonLayer;
using s2builderutil::S2PolylineLayer;

// TODO(user): Change the s2textformat API.
// TODO(ericv): Convert this file to use this function.
static unique_ptr<S2Polyline> MakePolyline(string const& str) {
  return unique_ptr<S2Polyline>(s2textformat::MakePolyline(str));
}

static void ExpectPolygonsEqual(S2Polygon const& expected,
                                S2Polygon const& actual) {
  EXPECT_TRUE(expected.Equals(&actual))
      << "\nExpected:\n" << s2textformat::ToString(&expected)
      << "\nActual:\n" << s2textformat::ToString(&actual);
}

static void ExpectPolygonsApproxEqual(S2Polygon const& expected,
                                      S2Polygon const& actual,
                                      S1Angle tolerance) {
  EXPECT_TRUE(expected.BoundaryApproxEquals(&actual, tolerance.radians()))
      << "\nExpected:  " << s2textformat::ToString(&expected)
      << "\nActual:    " << s2textformat::ToString(&actual)
      << "\nTolerance: " << tolerance.degrees();
}

static void ExpectPolylinesEqual(S2Polyline const& expected,
                                 S2Polyline const& actual) {
  EXPECT_TRUE(expected.Equals(&actual))
      << "\nExpected:\n" << s2textformat::ToString(&expected)
      << "\nActual:\n" << s2textformat::ToString(&actual);
}

TEST(S2Builder, SimpleVertexMerging) {
  // When IdentitySnapFunction is used (i.e., no special requirements on
  // vertex locations), check that vertices closer together than the snap
  // radius are merged together.

  S1Angle tolerance = S1Angle::Degrees(0.5);
  S2Builder builder((S2Builder::Options(IdentitySnapFunction(tolerance))));
  S2Polygon output;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output));
  unique_ptr<S2Polygon> input(s2textformat::MakePolygon(
      "0:0, 0.2:0.2, 0.1:0.2, 0.1:0.9, 0:1, 0.1:1.1, 0.9:1, 1:1, 1:0.9"));
  builder.AddPolygon(*input);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  unique_ptr<S2Polygon> expected(s2textformat::MakePolygon(
      "0:0, 0:1, 1:0.9"));
  ExpectPolygonsApproxEqual(*expected, output, tolerance);
}

TEST(S2Builder, SimpleS2CellIdSnapping) {
  // When S2CellIdSnapFunction is used, check that all output vertices are the
  // centers of S2CellIds at the specified level level.

  int level = S2CellIdSnapFunction::LevelForMaxSnapRadius(S1Angle::Degrees(1));
  S2CellIdSnapFunction snap_function(level);
  S2Builder builder((S2Builder::Options(snap_function)));
  S2Polygon output;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output));
  unique_ptr<S2Polygon> input(s2textformat::MakePolygon(
      "2:2, 3:4, 2:6, 4:5, 6:6, 5:4, 6:2, 4:3"));
  builder.AddPolygon(*input);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  ASSERT_EQ(1, output.num_loops());
  S2Loop const* loop = output.loop(0);
  for (int i = 0; i < loop->num_vertices(); ++i) {
    EXPECT_EQ(S2CellId::FromPoint(loop->vertex(i)).parent(level).ToPoint(),
              loop->vertex(i));
  }
  ExpectPolygonsApproxEqual(*input, output, snap_function.snap_radius());
}

TEST(S2Builder, SimpleIntLatLngSnapping) {
  S2Builder builder(S2Builder::Options(IntLatLngSnapFunction(0)));  // E0 coords
  S2Polygon output;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output));
  unique_ptr<S2Polygon> input(s2textformat::MakePolygon(
      "2.01:2.09, 3.24:4.49, 1.78:6.25, 3.51:5.49, 6.11:6.11, "
      "5.22:3.88, 5.55:2.49, 4.49:2.51"));
  unique_ptr<S2Polygon> expected(s2textformat::MakePolygon(
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

  S1Angle tolerance = S1Angle::Degrees(1);
  S2Builder builder((S2Builder::Options(IdentitySnapFunction(tolerance))));
  S2Polygon output;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output));
  // Spacing between vertices: approximately 2*pi*20/1000 = 0.125 degrees.
  S2Polygon input(S2Loop::MakeRegularLoop(S2Point(1, 0, 0),
                                          S1Angle::Degrees(20), 1000));
  builder.AddPolygon(input);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  ASSERT_EQ(1, output.num_loops());
  EXPECT_GE(output.loop(0)->num_vertices(), 70);
  EXPECT_LE(output.loop(0)->num_vertices(), 125);
  EXPECT_TRUE(output.BoundaryNear(&input, tolerance.radians()));
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
  unique_ptr<S2Polygon> input(s2textformat::MakePolygon(
      "0:0, 0:1, 1:.9, 2:.8, 3:.7, 4:.6, 5:.5, 6:.4, 7:.3, 8:.2, 9:.1, 10:0"));
  unique_ptr<S2Polygon> expected(s2textformat::MakePolygon(
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

TEST(S2Builder, SnappingIdempotence) {
  // Check that re-snapping a polygon does not change it.

  unique_ptr<S2Polygon> input(s2textformat::MakePolygon(
      "0:0, 0:10, 5:10, 1:5, 5:0"));
  unique_ptr<S2Polygon> expected(s2textformat::MakePolygon(
      "0.5:0, 0.5:10, 5:10, 1:5, 5:0"));
  S2Builder::Options options(IdentitySnapFunction(S1Angle::Degrees(0.6)));
  S2Builder builder(options);
  S2Polygon output1, output2;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output1));
  builder.AddPolygon(*input);
  builder.ForceVertex(S2LatLng::FromDegrees(0.5, 0).ToPoint());
  builder.ForceVertex(S2LatLng::FromDegrees(0.5, 10).ToPoint());
  S2Error error;
  ASSERT_TRUE(builder.Build(&error));
  ExpectPolygonsEqual(*expected, output1);

  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output2));
  builder.AddPolygon(output1);
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  ExpectPolygonsEqual(output1, output2);
}

TEST(S2Builder, kMaxSnapRadius) {
  // Verify that kMaxSnapRadius will allow snapping at S2CellId level 2.
  EXPECT_LE(S2CellIdSnapFunction::MinSnapRadiusForLevel(2),
            S2Builder::SnapFunction::kMaxSnapRadius());
}

TEST(S2Builder, S2CellIdSnappingAtAllLevels) {
  unique_ptr<const S2Polygon> input(
      s2textformat::MakePolygon("0:0, 0:2, 2:0; 0:0, 0:-2, -2:-2, -2:0"));
  for (int level = 0; level <= S2CellId::kMaxLevel; ++level) {
    S2CellIdSnapFunction snap_function(level);
    S2Builder builder((S2Builder::Options(snap_function)));
    S2Polygon output;
    builder.StartLayer(MakeUnique<S2PolygonLayer>(&output));
    builder.AddPolygon(*input);
    S2Error error;
    ASSERT_TRUE(builder.Build(&error)) << error.text();
    EXPECT_TRUE(output.IsValid());
    // S2PolygonBuilder uses a different snapping model, so we need to expand
    // the tolerance slightly to ensure that one polygon contains the other.
    S1Angle tolerance = (snap_function.snap_radius() /
                         S2PolygonBuilderOptions().edge_splice_fraction());
    EXPECT_TRUE(output.ApproxContains(input.get(), tolerance));
    EXPECT_TRUE(input->ApproxContains(&output, tolerance));
  }
}

TEST(S2Builder, SnappingDoesNotRotateVertices) {
  // This is already tested extensively elsewhere.
  unique_ptr<S2Polygon> input(s2textformat::MakePolygon(
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
  unique_ptr<S2Polyline> input(s2textformat::MakePolyline(
      "3:1, 1:3, 1:1, 3:3"));
  unique_ptr<S2Polyline> expected(s2textformat::MakePolyline(
      "3:1, 2:2, 1:3, 1:1, 2:2, 3:3"));
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
      &output, S2PolygonLayer::Options(S2Builder::EdgeType::UNDIRECTED)));
  unique_ptr<S2Polyline> input(s2textformat::MakePolyline(
      "3:1, 1:3, 1:1, 3:3, 3:1"));
  unique_ptr<S2Polygon> expected(s2textformat::MakePolygon(
      "1:1, 1:3, 2:2; 3:1, 2:2, 3:3"));
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

TEST(S2Builder, SelfIntersectionStressTest) {
  for (int iter = 0; iter < 50; ++iter) {
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
        &output, S2PolygonLayer::Options(S2Builder::EdgeType::UNDIRECTED)));
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
      cout << "S2Polyline: " << s2textformat::ToString(&input) << "\n";
      cout << "S2Polygon: " << s2textformat::ToString(&output) << "\n";
    }
    printf("iter=%4d: ms=%4lld, radius=%8.3g, loops=%d, vertices=%d\n",
           iter, timer.GetInMs(), cap.GetRadius().radians(),
           output.num_loops(), output.num_vertices());
  }
}

TEST(S2Builder, FractalStressTest) {
  for (int iter = 0; iter < (google::DEBUG_MODE ? 100 : 1000); ++iter) {
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
      cout << "S2Polygon: " << s2textformat::ToString(&input) << "\n";
      cout << "S2Polygon: " << s2textformat::ToString(&output) << "\n";
    }
    printf("iter=%4d: in_vertices=%d, out_vertices=%d\n",
           iter, input.num_vertices(), output.num_vertices());
  }
}

void TestSnappingWithForcedVertices(char const* input_str,
                                    S1Angle snap_radius,
                                    char const* vertices_str,
                                    char const* expected_str) {
  S2Builder builder((S2Builder::Options(IdentitySnapFunction(snap_radius))));
  vector<S2Point> vertices;
  s2textformat::ParsePoints(vertices_str, &vertices);
  for (auto const& vertex : vertices) {
    builder.ForceVertex(vertex);
  }
  S2Polyline output;
  builder.StartLayer(MakeUnique<S2PolylineLayer>(&output));
  builder.AddPolyline(*MakePolyline(input_str));
  S2Error error;
  EXPECT_TRUE(builder.Build(&error)) << error.text();
  EXPECT_EQ(expected_str, s2textformat::ToString(&output));
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
  unique_ptr<S2Polygon> input(
      s2textformat::MakePolygon(
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

  S1Angle tolerance = S2Testing::MetersToAngle(20/0.866);
  S2Builder builder((S2Builder::Options(IdentitySnapFunction(tolerance))));
  S2Polygon output;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output));
  builder.AddPolygon(*input);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error)) << error.text();
  EXPECT_TRUE(output.IsValid());
  unique_ptr<S2Polygon> expected(s2textformat::MakePolygon(
      "32.2991552:72.3429416, 32.2991881:72.3433077, 32.2996075:72.3436269; "
      "32.2988015:72.3421839, 32.2985465:72.3413832, 32.2983764:72.3416059, "
      "32.2985238:72.3423743, 32.2987176:72.3427807"));
  ExpectPolygonsEqual(*expected, output);
}

#if 0
TEST(S2PolygonBuilder, FractalStressTest) {
  // This test fails -- its purpose is to check whether S2Builder is in fact
  // more robust than the code it will replace.
  for (int iter = 0; iter < 1000; ++iter) {
    S2Testing::rnd.Reset(iter + 1);  // Easier to reproduce a specific case.
    S2Testing::Fractal fractal;
    fractal.SetLevelForApproxMaxEdges(12800);
    fractal.SetLevelForApproxMinEdges(12);
    fractal.set_fractal_dimension(1.5 + 0.5 * S2Testing::rnd.RandDouble());
    S2Polygon input(fractal.MakeLoop(S2Testing::GetRandomFrame(),
                                     S1Angle::Degrees(20)));
    int level = 2 + S2Testing::rnd.Uniform(20);
    S2PolygonBuilderOptions options;
    options.SetRobustnessRadius(
        S1Angle::Radians(S2::kMaxDiag.GetValue(level) / 2 + 1e-15));
    options.set_snap_to_cell_centers(false);
    S2PolygonBuilder builder(options);
    builder.AddPolygon(&input);
    S2Polygon output;
    ASSERT_TRUE(builder.AssemblePolygon(&output, nullptr));
    EXPECT_TRUE(output.IsValid());
  }
}
#endif
