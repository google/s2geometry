// Copyright Google Inc. All Rights Reserved.
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

#include "s2/s2hausdorff_distance_query.h"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include "absl/types/optional.h"
#include "s2/mutable_s2shape_index.h"
#include "s2/s1chord_angle.h"
#include "s2/s2latlng.h"
#include "s2/s2lax_polygon_shape.h"
#include "s2/s2lax_polyline_shape.h"
#include "s2/s2point.h"
#include "s2/s2point_vector_shape.h"
#include "s2/s2shape.h"
#include "s2/s2shape_index.h"
#include "s2/s2text_format.h"

using s2textformat::ParsePointsOrDie;
using std::make_unique;
using std::vector;
using DirectedResult = S2HausdorffDistanceQuery::DirectedResult;
using Result = S2HausdorffDistanceQuery::Result;
using Options = S2HausdorffDistanceQuery::Options;

// Test the constructors and accessors of Result and DirectedResult.
TEST(S2HausdorffDistanceQueryTest, ResultConstructorsAndAccessorsWork) {
  S2Point point1 = S2LatLng::FromDegrees(3, 4).ToPoint();
  S2Point point2 = S2LatLng::FromDegrees(5, 6).ToPoint();
  S1ChordAngle distance1 = S1ChordAngle::Degrees(5);
  S1ChordAngle distance2 = S1ChordAngle::Degrees(5);

  DirectedResult directed_result1(distance1, point1);
  DirectedResult directed_result2(distance2, point2);
  Result result12(directed_result1, directed_result2);

  EXPECT_EQ(directed_result1.target_point(), point1);
  EXPECT_EQ(directed_result1.distance(), distance1);
  EXPECT_EQ(directed_result2.target_point(), point2);
  EXPECT_EQ(directed_result2.distance(), distance2);

  EXPECT_EQ(result12.target_to_source().target_point(), point1);
  EXPECT_EQ(result12.source_to_target().target_point(), point2);
  EXPECT_EQ(result12.distance(), directed_result2.distance());
}

// Test the constructors and accessors of the Options.
TEST(S2HausdorffDistanceQueryTest, OptionsConstructorsAndAccessorsWork) {
  Options default_options;
  Options options;
  options.set_include_interiors(!default_options.include_interiors());

  EXPECT_TRUE(default_options.include_interiors());
  EXPECT_FALSE(options.include_interiors());
}

// Test the constructors and accessors of the Options.
TEST(S2HausdorffDistanceQueryTest, QueryOptionsAccessorsWorks) {
  S2HausdorffDistanceQuery query;
  const bool default_include_interiors = query.options().include_interiors();

  query.mutable_options()->set_include_interiors(!default_include_interiors);
  const bool modified_include_interiors = query.options().include_interiors();

  EXPECT_TRUE(default_include_interiors);
  EXPECT_FALSE(modified_include_interiors);
}

// Test involving 2 simple polyline shape indexes.
TEST(S2HausdorffDistanceQueryTest, SimplePolylineQueriesSucceed) {
  const vector<S2Point> a0 = ParsePointsOrDie("0:0, 0:1, 0:1.5");
  const vector<S2Point> a1 = ParsePointsOrDie("0:2, 0:1.5, -10:1");
  const vector<S2Point> b0 = ParsePointsOrDie("1:0, 1:1, 3:2");

  // Setup the shape indexes.
  const MutableS2ShapeIndex empty_index;

  // Shape index a consists of 2 polylines, a0 and a1.
  MutableS2ShapeIndex a;
  a.Add(make_unique<S2LaxPolylineShape>(a0));
  a.Add(make_unique<S2LaxPolylineShape>(a1));

  // Shape index b consists of 1 polylines: b0.
  MutableS2ShapeIndex b;
  b.Add(make_unique<S2LaxPolylineShape>(b0));

  // Calculate expected distances.
  // Directed a to b HD is achieved at the vertex 2 of a1 and vertex 1 of b0.
  const S1ChordAngle expected_a_to_b(a1[2], b0[1]);
  // Directed b to a HD is achieved at the vertex 2 of b0 and vertex 0 of a1.
  const S1ChordAngle expected_b_to_a(b0[2], a1[0]);

  Options options;
  S2HausdorffDistanceQuery query(options);

  absl::optional<DirectedResult> directed_empty_to_a =
      query.GetDirectedResult(&empty_index, &a);

  absl::optional<DirectedResult> directed_a_to_empty =
      query.GetDirectedResult(&a, &empty_index);

  S1ChordAngle directed_a_to_empty_distance =
      query.GetDirectedDistance(&a, &empty_index);
  // These two should be false since an empty set is an infinite distance away.
  bool empty_to_a_directed_distance_less = query.IsDirectedDistanceLess(
      &empty_index, &a, S1ChordAngle::Degrees(360));
  bool a_to_empty_directed_distance_less = query.IsDirectedDistanceLess(
      &a, &empty_index, S1ChordAngle::Degrees(360));

  EXPECT_FALSE(directed_empty_to_a);
  EXPECT_FALSE(directed_a_to_empty);
  EXPECT_FALSE(directed_a_to_empty);
  EXPECT_TRUE(directed_a_to_empty_distance.is_infinity());
  EXPECT_FALSE(empty_to_a_directed_distance_less);
  EXPECT_FALSE(a_to_empty_directed_distance_less);

  absl::optional<DirectedResult> directed_a_to_b =
      query.GetDirectedResult(&a, &b);
  absl::optional<DirectedResult> directed_b_to_a =
      query.GetDirectedResult(&b, &a);
  S1ChordAngle directed_a_to_b_distance = query.GetDirectedDistance(&a, &b);

  // Tests for IsDirectedDistanceLess with limits near the Hausdorff
  // distance.
  bool a_to_b_directed_distance_less_than_distance_plus =
      query.IsDirectedDistanceLess(
          &a, &b, directed_a_to_b_distance + S1ChordAngle::Degrees(1.0));
  bool a_to_b_directed_distance_less_than_distance_minus =
      query.IsDirectedDistanceLess(
          &a, &b, directed_a_to_b_distance - S1ChordAngle::Degrees(1.0));

  EXPECT_TRUE(directed_a_to_b);
  EXPECT_TRUE(directed_b_to_a);

  EXPECT_DOUBLE_EQ(directed_a_to_b->distance().degrees(),
                   expected_a_to_b.degrees());
  EXPECT_DOUBLE_EQ(directed_a_to_b_distance.degrees(),
                   expected_a_to_b.degrees());
  EXPECT_DOUBLE_EQ(directed_b_to_a->distance().degrees(),
                   expected_b_to_a.degrees());
  EXPECT_TRUE(a_to_b_directed_distance_less_than_distance_plus);
  EXPECT_FALSE(a_to_b_directed_distance_less_than_distance_minus);

  // Tests for undirected cases.
  absl::optional<Result> a_to_b = query.GetResult(&a, &b);
  absl::optional<Result> b_to_a = query.GetResult(&b, &a);
  S1ChordAngle b_to_a_distance = query.GetDistance(&b, &a);
  absl::optional<Result> bb = query.GetResult(&b, &b);

  EXPECT_TRUE(a_to_b);
  EXPECT_TRUE(b_to_a);
  EXPECT_TRUE(bb);

  // Tests for IsDistanceLess with limits near the Hausdorff distance and
  // the average of the two directed Hausdorff distances.
  double larger_a_and_b_distance =
      std::max(directed_a_to_b->distance().radians(),
               directed_b_to_a->distance().radians());
  double smaller_a_and_b_distance =
      std::min(directed_a_to_b->distance().radians(),
               directed_b_to_a->distance().radians());
  double average_a_and_b_distance =
      (larger_a_and_b_distance + smaller_a_and_b_distance) / 2.0;

  // THis should be true if we add a small epsilon upwards to account for any
  // floating point error.
  bool distance_less_larger_distance = query.IsDistanceLess(
      &a, &b, S1ChordAngle::Radians(larger_a_and_b_distance + 0.001));
  // The average should cause one direction to succeed and the other to fail so
  // overall this should return false.
  bool distance_less_average_distance = query.IsDistanceLess(
      &a, &b, S1ChordAngle::Radians(average_a_and_b_distance));
  // The IsWithin(Directed)DistanceLimit methods are inclusive so subtract some
  // small epsilon to cause the method to return false.
  bool distance_less_smaller_distance = query.IsDistanceLess(
      &a, &b, S1ChordAngle::Radians(smaller_a_and_b_distance - 0.001));

  bool bb_always_within =
      query.IsDistanceLess(&b, &b, S1ChordAngle::Degrees(0));

  EXPECT_DOUBLE_EQ(a_to_b->distance().degrees(), b_to_a->distance().degrees());
  EXPECT_DOUBLE_EQ(bb->distance().degrees(), 0);
  EXPECT_EQ(
      a_to_b->distance().degrees(),
      std::max(a_to_b->distance().degrees(), b_to_a->distance().degrees()));
  EXPECT_EQ(b_to_a_distance.degrees(), b_to_a->distance().degrees());

  EXPECT_TRUE(distance_less_larger_distance);
  EXPECT_FALSE(distance_less_average_distance);
  EXPECT_FALSE(distance_less_smaller_distance);
  EXPECT_TRUE(bb_always_within);
}

// Test involving a polyline shape (dimension == 1) and a point shape (dimension
// == 0).
TEST(S2HausdorffDistanceQueryTest, PointVectorShapeQueriesSucceed) {
  // Points for the polyline shape.
  const vector<S2Point> a_points = ParsePointsOrDie("2:0, 0:1, 1:2, 0:3, 0:4");
  // Points for the point vector shape.
  const vector<S2Point> b_points = ParsePointsOrDie("-1:2, -0.5:0.5, -0.5:3.5");

  MutableS2ShapeIndex a;
  a.Add(make_unique<S2LaxPolylineShape>(a_points));

  MutableS2ShapeIndex b;
  b.Add(make_unique<S2PointVectorShape>(b_points));

  Options options;
  S2HausdorffDistanceQuery query(options);

  // Directed Hausdorff distance from a to b is achieved at the vertex 0 of a
  // and vertex 1 of b.
  const S1ChordAngle expected_a_to_b(a_points[0], b_points[1]);

  // Directed Hausdorff distance from b to a is achieved at the vertex 0 of b
  // and vertex 3 of a.
  const S1ChordAngle expected_b_to_a(b_points[0], a_points[3]);

  // Undirected Hausdorff distance between a and b is the maximum of the two
  // directed Hausdorff distances.
  const S1ChordAngle expected_a_b = std::max(expected_a_to_b, expected_b_to_a);

  absl::optional<DirectedResult> directed_a_to_b =
      query.GetDirectedResult(&a, &b);
  absl::optional<DirectedResult> directed_b_to_a =
      query.GetDirectedResult(&b, &a);
  S1ChordAngle undirected_a_b = query.GetDistance(&a, &b);

  EXPECT_TRUE(directed_a_to_b);
  EXPECT_TRUE(directed_b_to_a);
  EXPECT_FALSE(undirected_a_b.is_infinity());
  EXPECT_EQ(undirected_a_b.degrees(), expected_a_b.degrees());
  EXPECT_DOUBLE_EQ(directed_a_to_b->distance().degrees(),
                   expected_a_to_b.degrees());
  EXPECT_EQ(directed_a_to_b->target_point(), a_points[0]);
  EXPECT_DOUBLE_EQ(directed_b_to_a->distance().degrees(),
                   expected_b_to_a.degrees());
  EXPECT_EQ(directed_b_to_a->target_point(), b_points[0]);

  bool a_to_b_directed_distance_less_plus = query.IsDirectedDistanceLess(
      &a, &b, S1ChordAngle::Degrees(expected_a_to_b.degrees() + 0.01));
  bool b_to_a_directed_distance_less_plus = query.IsDirectedDistanceLess(
      &b, &a, S1ChordAngle::Degrees(expected_b_to_a.degrees() + 0.01));
  bool a_to_b_directed_distance_less_minus = query.IsDirectedDistanceLess(
      &a, &b, S1ChordAngle::Degrees(expected_a_to_b.degrees() - 0.01));
  bool b_to_a_directed_distance_less_minus = query.IsDirectedDistanceLess(
      &b, &a, S1ChordAngle::Degrees(expected_b_to_a.degrees() - 0.01));

  bool a_b_distance_less_plus = query.IsDistanceLess(
      &a, &b, S1ChordAngle::Degrees(expected_a_b.degrees() + 0.01));
  bool b_a_distance_less_minus = query.IsDistanceLess(
      &b, &a, S1ChordAngle::Degrees(expected_b_to_a.degrees() - 0.01));

  EXPECT_TRUE(a_to_b_directed_distance_less_plus);
  EXPECT_TRUE(b_to_a_directed_distance_less_plus);
  EXPECT_FALSE(a_to_b_directed_distance_less_minus);
  EXPECT_FALSE(b_to_a_directed_distance_less_minus);
  EXPECT_TRUE(a_b_distance_less_plus);
  EXPECT_FALSE(b_a_distance_less_minus);
}

// Test involving partially overlapping polygons.
TEST(S2HausdorffDistanceQueryTest, OverlappingPolygons) {
  // The first polygon is a triangle. It's first two vertices are inside the
  // quadrangle b (defined below), and the last vertex is outside of b.
  MutableS2ShapeIndex a;
  a.Add(s2textformat::MakeLaxPolygonOrDie("1:1, 1:2, 3.5:1.5"));

  // The other polygon is a quadrangle.
  MutableS2ShapeIndex b;
  b.Add(s2textformat::MakeLaxPolygonOrDie("0:0, 0:3, 3:3, 3:0"));

  // A triangle.
  MutableS2ShapeIndex c;
  c.Add(s2textformat::MakeLaxPolygonOrDie("0:0, 0:2, 3:0"));

  // Error tolerance to account for the difference between the northern edge of
  // the quadrangle, which is a geodesic line, and the parallel lat=3 connecting
  // the vertices of that edge.
  static constexpr double kEpsilon = 3.0e-3;

  // The first query does not include the interiors.
  Options options;
  options.set_include_interiors(false);
  S2HausdorffDistanceQuery query_1(options);

  // The directed Hausdorff distance from the first query is achieved on the
  // vertex of the triangle that is inside the quadrangle, and is approximately
  // 1 degree away from the nearest edge of the quadrangle.
  S2Point expected_target_point_1 = S2LatLng::FromDegrees(1, 2).ToPoint();

  absl::optional<DirectedResult> a_to_b_1 = query_1.GetDirectedResult(&a, &b);

  bool c_to_b_less_than = query_1.IsDirectedDistanceLess(
      &c, &b, S1ChordAngle::Degrees(1.0 + kEpsilon));

  EXPECT_TRUE(a_to_b_1);
  EXPECT_NEAR(a_to_b_1->distance().degrees(), 1, kEpsilon);
  EXPECT_EQ(a_to_b_1->target_point(), expected_target_point_1);
  EXPECT_TRUE(c_to_b_less_than);

  // The second query has include_interiors set to true.
  options.set_include_interiors(true);
  S2HausdorffDistanceQuery query_2(options);

  // The directed Hausdorff distance from the second query is achieved on the
  // last vertex of the triangle that is outside the quadrangle, and is about
  // 0.5 degrees away from the nearest edge of the quadrangle.
  S2Point expected_target_point_2 = S2LatLng::FromDegrees(3.5, 1.5).ToPoint();

  absl::optional<DirectedResult> a_to_b_2 = query_2.GetDirectedResult(&a, &b);
  // C is fully contained in B so all points are 0 distance to B.
  bool c_to_b_less_than_2 =
      query_2.IsDirectedDistanceLess(&c, &b, S1ChordAngle::Degrees(kEpsilon));

  EXPECT_TRUE(a_to_b_2);
  EXPECT_NEAR(a_to_b_2->distance().degrees(), 0.5, kEpsilon);
  EXPECT_EQ(a_to_b_2->target_point(), expected_target_point_2);
  EXPECT_TRUE(c_to_b_less_than_2);
}

// Test involving full geometries.
TEST(S2HausdorffDistanceQueryTest, WholeWorld) {
  MutableS2ShapeIndex a;
  a.Add(std::make_unique<S2PointVectorShape>(
      s2textformat::ParsePointsOrDie("1:1")));

  std::unique_ptr<S2ShapeIndex> b = s2textformat::MakeIndexOrDie("# # full");

  Options options;
  options.set_include_interiors(true);
  S2HausdorffDistanceQuery query_1(options);
  absl::optional<DirectedResult> a_to_b_1 =
      query_1.GetDirectedResult(&a, b.get());

  EXPECT_TRUE(a_to_b_1);
  EXPECT_EQ(a_to_b_1->distance().degrees(), 0.0);

  // Going from full geometry to non-full geometry should return empty option.
  absl::optional<DirectedResult> b_to_a_1 =
      query_1.GetDirectedResult(b.get(), &a);

  EXPECT_FALSE(b_to_a_1);

  absl::optional<Result> undirected_a_to_b = query_1.GetResult(b.get(), &a);
  absl::optional<Result> undirected_b_to_a = query_1.GetResult(&a, b.get());

  EXPECT_FALSE(undirected_a_to_b);
  EXPECT_FALSE(undirected_b_to_a);

  // A point to the whole world should always work.
  bool a_to_b_directed_distance_less_zero =
      query_1.IsDirectedDistanceLess(&a, b.get(), S1ChordAngle::Zero());
  // The whole world to a point should always fail.
  bool b_to_a_directed_distance_less_inf =
      query_1.IsDirectedDistanceLess(b.get(), &a, S1ChordAngle::Infinity());
  // In undirected case, must consider distance from full geometry to single
  // point which is infinite.
  bool a_b_distance_less_zero =
      query_1.IsDistanceLess(&a, b.get(), S1ChordAngle::Infinity());

  EXPECT_TRUE(a_to_b_directed_distance_less_zero);
  EXPECT_FALSE(b_to_a_directed_distance_less_inf);
  EXPECT_FALSE(a_b_distance_less_zero);
}

TEST(S2HausdorffDistanceQueryTest, WholeWorldSameReference) {
  std::unique_ptr<S2ShapeIndex> a = s2textformat::MakeIndexOrDie("# # full");
  std::unique_ptr<S2ShapeIndex> b = s2textformat::MakeIndexOrDie("# # full");

  Options options;
  options.set_include_interiors(true);
  S2HausdorffDistanceQuery query_1(options);
  absl::optional<Result> a_to_b = query_1.GetResult(a.get(), b.get());
  EXPECT_FALSE(a_to_b);

  S2HausdorffDistanceQuery query_2(options);
  absl::optional<Result> a_to_a = query_1.GetResult(a.get(), a.get());
  EXPECT_FALSE(a_to_a);

  bool a_to_b_distance_less_inf =
      query_2.IsDistanceLess(a.get(), b.get(), S1ChordAngle::Infinity());
  bool a_to_a_distance_less_inf =
      query_2.IsDistanceLess(a.get(), a.get(), S1ChordAngle::Infinity());

  EXPECT_FALSE(a_to_b_distance_less_inf);
  EXPECT_FALSE(a_to_a_distance_less_inf);
}
