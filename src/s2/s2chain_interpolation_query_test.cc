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

#include "s2/s2chain_interpolation_query.h"

#include <algorithm>
#include <cstdint>
#include <vector>

#include <gtest/gtest.h>
#include "absl/types/span.h"
#include "s2/s1angle.h"
#include "s2/s2latlng.h"
#include "s2/s2lax_polygon_shape.h"
#include "s2/s2lax_polyline_shape.h"
#include "s2/s2point.h"
#include "s2/s2polyline.h"
#include "s2/s2shape.h"
#include "s2/s2text_format.h"

using std::vector;

static constexpr double kEpsilon = 1.e-8;
static constexpr S1Angle kEpsilonAngle = S1Angle::Radians(kEpsilon);

TEST(S2ChainInterpolationQueryTest, SimplePolylines) {
  // Setup the test inputs.
  constexpr double kLatitudeB = 1;
  constexpr double kLatitudeC = 2.5;
  constexpr double kTotalLengthAbc = kLatitudeC;

  S2Point a = S2LatLng::FromDegrees(0, 0).ToPoint();
  S2Point b = S2LatLng::FromDegrees(kLatitudeB, 0).ToPoint();
  S2Point c = S2LatLng::FromDegrees(kLatitudeC, 0).ToPoint();

  S2LaxPolylineShape empty_shape;
  S2LaxPolylineShape shape_ac({a, c});
  S2LaxPolylineShape shape_abc({a, b, c});
  S2LaxPolylineShape shape_bb({b, b});
  S2Polyline polyline_cc({c});
  S2Polyline::Shape shape_cc(&polyline_cc);

  S2ChainInterpolationQuery uninitialized_query;
  S2ChainInterpolationQuery query_empty (&empty_shape);
  S2ChainInterpolationQuery query_ac (&shape_ac);
  S2ChainInterpolationQuery query_abc (&shape_abc);
  S2ChainInterpolationQuery query_bb (&shape_bb);
  S2ChainInterpolationQuery query_cc (&shape_cc);

  vector<S2ChainInterpolationQuery::Result> ground_truth;
  const vector<double> distances = {-1.0,
                                         0.0,
                                         1.0e-8,
                                         kLatitudeB / 2,
                                         kLatitudeB - 1.0e-7,
                                         kLatitudeB,
                                         kLatitudeB + 1.0e-5,
                                         kLatitudeB + 0.5,
                                         kLatitudeC - 10.e-7,
                                         kLatitudeC,
                                         kLatitudeC + 10.e-16,
                                         1.e6};

  for (const double& distance : distances) {
    double lat = std::clamp(distance, 0.0, kTotalLengthAbc);
    const S2Point point = S2LatLng::FromDegrees(lat, 0).ToPoint();
    int32_t edge_id = distance < kLatitudeB ? 0 : 1;
    ground_truth.emplace_back(point, edge_id, S1Angle::Degrees(lat));
  }

  // Run the tests
  const double length_empty = query_empty.GetLength().degrees();
  const double length_abc = query_abc.GetLength().degrees();
  const double length_ac = query_ac.GetLength().degrees();
  const double length_bb = query_bb.GetLength().degrees();
  const double length_cc = query_cc.GetLength().degrees();

  bool empty_query_valid = false;
  auto uninitialized_query_result = uninitialized_query.AtFraction(0);
  auto ac_result_at_infinity = query_ac.AtDistance(S1Angle::Infinity());
  const S2Point ac_point_at_infinity = ac_result_at_infinity.point();

  vector<S2ChainInterpolationQuery::Result> ac, abc, bb, empty, cc;
  for (const double& distance : distances) {
    const double total_fraction = distance / kTotalLengthAbc;

    auto empty_query_result = query_empty.AtFraction(total_fraction);
    empty_query_valid |= empty_query_result.is_valid();

    ac.push_back(query_ac.AtFraction(total_fraction));
    abc.push_back(query_abc.AtFraction(total_fraction));
    bb.push_back(query_bb.AtFraction(total_fraction));
    cc.push_back(query_cc.AtFraction(total_fraction));
  }

  // Check the test results.
  EXPECT_FALSE(uninitialized_query_result.is_valid());
  EXPECT_FALSE(empty_query_valid);
  EXPECT_TRUE(ac_result_at_infinity.is_valid());

  EXPECT_LE(length_empty, kEpsilon);
  EXPECT_NEAR(length_ac, kTotalLengthAbc, kEpsilon);
  EXPECT_NEAR(length_abc, kTotalLengthAbc, kEpsilon);
  EXPECT_LE(length_bb, kEpsilon);
  EXPECT_LE(length_cc, kEpsilon);
  EXPECT_LE(S1Angle(ac_point_at_infinity, c).degrees(),  kEpsilon);

  for (int32_t i = 0; i < ground_truth.size(); ++i) {
    EXPECT_TRUE(ac[i].is_valid());
    EXPECT_TRUE(abc[i].is_valid());
    EXPECT_TRUE(bb[i].is_valid());
    EXPECT_FALSE(cc[i].is_valid());

    EXPECT_LE(S1Angle(ac[i].point(), ground_truth[i].point()), kEpsilonAngle);
    EXPECT_LE(S1Angle(abc[i].point(), ground_truth[i].point()), kEpsilonAngle);
    EXPECT_LE(S1Angle(bb[i].point(), shape_bb.vertex(0)), kEpsilonAngle);

    EXPECT_EQ(ac[i].edge_id(), 0);
    EXPECT_EQ(bb[i].edge_id(), 0);
    EXPECT_EQ(abc[i].edge_id(), ground_truth[i].edge_id());
  }
}

TEST(S2ChainInterpolationQueryTest, Distance) {
  // Setup the test inputs.
  const vector<double> distances = {
      -1.0,         -1.0e-8,  0.0,         1.0e-8,     0.2, 0.5,
      1.0 - 1.0e-8, 1.0,      1.0 + 1.e-8, 1.2,        1.2, 1.2 + 1.0e-10,
      1.5,          1.999999, 2.0,         2.00000001, 1.e6};
  const vector<S2Point> vertices = s2textformat::ParsePointsOrDie(
      "0:0, 0:0, 1.0e-7:0, 0.1:0, 0.2:0, 0.2:0, 0.6:0, 0.999999:0, 0.999999:0, "
      "1:0, 1:0, 1.000001:0, 1.000001:0, 1.1:0, 1.2:0, 1.2000001:0, 1.7:0, "
      "1.99999999:0, 2:0");
  const double kTotalLength =
      S1Angle(vertices.front(), vertices.back()).degrees();

  S2LaxPolylineShape shape(vertices);
  S2ChainInterpolationQuery query(&shape);

  // Run the tests
  const double length = query.GetLength().degrees();
  vector<S2ChainInterpolationQuery::Result> results;
  for (const double& d : distances) {
    results.push_back(query.AtDistance(S1Angle::Degrees(d)));
  }

  // Check the test results.
  EXPECT_NEAR(length, kTotalLength, kEpsilon);

  for (int32_t i = 0; i < distances.size(); ++i) {
    EXPECT_TRUE(results[i].is_valid());

    const double d = distances[i];
    const double lat = S2LatLng(results[i].point()).lat().degrees();
    const int edge_id = results[i].edge_id();
    const S1Angle distance = results[i].distance();


    if (d < 0) {
      EXPECT_DOUBLE_EQ(lat, 0);
      EXPECT_EQ(edge_id, 0);
      EXPECT_DOUBLE_EQ(distance.degrees(), 0.0);
    } else if (d > 2) {
      EXPECT_NEAR(lat, 2, kEpsilon);
      EXPECT_EQ(edge_id, shape.num_edges() - 1);
      EXPECT_DOUBLE_EQ(distance.degrees(), kTotalLength);
    } else {
      EXPECT_NEAR(lat, d, kEpsilon);
      EXPECT_GE(edge_id, 0);
      EXPECT_LT(edge_id, shape.num_edges());
      const auto& edge = shape.edge(edge_id);
      EXPECT_GE(lat, S2LatLng(edge.v0).lat().degrees());
      EXPECT_LE(lat, S2LatLng(edge.v1).lat().degrees());
      EXPECT_NEAR(distance.degrees(), d,  kEpsilon);
    }
  }
}

TEST(S2ChainInterpolationQueryTest, Chains) {
  // Setup the test inputs.
  vector<S2LaxPolygonShape::Loop> loops;
  loops.push_back(s2textformat::ParsePointsOrDie("0:0, 1:0"));
  loops.push_back(s2textformat::ParsePointsOrDie("2:0, 3:0"));

  S2LaxPolygonShape shape(loops);

  S2ChainInterpolationQuery query(&shape);
  S2ChainInterpolationQuery query0(&shape, 0);
  S2ChainInterpolationQuery query1(&shape, 1);

  // Run the tests.
  auto query_result = query.AtFraction(0.25);
  auto query0_result = query0.AtFraction(0.25);
  auto query1_result = query1.AtFraction(0.25);

  // Check the test results.
  EXPECT_TRUE(query_result.is_valid());
  EXPECT_TRUE(query0_result.is_valid());
  EXPECT_TRUE(query1_result.is_valid());

  EXPECT_NEAR(S2LatLng(query_result.point()).lat().degrees(), 1, kEpsilon);
  EXPECT_NEAR(S2LatLng(query0_result.point()).lat().degrees(), 0.5, kEpsilon);
  EXPECT_NEAR(S2LatLng(query1_result.point()).lat().degrees(), 2.5, kEpsilon);
}

TEST(S2ChainInterpolationQueryTest, GetLengthAtEdgeEmpty) {
  const S2LaxPolylineShape shape;
  const S2ChainInterpolationQuery query(&shape);
  EXPECT_EQ(query.GetLengthAtEdgeEnd(0), S1Angle::Zero());
}

TEST(S2ChainInterpolationQueryTest, GetLengthAtEdgePolyline) {
  const S2Point vertex0 = S2LatLng::FromDegrees(0.0, 0.0).ToPoint();
  const S2Point vertex1 = S2LatLng::FromDegrees(0.0, 1.0).ToPoint();
  const S2Point vertex2 = S2LatLng::FromDegrees(0.0, 3.0).ToPoint();
  const S2Point vertex3 = S2LatLng::FromDegrees(0.0, 6.0).ToPoint();

  const S2LaxPolylineShape shape{{vertex0, vertex1, vertex2, vertex3}};

  const S2ChainInterpolationQuery query(&shape);

  EXPECT_EQ(query.GetLength(), S1Angle::Degrees(6.0));
  EXPECT_EQ(query.GetLengthAtEdgeEnd(-100), S1Angle::Infinity());
  EXPECT_EQ(query.GetLengthAtEdgeEnd(0), S1Angle::Degrees(1.0));
  EXPECT_EQ(query.GetLengthAtEdgeEnd(1), S1Angle::Degrees(3.0));
  EXPECT_EQ(query.GetLengthAtEdgeEnd(2), S1Angle::Degrees(6.0));
  EXPECT_EQ(query.GetLengthAtEdgeEnd(100), S1Angle::Infinity());
}

TEST(S2ChainInterpolationQueryTest, GetLengthAtEdgePolygon) {
  const S2Point vertex0 = S2LatLng::FromDegrees(1.0, 1.0).ToPoint();
  const S2Point vertex1 = S2LatLng::FromDegrees(2.0, 1.0).ToPoint();
  const S2Point vertex2 = S2LatLng::FromDegrees(2.0, 3.0).ToPoint();
  const S2Point vertex3 = S2LatLng::FromDegrees(1.0, 3.0).ToPoint();
  const S2Point vertex4 = S2LatLng::FromDegrees(0.0, 0.0).ToPoint();
  const S2Point vertex5 = S2LatLng::FromDegrees(0.0, 4.0).ToPoint();
  const S2Point vertex6 = S2LatLng::FromDegrees(3.0, 4.0).ToPoint();
  const S2Point vertex7 = S2LatLng::FromDegrees(3.0, 0.0).ToPoint();

  vector<S2LaxPolygonShape::Loop> loops;
  loops.push_back({vertex0, vertex1, vertex2, vertex3});
  loops.push_back({vertex4, vertex5, vertex6, vertex7});

  const S2LaxPolygonShape shape(loops);

  constexpr S1Angle kTolerance = S1Angle::Degrees(0.01);

  const S2ChainInterpolationQuery query0(&shape, 0);

  EXPECT_NEAR(query0.GetLength().degrees(), 6.0, kTolerance.degrees());
  EXPECT_EQ(query0.GetLengthAtEdgeEnd(-100), S1Angle::Infinity());
  EXPECT_NEAR(query0.GetLengthAtEdgeEnd(0).degrees(), 1.0,
              kTolerance.degrees());
  EXPECT_NEAR(query0.GetLengthAtEdgeEnd(1).degrees(), 3.0,
              kTolerance.degrees());
  EXPECT_NEAR(query0.GetLengthAtEdgeEnd(2).degrees(), 4.0,
              kTolerance.degrees());
  EXPECT_NEAR(query0.GetLengthAtEdgeEnd(3).degrees(), 6.0,
              kTolerance.degrees());
  EXPECT_EQ(query0.GetLengthAtEdgeEnd(4), S1Angle::Infinity());
  EXPECT_EQ(query0.GetLengthAtEdgeEnd(5), S1Angle::Infinity());
  EXPECT_EQ(query0.GetLengthAtEdgeEnd(6), S1Angle::Infinity());
  EXPECT_EQ(query0.GetLengthAtEdgeEnd(7), S1Angle::Infinity());
  EXPECT_EQ(query0.GetLengthAtEdgeEnd(100), S1Angle::Infinity());

  const S2ChainInterpolationQuery query1(&shape, 1);

  EXPECT_NEAR(query1.GetLength().degrees(), 14.0, kTolerance.degrees());
  EXPECT_EQ(query1.GetLengthAtEdgeEnd(-100), S1Angle::Infinity());
  EXPECT_EQ(query1.GetLengthAtEdgeEnd(0), S1Angle::Infinity());
  EXPECT_EQ(query1.GetLengthAtEdgeEnd(1), S1Angle::Infinity());
  EXPECT_EQ(query1.GetLengthAtEdgeEnd(2), S1Angle::Infinity());
  EXPECT_EQ(query1.GetLengthAtEdgeEnd(3), S1Angle::Infinity());
  EXPECT_NEAR(query1.GetLengthAtEdgeEnd(4).degrees(), 4.0,
              kTolerance.degrees());
  EXPECT_NEAR(query1.GetLengthAtEdgeEnd(5).degrees(), 7.0,
              kTolerance.degrees());
  EXPECT_NEAR(query1.GetLengthAtEdgeEnd(6).degrees(), 11.0,
              kTolerance.degrees());
  EXPECT_NEAR(query1.GetLengthAtEdgeEnd(7).degrees(), 14.0,
              kTolerance.degrees());
  EXPECT_EQ(query1.GetLengthAtEdgeEnd(100), S1Angle::Infinity());
}

TEST(S2ChainInterpolationQueryTest, Slice) {
  S2ChainInterpolationQuery empty_query;
  EXPECT_EQ(s2textformat::ToString(empty_query.Slice(0, 1)), "");

  auto polyline = s2textformat::MakeLaxPolylineOrDie("0:0, 0:1, 0:2");
  S2ChainInterpolationQuery query(polyline.get());
  EXPECT_EQ(s2textformat::ToString(query.Slice(0, 1)), "0:0, 0:1, 0:2");
  EXPECT_EQ(s2textformat::ToString(query.Slice(0, 0.5)), "0:0, 0:1");
  EXPECT_EQ(s2textformat::ToString(query.Slice(1, 0.5)), "0:2, 0:1");
  EXPECT_EQ(s2textformat::ToString(query.Slice(0.25, 0.75)),
            "0:0.5, 0:1, 0:1.5");
}
