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

#include "s2/s2closest_edge_query.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "absl/container/flat_hash_set.h"
#include "absl/flags/flag.h"
#include "absl/log/absl_check.h"
#include "absl/log/absl_log.h"
#include "absl/log/log_streamer.h"
#include "absl/random/bit_gen_ref.h"
#include "absl/random/random.h"
#include "absl/types/span.h"

#include "s2/encoded_s2shape_index.h"
#include "s2/mutable_s2shape_index.h"
#include "s2/s1angle.h"
#include "s2/s1chord_angle.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2closest_edge_query_base.h"
#include "s2/s2closest_edge_query_testing.h"
#include "s2/s2edge_crossings.h"
#include "s2/s2edge_distances.h"
#include "s2/s2fractal.h"
#include "s2/s2latlng.h"
#include "s2/s2loop.h"
#include "s2/s2metrics.h"
#include "s2/s2min_distance_targets.h"
#include "s2/s2point.h"
#include "s2/s2point_vector_shape.h"
#include "s2/s2pointutil.h"
#include "s2/s2polygon.h"
#include "s2/s2predicates.h"
#include "s2/s2random.h"
#include "s2/s2shape.h"
#include "s2/s2shapeutil_coding.h"
#include "s2/s2shapeutil_count_edges.h"
#include "s2/s2shapeutil_shape_edge_id.h"
#include "s2/s2testing.h"
#include "s2/s2text_format.h"
#include "s2/util/math/matrix3x3.h"

using s2shapeutil::ShapeEdgeId;
using s2textformat::MakeIndexOrDie;
using s2textformat::MakePointOrDie;
using std::make_unique;
using std::min;
using std::pair;
using std::string;
using std::unique_ptr;
using std::vector;

TEST(S2ClosestEdgeQuery, NoEdges) {
  MutableS2ShapeIndex index;
  S2ClosestEdgeQuery query(&index);
  S2ClosestEdgeQuery::PointTarget target(S2Point(1, 0, 0));
  const auto edge = query.FindClosestEdge(&target);
  EXPECT_EQ(S1ChordAngle::Infinity(), edge.distance());
  EXPECT_EQ(-1, edge.shape_id());
  EXPECT_EQ(-1, edge.edge_id());
  EXPECT_FALSE(edge.is_interior());
  EXPECT_TRUE(edge.is_empty());
  EXPECT_EQ(S1ChordAngle::Infinity(), query.GetDistance(&target));
}

TEST(S2ClosestEdgeQuery, OptionsNotModified) {
  // Tests that FindClosestEdge(), GetDistance(), and IsDistanceLess() do not
  // modify query.options(), even though all of these methods have their own
  // specific options requirements.
  S2ClosestEdgeQuery::Options options;
  options.set_max_results(3);
  options.set_max_distance(S1ChordAngle::Degrees(3));
  options.set_max_error(S1ChordAngle::Degrees(0.001));
  auto index = MakeIndexOrDie("1:1 | 1:2 | 1:3 # #");
  S2ClosestEdgeQuery query(index.get(), options);
  S2ClosestEdgeQuery::PointTarget target(MakePointOrDie("2:2"));
  EXPECT_EQ(1, query.FindClosestEdge(&target).edge_id());
  EXPECT_NEAR(1.0, query.GetDistance(&target).degrees(), 1e-15);
  EXPECT_TRUE(query.IsDistanceLess(&target, S1ChordAngle::Degrees(1.5)));

  // Verify that none of the options above were modified.
  EXPECT_EQ(options.max_results(), query.options().max_results());
  EXPECT_EQ(options.max_distance(), query.options().max_distance());
  EXPECT_EQ(options.max_error(), query.options().max_error());
}

TEST(S2ClosestEdgeQuery, OptionsS1AngleSetters) {
  // Verify that the S1Angle and S1ChordAngle versions do the same thing.
  // This is mainly to prevent the (so far unused) S1Angle versions from
  // being detected as dead code.
  S2ClosestEdgeQuery::Options angle_options, chord_angle_options;
  angle_options.set_max_distance(S1Angle::Degrees(1));
  chord_angle_options.set_max_distance(S1ChordAngle::Degrees(1));
  EXPECT_EQ(chord_angle_options.max_distance(), angle_options.max_distance());

  angle_options.set_inclusive_max_distance(S1Angle::Degrees(1));
  chord_angle_options.set_inclusive_max_distance(S1ChordAngle::Degrees(1));
  EXPECT_EQ(chord_angle_options.max_distance(), angle_options.max_distance());

  angle_options.set_conservative_max_distance(S1Angle::Degrees(1));
  chord_angle_options.set_conservative_max_distance(S1ChordAngle::Degrees(1));
  EXPECT_EQ(chord_angle_options.max_distance(), angle_options.max_distance());
}

TEST(S2ClosestEdgeQuery, DistanceEqualToLimit) {
  // Tests the behavior of IsDistanceLess, IsDistanceLessOrEqual, and
  // IsConservativeDistanceLessOrEqual (and the corresponding Options) when
  // the distance to the target exactly equals the chosen limit.
  S2Point p0(MakePointOrDie("23:12")), p1(MakePointOrDie("47:11"));
  vector<S2Point> index_points{p0};
  MutableS2ShapeIndex index;
  index.Add(make_unique<S2PointVectorShape>(index_points));
  S2ClosestEdgeQuery query(&index);

  // Start with two identical points and a zero distance.
  S2ClosestEdgeQuery::PointTarget target0(p0);
  S1ChordAngle dist0 = S1ChordAngle::Zero();
  EXPECT_FALSE(query.IsDistanceLess(&target0, dist0));
  EXPECT_TRUE(query.IsDistanceLessOrEqual(&target0, dist0));
  EXPECT_TRUE(query.IsConservativeDistanceLessOrEqual(&target0, dist0));

  // Now try two points separated by a non-zero distance.
  S2ClosestEdgeQuery::PointTarget target1(p1);
  S1ChordAngle dist1(p0, p1);
  EXPECT_FALSE(query.IsDistanceLess(&target1, dist1));
  EXPECT_TRUE(query.IsDistanceLessOrEqual(&target1, dist1));
  EXPECT_TRUE(query.IsConservativeDistanceLessOrEqual(&target1, dist1));
}

TEST(S2ClosestEdgeQuery, TrueDistanceLessThanS1ChordAngleDistance) {
  // Tests that IsConservativeDistanceLessOrEqual returns points where the
  // true distance is slightly less than the one computed by S1ChordAngle.
  //
  // The points below had the worst error from among 100,000 random pairs.
  S2Point p0(0.78516762584829192, -0.50200400690845970, -0.36263449417782678);
  S2Point p1(0.78563011732429433, -0.50187655940493503, -0.36180828883938054);

  // The S1ChordAngle distance is ~4 ulps greater than the true distance.
  S1ChordAngle dist1(p0, p1);
  auto limit = dist1.Predecessor().Predecessor().Predecessor().Predecessor();
  ASSERT_LT(s2pred::CompareDistance(p0, p1, limit), 0);

  // Verify that IsConservativeDistanceLessOrEqual() still returns "p1".
  vector<S2Point> index_points{p0};
  MutableS2ShapeIndex index;
  index.Add(make_unique<S2PointVectorShape>(index_points));
  S2ClosestEdgeQuery query(&index);
  S2ClosestEdgeQuery::PointTarget target1(p1);
  EXPECT_FALSE(query.IsDistanceLess(&target1, limit));
  EXPECT_FALSE(query.IsDistanceLessOrEqual(&target1, limit));
  EXPECT_TRUE(query.IsConservativeDistanceLessOrEqual(&target1, limit));
}

TEST(S2ClosestEdgeQuery, TestReuseOfQuery) {
  // Tests that between queries, the internal mechanism for de-duplicating
  // results is re-set.  See b/71646017.
  auto index = MakeIndexOrDie("2:2 # #");
  S2ClosestEdgeQuery query(index.get());
  query.mutable_options()->set_max_error(S1Angle::Degrees(1));
  auto target_index = MakeIndexOrDie("## 0:0, 0:5, 5:5, 5:0");
  S2ClosestEdgeQuery::ShapeIndexTarget target(target_index.get());
  auto results1 = query.FindClosestEdges(&target);
  auto results2 = query.FindClosestEdges(&target);
  EXPECT_EQ(results1.size(), results2.size());
}

TEST(S2ClosestEdgeQuery, TargetPointInsideIndexedPolygon) {
  // Tests a target point in the interior of an indexed polygon.
  // (The index also includes a polyline loop with no interior.)
  auto index = MakeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
  S2ClosestEdgeQuery::Options options;
  options.set_include_interiors(true);
  options.set_max_distance(S1Angle::Degrees(1));
  S2ClosestEdgeQuery query(index.get(), options);
  S2ClosestEdgeQuery::PointTarget target(MakePointOrDie("2:12"));
  auto results = query.FindClosestEdges(&target);
  ASSERT_EQ(1, results.size());
  EXPECT_EQ(S1ChordAngle::Zero(), results[0].distance());
  EXPECT_EQ(1, results[0].shape_id());
  EXPECT_EQ(-1, results[0].edge_id());
  EXPECT_TRUE(results[0].is_interior());
  EXPECT_FALSE(results[0].is_empty());
}

TEST(S2ClosestEdgeQuery, ShapeFilteringWorks) {
  // Tests a target point in the interior of an indexed polygon.
  // (The index also includes a polyline loop with no interior.)
  auto index =
      MakeIndexOrDie("## 1:1, 1:-1, -1:-1, -1:1 | 2:2, 2:-2, -2:-2, -2:2");
  EXPECT_EQ(index->num_shape_ids(), 2);

  S2ClosestEdgeQuery::Options options;
  options.set_include_interiors(true);

  // Considering both shapes we should be within 0.1 degrees of the geometry.
  {
    S2ClosestEdgeQuery query(index.get(), options);

    S2ClosestEdgeQuery::PointTarget target(MakePointOrDie("0:1.5"));
    EXPECT_TRUE(query.IsDistanceLess(&target, S1ChordAngle::Degrees(0.1)));
  }

  // If we only consider shape 0, then we're 0.5 degrees away, so this should
  // be false now.
  {
    S2ClosestEdgeQuery query(index.get(), options);

    S2ClosestEdgeQuery::PointTarget target(MakePointOrDie("0:1.5"));
    EXPECT_FALSE(
        query.IsDistanceLess(&target, S1ChordAngle::Degrees(0.1),
                             [](int shape_id) { return shape_id == 0; }));
  }
}

class VisitClosestEdgesTest : public ::testing::Test {
 public:
  using Options = S2ClosestEdgeQuery::Options;
  using PointTarget = S2ClosestEdgeQuery::PointTarget;
  using Result = S2ClosestEdgeQuery::Result;
  using ResultVisitor = S2ClosestEdgeQuery::ResultVisitor;
  using ShapeFilter = S2ClosestEdgeQuery::ShapeFilter;

  VisitClosestEdgesTest() {
    // Construct a query with an index of simple geometry.
    index_ =
        MakeIndexOrDie("## 1:1, 1:-1, -1:-1, -1:1 | 2:2, 2:-2, -2:-2, -2:2");
    EXPECT_EQ(index_->num_shape_ids(), 2);

    query_.Init(index_.get());
  }

  // Generate a large fractal at (0, 0) and set the query to use it instead.
  int FractalQuery(absl::BitGenRef bitgen,
                   S1Angle radius = S1Angle::Degrees(10)) {
    S2Point z = S2LatLng::FromDegrees(0, 0).ToPoint();
    S2Point x = S2::RobustCrossProd(z, S2Point(0, 0, 1)).Normalize();
    S2Point y = S2::RobustCrossProd(z, x).Normalize();
    auto frame = Matrix3x3_d::FromCols(x, y, z);

    S2Fractal fractal(bitgen);
    fractal.SetLevelForApproxMaxEdges(10000);
    polygon_ = std::make_unique<S2Polygon>(fractal.MakeLoop(frame, radius));
    query_.Init(&polygon_->index());

    return polygon_->num_vertices();
  }

  // Returns the number of edges visited.  If a visitor is given, results are
  // passed to it as well.  The given shape filter (if any) is passed to the
  // query.
  int Visit(S2MinDistanceTarget* target, const Options& options = Options(),
            std::optional<ResultVisitor> visitor = {},
            ShapeFilter filter = {}) {
    int count = 0;
    query_.VisitClosestEdges(
        target, options,
        [&](const Result& result) {
          ++count;
          return !visitor || (*visitor)(result);
        },
        filter);
    return count;
  }

  // A function to use as a visitor that always returns false.
  static bool FalseVisitor(const Result&) { return false; }

 private:
  std::unique_ptr<S2ShapeIndex> index_;
  std::unique_ptr<S2Polygon> polygon_;
  S2ClosestEdgeQuery query_;
};

TEST_F(VisitClosestEdgesTest, CanVisitClosestEdges) {
  // The target point is contained by the second shape but not the first and
  // then there are 8 edges total so we should see 1 + 8 = 9 total results.
  PointTarget target(MakePointOrDie("0:1.5"));
  EXPECT_EQ(Visit(&target), 9);
}

TEST_F(VisitClosestEdgesTest, CanFilterShapes) {
  // Check that we can filter out individual shape ids, or all shapes.
  PointTarget target(MakePointOrDie("0:1.5"));
  EXPECT_EQ(Visit(&target, {}, {}, [](int id) { return id == 0; }), 4);
  EXPECT_EQ(Visit(&target, {}, {}, [](int id) { return id == 1; }), 5);
  EXPECT_EQ(Visit(&target, {}, {}, [](int) { return false; }), 0);
}

TEST_F(VisitClosestEdgesTest, UpdatingShapeFilterWorks) {
  absl::flat_hash_set<int> seen;
  const auto filter = [&](int shape_id) { return !seen.contains(shape_id); };

  // We should be able to filter shapes even while we're visiting.
  PointTarget target(MakePointOrDie("2.5:1.5"));
  EXPECT_EQ(Visit(
                &target, {},
                [&](const Result& result) {
                  seen.insert(result.shape_id());
                  return true;
                },
                filter),
            2);
  EXPECT_EQ(seen.size(), 2);
}

TEST_F(VisitClosestEdgesTest, CanBreakFromShapeIteration) {
  // If we return false immediately we should only see one result.
  PointTarget target(MakePointOrDie("0:0"));
  EXPECT_EQ(Visit(&target, {}, FalseVisitor), 1);
}

TEST_F(VisitClosestEdgesTest, CanBreakFromBruteForce) {
  using Result = S2ClosestEdgeQuery::Result;

  S2ClosestEdgeQuery::Options options;
  options.set_use_brute_force(true);
  options.set_include_interiors(false);

  // If we return false immediately we should only see one result.
  PointTarget target(MakePointOrDie("0:0"));
  EXPECT_EQ(Visit(&target, options, FalseVisitor), 1);
}

TEST_F(VisitClosestEdgesTest, CanBreakFromNormalIteration) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "CAN_BREAK_FROM_NORMAL_ITERATION",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  FractalQuery(bitgen);

  S2ClosestEdgeQuery::Options options;
  options.set_include_interiors(false);

  // If we return false immediately we should only see one result.
  PointTarget target(MakePointOrDie("0:0"));
  EXPECT_EQ(Visit(&target, options, FalseVisitor), 1);
}

TEST_F(VisitClosestEdgesTest, DistanceIsMonotonic) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "DISTANCE_IS_MONOTONIC",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  int num_vertices = FractalQuery(bitgen);

  S2ClosestEdgeQuery::Options options;
  options.set_include_interiors(false);

  PointTarget target(MakePointOrDie("3.14:15.962"));

  // Edge distance should increase monotonically.
  S1ChordAngle last_edge_distance = S1ChordAngle::Zero();
  const int results = Visit(&target, options, [&](const Result& result) {
    EXPECT_GE(result.distance(), last_edge_distance);
    last_edge_distance = result.distance();
    return true;
  });

  // And we should have seen a result for every edge of the fractal.
  EXPECT_EQ(results, num_vertices);
}

TEST_F(VisitClosestEdgesTest, CanLimitByDistance) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "CAN_LIMIT_BY_DISTANCE",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  int num_vertices = FractalQuery(bitgen);

  const S1ChordAngle kDistanceLimit = S1ChordAngle::Degrees(12);

  S2ClosestEdgeQuery::Options options;
  options.set_include_interiors(false);
  options.set_max_distance(kDistanceLimit);

  S1ChordAngle max_edge_distance = S1ChordAngle::Zero();

  PointTarget target(MakePointOrDie("3.14:15.962"));
  const int results = Visit(&target, options, [&](const Result& result) {
    if (result.distance() > max_edge_distance) {
      max_edge_distance = result.distance();
    }
    return true;
  });

  // We shouldn't see every edge of the polygon since we limited by distance.
  EXPECT_LT(results, num_vertices);

  // The maximum result distance we saw should be under the limit.
  EXPECT_LT(max_edge_distance, kDistanceLimit);
}

TEST_F(VisitClosestEdgesTest, CanLimitByNumResults) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "CAN_LIMIT_BY_NUM_RESULTS",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  FractalQuery(bitgen);

  constexpr int kResultLimit = 3141;

  S2ClosestEdgeQuery::Options options;
  options.set_include_interiors(false);
  options.set_max_results(kResultLimit);

  PointTarget target(MakePointOrDie("3.14:15.962"));
  EXPECT_EQ(Visit(&target, options), kResultLimit);
}

TEST(S2ClosestEdgeQuery, TargetPointOutsideIndexedPolygon) {
  // Tests a target point in the interior of a polyline loop with no
  // interior.  (The index also includes a nearby polygon.)
  auto index = MakeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
  S2ClosestEdgeQuery::Options options;
  options.set_include_interiors(true);
  options.set_max_distance(S1Angle::Degrees(1));
  S2ClosestEdgeQuery query(index.get(), options);
  S2ClosestEdgeQuery::PointTarget target(MakePointOrDie("2:2"));
  auto results = query.FindClosestEdges(&target);
  EXPECT_EQ(0, results.size());
}

TEST(S2ClosestEdgeQuery, TargetPolygonContainingIndexedPoints) {
  // Two points are contained within a polyline loop (no interior) and two
  // points are contained within a polygon.
  auto index = MakeIndexOrDie("2:2 | 3:3 | 1:11 | 3:13 # #");
  S2ClosestEdgeQuery query(index.get());
  query.mutable_options()->set_max_distance(S1Angle::Degrees(1));
  auto target_index = MakeIndexOrDie(
      "# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
  S2ClosestEdgeQuery::ShapeIndexTarget target(target_index.get());
  target.set_include_interiors(true);
  auto results = query.FindClosestEdges(&target);
  ASSERT_EQ(2, results.size());

  EXPECT_EQ(S1ChordAngle::Zero(), results[0].distance());
  EXPECT_EQ(0, results[0].shape_id());
  EXPECT_EQ(2, results[0].edge_id());  // 1:11
  EXPECT_FALSE(results[0].is_interior());
  S2Shape::Edge e0 = query.GetEdge(results[0]);
  EXPECT_TRUE(S2::ApproxEquals(e0.v0, S2LatLng::FromDegrees(1, 11).ToPoint()))
      << S2LatLng(e0.v0);
  EXPECT_TRUE(S2::ApproxEquals(e0.v1, S2LatLng::FromDegrees(1, 11).ToPoint()))
      << S2LatLng(e0.v1);

  EXPECT_EQ(S1ChordAngle::Zero(), results[1].distance());
  EXPECT_EQ(0, results[1].shape_id());
  EXPECT_EQ(3, results[1].edge_id());  // 3:13
  EXPECT_FALSE(results[1].is_interior());
  S2Shape::Edge e1 = query.GetEdge(results[1]);
  EXPECT_TRUE(S2::ApproxEquals(e1.v0, S2LatLng::FromDegrees(3, 13).ToPoint()))
      << S2LatLng(e1.v0);
  EXPECT_TRUE(S2::ApproxEquals(e1.v1, S2LatLng::FromDegrees(3, 13).ToPoint()))
      << S2LatLng(e1.v1);
}

TEST(S2ClosestEdgeQuery, EmptyTargetOptimized) {
  // Ensure that the optimized algorithm handles empty targets when a distance
  // limit is specified.
  MutableS2ShapeIndex index;
  index.Add(make_unique<S2Polygon::OwningShape>(make_unique<S2Polygon>(
      S2Loop::MakeRegularLoop(S2Point(1, 0, 0), S1Angle::Radians(0.1), 1000))));
  S2ClosestEdgeQuery query(&index);
  query.mutable_options()->set_max_distance(S1Angle::Radians(1e-5));
  MutableS2ShapeIndex target_index;
  S2ClosestEdgeQuery::ShapeIndexTarget target(&target_index);
  EXPECT_EQ(0, query.FindClosestEdges(&target).size());
}

TEST(S2ClosestEdgeQuery, EmptyPolygonTarget) {
  // Verifies that distances are measured correctly to empty polygon targets.
  auto empty_polygon_index = MakeIndexOrDie("# # empty");
  auto point_index = MakeIndexOrDie("1:1 # #");
  auto full_polygon_index = MakeIndexOrDie("# # full");
  S2ClosestEdgeQuery::ShapeIndexTarget target(empty_polygon_index.get());
  target.set_include_interiors(true);

  S2ClosestEdgeQuery empty_query(empty_polygon_index.get());
  empty_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Infinity(), empty_query.GetDistance(&target));

  S2ClosestEdgeQuery point_query(point_index.get());
  point_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Infinity(), point_query.GetDistance(&target));

  S2ClosestEdgeQuery full_query(full_polygon_index.get());
  full_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Infinity(), full_query.GetDistance(&target));
}

TEST(S2ClosestEdgeQuery, FullLaxPolygonTarget) {
  // Verifies that distances are measured correctly to full LaxPolygon targets.
  auto empty_polygon_index = MakeIndexOrDie("# # empty");
  auto point_index = MakeIndexOrDie("1:1 # #");
  auto full_polygon_index = MakeIndexOrDie("# # full");
  S2ClosestEdgeQuery::ShapeIndexTarget target(full_polygon_index.get());
  target.set_include_interiors(true);

  S2ClosestEdgeQuery empty_query(empty_polygon_index.get());
  empty_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Infinity(), empty_query.GetDistance(&target));

  S2ClosestEdgeQuery point_query(point_index.get());
  point_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Zero(), point_query.GetDistance(&target));

  S2ClosestEdgeQuery full_query(full_polygon_index.get());
  full_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Zero(), full_query.GetDistance(&target));
}

TEST(S2ClosestEdgeQuery, FullS2PolygonTarget) {
  // Verifies that distances are measured correctly to full S2Polygon targets
  // (which use a different representation of "full" than LaxPolygon does).
  auto empty_polygon_index = MakeIndexOrDie("# # empty");
  auto point_index = MakeIndexOrDie("1:1 # #");
  auto full_polygon_index = MakeIndexOrDie("# #");
  full_polygon_index->Add(make_unique<S2Polygon::OwningShape>(
      s2textformat::MakePolygonOrDie("full")));

  S2ClosestEdgeQuery::ShapeIndexTarget target(full_polygon_index.get());
  target.set_include_interiors(true);

  S2ClosestEdgeQuery empty_query(empty_polygon_index.get());
  empty_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Infinity(), empty_query.GetDistance(&target));

  S2ClosestEdgeQuery point_query(point_index.get());
  point_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Zero(), point_query.GetDistance(&target));

  S2ClosestEdgeQuery full_query(full_polygon_index.get());
  full_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Zero(), full_query.GetDistance(&target));
}

// Returns the `z`-sigma confidence bound for `num_trials` binomial trials
// with success probability `p`.  Uses Gaussian approximation.
static double ConfidenceBound(double p, int num_trials, double z) {
  return p * num_trials + z * std::sqrt(p * (1 - p) * num_trials);
}

TEST(S2ClosestEdgeQuery, IsConservativeDistanceLessOrEqual) {
  // Test
  int num_tested = 0;
  int num_conservative_needed = 0;
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "IS_CONSERVATIVE_DISTANCE_LESS_OR_EQUAL",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  constexpr int kNumIters = 10'000;
  for (int iter = 0; iter < kNumIters; ++iter) {
    S2Point x = s2random::Point(bitgen);
    S2Point dir = s2random::Point(bitgen);
    S1Angle r =
        S1Angle::Radians(M_PI * s2random::LogUniform(bitgen, 1e-30, 1.0));
    S2Point y = S2::GetPointOnLine(x, dir, r);
    S1ChordAngle limit(r);
    if (s2pred::CompareDistance(x, y, limit) <= 0) {
      MutableS2ShapeIndex index;
      index.Add(make_unique<S2PointVectorShape>(vector<S2Point>({x})));
      S2ClosestEdgeQuery query(&index);
      S2ClosestEdgeQuery::PointTarget target(y);
      EXPECT_TRUE(query.IsConservativeDistanceLessOrEqual(&target, limit));
      ++num_tested;
      if (!query.IsDistanceLess(&target, limit)) ++num_conservative_needed;
    }
  }
  // Verify that the observed values are within the expected range.
  // The success probabilities were obtained by running 10M iterations;
  // they are human-verified to be "reasonable".  With a 3-sigma threshold,
  // this test should be ~0.6% flaky.
  constexpr double kExpectedTestedFrac = 0.557;
  constexpr double kExpectedConservativeNeededFrac = 0.0280;
  EXPECT_GE(num_tested,
            ConfidenceBound(kExpectedTestedFrac, kNumIters, /*z=*/-3));
  EXPECT_LE(num_tested,
            ConfidenceBound(kExpectedTestedFrac, kNumIters, /*z=*/3));
  EXPECT_GE(
      num_conservative_needed,
      ConfidenceBound(kExpectedConservativeNeededFrac, kNumIters, /*z=*/-3));
  EXPECT_LE(
      num_conservative_needed,
      ConfidenceBound(kExpectedConservativeNeededFrac, kNumIters, /*z=*/3));
}

// The approximate radius of S2Cap from which query edges are chosen.
static const S1Angle kTestCapRadius = S2Testing::KmToAngle(10);

// An approximate bound on the distance measurement error for "reasonable"
// distances (say, less than Pi/2) due to using S1ChordAngle.
static const double kTestChordAngleError = 1e-15;

using TestingResult = pair<S2MinDistance, ShapeEdgeId>;

// Converts to the format required by CheckDistanceResults() in s2testing.h.
vector<TestingResult> ConvertResults(
    absl::Span<const S2ClosestEdgeQuery::Result> results) {
  vector<TestingResult> testing_results;
  for (const auto& result : results) {
    testing_results.push_back(
        TestingResult(result.distance(),
                      ShapeEdgeId(result.shape_id(), result.edge_id())));
  }
  return testing_results;
}

// Use "query" to find the closest edge(s) to the given target.  Verify that
// the results satisfy the search criteria.
static void GetClosestEdges(S2ClosestEdgeQuery::Target* target,
                            S2ClosestEdgeQuery* query,
                            vector<S2ClosestEdgeQuery::Result>* edges,
                            const absl::flat_hash_set<int>& allowed_shapes) {
  if (allowed_shapes.empty()) {
    query->FindClosestEdges(target, edges);
  } else {
    query->FindClosestEdges(target, edges, [&](int shape_id) {
      return allowed_shapes.contains(shape_id);
    });
  }

  EXPECT_LE(edges->size(), query->options().max_results());
  if (query->options().max_distance() ==
      S2ClosestEdgeQuery::Distance::Infinity()) {
    int max_edges = s2shapeutil::CountEdges(query->index());

    // If we have an allowed shape set, we should only see those edges.
    if (!allowed_shapes.empty()) {
      max_edges = 0;
      for (const int id : allowed_shapes) {
        max_edges += query->index().shape(id)->num_edges();
      }
    }

    int min_expected = min(query->options().max_results(), max_edges);

    if (!query->options().include_interiors()) {
      // We can predict exactly how many edges should be returned.
      EXPECT_EQ(min_expected, edges->size());
    } else {
      // All edges should be returned, and possibly some shape interiors.
      EXPECT_LE(min_expected, edges->size());
    }
  }
  for (const auto& edge : *edges) {
    // Check that the edge satisfies the max_distance() condition.
    EXPECT_LT(edge.distance(), query->options().max_distance());
  }
}

static S2ClosestEdgeQuery::Result TestFindClosestEdges(
    S2ClosestEdgeQuery::Target* target, S2ClosestEdgeQuery* query,
    absl::flat_hash_set<int> allowed_shapes = {}) {
  vector<S2ClosestEdgeQuery::Result> expected, actual;
  query->mutable_options()->set_use_brute_force(true);
  GetClosestEdges(target, query, &expected, allowed_shapes);

  // If we were restricted to some set of shape ids, make sure that the results
  // are only from that set.
  if (!allowed_shapes.empty()) {
    for (const auto& result : expected) {
      EXPECT_TRUE(allowed_shapes.contains(result.shape_id()));
    }
  }

  query->mutable_options()->set_use_brute_force(false);
  GetClosestEdges(target, query, &actual, allowed_shapes);
  EXPECT_TRUE(CheckDistanceResults(ConvertResults(expected),
                                   ConvertResults(actual),
                                   query->options().max_results(),
                                   query->options().max_distance(),
                                   query->options().max_error()))
      << "max_results=" << query->options().max_results()
      << ", max_distance=" << query->options().max_distance()
      << ", max_error=" << query->options().max_error();

  if (expected.empty()) return S2ClosestEdgeQuery::Result();

  // We need this assigned to an lvalue so that it sticks around when we call
  // the queries below.
  const auto Filter = [&](int shape_id) {
    return allowed_shapes.contains(shape_id);
  };

  const auto filter =
      allowed_shapes.empty() ? S2ClosestEdgeQuery::ShapeFilter() : Filter;

  // Note that when options.max_error() > 0, expected[0].distance() may not
  // be the minimum distance.  It is never larger by more than max_error(),
  // but the actual value also depends on max_results().
  //
  // Here we verify that GetDistance() and IsDistanceLess() return results
  // that are consistent with the max_error() setting.
  S1ChordAngle max_error = query->options().max_error();
  S1ChordAngle min_distance = expected[0].distance();
  EXPECT_LE(query->GetDistance(target, filter), min_distance + max_error);

  // Test IsDistanceLess().
  EXPECT_FALSE(query->IsDistanceLess(target, min_distance - max_error, filter));
  EXPECT_TRUE(
      query->IsConservativeDistanceLessOrEqual(target, min_distance, filter));

  // Return the closest edge result so that we can also test Project.
  return expected[0];
}

// The running time of this test is proportional to
//    (num_indexes + num_queries) * num_edges.
// (Note that every query is checked using the brute force algorithm.)
static void TestWithIndexFactory(const s2testing::ShapeIndexFactory& factory,
                                 int num_indexes, int num_edges,
                                 int num_queries,
                                 const absl::flat_hash_set<int>& allowed_shapes,
                                 absl::BitGenRef bitgen) {
  // Build a set of MutableS2ShapeIndexes containing the desired geometry.
  vector<S2Cap> index_caps;
  vector<unique_ptr<MutableS2ShapeIndex>> indexes;
  for (int i = 0; i < num_indexes; ++i) {
    index_caps.push_back(S2Cap(s2random::Point(bitgen), kTestCapRadius));
    indexes.push_back(make_unique<MutableS2ShapeIndex>());

    // Add at least two shapes.
    factory.AddEdges(index_caps.back(), num_edges, indexes.back().get());
    factory.AddEdges(index_caps.back(), num_edges, indexes.back().get());
  }
  for (int i = 0; i < num_queries; ++i) {
    int i_index = absl::Uniform(bitgen, 0, num_indexes);
    const S2Cap& index_cap = index_caps[i_index];

    // Choose query points from an area approximately 4x larger than the
    // geometry being tested.
    S1Angle query_radius = 2 * index_cap.GetRadius();
    S2Cap query_cap(index_cap.center(), query_radius);
    S2ClosestEdgeQuery query(indexes[i_index].get());

    // Occasionally we don't set any limit on the number of result edges.
    // (This may return all edges if we also don't set a distance limit.)
    if (absl::Bernoulli(bitgen, 0.8)) {
      query.mutable_options()->set_max_results(absl::Uniform(bitgen, 1, 11));
    }
    // We set a distance limit 2/3 of the time.
    if (absl::Bernoulli(bitgen, 2.0 / 3)) {
      query.mutable_options()->set_max_distance(
          absl::Uniform(bitgen, 0.0, 1.0) * query_radius);
    }
    if (absl::Bernoulli(bitgen, 0.5)) {
      // Choose a maximum error whose logarithm is uniformly distributed over
      // a reasonable range, except that it is sometimes zero.
      query.mutable_options()->set_max_error(S1Angle::Radians(
          s2random::LogUniform(bitgen, 1e-4, 1.0) * query_radius.radians()));
    }
    query.mutable_options()->set_include_interiors(
        absl::Bernoulli(bitgen, 0.5));
    int target_type = absl::Uniform(bitgen, 0, 4);
    if (target_type == 0) {
      // Find the edges closest to a given point.
      S2Point point = s2random::SamplePoint(bitgen, query_cap);
      S2ClosestEdgeQuery::PointTarget target(point);
      auto closest = TestFindClosestEdges(&target, &query, allowed_shapes);
      if (!closest.distance().is_infinity()) {
        // Also test the Project method.
        EXPECT_NEAR(
            closest.distance().ToAngle().radians(),
            S1Angle(point, query.Project(point, closest)).radians(),
            kTestChordAngleError);
      }
    } else if (target_type == 1) {
      // Find the edges closest to a given edge.
      S2Point a = s2random::SamplePoint(bitgen, query_cap);
      S2Point b = s2random::SamplePoint(
          bitgen,
          S2Cap(a, s2random::LogUniform(bitgen, 1e-4, 1.0) * query_radius));
      S2ClosestEdgeQuery::EdgeTarget target(a, b);
      TestFindClosestEdges(&target, &query, allowed_shapes);
    } else if (target_type == 2) {
      // Find the edges closest to a given cell.
      int min_level = S2::kMaxDiag.GetLevelForMaxValue(query_radius.radians());
      int level = absl::Uniform(absl::IntervalClosedClosed, bitgen, min_level,
                                S2CellId::kMaxLevel);
      S2Point a = s2random::SamplePoint(bitgen, query_cap);
      S2Cell cell(S2CellId(a).parent(level));
      S2ClosestEdgeQuery::CellTarget target(cell);
      TestFindClosestEdges(&target, &query, allowed_shapes);
    } else {
      ABSL_DCHECK_EQ(3, target_type);
      // Use another one of the pre-built indexes as the target.
      int j_index = absl::Uniform(bitgen, 0, num_indexes);
      S2ClosestEdgeQuery::ShapeIndexTarget target(indexes[j_index].get());
      target.set_include_interiors(absl::Bernoulli(bitgen, 0.5));
      TestFindClosestEdges(&target, &query, allowed_shapes);
    }
  }
}

static constexpr int kNumIndexes = 50;
static constexpr int kNumEdges = 100;
static constexpr int kNumQueries = 200;

class S2ClosestEdgeQueryShapeTest : public ::testing::TestWithParam<int> {
 public:
  absl::flat_hash_set<int> allowed_shapes() const {
    int shape_id = GetParam();
    if (shape_id < 0) {
      return {};
    }
    return {shape_id};
  }
};

TEST_P(S2ClosestEdgeQueryShapeTest, CircleEdges) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "CIRCLE_EDGES",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  TestWithIndexFactory(s2testing::RegularLoopShapeIndexFactory(), kNumIndexes,
                       kNumEdges, kNumQueries, allowed_shapes(), bitgen);
}

TEST_P(S2ClosestEdgeQueryShapeTest, FractalEdges) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRACTAL_EDGES",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  TestWithIndexFactory(s2testing::FractalLoopShapeIndexFactory(bitgen),
                       kNumIndexes, kNumEdges, kNumQueries, allowed_shapes(),
                       bitgen);
}

TEST_P(S2ClosestEdgeQueryShapeTest, PointCloudEdges) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "POINT_CLOUD_EDGES",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  TestWithIndexFactory(s2testing::PointCloudShapeIndexFactory(bitgen),
                       kNumIndexes, kNumEdges, kNumQueries, allowed_shapes(),
                       bitgen);
}

TEST_P(S2ClosestEdgeQueryShapeTest, ConservativeCellDistanceIsUsed) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "CONSERVATIVE_CELL_DISTANCE_IS_USED",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  // This test will be flaky if max_error() is not properly taken into
  // account when measuring distances to S2ShapeIndex cells.
  TestWithIndexFactory(s2testing::FractalLoopShapeIndexFactory(bitgen), 5, 100,
                       10, allowed_shapes(), bitgen);
}

INSTANTIATE_TEST_SUITE_P(AllowedShapeTests, S2ClosestEdgeQueryShapeTest,
                         testing::Values(-1, 0));

