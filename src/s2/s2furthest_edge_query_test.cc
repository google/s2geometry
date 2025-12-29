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

#include "s2/s2furthest_edge_query.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <benchmark/benchmark.h>
#include <gtest/gtest.h>

#include "absl/flags/flag.h"
#include "absl/log/absl_check.h"
#include "absl/log/log_streamer.h"
#include "absl/random/bit_gen_ref.h"
#include "absl/random/random.h"
#include "absl/strings/str_cat.h"
#include "absl/types/span.h"

#include "s2/mutable_s2shape_index.h"
#include "s2/s1angle.h"
#include "s2/s1chord_angle.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2closest_edge_query_testing.h"
#include "s2/s2edge_distances.h"
#include "s2/s2latlng.h"
#include "s2/s2max_distance_targets.h"
#include "s2/s2metrics.h"
#include "s2/s2point.h"
#include "s2/s2point_vector_shape.h"
#include "s2/s2pointutil.h"
#include "s2/s2polygon.h"
#include "s2/s2predicates.h"
#include "s2/s2random.h"
#include "s2/s2shape.h"
#include "s2/s2shapeutil_count_edges.h"
#include "s2/s2shapeutil_shape_edge_id.h"
#include "s2/s2testing.h"
#include "s2/s2text_format.h"

using absl::StrCat;
using s2textformat::MakeIndexOrDie;
using s2textformat::MakePointOrDie;
using s2textformat::ParsePointsOrDie;
using std::make_pair;
using std::make_unique;
using std::min;
using std::pair;
using std::string;
using std::unique_ptr;
using std::vector;

TEST(S2FurthestEdgeQuery, NoEdges) {
  MutableS2ShapeIndex index;
  S2FurthestEdgeQuery query(&index);
  S2FurthestEdgeQuery::PointTarget target(S2Point(1, 0, 0));
  const auto edge = query.FindFurthestEdge(&target);
  EXPECT_EQ(S1ChordAngle::Negative(), edge.distance());
  EXPECT_EQ(-1, edge.edge_id());
  EXPECT_EQ(-1, edge.shape_id());
  EXPECT_FALSE(edge.is_interior());
  EXPECT_TRUE(edge.is_empty());
  EXPECT_EQ(S1ChordAngle::Negative(), query.GetDistance(&target));
}

TEST(S2FurthestEdgeQuery, OptionsNotModified) {
  // Tests that FindFurthestEdge(), GetDistance(), and IsDistanceGreater() do
  // not modify query.options(), even though all of these methods have their
  // own specific options requirements.
  S2FurthestEdgeQuery::Options options;
  options.set_max_results(3);
  options.set_min_distance(S1ChordAngle::Degrees(1));
  options.set_max_error(S1ChordAngle::Degrees(0.001));
  auto index = MakeIndexOrDie("0:1 | 0:2 | 0:3 # #");
  S2FurthestEdgeQuery query(index.get(), options);
  S2FurthestEdgeQuery::PointTarget target(MakePointOrDie("0:4"));
  EXPECT_EQ(0, query.FindFurthestEdge(&target).edge_id());
  EXPECT_NEAR(3.0, query.GetDistance(&target).degrees(), 1e-15);
  EXPECT_TRUE(query.IsDistanceGreater(&target, S1ChordAngle::Degrees(1.5)));

  // Verify that none of the options above were modified.
  EXPECT_EQ(options.max_results(), query.options().max_results());
  EXPECT_EQ(options.min_distance(), query.options().min_distance());
  EXPECT_EQ(options.max_error(), query.options().max_error());
}

TEST(S2FurthestEdgeQuery, OptionsS1AngleSetters) {
  // Verify that the S1Angle and S1ChordAngle versions do the same thing.
  // This is mainly to prevent the (so far unused) S1Angle versions from
  // being detected as dead code.
  S2FurthestEdgeQuery::Options angle_options, chord_angle_options;
  angle_options.set_min_distance(S1Angle::Degrees(1));
  chord_angle_options.set_min_distance(S1ChordAngle::Degrees(1));
  EXPECT_EQ(chord_angle_options.min_distance(), angle_options.min_distance());

  angle_options.set_inclusive_min_distance(S1Angle::Degrees(1));
  chord_angle_options.set_inclusive_min_distance(S1ChordAngle::Degrees(1));
  EXPECT_EQ(chord_angle_options.min_distance(), angle_options.min_distance());

  angle_options.set_conservative_min_distance(S1Angle::Degrees(1));
  chord_angle_options.set_conservative_min_distance(S1ChordAngle::Degrees(1));
  EXPECT_EQ(chord_angle_options.min_distance(), angle_options.min_distance());
}

// In furthest edge queries, the following distance computation is used when
// updating max distances.
S1ChordAngle GetMaxDistanceToEdge(
    const S2Point& x, const S2Point& y0, const S2Point& y1) {
  S1ChordAngle dist = S1ChordAngle::Negative();
  S2::UpdateMaxDistance(x, y0, y1, &dist);
  return dist;
}

TEST(S2FurthestEdgeQuery, DistanceEqualToLimit) {
  // Tests the behavior of IsDistanceGreater, IsDistanceGreaterOrEqual, and
  // IsConservativeDistanceGreaterOrEqual (and the corresponding Options) when
  // the distance to the target exactly equals the chosen limit.
  S2Point p0(MakePointOrDie("23:12"));
  S2Point p1(MakePointOrDie("47:11"));
  vector<S2Point> index_points{p0};
  MutableS2ShapeIndex index;
  index.Add(make_unique<S2PointVectorShape>(index_points));
  S2FurthestEdgeQuery query(&index);

  // Start with antipodal points and a maximum (180 degrees) distance.
  S2FurthestEdgeQuery::PointTarget target0(-p0);
  S1ChordAngle dist_max = S1ChordAngle::Straight();
  EXPECT_FALSE(query.IsDistanceGreater(&target0, dist_max));
  EXPECT_TRUE(query.IsDistanceGreaterOrEqual(&target0, dist_max));
  EXPECT_TRUE(query.IsConservativeDistanceGreaterOrEqual(&target0, dist_max));

  // Now try two points separated by a non-maximal distance.
  S2FurthestEdgeQuery::PointTarget target1(-p1);
  S1ChordAngle dist1 = GetMaxDistanceToEdge(p0, -p1, -p1);
  EXPECT_FALSE(query.IsDistanceGreater(&target1, dist1));
  EXPECT_TRUE(query.IsDistanceGreaterOrEqual(&target1, dist1));
  EXPECT_TRUE(query.IsConservativeDistanceGreaterOrEqual(&target1, dist1));
}

TEST(S2FurthestEdgeQuery, TrueDistanceGreaterThanS1ChordAngleDistance) {
  // Tests that IsConservativeDistanceGreaterOrEqual returns points where the
  // true distance is slightly greater than the one computed by S1ChordAngle.
  //
  // The points below had the worst error from among 1x10^6 random pairs.
  S2Point p0(0.72362949088190598, -0.39019820403414807, -0.56930283812266336);
  S2Point p1(0.54383822931548842, 0.758981734255934404, 0.35803171284238039);

  // The S1ChordAngle distance is ~3 ulps greater than the true distance.
  S1ChordAngle dist1 = GetMaxDistanceToEdge(p0, p1, p1);
  auto limit = dist1.Successor().Successor().Successor();
  ASSERT_GT(s2pred::CompareDistance(p0, p1, limit), 0);

  // Verify that IsConservativeDistanceGreaterOrEqual() still returns "p1".
  vector<S2Point> index_points{p0};
  MutableS2ShapeIndex index;
  index.Add(make_unique<S2PointVectorShape>(index_points));
  S2FurthestEdgeQuery query(&index);
  S2FurthestEdgeQuery::PointTarget target1(p1);
  EXPECT_FALSE(query.IsDistanceGreater(&target1, limit));
  EXPECT_FALSE(query.IsDistanceGreaterOrEqual(&target1, limit));
  EXPECT_TRUE(query.IsConservativeDistanceGreaterOrEqual(&target1, limit));
}

TEST(S2FurthestEdgeQuery, AntipodalPointInsideIndexedPolygon) {
  // Tests a target point antipodal to the interior of an indexed polygon.
  // (The index also includes a polyline loop with no interior.)
  auto index = MakeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
  S2FurthestEdgeQuery::Options options;

  // First check that with include_interiors set to true, the distance is 180.
  options.set_include_interiors(true);
  options.set_min_distance(S1Angle::Degrees(178));
  S2FurthestEdgeQuery query(index.get(), options);
  S2FurthestEdgeQuery::PointTarget target(-MakePointOrDie("2:12"));
  auto results = query.FindFurthestEdges(&target);
  ASSERT_GT(results.size(), 0);
  EXPECT_EQ(S1ChordAngle::Straight(), results[0].distance());
  // Should find the polygon shape (id = 1).
  EXPECT_EQ(1, results[0].shape_id());
  // Should find the interior, so no specific edge id.
  EXPECT_EQ(-1, results[0].edge_id());
  EXPECT_TRUE(results[0].is_interior());
  EXPECT_FALSE(results[0].is_empty());

  // Next check that with include_interiors set to false, the distance is less
  // than 180 for the same target and index.
  query.mutable_options()->set_include_interiors(false);
  results = query.FindFurthestEdges(&target);
  ASSERT_GT(results.size(), 0);
  EXPECT_LE(results[0].distance(), S1ChordAngle::Straight());
  EXPECT_EQ(1, results[0].shape_id());
  // Found a specific edge, so id should be positive.
  EXPECT_EQ(3, results[0].edge_id());
  EXPECT_FALSE(results[0].is_interior());
  EXPECT_FALSE(results[0].is_empty());
  S2Shape::Edge e0 = query.GetEdge(results[0]);
  EXPECT_TRUE(S2::ApproxEquals(e0.v0, S2LatLng::FromDegrees(5, 10).ToPoint()))
      << S2LatLng(e0.v0);
  EXPECT_TRUE(S2::ApproxEquals(e0.v1, S2LatLng::FromDegrees(0, 10).ToPoint()))
      << S2LatLng(e0.v1);
}

TEST(S2FurthestEdgeQuery, AntipodalPointOutsideIndexedPolygon) {
  // Tests a target point antipodal to the interior of a polyline loop with no
  // interior.  The index also includes a polygon almost antipodal to the
  // target, but with all edges closer than the min_distance threshold.
  auto index = MakeIndexOrDie("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
  S2FurthestEdgeQuery::Options options;
  options.set_include_interiors(true);
  options.set_min_distance(S1Angle::Degrees(179));
  S2FurthestEdgeQuery query(index.get(), options);
  S2FurthestEdgeQuery::PointTarget target(-MakePointOrDie("2:2"));
  auto results = query.FindFurthestEdges(&target);
  EXPECT_EQ(0, results.size());
}

TEST(S2FurthestEdgeQuery, TargetPolygonContainingIndexedPoints) {
  // Two points are contained within a polyline loop (no interior) and two
  // points are contained within a polygon.
  auto index = MakeIndexOrDie("2:2 | 4:4 | 1:11 | 3:12 # #");
  S2FurthestEdgeQuery query(index.get());
  query.mutable_options()->set_use_brute_force(false);
  auto target_index = MakeIndexOrDie(
      "# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
  S2FurthestEdgeQuery::ShapeIndexTarget target(target_index.get());
  target.set_include_interiors(true);
  target.set_use_brute_force(true);
  auto results1 = query.FindFurthestEdges(&target);
  // All points should be returned since we did not specify max_results.
  ASSERT_EQ(4, results1.size());
  EXPECT_NE(S1ChordAngle::Zero(), results1[0].distance());
  EXPECT_EQ(0, results1[0].shape_id());
  EXPECT_EQ(0, results1[0].edge_id());  // 2:2 (to 5:15)
  EXPECT_NE(S1ChordAngle::Zero(), results1[1].distance());
  EXPECT_EQ(0, results1[1].shape_id());
  EXPECT_EQ(3, results1[1].edge_id());  // 3:12 (to 0:0)
}

TEST(S2FurthestEdgeQuery, AntipodalPolygonContainingIndexedPoints) {
  // Two antipodal points are contained within a polyline loop (no interior)
  // and two antipodal points are contained within a polygon.
  auto points = ParsePointsOrDie("2:2, 3:3, 1:11, 3:13");
  auto index = make_unique<MutableS2ShapeIndex>();
  vector<S2Point> antipodal_points;
  for (const auto& p : points) {
    antipodal_points.push_back(-p);
  }
  index->Add(make_unique<S2PointVectorShape>(std::move(antipodal_points)));

  S2FurthestEdgeQuery query(index.get());
  query.mutable_options()->set_min_distance(S1Angle::Degrees(179));
  auto target_index = MakeIndexOrDie(
      "# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
  S2FurthestEdgeQuery::ShapeIndexTarget target(target_index.get());
  target.set_include_interiors(true);
  auto results = query.FindFurthestEdges(&target);
  ASSERT_EQ(2, results.size());
  EXPECT_EQ(S1ChordAngle::Straight(), results[0].distance());
  EXPECT_EQ(0, results[0].shape_id());
  EXPECT_EQ(2, results[0].edge_id());  // 1:11
  EXPECT_EQ(S1ChordAngle::Straight(), results[1].distance());
  EXPECT_EQ(0, results[1].shape_id());
  EXPECT_EQ(3, results[1].edge_id());  // 3:13
}

TEST(S2FurthestEdgeQuery, EmptyPolygonTarget) {
  // Verifies that distances are measured correctly to empty polygon targets.
  auto empty_polygon_index = MakeIndexOrDie("# # empty");
  auto point_index = MakeIndexOrDie("1:1 # #");
  auto full_polygon_index = MakeIndexOrDie("# # full");
  S2FurthestEdgeQuery::ShapeIndexTarget target(empty_polygon_index.get());
  target.set_include_interiors(true);

  S2FurthestEdgeQuery empty_query(empty_polygon_index.get());
  empty_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Negative(), empty_query.GetDistance(&target));

  S2FurthestEdgeQuery point_query(point_index.get());
  point_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Negative(), point_query.GetDistance(&target));

  S2FurthestEdgeQuery full_query(full_polygon_index.get());
  full_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Negative(), full_query.GetDistance(&target));
}

TEST(S2FurthestEdgeQuery, FullLaxPolygonTarget) {
  // Verifies that distances are measured correctly to full LaxPolygon targets.
  auto empty_polygon_index = MakeIndexOrDie("# # empty");
  auto point_index = MakeIndexOrDie("1:1 # #");
  auto full_polygon_index = MakeIndexOrDie("# # full");
  S2FurthestEdgeQuery::ShapeIndexTarget target(full_polygon_index.get());
  target.set_include_interiors(true);

  S2FurthestEdgeQuery empty_query(empty_polygon_index.get());
  empty_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Negative(), empty_query.GetDistance(&target));

  S2FurthestEdgeQuery point_query(point_index.get());
  point_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Straight(), point_query.GetDistance(&target));

  S2FurthestEdgeQuery full_query(full_polygon_index.get());
  full_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Straight(), full_query.GetDistance(&target));
}

TEST(S2FurthestEdgeQuery, FullS2PolygonTarget) {
  // Verifies that distances are measured correctly to full S2Polygon targets
  // (which use a different representation of "full" than LaxPolygon does).
  auto empty_polygon_index = MakeIndexOrDie("# # empty");
  auto point_index = MakeIndexOrDie("1:1 # #");
  auto full_polygon_index = MakeIndexOrDie("# #");
  full_polygon_index->Add(make_unique<S2Polygon::OwningShape>(
      s2textformat::MakePolygonOrDie("full")));

  S2FurthestEdgeQuery::ShapeIndexTarget target(full_polygon_index.get());
  target.set_include_interiors(true);

  S2FurthestEdgeQuery empty_query(empty_polygon_index.get());
  empty_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Negative(), empty_query.GetDistance(&target));

  S2FurthestEdgeQuery point_query(point_index.get());
  point_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Straight(), point_query.GetDistance(&target));

  S2FurthestEdgeQuery full_query(full_polygon_index.get());
  full_query.mutable_options()->set_include_interiors(true);
  EXPECT_EQ(S1ChordAngle::Straight(), full_query.GetDistance(&target));
}

TEST(S2FurthestEdgeQuery, CheckSettings) {
  auto full_polygon_index = MakeIndexOrDie("# #");
  full_polygon_index->Add(make_unique<S2Polygon::OwningShape>(
      s2textformat::MakePolygonOrDie("full")));

  S2FurthestEdgeQuery::ShapeIndexTarget target(full_polygon_index.get());
  target.set_include_interiors(true);
  target.set_use_brute_force(true);

  EXPECT_TRUE(target.include_interiors());
  EXPECT_TRUE(target.use_brute_force());
}

//////////////////////////////////////////////////////////////////////////////
//  General query testing by comparing with brute force method.
//////////////////////////////////////////////////////////////////////////////

// The approximate radius of S2Cap from which query edges are chosen.
static const S1Angle kTestCapRadius = S2Testing::KmToAngle(10);

using Result = pair<S2MaxDistance, s2shapeutil::ShapeEdgeId>;

// Converts to the format required by CheckDistanceResults() in s2testing.h
// TODO(user): When S2ClosestEdgeQuery::Result is made into a class, some
// of the following code may become redundant with that in
// s2closest_edge_query.cc.
vector<Result> ConvertResults(
    absl::Span<const S2FurthestEdgeQuery::Result> edges) {
  vector<Result> results;
  for (const auto& edge : edges) {
    results.push_back(
        make_pair(S2MaxDistance(edge.distance()),
                  s2shapeutil::ShapeEdgeId(edge.shape_id(), edge.edge_id())));
  }

  return results;
}

// Use "query" to find the furthest edge(s) to the given target.  Verify that
// the results satisfy the search criteria.
static void GetFurthestEdges(S2FurthestEdgeQuery::Target* target,
                            S2FurthestEdgeQuery *query,
                            vector<S2FurthestEdgeQuery::Result>* edges) {
  query->FindFurthestEdges(target, edges);
  EXPECT_LE(edges->size(), query->options().max_results());
  if (query->options().min_distance() == S1ChordAngle::Negative()) {
    int min_expected = min(query->options().max_results(),
                           s2shapeutil::CountEdges(query->index()));
    if (!query->options().include_interiors()) {
      // We can predict exactly how many edges should be returned.
      EXPECT_EQ(min_expected, edges->size());
    } else {
      // All edges should be returned, and possibly some shape interiors.
      EXPECT_LE(min_expected, edges->size());
    }
  }
  for (const auto& edge : *edges) {
    // Check that the edge satisfies the min_distance() condition.
    EXPECT_GE(edge.distance(), S1ChordAngle(query->options().min_distance()));
  }
}

static void TestFindFurthestEdges(
    S2FurthestEdgeQuery::Target* target, S2FurthestEdgeQuery *query) {
  vector<S2FurthestEdgeQuery::Result> expected, actual;
  query->mutable_options()->set_use_brute_force(true);
  GetFurthestEdges(target, query, &expected);
  query->mutable_options()->set_use_brute_force(false);
  GetFurthestEdges(target, query, &actual);

  S1ChordAngle min_distance = query->options().min_distance();
  S1ChordAngle max_error = query->options().max_error();
  EXPECT_TRUE(CheckDistanceResults(
      ConvertResults(expected),
      ConvertResults(actual),
      query->options().max_results(),
      S2MaxDistance(min_distance),
      max_error))
      << "max_results=" << query->options().max_results()
      << ", max_distance=" << S1ChordAngle(min_distance)
      << ", max_error=" << max_error;

  if (expected.empty()) {
    return;
  }

  // Note that when options.max_error() > 0, expected[0].distance may not be
  // the maximum distance.  It is never smaller by more than max_error(), but
  // the actual value also depends on max_results().
  //
  // Here we verify that GetDistance() and IsDistanceGreater() return results
  // that are consistent with the max_error() setting.
  S1ChordAngle expected_distance = expected[0].distance();
  EXPECT_GE(query->GetDistance(target), expected_distance - max_error);

  // Test IsDistanceGreater().
  EXPECT_FALSE(query->IsDistanceGreater(
      target, expected_distance + max_error));
  EXPECT_TRUE(query->IsDistanceGreater(
      target, expected[0].distance().Predecessor()));
}

// The running time of this test is proportional to
//    (num_indexes + num_queries) * num_edges.
// (Note that every query is checked using the brute force algorithm.)
static void TestWithIndexFactory(const s2testing::ShapeIndexFactory& factory,
                                 int num_indexes, int num_edges,
                                 int num_queries, absl::BitGenRef bitgen) {
  // Build a set of MutableS2ShapeIndexes containing the desired geometry.
  vector<S2Cap> index_caps;
  vector<unique_ptr<MutableS2ShapeIndex>> indexes;
  for (int i = 0; i < num_indexes; ++i) {
    index_caps.push_back(S2Cap(s2random::Point(bitgen), kTestCapRadius));
    indexes.emplace_back(new MutableS2ShapeIndex);
    factory.AddEdges(index_caps.back(), num_edges, indexes.back().get());
  }

  for (int i = 0; i < num_queries; ++i) {
    int i_index = absl::Uniform(bitgen, 0, num_indexes);
    const S2Cap& index_cap = index_caps[i_index];

    // Choose query points from an area approximately 4x larger than the
    // geometry being tested.
    S1Angle query_radius = 2 * index_cap.GetRadius();
    S2FurthestEdgeQuery query(indexes[i_index].get());

    // Exercise the opposite-hemisphere code 1/5 of the time.
    int antipodal = absl::Bernoulli(bitgen, 0.2) ? -1 : 1;
    S2Cap query_cap(antipodal * index_cap.center(), query_radius);

    // Occasionally we don't set any limit on the number of result edges.
    // (This may return all edges if we also don't set a distance limit.)
    if (absl::Bernoulli(bitgen, 0.8)) {
      query.mutable_options()->set_max_results(absl::Uniform(bitgen, 1, 11));
    }
    // We set a distance limit 2/3 of the time.
    if (!absl::Bernoulli(bitgen, 1.0 / 3)) {
      query.mutable_options()->set_min_distance(
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
      // Find the edges furthest from a given point.
      S2Point point = s2random::SamplePoint(bitgen, query_cap);
      S2FurthestEdgeQuery::PointTarget target(point);
      TestFindFurthestEdges(&target, &query);
    } else if (target_type == 1) {
      // Find the edges furthest from a given edge.
      S2Point a = s2random::SamplePoint(bitgen, query_cap);
      S2Point b = s2random::SamplePoint(
          bitgen,
          S2Cap(a, s2random::LogUniform(bitgen, 1e-4, 1.0) * query_radius));
      S2FurthestEdgeQuery::EdgeTarget target(a, b);
      TestFindFurthestEdges(&target, &query);
    } else if (target_type == 2) {
      // Find the edges furthest from a given cell.
      int min_level = S2::kMaxDiag.GetLevelForMaxValue(query_radius.radians());
      int level = absl::Uniform(absl::IntervalClosedClosed, bitgen, min_level,
                                S2CellId::kMaxLevel);
      S2Point a = s2random::SamplePoint(bitgen, query_cap);
      S2Cell cell(S2CellId(a).parent(level));
      S2FurthestEdgeQuery::CellTarget target(cell);
      TestFindFurthestEdges(&target, &query);
    } else {
      ABSL_DCHECK_EQ(3, target_type);
      // Use another one of the pre-built indexes as the target.
      int j_index = absl::Uniform(bitgen, 0, num_indexes);
      S2FurthestEdgeQuery::ShapeIndexTarget target(indexes[j_index].get());
      target.set_include_interiors(absl::Bernoulli(bitgen, 0.5));
      TestFindFurthestEdges(&target, &query);
    }
  }
}

static constexpr int kNumIndexes = 50;
static constexpr int kNumEdges = 100;
static constexpr int kNumQueries = 200;

TEST(S2FurthestEdgeQuery, CircleEdges) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "CIRCLE_EDGES", absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  TestWithIndexFactory(s2testing::RegularLoopShapeIndexFactory(), kNumIndexes,
                       kNumEdges, kNumQueries, bitgen);
}

TEST(S2FurthestEdgeQuery, FractalEdges) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRACTAL_EDGES", absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  TestWithIndexFactory(s2testing::FractalLoopShapeIndexFactory(bitgen),
                       kNumIndexes, kNumEdges, kNumQueries, bitgen);
}

TEST(S2FurthestEdgeQuery, PointCloudEdges) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "POINT_CLOUD_EDGES", absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  TestWithIndexFactory(s2testing::PointCloudShapeIndexFactory(bitgen),
                       kNumIndexes, kNumEdges, kNumQueries, bitgen);
}


ABSL_FLAG(bool, bm_use_brute_force, false,
          "Benchmarks: Use the brute force implementation");

ABSL_FLAG(double, bm_radius_km, 1000,
          "Benchmarks: Default radius for indexed geometry");

// Calls FindFurthestEdges() the given number of times on a MutableS2ShapeIndex
// with approximately "num_index_edges" edges generated by "factory".  The
// geometry is generated within an S2Cap of the radius specified by
// "FLAGS_bm_radius_km" (the "index radius").
//
// Each query uses a target of the given "target_type".  If "target_type" is
// INDEX, then the target will have approximately "num_target_edges" edges.
//
//   - If min_distance_fraction > 0, then min_distance() is set to the given
//     fraction of the index radius.

//   - If max_error_fraction > 0, then max_error() is set to the given
//     fraction of the index radius.
//
// The remaining parameters are passed to S2ClosestEdgeQueryBenchmarkHarness
// where the geometry for the index and targets is generated.  See there for
// further explanation.
//
static void BenchmarkFindFurthest(
    absl::BitGenRef bitgen, benchmark::State& state,
    const s2testing::ShapeIndexFactory& factory, int num_index_edges,
    bool include_interiors, s2testing::TargetType target_type,
    int num_target_edges, double min_distance_fraction,
    double max_error_fraction, bool choose_target_from_index,
    double target_radius_fraction, double center_separation_fraction) {
  // We execute "state.max_iterations" queries spread out over a total of
  // kNumIndexSamples different geometry samples.  (Performance is affected by
  // how the shapes are positioned relative to the S2Cell hierarchy.)
  static constexpr int kNumIndexSamples = 8;

  MutableS2ShapeIndex index;
  S2FurthestEdgeQuery query(&index);
  query.mutable_options()->set_max_results(1);
  query.mutable_options()->set_include_interiors(include_interiors);
  const S1Angle radius =
      S2Testing::KmToAngle(absl::GetFlag(FLAGS_bm_radius_km));
  if (min_distance_fraction > 0) {
    query.mutable_options()->set_min_distance(min_distance_fraction * radius);
  }
  if (max_error_fraction > 0) {
    query.mutable_options()->set_max_error(max_error_fraction * radius);
  }
  query.mutable_options()->set_use_brute_force(
      absl::GetFlag(FLAGS_bm_use_brute_force));

  int delta = 0;  // Bresenham-type algorithm for geometry sampling.
  vector<unique_ptr<S2FurthestEdgeQuery::Target>> targets;
  vector<unique_ptr<MutableS2ShapeIndex>> target_indexes;

  s2testing::S2ClosestEdgeQueryBenchmarkHarness harness(
      factory, num_index_edges, include_interiors, target_type,
      num_target_edges, choose_target_from_index,
      absl::GetFlag(FLAGS_bm_radius_km), target_radius_fraction,
      center_separation_fraction, bitgen);

  int i_target = 0;
  for (auto s : state) {
    delta -= kNumIndexSamples;
    if (delta < 0) {
      // Generate a new index and a new set of targets to go with it.
      // Reset the random number seed so that we use the same sequence of
      // indexed shapes no matter how many iterations are specified.
      state.PauseTiming();
      delta += state.max_iterations;
      harness.GenerateEdgeQueryWithTargets<S2FurthestEdgeQuery>(
          &query, &index, &targets, &target_indexes);
      state.ResumeTiming();
    }
    query.FindFurthestEdge(targets[i_target].get());
    if (++i_target == targets.size()) i_target = 0;
  }
}

// The maximum number of edges to use in the benchmarks.
// Some benchmarks are run only with this number of edges.
static constexpr int kMaxEdges = 3 * 16384;

namespace {
// Define shorthand names to make the benchmark names less verbose.
using Fractal = s2testing::FractalLoopShapeIndexFactory;
using Regular = s2testing::RegularLoopShapeIndexFactory;
using PointCloud = s2testing::PointCloudShapeIndexFactory;
using TargetType = s2testing::TargetType;
}  // namespace

// Test searching within the general vicinity of the indexed shapes.
template <class Factory>
static void BM_FindFurthest(benchmark::State& state) {
  const string seed_str = StrCat(__func__, absl::GetFlag(FLAGS_s2_random_seed));
  std::seed_seq seed(seed_str.begin(), seed_str.end());
  std::mt19937_64 bitgen(seed);
  int num_edges = state.range(0);
  BenchmarkFindFurthest(bitgen, state, Factory(bitgen), num_edges, false,
                        TargetType::POINT, 0, -1, -1, false, 0.0, -2.0);
}

// As above, but include interiors in the query.
template <class Factory>
static void BM_FindFurthestInterior(benchmark::State& state) {
  const string seed_str = StrCat(__func__, absl::GetFlag(FLAGS_s2_random_seed));
  std::seed_seq seed(seed_str.begin(), seed_str.end());
  std::mt19937_64 bitgen(seed);
  int num_edges = state.range(0);
  BenchmarkFindFurthest(bitgen, state, Factory(bitgen), num_edges, true,
                        TargetType::POINT, 0, -1, -1, false, 0.0, -2.0);
}

// Test searching with an error tolerance.  Allowing 1% error makes searches
// faster.
template <class Factory>
static void BM_FindFurthestMaxErrorPct(benchmark::State& state) {
  const string seed_str = StrCat(__func__, absl::GetFlag(FLAGS_s2_random_seed));
  std::seed_seq seed(seed_str.begin(), seed_str.end());
  std::mt19937_64 bitgen(seed);
  int num_edges = state.range(0);
  int max_error_percent = state.range(1);
  BenchmarkFindFurthest(bitgen, state, Factory(bitgen), num_edges, false,
                        TargetType::POINT, 0, -1, 0.01 * max_error_percent,
                        false, 0.0, -2.0);
}

// Repeat the benchmarks for edge targets rather than point targets.
// (We don't measure searching with a small distance limit because it
// is only beneficial when the edge is small compared to the spacing of the
// index points, which is not the case here.)

template <class Factory>
static void BM_FindFurthestToEdge(benchmark::State& state) {
  const string seed_str = StrCat(__func__, absl::GetFlag(FLAGS_s2_random_seed));
  std::seed_seq seed(seed_str.begin(), seed_str.end());
  std::mt19937_64 bitgen(seed);
  int num_edges = state.range(0);
  BenchmarkFindFurthest(bitgen, state, Factory(bitgen), num_edges, false,
                        TargetType::EDGE, 0, -1, -1, false, -1.0, -2.0);
}

template <class Factory>
static void BM_FindFurthestToEdgeInterior(benchmark::State& state) {
  const string seed_str = StrCat(__func__, absl::GetFlag(FLAGS_s2_random_seed));
  std::seed_seq seed(seed_str.begin(), seed_str.end());
  std::mt19937_64 bitgen(seed);
  int num_edges = state.range(0);
  BenchmarkFindFurthest(bitgen, state, Factory(bitgen), num_edges, true,
                        TargetType::EDGE, 0, -1, -1, false, -1.0, -2.0);
}

// Repeat the benchmarks for S2Cell targets.

template <class Factory>
static void BM_FindFurthestToCell(benchmark::State& state) {
  const string seed_str = StrCat(__func__, absl::GetFlag(FLAGS_s2_random_seed));
  std::seed_seq seed(seed_str.begin(), seed_str.end());
  std::mt19937_64 bitgen(seed);
  int num_edges = state.range(0);
  BenchmarkFindFurthest(bitgen, state, Factory(bitgen), num_edges, false,
                        TargetType::CELL, 0, -1, -1, false, -1.0, -2.0);
}

template <class Factory>
static void BM_FindFurthestToCellInterior(benchmark::State& state) {
  const string seed_str = StrCat(__func__, absl::GetFlag(FLAGS_s2_random_seed));
  std::seed_seq seed(seed_str.begin(), seed_str.end());
  std::mt19937_64 bitgen(seed);
  int num_edges = state.range(0);
  BenchmarkFindFurthest(bitgen, state, Factory(bitgen), num_edges, true,
                        TargetType::CELL, 0, -1, -1, false, -1.0, -2.0);
}

// Repeat the benchmarks for MutableS2ShapeIndex targets.

template <class Factory>
static void BM_FindFurthestToSmallAbuttingIndex(benchmark::State& state) {
  const string seed_str = StrCat(__func__, absl::GetFlag(FLAGS_s2_random_seed));
  std::seed_seq seed(seed_str.begin(), seed_str.end());
  std::mt19937_64 bitgen(seed);
  int num_edges = state.range(0);
  // Measure the distance from an S2ShapeIndex with "num_edges" edges to a
  // small S2ShapeIndex target such that their bounding S2Caps touch.
  BenchmarkFindFurthest(bitgen, state, Factory(bitgen), num_edges, false,
                        TargetType::INDEX, 3 * 4, -1, -1, false, 1.0, 2.0);
}

template <class Factory>
static void BM_FindFurthestFromSmallAbuttingIndex(benchmark::State& state) {
  const string seed_str = StrCat(__func__, absl::GetFlag(FLAGS_s2_random_seed));
  std::seed_seq seed(seed_str.begin(), seed_str.end());
  std::mt19937_64 bitgen(seed);
  int num_edges = state.range(0);
  // Like the benchmark above, except that the two indexes are reversed.  (It
  // turns out that it is somewhat faster to use the bigger of the two
  // S2ShapeIndexes as the target.)
  BenchmarkFindFurthest(bitgen, state, Factory(bitgen), 3 * 4, false,
                        TargetType::INDEX, num_edges, -1, -1, false, 1.0, 2.0);
}

template <class Factory>
static void BM_FindFurthestToSameSizeAbuttingIndex(benchmark::State& state) {
  const string seed_str = StrCat(__func__, absl::GetFlag(FLAGS_s2_random_seed));
  std::seed_seq seed(seed_str.begin(), seed_str.end());
  std::mt19937_64 bitgen(seed);
  int num_edges = state.range(0);
  // Measure the distance between two indexes with the same number of edges
  // and whose bounding S2Caps touch each other.
  BenchmarkFindFurthest(bitgen, state, Factory(bitgen), num_edges, false,
                        TargetType::INDEX, num_edges, -1, -1, false, 1.0, 2.0);
}

template <class Factory>
static void BM_FindFurthestToSameSizeDistantIndex(benchmark::State& state) {
  const string seed_str = StrCat(__func__, absl::GetFlag(FLAGS_s2_random_seed));
  std::seed_seq seed(seed_str.begin(), seed_str.end());
  std::mt19937_64 bitgen(seed);
  int num_edges = state.range(0);
  // Measure the distance between two indexes with the same number of edges
  // and whose bounding S2Caps are distant relative to their radius.
  BenchmarkFindFurthest(bitgen, state, Factory(bitgen), num_edges, false,
                        TargetType::INDEX, num_edges, -1, -1, false, 1.0, 10.0);
}

template <class Factory>
static void BM_IsDistanceGreaterSameSizeDistantIndexFalse(
    benchmark::State& state) {
  const string seed_str = StrCat(__func__, absl::GetFlag(FLAGS_s2_random_seed));
  std::seed_seq seed(seed_str.begin(), seed_str.end());
  std::mt19937_64 bitgen(seed);
  int num_edges = state.range(0);
  // Like the test above, but instead just compare the maximum distance to a
  // threshold that is slightly greater than the maximum possible value.  (The
  // centers of the two geometries are 10 radii apart, so the distance is at
  // most 12 radii.  We set the threshold slightly higher than this.)
  BenchmarkFindFurthest(bitgen, state, Factory(bitgen), num_edges, false,
                        TargetType::INDEX, num_edges, 12.4, 12.4, false, 1.0,
                        10.0);
}

template <class Factory>
static void BM_IsDistanceGreaterSameSizeDistantIndexTrue(
    benchmark::State& state) {
  const string seed_str = StrCat(__func__, absl::GetFlag(FLAGS_s2_random_seed));
  std::seed_seq seed(seed_str.begin(), seed_str.end());
  std::mt19937_64 bitgen(seed);
  int num_edges = state.range(0);
  // Like the test above, but instead just compare the maximum distance to a
  // threshold that is usually slightly less than the true value.  (The
  // centers of the two geometries are 10 radii apart, and by experimentation
  // the maximum distance between the two fractals is at most about 12 radii.
  // We set the threshold slightly lower than this.)
  BenchmarkFindFurthest(bitgen, state, Factory(bitgen), num_edges, false,
                        TargetType::INDEX, num_edges, 11.6, 11.6, false, 1.0,
                        10.0);
}

BENCHMARK_TEMPLATE(BM_FindFurthest, Fractal)
->Arg(3*4)->Arg(3*16)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestInterior, Fractal)
->Arg(3*4)->Arg(3*16)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestMaxErrorPct, Fractal)
->ArgPair(kMaxEdges, 1)->ArgPair(kMaxEdges, 10);

BENCHMARK_TEMPLATE(BM_FindFurthestToEdge, Fractal)
->Arg(3*4)->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToEdgeInterior, Fractal)
->Arg(3*4)->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToCell, Fractal)
->Arg(3*4)->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToCellInterior, Fractal)
->Arg(3*4)->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToSmallAbuttingIndex, Fractal)
->Arg(3*16)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestFromSmallAbuttingIndex, Fractal)
->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToSameSizeAbuttingIndex, Fractal)
->Arg(3*4)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToSameSizeDistantIndex, Fractal)
->Arg(3*4)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_IsDistanceGreaterSameSizeDistantIndexFalse, Fractal)
->Arg(3*4)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_IsDistanceGreaterSameSizeDistantIndexTrue, Fractal)
->Arg(3*4)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

// Now repeat all the benchmarks for regular loops rather than fractal loops.
// We group the benchmarks together by the type of geometry so that it's
// easier to see what effect the various options have (max_distance, etc).

BENCHMARK_TEMPLATE(BM_FindFurthest, Regular)
->Arg(3*4)->Arg(3*16)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestInterior, Regular)
->Arg(3*4)->Arg(3*16)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestMaxErrorPct, Regular)
->ArgPair(kMaxEdges, 1)->ArgPair(kMaxEdges, 10);

BENCHMARK_TEMPLATE(BM_FindFurthestToEdge, Regular)
->Arg(3*4)->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToEdgeInterior, Regular)
->Arg(3*4)->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToCell, Regular)
->Arg(3*4)->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToCellInterior, Regular)
->Arg(3*4)->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToSmallAbuttingIndex, Regular)
->Arg(3*16)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestFromSmallAbuttingIndex, Regular)
->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToSameSizeAbuttingIndex, Regular)
->Arg(3*4)->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToSameSizeDistantIndex, Regular)
->Arg(3*4)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_IsDistanceGreaterSameSizeDistantIndexFalse, Regular)
->Arg(3*4)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_IsDistanceGreaterSameSizeDistantIndexTrue, Regular)
->Arg(3*4)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

// And again with point clouds.

BENCHMARK_TEMPLATE(BM_FindFurthest, PointCloud)
->Arg(3*4)->Arg(3*16)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestInterior, PointCloud)
->Arg(3*4)->Arg(3*16)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestMaxErrorPct, PointCloud)
->ArgPair(kMaxEdges, 1)->ArgPair(kMaxEdges, 10);

BENCHMARK_TEMPLATE(BM_FindFurthestToEdge, PointCloud)
->Arg(3*4)->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToEdgeInterior, PointCloud)
->Arg(3*4)->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToCell, PointCloud)
->Arg(3*4)->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToCellInterior, PointCloud)
->Arg(3*4)->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToSmallAbuttingIndex, PointCloud)
->Arg(3*16)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestFromSmallAbuttingIndex, PointCloud)
->Arg(3*256)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToSameSizeAbuttingIndex, PointCloud)
->Arg(3*4)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_FindFurthestToSameSizeDistantIndex, PointCloud)
->Arg(3*4)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_IsDistanceGreaterSameSizeDistantIndexFalse, PointCloud)
->Arg(3*4)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);

BENCHMARK_TEMPLATE(BM_IsDistanceGreaterSameSizeDistantIndexTrue, PointCloud)
->Arg(3*4)->Arg(3*256)->Arg(3*4096)->Arg(kMaxEdges);
