// Copyright 2015 Google Inc. All Rights Reserved.
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

#include "s2/s2closest_point_query.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "absl/flags/flag.h"
#include "absl/log/absl_check.h"
#include "absl/log/log_streamer.h"
#include "absl/random/bit_gen_ref.h"
#include "absl/random/random.h"
#include "absl/strings/str_cat.h"

#include "s2/mutable_s2shape_index.h"
#include "s2/s1angle.h"
#include "s2/s1chord_angle.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2closest_edge_query_testing.h"
#include "s2/s2closest_point_query_base.h"
#include "s2/s2fractal.h"
#include "s2/s2latlng.h"
#include "s2/s2latlng_rect.h"
#include "s2/s2loop.h"
#include "s2/s2metrics.h"
#include "s2/s2min_distance_targets.h"
#include "s2/s2point.h"
#include "s2/s2point_index.h"
#include "s2/s2pointutil.h"
#include "s2/s2random.h"
#include "s2/s2region.h"
#include "s2/s2testing.h"
#include "s2/util/math/matrix3x3.h"

using std::make_unique;
using std::pair;
using std::unique_ptr;
using std::vector;

using TestIndex = S2PointIndex<int>;
using TestQuery = S2ClosestPointQuery<int>;

TEST(S2ClosestPointQuery, NoPoints) {
  TestIndex index;
  TestQuery query(&index);
  S2ClosestPointQueryPointTarget target(S2Point(1, 0, 0));
  const auto results = query.FindClosestPoints(&target);
  EXPECT_EQ(0, results.size());
}

TEST(S2ClosestPointQuery, ManyDuplicatePoints) {
  static constexpr int kNumPoints = 10000;
  const S2Point kTestPoint(1, 0, 0);
  TestIndex index;
  for (int i = 0; i < kNumPoints; ++i) {
    index.Add(kTestPoint, i);
  }
  TestQuery query(&index);
  S2ClosestPointQueryPointTarget target(kTestPoint);
  const auto results = query.FindClosestPoints(&target);
  EXPECT_EQ(kNumPoints, results.size());
}

TEST(S2ClosestPointQuery, EmptyTargetOptimized) {
  // Ensure that the optimized algorithm handles empty targets when a distance
  // limit is specified.
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "EMPTY_TARGET_OPTIMIZED",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  TestIndex index;
  for (int i = 0; i < 1000; ++i) {
    index.Add(s2random::Point(bitgen), i);
  }
  TestQuery query(&index);
  query.mutable_options()->set_max_distance(S1Angle::Radians(1e-5));
  MutableS2ShapeIndex target_index;
  S2ClosestPointQueryShapeIndexTarget target(&target_index);
  EXPECT_EQ(0, query.FindClosestPoints(&target).size());
}

// An abstract class that adds points to an S2PointIndex for benchmarking.
struct PointIndexFactory {
 public:
  virtual ~PointIndexFactory() = default;

  // Requests that approximately "num_points" points located within the given
  // S2Cap bound should be added to "index".
  virtual void AddPoints(const S2Cap& index_cap, int num_points,
                         absl::BitGenRef bitgen, TestIndex* index) const = 0;
};

// Generates points that are regularly spaced along a circle.  Points along a
// circle are nearly the worst case for distance calculations, since many
// points are nearly equidistant from any query point that is not immediately
// adjacent to the circle.
struct CirclePointIndexFactory : public PointIndexFactory {
  void AddPoints(const S2Cap& index_cap, int num_points, absl::BitGenRef bitgen,
                 TestIndex* index) const override {
    vector<S2Point> points = S2Testing::MakeRegularPoints(
        index_cap.center(), index_cap.GetRadius(), num_points);
    for (int i = 0; i < points.size(); ++i) {
      index->Add(points[i], i);
    }
  }
};

// Generates the vertices of a fractal whose convex hull approximately
// matches the given cap.
struct FractalPointIndexFactory : public PointIndexFactory {
  void AddPoints(const S2Cap& index_cap, int num_points, absl::BitGenRef bitgen,
                 TestIndex* index) const override {
    S2Fractal fractal(bitgen);
    fractal.SetLevelForApproxMaxEdges(num_points);
    fractal.set_fractal_dimension(1.5);
    unique_ptr<S2Loop> loop(fractal.MakeLoop(
        s2random::FrameAt(bitgen, index_cap.center()), index_cap.GetRadius()));
    for (int i = 0; i < loop->num_vertices(); ++i) {
      index->Add(loop->vertex(i), i);
    }
  }
};

// Generates vertices on a square grid that includes the entire given cap.
struct GridPointIndexFactory : public PointIndexFactory {
  void AddPoints(const S2Cap& index_cap, int num_points, absl::BitGenRef bitgen,
                 TestIndex* index) const override {
    int sqrt_num_points = ceil(sqrt(num_points));
    Matrix3x3_d frame = s2random::FrameAt(bitgen, index_cap.center());
    double radius = index_cap.GetRadius().radians();
    double spacing = 2 * radius / sqrt_num_points;
    for (int i = 0; i < sqrt_num_points; ++i) {
      for (int j = 0; j < sqrt_num_points; ++j) {
        S2Point point(tan((i + 0.5) * spacing - radius),
                      tan((j + 0.5) * spacing - radius), 1.0);
        index->Add(S2::FromFrame(frame, point.Normalize()),
                   i * sqrt_num_points + j);
      }
    }
  }
};

// The approximate radius of S2Cap from which query points are chosen.
static const S1Angle kTestCapRadius = S2Testing::KmToAngle(10);

// The result format required by CheckDistanceResults() in s2testing.h.
using TestingResult = pair<S2MinDistance, int>;

// Use "query" to find the closest point(s) to the given target, and extract
// the query results into the given vector.  Also verify that the results
// satisfy the search criteria.
static void GetClosestPoints(TestQuery::Target* target, TestQuery* query,
                             vector<TestingResult>* results) {
  const auto query_results = query->FindClosestPoints(target);
  EXPECT_LE(query_results.size(), query->options().max_results());
  const S2Region* region = query->options().region();
  if (!region && query->options().max_distance() == S1ChordAngle::Infinity()) {
    // We can predict exactly how many points should be returned.
    EXPECT_EQ(std::min(query->options().max_results(),
                       query->index().num_points()),
              query_results.size());
  }
  for (const auto& result : query_results) {
    // Check that the point satisfies the region() condition.
    if (region) EXPECT_TRUE(region->Contains(result.point()));

    // Check that it satisfies the max_distance() condition.
    EXPECT_LT(result.distance(), query->options().max_distance());
    results->push_back(TestingResult(result.distance(), result.data()));
  }
}

static void TestFindClosestPoints(TestQuery::Target* target, TestQuery *query) {
  vector<TestingResult> expected, actual;
  query->mutable_options()->set_use_brute_force(true);
  GetClosestPoints(target, query, &expected);
  query->mutable_options()->set_use_brute_force(false);
  GetClosestPoints(target, query, &actual);
  EXPECT_TRUE(CheckDistanceResults(expected, actual,
                                   query->options().max_results(),
                                   query->options().max_distance(),
                                   query->options().max_error()))
      << "max_results=" << query->options().max_results()
      << ", max_distance=" << query->options().max_distance()
      << ", max_error=" << query->options().max_error();

  if (expected.empty()) return;

  // Note that when options.max_error() > 0, expected[0].distance may not be
  // the minimum distance.  It is never larger by more than max_error(), but
  // the actual value also depends on max_results().
  //
  // Here we verify that GetDistance() and IsDistanceLess() return results
  // that are consistent with the max_error() setting.
  S1ChordAngle max_error = query->options().max_error();
  S1ChordAngle min_distance = expected[0].first;
  EXPECT_LE(query->GetDistance(target), min_distance + max_error);

  // Test IsDistanceLess().
  EXPECT_FALSE(query->IsDistanceLess(target, min_distance - max_error));
  EXPECT_TRUE(query->IsDistanceLessOrEqual(target, min_distance));
  EXPECT_TRUE(query->IsConservativeDistanceLessOrEqual(target, min_distance));
}

// (Note that every query is checked using the brute force algorithm.)
static void TestWithIndexFactory(const PointIndexFactory& factory,
                                 int num_indexes, int num_points,
                                 int num_queries, absl::BitGenRef bitgen) {
  // Build a set of S2PointIndexes containing the desired geometry.
  vector<S2Cap> index_caps;
  vector<unique_ptr<TestIndex>> indexes;
  for (int i = 0; i < num_indexes; ++i) {
    index_caps.push_back(S2Cap(s2random::Point(bitgen), kTestCapRadius));
    indexes.push_back(make_unique<TestIndex>());
    factory.AddPoints(index_caps.back(), num_points, bitgen,
                      indexes.back().get());
  }
  for (int i = 0; i < num_queries; ++i) {
    int i_index = absl::Uniform<size_t>(bitgen, 0, num_indexes);
    const S2Cap& index_cap = index_caps[i_index];

    // Choose query points from an area approximately 4x larger than the
    // geometry being tested.
    S1Angle query_radius = 2 * index_cap.GetRadius();
    S2Cap query_cap(index_cap.center(), query_radius);
    TestQuery query(indexes[i_index].get());

    // Occasionally we don't set any limit on the number of result points.
    // (This may return all points if we also don't set a distance limit.)
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
    S2LatLngRect filter_rect = S2LatLngRect::FromCenterSize(
        S2LatLng(s2random::SamplePoint(bitgen, query_cap)),
        S2LatLng(absl::Uniform(bitgen, 0.0, 1.0) * kTestCapRadius,
                 absl::Uniform(bitgen, 0.0, 1.0) * kTestCapRadius));
    if (absl::Bernoulli(bitgen, 0.2)) {
      query.mutable_options()->set_region(&filter_rect);
    }
    int target_type = absl::Uniform(bitgen, 0, 4);
    if (target_type == 0) {
      // Find the points closest to a given point.
      S2Point point = s2random::SamplePoint(bitgen, query_cap);
      TestQuery::PointTarget target(point);
      TestFindClosestPoints(&target, &query);
    } else if (target_type == 1) {
      // Find the points closest to a given edge.
      S2Point a = s2random::SamplePoint(bitgen, query_cap);
      S2Point b = s2random::SamplePoint(
          bitgen,
          S2Cap(a, s2random::LogUniform(bitgen, 1e-4, 1.0) * query_radius));
      TestQuery::EdgeTarget target(a, b);
      TestFindClosestPoints(&target, &query);
    } else if (target_type == 2) {
      // Find the points closest to a given cell.
      int min_level = S2::kMaxDiag.GetLevelForMaxValue(query_radius.radians());
      int level = absl::Uniform(absl::IntervalClosedClosed, bitgen, min_level,
                                S2CellId::kMaxLevel);
      S2Point a = s2random::SamplePoint(bitgen, query_cap);
      S2Cell cell(S2CellId(a).parent(level));
      TestQuery::CellTarget target(cell);
      TestFindClosestPoints(&target, &query);
    } else {
      ABSL_DCHECK_EQ(3, target_type);
      MutableS2ShapeIndex target_index;
      s2testing::FractalLoopShapeIndexFactory(bitgen).AddEdges(index_cap, 100,
                                                               &target_index);
      TestQuery::ShapeIndexTarget target(&target_index);
      target.set_include_interiors(absl::Bernoulli(bitgen, 0.5));
      TestFindClosestPoints(&target, &query);
    }
  }
}

static constexpr int kNumIndexes = 10;
static constexpr int kNumPoints = 1000;
static constexpr int kNumQueries = 50;

TEST(S2ClosestPointQueryTest, CirclePoints) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "CIRCLE_POINTS",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  TestWithIndexFactory(CirclePointIndexFactory(), kNumIndexes, kNumPoints,
                       kNumQueries, bitgen);
}

TEST(S2ClosestPointQueryTest, FractalPoints) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRACTAL_POINTS",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  TestWithIndexFactory(FractalPointIndexFactory(), kNumIndexes, kNumPoints,
                       kNumQueries, bitgen);
}

TEST(S2ClosestPointQueryTest, GridPoints) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "GRID_POINT",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  TestWithIndexFactory(GridPointIndexFactory(), kNumIndexes, kNumPoints,
                       kNumQueries, bitgen);
}

TEST(S2ClosestPointQueryTest, ConservativeCellDistanceIsUsed) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "CONSERVATIVE_CELL_DISTANCE_IS_USED",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  TestWithIndexFactory(FractalPointIndexFactory(), 5, 100, 10, bitgen);
}

