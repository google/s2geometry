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

#include "s2/s2closestedgequery.h"

#include <memory>
#include <set>
#include <vector>

#include <gtest/gtest.h>
#include "s2/third_party/absl/memory/memory.h"
#include "s2/s1angle.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cellid.h"
#include "s2/s2edge_distances.h"
#include "s2/s2loop.h"
#include "s2/s2metrics.h"
#include "s2/s2shapeutil.h"
#include "s2/s2testing.h"
#include "s2/s2textformat.h"

using absl::make_unique;
using s2textformat::MakeIndex;
using s2textformat::MakePoint;
using std::fabs;
using std::make_pair;
using std::min;
using std::ostream;
using std::pair;
using std::unique_ptr;
using std::vector;

TEST(PointTarget, GetContainingShapes) {
  // Only shapes 2 and 4 should contain the target point.
  auto index = MakeIndex(
      "1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | 0:0, 0:4, 4:0");
  S2ClosestEdgeQuery::PointTarget target(MakePoint("1:1"));
  EXPECT_EQ((vector<int>{2}), target.GetContainingShapes(*index, 1));
  EXPECT_EQ((vector<int>{2, 4}), target.GetContainingShapes(*index, 5));
}

TEST(EdgeTarget, GetContainingShapes) {
  // Only shapes 2 and 4 should contain the target point.
  auto index = MakeIndex(
      "1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | 0:0, 0:4, 4:0");
  S2ClosestEdgeQuery::EdgeTarget target(MakePoint("1:2"), MakePoint("2:1"));
  EXPECT_EQ((vector<int>{2}), target.GetContainingShapes(*index, 1));
  EXPECT_EQ((vector<int>{2, 4}), target.GetContainingShapes(*index, 5));
}

TEST(CellTarget, GetContainingShapes) {
  auto index = MakeIndex(
      "1:1 # 1:1, 2:2 # 0:0, 0:3, 3:0 | 6:6, 6:9, 9:6 | -1:-1, -1:5, 5:-1");
  // Only shapes 2 and 4 should contain a very small cell near 1:1.
  S2CellId cellid1(MakePoint("1:1"));
  S2ClosestEdgeQuery::CellTarget target1{S2Cell(cellid1)};
  EXPECT_EQ((vector<int>{2}), target1.GetContainingShapes(*index, 1));
  EXPECT_EQ((vector<int>{2, 4}), target1.GetContainingShapes(*index, 5));

  // For a larger cell that properly contains one or more index cells, all
  // shapes that intersect the first such cell in S2CellId order are returned.
  // In the test below, this happens to again be the 1st and 3rd polygons
  // (whose shape_ids are 2 and 4).
  S2CellId cellid2 = cellid1.parent(5);
  S2ClosestEdgeQuery::CellTarget target2{S2Cell(cellid2)};
  EXPECT_EQ((vector<int>{2, 4}), target2.GetContainingShapes(*index, 5));
}

TEST(ShapeIndexTarget, GetContainingShapes) {
  // Create an index containing a repeated grouping of one point, one
  // polyline, and one polygon.
  auto index = MakeIndex(
      "1:1 | 4:4 | 7:7 | 10:10 # "
      "1:1, 1:2 | 4:4, 4:5 | 7:7, 7:8 | 10:10, 10:11 # "
      "0:0, 0:3, 3:0 | 3:3, 3:6, 6:3 | 6:6, 6:9, 9:6 | 9:9, 9:12, 12:9");

  // Construct a target consisting of one point, one polyline, and one polygon
  // with two loops where only the second loop is contained by a polygon in
  // the index above.
  auto target_index = MakeIndex(
      "1:1 # 4:5, 5:4 # 20:20, 20:21, 21:20; 10:10, 10:11, 11:10");

  S2ClosestEdgeQuery::ShapeIndexTarget target(target_index.get());
  // These are the shape_ids of the 1st, 2nd, and 4th polygons of "index"
  // (noting that the 4 points are represented by one PointVectorShape).
  EXPECT_EQ((vector<int>{5, 6, 8}), target.GetContainingShapes(*index, 5));
}

TEST(S2ClosestEdgeQuery, NoEdges) {
  S2ShapeIndex index;
  S2ClosestEdgeQuery query(&index);
  query.mutable_options()->set_max_edges(1);
  S2ClosestEdgeQuery::PointTarget target(S2Point(1, 0, 0));
  auto const edge = query.FindClosestEdge(&target);
  EXPECT_EQ(S1ChordAngle::Infinity(), edge.distance);
  EXPECT_EQ(-1, edge.edge_id);
  EXPECT_EQ(-1, edge.shape_id);
  EXPECT_EQ(S1ChordAngle::Infinity(), query.GetDistance(&target));
}

TEST(S2ClosestEdgeQuery, TargetPointInsideIndexedPolygon) {
  // Tests a target point in the interior of an indexed polygon.
  // (The index also includes a polyline loop with no interior.)
  auto index = MakeIndex("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
  S2ClosestEdgeQuery::Options options;
  options.set_include_interiors(true);
  options.set_max_distance(S1Angle::Degrees(1));
  S2ClosestEdgeQuery query(index.get(), options);
  S2ClosestEdgeQuery::PointTarget target(MakePoint("2:12"));
  auto results = query.FindClosestEdges(&target);
  ASSERT_EQ(1, results.size());
  EXPECT_EQ(S1ChordAngle::Zero(), results[0].distance);
  EXPECT_EQ(1, results[0].shape_id);
  EXPECT_EQ(-1, results[0].edge_id);
}

TEST(S2ClosestEdgeQuery, TargetPointOutsideIndexedPolygon) {
  // Tests a target point in the interior of a polyline loop with no
  // interior.  (The index also includes a nearby polygon.)
  auto index = MakeIndex("# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
  S2ClosestEdgeQuery::Options options;
  options.set_include_interiors(true);
  options.set_max_distance(S1Angle::Degrees(1));
  S2ClosestEdgeQuery query(index.get(), options);
  S2ClosestEdgeQuery::PointTarget target(MakePoint("2:2"));
  auto results = query.FindClosestEdges(&target);
  EXPECT_EQ(0, results.size());
}

TEST(S2ClosestEdgeQuery, TargetPolygonContainingIndexedPoints) {
  // Two points are contained within a polyline loop (no interior) and two
  // points are contained within a polygon.
  auto index = MakeIndex("2:2 | 3:3 | 1:11 | 3:13 # #");
  S2ClosestEdgeQuery query(index.get());
  query.mutable_options()->set_max_distance(S1Angle::Degrees(1));
  auto target_index = MakeIndex(
      "# 0:0, 0:5, 5:5, 5:0 # 0:10, 0:15, 5:15, 5:10");
  S2ClosestEdgeQuery::ShapeIndexTarget target(target_index.get());
  target.set_include_interiors(true);
  auto results = query.FindClosestEdges(&target);
  ASSERT_EQ(2, results.size());
  EXPECT_EQ(S1ChordAngle::Zero(), results[0].distance);
  EXPECT_EQ(0, results[0].shape_id);
  EXPECT_EQ(2, results[0].edge_id);  // 1:11
  EXPECT_EQ(S1ChordAngle::Zero(), results[1].distance);
  EXPECT_EQ(0, results[1].shape_id);
  EXPECT_EQ(3, results[1].edge_id);  // 3:13
}

// An abstract class that adds edges to an S2ShapeIndex for benchmarking.
class ShapeIndexFactory {
 public:
  virtual ~ShapeIndexFactory() {}

  // Requests that approximately "num_edges" edges located within the given
  // S2Cap bound should be added to "index".
  virtual void AddEdges(S2Cap const& index_cap, int num_edges,
                        S2ShapeIndex* index) const = 0;
};

// Generates a regular loop that approximately fills the given S2Cap.
//
// Regular loops are nearly the worst case for distance calculations, since
// many edges are nearly equidistant from any query point that is not
// immediately adjacent to the loop.
class RegularLoopShapeIndexFactory : public ShapeIndexFactory {
 public:
  void AddEdges(S2Cap const& index_cap, int num_edges,
                S2ShapeIndex* index) const override {
    index->Add(make_unique<S2Loop::OwningShape>(S2Loop::MakeRegularLoop(
        index_cap.center(), index_cap.GetRadius(), num_edges)));
  }
};

// Generates a fractal loop that approximately fills the given S2Cap.
class FractalLoopShapeIndexFactory : public ShapeIndexFactory {
 public:
  void AddEdges(S2Cap const& index_cap, int num_edges,
                S2ShapeIndex* index) const override {
    S2Testing::Fractal fractal;
    fractal.SetLevelForApproxMaxEdges(num_edges);
    index->Add(make_unique<S2Loop::OwningShape>(
        fractal.MakeLoop(S2Testing::GetRandomFrameAt(index_cap.center()),
                         index_cap.GetRadius())));
  }
};

// Generates a cloud of points that approximately fills the given S2Cap.
class PointCloudShapeIndexFactory : public ShapeIndexFactory {
 public:
  void AddEdges(S2Cap const& index_cap, int num_edges,
                S2ShapeIndex* index) const override {
    vector<S2Point> points;
    for (int i = 0; i < num_edges; ++i) {
      points.push_back(S2Testing::SamplePoint(index_cap));
    }
    index->Add(make_unique<s2shapeutil::PointVectorShape>(&points));
  }
};

// The approximate radius of S2Cap from which query edges are chosen.
static S1Angle const kRadius = S2Testing::KmToAngle(10);

// An approximate bound on the distance measurement error for "reasonable"
// distances (say, less than Pi/2) due to using S1ChordAngle.
static double kChordAngleError = 1e-15;

using Result = pair<S1Angle, s2shapeutil::ShapeEdgeId>;

// Converts to the format required by CheckDistanceResults() in s2testing.h
vector<Result> ConvertResults(vector<S2ClosestEdgeQuery::Result> const& edges) {
  vector<Result> results;
  for (auto const& edge : edges) {
    results.push_back(
        make_pair(edge.distance.ToAngle(),
                  s2shapeutil::ShapeEdgeId(edge.shape_id, edge.edge_id)));
  }
  return results;
}

// Use "query" to find the closest edge(s) to the given target, then convert
// the query results into two parallel vectors, one for distances and one for
// (shape_id, edge_id) pairs.  Also verify that the results satisfy the search
// criteria.
static void GetClosestEdges(S2ClosestEdgeQuery::Target* target,
                            S2ClosestEdgeQuery *query,
                            vector<S2ClosestEdgeQuery::Result>* edges) {
  query->FindClosestEdges(target, edges);
  EXPECT_LE(edges->size(), query->options().max_edges());
  if (query->options().max_distance() ==
      S2ClosestEdgeQuery::Distance::Infinity()) {
    int min_expected = min(query->options().max_edges(),
                           s2shapeutil::GetNumEdges(query->index()));
    if (!query->options().include_interiors()) {
      // We can predict exactly how many edges should be returned.
      EXPECT_EQ(min_expected, edges->size());
    } else {
      // All edges should be returned, and possibly some shape interiors.
      EXPECT_LE(min_expected, edges->size());
    }
  }
  for (auto const& edge : *edges) {
    // Check that the edge satisfies the max_distance() condition.
    EXPECT_LE(edge.distance, query->options().max_distance());
  }
}

static S2ClosestEdgeQuery::Result TestFindClosestEdges(
    S2ClosestEdgeQuery::Target* target, S2ClosestEdgeQuery *query) {
  vector<S2ClosestEdgeQuery::Result> expected, actual;
  query->mutable_options()->set_use_brute_force(true);
  GetClosestEdges(target, query, &expected);
  query->mutable_options()->set_use_brute_force(false);
  GetClosestEdges(target, query, &actual);
  EXPECT_TRUE(CheckDistanceResults(ConvertResults(expected),
                                   ConvertResults(actual),
                                   query->options().max_edges(),
                                   query->options().max_distance().ToAngle(),
                                   query->options().max_error().ToAngle()))
      << "max_edges=" << query->options().max_edges()
      << ", max_distance=" << query->options().max_distance()
      << ", max_error=" << query->options().max_error();

  if (expected.empty()) return S2ClosestEdgeQuery::Result();

  // Test GetDistance(), which should yield the same minimum distance as
  // GetClosestEdges() when max_error() == 0, but can yield a different result
  // when max_error() > 0 because it also sets max_edges() == 1.
  S1ChordAngle min_distance = expected[0].distance;
  S1ChordAngle get_distance = query->GetDistance(target);
  EXPECT_LE(get_distance, min_distance + query->options().max_error());

  // Return the closest edge result so that we can also test Project.
  return expected[0];
}

// The running time of this test is proportional to
//    (num_indexes + num_queries) * num_edges.
// (Note that every query is checked using the brute force algorithm.)
static void TestWithIndexFactory(ShapeIndexFactory const& factory,
                                 int num_indexes, int num_edges,
                                 int num_queries) {
  // Build a set of S2ShapeIndexes containing the desired geometry.
  vector<S2Cap> index_caps;
  vector<unique_ptr<S2ShapeIndex>> indexes;
  for (int i = 0; i < num_indexes; ++i) {
    S2Testing::rnd.Reset(FLAGS_s2_random_seed + i);
    index_caps.push_back(S2Cap(S2Testing::RandomPoint(), kRadius));
    indexes.emplace_back(new S2ShapeIndex);
    factory.AddEdges(index_caps.back(), num_edges, indexes.back().get());
  }
  for (int i = 0; i < num_queries; ++i) {
    S2Testing::rnd.Reset(FLAGS_s2_random_seed + i);
    int i_index = S2Testing::rnd.Uniform(num_indexes);
    S2Cap const& index_cap = index_caps[i_index];

    // Choose query points from an area approximately 4x larger than the
    // geometry being tested.
    S1Angle query_radius = 2 * index_cap.GetRadius();
    S2Cap query_cap(index_cap.center(), query_radius);
    S2ClosestEdgeQuery query(indexes[i_index].get());

    // Occasionally we don't set any limit on the number of result edges.
    // (This may return all edges if we also don't set a distance limit.)
    if (!S2Testing::rnd.OneIn(5)) {
      query.mutable_options()->set_max_edges(1 + S2Testing::rnd.Uniform(10));
    }
    // We set a distance limit 2/3 of the time.
    if (!S2Testing::rnd.OneIn(3)) {
      query.mutable_options()->set_max_distance(
          S2Testing::rnd.RandDouble() * query_radius);
    }
    if (S2Testing::rnd.OneIn(2)) {
      // Choose a maximum error whose logarithm is uniformly distributed over
      // a reasonable range, except that it is sometimes zero.
      query.mutable_options()->set_max_error(S1Angle::Radians(
          pow(1e-4, S2Testing::rnd.RandDouble()) * query_radius.radians()));
    }
    query.mutable_options()->set_include_interiors(S2Testing::rnd.OneIn(2));
    int target_type = S2Testing::rnd.Uniform(4);
    if (target_type == 0) {
      // Find the edges closest to a given point.
      S2Point point = S2Testing::SamplePoint(query_cap);
      S2ClosestEdgeQuery::PointTarget target(point);
      auto closest = TestFindClosestEdges(&target, &query);
      if (!closest.distance.is_infinity()) {
        // Also test the Project method.
        EXPECT_NEAR(
            closest.distance.ToAngle().radians(),
            S1Angle(point, query.Project(point, closest)).radians(),
            kChordAngleError);
      }
    } else if (target_type == 1) {
      // Find the edges closest to a given edge.
      S2Point a = S2Testing::SamplePoint(query_cap);
      S2Point b = S2Testing::SamplePoint(
          S2Cap(a, pow(1e-4, S2Testing::rnd.RandDouble()) * query_radius));
      S2ClosestEdgeQuery::EdgeTarget target(a, b);
      TestFindClosestEdges(&target, &query);
    } else if (target_type == 2) {
      // Find the edges closest to a given cell.
      int min_level = S2::kMaxDiag.GetLevelForMaxValue(query_radius.radians());
      int level = min_level + S2Testing::rnd.Uniform(
          S2CellId::kMaxLevel - min_level + 1);
      S2Point a = S2Testing::SamplePoint(query_cap);
      S2Cell cell(S2CellId(a).parent(level));
      S2ClosestEdgeQuery::CellTarget target(cell);
      TestFindClosestEdges(&target, &query);
    } else {
      DCHECK_EQ(3, target_type);
      // Use another one of the pre-built indexes as the target.
      int j_index = S2Testing::rnd.Uniform(num_indexes);
      S2ClosestEdgeQuery::ShapeIndexTarget target(indexes[j_index].get());
      target.set_include_interiors(S2Testing::rnd.OneIn(2));
      TestFindClosestEdges(&target, &query);
    }
  }
}

static int const kNumIndexes = 50;
static int const kNumEdges = 100;
static int const kNumQueries = 200;

TEST(S2ClosestEdgeQuery, CircleEdges) {
  TestWithIndexFactory(RegularLoopShapeIndexFactory(),
                       kNumIndexes, kNumEdges, kNumQueries);
}

TEST(S2ClosestEdgeQuery, FractalEdges) {
  TestWithIndexFactory(FractalLoopShapeIndexFactory(),
                       kNumIndexes, kNumEdges, kNumQueries);
}

TEST(S2ClosestEdgeQuery, PointCloudEdges) {
  TestWithIndexFactory(PointCloudShapeIndexFactory(),
                       kNumIndexes, kNumEdges, kNumQueries);
}

