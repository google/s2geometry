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
#include "s2/s2edge_distances.h"
#include "s2/s2edgeutil.h"
#include "s2/s2loop.h"
#include "s2/s2shapeutil.h"
#include "s2/s2testing.h"

using absl::MakeUnique;
using std::make_pair;
using std::min;
using std::ostream;
using std::pair;
using std::vector;

TEST(S2ClosestEdgeQuery, NoEdges) {
  S2ShapeIndex index;
  S2ClosestEdgeQuery query(&index);
  query.mutable_options()->set_max_edges(1);
  S2ClosestEdgeQuery::PointTarget const target(S2Point(1, 0, 0));
  auto const edge = query.FindClosestEdge(target);
  EXPECT_EQ(S1ChordAngle::Infinity(), edge.distance);
  EXPECT_EQ(-1, edge.edge_id);
  EXPECT_EQ(-1, edge.shape_id);
  EXPECT_EQ(S1ChordAngle::Infinity(), query.GetDistance(target));
}

// An abstract class that adds edges to an S2ShapeIndex for benchmarking.
class ShapeIndexFactory {
 public:
  virtual ~ShapeIndexFactory() {}

  // Given an index that will be queried using random points from "query_cap",
  // add approximately "num_edges" edges to "index".  (Typically the indexed
  // edges will occupy some fraction of this cap.)
  virtual void AddEdges(S2Cap const& query_cap, int num_edges,
                        S2ShapeIndex* index) const = 0;
};

// Generate a regular loop that occupies approximately 25% of the query cap,
// i.e. random query points have a 25% chance of being inside the loop.
//
// Regular loops are nearly the worst case for distance calculations, since
// many edges are nearly equidistant from any query point that is not
// immediately adjacent to the loop.
class RegularLoopShapeIndexFactory : public ShapeIndexFactory {
 public:
  void AddEdges(S2Cap const& query_cap, int num_edges,
                S2ShapeIndex* index) const override {
    index->Add(MakeUnique<S2Loop::OwningShape>(S2Loop::MakeRegularLoop(
        query_cap.center(), 0.5 * query_cap.GetRadius(), num_edges)));
  }
};

// Generate a fractal loop that occupies approximately 25% of the query cap,
// i.e. random query points have a 25% chance of being inside the loop.
class FractalLoopShapeIndexFactory : public ShapeIndexFactory {
 public:
  void AddEdges(S2Cap const& query_cap, int num_edges,
                S2ShapeIndex* index) const override {
    S2Testing::Fractal fractal;
    fractal.SetLevelForApproxMaxEdges(num_edges);
    index->Add(MakeUnique<S2Loop::OwningShape>(
        fractal.MakeLoop(S2Testing::GetRandomFrameAt(query_cap.center()),
                         0.5 * query_cap.GetRadius())));
  }
};

// The approximate radius of S2Cap from which query edges are chosen.
static S1Angle const kRadius = S2Testing::KmToAngle(100);

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
static void GetClosestEdges(S2ClosestEdgeQuery::Target const& target,
                            S2ClosestEdgeQuery *query,
                            vector<S2ClosestEdgeQuery::Result>* edges) {
  query->FindClosestEdges(target, edges);
  EXPECT_LE(edges->size(), query->options().max_edges());
  if (query->options().max_distance() ==
      S2ClosestEdgeQuery::Distance::Infinity()) {
    // We can predict exactly how many edges should be returned.
    EXPECT_EQ(min(query->options().max_edges(),
                  s2shapeutil::GetNumEdges(query->index())),
              edges->size());
  }
  for (auto const& edge : *edges) {
    // Check that the edge satisfies the max_distance() condition.
    EXPECT_LE(edge.distance, query->options().max_distance());
  }
}

static int TestFindClosestEdges(S2ClosestEdgeQuery::Target const& target,
                                S2ClosestEdgeQuery *query,
                                S2ClosestEdgeQuery::Result *closest_edge) {
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

  if (closest_edge != nullptr && !actual.empty())
    *closest_edge = actual.at(0);

  return actual.size();
}

// The running time of this test is proportional to
//    num_indexes * num_edges * num_queries.
static void TestWithIndexFactory(ShapeIndexFactory const& factory,
                                 int num_indexes, int num_edges,
                                 int num_queries) {
  S2ShapeIndex index;
  for (int i_index = 0; i_index < num_indexes; ++i_index) {
    S2Cap query_cap = S2Cap(S2Testing::RandomPoint(), kRadius);
    index.Reset();
    factory.AddEdges(query_cap, num_edges, &index);
    for (int i_query = 0; i_query < num_queries; ++i_query) {
      // Use a new query each time to avoid resetting default parameters.
      S2ClosestEdgeQuery query(&index);
      query.mutable_options()->set_max_edges(1 + S2Testing::rnd.Uniform(100));
      if (S2Testing::rnd.OneIn(2)) {
        query.mutable_options()->set_max_distance(
            S2Testing::rnd.RandDouble() * kRadius);
      }
      if (S2Testing::rnd.OneIn(2)) {
        // Choose a maximum error whose logarithm is uniformly distributed over
        // a reasonable range, except that it is sometimes zero.
        query.mutable_options()->set_max_error(S1Angle::Radians(
            pow(1e-4, S2Testing::rnd.RandDouble()) * kRadius.radians()));
      }
      if (S2Testing::rnd.OneIn(2)) {
        // Find the closest edges to a given query edge.
        S2Point a = S2Testing::SamplePoint(query_cap);
        S2Point b = S2Testing::SamplePoint(
            S2Cap(a, pow(1e-4, S2Testing::rnd.RandDouble()) * kRadius));
        S2ClosestEdgeQuery::EdgeTarget target(a, b);
        TestFindClosestEdges(target, &query, nullptr);
      } else {
        S2Point point = S2Testing::SamplePoint(query_cap);
        S2ClosestEdgeQuery::PointTarget target(point);
        // Find the closest edges to a given query point.
        S2ClosestEdgeQuery::Result closest_edge;
        int num_edges = TestFindClosestEdges(target, &query, &closest_edge);
        // Also test the GetDistance() and Project() query methods.
        if (num_edges == 0) {
          EXPECT_EQ(S1ChordAngle::Infinity(), query.GetDistance(target));
        } else {
          S1ChordAngle expected_min_distance = closest_edge.distance;
          S1ChordAngle actual_min_distance = query.GetDistance(target);
          EXPECT_LE(actual_min_distance,
                    expected_min_distance + query.options().max_error());
          EXPECT_NEAR(actual_min_distance.ToAngle().radians(),
                      S1Angle(point, query.Project(point)).radians(),
                      kChordAngleError);
        }
      }
    }
  }
}

static int const kNumIndexes = 5;
static int const kNumEdges = 500;
static int const kNumQueries = 1000;

TEST(S2ClosestEdgeQuery, CircleEdges) {
  TestWithIndexFactory(RegularLoopShapeIndexFactory(),
                       kNumIndexes, kNumEdges, kNumQueries);
}

TEST(S2ClosestEdgeQuery, FractalEdges) {
  TestWithIndexFactory(FractalLoopShapeIndexFactory(),
                       kNumIndexes, kNumEdges, kNumQueries);
}

