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
  query.FindClosestEdge(S2Point(1, 0, 0));
  EXPECT_EQ(0, query.num_edges());
  EXPECT_EQ(S1Angle::Infinity(), query.GetDistance(S2Point(1, 0, 0)));
  EXPECT_EQ(S2Point(1, 0, 0), query.Project(S2Point(1, 0, 0)));
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

// A (shape_id, edge_id) pair, used to identify a result edge.
struct EdgeId {
  EdgeId(int _shape_id, int _edge_id)
      : shape_id(_shape_id), edge_id(_edge_id) {
  }
  int shape_id, edge_id;
};
bool operator==(EdgeId const& x, EdgeId const& y) {
  return x.shape_id == y.shape_id && x.edge_id == y.edge_id;
}
bool operator<(EdgeId const& x, EdgeId const& y) {
  return (x.shape_id < y.shape_id ||
          (x.shape_id == y.shape_id && x.edge_id < y.edge_id));
}
ostream& operator<<(ostream& out, EdgeId const& x) {
  return out << '(' << x.shape_id << ", " << x.edge_id << ')';
}

using Result = pair<S1Angle, EdgeId>;

class PointTarget {
 public:
  explicit PointTarget(S2Point const& point) : point_(point) {}
  S1Angle GetDistanceToPoint(S2Point const& x) const {
    return S1Angle(x, point_);
  }
  S1Angle GetDistanceToEdge(S2Point const& v0, S2Point const& v1) const {
    return S2::GetDistance(point_, v0, v1);
  }
  void FindClosestEdges(S2ClosestEdgeQuery* query) const {
    query->FindClosestEdges(point_);
  }
 private:
  S2Point point_;
};

class EdgeTarget {
 public:
  EdgeTarget(S2Point const& a, S2Point const& b) : a_(a), b_(b) {}
  S1Angle GetDistanceToPoint(S2Point const& x) const {
    return S2::GetDistance(x, a_, b_);
  }
  S1Angle GetDistanceToEdge(S2Point const& v0, S2Point const& v1) const {
    S1ChordAngle distance = S1ChordAngle::Infinity();
    S2::UpdateEdgePairMinDistance(a_, b_, v0, v1, &distance);
    return distance.ToAngle();
  }
  void FindClosestEdges(S2ClosestEdgeQuery* query) const {
    query->FindClosestEdgesToEdge(a_, b_);
  }
 private:
  S2Point a_, b_;
};

// Use "query" to find the closest edge(s) to the given target, then convert
// the query results into two parallel vectors, one for distances and one for
// (shape_id, edge_id) pairs.  Also verify that the results satisfy the search
// criteria.
template <class Target>
static void GetClosestEdges(Target const& target, S2ClosestEdgeQuery *query,
                            vector<Result>* results) {
  target.FindClosestEdges(query);
  EXPECT_LE(query->num_edges(), query->max_edges());
  if (query->max_distance() == S1Angle::Infinity()) {
    // We can predict exactly how many edges should be returned.
    EXPECT_EQ(min(query->max_edges(), s2shapeutil::GetNumEdges(query->index())),
              query->num_edges());
  }
  for (int i = 0; i < query->num_edges(); ++i) {
    // Check that query->distance() is approximately equal to the
    // S1Angle(edge, target) distance.  They may be slightly different
    // because query->distance() is computed using S1ChordAngle.  Note that
    // the error gets considerably larger (1e-7) as the angle approaches Pi.
    auto edge = query->edge(i);
    EXPECT_NEAR(target.GetDistanceToEdge(edge.v0, edge.v1).radians(),
                query->distance(i).radians(), kChordAngleError);

    // Check that the edge satisfies the max_distance() condition.
    EXPECT_LE(query->distance(i), query->max_distance());
    results->push_back(make_pair(query->distance(i),
                                 EdgeId(query->shape_id(i),
                                        query->edge_id(i))));

    // Find the closest point on the edge and check its distance as well.
    S2Point closest = query->GetClosestPointOnEdge(i);
    EXPECT_NEAR(target.GetDistanceToPoint(closest).radians(),
                query->distance(i).radians(), kChordAngleError);
  }
}

template <class Target>
static void TestFindClosestEdges(Target const& target,
                                 S2ClosestEdgeQuery *query) {
  vector<Result> expected, actual;
  query->mutable_options()->set_use_brute_force(true);
  GetClosestEdges(target, query, &expected);
  query->mutable_options()->set_use_brute_force(false);
  GetClosestEdges(target, query, &actual);
  EXPECT_TRUE(CheckDistanceResults(expected, actual, query->max_edges(),
                                   query->max_distance(), query->max_error()))
      << "max_edges=" << query->max_edges()
      << ", max_distance=" << query->max_distance()
      << ", max_error=" << query->max_error();
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
      query.set_max_edges(1 + S2Testing::rnd.Uniform(100));
      if (S2Testing::rnd.OneIn(2)) {
        query.set_max_distance(S2Testing::rnd.RandDouble() * kRadius);
      }
      if (S2Testing::rnd.OneIn(2)) {
        // Choose a maximum error whose logarithm is uniformly distributed over
        // a reasonable range, except that it is sometimes zero.
        query.set_max_error(S1Angle::Radians(
            pow(1e-4, S2Testing::rnd.RandDouble()) * kRadius.radians()));
      }
      if (S2Testing::rnd.OneIn(2)) {
        // Find the closest edges to a given query edge.
        S2Point a = S2Testing::SamplePoint(query_cap);
        S2Point b = S2Testing::SamplePoint(
            S2Cap(a, pow(1e-4, S2Testing::rnd.RandDouble()) * kRadius));
        TestFindClosestEdges(EdgeTarget(a, b), &query);
      } else {
        S2Point target = S2Testing::SamplePoint(query_cap);
        // Find the closest edges to a given query point.
        TestFindClosestEdges(PointTarget(target), &query);
        // Also test the GetDistance() and Project() query methods.
        if (query.num_edges() == 0) {
          EXPECT_EQ(S1Angle::Infinity(), query.GetDistance(target));
          EXPECT_EQ(target, query.Project(target));
        } else {
          S1Angle expected_min_distance = query.distance(0);
          S1Angle actual_min_distance = query.GetDistance(target);
          EXPECT_LE(actual_min_distance,
                    expected_min_distance + query.max_error());
          EXPECT_NEAR(actual_min_distance.radians(),
                      S1Angle(target, query.Project(target)).radians(),
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

