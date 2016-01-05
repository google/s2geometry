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

#include "s2closestpointquery.h"

#include <memory>
#include <vector>

#include <gtest/gtest.h>
#include "s1angle.h"
#include "s2cap.h"
#include "s2cellid.h"
#include "s2loop.h"
#include "s2testing.h"

using std::make_pair;
using std::pair;
using std::unique_ptr;

//////////////////////////////////////////////////////////////////////////
// TODO(ericv): Move this code to s2testing.h after s2testing is no longer
// used in production code and this restriction is enforced (testonly=1).

#include <utility>
#include <set>
#include <algorithm>

namespace {

template <typename Id>
struct CompareFirst {
  inline bool operator()(std::pair<S1Angle, Id> const& x,
                         std::pair<S1Angle, Id> const& y) {
    return x.first < y.first;
  }
};

// Check that result set "x" contains all the expected results from "y", and
// does not include any duplicate results.
template <typename Id>
bool CheckResultSet(std::vector<std::pair<S1Angle, Id>> const& x,
                    std::vector<std::pair<S1Angle, Id>> const& y,
                    int max_size, S1Angle max_distance, S1Angle max_error,
                    S1Angle max_pruning_error, string const& label) {
  // Results should be sorted by distance.
  CompareFirst<Id> cmp;
  EXPECT_TRUE(std::is_sorted(x.begin(), x.end(), cmp));

  // Make sure there are no duplicate values.
  std::set<std::pair<S1Angle, Id>> x_set(x.begin(), x.end());
  EXPECT_EQ(x_set.size(), x.size()) << "Result set contains duplicates";

  // Result set X should contain all the items from U whose distance is less
  // than "limit" computed below.
  S1Angle limit = S1Angle::Zero();
  if (x.size() < max_size) {
    // Result set X was not limited by "max_size", so it should contain all
    // the items up to "max_distance", except that a few items right near the
    // distance limit may be missed because the distance measurements used for
    // pruning S2Cells are not conservative.
    limit = max_distance - max_pruning_error;
  } else if (x.size() > 0) {
    // Result set X contains only the closest "max_size" items, to within a
    // tolerance of "max_error + max_pruning_error".
    limit = x.back().first - max_error - max_pruning_error;
  }
  bool result = true;
  for (int i = 0; i < y.size(); ++i) {
    if (y[i].first < limit && std::count(x.begin(), x.end(), y[i]) != 1) {
      result = false;
      std::cout << label << " distance = " << y[i].first
                << ", id = " << y[i].second << "\n";
    }
  }
  return result;
}

// Compare two sets of "closest" items, where "expected" is computed via brute
// force (i.e., considering every possible candidate) and "actual" is computed
// using a spatial data structure.  Here "max_size" is a bound on the maximum
// number of items, "max_distance" is a limit on the distance to any item, and
// "max_error" is the maximum error allowed when selecting which items are
// closest (see S2ClosestEdgeQuery::max_error).
template <typename Id>
bool CheckDistanceResults(
    std::vector<std::pair<S1Angle, Id>> const& expected,
    std::vector<std::pair<S1Angle, Id>> const& actual,
    int max_size, S1Angle max_distance, S1Angle max_error) {
  static S1Angle const kMaxPruningError = S1Angle::Radians(1e-15);
  return (CheckResultSet(actual, expected, max_size, max_distance,
                         max_error, kMaxPruningError, "Missing") & /*not &&*/
          CheckResultSet(expected, actual, max_size, max_distance,
                         max_error, S1Angle::Zero(), "Extra"));
}

}  // namespace
//////////////////////////////////////////////////////////////////////////

typedef S2PointIndex<int> TestIndex;
typedef S2ClosestPointQuery<int> TestQuery;

TEST(S2ClosestPointQuery, NoPoints) {
  TestIndex index;
  TestQuery query(index);
  query.FindClosestPoint(S2Point(1, 0, 0));
  EXPECT_EQ(0, query.num_points());
}

TEST(S2ClosestPointQuery, ManyDuplicatePoints) {
  static int const kNumPoints = 10000;
  S2Point const kTestPoint(1, 0, 0);
  TestIndex index;
  for (int i = 0; i < kNumPoints; ++i) {
    index.Add(kTestPoint, i);
  }
  TestQuery query(index);
  query.FindClosestPoints(kTestPoint);
  EXPECT_EQ(kNumPoints, query.num_points());
}

// An abstract class that adds points to an S2PointIndex for benchmarking.
struct PointIndexFactory {
 public:
  virtual ~PointIndexFactory() {}

  // Given an index that will be queried using random points from "query_cap",
  // add approximately "num_vertices" points to "index".  (Typically the
  // indexed points will occupy some fraction of this cap.)
  virtual void AddPoints(S2Cap const& query_cap, int num_vertices,
                         TestIndex* index) const = 0;
};

// Generate points that are regularly spaced along a circle.  The circle is
// centered within the query cap and occupies 25% of its area, so that random
// query points have a 25% chance of being inside the circle.
//
// Points along a circle are nearly the worst case for distance calculations,
// since many points are nearly equidistant from any query point that is not
// immediately adjacent to the circle.
struct CirclePointIndexFactory : public PointIndexFactory {
  virtual void AddPoints(S2Cap const& query_cap, int num_vertices,
                         TestIndex* index) const {
    std::vector<S2Point> points = S2Testing::MakeRegularPoints(
        query_cap.center(), 0.5 * query_cap.GetRadius(), num_vertices);
    for (int i = 0; i < points.size(); ++i) {
      index->Add(points[i], i);
    }
  }
};

// Generate the vertices of a fractal whose convex hull approximately
// matches the query cap.
struct FractalPointIndexFactory : public PointIndexFactory {
  virtual void AddPoints(S2Cap const& query_cap, int num_vertices,
                         TestIndex* index) const {
    S2Testing::Fractal fractal;
    fractal.SetLevelForApproxMaxEdges(num_vertices);
    fractal.set_fractal_dimension(1.5);
    unique_ptr<S2Loop> loop(
        fractal.MakeLoop(S2Testing::GetRandomFrameAt(query_cap.center()),
                         query_cap.GetRadius()));
    for (int i = 0; i < loop->num_vertices(); ++i) {
      index->Add(loop->vertex(i), i);
    }
  }
};

// Generate vertices on a square grid that includes the entire query cap.
struct GridPointIndexFactory : public PointIndexFactory {
  virtual void AddPoints(S2Cap const& query_cap, int num_vertices,
                         TestIndex* index) const {
    int sqrt_num_points = ceil(sqrt(num_vertices));
    Matrix3x3_d frame = S2Testing::GetRandomFrameAt(query_cap.center());
    double radius = query_cap.GetRadius().radians();
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
static S1Angle const kRadius = S2Testing::KmToAngle(10);

// An approximate bound on the distance measurement error for "reasonable"
// distances (say, less than Pi/2) due to using S1ChordAngle.
static double kChordAngleError = 1e-15;

typedef pair<S1Angle, int> Result;

class PointTarget {
 public:
  explicit PointTarget(S2Point const& point) : point_(point) {}
  S1Angle GetDistance(S2Point const& x) const {
    return S1Angle(x, point_);
  }
  void FindClosestPoints(TestQuery* query) const {
    query->FindClosestPoints(point_);
  }
 private:
  S2Point point_;
};

class EdgeTarget {
 public:
  EdgeTarget(S2Point const& a, S2Point const& b) : a_(a), b_(b) {}
  S1Angle GetDistance(S2Point const& x) const {
    return S2EdgeUtil::GetDistance(x, a_, b_);
  }
  void FindClosestPoints(TestQuery* query) const {
    query->FindClosestPointsToEdge(a_, b_);
  }
 private:
  S2Point a_, b_;
};

// Use "query" to find the closest point(s) to the given target, and extract
// the query results into the given vector.  Also verify that the results
// satisfy the search criteria.
template <class Target>
static void GetClosestPoints(Target const& target, TestQuery* query,
                             std::vector<Result>* results) {
  target.FindClosestPoints(query);
  EXPECT_LE(query->num_points(), query->max_points());
  if (!query->region() && query->max_distance() == S1Angle::Infinity()) {
    // We can predict exactly how many points should be returned.
    EXPECT_EQ(std::min(query->max_points(), query->index().num_points()),
              query->num_points());
  }
  for (int i = 0; i < query->num_points(); ++i) {
    // Check that query->distance() is approximately equal to the
    // S1Angle(point, target) distance.  They may be slightly different
    // because query->distance() is computed using S1ChordAngle.  Note that
    // the error gets considerably larger (1e-7) as the angle approaches Pi.
    EXPECT_NEAR(target.GetDistance(query->point(i)).radians(),
                query->distance(i).radians(), kChordAngleError);
    // Check that the point satisfies the region() condition.
    if (query->region()) {
      EXPECT_TRUE(query->region()->VirtualContainsPoint(query->point(i)));
    }
    // Check that it satisfies the max_distance() condition.
    EXPECT_LE(query->distance(i), query->max_distance());
    results->push_back(make_pair(query->distance(i), query->data(i)));
  }
}

template <class Target>
static void TestFindClosestPoints(Target const& target, TestQuery *query) {
  std::vector<Result> expected, actual;
  query->UseBruteForce(true);
  GetClosestPoints(target, query, &expected);
  query->UseBruteForce(false);
  GetClosestPoints(target, query, &actual);
  EXPECT_TRUE(CheckDistanceResults(expected, actual, query->max_points(),
                                   query->max_distance(), S1Angle::Zero()))
      << "max_points=" << query->max_points()
      << ", max_distance=" << query->max_distance();
}

// The running time of this test is proportional to
//    num_indexes * num_vertices * num_queries.
static void TestWithIndexFactory(PointIndexFactory const& factory,
                                 int num_indexes, int num_vertices,
                                 int num_queries) {
  TestIndex index;
  for (int i_index = 0; i_index < num_indexes; ++i_index) {
    // Generate a point set and index it.
    S2Cap query_cap = S2Cap(S2Testing::RandomPoint(), kRadius);
    index.Reset();
    factory.AddPoints(query_cap, num_vertices, &index);
    for (int i_query = 0; i_query < num_queries; ++i_query) {
      // Use a new query each time to avoid resetting default parameters.
      TestQuery query(index);
      query.set_max_points(1 + S2Testing::rnd.Uniform(100));
      if (S2Testing::rnd.OneIn(2)) {
        query.set_max_distance(S2Testing::rnd.RandDouble() * kRadius);
      }
      S2LatLngRect rect = S2LatLngRect::FromCenterSize(
          S2LatLng(S2Testing::SamplePoint(query_cap)),
          S2LatLng(S2Testing::rnd.RandDouble() * kRadius,
                   S2Testing::rnd.RandDouble() * kRadius));
      if (S2Testing::rnd.OneIn(5)) {
        query.set_region(&rect);
      }
      if (S2Testing::rnd.OneIn(2)) {
        PointTarget target(S2Testing::SamplePoint(query_cap));
        TestFindClosestPoints(target, &query);
      } else {
        S2Point a = S2Testing::SamplePoint(query_cap);
        S2Point b = S2Testing::SamplePoint(
            S2Cap(a, pow(1e-4, S2Testing::rnd.RandDouble()) * kRadius));
        TestFindClosestPoints(EdgeTarget(a, b), &query);
      }
    }
  }
}

static int const kNumIndexes = 10;
static int const kNumVertices = 1000;
static int const kNumQueries = 50;

TEST(S2ClosestPointQueryTest, CirclePoints) {
  TestWithIndexFactory(CirclePointIndexFactory(),
                       kNumIndexes, kNumVertices, kNumQueries);
}

TEST(S2ClosestPointQueryTest, FractalPoints) {
  TestWithIndexFactory(FractalPointIndexFactory(),
                       kNumIndexes, kNumVertices, kNumQueries);
}

TEST(S2ClosestPointQueryTest, GridPoints) {
  TestWithIndexFactory(GridPointIndexFactory(),
                       kNumIndexes, kNumVertices, kNumQueries);
}

