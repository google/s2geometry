// Copyright 2005 Google Inc. All Rights Reserved.
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

#include "s2edgeutil.h"

#include <float.h>
#include <math.h>
#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include <glog/logging.h>
#include "base/stringprintf.h"
#include <gtest/gtest.h>
#include "r1interval.h"
#include "r2rect.h"
#include "s1chordangle.h"
#include "s2polyline.h"
#include "s2testing.h"
#include "s2textformat.h"
#include "util/math/exactfloat/exactfloat.h"
#include "util/math/vector3.h"

using std::max;
using std::swap;
using std::unique_ptr;
using std::vector;

extern void PrintIntersectionStats();
extern int exact_calls;

const int kDegen = -2;
void CompareResult(int actual, int expected) {
  // HACK ALERT: RobustCrossing() is allowed to return 0 or -1 if either edge
  // is degenerate.  We use the value kDegen to represent this possibility.
  if (expected == kDegen) {
    EXPECT_LE(actual, 0);
  } else {
    EXPECT_EQ(expected, actual);
  }
}

void TestCrossing(S2Point a, S2Point b, S2Point c, S2Point d,
                  int robust, bool edge_or_vertex, bool simple) {
  a = a.Normalize();
  b = b.Normalize();
  c = c.Normalize();
  d = d.Normalize();
  CompareResult(S2EdgeUtil::RobustCrossing(a, b, c, d), robust);
  if (simple) {
    EXPECT_EQ(robust > 0, S2EdgeUtil::SimpleCrossing(a, b, c, d));
  }
  S2EdgeUtil::EdgeCrosser crosser(&a, &b, &c);
  CompareResult(crosser.RobustCrossing(&d), robust);
  CompareResult(crosser.RobustCrossing(&c), robust);

  EXPECT_EQ(edge_or_vertex, S2EdgeUtil::EdgeOrVertexCrossing(a, b, c, d));
  EXPECT_EQ(edge_or_vertex, crosser.EdgeOrVertexCrossing(&d));
  EXPECT_EQ(edge_or_vertex, crosser.EdgeOrVertexCrossing(&c));

  // Check that the crosser can be re-used.
  crosser.Init(&c, &d);
  crosser.RestartAt(&a);
  CompareResult(crosser.RobustCrossing(&b), robust);
  CompareResult(crosser.RobustCrossing(&a), robust);
}

void TestCrossings(S2Point a, S2Point b, S2Point c, S2Point d,
                   int robust, bool edge_or_vertex, bool simple) {
  TestCrossing(a, b, c, d, robust, edge_or_vertex, simple);
  TestCrossing(b, a, c, d, robust, edge_or_vertex, simple);
  TestCrossing(a, b, d, c, robust, edge_or_vertex, simple);
  TestCrossing(b, a, d, c, robust, edge_or_vertex, simple);
  TestCrossing(a, a, c, d, kDegen, 0, false);
  TestCrossing(a, b, c, c, kDegen, 0, false);
  TestCrossing(a, b, a, b, 0, 1, false);
  TestCrossing(c, d, a, b, robust, edge_or_vertex ^ (robust == 0), simple);
}

TEST(S2EdgeUtil, Crossings) {
  // The real tests of edge crossings are in s2{loop,polygon}_unittest,
  // but we do a few simple tests here.

  // Two regular edges that cross.
  TestCrossings(S2Point(1, 2, 1), S2Point(1, -3, 0.5),
                S2Point(1, -0.5, -3), S2Point(0.1, 0.5, 3), 1, true, true);

  // Two regular edges that cross antipodal points.
  TestCrossings(S2Point(1, 2, 1), S2Point(1, -3, 0.5),
                S2Point(-1, 0.5, 3), S2Point(-0.1, -0.5, -3), -1, false, true);

  // Two edges on the same great circle.
  TestCrossings(S2Point(0, 0, -1), S2Point(0, 1, 0),
                S2Point(0, 1, 1), S2Point(0, 0, 1), -1, false, true);

  // Two edges that cross where one vertex is S2::Origin().
  TestCrossings(S2Point(1, 0, 0), S2::Origin(),
                S2Point(1, -0.1, 1), S2Point(1, 1, -0.1), 1, true, true);

  // Two edges that cross antipodal points where one vertex is S2::Origin().
  TestCrossings(S2Point(1, 0, 0), S2Point(0, 1, 0),
                S2Point(0, 0, -1), S2Point(-1, -1, 1), -1, false, true);

  // Two edges that share an endpoint.  The Ortho() direction is (-4,0,2),
  // and edge CD is further CCW around (2,3,4) than AB.
  TestCrossings(S2Point(2, 3, 4), S2Point(-1, 2, 5),
                S2Point(7, -2, 3), S2Point(2, 3, 4), 0, false, true);

  // Two edges that barely cross each other near the middle of one edge.  The
  // edge AB is approximately in the x=y plane, while CD is approximately
  // perpendicular to it and ends exactly at the x=y plane.
  TestCrossings(S2Point(1, 1, 1), S2Point(1, nextafter(1, 0), -1),
                S2Point(11, -12, -1), S2Point(10, 10, 1), 1, true, false);

  // In this version, the edges are separated by a distance of about 1e-15.
  TestCrossings(S2Point(1, 1, 1), S2Point(1, nextafter(1, 2), -1),
                S2Point(1, -1, 0), S2Point(1, 1, 0), -1, false, false);

  // Two edges that barely cross each other near the end of both edges.  This
  // example cannot be handled using regular double-precision arithmetic due
  // to floating-point underflow.
  TestCrossings(S2Point(0, 0, 1), S2Point(2, -1e-323, 1),
                S2Point(1, -1, 1), S2Point(1e-323, 0, 1), 1, true, false);

  // In this version, the edges are separated by a distance of about 1e-640.
  TestCrossings(S2Point(0, 0, 1), S2Point(2, 1e-323, 1),
                S2Point(1, -1, 1), S2Point(1e-323, 0, 1), -1, false, false);

  // Two edges that barely cross each other near the middle of one edge.
  // Computing the exact determinant of some of the triangles in this test
  // requires more than 2000 bits of precision.
  TestCrossings(S2Point(1, -1e-323, -1e-323), S2Point(1e-323, 1, 1e-323),
                S2Point(1, -1, 1e-323), S2Point(1, 1, 0),
                1, true, false);

  // In this version, the edges are separated by a distance of about 1e-640.
  TestCrossings(S2Point(1, 1e-323, -1e-323), S2Point(-1e-323, 1, 1e-323),
                S2Point(1, -1, 1e-323), S2Point(1, 1, 0),
                -1, false, false);
}


S2LatLngRect GetEdgeBound(S2Point const& a, S2Point const& b) {
  S2EdgeUtil::RectBounder bounder;
  bounder.AddPoint(&a);
  bounder.AddPoint(&b);
  return bounder.GetBound();
}

S2LatLngRect GetEdgeBound(double x1, double y1, double z1,
                          double x2, double y2, double z2) {
  return GetEdgeBound(S2Point(x1, y1, z1).Normalize(),
                      S2Point(x2, y2, z2).Normalize());
}

S2LatLng const kRectError = S2EdgeUtil::RectBounder::MaxErrorForTests();

TEST(RectBounder, MaxLatitudeSimple) {
  // Check cases where the min/max latitude is attained at a vertex.
  static double const kCubeLat = asin(1 / sqrt(3));  // 35.26 degrees
  EXPECT_TRUE(GetEdgeBound(1,1,1, 1,-1,-1).ApproxEquals(  // NOLINT
      S2LatLngRect(R1Interval(-kCubeLat, kCubeLat),
                   S1Interval(-M_PI_4, M_PI_4)), kRectError));
  EXPECT_TRUE(GetEdgeBound(1,-1,1, 1,1,-1).ApproxEquals(  // NOLINT
      S2LatLngRect(R1Interval(-kCubeLat, kCubeLat),
                   S1Interval(-M_PI_4, M_PI_4)), kRectError));

  // Check cases where the min/max latitude occurs in the edge interior.
  // These tests expect the result to be pretty close to the middle of the
  // allowable error range (i.e., by adding 0.5 * kRectError).

  // Max latitude, CW edge
  EXPECT_DOUBLE_EQ(M_PI_4 + 0.5 * kRectError.lat().radians(),
                   GetEdgeBound(1,1,1, 1,-1,1).lat().hi());
  // Max latitude, CCW edge
  EXPECT_DOUBLE_EQ(M_PI_4 + 0.5 * kRectError.lat().radians(),
                   GetEdgeBound(1,-1,1, 1,1,1).lat().hi());  // NOLINT
  // Min latitude, CW edge
  EXPECT_DOUBLE_EQ(-M_PI_4 - 0.5 * kRectError.lat().radians(),
                   GetEdgeBound(1,-1,-1, -1,-1,-1).lat().lo());  // NOLINT
  // Min latitude, CCW edge
  EXPECT_DOUBLE_EQ(-M_PI_4 - 0.5 * kRectError.lat().radians(),
                   GetEdgeBound(-1,1,-1, -1,-1,-1).lat().lo());  // NOLINT

  // Check cases where the edge passes through one of the poles.
  EXPECT_EQ(M_PI_2, GetEdgeBound(.3,.4,1, -.3,-.4,1).lat().hi());  // NOLINT
  EXPECT_EQ(-M_PI_2, GetEdgeBound(.3,.4,-1, -.3,-.4,-1).lat().lo());  // NOLINT
}

TEST(RectBounder, MaxLatitudeRandom) {
  // Check that the maximum latitude of edges is computed accurately to within
  // 3 * DBL_EPSILON (the expected maximum error).  We concentrate on maximum
  // latitudes near the equator and north pole since these are the extremes.

  S2Testing::Random* rnd = &S2Testing::rnd;
  const int kIters = 100;
  for (int iter = 0; iter < kIters; ++iter) {
    // Construct a right-handed coordinate frame (U,V,W) such that U points
    // slightly above the equator, V points at the equator, and W is slightly
    // offset from the north pole.
    S2Point u = S2Testing::RandomPoint();
    u[2] = DBL_EPSILON * 1e-6 * pow(1e12, rnd->RandDouble());  // log is uniform
    u = u.Normalize();
    S2Point v = S2::RobustCrossProd(S2Point(0, 0, 1), u).Normalize();
    S2Point w = S2::RobustCrossProd(u, v).Normalize();

    // Construct a line segment AB that passes through U, and check that the
    // maximum latitude of this segment matches the latitude of U.
    S2Point a = (u - rnd->RandDouble() * v).Normalize();
    S2Point b = (u + rnd->RandDouble() * v).Normalize();
    S2LatLngRect ab_bound = GetEdgeBound(a, b);
    EXPECT_NEAR(S2LatLng::Latitude(u).radians(),
                ab_bound.lat().hi(), kRectError.lat().radians());

    // Construct a line segment CD that passes through W, and check that the
    // maximum latitude of this segment matches the latitude of W.
    S2Point c = (w - rnd->RandDouble() * v).Normalize();
    S2Point d = (w + rnd->RandDouble() * v).Normalize();
    S2LatLngRect cd_bound = GetEdgeBound(c, d);
    EXPECT_NEAR(S2LatLng::Latitude(w).radians(),
                cd_bound.lat().hi(), kRectError.lat().radians());
  }
}

S2Point PerturbATowardsB(S2Point const& a, S2Point const& b) {
  S2Testing::Random* rnd = &S2Testing::rnd;
  double choice = rnd->RandDouble();
  if (choice < 0.1) {
    return a;
  }
  if (choice < 0.3) {
    // Return a point that is exactly proportional to A and that still
    // satisfies S2::IsUnitLength().
    for (;;) {
      S2Point b = (2 - a.Norm() + 5*(rnd->RandDouble()-0.5) * DBL_EPSILON) * a;
      if (b != a && S2::IsUnitLength(b))
        return b;
    }
  }
  if (choice < 0.5) {
    // Return a point such that the distance squared to A will underflow.
    return S2EdgeUtil::InterpolateAtDistance(S1Angle::Radians(1e-300), a, b);
  }
  // Otherwise return a point whose distance from A is near DBL_EPSILON such
  // that the log of the pdf is uniformly distributed.
  double distance = DBL_EPSILON * 1e-5 * pow(1e6, rnd->RandDouble());
  return S2EdgeUtil::InterpolateAtDistance(S1Angle::Radians(distance), a, b);
}

S2Point RandomPole() {
  return S2Point(0, 0, S2Testing::rnd.OneIn(2) ? 1 : -1);
}

S2Point PointNearPole() {
  return PerturbATowardsB(RandomPole(), S2Testing::RandomPoint());
}

S2Point PointNearEquator() {
  return PerturbATowardsB(S2Point(S2Testing::rnd.RandDouble(),
                                  S2Testing::rnd.RandDouble(), 0).Normalize(),
                          RandomPole());
}

TEST(RectBounder, NearlyIdenticalOrAntipodalPoints) {
  // Test pairs of points that are either:
  //  - identical
  //  - nearly or exactly proportional, e.g. (1,0,0) vs. (1+2e-16, 0, 0)
  //  - very close to each other
  // Furthermore we want to test cases where the two points are:
  //  - on a nearly-polar great circle
  //  - on a nearly-equatorial great circle
  //  - near the poles, but on any great circle
  //  - near the equator, but on any great circle
  //  - positioned arbitrarily
  // Also test the corresponding situations for antipodal points, i.e. by
  // negating one of the points so that they are almost 180 degrees apart.

  S2Testing::Random* rnd = &S2Testing::rnd;
  const int kIters = 10000;
  for (int iter = 0; iter < kIters; ++iter) {
    SCOPED_TRACE(StringPrintf("Iteration %d", iter));
    S2Point a, b;
    switch (rnd->Uniform(5)) {
      case 0:
        // Two nearby points on a nearly-polar great circle.
        a = S2Testing::RandomPoint();
        b = PerturbATowardsB(a, PointNearPole());
        break;
      case 1:
        // Two nearby points on a nearly-equatorial great circle.
        a = PointNearEquator();
        b = PerturbATowardsB(a, PointNearEquator());
        break;
      case 2:
        // Two nearby points near a pole, but on any great circle.
        a = PointNearPole();
        b = PerturbATowardsB(a, S2Testing::RandomPoint());
        break;
      case 3:
        // Two nearby points near the equator, but on any great circle.
        a = PointNearEquator();
        b = PerturbATowardsB(a, S2Testing::RandomPoint());
        break;
      case 4:
        // Two nearby points anywhere on the sphere.
        a = S2Testing::RandomPoint();
        b = PerturbATowardsB(a, S2Testing::RandomPoint());
        break;
    }
    // The two points are chosen to be so close to each other that the min/max
    // latitudes are nearly always achieved at the edge endpoints.  The only
    // thing we need to watch out for is that the latitude error bound is
    // slightly larger if the min/max latitude occurs in the edge interior.
    S2LatLngRect expected_bound = S2LatLngRect::FromPointPair(S2LatLng(a),
                                                              S2LatLng(b));
    S2LatLngRect bound = GetEdgeBound(a, b);
    EXPECT_TRUE(bound.Contains(expected_bound));
    EXPECT_TRUE(expected_bound.Expanded(kRectError).PolarClosure().
                Contains(bound));

    // If the two points are close enough and one point is negated (antipodal
    // points), the bound should be the entire sphere.
    if ((a - b).CrossProd(a + b).Norm() <= 6.110 * DBL_EPSILON) {
      EXPECT_EQ(S2LatLngRect::Full(), GetEdgeBound(a, -b));
    }
  }
}

S2LatLngRect GetSubregionBound(double x_lat, double x_lng,
                               double y_lat, double y_lng) {
  S2LatLngRect in = S2LatLngRect::FromPointPair(
      S2LatLng::FromRadians(x_lat, x_lng),
      S2LatLng::FromRadians(y_lat, y_lng));
  S2LatLngRect out = S2EdgeUtil::RectBounder::ExpandForSubregions(in);

  // Test that the bound is actually expanded.
  EXPECT_TRUE(out.Contains(in));
  if (in.lat() == S2LatLngRect::FullLat()) {
    EXPECT_FALSE(in.lat().Contains(out.lat()));
  }
  return out;
}


TEST(RectBounder, ExpandForSubregions) {
  // First we check the various situations where the bound contains
  // nearly-antipodal points.  The tests are organized into pairs where the
  // two bounds are similar except that the first bound meets the
  // nearly-antipodal criteria while the second does not.

  // Cases where the bound does not straddle the equator (but almost does),
  // and spans nearly 180 degrees in longitude.
  EXPECT_TRUE(GetSubregionBound(3e-16, 0, 1e-14, M_PI).is_full());
  EXPECT_FALSE(GetSubregionBound(9e-16, 0, 1e-14, M_PI).is_full());
  EXPECT_TRUE(GetSubregionBound(1e-16, 7e-16, 1e-14, M_PI).is_full());
  EXPECT_FALSE(GetSubregionBound(3e-16, 14e-16, 1e-14, M_PI).is_full());
  EXPECT_TRUE(GetSubregionBound(1e-100, 14e-16, 1e-14, M_PI).is_full());
  EXPECT_FALSE(GetSubregionBound(1e-100, 22e-16, 1e-14, M_PI).is_full());

  // Cases where the bound spans at most 90 degrees in longitude, and almost
  // 180 degrees in latitude.  Note that DBL_EPSILON is about 2.22e-16, which
  // implies that the double-precision value just below Pi/2 can be written as
  // (M_PI_2 - 2e-16).
  EXPECT_TRUE(GetSubregionBound(-M_PI_2, -1e-15, M_PI_2 - 7e-16, 0).
              is_full());
  EXPECT_FALSE(GetSubregionBound(-M_PI_2, -1e-15, M_PI_2 - 30e-16, 0).
              is_full());
  EXPECT_TRUE(GetSubregionBound(-M_PI_2 + 4e-16, 0, M_PI_2 - 2e-16, 1e-7).
              is_full());
  EXPECT_FALSE(GetSubregionBound(-M_PI_2 + 30e-16, 0, M_PI_2, 1e-7).
              is_full());
  EXPECT_TRUE(GetSubregionBound(-M_PI_2 + 4e-16, 0, M_PI_2 - 4e-16, M_PI_2).
              is_full());
  EXPECT_FALSE(GetSubregionBound(-M_PI_2, 0, M_PI_2 - 30e-16, M_PI_2).
              is_full());

  // Cases where the bound straddles the equator and spans more than 90
  // degrees in longitude.  These are the cases where the critical distance is
  // between a corner of the bound and the opposite longitudinal edge.  Unlike
  // the cases above, here the bound may contain nearly-antipodal points (to
  // within 3.055 * DBL_EPSILON) even though the latitude and longitude ranges
  // are both significantly less than (Pi - 3.055 * DBL_EPSILON).
  EXPECT_TRUE(GetSubregionBound(-M_PI_2, 0, M_PI_2 - 1e-8, M_PI - 1e-7).
              is_full());
  EXPECT_FALSE(GetSubregionBound(-M_PI_2, 0, M_PI_2 - 1e-7, M_PI - 1e-7).
               is_full());
  EXPECT_TRUE(GetSubregionBound(-M_PI_2 + 1e-12, -M_PI + 1e-4, M_PI_2, 0).
              is_full());
  EXPECT_TRUE(GetSubregionBound(-M_PI_2 + 1e-11, -M_PI + 1e-4, M_PI_2, 0).
               is_full());

  // Now we test cases where the bound does not contain nearly-antipodal
  // points, but it does contain points that are approximately 180 degrees
  // apart in latitude.
  EXPECT_TRUE(GetSubregionBound(1.5, -M_PI_2, 1.5, M_PI_2 - 2e-16).
              ApproxEquals(S2LatLngRect(R1Interval(1.5, 1.5),
                                        S1Interval::Full()), kRectError));
  EXPECT_TRUE(GetSubregionBound(1.5, -M_PI_2, 1.5, M_PI_2 - 7e-16).
              ApproxEquals(S2LatLngRect(R1Interval(1.5, 1.5),
                                        S1Interval(-M_PI_2, M_PI_2 - 7e-16)),
                           kRectError));

  // Test the full and empty bounds.
  EXPECT_TRUE(S2EdgeUtil::RectBounder::ExpandForSubregions(
      S2LatLngRect::Full()).is_full());
  EXPECT_TRUE(S2EdgeUtil::RectBounder::ExpandForSubregions(
      S2LatLngRect::Empty()).is_empty());

  // Check for cases where the bound is expanded to include one of the poles.
  EXPECT_TRUE(GetSubregionBound(-M_PI_2 + 1e-15, 0, -M_PI_2 + 1e-15, 0).
              ApproxEquals(S2LatLngRect(R1Interval(-M_PI_2, -M_PI_2 + 1e-15),
                                        S1Interval::Full()), kRectError));
  EXPECT_TRUE(GetSubregionBound(M_PI_2 - 1e-15, 0, M_PI_2 - 1e-15, 0).
              ApproxEquals(S2LatLngRect(R1Interval(M_PI_2 - 1e-15, M_PI_2),
                                        S1Interval::Full()), kRectError));
}

TEST(S2EdgeUtil, LongitudePruner) {
  S2EdgeUtil::LongitudePruner pruner1(S1Interval(0.75 * M_PI, -0.75 * M_PI),
                                      S2Point(0, 1, 2));
  EXPECT_FALSE(pruner1.Intersects(S2Point(1, 1, 3)));
  EXPECT_TRUE(pruner1.Intersects(S2Point(-1 - 1e-15, -1, 0)));
  EXPECT_TRUE(pruner1.Intersects(S2Point(-1, 0, 0)));
  EXPECT_TRUE(pruner1.Intersects(S2Point(-1, 0, 0)));
  EXPECT_TRUE(pruner1.Intersects(S2Point(1, -1, 8)));
  EXPECT_FALSE(pruner1.Intersects(S2Point(1, 0, -2)));
  EXPECT_TRUE(pruner1.Intersects(S2Point(-1, -1e-15, 0)));

  S2EdgeUtil::LongitudePruner pruner2(S1Interval(0.25 * M_PI, 0.25 * M_PI),
                                      S2Point(1, 0, 0));
  EXPECT_FALSE(pruner2.Intersects(S2Point(2, 1, 2)));
  EXPECT_TRUE(pruner2.Intersects(S2Point(1, 2, 3)));
  EXPECT_FALSE(pruner2.Intersects(S2Point(0, 1, 4)));
  EXPECT_FALSE(pruner2.Intersects(S2Point(-1e-15, -1, -1)));
}

void TestWedge(S2Point a0, S2Point ab1, S2Point a2, S2Point b0, S2Point b2,
               bool contains, bool intersects,
               S2EdgeUtil::WedgeRelation wedge_relation) {
  a0 = a0.Normalize();
  ab1 = ab1.Normalize();
  a2 = a2.Normalize();
  b0 = b0.Normalize();
  b2 = b2.Normalize();
  EXPECT_EQ(contains, S2EdgeUtil::WedgeContains(a0, ab1, a2, b0, b2));
  EXPECT_EQ(intersects, S2EdgeUtil::WedgeIntersects(a0, ab1, a2, b0, b2));
  EXPECT_EQ(wedge_relation, S2EdgeUtil::GetWedgeRelation(a0, ab1, a2, b0, b2));
}

TEST(S2EdgeUtil, Wedges) {
  // For simplicity, all of these tests use an origin of (0, 0, 1).
  // This shouldn't matter as long as the lower-level primitives are
  // implemented correctly.

  // Intersection in one wedge.
  TestWedge(S2Point(-1, 0, 10), S2Point(0, 0, 1), S2Point(1, 2, 10),
            S2Point(0, 1, 10), S2Point(1, -2, 10),
            false, true, S2EdgeUtil::WEDGE_PROPERLY_OVERLAPS);
  // Intersection in two wedges.
  TestWedge(S2Point(-1, -1, 10), S2Point(0, 0, 1), S2Point(1, -1, 10),
            S2Point(1, 0, 10), S2Point(-1, 1, 10),
            false, true, S2EdgeUtil::WEDGE_PROPERLY_OVERLAPS);

  // Normal containment.
  TestWedge(S2Point(-1, -1, 10), S2Point(0, 0, 1), S2Point(1, -1, 10),
            S2Point(-1, 0, 10), S2Point(1, 0, 10),
            true, true, S2EdgeUtil::WEDGE_PROPERLY_CONTAINS);
  // Containment with equality on one side.
  TestWedge(S2Point(2, 1, 10), S2Point(0, 0, 1), S2Point(-1, -1, 10),
            S2Point(2, 1, 10), S2Point(1, -5, 10),
            true, true, S2EdgeUtil::WEDGE_PROPERLY_CONTAINS);
  // Containment with equality on the other side.
  TestWedge(S2Point(2, 1, 10), S2Point(0, 0, 1), S2Point(-1, -1, 10),
            S2Point(1, -2, 10), S2Point(-1, -1, 10),
            true, true, S2EdgeUtil::WEDGE_PROPERLY_CONTAINS);

  // Containment with equality on both sides.
  TestWedge(S2Point(-2, 3, 10), S2Point(0, 0, 1), S2Point(4, -5, 10),
            S2Point(-2, 3, 10), S2Point(4, -5, 10),
            true, true, S2EdgeUtil::WEDGE_EQUALS);

  // Disjoint with equality on one side.
  TestWedge(S2Point(-2, 3, 10), S2Point(0, 0, 1), S2Point(4, -5, 10),
            S2Point(4, -5, 10), S2Point(-2, -3, 10),
            false, false, S2EdgeUtil::WEDGE_IS_DISJOINT);
  // Disjoint with equality on the other side.
  TestWedge(S2Point(-2, 3, 10), S2Point(0, 0, 1), S2Point(0, 5, 10),
            S2Point(4, -5, 10), S2Point(-2, 3, 10),
            false, false, S2EdgeUtil::WEDGE_IS_DISJOINT);
  // Disjoint with equality on both sides.
  TestWedge(S2Point(-2, 3, 10), S2Point(0, 0, 1), S2Point(4, -5, 10),
            S2Point(4, -5, 10), S2Point(-2, 3, 10),
            false, false, S2EdgeUtil::WEDGE_IS_DISJOINT);

  // B contains A with equality on one side.
  TestWedge(S2Point(2, 1, 10), S2Point(0, 0, 1), S2Point(1, -5, 10),
            S2Point(2, 1, 10), S2Point(-1, -1, 10),
            false, true, S2EdgeUtil::WEDGE_IS_PROPERLY_CONTAINED);
  // B contains A with equality on the other side.
  TestWedge(S2Point(2, 1, 10), S2Point(0, 0, 1), S2Point(1, -5, 10),
            S2Point(-2, 1, 10), S2Point(1, -5, 10),
            false, true, S2EdgeUtil::WEDGE_IS_PROPERLY_CONTAINED);
}

// Given a point X and an edge AB, check that the distance from X to AB is
// "distance_radians" and the closest point on AB is "expected_closest".
void CheckDistance(S2Point x, S2Point a, S2Point b,
                   double distance_radians, S2Point expected_closest) {
  x = x.Normalize();
  a = a.Normalize();
  b = b.Normalize();
  expected_closest = expected_closest.Normalize();
  EXPECT_NEAR(distance_radians, S2EdgeUtil::GetDistance(x, a, b).radians(),
              1e-15);
  S2Point closest = S2EdgeUtil::GetClosestPoint(x, a, b);
  if (expected_closest == S2Point(0, 0, 0)) {
    // This special value says that the result should be A or B.
    EXPECT_TRUE(closest == a || closest == b);
  } else {
    EXPECT_TRUE(S2::ApproxEquals(closest, expected_closest));
  }
  S1ChordAngle min_distance = S1ChordAngle::Zero();
  EXPECT_FALSE(S2EdgeUtil::UpdateMinDistance(x, a, b, &min_distance));
  min_distance = S1ChordAngle::Infinity();
  EXPECT_TRUE(S2EdgeUtil::UpdateMinDistance(x, a, b, &min_distance));
  EXPECT_NEAR(distance_radians, min_distance.ToAngle().radians(), 1e-15);
}

TEST(S2EdgeUtil, Distance) {
  CheckDistance(S2Point(1, 0, 0), S2Point(1, 0, 0), S2Point(0, 1, 0),
                0, S2Point(1, 0, 0));
  CheckDistance(S2Point(0, 1, 0), S2Point(1, 0, 0), S2Point(0, 1, 0),
                0, S2Point(0, 1, 0));
  CheckDistance(S2Point(1, 3, 0), S2Point(1, 0, 0), S2Point(0, 1, 0),
                0, S2Point(1, 3, 0));
  CheckDistance(S2Point(0, 0, 1), S2Point(1, 0, 0), S2Point(0, 1, 0),
                M_PI_2, S2Point(1, 0, 0));
  CheckDistance(S2Point(0, 0, -1), S2Point(1, 0, 0), S2Point(0, 1, 0),
                M_PI_2, S2Point(1, 0, 0));
  CheckDistance(S2Point(-1, -1, 0), S2Point(1, 0, 0), S2Point(0, 1, 0),
                0.75 * M_PI, S2Point(0, 0, 0));

  CheckDistance(S2Point(0, 1, 0), S2Point(1, 0, 0), S2Point(1, 1, 0),
                M_PI_4, S2Point(1, 1, 0));
  CheckDistance(S2Point(0, -1, 0), S2Point(1, 0, 0), S2Point(1, 1, 0),
                M_PI_2, S2Point(1, 0, 0));

  CheckDistance(S2Point(0, -1, 0), S2Point(1, 0, 0), S2Point(-1, 1, 0),
                M_PI_2, S2Point(1, 0, 0));
  CheckDistance(S2Point(-1, -1, 0), S2Point(1, 0, 0), S2Point(-1, 1, 0),
                M_PI_2, S2Point(-1, 1, 0));

  CheckDistance(S2Point(1, 1, 1), S2Point(1, 0, 0), S2Point(0, 1, 0),
                asin(sqrt(1./3)), S2Point(1, 1, 0));
  CheckDistance(S2Point(1, 1, -1), S2Point(1, 0, 0), S2Point(0, 1, 0),
                asin(sqrt(1./3)), S2Point(1, 1, 0));

  CheckDistance(S2Point(-1, 0, 0), S2Point(1, 1, 0), S2Point(1, 1, 0),
                0.75 * M_PI, S2Point(1, 1, 0));
  CheckDistance(S2Point(0, 0, -1), S2Point(1, 1, 0), S2Point(1, 1, 0),
                M_PI_2, S2Point(1, 1, 0));
  CheckDistance(S2Point(-1, 0, 0), S2Point(1, 0, 0), S2Point(1, 0, 0),
                M_PI, S2Point(1, 0, 0));
}

void CheckInterpolate(double t, S2Point a, S2Point b, S2Point expected) {
  a = a.Normalize();
  b = b.Normalize();
  expected = expected.Normalize();
  S2Point actual = S2EdgeUtil::Interpolate(t, a, b);

  // We allow a bit more than the usual 1e-15 error tolerance because
  // Interpolate() uses trig functions.
  EXPECT_TRUE(S2::ApproxEquals(expected, actual, 3e-15))
      << "Expected: " << expected << ", actual: " << actual;
}

TEST(S2EdgeUtil, Interpolate) {
  // Choose test points designed to expose floating-point errors.
  S2Point p1 = S2Point(0.1, 1e-30, 0.3).Normalize();
  S2Point p2 = S2Point(-0.7, -0.55, -1e30).Normalize();

  // A zero-length edge.
  CheckInterpolate(0, p1, p1, p1);
  CheckInterpolate(1, p1, p1, p1);

  // Start, end, and middle of a medium-length edge.
  CheckInterpolate(0, p1, p2, p1);
  CheckInterpolate(1, p1, p2, p2);
  CheckInterpolate(0.5, p1, p2, 0.5 * (p1 + p2));

  // Test that interpolation is done using distances on the sphere rather than
  // linear distances.
  CheckInterpolate(1./3, S2Point(1, 0, 0), S2Point(0, 1, 0),
                   S2Point(sqrt(3), 1, 0));
  CheckInterpolate(2./3, S2Point(1, 0, 0), S2Point(0, 1, 0),
                   S2Point(1, sqrt(3), 0));

  // Test that interpolation is accurate on a long edge (but not so long that
  // the definition of the edge itself becomes too unstable).
  {
    double const kLng = M_PI - 1e-2;
    S2Point a = S2LatLng::FromRadians(0, 0).ToPoint();
    S2Point b = S2LatLng::FromRadians(0, kLng).ToPoint();
    for (double f = 0.4; f > 1e-15; f *= 0.1) {
      CheckInterpolate(f, a, b,
                       S2LatLng::FromRadians(0, f * kLng).ToPoint());
      CheckInterpolate(1 - f, a, b,
                       S2LatLng::FromRadians(0, (1 - f) * kLng).ToPoint());
    }
  }

  // Test that interpolation on a 180 degree edge (antipodal endpoints) yields
  // a result with the correct distance from each endpoint.
  for (double t = 0; t <= 1; t += 0.125) {
    S2Point actual = S2EdgeUtil::Interpolate(t, p1, -p1);
    EXPECT_NEAR(S1Angle(actual, p1).radians(), t * M_PI, 3e-15);
  }
}

TEST(S2EdgeUtil, InterpolateCanExtrapolate) {
  const S2Point i(1, 0, 0);
  const S2Point j(0, 1, 0);
  // Initial vectors at 90 degrees.
  CheckInterpolate(0, i, j, S2Point(1, 0, 0));
  CheckInterpolate(1, i, j, S2Point(0, 1, 0));
  CheckInterpolate(1.5, i, j, S2Point(-1, 1, 0));
  CheckInterpolate(2, i, j, S2Point(-1, 0, 0));
  CheckInterpolate(3, i, j, S2Point(0, -1, 0));
  CheckInterpolate(4, i, j, S2Point(1, 0, 0));

  // Negative values of t.
  CheckInterpolate(-1, i, j, S2Point(0, -1, 0));
  CheckInterpolate(-2, i, j, S2Point(-1, 0, 0));
  CheckInterpolate(-3, i, j, S2Point(0, 1, 0));
  CheckInterpolate(-4, i, j, S2Point(1, 0, 0));

  // Initial vectors at 45 degrees.
  CheckInterpolate(2, i, S2Point(1, 1, 0), S2Point(0, 1, 0));
  CheckInterpolate(3, i, S2Point(1, 1, 0), S2Point(-1, 1, 0));
  CheckInterpolate(4, i, S2Point(1, 1, 0), S2Point(-1, 0, 0));

  // Initial vectors at 135 degrees.
  CheckInterpolate(2, i, S2Point(-1, 1, 0), S2Point(0, -1, 0));

  // Take a small fraction along the curve.
  S2Point p(S2EdgeUtil::Interpolate(0.001, i, j));
  // We should get back where we started.
  CheckInterpolate(1000, i, p, j);
}


TEST(S2EdgeUtil, RepeatedInterpolation) {
  // Check that points do not drift away from unit length when repeated
  // interpolations are done.
  for (int i = 0; i < 100; ++i) {
    S2Point a = S2Testing::RandomPoint();
    S2Point b = S2Testing::RandomPoint();
    for (int j = 0; j < 1000; ++j) {
      a = S2EdgeUtil::Interpolate(0.01, a, b);
    }
    EXPECT_TRUE(S2::IsUnitLength(a));
  }
}

class S2EdgeUtilTesting {
 public:
  static S2Point GetIntersection(S2Point const& a, S2Point const& b,
                                 S2Point const& c, S2Point const& d) {
    S2Point result = S2EdgeUtil::GetIntersection(a, b, c, d);
    ++tally_[S2EdgeUtil::last_intersection_method_];
    return result;
  }

  static void PrintIntersectionStats() {
    int total = 0;
    int totals[kNumMethods];
    for (int i = kNumMethods; --i >= 0; ) {
      total += tally_[i];
      totals[i] = total;
    }
    printf("%10s %16s %16s  %6s\n",
           "Method", "Successes", "Attempts", "Rate");
    for (int i = 0; i < kNumMethods; ++i) {
      if (tally_[i] == 0) continue;
      printf("%10s %9d %5.1f%% %9d %5.1f%%  %5.1f%%\n",
             S2EdgeUtil::GetIntersectionMethodName(
                 static_cast<S2EdgeUtil::IntersectionMethod>(i)),
             tally_[i], 100.0 * tally_[i] / total,
             totals[i], 100.0 * totals[i] / total,
             100.0 * tally_[i] / totals[i]);
    }
    for (int i = 0; i < kNumMethods; ++i) tally_[i] = 0;
  }

  // This returns the true intersection point of two line segments (a0,a1) and
  // (b0,b1), with a relative error of at most DBL_EPSILON in each coordinate
  // (i.e., one ulp, or twice the double precision rounding error).
  static S2Point GetIntersectionExact(S2Point const& a0, S2Point const& a1,
                                      S2Point const& b0, S2Point const& b1) {
    S2Point x = S2EdgeUtil::GetIntersectionExact(a0, a1, b0, b1);
    if (x.DotProd((a0 + a1) + (b0 + b1)) < 0) x = -x;
    return x;
  }

 private:
  static int const kNumMethods = S2EdgeUtil::NUM_INTERSECTION_METHODS;
  static int tally_[kNumMethods];
};
int S2EdgeUtilTesting::tally_[kNumMethods];

// The approximate maximum error in GetDistance() for small distances.
S1Angle kGetDistanceAbsError = S1Angle::Radians(3 * DBL_EPSILON);

TEST(S2EdgeUtil, IntersectionError) {
  // We repeatedly construct two edges that cross near a random point "p", and
  // measure the distance from the actual intersection point "x" to the
  // exact intersection point and also to the edges.

  S1Angle max_point_dist, max_edge_dist;
  S2Testing::Random* rnd = &S2Testing::rnd;
  for (int iter = 0; iter < 5000; ++iter) {
    // We construct two edges AB and CD that intersect near "p".  The angle
    // between AB and CD (expressed as a slope) is chosen randomly between
    // 1e-15 and 1e15 such that its logarithm is uniformly distributed.
    // Similarly, two edge lengths approximately between 1e-15 and 1 are
    // chosen.  The edge endpoints are chosen such that they are often very
    // close to the other edge (i.e., barely crossing).  Taken together this
    // ensures that we test both long and very short edges that intersect at
    // both large and very small angles.
    //
    // Sometimes the edges we generate will not actually cross, in which case
    // we simply try again.
    Vector3_d p, d1, d2;
    S2Testing::GetRandomFrame(&p, &d1, &d2);
    double slope = 1e-15 * pow(1e30, rnd->RandDouble());
    d2 = (d1 + slope * d2).Normalize();
    S2Point a, b, c, d;
    do {
      double ab_len = pow(1e-15, rnd->RandDouble());
      double cd_len = pow(1e-15, rnd->RandDouble());
      double a_fraction = pow(1e-5, rnd->RandDouble());
      if (rnd->OneIn(2)) a_fraction = 1 - a_fraction;
      double c_fraction = pow(1e-5, rnd->RandDouble());
      if (rnd->OneIn(2)) c_fraction = 1 - c_fraction;
      a = (p - a_fraction * ab_len * d1).Normalize();
      b = (p + (1 - a_fraction) * ab_len * d1).Normalize();
      c = (p - c_fraction * cd_len * d2).Normalize();
      d = (p + (1 - c_fraction) * cd_len * d2).Normalize();
    } while (S2EdgeUtil::RobustCrossing(a, b, c, d) <= 0);

    // Each constructed edge should be at most 1.5 * DBL_EPSILON away from the
    // original point P.
    EXPECT_LE(S2EdgeUtil::GetDistance(p, a, b),
              S1Angle::Radians(1.5 * DBL_EPSILON) + kGetDistanceAbsError);
    EXPECT_LE(S2EdgeUtil::GetDistance(p, c, d),
              S1Angle::Radians(1.5 * DBL_EPSILON) + kGetDistanceAbsError);

    // Verify that the expected intersection point is close to both edges and
    // also close to the original point P.  (It might not be very close to P
    // if the angle between the edges is very small.)
    S2Point expected = S2EdgeUtilTesting::GetIntersectionExact(a, b, c, d);
    EXPECT_LE(S2EdgeUtil::GetDistance(expected, a, b),
              S1Angle::Radians(3 * DBL_EPSILON) + kGetDistanceAbsError);
    EXPECT_LE(S2EdgeUtil::GetDistance(expected, c, d),
              S1Angle::Radians(3 * DBL_EPSILON) + kGetDistanceAbsError);
    EXPECT_LE(S1Angle(expected, p),
              S1Angle::Radians(3 * DBL_EPSILON / slope) +
              S2EdgeUtil::kIntersectionError);

    // Now we actually test the GetIntersection() method.
    S2Point actual = S2EdgeUtilTesting::GetIntersection(a, b, c, d);
    S1Angle dist_ab = S2EdgeUtil::GetDistance(actual, a, b);
    S1Angle dist_cd = S2EdgeUtil::GetDistance(actual, c, d);
    EXPECT_LE(dist_ab, S2EdgeUtil::kIntersectionError + kGetDistanceAbsError);
    EXPECT_LE(dist_cd, S2EdgeUtil::kIntersectionError + kGetDistanceAbsError);
    max_edge_dist = max(max_edge_dist, max(dist_ab, dist_cd));
    S1Angle point_dist(expected, actual);
    EXPECT_LE(point_dist, S2EdgeUtil::kIntersectionError);
    max_point_dist = max(max_point_dist, point_dist);
  }
  S2EdgeUtilTesting::PrintIntersectionStats();
  LOG(INFO) << "Max distance to either edge being intersected: "
            << max_edge_dist.radians();
  LOG(INFO) << "Maximum distance to expected intersection point: "
            << max_point_dist.radians();
}

// Chooses a point in the XY plane that is separated from X by at least 1e-15
// (to avoid choosing too many duplicate points) and by at most Pi/2 - 1e-3
// (to avoid nearly-diametric edges, since the test below is not sophisticated
// enough to test such edges).
static S2Point ChooseSemicirclePoint(S2Point const& x, S2Point const& y) {
  S2Testing::Random* rnd = &S2Testing::rnd;
  double sign = (2 * rnd->Uniform(2)) - 1;
  return (x + sign * 1e3 * pow(1e-18, rnd->RandDouble()) * y).Normalize();
}

TEST(S2EdgeUtil, GrazingIntersections) {
  // This test choose 5 points along a great circle (i.e., as collinear as
  // possible), and uses them to construct an edge AB and a triangle CDE such
  // that CD and CE both cross AB.  It then checks that the intersection
  // points returned by GetIntersection() have the correct relative ordering
  // along AB (to within kIntersectionError).
  for (int iter = 0; iter < 1000; ++iter) {
    Vector3_d x, y, z;
    S2Testing::GetRandomFrame(&x, &y, &z);
    S2Point a, b, c, d, e, ab;
    do {
      a = ChooseSemicirclePoint(x, y);
      b = ChooseSemicirclePoint(x, y);
      c = ChooseSemicirclePoint(x, y);
      d = ChooseSemicirclePoint(x, y);
      e = ChooseSemicirclePoint(x, y);
      ab = (a - b).CrossProd(a + b);
    } while (ab.Norm() < 50 * DBL_EPSILON ||
             S2EdgeUtil::RobustCrossing(a, b, c, d) <= 0 ||
             S2EdgeUtil::RobustCrossing(a, b, c, e) <= 0);
    S2Point xcd = S2EdgeUtilTesting::GetIntersection(a, b, c, d);
    S2Point xce = S2EdgeUtilTesting::GetIntersection(a, b, c, e);
    // Essentially this says that if CDE and CAB have the same orientation,
    // then CD and CE should intersect along AB in that order.
    ab = ab.Normalize();
    if (S1Angle(xcd, xce) > 2 * S2EdgeUtil::kIntersectionError) {
      EXPECT_EQ(S2::RobustCCW(c, d, e) == S2::RobustCCW(c, a, b),
                S2::RobustCCW(ab, xcd, xce) > 0);
    }
  }
  S2EdgeUtilTesting::PrintIntersectionStats();
}

TEST(S2EdgeUtil, GetIntersectionInvariants) {
  // Test that the result of GetIntersection does not change when the edges
  // are swapped and/or reversed.  The number of iterations is high because it
  // is difficult to generate test cases that show that CompareEdges() is
  // necessary and correct, for example.
  const int kIters = google::DEBUG_MODE ? 5000 : 50000;
  for (int iter = 0; iter < kIters; ++iter) {
    S2Point a, b, c, d;
    do {
      // GetIntersectionStable() sorts the two edges by length, so construct
      // edges (a,b) and (c,d) that cross and have exactly the same length.
      // This can be done by swapping the "x" and "y" coordinates.
      // [Swapping other coordinate pairs doesn't work because it changes the
      // order of addition in Norm2() == (x**2 + y**2) + z**2.]
      a = c = S2Testing::RandomPoint();
      b = d = S2Testing::RandomPoint();
      swap(c[0], c[1]);
      swap(d[0], d[1]);
    } while (S2EdgeUtil::RobustCrossing(a, b, c, d) <= 0);
    EXPECT_EQ((a - b).Norm2(), (c - d).Norm2());

    // Now verify that GetIntersection returns exactly the same result when
    // the edges are swapped and/or reversed.
    S2Point result = S2EdgeUtil::GetIntersection(a, b, c, d);
    if (S2Testing::rnd.OneIn(2)) { swap(a, b); }
    if (S2Testing::rnd.OneIn(2)) { swap(c, d); }
    if (S2Testing::rnd.OneIn(2)) { swap(a, c); swap(b, d); }
    EXPECT_EQ(result, S2EdgeUtil::GetIntersection(a, b, c, d));
  }
}

bool IsEdgeBNearEdgeA(string const& a_str, const string& b_str,
                      double max_error_degrees) {
  unique_ptr<S2Polyline> a(s2textformat::MakePolyline(a_str));
  EXPECT_EQ(2, a->num_vertices());
  unique_ptr<S2Polyline> b(s2textformat::MakePolyline(b_str));
  EXPECT_EQ(2, b->num_vertices());
  return S2EdgeUtil::IsEdgeBNearEdgeA(a->vertex(0), a->vertex(1),
                                      b->vertex(0), b->vertex(1),
                                      S1Angle::Degrees(max_error_degrees));
}

TEST(S2EdgeUtil, EdgeBNearEdgeA) {
  // Edge is near itself.
  EXPECT_TRUE(IsEdgeBNearEdgeA("5:5, 10:-5", "5:5, 10:-5", 1e-6));

  // Edge is near its reverse
  EXPECT_TRUE(IsEdgeBNearEdgeA("5:5, 10:-5", "10:-5, 5:5", 1e-6));

  // Short edge is near long edge.
  EXPECT_TRUE(IsEdgeBNearEdgeA("10:0, -10:0", "2:1, -2:1", 1.0));

  // Long edges cannot be near shorter edges.
  EXPECT_FALSE(IsEdgeBNearEdgeA("2:1, -2:1", "10:0, -10:0", 1.0));

  // Orthogonal crossing edges are not near each other...
  EXPECT_FALSE(IsEdgeBNearEdgeA("10:0, -10:0", "0:1.5, 0:-1.5", 1.0));

  // ... unless all points on B are within tolerance of A.
  EXPECT_TRUE(IsEdgeBNearEdgeA("10:0, -10:0", "0:1.5, 0:-1.5", 2.0));

  // Very long edges whose endpoints are close may have interior points that are
  // far apart.  An implementation that only considers the vertices of polylines
  // will incorrectly consider such edges as "close" when they are not.
  // Consider, for example, two consecutive lines of longitude.  As they
  // approach the poles, they become arbitrarily close together, but along the
  // equator they bow apart.
  EXPECT_FALSE(IsEdgeBNearEdgeA("89:1, -89:1", "89:2, -89:2", 0.5));
  EXPECT_TRUE(IsEdgeBNearEdgeA("89:1, -89:1", "89:2, -89:2", 1.5));

  // The two arcs here are nearly as long as S2 edges can be (just shy of 180
  // degrees), and their endpoints are less than 1 degree apart.  Their
  // midpoints, however, are at opposite ends of the sphere along its equator.
  EXPECT_FALSE(IsEdgeBNearEdgeA(
                   "0:-179.75, 0:-0.25", "0:179.75, 0:0.25", 1.0));

  // At the equator, the second arc here is 9.75 degrees from the first, and
  // closer at all other points.  However, the southern point of the second arc
  // (-1, 9.75) is too far from the first arc for the short-circuiting logic in
  // IsEdgeBNearEdgeA to apply.
  EXPECT_TRUE(IsEdgeBNearEdgeA("40:0, -5:0", "39:0.975, -1:0.975", 1.0));

  // Same as above, but B's orientation is reversed, causing the angle between
  // the normal vectors of circ(B) and circ(A) to be (180-9.75) = 170.5 degrees,
  // not 9.75 degrees.  The greatest separation between the planes is still 9.75
  // degrees.
  EXPECT_TRUE(IsEdgeBNearEdgeA("10:0, -10:0", "-.4:0.975, 0.4:0.975", 1.0));

  // A and B are on the same great circle, A and B partially overlap, but the
  // only part of B that does not overlap A is shorter than tolerance.
  EXPECT_TRUE(IsEdgeBNearEdgeA("0:0, 1:0", "0.9:0, 1.1:0", 0.25));

  // A and B are on the same great circle, all points on B are close to A at its
  // second endpoint, (1,0).
  EXPECT_TRUE(IsEdgeBNearEdgeA("0:0, 1:0", "1.1:0, 1.2:0", 0.25));

  // Same as above, but B's orientation is reversed.  This case is special
  // because the projection of the normal defining A onto the plane containing B
  // is the null vector, and must be handled by a special case.
  EXPECT_TRUE(IsEdgeBNearEdgeA("0:0, 1:0", "1.2:0, 1.1:0", 0.25));
}


TEST(S2EdgeUtil, CollinearEdgesThatDontTouch) {
  const int kIters = 500;
  for (int iter = 0; iter < kIters; ++iter) {
    S2Point a = S2Testing::RandomPoint();
    S2Point d = S2Testing::RandomPoint();
    S2Point b = S2EdgeUtil::Interpolate(0.05, a, d);
    S2Point c = S2EdgeUtil::Interpolate(0.95, a, d);
    EXPECT_GT(0, S2EdgeUtil::RobustCrossing(a, b, c, d));
    EXPECT_GT(0, S2EdgeUtil::RobustCrossing(a, b, c, d));
    S2EdgeUtil::EdgeCrosser crosser(&a, &b, &c);
    EXPECT_GT(0, crosser.RobustCrossing(&d));
    EXPECT_GT(0, crosser.RobustCrossing(&c));
  }
}


TEST(S2EdgeUtil, CoincidentZeroLengthEdgesThatDontTouch) {
  // It is important that the edge primitives can handle vertices that exactly
  // exactly proportional to each other, i.e. that are not identical but are
  // nevertheless exactly coincident when projected onto the unit sphere.
  // There are various ways that such points can arise.  For example,
  // Normalize() itself is not idempotent: there exist distinct points A,B
  // such that Normalize(A) == B  and Normalize(B) == A.  Another issue is
  // that sometimes calls to Normalize() are skipped when the result of a
  // calculation "should" be unit length mathematically (e.g., when computing
  // the cross product of two orthonormal vectors).
  //
  // This test checks pairs of edges AB and CD where A,B,C,D are exactly
  // coincident on the sphere and the norms of A,B,C,D are monotonically
  // increasing.  Such edge pairs should never intersect.  (This is not
  // obvious, since it depends on the particular symbolic perturbations used
  // by S2::RobustCCW().  It would be better to replace this with a test that
  // says that the CCW results must be consistent with each other.)
  const int kIters = 1000;
  for (int iter = 0; iter < kIters; ++iter) {
    // Construct a point P where every component is zero or a power of 2.
    S2Point p;
    for (int i = 0; i < 3; ++i) {
      int binary_exp = S2Testing::rnd.Skewed(11);
      p[i] = (binary_exp > 1022) ? 0 : pow(2, -binary_exp);
    }
    // If all components were zero, try again.  Note that normalization may
    // convert a non-zero point into a zero one due to underflow (!)
    p = p.Normalize();
    if (p == S2Point(0, 0, 0)) { --iter; continue; }

    // Now every non-zero component should have exactly the same mantissa.
    // This implies that if we scale the point by an arbitrary factor, every
    // non-zero component will still have the same mantissa.  Scale the points
    // so that they are all distinct and are still very likely to satisfy
    // S2::IsUnitLength (which allows for a small amount of error in the norm).
    S2Point a = (1-3e-16) * p;
    S2Point b = (1-1e-16) * p;
    S2Point c = p;
    S2Point d = (1+2e-16) * p;
    if (!S2::IsUnitLength(a) || !S2::IsUnitLength(d)) {
      --iter;
      continue;
    }
    // Verify that the expected edges do not cross.
    EXPECT_GT(0, S2EdgeUtil::RobustCrossing(a, b, c, d));
    S2EdgeUtil::EdgeCrosser crosser(&a, &b, &c);
    EXPECT_GT(0, crosser.RobustCrossing(&d));
    EXPECT_GT(0, crosser.RobustCrossing(&c));
  }
}

void TestFaceClipping(S2Point const& a_raw, S2Point const& b_raw) {
  S2Point a = a_raw.Normalize();
  S2Point b = b_raw.Normalize();
  // TODO(ericv): Remove the following line once S2::RobustCrossProd is
  // extended to use simulation of simplicity.
  if (a == -b) return;

  // First we test GetFaceSegments.
  S2EdgeUtil::FaceSegmentVector segments;
  S2EdgeUtil::GetFaceSegments(a, b, &segments);
  int n = segments.size();
  EXPECT_GE(n, 1);

  ::testing::Message msg;
  msg << "\nA=" << a_raw << "\nB=" << b_raw;
  msg << "\nN=" << S2::RobustCrossProd(a, b) << "\nSegments:\n";
  for (int i = 0; i < n; ++i) {
    S2EdgeUtil::FaceSegment const& s = segments[i];
    msg << i << ": face=" << s.face << ", a=" << s.a << ", b=" << s.b << "\n";
  }
  SCOPED_TRACE(msg);

  R2Rect biunit(R1Interval(-1, 1), R1Interval(-1, 1));
  double const kErrorRadians = S2EdgeUtil::kFaceClipErrorRadians;

  // The first and last vertices should approximately equal A and B.
  EXPECT_LE(a.Angle(S2::FaceUVtoXYZ(segments[0].face, segments[0].a)),
            kErrorRadians);
  EXPECT_LE(b.Angle(S2::FaceUVtoXYZ(segments[n-1].face, segments[n-1].b)),
            kErrorRadians);

  S2Point norm = S2::RobustCrossProd(a, b).Normalize();
  S2Point a_tangent = norm.CrossProd(a);
  S2Point b_tangent = b.CrossProd(norm);
  for (int i = 0; i < n; ++i) {
    // Vertices may not protrude outside the biunit square.
    EXPECT_TRUE(biunit.Contains(segments[i].a));
    EXPECT_TRUE(biunit.Contains(segments[i].b));
    if (i == 0) continue;

    // The two representations of each interior vertex (on adjacent faces)
    // must correspond to exactly the same S2Point.
    EXPECT_NE(segments[i-1].face, segments[i].face);
    EXPECT_EQ(S2::FaceUVtoXYZ(segments[i-1].face, segments[i-1].b),
              S2::FaceUVtoXYZ(segments[i].face, segments[i].a));

    // Interior vertices should be in the plane containing A and B, and should
    // be contained in the wedge of angles between A and B (i.e., the dot
    // products with a_tangent and b_tangent should be non-negative).
    S2Point p = S2::FaceUVtoXYZ(segments[i].face, segments[i].a).Normalize();
    EXPECT_LE(fabs(p.DotProd(norm)), kErrorRadians);
    EXPECT_GE(p.DotProd(a_tangent), -kErrorRadians);
    EXPECT_GE(p.DotProd(b_tangent), -kErrorRadians);
  }

  // Now we test ClipToPaddedFace (sometimes with a padding of zero).  We do
  // this by defining an (x,y) coordinate system for the plane containing AB,
  // and converting points along the great circle AB to angles in the range
  // [-Pi, Pi].  We then accumulate the angle intervals spanned by each
  // clipped edge; the union over all 6 faces should approximately equal the
  // interval covered by the original edge.
  S2Testing::Random* rnd = &S2Testing::rnd;
  double padding = rnd->OneIn(10) ? 0.0 : 1e-10 * pow(1e-5, rnd->RandDouble());
  S2Point x_axis = a, y_axis = a_tangent;
  S1Interval expected_angles(0, a.Angle(b));
  S1Interval max_angles = expected_angles.Expanded(kErrorRadians);
  S1Interval actual_angles;
  for (int face = 0; face < 6; ++face) {
    R2Point a_uv, b_uv;
    if (S2EdgeUtil::ClipToPaddedFace(a, b, face, padding, &a_uv, &b_uv)) {
      S2Point a_clip = S2::FaceUVtoXYZ(face, a_uv).Normalize();
      S2Point b_clip = S2::FaceUVtoXYZ(face, b_uv).Normalize();
      EXPECT_LE(fabs(a_clip.DotProd(norm)), kErrorRadians);
      EXPECT_LE(fabs(b_clip.DotProd(norm)), kErrorRadians);
      if (a_clip.Angle(a) > kErrorRadians) {
        EXPECT_DOUBLE_EQ(1 + padding, max(fabs(a_uv[0]), fabs(a_uv[1])));
      }
      if (b_clip.Angle(b) > kErrorRadians) {
        EXPECT_DOUBLE_EQ(1 + padding, max(fabs(b_uv[0]), fabs(b_uv[1])));
      }
      double a_angle = atan2(a_clip.DotProd(y_axis), a_clip.DotProd(x_axis));
      double b_angle = atan2(b_clip.DotProd(y_axis), b_clip.DotProd(x_axis));
      // Rounding errors may cause b_angle to be slightly less than a_angle.
      // We handle this by constructing the interval with FromPointPair(),
      // which is okay since the interval length is much less than M_PI.
      S1Interval face_angles = S1Interval::FromPointPair(a_angle, b_angle);
      EXPECT_TRUE(max_angles.Contains(face_angles));
      actual_angles = actual_angles.Union(face_angles);
    }
  }
  EXPECT_TRUE(actual_angles.Expanded(kErrorRadians).Contains(expected_angles));
}

void TestFaceClippingEdgePair(S2Point const& a, S2Point const& b) {
  TestFaceClipping(a, b);
  TestFaceClipping(b, a);
}

// This function is designed to choose line segment endpoints that are
// difficult to handle correctly.  Given two adjacent cube vertices P and Q,
// it returns either an edge midpoint, face midpoint, or corner vertex that is
// in the plane of PQ and that has been perturbed slightly.  It also sometimes
// returns a random point from anywhere on the sphere.
S2Point PerturbedCornerOrMidpoint(S2Point const& p, S2Point const& q) {
  S2Testing::Random* rnd = &S2Testing::rnd;
  S2Point a = (rnd->Uniform(3) - 1) * p + (rnd->Uniform(3) - 1) * q;
  if (rnd->OneIn(10)) {
    // This perturbation often has no effect except on coordinates that are
    // zero, in which case the perturbed value is so small that operations on
    // it often result in underflow.
    a += pow(1e-300, rnd->RandDouble()) * S2Testing::RandomPoint();
  } else if (rnd->OneIn(2)) {
    // For coordinates near 1 (say > 0.5), this perturbation yields values
    // that are only a few representable values away from the initial value.
    a += 4 * DBL_EPSILON * S2Testing::RandomPoint();
  } else {
    // A perturbation whose magnitude is in the range [1e-25, 1e-10].
    a += 1e-10 * pow(1e-15, rnd->RandDouble()) * S2Testing::RandomPoint();
  }
  if (a.Norm2() < DBL_MIN) {
    // If a.Norm2() is denormalized, Normalize() loses too much precision.
    return PerturbedCornerOrMidpoint(p, q);
  }
  return a;
}

TEST(S2EdgeUtil, FaceClipping) {
  // Start with a few simple cases.
  // An edge that is entirely contained within one cube face:
  TestFaceClippingEdgePair(S2Point(1, -0.5, -0.5), S2Point(1, 0.5, 0.5));
  // An edge that crosses one cube edge:
  TestFaceClippingEdgePair(S2Point(1, 0, 0), S2Point(0, 1, 0));
  // An edge that crosses two opposite edges of face 0:
  TestFaceClippingEdgePair(S2Point(0.75, 0, -1), S2Point(0.75, 0, 1));
  // An edge that crosses two adjacent edges of face 2:
  TestFaceClippingEdgePair(S2Point(1, 0, 0.75), S2Point(0, 1, 0.75));
  // An edges that crosses three cube edges (four faces):
  TestFaceClippingEdgePair(S2Point(1, 0.9, 0.95), S2Point(-1, 0.95, 0.9));

  // Comprehensively test edges that are difficult to handle, especially those
  // that nearly follow one of the 12 cube edges.
  S2Testing::Random* rnd = &S2Testing::rnd;
  R2Rect biunit(R1Interval(-1, 1), R1Interval(-1, 1));
  const int kIters = 1000;  // Test passes with 1e6 iterations
  for (int iter = 0; iter < kIters; ++iter) {
    SCOPED_TRACE(StringPrintf("Iteration %d", iter));
    // Choose two adjacent cube corners P and Q.
    int face = rnd->Uniform(6);
    int i = rnd->Uniform(4);
    int j = (i + 1) & 3;
    S2Point p = S2::FaceUVtoXYZ(face, biunit.GetVertex(i));
    S2Point q = S2::FaceUVtoXYZ(face, biunit.GetVertex(j));

    // Now choose two points that are nearly in the plane of PQ, preferring
    // points that are near cube corners, face midpoints, or edge midpoints.
    S2Point a = PerturbedCornerOrMidpoint(p, q);
    S2Point b = PerturbedCornerOrMidpoint(p, q);
    TestFaceClipping(a, b);
  }
}

// Choose a random point in the rectangle defined by points A and B, sometimes
// returning a point on the edge AB or the points A and B themselves.
R2Point ChooseRectPoint(R2Point const& a, R2Point const& b) {
  S2Testing::Random* rnd = &S2Testing::rnd;
  if (rnd->OneIn(5)) {
    return rnd->OneIn(2) ? a : b;
  } else if (rnd->OneIn(3)) {
    return a + rnd->RandDouble() * (b - a);
  } else {
    return R2Point(rnd->UniformDouble(a[0], b[0]),
                   rnd->UniformDouble(a[1], b[1]));
  }
}

// Given a point X on the line AB (which is checked), return the fraction "t"
// such that x = (1-t)*a + t*b.  Return 0 if A = B.
double GetFraction(R2Point const& x, R2Point const& a, R2Point const& b) {
  // A bound for the error in edge clipping plus the error in the calculation
  // below (which is similar to IntersectsRect).
  double const kError = (S2EdgeUtil::kEdgeClipErrorUVDist +
                         S2EdgeUtil::kIntersectsRectErrorUVDist);
  if (a == b) return 0.0;
  R2Point dir = (b - a).Normalize();
  EXPECT_LE(fabs((x - a).DotProd(dir.Ortho())), kError);
  return (x - a).DotProd(dir);
}

// Given a point P representing a possibly clipped endpoint A of an edge AB,
// verify that "clip" contains P, and that if clipping occurred (i.e., P != A)
// then P is on the boundary of "clip".
void CheckPointOnBoundary(R2Point const& p, R2Point const& a,
                          R2Rect const& clip) {
  EXPECT_TRUE(clip.Contains(p));
  if (p != a) {
    EXPECT_FALSE(clip.Contains(R2Point(nextafter(p[0], a[0]),
                                       nextafter(p[1], a[1]))));
  }
}

// Given an edge AB and a rectangle "clip", verify that IntersectsRect(),
// ClipEdge(), and ClipEdgeBound() produce consistent results.
void TestClipEdge(R2Point const& a, R2Point const& b, R2Rect const& clip) {
  // A bound for the error in edge clipping plus the error in the
  // IntersectsRect calculation below.
  double const kError = (S2EdgeUtil::kEdgeClipErrorUVDist +
                         S2EdgeUtil::kIntersectsRectErrorUVDist);
  R2Point a_clipped, b_clipped;
  if (!S2EdgeUtil::ClipEdge(a, b, clip, &a_clipped, &b_clipped)) {
    EXPECT_FALSE(S2EdgeUtil::IntersectsRect(a, b, clip.Expanded(-kError)));
  } else {
    EXPECT_TRUE(S2EdgeUtil::IntersectsRect(a, b, clip.Expanded(kError)));
    // Check that the clipped points lie on the edge AB, and that the points
    // have the expected order along the segment AB.
    EXPECT_LE(GetFraction(a_clipped, a, b), GetFraction(b_clipped, a, b));
    // Check that the clipped portion of AB is as large as possible.
    CheckPointOnBoundary(a_clipped, a, clip);
    CheckPointOnBoundary(b_clipped, b, clip);
  }
  // Choose a random initial bound to pass to ClipEdgeBound.
  R2Rect initial_clip = R2Rect::FromPointPair(ChooseRectPoint(a, b),
                                              ChooseRectPoint(a, b));
  R2Rect bound = S2EdgeUtil::GetClippedEdgeBound(a, b, initial_clip);
  if (bound.is_empty()) return;  // Precondition of ClipEdgeBound not met
  R2Rect max_bound = bound.Intersection(clip);
  if (!S2EdgeUtil::ClipEdgeBound(a, b, clip, &bound)) {
    EXPECT_FALSE(S2EdgeUtil::IntersectsRect(a, b, max_bound.Expanded(-kError)));
  } else {
    EXPECT_TRUE(S2EdgeUtil::IntersectsRect(a, b, max_bound.Expanded(kError)));
    // Check that the bound is as large as possible.
    int ai = (a[0] > b[0]), aj = (a[1] > b[1]);
    CheckPointOnBoundary(bound.GetVertex(ai, aj), a, max_bound);
    CheckPointOnBoundary(bound.GetVertex(1-ai, 1-aj), b, max_bound);
  }
}

// Given an interval "clip", randomly choose either a value in the interval, a
// value outside the interval, or one of the two interval endpoints, ensuring
// that all cases have reasonable probability for any interval "clip".
double ChooseEndpoint(R1Interval const& clip) {
  S2Testing::Random* rnd = &S2Testing::rnd;
  if (rnd->OneIn(5)) {
    return rnd->OneIn(2) ? clip.lo() : clip.hi();
  } else {
    switch (rnd->Uniform(3)) {
      case 0:  return clip.lo() - rnd->RandDouble();
      case 1:  return clip.hi() + rnd->RandDouble();
      default: return clip.lo() + rnd->RandDouble() * clip.GetLength();
    }
  }
}

// Given a rectangle "clip", choose a point that may lie in the rectangle
// interior, along an extended edge, exactly at a vertex, or in one of the
// eight regions exterior to "clip" that are separated by its extended edges.
// Also sometimes return points that are exactly on one of the extended
// diagonals of "clip".  All cases are reasonably likely to occur for any
// given rectangle "clip".
R2Point ChooseEndpoint(R2Rect const& clip) {
  if (S2Testing::rnd.OneIn(10)) {
    // Return a point on one of the two extended diagonals.
    int diag = S2Testing::rnd.Uniform(2);
    double t = S2Testing::rnd.UniformDouble(-1, 2);
    return (1 - t) * clip.GetVertex(diag) + t * clip.GetVertex(diag + 2);
  } else {
    return R2Point(ChooseEndpoint(clip[0]), ChooseEndpoint(clip[1]));
  }
}

// Given a rectangle "clip", test the S2EdgeUtil edge clipping methods using
// many edges that are randomly constructed to trigger special cases.
void TestEdgeClipping(R2Rect const& clip) {
  const int kIters = 1000;  // Test passes with 1e6 iterations
  for (int iter = 0; iter < kIters; ++iter) {
    SCOPED_TRACE(StringPrintf("Iteration %d", iter));
    TestClipEdge(ChooseEndpoint(clip), ChooseEndpoint(clip), clip);
  }
}

TEST(S2EdgeUtil, EdgeClipping) {
  S2Testing::Random* rnd = &S2Testing::rnd;
  // Test clipping against random rectangles.
  for (int i = 0; i < 5; ++i) {
    TestEdgeClipping(R2Rect::FromPointPair(
        R2Point(rnd->UniformDouble(-1, 1), rnd->UniformDouble(-1, 1)),
        R2Point(rnd->UniformDouble(-1, 1), rnd->UniformDouble(-1, 1))));
  }
  // Also clip against one-dimensional, singleton, and empty rectangles.
  TestEdgeClipping(R2Rect(R1Interval(-0.7, -0.7), R1Interval(0.3, 0.35)));
  TestEdgeClipping(R2Rect(R1Interval(0.2, 0.5), R1Interval(0.3, 0.3)));
  TestEdgeClipping(R2Rect(R1Interval(-0.7, 0.3), R1Interval(0, 0)));
  TestEdgeClipping(R2Rect::FromPoint(R2Point(0.3, 0.8)));
  TestEdgeClipping(R2Rect::Empty());
}

