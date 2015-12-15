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
#include <algorithm>
#include <cmath>
#include <limits>

#include <glog/logging.h>

#include "r1interval.h"
#include "s1chordangle.h"
#include "util/math/exactfloat/exactfloat.h"
#include "util/math/vector3.h"

// Avoid "using std::abs" to avoid possible programmer confusion with the
// integer-only C function.
using std::max;
using std::min;
using std::numeric_limits;
using std::swap;

// Error constant definitions.  See the header file for details.
double const S2EdgeUtil::kFaceClipErrorRadians = 3 * DBL_EPSILON;
double const S2EdgeUtil::kFaceClipErrorUVDist = 9 * DBL_EPSILON;
double const S2EdgeUtil::kFaceClipErrorUVCoord = 9 * M_SQRT1_2 * DBL_EPSILON;
double const S2EdgeUtil::kIntersectsRectErrorUVDist = 3 * M_SQRT2 * DBL_EPSILON;
double const S2EdgeUtil::kEdgeClipErrorUVCoord = 2.25 * DBL_EPSILON;
double const S2EdgeUtil::kEdgeClipErrorUVDist = 2.25 * DBL_EPSILON;

// kIntersectionError can be set somewhat arbitrarily, because the algorithm
// uses more precision if necessary in order to achieve the specified error.
// The only strict requirement is that kIntersectionError >= DBL_EPSILON
// radians.  However, using a larger error tolerance makes the algorithm more
// efficient because it reduces the number of cases where exact arithmetic is
// needed.
S1Angle const S2EdgeUtil::kIntersectionError = S1Angle::Radians(4*DBL_EPSILON);

// This value can be supplied as the vertex_merge_radius() to S2PolygonBuilder
// to ensure that intersection points that are supposed to be coincident are
// merged back together into a single vertex.  This is required in order for
// various polygon operations (union, intersection, etc) to work correctly.
// It is twice the intersection error because two coincident intersection
// points might have errors in opposite directions.
S1Angle const S2EdgeUtil::kIntersectionMergeRadius = 2 * kIntersectionError;

// DEPRECATED.  This value is related to an obsolete intersection algorithm
// that made weaker accuracy guarantees.
S1Angle const S2EdgeUtil::kIntersectionTolerance = S1Angle::Radians(1.5e-15);

S1Angle const S2EdgeUtil::kIntersectionExactError =
    S1Angle::Radians(DBL_EPSILON);

S2EdgeUtil::IntersectionMethod S2EdgeUtil::last_intersection_method_;

char const* S2EdgeUtil::GetIntersectionMethodName(IntersectionMethod method) {
  switch (method) {
    case SIMPLE:    return "Simple";
    case SIMPLE_LD: return "Simple_ld";
    case STABLE:    return "Stable";
    case STABLE_LD: return "Stable_ld";
    case EXACT:     return "Exact";
    default:        return "Unknown";
  }
}

bool S2EdgeUtil::SimpleCrossing(S2Point const& a, S2Point const& b,
                                S2Point const& c, S2Point const& d) {
  // We compute SimpleCCW() for triangles ACB, CBD, BDA, and DAC.  All
  // of these triangles need to have the same orientation (CW or CCW)
  // for an intersection to exist.  Note that this is slightly more
  // restrictive than the corresponding definition for planar edges,
  // since we need to exclude pairs of line segments that would
  // otherwise "intersect" by crossing two antipodal points.

  Vector3_d ab = a.CrossProd(b);
  double acb = -(ab.DotProd(c));
  double bda = ab.DotProd(d);
  if (acb * bda <= 0) return false;

  Vector3_d cd = c.CrossProd(d);
  double cbd = -(cd.DotProd(b));
  double dac = cd.DotProd(a);
  return (acb * cbd > 0) && (acb * dac > 0);
}

int S2EdgeUtil::RobustCrossing(S2Point const& a, S2Point const& b,
                               S2Point const& c, S2Point const& d) {
  S2EdgeUtil::EdgeCrosser crosser(&a, &b, &c);
  return crosser.RobustCrossing(&d);
}

bool S2EdgeUtil::VertexCrossing(S2Point const& a, S2Point const& b,
                                S2Point const& c, S2Point const& d) {
  // If A == B or C == D there is no intersection.  We need to check this
  // case first in case 3 or more input points are identical.
  if (a == b || c == d) return false;

  // If any other pair of vertices is equal, there is a crossing if and only
  // if OrderedCCW() indicates that the edge AB is further CCW around the
  // shared vertex O (either A or B) than the edge CD, starting from an
  // arbitrary fixed reference point.
  if (a == d) return S2::OrderedCCW(S2::Ortho(a), c, b, a);
  if (b == c) return S2::OrderedCCW(S2::Ortho(b), d, a, b);
  if (a == c) return S2::OrderedCCW(S2::Ortho(a), d, b, a);
  if (b == d) return S2::OrderedCCW(S2::Ortho(b), c, a, b);

  LOG(DFATAL) << "VertexCrossing called with 4 distinct vertices";
  return false;
}

bool S2EdgeUtil::EdgeOrVertexCrossing(S2Point const& a, S2Point const& b,
                                      S2Point const& c, S2Point const& d) {
  int crossing = RobustCrossing(a, b, c, d);
  if (crossing < 0) return false;
  if (crossing > 0) return true;
  return VertexCrossing(a, b, c, d);
}

typedef Vector3<long double> Vector3_ld;
typedef Vector3<ExactFloat> Vector3_xf;

// Computes the cross product of "x" and "y", normalizes it to be unit length,
// and stores the result in "result".  Also returns the length of the cross
// product before normalization, which is useful for estimating the amount of
// error in the result.  For numerical stability, "x" and "y" should both be
// approximately unit length.
template <class T>
static T RobustNormalWithLength(Vector3<T> const& x, Vector3<T> const& y,
                                Vector3<T>* result) {
  // This computes 2 * (x.CrossProd(y)), but has much better numerical
  // stability when "x" and "y" are unit length.
  Vector3<T> tmp = (x - y).CrossProd(x + y);
  T length = tmp.Norm();
  if (length != 0) {
    *result = (1 / length) * tmp;
  }
  return 0.5 * length;  // Since tmp == 2 * (x.CrossProd(y))
}

// If the intersection point of the edges (a0,a1) and (b0,b1) can be computed
// to within an error of at most S2EdgeUtil::kIntersectionError by this
// function, then set "result" to the intersection point and return true.
//
// The intersection point is not guaranteed to have the correct sign
// (i.e., it may be either "result" or "-result").
template <class T>
static bool GetIntersectionSimple(Vector3<T> const& a0, Vector3<T> const& a1,
                                  Vector3<T> const& b0, Vector3<T> const& b1,
                                  Vector3<T>* result) {
  // The code below computes the intersection point as
  //
  //    (a0.CrossProd(a1)).CrossProd(b0.CrossProd(b1))
  //
  // except that it has better numerical stability and also computes a
  // guaranteed error bound.
  //
  // Each cross product is computed as (X-Y).CrossProd(X+Y) using unit-length
  // input vectors, which eliminates most of the cancellation error.  However
  // the error in the direction of the cross product can still become large if
  // the two points are extremely close together.  We can show that as long as
  // the length of the cross product is at least (8*sqrt(3)+12) * DBL_EPSILON
  // (about 6e-15), then the directional error is at most 2.5 *
  // numeric_limits<T>::epsilon() (about 3e-19 when T == "long double").
  // (DBL_EPSILON appears in the first formula because the inputs are assumed
  // to be normalized in double precision rather than in the given type T.)
  //
  // The third cross product is different because its inputs already have some
  // error.  Letting "result_len" be the length of the cross product, it can
  // be shown that the error is at most
  //
  //   (1 + sqrt(3) + 6 / result_len) * numeric_limits<T>::epsilon()
  //
  // We want this error to be at most kIntersectionError, which is true as
  // long as "result_len" is at least kMinResultLen defined below.

  static const T kMinNormalLength = (8 * sqrt(3) + 12) * DBL_EPSILON;
  static const T kMinResultLen =
      6 / (S2EdgeUtil::kIntersectionError.radians() /
           numeric_limits<T>::epsilon() - (1 + sqrt(3)));

  // On some platforms "long double" is the same as "double", and on these
  // platforms this method always returns false (e.g. ARM, Win32).  Rather
  // than testing this directly, instead we look at kMinResultLen since this
  // is a direct measure of whether "long double" has sufficient accuracy to
  // be useful.  If kMinResultLen > 0.5, it means that this method will fail
  // even for edges that meet at an angle of 30 degrees.  (On Intel platforms
  // kMinResultLen corresponds to an intersection angle of about 0.04
  // degrees.)
  DCHECK_LE(kMinResultLen, 0.5);

  Vector3<T> a_norm, b_norm;
  if (RobustNormalWithLength(a0, a1, &a_norm) >= kMinNormalLength &&
      RobustNormalWithLength(b0, b1, &b_norm) >= kMinNormalLength &&
      RobustNormalWithLength(a_norm, b_norm, result) >= kMinResultLen) {
    return true;
  }
  return false;
}

static bool GetIntersectionSimpleLD(S2Point const& a0, S2Point const& a1,
                                    S2Point const& b0, S2Point const& b1,
                                    S2Point* result) {
  Vector3_ld result_ld;
  if (GetIntersectionSimple(Vector3_ld::Cast(a0), Vector3_ld::Cast(a1),
                            Vector3_ld::Cast(b0), Vector3_ld::Cast(b1),
                            &result_ld)) {
    *result = S2Point::Cast(result_ld);
    return true;
  }
  return false;
}

// Given a point X and a vector "a_norm" (not necessarily unit length),
// compute x.DotProd(a_norm) and return a bound on the error in the result.
// The remaining parameters allow this dot product to be computed more
// accurately and efficiently.  They include the length of "a_norm"
// ("a_norm_len") and the edge endpoints "a0" and "a1".
template <class T>
static T GetProjection(Vector3<T> const& x,
                       Vector3<T> const& a_norm, T a_norm_len,
                       Vector3<T> const& a0, Vector3<T> const& a1,
                       T* error) {
  // The error in the dot product is proportional to the lengths of the input
  // vectors, so rather than using "x" itself (a unit-length vector) we use
  // the vectors from "x" to the closer of the two edge endpoints.  This
  // typically reduces the error by a huge factor.
  Vector3<T> x0 = x - a0;
  Vector3<T> x1 = x - a1;
  T x0_dist2 = x0.Norm2();
  T x1_dist2 = x1.Norm2();

  // If both distances are the same, we need to be careful to choose one
  // endpoint deterministically so that the result does not change if the
  // order of the endpoints is reversed.
  T dist, result;
  if (x0_dist2 < x1_dist2 || (x0_dist2 == x1_dist2 && x0 < x1)) {
    dist = sqrt(x0_dist2);
    result = x0.DotProd(a_norm);
  } else {
    dist = sqrt(x1_dist2);
    result = x1.DotProd(a_norm);
  }
  // This calculation bounds the error from all sources: the computation of
  // the normal, the subtraction of one endpoint, and the dot product itself.
  // (DBL_EPSILON appears because the input points are assumed to be
  // normalized in double precision rather than in the given type T.)
  *error = ((3.5 * a_norm_len + 8 * sqrt(3) * DBL_EPSILON) * dist
            + 0.75 * std::abs(result)) * numeric_limits<T>::epsilon();
  return result;
}

// Helper function for GetIntersectionStable().  It expects that the edges
// (a0,a1) and (b0,b1) have been sorted so that the first edge is longer.
template <class T>
static bool GetIntersectionStableSorted(
    Vector3<T> const& a0, Vector3<T> const& a1,
    Vector3<T> const& b0, Vector3<T> const& b1, Vector3<T>* result) {
  DCHECK_GE((a1 - a0).Norm2(), (b1 - b0).Norm2());

  // Compute the normal of the plane through (a0, a1) in a stable way.
  Vector3<T> a_norm = (a0 - a1).CrossProd(a0 + a1);
  T a_norm_len = a_norm.Norm();
  T b_len = (b1 - b0).Norm();

  // Compute the projection (i.e., signed distance) of b0 and b1 onto the
  // plane through (a0, a1).  Distances are scaled by the length of a_norm.
  T b0_error, b1_error;
  T b0_dist = GetProjection(b0, a_norm, a_norm_len, a0, a1, &b0_error);
  T b1_dist = GetProjection(b1, a_norm, a_norm_len, a0, a1, &b1_error);

  // The total distance from b0 to b1 measured perpendicularly to (a0,a1) is
  // |b0_dist - b1_dist|.  Note that b0_dist and b1_dist generally have
  // opposite signs because b0 and b1 are on opposite sides of (a0, a1).  The
  // code below finds the intersection point by interpolating along the edge
  // (b0, b1) to a fractional distance of b0_dist / (b0_dist - b1_dist).
  //
  // It can be shown that the maximum error in the interpolation fraction is
  //
  //     (b0_dist * b1_error - b1_dist * b0_error) /
  //        (dist_sum * (dist_sum - error_sum))
  //
  // We save ourselves some work by scaling the result and the error bound by
  // "dist_sum", since the result is normalized to be unit length anyway.
  T dist_sum = std::abs(b0_dist - b1_dist);
  T error_sum = b0_error + b1_error;
  if (dist_sum <= error_sum) {
    return false;  // Error is unbounded in this case.
  }
  Vector3<T> x = b0_dist * b1 - b1_dist * b0;
  T error = b_len * std::abs(b0_dist * b1_error - b1_dist * b0_error) /
      (dist_sum - error_sum) + dist_sum * numeric_limits<T>::epsilon();

  // Finally we normalize the result, compute the corresponding error, and
  // check whether the total error is acceptable.
  T x_len = x.Norm();
  T const kMaxError = S2EdgeUtil::kIntersectionError.radians();
  if (error > (kMaxError - 0.5 * numeric_limits<T>::epsilon()) * x_len) {
    return false;
  }
  *result = (1 / x_len) * x;
  return true;
}

// Returns whether (a0,a1) is less than (b0,b1) with respect to a total
// ordering on edges that is invariant under edge reversals.
template <class T>
static bool CompareEdges(Vector3<T> const& a0, Vector3<T> const& a1,
                         Vector3<T> const& b0, Vector3<T> const& b1) {
  Vector3<T> const *pa0 = &a0, *pa1 = &a1;
  Vector3<T> const *pb0 = &b0, *pb1 = &b1;
  if (*pa0 >= *pa1) swap(pa0, pa1);
  if (*pb0 >= *pb1) swap(pb0, pb1);
  return *pa0 < *pb0 || (*pa0 == *pb0 && *pb0 < *pb1);
}

// If the intersection point of the edges (a0,a1) and (b0,b1) can be computed
// to within an error of at most S2EdgeUtil::kIntersectionError by this
// function, then set "result" to the intersection point and return true.
//
// The intersection point is not guaranteed to have the correct sign
// (i.e., it may be either "result" or "-result").
template <class T>
static bool GetIntersectionStable(Vector3<T> const& a0, Vector3<T> const& a1,
                                  Vector3<T> const& b0, Vector3<T> const& b1,
                                  Vector3<T>* result) {
  // Sort the two edges so that (a0,a1) is longer, breaking ties in a
  // deterministic way that does not depend on the ordering of the endpoints.
  // This is desirable for two reasons:
  //  - So that the result doesn't change when edges are swapped or reversed.
  //  - It reduces error, since the first edge is used to compute the edge
  //    normal (where a longer edge means less error), and the second edge
  //    is used for interpolation (where a shorter edge means less error).
  T a_len2 = (a1 - a0).Norm2();
  T b_len2 = (b1 - b0).Norm2();
  if (a_len2 < b_len2 || (a_len2 == b_len2 && CompareEdges(a0, a1, b0, b1))) {
    return GetIntersectionStableSorted(b0, b1, a0, a1, result);
  } else {
    return GetIntersectionStableSorted(a0, a1, b0, b1, result);
  }
}

static bool GetIntersectionStableLD(S2Point const& a0, S2Point const& a1,
                                    S2Point const& b0, S2Point const& b1,
                                    S2Point* result) {
  Vector3_ld result_ld;
  if (GetIntersectionStable(Vector3_ld::Cast(a0), Vector3_ld::Cast(a1),
                            Vector3_ld::Cast(b0), Vector3_ld::Cast(b1),
                            &result_ld)) {
    *result = S2Point::Cast(result_ld);
    return true;
  }
  return false;
}

S2Point S2EdgeUtil::S2PointFromExact(Vector3_xf const& x) {
  // TODO(ericv): In theory, this function may return S2Point(0, 0, 0) even
  // though "x" is non-zero.  This happens when all components of "x" have
  // absolute value less than about 1e-154, since in that case x.Norm2() is
  // zero in double precision.  This could be fixed by scaling "x" by an
  // appropriate power of 2 before the conversion.
  return S2Point(x[0].ToDouble(), x[1].ToDouble(), x[2].ToDouble()).Normalize();
}

// Compute the intersection point of (a0, a1) and (b0, b1) using exact
// arithmetic.  Note that the result is not exact because it is rounded to
// double precision.  Also, the intersection point is not guaranteed to have
// the correct sign (i.e., the return value may need to be negated).
S2Point S2EdgeUtil::GetIntersectionExact(S2Point const& a0, S2Point const& a1,
                                         S2Point const& b0, S2Point const& b1) {
  // Since we are using exact arithmetic, we don't need to worry about
  // numerical stability.
  Vector3_xf a0_xf = Vector3_xf::Cast(a0);
  Vector3_xf a1_xf = Vector3_xf::Cast(a1);
  Vector3_xf b0_xf = Vector3_xf::Cast(b0);
  Vector3_xf b1_xf = Vector3_xf::Cast(b1);
  Vector3_xf a_norm_xf = a0_xf.CrossProd(a1_xf);
  Vector3_xf b_norm_xf = b0_xf.CrossProd(b1_xf);
  Vector3_xf x_xf = a_norm_xf.CrossProd(b_norm_xf);

  // The final Normalize() call is done in double precision, which creates a
  // directional error of up to DBL_EPSILON.  (ToDouble() and Normalize()
  // each contribute up to 0.5 * DBL_EPSILON of directional error.)
  S2Point x = S2EdgeUtil::S2PointFromExact(x_xf);

  if (x == S2Point(0, 0, 0)) {
    // The two edges are exactly collinear, but we still consider them to be
    // "crossing" because of simulation of simplicity.  Out of the four
    // endpoints, exactly two lie in the interior of the other edge.  Of
    // those two we return the one that is lexicographically smallest.
    x = S2Point(10, 10, 10);  // Greater than any valid S2Point
    S2Point a_norm = S2EdgeUtil::S2PointFromExact(a_norm_xf);
    S2Point b_norm = S2EdgeUtil::S2PointFromExact(b_norm_xf);
    if (S2::OrderedCCW(b0, a0, b1, b_norm) && a0 < x) x = a0;
    if (S2::OrderedCCW(b0, a1, b1, b_norm) && a1 < x) x = a1;
    if (S2::OrderedCCW(a0, b0, a1, a_norm) && b0 < x) x = b0;
    if (S2::OrderedCCW(a0, b1, a1, a_norm) && b1 < x) x = b1;
  }
  DCHECK(S2::IsUnitLength(x));
  return x;
}

// Given three points "a", "x", "b", returns true if these three points occur
// in the given order along the edge (a,b) to within the given tolerance.
// More precisely, either "x" must be within "tolerance" of "a" or "b", or
// when "x" is projected onto the great circle through "a" and "b" it must lie
// along the edge (a,b) (i.e., the shortest path from "a" to "b").
static bool ApproximatelyOrdered(S2Point const& a, S2Point const& x,
                                 S2Point const& b, double tolerance) {
  if ((x - a).Norm2() <= tolerance * tolerance) return true;
  if ((x - b).Norm2() <= tolerance * tolerance) return true;
  return S2::OrderedCCW(a, x, b, S2::RobustCrossProd(a, b).Normalize());
}

S2Point S2EdgeUtil::GetIntersection(S2Point const& a0, S2Point const& a1,
                                    S2Point const& b0, S2Point const& b1) {
  DCHECK_GT(RobustCrossing(a0, a1, b0, b1), 0);

  // It is difficult to compute the intersection point of two edges accurately
  // when the angle between the edges is very small.  Previously we handled
  // this by only guaranteeing that the returned intersection point is within
  // kIntersectionError of each edge.  However, this means that when the edges
  // cross at a very small angle, the computed result may be very far from the
  // true intersection point.
  //
  // Instead this function now guarantees that the result is always within
  // kIntersectionError of the true intersection.  This requires using more
  // sophisticated techniques and in some cases extended precision.
  //
  // Three different techniques are implemented, but only two are used:
  //
  //  - GetIntersectionSimple() computes the intersection point using
  //    numerically stable cross products in "long double" precision.
  //
  //  - GetIntersectionStable() computes the intersection point using
  //    projection and interpolation, taking care to minimize cancellation
  //    error.  This method exists in "double" and "long double" versions.
  //
  //  - GetIntersectionExact() computes the intersection point using exact
  //    arithmetic and converts the final result back to an S2Point.
  //
  // We don't actually use the first method (GetIntersectionSimple) because it
  // turns out that GetIntersectionStable() is twice as fast and also much
  // more accurate (even in double precision).  The "long double" version
  // (only available on Intel platforms) uses 80-bit precision and is about
  // twice as slow.  The exact arithmetic version is about 100x slower.
  //
  // So our strategy is to first call GetIntersectionStable() in double
  // precision; if that doesn't work and this platform supports "long double",
  // then we try again in "long double"; if that doesn't work then we fall
  // back to exact arithmetic.

  static bool const kUseSimpleMethod = false;
  static bool const kHasLongDouble = (numeric_limits<long double>::epsilon() <
                                      numeric_limits<double>::epsilon());
  S2Point result;
  if (kUseSimpleMethod && GetIntersectionSimple(a0, a1, b0, b1, &result)) {
    last_intersection_method_ = SIMPLE;
  } else if (kUseSimpleMethod && kHasLongDouble &&
             GetIntersectionSimpleLD(a0, a1, b0, b1, &result)) {
    last_intersection_method_ = SIMPLE_LD;
  } else if (GetIntersectionStable(a0, a1, b0, b1, &result)) {
    last_intersection_method_ = STABLE;
  } else if (kHasLongDouble &&
             GetIntersectionStableLD(a0, a1, b0, b1, &result)) {
    last_intersection_method_ = STABLE_LD;
  } else {
    result = GetIntersectionExact(a0, a1, b0, b1);
    last_intersection_method_ = EXACT;
  }

  // Make sure the intersection point is on the correct side of the sphere.
  // Since all vertices are unit length, and edges are less than 180 degrees,
  // (a0 + a1) and (b0 + b1) both have positive dot product with the
  // intersection point.  We use the sum of all vertices to make sure that the
  // result is unchanged when the edges are swapped or reversed.
  if (result.DotProd((a0 + a1) + (b0 + b1)) < 0) result = -result;

  // Make sure that the intersection point lies on both edges.
  DCHECK(ApproximatelyOrdered(a0, result, a1, kIntersectionError.radians()));
  DCHECK(ApproximatelyOrdered(b0, result, b1, kIntersectionError.radians()));

  return result;
}

double S2EdgeUtil::GetDistanceFraction(S2Point const& x,
                                       S2Point const& a0, S2Point const& a1) {
  DCHECK_NE(a0, a1);
  double d0 = x.Angle(a0);
  double d1 = x.Angle(a1);
  return d0 / (d0 + d1);
}

S2Point S2EdgeUtil::InterpolateAtDistance(S1Angle ax_angle,
                                          S2Point const& a, S2Point const& b) {
  // As of crosstool v14, gcc tries to calculate sin(ax) and cos(ax) using one
  // sincos() call.  However, for some inputs sincos() returns significantly
  // different values between AMD and Intel.
  //
  // As a temporary workaround, "ax" is declared as "volatile" to prohibit the
  // compiler from using sincos(), because sin() and cos() don't seem to have
  // the problem.  See b/3088321 for details.
  volatile double ax = ax_angle.radians();

  DCHECK(S2::IsUnitLength(a));
  DCHECK(S2::IsUnitLength(b));

  // Use RobustCrossProd() to compute the tangent vector at A towards B.  The
  // result is always perpendicular to A, even if A=B or A=-B, but it is not
  // necessarily unit length.  (We effectively normalize it below.)
  Vector3_d normal = S2::RobustCrossProd(a, b);
  Vector3_d tangent = normal.CrossProd(a);
  DCHECK(tangent != S2Point(0, 0, 0));

  // Now compute the appropriate linear combination of A and "tangent".  With
  // infinite precision the result would always be unit length, but we
  // normalize it anyway to ensure that the error is within acceptable bounds.
  // (Otherwise errors can build up when the result of one interpolation is
  // fed into another interpolation.)
  return (cos(ax) * a + (sin(ax) / tangent.Norm()) * tangent).Normalize();
}

S2Point S2EdgeUtil::Interpolate(double t, S2Point const& a, S2Point const& b) {
  if (t == 0) return a;
  if (t == 1) return b;
  S1Angle ab(a, b);
  return InterpolateAtDistance(t * ab, a, b);
}

// This function computes the distance from a point X to a line segment AB.
// If the distance is less than "min_dist" or "always_update" is true, it
// updates "min_dist" and returns true.  Otherwise it returns false.
//
// The "Always" in the function name refers to the template argument, i.e.
// AlwaysUpdateMinDistance<true> always updates the given distance, while
// AlwaysUpdateMinDistance<false> does not.  This optimization increases the
// speed of GetDistance() by about 10% without creating code duplication.
template <bool always_update>
inline bool AlwaysUpdateMinDistance(S2Point const& x,
                                    S2Point const& a, S2Point const& b,
                                    S1ChordAngle* min_dist) {
  DCHECK(S2::IsUnitLength(x) && S2::IsUnitLength(a) && S2::IsUnitLength(b));

  // We divide the problem into two cases, based on whether the closest point
  // on AB is one of the two vertices (the "vertex case") or in the interior
  // (the "interior case").  Let C = A x B.  If X is in the spherical wedge
  // extending from A to B around the axis through C, then we are in the
  // interior case.  Otherwise we are in the vertex case.

  // Check whether we might be in the interior case.  For this to be true, XAB
  // and XBA must both be acute angles.  Checking this condition exactly is
  // expensive, so instead we consider the planar triangle ABX (which passes
  // through the sphere's interior).  The planar angles XAB and XBA are always
  // less than the corresponding spherical angles, so if we are in the
  // interior case then both of these angles must be acute.
  //
  // We check this by computing the squared edge lengths of the planar
  // triangle ABX, and testing acuteness using the law of cosines:
  //
  //             max(XA^2, XB^2) < AB^2 + min(XA^2, XB^2)
  //
  double xa2 = (x-a).Norm2(), xb2 = (x-b).Norm2(), ab2 = (a-b).Norm2();
  double dist2 = min(xa2, xb2);
  if (max(xa2, xb2) < ab2 + dist2) {
    // The minimum distance might be to a point on the edge interior.  Let R
    // be closest point to X that lies on the great circle through AB.  Rather
    // than computing the geodesic distance along the surface of the sphere,
    // instead we compute the "chord length" through the sphere's interior.
    // If the squared chord length exceeds min_dist.length2() then we can
    // return "false" immediately.
    //
    // The squared chord length XR^2 can be expressed as XQ^2 + QR^2, where Q
    // is the point X projected onto the plane through the great circle AB.
    // The distance XQ^2 can be written as (X.C)^2 / |C|^2 where C = A x B.
    // We ignore the QR^2 term and instead use XQ^2 as a lower bound, since it
    // is faster and the corresponding distance on the Earth's surface is
    // accurate to within 1% for distances up to about 1800km.
    S2Point c = S2::RobustCrossProd(a, b);
    double c2 = c.Norm2();
    double x_dot_c = x.DotProd(c);
    double x_dot_c2 = x_dot_c * x_dot_c;
    if (!always_update && x_dot_c2 >= c2 * min_dist->length2()) {
      // The closest point on the great circle AB is too far away.
      return false;
    }
    // Otherwise we do the exact, more expensive test for the interior case.
    // This test is very likely to succeed because of the conservative planar
    // test we did initially.
    S2Point cx = c.CrossProd(x);
    if (a.DotProd(cx) < 0 && b.DotProd(cx) > 0) {
      // Compute the squared chord length XR^2 = XQ^2 + QR^2 (see above).
      // This calculation has good accuracy for all chord lengths since it
      // is based on both the dot product and cross product (rather than
      // deriving one from the other).  However, note that the chord length
      // representation itself loses accuracy as the angle approaches Pi.
      double qr = 1 - sqrt(cx.Norm2() / c2);
      dist2 = (x_dot_c2 / c2) + (qr * qr);
    }
  }
  if (!always_update && dist2 >= min_dist->length2()) return false;
  *min_dist = S1ChordAngle::FromLength2(dist2);
  return true;
}

S1Angle S2EdgeUtil::GetDistance(S2Point const& x,
                                S2Point const& a, S2Point const& b) {
  S1ChordAngle min_dist;
  AlwaysUpdateMinDistance<true>(x, a, b, &min_dist);
  return min_dist.ToAngle();
}

bool S2EdgeUtil::UpdateMinDistance(S2Point const& x,
                                   S2Point const& a, S2Point const& b,
                                   S1ChordAngle* min_dist) {
  return AlwaysUpdateMinDistance<false>(x, a, b, min_dist);
}

S2Point S2EdgeUtil::GetClosestPoint(S2Point const& x,
                                    S2Point const& a, S2Point const& b,
                                    Vector3_d const& a_cross_b) {
  DCHECK(S2::IsUnitLength(a));
  DCHECK(S2::IsUnitLength(b));
  DCHECK(S2::IsUnitLength(x));

  // Find the closest point to X along the great circle through AB.
  S2Point p = x - (x.DotProd(a_cross_b) / a_cross_b.Norm2()) * a_cross_b;

  // If this point is on the edge AB, then it's the closest point.
  if (S2::SimpleCCW(a_cross_b, a, p) && S2::SimpleCCW(p, b, a_cross_b))
    return p.Normalize();

  // Otherwise, the closest point is either A or B.
  return ((x - a).Norm2() <= (x - b).Norm2()) ? a : b;
}

S2Point S2EdgeUtil::GetClosestPoint(S2Point const& x,
                                    S2Point const& a, S2Point const& b) {
  return GetClosestPoint(x, a, b, S2::RobustCrossProd(a, b));
}

bool S2EdgeUtil::UpdateEdgePairMinDistance(
    S2Point const& a0, S2Point const& a1,
    S2Point const& b0, S2Point const& b1,
    S1ChordAngle* min_dist) {
  if (*min_dist == S1ChordAngle::Zero()) {
    return false;
  }
  if (RobustCrossing(a0, a1, b0, b1) > 0) {
    *min_dist = S1ChordAngle::Zero();
    return true;
  }
  // Otherwise, the minimum distance is achieved at an endpoint of at least
  // one of the two edges.  We use "|" rather than "||" below to ensure that
  // all four possibilities are always checked.
  //
  // The calculation below computes each of the six vertex-vertex distances
  // twice (this could be optimized).
  return (UpdateMinDistance(a0, b0, b1, min_dist) |
          UpdateMinDistance(a1, b0, b1, min_dist) |
          UpdateMinDistance(b0, a0, a1, min_dist) |
          UpdateMinDistance(b1, a0, a1, min_dist));
}

std::pair<S2Point, S2Point> S2EdgeUtil::GetEdgePairClosestPoints(
      S2Point const& a0, S2Point const& a1,
      S2Point const& b0, S2Point const& b1) {
  if (RobustCrossing(a0, a1, b0, b1) > 0) {
    S2Point x = GetIntersection(a0, a1, b0, b1);
    return std::make_pair(x, x);
  }
  // We save some work by first determining which vertex/edge pair achieves
  // the minimum distance, and then computing the closest point on that edge.
  S1ChordAngle min_dist;
  AlwaysUpdateMinDistance<true>(a0, b0, b1, &min_dist);
  enum { A0, A1, B0, B1 } closest_vertex = A0;
  if (UpdateMinDistance(a1, b0, b1, &min_dist)) { closest_vertex = A1; }
  if (UpdateMinDistance(b0, a0, a1, &min_dist)) { closest_vertex = B0; }
  if (UpdateMinDistance(b1, a0, a1, &min_dist)) { closest_vertex = B1; }
  switch (closest_vertex) {
    case A0: return std::make_pair(a0, GetClosestPoint(a0, b0, b1));
    case A1: return std::make_pair(a1, GetClosestPoint(a1, b0, b1));
    case B0: return std::make_pair(GetClosestPoint(b0, a0, a1), b0);
    case B1: return std::make_pair(GetClosestPoint(b1, a0, a1), b1);
  }
}

bool S2EdgeUtil::IsEdgeBNearEdgeA(S2Point const& a0, S2Point const& a1,
                                  S2Point const& b0, S2Point const& b1,
                                  S1Angle tolerance) {
  DCHECK_LT(tolerance.radians(), M_PI / 2);
  DCHECK_GT(tolerance.radians(), 0);
  // The point on edge B=b0b1 furthest from edge A=a0a1 is either b0, b1, or
  // some interior point on B.  If it is an interior point on B, then it must be
  // one of the two points where the great circle containing B (circ(B)) is
  // furthest from the great circle containing A (circ(A)).  At these points,
  // the distance between circ(B) and circ(A) is the angle between the planes
  // containing them.

  Vector3_d a_ortho = S2::RobustCrossProd(a0, a1).Normalize();
  S2Point const a_nearest_b0 = GetClosestPoint(b0, a0, a1, a_ortho);
  S2Point const a_nearest_b1 = GetClosestPoint(b1, a0, a1, a_ortho);
  // If a_nearest_b0 and a_nearest_b1 have opposite orientation from a0 and a1,
  // we invert a_ortho so that it points in the same direction as a_nearest_b0 x
  // a_nearest_b1.  This helps us handle the case where A and B are oppositely
  // oriented but otherwise might be near each other.  We check orientation and
  // invert rather than computing a_nearest_b0 x a_nearest_b1 because those two
  // points might be equal, and have an unhelpful cross product.
  if (S2::RobustCCW(a_ortho, a_nearest_b0, a_nearest_b1) < 0)
    a_ortho *= -1;

  // To check if all points on B are within tolerance of A, we first check to
  // see if the endpoints of B are near A.  If they are not, B is not near A.
  S1Angle const b0_distance(b0, a_nearest_b0);
  S1Angle const b1_distance(b1, a_nearest_b1);
  if (b0_distance > tolerance || b1_distance > tolerance)
    return false;

  // If b0 and b1 are both within tolerance of A, we check to see if the angle
  // between the planes containing B and A is greater than tolerance.  If it is
  // not, no point on B can be further than tolerance from A (recall that we
  // already know that b0 and b1 are close to A, and S2Edges are all shorter
  // than 180 degrees).  The angle between the planes containing circ(A) and
  // circ(B) is the angle between their normal vectors.
  Vector3_d const b_ortho = S2::RobustCrossProd(b0, b1).Normalize();
  S1Angle const planar_angle(a_ortho, b_ortho);
  if (planar_angle <= tolerance)
    return true;


  // As planar_angle approaches M_PI, the projection of a_ortho onto the plane
  // of B approaches the null vector, and normalizing it is numerically
  // unstable.  This makes it unreliable or impossible to identify pairs of
  // points where circ(A) is furthest from circ(B).  At this point in the
  // algorithm, this can only occur for two reasons:
  //
  //  1.) b0 and b1 are closest to A at distinct endpoints of A, in which case
  //      the opposite orientation of a_ortho and b_ortho means that A and B are
  //      in opposite hemispheres and hence not close to each other.
  //
  //  2.) b0 and b1 are closest to A at the same endpoint of A, in which case
  //      the orientation of a_ortho was chosen arbitrarily to be that of a0
  //      cross a1.  B must be shorter than 2*tolerance and all points in B are
  //      close to one endpoint of A, and hence to A.
  //
  // The logic applies when planar_angle is robustly greater than M_PI/2, but
  // may be more computationally expensive than the logic beyond, so we choose a
  // value close to M_PI.
  if (planar_angle >= S1Angle::Radians(M_PI - 0.01)) {
    return (S1Angle(b0, a0) < S1Angle(b0, a1)) ==
        (S1Angle(b1, a0) < S1Angle(b1, a1));
  }

  // Finally, if either of the two points on circ(B) where circ(B) is furthest
  // from circ(A) lie on edge B, edge B is not near edge A.
  //
  // The normalized projection of a_ortho onto the plane of circ(B) is one of
  // the two points along circ(B) where it is furthest from circ(A).  The other
  // is -1 times the normalized projection.
  S2Point furthest = (a_ortho - a_ortho.DotProd(b_ortho) * b_ortho).Normalize();
  DCHECK(S2::IsUnitLength(furthest));
  S2Point furthest_inv = -1 * furthest;

  // A point p lies on B if you can proceed from b_ortho to b0 to p to b1 and
  // back to b_ortho without ever turning right.  We test this for furthest and
  // furthest_inv, and return true if neither point lies on B.
  return !((S2::RobustCCW(b_ortho, b0, furthest) > 0 &&
            S2::RobustCCW(furthest, b1, b_ortho) > 0) ||
           (S2::RobustCCW(b_ortho, b0, furthest_inv) > 0 &&
            S2::RobustCCW(furthest_inv, b1, b_ortho) > 0));
}


bool S2EdgeUtil::WedgeContains(S2Point const& a0, S2Point const& ab1,
                               S2Point const& a2, S2Point const& b0,
                               S2Point const& b2) {
  // For A to contain B (where each loop interior is defined to be its left
  // side), the CCW edge order around ab1 must be a2 b2 b0 a0.  We split
  // this test into two parts that test three vertices each.
  return S2::OrderedCCW(a2, b2, b0, ab1) && S2::OrderedCCW(b0, a0, a2, ab1);
}

bool S2EdgeUtil::WedgeIntersects(S2Point const& a0, S2Point const& ab1,
                                 S2Point const& a2, S2Point const& b0,
                                 S2Point const& b2) {
  // For A not to intersect B (where each loop interior is defined to be
  // its left side), the CCW edge order around ab1 must be a0 b2 b0 a2.
  // Note that it's important to write these conditions as negatives
  // (!OrderedCCW(a,b,c,o) rather than Ordered(c,b,a,o)) to get correct
  // results when two vertices are the same.
  return !(S2::OrderedCCW(a0, b2, b0, ab1) && S2::OrderedCCW(b0, a2, a0, ab1));
}

S2EdgeUtil::WedgeRelation S2EdgeUtil::GetWedgeRelation(
    S2Point const& a0, S2Point const& ab1, S2Point const& a2,
    S2Point const& b0, S2Point const& b2) {
  // There are 6 possible edge orderings at a shared vertex (all
  // of these orderings are circular, i.e. abcd == bcda):
  //
  //  (1) a2 b2 b0 a0: A contains B
  //  (2) a2 a0 b0 b2: B contains A
  //  (3) a2 a0 b2 b0: A and B are disjoint
  //  (4) a2 b0 a0 b2: A and B intersect in one wedge
  //  (5) a2 b2 a0 b0: A and B intersect in one wedge
  //  (6) a2 b0 b2 a0: A and B intersect in two wedges
  //
  // We do not distinguish between 4, 5, and 6.
  // We pay extra attention when some of the edges overlap.  When edges
  // overlap, several of these orderings can be satisfied, and we take
  // the most specific.
  if (a0 == b0 && a2 == b2) return WEDGE_EQUALS;

  if (S2::OrderedCCW(a0, a2, b2, ab1)) {
    // The cases with this vertex ordering are 1, 5, and 6,
    // although case 2 is also possible if a2 == b2.
    if (S2::OrderedCCW(b2, b0, a0, ab1)) return WEDGE_PROPERLY_CONTAINS;

    // We are in case 5 or 6, or case 2 if a2 == b2.
    return (a2 == b2) ? WEDGE_IS_PROPERLY_CONTAINED : WEDGE_PROPERLY_OVERLAPS;
  }

  // We are in case 2, 3, or 4.
  if (S2::OrderedCCW(a0, b0, b2, ab1)) return WEDGE_IS_PROPERLY_CONTAINED;
  return S2::OrderedCCW(a0, b0, a2, ab1) ?
      WEDGE_IS_DISJOINT : WEDGE_PROPERLY_OVERLAPS;
}

int S2EdgeUtil::EdgeCrosser::RobustCrossingInternal(S2Point const* d) {
  // Compute the actual result, and then save the current vertex D as the next
  // vertex C, and save the orientation of the next triangle ACB (which is
  // opposite to the current triangle BDA).
  int result = RobustCrossingInternal2(d);
  c_ = d;
  acb_ = -bda_;
  return result;
}

inline int S2EdgeUtil::EdgeCrosser::RobustCrossingInternal2(S2Point const* d) {
  // RobustCCW is very expensive, so we avoid calling it if at all possible.
  // First eliminate the cases where two vertices are equal.
  if (*a_ == *c_ || *a_ == *d || *b_ == *c_ || *b_ == *d) return 0;

  // At this point, a very common situation is that A,B,C,D are four points on
  // a line such that AB does not overlap CD.  (For example, this happens when
  // a line or curve is sampled finely, or when geometry is constructed by
  // computing the union of S2CellIds.)  Most of the time, we can determine
  // that AB and CD do not intersect by computing the two outward-facing
  // tangents at A and B (parallel to AB) and testing whether AB and CD are on
  // opposite sides of the plane perpendicular to one of these tangents.  This
  // is moderately expensive but still much cheaper than S2::ExpensiveCCW().
  if (!have_tangents_) {
    S2Point norm = S2::RobustCrossProd(*a_, *b_).Normalize();
    a_tangent_ = a_->CrossProd(norm);
    b_tangent_ = norm.CrossProd(*b_);
    have_tangents_ = true;
  }
  // The error in RobustCrossProd() is insignificant.  The maximum error in
  // the call to CrossProd() (i.e., the maximum norm of the error vector) is
  // (0.5 + 1/sqrt(3)) * DBL_EPSILON.  The maximum error in each call to
  // DotProd() below is DBL_EPSILON.  (There is also a small relative error
  // term that is insignificant because we are comparing the result against a
  // constant that is very close to zero.)
  static const double kError = (1.5 + 1/sqrt(3)) * DBL_EPSILON;
  if ((c_->DotProd(a_tangent_) > kError && d->DotProd(a_tangent_) > kError) ||
      (c_->DotProd(b_tangent_) > kError && d->DotProd(b_tangent_) > kError)) {
    return -1;
  }

  // Otherwise it's time to break out the big guns.
  if (acb_ == 0) acb_ = -S2::ExpensiveCCW(*a_, *b_, *c_);
  if (bda_ == 0) bda_ = S2::ExpensiveCCW(*a_, *b_, *d);
  if (bda_ != acb_) return -1;

  Vector3_d c_cross_d = c_->CrossProd(*d);
  int cbd = -S2::RobustCCW(*c_, *d, *b_, c_cross_d);
  if (cbd != acb_) return -1;
  int dac = S2::RobustCCW(*c_, *d, *a_, c_cross_d);
  return (dac == acb_) ? 1 : -1;
}

void S2EdgeUtil::RectBounder::AddPoint(S2Point const* b) {
  DCHECK(S2::IsUnitLength(*b));
  S2LatLng b_latlng(*b);
  if (bound_.is_empty()) {
    bound_.AddPoint(b_latlng);
  } else {
    // First compute the cross product N = A x B robustly.  This is the normal
    // to the great circle through A and B.  We don't use S2::RobustCrossProd()
    // since that method returns an arbitrary vector orthogonal to A if the two
    // vectors are proportional, and we want the zero vector in that case.
    Vector3_d n = (*a_ - *b).CrossProd(*a_ + *b);  // N = 2 * (A x B)

    // The relative error in N gets large as its norm gets very small (i.e.,
    // when the two points are nearly identical or antipodal).  We handle this
    // by choosing a maximum allowable error, and if the error is greater than
    // this we fall back to a different technique.  Since it turns out that
    // the other sources of error add up to at most 1.16 * DBL_EPSILON, and it
    // is desirable to have the total error be a multiple of DBL_EPSILON, we
    // have chosen the maximum error threshold here to be 3.84 * DBL_EPSILON.
    // It is possible to show that the error is less than this when
    //
    //   n.Norm() >= 8 * sqrt(3) / (3.84 - 0.5 - sqrt(3)) * DBL_EPSILON
    //            = 1.91346e-15 (about 8.618 * DBL_EPSILON)
    double n_norm = n.Norm();
    if (n_norm < 1.91346e-15) {
      // A and B are either nearly identical or nearly antipodal (to within
      // 4.309 * DBL_EPSILON, or about 6 nanometers on the earth's surface).
      if (a_->DotProd(*b) < 0) {
        // The two points are nearly antipodal.  The easiest solution is to
        // assume that the edge between A and B could go in any direction
        // around the sphere.
        bound_ = S2LatLngRect::Full();
      } else {
        // The two points are nearly identical (to within 4.309 * DBL_EPSILON).
        // In this case we can just use the bounding rectangle of the points,
        // since after the expansion done by GetBound() this rectangle is
        // guaranteed to include the (lat,lng) values of all points along AB.
        bound_ = bound_.Union(S2LatLngRect::FromPointPair(a_latlng_, b_latlng));
      }
    } else {
      // Compute the longitude range spanned by AB.
      S1Interval lng_ab = S1Interval::FromPointPair(a_latlng_.lng().radians(),
                                                    b_latlng.lng().radians());
      if (lng_ab.GetLength() >= M_PI - 2 * DBL_EPSILON) {
        // The points lie on nearly opposite lines of longitude to within the
        // maximum error of the calculation.  (Note that this test relies on
        // the fact that M_PI is slightly less than the true value of Pi, and
        // that representable values near M_PI are 2 * DBL_EPSILON apart.)
        // The easiest solution is to assume that AB could go on either side
        // of the pole.
        lng_ab = S1Interval::Full();
      }

      // Next we compute the latitude range spanned by the edge AB.  We start
      // with the range spanning the two endpoints of the edge:
      R1Interval lat_ab = R1Interval::FromPointPair(a_latlng_.lat().radians(),
                                                    b_latlng.lat().radians());

      // This is the desired range unless the edge AB crosses the plane
      // through N and the Z-axis (which is where the great circle through A
      // and B attains its minimum and maximum latitudes).  To test whether AB
      // crosses this plane, we compute a vector M perpendicular to this
      // plane and then project A and B onto it.
      Vector3_d m = n.CrossProd(S2Point(0, 0, 1));
      double m_a = m.DotProd(*a_);
      double m_b = m.DotProd(*b);

      // We want to test the signs of "m_a" and "m_b", so we need to bound
      // the error in these calculations.  It is possible to show that the
      // total error is bounded by
      //
      //  (1 + sqrt(3)) * DBL_EPSILON * n_norm + 8 * sqrt(3) * (DBL_EPSILON**2)
      //    = 6.06638e-16 * n_norm + 6.83174e-31

      double m_error = 6.06638e-16 * n_norm + 6.83174e-31;
      if (m_a * m_b < 0 || fabs(m_a) <= m_error || fabs(m_b) <= m_error) {
        // Minimum/maximum latitude *may* occur in the edge interior.
        //
        // The maximum latitude is 90 degrees minus the latitude of N.  We
        // compute this directly using atan2 in order to get maximum accuracy
        // near the poles.
        //
        // Our goal is compute a bound that contains the computed latitudes of
        // all S2Points P that pass the point-in-polygon containment test.
        // There are three sources of error we need to consider:
        //  - the directional error in N (at most 3.84 * DBL_EPSILON)
        //  - converting N to a maximum latitude
        //  - computing the latitude of the test point P
        // The latter two sources of error are at most 0.955 * DBL_EPSILON
        // individually, but it is possible to show by a more complex analysis
        // that together they can add up to at most 1.16 * DBL_EPSILON, for a
        // total error of 5 * DBL_EPSILON.
        //
        // We add 3 * DBL_EPSILON to the bound here, and GetBound() will pad
        // the bound by another 2 * DBL_EPSILON.
        double max_lat = min(
            atan2(sqrt(n[0]*n[0] + n[1]*n[1]), fabs(n[2])) + 3 * DBL_EPSILON,
            M_PI_2);

        // In order to get tight bounds when the two points are close together,
        // we also bound the min/max latitude relative to the latitudes of the
        // endpoints A and B.  First we compute the distance between A and B,
        // and then we compute the maximum change in latitude between any two
        // points along the great circle that are separated by this distance.
        // This gives us a latitude change "budget".  Some of this budget must
        // be spent getting from A to B; the remainder bounds the round-trip
        // distance (in latitude) from A or B to the min or max latitude
        // attained along the edge AB.
        double lat_budget = 2 * asin(0.5 * (*a_ - *b).Norm() * sin(max_lat));
        double max_delta = 0.5*(lat_budget - lat_ab.GetLength()) + DBL_EPSILON;

        // Test whether AB passes through the point of maximum latitude or
        // minimum latitude.  If the dot product(s) are small enough then the
        // result may be ambiguous.
        if (m_a <= m_error && m_b >= -m_error) {
          lat_ab.set_hi(min(max_lat, lat_ab.hi() + max_delta));
        }
        if (m_b <= m_error && m_a >= -m_error) {
          lat_ab.set_lo(max(-max_lat, lat_ab.lo() - max_delta));
        }
      }
      bound_ = bound_.Union(S2LatLngRect(lat_ab, lng_ab));
    }
  }
  a_ = b;
  a_latlng_ = b_latlng;
}

S2LatLngRect S2EdgeUtil::RectBounder::GetBound() const {
  // To save time, we ignore numerical errors in the computed S2LatLngs while
  // accumulating the bounds and then account for them here.
  //
  // S2LatLng(S2Point) has a maximum error of 0.955 * DBL_EPSILON in latitude.
  // In the worst case, we might have rounded "inwards" when computing the
  // bound and "outwards" when computing the latitude of a contained point P,
  // therefore we expand the latitude bounds by 2 * DBL_EPSILON in each
  // direction.  (A more complex analysis shows that 1.5 * DBL_EPSILON is
  // enough, but the expansion amount should be a multiple of DBL_EPSILON in
  // order to avoid rounding errors during the expansion itself.)
  //
  // S2LatLng(S2Point) has a maximum error of DBL_EPSILON in longitude, which
  // is simply the maximum rounding error for results in the range [-Pi, Pi].
  // This is true because the Gnu implementation of atan2() comes from the IBM
  // Accurate Mathematical Library, which implements correct rounding for this
  // instrinsic (i.e., it returns the infinite precision result rounded to the
  // nearest representable value, with ties rounded to even values).  This
  // implies that we don't need to expand the longitude bounds at all, since
  // we only guarantee that the bound contains the *rounded* latitudes of
  // contained points.  The *true* latitudes of contained points may lie up to
  // DBL_EPSILON outside of the returned bound.

  const S2LatLng kExpansion = S2LatLng::FromRadians(2 * DBL_EPSILON, 0);
  return bound_.Expanded(kExpansion).PolarClosure();
}

S2LatLngRect S2EdgeUtil::RectBounder::ExpandForSubregions(
    S2LatLngRect const& bound) {
  // Empty bounds don't need expansion.
  if (bound.is_empty()) return bound;

  // First we need to check whether the bound B contains any nearly-antipodal
  // points (to within 4.309 * DBL_EPSILON).  If so then we need to return
  // S2LatLngRect::Full(), since the subregion might have an edge between two
  // such points, and AddPoint() returns Full() for such edges.  Note that
  // this can happen even if B is not Full(); for example, consider a loop
  // that defines a 10km strip straddling the equator extending from
  // longitudes -100 to +100 degrees.
  //
  // It is easy to check whether B contains any antipodal points, but checking
  // for nearly-antipodal points is trickier.  Essentially we consider the
  // original bound B and its reflection through the origin B', and then test
  // whether the minimum distance between B and B' is less than 4.309 *
  // DBL_EPSILON.

  // "lng_gap" is a lower bound on the longitudinal distance between B and its
  // reflection B'.  (2.5 * DBL_EPSILON is the maximum combined error of the
  // endpoint longitude calculations and the GetLength() call.)
  double lng_gap = max(0.0, M_PI - bound.lng().GetLength() - 2.5 * DBL_EPSILON);

  // "min_abs_lat" is the minimum distance from B to the equator (if zero or
  // negative, then B straddles the equator).
  double min_abs_lat = max(bound.lat().lo(), -bound.lat().hi());

  // "lat_gap1" and "lat_gap2" measure the minimum distance from B to the
  // south and north poles respectively.
  double lat_gap1 = M_PI_2 + bound.lat().lo();
  double lat_gap2 = M_PI_2 - bound.lat().hi();

  if (min_abs_lat >= 0) {
    // The bound B does not straddle the equator.  In this case the minimum
    // distance is between one endpoint of the latitude edge in B closest to
    // the equator and the other endpoint of that edge in B'.  The latitude
    // distance between these two points is 2*min_abs_lat, and the longitude
    // distance is lng_gap.  We could compute the distance exactly using the
    // Haversine formula, but then we would need to bound the errors in that
    // calculation.  Since we only need accuracy when the distance is very
    // small (close to 4.309 * DBL_EPSILON), we substitute the Euclidean
    // distance instead.  This gives us a right triangle XYZ with two edges of
    // length x = 2*min_abs_lat and y ~= lng_gap.  The desired distance is the
    // length of the third edge "z", and we have
    //
    //         z  ~=  sqrt(x^2 + y^2)  >=  (x + y) / sqrt(2)
    //
    // Therefore the region may contain nearly antipodal points only if
    //
    //  2*min_abs_lat + lng_gap  <  sqrt(2) * 4.309 * DBL_EPSILON
    //                           ~= 1.354e-15
    //
    // Note that because the given bound B is conservative, "min_abs_lat" and
    // "lng_gap" are both lower bounds on their true values so we do not need
    // to make any adjustments for their errors.
    if (2 * min_abs_lat + lng_gap < 1.354e-15) {
      return S2LatLngRect::Full();
    }
  } else if (lng_gap >= M_PI_2) {
    // B spans at most Pi/2 in longitude.  The minimum distance is always
    // between one corner of B and the diagonally opposite corner of B'.  We
    // use the same distance approximation that we used above; in this case
    // we have an obtuse triangle XYZ with two edges of length x = lat_gap1
    // and y = lat_gap2, and angle Z >= Pi/2 between them.  We then have
    //
    //         z  >=  sqrt(x^2 + y^2)  >=  (x + y) / sqrt(2)
    //
    // Unlike the case above, "lat_gap1" and "lat_gap2" are not lower bounds
    // (because of the extra addition operation, and because M_PI_2 is not
    // exactly equal to Pi/2); they can exceed their true values by up to
    // 0.75 * DBL_EPSILON.  Putting this all together, the region may
    // contain nearly antipodal points only if
    //
    //   lat_gap1 + lat_gap2  <  (sqrt(2) * 4.309 + 1.5) * DBL_EPSILON
    //                        ~= 1.687e-15
    if (lat_gap1 + lat_gap2 < 1.687e-15) {
      return S2LatLngRect::Full();
    }
  } else {
    // Otherwise we know that (1) the bound straddles the equator and (2) its
    // width in longitude is at least Pi/2.  In this case the minimum
    // distance can occur either between a corner of B and the diagonally
    // opposite corner of B' (as in the case above), or between a corner of B
    // and the opposite longitudinal edge reflected in B'.  It is sufficient
    // to only consider the corner-edge case, since this distance is also a
    // lower bound on the corner-corner distance when that case applies.

    // Consider the spherical triangle XYZ where X is a corner of B with
    // minimum absolute latitude, Y is the closest pole to X, and Z is the
    // point closest to X on the opposite longitudinal edge of B'.  This is a
    // right triangle (Z = Pi/2), and from the spherical law of sines we have
    //
    //     sin(z) / sin(Z)  =  sin(y) / sin(Y)
    //     sin(max_lat_gap) / 1  =  sin(d_min) / sin(lng_gap)
    //     sin(d_min)  =  sin(max_lat_gap) * sin(lng_gap)
    //
    // where "max_lat_gap" = max(lat_gap1, lat_gap2) and "d_min" is the
    // desired minimum distance.  Now using the facts that sin(t) >= (2/Pi)*t
    // for 0 <= t <= Pi/2, that we only need an accurate approximation when
    // at least one of "max_lat_gap" or "lng_gap" is extremely small (in
    // which case sin(t) ~= t), and recalling that "max_lat_gap" has an error
    // of up to 0.75 * DBL_EPSILON, we want to test whether
    //
    //   max_lat_gap * lng_gap  <  (4.309 + 0.75) * (Pi/2) * DBL_EPSILON
    //                          ~= 1.765e-15
    if (max(lat_gap1, lat_gap2) * lng_gap < 1.765e-15) {
      return S2LatLngRect::Full();
    }
  }
  // Next we need to check whether the subregion might contain any edges that
  // span (M_PI - 2 * DBL_EPSILON) radians or more in longitude, since AddPoint
  // sets the longitude bound to Full() in that case.  This corresponds to
  // testing whether (lng_gap <= 0) in "lng_expansion" below.

  // Otherwise, the maximum latitude error in AddPoint is 4.8 * DBL_EPSILON.
  // In the worst case, the errors when computing the latitude bound for a
  // subregion could go in the opposite direction as the errors when computing
  // the bound for the original region, so we need to double this value.
  // (More analysis shows that it's okay to round down to a multiple of
  // DBL_EPSILON.)
  //
  // For longitude, we rely on the fact that atan2 is correctly rounded and
  // therefore no additional bounds expansion is necessary.

  double lat_expansion = 9 * DBL_EPSILON;
  double lng_expansion = (lng_gap <= 0) ? M_PI : 0;
  return bound.Expanded(S2LatLng::FromRadians(lat_expansion,
                                              lng_expansion)).PolarClosure();
}

S2LatLng S2EdgeUtil::RectBounder::MaxErrorForTests() {
  // The maximum error in the latitude calculation is
  //    3.84 * DBL_EPSILON   for the RobustCrossProd calculation
  //    0.96 * DBL_EPSILON   for the Latitude() calculation
  //    5    * DBL_EPSILON   added by AddPoint/GetBound to compensate for error
  //    ------------------
  //    9.80 * DBL_EPSILON   maximum error in result
  //
  // The maximum error in the longitude calculation is DBL_EPSILON.  GetBound
  // does not do any expansion because this isn't necessary in order to
  // bound the *rounded* longitudes of contained points.
  return S2LatLng::FromRadians(10 * DBL_EPSILON, 1 * DBL_EPSILON);
}

S2EdgeUtil::LongitudePruner::LongitudePruner(S1Interval const& interval,
                                             S2Point const& v0)
  : interval_(interval), lng0_(S2LatLng::Longitude(v0).radians()) {
}

// S2PointUVW is used to document that a given S2Point is expressed in the
// (u,v,w) coordinates of some cube face.
typedef S2Point S2PointUVW;

// The three functions below all compare a sum (u + v) to a third value w.
// They are implemented in such a way that they produce an exact result even
// though all calculations are done with ordinary floating-point operations.
// Here are the principles on which these functions are based:
//
// A. If u + v < w in floating-point, then u + v < w in exact arithmetic.
//
// B. If u + v < w in exact arithmetic, then at least one of the following
//    expressions is true in floating-point:
//       u + v < w
//       u < w - v
//       v < w - u
//
//    Proof: By rearranging terms and substituting ">" for "<", we can assume
//    that all values are non-negative.  Now clearly "w" is not the smallest
//    value, so assume WLOG that "u" is the smallest.  We want to show that
//    u < w - v in floating-point.  If v >= w/2, the calculation of w - v is
//    exact since the result is smaller in magnitude than either input value,
//    so the result holds.  Otherwise we have u <= v < w/2 and w - v >= w/2
//    (even in floating point), so the result also holds.

// Return true if u + v == w exactly.
inline static bool SumEquals(double u, double v, double w) {
  return (u + v == w) && (u == w - v) && (v == w - u);
}

// Return true if a given directed line L intersects the cube face F.  The
// line L is defined by its normal N in the (u,v,w) coordinates of F.
inline static bool IntersectsFace(S2PointUVW const& n) {
  // L intersects the [-1,1]x[-1,1] square in (u,v) if and only if the dot
  // products of N with the four corner vertices (-1,-1,1), (1,-1,1), (1,1,1),
  // and (-1,1,1) do not all have the same sign.  This is true exactly when
  // |Nu| + |Nv| >= |Nw|.  The code below evaluates this expression exactly
  // (see comments above).
  double u = fabs(n[0]), v = fabs(n[1]), w = fabs(n[2]);
  // We only need to consider the cases where u or v is the smallest value,
  // since if w is the smallest then both expressions below will have a
  // positive LHS and a negative RHS.
  return (v >= w - u) && (u >= w - v);
}

// Given a directed line L intersecting a cube face F, return true if L
// intersects two opposite edges of F (including the case where L passes
// exactly through a corner vertex of F).  The line L is defined by its
// normal N in the (u,v,w) coordinates of F.
inline static bool IntersectsOppositeEdges(S2PointUVW const& n) {
  // The line L intersects opposite edges of the [-1,1]x[-1,1] (u,v) square if
  // and only exactly two of the corner vertices lie on each side of L.  This
  // is true exactly when ||Nu| - |Nv|| >= |Nw|.  The code below evaluates this
  // expression exactly (see comments above).
  double u = fabs(n[0]), v = fabs(n[1]), w = fabs(n[2]);
  // If w is the smallest, the following line returns an exact result.
  if (fabs(u - v) != w) return fabs(u - v) >= w;
  // Otherwise u - v = w exactly, or w is not the smallest value.  In either
  // case the following line returns the correct result.
  return (u >= v) ? (u - w >= v) : (v - w >= u);
}

// Given cube face F and a directed line L (represented by its CCW normal N in
// the (u,v,w) coordinates of F), compute the axis of the cube face edge where
// L exits the face: return 0 if L exits through the u=-1 or u=+1 edge, and 1
// if L exits through the v=-1 or v=+1 edge.  Either result is acceptable if L
// exits exactly through a corner vertex of the cube face.
static int GetExitAxis(S2PointUVW const& n) {
  DCHECK(IntersectsFace(n));
  if (IntersectsOppositeEdges(n)) {
    // The line passes through through opposite edges of the face.
    // It exits through the v=+1 or v=-1 edge if the u-component of N has a
    // larger absolute magnitude than the v-component.
    return (fabs(n[0]) >= fabs(n[1])) ? 1 : 0;
  } else {
    // The line passes through through two adjacent edges of the face.
    // It exits the v=+1 or v=-1 edge if an even number of the components of N
    // are negative.  We test this using signbit() rather than multiplication
    // to avoid the possibility of underflow.
    DCHECK(n[0] != 0 && n[1] != 0  && n[2] != 0);
    using std::signbit;
    return ((signbit(n[0]) ^ signbit(n[1]) ^ signbit(n[2])) == 0) ? 1 : 0;
  }
}

// Given a cube face F, a directed line L (represented by its CCW normal N in
// the (u,v,w) coordinates of F), and result of GetExitAxis(N), return the
// (u,v) coordinates of the point where L exits the cube face.
static R2Point GetExitPoint(S2PointUVW const& n, int axis) {
  if (axis == 0) {
    double u = (n[1] > 0) ? 1.0 : -1.0;
    return R2Point(u, (-u * n[0] - n[2]) / n[1]);
  } else {
    double v = (n[0] < 0) ? 1.0 : -1.0;
    return R2Point((-v * n[1] - n[2]) / n[0], v);
  }
}

// Given a line segment AB whose origin A has been projected onto a given cube
// face, determine whether it is necessary to project A onto a different face
// instead.  This can happen because the normal of the line AB is not computed
// exactly, so that the line AB (defined as the set of points perpendicular to
// the normal) may not intersect the cube face containing A.  Even if it does
// intersect the face, the "exit point" of the line from that face may be on
// the wrong side of A (i.e., in the direction away from B).  If this happens,
// we reproject A onto the adjacent face where the line AB approaches A most
// closely.  This moves the origin by a small amount, but never more than the
// error tolerances documented in the header file.
static int MoveOriginToValidFace(int face, S2Point const& a,
                                 S2Point const& ab, R2Point* a_uv) {
  // Fast path: if the origin is sufficiently far inside the face, it is
  // always safe to use it.
  double const kMaxSafeUVCoord = 1 - S2EdgeUtil::kFaceClipErrorUVCoord;
  if (max(fabs((*a_uv)[0]), fabs((*a_uv)[1])) <= kMaxSafeUVCoord) {
    return face;
  }
  // Otherwise check whether the normal AB even intersects this face.
  S2PointUVW n = S2::FaceXYZtoUVW(face, ab);
  if (IntersectsFace(n)) {
    // Check whether the point where the line AB exits this face is on the
    // wrong side of A (by more than the acceptable error tolerance).
    S2Point exit = S2::FaceUVtoXYZ(face, GetExitPoint(n, GetExitAxis(n)));
    S2Point a_tangent = ab.Normalize().CrossProd(a);
    if ((exit - a).DotProd(a_tangent) >= -S2EdgeUtil::kFaceClipErrorRadians) {
      return face;  // We can use the given face.
    }
  }
  // Otherwise we reproject A to the nearest adjacent face.  (If line AB does
  // not pass through a given face, it must pass through all adjacent faces.)
  if (fabs((*a_uv)[0]) >= fabs((*a_uv)[1])) {
    face = S2::GetUVWFace(face, 0 /*U axis*/, (*a_uv)[0] > 0);
  } else {
    face = S2::GetUVWFace(face, 1 /*V axis*/, (*a_uv)[1] > 0);
  }
  DCHECK(IntersectsFace(S2::FaceXYZtoUVW(face, ab)));
  S2::ValidFaceXYZtoUV(face, a, a_uv);
  (*a_uv)[0] = max(-1.0, min(1.0, (*a_uv)[0]));
  (*a_uv)[1] = max(-1.0, min(1.0, (*a_uv)[1]));
  return face;
}

// Return the next face that should be visited by GetFaceSegments, given that
// we have just visited "face" and we are following the line AB (represented
// by its normal N in the (u,v,w) coordinates of that face).  The other
// arguments include the point where AB exits "face", the corresponding
// exit axis, and the "target face" containing the destination point B.
static int GetNextFace(int face, R2Point const& exit, int axis,
                       S2PointUVW const& n, int target_face) {
  // We return the face that is adjacent to the exit point along the given
  // axis.  If line AB exits *exactly* through a corner of the face, there are
  // two possible next faces.  If one is the "target face" containing B, then
  // we guarantee that we advance to that face directly.
  //
  // The three conditions below check that (1) AB exits approximately through
  // a corner, (2) the adjacent face along the non-exit axis is the target
  // face, and (3) AB exits *exactly* through the corner.  (The SumEquals()
  // code checks whether the dot product of (u,v,1) and "n" is exactly zero.)
  if (fabs(exit[1 - axis]) == 1 &&
      S2::GetUVWFace(face, 1 - axis, exit[1 - axis] > 0) == target_face &&
      SumEquals(exit[0] * n[0], exit[1] * n[1], -n[2])) {
    return target_face;
  }
  // Otherwise return the face that is adjacent to the exit point in the
  // direction of the exit axis.
  return S2::GetUVWFace(face, axis, exit[axis] > 0);
}

void S2EdgeUtil::GetFaceSegments(S2Point const& a, S2Point const& b,
                                 FaceSegmentVector* segments) {
  DCHECK(S2::IsUnitLength(a));
  DCHECK(S2::IsUnitLength(b));
  segments->clear();

  // Fast path: both endpoints are on the same face.
  FaceSegment segment;
  int a_face = S2::XYZtoFaceUV(a, &segment.a);
  int b_face = S2::XYZtoFaceUV(b, &segment.b);
  if (a_face == b_face) {
    segment.face = a_face;
    segments->push_back(segment);
    return;
  }
  // Starting at A, we follow AB from face to face until we reach the face
  // containing B.  The following code is designed to ensure that we always
  // reach B, even in the presence of numerical errors.
  //
  // First we compute the normal to the plane containing A and B.  This normal
  // becomes the ultimate definition of the line AB; it is used to resolve all
  // questions regarding where exactly the line goes.  Unfortunately due to
  // numerical errors, the line may not quite intersect the faces containing
  // the original endpoints.  We handle this by moving A and/or B slightly if
  // necessary so that they are on faces intersected by the line AB.
  S2Point ab = S2::RobustCrossProd(a, b);
  a_face = MoveOriginToValidFace(a_face, a, ab, &segment.a);
  b_face = MoveOriginToValidFace(b_face, b, -ab, &segment.b);

  // Now we simply follow AB from face to face until we reach B.
  segment.face = a_face;
  R2Point b_saved = segment.b;
  for (int face = a_face; face != b_face; ) {
    // Complete the current segment by finding the point where AB exits the
    // current face.
    S2PointUVW n = S2::FaceXYZtoUVW(face, ab);
    int exit_axis = GetExitAxis(n);
    segment.b = GetExitPoint(n, exit_axis);
    segments->push_back(segment);

    // Compute the next face intersected by AB, and translate the exit point
    // of the current segment into the (u,v) coordinates of the next face.
    // This becomes the first point of the next segment.
    S2Point exit_xyz = S2::FaceUVtoXYZ(face, segment.b);
    face = GetNextFace(face, segment.b, exit_axis, n, b_face);
    S2PointUVW exit_uvw = S2::FaceXYZtoUVW(face, exit_xyz);
    segment.face = face;
    segment.a = R2Point(exit_uvw[0], exit_uvw[1]);
  }
  // Finish the last segment.
  segment.b = b_saved;
  segments->push_back(segment);
}

// This helper function does two things.  First, it clips the line segment AB
// to find the clipped destination B' on a given face.  (The face is specified
// implicitly by expressing *all arguments* in the (u,v,w) coordinates of that
// face.)  Second, it partially computes whether the segment AB intersects
// this face at all.  The actual condition is fairly complicated, but it turns
// out that it can be expressed as a "score" that can be computed
// independently when clipping the two endpoints A and B.  This function
// returns the score for the given endpoint, which is an integer ranging from
// 0 to 3.  If the sum of the two scores is 3 or more, then AB does not
// intersect this face.  See the calling function for the meaning of the
// various parameters.
static int ClipDestination(
    S2PointUVW const& a, S2PointUVW const& b, S2PointUVW const& scaled_n,
    S2PointUVW const& a_tangent, S2PointUVW const& b_tangent, double scale_uv,
    R2Point* uv) {
  DCHECK(IntersectsFace(scaled_n));

  // Optimization: if B is within the safe region of the face, use it.
  double const kMaxSafeUVCoord = 1 - S2EdgeUtil::kFaceClipErrorUVCoord;
  if (b[2] > 0) {
    *uv = R2Point(b[0] / b[2], b[1] / b[2]);
    if (max(fabs((*uv)[0]), fabs((*uv)[1])) <= kMaxSafeUVCoord)
      return 0;
  }
  // Otherwise find the point B' where the line AB exits the face.
  *uv = scale_uv * GetExitPoint(scaled_n, GetExitAxis(scaled_n));
  S2PointUVW p((*uv)[0], (*uv)[1], 1.0);

  // Determine if the exit point B' is contained within the segment.  We do this
  // by computing the dot products with two inward-facing tangent vectors at A
  // and B.  If either dot product is negative, we say that B' is on the "wrong
  // side" of that point.  As the point B' moves around the great circle AB past
  // the segment endpoint B, it is initially on the wrong side of B only; as it
  // moves further it is on the wrong side of both endpoints; and then it is on
  // the wrong side of A only.  If the exit point B' is on the wrong side of
  // either endpoint, we can't use it; instead the segment is clipped at the
  // original endpoint B.
  //
  // We reject the segment if the sum of the scores of the two endpoints is 3
  // or more.  Here is what that rule encodes:
  //  - If B' is on the wrong side of A, then the other clipped endpoint A'
  //    must be in the interior of AB (otherwise AB' would go the wrong way
  //    around the circle).  There is a similar rule for A'.
  //  - If B' is on the wrong side of either endpoint (and therefore we must
  //    use the original endpoint B instead), then it must be possible to
  //    project B onto this face (i.e., its w-coordinate must be positive).
  //    This rule is only necessary to handle certain zero-length edges (A=B).
  int score = 0;
  if ((p - a).DotProd(a_tangent) < 0) {
    score = 2;  // B' is on wrong side of A.
  } else if ((p - b).DotProd(b_tangent) < 0) {
    score = 1;  // B' is on wrong side of B.
  }
  if (score > 0) {  // B' is not in the interior of AB.
    if (b[2] <= 0) {
      score = 3;    // B cannot be projected onto this face.
    } else {
      *uv = R2Point(b[0] / b[2], b[1] / b[2]);
    }
  }
  return score;
}

bool S2EdgeUtil::ClipToPaddedFace(S2Point const& a_xyz, S2Point const& b_xyz,
                                  int face, double padding,
                                  R2Point* a_uv, R2Point* b_uv) {
  DCHECK_GE(padding, 0);
  // Fast path: both endpoints are on the given face.
  if (S2::GetFace(a_xyz) == face && S2::GetFace(b_xyz) == face) {
    S2::ValidFaceXYZtoUV(face, a_xyz, a_uv);
    S2::ValidFaceXYZtoUV(face, b_xyz, b_uv);
    return true;
  }
  // Convert everything into the (u,v,w) coordinates of the given face.  Note
  // that the cross product *must* be computed in the original (x,y,z)
  // coordinate system because RobustCrossProd (unlike the mathematical cross
  // product) can produce different results in different coordinate systems
  // when one argument is a linear multiple of the other, due to the use of
  // symbolic perturbations.
  S2PointUVW n = S2::FaceXYZtoUVW(face, S2::RobustCrossProd(a_xyz, b_xyz));
  S2PointUVW a = S2::FaceXYZtoUVW(face, a_xyz);
  S2PointUVW b = S2::FaceXYZtoUVW(face, b_xyz);

  // Padding is handled by scaling the u- and v-components of the normal.
  // Letting R=1+padding, this means that when we compute the dot product of
  // the normal with a cube face vertex (such as (-1,-1,1)), we will actually
  // compute the dot product with the scaled vertex (-R,-R,1).  This allows
  // methods such as IntersectsFace(), GetExitAxis(), etc, to handle padding
  // with no further modifications.
  double const scale_uv = 1 + padding;
  S2PointUVW scaled_n(scale_uv * n[0], scale_uv * n[1], n[2]);
  if (!IntersectsFace(scaled_n)) return false;

  // TODO(ericv): This is a temporary hack until I rewrite S2::RobustCrossProd;
  // it avoids loss of precision in Normalize() when the vector is so small
  // that it underflows.
  if (max(fabs(n[0]), max(fabs(n[1]), fabs(n[2]))) < ldexp(1, -511)) {
    n *= ldexp(1, 563);
  }  // END OF HACK
  n = n.Normalize();
  S2PointUVW a_tangent = n.CrossProd(a);
  S2PointUVW b_tangent = b.CrossProd(n);
  // As described above, if the sum of the scores from clipping the two
  // endpoints is 3 or more, then the segment does not intersect this face.
  int a_score = ClipDestination(b, a, -scaled_n, b_tangent, a_tangent,
                                scale_uv, a_uv);
  int b_score = ClipDestination(a, b, scaled_n, a_tangent, b_tangent,
                                scale_uv, b_uv);
  return a_score + b_score < 3;
}

bool S2EdgeUtil::IntersectsRect(R2Point const& a, R2Point const& b,
                                R2Rect const& rect) {
  // First check whether the bound of AB intersects "rect".
  R2Rect bound = R2Rect::FromPointPair(a, b);
  if (!rect.Intersects(bound)) return false;

  // Otherwise AB intersects "rect" if and only if all four vertices of "rect"
  // do not lie on the same side of the extended line AB.  We test this by
  // finding the two vertices of "rect" with minimum and maximum projections
  // onto the normal of AB, and computing their dot products with the edge
  // normal.
  R2Point n = (b - a).Ortho();
  int i = (n[0] >= 0) ? 1 : 0;
  int j = (n[1] >= 0) ? 1 : 0;
  double max = n.DotProd(rect.GetVertex(i, j) - a);
  double min = n.DotProd(rect.GetVertex(1-i, 1-j) - a);
  return (max >= 0) && (min <= 0);
}

inline static bool UpdateEndpoint(R1Interval* bound, int end, double value) {
  if (end == 0) {
    if (bound->hi() < value) return false;
    if (bound->lo() < value) bound->set_lo(value);
  } else {
    if (bound->lo() > value) return false;
    if (bound->hi() > value) bound->set_hi(value);
  }
  return true;
}

// Given a line segment from (a0,a1) to (b0,b1) and a bounding interval for
// each axis, clip the segment further if necessary so that "bound0" does not
// extend outside the given interval "clip".  "diag" is a a precomputed helper
// variable that indicates which diagonal of the bounding box is spanned by AB:
// it is 0 if AB has positive slope, and 1 if AB has negative slope.
inline static bool ClipBoundAxis(double a0, double b0, R1Interval* bound0,
                                 double a1, double b1, R1Interval* bound1,
                                 int diag, R1Interval const& clip0) {
  if (bound0->lo() < clip0.lo()) {
    if (bound0->hi() < clip0.lo()) return false;
    (*bound0)[0] = clip0.lo();
    if (!UpdateEndpoint(bound1, diag, S2EdgeUtil::InterpolateDouble(
            clip0.lo(), a0, b0, a1, b1)))
      return false;
  }
  if (bound0->hi() > clip0.hi()) {
    if (bound0->lo() > clip0.hi()) return false;
    (*bound0)[1] = clip0.hi();
    if (!UpdateEndpoint(bound1, 1-diag, S2EdgeUtil::InterpolateDouble(
            clip0.hi(), a0, b0, a1, b1)))
      return false;
  }
  return true;
}

R2Rect S2EdgeUtil::GetClippedEdgeBound(R2Point const& a, R2Point const& b,
                                       R2Rect const& clip) {
  R2Rect bound = R2Rect::FromPointPair(a, b);
  if (ClipEdgeBound(a, b, clip, &bound)) return bound;
  return R2Rect::Empty();
}

bool S2EdgeUtil::ClipEdgeBound(R2Point const& a, R2Point const& b,
                               R2Rect const& clip, R2Rect* bound) {
  // "diag" indicates which diagonal of the bounding box is spanned by AB: it
  // is 0 if AB has positive slope, and 1 if AB has negative slope.  This is
  // used to determine which interval endpoints need to be updated each time
  // the edge is clipped.
  int diag = (a[0] > b[0]) != (a[1] > b[1]);
  return (ClipBoundAxis(a[0], b[0], &(*bound)[0], a[1], b[1], &(*bound)[1],
                        diag, clip[0]) &&
          ClipBoundAxis(a[1], b[1], &(*bound)[1], a[0], b[0], &(*bound)[0],
                        diag, clip[1]));
}

bool S2EdgeUtil::ClipEdge(R2Point const& a, R2Point const& b,
                          R2Rect const& clip,
                          R2Point* a_clipped, R2Point* b_clipped) {
  // Compute the bounding rectangle of AB, clip it, and then extract the new
  // endpoints from the clipped bound.
  R2Rect bound = R2Rect::FromPointPair(a, b);
  if (ClipEdgeBound(a, b, clip, &bound)) {
    int ai = (a[0] > b[0]), aj = (a[1] > b[1]);
    *a_clipped = bound.GetVertex(ai, aj);
    *b_clipped = bound.GetVertex(1-ai, 1-aj);
    return true;
  }
  return false;
}
