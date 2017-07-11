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

#include "s2/s2edge_crossings.h"
#include "s2/s2edge_crossings-internal.h"

#include <cmath>

#include <glog/logging.h>
#include "s2/s1angle.h"
#include "s2/s2edge_crosser.h"
#include "s2/s2pointutil.h"
#include "s2/s2predicates.h"
#include "s2/s2predicates-internal.h"
#include "s2/util/math/exactfloat/exactfloat.h"

namespace S2 {

using internal::GetIntersectionExact;
using internal::IntersectionMethod;
using internal::intersection_method_tally_;

// All error bounds in this file are expressed in terms of the maximum
// rounding error for a floating-point type.  The rounding error is half of
// the numeric_limits<T>::epsilon() value.
static constexpr double DBL_ERR = s2pred::rounding_epsilon<double>();

// kIntersectionError can be set somewhat arbitrarily, because the algorithm
// uses more precision when necessary in order to achieve the specified error.
// The only strict requirement is that kIntersectionError >= 2 * DBL_ERR
// radians.  However, using a larger error tolerance makes the algorithm more
// efficient because it reduces the number of cases where exact arithmetic is
// needed.
S1Angle const kIntersectionError = S1Angle::Radians(8 * DBL_ERR);

S1Angle const kIntersectionMergeRadius = 2 * kIntersectionError;

namespace internal {

S1Angle const kIntersectionExactError = S1Angle::Radians(2 * DBL_ERR);

int* intersection_method_tally_ = nullptr;

char const* GetIntersectionMethodName(IntersectionMethod method) {
  switch (method) {
    case IntersectionMethod::SIMPLE:    return "Simple";
    case IntersectionMethod::SIMPLE_LD: return "Simple_ld";
    case IntersectionMethod::STABLE:    return "Stable";
    case IntersectionMethod::STABLE_LD: return "Stable_ld";
    case IntersectionMethod::EXACT:     return "Exact";
    default:                            return "Unknown";
  }
}

}  // namespace internal

int CrossingSign(S2Point const& a, S2Point const& b,
                 S2Point const& c, S2Point const& d) {
  S2EdgeCrosser crosser(&a, &b, &c);
  return crosser.CrossingSign(&d);
}

bool VertexCrossing(S2Point const& a, S2Point const& b,
                    S2Point const& c, S2Point const& d) {
  // If A == B or C == D there is no intersection.  We need to check this
  // case first in case 3 or more input points are identical.
  if (a == b || c == d) return false;

  // If any other pair of vertices is equal, there is a crossing if and only
  // if OrderedCCW() indicates that the edge AB is further CCW around the
  // shared vertex O (either A or B) than the edge CD, starting from an
  // arbitrary fixed reference point.
  if (a == d) return s2pred::OrderedCCW(S2::Ortho(a), c, b, a);
  if (b == c) return s2pred::OrderedCCW(S2::Ortho(b), d, a, b);
  if (a == c) return s2pred::OrderedCCW(S2::Ortho(a), d, b, a);
  if (b == d) return s2pred::OrderedCCW(S2::Ortho(b), c, a, b);

  LOG(DFATAL) << "VertexCrossing called with 4 distinct vertices";
  return false;
}

bool EdgeOrVertexCrossing(S2Point const& a, S2Point const& b,
                          S2Point const& c, S2Point const& d) {
  int crossing = CrossingSign(a, b, c, d);
  if (crossing < 0) return false;
  if (crossing > 0) return true;
  return VertexCrossing(a, b, c, d);
}

using Vector3_ld = Vector3<long double>;
using Vector3_xf = Vector3<ExactFloat>;

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
// to within an error of at most kIntersectionError by this function, then set
// "result" to the intersection point and return true.
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
  // the length of the cross product is at least (16 * sqrt(3) + 24) * DBL_ERR
  // (about 6e-15), then the directional error is at most 5 * T_ERR (about
  // 3e-19 when T == "long double").  (DBL_ERR appears in the first formula
  // because the inputs are assumed to be normalized in double precision
  // rather than in the given type T.)
  //
  // The third cross product is different because its inputs already have some
  // error.  Letting "result_len" be the length of the cross product, it can
  // be shown that the error is at most
  //
  //   (2 + 2 * sqrt(3) + 12 / result_len) * T_ERR
  //
  // We want this error to be at most kIntersectionError, which is true as
  // long as "result_len" is at least kMinResultLen defined below.

  constexpr T T_ERR = s2pred::rounding_epsilon<T>();
  static const T kMinNormalLength = (16 * sqrt(3) + 24) * DBL_ERR;
  static const T kMinResultLen =
      12 / (kIntersectionError.radians() / T_ERR - (2 + 2 * sqrt(3)));

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
  // (DBL_ERR appears because the input points are assumed to be normalized in
  // double precision rather than in the given type T.)
  //
  // For reference, the bounds that went into this calculation are:
  // ||N'-N|| <= ((1 + 2 * sqrt(3))||N|| + 32 * sqrt(3) * DBL_ERR) * T_ERR
  // |(A.B)'-(A.B)| <= (1.5 * (A.B) + 1.5 * ||A|| * ||B||) * T_ERR
  // ||(X-Y)'-(X-Y)|| <= ||X-Y|| * T_ERR
  constexpr T T_ERR = s2pred::rounding_epsilon<T>();
  *error = (((3.5 + 2 * sqrt(3)) * a_norm_len + 32 * sqrt(3) * DBL_ERR)
            * dist + 1.5 * std::abs(result)) * T_ERR;
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
  constexpr T T_ERR = s2pred::rounding_epsilon<T>();
  T error = b_len * std::abs(b0_dist * b1_error - b1_dist * b0_error) /
      (dist_sum - error_sum) + 2 * T_ERR * dist_sum;

  // Finally we normalize the result, compute the corresponding error, and
  // check whether the total error is acceptable.
  T x_len = x.Norm();
  T const kMaxError = kIntersectionError.radians();
  if (error > (kMaxError - T_ERR) * x_len) {
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
  if (*pa0 >= *pa1) std::swap(pa0, pa1);
  if (*pb0 >= *pb1) std::swap(pb0, pb1);
  return *pa0 < *pb0 || (*pa0 == *pb0 && *pb0 < *pb1);
}

// If the intersection point of the edges (a0,a1) and (b0,b1) can be computed
// to within an error of at most kIntersectionError by this function, then set
// "result" to the intersection point and return true.
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

static S2Point S2PointFromExact(Vector3_xf const& x) {
  // TODO(ericv): In theory, this function may return S2Point(0, 0, 0) even
  // though "x" is non-zero.  This happens when all components of "x" have
  // absolute value less than about 1e-154, since in that case x.Norm2() is
  // zero in double precision.  This could be fixed by scaling "x" by an
  // appropriate power of 2 before the conversion.
  return S2Point(x[0].ToDouble(), x[1].ToDouble(), x[2].ToDouble()).Normalize();
}

namespace internal {

// Compute the intersection point of (a0, a1) and (b0, b1) using exact
// arithmetic.  Note that the result is not exact because it is rounded to
// double precision.  Also, the intersection point is not guaranteed to have
// the correct sign (i.e., the return value may need to be negated).
S2Point GetIntersectionExact(S2Point const& a0, S2Point const& a1,
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
  // directional error of up to 2 * DBL_ERR.  (ToDouble() and Normalize() each
  // contribute up to DBL_ERR of directional error.)
  S2Point x = S2PointFromExact(x_xf);

  if (x == S2Point(0, 0, 0)) {
    // The two edges are exactly collinear, but we still consider them to be
    // "crossing" because of simulation of simplicity.  Out of the four
    // endpoints, exactly two lie in the interior of the other edge.  Of
    // those two we return the one that is lexicographically smallest.
    x = S2Point(10, 10, 10);  // Greater than any valid S2Point
    S2Point a_norm = S2PointFromExact(a_norm_xf);
    S2Point b_norm = S2PointFromExact(b_norm_xf);
    if (s2pred::OrderedCCW(b0, a0, b1, b_norm) && a0 < x) x = a0;
    if (s2pred::OrderedCCW(b0, a1, b1, b_norm) && a1 < x) x = a1;
    if (s2pred::OrderedCCW(a0, b0, a1, a_norm) && b0 < x) x = b0;
    if (s2pred::OrderedCCW(a0, b1, a1, a_norm) && b1 < x) x = b1;
  }
  DCHECK(S2::IsUnitLength(x));
  return x;
}

}  // namespace internal

// Given three points "a", "x", "b", returns true if these three points occur
// in the given order along the edge (a,b) to within the given tolerance.
// More precisely, either "x" must be within "tolerance" of "a" or "b", or
// when "x" is projected onto the great circle through "a" and "b" it must lie
// along the edge (a,b) (i.e., the shortest path from "a" to "b").
static bool ApproximatelyOrdered(S2Point const& a, S2Point const& x,
                                 S2Point const& b, double tolerance) {
  if ((x - a).Norm2() <= tolerance * tolerance) return true;
  if ((x - b).Norm2() <= tolerance * tolerance) return true;
  return s2pred::OrderedCCW(a, x, b, S2::RobustCrossProd(a, b).Normalize());
}

S2Point GetIntersection(S2Point const& a0, S2Point const& a1,
                        S2Point const& b0, S2Point const& b1) {
  DCHECK_GT(CrossingSign(a0, a1, b0, b1), 0);

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
  static bool const kHasLongDouble = (s2pred::rounding_epsilon<long double>() <
                                      s2pred::rounding_epsilon<double>());
  S2Point result;
  IntersectionMethod method;
  if (kUseSimpleMethod && GetIntersectionSimple(a0, a1, b0, b1, &result)) {
    method = IntersectionMethod::SIMPLE;
  } else if (kUseSimpleMethod && kHasLongDouble &&
             GetIntersectionSimpleLD(a0, a1, b0, b1, &result)) {
    method = IntersectionMethod::SIMPLE_LD;
  } else if (GetIntersectionStable(a0, a1, b0, b1, &result)) {
    method = IntersectionMethod::STABLE;
  } else if (kHasLongDouble &&
             GetIntersectionStableLD(a0, a1, b0, b1, &result)) {
    method = IntersectionMethod::STABLE_LD;
  } else {
    result = GetIntersectionExact(a0, a1, b0, b1);
    method = IntersectionMethod::EXACT;
  }
  if (intersection_method_tally_) {
    ++intersection_method_tally_[static_cast<int>(method)];
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

bool SimpleCrossing(S2Point const& a, S2Point const& b,
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

}  // namespace S2
