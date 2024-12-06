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

#ifndef S2_S1CHORD_ANGLE_H_
#define S2_S1CHORD_ANGLE_H_

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <limits>
#include <ostream>
#include <type_traits>

#include "absl/log/absl_check.h"
#include "s2/_fp_contract_off.h"  // IWYU pragma: keep
#include "s2/s1angle.h"
#include "s2/s2point.h"
#include "s2/s2pointutil.h"

// S1ChordAngle represents the angle subtended by a chord (i.e., the straight
// line segment connecting two points on the sphere).  Its representation
// makes it very efficient for computing and comparing distances, but unlike
// S1Angle it is only capable of representing angles between 0 and Pi radians.
// S1ChordAngle is intended for applications where many angles need to be
// computed and compared, otherwise it is simpler to use S1Angle.
//
// S1ChordAngle also loses some accuracy as the angle approaches Pi radians.
// There are several different ways to measure this error, including the
// representational error (i.e., how accurately S1ChordAngle can represent
// angles near Pi radians), the conversion error (i.e., how much precision is
// lost when an S1Angle is converted to an S1ChordAngle), and the measurement
// error (i.e., how accurate the S1ChordAngle(a, b) constructor is when the
// points A and B are separated by angles close to Pi radians).  All of these
// errors differ by a small constant factor.
//
// For the measurement error (which is the largest of these errors and also
// the most important in practice), let the angle between A and B be (Pi - x)
// radians, i.e. A and B are within "x" radians of being antipodal.  The
// corresponding chord length is
//
//    r = 2 * sin((Pi - x) / 2) = 2 * cos(x / 2) .
//
// For values of x not close to Pi the relative error in the squared chord
// length is at most 4.5 * DBL_EPSILON (see GetS2PointConstructorMaxError).
// The relative error in "r" is thus at most 2.25 * DBL_EPSILON ~= 5e-16.  To
// convert this error into an equivalent angle, we have
//
//    |dr / dx| = sin(x / 2)
//
// and therefore
//
//    |dx| = dr / sin(x / 2)
//         = 5e-16 * (2 * cos(x / 2)) / sin(x / 2)
//         = 1e-15 / tan(x / 2)
//
// The maximum error is attained when
//
//    x  = |dx|
//       = 1e-15 / tan(x / 2)
//      ~= 1e-15 / (x / 2)
//      ~= sqrt(2e-15)
//
// In summary, the measurement error for an angle (Pi - x) is at most
//
//    dx  = min(1e-15 / tan(x / 2), sqrt(2e-15))
//      (~= min(2e-15 / x, sqrt(2e-15)) when x is small).
//
// On the Earth's surface (assuming a radius of 6371km), this corresponds to
// the following worst-case measurement errors:
//
//     Accuracy:             Unless antipodal to within:
//     ---------             ---------------------------
//     6.4 nanometers        10,000 km (90 degrees)
//     1 micrometer          81.2 kilometers
//     1 millimeter          81.2 meters
//     1 centimeter          8.12 meters
//     28.5 centimeters      28.5 centimeters
//
// The representational and conversion errors referred to earlier are somewhat
// smaller than this.  For example, maximum distance between adjacent
// representable S1ChordAngle values is only 13.5 cm rather than 28.5 cm.  To
// see this, observe that the closest representable value to r^2 = 4 is
// r^2 =  4 * (1 - DBL_EPSILON / 2).  Thus r = 2 * (1 - DBL_EPSILON / 4) and
// the angle between these two representable values is
//
//    x  = 2 * acos(r / 2)
//       = 2 * acos(1 - DBL_EPSILON / 4)
//      ~= 2 * asin(sqrt(DBL_EPSILON / 2)
//      ~= sqrt(2 * DBL_EPSILON)
//      ~= 2.1e-8
//
// which is 13.5 cm on the Earth's surface.
//
// The worst case rounding error occurs when the value halfway between these
// two representable values is rounded up to 4.  This halfway value is
// r^2 = (4 * (1 - DBL_EPSILON / 4)), thus r = 2 * (1 - DBL_EPSILON / 8) and
// the worst case rounding error is
//
//    x  = 2 * acos(r / 2)
//       = 2 * acos(1 - DBL_EPSILON / 8)
//      ~= 2 * asin(sqrt(DBL_EPSILON / 4)
//      ~= sqrt(DBL_EPSILON)
//      ~= 1.5e-8
//
// which is 9.5 cm on the Earth's surface.
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.
class S1ChordAngle {
 public:
  // Maximum relative error when summing two S1ChordAngles together.
  // Absolute error is length2() of the summed S1ChordAngle() times this value.
  // See the definition of operator+() for the error analysis.
  static constexpr double kRelativeSumError = 2.02 * DBL_EPSILON;

  // The default constructor yields a zero angle.  This is useful for STL
  // containers and class methods with output arguments.
  constexpr S1ChordAngle() = default;

  // Construct the S1ChordAngle corresponding to the distance between the two
  // given points.  The points must be unit length.
  S1ChordAngle(const S2Point& x, const S2Point& y);

  // Return the zero chord angle.
  static constexpr S1ChordAngle Zero();

  // Return a chord angle of 90 degrees (a "right angle").
  static constexpr S1ChordAngle Right();

  // Return a chord angle of 180 degrees (a "straight angle").  This is the
  // maximum finite chord angle.
  static constexpr S1ChordAngle Straight();

  // Return a chord angle larger than any finite chord angle.  The only valid
  // operations on Infinity() are comparisons, S1Angle conversions, and
  // Successor() / Predecessor().
  static constexpr S1ChordAngle Infinity();

  // Return a chord angle smaller than Zero().  The only valid operations on
  // Negative() are comparisons, S1Angle conversions, and Successor() /
  // Predecessor().
  static constexpr S1ChordAngle Negative();

  // Conversion from an S1Angle.  Angles outside the range [0, Pi] are handled
  // as follows: Infinity() is mapped to Infinity(), negative angles are
  // mapped to Negative(), and finite angles larger than Pi are mapped to
  // Straight().
  //
  // Note that this operation is relatively expensive and should be avoided.
  // To use S1ChordAngle effectively, you should structure your code so that
  // input arguments are converted to S1ChordAngles at the beginning of your
  // algorithm, and results are converted back to S1Angles only at the end.
  explicit S1ChordAngle(S1Angle angle);

  // Convenience methods implemented by converting from an S1Angle.
  static S1ChordAngle Radians(double radians);
  static S1ChordAngle Degrees(double degrees);
  static S1ChordAngle E5(int32_t e5);
  static S1ChordAngle E6(int32_t e6);
  static S1ChordAngle E7(int32_t e7);

  // Construct an S1ChordAngle that is an upper bound on the given S1Angle,
  // i.e. such that FastUpperBoundFrom(x).ToAngle() >= x.  Unlike the S1Angle
  // constructor above, this method is very fast, and the bound is accurate to
  // within 1% for distances up to about 3100km on the Earth's surface.
  static constexpr S1ChordAngle FastUpperBoundFrom(S1Angle angle);

  // Construct an S1ChordAngle from the squared chord length.  Note that the
  // argument is automatically clamped to a maximum of 4.0 to handle possible
  // roundoff errors.  The argument must be non-negative.
  static constexpr S1ChordAngle FromLength2(double length2);

  // Converts to an S1Angle.  Can be used just like an S1Angle constructor:
  //
  //   S1ChordAngle x = ...;
  //   return S1Angle(x);
  //
  // Infinity() is converted to S1Angle::Infinity(), and Negative() is
  // converted to an unspecified negative S1Angle.
  //
  // Note that the conversion uses trigonometric functions and therefore
  // should be avoided in inner loops.
  explicit operator S1Angle() const;

  // Converts to an S1Angle (equivalent to the operator above).
  S1Angle ToAngle() const;

  // Convenience methods implemented by calling ToAngle() first.  Note that
  // because of the S1Angle conversion these methods are relatively expensive
  // (despite their lowercase names), so the results should be cached if they
  // are needed inside loops.
  double radians() const;
  double degrees() const;
  int32_t e5() const;
  int32_t e6() const;
  int32_t e7() const;

  // All operators and functions are declared here so that we can put them all
  // in one place.  (The compound assignment operators must be put here.)

  // Comparison operators.
  friend bool operator==(S1ChordAngle x, S1ChordAngle y);
  friend bool operator!=(S1ChordAngle x, S1ChordAngle y);
  friend bool operator<(S1ChordAngle x, S1ChordAngle y);
  friend bool operator>(S1ChordAngle x, S1ChordAngle y);
  friend bool operator<=(S1ChordAngle x, S1ChordAngle y);
  friend bool operator>=(S1ChordAngle x, S1ChordAngle y);

  // Comparison predicates.
  constexpr bool is_zero() const;
  constexpr bool is_negative() const;
  constexpr bool is_infinity() const;
  constexpr bool is_special() const;  // Negative or infinity.

  // Only addition and subtraction of S1ChordAngles is supported.  These
  // methods add or subtract the corresponding S1Angles, and clamp the result
  // to the range [0, Pi].  Both arguments must be non-negative and
  // non-infinite.
  //
  // REQUIRES: !a.is_special() && !b.is_special()
  friend S1ChordAngle operator+(S1ChordAngle a, S1ChordAngle b);
  friend S1ChordAngle operator-(S1ChordAngle a, S1ChordAngle b);
  S1ChordAngle& operator+=(S1ChordAngle a);
  S1ChordAngle& operator-=(S1ChordAngle a);

  // Trigonmetric functions.  It is more accurate and efficient to call these
  // rather than first converting to an S1Angle.
  friend double sin(S1ChordAngle a);
  friend double cos(S1ChordAngle a);
  friend double tan(S1ChordAngle a);

  // Returns sin(a)**2, but computed more efficiently.
  friend double sin2(S1ChordAngle a);

  // The squared length of the chord.  (Most clients will not need this.)
  double length2() const { return length2_; }

  // Returns the smallest representable S1ChordAngle larger than this object.
  // This can be used to convert a "<" comparison to a "<=" comparison.  For
  // example:
  //
  //   S2ClosestEdgeQuery query(...);
  //   S1ChordAngle limit = ...;
  //   if (query.IsDistanceLess(target, limit.Successor())) {
  //     // Distance to "target" is less than or equal to "limit".
  //   }
  //
  // Note the following special cases:
  //   Negative().Successor() == Zero()
  //   Straight().Successor() == Infinity()
  //   Infinity().Successor() == Infinity()
  S1ChordAngle Successor() const;

  // Like Successor(), but returns the largest representable S1ChordAngle less
  // than this object.
  //
  // Note the following special cases:
  //   Infinity().Predecessor() == Straight()
  //   Zero().Predecessor() == Negative()
  //   Negative().Predecessor() == Negative()
  S1ChordAngle Predecessor() const;

  // Returns a new S1ChordAngle that has been adjusted by the given error
  // bound (which can be positive or negative).  "error" should be the value
  // returned by one of the error bound methods below.  For example:
  //    S1ChordAngle a(x, y);
  //    S1ChordAngle a1 = a.PlusError(a.GetS2PointConstructorMaxError());
  S1ChordAngle PlusError(double error) const;

  // Return the maximum error in length2() for the S1ChordAngle(x, y)
  // constructor, assuming that "x" and "y" are normalized to within the
  // bounds guaranteed by S2Point::Normalize().  (The error is defined with
  // respect to the true distance after the points are projected to lie
  // exactly on the sphere.)
  double GetS2PointConstructorMaxError() const;

  // Return the maximum error in length2() for the S1Angle constructor.
  double GetS1AngleConstructorMaxError() const;

  // Return true if the internal representation is valid.  Negative() and
  // Infinity() are both considered valid.
  constexpr bool is_valid() const;

  // When S1ChordAngle is used as a key in one of the absl::btree container
  // types, indicate that linear rather than binary search should be used.
  // This is much faster when the comparison function is cheap.
  typedef std::true_type absl_btree_prefer_linear_node_search;

 private:
  // S1ChordAngles are represented by the squared chord length, which can
  // range from 0 to 4.  Infinity() uses an infinite squared length.
  explicit constexpr S1ChordAngle(double length2) : length2_(length2) {
    ABSL_DCHECK(is_valid());
  }
  double length2_ = 0.0;
};


//////////////////   Implementation details follow   ////////////////////

inline S1ChordAngle::S1ChordAngle(S1Angle angle) {
  if (angle.radians() < 0) {
    *this = Negative();
  } else if (angle == S1Angle::Infinity()) {
    *this = Infinity();
  } else {
    // The chord length is 2 * sin(angle / 2).
    double length = 2 * sin(0.5 * std::min(M_PI, angle.radians()));
    length2_ = length * length;
  }
  ABSL_DCHECK(is_valid());
}

inline S1ChordAngle::S1ChordAngle(const S2Point& x, const S2Point& y)
    : length2_(std::min(4.0, (x - y).Norm2())) {
  ABSL_DCHECK(S2::IsUnitLength(x));
  ABSL_DCHECK(S2::IsUnitLength(y));
  // The squared distance may slightly exceed 4.0 due to roundoff errors.
  // The maximum error in the result is 2 * DBL_EPSILON * length2_.
  ABSL_DCHECK(is_valid());
}

inline constexpr bool S1ChordAngle::is_valid() const {
  return (length2_ >= 0 && length2_ <= 4.0) || is_special();
}

inline constexpr S1ChordAngle S1ChordAngle::FromLength2(double length2) {
  return S1ChordAngle(std::min(4.0, length2));
}

inline constexpr S1ChordAngle S1ChordAngle::Zero() {
  return S1ChordAngle(0);
}

inline constexpr S1ChordAngle S1ChordAngle::Right() {
  return S1ChordAngle(2);
}

inline constexpr S1ChordAngle S1ChordAngle::Straight() {
  return S1ChordAngle(4);
}

inline constexpr S1ChordAngle S1ChordAngle::Infinity() {
  return S1ChordAngle(std::numeric_limits<double>::infinity());
}

inline constexpr S1ChordAngle S1ChordAngle::Negative() {
  return S1ChordAngle(-1);
}

inline S1ChordAngle S1ChordAngle::Radians(double radians) {
  return S1ChordAngle(S1Angle::Radians(radians));
}

inline S1ChordAngle S1ChordAngle::Degrees(double degrees) {
  return S1ChordAngle(S1Angle::Degrees(degrees));
}

inline S1ChordAngle S1ChordAngle::E5(int32_t e5) {
  return S1ChordAngle(S1Angle::E5(e5));
}

inline S1ChordAngle S1ChordAngle::E6(int32_t e6) {
  return S1ChordAngle(S1Angle::E6(e6));
}

inline S1ChordAngle S1ChordAngle::E7(int32_t e7) {
  return S1ChordAngle(S1Angle::E7(e7));
}

inline constexpr S1ChordAngle S1ChordAngle::FastUpperBoundFrom(S1Angle angle) {
  // This method uses the distance along the surface of the sphere as an upper
  // bound on the distance through the sphere's interior.
  return S1ChordAngle::FromLength2(angle.radians() * angle.radians());
}

inline S1ChordAngle::operator S1Angle() const {
  return ToAngle();
}

inline double S1ChordAngle::radians() const {
  return ToAngle().radians();
}

inline double S1ChordAngle::degrees() const {
  return ToAngle().degrees();
}

inline int32_t S1ChordAngle::e5() const { return ToAngle().e5(); }

inline int32_t S1ChordAngle::e6() const { return ToAngle().e6(); }

inline int32_t S1ChordAngle::e7() const { return ToAngle().e7(); }

inline constexpr bool S1ChordAngle::is_zero() const { return length2_ == 0; }

inline constexpr bool S1ChordAngle::is_negative() const {
  // TODO(ericv): Consider stricter check here -- only allow Negative().
  return length2_ < 0;
}

inline constexpr bool S1ChordAngle::is_infinity() const {
  return length2_ == std::numeric_limits<double>::infinity();
}

inline constexpr bool S1ChordAngle::is_special() const {
  return is_negative() || is_infinity();
}

inline bool operator==(S1ChordAngle x, S1ChordAngle y) {
  return x.length2() == y.length2();
}

inline bool operator!=(S1ChordAngle x, S1ChordAngle y) {
  return x.length2() != y.length2();
}

inline bool operator<(S1ChordAngle x, S1ChordAngle y) {
  return x.length2() < y.length2();
}

inline bool operator>(S1ChordAngle x, S1ChordAngle y) {
  return x.length2() > y.length2();
}

inline bool operator<=(S1ChordAngle x, S1ChordAngle y) {
  return x.length2() <= y.length2();
}

inline bool operator>=(S1ChordAngle x, S1ChordAngle y) {
  return x.length2() >= y.length2();
}

inline S1ChordAngle& S1ChordAngle::operator+=(S1ChordAngle a) {
  return (*this = *this + a);
}

inline S1ChordAngle& S1ChordAngle::operator-=(S1ChordAngle a) {
  return (*this = *this - a);
}

// Outputs the chord angle as the equivalent S1Angle.
std::ostream& operator<<(std::ostream& os, S1ChordAngle a);

#endif  // S2_S1CHORD_ANGLE_H_
