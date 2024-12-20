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

#include "s2/s1chord_angle.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <ostream>

#include "absl/log/absl_check.h"
#include "s2/s1angle.h"

using std::max;
using std::min;
// Android with gnustl has ::nextafter but not std::nextafter.
// https://github.com/android-ndk/ndk/issues/82
// Check for gnustl with _GLIBCXX_CMATH, which is its cmath include
// guard.
#if !defined(__ANDROID__) || !defined(_GLIBCXX_CMATH)
using std::nextafter;
#endif

static constexpr double kMaxLength2 = 4.0;

S1Angle S1ChordAngle::ToAngle() const {
  if (is_negative()) return S1Angle::Radians(-1);
  if (is_infinity()) return S1Angle::Infinity();
  return S1Angle::Radians(2 * asin(0.5 * sqrt(length2_)));
}

S1ChordAngle S1ChordAngle::Successor() const {
  if (length2_ >= kMaxLength2) return Infinity();
  if (length2_ < 0.0) return Zero();
  return S1ChordAngle(nextafter(length2_, 10.0));
}

S1ChordAngle S1ChordAngle::Predecessor() const {
  if (length2_ <= 0.0) return Negative();
  if (length2_ > kMaxLength2) return Straight();
  return S1ChordAngle(nextafter(length2_, -10.0));
}

S1ChordAngle S1ChordAngle::PlusError(double error) const {
  // If angle is Negative() or Infinity(), don't change it.
  // Otherwise clamp it to the valid range.
  if (is_special()) return *this;
  return S1ChordAngle(max(0.0, min(kMaxLength2, length2_ + error)));
}

double S1ChordAngle::GetS2PointConstructorMaxError() const {
  // There is a relative error of 2.5 * DBL_EPSILON when computing the squared
  // distance, plus a relative error of 2 * DBL_EPSILON and an absolute error
  // of (16 * DBL_EPSILON**2) because the lengths of the input points may
  // differ from 1 by up to (2 * DBL_EPSILON) each.  (This is the maximum
  // length error in S2Point::Normalize.)
  return 4.5 * DBL_EPSILON * length2_ + 16 * DBL_EPSILON * DBL_EPSILON;
}

double S1ChordAngle::GetS1AngleConstructorMaxError() const {
  // Assuming that an accurate math library is being used, the sin() call and
  // the multiply each have a relative error of 0.5 * DBL_EPSILON.  However
  // the sin() error is squared.
  return 1.5 * DBL_EPSILON * length2_;
}

S1ChordAngle operator+(S1ChordAngle a, S1ChordAngle b) {
  // Error Analysis:
  //
  //   u is the unit round-off, equal to ε/2 = 2^-53
  //   x and y below are both computed with (1+u)² error.  So we have
  //     length2 = x + y + 2*sqrt(x * y)
  //     length2 = ((x + y)(1+u)³ + 2*sqrt((x * y)((1+u)³))(1+u))(1+u)
  //
  //   Taking (1+u)^1.5 out of the square root and rounding (1+u)^2.5 to (1+u)³:
  //     length2 = (x + y + 2*sqrt(x * y))(1+u)⁴
  //
  //   Bounding (1+u)⁴ as 4*1.01u = 4.04*ε/2 = 2.02ε, total relative error is:
  //     length2()*2.02ε.

  // Note that this method is much more efficient than converting the chord
  // angles to S1Angles and adding those.  It requires only one square root
  // plus a few additions and multiplications.
  ABSL_DCHECK(!a.is_special()) << a;
  ABSL_DCHECK(!b.is_special()) << b;

  // Optimization for the common case where "b" is an error tolerance
  // parameter that happens to be set to zero.
  double a2 = a.length2(), b2 = b.length2();
  if (b2 == 0) return a;

  // Clamp the angle sum to at most 180 degrees.
  if (a2 + b2 >= kMaxLength2) return S1ChordAngle::Straight();

  // Let "a" and "b" be the (non-squared) chord lengths, and let c = a+b.
  // Let A, B, and C be the corresponding half-angles (a = 2*sin(A), etc).
  // Then the formula below can be derived from c = 2 * sin(A+B) and the
  // relationships   sin(A+B) = sin(A)*cos(B) + sin(B)*cos(A)
  //                 cos(X) = sqrt(1 - sin^2(X)) .

  double x = a2 * (1 - 0.25 * b2);  // is_valid() => non-negative
  double y = b2 * (1 - 0.25 * a2);  // is_valid() => non-negative
  return S1ChordAngle(min(kMaxLength2, x + y + 2 * sqrt(x * y)));
}

S1ChordAngle operator-(S1ChordAngle a, S1ChordAngle b) {
  // See comments in operator+().
  ABSL_DCHECK(!a.is_special()) << a;
  ABSL_DCHECK(!b.is_special()) << b;
  double a2 = a.length2(), b2 = b.length2();
  if (b2 == 0) return a;
  if (a2 <= b2) return S1ChordAngle::Zero();
  double x = a2 * (1 - 0.25 * b2);
  double y = b2 * (1 - 0.25 * a2);

  // The calculation below is formulated differently (with two square roots
  // rather than one) to avoid excessive cancellation error when two nearly
  // equal values are subtracted.
  double c = max(0.0, sqrt(x) - sqrt(y));
  return S1ChordAngle(c * c);
}

double sin2(S1ChordAngle a) {
  ABSL_DCHECK(!a.is_special());
  // Let "a" be the (non-squared) chord length, and let A be the corresponding
  // half-angle (a = 2*sin(A)).  The formula below can be derived from:
  //   sin(2*A) = 2 * sin(A) * cos(A)
  //   cos^2(A) = 1 - sin^2(A)
  // This is much faster than converting to an angle and computing its sine.
  return a.length2() * (1 - 0.25 * a.length2());
}

double sin(S1ChordAngle a) {
  return sqrt(sin2(a));
}

double cos(S1ChordAngle a) {
  // cos(2*A) = cos^2(A) - sin^2(A) = 1 - 2*sin^2(A)
  ABSL_DCHECK(!a.is_special());
  return 1 - 0.5 * a.length2();
}

double tan(S1ChordAngle a) {
  return sin(a) / cos(a);
}

std::ostream& operator<<(std::ostream& os, S1ChordAngle a) {
  return os << a.ToAngle();
}
