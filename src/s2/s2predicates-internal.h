// Copyright 2016 Google Inc. All Rights Reserved.
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
//
// The following functions are not part of the public API.  Currently they are
// only used internally for testing purposes.

#ifndef S2_S2PREDICATES_INTERNAL_H_
#define S2_S2PREDICATES_INTERNAL_H_

#include <limits>

#include "s2/base/casts.h"
#include "s2/s1chordangle.h"
#include "s2/s2predicates.h"
#include "s2/util/math/exactfloat/exactfloat.h"
#include "s2/util/math/vector.h"

namespace s2pred {

// Returns 2 ** (-digits).  This could be implemented using "ldexp" except
// that std::ldexp is not constexpr in C++11.
constexpr double epsilon_for_digits(int digits) {
  return (digits < 64 ? 1.0 / (1ull << digits) :
          epsilon_for_digits(digits - 63) / (1ull << 63));
}

// Returns the maximum rounding error for arithmetic operations in type T.
// We could simply return 0.5 * numeric_limits<T>::epsilon(), except that some
// platforms implement "long double" using "double-double" arithmetic, and for
// those platforms we need to compute the rounding error manually based on
// numeric_limits<T>::digits (the number of bits of mantissa precision).
template <typename T> constexpr T rounding_epsilon() {
  return epsilon_for_digits(std::numeric_limits<T>::digits);
}

using Vector3_ld = Vector3<long double>;
using Vector3_xf = Vector3<ExactFloat>;

inline static Vector3_ld ToLD(S2Point const& x) {
  return Vector3_ld::Cast(x);
}

inline static long double ToLD(double x) {
  return implicit_cast<long double>(x);
}

inline static Vector3_xf ToExact(S2Point const& x) {
  return Vector3_xf::Cast(x);
}

int StableSign(S2Point const& a, S2Point const& b, S2Point const& c);

int ExactSign(S2Point const& a, S2Point const& b, S2Point const& c,
              bool perturb);

int SymbolicallyPerturbedSign(
    Vector3_xf const& a, Vector3_xf const& b,
    Vector3_xf const& c, Vector3_xf const& b_cross_c);

template <class T>
int TriageCompareCosDistances(Vector3<T> const& x,
                              Vector3<T> const& a, Vector3<T> const& b);

template <class T>
int TriageCompareSin2Distances(Vector3<T> const& x,
                               Vector3<T> const& a, Vector3<T> const& b);

int ExactCompareDistances(Vector3_xf const& x,
                          Vector3_xf const& a, Vector3_xf const& b);

int SymbolicCompareDistances(S2Point const& x,
                             S2Point const& a, S2Point const& b);

template <class T>
int TriageCompareSin2Distance(Vector3<T> const& x, Vector3<T> const& y, T r2);

template <class T>
int TriageCompareCosDistance(Vector3<T> const& x, Vector3<T> const& y, T r2);

int ExactCompareDistance(Vector3_xf const& x, Vector3_xf const& y,
                         ExactFloat const& r2);

template <class T>
int TriageCompareEdgeDistance(Vector3<T> const& x, Vector3<T> const& a0,
                              Vector3<T> const& a1, T r2);

int ExactCompareEdgeDistance(S2Point const& x, S2Point const& a0,
                             S2Point const& a1, S1ChordAngle r);

template <class T>
int TriageCompareEdgeDirections(
    Vector3<T> const& a0, Vector3<T> const& a1,
    Vector3<T> const& b0, Vector3<T> const& b1);

int ExactCompareEdgeDirections(Vector3_xf const& a0, Vector3_xf const& a1,
                               Vector3_xf const& b0, Vector3_xf const& b1);

template <class T>
int TriageEdgeCircumcenterSign(Vector3<T> const& x0, Vector3<T> const& x1,
                               Vector3<T> const& a, Vector3<T> const& b,
                               Vector3<T> const& c, int abc_sign);

int ExactEdgeCircumcenterSign(Vector3_xf const& x0, Vector3_xf const& x1,
                              Vector3_xf const& a, Vector3_xf const& b,
                              Vector3_xf const& c, int abc_sign);

int SymbolicEdgeCircumcenterSign(
    S2Point const& x0, S2Point const& x1,
    S2Point const& a_arg, S2Point const& b_arg, S2Point const& c_arg);

template <class T>
Excluded TriageVoronoiSiteExclusion(Vector3<T> const& a, Vector3<T> const& b,
                                    Vector3<T> const& x0, Vector3<T> const& x1,
                                    T r2);

Excluded ExactVoronoiSiteExclusion(Vector3_xf const& a, Vector3_xf const& b,
                                   Vector3_xf const& x0, Vector3_xf const& x1,
                                   ExactFloat const& r2);
}  // namespace s2pred

#endif  // S2_S2PREDICATES_INTERNAL_H_
