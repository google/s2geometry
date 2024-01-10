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

#ifndef S2_S2POINT_H_
#define S2_S2POINT_H_

#include <utility>

#include "absl/base/attributes.h"
#include "absl/hash/hash.h"
#include "s2/util/coding/coder.h"
#include "s2/_fp_contract_off.h"
#include "s2/s2coder.h"
#include "s2/s2error.h"
#include "s2/util/math/vector.h"  // IWYU pragma: export

// An S2Point represents a point on the unit sphere as a 3D vector.  Usually
// points are normalized to be unit length, but some methods do not require
// this.  See util/math/vector.h for the methods available.  Among other
// things, there are overloaded operators that make it convenient to write
// arithmetic expressions (e.g. (1-x)*p1 + x*p2).
class S2Point : public Vector3_d {
  using ValType = double;

 public:
  typedef s2coding::S2BasicCoder<S2Point> Coder;

  // Inherit base class constructors.
  using Base = Vector3_d;
  using Base::Base;

  // Due to an ambiguity in original C++11 specificiation, it was unclear
  // whether imported base class default constructors should be considered
  // when deciding to delete the default constructor of a class.  GCC and
  // Clang both accept the base class default ctor, while MSVC 2017 and
  // later do not.
  // The ambiguity was resolved in
  // https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2015/p0136r1.html
  // and accepted as part of the C++17 approval process, with the intent of
  // being retroactively applied to C++11.
  // However, while MSVC accepted the change for C++17 and forward, it did
  // not implement the change for prior versions. See their conformance page:
  // https://learn.microsoft.com/en-us/cpp/overview/visual-cpp-language-conformance?view=msvc-170
  // (search for P0136R1).
  // The explicit declaration here of a default ctor is a workaround until
  // this codebase is targeted at C++17 and above, at which point MSVC, GCC
  // and Clang will all have identical behavior and won't require this
  // explicit declaration of a default ctor.
  S2Point() = default;

  // When S2Point was defined as a Vector3_d we could mix and match the two
  // names.  With inheritance upcasting to a Vector3_d is easy, but we need to
  // explicitly allow the other direction, even though there's no data to
  // modify.  These are not marked explicit because this translation wasn't
  // explicit before.

  // NOLINTNEXTLINE(google-explicit-constructor)
  S2Point(const Base& base) : Base(base) {}

  // NOLINTNEXTLINE(google-explicit-constructor)
  S2Point(Base&& base) : Base(std::move(base)) {}

  // Initialize S2Point from a Decoder instance.
  bool Init(Decoder* decoder, S2Error& error) {
    if (decoder->avail() < sizeof(S2Point)) {
      error.Init(S2Error::DATA_LOSS, "Not enough data to decode S2Point");
      return false;
    }

    x(decoder->getdouble());
    y(decoder->getdouble());
    z(decoder->getdouble());
    return true;
  }

  S2Point& operator=(const Base& base) {
    Base::operator=(base);
    return *this;
  }

  S2Point& operator=(Base&& base) {
    Base::operator=(std::move(base));
    return *this;
  }

  // We can freely convert between S2Point and Vector3_d with no cost, but there
  // are a few corner cases where returning a Vector3_d can cause problems.
  // Notably the type of a ternary operator is evaluated independently of the
  // type being assigned to, so something like:
  //
  //   S2Point pnt = (x < 0) ? start_pnt : start_pnt + step;
  // (where start_pnt and step are both S2Point)
  //
  // Would fail to compile because start_pnt is an S2Point but start_pnt + step
  // is a Vector3_d, which is likely surprising to people.  So add overloads for
  // functions that return a Vector3_d to force return types to be covariant.
  S2Point& operator+=(const S2Point& b) {
    Base::operator+=(b);
    return *this;
  }
  S2Point& operator-=(const S2Point& b) {
    Base::operator-=(b);
    return *this;
  }
  S2Point& operator*=(const ValType& v) {
    Base::operator*=(v);
    return *this;
  }
  S2Point& operator/=(const ValType& v) {
    Base::operator/=(v);
    return *this;
  }

  S2Point operator+(const S2Point& b) const { return Base::operator+(b); }
  S2Point operator-(const S2Point& b) const { return Base::operator-(b); }
  S2Point operator*(const ValType& v) const { return Base::operator*(v); }
  S2Point operator/(const ValType& v) const { return Base::operator/(v); }

  friend S2Point operator-(const S2Point& pnt) {
    return -static_cast<const Base&>(pnt);
  }

  template <typename T>
  static S2Point Cast(const Vector3<T>& b) {
    return Base::Cast(b);
  }

  S2Point MulComponents(const S2Point& b) const {
    return Base::MulComponents(b);
  }
  S2Point DivComponents(const S2Point& b) const {
    return Base::DivComponents(b);
  }

  friend S2Point Max(const S2Point& a, const S2Point& b) {
    return Max(static_cast<const Base&>(a), static_cast<const Base&>(b));
  }

  friend S2Point Min(const S2Point& a, const S2Point& b) {
    return Min(static_cast<const Base&>(a), static_cast<const Base&>(b));
  }

  S2Point Normalize() const { return Base::Normalize(); }
  S2Point Sqrt() const { return Base::Sqrt(); }
  S2Point Floor() const { return Base::Floor(); }
  S2Point Ceil() const { return Base::Ceil(); }
  S2Point FRound() const { return Base::FRound(); }
  static S2Point NaN() { return Base::NaN(); }

  void Encode(Encoder* encoder) const {
    encoder->Ensure(sizeof(S2Point));
    encoder->putn(Data(), sizeof(S2Point));
  }
};

// S2PointHash can be used with standard containers (e.g., unordered_set) or
// nonstandard extensions (e.g., hash_map).  It is defined such that if two
// S2Points compare equal to each other, they have the same hash.  (This
// requires that positive and negative zero hash to the same value.)
using S2PointHash = absl::Hash<S2Point>;

#endif  // S2_S2POINT_H_
