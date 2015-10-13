// Copyright 2011 Google Inc. All Rights Reserved.
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


// This file provides hash functions for Vector3<T>, assuming T is a POD.
//
// WARNING: hashing floats and doubles is a tricky business due to multiple
// bit patterns having the same or similar meanings.  E.g., there are multiple
// bit patterns that mean "NaN."  The code here hashes +0 and -0 to the same
// hash value, but still may not be ideal for all situations.

#ifndef UTIL_MATH_VECTOR3_HASH_H_
#define UTIL_MATH_VECTOR3_HASH_H_

// IWYU sees __gnu_cxx::size_t and reaches for a non-standard gcc header
// IWYU pragma: no_include <ext/new_allocator.h>
#include <stddef.h>

#include "base/int128.h"
#include "base/integral_types.h"
#include <glog/logging.h>
#include "base/macros.h"
#include "base/port.h"
#include "base/type_traits.h"
#include "util/hash/hash.h"
#include "util/hash/string_hash.h"
#include "util/math/vector3.h"

namespace vector3_hash_internal {

inline uint64 CollapseZero(uint64 bits) {
  // IEEE 754 has two representations for zero, positive zero and negative
  // zero.  These two values compare as equal, and therefore we need them to
  // hash to the same value.
  //
  // We handle this by simply clearing the top bit of every value,
  // which clears the sign bit on both big-endian and little-endian
  // architectures.  This creates some additional hash collisions between
  // points that differ only in the sign of their components, but this is
  // rarely a problem with real data.
  //
  // The obvious alternative is to explicitly map all occurrences of positive
  // zero to negative zero (or vice versa), but this is more expensive and
  // makes the average case slower.
  //
  // We also mask off the low-bit because we've seen differences in
  // some floating point operations (specifically 'fcos' on i386)
  // between different implementations of the same architecure
  // (e.g. 'Xeon 5345' vs. 'Opteron 270').  It's unknown how many bits
  // of mask are sufficient to cover real world cases, but the intent
  // is to be as conservative as possible in discarding bits.
  return bits & GG_LONGLONG(0x7ffffffffffffffe);
}

// See CollapseZero().  This is the same idea, but in this case the
// input is to be treated as two 32-bit floats packed into a single
// 64-bit value.
inline uint64 CollapseZeroes(uint64 bits) {
  return bits & GG_LONGLONG(0x7ffffffe7ffffffe);
}

// A hash for vectors of floats or doubles.  We call CollapseZero() on each
// element of v because of +0 vs -0 etc.  See above.
template <class Vec> size_t FloatHash(const Vec& v) {
  const bool kIsFloat = std::is_same<const Vec&, const Vector3_f&>::value;
  const bool kIsDouble = std::is_same<const Vec&, const Vector3_d&>::value;
  CHECK(kIsFloat || kIsDouble);
  if (kIsFloat) {
    const char* data = reinterpret_cast<const char*>(v.Data());
    uint128 u(CollapseZeroes(UNALIGNED_LOAD64(data)),
              CollapseZeroes(UNALIGNED_LOAD64(data + 4)) + MIX64);
    return HASH_NAMESPACE::hash<uint128>()(u);
  } else {
    const char* data = reinterpret_cast<const char*>(v.Data());
    uint128 u(CollapseZero(UNALIGNED_LOAD64(data)),
              CollapseZero(UNALIGNED_LOAD64(data + 8)) + MIX64);
    u = uint128(HASH_NAMESPACE::hash<uint128>()(u),
                CollapseZero(UNALIGNED_LOAD64(data + 16)));
    return HASH_NAMESPACE::hash<uint128>()(u);
  }
}

};  // namespace vector3_hash_internal

template <class T> struct GoodFastHash;

// This hash function may change from time to time.
template <class VType> struct GoodFastHash<Vector3<VType> > {
  size_t operator()(const Vector3<VType>& v) const {
    COMPILE_ASSERT(std::is_pod<VType>::value, POD_expected);
    if (std::is_floating_point<VType>::value)
      return vector3_hash_internal::FloatHash(v);
    else
      return HashStringThoroughly(reinterpret_cast<const char*>(v.Data()),
                                  v.Size() * sizeof(v.Data()[0]));
  }
};

HASH_NAMESPACE_DECLARATION_START
template <class T> struct hash;

// This hash function may change from time to time.
template <class VType> struct hash<Vector3<VType> > {
  size_t operator()(const Vector3<VType>& v) const {
    return GoodFastHash<Vector3<VType> >()(v);
  }
};
HASH_NAMESPACE_DECLARATION_END

#endif  // UTIL_MATH_VECTOR3_HASH_H_
