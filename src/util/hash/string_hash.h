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

//
// These are the core hashing routines which operate on strings. We define
// strings loosely as a sequence of bytes, and these routines are designed to
// work with the most fundamental representations of a string of bytes.
//
// These routines provide "good" hash functions in terms of both quality and
// speed. Their values can and will change as their implementations change and
// evolve.

#ifndef UTIL_HASH_STRING_HASH_H_
#define UTIL_HASH_STRING_HASH_H_

#include <limits.h>
#include <stddef.h>

#include "base/port.h"
#include "base/integral_types.h"
#include "util/hash/city.h"
#include "util/hash/jenkins.h"
#include "util/hash/jenkins_lookup2.h"

namespace hash_internal {
// Arbitrary mix constants (pi).
static const uint32 kMix32 = 0x12b9b0a1UL;
static const uint64 kMix64 = GG_ULONGLONG(0x2b992ddfa23249d6);

template <size_t Bits = sizeof(size_t) * CHAR_BIT> struct Thoroughly;

template <>
struct Thoroughly<64> {
  static size_t Hash(const char* s, size_t len, size_t seed) {
    return static_cast<size_t>(util_hash::CityHash64WithSeed(s, len, seed));
  }
  static size_t Hash(const char* s, size_t len) {
    return static_cast<size_t>(util_hash::CityHash64(s, len));
  }
  static size_t Hash(const char* s, size_t len, size_t seed0, size_t seed1) {
    return static_cast<size_t>(util_hash::CityHash64WithSeeds(s, len,
                                                              seed0, seed1));
  }
};

template <>
struct Thoroughly<32> {
  static size_t Hash(const char* s, size_t len, size_t seed) {
    return static_cast<size_t>(
        util_hash::CityHash32WithSeed(s, len, static_cast<uint32>(seed)));
  }
  static size_t Hash(const char* s, size_t len) {
    return static_cast<size_t>(util_hash::CityHash32(s, len));
  }
  static size_t Hash(const char* s, size_t len, size_t seed0, size_t seed1) {
    seed0 += Hash(s, len / 2, seed1);
    s += len / 2;
    len -= len / 2;
    return Hash(s, len, seed0);
  }
};
}  // namespace hash_internal

// We use different algorithms depending on the size of size_t.
inline size_t HashStringThoroughlyWithSeed(const char* s, size_t len,
                                           size_t seed) {
  return hash_internal::Thoroughly<>::Hash(s, len, seed);
}

inline size_t HashStringThoroughly(const char* s, size_t len) {
  return hash_internal::Thoroughly<>::Hash(s, len);
}

inline size_t HashStringThoroughlyWithSeeds(const char* s, size_t len,
                                            size_t seed0, size_t seed1) {
  return hash_internal::Thoroughly<>::Hash(s, len, seed0, seed1);
}

#endif  // UTIL_HASH_STRING_HASH_H_
