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

#ifndef UTIL_HASH_HASH_H_
#define UTIL_HASH_HASH_H_

#include "base/int128.h"
#include "base/integral_types.h"
#include "util/hash/hash128to64.h"
#include "util/hash/jenkins.h"
#include "util/hash/jenkins_lookup2.h"
#include "util/hash/string_hash.h"

inline uint32 HashTo32(const char *s, size_t slen) {
  uint32 retval = Hash32StringWithSeed(s, slen, MIX32);
  return retval == kIllegalHash32 ? static_cast<uint32>(retval-1) : retval;
}

HASH_NAMESPACE_DECLARATION_START

template<> struct hash<uint128> {
  size_t operator()(const uint128& x) const {
    if (sizeof(&x) == 8) {  // 64-bit systems have 8-byte pointers.
      return Hash128to64(x);
    } else {
      uint32 a = static_cast<uint32>(Uint128Low64(x)) +
          static_cast<uint32>(0xba79e1bd);
      uint32 b = static_cast<uint32>(Uint128Low64(x) >> 32) +
          static_cast<uint32>(0xba79e1bd);
      uint32 c = static_cast<uint32>(Uint128High64(x)) + MIX32;
      mix(a, b, c);
      a += static_cast<uint32>(Uint128High64(x) >> 32);
      mix(a, b, c);
      return c;
    }
  }
};

HASH_NAMESPACE_DECLARATION_END

// Various files define GoodFastHash specializations.  Declare it here
// so they don't have to.
template<class X> struct GoodFastHash;

#endif  // UTIL_HASH_HASH_H_
