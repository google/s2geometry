// Copyright 2010 Google Inc. All Rights Reserved.
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
// This file provides CityHash64() and related functions.
//
// The externally visible functions follow the naming conventions of
// hash.h, where the size of the output is part of the name.  For
// example, CityHash64 returns a 64-bit hash.  The internal helpers do
// not have the return type in their name, but instead have names like
// HashLenXX or HashLenXXtoYY, where XX and YY are input string lengths.
//
// Most of the constants and tricks here were copied from murmur.cc or
// hash.h, or discovered by trial and error.  It's probably possible to further
// optimize the code here by writing a program that systematically explores
// more of the space of possible hash functions, or by using SIMD instructions.

#include "util/hash/city.h"


#include "util/hash/farmhash.h"
#include "util/hash/hash.h"  // for HashSeed()

namespace util_hash {

uint64 CityHash64(const char *s, size_t len) {
  return farmhash::Hash64(s, len) ^ HashSeed();
}

uint64 CityHash64WithSeed(const char *s, size_t len, uint64 seed) {
  return farmhash::Hash64WithSeed(s, len, seed + HashSeed());
}

uint64 CityHash64WithSeeds(const char *s, size_t len,
                           uint64 seed0, uint64 seed1) {
  return farmhash::Hash64WithSeeds(s, len, seed0, seed1 + HashSeed());
}

uint128 CityHash128(const char *s, size_t len) {
  return farmhashcc::CityHash128WithSeed(s, len, uint128(9, 999));
}

uint128 CityHash128WithSeed(const char *s, size_t len, uint128 seed) {
  return farmhashcc::CityHash128WithSeed(s, len, seed);
}

uint32 CityHash32(const char *s, size_t len) {
  return farmhash::Hash32(s, len) ^ HashSeed();
}

uint32 CityHash32WithSeed(const char *s, size_t len, uint32 seed) {
  return farmhash::Hash32WithSeed(s, len, seed - HashSeed());
}

}  // namespace util_hash
