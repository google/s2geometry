// Copyright 2009 Google Inc. All Rights Reserved.
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
// Author: jyrki@google.com (Jyrki Alakuijala)
//
// Tests for the fast hashing algorithm based on Austin Appleby's
// MurmurHash 2.0 algorithm. See http://murmurhash.googlepages.com/
//
// Testing is split into small and large parts. The small test runs in
// 10 seconds and the large takes around 30 seconds in optimized mode, and
// substantially longer in debug.

#include <stdint.h>

#include "base/integral_types.h"
#include "gtest/gtest.h"
#include "util/endian/endian.h"
#include "util/hash/hasheval/string_hash_eval.h"
#include "util/hash/murmur.h"
#include "util/random/acmrandom.h"


namespace {

// Reference public domain implementation from
// http://murmurhash.googlepages.com/
// If you need this code or parts of it in production code, move it first to
// third_party.
uint64_t MurmurHash64Reference(const void * key, int len) {
  unsigned int seed = 0;
  const uint64_t m = 0xc6a4a7935bd1e995ULL;
  const int r = 47;

  uint64_t h = seed ^ (len * m);

  const uint64_t * data = (const uint64_t *)key;
  const uint64_t * end = data + (len/8);

  while (data != end) {
    uint64_t k = LittleEndian::Load64(data++);

    k *= m;
    k ^= k >> r;
    k *= m;

    h ^= k;
    h *= m;
  }

  const unsigned char * data2 = (const unsigned char*)data;

  switch (len & 7) {
    case 7:
      h ^= uint64_t(data2[6]) << 48;
      FALLTHROUGH_INTENDED;
    case 6:
      h ^= uint64_t(data2[5]) << 40;
      FALLTHROUGH_INTENDED;
    case 5:
      h ^= uint64_t(data2[4]) << 32;
      FALLTHROUGH_INTENDED;
    case 4:
      h ^= uint64_t(data2[3]) << 24;
      FALLTHROUGH_INTENDED;
    case 3:
      h ^= uint64_t(data2[2]) << 16;
      FALLTHROUGH_INTENDED;
    case 2:
      h ^= uint64_t(data2[1]) << 8;
      FALLTHROUGH_INTENDED;
    case 1:
      h ^= uint64_t(data2[0]);
      h *= m;
  }

  h ^= h >> r;
  h *= m;
  h ^= h >> r;

  return h;
}

}  // namespace

namespace util_hash {

TEST(RandomStrings, AgainstReferenceImplementation) {
  char data[1024] = { 0 };
  const int kIters = 1000000;
  ACMRandom rnd(ACMRandom::DeterministicSeed());
  for (int i = 0; i < kIters; ++i) {
    const int min_length = 1;
    const int len = min_length + rnd.Uniform(sizeof(data) - min_length);
    InitializeRandomString(len, &data[0], &rnd);
    const uint64 hash = MurmurHash64(&data[0], len);
    uint64 reference_hash = MurmurHash64Reference(&data[0], len);
    if (reference_hash >> 32 == 0) {
      reference_hash |= 1ULL << 32;
    }
    EXPECT_EQ(reference_hash, hash);
  }
}

}  // namespace util_hash