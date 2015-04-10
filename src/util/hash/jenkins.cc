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

#include "util/hash/jenkins.h"


#include <glog/logging.h>
#include "util/hash/jenkins_lookup2.h"
// These functions can be updated using farmhash later.
// The current version is copied from mobile/util/hash/hash.cc.
static const uint32 kPrimes32[16] ={
  65537, 65539, 65543, 65551, 65557, 65563, 65579, 65581,
  65587, 65599, 65609, 65617, 65629, 65633, 65647, 65651,
};

static const uint64 kPrimes64[] ={
  GG_ULONGLONG(4294967311), GG_ULONGLONG(4294967357),
  GG_ULONGLONG(4294967371), GG_ULONGLONG(4294967377),
  GG_ULONGLONG(4294967387), GG_ULONGLONG(4294967389),
  GG_ULONGLONG(4294967459), GG_ULONGLONG(4294967477),
  GG_ULONGLONG(4294967497), GG_ULONGLONG(4294967513),
  GG_ULONGLONG(4294967539), GG_ULONGLONG(4294967543),
  GG_ULONGLONG(4294967549), GG_ULONGLONG(4294967561),
  GG_ULONGLONG(4294967563), GG_ULONGLONG(4294967569)
};

uint32 Hash32StringWithSeedReferenceImplementation(const char *s,
                                                   size_t len,
                                                   uint32 seed) {
  uint32 n = seed;
  size_t prime1 = 0, prime2 = 8;  // Indices into kPrimes32
  union {
    uint16 n;
    char bytes[sizeof(uint16)];
  } chunk;
  for (const char *i = s, *const end = s + len; i != end; ) {
    chunk.bytes[0] = *i++;
    chunk.bytes[1] = i == end ? 0 : *i++;
    n = n * kPrimes32[prime1++] ^ chunk.n * kPrimes32[prime2++];
    prime1 &= 0x0F;
    prime2 &= 0x0F;
  }
  return n;
}

uint32 Hash32StringWithSeed(const char *s, size_t len, uint32 seed) {
  uint32 a, b;
  b = 0x12f905ffUL;
  uint32 c = seed;
  while (len > 12) {
    a = Hash32StringWithSeedReferenceImplementation(s, 12, c);
    mix(a, b, c);
    s += 12;
    len -= 12;
  }
  if (len > 0) {
    a = Hash32StringWithSeedReferenceImplementation(s, len, c);
    mix(a, b, c);
  }
  return c;
}

uint64 Hash64StringWithSeed(const char *s, size_t len, uint64 seed) {
  uint64 n = seed;
  size_t prime1 = 0, prime2 = 8;  // Indices into kPrimes64
  union {
    uint32 n;
    char bytes[sizeof(uint32)];
  } chunk;
  for (const char *i = s, *const end = s + len; i != end; ) {
    chunk.bytes[0] = *i++;
    chunk.bytes[1] = i == end ? 0 : *i++;
    chunk.bytes[2] = i == end ? 0 : *i++;
    chunk.bytes[3] = i == end ? 0 : *i++;
    n = n * kPrimes64[prime1++] ^ chunk.n * kPrimes64[prime2++];
    prime1 &= 0x0F;
    prime2 &= 0x0F;
  }
  return n;
}
