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

//         jyrki@google.com (Jyrki Alakuijala)
//
// MurmurHash is a fast multiplication and shifting based algorithm,
// based on Austin Appleby's MurmurHash 2.0 algorithm.
// See http://murmurhash.googlepages.com/ for more details.
// Extension from 64 to 128 bits is by Jyrki Alakuijala.
//
// Its speed is about 3x slower than our fast Adler32 implementation,
// but still around 7x faster than Fingerprint96, and 2x faster than
// normal 64 bit Fingerprinting. Multiplication with a constant hashes
// bits efficiently to the higher bits. A shift by 47 positions and
// xor mixes the lower 17 bits.
//
// Here are the benchmark results. The number after the slash is
// log2 of the used size (10 means 1024 bytes, 25 means 32 MB).
//
// CPU: Intel Core2 (2 cores) dL1:32KB dL2:4096KB
// Benchmark                      Time(ns)    CPU(ns) Iterations
// -------------------------------------------------
// BM_MurmurHash64/2                    14         13  317200898
// BM_MurmurHash64/4                    31         27  159564526
// BM_MurmurHash64/7                    77         75   55847573
// BM_MurmurHash64/10                  499        457    8973803
// BM_MurmurHash64/25             16004739   15197974        269
// BM_Fprint64/2                        27         24  173256462
// BM_Fprint64/4                        28         26  162026910
// BM_Fprint64/7                       139        130   34088779
// BM_Fprint64/10                      969        874    4666371
// BM_Fprint64/25                 31474522   29446278        144
// BM_Fprint96/2                       339        303   10000000
// BM_Fprint96/4                       296        287   14521913
// BM_Fprint96/7                       686        660    6363241
// BM_Fprint96/10                     3508       3240    1000000
// BM_Fprint96/25                107542968   97926110        100
//
// There are no special reserved values in the value space.
//
// Murmur is not only faster than Fingerprint, but better too.
// The quality of the hash is much better than FP64 (Jenkins' hash),
// with less collisions found even in practical applications.

#ifndef S2GEOMETRY_UTIL_HASH_MURMUR_H_
#define S2GEOMETRY_UTIL_HASH_MURMUR_H_

#include <stddef.h>
#include <stdlib.h>  // for size_t.

#include "base/integral_types.h"

namespace util_hash {

// Like MurmurHash64 but with a seed which allows this hash function to be used
// in algorithms that need a family of parameterized hash functions
// e.g. Minhash.
uint64 MurmurHash64WithSeed(const char *buf, size_t len, uint64 seed);

// Hash function for a byte array.
// The hash mapping of MurmurHash64 will never change.
inline uint64 MurmurHash64(const char *buf, const size_t len) {
  return MurmurHash64WithSeed(buf, len, 0ULL);
}

// Class to compute MurmurHash64 compatible hash values from substrings.
class MurmurCat {
 public:
  MurmurCat() : hash_(0), last_val_(0), offset_(0) {}

  // Initializer. For MurmurHash64 compatibility, use the total length
  // of all the substrings. Use a non-zero fixed constant (like 1) if
  // it is expensive to compute the actual length and if the compatibility
  // with MurmurHash64(...) is not important.
  void Init(uint64 seed, size_t total_len);

  // Append all substrings by this method.
  void Append(const char *buf, const size_t len);

  // Get the final hash after all Append(...) calls.
  uint64 GetHash() const;

  // Convenience function
  static uint64 Hash64WithSeed(const char *buf, size_t len, uint64 seed) {
    return MurmurHash64WithSeed(buf, len, seed);
  }

 private:
  uint64 hash_;
  uint64 last_val_;
  int offset_;
};

}  // namespace util_hash

#endif  // S2GEOMETRY_UTIL_HASH_MURMUR_H_
