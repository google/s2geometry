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

#include "util/hash/murmur.h"

#include "base/integral_types.h"
#include <glog/logging.h>
#include "util/endian/endian.h"

namespace util_hash {

namespace {

// Load the last remaining bytes in a little-endian manner (the byte that
// with the highest index ends up in the highest bits in the return value)
// into a uint64 to be used for checksumming those bytes that could not be
// directly loaded as a uint64.
inline uint64 LoadBytes(const char * const buf, int len) {
  DCHECK_LT(len, 9);
  uint64 val = 0;
  --len;
  do {
    val = (val << 8) | reinterpret_cast<const uint8 *>(buf)[len];
  } while (--len >= 0);
  // (--len >= 0) is about 10 % faster in the small string ubenchmarks
  // than (len--).
  return val;
}

// We need to mix some of the bits that get propagated and mixed into the
// high bits by multiplication back into the low bits. 17 last bits get
// a more efficiently mixed with this.
inline uint64 ShiftMix(uint64 val) {
  return val ^ (val >> 47);
}

}  // namespace

void MurmurCat::Init(uint64 seed, size_t total_len) {
  static const uint64 mul = 0xc6a4a7935bd1e995ULL;
  hash_ = seed ^ (total_len * mul);
  offset_ = 0;
}

void MurmurCat::Append(const char *buf, size_t len) {
  static const uint64 mul = 0xc6a4a7935bd1e995ULL;
  if (offset_ + len < 8) {
    if (len == 0) {
      return;
    }
    // Not enough data to murmur, accumulate more.
    last_val_ |= LoadBytes(buf, len) << (offset_ * 8);
    offset_ += len;
    return;  // We don't have enough to hash yet.
  }
  if (offset_ != 0) {
    // We have enough data to murmur, and there is old data
    // in last_val_.
    int num_bytes = 8 - offset_;
    last_val_ |= LoadBytes(buf, num_bytes) << (offset_ * 8);

    const uint64 data = ShiftMix(last_val_ * mul) * mul;
    hash_ ^= data;
    hash_ *= mul;
    buf += num_bytes;
    len -= num_bytes;
  }
  const int len_aligned = len & ~0x7;
  const char * const end = buf + len_aligned;
  for (const char *p = buf; p != end; p += 8) {
    const uint64 data = ShiftMix(LittleEndian::Load64(p) * mul) * mul;
    hash_ ^= data;
    hash_ *= mul;
  }
  if ((len & 0x7) != 0) {
    offset_ = len & 0x7;
    last_val_ = LoadBytes(end, len & 0x7);
  } else {
    last_val_ = 0;
    offset_ = 0;
  }
}

uint64 MurmurCat::GetHash() const {
  static const uint64 mul = 0xc6a4a7935bd1e995ULL;
  uint64 hash = hash_;
  if ((offset_ & 0x7) != 0) {
    const uint64 data = last_val_;
    hash ^= data;
    hash *= mul;
  }
  hash = ShiftMix(hash) * mul;
  hash = ShiftMix(hash);
  return hash;
}

uint64 MurmurHash64WithSeed(const char *buf,
                            const size_t len,
                            const uint64 seed) {
  static const uint64 mul = 0xc6a4a7935bd1e995ULL;
  // Let's remove the bytes not divisible by the sizeof(uint64).
  // This allows the inner loop to process the data as 64 bit integers.
  const int len_aligned = len & ~0x7;
  const char * const end = buf + len_aligned;
  uint64 hash = seed ^ (len * mul);
  for (const char *p = buf; p != end; p += 8) {
    const uint64 data = ShiftMix(LittleEndian::Load64(p) * mul) * mul;
    hash ^= data;
    hash *= mul;
  }
  if ((len & 0x7) != 0) {
    const uint64 data = LoadBytes(end, len & 0x7);
    hash ^= data;
    hash *= mul;
  }
  hash = ShiftMix(hash) * mul;
  hash = ShiftMix(hash);
  return hash;
}

}  // namespace util_hash