// Copyright Google Inc. All Rights Reserved.
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

#ifndef S2_UTIL_ENDIAN_ENDIAN_H_
#define S2_UTIL_ENDIAN_ENDIAN_H_

#include <cstdint>

#include "s2/base/port.h"
#include "absl/numeric/bits.h"
#include "absl/numeric/int128.h"

namespace s2endian {
// TODO(user): Trim unused Encoder functions and drop int128 API.
inline absl::uint128 byteswap(absl::uint128 host_int) {
  return absl::MakeUint128(absl::byteswap(absl::Uint128Low64(host_int)),
                           absl::byteswap(absl::Uint128High64(host_int)));
}
}  // namespace s2endian

// Utilities to convert numbers between the current hosts's native byte
// order and little-endian byte order
//
// Load/Store methods are alignment safe
class LittleEndian {
 public:
  // Conversion functions.
  static uint16_t FromHost16(uint16_t x) { return byteswap_if_big_endian(x); }
  static uint16_t ToHost16(uint16_t x) { return byteswap_if_big_endian(x); }

  static uint32_t FromHost32(uint32_t x) { return byteswap_if_big_endian(x); }
  static uint32_t ToHost32(uint32_t x) { return byteswap_if_big_endian(x); }

  static uint64_t FromHost64(uint64_t x) { return byteswap_if_big_endian(x); }
  static uint64_t ToHost64(uint64_t x) { return byteswap_if_big_endian(x); }

  static absl::uint128 FromHost128(absl::uint128 x) {
    return byteswap_if_big_endian(x);
  }
  static absl::uint128 ToHost128(absl::uint128 x) {
    return byteswap_if_big_endian(x);
  }

  static constexpr bool IsLittleEndian() {
    return absl::endian::native == absl::endian::little;
  }

  // Functions to do unaligned loads and stores in little-endian order.
  static uint16_t Load16(const void* p) {
    return ToHost16(UNALIGNED_LOAD16(p));
  }

  static void Store16(void* p, uint16_t v) {
    UNALIGNED_STORE16(p, FromHost16(v));
  }

  static uint32_t Load32(const void* p) {
    return ToHost32(UNALIGNED_LOAD32(p));
  }

  static void Store32(void* p, uint32_t v) {
    UNALIGNED_STORE32(p, FromHost32(v));
  }

  static uint64_t Load64(const void* p) {
    return ToHost64(UNALIGNED_LOAD64(p));
  }

  static void Store64(void* p, uint64_t v) {
    UNALIGNED_STORE64(p, FromHost64(v));
  }

  static absl::uint128 Load128(const void* p) {
    return absl::MakeUint128(
        ToHost64(UNALIGNED_LOAD64(reinterpret_cast<const uint64_t*>(p) + 1)),
        ToHost64(UNALIGNED_LOAD64(p)));
  }

  static void Store128(void* p, const absl::uint128 v) {
    UNALIGNED_STORE64(p, FromHost64(absl::Uint128Low64(v)));
    UNALIGNED_STORE64(reinterpret_cast<uint64_t*>(p) + 1,
                      FromHost64(absl::Uint128High64(v)));
  }

 private:
  template <typename T>
  static T byteswap_if_big_endian(T x) {
    static_assert(absl::endian::native == absl::endian::big ||
                  absl::endian::native == absl::endian::little,
                  "Unsupported endianness.");
    if constexpr (absl::endian::native == absl::endian::big) {
      using absl::byteswap;
      using s2endian::byteswap;
      x = byteswap(x);
    }
    return x;
  }
};

#endif  // S2_UTIL_ENDIAN_ENDIAN_H_
