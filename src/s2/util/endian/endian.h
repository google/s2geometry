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

#include "absl/numeric/int128.h"
#include "s2/util/gtl/unaligned.h"

// Check for absl::byteswap availability (added in abseil LTS 2025.01.27)
#if __has_include("absl/numeric/bits.h")
#include "absl/numeric/bits.h"
#endif

// Provide fallback byteswap implementations for older abseil versions
namespace s2endian {

// Fallback byteswap for uint16_t
inline uint16_t byteswap_u16(uint16_t x) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_bswap16(x);
#else
  return static_cast<uint16_t>((x >> 8) | (x << 8));
#endif
}

// Fallback byteswap for uint32_t
inline uint32_t byteswap_u32(uint32_t x) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_bswap32(x);
#else
  return ((x & 0xFF000000U) >> 24) |
         ((x & 0x00FF0000U) >> 8) |
         ((x & 0x0000FF00U) << 8) |
         ((x & 0x000000FFU) << 24);
#endif
}

// Fallback byteswap for uint64_t
inline uint64_t byteswap_u64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_bswap64(x);
#else
  return ((x & 0xFF00000000000000ULL) >> 56) |
         ((x & 0x00FF000000000000ULL) >> 40) |
         ((x & 0x0000FF0000000000ULL) >> 24) |
         ((x & 0x000000FF00000000ULL) >> 8) |
         ((x & 0x00000000FF000000ULL) << 8) |
         ((x & 0x0000000000FF0000ULL) << 24) |
         ((x & 0x000000000000FF00ULL) << 40) |
         ((x & 0x00000000000000FFULL) << 56);
#endif
}

// TODO(user): Trim unused Encoder functions and drop int128 API.
inline absl::uint128 byteswap(absl::uint128 host_int) {
  return absl::MakeUint128(byteswap_u64(absl::Uint128Low64(host_int)),
                           byteswap_u64(absl::Uint128High64(host_int)));
}

}  // namespace s2endian

// Detect endianness at compile time
namespace s2endian_detail {
#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__) && \
    defined(__ORDER_BIG_ENDIAN__)
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
constexpr bool kIsLittleEndian = true;
#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
constexpr bool kIsLittleEndian = false;
#else
#error "Unknown endianness"
#endif
#elif defined(_WIN32) || defined(__x86_64__) || defined(__i386__) || \
      defined(__aarch64__) || defined(__arm__)
// Most common platforms are little-endian
constexpr bool kIsLittleEndian = true;
#else
// Fallback: assume little-endian (safe for most platforms)
constexpr bool kIsLittleEndian = true;
#endif
}  // namespace s2endian_detail

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
    return s2endian_detail::kIsLittleEndian;
  }

  // Functions to do unaligned loads and stores in little-endian order.
  template <typename T>
  static T Load(const void* p) {
    const char* pc = reinterpret_cast<const char*>(p);
    return byteswap_if_big_endian(gtl::UnalignedLoad<T>(pc));
  }

  static uint16_t Load16(const void* p) {
    const char* pc = reinterpret_cast<const char*>(p);
    return ToHost16(gtl::UnalignedLoad<uint16_t>(pc));
  }

  static void Store16(void* p, uint16_t v) {
    char* pc = reinterpret_cast<char*>(p);
    gtl::UnalignedStore<uint16_t>(FromHost16(v), pc);
  }

  static uint32_t Load32(const void* p) {
    const char* pc = reinterpret_cast<const char*>(p);
    return ToHost32(gtl::UnalignedLoad<uint32_t>(pc));
  }

  static void Store32(void* p, uint32_t v) {
    char* pc = reinterpret_cast<char*>(p);
    gtl::UnalignedStore<uint32_t>(FromHost32(v), pc);
  }

  static uint64_t Load64(const void* p) {
    const char* pc = reinterpret_cast<const char*>(p);
    return ToHost64(gtl::UnalignedLoad<uint64_t>(pc));
  }

  static void Store64(void* p, uint64_t v) {
    char* pc = reinterpret_cast<char*>(p);
    gtl::UnalignedStore<uint64_t>(FromHost64(v), pc);
  }

  static absl::uint128 Load128(const void* p) {
    const char* pc = reinterpret_cast<const char*>(p);
    return absl::MakeUint128(
        ToHost64(gtl::UnalignedLoad<uint64_t>(pc + sizeof(uint64_t))),
        ToHost64(gtl::UnalignedLoad<uint64_t>(pc)));
  }

  static void Store128(void* p, const absl::uint128 v) {
    char* pc = reinterpret_cast<char*>(p);
    gtl::UnalignedStore<uint64_t>(FromHost64(absl::Uint128Low64(v)), pc);
    gtl::UnalignedStore<uint64_t>(FromHost64(absl::Uint128High64(v)),
                                  pc + sizeof(uint64_t));
  }

 private:
  // Byteswap implementations for each type
  static uint16_t do_byteswap(uint16_t x) { return s2endian::byteswap_u16(x); }
  static uint32_t do_byteswap(uint32_t x) { return s2endian::byteswap_u32(x); }
  static uint64_t do_byteswap(uint64_t x) { return s2endian::byteswap_u64(x); }
  static absl::uint128 do_byteswap(absl::uint128 x) { return s2endian::byteswap(x); }

  template <typename T>
  static T byteswap_if_big_endian(T x) {
    if constexpr (!s2endian_detail::kIsLittleEndian) {
      return do_byteswap(x);
    }
    return x;
  }
};

#endif  // S2_UTIL_ENDIAN_ENDIAN_H_
