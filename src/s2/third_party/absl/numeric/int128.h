// Copyright 2004 Google Inc. All Rights Reserved.
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

// All Rights Reserved.
//

#ifndef S2_THIRD_PARTY_ABSL_NUMERIC_INT128_H_
#define S2_THIRD_PARTY_ABSL_NUMERIC_INT128_H_

#include <cmath>
#include <cstring>
#include <iosfwd>
#include <limits>

#include "s2/third_party/absl/base/config.h"
#include "s2/third_party/absl/base/integral_types.h"
#include "s2/third_party/absl/base/macros.h"
#include "s2/third_party/absl/base/port.h"

// An unsigned 128-bit integer type. The API is meant to mimic an intrinsic type
// as closely as practical including cases exhibiting undefined behavior (e.g.
// division by zero).
//
// Supports:
//   - Implicit constexpr construction from integral types.
//   - Explicit constexpr conversion to integral types.
//   - Explicit construction from and conversion to floating point types.
//   - Interaction with __int128 where the type is available.
//
// When constructing with values greater than 2**64 use the constexpr
// absl::MakeUint128 factory function declared below. e.g.
//
//     uint128 big = absl::MakeUint128(1, 0);
//
// Ways in which this API differs from a built-in integer:
//   - Implicit conversion that does not preserve value is valid for built-in
//     types except when treated as an error by the compiler. Our type behaves
//     as if these compiler errors (e.g. Wnarrowing, Wconversion) are always on.
//   - A built-in integer will implicitly convert to a floating point type when
//     interacting through arithmetic operators. Our type must first be
//     explicitly cast to a floating point type.
// TODO(user) Remove zeroing behavior when shifting by >= 128; should be UB.
//
// This type will eventually be removed and uint128 will be an alias for a
// standard uint128_t once it exists. Code written with this type will continue
// to compile when that happens, provided replacement helper functions
// Uint128(Low|High)64 and MakeUint128 are made. The leftover cruft to be
// removed will be unnecessary static_cast's when doing arithmetic operations
// between uint128 and floating point types.
//
// Note that the alignment requirement of uint128 is due to change, so users
// should take care to avoid depending on the current 8 byte alignment.
// TODO(user) Remove alignment note above once alignof(uint128) becomes 16.
class uint128;

namespace absl {

constexpr uint128 MakeUint128(uint64 top, uint64 bottom);

}  // namespace absl

class uint128 {
 public:
  uint128() = default;

  // Constructors from arithmetic types
  constexpr uint128(int v);                 // NOLINT(runtime/explicit)
  constexpr uint128(unsigned int v);        // NOLINT(runtime/explicit)
  constexpr uint128(long v);                // NOLINT(runtime/int)
  constexpr uint128(unsigned long v);       // NOLINT(runtime/int)
  constexpr uint128(long long v);           // NOLINT(runtime/int)
  constexpr uint128(unsigned long long v);  // NOLINT(runtime/int)
#ifdef ABSL_HAVE_INTRINSIC_INT128
  constexpr uint128(__int128 v);           // NOLINT(runtime/explicit)
  constexpr uint128(unsigned __int128 v);  // NOLINT(runtime/explicit)
#endif  // ABSL_HAVE_INTRINSIC_INT128
  explicit uint128(float v);        // NOLINT(runtime/explicit)
  explicit uint128(double v);       // NOLINT(runtime/explicit)
  explicit uint128(long double v);  // NOLINT(runtime/explicit)

  // Assignment operators from arithmetic types
  uint128& operator=(int v);
  uint128& operator=(unsigned int v);
  uint128& operator=(long v);                // NOLINT(runtime/int)
  uint128& operator=(unsigned long v);       // NOLINT(runtime/int)
  uint128& operator=(long long v);           // NOLINT(runtime/int)
  uint128& operator=(unsigned long long v);  // NOLINT(runtime/int)
#ifdef ABSL_HAVE_INTRINSIC_INT128
  uint128& operator=(__int128 v);
  uint128& operator=(unsigned __int128 v);
#endif  // ABSL_HAVE_INTRINSIC_INT128

  // Conversion operators to other arithmetic types
  constexpr explicit operator bool() const;
  constexpr explicit operator char() const;
  constexpr explicit operator signed char() const;
  constexpr explicit operator unsigned char() const;
  constexpr explicit operator char16_t() const;
  constexpr explicit operator char32_t() const;
  constexpr explicit operator wchar_t() const;
  constexpr explicit operator short() const;  // NOLINT(runtime/int)
  // NOLINTNEXTLINE(runtime/int)
  constexpr explicit operator unsigned short() const;
  constexpr explicit operator int() const;
  constexpr explicit operator unsigned int() const;
  constexpr explicit operator long() const;  // NOLINT(runtime/int)
  // NOLINTNEXTLINE(runtime/int)
  constexpr explicit operator unsigned long() const;
  // NOLINTNEXTLINE(runtime/int)
  constexpr explicit operator long long() const;
  // NOLINTNEXTLINE(runtime/int)
  constexpr explicit operator unsigned long long() const;
#ifdef ABSL_HAVE_INTRINSIC_INT128
  constexpr explicit operator __int128() const;
  constexpr explicit operator unsigned __int128() const;
#endif  // ABSL_HAVE_INTRINSIC_INT128
  explicit operator float() const;
  explicit operator double() const;
  explicit operator long double() const;

  // Trivial copy constructor, assignment operator and destructor.

  // Arithmetic operators.
  uint128& operator+=(const uint128& b);
  uint128& operator-=(const uint128& b);
  uint128& operator*=(const uint128& b);
  // Long division/modulo for uint128.
  uint128& operator/=(const uint128& b);
  uint128& operator%=(const uint128& b);
  uint128 operator++(int);
  uint128 operator--(int);
  uint128& operator<<=(int);
  uint128& operator>>=(int);
  uint128& operator&=(const uint128& b);
  uint128& operator|=(const uint128& b);
  uint128& operator^=(const uint128& b);
  uint128& operator++();
  uint128& operator--();

  friend uint64 Uint128Low64(const uint128& v);
  friend uint64 Uint128High64(const uint128& v);
  friend constexpr uint128 absl::MakeUint128(uint64 top, uint64 bottom);

 private:
  constexpr uint128(uint64 top, uint64 bottom);

  // TODO(user) Update implementation to use __int128 once all users of
  // uint128 are fixed to not depend on alignof(uint128) == 8. Also add
  // alignas(16) to class definition to keep alignment consistent across
  // platforms.
#if defined(ABSL_IS_LITTLE_ENDIAN)
  uint64 lo_;
  uint64 hi_;
#elif defined(ABSL_IS_BIG_ENDIAN)
  uint64 hi_;
  uint64 lo_;
#else  // byte order
#error "Unsupported byte order: must be little-endian or big-endian."
#endif  // byte order
};

extern const uint128 kuint128max;

// allow uint128 to be logged
extern std::ostream& operator<<(std::ostream& o, const uint128& b);

// TODO(user) add operator>>(std::istream&, uint128&)

// Methods to access low and high pieces of 128-bit value.
uint64 Uint128Low64(const uint128& v);
uint64 Uint128High64(const uint128& v);

// TODO(b/31950287): Implement signed 128-bit type

// --------------------------------------------------------------------------
//                      Implementation details follow
// --------------------------------------------------------------------------

namespace absl {

inline constexpr uint128 MakeUint128(uint64 top, uint64 bottom) {
  return uint128(top, bottom);
}

}  // namespace absl

// Assignment from integer types.

inline uint128& uint128::operator=(int v) {
  return *this = uint128(v);
}

inline uint128& uint128::operator=(unsigned int v) {
  return *this = uint128(v);
}

inline uint128& uint128::operator=(long v) {  // NOLINT(runtime/int)
  return *this = uint128(v);
}

// NOLINTNEXTLINE(runtime/int)
inline uint128& uint128::operator=(unsigned long v) {
  return *this = uint128(v);
}

// NOLINTNEXTLINE(runtime/int)
inline uint128& uint128::operator=(long long v) {
  return *this = uint128(v);
}

// NOLINTNEXTLINE(runtime/int)
inline uint128& uint128::operator=(unsigned long long v) {
  return *this = uint128(v);
}

#ifdef ABSL_HAVE_INTRINSIC_INT128
inline uint128& uint128::operator=(__int128 v) {
  return *this = uint128(v);
}

inline uint128& uint128::operator=(unsigned __int128 v) {
  return *this = uint128(v);
}
#endif  // ABSL_HAVE_INTRINSIC_INT128

// Shift and arithmetic operators.

inline uint128 operator<<(const uint128& lhs, int amount) {
  return uint128(lhs) <<= amount;
}

inline uint128 operator>>(const uint128& lhs, int amount) {
  return uint128(lhs) >>= amount;
}

inline uint128 operator+(const uint128& lhs, const uint128& rhs) {
  return uint128(lhs) += rhs;
}

inline uint128 operator-(const uint128& lhs, const uint128& rhs) {
  return uint128(lhs) -= rhs;
}

inline uint128 operator*(const uint128& lhs, const uint128& rhs) {
  return uint128(lhs) *= rhs;
}

inline uint128 operator/(const uint128& lhs, const uint128& rhs) {
  return uint128(lhs) /= rhs;
}

inline uint128 operator%(const uint128& lhs, const uint128& rhs) {
  return uint128(lhs) %= rhs;
}

inline uint64 Uint128Low64(const uint128& v) { return v.lo_; }

inline uint64 Uint128High64(const uint128& v) { return v.hi_; }

// Constructors from integer types.

#if defined(ABSL_IS_LITTLE_ENDIAN)

inline constexpr uint128::uint128(uint64 top, uint64 bottom)
    : lo_(bottom), hi_(top) {}

inline constexpr uint128::uint128(int v)
    : lo_(v), hi_(v < 0 ? std::numeric_limits<uint64>::max() : 0) {}
inline constexpr uint128::uint128(long v)  // NOLINT(runtime/int)
    : lo_(v), hi_(v < 0 ? std::numeric_limits<uint64>::max() : 0) {}
inline constexpr uint128::uint128(long long v)  // NOLINT(runtime/int)
    : lo_(v), hi_(v < 0 ? std::numeric_limits<uint64>::max() : 0) {}

inline constexpr uint128::uint128(unsigned int v) : lo_(v), hi_(0) {}
// NOLINTNEXTLINE(runtime/int)
inline constexpr uint128::uint128(unsigned long v) : lo_(v), hi_(0) {}
// NOLINTNEXTLINE(runtime/int)
inline constexpr uint128::uint128(unsigned long long v)
    : lo_(v), hi_(0) {}

#ifdef ABSL_HAVE_INTRINSIC_INT128
inline constexpr uint128::uint128(__int128 v)
    : lo_(static_cast<uint64>(v & ~uint64{0})),
      hi_(static_cast<uint64>(static_cast<unsigned __int128>(v) >> 64)) {}
inline constexpr uint128::uint128(unsigned __int128 v)
    : lo_(static_cast<uint64>(v & ~uint64{0})),
      hi_(static_cast<uint64>(v >> 64)) {}
#endif  // ABSL_HAVE_INTRINSIC_INT128

#elif defined(ABSL_IS_BIG_ENDIAN)

inline constexpr uint128::uint128(uint64 top, uint64 bottom)
    : hi_(top), lo_(bottom) {}

inline constexpr uint128::uint128(int v)
    : hi_(v < 0 ? std::numeric_limits<uint64>::max() : 0), lo_(v) {}
inline constexpr uint128::uint128(long v)  // NOLINT(runtime/int)
    : hi_(v < 0 ? std::numeric_limits<uint64>::max() : 0), lo_(v) {}
inline constexpr uint128::uint128(long long v)  // NOLINT(runtime/int)
    : hi_(v < 0 ? std::numeric_limits<uint64>::max() : 0), lo_(v) {}

inline constexpr uint128::uint128(unsigned int v) : hi_(0), lo_(v) {}
// NOLINTNEXTLINE(runtime/int)
inline constexpr uint128::uint128(unsigned long v) : hi_(0), lo_(v) {}
// NOLINTNEXTLINE(runtime/int)
inline constexpr uint128::uint128(unsigned long long v)
    : hi_(0), lo_(v) {}

#ifdef ABSL_HAVE_INTRINSIC_INT128
inline constexpr uint128::uint128(__int128 v)
    : hi_(static_cast<uint64>(static_cast<unsigned __int128>(v) >> 64)),
      lo_(static_cast<uint64>(v & ~uint64{0})) {}
inline constexpr uint128::uint128(unsigned __int128 v)
    : hi_(static_cast<uint64>(v >> 64)),
      lo_(static_cast<uint64>(v & ~uint64{0})) {}
#endif  // ABSL_HAVE_INTRINSIC_INT128

#else  // byte order
#error "Unsupported byte order: must be little-endian or big-endian."
#endif  // byte order

// Conversion operators to integer types.

inline constexpr uint128::operator bool() const {
  return lo_ || hi_;
}

inline constexpr uint128::operator char() const {
  return static_cast<char>(lo_);
}

inline constexpr uint128::operator signed char() const {
  return static_cast<signed char>(lo_);
}

inline constexpr uint128::operator unsigned char() const {
  return static_cast<unsigned char>(lo_);
}

inline constexpr uint128::operator char16_t() const {
  return static_cast<char16_t>(lo_);
}

inline constexpr uint128::operator char32_t() const {
  return static_cast<char32_t>(lo_);
}

inline constexpr uint128::operator wchar_t() const {
  return static_cast<wchar_t>(lo_);
}

// NOLINTNEXTLINE(runtime/int)
inline constexpr uint128::operator short() const {
  return static_cast<short>(lo_);  // NOLINT(runtime/int)
}

// NOLINTNEXTLINE(runtime/int)
inline constexpr uint128::operator unsigned short() const {
  return static_cast<unsigned short>(lo_);  // NOLINT(runtime/int)
}

inline constexpr uint128::operator int() const {
  return static_cast<int>(lo_);
}

inline constexpr uint128::operator unsigned int() const {
  return static_cast<unsigned int>(lo_);
}

// NOLINTNEXTLINE(runtime/int)
inline constexpr uint128::operator long() const {
  return static_cast<long>(lo_);  // NOLINT(runtime/int)
}

// NOLINTNEXTLINE(runtime/int)
inline constexpr uint128::operator unsigned long() const {
  return static_cast<unsigned long>(lo_);  // NOLINT(runtime/int)
}

// NOLINTNEXTLINE(runtime/int)
inline constexpr uint128::operator long long() const {
  return static_cast<long long>(lo_);  // NOLINT(runtime/int)
}

// NOLINTNEXTLINE(runtime/int)
inline constexpr uint128::operator unsigned long long() const {
  return static_cast<unsigned long long>(lo_);  // NOLINT(runtime/int)
}

#ifdef ABSL_HAVE_INTRINSIC_INT128
inline constexpr uint128::operator __int128() const {
  return (static_cast<__int128>(hi_) << 64) + lo_;
}

inline constexpr uint128::operator unsigned __int128() const {
  return (static_cast<unsigned __int128>(hi_) << 64) + lo_;
}
#endif  // ABSL_HAVE_INTRINSIC_INT128

// Conversion operators to floating point types.

inline uint128::operator float() const {
  return static_cast<float>(lo_) + std::ldexp(static_cast<float>(hi_), 64);
}

inline uint128::operator double() const {
  return static_cast<double>(lo_) + std::ldexp(static_cast<double>(hi_), 64);
}

inline uint128::operator long double() const {
  return static_cast<long double>(lo_) +
         std::ldexp(static_cast<long double>(hi_), 64);
}

// Comparison operators.

inline bool operator==(const uint128& lhs, const uint128& rhs) {
  return (Uint128Low64(lhs) == Uint128Low64(rhs) &&
          Uint128High64(lhs) == Uint128High64(rhs));
}

inline bool operator!=(const uint128& lhs, const uint128& rhs) {
  return !(lhs == rhs);
}

inline bool operator<(const uint128& lhs, const uint128& rhs) {
  return (Uint128High64(lhs) == Uint128High64(rhs))
             ? (Uint128Low64(lhs) < Uint128Low64(rhs))
             : (Uint128High64(lhs) < Uint128High64(rhs));
}

inline bool operator>(const uint128& lhs, const uint128& rhs) {
  return (Uint128High64(lhs) == Uint128High64(rhs))
             ? (Uint128Low64(lhs) > Uint128Low64(rhs))
             : (Uint128High64(lhs) > Uint128High64(rhs));
}

inline bool operator<=(const uint128& lhs, const uint128& rhs) {
  return (Uint128High64(lhs) == Uint128High64(rhs))
             ? (Uint128Low64(lhs) <= Uint128Low64(rhs))
             : (Uint128High64(lhs) <= Uint128High64(rhs));
}

inline bool operator>=(const uint128& lhs, const uint128& rhs) {
  return (Uint128High64(lhs) == Uint128High64(rhs))
             ? (Uint128Low64(lhs) >= Uint128Low64(rhs))
             : (Uint128High64(lhs) >= Uint128High64(rhs));
}

// Unary operators.

inline uint128 operator-(const uint128& val) {
  const uint64 hi_flip = ~Uint128High64(val);
  const uint64 lo_flip = ~Uint128Low64(val);
  const uint64 lo_add = lo_flip + 1;
  if (lo_add < lo_flip) {
    return absl::MakeUint128(hi_flip + 1, lo_add);
  }
  return absl::MakeUint128(hi_flip, lo_add);
}

inline bool operator!(const uint128& val) {
  return !Uint128High64(val) && !Uint128Low64(val);
}

// Logical operators.

inline uint128 operator~(const uint128& val) {
  return absl::MakeUint128(~Uint128High64(val), ~Uint128Low64(val));
}

inline uint128 operator|(const uint128& lhs, const uint128& rhs) {
  return absl::MakeUint128(Uint128High64(lhs) | Uint128High64(rhs),
                           Uint128Low64(lhs) | Uint128Low64(rhs));
}

inline uint128 operator&(const uint128& lhs, const uint128& rhs) {
  return absl::MakeUint128(Uint128High64(lhs) & Uint128High64(rhs),
                           Uint128Low64(lhs) & Uint128Low64(rhs));
}

inline uint128 operator^(const uint128& lhs, const uint128& rhs) {
  return absl::MakeUint128(Uint128High64(lhs) ^ Uint128High64(rhs),
                           Uint128Low64(lhs) ^ Uint128Low64(rhs));
}

inline uint128& uint128::operator|=(const uint128& other) {
  hi_ |= other.hi_;
  lo_ |= other.lo_;
  return *this;
}

inline uint128& uint128::operator&=(const uint128& other) {
  hi_ &= other.hi_;
  lo_ &= other.lo_;
  return *this;
}

inline uint128& uint128::operator^=(const uint128& other) {
  hi_ ^= other.hi_;
  lo_ ^= other.lo_;
  return *this;
}

// Shift and arithmetic assign operators.

inline uint128& uint128::operator<<=(int amount) {
  // uint64 shifts of >= 64 are undefined, so we will need some special-casing.
  if (amount < 64) {
    if (amount != 0) {
      hi_ = (hi_ << amount) | (lo_ >> (64 - amount));
      lo_ = lo_ << amount;
    }
  } else if (amount < 128) {
    hi_ = lo_ << (amount - 64);
    lo_ = 0;
  } else {
    hi_ = 0;
    lo_ = 0;
  }
  return *this;
}

inline uint128& uint128::operator>>=(int amount) {
  // uint64 shifts of >= 64 are undefined, so we will need some special-casing.
  if (amount < 64) {
    if (amount != 0) {
      lo_ = (lo_ >> amount) | (hi_ << (64 - amount));
      hi_ = hi_ >> amount;
    }
  } else if (amount < 128) {
    lo_ = hi_ >> (amount - 64);
    hi_ = 0;
  } else {
    lo_ = 0;
    hi_ = 0;
  }
  return *this;
}

inline uint128& uint128::operator+=(const uint128& b) {
  hi_ += b.hi_;
  uint64 lolo = lo_ + b.lo_;
  if (lolo < lo_)
    ++hi_;
  lo_ = lolo;
  return *this;
}

inline uint128& uint128::operator-=(const uint128& b) {
  hi_ -= b.hi_;
  if (b.lo_ > lo_)
    --hi_;
  lo_ -= b.lo_;
  return *this;
}

inline uint128& uint128::operator*=(const uint128& b) {
#if defined(ABSL_HAVE_INTRINSIC_INT128)
  // TODO(user) Remove once alignment issues are resolved and unsigned __int128
  // can be used for uint128 storage.
  *this =
      static_cast<unsigned __int128>(*this) * static_cast<unsigned __int128>(b);
  return *this;
#else   // ABSL_HAVE_INTRINSIC128
  uint64 a96 = hi_ >> 32;
  uint64 a64 = hi_ & 0xffffffff;
  uint64 a32 = lo_ >> 32;
  uint64 a00 = lo_ & 0xffffffff;
  uint64 b96 = b.hi_ >> 32;
  uint64 b64 = b.hi_ & 0xffffffff;
  uint64 b32 = b.lo_ >> 32;
  uint64 b00 = b.lo_ & 0xffffffff;
  // multiply [a96 .. a00] x [b96 .. b00]
  // terms higher than c96 disappear off the high side
  // terms c96 and c64 are safe to ignore carry bit
  uint64 c96 = a96 * b00 + a64 * b32 + a32 * b64 + a00 * b96;
  uint64 c64 = a64 * b00 + a32 * b32 + a00 * b64;
  this->hi_ = (c96 << 32) + c64;
  this->lo_ = 0;
  // add terms after this one at a time to capture carry
  *this += uint128(a32 * b00) << 32;
  *this += uint128(a00 * b32) << 32;
  *this += a00 * b00;
  return *this;
#endif  // ABSL_HAVE_INTRINSIC128
}

// Increment/decrement operators.

inline uint128 uint128::operator++(int) {
  uint128 tmp(*this);
  *this += 1;
  return tmp;
}

inline uint128 uint128::operator--(int) {
  uint128 tmp(*this);
  *this -= 1;
  return tmp;
}

inline uint128& uint128::operator++() {
  *this += 1;
  return *this;
}

inline uint128& uint128::operator--() {
  *this -= 1;
  return *this;
}

#endif  // S2_THIRD_PARTY_ABSL_NUMERIC_INT128_H_
