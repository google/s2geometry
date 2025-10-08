// Copyright 2025 Google LLC
// Author: smcallis@google.com (Sean McAllister)
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef S2_UTIL_MATH_EXACTFLOAT_BIGNUM_H_
#define S2_UTIL_MATH_EXACTFLOAT_BIGNUM_H_

#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <type_traits>

#include "absl/algorithm/container.h"
#include "absl/base/attributes.h"
#include "absl/container/inlined_vector.h"
#include "absl/log/absl_check.h"
#include "absl/numeric/bits.h"
#include "absl/numeric/int128.h"
#include "absl/strings/ascii.h"
#include "absl/strings/str_format.h"

namespace exactfloat_internal {

// A digit of a bignum. A contraction of "big digit" (rhymes with the latter).
using Bigit = uint64_t;

// A class to support arithmetic on large, arbitrary precision integers.
//
// Large integers are represented as an array of uint64_t values.
class Bignum {
 public:
  // The most common use of ExactFloat involves evaluating a 3x3 determinant to
  // determine whether 3 points are oriented clockwise or counter-clockwise.
  //
  // The typical number of mantissa bits in the result is probably about 170, so
  // we allocate 4 bigits (256 bits) inline.
  using BigitVector = absl::InlinedVector<Bigit, 4>;

  static constexpr int kBigitBits = std::numeric_limits<Bigit>::digits;

  Bignum() = default;

  // Constructs a bignum from an integral value (signed or unsigned).
  template <typename T,
            typename = std::enable_if_t<std::numeric_limits<T>::is_integer>>
  explicit Bignum(T value);

  //--------------------------------------
  // String formatting and parsing.
  //--------------------------------------

  // Constructs a bignum from an ASCII string containing decimal digits.
  //
  // The input string must only have an optional leading +/- and decimal digits.
  // Any other characters will yield std::nullopt.
  static std::optional<Bignum> FromString(absl::string_view s);

  // Formats the bignum as a decimal integer into an abseil sink.
  template <typename Sink>
  friend void AbslStringify(Sink& sink, const Bignum& b);

  friend std::ostream& operator<<(std::ostream& os, const Bignum& b) {
    return os << absl::StrFormat("%v", b);
  }

  friend std::ostream& operator<<(std::ostream& os,
                                  const std::optional<Bignum>& b) {
    if (!b) {
      return os << "[nullopt]";
    }
    return os << *b;
  }

  //--------------------------------------
  // Casting and coercion.
  //--------------------------------------

  // Returns true if bignum can be stored in T without truncation or overflow.
  template <typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
  bool FitsIn() const;

  // Cast to an integral type T unconditionally. Use FitsIn<T>() or
  // ConvertTo<T> to perform conversion with bounds checking.
  template <typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
  T Cast() const;

  // Casts the value to the given type if it fits, otherwise std::nullopt.
  template <typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
  std::optional<T> ConvertTo() const {
    if (!FitsIn<T>()) {
      return std::nullopt;
    }
    return Cast<T>();
  }

  //--------------------------------------
  // General accessors.
  //--------------------------------------

  // Returns the number of bits required for the magnitude of the value.
  //
  // Named to match std::bit_width.
  friend int bit_width(const Bignum& a);

  // Returns the number of consecutive 0 bits in the value, starting from the
  // least significant bit.
  //
  // Named to match std::countr_zero.
  friend int countr_zero(const Bignum& a);

  // Returns true if the n-th bit of the number's magnitude is 1.
  bool is_bit_set(int nbit) const;

  // Clears this bignum and sets it to zero.
  Bignum& set_zero() {
    negative_ = false;
    bigits_.clear();
    return *this;
  }

  // Sets the negative flag on the value. If the value is zero, has no effect.
  Bignum& set_negative(bool negative = true) {
    negative_ = !bigits_.empty() && negative;
    return *this;
  }

  // Returns true if the number is zero.
  bool is_zero() const { return bigits_.empty(); }

  // Returns true if the number is less than zero.
  bool is_negative() const { return negative_; }

  // Returns true if the number is odd (least significant bit is 1).
  bool is_odd() const { return is_bit_set(0); }

  // Returns true if the number is even (least significant bit is 0).
  bool is_even() const { return !is_odd(); }

  //--------------------------------------
  // Comparisons.
  //--------------------------------------

  bool operator==(const Bignum& b) const {
    return negative_ == b.negative_ && bigits_ == b.bigits_;
  }

  bool operator!=(const Bignum& b) const { return !(*this == b); }
  bool operator<(const Bignum& b) const { return Compare(b) < 0; }
  bool operator<=(const Bignum& b) const { return Compare(b) <= 0; }
  bool operator>(const Bignum& b) const { return Compare(b) > 0; }
  bool operator>=(const Bignum& b) const { return Compare(b) >= 0; }

  //--------------------------------------
  // Arithmetic operators.
  //--------------------------------------

  Bignum operator+() const { return *this; }
  Bignum operator-() const;

  // Raise this value to the given power, which must be non-negative.
  Bignum Pow(int32_t pow) const;

  Bignum& operator+=(const Bignum& b);
  Bignum& operator-=(const Bignum& b);
  Bignum& operator*=(const Bignum& b);
  Bignum& operator<<=(int nbit);
  Bignum& operator>>=(int nbit);

  friend Bignum operator+(Bignum a, const Bignum& b) { return a += b; }
  friend Bignum operator-(Bignum a, const Bignum& b) { return a -= b; }
  friend Bignum operator*(Bignum a, const Bignum& b) { return a *= b; }
  friend Bignum operator<<(Bignum a, int nbit) { return a <<= nbit; }
  friend Bignum operator>>(Bignum a, int nbit) { return a >>= nbit; }

  // Negates this bignum in place.
  void negate() {
    negative_ = !negative_;
    if (bigits_.empty()) {
      negative_ = false;
    }
  }

  // Compares to another bignum, returning -1, 0, +1.
  int Compare(const Bignum& b) const;

 private:
  // Constructs a Bignum from bigits and an optional sign bit.
  explicit Bignum(BigitVector bigits, bool negative = false)
      : bigits_(std::move(bigits)), negative_(negative) {
    Normalize();
  }

  // Drop leading zero bigits, and ensure sign is positive if result is zero.
  void Normalize() {
    while (!bigits_.empty() && bigits_.back() == 0) {
      bigits_.pop_back();
    }

    if (bigits_.empty()) {
      negative_ = false;
    }
  }

  // We store bignums in sign-magnitude form. bigits_ contains the individual
  // 64-bit digits of the bignum. If bigits_ is non-empty, then the last element
  // must be non-zero and when it is empty (representing a zero value),
  // negative_ must be false.
  BigitVector bigits_;
  bool negative_ = false;
};

////////////////////////////////////////////////////////////////////////////////
//                           Implementation Details
////////////////////////////////////////////////////////////////////////////////

// Constructs a bignum from an integral value (signed or unsigned).
template <typename T, typename>
Bignum::Bignum(T value) {
  using UT = std::make_unsigned_t<T>;

  if (value == 0) {
    return;
  }

  negative_ = false;
  if constexpr (std::is_signed_v<T>) {
    // Put into constexpr if to avoid warnings when T is unsigned.
    negative_ = (value < 0);
  }

  // Get magnitude of value, handle minimum value of T cleanly.
  UT mag = static_cast<UT>(value);
  if constexpr (std::is_signed_v<T>) {
    if (value < 0) {
      mag = UT(0) - mag;
    }
  }

  // Pack the magnitude into bigits.
  if constexpr (std::numeric_limits<UT>::digits <= kBigitBits) {
    bigits_.push_back(static_cast<Bigit>(mag));
  } else {
    while (mag) {
      bigits_.push_back(static_cast<Bigit>(mag));
      mag >>= kBigitBits;
    }
  }
}

// Formats the bignum as a decimal integer into an abseil sink.
template <typename Sink>
void AbslStringify(Sink& sink, const Bignum& b) {
  if (b.is_zero()) {
    sink.Append("0");
    return;
  }

  // Sign
  if (b.is_negative()) {
    sink.Append("-");
  }

  // Work on a copy of the magnitude.
  Bignum copy = b;
  copy.negative_ = false;

  // 10**19 is the largest power of 10 that fits in 64-bits. So we can
  // repeatedly divide and modulo the bignum to get uint64_t values we can
  // format as 19 decimal digits.
  static_assert(sizeof(Bigit) * 8 == 64);
  static constexpr uint64_t kChunkDivisor = 10'000'000'000'000'000'000u;
  Bignum::BigitVector chunks;

  while (!copy.is_zero()) {
    absl::uint128 rem = 0;
    for (int i = static_cast<int>(copy.bigits_.size()) - 1; i >= 0; --i) {
      absl::uint128 acc = (rem << 64) + copy.bigits_[i];
      Bigit quot = static_cast<Bigit>(acc / kChunkDivisor);
      rem = acc - absl::uint128(quot) * kChunkDivisor;
      copy.bigits_[i] = quot;
    }

    copy.Normalize();
    chunks.push_back(static_cast<Bigit>(rem));
  }
  ABSL_DCHECK(!chunks.empty());

  // Emit most significant chunk without zero padding.
  absl::Format(&sink, "%d", chunks.back());

  // Emit remaining chunks as fixed-width 19-digit zero-padded blocks.
  for (int i = static_cast<int>(chunks.size()) - 2; i >= 0; --i) {
    absl::Format(&sink, "%019d", chunks[i]);
  }
}

template <typename T, typename>
inline bool Bignum::FitsIn() const {
  using UT = std::make_unsigned_t<T>;

  if (is_zero()) {
    return true;
  }

  // Maximum number of bits that could fit in the output type.
  constexpr int kTBitWidth = std::numeric_limits<UT>::digits;
  constexpr int kMaxBigits = (kTBitWidth + (kBigitBits - 1)) / kBigitBits;

  // Fast reject if the bignum couldn't conceivably fit.
  if (bigits_.size() > kMaxBigits) {
    return false;
  }

  // Unsigned type T can hold the value iff the value is non-negative and
  // the bitwidth is <= the maximum bit width of the type.
  if constexpr (!std::is_signed_v<T>) {
    if (is_negative()) {
      return false;
    }
    return bit_width(*this) <= kTBitWidth;
  }

  // T is signed and our bignum isn't zero.
  ABSL_DCHECK(std::is_signed_v<T> && !is_zero());

  if (is_negative()) {
    // Magnitude must fit in negative value. If the value is negative and
    // the same bit width as the output type, the only valid value is
    // -2^(k-1).
    if (bit_width(*this) == kTBitWidth) {
      return countr_zero(*this) == kTBitWidth - 1;
    }
    return bit_width(*this) < kTBitWidth;
  } else /* positive */ {
    return bit_width(*this) <= (kTBitWidth - 1);
  }
}

template <typename T, typename>
T Bignum::Cast() const {
  using UT = std::make_unsigned_t<T>;

  constexpr int kTBitWidth = std::numeric_limits<UT>::digits;

  if (bigits_.empty()) {
    return 0;
  }

  // T fits in a Bigit, so just cast to truncate.
  UT residue = 0;
  if (kTBitWidth <= kBigitBits) {
    residue = static_cast<UT>(bigits_[0]);
  } else {
    ABSL_DCHECK_EQ(kTBitWidth % kBigitBits, 0);
    for (int i = 0; i < kTBitWidth / kBigitBits; ++i) {
      UT chunk = static_cast<UT>(bigits_[i]);
      residue |= chunk << (i * kBigitBits);
    }
  }

  // Compute two's complement of the residue if value is negative.
  if (is_negative()) {
    residue = UT(0) - residue;
  }

  return static_cast<T>(residue);
}

}  // namespace exactfloat_internal

#endif  // S2_UTIL_MATH_EXACTFLOAT_BIGNUM_H_
