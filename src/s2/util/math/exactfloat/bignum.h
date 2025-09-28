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

// A class to support arithmetic on large, arbitrary precision integers.
//
// Large integers are represented as an array of uint64_t values.
class Bignum {
 public:
  // Wrap uint64_t in a struct so we can make value-initialization a noop.
  //
  // Avoiding value-initialization overhead saves us 50% on some benchmarks.
  struct Bigit {
    static constexpr int kBits = std::numeric_limits<uint64_t>::digits;

    Bigit() {}
    constexpr Bigit(uint64_t value) : value_(value) {}
    explicit Bigit(absl::uint128 value) : value_(absl::Uint128Low64(value)) {}

    constexpr operator uint64_t() const { return value_; }
    constexpr Bigit& operator=(uint64_t value) {
      value_ = value;
      return *this;
    }

    ABSL_ATTRIBUTE_ALWAYS_INLINE constexpr Bigit& operator--(int) {
      value_--;
      return *this;
    }

    ABSL_ATTRIBUTE_ALWAYS_INLINE constexpr friend Bigit operator*(  //
        int a, Bigit b) {
      return a * b.value_;
    }

    ABSL_ATTRIBUTE_ALWAYS_INLINE friend absl::uint128 operator*(  //
        absl::uint128 a, Bigit b) {
      return a * b.value_;
    }

    ABSL_ATTRIBUTE_ALWAYS_INLINE friend absl::uint128 operator+(  //
        absl::uint128 a, Bigit b) {
      return a + b.value_;
    }

    uint64_t value_;
  };

  using BigitVector = absl::InlinedVector<Bigit, 2>;

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

  friend std::ostream& operator<<(  //
      std::ostream& os, const std::optional<Bignum>& b) {
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
  friend int bit_width(const Bignum& a);

  // Returns the number of consecutive 0 bits in the value, starting from the
  // least significant bit.
  friend int countr_zero(const Bignum& a);

  // Returns true if the n-th bit of the number's magnitude is set.
  bool Bit(int nbit) const;

  // Clears this bignum and sets it to zero.
  Bignum& SetZero() {
    sign_ = 0;
    bigits_.clear();
    return *this;
  }

  // Unconditionally makes the sign of this bignum negative.
  Bignum& SetNegative() {
    sign_ = -1;
    return *this;
  }

  // Unconditionally makes the sign of this bignum positive.
  Bignum& SetPositive() {
    sign_ = +1;
    return *this;
  }

  // Unconditionally set the sign of this bignum to match the sign of the
  // argument. If the argument is zero, set the bignum to zero.
  Bignum& SetSign(int sign) {
    if (sign == 0) {
      return SetZero();
    }

    if (sign < 0) {
      return SetNegative();
    }
    return SetPositive();
  }

  // Returns true if the number is zero.
  bool is_zero() const {  //
    return sign_ == 0;
  }

  // Returns true if the number is greater than zero.
  bool positive() const {  //
    return sign_ > 0;
  }

  // Returns true if the number is less than zero.
  bool negative() const {  //
    return sign_ < 0;
  }

  // Returns true if the number is odd (least significant bit is 1).
  bool is_odd() const { return Bit(0); }

  // Returns true if the number is even (least significant bit is 0).
  bool is_even() const { return !is_odd(); }

  //--------------------------------------
  // Comparisons.
  //--------------------------------------

  bool operator==(const Bignum& b) const {
    return sign_ == b.sign_ && bigits_ == b.bigits_;
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

 private:
  // Constructs a Bignum from bigits and an optional sign bit.
  explicit Bignum(BigitVector bigits, int sign = +1)
      : bigits_(std::move(bigits)) {
    NormalizeSign(sign);
  }

  // Returns the number of bigits in this bignum.
  size_t size() const {  //
    return bigits_.size();
  }

  // Returns true if this value has no digits.
  bool empty() const {  //
    return bigits_.empty();
  }

  // Compare to another bignum, returns -1, 0, +1.
  int Compare(const Bignum& b) const {
    if (sign_ != b.sign_) {
      return sign_ < b.sign_ ? -1 : 1;
    }

    // Signs are equal, are they both zero?
    if (sign_ == 0) {
      return 0;
    }

    // Signs are equal and non-zero, compare magnitude.
    return positive() ? CmpAbs(b) : -CmpAbs(b);
  }

  // Multiplies two unsigned bigit vectors together using Karatsuba's algorithm.
  static BigitVector KaratsubaMul(  //
      absl::Span<const Bigit> a, absl::Span<const Bigit> b);

  // Drop leading zero bigits.
  void Normalize() {
    while (!empty() && bigits_.back() == 0) {
      bigits_.pop_back();
    }

    if (empty()) {
      sign_ = 0;
    }
  }

  // Drop leading zero bigits and canonicalize sign.
  void NormalizeSign(int sign) {
    Normalize();
    sign_ = empty() ? 0 : sign;
  }

  // Returns true if the bignum is in normal form (no extra leading zeros).
  bool Normalized() const {  //
    return bigits_.empty() || bigits_.back() != 0;
  }

  // Returns true if the bignum magnitude is the given power of two.
  bool IsPow2(int pow2) const {
    const int bigits = pow2 / Bigit::kBits;
    if (bigits_.size() != bigits + 1) {
      return false;
    }

    // Verify lower words are zero.
    for (int i = 0; i < bigits; ++i) {
      if (bigits_[i] != 0) {
        return false;
      }
    }

    // Check final word is power of two.
    pow2 -= bigits * Bigit::kBits;
    ABSL_DCHECK_LT(pow2, Bigit::kBits);
    return bigits_.back() == (Bigit(1) << pow2);
  }

  // Compares magnitude with another bignum, returning -1, 0, or +1.
  //
  // Magnitudes are compared lexicographically from the most significant bigit
  // (bigits_.back()) to the least significant (bigits_[0]).
  int CmpAbs(const Bignum& b) const;

  BigitVector bigits_;
  int sign_ = 0;
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

  sign_ = +1;
  if constexpr (std::is_signed_v<T>) {
    sign_ = (value < 0) ? -1 : +1;
  }

  // Get magnitude of value, handle minimum value of T cleanly.
  UT mag = static_cast<UT>(value);
  if constexpr (std::is_signed_v<T>) {
    if (value < 0) {
      mag = UT(0) - mag;
    }
  }

  // Pack the magnitude into bigits.
  if constexpr (std::numeric_limits<UT>::digits <= Bigit::kBits) {
    bigits_.push_back(static_cast<Bigit>(mag));
  } else {
    while (mag) {
      bigits_.push_back(static_cast<Bigit>(mag));
      mag >>= Bigit::kBits;
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
  if (b.negative()) {
    sink.Append("-");
  }

  // Work on a copy of the magnitude.
  Bignum copy = b;
  copy.sign_ = 1;

  // Repeatedly divide and modulo by 10^19 to get decimal chunks.
  static constexpr uint64_t kBase = 10'000'000'000'000'000'000u;
  Bignum::BigitVector chunks;

  while (!copy.is_zero()) {
    absl::uint128 rem = 0;
    for (int i = static_cast<int>(copy.bigits_.size()) - 1; i >= 0; --i) {
      absl::uint128 acc = (rem << 64) + copy.bigits_[i];
      Bignum::Bigit quot = static_cast<Bignum::Bigit>(acc / kBase);
      rem = acc - absl::uint128(quot) * kBase;
      copy.bigits_[i] = quot;
    }

    copy.Normalize();
    chunks.push_back(static_cast<Bignum::Bigit>(rem));
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

  if (sign_ == 0) {
    return true;
  }

  // Maximum number of bits that could fit in the output type.
  constexpr int kTBitWidth = std::numeric_limits<UT>::digits;
  constexpr int kMaxBigits = (kTBitWidth + (Bigit::kBits - 1)) / Bigit::kBits;

  // Fast reject if the bignum couldn't conceivably fit.
  if (bigits_.size() > kMaxBigits) {
    return false;
  }

  // Unsigned type T can hold the value iff the value is non-negative and
  // the bitwidth is <= the maximum bit width of the type.
  if constexpr (!std::is_signed_v<T>) {
    if (negative()) {
      return false;
    }
    return bit_width(*this) <= kTBitWidth;
  }

  // T is signed and our bignum isn't zero.
  ABSL_DCHECK(std::is_signed_v<T> && !is_zero());

  if (positive()) {
    return bit_width(*this) <= (kTBitWidth - 1);
  } else /* negative() */ {
    // Magnitude must fit in negative value. If the value is negative and
    // the same bit width as the output type, the only valid value is
    // -2^(k-1).
    if (bit_width(*this) == kTBitWidth) {
      return IsPow2(kTBitWidth - 1);
    }
    return bit_width(*this) < kTBitWidth;
  }
}

template <typename T, typename>
T Bignum::Cast() const {
  using UT = std::make_unsigned_t<T>;

  constexpr int kTBitWidth = std::numeric_limits<UT>::digits;

  if (empty()) {
    return 0;
  }

  // Grab the bottom bits into an unsigned value.
  UT residue = 0;
  for (size_t i = 0; i < bigits_.size(); ++i) {
    const int shift = i * Bigit::kBits;
    if (shift >= kTBitWidth) {
      break;
    }

    const int room = kTBitWidth - shift;
    UT chunk = static_cast<UT>(bigits_[i]);
    if (room < Bigit::kBits && room < std::numeric_limits<UT>::digits) {
      chunk &= (UT(1) << room) - UT(1);
    }
    residue |= (chunk << shift);
  }

  // Compute two's complement of the residue if value is negative.
  if (negative()) {
    residue = UT(0) - residue;
  }

  return static_cast<T>(residue);
}

inline int Bignum::CmpAbs(const Bignum& b) const {
  if (size() != b.size()) {
    return size() < b.size() ? -1 : +1;
  }

  for (int i = size() - 1; i >= 0; --i) {
    if (bigits_[i] != b.bigits_[i]) {
      return bigits_[i] < b.bigits_[i] ? -1 : +1;
    }
  }

  return 0;
}

}  // namespace exactfloat_internal
