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

#include <cstdint>
#include <limits>
#include <optional>
#include <type_traits>

#include "absl/container/inlined_vector.h"
#include "absl/log/absl_check.h"
#include "absl/numeric/bits.h"
#include "absl/numeric/int128.h"
#include "absl/strings/ascii.h"
#include "absl/strings/str_format.h"

namespace internal {

// Most of the STL cannot be overloaded per the spec, so we need to roll our own
// wrappers that will also work with absl::int128.

template <typename T>
constexpr bool IsInt = std::numeric_limits<T>::is_integer;

template <typename T>
constexpr bool IsSigned() {
  if constexpr (std::is_same_v<T, absl::uint128>) {
    return false;
  }

  if constexpr (std::is_same_v<T, absl::int128>) {
    return true;
  }

  return std::is_signed_v<T>;
}

template <typename T>
constexpr auto InferUnsigned() {
  if constexpr (std::is_same_v<T, absl::int128> ||
                std::is_same_v<T, absl::uint128>) {
    return absl::uint128{};
  } else {
    return std::make_unsigned_t<T>{};
  }
}

template <typename T>
using MakeUnsigned = decltype(InferUnsigned<T>());

}  // namespace internal

class Bignum {
 private:
  using Bigit = uint64_t;

  static constexpr int kKaratsubaThreshold = 32;
  static constexpr int kBigitBits = std::numeric_limits<Bigit>::digits;

 public:
  Bignum() = default;

  // Constructs a bignum from an integral value (signed or unsigned).
  template <typename T, typename = std::enable_if_t<internal::IsInt<T>>>
  explicit Bignum(T value) {
    using UT = internal::MakeUnsigned<T>;

    if (value == 0) {
      return;
    }

    sign_ = +1;
    if constexpr (internal::IsSigned<T>()) {
      sign_ = (value < 0) ? -1 : +1;
    }

    // Get magnitude of value, handle minimum value of T cleanly.
    UT mag = static_cast<UT>(value);
    if constexpr (internal::IsSigned<T>()) {
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

  // Constructs a bignum from an ASCII string containing decimal digits.
  //
  // The input string must only have an optional leading +/- and decimal digits.
  // Any other characters will yield std::nullopt.
  static std::optional<Bignum> FromString(absl::string_view s) {
    // We can fit ~10^19 into a uint64_t.
    constexpr int kMaxChunkDigits = 19;

    // NOTE: We use a simple multiply-and-add (aka Horner's) method here for the
    // sake of simplicity. This isn't the fastest algorithm, being quadratic in
    // the number of chunks the input has. If we use divide and conquer approach
    // or an FFT based multiply we could probably make this ~O(n^1.5) or
    // semi-linear.

    // Precomputed powers of 10.
    static constexpr uint64_t kPow10[20] = {1ull,
                                            10ull,
                                            100ull,
                                            1000ull,
                                            10000ull,
                                            100000ull,
                                            1000000ull,
                                            10000000ull,
                                            100000000ull,
                                            1000000000ull,
                                            10000000000ull,
                                            100000000000ull,
                                            1000000000000ull,
                                            10000000000000ull,
                                            100000000000000ull,
                                            1000000000000000ull,
                                            10000000000000000ull,
                                            100000000000000000ull,
                                            1000000000000000000ull,
                                            10000000000000000000ull};

    Bignum out;
    if (s.empty()) {
      return out;
    }

    // Reserve space for bigits.
    out.bigits_.reserve((s.size() + kMaxChunkDigits - 1) / kMaxChunkDigits);

    int sign = +1;
    uint64_t chunk = 0;
    int clen = 0;

    // Finish processing the current chunk.
    auto FlushChunk = [&]() {
      if (clen) {
        out.MulAddSmall(kPow10[clen], chunk);
        chunk = 0;
        clen = 0;
      }
    };

    // Consume optional +/- at the front.
    int start = 0;
    if ((s[0] == '+' || s[0] == '-')) {
      sign = (s[0] == '-') ? -1 : +1;
      ++start;
    }

    bool seen_digit = false;
    for (char c : s.substr(start)) {
      if (!absl::ascii_isdigit(c)) {
        return std::nullopt;
      }

      // Accumulate digit into the local 64-bit chunk.  Skip leading zeros.
      uint64_t digit = static_cast<uint64_t>(c - '0');
      if (!seen_digit && digit == 0) {
        continue;
      }
      seen_digit = true;

      chunk = 10 * chunk + digit;
      ++clen;

      if (clen == kMaxChunkDigits) {
        FlushChunk();
      }
    }
    FlushChunk();

    out.NormalizeSign(sign);
    return out;
  }

  // Formats the bignum as a decimal integer into an abseil sink.
  template <typename Sink>
  friend void AbslStringify(Sink& sink, const Bignum& b) {
    if (b.zero()) {
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
    static constexpr uint64_t kBase = 10000000000000000000ull;
    absl::InlinedVector<uint64_t, 4> chunks;

    while (!copy.zero()) {
      absl::uint128 rem = 0;
      for (int i = static_cast<int>(copy.bigits_.size()) - 1; i >= 0; --i) {
        absl::uint128 acc = (rem << 64) + copy.bigits_[i];
        uint64_t quot = static_cast<uint64_t>(acc / kBase);
        rem = acc - absl::uint128(quot) * kBase;
        copy.bigits_[i] = quot;
      }

      copy.Normalize();
      chunks.push_back(static_cast<uint64_t>(rem));
    }
    ABSL_DCHECK(!chunks.empty());

    // Emit most significant chunk without zero padding.
    absl::Format(&sink, "%d", chunks.back());

    // Emit remaining chunks as fixed-width 19-digit zero-padded blocks.
    for (int i = static_cast<int>(chunks.size()) - 2; i >= 0; --i) {
      absl::Format(&sink, "%019d", chunks[i]);
    }
  }

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

  // Returns true if bignum can be stored in T without truncation or overflow.
  template <typename T, typename = std::enable_if_t<internal::IsInt<T>>>
  bool Compatible() const {
    using UT = internal::MakeUnsigned<T>;

    if (sign_ == 0) {
      return true;
    }

    // Maximum number of bits that could fit in the output type.
    constexpr int kTBitWidth = std::numeric_limits<UT>::digits;
    constexpr int kMaxBigits = (kTBitWidth + (kBigitBits - 1)) / kBigitBits;

    // Fast reject if the bignum couldn't conceivably fit.
    if (bigits_.size() > kMaxBigits) {
      return false;
    }

    // Unsigned type T can hold the value iff the value is non-negative and the
    // bitwidth is <= the maximum bit width of the type.
    if constexpr (!internal::IsSigned<T>()) {
      if (negative()) {
        return false;
      }
      return BitWidth() <= kTBitWidth;
    }

    // T is signed and our bignum isn't zero.
    ABSL_DCHECK(internal::IsSigned<T>() && !zero());

    if (positive()) {
      return BitWidth() <= (kTBitWidth - 1);
    } else /* negative() */ {
      // Magnitude must fit in negative value. If the value is negative and the
      // same bit width as the output type, the only valid value is -2^(k-1).
      if (BitWidth() == kTBitWidth) {
        return IsPow2(kTBitWidth - 1);
      }
      return BitWidth() < kTBitWidth;
    }
  }

  // Cast to an integral type T unconditionally. Use Compatible<T>() or
  // Convert<T> to perform conversion with bounds checking.
  template <typename T, typename = std::enable_if_t<internal::IsInt<T>>>
  T Cast() const {
    using UT = internal::MakeUnsigned<T>;

    constexpr int kTBitWidth = std::numeric_limits<UT>::digits;

    if (empty()) {
      return 0;
    }

    // Grab the bottom bits into an unsigned value.
    UT residue = 0;
    for (size_t i = 0; i < bigits_.size(); ++i) {
      const int shift = i * kBigitBits;
      if (shift >= kTBitWidth) {
        break;
      }

      const int room = kTBitWidth - shift;
      UT chunk = static_cast<UT>(bigits_[i]);
      if (room < kBigitBits && room < std::numeric_limits<UT>::digits) {
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

  // Casts the value to the given type if it fits, otherwise std::nullopt.
  template <typename T, typename = std::enable_if_t<internal::IsInt<T>>>
  std::optional<T> Convert() const {
    if (!Compatible<T>()) {
      return std::nullopt;
    }
    return Cast<T>();
  }

  // // Creates a Bignum by parsing a decimal representation from a string.
  // explicit Bignum(absl::string_view dec) { ParseDecimal_(dec); }

  // Returns the number of bits required to represent the bignum.
  int BitWidth() const {
    ABSL_DCHECK(Normalized());
    if (empty()) {
      return 0;
    }

    // Bit width is the bits in the least significant bigits + bit width of the
    // most significant word.
    const int msw_width = (kBigitBits - absl::countl_zero(bigits_.back()));
    const int lsw_width = (bigits_.size() - 1) * kBigitBits;
    return msw_width + lsw_width;
  }

  // Returns the number of consecutive 0 bits in the value, starting from the
  // least significant bit.
  int CountrZero() const {
    if (zero()) {
      return 0;
    }

    int nzero = 0;
    for (Bigit bigit : bigits_) {
      if (bigit == 0) {
        nzero += kBigitBits;
      } else {
        nzero += absl::countr_zero(bigit);
        break;
      }
    }
    return nzero;
  }

  // Returns true if the n-th bit of the number's magnitude is set.
  bool Bit(int nbit) const {
    ABSL_DCHECK_GE(nbit, 0);
    if (zero()) {
      return false;
    }

    const int digit = nbit / kBigitBits;
    const int shift = nbit % kBigitBits;

    if (digit >= size()) {
      return false;
    }

    return ((bigits_[digit] >> shift) & 0x1) != 0;
  }

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
  bool zero() const {  //
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
  bool odd() const { return Bit(0); }

  // Returns true if the number is even (least significant bit is 0).
  bool even() const { return !odd(); }

  bool operator==(const Bignum& b) const {
    return sign_ == b.sign_ && bigits_ == b.bigits_;
  }

  bool operator!=(const Bignum& b) const { return !(*this == b); }

  bool operator<(const Bignum& b) const { return Compare(b) < 0; }

  bool operator<=(const Bignum& b) const { return Compare(b) <= 0; }

  bool operator>(const Bignum& b) const { return Compare(b) > 0; }

  bool operator>=(const Bignum& b) const { return Compare(b) >= 0; }

  Bignum operator+() const { return *this; }

  Bignum operator-() const {
    Bignum result = *this;
    result.sign_ = -result.sign_;
    return result;
  }

  Bignum& operator+=(const Bignum& b) {
    if (b.zero()) {
      return *this;
    }

    if (zero()) {
      *this = b;
      return *this;
    }

    if (sign_ == b.sign_) {
      // Same sign:
      //   (+a) + (+b) == +(a + b)
      //   (-a) + (-b) == -(a + b)
      AddAbs(b);
    } else {
      if (CmpAbs(b) >= 0) {
        // |a| >= |b|, so a - b is same sign as a.
        SubAbsGe(b);
        NormalizeSign(sign_);
      } else {
        // |a| < |b|, so a - b is same sign as b.
        SubAbsLt(b);
        NormalizeSign(b.sign_);
      }
    }

    return *this;
  }

  Bignum& operator-=(const Bignum& b) {
    if (this == &b) {
      bigits_.clear();
      sign_ = 0;
      return *this;
    }

    if (b.zero()) {
      return *this;
    }

    if (zero()) {
      return *this = -b;
    }

    if (sign_ != b.sign_) {
      AddAbs(b);
    } else {
      if (CmpAbs(b) >= 0) {
        SubAbsGe(b);
        NormalizeSign(sign_);
      } else {
        SubAbsLt(b);
        NormalizeSign(-sign_);
      }
    }

    return *this;
  }

  // Left-shift the bignum by nbit.
  Bignum& operator<<=(int nbit) {
    ABSL_DCHECK_GE(nbit, 0);
    if (zero() || nbit == 0) {
      return *this;
    }

    const int nbigit = nbit / kBigitBits;
    const int nrem = nbit % kBigitBits;

    // First, handle the whole-bigit shift by inserting zeros.
    bigits_.insert(bigits_.begin(), nbigit, 0);

    // Then, handle the within-bigit shift, if any.
    if (nrem != 0) {
      Bigit carry = 0;
      for (size_t i = 0; i < bigits_.size(); ++i) {
        const Bigit old_val = bigits_[i];
        bigits_[i] = (old_val << nrem) | carry;
        carry = old_val >> (kBigitBits - nrem);
      }

      if (carry) {
        bigits_.push_back(carry);
      }
    }

    return *this;
  }

  // Right-shift the bignum by nbit.
  Bignum& operator>>=(int nbit) {
    ABSL_DCHECK_GE(nbit, 0);
    if (zero() || nbit == 0) {
      return *this;
    }

    // Shifting by more than the bit width results in zero.
    if (nbit >= BitWidth()) {
      bigits_.clear();
      sign_ = 0;
      return *this;
    }

    const int nbigit = nbit / kBigitBits;
    const int nrem = nbit % kBigitBits;

    // First, handle the whole-bigit shift by removing bigits.
    bigits_.erase(bigits_.begin(), bigits_.begin() + nbigit);

    // Then, handle the within-bigit shift, if any.
    if (nrem != 0) {
      Bigit carry = 0;
      for (int i = static_cast<int>(bigits_.size()) - 1; i >= 0; --i) {
        const Bigit old_val = bigits_[i];
        bigits_[i] = (old_val >> nrem) | carry;
        carry = old_val << (kBigitBits - nrem);
      }
    }

    // Result might be smaller or zero, so normalize.
    NormalizeSign(sign_);
    return *this;
  }

  Bignum& operator*=(const Bignum& b) {
    if (zero() || b.zero()) {
      bigits_.clear();
      sign_ = 0;
      return *this;
    }

    const int new_sign = sign_ * b.sign_;
    bigits_ = MulAbs(bigits_, b.bigits_);
    NormalizeSign(new_sign);
    return *this;
  }

  // Raise this value to the given power, which must be non-negative.
  Bignum Pow(int32_t pow) const {
    ABSL_DCHECK_GE(pow, 0);

    // Anything to the zero-th power is 1 (including zero).
    if (pow == 0) {
      return Bignum(1);
    }

    if (zero()) {
      return Bignum(0);
    }

    if (*this == Bignum(1)) {
      return Bignum(1);
    }

    if (*this == Bignum(-1)) {
      return (pow % 2 != 0) ? Bignum(-1) : Bignum(1);
    }

    // Core algorithm: Exponentiation by squaring.
    Bignum result(1);
    Bignum base = *this;  // A mutable copy of the base.
    uint32_t upow = static_cast<uint32_t>(pow);

    while (upow > 0) {
      if (upow & 1) {  // If current exponent bit is 1, multiply into result.
        result *= base;
      }
      base *= base;
      upow >>= 1;
    }

    return result;
  }

  friend Bignum operator*(Bignum a, const Bignum& b) { return a *= b; }

  friend Bignum operator+(Bignum a, const Bignum& b) { return a += b; }

  friend Bignum operator-(Bignum a, const Bignum& b) { return a -= b; }

  friend Bignum operator<<(Bignum a, int nbit) { return a <<= nbit; }

  friend Bignum operator>>(Bignum a, int nbit) { return a >>= nbit; }

 private:
  // Construct a Bignum from bigits and an optional sign bit.
  explicit Bignum(absl::Span<const Bigit> bigits, int sign = +1) {
    if (bigits.empty()) {
      return;
    }

    bigits_.assign(bigits.begin(), bigits.end());
    NormalizeSign(sign);
  }

  // Returns the number of bigits in this bignum.
  int size() const {  //
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

  // Compute value = value * mul + add where mul, add ≤ 10^19. We can accumulate
  // using a 128 bit integer in a single pass over the bigits for small terms.
  void MulAddSmall(uint64_t mul, uint64_t add) {
    absl::uint128 carry = add;
    for (size_t i = 0, n = bigits_.size(); i < n; ++i) {
      absl::uint128 prod = absl::uint128(bigits_[i]) * mul + carry;
      bigits_[i] = absl::Uint128Low64(prod);
      carry = absl::Uint128High64(prod);
    }

    if (carry != 0) {
      bigits_.push_back(absl::Uint128Low64(carry));
    }
  }

  // Multiplies two bigit operands. Uses either simple quadratic multiplication
  // (SimpleMulAbs) or divide-and-conquer multiplication (KaratsubaMulAbs) based
  // on the size of the operands.
  static absl::InlinedVector<Bignum::Bigit, 2> MulAbs(
      absl::Span<const Bigit> a, absl::Span<const Bigit> b) {
    // Fast path for single-bigit multiplication.
    if (a.size() == 1 && b.size() == 1) {
      absl::uint128 prod = absl::uint128(a[0]) * b[0];
      const uint64_t lo = absl::Uint128Low64(prod);
      const uint64_t hi = absl::Uint128High64(prod);
      if (hi == 0) {
        return {lo};
      }
      return {lo, hi};
    }

    if (a.size() < kKaratsubaThreshold || b.size() < kKaratsubaThreshold) {
      return SimpleMulAbs(a, b);
    }
    return KaratsubaMulAbs(a, b);
  }

  // Performs simple quadratic long multiplication between two sets of bigits.
  static absl::InlinedVector<Bignum::Bigit, 2> SimpleMulAbs(
      absl::Span<const Bigit> a, absl::Span<const Bigit> b) {
    if (a.empty() || b.empty()) {
      return {};
    }

    absl::InlinedVector<Bigit, 2> result(a.size() + b.size(), 0);
    for (size_t i = 0; i < a.size(); ++i) {
      if (a[i] == 0) {
        continue;
      }

      absl::uint128 carry = 0;
      for (size_t j = 0; j < b.size(); ++j) {
        absl::uint128 prod = absl::uint128(a[i]) * b[j] + result[i + j] + carry;
        result[i + j] = absl::Uint128Low64(prod);
        carry = absl::Uint128High64(prod);
      }

      // Propagate final carry. This can ripple through multiple bigits.
      for (size_t k = i + b.size(); carry != 0 && k < result.size(); ++k) {
        absl::uint128 sum = absl::uint128(result[k]) + carry;
        result[k] = absl::Uint128Low64(sum);
        carry = absl::Uint128High64(sum);
      }
      ABSL_DCHECK_EQ(carry, 0);  // Result vector must be large enough.
    }

    // Normalize vector by removing leading zeros.
    while (!result.empty() && result.back() == 0) {
      result.pop_back();
    }
    return result;
  }

  // Repeatedly divides a multiplication in half and recurses, stitching
  // results back together to get the final result.
  static absl::InlinedVector<Bignum::Bigit, 2> KaratsubaMulAbs(
      absl::Span<const Bigit> a, absl::Span<const Bigit> b) {
    // Base case is handled by MulAbs dispatcher, so we only handle recursion.
    const size_t n = std::max(a.size(), b.size());
    const size_t m = (n + 1) / 2;

    const Bignum a0(a.subspan(0, std::min(m, a.size())));
    const Bignum a1(a.size() > m ? a.subspan(m) : absl::Span<const Bigit>());
    const Bignum b0(b.subspan(0, std::min(m, b.size())));
    const Bignum b1(b.size() > m ? b.subspan(m) : absl::Span<const Bigit>());

    Bignum z2;
    z2.bigits_ = MulAbs(a1.bigits_, b1.bigits_);
    z2.NormalizeSign(+1);

    Bignum z0;
    z0.bigits_ = MulAbs(a0.bigits_, b0.bigits_);
    z0.NormalizeSign(+1);

    Bignum z1;
    z1.bigits_ = MulAbs((a0 + a1).bigits_, (b0 + b1).bigits_);
    z1.NormalizeSign(+1);

    z1 -= z2;
    z1 -= z0;

    // Recombine: result = (z2 << 2*m) + (z1 << m) + z0
    z2 <<= (2 * m * kBigitBits);
    z1 <<= (m * kBigitBits);

    Bignum result = z2 + z1 + z0;
    return result.bigits_;
  }

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

  // Returns true if the bignum is the given power of two.
  bool IsPow2(int pow2) const {
    const int bigits = pow2 / kBigitBits;
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
    pow2 -= bigits * kBigitBits;
    ABSL_DCHECK_LT(pow2, kBigitBits);
    return bigits_.back() == (Bigit(1) << pow2);
  }

  // Compares magnitude with another bignum, returning -1, 0, or +1.
  int CmpAbs(const Bignum& b) const;

  // Adds another bignum to this bignum in place.
  void AddAbs(const Bignum& b);

  // In-place subtraction: *this = |*this| - |b|, assuming |*this| >= |b|.
  void SubAbsGe(const Bignum& b);

  // In-place subtraction: *this = |b| - |*this|, assuming |*this| < |b|.
  void SubAbsLt(const Bignum& b);

  absl::InlinedVector<Bigit, 2> bigits_;
  char sign_ = 0;
};

////////////////////////////////////////////////////////////////////////////////
//                           Implementation Details
////////////////////////////////////////////////////////////////////////////////

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

inline void Bignum::AddAbs(const Bignum& b) {
  // Grow if needed.
  const bool a_longer = size() > b.size();
  const size_t min_size = std::min(size(), b.size());
  const size_t max_size = std::max(size(), b.size());
  bigits_.resize(max_size, 0);

  // Add common parts.
  absl::uint128 sum;
  absl::uint128 carry = 0;
  for (size_t i = 0; i < min_size; ++i) {
    sum = absl::uint128(bigits_[i]) + b.bigits_[i] + carry;
    bigits_[i] = absl::Uint128Low64(sum);
    carry = absl::Uint128High64(sum);
  }

  // Propagate carry through the longer operand.
  const auto* longer = a_longer ? this : &b;
  for (size_t i = min_size; i < max_size; ++i) {
    sum = absl::uint128(longer->bigits_[i]) + carry;
    bigits_[i] = absl::Uint128Low64(sum);
    carry = absl::Uint128High64(sum);
  }

  if (carry) {
    bigits_.push_back(absl::Uint128Low64(carry));
  }
}

inline void Bignum::SubAbsGe(const Bignum& b) {
  ABSL_DCHECK_GE(CmpAbs(b), 0);
  uint64_t borrow = 0;

  size_t i = 0;
  for (; i < b.size(); ++i) {
    const uint64_t d1 = bigits_[i];
    const uint64_t d2 = b.bigits_[i];
    const uint64_t diff = d1 - d2 - borrow;
    borrow = (d1 < d2) || (borrow && d1 == d2);
    bigits_[i] = diff;
  }

  for (; borrow && i < bigits_.size(); ++i) {
    borrow = (bigits_[i] == 0);
    bigits_[i]--;
  }
  ABSL_DCHECK(!borrow);
  Normalize();
}

inline void Bignum::SubAbsLt(const Bignum& b) {
  ABSL_DCHECK_LT(CmpAbs(b), 0);
  uint64_t borrow = 0;
  const size_t n_this = bigits_.size();
  const size_t n_b = b.size();
  bigits_.resize(n_b);

  size_t i = 0;
  for (; i < n_this; ++i) {
    const uint64_t d1 = b.bigits_[i];
    const uint64_t d2 = bigits_[i];
    const uint64_t diff = d1 - d2 - borrow;
    borrow = (d1 < d2) || (borrow && d1 == d2);
    bigits_[i] = diff;
  }

  for (; i < n_b; ++i) {
    const uint64_t d1 = b.bigits_[i];
    const uint64_t diff = d1 - borrow;
    borrow = (borrow && d1 == 0);
    bigits_[i] = diff;
  }
  ABSL_DCHECK(!borrow);
  Normalize();
}
