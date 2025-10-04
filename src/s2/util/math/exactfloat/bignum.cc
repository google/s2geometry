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

#include "s2/util/math/exactfloat/bignum.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <utility>
#include <vector>

#include "absl/algorithm/container.h"
#include "absl/base/nullability.h"
#include "absl/log/absl_check.h"
#include "absl/numeric/bits.h"
#include "absl/numeric/int128.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"

namespace exactfloat_internal {

// Number of bigits in the result of a multiplication before we fall back to
// simple multiplication in the Karatsuba recursion. Determined empirically.
static constexpr int kSimpleMulThreshold = 64;

// Computes out[i] = a[i]*b + c
//
// Returns the final carry, if any.
inline Bigit MulAdd(absl::Span<Bigit> out, absl::Span<const Bigit> a, Bigit b,
                    Bigit c);

// Compares magnitude magnitude of two bigit vectors, returning -1, 0, or +1.
//
// Magnitudes are compared lexicographically from the most significant bigit
// to the least significant.
int CmpAbs(absl::Span<const Bigit> a, absl::Span<const Bigit> b) {
  if (a.size() != b.size()) {
    return a.size() < b.size() ? -1 : +1;
  }

  for (int i = a.size() - 1; i >= 0; --i) {
    if (a[i] != b[i]) {
      return a[i] < b[i] ? -1 : +1;
    }
  }

  return 0;
}

int Bignum::Compare(const Bignum& b) const {
  if (is_negative() != b.is_negative()) {
    return is_negative() ? -1 : +1;
  }

  // Signs are equal, are they both zero?
  if (is_zero() && b.is_zero()) {
    return 0;
  }

  // Signs are equal and non-zero, compare magnitude.
  const int compare = CmpAbs(bigits_, b.bigits_);
  return is_negative() ? -compare : compare;
}

std::optional<Bignum> Bignum::FromString(absl::string_view s) {
  // A chunk is up to 19 decimal digits, which can always fit into a Bigit.
  constexpr int64_t kMaxChunkDigits = std::numeric_limits<Bigit>::digits10;

  // NOTE: We use a simple multiply-and-add (aka Horner's) method here for the
  // sake of simplicity. This isn't the fastest algorithm, being quadratic in
  // the number of chunks the input has. If we use divide and conquer approach
  // or an FFT based multiply we could probably make this ~O(n^1.5) or
  // semi-linear.

  // Precomputed powers of 10.
  static constexpr auto kPow10 = []() {
    std::array<Bigit, kMaxChunkDigits + 1> out = {1};
    for (size_t i = 1; i < out.size(); ++i) {
      out[i] = 10 * out[i - 1];
    }
    return out;
  }();

  Bignum out;
  if (s.empty()) {
    return out;
  }

  out.bigits_.reserve((s.size() + kMaxChunkDigits - 1) / kMaxChunkDigits);

  bool negative = false;

  // Consume optional +/- at the front.
  auto begin = s.cbegin();
  if ((*begin == '+' || *begin == '-')) {
    negative = (s[0] == '-');
    ++begin;
  }

  const auto end = s.cend();
  while (begin < end) {
    int64_t chunk_len = std::min(
        static_cast<int64_t>(std::distance(begin, end)), kMaxChunkDigits);

    Bigit chunk = 0;
    auto result = std::from_chars(begin, begin + chunk_len, chunk);
    if (result.ec != std::errc() || (result.ptr - begin) != chunk_len) {
      return std::nullopt;
    }
    begin += chunk_len;

    // Shift out up by chunk_len digits and add the chunk to it.
    auto outspan = absl::MakeSpan(out.bigits_);
    Bigit carry = MulAdd(outspan, outspan, kPow10[chunk_len], chunk);
    if (carry) {
      out.bigits_.emplace_back(carry);
    }
  }

  out.negative_ = negative;
  out.Normalize();
  return out;
}

int bit_width(const Bignum& a) {
  ABSL_DCHECK(a.is_normalized());
  if (a.bigits_.empty()) {
    return 0;
  }

  // Bit width is the bits in the least significant bigits + bit width of
  // the most significant word.
  const int msw_width =
      (Bignum::kBigitBits - absl::countl_zero(a.bigits_.back()));
  const int lsw_width = (a.bigits_.size() - 1) * Bignum::kBigitBits;
  return msw_width + lsw_width;
}

int countr_zero(const Bignum& a) {
  if (a.is_zero()) {
    return 0;
  }

  int nzero = 0;
  for (Bigit bigit : a.bigits_) {
    if (bigit == 0) {
      nzero += Bignum::kBigitBits;
    } else {
      nzero += absl::countr_zero(bigit);
      break;
    }
  }
  return nzero;
}

bool Bignum::is_bit_set(int nbit) const {
  ABSL_DCHECK_GE(nbit, 0);
  if (is_zero()) {
    return false;
  }

  const size_t digit = nbit / kBigitBits;
  const size_t shift = nbit % kBigitBits;

  if (digit >= bigits_.size()) {
    return false;
  }

  return ((bigits_[digit] >> shift) & 0x1) != 0;
}

Bignum Bignum::operator-() const {
  Bignum result = *this;
  result.set_negative(!result.negative_);
  return result;
}

Bignum& Bignum::operator<<=(int nbit) {
  ABSL_DCHECK_GE(nbit, 0);
  if (is_zero() || nbit == 0) {
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

Bignum& Bignum::operator>>=(int nbit) {
  ABSL_DCHECK_GE(nbit, 0);
  if (is_zero() || nbit == 0) {
    return *this;
  }

  // Shifting by more than the bit width results in zero.
  if (nbit >= bit_width(*this)) {
    return set_zero();
  }

  const int nbigit = nbit / kBigitBits;
  const int nrem = nbit % kBigitBits;

  // First, handle the whole-bigit shift by removing bigits.
  bigits_.erase(bigits_.begin(), bigits_.begin() + nbigit);

  // Then, handle the within-bigit shift, if any.
  if (nrem != 0) {
    Bigit carry = 0;
    for (auto i = static_cast<int64_t>(bigits_.size()) - 1; i >= 0; --i) {
      const Bigit old_val = bigits_[i];
      bigits_[i] = (old_val >> nrem) | carry;
      carry = old_val << (kBigitBits - nrem);
    }
  }

  // Result might be smaller or zero, so normalize.
  Normalize();
  return *this;
}

// Raise this value to the given power, which must be non-negative.
Bignum Bignum::Pow(int32_t pow) const {
  ABSL_DCHECK_GE(pow, 0);

  // Anything to the zero-th power is 1 (including zero).
  if (pow == 0) {
    return Bignum(1);
  }

  if (is_zero()) {
    return Bignum(0);
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

// Computes a + b + carry and updates the carry.
inline Bigit AddCarry(Bigit a, Bigit b, Bigit* absl_nonnull carry) {
  auto sum = absl::uint128(a) + b + *carry;
  *carry = absl::Uint128High64(sum);
  return static_cast<Bigit>(sum);
}

// Computes a - b - borrow  and updates the borrow.
//
// NOTE: Borrow must be one or zero.
inline Bigit SubBorrow(Bigit a, Bigit b, Bigit* absl_nonnull borrow) {
  ABSL_DCHECK_LE(*borrow, Bigit(1));
  Bigit diff = a - b - *borrow;
  *borrow = (a < b) || (*borrow && (a == b));
  return diff;
}

// Computes a * b + carry and updates the carry.
inline Bigit MulCarry(Bigit a, Bigit b, Bigit* absl_nonnull carry) {
  auto sum = absl::uint128(a) * b + *carry;
  *carry = absl::Uint128High64(sum);
  return static_cast<Bigit>(sum);
}

// Computes out += a * b + carry and updates the carry.
//
// NOTE: Will not overflow even if a, b, and c are their maximum values.
inline void MulAddCarry(Bigit& out, Bigit a, Bigit b,
                        Bigit* absl_nonnull carry) {
  auto sum = absl::uint128(a) * b + *carry + out;
  *carry = absl::Uint128High64(sum);
  out = static_cast<Bigit>(sum);
}

// Computes a += b in place. Returns the final carry (if any).
//
// A operand must be at least as large as B. When adding two same-sized values,
// the result may overflow and be larger than either of them, in which case we
// will return the final carry value.
//
// This allows a work flow like this:
//    Bigit carry = AddInPlace(a, b);
//    if (carry) {
//      a.bigits_.emplace_back(carry);
//    }
//
// Rather than having to expand A to B.bigits_.size() + 1, and popping off the
// top bigit if it's unused (which is the most common case).
inline Bigit AddInPlace(absl::Span<Bigit> a, absl::Span<const Bigit> b) {
  ABSL_DCHECK_GE(a.size(), b.size());

  Bigit carry = 0;

  // Dispatch four at a time to help loop unrolling.
  size_t i = 0;
  while (i + 4 <= b.size()) {
    for (int j = 0; j < 4; ++j, ++i) {
      a[i] = AddCarry(a[i], b[i], &carry);
    }
  }

  // Finish remainder.
  for (; i < b.size(); ++i) {
    a[i] = AddCarry(a[i], b[i], &carry);
  }

  // Propagate carry through the rest of a.
  for (; carry && i < a.size(); ++i) {
    a[i] = AddCarry(a[i], 0, &carry);
  }

  return carry;
}

// Computes dst = a + b out of place. Returns the number of bigits actually
// written into dst.
//
// NOTE: dst must be sized to be larger than max(a.size(), b.size()) + 1 (i.e.
// it must be able to hold the carry bigit, if any.
//
// This allows for using a pre-allocated buffer to store the result of an
// addition followed by trimming down to size:
//
//     auto out = arena.Alloc(std::max(a.size(), b.size()) + 1);
//     out = out.first(AddInto(out, a, b));
//
// Which is used in the Karatsuba multiplication, where we don't have the option
// to expand the allocate space on demand.
inline size_t AddOutOfPlace(absl::Span<Bigit> dst, absl::Span<const Bigit> a,
                            absl::Span<const Bigit> b) {
  const size_t max_size = std::max(a.size(), b.size());
  const size_t min_size = std::min(a.size(), b.size());
  ABSL_DCHECK_GE(dst.size(), max_size + 1);

  // Add common parts.
  Bigit carry = 0;

  // Dispatch four at a time to help loop unrolling.
  size_t i = 0;
  while (i + 4 < min_size) {
    for (int j = 0; j < 4; ++j, ++i) {
      dst[i] = AddCarry(a[i], b[i], &carry);
    }
  }

  // Finish remainder of the parts common to A and B.
  for (; i < min_size; ++i) {
    dst[i] = AddCarry(a[i], b[i], &carry);
  }

  // Copy remaining digits from the longer operand and propagate carry.
  auto longer = (a.size() > b.size()) ? a : b;

  // Dispatch four at a time for the remaining part.
  const size_t size = longer.size();
  while (i + 4 < size) {
    for (int j = 0; j < 4; ++j, ++i) {
      dst[i] = AddCarry(longer[i], 0, &carry);
    }
  }

  // Propagate carry through the longer operand.
  for (; i < size; ++i) {
    dst[i] = AddCarry(longer[i], 0, &carry);
  }

  if (carry) {
    dst[i++] = carry;
    return max_size + 1;
  }

  return max_size;
}

// Computes a -= b.
//
// REQUIRES: |a| >= |b|.
inline void SubInPlace(absl::Span<Bigit> a, absl::Span<const Bigit> b) {
  ABSL_DCHECK_GE(a.size(), b.size());
  ABSL_DCHECK_GE(CmpAbs(a, b), 0);

  Bigit borrow = 0;

  // Dispatch four at a time to help loop unrolling.
  size_t size = b.size();
  size_t i = 0;
  while (i + 4 <= size) {
    for (int j = 0; j < 4; ++j, ++i) {
      a[i] = SubBorrow(a[i], b[i], &borrow);
    }
  }

  // Finish remainder of subtraction.
  for (; i < size; ++i) {
    a[i] = SubBorrow(a[i], b[i], &borrow);
  }

  // Propagate the borrow through a.
  for (; borrow && i < a.size(); ++i) {
    borrow = (a[i] == 0);
    a[i]--;
  }
}

// Computes dst = a - b.
//
// Requires |a| >= |b| and dst is thus the same size as a.
// A must be expanded to match the size of B and the total number of digits
// actually set in A must be passed in via a_digits.
inline Bigit SubOutOfPlace(absl::Span<Bigit> dst, absl::Span<const Bigit> a,
                           absl::Span<const Bigit> b, size_t digits) {
  ABSL_DCHECK_EQ(dst.size(), a.size());
  ABSL_DCHECK_GE(CmpAbs(a, b), 1);

  Bigit borrow = 0;

  // Dispatch four at a time to help loop unrolling.
  size_t size = digits;
  size_t i = 0;
  while (i + 4 < size) {
    for (int j = 0; j < 4; ++j, ++i) {
      dst[i] = SubBorrow(a[i], b[i], &borrow);
    }
  }

  // Finish remainder.
  for (; i < digits; ++i) {
    dst[i] = SubBorrow(a[i], b[i], &borrow);
  }

  // Propagate borrow through the rest of a.
  for (; borrow && i < a.size(); ++i) {
    dst[i] = SubBorrow(a[i], 0, &borrow);
  }

  return borrow;
}

Bigit MulAdd(absl::Span<Bigit> out, absl::Span<const Bigit> a, Bigit b,
             Bigit carry) {
  ABSL_DCHECK_GE(out.size(), a.size());

  int left = a.size();

  // Dispatch four at a time to help loop unrolling.
  size_t i = 0;
  while (i + 4 <= a.size()) {
    for (int j = 0; j < 4; ++j, ++i) {
      out[i] = MulCarry(a[i], b, &carry);
      --left;
    }
  }

  for (; i < a.size(); ++i) {
    out[i] = MulCarry(a[i], b, &carry);
  }

  return carry;
}

// Computes out[i] += a[i]*b in place.
//
// Returns the final carry, if any.
inline Bigit MulAddInPlace(absl::Span<Bigit> out, absl::Span<const Bigit> a,
                           Bigit b) {
  // Dispatch four at a time to help loop unrolling.
  Bigit carry = 0;
  size_t i = 0;
  while (i + 4 <= a.size()) {
    for (int j = 0; j < 4; ++j, ++i) {
      MulAddCarry(out[i], a[i], b, &carry);
    }
  }

  // Finish remainder.
  for (; i < a.size(); ++i) {
    MulAddCarry(out[i], a[i], b, &carry);
  }

  return carry;
}

// Implements the standard grade school long multiplication algorithm. The
// output is computed by multiplying A by each digit of B and summing the
// results as we go. This is a quadratic algorithm and only serves as the base
// case for the recursive Karatsuba algorithm below.
//
// NOTE: out must be at least as large as the sums of the sizes of A and B.
inline void MulQuadratic(absl::Span<Bigit> out, absl::Span<const Bigit> a,
                         absl::Span<const Bigit> b) {
  ABSL_DCHECK_GE(out.size(), a.size() + b.size());

  // Make sure A is the longer of the two arguments.
  if (a.size() < b.size()) {
    using std::swap;
    swap(a, b);
  }

  if (b.empty()) {
    absl::c_fill(out, 0);
    return;
  }

  // Each call to MulAdd and MulAddInPlace only updates a.size() elements of out
  // so we manually set the carries as we go. We grab a span to the upper half
  // of out starting at a.size() to facilitate this.
  auto upper = out.subspan(a.size());
  upper[0] = MulAdd(out, a, b[0], 0);

  const size_t size = b.size();
  size_t i = 1;
  for (; i < size; ++i) {
    upper[i] = MulAddInPlace(out.subspan(i), a, b[i]);
  }

  // Finish zeroing out the upper half.
  for (; i < upper.size(); ++i) {
    upper[i] = 0;
  }
}

// Split a span into at most two contiguous spans of length a and b.
//
// If a + b < span.size() then the two spans only cover part of the input.
// If span.size() <= a, then the second span is empty.
template <typename T>
inline std::pair<absl::Span<T>, absl::Span<T>> Split(absl::Span<T> span,
                                                     size_t a, size_t b) {
  if (a < span.size()) {
    return {span.subspan(0, a), span.subspan(a, b)};
  }
  return {span.subspan(0, a), {}};
};

// A simple bump allocator to allow us to very efficiently allocate temporary
// space when recursing in the Karatsuba multiply. The arena is pre-sized and
// returns spans of memory via Alloc() which are then returned to the arena via
// Release.
class Arena {
 public:
  explicit Arena(size_t size) { data_.reserve(size); }

  // Allocates a span of length n from the arena.
  absl::Span<Bigit> Alloc(size_t n) {
    ABSL_DCHECK_LE(used_ + n, data_.capacity());
    size_t start = used_;
    used_ += n;
    return absl::Span<Bigit>(data_.data() + start, n);
  }

  size_t Used() const { return used_; }

  // Resets the arena to the given position which must be < Used().
  void Reset(size_t to) {
    ABSL_DCHECK_LE(to, used_);
    used_ = to;
  }

 private:
  size_t used_ = 0;
  std::vector<Bigit> data_;
};

inline void KaratsubaMulRecursive(absl::Span<Bigit> dst,
                                  absl::Span<const Bigit> a,
                                  absl::Span<const Bigit> b,
                                  Arena* absl_nonnull arena) {
  ABSL_DCHECK_GE(dst.size(), a.size() + b.size());
  if (a.empty() || b.empty()) {
    absl::c_fill(dst, 0);
    return;
  }

  int arena_start = arena->Used();

  // Karatsuba lets us represent two numbers of M bigits each, A and B, as:
  //
  //   A = a1*10^(M/2) + a0
  //   B = b1*10^(M/2) + b0
  //
  // Which we can multiply out:
  //   AB = (a1*10^(M/2) + a0)*(b1*10^(M/2) + b0);
  //      = a1*b1*10^M + (a1*b0 + a0*b1)*10^(M/2) + a0*b0
  //      = z2 * 10^M  + z1*10^(M/2)            + z0
  //
  // Where:
  //   z0 = a0*b0
  //   z1 = a1*b0 + a0*b1
  //   z2 = a1*b1
  //
  // We can replace the multiplications in z1 by computing:
  //
  //   z3 = (a0 + a1)*(b0 + b1)
  //
  // And noting z1 = z3 - z2 - z0
  //
  // This lets us compute a 2M digit multiply with three M digit multiplies,
  // with those individual multiplies able to be recursively divided.

  // Fall back to long multiplication when we're small enough.
  if (std::min(a.size(), b.size()) <= kSimpleMulThreshold) {
    MulQuadratic(dst, a, b);
    return;
  }

  const int half = (std::max(a.size(), b.size()) + 1) / 2;

  // Split the inputs into contiguous subspans.
  auto [a0, a1] = Split(a, half, half);
  auto [b0, b1] = Split(b, half, half);

  // Make space to hold results in the output and multiply sub-terms.
  //   z0 = a0 * b0
  //   z2 = a1 * b1
  auto [z0, z2] = Split(dst, a0.size() + b0.size(), a1.size() + b1.size());
  KaratsubaMulRecursive(z0, a0, b0, arena);
  KaratsubaMulRecursive(z2, a1, b1, arena);

  // Compute (a0 + a1) and (b0 + b1)
  //
  // If the upper terms are zero we can just re-use the terms we have, otherwise
  // we compute the sum and pop off the MSB bigit if no carry occurred.
  absl::Span<const Bigit> asum = a0;
  if (!a1.empty()) {
    absl::Span<Bigit> tmp = arena->Alloc(half + 1);
    asum = tmp.first(AddOutOfPlace(tmp, a0, a1));
  }

  absl::Span<const Bigit> bsum = b0;
  if (!b1.empty()) {
    absl::Span<Bigit> tmp = arena->Alloc(half + 1);
    bsum = tmp.first(AddOutOfPlace(tmp, b0, b1));
  }

  // Compute z1 = asum*bsum - z0 - z2 = (a0 + a1)*(b0 + b1) - z0 - z2
  auto z1 = arena->Alloc(asum.size() + bsum.size());

  // Compute asum * bsum into the beginning of z1
  KaratsubaMulRecursive(z1, asum, bsum, arena);

  // NOTE: (a0 + a1) * (b0 + b1) >= a0*b0 + a1*b1 so this never underflows.
  SubInPlace(z1, z0);
  if (!a1.empty() && !b1.empty()) {
    SubInPlace(z1, z2);
  }

  // Z1 may overflow because of a carry in (a0 + b0) or (a1 + b1) but
  // subtracting z0 and z2 will always bring it back in range, trim any leading
  // zeros to shorten the value if needed.
  while (z1.back() == 0) {
    z1 = z1.first(z1.size() - 1);
  }

  // We need to add z1*10^half which we can do by adding it at an offset.
  AddInPlace(dst.subspan(half), z1);

  // Release temporary memory we used.
  arena->Reset(arena_start);
}

// Multiplies two unsigned bigit vectors together using Karatsuba's algorithm.
//
// This algorithm recursively subdivides the inputs until one or both is below
// some threshold, and then falls back to standard long multiplication.
void KaratsubaMul(absl::Span<Bigit> out, absl::Span<const Bigit> a,
                  absl::Span<const Bigit> b) {
  ABSL_DCHECK_GE(out.size(), a.size() + b.size());
  if (a.empty() || b.empty()) {
    absl::c_fill(out, 0);
    return;
  }

  // Each step of Karatsuba splits at:
  //   N = (std::max(a.size() + b.size() + 1) / 2
  //
  // We have to hold a total of 4*(N + 1) bigits as temporaries at each step.
  //
  // Simulate the recursion (log(n) steps) and compute the arena size.
  int a_size = a.size();
  int b_size = b.size();
  int peak = 0;
  while (std::min(a_size, b_size) > kSimpleMulThreshold) {
    int half = (std::max(a_size, b_size) + 1) / 2;
    int next = half + 1;
    peak += 4 * next;
    a_size = next;
    b_size = next;
  };

  Arena arena(peak);
  KaratsubaMulRecursive(out, a, b, &arena);
}

Bignum& Bignum::operator+=(const Bignum& b) {
  if (b.is_zero()) {
    return *this;
  }

  if (is_zero()) {
    *this = b;
    return *this;
  }

  if (is_negative() == b.is_negative()) {
    // Same sign:
    //   +|a| + +|b| == +(|a| + |b|)
    //   -|a| + -|b| == -(|a| + |b|)
    //
    // So we can just sum magnitudes, final sign is the same as A.
    bigits_.resize(std::max(bigits_.size(), b.bigits_.size()), 0);
    Bigit carry = AddInPlace(absl::MakeSpan(bigits_), b.bigits_);
    if (carry) {
      bigits_.emplace_back(carry);
    }
  } else {
    // We know the signs are different, so there's two options:
    //   -|a| + +|b| = ?(|b| - |a|)
    //   +|a| + -|b| = ?(|a| - |b|)
    //
    // With the final sign being dependent on how |a| and |b| relate.
    if (CmpAbs(bigits_, b.bigits_) >= 0) {
      // |a| >= |b|
      //   -|a| + +|b| --> -(|a| - |b|)
      //   +|a| + -|b| --> +(|a| - |b|)
      //
      // So we can subtract magnitudes, final sign is the same as A.
      SubInPlace(absl::MakeSpan(bigits_), b.bigits_);
    } else {
      // |a| < |b|
      //   -|a| + +|b| --> +(|b| - |a|)
      //   +|a| + -|b| --> -(|b| - |a|)
      //
      // So we can compute |b| - |a| and the final sign is the same as B.
      size_t prev_size = bigits_.size();
      bigits_.resize(b.bigits_.size());
      SubOutOfPlace(absl::MakeSpan(bigits_), b.bigits_, bigits_, prev_size);
      negative_ = b.is_negative();
    }
  }

  Normalize();
  return *this;
}

Bignum& Bignum::operator-=(const Bignum& b) {
  if (this == &b) {
    set_zero();
    return *this;
  }

  // Compute -(-a + b) == a - b
  negate();
  *this += b;
  negate();
  return *this;
}

Bignum& Bignum::operator*=(const Bignum& b) {
  if (is_zero() || b.is_zero()) {
    return set_zero();
  }

  // Result is only negative if signs are different.
  const bool negative = (is_negative() != b.is_negative());

  // Fast path for single-bigit multiplication.
  if (bigits_.size() == 1 && b.bigits_.size() == 1) {
    absl::uint128 prod = absl::uint128(bigits_[0]) * b.bigits_[0];
    const uint64_t lo = absl::Uint128Low64(prod);
    const uint64_t hi = absl::Uint128High64(prod);
    if (hi == 0) {
      bigits_ = {lo};
    } else {
      bigits_ = {lo, hi};
    }
    set_negative(negative);
    return *this;
  }

  // Use Karatsuba multiplication.
  // If the inputs are small enough this will just do long multiplication.
  BigitVector result;
  result.resize(bigits_.size() + b.bigits_.size());
  KaratsubaMul(absl::MakeSpan(result), bigits_, b.bigits_);
  bigits_ = std::move(result);

  negative_ = negative;
  Normalize();
  return *this;
}

}  // namespace exactfloat_internal
