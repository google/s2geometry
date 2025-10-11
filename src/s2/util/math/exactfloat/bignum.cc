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

#ifdef __x86_64__
#include <immintrin.h>
#endif

#include <algorithm>
#include <array>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <utility>

#include "absl/algorithm/container.h"
#include "absl/base/nullability.h"
#include "absl/log/absl_check.h"
#include "absl/numeric/bits.h"
#include "absl/numeric/int128.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"

namespace exactfloat_internal {

// Number of bigits in smaller of the two operands before we fall back to simple
// multiplication in the Karatsuba recursion. Determined empirically.
static constexpr int kSimpleMulThreshold = 24;

// Computes dst[i] = a[i]*b + c
//
// Returns the final carry, if any.
inline Bigit MulWithCarry(absl::Span<Bigit> dst, absl::Span<const Bigit> a,
                          Bigit b, Bigit carry);

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
  constexpr size_t kMaxChunkDigits = std::numeric_limits<Bigit>::digits10;

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
    int64_t chunk_len = std::min(static_cast<size_t>(std::distance(begin, end)),
                                 kMaxChunkDigits);

    Bigit chunk = 0;
    auto result = std::from_chars(begin, begin + chunk_len, chunk);
    if (result.ec != std::errc() || (result.ptr - begin) != chunk_len) {
      return std::nullopt;
    }
    begin += chunk_len;

    // Shift left by chunk_len digits and add the chunk to it.
    auto outspan = absl::MakeSpan(out.bigits_);
    Bigit carry = MulWithCarry(outspan, outspan, kPow10[chunk_len], chunk);
    if (carry) {
      out.bigits_.emplace_back(carry);
    }
  }

  out.negative_ = negative;
  out.Normalize();
  return out;
}

int bit_width(const Bignum& a) {
  ABSL_DCHECK(a.bigits_.empty() || a.bigits_.back() != 0);
  if (a.is_zero()) {
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
  const size_t digit = nbit / kBigitBits;
  const size_t shift = nbit % kBigitBits;

  if (digit >= bigits_.size()) {
    return false;
  }

  return ((bigits_[digit] >> shift) & 0x1) != 0;
}

Bignum Bignum::operator-() const {
  Bignum result = *this;
  result.negate();
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
inline Bigit AddBigit(Bigit a, Bigit b, Bigit* absl_nonnull carry) {
#ifdef __x86_64__
  Bigit out;
  *carry =
      _addcarry_u64(*carry, a, b, reinterpret_cast<unsigned long long*>(&out));
  return out;
#else
  auto sum = absl::uint128(a) + b + *carry;
  *carry = absl::Uint128High64(sum);
  return static_cast<Bigit>(sum);
#endif
}

// Computes a - b - borrow  and updates the borrow.
//
// NOTE: Borrow must be one or zero.
inline Bigit SubBigit(Bigit a, Bigit b, Bigit* absl_nonnull borrow) {
  ABSL_DCHECK_LE(*borrow, Bigit(1));
#ifdef __x86_64__
  Bigit out;
  *borrow = _subborrow_u64(*borrow, a, b,
                           reinterpret_cast<unsigned long long*>(&out));
  return out;
#else
  Bigit diff = a - b - *borrow;
  *borrow = (a < b) || (*borrow && (a == b));
  return diff;
#endif
}

// Computes a * b + carry and updates the carry.
inline Bigit MulBigit(Bigit a, Bigit b, Bigit* absl_nonnull carry) {
  auto sum = absl::uint128(a) * b + *carry;
  *carry = absl::Uint128High64(sum);
  return static_cast<Bigit>(sum);
}

// Computes sum += a * b + carry and updates the carry.
//
// NOTE: Will not overflow even if a, b, and c are their maximum values.
inline void MulAddBigit(Bigit* absl_nonnull sum, Bigit a, Bigit b,
                        Bigit* absl_nonnull carry) {
  auto term = absl::uint128(a) * b + *carry + *sum;
  *carry = absl::Uint128High64(term);
  *sum = static_cast<Bigit>(term);
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
      a[i] = AddBigit(a[i], b[i], &carry);
    }
  }

  // Finish remainder.
  for (; i < b.size(); ++i) {
    a[i] = AddBigit(a[i], b[i], &carry);
  }

  // Propagate carry through the rest of a.
  for (; carry && i < a.size(); ++i) {
    a[i] = AddBigit(a[i], 0, &carry);
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
inline size_t Add(absl::Span<Bigit> dst, absl::Span<const Bigit> a,
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
      dst[i] = AddBigit(a[i], b[i], &carry);
    }
  }

  // Finish remainder of the parts common to A and B.
  for (; i < min_size; ++i) {
    dst[i] = AddBigit(a[i], b[i], &carry);
  }

  // Copy remaining digits from the longer operand and propagate carry.
  auto longer = (a.size() > b.size()) ? a : b;

  // Dispatch four at a time for the remaining part.
  const size_t size = longer.size();
  while (i + 4 < size) {
    for (int j = 0; j < 4; ++j, ++i) {
      dst[i] = AddBigit(longer[i], 0, &carry);
    }
  }

  // Propagate carry through the longer operand.
  for (; i < size; ++i) {
    dst[i] = AddBigit(longer[i], 0, &carry);
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
      a[i] = SubBigit(a[i], b[i], &borrow);
    }
  }

  // Finish remainder of subtraction.
  for (; i < size; ++i) {
    a[i] = SubBigit(a[i], b[i], &borrow);
  }

  // Propagate the borrow through a.
  for (; borrow && i < a.size(); ++i) {
    borrow = (a[i] == 0);
    a[i]--;
  }
}

// Computes a = b - a.
//
// NOTE: Requires |b| >= |a|.
//
// Since we write the result to a, but b is larger, a must be expanded with
// enough leading zeros to fit the result.
inline void SubReverseInPlace(absl::Span<Bigit> a, absl::Span<const Bigit> b) {
  ABSL_DCHECK_GE(a.size(), b.size());
  ABSL_DCHECK_GE(CmpAbs(b, a), 0);

  Bigit borrow = 0;

  // Dispatch four at a time to help loop unrolling.
  size_t size = a.size();
  size_t i = 0;
  while (i + 4 <= size) {
    for (int j = 0; j < 4; ++j, ++i) {
      a[i] = SubBigit(b[i], a[i], &borrow);
    }
  }

  // Finish remainder.
  for (; i < size; ++i) {
    a[i] = SubBigit(b[i], a[i], &borrow);
  }

  // Propagate borrow through the rest of a.
  for (; borrow && i < a.size(); ++i) {
    a[i] = SubBigit(a[i], 0, &borrow);
  }
}

inline Bigit MulWithCarry(absl::Span<Bigit> dst, absl::Span<const Bigit> a,
                          Bigit b, Bigit carry) {
  ABSL_DCHECK_GE(dst.size(), a.size());

  // Dispatch four at a time to help loop unrolling.
  size_t i = 0;
  while (i + 4 <= a.size()) {
    for (int j = 0; j < 4; ++j, ++i) {
      dst[i] = MulBigit(a[i], b, &carry);
    }
  }

  for (; i < a.size(); ++i) {
    dst[i] = MulBigit(a[i], b, &carry);
  }

  return carry;
}

// Computes sum[i] += a[i]*b in place.
//
// Returns the final carry, if any.
inline Bigit MulAddInPlace(absl::Span<Bigit> sum, absl::Span<const Bigit> a,
                           Bigit b) {
  // Dispatch four at a time to help loop unrolling.
  Bigit carry = 0;
  size_t i = 0;
  while (i + 4 <= a.size()) {
    for (int j = 0; j < 4; ++j, ++i) {
      MulAddBigit(&sum[i], a[i], b, &carry);
    }
  }

  // Finish remainder.
  for (; i < a.size(); ++i) {
    MulAddBigit(&sum[i], a[i], b, &carry);
  }

  return carry;
}

// Implements the standard grade school long multiplication algorithm. The
// output is computed by multiplying A by each digit of B and summing the
// results as we go. This is a quadratic algorithm and only serves as the base
// case for the recursive Karatsuba algorithm below.
//
// NOTE: out must be at least as large as the sums of the sizes of A and B.
inline void MulQuadratic(absl::Span<Bigit> dst, absl::Span<const Bigit> a,
                         absl::Span<const Bigit> b) {
  ABSL_DCHECK_GE(dst.size(), a.size() + b.size());

  // Make sure A is the longer of the two arguments.
  if (a.size() < b.size()) {
    using std::swap;
    swap(a, b);
  }

  if (b.empty()) {
    absl::c_fill(dst, 0);
    return;
  }

  // Each call to MulAdd and MulAddInPlace only updates a.size() elements of out
  // so we manually set the carries as we go. We grab a span to the upper half
  // of out starting at a.size() to facilitate this.
  auto upper = dst.subspan(a.size());
  upper[0] = MulWithCarry(dst, a, b[0], 0);

  const size_t size = b.size();
  size_t i = 1;
  for (; i < size; ++i) {
    upper[i] = MulAddInPlace(dst.subspan(i), a, b[i]);
  }

  // Finish zeroing out the upper half.
  for (; i < upper.size(); ++i) {
    upper[i] = 0;
  }
}

// Split a span into two contiguous spans of length at most a and b.
//
// If span.size() <= a, the second span is empty, otherwise the second span
// has length at most b. If span.size() > a + b, then the two spans only cover
// part of the input span.
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
//
// NOTE: We use std::unique_ptr here instead of std::vector because we don't
// want to initialize the memory unnecessarily and using std::vector without
// resizing (and thus initializing) the container leads to false positives with
// ASAN.
class Arena {
 public:
  // TODO: Use make_unique_for_overwrite when on C++20.
  explicit Arena(size_t size) : size_(size), data_(new Bigit[size]) {}

  // Allocates a span of length n from the arena.
  absl::Span<Bigit> Alloc(size_t n) {
    ABSL_DCHECK_LE(used_ + n, size_);
    size_t start = used_;
    used_ += n;
    return absl::Span<Bigit>(data_.get() + start, n);
  }

  size_t Available() const { return size_ - used_; }

  size_t Used() const { return used_; }

  // Resets the arena to the given position which must be < Used().
  void Reset(size_t to) {
    ABSL_DCHECK_LE(to, used_);
    used_ = to;
  }

 private:
  size_t size_ = 0;
  size_t used_ = 0;
  std::unique_ptr<Bigit[]> data_;
};

// Returns the total arena size needed to multiply two number of a_size and
// b_size bigits using the recursive Karatsuba implementation.
inline size_t ArenaSize(size_t a_size, size_t b_size) {
  // Each step of Karatsuba splits at:
  //   N = (std::max(a.size() + b.size() + 1) / 2
  //
  // We have to hold a total of 4*(N + 1) bigits as temporaries at each step.
  //
  // Simulate the recursion (log(n) steps) and compute the arena size.
  int peak = 0;
  while (std::min(a_size, b_size) > kSimpleMulThreshold) {
    int half = (std::max(a_size, b_size) + 1) / 2;
    int next = half + 1;
    peak += 4 * next;
    a_size = next;
    b_size = next;
  };
  return peak;
}

// Recursive step in the Karatsuba multiplication. dst must be large enough to
// hold the product of a and b (i.e. it must be at least as large as a.size() +
// b.size()). The product is computed and stored in-place in dst.
//
// Additionally an arena must be provided for temporary storage for intermediate
// products. The arena must have at least ArenaSize(a.size(), b.size()) space
// available.
inline void KaratsubaMulRecursive(absl::Span<Bigit> dst,
                                  absl::Span<const Bigit> a,
                                  absl::Span<const Bigit> b,
                                  Arena* absl_nonnull arena) {
  ABSL_DCHECK_GE(dst.size(), a.size() + b.size());
  ABSL_DCHECK_GE(arena->Available(), ArenaSize(a.size(), b.size()));
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

  const size_t half = (std::max(a.size(), b.size()) + 1) / 2;

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
    asum = tmp.first(Add(tmp, a0, a1));
  }

  absl::Span<const Bigit> bsum = b0;
  if (!b1.empty()) {
    absl::Span<Bigit> tmp = arena->Alloc(half + 1);
    bsum = tmp.first(Add(tmp, b0, b1));
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

  // We need to add z1*10^half, which we can do by simply adding z1 at a shifted
  // position in the output.
  auto dst_z1 = dst.subspan(half);

  // Although the value of z1 is guaranteed to fit in the available space of
  // dst, it may have one or more high-order zero bigits because it was sized
  // conservatively to hold the intermediate result (asum * bsum). We trim these
  // leading zeros if necessary to ensure that the Add() operation below does
  // not attempt to write zero bigits past the end of dst.
  AddInPlace(dst_z1, z1.first(std::min(z1.size(), dst_z1.size())));

  // Release temporary memory we used.
  arena->Reset(arena_start);
}

// Multiplies two unsigned bigit vectors together using Karatsuba's algorithm.
//
// This algorithm recursively subdivides the inputs until one or both is below
// some threshold, and then falls back to standard long multiplication.
void KaratsubaMul(absl::Span<Bigit> dst, absl::Span<const Bigit> a,
                  absl::Span<const Bigit> b) {
  ABSL_DCHECK_GE(dst.size(), a.size() + b.size());
  if (a.empty() || b.empty()) {
    absl::c_fill(dst, 0);
    return;
  }

  Arena arena(ArenaSize(a.size(), b.size()));
  KaratsubaMulRecursive(dst, a, b, &arena);
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
      bigits_.resize(b.bigits_.size());
      SubReverseInPlace(absl::MakeSpan(bigits_), b.bigits_);
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
