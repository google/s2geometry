#include "s2/util/math/exactfloat/exactfloat_internal.h"

// Threshold for fallback to simple multiplication, determined empirically.
static constexpr int kKaratsubaThreshold = 64;

// Avoid the dependent name clutter.
using Bigit = typename Bignum::Bigit;

static Bigit MulAdd(  //
    absl::Span<Bigit> out, absl::Span<const Bigit> a, Bigit b, Bigit c);

std::optional<Bignum> Bignum::FromString(absl::string_view s) {
  // A chunk is up to 19 decimal digits, which can always fit into a Bigit.
  constexpr int kMaxChunkDigits = std::numeric_limits<uint64_t>::digits10;

  // NOTE: We use a simple multiply-and-add (aka Horner's) method here for the
  // sake of simplicity. This isn't the fastest algorithm, being quadratic in
  // the number of chunks the input has. If we use divide and conquer approach
  // or an FFT based multiply we could probably make this ~O(n^1.5) or
  // semi-linear.

  // Precomputed powers of 10.
  static const auto kPow10 = []() {
    std::array<Bigit, 20> out;

    Bigit value = 1;
    for (int i = 0; i < out.size(); ++i) {
      out[i] = value;
      value = value * 10;
    }
    return out;
  }();

  Bignum out;
  if (s.empty()) {
    return out;
  }

  // Reserve space for bigits.
  out.bigits_.reserve((s.size() + kMaxChunkDigits - 1) / kMaxChunkDigits);

  int sign = +1;
  Bigit chunk = 0;
  int clen = 0;

  // Finish processing the current chunk.
  auto FlushChunk = [&]() {
    if (clen) {
      auto outspan = absl::MakeSpan(out.bigits_);
      if (Bigit carry = MulAdd(outspan, outspan, kPow10[clen], chunk)) {
        out.bigits_.emplace_back(carry);
      }
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

    // Accumulate digit into the local 64-bit chunk.  Skip leading
    // zeros.
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

int bit_width(const Bignum& a) {
  ABSL_DCHECK(a.Normalized());
  if (a.empty()) {
    return 0;
  }

  // Bit width is the bits in the least significant bigits + bit width of
  // the most significant word.
  const int msw_width =
      (Bigit::kBits - absl::countl_zero(a.bigits_.back().value_));
  const int lsw_width = (a.bigits_.size() - 1) * Bigit::kBits;
  return msw_width + lsw_width;
}

int countr_zero(const Bignum& a) {
  if (a.is_zero()) {
    return 0;
  }

  int nzero = 0;
  for (Bigit bigit : a.bigits_) {
    if (bigit == 0) {
      nzero += Bigit::kBits;
    } else {
      nzero += absl::countr_zero(static_cast<uint64_t>(bigit));
      break;
    }
  }
  return nzero;
}

bool Bignum::Bit(int nbit) const {
  ABSL_DCHECK_GE(nbit, 0);
  if (is_zero()) {
    return false;
  }

  const int digit = nbit / Bigit::kBits;
  const int shift = nbit % Bigit::kBits;

  if (digit >= size()) {
    return false;
  }

  return ((bigits_[digit] >> shift) & 0x1) != 0;
}

Bignum Bignum::operator-() const {
  Bignum result = *this;
  result.sign_ = -result.sign_;
  return result;
}

Bignum& Bignum::operator<<=(int nbit) {
  ABSL_DCHECK_GE(nbit, 0);
  if (is_zero() || nbit == 0) {
    return *this;
  }

  const int nbigit = nbit / Bigit::kBits;
  const int nrem = nbit % Bigit::kBits;

  // First, handle the whole-bigit shift by inserting zeros.
  bigits_.insert(bigits_.begin(), nbigit, 0);

  // Then, handle the within-bigit shift, if any.
  if (nrem != 0) {
    Bigit carry = 0;
    for (size_t i = 0; i < bigits_.size(); ++i) {
      const Bigit old_val = bigits_[i];
      bigits_[i] = (old_val << nrem) | carry;
      carry = old_val >> (Bigit::kBits - nrem);
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
    return SetZero();
  }

  const int nbigit = nbit / Bigit::kBits;
  const int nrem = nbit % Bigit::kBits;

  // First, handle the whole-bigit shift by removing bigits.
  bigits_.erase(bigits_.begin(), bigits_.begin() + nbigit);

  // Then, handle the within-bigit shift, if any.
  if (nrem != 0) {
    Bigit carry = 0;
    for (int i = static_cast<int>(bigits_.size()) - 1; i >= 0; --i) {
      const Bigit old_val = bigits_[i];
      bigits_[i] = (old_val >> nrem) | carry;
      carry = old_val << (Bigit::kBits - nrem);
    }
  }

  // Result might be smaller or zero, so normalize.
  NormalizeSign(sign_);
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

// Computes a + b + c and updates the carry.
static Bigit AddCarry(Bigit a, Bigit b, Bigit& c) {
  auto sum = absl::uint128(a) + b + c;
  c = absl::Uint128High64(sum);
  return static_cast<Bigit>(sum);
}

// Computes a - b - c and updates the borrow.
static Bigit SubBorrow(Bigit a, Bigit b, Bigit& borrow) {
  Bigit diff = a - b - borrow;
  borrow = (a < b) || (borrow && (a == b));
  return diff;
}

// Computes a * b + c and updates the carry.
static Bigit MulCarry(Bigit a, Bigit b, Bigit& c) {
  auto sum = absl::uint128(a) * b + c;
  c = absl::Uint128High64(sum);
  return static_cast<Bigit>(sum);
}

// Computes out += a * b + c and updates the carry.
static void MulAddCarry(Bigit& out, Bigit a, Bigit b, Bigit& c) {
  auto sum = absl::uint128(a) * b + c + out;
  c = absl::Uint128High64(sum);
  out = static_cast<Bigit>(sum);
}

// Computes a += b in place. Returns the final carry (if any).
// NOTE: the a operand must be pre-expanded to fit b.
static Bigit AddInPlace(absl::Span<Bigit> a, absl::Span<const Bigit> b) {
  ABSL_DCHECK_GE(a.size(), b.size());

  Bigit* pa = a.data();
  const Bigit* pb = b.data();

  int left = b.size();
  Bigit carry = 0;

  // Dispatch four at a time to help loop unrolling.
  while (left >= 4) {
    for (int i = 0; i < 4; ++i) {
      *pa = AddCarry(*pa, *pb++, carry);
      ++pa;
      --left;
    }
  }

  // Finish remainder.
  while (left--) {
    *pa = AddCarry(*pa, *pb++, carry);
    ++pa;
  }

  // Propagate carry through the rest of a.
  int remaining = a.size() - b.size();
  while (carry && remaining--) {
    *pa = AddCarry(*pa, 0, carry);
    ++pa;
  }

  return carry;
}

static ssize_t AddInto(  //
    absl::Span<Bigit> dst, absl::Span<const Bigit> a,
    absl::Span<const Bigit> b) {
  const size_t max_size = std::max(a.size(), b.size());
  const size_t min_size = std::min(a.size(), b.size());
  ABSL_DCHECK_GE(dst.size(), max_size + 1);

  Bigit* pdst = dst.data();
  const Bigit* pa = a.data();
  const Bigit* pb = b.data();

  // Add common parts.
  Bigit carry = 0;

  // Dispatch four at a time to help loop unrolling.
  int size = min_size;
  int i = 0;
  while (size >= i + 4) {
    for (int j = 0; j < 4; ++j) {
      pdst[i] = AddCarry(pa[i], pb[i], carry);
      ++i;
    }
  }

  // Finish remainder of common parts.
  for (; i < size; ++i) {
    pdst[i] = AddCarry(pa[i], pb[i], carry);
  }

  // Copy remaining digits from the longer operand and propagate carry.
  auto longer = (a.size() > b.size()) ? a : b;
  const Bigit* plonger = (a.size() > b.size()) ? pa : pb;

  // Dispatch four at a time for the remaining part.
  size = longer.size();
  while (size >= i + 4) {
    for (int j = 0; j < 4; ++j) {
      pdst[i] = AddCarry(plonger[i], 0, carry);
      ++i;
    }
  }

  // Finish remainder.
  for (; i < size; ++i) {
    pdst[i] = AddCarry(plonger[i], 0, carry);
  }

  if (carry) {
    pdst[i++] = carry;
    return max_size + 1;
  }
  return max_size;
}

// Computes a -= b. Returns the final borrow (if any).
//
// REQUIRES: |a| < |b|.
// NOTE: A must be pre-expanded to match the size of b.
static Bigit SubLtIp(  //
    absl::Span<Bigit> a, absl::Span<const Bigit> b, ssize_t na) {
  ABSL_DCHECK_EQ(a.size(), b.size());

  Bigit* pa = a.data();
  const Bigit* pb = b.data();
  Bigit borrow = 0;

  // Dispatch four at a time to help loop unrolling.
  int size = na;
  int i = 0;
  while (size >= i + 4) {
    for (int j = 0; j < 4; ++j) {
      pa[i] = SubBorrow(pb[i], pa[i], borrow);
      ++i;
    }
  }

  // Finish remainder.
  for (; i < na; ++i) {
    pa[i] = SubBorrow(pb[i], pa[i], borrow);
  }

  // Propagate borrow through the rest of b.
  for (; borrow && i < b.size(); ++i) {
    pa[i] = SubBorrow(pb[i], 0, borrow);
  }
  return borrow;
}

// Computes a -= b. Returns the final borrow (if any).
//
// REQUIRES: |a| >= |b|.
static Bigit SubGeIp(absl::Span<Bigit> a, absl::Span<const Bigit> b) {
  ABSL_DCHECK_GE(a.size(), b.size());

  Bigit borrow = 0;

  Bigit* pa = a.data();
  const Bigit* pb = b.data();

  // Dispatch four at a time to help loop unrolling.
  int size = b.size();
  int done = 0;
  while (size >= done + 4) {
    for (int i = 0; i < 4; ++i) {
      pa[done] = SubBorrow(pa[done], pb[done], borrow);
      ++done;
    }
  }

  // Finish remainder of subtraction.
  for (; done < size; ++done) {
    pa[done] = SubBorrow(pa[done], pb[done], borrow);
  }

  // Propagate the borrow through a.
  for (; borrow && done < a.size(); ++done) {
    borrow = (a[done] == 0);
    a[done]--;
  }
  return borrow;
}

// Computes out[i] = a[i]*b + c
//
// Returns the final carry, if any.
Bigit MulAdd(  //
    absl::Span<Bigit> out, absl::Span<const Bigit> a, Bigit b, Bigit c = 0) {
  ABSL_DCHECK_GE(out.size(), a.size());

  Bigit* pout = out.data();
  const Bigit* pa = a.data();

  int left = a.size();

  // Dispatch four at a time to help loop unrolling.
  while (left >= 4) {
    for (int i = 0; i < 4; ++i) {
      *pout++ = MulCarry(*pa++, b, c);
      --left;
    }
  }

  while (left--) {
    *pout++ = MulCarry(*pa++, b, c);
  }
  return c;
}

// Computes out[i] += a[i]*b in place.
//
// Returns the final carry, if any.
static Bigit MulAddIp(  //
    absl::Span<Bigit> out, absl::Span<const Bigit> a, Bigit b) {
  Bigit* pout = out.data();
  const Bigit* pa = a.data();

  int left = a.size();

  // Dispatch four at a time to help loop unrolling.
  Bigit carry = 0;
  while (left >= 4) {
    for (int i = 0; i < 4; ++i) {
      MulAddCarry(*pout++, *pa++, b, carry);
      --left;
    }
  }

  // Finish remainder.
  while (left--) {
    MulAddCarry(*pout++, *pa++, b, carry);
  }

  return carry;
}

static void MulQuadratic(   //
    absl::Span<Bigit> out,  //
    absl::Span<const Bigit> a, absl::Span<const Bigit> b) {
  ABSL_DCHECK_EQ(out.size(), a.size() + b.size());

  // Make sure A is the longer of the two arguments.
  if (a.size() < b.size()) {
    using std::swap;
    swap(a, b);
  }

  if (b.empty()) {
    absl::c_fill(out, 0);
    return;
  }

  auto upper = out.subspan(a.size());
  upper[0] = MulAdd(out, a, b[0]);

  const int size = b.size();
  int i = 1;
  while (size >= i + 4) {
    for (int j = 0; j < 4; ++j) {
      upper[i] = MulAddIp(out.subspan(i), a, b[i]);
      ++i;
    }
  }

  // Finish remainder (if any).
  for (; i < size; ++i) {
    upper[i] = MulAddIp(out.subspan(i), a, b[i]);
  }

  // Finish zeroing out upper half.
  for (; i < upper.size(); ++i) {
    upper[i] = 0;
  }
}

// Split a span into two contiguous pieces of length a and b, respectively.
template <typename T>
static std::pair<absl::Span<T>, absl::Span<T>> Split(  //
    absl::Span<T> span, int a, int b) {
  return {span.subspan(0, a), span.subspan(a, b)};
};

// A simple bump allocator to avoid allocating memory during recursion.
class Arena {
 public:
  explicit Arena(ssize_t size) { data_.reserve(size); }

  // Allocates a span of length n from the arena.
  absl::Span<Bigit> Alloc(ssize_t n) {
    ABSL_DCHECK_LE(used_ + n, data_.capacity());
    size_t start = used_;
    used_ += n;
    return absl::Span<Bigit>(data_.data() + start, n);
  }

  void Release(ssize_t n) {
    ABSL_DCHECK_LE(n, used_);
    used_ -= n;
  }

 private:
  ssize_t used_ = 0;
  absl::InlinedVector<Bigit, 1024> data_;
};

static void KaratsubaMulRec(  //
    absl::Span<Bigit> dst,    //
    absl::Span<const Bigit> a, absl::Span<const Bigit> b, Arena& arena) {
  ABSL_DCHECK_EQ(dst.size(), a.size() + b.size());
  if (a.empty() || b.empty()) {
    return;
  }

  // Karatsuba lets us represent two numbers, A and B thusly:
  //   A = a1*10^M + a0
  //   B = b1*10^M + b0
  //
  // Which we can multiply out:
  //   AB = (a1*10^M + a0)*(b1*10^M + b0);
  //      = a1*b1*10^(2M) + (a1*b0 + a0*b1)*10^M + a0*b0
  //      = z2 * 10^2M    + z1*10^M              + z0
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
  if (dst.size() < kKaratsubaThreshold) {
    MulQuadratic(dst, a, b);
    return;
  }

  const int half = (std::min(a.size(), b.size()) + 1) / 2;

  // Split the inputs into contiguous subspans.
  auto [a0, a1] = Split(a, half, half);
  auto [b0, b1] = Split(b, half, half);

  // Make space to hold results in the output and multiply sub-terms.
  //   z0 = a0 * b0
  //   z2 = a1 * b1
  auto [z0, z2] = Split(dst, 2 * half, 2 * half);

  KaratsubaMulRec(z0, a0, b0, arena);
  KaratsubaMulRec(z2, a1, b1, arena);

  // Compute (a0 + a1) and (b0 + b1) using space from the arena.
  //
  // The sums may or may not carry. We pop the extra bigit off if they
  // don't.
  auto sa = arena.Alloc(half + 1);
  auto sb = arena.Alloc(half + 1);
  sa = sa.first(AddInto(sa, a0, a1));
  sb = sb.first(AddInto(sb, b0, b1));

  // Compute z1 = sa*sb - z0 - z2 = (a0 + a1)*(b0 + b1) - z0 - z2
  auto z1 = arena.Alloc(sa.size() + sb.size());

  // Compute sa * sb into the beginning of z1
  KaratsubaMulRec(z1, sa, sb, arena);

  // NOTE: (a0 + a1) * (b0 + b1) >= a0*b0 + a1*b1 so this never underflows.
  SubGeIp(z1, z0);
  SubGeIp(z1, z2);

  // We need to add z1*10^half which we can do by adding it offset.
  AddInPlace(dst.subspan(half), z1);

  // Release temporary memory we used.
  arena.Release(z1.size() + sb.size() + sa.size());
}

Bignum::BigitVector Bignum::KaratsubaMul(  //
    absl::Span<const Bigit> a, absl::Span<const Bigit> b) {
  if (a.empty() || b.empty()) {
    return {};
  }

  // Each step of Karatsuba splits at:
  //   N = std::ceil(std::min(a.size(), b.size())/2)
  //
  // We have to hold a total of 4*(N + 1) bigits as temporaries at each step.
  //
  // Simulate the recursion (log(n) steps) and compute the arena size.
  int size = a.size() + b.size();
  int peak = 0;
  do {
    int half = (size + 1) / 2;
    int next = half + 1;
    peak += 4 * next;
    size = next;
  } while (size > kKaratsubaThreshold);

  Arena arena(peak);
  BigitVector out(a.size() + b.size(), 0);
  KaratsubaMulRec(absl::MakeSpan(out), a, b, arena);
  return out;
}

Bignum& Bignum::operator+=(const Bignum& b) {
  if (b.is_zero()) {
    return *this;
  }

  if (is_zero()) {
    *this = b;
    return *this;
  }

  if (sign_ == b.sign_) {
    // Same sign:
    //   (+a) + (+b) == +(a + b)
    //   (-a) + (-b) == -(a + b)
    bigits_.resize(std::max(size(), b.size()), 0);
    Bigit carry = AddInPlace(absl::MakeSpan(bigits_), b.bigits_);
    if (carry) {
      bigits_.emplace_back(carry);
    }
    Normalize();
  } else {
    if (CmpAbs(b) >= 0) {
      // |a| >= |b|, so a - b is same sign as a.
      SubGeIp(absl::MakeSpan(bigits_), b.bigits_);
      NormalizeSign(sign_);
    } else {
      // |a| < |b|, so a - b is same sign as b.
      const int prev_size = size();
      bigits_.resize(b.size());
      SubLtIp(absl::MakeSpan(bigits_), b.bigits_, prev_size);
      NormalizeSign(b.sign_);
    }
  }

  return *this;
}

Bignum& Bignum::operator-=(const Bignum& b) {
  if (this == &b) {
    return SetZero();
  }

  if (b.is_zero()) {
    return *this;
  }

  if (is_zero()) {
    return *this = -b;
  }

  if (sign_ != b.sign_) {
    bigits_.resize(std::max(size(), b.size()), 0);
    uint64_t carry = AddInPlace(absl::MakeSpan(bigits_), b.bigits_);
    if (carry) {
      bigits_.emplace_back(carry);
    }
    Normalize();
  } else {
    if (CmpAbs(b) >= 0) {
      SubGeIp(absl::MakeSpan(bigits_), b.bigits_);
      NormalizeSign(sign_);
    } else {
      const int prev_size = size();
      bigits_.resize(b.size());
      SubLtIp(absl::MakeSpan(bigits_), b.bigits_, prev_size);
      NormalizeSign(-sign_);
    }
  }

  return *this;
}

Bignum& Bignum::operator*=(const Bignum& b) {
  if (is_zero() || b.is_zero()) {
    return SetZero();
  }

  const int new_sign = sign_ * b.sign_;

  // Fast path for single-bigit multiplication.
  if (size() == 1 && b.size() == 1) {
    absl::uint128 prod = absl::uint128(bigits_[0]) * b.bigits_[0];
    const uint64_t lo = absl::Uint128Low64(prod);
    const uint64_t hi = absl::Uint128High64(prod);
    if (hi == 0) {
      bigits_ = {lo};
    } else {
      bigits_ = {lo, hi};
    }
    sign_ = new_sign;
    return *this;
  }

  // Use Karatsuba multiplication.
  // If the inputs are small enough this will just do long multiplication.
  bigits_ = KaratsubaMul(bigits_, b.bigits_);
  NormalizeSign(new_sign);
  return *this;
}
