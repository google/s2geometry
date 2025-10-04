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

#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <random>
#include <string>
#include <vector>

// TODO: remove once benchmarks are available
#if 0
#include "benchmark/benchmark.h"
#endif

#include "absl/random/bit_gen_ref.h"
#include "absl/random/random.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_format.h"
#include "absl/strings/string_view.h"
#include "gtest/gtest.h"
#include "openssl/bn.h"
#include "openssl/crypto.h"

namespace exactfloat_internal {

using ::testing::TestWithParam;

constexpr uint64_t kU8max = std::numeric_limits<uint8_t>::max();
constexpr uint64_t kU16max = std::numeric_limits<uint16_t>::max();
constexpr uint64_t kU32max = std::numeric_limits<uint32_t>::max();
constexpr uint64_t kU64max = std::numeric_limits<uint64_t>::max();

constexpr int64_t kI8max = std::numeric_limits<int8_t>::max();
constexpr int64_t kI16max = std::numeric_limits<int16_t>::max();
constexpr int64_t kI32max = std::numeric_limits<int32_t>::max();
constexpr int64_t kI64max = std::numeric_limits<int64_t>::max();

constexpr int64_t kI8min = std::numeric_limits<int8_t>::min();
constexpr int64_t kI16min = std::numeric_limits<int16_t>::min();
constexpr int64_t kI32min = std::numeric_limits<int32_t>::min();
constexpr int64_t kI64min = std::numeric_limits<int64_t>::min();

// To reduce duplication.
inline auto Bn(absl::string_view str) { return Bignum::FromString(str); };

TEST(BignumTest, ZeroInputsNormalizeToZero) {
  EXPECT_EQ(Bn(""), Bignum(0));
  EXPECT_EQ(Bn("+"), Bignum(0));
  EXPECT_EQ(Bn("-"), Bignum(0));
  EXPECT_EQ(Bn("0"), Bignum(0));
  EXPECT_EQ(Bn("-0"), Bignum(0));
  EXPECT_EQ(Bn("000000"), Bignum(0));
  EXPECT_EQ(Bn("-000000"), Bignum(0));
  EXPECT_EQ(Bn("-00000000000000000000"), Bignum(0));
  EXPECT_EQ(Bn("+00000000000000000000"), Bignum(0));
}

TEST(BignumTest, BasicSmallNumbers) {
  EXPECT_EQ(Bn("42"), Bignum(42));
  EXPECT_EQ(Bn("-17"), Bignum(-17));
  EXPECT_EQ(Bn("+123"), Bignum(123));
}

TEST(BignumTest, LeadingZeros) {
  EXPECT_EQ(Bn("0000042"), Bignum(42));
  EXPECT_EQ(Bn("-00017"), Bignum(-17));
  EXPECT_EQ(Bn("+00007"), Bignum(7));

  // Larger bignums.
  // 10^19 and 10^19 - 1, near chunk boundaries.
  EXPECT_EQ(Bn("100000000000000000"), Bn("000100000000000000000"));
  EXPECT_EQ(Bn("99999999999999999"), Bn("000099999999999999999"));

  // 2^64 - 1 cross-check against integral constructor.
  EXPECT_EQ(Bn("18446744073709551615"), Bignum(18446744073709551615ull));

  // 2^64 cross-check with leading zeros variant via string constructor.
  EXPECT_EQ(Bn("18446744073709551616"), Bn("00018446744073709551616"));
  EXPECT_EQ(Bn("-18446744073709551616"), Bn("-00018446744073709551616"));

  // Multi-digit bignums.
  // 2^80
  EXPECT_EQ(Bn("1208925819614629174706176"),
            Bn("0001208925819614629174706176"));
  EXPECT_EQ(Bn("-1208925819614629174706176"),
            Bn("-0001208925819614629174706176"));

  // 2^128
  EXPECT_EQ(Bn("340282366920938463463374607431768211456"),
            Bn("000340282366920938463463374607431768211456"));
  EXPECT_EQ(Bn("-340282366920938463463374607431768211456"),
            Bn("-000340282366920938463463374607431768211456"));
}

TEST(BignumTest, MultipleSignsBeforeDigitsCausesFailure) {
  EXPECT_EQ(Bn("++123"), std::nullopt);
  EXPECT_EQ(Bn("--42"), std::nullopt);
  EXPECT_EQ(Bn("+-9"), std::nullopt);
  EXPECT_EQ(Bn("-+9"), std::nullopt);
}

TEST(BignumTest, SignInWrongPlaceCausesFailure) {
  EXPECT_EQ(Bn("123-"), std::nullopt);
  EXPECT_EQ(Bn("456+"), std::nullopt);
  EXPECT_EQ(Bn("-789+"), std::nullopt);
  EXPECT_EQ(Bn("+314-"), std::nullopt);
}

TEST(BignumTest, ZeroAlwaysFitsIn) {
  const Bignum zero(0);
  EXPECT_TRUE(zero.FitsIn<int8_t>());
  EXPECT_TRUE(zero.FitsIn<uint8_t>());
  EXPECT_TRUE(zero.FitsIn<int16_t>());
  EXPECT_TRUE(zero.FitsIn<uint16_t>());
  EXPECT_TRUE(zero.FitsIn<int32_t>());
  EXPECT_TRUE(zero.FitsIn<uint32_t>());
  EXPECT_TRUE(zero.FitsIn<int64_t>());
  EXPECT_TRUE(zero.FitsIn<uint64_t>());
}

TEST(BignumTest, ZeroAlwaysCastsToZero) {
  Bignum zero;
  EXPECT_EQ(zero.Cast<int8_t>(), 0);
  EXPECT_EQ(zero.Cast<uint8_t>(), 0);
  EXPECT_EQ(zero.Cast<int16_t>(), 0);
  EXPECT_EQ(zero.Cast<uint16_t>(), 0);
  EXPECT_EQ(zero.Cast<int32_t>(), 0);
  EXPECT_EQ(zero.Cast<uint32_t>(), 0);
  EXPECT_EQ(zero.Cast<int64_t>(), 0);
  EXPECT_EQ(zero.Cast<uint64_t>(), 0);
}

TEST(BignumTest, NegativeOnlyFitsInSigned) {
  const Bignum small_neg(-1);
  EXPECT_FALSE(small_neg.FitsIn<uint8_t>());
  EXPECT_FALSE(small_neg.FitsIn<uint16_t>());
  EXPECT_FALSE(small_neg.FitsIn<uint32_t>());
  EXPECT_FALSE(small_neg.FitsIn<uint64_t>());

  EXPECT_TRUE(small_neg.FitsIn<int8_t>());
  EXPECT_TRUE(small_neg.FitsIn<int16_t>());
  EXPECT_TRUE(small_neg.FitsIn<int32_t>());
  EXPECT_TRUE(small_neg.FitsIn<int64_t>());
}

TEST(BignumTest, FitsInUnsignedBoundsChecks) {
  const Bignum bn_u8max(kU8max);
  const Bignum bn_u8over(kU8max + 1);
  EXPECT_TRUE(bn_u8max.FitsIn<uint8_t>());
  EXPECT_TRUE(bn_u8max.FitsIn<uint16_t>());
  EXPECT_TRUE(bn_u8max.FitsIn<uint32_t>());
  EXPECT_FALSE(bn_u8over.FitsIn<uint8_t>());
  EXPECT_TRUE(bn_u8over.FitsIn<uint16_t>());
  EXPECT_TRUE(bn_u8over.FitsIn<uint32_t>());

  const Bignum bn_u16max(kU16max);
  const Bignum bn_u16over(kU16max + 1);
  EXPECT_FALSE(bn_u16max.FitsIn<uint8_t>());
  EXPECT_TRUE(bn_u16max.FitsIn<uint16_t>());
  EXPECT_TRUE(bn_u16max.FitsIn<uint32_t>());
  EXPECT_FALSE(bn_u16over.FitsIn<uint8_t>());
  EXPECT_FALSE(bn_u16over.FitsIn<uint16_t>());
  EXPECT_TRUE(bn_u16over.FitsIn<uint32_t>());

  const Bignum bn_u32max(kU32max);
  const Bignum bn_u32over(kU32max + 1);
  EXPECT_FALSE(bn_u32max.FitsIn<uint8_t>());
  EXPECT_FALSE(bn_u32max.FitsIn<uint16_t>());
  EXPECT_TRUE(bn_u32max.FitsIn<uint32_t>());
  EXPECT_FALSE(bn_u32over.FitsIn<uint8_t>());
  EXPECT_FALSE(bn_u32over.FitsIn<uint16_t>());
  EXPECT_FALSE(bn_u32over.FitsIn<uint32_t>());

  const Bignum bn_u64max(kU64max);
  EXPECT_TRUE(bn_u64max.FitsIn<uint64_t>());

  // 2^64, need to use string constructor.
  Bignum bn0 = *Bn("18446744073709551616");
  EXPECT_FALSE(bn0.FitsIn<uint64_t>());
}

TEST(BignumTest, FitsInSignedBoundsChecks) {
  const Bignum bn_i8max(kI8max);
  const Bignum bn_i8over(kI8max + 1);
  EXPECT_TRUE(bn_i8max.FitsIn<int8_t>());
  EXPECT_TRUE(bn_i8max.FitsIn<int16_t>());
  EXPECT_TRUE(bn_i8max.FitsIn<int32_t>());
  EXPECT_FALSE(bn_i8over.FitsIn<int8_t>());
  EXPECT_TRUE(bn_i8over.FitsIn<int16_t>());
  EXPECT_TRUE(bn_i8over.FitsIn<int32_t>());

  const Bignum bn_i16max(kI16max);
  const Bignum bn_i16over(kI16max + 1);
  EXPECT_FALSE(bn_i16max.FitsIn<int8_t>());
  EXPECT_TRUE(bn_i16max.FitsIn<int16_t>());
  EXPECT_TRUE(bn_i16max.FitsIn<int32_t>());
  EXPECT_FALSE(bn_i16over.FitsIn<int8_t>());
  EXPECT_FALSE(bn_i16over.FitsIn<int16_t>());
  EXPECT_TRUE(bn_i16over.FitsIn<int32_t>());

  const Bignum bn_i32max(kI32max);
  const Bignum bn_i32over(kI32max + 1);
  EXPECT_FALSE(bn_i32max.FitsIn<int8_t>());
  EXPECT_FALSE(bn_i32max.FitsIn<int16_t>());
  EXPECT_TRUE(bn_i32max.FitsIn<int32_t>());
  EXPECT_FALSE(bn_i32over.FitsIn<int8_t>());
  EXPECT_FALSE(bn_i32over.FitsIn<int16_t>());
  EXPECT_FALSE(bn_i32over.FitsIn<int32_t>());

  Bignum bn_i64max(kI64max);
  EXPECT_TRUE(bn_i64max.FitsIn<int64_t>());

  // 2^63, need to use string constructor.
  Bignum bn0 = *Bn("9223372036854775808");
  EXPECT_FALSE(bn0.FitsIn<int64_t>());

  const Bignum bn_i8min(kI8min);
  const Bignum bn_i8under(kI8min - 1);
  EXPECT_TRUE(bn_i8min.FitsIn<int8_t>());
  EXPECT_TRUE(bn_i8min.FitsIn<int16_t>());
  EXPECT_TRUE(bn_i8min.FitsIn<int32_t>());
  EXPECT_FALSE(bn_i8under.FitsIn<int8_t>());
  EXPECT_TRUE(bn_i8under.FitsIn<int16_t>());
  EXPECT_TRUE(bn_i8under.FitsIn<int32_t>());

  const Bignum bn_i16min(kI16min);
  const Bignum bn_i16under(kI16min - 1);
  EXPECT_FALSE(bn_i16min.FitsIn<int8_t>());
  EXPECT_TRUE(bn_i16min.FitsIn<int16_t>());
  EXPECT_TRUE(bn_i16min.FitsIn<int32_t>());
  EXPECT_FALSE(bn_i16under.FitsIn<int8_t>());
  EXPECT_FALSE(bn_i16under.FitsIn<int16_t>());
  EXPECT_TRUE(bn_i16under.FitsIn<int32_t>());

  const Bignum bn_i32min(kI32min);
  const Bignum bn_i32under(kI32min - 1);
  EXPECT_FALSE(bn_i32min.FitsIn<int8_t>());
  EXPECT_FALSE(bn_i32min.FitsIn<int16_t>());
  EXPECT_TRUE(bn_i32min.FitsIn<int32_t>());
  EXPECT_FALSE(bn_i32under.FitsIn<int8_t>());
  EXPECT_FALSE(bn_i32under.FitsIn<int16_t>());
  EXPECT_FALSE(bn_i32under.FitsIn<int32_t>());

  Bignum bn_i64min(kI64min);
  EXPECT_TRUE(bn_i64min.FitsIn<int64_t>());
}

TEST(BignumTest, FitsInBasicSanityChecks) {
  Bignum pos42(42);
  EXPECT_TRUE(pos42.FitsIn<int8_t>());
  EXPECT_TRUE(pos42.FitsIn<uint8_t>());
  EXPECT_TRUE(pos42.FitsIn<int16_t>());
  EXPECT_TRUE(pos42.FitsIn<uint16_t>());
  EXPECT_TRUE(pos42.FitsIn<int32_t>());
  EXPECT_TRUE(pos42.FitsIn<uint32_t>());
  EXPECT_TRUE(pos42.FitsIn<int64_t>());
  EXPECT_TRUE(pos42.FitsIn<uint64_t>());

  Bignum neg42(-42);
  EXPECT_TRUE(neg42.FitsIn<int8_t>());
  EXPECT_FALSE(neg42.FitsIn<uint8_t>());
  EXPECT_TRUE(neg42.FitsIn<int16_t>());
  EXPECT_FALSE(neg42.FitsIn<uint16_t>());
  EXPECT_TRUE(neg42.FitsIn<int32_t>());
  EXPECT_FALSE(neg42.FitsIn<uint32_t>());
  EXPECT_TRUE(neg42.FitsIn<int64_t>());
  EXPECT_FALSE(neg42.FitsIn<uint64_t>());
}

TEST(BignumTest, UnsignedCasting) {
  EXPECT_EQ(Bignum(300).Cast<uint8_t>(), static_cast<uint8_t>(300));

  // Negative to unsigned -> maximum value.
  Bignum bn1(-1);
  EXPECT_EQ(bn1.Cast<uint8_t>(), std::numeric_limits<uint8_t>::max());
  EXPECT_EQ(bn1.Cast<uint16_t>(), std::numeric_limits<uint16_t>::max());
  EXPECT_EQ(bn1.Cast<uint32_t>(), std::numeric_limits<uint32_t>::max());
  EXPECT_EQ(bn1.Cast<uint64_t>(), std::numeric_limits<uint64_t>::max());

  // Big values via decimal: 2^64, 2^128-1, 2^128
  Bignum bn2 = *Bn("18446744073709551616");
  EXPECT_EQ(bn2.Cast<uint64_t>(), 0);

  // 2^128 - 1 -> lower 64 bits = 2^64 - 1
  Bignum bn3 = *Bn("340282366920938463463374607431768211455");
  EXPECT_EQ(bn3.Cast<uint64_t>(), std::numeric_limits<uint64_t>::max());

  // 2^128 -> lower 64 bits = 0
  Bignum bn4 = *Bn("340282366920938463463374607431768211456");
  EXPECT_EQ(bn4.Cast<uint64_t>(), 0);
}

TEST(BignumTest, SignedCasting) {
  // In-range positives stay the same.
  Bignum bn0(127);
  EXPECT_EQ(bn0.Cast<int8_t>(), 127);

  // Positive overflow wraps into negative range.
  Bignum bn1(128);
  EXPECT_EQ(bn1.Cast<int8_t>(), -128);

  // In-range negatives stay the same.
  Bignum bn2(-128);
  EXPECT_EQ(bn2.Cast<int8_t>(), -128);

  // Negative overflow wraps into positive range.
  Bignum bn3(-129);
  EXPECT_EQ(bn3.Cast<int8_t>(), 127);

  // +2^63 over to -2^63 in signed int64.
  Bignum bn4 = *Bn("9223372036854775808");
  EXPECT_EQ(bn4.Cast<int64_t>(), std::numeric_limits<int64_t>::min());

  // +2^64 - 1 casts to -1.
  Bignum bn5 = *Bn("18446744073709551615");
  EXPECT_EQ(bn5.Cast<int64_t>(), -1);

  // -(2^63) - 1 casts to 2^63 - 3
  Bignum bn6 = *Bn("-9223372036854775809");
  EXPECT_EQ(bn6.Cast<int64_t>(), std::numeric_limits<int64_t>::max());
}

TEST(BignumTest, CastingLargeResidues) {
  // 2^80 + 0x1234 -> low 64 bits should be 0x1234.
  Bignum bn0 = *Bn("1208925819614629174710836");
  EXPECT_EQ(bn0.Cast<uint64_t>(), 0x1234);
  EXPECT_EQ(bn0.Cast<int64_t>(), 0x1234);

  // -(2^80 + 1) -> low 64 bits = 0xFFFFFFFFFFFFFFFF; signed = -1
  Bignum bn1 = *Bn("-1208925819614629174706177");
  EXPECT_EQ(bn1.Cast<uint64_t>(), std::numeric_limits<uint64_t>::max());
  EXPECT_EQ(bn1.Cast<int64_t>(), -1);
}

TEST(BignumTest, UnaryOperators) {
  EXPECT_EQ(+Bignum(0), Bignum(0));
  EXPECT_EQ(-Bignum(0), Bignum(0));
  EXPECT_EQ(+Bignum(+42), Bignum(+42));
  EXPECT_EQ(+Bignum(-42), Bignum(-42));
  EXPECT_EQ(-Bignum(+42), Bignum(-42));
  EXPECT_EQ(-Bignum(-17), Bignum(+17));
}

TEST(BignumTest, Addition) {
  // Basic combinations of signs
  EXPECT_EQ(Bignum(+5) + Bignum(+3), Bignum(+8));
  EXPECT_EQ(Bignum(-5) + Bignum(-3), Bignum(-8));
  EXPECT_EQ(Bignum(+5) + Bignum(-3), Bignum(+2));
  EXPECT_EQ(Bignum(-5) + Bignum(+3), Bignum(-2));

  // Identity and additive inverse
  EXPECT_EQ(Bignum(42) + Bignum(0), Bignum(42));
  EXPECT_EQ(Bignum(0) + Bignum(42), Bignum(42));
  EXPECT_EQ(Bignum(5) + Bignum(-5), Bignum(0));

  // Carry propagation
  const auto bn_u64max = Bignum(kU64max);
  EXPECT_EQ(bn_u64max + Bignum(1), *Bn("18446744073709551616"));
  EXPECT_EQ(bn_u64max + bn_u64max, *Bn("36893488147419103230"));

  // Aliasing (x += x)
  Bignum a = bn_u64max;
  a += a;
  EXPECT_EQ(a, *Bn("36893488147419103230"));
}

TEST(BignumTest, Subtraction) {
  // Basic combinations of signs
  EXPECT_EQ(Bignum(+5) - Bignum(+3), Bignum(+2));
  EXPECT_EQ(Bignum(+3) - Bignum(+5), Bignum(-2));
  EXPECT_EQ(Bignum(-5) - Bignum(-3), Bignum(-2));
  EXPECT_EQ(Bignum(+5) - Bignum(-3), Bignum(+8));
  EXPECT_EQ(Bignum(-5) - Bignum(+3), Bignum(-8));

  // Identity and subtracting to zero
  EXPECT_EQ(Bignum(42) - Bignum(0), Bignum(42));
  EXPECT_EQ(Bignum(0) - Bignum(42), Bignum(-42));
  EXPECT_EQ(Bignum(42) - Bignum(42), Bignum(0));

  // Borrow propagation
  const auto bn_u64max = Bignum(kU64max);
  const auto two_pow_64 = *Bn("18446744073709551616");
  EXPECT_EQ(two_pow_64 - Bignum(1), bn_u64max);

  // Aliasing (x -= x)
  Bignum a(100);
  a -= a;
  EXPECT_EQ(a, Bignum(0));
}

TEST(BignumTest, MixedOperations) {
  Bignum a(10), b(20), c(-5);
  EXPECT_EQ((a + b) - a, b);
  EXPECT_EQ((b - a) + a, b);
  EXPECT_EQ(a + c, Bignum(5));
  EXPECT_EQ(c - a, Bignum(-15));
}

TEST(BignumTest, LargeNumberArithmetic) {
  const auto two_pow_128_minus_1 =
      *Bn("340282366920938463463374607431768211455");
  const auto two_pow_128 = *Bn("340282366920938463463374607431768211456");
  const auto two_pow_64 = *Bn("18446744073709551616");

  // Test multi-bigit carry propagation: (2^128 - 1) + 1 = 2^128
  EXPECT_EQ(two_pow_128_minus_1 + Bignum(1), two_pow_128);

  // Test multi-bigit borrow propagation: 2^128 - 1 = (2^128 - 1)
  EXPECT_EQ(two_pow_128 - Bignum(1), two_pow_128_minus_1);

  // Subtraction resulting in a sign change with large numbers.
  const auto neg_two_pow_128_minus_1 =
      *Bn("-340282366920938463463374607431768211455");
  EXPECT_EQ(Bignum(1) - two_pow_128, neg_two_pow_128_minus_1);

  // Addition of large numbers with different signs (triggers subtraction).
  EXPECT_EQ(two_pow_128 + Bignum(-1), two_pow_128_minus_1);

  // Subtraction of large numbers with different signs (triggers addition).
  EXPECT_EQ(two_pow_128_minus_1 - Bignum(-1), two_pow_128);

  // Add two different large positive numbers: 2^128 + 2^64
  const auto sum_128_64 = *Bn("340282366920938463481821351505477763072");
  EXPECT_EQ(two_pow_128 + two_pow_64, sum_128_64);

  // Subtract two different large positive numbers: 2^128 - 2^64
  const auto diff_128_64 = *Bn("340282366920938463444927863358058659840");
  EXPECT_EQ(two_pow_128 - two_pow_64, diff_128_64);
}

TEST(BignumTest, LeftShift) {
  EXPECT_EQ((Bignum(1) << 0), Bignum(1));
  EXPECT_EQ((Bignum(1) << 1), Bignum(2));
  EXPECT_EQ((Bignum(1) << 63), *Bn("9223372036854775808"));
  EXPECT_EQ((Bignum(1) << 64), *Bn("18446744073709551616"));
  EXPECT_EQ((Bignum(-1) << 64), *Bn("-18446744073709551616"));

  const auto bn_u64max = Bignum(kU64max);
  const auto two_pow_128_minus_two_pow_64 =
      *Bn("340282366920938463444927863358058659840");
  EXPECT_EQ((bn_u64max << 64), two_pow_128_minus_two_pow_64);

  Bignum a(5);
  a <<= 2;
  EXPECT_EQ(a, Bignum(20));

  // Shifting zero or by zero amount.
  EXPECT_EQ((Bignum(0) << 100), Bignum(0));
  EXPECT_EQ((Bignum(123) << 0), Bignum(123));
}

TEST(BignumTest, RightShift) {
  EXPECT_EQ((Bignum(8) >> 0), Bignum(8));
  EXPECT_EQ((Bignum(8) >> 3), Bignum(1));
  EXPECT_EQ((Bignum(8) >> 4), Bignum(0));
  EXPECT_EQ((Bignum(7) >> 2), Bignum(1));
  EXPECT_EQ((Bignum(-7) >> 2), Bignum(-1));

  const auto two_pow_64 = *Bn("18446744073709551616");
  EXPECT_EQ((two_pow_64 >> 1), *Bn("9223372036854775808"));
  EXPECT_EQ((two_pow_64 >> 64), Bignum(1));
  EXPECT_EQ((two_pow_64 >> 65), Bignum(0));

  const auto u64max = Bignum(std::numeric_limits<uint64_t>::max());
  const auto two_pow_128_minus_1 =
      *Bn("340282366920938463463374607431768211455");
  EXPECT_EQ((two_pow_128_minus_1 >> 64), u64max);

  const auto two_pow_65_minus_1 = *Bn("36893488147419103231");
  EXPECT_EQ((two_pow_65_minus_1 >> 63), Bignum(3));

  Bignum b(20);
  b >>= 2;
  EXPECT_EQ(b, Bignum(5));

  // Shifting zero or by zero amount.
  EXPECT_EQ((Bignum(0) >> 100), Bignum(0));
  EXPECT_EQ((Bignum(123) >> 0), Bignum(123));

  // Shifting to zero.
  EXPECT_EQ((Bignum(100) >> 100), Bignum(0));
}

TEST(BignumTest, Multiplication) {
  // Zero
  EXPECT_EQ(Bignum(123) * Bignum(0), Bignum(0));
  EXPECT_EQ(Bignum(0) * Bignum(456), Bignum(0));

  // Identity
  EXPECT_EQ(Bignum(123) * Bignum(1), Bignum(123));
  EXPECT_EQ(Bignum(1) * Bignum(456), Bignum(456));
  EXPECT_EQ(Bignum(-123) * Bignum(1), Bignum(-123));

  // Signs
  EXPECT_EQ(Bignum(10) * Bignum(20), Bignum(200));
  EXPECT_EQ(Bignum(-10) * Bignum(20), Bignum(-200));
  EXPECT_EQ(Bignum(10) * Bignum(-20), Bignum(-200));
  EXPECT_EQ(Bignum(-10) * Bignum(-20), Bignum(200));

  // Simple carry
  const auto bn_u32max = Bignum(kU32max);
  EXPECT_EQ(bn_u32max * Bignum(2), *Bn("8589934590"));

  // 1x1 bigit fast path
  const auto bn_u64max = Bignum(kU64max);
  EXPECT_EQ(Bignum(2) * bn_u64max, *Bn("36893488147419103230"));

  // 1xN bigit multiplication
  const auto two_pow_128_minus_1 =
      *Bn("340282366920938463463374607431768211455");
  const auto res_1xN = *Bn("680564733841876926926749214863536422910");
  EXPECT_EQ(two_pow_128_minus_1 * Bignum(2), res_1xN);

  // Check that aliasing doesn't cause problems.
  Bignum a(100);
  a *= a;
  EXPECT_EQ(a, Bignum(10000));

  Bignum b = *Bn("10000000000000000000");  // > 64 bits
  b *= b;
  EXPECT_EQ(b, *Bn("100000000000000000000000000000000000000"));

  // Karatsuba threshold test
  // (2^128 - 1) * (2^128 - 1) = 2^256 - 2*2^128 + 1
  // This should trigger karatsuba if the threshold is low enough.
  EXPECT_EQ(two_pow_128_minus_1 * two_pow_128_minus_1,
            *Bn("115792089237316195423570985008687907852589419931798687112530"
                "834793049593217025"));

  // Karatsuba with uneven operands
  EXPECT_EQ(two_pow_128_minus_1 * bn_u64max,
            *Bn("6277101735386680763495507056286727952620534092958556749825"));
}

TEST(BignumTest, CountrZero) {
  EXPECT_EQ(countr_zero(Bignum(0)), 0);
  EXPECT_EQ(countr_zero(Bignum(1)), 0);
  EXPECT_EQ(countr_zero(Bignum(7)), 0);
  EXPECT_EQ(countr_zero(Bignum(-7)), 0);

  EXPECT_EQ(countr_zero(Bignum(2)), 1);
  EXPECT_EQ(countr_zero(Bignum(8)), 3);
  EXPECT_EQ(countr_zero(Bignum(10)), 1);  // 0b1010
  EXPECT_EQ(countr_zero(Bignum(12)), 2);  // 0b1100

  auto two_pow_64 = Bignum(1) << 64;
  EXPECT_EQ(countr_zero(two_pow_64), 64);

  auto large_shifted = Bignum(6) << 100;  // 0b110 << 100
  EXPECT_EQ(countr_zero(large_shifted), 101);

  auto neg_large_shifted = Bignum(-5) << 200;
  EXPECT_EQ(countr_zero(neg_large_shifted), 200);
}

TEST(BignumTest, is_bit_set) {
  EXPECT_FALSE(Bignum(0).is_bit_set(0));
  EXPECT_FALSE(Bignum(0).is_bit_set(100));

  // 5 = 0b101
  Bignum five(5);
  EXPECT_TRUE(five.is_bit_set(0));
  EXPECT_FALSE(five.is_bit_set(1));
  EXPECT_TRUE(five.is_bit_set(2));
  EXPECT_FALSE(five.is_bit_set(3));

  // Negative numbers should test the magnitude.
  Bignum neg_five(-5);
  EXPECT_TRUE(neg_five.is_bit_set(0));
  EXPECT_FALSE(neg_five.is_bit_set(1));
  EXPECT_TRUE(neg_five.is_bit_set(2));

  // Test edges of and across bigits.
  Bignum high_is_bit_set_63 = Bignum(1) << 63;
  EXPECT_FALSE(high_is_bit_set_63.is_bit_set(62));
  EXPECT_TRUE(high_is_bit_set_63.is_bit_set(63));
  EXPECT_FALSE(high_is_bit_set_63.is_bit_set(64));

  Bignum cross_bigit = (Bignum(1) << 100) + Bignum(1);
  EXPECT_TRUE(cross_bigit.is_bit_set(0));
  EXPECT_TRUE(cross_bigit.is_bit_set(100));
  EXPECT_FALSE(cross_bigit.is_bit_set(50));
  EXPECT_FALSE(cross_bigit.is_bit_set(1000));
}

TEST(BignumTest, Pow) {
  // Edge cases
  EXPECT_EQ(Bignum(0).Pow(0), Bignum(1));
  EXPECT_EQ(Bignum(123).Pow(0), Bignum(1));
  EXPECT_EQ(Bignum(0).Pow(123), Bignum(0));
  EXPECT_EQ(Bignum(1).Pow(12345), Bignum(1));

  // Negative base
  EXPECT_EQ(Bignum(-1).Pow(2), Bignum(1));
  EXPECT_EQ(Bignum(-1).Pow(3), Bignum(-1));
  EXPECT_EQ(Bignum(-2).Pow(2), Bignum(4));
  EXPECT_EQ(Bignum(-2).Pow(3), Bignum(-8));

  // Basic powers
  EXPECT_EQ(Bignum(2).Pow(10), Bignum(1024));
  EXPECT_EQ(Bignum(3).Pow(5), Bignum(243));
  EXPECT_EQ(Bignum(10).Pow(18), *Bn("1000000000000000000"));

  // Large exponent
  Bignum two_pow_100 = Bignum(1) << 100;
  EXPECT_EQ(Bignum(2).Pow(100), two_pow_100);

  // Large base
  Bignum ten_pow_19 = *Bn("10000000000000000000");
  Bignum ten_pow_38 = *Bn("100000000000000000000000000000000000000");
  EXPECT_EQ(ten_pow_19.Pow(2), ten_pow_38);
}

TEST(BignumTest, SetZero) {
  Bignum a(123);
  a.set_zero();
  EXPECT_TRUE(a.is_zero());

  Bignum b(-456);
  b.set_zero();
  EXPECT_EQ(b, Bignum(0));
}

TEST(BignumTest, SetNegativeSetPositive) {
  Bignum a(42);
  a.set_negative();
  EXPECT_TRUE(a.is_negative());
  EXPECT_EQ(a, Bignum(-42));

  a.set_negative(false);
  EXPECT_FALSE(a.is_negative());
  EXPECT_EQ(a, Bignum(42));

  // set_negative() has no effect on zero.
  Bignum b(0);
  b.set_negative();
  EXPECT_FALSE(b.is_negative());
  EXPECT_EQ(b, Bignum(0));
}

TEST(BignumTest, Comparisons) {
  EXPECT_EQ(*Bn("123"), Bignum(123));
  EXPECT_EQ(*Bn("-123"), Bignum(-123));
  EXPECT_NE(*Bn("123"), Bignum(-123));
  EXPECT_EQ(Bignum(0), Bignum(0));
  EXPECT_NE(Bignum(0), Bignum(1));

  // Positive vs Positive
  EXPECT_LT(Bignum(100), Bignum(200));
  EXPECT_GT(Bignum(200), Bignum(100));
  EXPECT_LE(Bignum(100), Bignum(200));
  EXPECT_GE(Bignum(200), Bignum(100));

  // Negative vs Negative
  EXPECT_LT(Bignum(-200), Bignum(-100));
  EXPECT_GT(Bignum(-100), Bignum(-200));
  EXPECT_GE(Bignum(-100), Bignum(-200));
  EXPECT_LE(Bignum(-200), Bignum(-100));

  // Positive vs Negative
  EXPECT_LT(Bignum(-10), Bignum(10));
  EXPECT_GT(Bignum(10), Bignum(-10));

  // Zero
  EXPECT_LT(Bignum(-1), Bignum(0));
  EXPECT_LT(Bignum(0), Bignum(1));
  EXPECT_GT(Bignum(0), Bignum(-1));
  EXPECT_GT(Bignum(1), Bignum(0));

  // Multi-bigit
  const auto two_pow_64 = *Bn("18446744073709551616");
  EXPECT_LT(Bignum(0), two_pow_64);
  EXPECT_GT(two_pow_64, Bignum(0));
  EXPECT_LT(Bignum(-1), two_pow_64);
  EXPECT_GT(two_pow_64, Bignum(-1));

  EXPECT_LE(Bignum(100), Bignum(200));
  EXPECT_LE(Bignum(100), Bignum(100));
  EXPECT_LE(Bignum(-200), Bignum(-100));
  EXPECT_LE(Bignum(-100), Bignum(-100));
  EXPECT_GT(Bignum(200), Bignum(100));

  EXPECT_GE(Bignum(200), Bignum(100));
  EXPECT_GE(Bignum(100), Bignum(100));
  EXPECT_GE(Bignum(-100), Bignum(-200));
  EXPECT_GE(Bignum(-100), Bignum(-100));
  EXPECT_LT(Bignum(100), Bignum(200));

  EXPECT_LE(Bignum(0), Bignum(0));
  EXPECT_GE(Bignum(0), Bignum(0));
}

// RAII wrapper for OpenSSL BIGNUM
class OpenSSLBignum {
 public:
  OpenSSLBignum() : bn_(BN_new()) {}

  // Construct from a decimal number in a string.
  //
  // We take decimal as a string so that it's explicitly zero-terminated.
  explicit OpenSSLBignum(const std::string& decimal) : bn_(BN_new()) {
    BN_dec2bn(&bn_, decimal.data());
  }

  explicit OpenSSLBignum(uint64_t value) : bn_(BN_new()) {
    BN_set_word(bn_, value);
  }

  ~OpenSSLBignum() { BN_free(bn_); }

  OpenSSLBignum(OpenSSLBignum&& other) noexcept : bn_(other.bn_) {
    other.bn_ = nullptr;
  }

  OpenSSLBignum& operator=(OpenSSLBignum&& other) noexcept {
    if (this != &other) {
      BN_free(bn_);
      bn_ = other.bn_;
      other.bn_ = nullptr;
    }
    return *this;
  }

  OpenSSLBignum(const OpenSSLBignum& other) : bn_(BN_dup(other.bn_)) {}

  OpenSSLBignum& operator=(const OpenSSLBignum& other) {
    if (this != &other) {
      BN_copy(bn_, other.bn_);
    }
    return *this;
  }

  BIGNUM* get() const { return bn_; }

 private:
  BIGNUM* bn_;
};

// Power of two for fast modulo.
const int kRandomBignumCount = 128;

static std::vector<std::string> GenerateRandomNumberStrings(
    absl::BitGenRef bitgen, int bits) {
  std::vector<std::string> numbers;
  numbers.reserve(kRandomBignumCount);

  for (int i = 0; i < kRandomBignumCount; ++i) {
    std::string num;

    // Generate approximately `bits` worth of decimal digits
    int decimal_digits = (bits * 3) / 10;  // log10(2^bits) ≈ bits * 0.301

    // First digit can't be zero
    absl::StrAppend(&num, absl::StrFormat("%d", absl::Uniform(bitgen, 1, 9)));
    for (int j = 1; j < decimal_digits; ++j) {
      num += absl::Uniform(bitgen, '0', '9');
    }

    numbers.push_back(num);
  }

  return numbers;
}

// Basic correctness test to ensure OpenSSL integration is working
TEST(BignumTest, OpenSSLIntegration) {
  OpenSSLBignum a(123);
  OpenSSLBignum b(456);
  OpenSSLBignum result;

  BN_add(result.get(), a.get(), b.get());

  char* str = BN_bn2dec(result.get());
  EXPECT_STREQ(str, "579");
  OPENSSL_free(str);
}

TEST(BignumTest, ResultsMatch) {
  // Test that and OpenSSL produce the same results
  const Bignum w_a(12345);
  const Bignum w_b(67890);
  Bignum w_result = w_a + w_b;

  const OpenSSLBignum ssl_a(12345);
  const OpenSSLBignum ssl_b(67890);
  OpenSSLBignum ssl_result;
  BN_add(ssl_result.get(), ssl_a.get(), ssl_b.get());

  char* ssl_str = BN_bn2dec(ssl_result.get());
  std::string w_str = absl::StrFormat("%v", w_result);

  EXPECT_EQ(w_str, std::string(ssl_str));
  OPENSSL_free(ssl_str);
}

// Different number sizes for benchmarking.
enum class NumberSizeClass : uint32_t {
  kSmall = 64,
  kMedium = 256,
  kLarge = 1024,
  kHuge = 4096,
  kMega = 18000
};

std::vector<std::string> RandomNumberStrings(absl::BitGenRef bitgen,
                                             NumberSizeClass size_class) {
  return GenerateRandomNumberStrings(bitgen, static_cast<int>(size_class));
}

class VsOpenSSLTest
    : public TestWithParam<std::pair<NumberSizeClass, NumberSizeClass>> {
 protected:
  std::vector<std::pair<std::string, std::string>> Numbers() {
    auto numbers0 = RandomNumberStrings(bitgen_, GetParam().first);
    auto numbers1 = RandomNumberStrings(bitgen_, GetParam().second);
    ABSL_CHECK_EQ(numbers0.size(), numbers1.size());

    std::vector<std::pair<std::string, std::string>> numbers;
    numbers.reserve(numbers0.size());

    for (size_t i = 0; i < numbers0.size(); ++i) {
      numbers.emplace_back(numbers0[i], numbers1[i]);
    }
    return numbers;
  }

 private:
  absl::BitGen bitgen_;
};

TEST_P(VsOpenSSLTest, MultiplyCorrect) {
  // Test that multiplication produces the same results as OpenSSL.
  BN_CTX* ctx = BN_CTX_new();
  for (const auto& [a, b] : Numbers()) {
    const Bignum bn_a = *Bignum::FromString(a);
    const Bignum bn_b = *Bignum::FromString(b);
    const Bignum bn_result = bn_a * bn_b;

    const OpenSSLBignum ssl_a(a);
    const OpenSSLBignum ssl_b(b);
    OpenSSLBignum ssl_result;
    BN_mul(ssl_result.get(), ssl_a.get(), ssl_b.get(), ctx);

    // Compare string representations
    char* ssl_str = BN_bn2dec(ssl_result.get());
    std::string bn_str = absl::StrFormat("%v", bn_result);

    EXPECT_EQ(bn_str, std::string(ssl_str))
        << "Mismatch for multiplication"
        << "\nBignum result: " << bn_str.substr(0, 100) << "..."
        << "\nOpenSSL result: " << std::string(ssl_str).substr(0, 100) << "...";
    OPENSSL_free(ssl_str);
  }
  BN_CTX_free(ctx);
}

TEST_P(VsOpenSSLTest, AdditionCorrect) {
  // Test that addition produces correct results by comparing to OpenSSL.
  for (const auto& [a, b] : Numbers()) {
    const Bignum bn_a = *Bignum::FromString(a);
    const Bignum bn_b = *Bignum::FromString(b);

    const Bignum bn_result = bn_a + bn_b;

    const OpenSSLBignum ssl_a(a);
    const OpenSSLBignum ssl_b(b);
    OpenSSLBignum ssl_result;
    BN_add(ssl_result.get(), ssl_a.get(), ssl_b.get());

    // Compare string representations
    char* ssl_str = BN_bn2dec(ssl_result.get());
    std::string bn_str = absl::StrFormat("%v", bn_result);

    EXPECT_EQ(bn_str, std::string(ssl_str))
        << "Mismatch for addition"
        << "\nBignum result: " << bn_str.substr(0, 100) << "..."
        << "\nOpenSSL result: " << std::string(ssl_str).substr(0, 100) << "...";
    OPENSSL_free(ssl_str);
  }
}

TEST_P(VsOpenSSLTest, SubtractionCorrect) {
  // Test that subtraction produces correct results by comparing to
  // OpenSSL.
  for (const auto& [a, b] : Numbers()) {
    const Bignum bn_a = *Bignum::FromString(a);
    const Bignum bn_b = *Bignum::FromString(b);

    const Bignum bn_result = bn_a - bn_b;

    const OpenSSLBignum ssl_a(a);
    const OpenSSLBignum ssl_b(b);
    OpenSSLBignum ssl_result;
    BN_sub(ssl_result.get(), ssl_a.get(), ssl_b.get());

    // Compare string representations
    char* ssl_str = BN_bn2dec(ssl_result.get());
    std::string bn_str = absl::StrFormat("%v", bn_result);

    EXPECT_EQ(bn_str, std::string(ssl_str))
        << "Mismatch for addition"
        << "\nBignum result: " << bn_str.substr(0, 100) << "..."
        << "\nOpenSSL result: " << std::string(ssl_str).substr(0, 100) << "...";
    OPENSSL_free(ssl_str);
  }
}

// clang-format off
INSTANTIATE_TEST_SUITE_P(
    VsOpenSSL, VsOpenSSLTest, ::testing::Values(
      std::make_pair(NumberSizeClass::kSmall, NumberSizeClass::kSmall),
      std::make_pair(NumberSizeClass::kSmall, NumberSizeClass::kHuge),
      std::make_pair(NumberSizeClass::kHuge, NumberSizeClass::kSmall),
      std::make_pair(NumberSizeClass::kMedium, NumberSizeClass::kMedium),
      std::make_pair(NumberSizeClass::kLarge, NumberSizeClass::kLarge),
      std::make_pair(NumberSizeClass::kHuge, NumberSizeClass::kHuge),
      std::make_pair(NumberSizeClass::kMega, NumberSizeClass::kMega)));
// clang-format on

// TODO: Enable once benchmark is integrated.
#if 0

template <typename BinaryOp>
void BignumBinaryOpBenchmark(benchmark::State& state,
                             const std::vector<std::string>& number_strings,
                             BinaryOp op) {
  std::vector<Bignum> numbers;
  for (const auto& str : number_strings) {
    numbers.push_back(*Bignum::FromString(str));
  }

  Bignum result;
  size_t idx = 0;
  for (auto _ : state) {
    const Bignum& a = numbers[(idx + 0) % kRandomBignumCount];
    const Bignum& b = numbers[(idx + 1) % kRandomBignumCount];
    result = op(a, b);
    benchmark::DoNotOptimize(result);
    ++idx;
  }
}

void BignumPowBenchmark(benchmark::State& state,
                        const std::vector<std::string>& number_strings,
                        int exponent) {
  std::vector<Bignum> numbers;
  for (const auto& str : number_strings) {
    numbers.push_back(*Bignum::FromString(str));
  }

  Bignum result;
  size_t idx = 0;
  for (auto _ : state) {
    const Bignum& base = numbers[(idx + 0) % kRandomBignumCount];
    result = base.Pow(exponent);
    benchmark::DoNotOptimize(result);
    ++idx;
  }
}

template <typename BinaryOp>
void OpenSSLBinaryOpBenchmark(benchmark::State& state,
                              const std::vector<std::string>& number_strings,
                              BinaryOp op) {
  std::vector<OpenSSLBignum> numbers;
  for (const auto& str : number_strings) {
    numbers.emplace_back(str);
  }

  size_t idx = 0;

  for (auto _ : state) {
    OpenSSLBignum result;
    const OpenSSLBignum& a = numbers[(idx + 0) % kRandomBignumCount];
    const OpenSSLBignum& b = numbers[(idx + 1) % kRandomBignumCount];
    op(result.get(), a.get(), b.get());
    benchmark::DoNotOptimize(result.get());
    ++idx;
  }
}

template <typename MulOp>
void OpenSSLMulOpBenchmark(benchmark::State& state,
                           const std::vector<std::string>& number_strings,
                           MulOp op) {
  std::vector<OpenSSLBignum> numbers;
  for (const auto& str : number_strings) {
    numbers.emplace_back(str);
  }

  BN_CTX* ctx = BN_CTX_new();
  size_t idx = 0;

  for (auto _ : state) {
    OpenSSLBignum result;
    const OpenSSLBignum& a = numbers[(idx + 0) % kRandomBignumCount];
    const OpenSSLBignum& b = numbers[(idx + 1) % kRandomBignumCount];
    op(result.get(), a.get(), b.get(), ctx);
    benchmark::DoNotOptimize(result.get());
    ++idx;
  }

  BN_CTX_free(ctx);
}

void OpenSSLPowBenchmark(benchmark::State& state,
                         const std::vector<std::string>& number_strings,
                         int exponent) {
  std::vector<OpenSSLBignum> numbers;
  for (const auto& str : number_strings) {
    numbers.emplace_back(str);
  }

  const OpenSSLBignum exp(exponent);
  BN_CTX* ctx = BN_CTX_new();
  size_t idx = 0;

  for (auto _ : state) {
    OpenSSLBignum result;
    const OpenSSLBignum& base = numbers[(idx + 0) % kRandomBignumCount];
    BN_exp(result.get(), base.get(), exp.get(), ctx);
    benchmark::DoNotOptimize(result.get());
    ++idx;
  }

  BN_CTX_free(ctx);
}

std::vector<std::string> SmallNumbers(absl::BitGenRef bitgen) {
  return RandomNumberStrings(bitgen, NumberSizeClass::kSmall);
}

std::vector<std::string> MediumNumbers(absl::BitGenRef bitgen) {
  return RandomNumberStrings(bitgen, NumberSizeClass::kMedium);
}

std::vector<std::string> LargeNumbers(absl::BitGenRef bitgen) {
  return RandomNumberStrings(bitgen, NumberSizeClass::kLarge);
}

std::vector<std::string> HugeNumbers(absl::BitGenRef bitgen) {
  return RandomNumberStrings(bitgen, NumberSizeClass::kHuge);
}

std::vector<std::string> MegaNumbers(absl::BitGenRef bitgen) {
  return RandomNumberStrings(bitgen, NumberSizeClass::kMega);
}

void BM_Bignum_AddSmall(benchmark::State& state) {
  std::mt19937_64 bitgen;
  BignumBinaryOpBenchmark(state, SmallNumbers(bitgen), std::plus<Bignum>{});
}
BENCHMARK(BM_Bignum_AddSmall);

void BM_Bignum_AddMedium(benchmark::State& state) {
  std::mt19937_64 bitgen;
  BignumBinaryOpBenchmark(state, MediumNumbers(bitgen), std::plus<Bignum>{});
}
BENCHMARK(BM_Bignum_AddMedium);

void BM_Bignum_AddLarge(benchmark::State& state) {
  std::mt19937_64 bitgen;
  BignumBinaryOpBenchmark(state, LargeNumbers(bitgen), std::plus<Bignum>{});
}
BENCHMARK(BM_Bignum_AddLarge);

void BM_Bignum_AddHuge(benchmark::State& state) {
  std::mt19937_64 bitgen;
  BignumBinaryOpBenchmark(state, HugeNumbers(bitgen), std::plus<Bignum>{});
}
BENCHMARK(BM_Bignum_AddHuge);

void BM_Bignum_AddMega(benchmark::State& state) {
  std::mt19937_64 bitgen;
  BignumBinaryOpBenchmark(state, MegaNumbers(bitgen), std::plus<Bignum>{});
}
BENCHMARK(BM_Bignum_AddMega);

void BM_OpenSSL_AddSmall(benchmark::State& state) {
  std::mt19937_64 bitgen;
  OpenSSLBinaryOpBenchmark(state, SmallNumbers(bitgen), BN_add);
}
BENCHMARK(BM_OpenSSL_AddSmall);

void BM_OpenSSL_AddMedium(benchmark::State& state) {
  std::mt19937_64 bitgen;
  OpenSSLBinaryOpBenchmark(state, MediumNumbers(bitgen), BN_add);
}
BENCHMARK(BM_OpenSSL_AddMedium);

void BM_OpenSSL_AddLarge(benchmark::State& state) {
  std::mt19937_64 bitgen;
  OpenSSLBinaryOpBenchmark(state, LargeNumbers(bitgen), BN_add);
}
BENCHMARK(BM_OpenSSL_AddLarge);

void BM_OpenSSL_AddHuge(benchmark::State& state) {
  std::mt19937_64 bitgen;
  OpenSSLBinaryOpBenchmark(state, HugeNumbers(bitgen), BN_add);
}
BENCHMARK(BM_OpenSSL_AddHuge);

void BM_OpenSSL_AddMega(benchmark::State& state) {
  std::mt19937_64 bitgen;
  OpenSSLBinaryOpBenchmark(state, MegaNumbers(bitgen), BN_add);
}
BENCHMARK(BM_OpenSSL_AddMega);

void BM_Bignum_MulSmall(benchmark::State& state) {
  std::mt19937_64 bitgen;
  BignumBinaryOpBenchmark(state, SmallNumbers(bitgen),
                          std::multiplies<Bignum>{});
}
BENCHMARK(BM_Bignum_MulSmall);

void BM_Bignum_MulMedium(benchmark::State& state) {
  std::mt19937_64 bitgen;
  BignumBinaryOpBenchmark(state, MediumNumbers(bitgen),
                          std::multiplies<Bignum>{});
}
BENCHMARK(BM_Bignum_MulMedium);

void BM_Bignum_MulLarge(benchmark::State& state) {
  std::mt19937_64 bitgen;
  BignumBinaryOpBenchmark(state, LargeNumbers(bitgen),
                          std::multiplies<Bignum>{});
}
BENCHMARK(BM_Bignum_MulLarge);

void BM_Bignum_MulHuge(benchmark::State& state) {
  std::mt19937_64 bitgen;
  BignumBinaryOpBenchmark(state, HugeNumbers(bitgen),
                          std::multiplies<Bignum>{});
}
BENCHMARK(BM_Bignum_MulHuge);

void BM_Bignum_MulMega(benchmark::State& state) {
  std::mt19937_64 bitgen;
  BignumBinaryOpBenchmark(state, MegaNumbers(bitgen),
                          std::multiplies<Bignum>{});
}
BENCHMARK(BM_Bignum_MulMega);

void BM_OpenSSL_MulSmall(benchmark::State& state) {
  std::mt19937_64 bitgen;
  OpenSSLMulOpBenchmark(state, SmallNumbers(bitgen), BN_mul);
}
BENCHMARK(BM_OpenSSL_MulSmall);

void BM_OpenSSL_MulMedium(benchmark::State& state) {
  std::mt19937_64 bitgen;
  OpenSSLMulOpBenchmark(state, MediumNumbers(bitgen), BN_mul);
}
BENCHMARK(BM_OpenSSL_MulMedium);

void BM_OpenSSL_MulLarge(benchmark::State& state) {
  std::mt19937_64 bitgen;
  OpenSSLMulOpBenchmark(state, LargeNumbers(bitgen), BN_mul);
}
BENCHMARK(BM_OpenSSL_MulLarge);

void BM_OpenSSL_MulHuge(benchmark::State& state) {
  std::mt19937_64 bitgen;
  OpenSSLMulOpBenchmark(state, HugeNumbers(bitgen), BN_mul);
}
BENCHMARK(BM_OpenSSL_MulHuge);

void BM_OpenSSL_MulMega(benchmark::State& state) {
  std::mt19937_64 bitgen;
  OpenSSLMulOpBenchmark(state, MegaNumbers(bitgen), BN_mul);
}
BENCHMARK(BM_OpenSSL_MulMega);

void BM_Bignum_PowSmall(benchmark::State& state) {
  std::mt19937_64 bitgen;
  BignumPowBenchmark(state, SmallNumbers(bitgen), 20);
}
BENCHMARK(BM_Bignum_PowSmall);

void BM_Bignum_PowMedium(benchmark::State& state) {
  std::mt19937_64 bitgen;
  BignumPowBenchmark(state, MediumNumbers(bitgen), 10);
}
BENCHMARK(BM_Bignum_PowMedium);

void BM_OpenSSL_PowSmall(benchmark::State& state) {
  std::mt19937_64 bitgen;
  OpenSSLPowBenchmark(state, SmallNumbers(bitgen), 20);
}
BENCHMARK(BM_OpenSSL_PowSmall);

void BM_OpenSSL_PowMedium(benchmark::State& state) {
  std::mt19937_64 bitgen;
  OpenSSLPowBenchmark(state, MediumNumbers(bitgen), 10);
}
BENCHMARK(BM_OpenSSL_PowMedium);
#endif

}  // namespace exactfloat_internal
