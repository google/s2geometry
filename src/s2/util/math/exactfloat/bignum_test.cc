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

#include <memory>
#include <random>
#include <string>
#include <vector>

#if 0
#include "absl/base/no_destructor.h"
#include "absl/strings/string_view.h"
#include "benchmark/benchmark.h"
#include "openssl/bn.h"
#include "openssl/crypto.h"
#endif

#include "gtest/gtest.h"

const uint64_t u8max = std::numeric_limits<uint8_t>::max();
const uint64_t u16max = std::numeric_limits<uint16_t>::max();
const uint64_t u32max = std::numeric_limits<uint32_t>::max();
const uint64_t u64max = std::numeric_limits<uint64_t>::max();

const int64_t i8max = std::numeric_limits<int8_t>::max();
const int64_t i16max = std::numeric_limits<int16_t>::max();
const int64_t i32max = std::numeric_limits<int32_t>::max();
const int64_t i64max = std::numeric_limits<int64_t>::max();

const int64_t i8min = std::numeric_limits<int8_t>::min();
const int64_t i16min = std::numeric_limits<int16_t>::min();
const int64_t i32min = std::numeric_limits<int32_t>::min();
const int64_t i64min = std::numeric_limits<int64_t>::min();

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

TEST(BignumTest, ZeroAlwaysCompatible) {
  const Bignum zero(0);
  EXPECT_TRUE(zero.Compatible<int8_t>());
  EXPECT_TRUE(zero.Compatible<uint8_t>());
  EXPECT_TRUE(zero.Compatible<int16_t>());
  EXPECT_TRUE(zero.Compatible<uint16_t>());
  EXPECT_TRUE(zero.Compatible<int32_t>());
  EXPECT_TRUE(zero.Compatible<uint32_t>());
  EXPECT_TRUE(zero.Compatible<int64_t>());
  EXPECT_TRUE(zero.Compatible<uint64_t>());
  EXPECT_TRUE(zero.Compatible<absl::int128>());
  EXPECT_TRUE(zero.Compatible<absl::uint128>());
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

TEST(BignumTest, NegativeOnlyCompatibleSigned) {
  const Bignum small_neg(-1);
  EXPECT_FALSE(small_neg.Compatible<uint8_t>());
  EXPECT_FALSE(small_neg.Compatible<uint16_t>());
  EXPECT_FALSE(small_neg.Compatible<uint32_t>());
  EXPECT_FALSE(small_neg.Compatible<uint64_t>());
  EXPECT_FALSE(small_neg.Compatible<absl::uint128>());

  EXPECT_TRUE(small_neg.Compatible<int8_t>());
  EXPECT_TRUE(small_neg.Compatible<int16_t>());
  EXPECT_TRUE(small_neg.Compatible<int32_t>());
  EXPECT_TRUE(small_neg.Compatible<int64_t>());
  EXPECT_TRUE(small_neg.Compatible<absl::int128>());
}

TEST(BignumTest, CompatibleUnsignedBoundsChecks) {
  const Bignum bn_u8max(u8max);
  const Bignum bn_u8over(u8max + 1);
  EXPECT_TRUE(bn_u8max.Compatible<uint8_t>());
  EXPECT_TRUE(bn_u8max.Compatible<uint16_t>());
  EXPECT_TRUE(bn_u8max.Compatible<uint32_t>());
  EXPECT_FALSE(bn_u8over.Compatible<uint8_t>());
  EXPECT_TRUE(bn_u8over.Compatible<uint16_t>());
  EXPECT_TRUE(bn_u8over.Compatible<uint32_t>());

  const Bignum bn_u16max(u16max);
  const Bignum bn_u16over(u16max + 1);
  EXPECT_FALSE(bn_u16max.Compatible<uint8_t>());
  EXPECT_TRUE(bn_u16max.Compatible<uint16_t>());
  EXPECT_TRUE(bn_u16max.Compatible<uint32_t>());
  EXPECT_FALSE(bn_u16over.Compatible<uint8_t>());
  EXPECT_FALSE(bn_u16over.Compatible<uint16_t>());
  EXPECT_TRUE(bn_u16over.Compatible<uint32_t>());

  const Bignum bn_u32max(u32max);
  const Bignum bn_u32over(u32max + 1);
  EXPECT_FALSE(bn_u32max.Compatible<uint8_t>());
  EXPECT_FALSE(bn_u32max.Compatible<uint16_t>());
  EXPECT_TRUE(bn_u32max.Compatible<uint32_t>());
  EXPECT_FALSE(bn_u32over.Compatible<uint8_t>());
  EXPECT_FALSE(bn_u32over.Compatible<uint16_t>());
  EXPECT_FALSE(bn_u32over.Compatible<uint32_t>());

  const Bignum bn_u64max(u64max);
  EXPECT_TRUE(bn_u64max.Compatible<uint64_t>());

  // 2^64, need to use string constructor.
  Bignum bn0 = *Bn("18446744073709551616");
  EXPECT_FALSE(bn0.Compatible<uint64_t>());
  EXPECT_TRUE(bn0.Compatible<absl::uint128>());

  // (2^128 - 1) fits in absl::uint128.
  Bignum bn1 = *Bn("340282366920938463463374607431768211455");
  EXPECT_TRUE(bn1.Compatible<absl::uint128>());

  // 2^128 does not fit in absl::uint128.
  Bignum bn2 = *Bn("340282366920938463463374607431768211456");
  EXPECT_FALSE(bn2.Compatible<absl::uint128>());
}

TEST(BignumTest, CompatibleSignedBoundsChecks) {
  const Bignum bn_i8max(i8max);
  const Bignum bn_i8over(i8max + 1);
  EXPECT_TRUE(bn_i8max.Compatible<int8_t>());
  EXPECT_TRUE(bn_i8max.Compatible<int16_t>());
  EXPECT_TRUE(bn_i8max.Compatible<int32_t>());
  EXPECT_FALSE(bn_i8over.Compatible<int8_t>());
  EXPECT_TRUE(bn_i8over.Compatible<int16_t>());
  EXPECT_TRUE(bn_i8over.Compatible<int32_t>());

  const Bignum bn_i16max(i16max);
  const Bignum bn_i16over(i16max + 1);
  EXPECT_FALSE(bn_i16max.Compatible<int8_t>());
  EXPECT_TRUE(bn_i16max.Compatible<int16_t>());
  EXPECT_TRUE(bn_i16max.Compatible<int32_t>());
  EXPECT_FALSE(bn_i16over.Compatible<int8_t>());
  EXPECT_FALSE(bn_i16over.Compatible<int16_t>());
  EXPECT_TRUE(bn_i16over.Compatible<int32_t>());

  const Bignum bn_i32max(i32max);
  const Bignum bn_i32over(i32max + 1);
  EXPECT_FALSE(bn_i32max.Compatible<int8_t>());
  EXPECT_FALSE(bn_i32max.Compatible<int16_t>());
  EXPECT_TRUE(bn_i32max.Compatible<int32_t>());
  EXPECT_FALSE(bn_i32over.Compatible<int8_t>());
  EXPECT_FALSE(bn_i32over.Compatible<int16_t>());
  EXPECT_FALSE(bn_i32over.Compatible<int32_t>());

  Bignum bn_i64max(i64max);
  EXPECT_TRUE(bn_i64max.Compatible<int64_t>());

  // 2^63, need to use string constructor.
  Bignum bn0 = *Bn("9223372036854775808");
  EXPECT_FALSE(bn0.Compatible<int64_t>());

  const Bignum bn_i8min(i8min);
  const Bignum bn_i8under(i8min - 1);
  EXPECT_TRUE(bn_i8min.Compatible<int8_t>());
  EXPECT_TRUE(bn_i8min.Compatible<int16_t>());
  EXPECT_TRUE(bn_i8min.Compatible<int32_t>());
  EXPECT_FALSE(bn_i8under.Compatible<int8_t>());
  EXPECT_TRUE(bn_i8under.Compatible<int16_t>());
  EXPECT_TRUE(bn_i8under.Compatible<int32_t>());

  const Bignum bn_i16min(i16min);
  const Bignum bn_i16under(i16min - 1);
  EXPECT_FALSE(bn_i16min.Compatible<int8_t>());
  EXPECT_TRUE(bn_i16min.Compatible<int16_t>());
  EXPECT_TRUE(bn_i16min.Compatible<int32_t>());
  EXPECT_FALSE(bn_i16under.Compatible<int8_t>());
  EXPECT_FALSE(bn_i16under.Compatible<int16_t>());
  EXPECT_TRUE(bn_i16under.Compatible<int32_t>());

  const Bignum bn_i32min(i32min);
  const Bignum bn_i32under(i32min - 1);
  EXPECT_FALSE(bn_i32min.Compatible<int8_t>());
  EXPECT_FALSE(bn_i32min.Compatible<int16_t>());
  EXPECT_TRUE(bn_i32min.Compatible<int32_t>());
  EXPECT_FALSE(bn_i32under.Compatible<int8_t>());
  EXPECT_FALSE(bn_i32under.Compatible<int16_t>());
  EXPECT_FALSE(bn_i32under.Compatible<int32_t>());

  Bignum bn_i64min(i64min);
  EXPECT_TRUE(bn_i64min.Compatible<int64_t>());

  // -(2^63) - 1 doesn't fit in int64_t.
  Bignum b0 = *Bn("-9223372036854775809");
  EXPECT_FALSE(b0.Compatible<int64_t>());

  // Exact min and max of signed 128.
  Bignum bn_s128min = *Bn("-170141183460469231731687303715884105728");
  Bignum bn_s128max = *Bn("170141183460469231731687303715884105727");
  EXPECT_TRUE(bn_s128min.Compatible<absl::int128>());
  EXPECT_TRUE(bn_s128max.Compatible<absl::int128>());

  // +2^127 does not fit in signed 128, but does in unsigned 128.
  Bignum bn1 = *Bn("170141183460469231731687303715884105728");
  EXPECT_FALSE(bn1.Compatible<absl::int128>());
  EXPECT_TRUE(bn1.Compatible<absl::uint128>());

  // Below min: -(2^127) - 1 should not fit.
  Bignum bn2 = *Bn("-170141183460469231731687303715884105729");
  EXPECT_FALSE(bn2.Compatible<absl::int128>());
}

TEST(BignumTest, CompatibleBasicSanityChecks) {
  Bignum pos42(42);
  EXPECT_TRUE(pos42.Compatible<int8_t>());
  EXPECT_TRUE(pos42.Compatible<uint8_t>());
  EXPECT_TRUE(pos42.Compatible<int16_t>());
  EXPECT_TRUE(pos42.Compatible<uint16_t>());
  EXPECT_TRUE(pos42.Compatible<int32_t>());
  EXPECT_TRUE(pos42.Compatible<uint32_t>());
  EXPECT_TRUE(pos42.Compatible<int64_t>());
  EXPECT_TRUE(pos42.Compatible<uint64_t>());
  EXPECT_TRUE(pos42.Compatible<absl::int128>());
  EXPECT_TRUE(pos42.Compatible<absl::uint128>());

  Bignum neg42(-42);
  EXPECT_TRUE(neg42.Compatible<int8_t>());
  EXPECT_FALSE(neg42.Compatible<uint8_t>());
  EXPECT_TRUE(neg42.Compatible<int16_t>());
  EXPECT_FALSE(neg42.Compatible<uint16_t>());
  EXPECT_TRUE(neg42.Compatible<int32_t>());
  EXPECT_FALSE(neg42.Compatible<uint32_t>());
  EXPECT_TRUE(neg42.Compatible<int64_t>());
  EXPECT_FALSE(neg42.Compatible<uint64_t>());
  EXPECT_TRUE(neg42.Compatible<absl::int128>());
  EXPECT_FALSE(neg42.Compatible<absl::uint128>());
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

TEST(BignumTest, AbslUint128Casting) {
  Bignum neg1(-1);
  EXPECT_EQ(neg1.Cast<absl::uint128>(), ~absl::uint128(0));

  // 2^128 -> low 128 bits == 0
  Bignum bn1 = *Bn("340282366920938463463374607431768211456");
  EXPECT_EQ(bn1.Cast<absl::uint128>(), absl::uint128(0));

  // 2^200 + 5 -> low 128 bits == 5
  Bignum bn2 =
      *Bn("1606938044258990275541962092341162602522202993782792835301381");
  EXPECT_EQ(bn2.Cast<absl::uint128>(), absl::uint128(5));
}

TEST(BignumTest, AbslInt128Casting) {
  const absl::int128 two127 = absl::int128(1) << 127;

  // +2^127 -> wraps to -2^127
  Bignum bn0 = *Bn("170141183460469231731687303715884105728");
  EXPECT_EQ(bn0.Cast<absl::int128>(), 0 - two127);

  // -(2^127) - 1 -> wraps to +2^127 - 1
  Bignum bn1 = *Bn("-170141183460469231731687303715884105729");
  EXPECT_EQ(bn1.Cast<absl::int128>(), two127 - 1);

  // 2^200 + 5 -> low 128 bits == 5
  Bignum bn2 =
      *Bn("1606938044258990275541962092341162602522202993782792835301381");
  EXPECT_EQ(bn2.Cast<absl::int128>(), absl::int128(5));

  // -(2^200 + 5) -> low 128 bits == -5
  Bignum bn3 =
      *Bn("-1606938044258990275541962092341162602522202993782792835301381");
  EXPECT_EQ(bn3.Cast<absl::int128>(), absl::int128(-5));
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
  const auto bn_u64max = Bignum(u64max);
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
  const auto bn_u64max = Bignum(u64max);
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

  const auto bn_u64max = Bignum(u64max);
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
  const auto bn_u32max = Bignum(u32max);
  EXPECT_EQ(bn_u32max * Bignum(2), *Bn("8589934590"));

  // 1x1 bigit fast path
  const auto bn_u64max = Bignum(u64max);
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
  EXPECT_EQ(Bignum(0).CountrZero(), 0);
  EXPECT_EQ(Bignum(1).CountrZero(), 0);
  EXPECT_EQ(Bignum(7).CountrZero(), 0);
  EXPECT_EQ(Bignum(-7).CountrZero(), 0);

  EXPECT_EQ(Bignum(2).CountrZero(), 1);
  EXPECT_EQ(Bignum(8).CountrZero(), 3);
  EXPECT_EQ(Bignum(10).CountrZero(), 1);  // 0b1010
  EXPECT_EQ(Bignum(12).CountrZero(), 2);  // 0b1100

  auto two_pow_64 = Bignum(1) << 64;
  EXPECT_EQ(two_pow_64.CountrZero(), 64);

  auto large_shifted = Bignum(6) << 100;  // 0b110 << 100
  EXPECT_EQ(large_shifted.CountrZero(), 101);

  auto neg_large_shifted = Bignum(-5) << 200;
  EXPECT_EQ(neg_large_shifted.CountrZero(), 200);
}

TEST(BignumTest, Bit) {
  EXPECT_FALSE(Bignum(0).Bit(0));
  EXPECT_FALSE(Bignum(0).Bit(100));

  // 5 = 0b101
  Bignum five(5);
  EXPECT_TRUE(five.Bit(0));
  EXPECT_FALSE(five.Bit(1));
  EXPECT_TRUE(five.Bit(2));
  EXPECT_FALSE(five.Bit(3));

  // Negative numbers should test the magnitude.
  Bignum neg_five(-5);
  EXPECT_TRUE(neg_five.Bit(0));
  EXPECT_FALSE(neg_five.Bit(1));
  EXPECT_TRUE(neg_five.Bit(2));

  // Test edges of and across bigits.
  Bignum high_bit_63 = Bignum(1) << 63;
  EXPECT_FALSE(high_bit_63.Bit(62));
  EXPECT_TRUE(high_bit_63.Bit(63));
  EXPECT_FALSE(high_bit_63.Bit(64));

  Bignum cross_bigit = (Bignum(1) << 100) + Bignum(1);
  EXPECT_TRUE(cross_bigit.Bit(0));
  EXPECT_TRUE(cross_bigit.Bit(100));
  EXPECT_FALSE(cross_bigit.Bit(50));
  EXPECT_FALSE(cross_bigit.Bit(1000));
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
  a.SetZero();
  EXPECT_TRUE(a.zero());

  Bignum b(-456);
  b.SetZero();
  EXPECT_EQ(b, Bignum(0));
}

TEST(BignumTest, SetNegativeSetPositive) {
  Bignum a(42);
  a.SetNegative();
  EXPECT_TRUE(a.negative());
  EXPECT_EQ(a, Bignum(-42));

  a.SetPositive();
  EXPECT_TRUE(a.positive());
  EXPECT_EQ(a, Bignum(42));
}

TEST(BignumTest, SetSign) {
  Bignum a(99);
  a.SetSign(-10);  // any negative
  EXPECT_EQ(a, Bignum(-99));

  a.SetSign(5);  // any positive
  EXPECT_EQ(a, Bignum(99));

  a.SetSign(0);
  EXPECT_TRUE(a.zero());
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

// TODO: Enable once benchmark is integrated.
#if 0

// RAII wrapper for OpenSSL BIGNUM
class OpenSSLBignum {
 public:
  OpenSSLBignum() : bn_(BN_new()) {}

  // Construct from a decimal number in a string.
  explicit OpenSSLBignum(const absl::string_view& decimal) : bn_(BN_new()) {
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

static std::vector<std::string> GenerateRandomNumbers(int bits) {
  std::vector<std::string> numbers;
  std::mt19937_64 rng(42);  // Fixed seed for reproducibility

  for (int i = 0; i < kRandomBignumCount; ++i) {
    std::string num;

    // Generate approximately `bits` worth of decimal digits
    int decimal_digits = (bits * 3) / 10;  // log10(2^bits) ≈ bits * 0.301

    // First digit can't be zero
    std::uniform_int_distribution<int> first_digit(1, 9);
    num += std::to_string(first_digit(rng));

    std::uniform_int_distribution<int> digit(0, 9);
    for (int j = 1; j < decimal_digits; ++j) {
      num += std::to_string(digit(rng));
    }

    numbers.push_back(num);
  }

  return numbers;
}

// Basic correctness test to ensure OpenSSL integration is working
TEST(BignumTestBenchmarkTest, OpenSSLIntegration) {
  OpenSSLBignum a(123);
  OpenSSLBignum b(456);
  OpenSSLBignum result;

  BN_add(result.get(), a.get(), b.get());

  char* str = BN_bn2dec(result.get());
  EXPECT_STREQ(str, "579");
  OPENSSL_free(str);
}

TEST(BignumTestBenchmarkTest, ResultsMatch) {
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

const std::vector<std::string>& SmallNumbers() {
  static absl::NoDestructor<std::vector<std::string>> numbers(  //
      GenerateRandomNumbers(64));
  return *numbers;
}

const std::vector<std::string>& MediumNumbers() {
  static absl::NoDestructor<std::vector<std::string>> numbers(  //
      GenerateRandomNumbers(256));
  return *numbers;
}

const std::vector<std::string>& LargeNumbers() {
  static absl::NoDestructor<std::vector<std::string>> numbers(  //
      GenerateRandomNumbers(1024));
  return *numbers;
}

const std::vector<std::string>& HugeNumbers() {
  static absl::NoDestructor<std::vector<std::string>> numbers(  //
      GenerateRandomNumbers(4096));
  return *numbers;
}

const std::vector<std::string>& MegaNumbers() {
  static absl::NoDestructor<std::vector<std::string>> numbers(  //
      GenerateRandomNumbers(18000));
  return *numbers;
}

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

void BM_Bignum_AddSmall(benchmark::State& state) {
  BignumBinaryOpBenchmark(state, SmallNumbers(), std::plus<Bignum>{});
}
BENCHMARK(BM_Bignum_AddSmall);

void BM_Bignum_AddMedium(benchmark::State& state) {
  BignumBinaryOpBenchmark(state, MediumNumbers(), std::plus<Bignum>{});
}
BENCHMARK(BM_Bignum_AddMedium);

void BM_Bignum_AddLarge(benchmark::State& state) {
  BignumBinaryOpBenchmark(state, LargeNumbers(), std::plus<Bignum>{});
}
BENCHMARK(BM_Bignum_AddLarge);

void BM_Bignum_AddHuge(benchmark::State& state) {
  BignumBinaryOpBenchmark(state, HugeNumbers(), std::plus<Bignum>{});
}
BENCHMARK(BM_Bignum_AddHuge);

void BM_Bignum_AddMega(benchmark::State& state) {
  BignumBinaryOpBenchmark(state, MegaNumbers(), std::plus<Bignum>{});
}
BENCHMARK(BM_Bignum_AddMega);

void BM_OpenSSL_AddSmall(benchmark::State& state) {
  OpenSSLBinaryOpBenchmark(state, SmallNumbers(), BN_add);
}
BENCHMARK(BM_OpenSSL_AddSmall);

void BM_OpenSSL_AddMedium(benchmark::State& state) {
  OpenSSLBinaryOpBenchmark(state, MediumNumbers(), BN_add);
}
BENCHMARK(BM_OpenSSL_AddMedium);

void BM_OpenSSL_AddLarge(benchmark::State& state) {
  OpenSSLBinaryOpBenchmark(state, LargeNumbers(), BN_add);
}
BENCHMARK(BM_OpenSSL_AddLarge);

void BM_OpenSSL_AddHuge(benchmark::State& state) {
  OpenSSLBinaryOpBenchmark(state, HugeNumbers(), BN_add);
}
BENCHMARK(BM_OpenSSL_AddHuge);

void BM_OpenSSL_AddMega(benchmark::State& state) {
  OpenSSLBinaryOpBenchmark(state, MegaNumbers(), BN_add);
}
BENCHMARK(BM_OpenSSL_AddMega);

void BM_Bignum_MulSmall(benchmark::State& state) {
  BignumBinaryOpBenchmark(state, SmallNumbers(), std::multiplies<Bignum>{});
}
BENCHMARK(BM_Bignum_MulSmall);

void BM_Bignum_MulMedium(benchmark::State& state) {
  BignumBinaryOpBenchmark(state, MediumNumbers(), std::multiplies<Bignum>{});
}
BENCHMARK(BM_Bignum_MulMedium);

void BM_Bignum_MulLarge(benchmark::State& state) {
  BignumBinaryOpBenchmark(state, LargeNumbers(), std::multiplies<Bignum>{});
}
BENCHMARK(BM_Bignum_MulLarge);

void BM_Bignum_MulHuge(benchmark::State& state) {
  BignumBinaryOpBenchmark(state, HugeNumbers(), std::multiplies<Bignum>{});
}
BENCHMARK(BM_Bignum_MulHuge);

void BM_Bignum_MulMega(benchmark::State& state) {
  BignumBinaryOpBenchmark(state, MegaNumbers(), std::multiplies<Bignum>{});
}
BENCHMARK(BM_Bignum_MulMega);

void BM_OpenSSL_MulSmall(benchmark::State& state) {
  OpenSSLMulOpBenchmark(state, SmallNumbers(), BN_mul);
}
BENCHMARK(BM_OpenSSL_MulSmall);

void BM_OpenSSL_MulMedium(benchmark::State& state) {
  OpenSSLMulOpBenchmark(state, MediumNumbers(), BN_mul);
}
BENCHMARK(BM_OpenSSL_MulMedium);

void BM_OpenSSL_MulLarge(benchmark::State& state) {
  OpenSSLMulOpBenchmark(state, LargeNumbers(), BN_mul);
}
BENCHMARK(BM_OpenSSL_MulLarge);

void BM_OpenSSL_MulHuge(benchmark::State& state) {
  OpenSSLMulOpBenchmark(state, HugeNumbers(), BN_mul);
}
BENCHMARK(BM_OpenSSL_MulHuge);

void BM_OpenSSL_MulMega(benchmark::State& state) {
  OpenSSLMulOpBenchmark(state, MegaNumbers(), BN_mul);
}
BENCHMARK(BM_OpenSSL_MulMega);

void BM_Bignum_PowSmall(benchmark::State& state) {
  BignumPowBenchmark(state, SmallNumbers(), 20);
}
BENCHMARK(BM_Bignum_PowSmall);

void BM_Bignum_PowMedium(benchmark::State& state) {
  BignumPowBenchmark(state, MediumNumbers(), 10);
}
BENCHMARK(BM_Bignum_PowMedium);

void BM_OpenSSL_PowSmall(benchmark::State& state) {
  OpenSSLPowBenchmark(state, SmallNumbers(), 20);
}
BENCHMARK(BM_OpenSSL_PowSmall);

void BM_OpenSSL_PowMedium(benchmark::State& state) {
  OpenSSLPowBenchmark(state, MediumNumbers(), 10);
}
BENCHMARK(BM_OpenSSL_PowMedium);
#endif
