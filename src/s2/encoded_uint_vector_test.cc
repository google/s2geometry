// Copyright 2018 Google Inc. All Rights Reserved.
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

// Author: ericv@google.com (Eric Veach)

#include "s2/encoded_uint_vector.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <benchmark/benchmark.h>
#include <gtest/gtest.h>
#include "absl/log/absl_check.h"
#include "s2/util/coding/coder.h"

using std::vector;

namespace s2coding {

// Make sure that this class is compact since it is extensively used.
// 16 for 64-bit, 12 for 32-bit.
static_assert(sizeof(EncodedUintVector<uint64_t>) <= 16, "too big");

template <class T>
void TestEncodedUintVector(const vector<T>& expected, size_t expected_bytes) {
  Encoder encoder;
  EncodeUintVector<T>(expected, &encoder);
  EXPECT_EQ(expected_bytes, encoder.length());
  Decoder decoder(encoder.base(), encoder.length());
  EncodedUintVector<T> actual;
  ASSERT_TRUE(actual.Init(&decoder));
  EXPECT_EQ(actual.Decode(), expected);
}

TEST(EncodedUintVectorTest, Empty) {
  TestEncodedUintVector(vector<uint32_t>{}, 1);
}

TEST(EncodedUintVectorTest, Zero) {
  TestEncodedUintVector(vector<uint64_t>{0}, 2);
}

TEST(EncodedUintVectorTest, RepeatedZeros) {
  TestEncodedUintVector(vector<uint16_t>{0, 0, 0}, 4);
}

TEST(EncodedUintVectorTest, MaxInt) {
  TestEncodedUintVector(vector<uint64_t>{~0ULL}, 9);
}

TEST(EncodedUintVectorTest, OneByte) {
  TestEncodedUintVector(vector<uint64_t>{0, 255, 1, 254}, 5);
}

TEST(EncodedUintVectorTest, TwoBytes) {
  TestEncodedUintVector(vector<uint64_t>{0, 255, 256, 254}, 9);
}

TEST(EncodedUintVectorTest, ThreeBytes) {
  TestEncodedUintVector(vector<uint64_t>{0xffffff, 0x0102, 0, 0x050403}, 13);
}

TEST(EncodedUintVectorTest, EightBytes) {
  TestEncodedUintVector(vector<uint64_t>{~0ULL, 0, 0x0102030405060708}, 25);
}

template <class T>
vector<T> MakeSortedTestVector(int bytes_per_value, int num_values) {
  ABSL_DCHECK_LE(bytes_per_value, sizeof(T));
  T limit_value = ~T{0} >> (8 * (sizeof(T) - bytes_per_value));
  vector<T> values;
  for (int i = 0; i + 1 < num_values; ++i) {
    values.push_back(limit_value * (static_cast<double>(i) / (num_values - 1)));
  }
  // The last value needs special handling since casting it to "double" loses
  // precision when T == uint64_t.
  values.push_back(limit_value);
  ABSL_CHECK(std::is_sorted(values.begin(), values.end()));
  return values;
}

template <class T>
EncodedUintVector<T> MakeEncodedVector(const vector<T>& values,
                                       Encoder* encoder) {
  EncodeUintVector<T>(values, encoder);
  Decoder decoder(encoder->base(), encoder->length());
  EncodedUintVector<T> actual;
  ABSL_CHECK(actual.Init(&decoder));
  return actual;
}

template <class T>
void TestLowerBound(int bytes_per_value, int num_values) {
  auto v = MakeSortedTestVector<T>(bytes_per_value, num_values);
  Encoder encoder;
  auto actual = MakeEncodedVector(v, &encoder);
  for (T x : v) {
    EXPECT_EQ(std::lower_bound(v.begin(), v.end(), x) - v.begin(),
              actual.lower_bound(x));
    if (x > 0) {
      EXPECT_EQ(std::lower_bound(v.begin(), v.end(), x - 1) - v.begin(),
                actual.lower_bound(x - 1));
    }
  }
}

TEST(EncodedUintVector, LowerBound) {
  for (int bytes_per_value = 8; bytes_per_value <= 8; ++bytes_per_value) {
    TestLowerBound<uint64_t>(bytes_per_value, 10);
    if (bytes_per_value <= 4) {
      TestLowerBound<uint32_t>(bytes_per_value, 500);
      if (bytes_per_value <= 2) {
        TestLowerBound<uint16_t>(bytes_per_value, 100);
      }
    }
  }
}

TEST(EncodedUintVectorTest, RoundtripEncoding) {
  vector<uint64_t> values{10, 20, 30, 40};

  Encoder a_encoder;
  auto a = MakeEncodedVector<uint64_t>(values, &a_encoder);
  ASSERT_EQ(a.Decode(), values);

  Encoder b_encoder;
  a.Encode(&b_encoder);
  Decoder decoder(b_encoder.base(), b_encoder.length());

  EncodedUintVector<uint64_t> v2;
  ASSERT_TRUE(v2.Init(&decoder));

  EXPECT_EQ(v2.Decode(), values);
}

template <class T>
static void BM_DecodeValue(benchmark::State& state) {
  const int bytes_per_value = state.range(0);
  constexpr int kNumValues = 1 << 10;  // Must be power of 2.
  ABSL_DCHECK_EQ(kNumValues & (kNumValues - 1), 0);

  vector<T> values(kNumValues);
  for (T& x : values) {
    for (int i = 0; i < bytes_per_value; ++i) {
      x = (x << 8) | 0xa5;  // Immaterial.
    }
  }
  Encoder encoder;
  EncodeUintVector<T>(values, &encoder);
  Decoder decoder(encoder.base(), encoder.length());
  EncodedUintVector<T> actual;
  ABSL_CHECK(actual.Init(&decoder));
  int i = 0;
  for (auto _ : state) {
    benchmark::DoNotOptimize(actual[i]);
    i = (i + 1) & (kNumValues - 1);
  }
}
BENCHMARK_TEMPLATE(BM_DecodeValue, uint32_t)
    ->Arg(0)
    ->Arg(1)
    ->Arg(2)
    ->Arg(3)
    ->Arg(4);
BENCHMARK_TEMPLATE(BM_DecodeValue, uint64_t)
    ->Arg(0)
    ->Arg(1)
    ->Arg(2)
    ->Arg(3)
    ->Arg(4)
    ->Arg(5)
    ->Arg(6)
    ->Arg(7)
    ->Arg(8);

template <class T>
static void BM_LowerBound(benchmark::State& state) {
  const int num_values = state.range(0);
  const int bytes_per_value = state.range(1);
  auto values = MakeSortedTestVector<T>(bytes_per_value, num_values);
  T limit_value = values.back();
  Encoder encoder;
  auto actual = MakeEncodedVector(values, &encoder);
  T value = 0, delta = 3 * (limit_value / num_values);
  for (auto _ : state) {
    benchmark::DoNotOptimize(actual.lower_bound(value));
    if ((value += delta) > limit_value) value -= limit_value;
  }
}
// Use case: S2LaxPolygonShape::cumulative_vertices.
// Stores one cumulative vertex count per loop.
BENCHMARK_TEMPLATE(BM_LowerBound, uint32_t)
    ->ArgPair(10, 2)
    ->ArgPair(1000, 2)
    ->ArgPair(1000000, 3)
    ->ArgPair(1000000, 4);

// Use case: EncodedS2ShapeIndex::cell_ids_.
// Stores one S2CellId value per cell.
BENCHMARK_TEMPLATE(BM_LowerBound, uint64_t)
    ->ArgPair(10, 1)
    ->ArgPair(1000, 2)
    ->ArgPair(100000, 3)
    ->ArgPair(1000000, 4)
    ->ArgPair(1000000, 5)
    ->ArgPair(1000000, 6)
    ->ArgPair(1000000, 7)
    ->ArgPair(1000000, 8);

}  // namespace s2coding
