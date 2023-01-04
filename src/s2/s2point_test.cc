// Copyright 2005 Google Inc. All Rights Reserved.
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

#include "s2/s2point.h"

#include <cstddef>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/container/flat_hash_set.h"
#include "absl/flags/flag.h"
#include "absl/hash/hash.h"
#include "s2/s2coder_testing.h"
#include "s2/s2error.h"
#include "s2/s2testing.h"

namespace {

using ::testing::Eq;

TEST(S2Point, HashSpreads) {
  int kTestPoints = 1 << 16;
  absl::flat_hash_set<size_t> set;
  absl::flat_hash_set<S2Point, S2PointHash> points;
  S2PointHash hasher;
  S2Point base = S2Point(1, 1, 1);
  for (int i = 0; i < kTestPoints; ++i) {
    // All points in a tiny cap to test avalanche property of hash
    // function (the cap would be of radius 1mm on Earth (4*10^9/2^35).
    S2Point perturbed = base + S2Testing::RandomPoint() / (1ULL << 35);
    perturbed = perturbed.Normalize();
    set.insert(hasher(perturbed));
    points.insert(perturbed);
  }
  // A real collision is extremely unlikely.
  EXPECT_EQ(0, kTestPoints - points.size());
  // Allow a few for the hash.
  EXPECT_GE(10, kTestPoints - set.size());
}

TEST(S2Point, IsAVector) {
  // Check that we haven't accidentally increased the size or modified
  // the alignment of S2Point.
  EXPECT_EQ(sizeof(S2Point), sizeof(Vector3_d));
  EXPECT_EQ(alignof(S2Point), alignof(Vector3_d));
}

TEST(S2Point, CoderWorks) {
  S2Testing::rnd.Reset(absl::GetFlag(FLAGS_s2_random_seed));

  S2Point point = S2Testing::RandomPoint();
  S2Error error;
  S2Point decoded = s2coding::RoundTrip(S2Point::Coder(), point, error);
  EXPECT_TRUE(error.ok());
  EXPECT_EQ(decoded, point);
}

TEST(S2Point, SubtractionWorks) {
  S2Point a(1, 2, 3);
  S2Point b(1, 1, 1);
  a -= b;
  EXPECT_THAT(a, Eq(S2Point(0, 1, 2)));
}

TEST(S2Point, ElementWiseDivisionWorks) {
  S2Point a(4, 8, 16);
  S2Point b(2, 2, 2);
  EXPECT_THAT(a.DivComponents(b), Eq(S2Point(2, 4, 8)));
}

TEST(S2Point, SqrtWorks) {
  S2Point a(4, 9, 16);
  EXPECT_THAT(a.Sqrt(), Eq(S2Point(2, 3, 4)));
}

TEST(S2Point, FloorWorks) {
  S2Point a(1.4, 1.5, 1.6);
  EXPECT_THAT(a.Floor(), Eq(S2Point(1, 1, 1)));
}

TEST(S2Point, CeilWorks) {
  S2Point a(1.4, 1.5, 1.6);
  EXPECT_THAT(a.Ceil(), Eq(S2Point(2, 2, 2)));
}

TEST(S2Point, FRoundWorks) {
  S2Point a(1.4, 1.5, 1.6);
  EXPECT_THAT(a.FRound(), Eq(S2Point(1, 2, 2)));
}

}  // namespace
