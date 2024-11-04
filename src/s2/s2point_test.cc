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

#include <algorithm>
#include <cstddef>
#include <string>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/algorithm/container.h"
#include "absl/flags/flag.h"
#include "absl/hash/hash.h"
#include "absl/log/log_streamer.h"
#include "absl/random/random.h"
#include "absl/strings/str_format.h"
#include "s2/s2coder_testing.h"
#include "s2/s2error.h"
#include "s2/s2random.h"
#include "s2/s2testing.h"

namespace {

using ::testing::Eq;

// A matcher that applies DoubleNear to each field of an S2Point.
MATCHER_P2(NearPoint, P, tol,
           absl::StrFormat("Near point (%f, %f, %f) (tol: %.3e)", P.x(), P.y(),
                           P.z(), tol)) {
  return ::testing::Value(arg.x(), ::testing::DoubleNear(P.x(), tol)) &&
         ::testing::Value(arg.y(), ::testing::DoubleNear(P.y(), tol)) &&
         ::testing::Value(arg.z(), ::testing::DoubleNear(P.z(), tol));
}

TEST(S2Point, HashSpreads) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "HASH_SPREADS",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  int kTestPoints = 1 << 16;
  std::vector<size_t> hashes;
  hashes.reserve(kTestPoints);
  std::vector<S2Point> points;
  points.reserve(kTestPoints);
  S2PointHash hasher;
  S2Point base = S2Point(1, 1, 1);
  for (int i = 0; i < kTestPoints; ++i) {
    // All points in a tiny cap to test avalanche property of hash
    // function (the cap would be of radius 1mm on Earth (4*10^9/2^35).
    S2Point perturbed = base + s2random::Point(bitgen) / (1ULL << 35);
    perturbed = perturbed.Normalize();
    hashes.push_back(hasher(perturbed));
    points.push_back(perturbed);
  }

  absl::c_sort(hashes);
  hashes.erase(std::unique(hashes.begin(), hashes.end()), hashes.end());

  absl::c_sort(points);
  points.erase(std::unique(points.begin(), points.end()), points.end());

  // Point collisions are somewhat unlikely.  This test is ~5% flaky with a
  // value of 0 here and ~0.4% flaky with a value of 1.  Using 2 is <0.05%
  // flaky.  A collision analysis would be nice to have, rather than using
  // empirically determined thresholds.
  const int num_point_collisions = kTestPoints - points.size();
  const int num_hash_collisions = kTestPoints - hashes.size();
  EXPECT_LE(num_point_collisions, 2);
  // Hashes should add no additional collisions.
  EXPECT_EQ(num_point_collisions, num_hash_collisions);
}

TEST(S2Point, IsAVector) {
  // Check that we haven't accidentally increased the size or modified
  // the alignment of S2Point.
  EXPECT_EQ(sizeof(S2Point), sizeof(Vector3_d));
  EXPECT_EQ(alignof(S2Point), alignof(Vector3_d));
}

TEST(S2Point, CoderWorks) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "CODER_WORKS",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  S2Point point = s2random::Point(bitgen);
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
