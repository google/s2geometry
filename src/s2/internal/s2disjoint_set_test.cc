// Copyright 2022 Google Inc. All Rights Reserved.
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


#include "s2/internal/s2disjoint_set.h"

#include <optional>

#include <gtest/gtest.h>
#include "absl/log/absl_check.h"
#include "s2/s2point.h"

namespace internal {
namespace {

TEST(DisjointSetTest, S2PointSetCompiles) {
  DisjointSet<S2Point> set;
  (void)set;
}

TEST(DisjointSetTest, InsertMoreThanOnceFails) {
  DisjointSet<int> set;
  EXPECT_TRUE(set.Add(1));
  EXPECT_FALSE(set.Add(1));
  EXPECT_FALSE(set.Add(1));
}

TEST(DisjointSetTest, FindRootWorks) {
  DisjointSet<int> set;
  set.Add(1);
  std::optional<int> root = set.FindRoot(1);
  ABSL_CHECK(root.has_value());
  EXPECT_EQ(root.value(), 1);
  EXPECT_FALSE(set.FindRoot(2).has_value());
}

TEST(DisjointSetTest, UnionWorks) {
  DisjointSet<int> set;
  for (int i = 0; i < 10; ++i) {
    ABSL_CHECK(set.Add(i));
  }

  // Union into two disjoint sets
  for (int i = 0; i < 4; ++i) {
    ABSL_CHECK(set.Union(i + 0, i + 1));
    ABSL_CHECK(set.Union(i + 5, i + 6));
  }

  // The values in the two sets should all have the same root now.
  for (int i = 0; i < 5; ++i) {
    std::optional<int> root = set.FindRoot(i);
    ABSL_CHECK(root.has_value());
    EXPECT_EQ(root.value(), 0);

    root = set.FindRoot(i + 5);
    ABSL_CHECK(root.has_value());
    EXPECT_EQ(root.value(), 5);
  }

  // Check that attempting to union with something not in the set fails.
  EXPECT_FALSE(set.Union(0, 13));
  EXPECT_FALSE(set.Union(13, 0));
  EXPECT_FALSE(set.Union(12, 13));

  // Union of two elements from two sets should unify the sets
  ABSL_CHECK(set.Union(3, 7));
  for (int i = 0; i < 10; ++i) {
    std::optional<int> root = set.FindRoot(i);
    ABSL_CHECK(root.has_value());
    EXPECT_EQ(root.value(), 0);
  }
}

TEST(DisjointSetTest, SizeAndClearWorks) {
  DisjointSet<int> set;
  set.Reserve(10);
  for (int i = 0; i < 10; ++i) {
    ABSL_CHECK(set.Add(i));
  }

  EXPECT_EQ(set.Size(), 10);
  for (int i = 0; i < set.Size() - 1; ++i) {
    ABSL_CHECK(set.Union(i, i + 1));
  }

  // Unioning doesn't change size.
  EXPECT_EQ(set.Size(), 10);
  set.Clear();
  EXPECT_EQ(set.Size(), 0);
}

}  // namespace
}  // namespace internal
