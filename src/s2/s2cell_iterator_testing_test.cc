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

#include "s2/s2cell_iterator_testing.h"

#include <utility>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/container/btree_map.h"
#include "s2/s2cell_id.h"

namespace {

using ::testing::Eq;
using ::testing::Ge;
using ::testing::IsFalse;
using ::testing::IsTrue;

S2CellId MakeId(absl::string_view token) { return S2CellId::FromToken(token); }

TEST(MockIterator, BasicFunctionality) {
  absl::btree_map<S2CellId, int> map = {
      std::make_pair(MakeId("89c259"), 0),
      std::make_pair(MakeId("89c2598c"), 1),
      std::make_pair(MakeId("89c2599"), 2),
      std::make_pair(MakeId("89c259c"), 3),
      std::make_pair(MakeId("89c259d2"), 4),
      std::make_pair(MakeId("89c25c"), 5),
  };

  auto iter = MakeMockS2CellIterator(&map);

  // Iterator should be positioned at the beginning after construction.
  EXPECT_THAT(iter.id(), Eq(MakeId("89c259")));
  EXPECT_THAT(iter.value(), Eq(0));

  // Finish should make the iterator done.
  iter.Finish();
  EXPECT_THAT(iter.done(), IsTrue());

  // Move to beginning, iterator should no longer test as done.
  iter.Begin();
  EXPECT_THAT(iter.done(), IsFalse());

  // Move to next element, it should be the second one.
  iter.Next();
  EXPECT_THAT(iter.id(), MakeId("89c2598c"));
  EXPECT_THAT(iter.value(), Eq(1));

  // Move past the second element, is it greater?
  iter.Next();
  EXPECT_THAT(iter.id(), Ge(MakeId("89c2598c")));
  EXPECT_THAT(iter.value(), Eq(2));

  // Move back, should be on second element again.
  iter.Prev();
  EXPECT_THAT(iter.id(), Eq(MakeId("89c2598c")));
  EXPECT_THAT(iter.value(), Eq(1));

  // Seek to parent of nodes, should be positioned on first child.
  iter.Seek(MakeId("89c259c"));
  EXPECT_THAT(iter.id(), Eq(MakeId("89c259c")));
  iter.Next();
  EXPECT_THAT(iter.id(), Eq(MakeId("89c259d2")));
  iter.Next();
  iter.Next();
  EXPECT_THAT(iter.done(), IsTrue());
}

}  // namespace
