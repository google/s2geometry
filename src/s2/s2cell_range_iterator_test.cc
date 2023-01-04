// Copyright 2013 Google Inc. All Rights Reserved.
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

#include "s2/s2cell_range_iterator.h"

#include <memory>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "s2/mutable_s2shape_index.h"
#include "s2/s2cell_id.h"
#include "s2/s2text_format.h"

namespace {

using ::testing::Eq;
using ::testing::IsFalse;
using ::testing::IsTrue;

TEST(S2CellRangeIterator, Relation) {
  // Create an index with one point each on S2CellId faces 0, 1, and 2.
  auto index = s2textformat::MakeIndexOrDie("0:0 | 0:90 | 90:0 # #");
  auto it0 = MakeS2CellRangeIterator(index.get());
  auto it1 = MakeS2CellRangeIterator(index.get());
  it1.Next();
  EXPECT_THAT(it0.Relation(it1), Eq(-1));
  EXPECT_THAT(it1.Relation(it0), Eq(+1));
  it1.Prev();
  EXPECT_THAT(it1.Relation(it0), Eq(0));
  it1.Finish();
  EXPECT_THAT(it1.Relation(it0), Eq(+1));
}

TEST(S2CellRangeIterator, Next) {
  // Create an index with one point each on S2CellId faces 0, 1, and 2.
  auto index = s2textformat::MakeIndexOrDie("0:0 | 0:90 | 90:0 # #");
  auto it = MakeS2CellRangeIterator(index.get());
  EXPECT_THAT(0, Eq(it.id().face()));
  it.Next();
  EXPECT_THAT(1, Eq(it.id().face()));
  it.Next();
  EXPECT_THAT(2, Eq(it.id().face()));
  it.Next();
  EXPECT_THAT(S2CellId::Sentinel(), Eq(it.id()));
  EXPECT_THAT(it.done(), IsTrue());
}

TEST(S2CellRangeIterator, Locate) {
  // Create two indices with one point each on S2CellId faces 0, 1, and 2.
  auto index0 = s2textformat::MakeIndexOrDie("0:0 | 0:90 | 90:0 # #");
  auto index1 = s2textformat::MakeIndexOrDie("0:0 | 0:90 | 90:0 # #");
  auto it0 = MakeS2CellRangeIterator(index0.get());
  auto it1 = MakeS2CellRangeIterator(index1.get());
  it0.Next();
  it1.Locate(it0);
  EXPECT_THAT(it1.id(), Eq(it0.id()));
}

TEST(S2CellRangeIterator, EmptyIndex) {
  auto empty = s2textformat::MakeIndexOrDie("# #");
  auto non_empty = s2textformat::MakeIndexOrDie("0:0 # #");
  auto empty_it = MakeS2CellRangeIterator(empty.get());
  auto non_empty_it = MakeS2CellRangeIterator(non_empty.get());
  EXPECT_THAT(non_empty_it.done(), IsFalse());
  EXPECT_THAT(empty_it.done(), IsTrue());

  empty_it.SeekTo(non_empty_it);
  EXPECT_THAT(empty_it.done(), IsTrue());

  empty_it.SeekBeyond(non_empty_it);
  EXPECT_THAT(empty_it.done(), IsTrue());

  empty_it.SeekTo(empty_it);
  EXPECT_THAT(empty_it.done(), IsTrue());

  empty_it.SeekBeyond(empty_it);
  EXPECT_THAT(empty_it.done(), IsTrue());
}

}  // namespace
