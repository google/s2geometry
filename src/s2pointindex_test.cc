// Copyright 2015 Google Inc. All Rights Reserved.
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

#include "s2pointindex.h"

#include <set>

#include <gtest/gtest.h>
#include "s2cellid.h"
#include "s2cellunion.h"
#include "s2testing.h"

class S2PointIndexTest : public ::testing::Test {
 private:
  typedef S2PointIndex<int> Index;
  typedef Index::PointData PointData;
  typedef std::multiset<PointData> Contents;
  Index index_;
  Contents contents_;

 public:
  void Add(S2Point const& point, int data) {
    index_.Add(point, data);
    contents_.insert(PointData(point, data));
  }

  void Remove(S2Point const& point, int data) {
    index_.Remove(point, data);
    // If there are multiple copies, remove only one.
    contents_.erase(contents_.find(PointData(point, data)));
  }

  void Verify() {
    Contents remaining = contents_;
    for (Index::Iterator it(index_); !it.Done(); it.Next()) {
      Contents::iterator element = remaining.find(it.point_data());
      EXPECT_TRUE(element != remaining.end());
      remaining.erase(element);
    }
    EXPECT_TRUE(remaining.empty());
  }

  void TestIteratorMethods() {
    Index::Iterator it(index_);
    EXPECT_TRUE(it.AtBegin());
    it.Finish();
    EXPECT_TRUE(it.Done());

    // Iterate through all the cells in the index.
    S2CellId prev_cellid = S2CellId::None();
    S2CellId min_cellid = S2CellId::Begin(S2CellId::kMaxLevel);
    for (it.Reset(); !it.Done(); it.Next()) {
      S2CellId cellid = it.id();
      EXPECT_EQ(cellid, S2CellId::FromPoint(it.point()));

      typename Index::Iterator it2(index_);
      if (cellid == prev_cellid) {
        it2.Seek(cellid);
      }

      // Generate a cellunion that covers the range of empty leaf cells between
      // the last cell and this one.  Then make sure that seeking to any of
      // those cells takes us to the immediately following cell.
      S2CellUnion skipped;
      skipped.InitFromBeginEnd(min_cellid, cellid.range_min());
      for (int i = 0; i < skipped.num_cells(); ++i) {
        it2.Seek(skipped.cell_id(i));
        EXPECT_EQ(cellid, it2.id());
      }
      // Test Prev(), Next(), Seek(), and SeekForward().
      if (prev_cellid.is_valid()) {
        EXPECT_FALSE(it.AtBegin());
        it2 = it;
        it2.Prev();
        EXPECT_EQ(prev_cellid, it2.id());
        it2.Next();
        EXPECT_EQ(cellid, it2.id());
        it2.Seek(prev_cellid);
        EXPECT_EQ(prev_cellid, it2.id());
        it2.SeekForward(cellid);
        EXPECT_EQ(cellid, it2.id());
        it2.SeekForward(prev_cellid);
        EXPECT_EQ(cellid, it2.id());
      }
      prev_cellid = cellid;
      min_cellid = cellid.range_max().next();
    }
  }
};

TEST_F(S2PointIndexTest, NoPoints) {
  TestIteratorMethods();
}

TEST_F(S2PointIndexTest, RandomPoints) {
  for (int i = 0; i < 1000; ++i) {
    Add(S2Testing::RandomPoint(), S2Testing::rnd.Uniform(100));
  }
  Verify();
  TestIteratorMethods();
}
