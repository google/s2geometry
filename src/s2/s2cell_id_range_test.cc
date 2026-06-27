// Copyright 2026 Google Inc. All Rights Reserved.
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

#include "s2/s2cell_id_range.h"

#include <cstdint>
#include <vector>

#include <gtest/gtest.h>

#include "s2/s2cell_id.h"

namespace {

// Convenience: range of all cells at a given level on a single face.
S2CellIdRange FaceRange(int face, int level) {
  S2CellId f = S2CellId::FromFace(face);
  return S2CellIdRange{f.child_begin(level), f.child_end(level)};
}

// Collect all cells produced by iterating a range forward.
std::vector<S2CellId> Collect(S2CellIdRange r) {
  std::vector<S2CellId> result;
  S2CellIdForwardIterator it{r.begin, r.end};
  while (auto cell = it.next()) {
    result.push_back(*cell);
  }
  return result;
}

// Collect all cells produced by iterating a range in reverse.
std::vector<S2CellId> CollectReverse(S2CellIdRange r) {
  std::vector<S2CellId> result;
  S2CellIdReverseIterator it{r.end, r.begin};
  while (auto cell = it.next()) {
    result.push_back(*cell);
  }
  return result;
}

// ---------------------------------------------------------------------------
// S2CellIdRange::size
// ---------------------------------------------------------------------------

TEST(S2CellIdRangeTest, SizeEmptyRange) {
  S2CellId face = S2CellId::FromFace(0);
  S2CellIdRange r{face.child_begin(1), face.child_begin(1)};
  EXPECT_EQ(r.size(), 0);
}

TEST(S2CellIdRangeTest, SizeFaceLevel1) {
  // Each face has 4 level-1 children.
  S2CellIdRange r = FaceRange(0, 1);
  EXPECT_EQ(r.size(), 4);
}

TEST(S2CellIdRangeTest, SizeFaceLevel2) {
  // Each face has 16 level-2 children.
  S2CellIdRange r = FaceRange(0, 2);
  EXPECT_EQ(r.size(), 16);
}

TEST(S2CellIdRangeTest, SizeAllFacesLevel0) {
  S2CellIdRange r{S2CellId::Begin(0), S2CellId::End(0)};
  EXPECT_EQ(r.size(), 6);
}

// ---------------------------------------------------------------------------
// S2CellIdRange::at
// ---------------------------------------------------------------------------

TEST(S2CellIdRangeTest, AtFirstCell) {
  S2CellIdRange r = FaceRange(0, 1);
  EXPECT_EQ(r.at(0), r.begin);
}

TEST(S2CellIdRangeTest, AtLastCell) {
  S2CellIdRange r = FaceRange(0, 1);
  EXPECT_EQ(r.at(3), r.begin.advance(3));
}

TEST(S2CellIdRangeTest, AtMiddleCell) {
  S2CellIdRange r = FaceRange(0, 2);
  EXPECT_EQ(r.at(7), r.begin.advance(7));
}

// ---------------------------------------------------------------------------
// S2CellIdRange::slice
// ---------------------------------------------------------------------------

TEST(S2CellIdRangeTest, SliceFullRange) {
  S2CellIdRange r = FaceRange(0, 1);
  S2CellIdRange s = r.slice(0, 4);
  EXPECT_EQ(s.begin, r.begin);
  EXPECT_EQ(s.end, r.end);
}

TEST(S2CellIdRangeTest, SliceFirstTwo) {
  S2CellIdRange r = FaceRange(0, 1);
  S2CellIdRange s = r.slice(0, 2);
  EXPECT_EQ(s.size(), 2);
  EXPECT_EQ(s.begin, r.begin);
  EXPECT_EQ(s.end, r.begin.advance(2));
}

TEST(S2CellIdRangeTest, SliceLastTwo) {
  S2CellIdRange r = FaceRange(0, 1);
  S2CellIdRange s = r.slice(2, 4);
  EXPECT_EQ(s.size(), 2);
  EXPECT_EQ(s.begin, r.begin.advance(2));
  EXPECT_EQ(s.end, r.end);
}

TEST(S2CellIdRangeTest, SliceEmpty) {
  S2CellIdRange r = FaceRange(0, 1);
  S2CellIdRange s = r.slice(2, 2);
  EXPECT_EQ(s.size(), 0);
  EXPECT_EQ(s.begin, s.end);
}

TEST(S2CellIdRangeTest, SliceSingleCell) {
  S2CellIdRange r = FaceRange(0, 1);
  S2CellIdRange s = r.slice(1, 2);
  EXPECT_EQ(s.size(), 1);
  EXPECT_EQ(s.at(0), r.begin.advance(1));
}

// ---------------------------------------------------------------------------
// S2CellIdRange::contains
// ---------------------------------------------------------------------------

TEST(S2CellIdRangeTest, ContainsCellInRange) {
  S2CellIdRange r = FaceRange(0, 1);
  EXPECT_TRUE(r.contains(r.begin));
  EXPECT_TRUE(r.contains(r.begin.advance(1)));
  EXPECT_TRUE(r.contains(r.begin.advance(3)));
}

TEST(S2CellIdRangeTest, ContainsExcludesEnd) {
  S2CellIdRange r = FaceRange(0, 1);
  EXPECT_FALSE(r.contains(r.end));
}

TEST(S2CellIdRangeTest, ContainsCellFromDifferentFace) {
  S2CellIdRange r = FaceRange(0, 1);
  S2CellId other = S2CellId::FromFace(1).child_begin(1);
  EXPECT_FALSE(r.contains(other));
}

TEST(S2CellIdRangeTest, ContainsCellAtWrongLevel) {
  S2CellIdRange r = FaceRange(0, 1);
  // A level-2 cell that falls within the level-1 range bounds should not match.
  S2CellId level2 = r.begin.child_begin();
  EXPECT_FALSE(r.contains(level2));
}

TEST(S2CellIdRangeTest, ContainsEmptyRange) {
  S2CellId face = S2CellId::FromFace(0);
  S2CellIdRange r{face.child_begin(1), face.child_begin(1)};
  EXPECT_FALSE(r.contains(face.child_begin(1)));
}

// ---------------------------------------------------------------------------
// S2CellIdForwardIterator
// ---------------------------------------------------------------------------

TEST(S2CellIdForwardIteratorTest, IteratesAllCells) {
  S2CellIdRange r = FaceRange(0, 1);
  std::vector<S2CellId> cells = Collect(r);
  ASSERT_EQ(cells.size(), 4u);
  for (int i = 0; i < 4; ++i) {
    EXPECT_EQ(cells[i], r.begin.advance(i));
  }
}

TEST(S2CellIdForwardIteratorTest, EmptyRangeProducesNoItems) {
  S2CellId face = S2CellId::FromFace(0);
  S2CellIdRange r{face.child_begin(1), face.child_begin(1)};
  EXPECT_TRUE(Collect(r).empty());
}

TEST(S2CellIdForwardIteratorTest, ReturnsNulloptWhenExhausted) {
  S2CellIdRange r = FaceRange(0, 1);
  S2CellIdForwardIterator it{r.begin, r.end};
  for (int i = 0; i < 4; ++i) it.next();
  EXPECT_EQ(it.next(), std::nullopt);
  // Calling next() again is safe.
  EXPECT_EQ(it.next(), std::nullopt);
}

// ---------------------------------------------------------------------------
// S2CellIdReverseIterator
// ---------------------------------------------------------------------------

TEST(S2CellIdReverseIteratorTest, IteratesAllCellsInReverse) {
  S2CellIdRange r = FaceRange(0, 1);
  std::vector<S2CellId> cells = CollectReverse(r);
  ASSERT_EQ(cells.size(), 4u);
  for (int i = 0; i < 4; ++i) {
    EXPECT_EQ(cells[i], r.begin.advance(3 - i));
  }
}

TEST(S2CellIdReverseIteratorTest, EmptyRangeProducesNoItems) {
  S2CellId face = S2CellId::FromFace(0);
  S2CellIdRange r{face.child_begin(1), face.child_begin(1)};
  EXPECT_TRUE(CollectReverse(r).empty());
}

TEST(S2CellIdReverseIteratorTest, SingleCellRange) {
  S2CellIdRange r = FaceRange(0, 1);
  S2CellIdRange single = r.slice(2, 3);
  std::vector<S2CellId> cells = CollectReverse(single);
  ASSERT_EQ(cells.size(), 1u);
  EXPECT_EQ(cells[0], single.at(0));
}

TEST(S2CellIdReverseIteratorTest, ReturnsNulloptWhenExhausted) {
  S2CellIdRange r = FaceRange(0, 1);
  S2CellIdReverseIterator it{r.end, r.begin};
  for (int i = 0; i < 4; ++i) it.next();
  EXPECT_EQ(it.next(), std::nullopt);
  EXPECT_EQ(it.next(), std::nullopt);
}

// ---------------------------------------------------------------------------
// Round-trip: forward then reverse produces symmetric sequences
// ---------------------------------------------------------------------------

TEST(S2CellIdRangeTest, ForwardAndReverseAreSymmetric) {
  S2CellIdRange r = FaceRange(2, 2);  // 16 cells on face 2
  std::vector<S2CellId> fwd = Collect(r);
  std::vector<S2CellId> rev = CollectReverse(r);
  ASSERT_EQ(fwd.size(), rev.size());
  for (size_t i = 0; i < fwd.size(); ++i) {
    EXPECT_EQ(fwd[i], rev[fwd.size() - 1 - i]);
  }
}

}  // namespace
