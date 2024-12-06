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


#include "s2/internal/s2index_cell_data.h"

#include <memory>

#include <gmock/gmock.h>
#include "absl/log/absl_check.h"
#include "s2/mutable_s2shape_index.h"
#include "s2/s2cell_range_iterator.h"
#include "s2/s2testing.h"
#include "s2/s2text_format.h"

namespace internal {
namespace {

using ::std::unique_ptr;
using ::testing::Eq;
using ::testing::IsEmpty;
using ::testing::Ne;
using ::testing::Not;

TEST(S2IndexCellData, DimensionFilteringWorks) {
  unique_ptr<MutableS2ShapeIndex> index = s2textformat::MakeIndexOrDie(
      "0:0"
      "#1:1, 2:2"
      "#1:0, 0:1, -1:0, 0:-1");

  MutableS2ShapeIndex::Iterator iter(index.get(), S2ShapeIndex::BEGIN);

  {
    // Check that we get all dimensions by default.
    S2IndexCellData data;
    data.LoadCell(index.get(), iter.id(), &iter.cell());
    EXPECT_THAT(data.dim_edges(0), Not(IsEmpty()));
    EXPECT_THAT(data.dim_edges(1), Not(IsEmpty()));
    EXPECT_THAT(data.dim_edges(2), Not(IsEmpty()));
  }

  {
    // No dimensions should work too, we just don't decode edges.
    S2IndexCellData data;
    data.set_dim_wanted(0, false);
    data.set_dim_wanted(1, false);
    data.set_dim_wanted(2, false);
    data.LoadCell(index.get(), iter.id(), &iter.cell());
    EXPECT_THAT(data.dim_edges(0), IsEmpty());
    EXPECT_THAT(data.dim_edges(1), IsEmpty());
    EXPECT_THAT(data.dim_edges(2), IsEmpty());
  }

  // Should be able to get ranges even if a dimension is turned off.
  {
    S2IndexCellData data;
    data.set_dim_wanted(0, false);
    data.set_dim_wanted(1, true);
    data.set_dim_wanted(2, true);
    data.LoadCell(index.get(), iter.id(), &iter.cell());
    EXPECT_THAT(data.dim_range_edges(0, 0), IsEmpty());
    EXPECT_THAT(data.dim_range_edges(0, 2), Not(IsEmpty()));
  }

  {
    S2IndexCellData data;
    data.set_dim_wanted(0, false);
    data.set_dim_wanted(1, true);
    data.set_dim_wanted(2, false);
    data.LoadCell(index.get(), iter.id(), &iter.cell());
    EXPECT_THAT(data.dim_edges(0), IsEmpty());
    EXPECT_THAT(data.dim_edges(1), Not(IsEmpty()));
    EXPECT_THAT(data.dim_edges(2), IsEmpty());
  }

  {
    S2IndexCellData data;
    data.set_dim_wanted(0, true);
    data.set_dim_wanted(1, false);
    data.set_dim_wanted(2, true);
    data.LoadCell(index.get(), iter.id(), &iter.cell());
    EXPECT_THAT(data.dim_edges(0), Not(IsEmpty()));
    EXPECT_THAT(data.dim_edges(1), IsEmpty());
    EXPECT_THAT(data.dim_edges(2), Not(IsEmpty()));
  }

  {
    S2IndexCellData data;
    data.set_dim_wanted(0, true);
    data.set_dim_wanted(1, false);
    data.set_dim_wanted(2, false);
    data.LoadCell(index.get(), iter.id(), &iter.cell());
    EXPECT_THAT(data.dim_edges(0), Not(IsEmpty()));
    EXPECT_THAT(data.dim_edges(1), IsEmpty());
    EXPECT_THAT(data.dim_edges(2), IsEmpty());
  }

  {
    S2IndexCellData data;
    data.set_dim_wanted(0, false);
    data.set_dim_wanted(1, false);
    data.set_dim_wanted(2, true);
    data.LoadCell(index.get(), iter.id(), &iter.cell());
    EXPECT_THAT(data.dim_edges(0), IsEmpty());
    EXPECT_THAT(data.dim_edges(1), IsEmpty());
    EXPECT_THAT(data.dim_edges(2), Not(IsEmpty()));
  }
}

// We cache cell centers and S2Cell instances because they're expensive to
// compute if they're not needed.  Make sure that when we load a new unique
// cell, those values are updated if we access them.
TEST(S2IndexCellData, CellAndCenterRecomputed) {
  // A line between two faces will guarantee we get at least two cells.
  unique_ptr<MutableS2ShapeIndex> index =
      s2textformat::MakeIndexOrDie("# 0:0, 0:-90 #");

  S2IndexCellData data;
  auto iter = MakeS2CellRangeIterator(index.get());

  data.LoadCell(index.get(), iter.id(), &iter.iterator().cell());
  const S2Point center0 = data.center();
  const S2Cell cell0 = data.cell();

  iter.Next();
  ABSL_CHECK(!iter.done());

  data.LoadCell(index.get(), iter.id(), &iter.iterator().cell());
  const S2Point center1 = data.center();
  const S2Cell cell1 = data.cell();

  EXPECT_THAT(cell0, Ne(cell1));
  EXPECT_THAT(center0, Ne(center1));

  // Load the same cell again, nothing should change.
  data.LoadCell(index.get(), iter.id(), &iter.iterator().cell());
  const S2Point center2 = data.center();
  const S2Cell cell2 = data.cell();

  EXPECT_THAT(cell1, Eq(cell2));
  EXPECT_THAT(center1, Eq(center2));
}

}  // namespace
}  // namespace internal
