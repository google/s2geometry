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

// Author: ericv@google.com (Eric Veach)

#include "s2/s2padded_cell.h"

#include <algorithm>
#include <cmath>
#include <string>

#include <gtest/gtest.h>
#include "absl/log/log_streamer.h"
#include "absl/random/bit_gen_ref.h"
#include "absl/random/random.h"
#include "s2/r1interval.h"
#include "s2/r2.h"
#include "s2/r2rect.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2point.h"
#include "s2/s2random.h"
#include "s2/s2testing.h"

using std::min;

namespace {

void CompareS2CellToPadded(const S2Cell& cell, const S2PaddedCell& pcell,
                           double padding) {
  EXPECT_EQ(cell.id(), pcell.id());
  EXPECT_EQ(cell.level(), pcell.level());
  EXPECT_EQ(padding, pcell.padding());
  EXPECT_EQ(cell.GetBoundUV().Expanded(padding), pcell.bound());
  R2Point center_uv = cell.id().GetCenterUV();
  EXPECT_EQ(R2Rect::FromPoint(center_uv).Expanded(padding), pcell.middle());
  EXPECT_EQ(cell.GetCenter(), pcell.GetCenter());
}

TEST(S2PaddedCell, S2CellMethods) {
  // Test the S2PaddedCell methods that have approximate S2Cell equivalents.
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "S2CELL_METHODS",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  constexpr int kIters = 1000;
  for (int iter = 0; iter < kIters; ++iter) {
    S2CellId id = s2random::CellId(bitgen);
    double padding = s2random::LogUniform(bitgen, 1e-15, 1.0);
    S2Cell cell(id);
    S2PaddedCell pcell(id, padding);
    CompareS2CellToPadded(cell, pcell, padding);
    if (!id.is_leaf()) {
      S2Cell children[4];
      EXPECT_TRUE(cell.Subdivide(children));
      for (int pos = 0; pos < 4; ++pos) {
        int i, j;
        pcell.GetChildIJ(pos, &i, &j);
        CompareS2CellToPadded(children[pos], S2PaddedCell(pcell, i, j),
                              padding);
      }
    }
  }
}

TEST(S2PaddedCell, GetEntryExitVertices) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "ENTRY_EXIT_VERTICES",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  constexpr int kIters = 1000;
  for (int iter = 0; iter < kIters; ++iter) {
    S2CellId id = s2random::CellId(bitgen);
    // Check that entry/exit vertices do not depend on padding.
    EXPECT_EQ(S2PaddedCell(id, 0).GetEntryVertex(),
              S2PaddedCell(id, 0.5).GetEntryVertex());
    EXPECT_EQ(S2PaddedCell(id, 0).GetExitVertex(),
              S2PaddedCell(id, 0.5).GetExitVertex());

    // Check that the exit vertex of one cell is the same as the entry vertex
    // of the immediately following cell.  (This also tests wrapping from the
    // end to the start of the S2CellId curve with high probability.)
    EXPECT_EQ(S2PaddedCell(id, 0).GetExitVertex(),
              S2PaddedCell(id.next_wrap(), 0).GetEntryVertex());

    // Check that the entry vertex of a cell is the same as the entry vertex
    // of its first child, and similarly for the exit vertex.
    if (!id.is_leaf()) {
      EXPECT_EQ(S2PaddedCell(id, 0).GetEntryVertex(),
                S2PaddedCell(id.child(0), 0).GetEntryVertex());
      EXPECT_EQ(S2PaddedCell(id, 0).GetExitVertex(),
                S2PaddedCell(id.child(3), 0).GetExitVertex());
    }
  }
}

static double SampleInterval(absl::BitGenRef bitgen, const R1Interval& x) {
  return absl::Uniform<double>(bitgen, x.lo(), x.hi());
}

TEST(S2PaddedCell, ShrinkToFit) {
  constexpr int kIters = 1000;
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "SHRINK_TO_FIT",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int iter = 0; iter < kIters; ++iter) {
    // Start with the desired result and work backwards.
    S2CellId result = s2random::CellId(bitgen);
    R2Rect result_uv = result.GetBoundUV();
    R2Point size_uv = result_uv.GetSize();

    // Find the biggest rectangle that fits in "result" after padding.
    // (These calculations ignore numerical errors.)
    double max_padding = 0.5 * min(size_uv[0], size_uv[1]);
    double padding = absl::Uniform(bitgen, 0.0, max_padding);
    R2Rect max_rect = result_uv.Expanded(-padding);

    // Start with a random subset of the maximum rectangle.
    double ax = SampleInterval(bitgen, max_rect[0]);
    double ay = SampleInterval(bitgen, max_rect[1]);
    double bx = SampleInterval(bitgen, max_rect[0]);
    double by = SampleInterval(bitgen, max_rect[1]);
    R2Point a(ax, ay);
    R2Point b(bx, by);
    if (!result.is_leaf()) {
      // If the result is not a leaf cell, we must ensure that no child of
      // "result" also satisfies the conditions of ShrinkToFit().  We do this
      // by ensuring that "rect" intersects at least two children of "result"
      // (after padding).
      int axis = absl::Uniform(bitgen, 0, 2);
      double center = result.GetCenterUV()[axis];

      // Find the range of coordinates that are shared between child cells
      // along that axis.
      R1Interval shared(center - padding, center + padding);
      double mid = SampleInterval(bitgen, shared.Intersection(max_rect[axis]));
      a[axis] = SampleInterval(bitgen, R1Interval(max_rect[axis].lo(), mid));
      b[axis] = SampleInterval(bitgen, R1Interval(mid, max_rect[axis].hi()));
    }
    R2Rect rect = R2Rect::FromPointPair(a, b);

    // Choose an arbitrary ancestor as the S2PaddedCell.
    S2CellId initial_id =
        result.parent(absl::Uniform(bitgen, 0, result.level() + 1));
    EXPECT_EQ(result, S2PaddedCell(initial_id, padding).ShrinkToFit(rect));
  }
}

}  // namespace
