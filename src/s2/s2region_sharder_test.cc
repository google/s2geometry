// Copyright Google Inc. All Rights Reserved.
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

#include "s2/s2region_sharder.h"

#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/types/span.h"
#include "s2/s2cell_id.h"
#include "s2/s2cell_index.h"
#include "s2/s2cell_union.h"

namespace {

using std::vector;

class S2RegionSharderTest : public ::testing::Test {
 public:
  S2CellIndex IndexFromCoverings(absl::Span<const S2CellUnion> coverings) {
    S2CellIndex index;
    for (int i = 0; i < coverings.size(); ++i) {
      index.Add(coverings[i], i);
    }
    index.Build();
    return index;
  }
};

TEST_F(S2RegionSharderTest, GetMostIntersectingShard) {
  vector<S2CellUnion> coverings{
      S2CellUnion({S2CellId::FromFacePosLevel(0, 0, 10)}),
      S2CellUnion({
          S2CellId::FromFacePosLevel(1, 1, 9),
          S2CellId::FromFacePosLevel(3, 0, 8),
      }),
      S2CellUnion({S2CellId::FromFacePosLevel(5, 0, 10)}),
  };

  const auto Run = [](const S2RegionSharder& sharder) {
    // Overlap with only 1 shard
    EXPECT_EQ(0, sharder.GetMostIntersectingShard(
                     S2CellUnion({
                         S2CellId::FromFacePosLevel(0, 0, 11),
                     }),
                     42));

    // Overlap with multiple shards, picks the shard with more overlap.
    EXPECT_EQ(1, sharder.GetMostIntersectingShard(
                     S2CellUnion({
                         S2CellId::FromFacePosLevel(0, 0, 10),
                         S2CellId::FromFacePosLevel(3, 0, 9),
                         S2CellId::FromFacePosLevel(3, 1, 9),
                     }),
                     42));

    // Overlap with no shards.
    EXPECT_EQ(42, sharder.GetMostIntersectingShard(
                      S2CellUnion({S2CellId::FromFacePosLevel(4, 0, 10)}), 42));
  };

  // Run with an internal and external index.
  S2CellIndex index = IndexFromCoverings(coverings);
  Run(S2RegionSharder(&index));
  Run(S2RegionSharder(coverings));
}

TEST_F(S2RegionSharderTest, GetIntersectingShards) {
  vector<S2CellUnion> coverings{
      S2CellUnion({S2CellId::FromFacePosLevel(0, 0, 10)}),
      S2CellUnion({
          S2CellId::FromFacePosLevel(1, 1, 9),
          S2CellId::FromFacePosLevel(3, 0, 8),
      }),
      S2CellUnion({S2CellId::FromFacePosLevel(5, 0, 10)}),
  };

  const auto Run = [](const S2RegionSharder& sharder) {
    // Overlap with only 1 shard
    EXPECT_THAT(sharder.GetIntersectingShards(S2CellUnion({
                    S2CellId::FromFacePosLevel(0, 0, 11),
                })),
                testing::UnorderedElementsAre(0));

    // Overlap with multiple shards, picks the shard with more overlap.
    EXPECT_THAT(sharder.GetIntersectingShards(S2CellUnion({
                    S2CellId::FromFacePosLevel(0, 0, 10),
                    S2CellId::FromFacePosLevel(3, 0, 9),
                    S2CellId::FromFacePosLevel(3, 1, 9),
                })),
                testing::UnorderedElementsAre(0, 1));

    // Overlap with no shards.
    EXPECT_THAT(sharder.GetIntersectingShards(
                    S2CellUnion({S2CellId::FromFacePosLevel(4, 0, 10)})),
                testing::IsEmpty());
  };

  // Run with an internal and external index.
  S2CellIndex index = IndexFromCoverings(coverings);
  Run(S2RegionSharder(&index));
  Run(S2RegionSharder(coverings));
}

}  // namespace
