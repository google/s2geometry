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


#include "s2/s2cell_iterator_join.h"

#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/container/flat_hash_set.h"
#include "s2/mutable_s2shape_index.h"
#include "s2/s2cell.h"
#include "s2/s2cell_iterator_testing.h"
#include "s2/s2point_index.h"
#include "s2/s2polygon.h"
#include "s2/s2testing.h"

namespace {

using ::testing::Contains;
using ::testing::Eq;
using ::testing::Le;

// Cell IDs covering Central Park in New York City
constexpr static absl::string_view kCentralParkATokens[] = {
    "89c2589",  "89c258a1", "89c258a3", "89c258bc",
    "89c258c1", "89c258ec", "89c258f4"};

// Cell Ids also covering Central Park but a subset of CentralParkA.
constexpr static absl::string_view kCentralParkBTokens[] = {
    "89c2589", "89c258a03", "89c258a1c", "89c258a3", "89c258bd", "89c258be1"};

// Builds an S2JoinRow from S2CellId tokens.
std::pair<S2CellId, S2CellId> PairFromTokens(std::string token_a,
                                             std::string token_b) {
  return std::make_pair(S2CellId::FromToken(token_a),
                        S2CellId::FromToken(token_b));
}

// Builds a btree_map from a list of token to use with MockS2CellIterator.
absl::btree_map<S2CellId, int> TokenMap(
    absl::Span<const absl::string_view> tokens) {
  absl::btree_map<S2CellId, int> map;
  int count = 0;
  for (const auto& token : tokens) {
    map.emplace(S2CellId::FromToken(token), ++count);
  }
  return map;
}

TEST(S2CellJoinIterator, HeterogeneousJoinCompiles) {
  MutableS2ShapeIndex shape_index;
  S2PointIndex<int> point_index;
  MakeS2CellIteratorJoin(&shape_index, &point_index)
      .Join([](const MutableS2ShapeIndex::Iterator& itera,
               const S2PointIndex<int>::Iterator& iterb) { return true; });
}

TEST(S2CellJoinIterator, ExactJoinWorks) {
  auto cpa_map = TokenMap(absl::MakeSpan(kCentralParkATokens));
  auto cpb_map = TokenMap(absl::MakeSpan(kCentralParkBTokens));

  std::vector<std::pair<S2CellId, S2CellId>> rows;
  MakeS2CellIteratorJoin(MakeMockS2CellIterator(&cpa_map),
                         MakeMockS2CellIterator(&cpb_map))
      .Join([&rows](const MockS2CellIterator<int>& iter_a,
                    const MockS2CellIterator<int>& iter_b) {
        rows.emplace_back(iter_a.id(), iter_b.id());
        EXPECT_TRUE(rows.back().first.contains(rows.back().second));
        return true;
      });

  EXPECT_THAT(rows.size(), Eq(std::min(cpa_map.size(), cpb_map.size())));

  const std::vector<std::pair<S2CellId, S2CellId>> truth = {
      PairFromTokens("89c2589", "89c2589"),
      PairFromTokens("89c258a1", "89c258a03"),
      PairFromTokens("89c258a1", "89c258a1c"),
      PairFromTokens("89c258a3", "89c258a3"),
      PairFromTokens("89c258bc", "89c258bd"),
      PairFromTokens("89c258bc", "89c258be1")};

  EXPECT_THAT(rows.size(), Eq(truth.size()));
  for (int i = 0; i < truth.size(); ++i) {
    EXPECT_THAT(rows[i], Eq(truth[i]));
  }
}

TEST(S2CellJoinIterator, ExactJoinSeekingWorks) {
  // Sometimes we have to seek more than once to find a pair of cells that
  // overlap.  2d5e3 below doesn't overlap anything in map_b, and so shouldn't
  // be reported.  Instead the join should seek twice and skip it, make sure it
  // does.
  auto map_a = TokenMap({"2d5dd7", "2d5ddc", "2d5e3", "2d5e801", "2d5e803"});
  auto map_b = TokenMap({"2d5d", "2d5e84"});

  const std::vector<std::pair<S2CellId, S2CellId>> truth = {
      PairFromTokens("2d5dd7", "2d5d"),     //
      PairFromTokens("2d5ddc", "2d5d"),     //
      PairFromTokens("2d5e801", "2d5e84"),  //
      PairFromTokens("2d5e803", "2d5e84")   //
  };

  std::vector<std::pair<S2CellId, S2CellId>> rows;
  MakeS2CellIteratorJoin(MakeMockS2CellIterator(&map_a),
                         MakeMockS2CellIterator(&map_b))
      .Join([&rows](const MockS2CellIterator<int>& iter_a,
                    const MockS2CellIterator<int>& iter_b) {
        rows.emplace_back(iter_a.id(), iter_b.id());
        return true;
      });

  EXPECT_THAT(rows.size(), Eq(truth.size()));
  for (int i = 0; i < truth.size(); ++i) {
    EXPECT_THAT(rows[i], Eq(truth[i]));
  }
}

TEST(S2CellJoinIterator, NearJoinWorks) {
  auto cpa_map = TokenMap(absl::MakeSpan(kCentralParkATokens));
  auto cpb_map = TokenMap(absl::MakeSpan(kCentralParkBTokens));

  const S1ChordAngle kTolerance = S1ChordAngle::Degrees(1);

  absl::flat_hash_set<std::pair<S2CellId, S2CellId>> rows;
  MakeS2CellIteratorJoin(MakeMockS2CellIterator(&cpa_map),
                         MakeMockS2CellIterator(&cpb_map), kTolerance)
      .Join([&rows](const MockS2CellIterator<int>& iter_a,
                    const MockS2CellIterator<int>& iter_b) {
        rows.emplace(iter_a.id(), iter_b.id());
        return true;
      });

  const std::vector<std::pair<S2CellId, S2CellId>> truth = {
      PairFromTokens("89c2589", "89c2589"),
      PairFromTokens("89c258a1", "89c258a03"),
      PairFromTokens("89c258a1", "89c258a1c"),
      PairFromTokens("89c258a3", "89c258a3"),
      PairFromTokens("89c258bc", "89c258bd"),
      PairFromTokens("89c258bc", "89c258be1")};

  // The contents of the exact join should be a subset of the rough join.
  for (const auto& row : truth) {
    EXPECT_THAT(rows, Contains(row));
  }

  // All the pairs should be within tolerance of each other.
  for (const auto& row : truth) {
    EXPECT_THAT(S2Cell(row.first).GetDistance(S2Cell(row.second)),
                Le(kTolerance));
  }

  // These cells are more than 0 degrees apart (i.e. not touching) and not part
  // of the exact results, but should be included in the tolerant result.
  const std::vector<std::pair<S2CellId, S2CellId>> tolerant_truth = {
      PairFromTokens("89c258a1", "89c258bd"),
      PairFromTokens("89c258a1", "89c258be1"),
      PairFromTokens("89c258a3", "89c258a03"),
      PairFromTokens("89c258a3", "89c258be1"),
      PairFromTokens("89c258bc", "89c258a03"),
      PairFromTokens("89c258bc", "89c258a1c"),
      PairFromTokens("89c258c1", "89c258a03"),
      PairFromTokens("89c258c1", "89c258a1c"),
      PairFromTokens("89c258c1", "89c258a3"),
      PairFromTokens("89c258c1", "89c258bd"),
      PairFromTokens("89c258c1", "89c258be1"),
      PairFromTokens("89c258ec", "89c258a03"),
      PairFromTokens("89c258ec", "89c258a1c"),
      PairFromTokens("89c258ec", "89c258a3"),
      PairFromTokens("89c258ec", "89c258bd"),
      PairFromTokens("89c258ec", "89c258be1"),
      PairFromTokens("89c258f4", "89c258a03"),
      PairFromTokens("89c258f4", "89c258a1c"),
      PairFromTokens("89c258f4", "89c258a3"),
      PairFromTokens("89c258f4", "89c258bd"),
      PairFromTokens("89c258f4", "89c258be1")};

  for (const auto& row : tolerant_truth) {
    EXPECT_THAT(rows, Contains(row));
  }
}


}  // namespace
