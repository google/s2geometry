// Copyright 2017 Google Inc. All Rights Reserved.
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

#include "s2/s2region_term_indexer.h"

#include <cstdio>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "absl/container/flat_hash_map.h"
#include "absl/flags/flag.h"

#include "s2/base/commandlineflags.h"
#include "s2/base/logging.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2cell_union.h"
#include "s2/s2latlng.h"
#include "s2/s2testing.h"

using std::string;
using std::vector;

S2_DEFINE_int32(iters, 400, "number of iterations for testing");

namespace {

enum DataType { POINT, CAP };

void TestRandomCaps(const S2RegionTermIndexer::Options& options,
                    DataType index_type, DataType query_type) {
  // This function creates an index consisting either of points (if
  // options.index_contains_points_only() is true) or S2Caps of random size.
  // It then executes queries consisting of points (if query_type == POINT)
  // or S2Caps of random size (if query_type == CAP).
  S2RegionTermIndexer indexer(options);
  S2RegionCoverer coverer(options);
  vector<S2Cap> caps;
  vector<S2CellUnion> coverings;
  absl::flat_hash_map<string, vector<int>> index;
  int index_terms = 0, query_terms = 0;
  for (int i = 0; i < absl::GetFlag(FLAGS_iters); ++i) {
    // Choose the region to be indexed: either a single point or a cap
    // of random size (up to a full sphere).
    S2Cap cap;
    vector<string> terms;
    if (index_type == DataType::POINT) {
      cap = S2Cap::FromPoint(S2Testing::RandomPoint());
      terms = indexer.GetIndexTerms(cap.center(), "");
    } else {
      cap = S2Testing::GetRandomCap(
          0.3 * S2Cell::AverageArea(options.max_level()),
          4.0 * S2Cell::AverageArea(options.min_level()));
      terms = indexer.GetIndexTerms(cap, "");
    }
    caps.push_back(cap);
    coverings.push_back(coverer.GetCovering(cap));
    for (const string& term : terms) {
      index[term].push_back(i);
    }
    index_terms += terms.size();
  }
  for (int i = 0; i < absl::GetFlag(FLAGS_iters); ++i) {
    // Choose the region to be queried: either a random point or a cap of
    // random size.
    S2Cap cap;
    vector<string> terms;
    if (query_type == DataType::POINT) {
      cap = S2Cap::FromPoint(S2Testing::RandomPoint());
      terms = indexer.GetQueryTerms(cap.center(), "");
    } else {
      cap = S2Testing::GetRandomCap(
          0.3 * S2Cell::AverageArea(options.max_level()),
          4.0 * S2Cell::AverageArea(options.min_level()));
      terms = indexer.GetQueryTerms(cap, "");
    }
    // Compute the expected results of the S2Cell query by brute force.
    S2CellUnion covering = coverer.GetCovering(cap);
    std::set<int> expected, actual;
    for (int j = 0; j < caps.size(); ++j) {
      if (covering.Intersects(coverings[j])) {
        expected.insert(j);
      }
    }
    for (const string& term : terms) {
      actual.insert(index[term].begin(), index[term].end());
    }
    EXPECT_EQ(expected, actual);
    query_terms += terms.size();
  }
  printf("Index terms/doc: %.2f,  Query terms/doc: %.2f\n",
         static_cast<double>(index_terms) / absl::GetFlag(FLAGS_iters),
         static_cast<double>(query_terms) / absl::GetFlag(FLAGS_iters));
}

using TestCase = std::tuple<DataType, DataType, bool, bool, bool>;

class S2RegionTermIndexerTest : public testing::TestWithParam<TestCase> {
protected:
  void SetUp() override {
    index_type = std::get<0>(GetParam());
    query_type = std::get<1>(GetParam());
    options.set_optimize_for_space(std::get<2>(GetParam()));
    options.set_index_contains_points_only(std::get<3>(GetParam()));
    options.set_query_contains_points_only(std::get<4>(GetParam()));
    if (index_type != DataType::POINT && options.index_contains_points_only()) {
      GTEST_SKIP() << "Case index_type != DataType::POINT && "
                      "options.index_contains_points_only() is invalid.";
    }
    if (query_type != DataType::POINT && options.query_contains_points_only()) {
      GTEST_SKIP() << "Case query_type != DataType::POINT && "
                      "options.query_contains_points_only() is invalid.";
    }
  }

  S2RegionTermIndexer::Options options;
  DataType index_type{};
  DataType query_type{};
};

// We run one test case for each combination:
// index_type: POINT, CAP
// query_type: POINT, CAP
// optimize_for_space: false, true
// index_contains_points_only: false, true
// query_contains_points_only: false, true
INSTANTIATE_TEST_CASE_P(
    S2RegionTermIndexerTests, S2RegionTermIndexerTest,
    testing::Combine(testing::Values(DataType::POINT, DataType::CAP),
                     testing::Values(DataType::POINT, DataType::CAP),
                     testing::Bool(), testing::Bool(), testing::Bool()));

TEST_P(S2RegionTermIndexerTest, DefaultParametersValues) {
  TestRandomCaps(options, index_type, query_type);
}

TEST_P(S2RegionTermIndexerTest, UseFaceCells) {
  options.set_min_level(0);
  options.set_max_level(16);
  options.set_max_cells(20);
  TestRandomCaps(options, index_type, query_type);
}

TEST_P(S2RegionTermIndexerTest, ConstrainMinMaxLevels) {
  options.set_min_level(6);
  options.set_max_level(12);
  options.set_level_mod(3);
  TestRandomCaps(options, index_type, query_type);
}

TEST_P(S2RegionTermIndexerTest, UseLeafCells) {
  options.set_min_level(4);
  options.set_max_level(S2CellId::kMaxLevel);
  options.set_max_cells(8);
  TestRandomCaps(options, index_type, query_type);
}

TEST_P(S2RegionTermIndexerTest, UseFaceCells2) {
  options.set_min_level(0);
  options.set_max_level(S2CellId::kMaxLevel);
  options.set_level_mod(2);
  options.set_max_cells(20);
  TestRandomCaps(options, index_type, query_type);
}

TEST(S2RegionTermIndexer, MarkerCharacter) {
  S2RegionTermIndexer::Options options;
  options.set_min_level(20);
  options.set_max_level(20);

  S2RegionTermIndexer indexer(options);
  S2Point point = S2LatLng::FromDegrees(10, 20).ToPoint();
  EXPECT_EQ(indexer.options().marker_character(), '$');
  EXPECT_EQ(indexer.GetQueryTerms(point, ""),
            vector<string>({"11282087039", "$11282087039"}));

  indexer.mutable_options()->set_marker_character(':');
  EXPECT_EQ(indexer.options().marker_character(), ':');
  EXPECT_EQ(indexer.GetQueryTerms(point, ""),
            vector<string>({"11282087039", ":11282087039"}));
}

TEST(S2RegionTermIndexer, MaxLevelSetLoosely) {
  // Test that correct terms are generated even when (max_level - min_level)
  // is not a multiple of level_mod.
  S2RegionTermIndexer::Options options;
  options.set_min_level(1);
  options.set_level_mod(2);
  options.set_max_level(19);
  S2RegionTermIndexer indexer1(options);
  options.set_max_level(20);
  S2RegionTermIndexer indexer2(options);

  S2Point point = S2Testing::RandomPoint();
  EXPECT_EQ(indexer1.GetIndexTerms(point, ""),
            indexer2.GetIndexTerms(point, ""));
  EXPECT_EQ(indexer1.GetQueryTerms(point, ""),
            indexer2.GetQueryTerms(point, ""));

  S2Cap cap = S2Testing::GetRandomCap(0.0, 1.0);  // Area range.
  EXPECT_EQ(indexer1.GetIndexTerms(cap, ""),
            indexer2.GetIndexTerms(cap, ""));
  EXPECT_EQ(indexer1.GetQueryTerms(cap, ""),
            indexer2.GetQueryTerms(cap, ""));
}

TEST(S2RegionTermIndexer, MoveConstructor) {
  S2RegionTermIndexer x;
  x.mutable_options()->set_max_cells(12345);
  S2RegionTermIndexer y = std::move(x);
  EXPECT_EQ(12345, y.options().max_cells());
}

TEST(S2RegionTermIndexer, MoveAssignmentOperator) {
  S2RegionTermIndexer x;
  x.mutable_options()->set_max_cells(12345);
  S2RegionTermIndexer y;
  y.mutable_options()->set_max_cells(0);
  y = std::move(x);
  EXPECT_EQ(12345, y.options().max_cells());
}

}  // namespace
