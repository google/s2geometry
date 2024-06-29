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

#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "absl/container/flat_hash_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/flags/flag.h"
#include "absl/log/log_streamer.h"
#include "absl/random/bit_gen_ref.h"
#include "absl/random/random.h"
#include "absl/strings/str_format.h"
#include "absl/strings/string_view.h"

#include "s2/base/commandlineflags.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2cell_union.h"
#include "s2/s2latlng.h"
#include "s2/s2point.h"
#include "s2/s2random.h"
#include "s2/s2region_coverer.h"
#include "s2/s2testing.h"

using absl::string_view;
using std::string;
using std::vector;

S2_DEFINE_int32(iters, 400, "number of iterations for testing");

namespace {

enum class QueryType { POINT, CAP };

void TestRandomCaps(absl::BitGenRef bitgen,
                    const S2RegionTermIndexer::Options& options,
                    QueryType query_type) {
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
    if (options.index_contains_points_only()) {
      cap = S2Cap::FromPoint(s2random::Point(bitgen));
      terms = indexer.GetIndexTerms(cap.center(), "");
    } else {
      cap =
          s2random::Cap(bitgen, 0.3 * S2Cell::AverageArea(options.max_level()),
                        4.0 * S2Cell::AverageArea(options.min_level()));
      terms = indexer.GetIndexTerms(cap, "");
    }
    caps.push_back(cap);
    coverings.push_back(coverer.GetCovering(cap));
    for (string_view term : terms) {
      index[term].push_back(i);
    }
    index_terms += terms.size();
  }
  for (int i = 0; i < absl::GetFlag(FLAGS_iters); ++i) {
    // Choose the region to be queried: either a random point or a cap of
    // random size.
    S2Cap cap;
    vector<string> terms;
    if (query_type == QueryType::POINT) {
      cap = S2Cap::FromPoint(s2random::Point(bitgen));
      terms = indexer.GetQueryTerms(cap.center(), "");
    } else {
      cap =
          s2random::Cap(bitgen, 0.3 * S2Cell::AverageArea(options.max_level()),
                        4.0 * S2Cell::AverageArea(options.min_level()));
      terms = indexer.GetQueryTerms(cap, "");
    }
    // Compute the expected results of the S2Cell query by brute force.
    S2CellUnion covering = coverer.GetCovering(cap);
    absl::flat_hash_set<int> expected, actual;
    for (int j = 0; j < caps.size(); ++j) {
      if (covering.Intersects(coverings[j])) {
        expected.insert(j);
      }
    }
    for (string_view term : terms) {
      actual.insert(index[term].begin(), index[term].end());
    }
    EXPECT_EQ(expected, actual);
    query_terms += terms.size();
  }
  absl::PrintF("Index terms/doc: %.2f,  Query terms/doc: %.2f\n",
               static_cast<double>(index_terms) / absl::GetFlag(FLAGS_iters),
               static_cast<double>(query_terms) / absl::GetFlag(FLAGS_iters));
}

// We run one test case for each combination of space vs. time optimization,
// and indexing regions vs. only points.

TEST(S2RegionTermIndexer, IndexRegionsQueryRegionsOptimizeTime) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "INDEX_REGIONS_QUERY_REGIONS_OPTIMIZE_TIME",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  S2RegionTermIndexer::Options options;
  options.set_optimize_for_space(false);       // Optimize for time.
  options.set_min_level(0);                    // Use face cells.
  options.set_max_level(16);
  options.set_max_cells(20);
  TestRandomCaps(bitgen, options, QueryType::CAP);
}

TEST(S2RegionTermIndexer, IndexRegionsQueryPointsOptimizeTime) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "INDEX_REGIONS_QUERY_POINTS_OPTIMIZE_TIME",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  S2RegionTermIndexer::Options options;
  options.set_optimize_for_space(false);       // Optimize for time.
  options.set_min_level(0);                    // Use face cells.
  options.set_max_level(16);
  options.set_max_cells(20);
  TestRandomCaps(bitgen, options, QueryType::POINT);
}

TEST(S2RegionTermIndexer, IndexRegionsQueryRegionsOptimizeTimeWithLevelMod) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "INDEX_REGIONS_QUERY_REGIONS_OPTIMIZE_TIME_WITH_LEVEL_MOD",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  S2RegionTermIndexer::Options options;
  options.set_optimize_for_space(false);       // Optimize for time.
  options.set_min_level(6);                    // Constrain min/max levels.
  options.set_max_level(12);
  options.set_level_mod(3);
  TestRandomCaps(bitgen, options, QueryType::CAP);
}

TEST(S2RegionTermIndexer, IndexRegionsQueryRegionsOptimizeSpace) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "INDEX_REGIONS_QUERY_REGIONS_OPTIMIZE_SPACE",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  S2RegionTermIndexer::Options options;
  options.set_optimize_for_space(true);        // Optimize for space.
  options.set_min_level(4);
  options.set_max_level(S2CellId::kMaxLevel);  // Use leaf cells.
  options.set_max_cells(8);
  TestRandomCaps(bitgen, options, QueryType::CAP);
}

TEST(S2RegionTermIndexer, IndexPointsQueryRegionsOptimizeTime) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "INDEX_POINTS_QUERY_REGIONS_OPTIMIZE_TIME",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  S2RegionTermIndexer::Options options;
  options.set_optimize_for_space(false);       // Optimize for time.
  options.set_min_level(0);                    // Use face cells.
  options.set_max_level(S2CellId::kMaxLevel);
  options.set_level_mod(2);
  options.set_max_cells(20);
  options.set_index_contains_points_only(true);
  TestRandomCaps(bitgen, options, QueryType::CAP);
}

TEST(S2RegionTermIndexer, IndexPointsQueryRegionsOptimizeSpace) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "INDEX_POINTS_QUERY_REGIONS_OPTIMIZE_SPACE",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  S2RegionTermIndexer::Options options;
  options.set_optimize_for_space(true);        // Optimize for space.
  options.set_index_contains_points_only(true);
  // Use default parameter values.
  TestRandomCaps(bitgen, options, QueryType::CAP);
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

  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "MAX_LEVEL_SET_LOOSELY",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  S2Point point = s2random::Point(bitgen);
  EXPECT_EQ(indexer1.GetIndexTerms(point, ""),
            indexer2.GetIndexTerms(point, ""));
  EXPECT_EQ(indexer1.GetQueryTerms(point, ""),
            indexer2.GetQueryTerms(point, ""));

  S2Cap cap = s2random::Cap(bitgen, 0.0, 1.0);  // Area range.
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
