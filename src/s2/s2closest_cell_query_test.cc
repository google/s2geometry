// Copyright 2018 Google Inc. All Rights Reserved.
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

#include "s2/s2closest_cell_query.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "absl/flags/flag.h"
#include "absl/log/absl_check.h"
#include "absl/log/log_streamer.h"
#include "absl/random/bit_gen_ref.h"
#include "absl/random/random.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_format.h"

#include "s2/mutable_s2shape_index.h"
#include "s2/s1angle.h"
#include "s2/s1chord_angle.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2cell_index.h"
#include "s2/s2cell_union.h"
#include "s2/s2closest_cell_query_base.h"
#include "s2/s2closest_edge_query_testing.h"
#include "s2/s2edge_distances.h"
#include "s2/s2latlng.h"
#include "s2/s2latlng_rect.h"
#include "s2/s2metrics.h"
#include "s2/s2min_distance_targets.h"
#include "s2/s2point.h"
#include "s2/s2random.h"
#include "s2/s2region.h"
#include "s2/s2region_coverer.h"
#include "s2/s2testing.h"
#include "s2/s2text_format.h"

namespace {

using s2testing::FractalLoopShapeIndexFactory;
using s2textformat::MakeCellIdOrDie;
using s2textformat::MakePointOrDie;
using std::make_unique;
using std::pair;
using std::unique_ptr;
using std::vector;

using LabelledCell = S2CellIndex::LabelledCell;

TEST(S2ClosestCellQuery, NoCells) {
  S2CellIndex index;
  index.Build();
  S2ClosestCellQuery query(&index);
  S2ClosestCellQuery::PointTarget target(S2Point(1, 0, 0));
  const auto result = query.FindClosestCell(&target);
  EXPECT_EQ(S1ChordAngle::Infinity(), result.distance());
  EXPECT_EQ(S2CellId::None(), result.cell_id());
  EXPECT_EQ(-1, result.label());
  EXPECT_TRUE(result.is_empty());
  EXPECT_EQ(S1ChordAngle::Infinity(), query.GetDistance(&target));
}

TEST(S2ClosestCellQuery, OptionsNotModified) {
  // Tests that FindClosestCell(), GetDistance(), and IsDistanceLess() do not
  // modify query.options(), even though all of these methods have their own
  // specific options requirements.
  S2ClosestCellQuery::Options options;
  options.set_max_results(3);
  options.set_max_distance(S1ChordAngle::Degrees(3));
  options.set_max_error(S1ChordAngle::Degrees(0.001));
  S2CellIndex index;
  index.Add(S2CellId(MakePointOrDie("1:1")), 1);
  index.Add(S2CellId(MakePointOrDie("1:2")), 2);
  index.Add(S2CellId(MakePointOrDie("1:3")), 3);
  index.Build();
  S2ClosestCellQuery query(&index, options);
  S2ClosestCellQuery::PointTarget target(MakePointOrDie("2:2"));
  EXPECT_EQ(2, query.FindClosestCell(&target).label());
  EXPECT_NEAR(1.0, query.GetDistance(&target).degrees(), 1e-7);
  EXPECT_TRUE(query.IsDistanceLess(&target, S1ChordAngle::Degrees(1.5)));

  // Verify that none of the options above were modified.
  EXPECT_EQ(options.max_results(), query.options().max_results());
  EXPECT_EQ(options.max_distance(), query.options().max_distance());
  EXPECT_EQ(options.max_error(), query.options().max_error());
}

TEST(S2ClosestCellQuery, OptionsS1AngleSetters) {
  // Verify that the S1Angle and S1ChordAngle versions do the same thing.
  // This is mainly to prevent the (so far unused) S1Angle versions from
  // being detected as dead code.
  S2ClosestCellQuery::Options angle_options, chord_angle_options;
  angle_options.set_max_distance(S1Angle::Degrees(1));
  chord_angle_options.set_max_distance(S1ChordAngle::Degrees(1));
  EXPECT_EQ(chord_angle_options.max_distance(), angle_options.max_distance());

  angle_options.set_inclusive_max_distance(S1Angle::Degrees(1));
  chord_angle_options.set_inclusive_max_distance(S1ChordAngle::Degrees(1));
  EXPECT_EQ(chord_angle_options.max_distance(), angle_options.max_distance());

  angle_options.set_conservative_max_distance(S1Angle::Degrees(1));
  chord_angle_options.set_conservative_max_distance(S1ChordAngle::Degrees(1));
  EXPECT_EQ(chord_angle_options.max_distance(), angle_options.max_distance());
}

TEST(S2ClosestCellQuery, DistanceEqualToLimit) {
  // Tests the behavior of IsDistanceLess, IsDistanceLessOrEqual, and
  // IsConservativeDistanceLessOrEqual (and the corresponding Options) when
  // the distance to the target exactly equals the chosen limit.
  S2CellId id0(MakePointOrDie("23:12")), id1(MakePointOrDie("47:11"));
  S2CellIndex index;
  index.Add(id0, 0);
  index.Build();
  S2ClosestCellQuery query(&index);

  // Start with two identical cells and a zero distance.
  S2ClosestCellQuery::CellTarget target0(S2Cell{id0});
  S1ChordAngle dist0 = S1ChordAngle::Zero();
  EXPECT_FALSE(query.IsDistanceLess(&target0, dist0));
  EXPECT_TRUE(query.IsDistanceLessOrEqual(&target0, dist0));
  EXPECT_TRUE(query.IsConservativeDistanceLessOrEqual(&target0, dist0));

  // Now try two cells separated by a non-zero distance.
  S2ClosestCellQuery::CellTarget target1(S2Cell{id1});
  S1ChordAngle dist1 = S2Cell(id0).GetDistance(S2Cell(id1));
  EXPECT_FALSE(query.IsDistanceLess(&target1, dist1));
  EXPECT_TRUE(query.IsDistanceLessOrEqual(&target1, dist1));
  EXPECT_TRUE(query.IsConservativeDistanceLessOrEqual(&target1, dist1));
}

TEST(S2ClosestCellQuery, TargetPointInsideIndexedCell) {
  // Tests a target point in the interior of an indexed cell.
  S2CellId cell_id = MakeCellIdOrDie("4/012");
  S2CellIndex index;
  index.Add(cell_id, 1);
  index.Build();
  S2ClosestCellQuery query(&index);
  S2ClosestCellQuery::PointTarget target(cell_id.ToPoint());
  auto result = query.FindClosestCell(&target);
  EXPECT_EQ(S1ChordAngle::Zero(), result.distance());
  EXPECT_EQ(cell_id, result.cell_id());
  EXPECT_EQ(1, result.label());
}

TEST(S2ClosestCellQuery, EmptyTargetOptimized) {
  // Ensure that the optimized algorithm handles empty targets when a distance
  // limit is specified.
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "EMPTY_TARGET_OPTIMIZED", absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  S2CellIndex index;
  for (int i = 0; i < 1000; ++i) {
    index.Add(s2random::CellId(bitgen), i);
  }
  index.Build();
  S2ClosestCellQuery query(&index);
  query.mutable_options()->set_max_distance(S1Angle::Radians(1e-5));
  MutableS2ShapeIndex target_index;
  S2ClosestCellQuery::ShapeIndexTarget target(&target_index);
  EXPECT_EQ(0, query.FindClosestCells(&target).size());
}

TEST(S2ClosestCellQuery, EmptyCellUnionTarget) {
  // Verifies that distances are measured correctly to empty S2CellUnion
  // targets.
  S2ClosestCellQuery::CellUnionTarget target(S2CellUnion{});

  S2CellIndex empty_index;
  empty_index.Build();
  S2ClosestCellQuery empty_query(&empty_index);
  EXPECT_EQ(S1ChordAngle::Infinity(), empty_query.GetDistance(&target));

  S2CellIndex one_cell_index;
  one_cell_index.Add(MakeCellIdOrDie("1/123123"), 1);
  one_cell_index.Build();
  S2ClosestCellQuery one_cell_query(&one_cell_index);
  EXPECT_EQ(S1ChordAngle::Infinity(), one_cell_query.GetDistance(&target));
}

// An abstract class that adds cells to an S2CellIndex for benchmarking.
struct CellIndexFactory {
 public:
  virtual ~CellIndexFactory() = default;

  // Requests that approximately "num_cells" cells located within the given
  // S2Cap bound should be added to "index".  Uses "bitgen" for randomness.
  virtual void AddCells(const S2Cap& index_cap, int num_cells,
                        absl::BitGenRef bitgen, S2CellIndex* index) const = 0;
};

// Generates a cloud of points that approximately fills the given S2Cap, and
// adds a leaf S2CellId for each one.
class PointCloudCellIndexFactory : public CellIndexFactory {
 public:
  void AddCells(const S2Cap& index_cap, int num_cells, absl::BitGenRef bitgen,
                S2CellIndex* index) const override {
    for (int i = 0; i < num_cells; ++i) {
      index->Add(S2CellId(s2random::SamplePoint(bitgen, index_cap)), i);
    }
  }
};

// Generates a collection of S2Caps that are approximately within the given
// "index_cap", generates a covering with "max_cells_per_cap" for each one,
// and adds the coverings to the index.  The radius of each cap is chosen
// randomly such that the total area of the coverings is approximately
// "cap_density" times the area of "index_cap".  In other words, a random
// point inside "index_cap" is likely to intersect about "cap_density"
// coverings (within a factor of 2 or so).
class CapsCellIndexFactory : public CellIndexFactory {
 public:
  CapsCellIndexFactory(int max_cells_per_cap, double cap_density)
      : max_cells_per_cap_(max_cells_per_cap),
        cap_density_(cap_density) {}

  void AddCells(const S2Cap& index_cap, int num_cells, absl::BitGenRef bitgen,
                S2CellIndex* index) const override {
    // All of this math is fairly approximate, since the coverings don't have
    // exactly the given number of cells, etc.
    int num_caps = (num_cells - 1) / max_cells_per_cap_ + 1;
    double max_area = index_cap.GetArea() * cap_density_ / num_caps;
    for (int i = 0; i < num_caps; ++i) {
      // The coverings are bigger than the caps, so we compensate for this by
      // choosing the cap area randomly up to the limit value.
      auto cap = S2Cap::FromCenterArea(s2random::SamplePoint(bitgen, index_cap),
                                       absl::Uniform(bitgen, 0.0, max_area));
      S2RegionCoverer coverer;
      coverer.mutable_options()->set_max_cells(max_cells_per_cap_);
      index->Add(coverer.GetCovering(cap), i);
    }
  }

 private:
  int max_cells_per_cap_;
  double cap_density_;
};

// The approximate radius of S2Cap from which query cells are chosen.
static const S1Angle kTestCapRadius = S2Testing::KmToAngle(10);

using TestingResult = pair<S2MinDistance, LabelledCell>;

// Use "query" to find the closest cell(s) to the given target, and extract
// the query results into the given vector.  Also verify that the results
// satisfy the search criteria.
static void GetClosestCells(S2ClosestCellQuery::Target* target,
                            S2ClosestCellQuery* query,
                            vector<TestingResult>* results) {
  const auto query_results = query->FindClosestCells(target);
  EXPECT_LE(query_results.size(), query->options().max_results());
  const S2Region* region = query->options().region();
  if (!region && query->options().max_distance() == S1ChordAngle::Infinity()) {
    // We can predict exactly how many cells should be returned.
    EXPECT_EQ(std::min(query->options().max_results(),
                       query->index().num_cells()),
              query_results.size());
  }
  for (const auto& result : query_results) {
    // Check that the cell satisfies the region() condition.
    if (region) EXPECT_TRUE(region->MayIntersect(S2Cell(result.cell_id())));

    // Check that it satisfies the max_distance() condition.
    EXPECT_LT(result.distance(), query->options().max_distance());
    results->push_back(TestingResult(
        result.distance(), LabelledCell{result.cell_id(), result.label()}));
  }
}

static void TestFindClosestCells(S2ClosestCellQuery::Target* target,
                                 S2ClosestCellQuery *query) {
  vector<TestingResult> expected, actual;
  query->mutable_options()->set_use_brute_force(true);
  GetClosestCells(target, query, &expected);
  query->mutable_options()->set_use_brute_force(false);
  GetClosestCells(target, query, &actual);
  EXPECT_TRUE(CheckDistanceResults(expected, actual,
                                   query->options().max_results(),
                                   query->options().max_distance(),
                                   query->options().max_error()))
      << "max_results=" << query->options().max_results()
      << ", max_distance=" << query->options().max_distance()
      << ", max_error=" << query->options().max_error();

  if (expected.empty()) return;

  // Note that when options.max_error() > 0, expected[0].distance() may not be
  // the minimum distance.  It is never larger by more than max_error(), but
  // the actual value also depends on max_results().
  //
  // Here we verify that GetDistance() and IsDistanceLess() return results
  // that are consistent with the max_error() setting.
  S1ChordAngle max_error = query->options().max_error();
  S1ChordAngle min_distance = expected[0].first;
  EXPECT_LE(query->GetDistance(target), min_distance + max_error);

  // Test IsDistanceLess().
  EXPECT_FALSE(query->IsDistanceLess(target, min_distance - max_error));
  EXPECT_TRUE(query->IsConservativeDistanceLessOrEqual(target, min_distance));
}

// The running time of this test is proportional to
//    (num_indexes + num_queries) * num_cells.
// (Note that every query is checked using the brute force algorithm.)
static void TestWithIndexFactory(const CellIndexFactory& factory,
                                 int num_indexes, int num_cells,
                                 int num_queries, absl::BitGenRef bitgen) {
  // Build a set of S2CellIndexes containing the desired geometry.
  vector<S2Cap> index_caps;
  vector<unique_ptr<S2CellIndex>> indexes;
  for (int i = 0; i < num_indexes; ++i) {
    index_caps.push_back(S2Cap(s2random::Point(bitgen), kTestCapRadius));
    auto index = make_unique<S2CellIndex>();
    factory.AddCells(index_caps.back(), num_cells, bitgen, index.get());
    index->Build();
    indexes.push_back(std::move(index));
  }
  for (int i = 0; i < num_queries; ++i) {
    int i_index = absl::Uniform(bitgen, 0, num_indexes);
    const S2Cap& index_cap = index_caps[i_index];

    // Choose query points from an area approximately 4x larger than the
    // geometry being tested.
    S1Angle query_radius = 2 * index_cap.GetRadius();
    S2Cap query_cap(index_cap.center(), query_radius);
    S2ClosestCellQuery query(indexes[i_index].get());

    // Occasionally we don't set any limit on the number of result cells.
    // (This may return all cells if we also don't set a distance limit.)
    if (absl::Bernoulli(bitgen, 0.9)) {
      query.mutable_options()->set_max_results(absl::Uniform(bitgen, 1, 11));
    }
    // We set a distance limit 2/3 of the time.
    if (absl::Bernoulli(bitgen, 2.0 / 3)) {
      query.mutable_options()->set_max_distance(
          absl::Uniform(bitgen, 0.0, 1.0) * query_radius);
    }
    if (absl::Bernoulli(bitgen, 0.5)) {
      // Choose a maximum error whose logarithm is uniformly distributed over
      // a reasonable range, except that it is sometimes zero.
      query.mutable_options()->set_max_error(S1Angle::Radians(
          s2random::LogUniform(bitgen, 1e-4, 1.0) * query_radius.radians()));
    }
    S2LatLngRect filter_rect = S2LatLngRect::FromCenterSize(
        S2LatLng(s2random::SamplePoint(bitgen, query_cap)),
        S2LatLng(absl::Uniform(bitgen, 0.0, 1.0) * kTestCapRadius,
                 absl::Uniform(bitgen, 0.0, 1.0) * kTestCapRadius));
    if (absl::Bernoulli(bitgen, 0.2)) {
      query.mutable_options()->set_region(&filter_rect);
    }
    int target_type = absl::Uniform(bitgen, 0, 5);
    if (target_type == 0) {
      // Find the cells closest to a given point.
      S2Point point = s2random::SamplePoint(bitgen, query_cap);
      S2ClosestCellQuery::PointTarget target(point);
      TestFindClosestCells(&target, &query);
    } else if (target_type == 1) {
      // Find the cells closest to a given edge.
      S2Point a = s2random::SamplePoint(bitgen, query_cap);
      S2Point b = s2random::SamplePoint(
          bitgen,
          S2Cap(a, s2random::LogUniform(bitgen, 1e-4, 1.0) * query_radius));
      S2ClosestCellQuery::EdgeTarget target(a, b);
      TestFindClosestCells(&target, &query);
    } else if (target_type == 2) {
      // Find the cells closest to a given cell.
      int min_level = S2::kMaxDiag.GetLevelForMaxValue(query_radius.radians());
      int level = absl::Uniform(absl::IntervalClosedClosed, bitgen, min_level,
                                S2CellId::kMaxLevel);
      S2Point a = s2random::SamplePoint(bitgen, query_cap);
      S2Cell cell(S2CellId(a).parent(level));
      S2ClosestCellQuery::CellTarget target(cell);
      TestFindClosestCells(&target, &query);
    } else if (target_type == 3) {
      // Find the cells closest to an S2Cap covering.
      S2Point center = s2random::SamplePoint(bitgen, query_cap);
      S2Cap cap(center, s2random::LogUniform(bitgen, 1e-5, 0.1) * query_radius);
      S2RegionCoverer coverer;
      coverer.mutable_options()->set_max_cells(16);
      S2ClosestCellQuery::CellUnionTarget target(coverer.GetCovering(cap));
      TestFindClosestCells(&target, &query);
    } else {
      ABSL_DCHECK_EQ(4, target_type);
      MutableS2ShapeIndex target_index;
      s2testing::FractalLoopShapeIndexFactory(bitgen).AddEdges(index_cap, 100,
                                                               &target_index);
      S2ClosestCellQuery::ShapeIndexTarget target(&target_index);
      target.set_include_interiors(absl::Bernoulli(bitgen, 0.5));
      TestFindClosestCells(&target, &query);
    }
  }
}

static constexpr int kNumIndexes = 20;
static constexpr int kNumCells = 100;
static constexpr int kNumQueries = 100;

TEST(S2ClosestCellQuery, PointCloudCells) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "POINT_CLOUD_CELLS", absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  TestWithIndexFactory(PointCloudCellIndexFactory(), kNumIndexes, kNumCells,
                       kNumQueries, bitgen);
}

TEST(S2ClosestCellQuery, CapsCells) {
  absl::BitGen bitgen(
      S2Testing::MakeTaggedSeedSeq("CAPS_CELLS", absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  TestWithIndexFactory(
      CapsCellIndexFactory(16 /*max_cells_per_cap*/, 0.1 /*density*/),
      kNumIndexes, kNumCells, kNumQueries, bitgen);
}

TEST(S2ClosestCellQuery, ConservativeCellDistanceIsUsed) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "CONSERVATIVE_CELL_DISTANCE_IS_USED", absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  TestWithIndexFactory(PointCloudCellIndexFactory(), 5, 100, 10, bitgen);
}

}  // namespace
