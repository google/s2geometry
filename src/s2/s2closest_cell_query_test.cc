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

#include <benchmark/benchmark.h>
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

// Description of the benchmarks:
//
//  - All benchmarks index a collection of cells chosen within a random
//    S2Cap (the "index cap") of configurable radius (default 100km).
//
//  - The "PointCloud" benchmarks index a collection of random leaf cells
//    within the index cap, whereas the "CoarseCaps" benchmarks use coarse
//    coverings (up to 8 cells) of randomly sized S2Caps that fill about
//    10% of the index cap.  In both cases, the benchmark parameter
//    immediately following the benchmark name (e.g. /48k) represents the
//    approximate total number of indexed cells.
//
//  - The "FindClosest/x" benchmarks find the closest cell to the target
//    geometry (equivalent to measuring the minimum distance).
//
//  - The "FindClosestWhenZero" benchmarks do the same thing, but they
//    always use target geometry that intersects the indexed cells.
//
//  - "IsDistanceLess/a/b" uses an index with "a" cells, and checks whether
//    the distance to the target geometry is less than a threshold radius
//    "r" chosen such a disc of radius "r" around any point of the target
//    geometry contains an average of approximately 2**b indexed cell
//    center points.  So for example, IsDistanceLess<PointCloud>/48k/-3
//    uses an index of 48k leaf cells, and queries random discs that
//    contain an average of approximately 0.125 leaf cell centers each.
//
//  - "FindCellsNearby/a/b" is like IsDistanceLess, but finds all cells
//    within the threshold radius rather than just testing whether such a
//    cell exists.  For example, BM_FindCellsNearby<PointCloud>/48k/3 uses
//    an index of 48k leaf cells, and queries random discs that contain an
//    average of approximately 8 leaf cell centers each.
//
//  - The target (a point, edge, cell, S2CellUnion, or S2ShapeIndex) is
//    chosen so that it approximately fills a "target cap".  The target cap
//    is randomly chosen such that the radius of "target_cap" and the
//    distance between the centers of "target_cap" and "index_cap" are each
//    some given fraction of the index cap radius.  Most of the benchmarks
//    choose the distance between the centers randomly up to twice the
//    index cap radius.  In other words, 50% of the time the target is
//    inside the index cap, and 50% of the time the target is outside that
//    area.  (Note that even when the target is inside the index cap it
//    does not necessarily intersect any indexed cells.)
//
//  - The benchmarked target geometries are as follows:
//
//    - [No suffix]: point target.
//    - LongEdge: Edges of random length up to 10% of the index cap diameter.
//    - SmallCell: Cells of random size up to 1% of the index cap radius.
//    - SmallCoarseCellUnion: an 8-cell covering of an S2Cap whose radius
//      is up to 1% of the index radius.
//    - SmallFineCellUnion: like the above, but up to a 100-cell covering.
//    - SmallCoarseShapeIndex: a fractal loop with 12 edges, filling an
//      S2Cap whose radius is up to 1% of the index radius.
//    - SmallFindShapeindex: like the above, but the loop has 768 edges.

ABSL_FLAG(bool, bm_use_brute_force, false,
          "Benchmarks: Use the brute force implementation");

ABSL_FLAG(double, bm_radius_km, 100,
          "Benchmarks: Default radius for indexed geometry");

namespace {

enum class TargetType { POINT, EDGE, CELL, CELL_UNION, SHAPE_INDEX };

static S1Angle FractionToRadius(absl::BitGenRef bitgen, double fraction) {
  if (fraction < 0) {
    fraction = absl::Uniform(bitgen, 0.0, -fraction);
  }
  return fraction * S2Testing::KmToAngle(absl::GetFlag(FLAGS_bm_radius_km));
}

static S2Point SampleBoundary(absl::BitGenRef bitgen, const S2Cap& cap) {
  return S2::GetPointOnLine(cap.center(), s2random::Point(bitgen),
                            cap.GetRadius());
}

static S2CellId SampleCell(absl::BitGenRef bitgen, const S2CellIndex& index) {
  int num_cells = 0;
  for (S2CellIndex::CellIterator it(&index); !it.done(); it.Next()) {
    ++num_cells;
  }
  ABSL_DCHECK_GT(num_cells, 0);
  S2CellIndex::CellIterator it(&index);
  for (int i = absl::Uniform(bitgen, 0, num_cells); --i >= 0; it.Next())
    continue;
  return it.cell_id();
}

// Calls FindClosestCells() the given number of times on an S2CellIndex with
// approximately "num_index_cells" cells generated by "factory".  The geometry
// is generated within an S2Cap of the radius specified by
// "FLAGS_bm_radius_km" (the "index radius").  Parameters with "fraction" in
// their names are expressed as a fraction of this radius.
//
// Each query uses a target of the given "target_type".  If "target_type" is
// CELL_UNION, then the union will have approximately "num_target_objects"
// cells, and if "target_type" is SHAPE_INDEX, then the shape index will have
// approximately "num_target_objects" edges.  The number of results is limited
// to "max_results" (usually 1).
//
//   - If max_distance_fraction > 0, then max_distance() is set to the given
//     fraction of the index radius.
//
//   - If max_error_fraction > 0, then max_error() is set to the given
//     fraction of the index radius.
//
//   - If "choose_target_from_index" is true, then the target will be chosen
//     from the geometry in the index itself, otherwise it will be chosen
//     randomly according to the parameters below:
//
//   - If target_radius_fraction > 0, the target radius will be approximately
//     the given fraction of the index radius; if target_radius_fraction < 0,
//     it will be chosen randomly up to corresponding positive fraction.
//
//   - If center_separation_fraction > 0, then the centers of index and target
//     bounding caps will be separated by the given fraction of the index
//     radius; if center_separation_fraction < 0, they will be separated by up
//     to the corresponding positive fraction.
static void BenchmarkFindClosest(
    benchmark::State& state, const CellIndexFactory& factory,
    int num_index_cells, int max_results, double max_distance_fraction,
    double max_error_fraction, TargetType target_type, int num_target_objects,
    bool choose_target_from_index, double target_radius_fraction,
    double center_separation_fraction) {
  const std::string seed_str = absl::StrCat(
      "BENCHMARK_FIND_CLOSEST", absl::GetFlag(FLAGS_s2_random_seed));
  std::seed_seq seed(seed_str.begin(), seed_str.end());
  std::mt19937_64 bitgen(seed);

  // We execute "state.max_iterations" queries spread out over a total of
  // kNumIndexSamples different geometry samples.
  static constexpr int kNumIndexSamples = 8;

  // To save time, we generate at most this many distinct targets per index.
  static constexpr int kMaxTargetsPerIndex = 100;

  S2CellIndex index;
  S2ClosestCellQuery query(&index);
  query.mutable_options()->set_max_results(max_results);
  const S1Angle radius =
      S2Testing::KmToAngle(absl::GetFlag(FLAGS_bm_radius_km));
  if (max_distance_fraction > 0) {
    query.mutable_options()->set_max_distance(max_distance_fraction * radius);
  }
  if (max_error_fraction > 0) {
    query.mutable_options()->set_max_error(max_error_fraction * radius);
  }
  query.mutable_options()->set_use_brute_force(
      absl::GetFlag(FLAGS_bm_use_brute_force));

  int delta = 0;  // Bresenham-type algorithm for geometry sampling.
  vector<unique_ptr<S2ClosestCellQuery::Target>> targets;
  vector<unique_ptr<MutableS2ShapeIndex>> target_shape_indexes;
  vector<S2ClosestCellQuery::Result> results;
  int64_t num_results = 0;
  int i_target = 0;
  for (auto _ : state) {
    delta -= kNumIndexSamples;
    if (delta < 0) {
      // Generate a new index and a new set of targets to go with it.
      // Reset the random number seed so that we use the same sequence of
      // indexed shapes no matter how many iterations are specified.
      state.PauseTiming();
      delta += state.max_iterations;
      index.Clear();
      S2Cap index_cap(s2random::Point(bitgen), radius);
      factory.AddCells(index_cap, num_index_cells, bitgen, &index);
      index.Build();
      query.ReInit();
      targets.clear();
      for (int i = 0; i < kMaxTargetsPerIndex; ++i) {
        // "target_cap" is chosen randomly such that the radius of
        // "target_cap" and the distance between the centers of "target_cap"
        // and "index_cap" are each some specified fraction of the index cap
        // radius.  (The index cap loosely bounds the indexed cells.)
        S1Angle target_dist =
            FractionToRadius(bitgen, center_separation_fraction);
        S2Cap target_cap(
            SampleBoundary(bitgen, S2Cap(index_cap.center(), target_dist)),
            FractionToRadius(bitgen, target_radius_fraction));

        if (target_type == TargetType::POINT) {
          S2Point v0;
          if (choose_target_from_index) {
            v0 = SampleCell(bitgen, index).ToPoint();
          } else {
            v0 = target_cap.center();
          }
          targets.push_back(make_unique<S2ClosestCellQuery::PointTarget>(v0));

        } else if (target_type == TargetType::EDGE) {
          S2Point v0, v1;
          if (choose_target_from_index) {
            v0 = SampleCell(bitgen, index).ToPoint();
            v1 = SampleCell(bitgen, index).ToPoint();
          } else {
            v0 = SampleBoundary(bitgen, target_cap);
            v1 = SampleBoundary(bitgen, target_cap);
          }
          targets.push_back(make_unique<S2ClosestCellQuery::EdgeTarget>(
              v0, v1));

        } else if (target_type == TargetType::CELL) {
          S2CellId cellid;
          if (choose_target_from_index) {
            cellid = SampleCell(bitgen, index);
          } else {
            cellid = S2CellId(target_cap.center()).parent(
                S2::kMaxDiag.GetClosestLevel(target_cap.GetRadius().radians()));
          }
          targets.push_back(make_unique<S2ClosestCellQuery::CellTarget>(
              S2Cell(cellid)));

        } else if (target_type == TargetType::CELL_UNION) {
          if (choose_target_from_index) {
            target_cap =
                S2Cap(SampleCell(bitgen, index).ToPoint(), target_cap.radius());
          }
          CapsCellIndexFactory target_factory(num_target_objects, 1.0);
          S2CellIndex target_index;
          factory.AddCells(target_cap, num_target_objects, bitgen,
                           &target_index);
          target_index.Build();
          vector<S2CellId> target_ids;
          for (S2CellIndex::CellIterator it(&target_index);
               !it.done(); it.Next()) {
              target_ids.push_back(it.cell_id());
          }
          targets.push_back(make_unique<S2ClosestCellQuery::CellUnionTarget>(
              S2CellUnion(std::move(target_ids))));

        } else {
          ABSL_DCHECK(target_type == TargetType::SHAPE_INDEX);
          if (choose_target_from_index) {
            target_cap =
                S2Cap(SampleCell(bitgen, index).ToPoint(), target_cap.radius());
          }
          auto target_index = make_unique<MutableS2ShapeIndex>();
          FractalLoopShapeIndexFactory(bitgen).AddEdges(
              target_cap, num_target_objects, &*target_index);
          target_index->ForceBuild();
          targets.push_back(make_unique<S2ClosestCellQuery::ShapeIndexTarget>(
              &*target_index));
          target_shape_indexes.push_back(std::move(target_index));
        }
      }
      state.ResumeTiming();
    }
    query.FindClosestCells(targets[i_target].get(), &results);
    num_results += results.size();
    if (++i_target == targets.size()) i_target = 0;
  }
  state.SetLabel(absl::StrFormat(
      "results:%.3f", static_cast<double>(num_results) / state.iterations()));
}

////////////////////////////////////////////////////////////////////////
// See overview of benchmarks above ("Description of the benchmarks").

// The maximum number of cells to use in the benchmarks.
// Some benchmarks are run only with this number of cells.
static constexpr int kMaxCells = 3 * 16384;

// Define shorthand names to make the benchmark names less verbose.
using PointCloud = PointCloudCellIndexFactory;

// Test finding the exact minimum distance for points in the general vicinity
// of the indexed cells (chosen randomly within a region 4x larger).
template <class Factory>
static void BM_FindClosest(benchmark::State& state) {
  const int num_cells = state.range(0);
  BenchmarkFindClosest(state, Factory(), num_cells, 1, -1, -1,
                       TargetType::POINT, 0, false, 0.0, -2.0);
}
BENCHMARK_TEMPLATE(BM_FindClosest, PointCloud)
// ->Arg(6)->Arg(9)->Arg(12)->Arg(15)->Arg(20)  // Tuning
->Arg(3*4)->Arg(3*16)->Arg(3*256)->Arg(3*4096)->Arg(kMaxCells);

// Test finding the exact minimum distance when that distance is zero.
template <class Factory>
static void BM_FindClosestWhenZero(benchmark::State& state) {
  const int num_cells = state.range(0);
  BenchmarkFindClosest(state, Factory(), num_cells, 1, -1, -1,
                       TargetType::POINT, 0, true, 0.0, -2.0);
}
BENCHMARK_TEMPLATE(BM_FindClosestWhenZero, PointCloud)
->Arg(3 * 256)->Arg(kMaxCells);

// Test checking whether the minimum distance is above or below a given
// threshold, where the threshold is chosen so that the query disc contains
// approximately 2**log2_num_results cells on average.
template <class Factory>
static void BM_IsDistanceLess(benchmark::State& state) {
  const int num_cells = state.range(0);
  const int log2_num_results = state.range(1);
  // The target is located within the index cap about 50% of the time, because
  // the distance to the cap center is chosen uniformly at random up to twice
  // the index cap radius.
  double max_dist_fraction = sqrt(2 * pow(2.0, log2_num_results) / num_cells);
  BenchmarkFindClosest(state, Factory(), num_cells, 1, max_dist_fraction, 1.0,
                       TargetType::POINT, 0, false, 0.0, -2.0);
}
// This operation gets faster when log2_num_results is small *or* large.
BENCHMARK_TEMPLATE(BM_IsDistanceLess, PointCloud)
// ->ArgPair(6, 0)->ArgPair(9, 0)->ArgPair(12, 0)  // Tuning
// ->ArgPair(15, 0)->ArgPair(20, 0)                // Tuning
->ArgPair(3*4, -3)->ArgPair(3*16, -3)->ArgPair(3*256, -3)->ArgPair(3*4096, -3)
->ArgPair(kMaxCells, -6)->ArgPair(kMaxCells, -3)
->ArgPair(kMaxCells, 0)->ArgPair(kMaxCells, 3);

// Like the test above, but measures the speed of finding all results within a
// given distance threshold.
template <class Factory>
static void BM_FindCellsNearby(benchmark::State& state) {
  const int num_cells = state.range(0);
  const int log2_num_results = state.range(1);
  double max_dist_fraction = sqrt(2 * pow(2.0, log2_num_results) / num_cells);
  int max_results = S2ClosestCellQuery::Options::kMaxMaxResults;
  BenchmarkFindClosest(state, Factory(), num_cells, max_results,
                       max_dist_fraction, 10.0, TargetType::POINT, 0, false,
                       0.0, -2.0);
}
// This operation gets faster when log2_num_results is small *or* large.
BENCHMARK_TEMPLATE(BM_FindCellsNearby, PointCloud)
->ArgPair(3*4, -3)->ArgPair(3*16, -3)->ArgPair(3*256, -3)->ArgPair(3*4096, -3)
->ArgPair(kMaxCells, -6)->ArgPair(kMaxCells, -3)
->ArgPair(kMaxCells, 0)->ArgPair(kMaxCells, 3);

// Repeat the benchmarks for edge targets rather than point targets.
// The edge length is chosen randomly up to about 10% of the index diameter.

template <class Factory>
static void BM_FindClosestToLongEdge(benchmark::State& state) {
  const int num_cells = state.range(0);
  BenchmarkFindClosest(state, Factory(), num_cells, 1, -1, -1, TargetType::EDGE,
                       0, false, -0.1, -2.0);
}
BENCHMARK_TEMPLATE(BM_FindClosestToLongEdge, PointCloud)
// ->Arg(6)->Arg(9)->Arg(12)->Arg(15)->Arg(20)  // Tuning
->Arg(3*4)->Arg(3*256)->Arg(kMaxCells);

template <class Factory>
static void BM_FindClosestWhenZeroToLongEdge(benchmark::State& state) {
  const int num_cells = state.range(0);
  BenchmarkFindClosest(state, Factory(), num_cells, 1, -1, -1, TargetType::EDGE,
                       0, true, -0.1, -2.0);
}
BENCHMARK_TEMPLATE(BM_FindClosestWhenZeroToLongEdge, PointCloud)
->Arg(3 * 256)->Arg(kMaxCells);

template <class Factory>
static void BM_IsDistanceLessToLongEdge(benchmark::State& state) {
  const int num_cells = state.range(0);
  const int log2_num_results = state.range(1);
  double max_dist_fraction = sqrt(2 * pow(2.0, log2_num_results) / num_cells);
  BenchmarkFindClosest(state, Factory(), num_cells, 1, max_dist_fraction, 1.0,
                       TargetType::EDGE, 0, false, -0.1, -2.0);
}
BENCHMARK_TEMPLATE(BM_IsDistanceLessToLongEdge, PointCloud)
// ->ArgPair(6, 0)->ArgPair(9, 0)->ArgPair(12, 0)  // Tuning
// ->ArgPair(15, 0)->ArgPair(20, 0)                // Tuning
->ArgPair(3*4, -3)->ArgPair(3*256, -3)
->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

template <class Factory>
static void BM_FindCellsNearbyLongEdge(benchmark::State& state) {
  const int num_cells = state.range(0);
  const int log2_num_results = state.range(1);
  double max_dist_fraction = sqrt(2 * pow(2.0, log2_num_results) / num_cells);
  int max_results = S2ClosestCellQuery::Options::kMaxMaxResults;
  BenchmarkFindClosest(state, Factory(), num_cells, max_results,
                       max_dist_fraction, 10.0, TargetType::EDGE, 0, false,
                       -0.1, -2.0);
}
BENCHMARK_TEMPLATE(BM_FindCellsNearbyLongEdge, PointCloud)
->ArgPair(3*4, -3)->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3);

// Repeat the benchmarks for S2Cell targets.  The S2Cell diagonal is chosen
// randomly up to about 1% of the index radius.

template <class Factory>
static void BM_FindClosestToSmallCell(benchmark::State& state) {
  const int num_cells = state.range(0);
  BenchmarkFindClosest(state, Factory(), num_cells, 1, -1, -1, TargetType::CELL,
                       0, false, -0.01, -2.0);
}
BENCHMARK_TEMPLATE(BM_FindClosestToSmallCell, PointCloud)
// ->Arg(6)->Arg(9)->Arg(12)->Arg(15)->Arg(20)  // Tuning
->Arg(3*4)->Arg(3*256)->Arg(kMaxCells);

template <class Factory>
static void BM_FindClosestWhenZeroToSmallCell(benchmark::State& state) {
  const int num_cells = state.range(0);
  BenchmarkFindClosest(state, Factory(), num_cells, 1, -1, -1, TargetType::CELL,
                       0, true, -0.01, -2.0);
}
BENCHMARK_TEMPLATE(BM_FindClosestWhenZeroToSmallCell, PointCloud)
->Arg(3 * 256)->Arg(kMaxCells);

template <class Factory>
static void BM_IsDistanceLessToSmallCell(benchmark::State& state) {
  const int num_cells = state.range(0);
  const int log2_num_results = state.range(1);
  double max_dist_fraction = sqrt(2 * pow(2.0, log2_num_results) / num_cells);
  BenchmarkFindClosest(state, Factory(), num_cells, 1, max_dist_fraction, 1.0,
                       TargetType::CELL, 0, false, -0.01, -2.0);
}
BENCHMARK_TEMPLATE(BM_IsDistanceLessToSmallCell, PointCloud)
// ->ArgPair(6, 0)->ArgPair(9, 0)->ArgPair(12, 0)  // Tuning
// ->ArgPair(15, 0)->ArgPair(20, 0)                // Tuning
->ArgPair(3*4, -3)->ArgPair(3*256, -3)
->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

template <class Factory>
static void BM_FindCellsNearbySmallCell(benchmark::State& state) {
  const int num_cells = state.range(0);
  const int log2_num_results = state.range(1);
  double max_dist_fraction = sqrt(2 * pow(2.0, log2_num_results) / num_cells);
  int max_results = S2ClosestCellQuery::Options::kMaxMaxResults;
  BenchmarkFindClosest(state, Factory(), num_cells, max_results,
                       max_dist_fraction, 10.0, TargetType::CELL, 0, false,
                       -0.01, -2.0);
}
BENCHMARK_TEMPLATE(BM_FindCellsNearbySmallCell, PointCloud)
->ArgPair(3*4, -3)->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3);


// Repeat the benchmarks for S2CellUnion targets.  Each S2CellUnion consists
// of up to 8 cells ("coarse") or up to 100 cells ("fine") that cover an S2Cap
// whose radius is chosen randomly up to 1% of the index radius.

template <class Factory>
static void BM_FindClosestToSmallCellUnion(benchmark::State& state,
                                           int num_target_cells) {
  const int num_cells = state.range(0);
  BenchmarkFindClosest(state, Factory(), num_cells, 1, -1, -1,
                       TargetType::CELL_UNION, num_target_cells, false, -0.01,
                       -2.0);
}

template <class Factory>
static void BM_FindClosestToSmallCoarseCellUnion(benchmark::State& state) {
  BM_FindClosestToSmallCellUnion<Factory>(state, 8);
}
BENCHMARK_TEMPLATE(BM_FindClosestToSmallCoarseCellUnion, PointCloud)
// ->Arg(6)->Arg(9)->Arg(12)->Arg(15)->Arg(20)  // Tuning
->Arg(3*4)->Arg(3*256)->Arg(kMaxCells);

template <class Factory>
static void BM_FindClosestToSmallFineCellUnion(benchmark::State& state) {
  BM_FindClosestToSmallCellUnion<Factory>(state, 100);
}
BENCHMARK_TEMPLATE(BM_FindClosestToSmallFineCellUnion, PointCloud)
// ->Arg(6)->Arg(9)->Arg(12)->Arg(15)->Arg(20)  // Tuning
->Arg(kMaxCells);

template <class Factory>
static void BM_FindClosestWhenZeroToSmallCellUnion(benchmark::State& state,
                                                   int num_target_cells) {
  const int num_cells = state.range(0);
  BenchmarkFindClosest(state, Factory(), num_cells, 1, -1, -1,
                       TargetType::CELL_UNION, num_target_cells, true, -0.01,
                       -2.0);
}

template <class Factory>
static void BM_FindClosestWhenZeroToSmallCoarseCellUnion(
    benchmark::State& state) {
  BM_FindClosestWhenZeroToSmallCellUnion<Factory>(state, 8);
}
BENCHMARK_TEMPLATE(BM_FindClosestWhenZeroToSmallCoarseCellUnion, PointCloud)
->Arg(3*256)->Arg(kMaxCells);

template <class Factory>
static void BM_FindClosestWhenZeroToSmallFineCellUnion(
    benchmark::State& state) {
  BM_FindClosestWhenZeroToSmallCellUnion<Factory>(state, 100);
}
BENCHMARK_TEMPLATE(BM_FindClosestWhenZeroToSmallFineCellUnion, PointCloud)
->Arg(kMaxCells);

template <class Factory>
static void BM_IsDistanceLessToSmallCellUnion(benchmark::State& state,
                                              int num_cells,
                                              int num_target_cells,
                                              int log2_num_results) {
  double max_dist_fraction = sqrt(2 * pow(2.0, log2_num_results) / num_cells);
  BenchmarkFindClosest(state, Factory(), num_cells, 1, max_dist_fraction, 1.0,
                       TargetType::CELL_UNION, num_target_cells, false, -0.01,
                       -2.0);
}

template <class Factory>
static void BM_IsDistanceLessToSmallCoarseCellUnion(benchmark::State& state) {
  const int num_cells = state.range(0);
  const int log2_num_results = state.range(1);
  BM_IsDistanceLessToSmallCellUnion<Factory>(state, num_cells, 8,
                                             log2_num_results);
}
BENCHMARK_TEMPLATE(BM_IsDistanceLessToSmallCoarseCellUnion, PointCloud)
// ->ArgPair(6, 0)->ArgPair(9, 0)->ArgPair(12, 0)  // Tuning
// ->ArgPair(15, 0)->ArgPair(20, 0)                // Tuning
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

template <class Factory>
static void BM_IsDistanceLessToSmallFineCellUnion(benchmark::State& state) {
  const int num_cells = state.range(0);
  const int log2_num_results = state.range(1);
  BM_IsDistanceLessToSmallCellUnion<Factory>(state, num_cells, 100,
                                             log2_num_results);
}
BENCHMARK_TEMPLATE(BM_IsDistanceLessToSmallFineCellUnion, PointCloud)
// ->ArgPair(6, 0)->ArgPair(9, 0)->ArgPair(12, 0)  // Tuning
// ->ArgPair(15, 0)->ArgPair(20, 0)                // Tuning
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

template <class Factory>
static void BM_FindCellsNearbySmallCellUnion(benchmark::State& state,
                                             int num_target_cells) {
  const int num_cells = state.range(0);
  const int log2_num_results = state.range(1);
  double max_dist_fraction = sqrt(2 * pow(2.0, log2_num_results) / num_cells);
  int max_results = S2ClosestCellQuery::Options::kMaxMaxResults;
  BenchmarkFindClosest(state, Factory(), num_cells, max_results,
                       max_dist_fraction, 10.0, TargetType::CELL_UNION,
                       num_target_cells, false, -0.01, -2.0);
}

template <class Factory>
static void BM_FindCellsNearbySmallCoarseCellUnion(benchmark::State& state) {
  BM_FindCellsNearbySmallCellUnion<Factory>(state, 8);
}
BENCHMARK_TEMPLATE(BM_FindCellsNearbySmallCoarseCellUnion, PointCloud)
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

template <class Factory>
static void BM_FindCellsNearbySmallFineCellUnion(benchmark::State& state) {
  BM_FindCellsNearbySmallCellUnion<Factory>(state, 100);
}
BENCHMARK_TEMPLATE(BM_FindCellsNearbySmallFineCellUnion, PointCloud)
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);


// Repeat the benchmarks for S2ShapeIndex targets.  Each ShapeIndex is a
// fractal loop with either 12 edges ("coarse") or 768 edges ("fine") that
// approximately fills an S2Cap whose radius is chosen randomly up to 1% of
// the index radius.

template <class Factory>
static void BM_FindClosestToSmallShapeIndex(benchmark::State& state,
                                            int num_target_edges) {
  const int num_cells = state.range(0);
  BenchmarkFindClosest(state, Factory(), num_cells, 1, -1, -1,
                       TargetType::SHAPE_INDEX, num_target_edges, false, -0.01,
                       -2.0);
}

template <class Factory>
static void BM_FindClosestToSmallCoarseShapeIndex(benchmark::State& state) {
  BM_FindClosestToSmallShapeIndex<Factory>(state, 3 * 4);
}
BENCHMARK_TEMPLATE(BM_FindClosestToSmallCoarseShapeIndex, PointCloud)
// ->Arg(6)->Arg(9)->Arg(12)->Arg(15)->Arg(20)  // Tuning
->Arg(3*4)->Arg(3*256)->Arg(kMaxCells);

template <class Factory>
static void BM_FindClosestToSmallFineShapeIndex(benchmark::State& state) {
  BM_FindClosestToSmallShapeIndex<Factory>(state, 3 * 256);
}
BENCHMARK_TEMPLATE(BM_FindClosestToSmallFineShapeIndex, PointCloud)
// ->Arg(6)->Arg(9)->Arg(12)->Arg(15)->Arg(20)  // Tuning
->Arg(kMaxCells);

template <class Factory>
static void BM_FindClosestWhenZeroToSmallShapeIndex(benchmark::State& state,
                                                    int num_target_edges) {
  const int num_cells = state.range(0);
  BenchmarkFindClosest(state, Factory(), num_cells, 1, -1, -1,
                       TargetType::SHAPE_INDEX, num_target_edges, true, -0.01,
                       -2.0);
}

template <class Factory>
static void BM_FindClosestWhenZeroToSmallCoarseShapeIndex(
    benchmark::State& state) {
  BM_FindClosestWhenZeroToSmallShapeIndex<Factory>(state, 3 * 4);
}
BENCHMARK_TEMPLATE(BM_FindClosestWhenZeroToSmallCoarseShapeIndex, PointCloud)
->Arg(3*4)->Arg(3*256)->Arg(kMaxCells);

template <class Factory>
static void BM_FindClosestWhenZeroToSmallFineShapeIndex(
    benchmark::State& state) {
  BM_FindClosestWhenZeroToSmallShapeIndex<Factory>(state, 3 * 256);
}
BENCHMARK_TEMPLATE(BM_FindClosestWhenZeroToSmallFineShapeIndex, PointCloud)
->Arg(kMaxCells);

template <class Factory>
static void BM_IsDistanceLessToSmallShapeIndex(benchmark::State& state,
                                               int num_target_edges) {
  const int num_cells = state.range(0);
  const int log2_num_results = state.range(1);
  double max_dist_fraction = sqrt(2 * pow(2.0, log2_num_results) / num_cells);
  BenchmarkFindClosest(state, Factory(), num_cells, 1, max_dist_fraction, 1.0,
                       TargetType::SHAPE_INDEX, num_target_edges, false, -0.01,
                       2.0);
}

template <class Factory>
static void BM_IsDistanceLessToSmallCoarseShapeIndex(benchmark::State& state) {
  BM_IsDistanceLessToSmallShapeIndex<Factory>(state, 3 * 4);
}
BENCHMARK_TEMPLATE(BM_IsDistanceLessToSmallCoarseShapeIndex, PointCloud)
// ->ArgPair(6, 0)->ArgPair(9, 0)->ArgPair(12, 0)  // Tuning
// ->ArgPair(15, 0)->ArgPair(20, 0)                // Tuning
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

template <class Factory>
static void BM_IsDistanceLessToSmallFineShapeIndex(benchmark::State& state) {
  BM_IsDistanceLessToSmallShapeIndex<Factory>(state, 3 * 256);
}
BENCHMARK_TEMPLATE(BM_IsDistanceLessToSmallFineShapeIndex, PointCloud)
// ->ArgPair(6, 0)->ArgPair(9, 0)->ArgPair(12, 0)  // Tuning
// ->ArgPair(15, 0)->ArgPair(20, 0)                // Tuning
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

template <class Factory>
static void BM_FindCellsNearbySmallShapeIndex(benchmark::State& state,
                                              int num_target_edges) {
  const int num_cells = state.range(0);
  const int log2_num_results = state.range(1);
  double max_dist_fraction = sqrt(2 * pow(2.0, log2_num_results) / num_cells);
  int max_results = S2ClosestCellQuery::Options::kMaxMaxResults;
  BenchmarkFindClosest(state, Factory(), num_cells, max_results,
                       max_dist_fraction, 10.0, TargetType::SHAPE_INDEX,
                       num_target_edges, false, -0.01, -2.0);
}

template <class Factory>
static void BM_FindCellsNearbySmallCoarseShapeIndex(benchmark::State& state) {
  BM_FindCellsNearbySmallShapeIndex<Factory>(state, 3 * 4);
}
BENCHMARK_TEMPLATE(BM_FindCellsNearbySmallCoarseShapeIndex, PointCloud)
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

template <class Factory>
static void BM_FindCellsNearbySmallFineShapeIndex(benchmark::State& state) {
  BM_FindCellsNearbySmallShapeIndex<Factory>(state, 3 * 256);
}
BENCHMARK_TEMPLATE(BM_FindCellsNearbySmallFineShapeIndex, PointCloud)
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);


// Now repeat all the benchmarks for S2Cap coverings rather than point clouds.
// We group the benchmarks together by the type of geometry so that it's
// easier to see what effect the various options have (max_distance, etc).

// There are lots of combinations that we could benchmark.  Here we assume
// that the index consists of relatively coarse coverings (up to 8 cells), and
// that the coverings generally don't overlap (they fill about 10% of the
// indexed region).
class CoarseCaps : public CapsCellIndexFactory {
 public:
  CoarseCaps() : CapsCellIndexFactory(8 /*max_cells*/, 0.1 /*density*/) {}
};

BENCHMARK_TEMPLATE(BM_FindClosest, CoarseCaps)
// ->Arg(8)->Arg(16)->Arg(24)->Arg(32)  // Tuning
->Arg(3*16)->Arg(3*256)->Arg(3*4096)->Arg(kMaxCells);

BENCHMARK_TEMPLATE(BM_FindClosestWhenZero, CoarseCaps)
->Arg(3 * 256)->Arg(kMaxCells);

BENCHMARK_TEMPLATE(BM_IsDistanceLess, CoarseCaps)
// ->ArgPair(8, 0)->ArgPair(16, 0)->ArgPair(24, 0)->ArgPair(32, 0)  // Tuning
->ArgPair(3*16, -3)->ArgPair(3*256, -3)->ArgPair(3*4096, -3)
->ArgPair(kMaxCells, -6)->ArgPair(kMaxCells, -3)
->ArgPair(kMaxCells, 0)->ArgPair(kMaxCells, 3);

BENCHMARK_TEMPLATE(BM_FindCellsNearby, CoarseCaps)
->ArgPair(3*16, -3)->ArgPair(3*256, -3)->ArgPair(3*4096, -3)
->ArgPair(kMaxCells, -6)->ArgPair(kMaxCells, -3)
->ArgPair(kMaxCells, 0)->ArgPair(kMaxCells, 3);

BENCHMARK_TEMPLATE(BM_FindClosestToLongEdge, CoarseCaps)
// ->Arg(8)->Arg(16)->Arg(24)->Arg(32)  // Tuning
->Arg(3*256)->Arg(kMaxCells);

BENCHMARK_TEMPLATE(BM_FindClosestWhenZeroToLongEdge, CoarseCaps)
->Arg(3 * 256)->Arg(kMaxCells);

BENCHMARK_TEMPLATE(BM_IsDistanceLessToLongEdge, CoarseCaps)
// ->ArgPair(8, 0)->ArgPair(16, 0)->ArgPair(24, 0)->ArgPair(32, 0)  // Tuning
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

BENCHMARK_TEMPLATE(BM_FindCellsNearbyLongEdge, CoarseCaps)
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3);

BENCHMARK_TEMPLATE(BM_FindClosestToSmallCell, CoarseCaps)
// ->Arg(8)->Arg(16)->Arg(24)->Arg(32)  // Tuning
->Arg(3*256)->Arg(kMaxCells);

BENCHMARK_TEMPLATE(BM_FindClosestWhenZeroToSmallCell, CoarseCaps)
->Arg(3 * 256)->Arg(kMaxCells);

BENCHMARK_TEMPLATE(BM_IsDistanceLessToSmallCell, CoarseCaps)
// ->ArgPair(8, 0)->ArgPair(16, 0)->ArgPair(24, 0)->ArgPair(32, 0)  // Tuning
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

BENCHMARK_TEMPLATE(BM_FindCellsNearbySmallCell, CoarseCaps)
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3);

BENCHMARK_TEMPLATE(BM_FindClosestToSmallCoarseCellUnion, CoarseCaps)
// ->Arg(8)->Arg(16)->Arg(24)->Arg(32)  // Tuning
->Arg(3*256)->Arg(kMaxCells);

BENCHMARK_TEMPLATE(BM_FindClosestToSmallFineCellUnion, CoarseCaps)
// ->Arg(8)->Arg(16)->Arg(24)->Arg(32)  // Tuning
->Arg(kMaxCells);

BENCHMARK_TEMPLATE(BM_FindClosestWhenZeroToSmallCoarseCellUnion, CoarseCaps)
->Arg(3*256)->Arg(kMaxCells);

BENCHMARK_TEMPLATE(BM_FindClosestWhenZeroToSmallFineCellUnion, CoarseCaps)
->Arg(kMaxCells);

BENCHMARK_TEMPLATE(BM_IsDistanceLessToSmallCoarseCellUnion, CoarseCaps)
// ->ArgPair(8, 0)->ArgPair(16, 0)->ArgPair(24, 0)->ArgPair(32, 0)  // Tuning
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

BENCHMARK_TEMPLATE(BM_IsDistanceLessToSmallFineCellUnion, CoarseCaps)
// ->ArgPair(8, 0)->ArgPair(16, 0)->ArgPair(24, 0)->ArgPair(32, 0)  // Tuning
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

BENCHMARK_TEMPLATE(BM_FindCellsNearbySmallCoarseCellUnion, CoarseCaps)
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

BENCHMARK_TEMPLATE(BM_FindCellsNearbySmallFineCellUnion, CoarseCaps)
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

BENCHMARK_TEMPLATE(BM_FindClosestToSmallCoarseShapeIndex, CoarseCaps)
// ->Arg(8)->Arg(16)->Arg(24)->Arg(32)  // Tuning
->Arg(3*256)->Arg(kMaxCells);

BENCHMARK_TEMPLATE(BM_FindClosestToSmallFineShapeIndex, CoarseCaps)
// ->Arg(8)->Arg(16)->Arg(24)->Arg(32)  // Tuning
->Arg(kMaxCells);

BENCHMARK_TEMPLATE(BM_FindClosestWhenZeroToSmallCoarseShapeIndex, CoarseCaps)
->Arg(3*256)->Arg(kMaxCells);

BENCHMARK_TEMPLATE(BM_FindClosestWhenZeroToSmallFineShapeIndex, CoarseCaps)
->Arg(kMaxCells);

BENCHMARK_TEMPLATE(BM_IsDistanceLessToSmallCoarseShapeIndex, CoarseCaps)
// ->ArgPair(8, 0)->ArgPair(16, 0)->ArgPair(24, 0)->ArgPair(32, 0)  // Tuning
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

BENCHMARK_TEMPLATE(BM_IsDistanceLessToSmallFineShapeIndex, CoarseCaps)
// ->ArgPair(8, 0)->ArgPair(16, 0)->ArgPair(24, 0)->ArgPair(32, 0)  // Tuning
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

BENCHMARK_TEMPLATE(BM_FindCellsNearbySmallCoarseShapeIndex, CoarseCaps)
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

BENCHMARK_TEMPLATE(BM_FindCellsNearbySmallFineShapeIndex, CoarseCaps)
->ArgPair(3*256, -3)->ArgPair(kMaxCells, -3)->ArgPair(kMaxCells, 6);

}  // namespace
