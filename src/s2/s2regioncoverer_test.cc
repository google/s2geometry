// Copyright 2005 Google Inc. All Rights Reserved.
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

#include "s2/s2regioncoverer.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtest/gtest.h>

#include "s2/base/stringprintf.h"
#include "s2/base/strtoint.h"
#include "s2/s1angle.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cellid.h"
#include "s2/s2cellunion.h"
#include "s2/s2latlng.h"
#include "s2/s2region.h"
#include "s2/s2testing.h"
#include "s2/third_party/absl/base/integral_types.h"
#include "s2/third_party/absl/strings/str_split.h"

using std::max;
using std::min;
using std::priority_queue;
using std::unordered_map;
using std::vector;

DEFINE_string(max_cells, "4,8",
              "Comma-separated list of values to use for 'max_cells'");

DEFINE_int32(iters, google::DEBUG_MODE ? 1000 : 100000,
             "Number of random caps to try for each max_cells value");

TEST(S2RegionCoverer, RandomCells) {
  S2RegionCoverer::Options options;
  options.set_max_cells(1);
  S2RegionCoverer coverer(options);

  // Test random cell ids at all levels.
  for (int i = 0; i < 10000; ++i) {
    S2CellId id = S2Testing::GetRandomCellId();
    SCOPED_TRACE(StringPrintf("Iteration %d, cell ID token %s",
                              i, id.ToToken().c_str()));
    vector<S2CellId> covering = coverer.GetCovering(S2Cell(id)).Release();
    EXPECT_EQ(1, covering.size());
    EXPECT_EQ(id, covering[0]);
  }
}

static void CheckCovering(const S2RegionCoverer::Options& options,
                          const S2Region& region,
                          const vector<S2CellId>& covering,
                          bool interior) {
  // Keep track of how many cells have the same options.min_level() ancestor.
  unordered_map<S2CellId, int, S2CellIdHash> min_level_cells;
  for (S2CellId cell_id : covering) {
    int level = cell_id.level();
    EXPECT_GE(level, options.min_level());
    EXPECT_LE(level, options.max_level());
    EXPECT_EQ((level - options.min_level()) % options.level_mod(), 0);
    min_level_cells[cell_id.parent(options.min_level())] += 1;
  }
  if (covering.size() > options.max_cells()) {
    // If the covering has more than the requested number of cells, then check
    // that the cell count cannot be reduced by using the parent of some cell.
    for (unordered_map<S2CellId, int, S2CellIdHash>::const_iterator i =
             min_level_cells.begin();
         i != min_level_cells.end(); ++i) {
      EXPECT_EQ(i->second, 1);
    }
  }
  if (interior) {
    for (S2CellId cell_id : covering) {
      EXPECT_TRUE(region.Contains(S2Cell(cell_id)));
    }
  } else {
    S2CellUnion cell_union(covering);
    S2Testing::CheckCovering(region, cell_union, true);
  }
}

TEST(S2RegionCoverer, RandomCaps) {
  static const int kMaxLevel = S2CellId::kMaxLevel;
  S2RegionCoverer::Options options;
  for (int i = 0; i < 1000; ++i) {
    do {
      options.set_min_level(S2Testing::rnd.Uniform(kMaxLevel + 1));
      options.set_max_level(S2Testing::rnd.Uniform(kMaxLevel + 1));
    } while (options.min_level() > options.max_level());
    options.set_max_cells(S2Testing::rnd.Skewed(10));
    options.set_level_mod(1 + S2Testing::rnd.Uniform(3));
    double max_area =  min(4 * M_PI, (3 * options.max_cells() + 1) *
                           S2Cell::AverageArea(options.min_level()));
    S2Cap cap = S2Testing::GetRandomCap(0.1 * S2Cell::AverageArea(kMaxLevel),
                                        max_area);
    S2RegionCoverer coverer(options);
    vector<S2CellId> covering, interior;
    coverer.GetCovering(cap, &covering);
    CheckCovering(options, cap, covering, false);
    coverer.GetInteriorCovering(cap, &interior);
    CheckCovering(options, cap, interior, true);

    // Check that GetCovering is deterministic.
    vector<S2CellId> covering2;
    coverer.GetCovering(cap, &covering2);
    EXPECT_EQ(covering, covering2);

    // Also check S2CellUnion::Denormalize().  The denormalized covering
    // may still be different and smaller than "covering" because
    // S2RegionCoverer does not guarantee that it will not output all four
    // children of the same parent.
    S2CellUnion cells(covering);
    vector<S2CellId> denormalized;
    cells.Denormalize(options.min_level(), options.level_mod(), &denormalized);
    CheckCovering(options, cap, denormalized, false);
  }
}

TEST(S2RegionCoverer, SimpleCoverings) {
  static const int kMaxLevel = S2CellId::kMaxLevel;
  S2RegionCoverer::Options options;
  options.set_max_cells(kint32max);
  for (int i = 0; i < 1000; ++i) {
    int level = S2Testing::rnd.Uniform(kMaxLevel + 1);
    options.set_min_level(level);
    options.set_max_level(level);
    double max_area =  min(4 * M_PI, 1000 * S2Cell::AverageArea(level));
    S2Cap cap = S2Testing::GetRandomCap(0.1 * S2Cell::AverageArea(kMaxLevel),
                                        max_area);
    vector<S2CellId> covering;
    S2RegionCoverer::GetSimpleCovering(cap, cap.center(), level, &covering);
    CheckCovering(options, cap, covering, false);
  }
}

// We keep a priority queue of the caps that had the worst approximation
// ratios so that we can print them at the end.
struct WorstCap {
  double ratio;
  S2Cap cap;
  int num_cells;
  bool operator<(const WorstCap& o) const { return ratio > o.ratio; }
  WorstCap(double r, const S2Cap& c, int n) : ratio(r), cap(c), num_cells(n) {}
};

static void TestAccuracy(int max_cells) {
  SCOPED_TRACE(StringPrintf("%d cells", max_cells));

  static const int kNumMethods = 1;
  // This code is designed to evaluate several approximation algorithms and
  // figure out which one works better.  The way to do this is to hack the
  // S2RegionCoverer interface to add a global variable to control which
  // algorithm (or variant of an algorithm) is selected, and then assign to
  // this variable in the "method" loop below.  The code below will then
  // collect statistics on all methods, including how often each one wins in
  // terms of cell count and approximation area.

  S2RegionCoverer coverer;
  coverer.mutable_options()->set_max_cells(max_cells);

  double ratio_total[kNumMethods] = {0};
  double min_ratio[kNumMethods];  // initialized in loop below
  double max_ratio[kNumMethods] = {0};
  vector<double> ratios[kNumMethods];
  int cell_total[kNumMethods] = {0};
  int area_winner_tally[kNumMethods] = {0};
  int cell_winner_tally[kNumMethods] = {0};
  static const int kMaxWorstCaps = 10;
  priority_queue<WorstCap> worst_caps[kNumMethods];

  for (int method = 0; method < kNumMethods; ++method) {
    min_ratio[method] = 1e20;
  }
  for (int i = 0; i < FLAGS_iters; ++i) {
    // Choose the log of the cap area to be uniformly distributed over
    // the allowable range.  Don't try to approximate regions that are so
    // small they can't use the given maximum number of cells efficiently.
    const double min_cap_area = S2Cell::AverageArea(S2CellId::kMaxLevel)
                                * max_cells * max_cells;
    // Coverings for huge caps are not interesting, so limit the max area too.
    S2Cap cap = S2Testing::GetRandomCap(min_cap_area, 0.1 * M_PI);
    double cap_area = cap.GetArea();

    double min_area = 1e30;
    int min_cells = 1 << 30;
    double area[kNumMethods];
    int cells[kNumMethods];
    for (int method = 0; method < kNumMethods; ++method) {
      // If you want to play with different methods, do this:
      // S2RegionCoverer::method_number = method;

      vector<S2CellId> covering;
      coverer.GetCovering(cap, &covering);

      double union_area = 0;
      for (S2CellId cell_id : covering) {
        union_area += S2Cell(cell_id).ExactArea();
      }
      cells[method] = covering.size();
      min_cells = min(cells[method], min_cells);
      area[method] = union_area;
      min_area = min(area[method], min_area);
      cell_total[method] += cells[method];
      double ratio = area[method] / cap_area;
      ratio_total[method] += ratio;
      min_ratio[method] = min(ratio, min_ratio[method]);
      max_ratio[method] = max(ratio, max_ratio[method]);
      ratios[method].push_back(ratio);
      if (worst_caps[method].size() < kMaxWorstCaps) {
        worst_caps[method].push(WorstCap(ratio, cap, cells[method]));
      } else if (ratio > worst_caps[method].top().ratio) {
        worst_caps[method].pop();
        worst_caps[method].push(WorstCap(ratio, cap, cells[method]));
      }
    }
    for (int method = 0; method < kNumMethods; ++method) {
      if (area[method] == min_area) ++area_winner_tally[method];
      if (cells[method] == min_cells) ++cell_winner_tally[method];
    }
  }
  for (int method = 0; method < kNumMethods; ++method) {
    printf("\nMax cells %d, method %d:\n", max_cells, method);
    printf("  Average cells: %.4f\n", cell_total[method] /
           static_cast<double>(FLAGS_iters));
    printf("  Average area ratio: %.4f\n", ratio_total[method] / FLAGS_iters);
    vector<double>& mratios = ratios[method];
    std::sort(mratios.begin(), mratios.end());
    printf("  Median ratio: %.4f\n", mratios[mratios.size() / 2]);
    printf("  Max ratio: %.4f\n", max_ratio[method]);
    printf("  Min ratio: %.4f\n", min_ratio[method]);
    if (kNumMethods > 1) {
      printf("  Cell winner probability: %.4f\n",
             cell_winner_tally[method] / static_cast<double>(FLAGS_iters));
      printf("  Area winner probability: %.4f\n",
             area_winner_tally[method] / static_cast<double>(FLAGS_iters));
    }
    printf("  Caps with the worst approximation ratios:\n");
    for (; !worst_caps[method].empty(); worst_caps[method].pop()) {
      const WorstCap& w = worst_caps[method].top();
      S2LatLng ll(w.cap.center());
      printf("    Ratio %.4f, Cells %d, "
             "Center (%.8f, %.8f), Km %.6f\n",
             w.ratio, w.num_cells,
             ll.lat().degrees(), ll.lng().degrees(),
             w.cap.GetRadius().radians() * 6367.0);
    }
  }
}

TEST(S2RegionCoverer, Accuracy) {
  vector<string> max_cells =
      absl::StrSplit(FLAGS_max_cells, ',', absl::SkipEmpty());
  for (int i = 0; i < max_cells.size(); ++i) {
    TestAccuracy(atoi32(max_cells[i].c_str()));
  }
}

TEST(S2RegionCoverer, InteriorCovering) {
  // We construct the region the following way. Start with S2 cell of level l.
  // Remove from it one of its grandchildren (level l+2). If we then set
  //   min_level < l + 1
  //   max_level > l + 2
  //   max_cells = 3
  // the best interior covering should contain 3 children of the initial cell,
  // that were not effected by removal of a grandchild.
  const int level = 12;
  S2CellId small_cell =
      S2CellId(S2Testing::RandomPoint()).parent(level + 2);
  S2CellId large_cell = small_cell.parent(level);
  S2CellUnion diff =
      S2CellUnion({large_cell}).Difference(S2CellUnion({small_cell}));
  S2RegionCoverer::Options options;
  options.set_max_cells(3);
  options.set_max_level(level + 3);
  options.set_min_level(level);
  S2RegionCoverer coverer(options);
  vector<S2CellId> interior;
  coverer.GetInteriorCovering(diff, &interior);
  ASSERT_EQ(interior.size(), 3);
  for (int i = 0; i < 3; ++i) {
    EXPECT_EQ(interior[i].level(), level + 1);
  }
}
