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

#include <cstdint>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/meta/type_traits.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2cell_index.h"
#include "s2/s2cell_union.h"
#include "s2/s2region.h"

using std::vector;

S2RegionSharder::S2RegionSharder(const vector<S2CellUnion>& shards) {
  for (int i = 0; i < shards.size(); ++i) {
    owned_index_.Add(shards[i], i);
  }
  owned_index_.Build();
}

absl::flat_hash_map<int, S2CellUnion> ToS2CellUnionMap(
    absl::flat_hash_map<int, vector<S2CellId>> input) {
  absl::flat_hash_map<int, S2CellUnion> output;

  for (auto iter = input.begin(); iter != input.end(); ++iter) {
    output.insert({iter->first, S2CellUnion(std::move(iter->second))});
  }

  return output;
}

absl::flat_hash_map<int, S2CellUnion> S2RegionSharder::GetIntersectionsByShard(
    const S2Region& region) const {
  vector<S2CellId> region_covering_cells;
  region.GetCellUnionBound(&region_covering_cells);
  const S2CellUnion region_covering(region_covering_cells);

  // Compute the intersection between the region covering and each shard
  // covering.
  absl::flat_hash_map<int, vector<S2CellId>> shards;
  index_->VisitIntersectingCells(
      region_covering,
      [&](const S2CellId cell_id, const S2CellIndex::Label label) {
        auto iter = shards.find(label);
        if (iter == shards.end()) {
          shards.insert({label, {cell_id}});
        } else {
          iter->second.push_back(cell_id);
        }
        return true;
      });

  // The fast covering is very loose, but often it only intersects one
  // shard.
  if (shards.size() == 1) {
    return ToS2CellUnionMap(std::move(shards));
  }

  // Now see if we can reduce shard sizes by tightening each covering.
  for (auto iter = shards.begin(); iter != shards.end();) {
    auto& cell_ids = iter->second;

    // Ensure the covering is normalized.
    const S2CellUnion covering(std::move(cell_ids));

    // Get the intersection since the shard's covering may have smaller cells
    // than the region's.  We know the region covering is normalized.
    vector<S2CellId> intersection;
    for (const S2CellId cell_id : covering.Intersection(region_covering)) {
      // Only add cells that intersect the region.  Since the fast covering
      // tends to intersect multiple shards for regions that are near a boundary
      // between shards, there is a high chance of testing the region against
      // the same cell more than once.  These tests can be expensive, but rather
      // than make this algorithm more complicated, we choose to push clients
      // toward regions that are fast, such as S2ShapeIndexRegion.
      if (region.MayIntersect(S2Cell(cell_id))) {
        intersection.push_back(cell_id);
      }
    }

    if (intersection.empty()) {
      // No intersection between the shard and the input region. Remove the
      // shard as a candidate.
      shards.erase(iter++);
    } else {
      // Replace the original intersecting shards with the finer-grained
      // intersection.
      cell_ids = intersection;
      ++iter;
    }
  }

  return ToS2CellUnionMap(std::move(shards));
}

int S2RegionSharder::GetMostIntersectingShard(const S2Region& region,
                                              const int default_shard) const {
  const absl::flat_hash_map<int, S2CellUnion> intersecting_shards =
      GetIntersectionsByShard(region);

  // Having clipped each shard covering down to its intersection with the
  // region, return the best intersection, or the default shard if there are
  // no intersecting candidates.
  int best_shard = default_shard;
  int64_t best_sum = 0;
  for (const auto& shard : intersecting_shards) {
    int64_t sum = 0;
    for (const S2CellId cell_id : shard.second) {
      sum += cell_id.lsb();
    }
    if (sum > best_sum) {
      best_shard = shard.first;
      best_sum = sum;
    }
  }

  return best_shard;
}

vector<int> S2RegionSharder::GetIntersectingShards(
    const S2Region& region) const {
  vector<int> shard_numbers;

  for (const auto& shard : GetIntersectionsByShard(region)) {
    if (!shard.second.empty()) {
      shard_numbers.push_back(shard.first);
    }
  }

  return shard_numbers;
}
