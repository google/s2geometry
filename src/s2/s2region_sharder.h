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

#ifndef S2_S2REGION_SHARDER_H_
#define S2_S2REGION_SHARDER_H_

#include <vector>

#include "absl/container/flat_hash_map.h"
#include "s2/s2cell_index.h"
#include "s2/s2cell_union.h"
#include "s2/s2region.h"

// S2RegionSharder implements a sharding function that provides shard IDs whose
// boundaries intersect with an S2Region. This class is especially suited to
// testing regions against shards that are usually very different in size from
// the regions; in that case, S2CellUnion coverings of the region tend to cover
// too complex (small regions are often contained by a single cell of the
// shard's covering).
class S2RegionSharder {
 public:
  // Construct a new S2RegionSharder from a pre-existing S2CellIndex.  The index
  // must remain alive at least as long as the S2RegionSharder.
  explicit S2RegionSharder(const S2CellIndex* index) : index_(index) {}

  // Construct a new S2RegionSharder from the given list of S2CellUnions.
  explicit S2RegionSharder(const std::vector<S2CellUnion>& shards);

  S2RegionSharder(const S2RegionSharder&) = delete;
  S2RegionSharder& operator=(const S2RegionSharder&) = delete;

  // Returns an index into the original vector of S2CellUnions indicating the
  // covering with the most overlap to the input region.  If no shards overlap
  // with the input region, returns 'default_shard'.
  int GetMostIntersectingShard(const S2Region& region, int default_shard) const;

  // Returns a list of shard numbers which intersect with the input 'region'.
  // Shard numbers are not guaranteed to be sorted in any particular order.  If
  // no shards overlap, returns an empty vector.
  std::vector<int> GetIntersectingShards(const S2Region& region) const;

 private:
  absl::flat_hash_map<int, S2CellUnion> GetIntersectionsByShard(
      const S2Region& region) const;

  S2CellIndex owned_index_;
  const S2CellIndex* index_ = &owned_index_;
};

#endif  // S2_S2REGION_SHARDER_H_
