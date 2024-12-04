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

#include "s2/internal/s2incident_edge_tracker.h"

#include <cstdint>
#include <utility>

#include "absl/log/absl_check.h"

namespace internal {

void S2IncidentEdgeTracker::StartShape(int shape_id) {
  nursery_.clear();
  current_shape_ = shape_id;
}

void S2IncidentEdgeTracker::AddEdge(int edge_id, const S2Shape::Edge& edge) {
  ABSL_DCHECK_GE(current_shape_, 0);

  // Push non-degenerate edges into the nursery.
  nursery_.push_back({edge.v0, edge_id});
  if (!edge.IsDegenerate()) {
    nursery_.push_back({edge.v1, edge_id});
  }
}

void S2IncidentEdgeTracker::FinishShape() {
  ABSL_DCHECK_GE(current_shape_, 0);

  using std::swap;

  // We want to keep any vertices with more than two incident edges.  We could
  // sort the array by vertex and remove any with fewer, but that would require
  // shifting the array and could turn quadratic quickly.
  //
  // Instead we'll scan forward from each vertex, swapping entries with the same
  // vertex into a contiguous range.  Once we've done all the swapping we can
  // just make sure that we have at least three edges in the range.
  const int nursery_size = nursery_.size();
  for (int start = 0; start < nursery_size;) {
    int end = start + 1;

    // Scan to the end of the array, swap entries so that entries with the same
    // vertex as the start are adjacent.
    int next = start;
    const S2Point& curr_vertex = nursery_[start].vertex;
    while (++next < nursery_size) {
      if (nursery_[next].vertex == curr_vertex) {
        swap(nursery_[next], nursery_[end++]);
      }
    }

    // Most vertices will have two incident edges (the incoming edge and the
    // outgoing edge), which aren't interesting, skip them.
    const int num_edges = end - start;
    if (ABSL_PREDICT_TRUE(num_edges <= 2)) {
      start = end;
      continue;
    }

    const IncidentEdgeKey key = {current_shape_, nursery_[start].vertex};

    // If we don't have this key yet, create it manually with a pre-sized
    // hash set to avoid rehashes.  Start with a size of 8, which is four
    // edges with a load factor of 50%.
    auto iter = incident_edge_map_.find(key);
    if (iter == incident_edge_map_.end()) {
      iter = incident_edge_map_.emplace(key, 8).first;
    }

    absl::flat_hash_set<int32_t>& edges = iter->second;
    edges.reserve(edges.size() + num_edges);
    while (start != end) {
      edges.insert({nursery_[start++].edge_id});
    }
  }
}

}  // namespace internal
