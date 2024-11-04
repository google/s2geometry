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

#ifndef S2_INTERNAL_S2INCIDENT_EDGE_TRACKER_H_
#define S2_INTERNAL_S2INCIDENT_EDGE_TRACKER_H_

#include <algorithm>
#include <cstdint>
#include <utility>
#include <vector>

#include "absl/container/btree_map.h"
#include "absl/container/flat_hash_set.h"
#include "s2/s2point.h"
#include "s2/s2shape.h"

namespace internal {

// A tuple of (shape id, vertex) that compares by shape id.
struct IncidentEdgeKey {
  IncidentEdgeKey() = default;

  explicit IncidentEdgeKey(int shape_id)
      : shape_id(shape_id) { /* Leave vertex undefined */
  }

  IncidentEdgeKey(int shape_id, S2Point vertex)
      : shape_id(shape_id), vertex(vertex) {}

  int32_t shape_id;
  S2Point vertex;

  // We need a strict ordering to be a valid key for an ordered container, but
  // we don't actually care about the ordering of the vertices (as long as
  // they're grouped by shape id).  Vertices are 3D points so they don't have a
  // natural ordering, so we'll just compare them lexicographically.
  bool operator<(const IncidentEdgeKey& b) const {
    if (shape_id < b.shape_id) return true;
    if (shape_id > b.shape_id) return false;
    for (int i = 0; i < 3; ++i) {
      if (vertex[i] < b.vertex[i]) return true;
      if (vertex[i] > b.vertex[i]) return false;
    }
    return false;  // Must be equal.
  }
};

// Map from keys of (shape id, vertex) to a set of edge_ids.  We can and do
// encounter the same edges multiple times, so we need to use a set to
// deduplicate edges as they're inserted.
using IncidentEdgeSet =
    absl::btree_map<IncidentEdgeKey, absl::flat_hash_set<int32_t>>;

// A class for detecting and tracking shape edges that are incident on the same
// vertex.  Edges of multiple shapes may be tracked, but lookup is by shape id
// and vertex: there is no facility to get all edges of all shapes at a vertex.
// Edge vertices must compare exactly equal to be considered the same vertex, no
// tolerance is applied as this isn't intended for e.g.: snapping shapes
// together, which S2Builder does better and more robustly.
//
// To use, instantiate and then add edges with one or more sequences of calls,
// where each sequence begins with StartShape(), followed by AddEdge() calls to
// add edges for that shape, and ends with FinishShape().  Those sequences do
// not need to visit shapes or edges in order. Then, call incidentEdges() to get
// the resulting map from IncidentEdgeKeys (which are shapeId, vertex pairs) to
// a set of edgeIds of the shape that are incident to that vertex..
//
// This works on a block of edges at a time, meaning that to detect incident
// edges on a particular vertex, we must have at least three edges incident
// at that vertex when FinishShape() is called. We don't maintain partial
// information between calls.  However, subject to this constraint, a single
// shape's edges may be defined with multiple sequences of StartShape(),
// AddEdge()... , FinishShape() calls.
//
// The reason for this is simple: most edges don't have more than two incident
// edges (the incoming and outgoing edge). If we had to maintain information
// between calls, we'd end up with a map that contains every vertex, to no
// benefit.  Instead, when FinishShape() is called, we discard vertices that
// contain two or fewer incident edges.
//
// In principle this isn't a real limitation because generally we process an
// S2ShapeIndex cell at a time, and, if a vertex has multiple edges, we'll see
// all the edges in the same cell as the vertex, and, in general, it's possible
// to aggregate edges before calling.
//
// The tracker maintains incident edges until it's cleared.  If you call it with
// each cell from an S2ShapeIndex, then at the end you will have all the
// incident edge information for the whole index.  If only a subset is needed,
// call Reset() when you're done.
class S2IncidentEdgeTracker {
 public:
  // Add edges to the edge tracker.  After calling, any vertices with multiple
  // (> 2) incident edges will appear in the incident edge map.
  void StartShape(int shape_id);
  void AddEdge(int edge_id, const S2Shape::Edge& edge);
  void FinishShape();

  // Clear any accumulated state.
  void Reset() { incident_edge_map_.clear(); }

  // Returns a const reference to the incident edge map.
  const IncidentEdgeSet& IncidentEdges() const { return incident_edge_map_; }

 private:
  // Scratch space to temporarily store incident edges in.  Having a separate
  // space for these lets us remove vertices with only one incident edge
  // easily.
  struct VertexEdge {
    S2Point vertex;
    int32_t edge_id;
  };
  std::vector<VertexEdge> nursery_;

  int current_shape_ = -1;

  // Actual incident edge map.
  IncidentEdgeSet incident_edge_map_;
};

}  // namespace internal

#endif  // S2_INTERNAL_S2INCIDENT_EDGE_TRACKER_H_
