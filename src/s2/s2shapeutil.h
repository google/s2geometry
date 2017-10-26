// Copyright 2013 Google Inc. All Rights Reserved.
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
//
// This file contains miscellaneous S2ShapeIndex utility functions and classes:
//
//  - RangeIterator: useful for merging two or more S2ShapeIndexes.
//  - ShapeEdgeId: represents a (shape_id, edge_id) pair.
//  - ShapeEdge: represents a ShapeEdgeId plus the actual edge vertices.
//  - VisitCrossings: finds all edge intersections between S2ShapeIndexes.
//
// And functions that are mainly useful internally in the S2 library:
//
//  - FindAnyCrossing: helper function for S2Loop and S2Polygon validation.

#ifndef S2_S2SHAPEUTIL_H_
#define S2_S2SHAPEUTIL_H_

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "s2/third_party/absl/types/span.h"
#include "s2/_fpcontractoff.h"
#include "s2/s2edge_vector_shape.h"
#include "s2/s2lax_loop_shape.h"
#include "s2/s2lax_polygon_shape.h"
#include "s2/s2lax_polyline_shape.h"
#include "s2/s2loop.h"
#include "s2/s2point_vector_shape.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"
#include "s2/s2shapeindex.h"

class S2Loop;
class S2Error;

namespace s2shapeutil {

/////////////////////////////////////////////////////////////////////////////
////////////////////// Utility Functions and Classes ////////////////////////

// RangeIterator is a wrapper over S2ShapeIndex::Iterator with extra methods
// that are useful for merging the contents of two or more S2ShapeIndexes.
class RangeIterator {
 public:
  // Construct a new RangeIterator positioned at the first cell of the index.
  explicit RangeIterator(const S2ShapeIndex& index);

  // The current S2CellId and cell contents.
  S2CellId id() const { return it_.id(); }
  const S2ShapeIndexCell& cell() const { return it_.cell(); }

  // The min and max leaf cell ids covered by the current cell.  If done() is
  // true, these methods return a value larger than any valid cell id.
  S2CellId range_min() const { return range_min_; }
  S2CellId range_max() const { return range_max_; }

  void Next();
  bool done() { return it_.done(); }

  // Position the iterator at the first cell that overlaps or follows
  // "target", i.e. such that range_max() >= target.range_min().
  void SeekTo(const RangeIterator& target);

  // Position the iterator at the first cell that follows "target", i.e. the
  // first cell such that range_min() > target.range_max().
  void SeekBeyond(const RangeIterator& target);

 private:
  // Updates internal state after the iterator has been repositioned.
  void Refresh();
  S2ShapeIndex::Iterator it_;
  S2CellId range_min_, range_max_;
};

// A parameter that controls the reporting of edge intersections.
//
//  - CrossingType::INTERIOR reports intersections that occur at a point
//    interior to both edges (i.e., not at a vertex).
//
//  - CrossingType::ALL reports all intersections, even those where two edges
//    intersect only because they share a common vertex.
//
//  - CrossingType::NON_ADJACENT reports all intersections except for pairs of
//    the form (AB, BC) where both edges are from the same S2ShapeIndex.
//    (This is an optimization for checking polygon validity.)
enum class CrossingType { INTERIOR, ALL, NON_ADJACENT };

// ShapeEdgeId is a unique identifier for an edge within an S2ShapeIndex,
// consisting of a (shape_id, edge_id) pair.  It is similar to
// std::pair<int32, int32> except that it has named fields.
// It should be passed and returned by value.
struct ShapeEdgeId {
 public:
  int32 shape_id, edge_id;
  ShapeEdgeId() : shape_id(-1), edge_id(-1) {}
  ShapeEdgeId(int32 _shape_id, int32 _edge_id);
  bool operator==(ShapeEdgeId other) const;
  bool operator!=(ShapeEdgeId other) const;
  bool operator<(ShapeEdgeId other) const;
  bool operator>(ShapeEdgeId other) const;
  bool operator<=(ShapeEdgeId other) const;
  bool operator>=(ShapeEdgeId other) const;
};
std::ostream& operator<<(std::ostream& os, ShapeEdgeId id);

// Hasher for ShapeEdgeId.
// Example use: std::unordered_set<ShapeEdgeId, ShapeEdgeIdHash>.
struct ShapeEdgeIdHash {
  size_t operator()(ShapeEdgeId id) const {
    // The following preserves all bits even when edge_id < 0.
    return std::hash<uint64>()((static_cast<uint64>(id.shape_id) << 32) |
                               static_cast<uint32>(id.edge_id));
  }
};

// A class representing a ShapeEdgeId together with the two endpoints of that
// edge.  It should be passed by reference.
struct ShapeEdge {
 public:
  ShapeEdge() {}
  ShapeEdge(const S2Shape& shape, int32 edge_id);
  ShapeEdge(int32 shape_id, int32 edge_id, const S2Shape::Edge& edge);
  ShapeEdgeId id() const { return id_; }
  const S2Point& v0() const { return edge_.v0; }
  const S2Point& v1() const { return edge_.v1; }

 private:
  ShapeEdgeId id_;
  S2Shape::Edge edge_;
};

// A function that is called with pairs of crossing edges.  The function may
// return false in order to request that the algorithm should be terminated,
// i.e. no further crossings are needed.
//
// "is_interior" indicates that the crossing is at a point interior to both
// edges (i.e., not at a vertex).  (The calling function already has this
// information and it is moderately expensive to recompute.)
using EdgePairVisitor = std::function<
  bool (const ShapeEdge& a, const ShapeEdge& b, bool is_interior)>;

// Visits all pairs of crossing edges in the given S2ShapeIndex, terminating
// early if the given EdgePairVisitor function returns false (in which case
// VisitCrossings returns false as well).  "type" indicates whether all
// crossings should be visited, or only interior crossings.
//
// CAVEAT: Crossings may be visited more than once.
bool VisitCrossings(const S2ShapeIndex& index, CrossingType type,
                    const EdgePairVisitor& visitor);

// Like the above, but visits all pairs of crossing edges where one edge comes
// from each S2ShapeIndex.
//
// REQUIRES: type != CrossingType::NON_ADJACENT (not supported)
bool VisitCrossings(const S2ShapeIndex& a, const S2ShapeIndex& b,
                    CrossingType type, const EdgePairVisitor& visitor);

// Returns the total number of edges in all indexed shapes.  This method
// takes time linear in the number of shapes.
template <class S2ShapeIndexType>
int GetNumEdges(const S2ShapeIndexType& index);


/////////////////////////////////////////////////////////////////////////////
///////////////// Methods used internally by the S2 library /////////////////


// Given an S2ShapeIndex containing a single polygonal shape (e.g., an
// S2Polygon or S2Loop), return true if any loop has a self-intersection
// (including duplicate vertices) or crosses any other loop (including vertex
// crossings and duplicate edges) and set "error" to a human-readable error
// message.  Otherwise return false and leave "error" unchanged.
//
// This method is used to implement the FindValidationError methods of S2Loop
// and S2Polygon.
//
// TODO(ericv): Add an option to support S2LaxPolygonShape rules (i.e.,
// duplicate vertices and edges are allowed, but loop crossings are not).
bool FindAnyCrossing(const S2ShapeIndex& index, S2Error* error);


//////////////////   Implementation details follow   ////////////////////


inline ShapeEdgeId::ShapeEdgeId(int32 _shape_id, int32 _edge_id)
    : shape_id(_shape_id), edge_id(_edge_id) {
}

inline bool ShapeEdgeId::operator==(ShapeEdgeId other) const {
  return shape_id == other.shape_id && edge_id == other.edge_id;
}

inline bool ShapeEdgeId::operator!=(ShapeEdgeId other) const {
  return !(*this == other);
}

inline bool ShapeEdgeId::operator<(ShapeEdgeId other) const {
  if (shape_id < other.shape_id) return true;
  if (other.shape_id < shape_id) return false;
  return edge_id < other.edge_id;
}

inline bool ShapeEdgeId::operator>(ShapeEdgeId other) const {
  return other < *this;
}

inline bool ShapeEdgeId::operator<=(ShapeEdgeId other) const {
  return !(other < *this);
}

inline bool ShapeEdgeId::operator>=(ShapeEdgeId other) const {
  return !(*this < other);
}

inline ShapeEdge::ShapeEdge(const S2Shape& shape, int32 edge_id)
    : ShapeEdge(shape.id(), edge_id, shape.edge(edge_id)) {
}

inline ShapeEdge::ShapeEdge(int32 shape_id, int32 edge_id,
                            const S2Shape::Edge& edge)
    : id_(shape_id, edge_id), edge_(edge) {
}

template <class S2ShapeIndexType>
int GetNumEdges(const S2ShapeIndexType& index) {
  // There is no need to apply updates before counting edges, since shapes are
  // added or removed from the "shapes_" vector immediately.
  int num_edges = 0;
  for (int s = index.num_shape_ids(); --s >= 0; ) {
    S2Shape* shape = index.shape(s);
    if (shape == nullptr) continue;
    num_edges += shape->num_edges();
  }
  return num_edges;
}

}  // namespace s2shapeutil

#endif  // S2_S2SHAPEUTIL_H_
