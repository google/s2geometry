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

#ifndef S2_GEOMETRY_S2CLOSESTEDGEQUERY_H_
#define S2_GEOMETRY_S2CLOSESTEDGEQUERY_H_

#include <algorithm>
#include "priority_queue_sequence.h"
#include "s1angle.h"
#include "s1chordangle.h"
#include "fpcontractoff.h"
#include "s2.h"
#include "s2cellid.h"
#include "s2shapeindex.h"
#include "util/gtl/inlined_vector.h"

// S2ClosestEdgeQuery is a helper class for finding the minimum distance from
// a given point to a set of edges (stored in an S2ShapeIndex).
class S2ClosestEdgeQuery {
 public:
  // Convenience constructor that calls Init().
  explicit S2ClosestEdgeQuery(S2ShapeIndex const& index);

  // Default constructor; requires Init() to be called.
  S2ClosestEdgeQuery() {}
  ~S2ClosestEdgeQuery();

  // REQUIRES: "index" is not modified after this method is called.
  void Init(S2ShapeIndex const& index);

  // Return the minimum distance from the given point to any edge of the
  // S2ShapeIndex.  If there are no edges, return S1Angle::Infinity().  The
  // edge that determined the minimum distance can be obtained by calling
  // shape_id() and edge_id() after calling this method.
  //
  // REQUIRES: "target" is unit length.
  S1Angle GetDistance(S2Point const& target);

  // Return the closest point on any edge of the S2ShapeIndex to the given
  // point.  If the index has no edges, return the input argument.  The edge
  // that determined the minimum distance can be obtained by calling
  // shape_id() and edge_id() after calling this method.
  //
  // REQUIRES: "target" is unit length.
  S2Point Project(S2Point const& target);

  // A faster approximate version of GetDistance().  The result is never
  // smaller than the true distance, but it may be larger by up to
  // "max_error".  shape_id() and edge_id() can be used to obtain the edge
  // that determined the approximate minimum distance.
  S1Angle ApproxGetDistance(S2Point const& target, S1Angle max_error);

  // A faster approximate version of Project().  The result will alway be as
  // close as possible to some edge of the index, but there may be some other
  // edge that is closer by up to "max_error".  shape_id() and edge_id() can
  // be used to obtain the edge on which the given point was projected.
  S2Point ApproxProject(S2Point const& target, S1Angle max_error);

  // Return the shape and edge id of the edge that determined the result of
  // the most recent query (i.e., the closest edge to the target point).
  // Return -1 if the shape index has no edges.
  int shape_id() const { return shape_id_; }
  int edge_id() const { return edge_id_; }

  // Manually determine whether distances are computed using "brute force"
  // (i.e., by examining every edge) rather than using the S2ShapeIndex.  This
  // is useful for testing, benchmarking, and debugging.
  //
  // Changing this setting may cause Project() to return a different result
  // (e.g., when two edges are nearly equidistant from the target).  Also note
  // that the brute force algorithm ignores "max_error".
  //
  // REQUIRES: Init() has been called.
  void UseBruteForce(bool use_brute_force);

 private:
  // The algorithm maintains a priority queue of S2CellIds that contain at
  // least one S2ShapeIndexCell, sorted in increasing order of distance from
  // the target point.
  struct QueueEntry {
    // A lower bound on the distance from the target point to any edge point
    // within "id".  This is the key of the priority queue.
    S1ChordAngle distance;

    // The cell being queued.
    S2CellId id;

    // If "id" belongs to the index, this field stores the corresponding
    // S2ShapeIndexCell.  Otherwise "id" is a proper ancestor of one or more
    // S2ShapeIndexCells and this field stores NULL.  The purpose of this
    // field is to avoid an extra Seek() when the queue entry is processed.
    S2ShapeIndexCell const* index_cell;

    QueueEntry(S1ChordAngle _distance, S2CellId _id,
               S2ShapeIndexCell const* _index_cell)
        : distance(_distance), id(_id), index_cell(_index_cell) {
    }
    bool operator<(QueueEntry const& other) const {
      // The priority queue returns the largest elements first, so we want the
      // "largest" entry to have the smallest distance.
      return distance > other.distance;
    }
  };

  // Private methods are documented in the .cc file.
  void FindClosestEdge(S2Point const& target, S1Angle max_error);
  void FindClosestEdgeBruteForce(S2Point const& target);
  void UpdateDistance(QueueEntry const& entry);
  void EnqueueCell(S2CellId id, S2ShapeIndexCell const* index_cell);
  void EnqueueCurrentCell(S2CellId id);
  void AddInitialRange(S2ShapeIndex::Iterator const& first,
                       S2ShapeIndex::Iterator const& last);

  ////////// Fields that are constant after Init() is called /////////////

  S2ShapeIndex const* index_;

  // If the index has few edges, it is cheaper to use a brute force algorithm.
  bool use_brute_force_;

  // During Init() we precompute the top-level S2CellIds that will be added to
  // the priority queue.  There can be at most 6 of these cells.
  typedef std::pair<S2CellId, S2ShapeIndexCell const*> TopCell;
  util::gtl::InlinedVector<TopCell, 6> top_cells_;

  ////////// Fields that are constant during a query /////////////

  S2Point target_;
  S1ChordAngle max_error_;

  ////////// Fields that are updated during a query /////////////

  // These fields store the tentative result of the query.
  S1ChordAngle min_distance_;
  int shape_id_, edge_id_;

  // We only update "min_distance_" if it would decrease by at least
  // "max_error_", since computing the actual distance is fairly expensive.
  // We do this by maintaining a value "goal_distance_" equal to
  // (min_distance_ - max_error_) and comparing distance bounds against it.
  S1ChordAngle goal_distance_;

  // A temporary declared here to avoid initializing iterators multiple times.
  S2ShapeIndex::Iterator iter_;

  // We use a priority queue that exposes the underlying vector because the
  // standard STL one does not support clearing the remaining queue entries.
  typedef priority_queue_sequence<QueueEntry> CellQueue;
  CellQueue queue_;
};


//////////////////   Implementation details follow   ////////////////////


inline S2ClosestEdgeQuery::S2ClosestEdgeQuery(S2ShapeIndex const& index) {
  Init(index);
}

inline S1Angle S2ClosestEdgeQuery::GetDistance(S2Point const& target) {
  return ApproxGetDistance(target, S1Angle::Zero());
}

inline S2Point S2ClosestEdgeQuery::Project(S2Point const& target) {
  return ApproxProject(target, S1Angle::Zero());
}

#endif  // S2_GEOMETRY_S2CLOSESTEDGEQUERY_H_
