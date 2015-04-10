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

#include "s2closestedgequery.h"

#include "s1angle.h"
#include "s2.h"
#include "s2cell.h"
#include "s2cellid.h"
#include "s2edgeutil.h"

void S2ClosestEdgeQuery::Init(S2ShapeIndex const& index) {
  // This constant was tuned using the benchmarks.
  static int const kMaxBruteForceEdges = 180;
  index_ = &index;
  int total_edges = 0;
  for (int id = 0; id < index.num_shape_ids(); ++id) {
    total_edges += index.shape(id)->num_edges();
  }
  UseBruteForce(total_edges <= kMaxBruteForceEdges);
}

void S2ClosestEdgeQuery::UseBruteForce(bool use_brute_force) {
  use_brute_force_ = use_brute_force;
  if (use_brute_force) return;

  // Find the range of S2Cells spanned by the index and choose a level such
  // that the entire index can be covered with just a few cells.  These are the
  // "top level" cells.  There are at most 4 top level cells unless the index
  // spans more than one face, in which case up to 6 cells may be needed.
  top_cells_.clear();
  iter_.Init(*index_);
  if (iter_.Done()) return;
  S2ShapeIndex::Iterator next = iter_, last = iter_;
  last.Finish();
  last.Prev();
  if (next.id() != last.id()) {
    int level = next.id().GetCommonAncestorLevel(last.id()) + 1;
    // For each top-level cell, we find the range of index cells within it and
    // shrink the top-level cell if necessary so that it just covers them.
    // "next" is always the next index cell that has not been covered yet.
    S2CellId last_id = last.id().parent(level);
    for (S2CellId id = next.id().parent(level); id != last_id; id = id.next()) {
      if (id.range_max() < next.id()) continue;
      // Compute the range of index cells within the top-level cell "id".
      S2ShapeIndex::Iterator cell_first = next;
      next.Seek(id.range_max().next());
      S2ShapeIndex::Iterator cell_last = next;
      cell_last.Prev();
      AddInitialRange(cell_first, cell_last);
    }
  }
  AddInitialRange(next, last);
}

// Add an entry to top_cells_ that covers the given inclusive range of cells.
// REQUIRES: "first" and "last" have a common ancestor.
void S2ClosestEdgeQuery::AddInitialRange(S2ShapeIndex::Iterator const& first,
                                         S2ShapeIndex::Iterator const& last) {
  if (first.id() == last.id()) {
    // The range consists of a single index cell.
    top_cells_.push_back(TopCell(first.id(), first.cell()));
  } else {
    // Add the lowest common ancestor of the given range.
    int level = first.id().GetCommonAncestorLevel(last.id());
    DCHECK_GE(level, 0);
    top_cells_.push_back(TopCell(first.id().parent(level), NULL));
  }
}

S1Angle S2ClosestEdgeQuery::ApproxGetDistance(S2Point const& target,
                                              S1Angle max_error) {
  FindClosestEdge(target, max_error);
  return min_distance_.ToAngle();
}

S2Point S2ClosestEdgeQuery::ApproxProject(S2Point const& target,
                                          S1Angle max_error) {
  FindClosestEdge(target, max_error);
  if (shape_id_ < 0) return target;  // There are no edges.
  S2Point const *v0, *v1;
  index_->shape(shape_id_)->GetEdge(edge_id_, &v0, &v1);
  return S2EdgeUtil::GetClosestPoint(target, *v0, *v1);
}

void S2ClosestEdgeQuery::FindClosestEdgeBruteForce(S2Point const& target) {
  for (int id = 0; id < index_->num_shape_ids(); ++id) {
    S2Shape const* shape = index_->shape(id);
    int num_edges = shape->num_edges();
    for (int e = 0; e < num_edges; ++e) {
      S2Point const *a, *b;
      shape->GetEdge(e, &a, &b);
      if (S2EdgeUtil::UpdateMinDistance(target, *a, *b, &min_distance_)) {
        shape_id_ = id;
        edge_id_ = e;
      }
    }
  }
}

void S2ClosestEdgeQuery::FindClosestEdge(S2Point const& target,
                                         S1Angle max_error) {
  target_ = target;
  max_error_ = S1ChordAngle(max_error);
  min_distance_ = goal_distance_ = S1ChordAngle::Infinity();
  shape_id_ = edge_id_ = -1;
  if (use_brute_force_) {
    return FindClosestEdgeBruteForce(target);
  }
  queue_.clear();
  for (int i = 0; i < top_cells_.size(); ++i) {
    EnqueueCell(top_cells_[i].first, top_cells_[i].second);
  }
  // Repeatedly find the closest S2Cell to "target" and either split it into
  // its four children or process all of its edges.
  while (!queue_.empty()) {
    // We need to copy the top entry before removing it, and we need to remove
    // it before adding any new entries to the queue.
    QueueEntry entry = queue_.top();
    queue_.pop();
    S2CellId id = entry.id;
    // If all remaining cells are too far away, we can stop now.
    if (entry.distance >= goal_distance_)
      break;

    // If this is already known to be an index cell, just process it.
    if (entry.index_cell != NULL) {
      UpdateDistance(entry);
    } else {
      // Otherwise split the cell into its four children.  Before adding a
      // child back to the queue, we first check whether it is empty.  We do
      // this in two seek operations rather than four by seeking to the key
      // between children 0 and 1 and to the key between children 2 and 3.
      iter_.Seek(id.child(1).range_min());
      if (!iter_.Done() && iter_.id() <= id.child(1).range_max()) {
        EnqueueCurrentCell(id.child(1));
      }
      if (!iter_.AtBegin()) {
        iter_.Prev();
        if (iter_.id() >= id.range_min()) {
          EnqueueCurrentCell(id.child(0));
        }
      }
      iter_.Seek(id.child(3).range_min());
      if (!iter_.Done() && iter_.id() <= id.range_max()) {
        EnqueueCurrentCell(id.child(3));
      }
      if (!iter_.AtBegin()) {
        iter_.Prev();
        if (iter_.id() >= id.child(2).range_min()) {
          EnqueueCurrentCell(id.child(2));
        }
      }
    }
  }
}

// Process all the edges of the given index cell.
void S2ClosestEdgeQuery::UpdateDistance(QueueEntry const& entry) {
  S2ShapeIndexCell const* index_cell = entry.index_cell;
  for (int s = 0; s < index_cell->num_shapes(); ++s) {
    S2ClippedShape const& clipped = index_cell->clipped(s);
    S2Shape const* shape = index_->shape(clipped.shape_id());
    for (int j = 0; j < clipped.num_edges(); ++j) {
      int edge_id = clipped.edge(j);
      S2Point const *a, *b;
      shape->GetEdge(edge_id, &a, &b);
      if (S2EdgeUtil::UpdateMinDistance(target_, *a, *b, &goal_distance_)) {
        min_distance_ = goal_distance_;
        goal_distance_ = min_distance_ - max_error_;
        shape_id_ = shape->id();
        edge_id_ = edge_id;
        // If no edge in this cell can further reduce the distance by at least
        // max_error_, we can stop now.
        if (entry.distance >= goal_distance_) return;
      }
    }
  }
}

// Enqueue the given cell id.
// REQUIRES: iter_ is positioned at a cell contained by "id".
inline void S2ClosestEdgeQuery::EnqueueCurrentCell(S2CellId id) {
  DCHECK(id.contains(iter_.id()));
  if (iter_.id() == id) {
    EnqueueCell(id, iter_.cell());
  } else {
    EnqueueCell(id, NULL);
  }
}

// Return the number of edges in the given index cell.
int CountEdges(S2ShapeIndexCell const* cell) {
  int count = 0;
  for (int s = 0; s < cell->num_shapes(); ++s) {
    count += cell->clipped(s).num_edges();
  }
  return count;
}

// Add the given cell id to the queue.  "index_cell" is the corresponding
// S2ShapeIndexCell, or NULL if "id" is not an index cell.
void S2ClosestEdgeQuery::EnqueueCell(S2CellId id,
                                     S2ShapeIndexCell const* index_cell) {
  if (index_cell) {
    // If this index cell has only a few edges, then it is faster to check
    // them directly rather than computing the minimum distance to the S2Cell
    // and inserting it into the queue.
    static int const kMinEdgesToEnqueue = 10;
    int num_edges = CountEdges(index_cell);
    if (num_edges == 0) return;
    if (num_edges < kMinEdgesToEnqueue) {
      // Set "distance" to zero in order to force all edges to be examined.
      UpdateDistance(QueueEntry(S1ChordAngle::Zero(), id, index_cell));
      return;
    }
  }
  // Otherwise compute the minimum distance to any point in the cell and add
  // it to the priority queue.
  S2Cell cell(id);
  S1ChordAngle distance(cell.GetDistance(target_));  // XXX FIX  XXX
  if (distance >= goal_distance_) return;
  queue_.push(QueueEntry(distance, id, index_cell));
}