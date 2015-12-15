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
#include "s2cap.h"
#include "s2cell.h"
#include "s2cellid.h"
#include "s2cellunion.h"
#include "s2edgeutil.h"
#include "s2regioncoverer.h"

using std::vector;

S2ClosestEdgeQuery::~S2ClosestEdgeQuery() {
  // Prevent inline destructor bloat by providing a definition.
}

void S2ClosestEdgeQuery::Init(S2ShapeIndex const& index) {
  index_ = &index;
  max_edges_ = std::numeric_limits<int>::max();
  max_distance_ = S1Angle::Infinity();
  max_error_arg_ = S1Angle::Zero();
  Reset();
}

void S2ClosestEdgeQuery::Reset() {
  // This constant was tuned using the benchmarks.
  static int const kMaxBruteForceEdges = 180;

  results_.clear();
  iter_.Init(*index_);
  UseBruteForce(index_->GetNumEdges() <= kMaxBruteForceEdges);
}

void S2ClosestEdgeQuery::UseBruteForce(bool use_brute_force) {
  use_brute_force_ = use_brute_force;
  if (use_brute_force) return;

  // Find the range of S2Cells spanned by the index and choose a level such
  // that the entire index can be covered with just a few cells.  These are
  // the "top-level" cells.  There are two cases:
  //
  //  - If the index spans more than one face, then there is one top-level cell
  // per spanned face, just big enough to cover the index cells on that face.
  //
  //  - If the index spans only one face, then we find the smallest cell "C"
  // that covers the index cells on that face (just like the case above).
  // Then for each of the 4 children of "C", if the child contains any index
  // cells then we create a top-level cell that is big enough to just fit
  // those index cells (i.e., shrinking the child as much as possible to fit
  // its contents).  This essentially replicates what would happen if we
  // started with "C" as the top-level cell, since "C" would immediately be
  // split, except that we take the time to prune the children further since
  // this will save work on every subsequent query.
  index_covering_.clear();
  index_cells_.clear();
  iter_.Init(*index_);
  if (iter_.Done()) return;  // Empty index.

  index_covering_.reserve(6);
  // Don't need to reserve index_cells_ since it is an InlinedVector.
  S2ShapeIndex::Iterator next = iter_, last = iter_;
  last.Finish();
  last.Prev();
  if (next.id() != last.id()) {
    // The index has at least two cells.  Choose a level such that the entire
    // index can be spanned with at most 6 cells (if the index spans multiple
    // faces) or 4 cells (it the index spans a single face).
    int level = next.id().GetCommonAncestorLevel(last.id()) + 1;

    // Visit each potential top-level cell except the last (handled below).
    S2CellId last_id = last.id().parent(level);
    for (S2CellId id = next.id().parent(level); id != last_id; id = id.next()) {
      // Skip any top-level cells that don't contain any index cells.
      if (id.range_max() < next.id()) continue;

      // Find the range of index cells contained by this top-level cell and
      // then shrink the cell if necessary so that it just covers them.
      S2ShapeIndex::Iterator cell_first = next;
      next.Seek(id.range_max().next());
      S2ShapeIndex::Iterator cell_last = next;
      cell_last.Prev();
      AddInitialRange(cell_first, cell_last);
    }
  }
  AddInitialRange(next, last);
}

// Add an entry to index_covering_ and index_cells_ that covers the given
// inclusive range of cells.
//
// REQUIRES: "first" and "last" have a common ancestor.
void S2ClosestEdgeQuery::AddInitialRange(S2ShapeIndex::Iterator const& first,
                                         S2ShapeIndex::Iterator const& last) {
  if (first.id() == last.id()) {
    // The range consists of a single index cell.
    index_covering_.push_back(first.id());
    index_cells_.push_back(first.cell());
  } else {
    // Add the lowest common ancestor of the given range.
    int level = first.id().GetCommonAncestorLevel(last.id());
    DCHECK_GE(level, 0);
    index_covering_.push_back(first.id().parent(level));
    index_cells_.push_back(nullptr);
  }
}

S2Point S2ClosestEdgeQuery::GetClosestPointOnEdge(int i) const {
  S2Point const *v0, *v1;
  GetEdge(i, &v0, &v1);
  return target_->GetClosestPointOnEdge(*v0, *v1);
}

S1Angle S2ClosestEdgeQuery::GetDistance(S2Point const& target) {
  FindClosestEdge(target);
  if (num_edges() == 0) return S1Angle::Infinity();
  return distance(0);
}

S2Point S2ClosestEdgeQuery::Project(S2Point const& target) {
  FindClosestEdge(target);
  if (num_edges() == 0) return target;
  return GetClosestPointOnEdge(0);
}

void S2ClosestEdgeQuery::FindClosestEdges(S2Point const& target) {
  target_.reset(new PointTarget(target));
  FindClosestEdgesToTarget();
}

void S2ClosestEdgeQuery::FindClosestEdgesToEdge(S2Point const& a,
                                                S2Point const& b) {
  target_.reset(new EdgeTarget(a, b));
  FindClosestEdgesToTarget();
}

void S2ClosestEdgeQuery::FindClosestEdgesToTarget() {
  max_distance_limit_ = S1ChordAngle(max_distance_);
  max_error_ = S1ChordAngle(max_error_arg_);
  DCHECK(tmp_results_.empty());
  tmp_result_singleton_ = Result(S1ChordAngle::Infinity(), -1, -1);
  results_.clear();

  if (use_brute_force_) {
    FindClosestEdgesBruteForce();
  } else {
    FindClosestEdgesOptimized();
  }
  if (max_edges_ > 1) {
    results_.reserve(tmp_results_.size());
    for (ResultSet::const_iterator it = tmp_results_.begin(),
             end = tmp_results_.end(); it != end; ++it) {
      results_.push_back(*it);
    }
    tmp_results_.clear();
  } else if (tmp_result_singleton_.shape_id >= 0) {
    results_.push_back(tmp_result_singleton_);
  }
}

void S2ClosestEdgeQuery::FindClosestEdgesBruteForce() {
  for (int id = 0; id < index_->num_shape_ids(); ++id) {
    S2Shape const* shape = index_->shape(id);
    if (shape == nullptr) continue;
    int num_edges = shape->num_edges();
    for (int e = 0; e < num_edges; ++e) {
      MaybeAddResult(shape, e);
    }
  }
}

void S2ClosestEdgeQuery::FindClosestEdgesOptimized() {
  InitQueue();
  // Repeatedly find the closest S2Cell to "target" and either split it into
  // its four children or process all of its edges.
  while (!queue_.empty()) {
    // We need to copy the top entry before removing it, and we need to
    // remove it before adding any new entries to the queue.
    QueueEntry entry = queue_.top();
    queue_.pop();
    if (entry.distance >= max_distance_limit_) {
      queue_ = CellQueue();  // Clear any remaining entries.
      break;
    }
    // If this is already known to be an index cell, just process it.
    if (entry.index_cell != nullptr) {
      ProcessEdges(entry);
      continue;
    }
    // Otherwise split the cell into its four children.  Before adding a
    // child back to the queue, we first check whether it is empty.  We do
    // this in two seek operations rather than four by seeking to the key
    // between children 0 and 1 and to the key between children 2 and 3.
    S2CellId id = entry.id;
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

void S2ClosestEdgeQuery::InitQueue() {
  DCHECK(queue_.empty());

  // Optimization: if the user is searching for just the closest edge, and the
  // target happens to be intersect an index cell, then we try to limit the
  // search region to a small disc by first processing the edges in that cell.
  // This sets max_distance_limit_ based on the closest edge in that cell,
  // which we can then use to limit the search area.  This does mean that the
  // cell containing "target" will be processed twice, but in general this is
  // still faster.
  if (max_edges_ == 1 && iter_.Locate(target_->center())) {
    ProcessEdges(QueueEntry(S1ChordAngle::Zero(), iter_.id(), iter_.cell()));
  }
  if (max_distance_limit_ == S1ChordAngle::Infinity()) {
    // Start with the precomputed index covering.
    for (int i = 0; i < index_covering_.size(); ++i) {
      EnqueueCell(index_covering_[i], index_cells_[i]);
    }
  } else {
    // Compute a covering of the search disc and intersect it with the
    // precomputed index covering.
    S2RegionCoverer coverer;
    coverer.set_max_cells(4);
    S2Cap search_cap(target_->center(),
                     target_->radius() + max_distance_limit_.ToAngle());
    coverer.GetFastCovering(search_cap, &max_distance_covering_);
    S2CellUnion::GetIntersection(index_covering_, max_distance_covering_,
                                 &initial_cells_);

    // Now we need to clean up the initial cells to ensure that they all
    // contain at least one cell of the S2ShapeIndex.  (Some may not intersect
    // the index at all, while other may be descendants of an index cell.)
    for (int i = 0, j = 0; i < initial_cells_.size(); ) {
      S2CellId id_i = initial_cells_[i];
      // Find the top-level cell that contains this initial cell.
      while (index_covering_[j].range_max() < id_i) ++j;
      S2CellId id_j = index_covering_[j];
      if (id_i == id_j) {
        // This initial cell is one of the top-level cells.  Use the
        // precomputed S2ShapeIndexCell pointer to avoid an index seek.
        EnqueueCell(id_j, index_cells_[j]);
        ++i, ++j;
      } else {
        // This initial cell is a proper descendant of a top-level cell.
        // Check how it is related to the cells of the S2ShapeIndex.
        S2ShapeIndex::CellRelation r = iter_.Locate(id_i);
        if (r == S2ShapeIndex::INDEXED) {
          // This cell is a descendant of an index cell.  Enqueue it and skip
          // any other initial cells that are also descendants of this cell.
          EnqueueCell(iter_.id(), iter_.cell());
          S2CellId const last_id = iter_.id().range_max();
          while (++i < initial_cells_.size() && initial_cells_[i] <= last_id)
            continue;
        } else {
          // Enqueue the cell only if it contains at least one index cell.
          if (r == S2ShapeIndex::SUBDIVIDED) EnqueueCell(id_i, nullptr);
          ++i;
        }
      }
    }
  }
}

void S2ClosestEdgeQuery::MaybeAddResult(S2Shape const* shape, int edge_id) {
  S2Point const *v0, *v1;
  shape->GetEdge(edge_id, &v0, &v1);

  S1ChordAngle distance = max_distance_limit_;
  if (!target_->UpdateMinDistance(*v0, *v1, &distance)) return;

  if (max_edges_ == 1) {
    // Optimization for the common case where only the closest edge is wanted.
    tmp_result_singleton_ = Result(distance, shape->id(), edge_id);
    max_distance_limit_ = distance - max_error_;
  } else {
    // Add this edge to tmp_results_.  Note that even if we already have
    // enough edges, we can't erase an element before insertion because the
    // "new" edge might in fact be a duplicate.
    tmp_results_.insert(Result(distance, shape->id(), edge_id));
    int size = tmp_results_.size();
    if (size >= max_edges_) {
      if (size > max_edges_) {
        tmp_results_.erase(--tmp_results_.end());
      }
      max_distance_limit_ = (--tmp_results_.end())->distance - max_error_;
    }
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

// Process all the edges of the given index cell.
void S2ClosestEdgeQuery::ProcessEdges(QueueEntry const& entry) {
  S2ShapeIndexCell const* index_cell = entry.index_cell;
  for (int s = 0; s < index_cell->num_shapes(); ++s) {
    S2ClippedShape const& clipped = index_cell->clipped(s);
    S2Shape const* shape = index_->shape(clipped.shape_id());
    for (int j = 0; j < clipped.num_edges(); ++j) {
      MaybeAddResult(shape, clipped.edge(j));
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
    EnqueueCell(id, nullptr);
  }
}

// Add the given cell id to the queue.  "index_cell" is the corresponding
// S2ShapeIndexCell, or nullptr if "id" is not an index cell.
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
      // Set "distance" to zero to avoid the expense of computing it.
      ProcessEdges(QueueEntry(S1ChordAngle::Zero(), id, index_cell));
      return;
    }
  }
  // Otherwise compute the minimum distance to any point in the cell and add
  // it to the priority queue.
  S2Cell cell(id);
  S1ChordAngle distance = target_->GetDistance(cell);
  if (distance >= max_distance_limit_) return;
  queue_.push(QueueEntry(distance, id, index_cell));
}
