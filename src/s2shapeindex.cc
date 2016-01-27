// Copyright 2012 Google Inc. All Rights Reserved.
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

#include "s2shapeindex.h"

#include <cmath>
#include <algorithm>

#include <gflags/gflags.h>

#include "base/atomicops.h"
#include "base/spinlock.h"
#include "r1interval.h"
#include "r2.h"
#include "r2rect.h"
#include "s2cellid.h"
#include "s2cellunion.h"
#include "s2edgeutil.h"
#include "s2paddedcell.h"
#include "util/gtl/stl_util.h"

using std::max;
using std::vector;

// FLAGS_s2shapeindex_default_max_edges_per_cell
//
// The default maximum number of edges per cell (not counting "long" edges).
// If a cell has more than this many edges, and it is not a leaf cell, then it
// is subdivided.  This flag can be overridden via S2ShapeIndexOptions.
// Reasonable values range from 10 to about 50 or so.
DEFINE_int32(
    s2shapeindex_default_max_edges_per_cell, 10,
    "Default maximum number of edges (not counting 'long' edges) per cell; "
    "reasonable values range from 10 to 50.  Small values makes queries "
    "faster, while large values make construction faster and use less "
    "memory.");

// FLAGS_s2shapeindex_tmp_memory_budget_mb
//
// Attempt to limit the amount of temporary memory allocated while building or
// updating an S2ShapeIndex to at most this value.  This is achieved by
// splitting the updates into multiple batches when necessary.  (The memory
// required is proportional to the number of edges being updated at once.)
//
// Note that this limit is not a hard guarantee, for several reasons:
//  (1) the memory estimates are only approximations;
//  (2) all edges in a given shape are added or removed at once, so shapes
//      with huge numbers of edges may exceed the budget;
//  (3) shapes being removed are always processed in a single batch.  (This
//      could be fixed, but it seems better to keep the code simpler for now.)
DEFINE_int32(
    s2shapeindex_tmp_memory_budget_mb, 100,
    "Attempt to limit the amount of temporary memory used by S2ShapeIndex "
    "when creating or updating very large indexes to at most this value.  "
    "If more memory than this is needed, updates will automatically "
    "be split into batches internally.");

// FLAGS_s2shapeindex_cell_size_to_long_edge_ratio
//
// The cell size relative to the length of an edge at which it is first
// considered to be "long".  Long edges do not contribute toward the decision
// to subdivide a cell further.  For example, a value of 2.0 means that the
// cell must be at least twice the size of the edge in order for that edge to
// be counted.  There are two reasons for not counting long edges: (1) such
// edges typically need to be propagated to several children, which increases
// time and memory costs without much benefit, and (2) in pathological cases,
// many long edges close together could force subdivision to continue all the
// way to the leaf cell level.
DEFINE_double(
    s2shapeindex_cell_size_to_long_edge_ratio, 1.0,
    "The cell size relative to the length of an edge at which it is first "
    "considered to be 'long'.  Long edges do not contribute to the decision "
    "to subdivide a cell further.  The size and speed of the index are "
    "typically not very sensitive to this parameter.  Reasonable values range "
    "from 0.1 to 10, with smaller values causing more aggressive subdivision "
    "of long edges grouped closely together.");

// The total error when clipping an edge comes from two sources:
// (1) Clipping the original spherical edge to a cube face (the "face edge").
//     The maximum error in this step is S2EdgeUtil::kFaceClipErrorUVCoord.
// (2) Clipping the face edge to the u- or v-coordinate of a cell boundary.
//     The maximum error in this step is S2EdgeUtil::kEdgeClipErrorUVCoord.
// Finally, since we encounter the same errors when clipping query edges, we
// double the total error so that we only need to pad edges during indexing
// and not at query time.
double const S2ShapeIndex::kCellPadding = 2 *
    (S2EdgeUtil::kFaceClipErrorUVCoord + S2EdgeUtil::kEdgeClipErrorUVCoord);

bool S2ClippedShape::ContainsEdge(int id) const {
  // Linear search is fast because the number of edges per shape is typically
  // very small (less than 10).
  for (int e = 0; e < num_edges(); ++e) {
    if (edge(e) == id) return true;
  }
  return false;
}

// Initialize an S2ClippedShape to hold the given number of edges.
inline void S2ClippedShape::Init(int32 shape_id, int32 num_edges) {
  shape_id_ = shape_id;
  num_edges_ = num_edges;
  contains_center_ = false;
  if (!is_inline()) {
    edges_ = new int32[num_edges];
  }
}

// Set "contains_center_" to indicate whether this clipped shape contains the
// center of the cell to which it belongs.
inline void S2ClippedShape::set_contains_center(bool contains_center) {
  contains_center_ = contains_center;
}

// Set the i-th edge of this clipped shape to be the given edge of the
// original shape.
inline void S2ClippedShape::set_edge(int i, int edge) {
  if (is_inline()) {
    inline_edges_[i] = edge;
  } else {
    edges_[i] = edge;
  }
}

// Free any memory allocated by this S2ClippedShape.  We don't do this in
// the destructor because S2ClippedShapes are copied by STL code, and we
// don't want to repeatedly copy and free the edge data.  Instead the data
// is owned by the containing S2ShapeIndexCell.
inline void S2ClippedShape::Destruct() {
  if (!is_inline()) delete[] edges_;
}

S2ShapeIndexCell::~S2ShapeIndexCell() {
  // Free memory for all shapes owned by this cell.
  for (int i = 0; i < shapes_.size(); ++i) {
    shapes_[i].Destruct();
  }
  shapes_.clear();
}

S2ClippedShape const*
S2ShapeIndexCell::find_clipped(int shape_id) const {
  // Linear search is fine because the number of shapes per cell is typically
  // very small (most often 1), and is large only for pathological inputs
  // (e.g. very deeply nested loops).
  for (int i = 0; i < shapes_.size(); ++i) {
    if (shapes_[i].shape_id() == shape_id) return &shapes_[i];
  }
  return nullptr;
}

// Allocate room for "n" additional clipped shapes in the cell, and return a
// pointer to the first new clipped shape.  Expects that all new clipped
// shapes will have a larger shape id than any current shape, and that shapes
// will be added in increasing shape id order.
S2ClippedShape* S2ShapeIndexCell::add_shapes(int n) {
  int size = shapes_.size();
  shapes_.resize(size + n);
  return &shapes_[size];
}

S2ShapeIndexOptions::S2ShapeIndexOptions()
    : max_edges_per_cell_(FLAGS_s2shapeindex_default_max_edges_per_cell) {
}

void S2ShapeIndexOptions::set_max_edges_per_cell(
    int max_edges_per_cell) {
  max_edges_per_cell_ = max_edges_per_cell;
}

S2Point S2ShapeIndex::Iterator::center() const {
  return id().ToPoint();
}

bool S2ShapeIndex::Iterator::Locate(S2Point const& target_point) {
  // Let I = cell_map_->lower_bound(T), where T is the leaf cell containing
  // "target_point".  Then if T is contained by an index cell, then the
  // containing cell is either I or I'.  We test for containment by comparing
  // the ranges of leaf cells spanned by T, I, and I'.

  S2CellId target = S2CellId::FromPoint(target_point);
  Seek(target);
  if (!Done() && id().range_min() <= target) return true;
  if (!AtBegin()) {
    Prev();
    if (id().range_max() >= target) return true;
  }
  return false;
}

S2ShapeIndex::CellRelation S2ShapeIndex::Iterator::Locate(S2CellId target) {
  // Let T be the target, let I = cell_map_->lower_bound(T.range_min()), and
  // let I' be the predecessor of I.  If T contains any index cells, then T
  // contains I.  Similarly, if T is contained by an index cell, then the
  // containing cell is either I or I'.  We test for containment by comparing
  // the ranges of leaf cells spanned by T, I, and I'.

  Seek(target.range_min());
  if (!Done()) {
    if (id() >= target && id().range_min() <= target) return INDEXED;
    if (id() <= target.range_max()) return SUBDIVIDED;
  }
  if (!AtBegin()) {
    Prev();
    if (id().range_max() >= target) return INDEXED;
  }
  return DISJOINT;
}

S2ShapeIndex::S2ShapeIndex()
    : pending_additions_begin_(0),
      index_status_(FRESH) {
}

S2ShapeIndex::S2ShapeIndex(S2ShapeIndexOptions const& options)
    : options_(options),
      pending_additions_begin_(0),
      index_status_(FRESH) {
}

void S2ShapeIndex::Init(S2ShapeIndexOptions const& options) {
  DCHECK(shapes_.empty());
  options_ = options;
}

S2ShapeIndex::~S2ShapeIndex() {
  Reset();
}

void S2ShapeIndex::Add(S2Shape* shape) {
  // Additions are processed lazily by ApplyUpdates().
  shape->id_ = shapes_.size();
  shapes_.push_back(shape);
  base::subtle::NoBarrier_Store(&index_status_, STALE);
}

void S2ShapeIndex::Remove(S2Shape* shape) {
  // This class updates itself lazily, because it is much more efficient to
  // process additions and removals in batches.  However this means that when
  // a shape is removed, we need to make a copy of all its edges, since the
  // client is free to delete "shape" once this call is finished.

  int32 shape_id = shape->id();
  DCHECK(shapes_[shape_id] != nullptr);
  shapes_[shape_id] = nullptr;
  if (shape_id >= pending_additions_begin_) {
    // We are removing a shape that has not yet been added to the index,
    // so there is nothing else to do.
  } else {
    if (!pending_removals_) {
      pending_removals_.reset(new vector<RemovedShape>);
    }
    // We build the new RemovedShape in place, since it includes a potentially
    // large vector of edges that might be expensive to copy.
    pending_removals_->push_back(RemovedShape());
    RemovedShape* removed = &pending_removals_->back();
    removed->shape_id = shape->id();
    removed->has_interior = shape->has_interior();
    removed->contains_origin = shape->contains_origin();
    int num_edges = shape->num_edges();
    removed->edges.reserve(num_edges);
    for (int e = 0; e < num_edges; ++e) {
      S2Point const *va, *vb;
      shape->GetEdge(e, &va, &vb);
      removed->edges.push_back(std::make_pair(*va, *vb));
    }
  }
  // Return ownership to the caller; do not delete the shape.
  base::subtle::NoBarrier_Store(&index_status_, STALE);
}

void S2ShapeIndex::RemoveAll() {
  // Note that vector::clear() does not actually free storage.
  vector<S2Shape*> empty_shapes;
  shapes_.swap(empty_shapes);

  Iterator it;
  for (it.InitStale(*this); !it.Done(); it.Next()) {
    delete it.cell();
  }
  cell_map_.clear();
  pending_additions_begin_ = 0;
  pending_removals_.reset();
  DCHECK(update_state_ == nullptr);
  base::subtle::NoBarrier_Store(&index_status_, FRESH);
}

void S2ShapeIndex::Reset() {
  for (S2Shape* shape : shapes_) {
    delete shape;
  }
  RemoveAll();
}

// Helper class for testing whether the point "p" is contained by one or more
// shapes that intersect a given index cell.
class PointContainmentTester {
 public:
  // "it" should be positioned at the cell to be tested.  "p" is the target.
  PointContainmentTester(S2ShapeIndex::Iterator const& it, S2Point const& p)
      : it_(it), point_(p), crosser_initialized_(false) {
  }
  // Return true if the given shape contains "p".  "clipped" should point to
  // the S2ClippedShape for "shape" in the given index cell.
  bool ContainedBy(S2Shape const* shape, S2ClippedShape const& clipped);
 private:
  S2ShapeIndex::Iterator const& it_;
  S2Point const& point_;
  bool crosser_initialized_;
  S2Point center_;
  S2EdgeUtil::EdgeCrosser crosser_;
};

bool PointContainmentTester::ContainedBy(S2Shape const* shape,
                                         S2ClippedShape const& clipped) {
  if (!shape->has_interior()) {
    return false;
  }
  bool inside = clipped.contains_center();
  int num_clipped = clipped.num_edges();
  if (num_clipped == 0) {
    return inside;
  }
  // We initialize the EdgeCrosser lazily.  This saves work when the cell
  // does not intersect any edges (which happens with interior cells).
  if (!crosser_initialized_) {
    center_ = it_.center();
    crosser_.Init(&center_, &point_);
    crosser_initialized_ = true;
  }
  // Test containment by drawing a line segment from the cell center to the
  // given point and counting edge crossings.
  for (int i = 0; i < num_clipped; ++i) {
    S2Point const *a0, *a1;
    shape->GetEdge(clipped.edge(i), &a0, &a1);
    inside ^= crosser_.EdgeOrVertexCrossing(a0, a1);
  }
  return inside;
}

bool S2ShapeIndex::ShapeContains(S2Shape const* shape, S2Point const& p) const {
  // Look up the S2ShapeIndex cell containing this point.
  S2ShapeIndex::Iterator it(*this);
  if (!it.Locate(p)) return false;

  S2ClippedShape const* clipped = it.cell()->find_clipped(shape);
  if (clipped == nullptr) return false;

  PointContainmentTester point_tester(it, p);
  return point_tester.ContainedBy(shape, *clipped);
}

bool S2ShapeIndex::GetContainingShapes(S2Point const& p,
                                       std::vector<S2Shape*>* shapes) const {
  // Look up the S2ShapeIndex cell containing this point.
  shapes->clear();
  S2ShapeIndex::Iterator it(*this);
  if (!it.Locate(p)) return false;

  int num_shapes = it.cell()->num_shapes();
  if (num_shapes == 0) return false;

  PointContainmentTester point_tester(it, p);
  for (int s = 0; s < num_shapes; ++s) {
    S2ClippedShape const& clipped = it.cell()->clipped(s);
    S2Shape* shape = shapes_[clipped.shape_id()];
    if (point_tester.ContainedBy(shape, clipped)) {
      shapes->push_back(shape);
    }
  }
  return !shapes->empty();
}

// FaceEdge and ClippedEdge store temporary edge data while the index is being
// updated.  FaceEdge represents an edge that has been projected onto a given
// face, while ClippedEdge represents the portion of that edge that has been
// clipped to a given S2Cell.
//
// While it would be possible to combine all the edge information into one
// structure, there are two good reasons for separating it:
//
//  - Memory usage.  Separating the two classes means that we only need to
//    store one copy of the per-face data no matter how many times an edge is
//    subdivided, and it also lets us delay computing bounding boxes until
//    they are needed for processing each face (when the dataset spans
//    multiple faces).
//
//  - Performance.  UpdateEdges is significantly faster on large polygons when
//    the data is separated, because it often only needs to access the data in
//    ClippedEdge and this data is cached more successfully.

struct S2ShapeIndex::FaceEdge {
  int32 shape_id;     // The shape that this edge belongs to
  int32 edge_id;      // Edge id within that shape
  int32 max_level;    // Not desirable to subdivide this edge beyond this level
  bool has_interior;  // Belongs to a shape that has an interior
  R2Point a, b;       // The edge endpoints, clipped to a given face
  S2Point const* va;  // The original S2Loop vertices
  S2Point const* vb;
};

struct S2ShapeIndex::ClippedEdge {
  FaceEdge const* face_edge;  // The original unclipped edge
  R2Rect bound;               // Bounding box for the clipped portion
};

// Given a set of shapes, InteriorTracker keeps track of which shapes contain
// a particular point (the "focus").  It provides an efficient way to move the
// focus from one point to another and incrementally update the set of shapes
// which contain it.  We use this to compute which shapes contain the center
// of every S2CellId in the index, by advancing the focus from one cell center
// to the next.
//
// Initially the focus is S2::Origin(), and therefore we can initialize the
// state of every shape to its contains_origin() value.  Next we advance the
// focus to the start of the S2CellId space-filling curve, by drawing a line
// segment between this point and S2::Origin() and testing whether every edge
// of every shape intersects it.  Then we visit all the cells that are being
// added to the S2ShapeIndex in increasing order of S2CellId.  For each cell,
// we draw two edges: one from the entry vertex to the center, and another
// from the center to the exit vertex (where "entry" and "exit" refer to the
// points where the space-filling curve enters and exits the cell).  By
// counting edge crossings we can incrementally compute which shapes contain
// the cell center.  Note that the same set of shapes will always contain the
// exit point of one cell and the entry point of the next cell in the index,
// because either (a) these two points are actually the same, or (b) the
// intervening cells in S2CellId order are all empty, and therefore there are
// no edge crossings if we follow this path from one cell to the other.
class S2ShapeIndex::InteriorTracker {
 public:
  // Initialize the InteriorTracker.  You must call AddShape() for each shape
  // that will be tracked before calling MoveTo() or DrawTo().
  InteriorTracker();

  // Return true if any shapes are being tracked.
  bool is_active() const { return is_active_; }

  // Add a shape whose interior should be tracked.  "is_inside" indicates
  // whether the current focus point is inside the shape.  Alternatively, if
  // the focus point is in the process of being moved (via MoveTo/DrawTo), you
  // can also specify "is_inside" at the old focus point and call TestEdge()
  // for every edge of the shape that might cross the current DrawTo() line.
  // This updates the state to correspond to the new focus point.
  //
  // REQUIRES: shape->has_interior()
  void AddShape(int32 shape_id, bool is_inside);

  // Move the focus to the given point.  This method should only be used when
  // it is known that there are no edge crossings between the old and new
  // focus locations; otherwise use DrawTo().
  void MoveTo(S2Point const& b) { b_ = b; }

  // Move the focus to the given point.  After this method is called,
  // TestEdge() should be called with all edges that may cross the line
  // segment between the old and new focus locations.
  void DrawTo(S2Point const& b);

  // Indicate that the given edge of the given shape may cross the line
  // segment between the old and new focus locations (see DrawTo).
  // REQUIRES: shape->has_interior()
  inline void TestEdge(int32 shape_id,
                       S2Point const* c, S2Point const* d);

  // The set of shape ids that contain the current focus.
  ShapeIdSet const& shape_ids() const { return shape_ids_; }

  // Indicate that the last argument to MoveTo() or DrawTo() was the entry
  // vertex of the given S2CellId, i.e. the tracker is positioned at the start
  // of this cell.  By using this method together with at_cellid(), the caller
  // can avoid calling MoveTo() in cases where the exit vertex of the previous
  // cell is the same as the entry vertex of the current cell.
  void set_next_cellid(S2CellId next_cellid) {
    next_cellid_ = next_cellid.range_min();
  }

  // Return true if the focus is already at the entry vertex of the given
  // S2CellId (provided that the caller calls set_next_cellid() as each cell
  // is processed).
  bool at_cellid(S2CellId cellid) const {
    return cellid.range_min() == next_cellid_;
  }

  // Make an internal copy of the state for shape ids below the given limit,
  // and then clear the state for those shapes.  This is used during
  // incremental updates to track the state of added and removed shapes
  // separately.
  void SaveAndClearStateBefore(int32 limit_shape_id);

  // Restore the state previously saved by SaveAndClearStateBefore().  This
  // only affects the state for shape_ids below "limit_shape_id".
  void RestoreStateBefore(int32 limit_shape_id);

 private:
  // Remove "shape_id" from shape_ids_ if it exists, otherwise insert it.
  void ToggleShape(int shape_id);

  // Return a pointer to the first entry "x" where x >= shape_id.
  ShapeIdSet::iterator lower_bound(int32 shape_id);

  bool is_active_;
  S2Point a_, b_;
  S2CellId next_cellid_;
  S2EdgeUtil::EdgeCrosser crosser_;
  ShapeIdSet shape_ids_;

  // Shape ids saved by SaveAndClearStateBefore().  The state is never saved
  // recursively so we don't need to worry about maintaining a stack.
  ShapeIdSet saved_ids_;
};

// As shapes are added, we compute which ones contain the start of the
// S2CellId space-filling curve by drawing an edge from S2::Origin() to this
// point and counting how many shape edges cross this edge.
S2ShapeIndex::InteriorTracker::InteriorTracker()
    : is_active_(false), b_(S2::Origin()),
      next_cellid_(S2CellId::Begin(S2CellId::kMaxLevel)) {
  DrawTo(S2::FaceUVtoXYZ(0, -1, -1).Normalize());  // S2CellId curve start
}

void S2ShapeIndex::InteriorTracker::AddShape(int32 shape_id,
                                             bool contains_origin) {
  is_active_ = true;
  if (contains_origin) {
    ToggleShape(shape_id);
  }
}

void S2ShapeIndex::InteriorTracker::ToggleShape(int shape_id) {
  // Since shape_ids_.size() is typically *very* small (0, 1, or 2), it turns
  // out to be significantly faster to maintain a sorted array rather than
  // using an STL set or btree_set.
  if (shape_ids_.empty()) {
    shape_ids_.push_back(shape_id);
  } else if (shape_ids_[0] == shape_id) {
    shape_ids_.erase(shape_ids_.begin());
  } else {
    ShapeIdSet::iterator pos = shape_ids_.begin();
    while (*pos < shape_id) {
      if (++pos == shape_ids_.end()) {
        shape_ids_.push_back(shape_id);
        return;
      }
    }
    if (*pos == shape_id) {
      shape_ids_.erase(pos);
    } else {
      shape_ids_.insert(pos, shape_id);
    }
  }
}

void S2ShapeIndex::InteriorTracker::DrawTo(S2Point const& b) {
  a_ = b_;
  b_ = b;
  crosser_.Init(&a_, &b_);
}

inline void S2ShapeIndex::InteriorTracker::TestEdge(int32 shape_id,
                                                    S2Point const* c,
                                                    S2Point const* d) {
  if (crosser_.EdgeOrVertexCrossing(c, d)) {
    ToggleShape(shape_id);
  }
}

// Like std::lower_bound(shape_ids_.begin(), shape_ids_.end(), shape_id), but
// implemented with linear rather than binary search because the number of
// shapes being tracked is typically very small.
inline S2ShapeIndex::ShapeIdSet::iterator
S2ShapeIndex::InteriorTracker::lower_bound(int32 shape_id) {
  ShapeIdSet::iterator pos = shape_ids_.begin();
  while (pos != shape_ids_.end() && *pos < shape_id) { ++pos; }
  return pos;
}

void S2ShapeIndex::InteriorTracker::SaveAndClearStateBefore(
    int32 limit_shape_id) {
  DCHECK(saved_ids_.empty());
  ShapeIdSet::iterator limit = lower_bound(limit_shape_id);
  saved_ids_.assign(shape_ids_.begin(), limit);
  shape_ids_.erase(shape_ids_.begin(), limit);
}

void S2ShapeIndex::InteriorTracker::RestoreStateBefore(int32 limit_shape_id) {
  shape_ids_.erase(shape_ids_.begin(), lower_bound(limit_shape_id));
  shape_ids_.insert(shape_ids_.begin(), saved_ids_.begin(), saved_ids_.end());
  saved_ids_.clear();
}

// Apply any pending updates in a thread-safe way.
void S2ShapeIndex::ApplyUpdatesThreadSafe() {
  lock_.Lock();
  if (base::subtle::NoBarrier_Load(&index_status_) == FRESH) {
    lock_.Unlock();
  } else if (base::subtle::NoBarrier_Load(&index_status_) == UPDATING) {
    // Wait until the updating thread is finished.  We do this by attempting
    // to lock a mutex that is held by the updating thread.  When this mutex
    // is unlocked the index_status_ is guaranteed to be FRESH.
    ++update_state_->num_waiting;
    lock_.Unlock();
    update_state_->wait_mutex.Lock();
    lock_.Lock();
    --update_state_->num_waiting;
    UnlockAndSignal();  // Notify other waiting threads.
  } else {
    DCHECK_EQ(STALE, index_status_);
    base::subtle::NoBarrier_Store(&index_status_, UPDATING);
    // Allocate the extra state needed for thread synchronization.  We keep
    // the spinlock held while doing this, because (1) memory allocation is
    // fast, so the chance of a context switch while holding the lock is low;
    // (2) by far the most common situation is that there is no contention,
    // and this saves an extra lock and unlock step; (3) even in the rare case
    // where there is contention, the main side effect is that some other
    // thread will burn a few CPU cycles rather than sleeping.
    update_state_.reset(new UpdateState);
    // lock_.Lock wait_mutex *before* calling Unlock() to ensure that all other
    // threads will block on it.
    update_state_->wait_mutex.Lock();
    // Release the spinlock before doing any real work.
    lock_.Unlock();
    ApplyUpdatesInternal();
    lock_.Lock();
    // index_status_ can be updated to FRESH only while locked *and* using
    // an atomic store operation, so that MaybeApplyUpdates() can check
    // whether the index is FRESH without acquiring the spinlock.
    base::subtle::Release_Store(&index_status_, FRESH);
    UnlockAndSignal();  // Notify any waiting threads.
  }
}

// Releases lock_ and wakes up any waiting threads by releasing wait_mutex.
// If this was the last waiting thread, also deletes update_state_.
// REQUIRES: lock_ is held.
// REQUIRES: wait_mutex is held.
inline void S2ShapeIndex::UnlockAndSignal() {
  DCHECK_EQ(FRESH, index_status_);
  int num_waiting = update_state_->num_waiting;
  lock_.Unlock();
  // Allow another waiting thread to proceed.  Note that no new threads can
  // start waiting because the index_status_ is now FRESH, and the caller is
  // required to prevent any new mutations from occurring while these const
  // methods are running.
  //
  // We need to unlock wait_mutex before destroying it even if there are no
  // waiting threads.
  update_state_->wait_mutex.Unlock();
  if (num_waiting == 0) {
    update_state_.reset();
  }
}

void S2ShapeIndex::ForceApplyUpdates() {
  // No locks required because this is not a const method.  It is the client's
  // responsibility to ensure correct thread synchronization.
  if (base::subtle::NoBarrier_Load(&index_status_) != FRESH) {
    ApplyUpdatesInternal();
    base::subtle::NoBarrier_Store(&index_status_, FRESH);
  }
}

// A BatchDescriptor represents a set of pending updates that will be applied
// at the same time.  The batch consists of all updates with shape ids between
// the current value of "ShapeIndex::pending_additions_begin_" (inclusive) and
// "additions_end" (exclusive).  The first batch to be processed also
// implicitly includes all shapes being removed.  "num_edges" is the total
// number of edges that will be added or removed in this batch.
struct S2ShapeIndex::BatchDescriptor {
  BatchDescriptor(int _additions_end, int _num_edges)
      : additions_end(_additions_end), num_edges(_num_edges) {
  }
  int additions_end;
  int num_edges;
};

// This method updates the index by applying all pending additions and
// removals.  It does *not* update index_status_ (see ApplyUpdatesThreadSafe).
void S2ShapeIndex::ApplyUpdatesInternal() {
  // Check whether we have so many edges to process that we should process
  // them in multiple batches to save memory.  Building the index can use up
  // to 20x as much memory (per edge) as the final index size.
  vector<BatchDescriptor> batches;
  GetUpdateBatches(&batches);
  for (int i = 0; i < batches.size(); ++i) {
    vector<FaceEdge> all_edges[6];
    BatchDescriptor const& batch = batches[i];
    VLOG(1) << "Batch " << i << ": shape_limit=" << batch.additions_end
               << ", edges=" << batch.num_edges;

    ReserveSpace(batch, all_edges);
    InteriorTracker tracker;
    if (pending_removals_) {
      // The first batch implicitly includes all shapes being removed.
      for (int j = 0; j < pending_removals_->size(); ++j) {
        RemoveShape((*pending_removals_)[j], all_edges, &tracker);
      }
    }
    for (int id = pending_additions_begin_; id < batch.additions_end; ++id) {
      AddShape(id, all_edges, &tracker);
    }
    for (int face = 0; face < 6; ++face) {
      UpdateFaceEdges(face, all_edges[face], &tracker);
      // Save memory by clearing vectors after we are done with them.
      vector<FaceEdge> empty;
      all_edges[face].swap(empty);
    }
    // We can't clear pending_removals_ until all updates are processed,
    // because FaceEdge stores pointers directly to the removed vertices.
    pending_removals_.reset(nullptr);
    pending_additions_begin_ = batch.additions_end;
  }
  // It is the caller's responsibility to update index_status_.
}

// Count the number of edges being updated, and break them into several
// batches if necessary to reduce the amount of memory needed.  (See the
// documentation for FLAGS_s2shapeindex_tmp_memory_budget_mb.)
void S2ShapeIndex::GetUpdateBatches(vector<BatchDescriptor>* batches) const {
  // Count the edges being removed and added.
  int num_edges_removed = 0;
  if (pending_removals_) {
    for (int i = 0; i < pending_removals_->size(); ++i) {
      num_edges_removed += (*pending_removals_)[i].edges.size();
    }
  }
  int num_edges_added = 0;
  for (int id = pending_additions_begin_; id < shapes_.size(); ++id) {
    S2Shape const* shape = shapes_[id];
    if (shape == nullptr) continue;
    num_edges_added += shape->num_edges();
  }
  int num_edges = num_edges_removed + num_edges_added;

  // The following memory estimates are based on heap profiling.
  //
  // The final size of an S2ShapeIndex depends mainly on how finely the index
  // is subdivided, as controlled by S2ShapeIndexOptions::max_edges_per_cell()
  // and --s2shapeindex_default_max_edges_per_cell. For realistic values of
  // max_edges_per_cell() and shapes with moderate numbers of edges, it is
  // difficult to get much below 8 bytes per edge.  [The minimum possible size
  // is 4 bytes per edge (to store a 32-bit edge id in an S2ClippedShape) plus
  // 24 bytes per shape (for the S2ClippedShape itself plus a pointer in the
  // shapes_ vector.]
  //
  // The temporary memory consists mainly of the FaceEdge and ClippedEdge
  // structures plus a ClippedEdge pointer for every level of recursive
  // subdivision.  For very large indexes this can be 200 bytes per edge.
  size_t const kFinalBytesPerEdge = 8;
  size_t const kTmpBytesPerEdge = 200;
  size_t const kTmpMemoryBudgetBytes =
      static_cast<size_t>(FLAGS_s2shapeindex_tmp_memory_budget_mb) << 20;

  // We arbitrarily limit the number of batches just as a safety measure.
  // With the current default memory budget of 100 MB, this limit is not
  // reached even when building an index of 350 million edges.
  int const kMaxUpdateBatches = 100;

  if (num_edges * kTmpBytesPerEdge <= kTmpMemoryBudgetBytes) {
    // We can update all edges at once without exceeding kTmpMemoryBudgetBytes.
    batches->push_back(BatchDescriptor(shapes_.size(), num_edges));
    return;
  }
  // Otherwise, break the updates into up to several batches, where the size
  // of each batch is chosen so that all batches use approximately the same
  // high-water memory.  GetBatchSizes() returns the recommended number of
  // edges in each batch.
  vector<int> batch_sizes;
  GetBatchSizes(num_edges, kMaxUpdateBatches, kFinalBytesPerEdge,
                kTmpBytesPerEdge, kTmpMemoryBudgetBytes, &batch_sizes);

  // We always process removed edges in a single batch, since (1) they already
  // take up a lot of memory because we have copied all their edges, and (2)
  // AbsorbIndexCell() uses (shapes_[id] == nullptr) to detect when a shape is
  // being removed, so in order to split the removals into batches we would
  // need a different approach (e.g., temporarily add fake entries to shapes_
  // and restore them back to nullptr as shapes are actually removed).
  num_edges = 0;
  if (pending_removals_) {
    num_edges += num_edges_removed;
    if (num_edges >= batch_sizes[0]) {
      batches->push_back(BatchDescriptor(pending_additions_begin_, num_edges));
      num_edges = 0;
    }
  }
  // Keep adding shapes to each batch until the recommended number of edges
  // for that batch is reached, then move on to the next batch.
  for (int id = pending_additions_begin_; id < shapes_.size(); ++id) {
    S2Shape const* shape = shapes_[id];
    if (shape == nullptr) continue;
    num_edges += shape->num_edges();
    if (num_edges >= batch_sizes[batches->size()]) {
      batches->push_back(BatchDescriptor(id + 1, num_edges));
      num_edges = 0;
    }
  }
  // Some shapes have no edges.  If a shape with no edges is the last shape to
  // be added or removed, then the final batch may not include it, so we fix
  // that problem here.
  batches->back().additions_end = shapes_.size();
  DCHECK_LE(batches->size(), kMaxUpdateBatches);
}

// Given "num_items" items, each of which uses "tmp_bytes_per_item" while it
// is being updated but only "final_bytes_per_item" in the end, divide the
// items into batches that have approximately the same *total* memory usage
// consisting of the temporary memory needed for the items in the current
// batch plus the final size of all the items that have already been
// processed.  Use the fewest number of batches (but never more than
// "max_batches") such that the total memory usage does not exceed the
// combined final size of all the items plus "tmp_memory_budget_bytes".
/* static */
void S2ShapeIndex::GetBatchSizes(int num_items, int max_batches,
                                 double final_bytes_per_item,
                                 double tmp_bytes_per_item,
                                 double tmp_memory_budget_bytes,
                                 vector<int>* batch_sizes) {
  // This code tries to fit all the data into the same memory space
  // ("total_budget_bytes") at every iteration.  The data consists of some
  // number of processed items (at "final_bytes_per_item" each), plus some
  // number being updated (at "tmp_bytes_per_item" each).  The space occupied
  // by the items being updated is the "free space".  At each iteration, the
  // free space is multiplied by (1 - final_bytes_per_item/tmp_bytes_per_item)
  // as the items are converted into their final form.
  double final_bytes = num_items * final_bytes_per_item;
  double final_bytes_ratio = final_bytes_per_item / tmp_bytes_per_item;
  double free_space_multiplier = 1 - final_bytes_ratio;

  // The total memory budget is the greater of the final size plus the allowed
  // temporary memory, or the minimum amount of memory required to limit the
  // number of batches to "max_batches".
  double total_budget_bytes = max(
      final_bytes + tmp_memory_budget_bytes,
      final_bytes / (1 - pow(free_space_multiplier, max_batches)));

  // "max_batch_items" is the number of items in the current batch.
  double max_batch_items = total_budget_bytes / tmp_bytes_per_item;
  batch_sizes->clear();
  for (int i = 0; i + 1 < max_batches && num_items > 0; ++i) {
    int batch_items =
        std::min(num_items, static_cast<int>(max_batch_items + 1));
    batch_sizes->push_back(batch_items);
    num_items -= batch_items;
    max_batch_items *= free_space_multiplier;
  }
  DCHECK_LE(batch_sizes->size(), max_batches);
}

// Reserve an appropriate amount of space for the top-level face edges in the
// current batch.  This data structure uses about half of the temporary memory
// needed during index construction.  Furthermore, if the arrays are grown via
// push_back() then up to 10% of the total run time consists of copying data
// as these arrays grow, so it is worthwhile to preallocate space for them.
void S2ShapeIndex::ReserveSpace(BatchDescriptor const& batch,
                                vector<FaceEdge> all_edges[6]) const {
  // If the number of edges is relatively small, then the fastest approach is
  // to simply reserve space on every face for the maximum possible number of
  // edges.  We use a different threshold for this calculation than for
  // deciding when to break updates into batches, because the cost/benefit
  // ratio is different.  (Here the only extra expense is that we need to
  // sample the edges to estimate how many edges per face there are.)
  size_t const kMaxCheapBytes = 30 << 20;  // 30 MB
  int const kMaxCheapEdges = kMaxCheapBytes / (6 * sizeof(FaceEdge));
  if (batch.num_edges <= kMaxCheapEdges) {
    for (int face = 0; face < 6; ++face) {
      all_edges[face].reserve(batch.num_edges);
    }
    return;
  }
  // Otherwise we estimate the number of edges on each face by taking a random
  // sample.  The goal is to come up with an estimate that is fast and
  // accurate for non-pathological geometry.  If our estimates happen to be
  // wrong, the vector will still grow automatically - the main side effects
  // are that memory usage will be larger (by up to a factor of 3), and
  // constructing the index will be about 10% slower.
  //
  // Given a desired sample size, we choose equally spaced edges from
  // throughout the entire data set.  We use a Bresenham-type algorithm to
  // choose the samples.
  int const kDesiredSampleSize = 10000;
  int const sample_interval = max(1, batch.num_edges / kDesiredSampleSize);

  // Initialize "edge_id" to be midway through the first sample interval.
  // Because samples are equally spaced the actual sample size may differ
  // slightly from the desired sample size.
  int edge_id = sample_interval / 2;
  int const actual_sample_size = (batch.num_edges + edge_id) / sample_interval;
  int face_count[6] = { 0, 0, 0, 0, 0, 0 };
  if (pending_removals_) {
    for (int i = 0; i < pending_removals_->size(); ++i) {
      RemovedShape const& removed = (*pending_removals_)[i];
      edge_id += removed.edges.size();
      while (edge_id >= sample_interval) {
        edge_id -= sample_interval;
        face_count[S2::GetFace(removed.edges[edge_id].first)] += 1;
      }
    }
  }
  for (int id = pending_additions_begin_; id < batch.additions_end; ++id) {
    S2Shape const* shape = shapes_[id];
    if (shape == nullptr) continue;
    edge_id += shape->num_edges();
    while (edge_id >= sample_interval) {
      edge_id -= sample_interval;
      S2Point const *a, *b;
      shape->GetEdge(edge_id, &a, &b);
      // For speed, we only count the face containing one endpoint of the
      // edge.  In general the edge could span all 6 faces (with padding), but
      // it's not worth the expense to compute this more accurately.
      face_count[S2::GetFace(*a)] += 1;
    }
  }
  // Now given the raw face counts, compute a confidence interval such that we
  // will be unlikely to allocate too little space.  Computing accurate
  // binomial confidence intervals is expensive and not really necessary.
  // Instead we use a simple approximation:
  //  - For any face with at least 1 sample, we use at least a 4-sigma
  //    confidence interval.  (The chosen width is adequate for the worst case
  //    accuracy, which occurs when the face contains approximately 50% of the
  //    edges.)  Assuming that our sample is representative, the probability of
  //    reserving too little space is approximately 1 in 30,000.
  //  - For faces with no samples at all, we don't bother reserving space.
  //    It is quite likely that such faces are truly empty, so we save time
  //    and memory this way.  If the face does contain some edges, there will
  //    only be a few so it is fine to let the vector grow automatically.
  // On average, we reserve 2% extra space for each face that has geometry.

  // kMaxSemiWidth is the maximum semi-width over all probabilities p of a
  // 4-sigma binomial confidence interval with a sample size of 10,000.
  double const kMaxSemiWidth = 0.02;
  double const sample_ratio = 1.0 / actual_sample_size;
  for (int face = 0; face < 6; ++face) {
    if (face_count[face] == 0) continue;
    double fraction = sample_ratio * face_count[face] + kMaxSemiWidth;
    all_edges[face].reserve(fraction * batch.num_edges);
  }
}

// Clip all edges of the given shape to the six cube faces, add the clipped
// edges to "all_edges", and start tracking its interior if necessary.
void S2ShapeIndex::AddShape(int id, vector<FaceEdge> all_edges[6],
                            InteriorTracker* tracker) const {
  S2Shape const* shape = shapes_[id];
  if (shape == nullptr) {
    return;  // This shape has already been removed.
  }
  FaceEdge edge;
  edge.shape_id = id;
  edge.has_interior = shape->has_interior();
  if (edge.has_interior) tracker->AddShape(id, shape->contains_origin());
  int num_edges = shape->num_edges();
  for (int e = 0; e < num_edges; ++e) {
    edge.edge_id = e;
    shape->GetEdge(e, &edge.va, &edge.vb);
    edge.max_level = GetEdgeMaxLevel(*edge.va, *edge.vb);
    if (edge.has_interior) tracker->TestEdge(id, edge.va, edge.vb);
    AddFaceEdge(&edge, all_edges);
  }
}

void S2ShapeIndex::RemoveShape(RemovedShape const& removed,
                               vector<FaceEdge> all_edges[6],
                               InteriorTracker* tracker) const {
  FaceEdge edge;
  edge.edge_id = -1;  // Not used or needed for removed edges.
  edge.shape_id = removed.shape_id;
  edge.has_interior = removed.has_interior;
  if (edge.has_interior) {
    tracker->AddShape(edge.shape_id, removed.contains_origin);
  }
  for (int j = 0; j < removed.edges.size(); ++j) {
    edge.va = &removed.edges[j].first;
    edge.vb = &removed.edges[j].second;
    edge.max_level = GetEdgeMaxLevel(*edge.va, *edge.vb);
    if (edge.has_interior) {
      tracker->TestEdge(edge.shape_id, edge.va, edge.vb);
    }
    AddFaceEdge(&edge, all_edges);
  }
}

inline void S2ShapeIndex::AddFaceEdge(FaceEdge* edge,
                                      vector<FaceEdge> all_edges[6]) const {
  // Fast path: both endpoints are on the same face, and are far enough from
  // the edge of the face that don't intersect any (padded) adjacent face.
  int a_face = S2::GetFace(*edge->va);
  if (a_face == S2::GetFace(*edge->vb)) {
    S2::ValidFaceXYZtoUV(a_face, *edge->va, &edge->a);
    S2::ValidFaceXYZtoUV(a_face, *edge->vb, &edge->b);
    double const kMaxUV = 1 - kCellPadding;
    if (fabs(edge->a[0]) <= kMaxUV && fabs(edge->a[1]) <= kMaxUV &&
        fabs(edge->b[0]) <= kMaxUV && fabs(edge->b[1]) <= kMaxUV) {
      all_edges[a_face].push_back(*edge);
      return;
    }
  }
  // Otherwise we simply clip the edge to all six faces.
  for (int face = 0; face < 6; ++face) {
    if (S2EdgeUtil::ClipToPaddedFace(*edge->va, *edge->vb, face, kCellPadding,
                                     &edge->a, &edge->b)) {
      all_edges[face].push_back(*edge);
    }
  }
}

// Return the first level at which the edge will *not* contribute towards
// the decision to subdivide.
int S2ShapeIndex::GetEdgeMaxLevel(S2Point const& a, S2Point const& b) const {
  // Compute the maximum cell size for which this edge is considered "long".
  // The calculation does not need to be perfectly accurate, so we use Norm()
  // rather than Angle() for speed.
  double cell_size = ((a - b).Norm() *
                      FLAGS_s2shapeindex_cell_size_to_long_edge_ratio);
  // Now return the first level encountered during subdivision where the
  // average cell size is at most "cell_size".
  return S2::kAvgEdge.GetMinLevel(cell_size);
}

// EdgeAllocator provides temporary storage for new ClippedEdges that are
// created during indexing.  It is essentially a stack model, where edges are
// allocated as the recursion does down and freed as it comes back up.
//
// It also provides a mutable vector of FaceEdges that is used when
// incrementally updating the index (see AbsorbIndexCell).
class S2ShapeIndex::EdgeAllocator {
 public:
  EdgeAllocator() : size_(0) {}

  ~EdgeAllocator() {
    STLDeleteElements(&clipped_edges_);
  }

  // Return a pointer to a newly allocated edge.
  ClippedEdge* NewClippedEdge() {
    if (size_ == clipped_edges_.size()) {
      clipped_edges_.push_back(new ClippedEdge);
    }
    return clipped_edges_[size_++];
  }
  // Return the number of allocated edges.
  size_t size() const { return size_; }

  // Reset the allocator to only contain the first "size" allocated edges.
  void Reset(size_t size) { size_ = size; }

  vector<FaceEdge>* mutable_face_edges() {
    return &face_edges_;
  }

 private:
  // We can't use vector<ClippedEdge> because edges are not allowed to move
  // once they have been allocated.  Instead we keep a pool of allocated edges
  // that are all deleted together at the end.
  size_t size_;
  vector<ClippedEdge*> clipped_edges_;

  // On the other hand, we can use vector<FaceEdge> because they are allocated
  // only at one level during the recursion (namely, the level at which we
  // absorb an existing index cell).
  vector<FaceEdge> face_edges_;

  EdgeAllocator(EdgeAllocator const&) = delete;
  void operator=(EdgeAllocator const&) = delete;
};

// Given a face and a vector of edges that intersect that face, add or remove
// all the edges from the index.  (An edge is added if shapes_[id] is not
// nullptr, and removed otherwise.)
void S2ShapeIndex::UpdateFaceEdges(int face,
                                   vector<FaceEdge> const& face_edges,
                                   InteriorTracker* tracker) {
  int num_edges = face_edges.size();
  if (num_edges == 0 && tracker->shape_ids().empty()) return;

  // Create the initial ClippedEdge for each FaceEdge.  Additional clipped
  // edges are created when edges are split between child cells.  We create
  // two arrays, one containing the edge data and another containing pointers
  // to those edges, so that during the recursion we only need to copy
  // pointers in order to propagate an edge to the correct child.
  vector<ClippedEdge> clipped_edge_storage;
  vector<ClippedEdge const*> clipped_edges;
  clipped_edge_storage.reserve(num_edges);
  clipped_edges.reserve(num_edges);
  R2Rect bound = R2Rect::Empty();
  for (int e = 0; e < num_edges; ++e) {
    ClippedEdge clipped;
    clipped.face_edge = &face_edges[e];
    clipped.bound = R2Rect::FromPointPair(face_edges[e].a, face_edges[e].b);
    clipped_edge_storage.push_back(clipped);
    clipped_edges.push_back(&clipped_edge_storage.back());
    bound.AddRect(clipped.bound);
  }
  // Construct the initial face cell containing all the edges, and then update
  // all the edges in the index recursively.
  EdgeAllocator alloc;
  S2CellId face_id = S2CellId::FromFace(face);
  S2PaddedCell pcell(face_id, kCellPadding);

  // "disjoint_from_index" means that the current cell being processed (and
  // all its descendants) are not already present in the index.
  bool disjoint_from_index = is_first_update();
  if (num_edges > 0) {
    S2CellId shrunk_id = ShrinkToFit(pcell, bound);
    if (shrunk_id != pcell.id()) {
      // All the edges are contained by some descendant of the face cell.  We
      // can save a lot of work by starting directly with that cell, but if we
      // are in the interior of at least one shape then we need to create
      // index entries for the cells we are skipping over.
      SkipCellRange(face_id.range_min(), shrunk_id.range_min(),
                    tracker, &alloc, disjoint_from_index);
      pcell = S2PaddedCell(shrunk_id, kCellPadding);
      UpdateEdges(pcell, &clipped_edges, tracker, &alloc, disjoint_from_index);
      SkipCellRange(shrunk_id.range_max().next(), face_id.range_max().next(),
                    tracker, &alloc, disjoint_from_index);
      return;
    }
  }
  // Otherwise (no edges, or no shrinking is possible), subdivide normally.
  UpdateEdges(pcell, &clipped_edges, tracker, &alloc, disjoint_from_index);
}

inline S2CellId S2ShapeIndex::ShrinkToFit(S2PaddedCell const& pcell,
                                          R2Rect const& bound) const {
  S2CellId shrunk_id = pcell.ShrinkToFit(bound);
  if (!is_first_update() && shrunk_id != pcell.id()) {
    // Don't shrink any smaller than the existing index cells, since we need
    // to combine the new edges with those cells.
    Iterator iter;
    iter.InitStale(*this);  // "Stale" avoids applying updates recursively.
    CellRelation r = iter.Locate(shrunk_id);
    if (r == INDEXED) { shrunk_id = iter.id(); }
  }
  return shrunk_id;
}

// Skip over the cells in the given range, creating index cells if we are
// currently in the interior of at least one shape.
void S2ShapeIndex::SkipCellRange(S2CellId begin, S2CellId end,
                                 InteriorTracker* tracker,
                                 EdgeAllocator* alloc,
                                 bool disjoint_from_index) {
  // If we aren't in the interior of a shape, then skipping over cells is easy.
  if (tracker->shape_ids().empty()) return;

  // Otherwise generate the list of cell ids that we need to visit, and create
  // an index entry for each one.
  S2CellUnion skipped;
  skipped.InitFromBeginEnd(begin, end);
  for (int i = 0; i < skipped.num_cells(); ++i) {
    vector<ClippedEdge const*> clipped_edges;
    UpdateEdges(S2PaddedCell(skipped.cell_id(i), kCellPadding),
                &clipped_edges, tracker, alloc, disjoint_from_index);
  }
}

// Given a cell and a set of ClippedEdges whose bounding boxes intersect that
// cell, add or remove all the edges from the index.  Temporary space for
// edges that need to be subdivided is allocated from the given EdgeAllocator.
// "disjoint_from_index" is an optimization hint indicating that cell_map_
// does not contain any entries that overlap the given cell.
void S2ShapeIndex::UpdateEdges(S2PaddedCell const& pcell,
                               vector<ClippedEdge const*>* edges,
                               InteriorTracker* tracker,
                               EdgeAllocator* alloc,
                               bool disjoint_from_index) {
  // Cases where an index cell is not needed should be detected before this.
  DCHECK(!edges->empty() || !tracker->shape_ids().empty());

  // This function is recursive with a maximum recursion depth of 30
  // (S2CellId::kMaxLevel).  Note that using an explicit stack does not seem
  // to be any faster based on profiling.

  // Incremental updates are handled as follows.  All edges being added or
  // removed are combined together in "edges", and all shapes with interiors
  // are tracked using "tracker".  We subdivide recursively as usual until we
  // encounter an existing index cell.  At this point we "absorb" the index
  // cell as follows:
  //
  //   - Edges and shapes that are being removed are deleted from "edges" and
  //     "tracker".
  //   - All remaining edges and shapes from the index cell are added to
  //     "edges" and "tracker".
  //   - Continue subdividing recursively, creating new index cells as needed.
  //   - When the recursion gets back to the cell that was absorbed, we
  //     restore "edges" and "tracker" to their previous state.
  //
  // Note that the only reason that we include removed shapes in the recursive
  // subdivision process is so that we can find all of the index cells that
  // contain those shapes efficiently, without maintaining an explicit list of
  // index cells for each shape (which would be expensive in terms of memory).
  bool index_cell_absorbed = false;
  if (!disjoint_from_index) {
    // There may be existing index cells contained inside "pcell".  If we
    // encounter such a cell, we need to combine the edges being updated with
    // the existing cell contents by "absorbing" the cell.
    Iterator iter;
    iter.InitStale(*this);  // "Stale" avoids applying updates recursively.
    CellRelation r = iter.Locate(pcell.id());
    if (r == DISJOINT) {
      disjoint_from_index = true;
    } else if (r == INDEXED) {
      // Absorb the index cell by transferring its contents to "edges" and
      // deleting it.  We also start tracking the interior of any new shapes.
      AbsorbIndexCell(pcell, iter, edges, tracker, alloc);
      index_cell_absorbed = true;
      disjoint_from_index = true;
    } else {
      DCHECK_EQ(SUBDIVIDED, r);
    }
  }

  // If there are existing index cells below us, then we need to keep
  // subdividing so that we can merge with those cells.  Otherwise,
  // MakeIndexCell checks if the number of edges is small enough, and creates
  // an index cell if possible (returning true when it does so).
  if (!disjoint_from_index || !MakeIndexCell(pcell, *edges, tracker)) {
    // Reserve space for the edges that will be passed to each child.  This is
    // important since otherwise the running time is dominated by the time
    // required to grow the vectors.  The amount of memory involved is
    // relatively small, so we simply reserve the maximum space for every child.
    vector<ClippedEdge const*> child_edges[2][2];  // [i][j]
    int num_edges = edges->size();
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        child_edges[i][j].reserve(num_edges);
      }
    }

    // Remember the current size of the EdgeAllocator so that we can free any
    // edges that are allocated during edge splitting.
    size_t alloc_size = alloc->size();

    // Compute the middle of the padded cell, defined as the rectangle in
    // (u,v)-space that belongs to all four (padded) children.  By comparing
    // against the four boundaries of "middle" we can determine which children
    // each edge needs to be propagated to.
    R2Rect const& middle = pcell.middle();

    // Build up a vector edges to be passed to each child cell.  The (i,j)
    // directions are left (i=0), right (i=1), lower (j=0), and upper (j=1).
    // Note that the vast majority of edges are propagated to a single child.
    // This case is very fast, consisting of between 2 and 4 floating-point
    // comparisons and copying one pointer.  (ClipVAxis is inline.)
    for (int e = 0; e < num_edges; ++e) {
      ClippedEdge const* edge = (*edges)[e];
      if (edge->bound[0].hi() <= middle[0].lo()) {
        // Edge is entirely contained in the two left children.
        ClipVAxis(edge, middle[1], child_edges[0], alloc);
      } else if (edge->bound[0].lo() >= middle[0].hi()) {
        // Edge is entirely contained in the two right children.
        ClipVAxis(edge, middle[1], child_edges[1], alloc);
      } else if (edge->bound[1].hi() <= middle[1].lo()) {
        // Edge is entirely contained in the two lower children.
        child_edges[0][0].push_back(ClipUBound(edge, 1, middle[0].hi(), alloc));
        child_edges[1][0].push_back(ClipUBound(edge, 0, middle[0].lo(), alloc));
      } else if (edge->bound[1].lo() >= middle[1].hi()) {
        // Edge is entirely contained in the two upper children.
        child_edges[0][1].push_back(ClipUBound(edge, 1, middle[0].hi(), alloc));
        child_edges[1][1].push_back(ClipUBound(edge, 0, middle[0].lo(), alloc));
      } else {
        // The edge bound spans all four children.  The edge itself intersects
        // either three or four (padded) children.
        ClippedEdge const* left = ClipUBound(edge, 1, middle[0].hi(), alloc);
        ClipVAxis(left, middle[1], child_edges[0], alloc);
        ClippedEdge const* right = ClipUBound(edge, 0, middle[0].lo(), alloc);
        ClipVAxis(right, middle[1], child_edges[1], alloc);
      }
    }
    // Free any memory reserved for children that turned out to be empty.  This
    // step is cheap and reduces peak memory usage by about 10% when building
    // large S2ShapeIndexes (> 10M edges).
    for (int i = 0; i < 2; ++i) {
      for (int j = 0; j < 2; ++j) {
        if (child_edges[i][j].empty()) {
          vector<ClippedEdge const*> empty;
          empty.swap(child_edges[i][j]);
        }
      }
    }
    // Now recursively update the edges in each child.  We call the children in
    // increasing order of S2CellId so that when the index is first constructed,
    // all insertions into cell_map_ are at the end (which is much faster).
    for (int pos = 0; pos < 4; ++pos) {
      int i, j;
      pcell.GetChildIJ(pos, &i, &j);
      if (!child_edges[i][j].empty() || !tracker->shape_ids().empty()) {
        UpdateEdges(S2PaddedCell(pcell, i, j), &child_edges[i][j],
                    tracker, alloc, disjoint_from_index);
      }
    }
    // Free any temporary edges that were allocated during clipping.
    alloc->Reset(alloc_size);
  }
  if (index_cell_absorbed) {
    // Restore the state for any edges being removed that we are tracking.
    tracker->RestoreStateBefore(pending_additions_begin_);
  }
}

// Given an edge and an interval "middle" along the v-axis, clip the edge
// against the boundaries of "middle" and add the edge to the corresponding
// children.
/* static */
inline void S2ShapeIndex::ClipVAxis(ClippedEdge const* edge,
                                    R1Interval const& middle,
                                    vector<ClippedEdge const*> child_edges[2],
                                    EdgeAllocator* alloc) {
  if (edge->bound[1].hi() <= middle.lo()) {
    // Edge is entirely contained in the lower child.
    child_edges[0].push_back(edge);
  } else if (edge->bound[1].lo() >= middle.hi()) {
    // Edge is entirely contained in the upper child.
    child_edges[1].push_back(edge);
  } else {
    // The edge bound spans both children.
    child_edges[0].push_back(ClipVBound(edge, 1, middle.hi(), alloc));
    child_edges[1].push_back(ClipVBound(edge, 0, middle.lo(), alloc));
  }
}

// Given an edge, clip the given endpoint (lo=0, hi=1) of the u-axis so that
// it does not extend past the given value.
/* static */
S2ShapeIndex::ClippedEdge const*
S2ShapeIndex::ClipUBound(ClippedEdge const* edge, int u_end, double u,
                         EdgeAllocator* alloc) {
  // First check whether the edge actually requires any clipping.  (Sometimes
  // this method is called when clipping is not necessary, e.g. when one edge
  // endpoint is in the overlap area between two padded child cells.)
  if (u_end == 0) {
    if (edge->bound[0].lo() >= u) return edge;
  } else {
    if (edge->bound[0].hi() <= u) return edge;
  }
  // We interpolate the new v-value from the endpoints of the original edge.
  // This has two advantages: (1) we don't need to store the clipped endpoints
  // at all, just their bounding box; and (2) it avoids the accumulation of
  // roundoff errors due to repeated interpolations.  The result needs to be
  // clamped to ensure that it is in the appropriate range.
  FaceEdge const& e = *edge->face_edge;
  double v = edge->bound[1].ClampPoint(
      S2EdgeUtil::InterpolateDouble(u, e.a[0], e.b[0], e.a[1], e.b[1]));

  // Determine which endpoint of the v-axis bound to update.  If the edge
  // slope is positive we update the same endpoint, otherwise we update the
  // opposite endpoint.
  int v_end = u_end ^ ((e.a[0] > e.b[0]) != (e.a[1] > e.b[1]));
  return UpdateBound(edge, u_end, u, v_end, v, alloc);
}

// Given an edge, clip the given endpoint (lo=0, hi=1) of the v-axis so that
// it does not extend past the given value.
/* static */
S2ShapeIndex::ClippedEdge const*
S2ShapeIndex::ClipVBound(ClippedEdge const* edge, int v_end, double v,
                         EdgeAllocator* alloc) {
  // See comments in ClipUBound.
  if (v_end == 0) {
    if (edge->bound[1].lo() >= v) return edge;
  } else {
    if (edge->bound[1].hi() <= v) return edge;
  }
  FaceEdge const& e = *edge->face_edge;
  double u = edge->bound[0].ClampPoint(
      S2EdgeUtil::InterpolateDouble(v, e.a[1], e.b[1], e.a[0], e.b[0]));
  int u_end = v_end ^ ((e.a[0] > e.b[0]) != (e.a[1] > e.b[1]));
  return UpdateBound(edge, u_end, u, v_end, v, alloc);
}

// Given an edge and two bound endpoints that need to be updated, allocate and
// return a new edge with the updated bound.
/* static */
inline S2ShapeIndex::ClippedEdge const*
S2ShapeIndex::UpdateBound(ClippedEdge const* edge, int u_end, double u,
                          int v_end, double v, EdgeAllocator* alloc) {
  ClippedEdge* clipped = alloc->NewClippedEdge();
  clipped->face_edge = edge->face_edge;
  clipped->bound[0][u_end] = u;
  clipped->bound[1][v_end] = v;
  clipped->bound[0][1-u_end] = edge->bound[0][1-u_end];
  clipped->bound[1][1-v_end] = edge->bound[1][1-v_end];
  DCHECK(!clipped->bound.is_empty());
  DCHECK(edge->bound.Contains(clipped->bound));
  return clipped;
}

// Absorb an index cell by transferring its contents to "edges" and/or
// "tracker", and then delete this cell from the index.  If "edges" includes
// any edges that are being removed, this method also updates their
// InteriorTracker state to correspond to the exit vertex of this cell, and
// saves the InteriorTracker state by calling SaveAndClearStateBefore().  It
// is the caller's responsibility to restore this state by calling
// RestoreStateBefore() when processing of this cell is finished.
void S2ShapeIndex::AbsorbIndexCell(S2PaddedCell const& pcell,
                                   Iterator const& iter,
                                   vector<ClippedEdge const*>* edges,
                                   InteriorTracker* tracker,
                                   EdgeAllocator* alloc) {
  DCHECK_EQ(pcell.id(), iter.id());

  // When we absorb a cell, we erase all the edges that are being removed.
  // However when we are finished with this cell, we want to restore the state
  // of those edges (since that is how we find all the index cells that need
  // to be updated).  The edges themselves are restored automatically when
  // UpdateEdges returns from its recursive call, but the InteriorTracker
  // state needs to be restored explicitly.
  //
  // Here we first update the InteriorTracker state for removed edges to
  // correspond to the exit vertex of this cell, and then save the
  // InteriorTracker state.  This state will be restored by UpdateEdges when
  // it is finished processing the contents of this cell.
  if (tracker->is_active() && !edges->empty() &&
      is_shape_being_removed((*edges)[0]->face_edge->shape_id)) {
    // We probably need to update the InteriorTracker.  ("Probably" because
    // it's possible that all shapes being removed do not have interiors.)
    if (!tracker->at_cellid(pcell.id())) {
      tracker->MoveTo(pcell.GetEntryVertex());
    }
    tracker->DrawTo(pcell.GetExitVertex());
    tracker->set_next_cellid(pcell.id().next());
    for (int e = 0; e < edges->size(); ++e) {
      FaceEdge const* face_edge = (*edges)[e]->face_edge;
      if (!is_shape_being_removed(face_edge->shape_id)) {
        break;  // All shapes being removed come first.
      }
      if (face_edge->has_interior) {
        tracker->TestEdge(face_edge->shape_id, face_edge->va, face_edge->vb);
      }
    }
  }
  // Save the state of the edges being removed, so that it can be restored
  // when we are finished processing this cell and its children.  We don't
  // need to save the state of the edges being added because they aren't being
  // removed from "edges" and will therefore be updated normally as we visit
  // this cell and its children.
  tracker->SaveAndClearStateBefore(pending_additions_begin_);

  // Create a FaceEdge for each edge in this cell that isn't being removed.
  vector<FaceEdge>* face_edges = alloc->mutable_face_edges();
  face_edges->clear();
  bool tracker_moved = false;
  S2ShapeIndexCell const* cell = iter.cell();
  for (int s = 0; s < cell->num_shapes(); ++s) {
    S2ClippedShape const& clipped = cell->clipped(s);
    int shape_id = clipped.shape_id();
    S2Shape const* shape = shapes_[shape_id];
    if (shape == nullptr) continue;  // This shape is being removed.
    int num_clipped = clipped.num_edges();

    // If this shape has an interior, start tracking whether we are inside the
    // shape.  UpdateEdges() wants to know whether the entry vertex of this
    // cell is inside the shape, but we only know whether the center of the
    // cell is inside the shape, so we need to test all the edges against the
    // line segment from the cell center to the entry vertex.
    FaceEdge edge;
    edge.shape_id = shape->id();
    edge.has_interior = shape->has_interior();
    if (edge.has_interior) {
      tracker->AddShape(shape_id, clipped.contains_center());
      // There might not be any edges in this entire cell (i.e., it might be
      // in the interior of all shapes), so we delay updating the tracker
      // until we see the first edge.
      if (!tracker_moved && num_clipped > 0) {
        tracker->MoveTo(pcell.GetCenter());
        tracker->DrawTo(pcell.GetEntryVertex());
        tracker->set_next_cellid(pcell.id());
        tracker_moved = true;
      }
    }
    for (int i = 0; i < num_clipped; ++i) {
      int e = clipped.edge(i);
      edge.edge_id = e;
      shape->GetEdge(e, &edge.va, &edge.vb);
      edge.max_level = GetEdgeMaxLevel(*edge.va, *edge.vb);
      if (edge.has_interior) tracker->TestEdge(shape_id, edge.va, edge.vb);
      if (!S2EdgeUtil::ClipToPaddedFace(*edge.va, *edge.vb, pcell.id().face(),
                                        kCellPadding, &edge.a, &edge.b)) {
        LOG(DFATAL) << "Invariant failure in S2ShapeIndex";
      }
      face_edges->push_back(edge);
    }
  }
  // Now create a ClippedEdge for each FaceEdge, and put them in "new_edges".
  vector<ClippedEdge const*> new_edges;
  for (int i = 0; i < face_edges->size(); ++i) {
    FaceEdge const& face_edge = (*face_edges)[i];
    ClippedEdge* clipped = alloc->NewClippedEdge();
    clipped->face_edge = &face_edge;
    clipped->bound = S2EdgeUtil::GetClippedEdgeBound(face_edge.a, face_edge.b,
                                                     pcell.bound());
    new_edges.push_back(clipped);
  }
  // Discard any edges from "edges" that are being removed, and append the
  // remainder to "new_edges".  (This keeps the edges sorted by shape id.)
  for (int i = 0; i < edges->size(); ++i) {
    ClippedEdge const* clipped = (*edges)[i];
    if (!is_shape_being_removed(clipped->face_edge->shape_id)) {
      new_edges.insert(new_edges.end(), edges->begin() + i, edges->end());
      break;
    }
  }
  // Update the edge list and delete this cell from the index.
  edges->swap(new_edges);
  cell_map_.erase(pcell.id());
  delete cell;
}

// Attempt to build an index cell containing the given edges, and return true
// if successful.  (Otherwise the edges should be subdivided further.)
bool S2ShapeIndex::MakeIndexCell(S2PaddedCell const& pcell,
                                 vector<ClippedEdge const*> const& edges,
                                 InteriorTracker* tracker) {
  if (edges.empty() && tracker->shape_ids().empty()) {
    // No index cell is needed.  (In most cases this situation is detected
    // before we get to this point, but this can happen when all shapes in a
    // cell are removed.)
    return true;
  }

  // Count the number of edges that have not reached their maximum level yet.
  // Return false if there are too many such edges.
  int count = 0;
  for (int e = 0; e < edges.size(); ++e) {
    count += (pcell.level() < edges[e]->face_edge->max_level);
    if (count > options_.max_edges_per_cell())
      return false;
  }

  // Possible optimization: Continue subdividing as long as exactly one child
  // of "pcell" intersects the given edges.  This can be done by finding the
  // bounding box of all the edges and calling ShrinkToFit():
  //
  // S2CellId cellid = pcell.ShrinkToFit(GetRectBound(edges));
  //
  // Currently this is not beneficial; it slows down construction by 4-25%
  // (mainly computing the union of the bounding rectangles) and also slows
  // down queries (since more recursive clipping is required to get down to
  // the level of a spatial index cell).  But it may be worth trying again
  // once "contains_center" is computed and all algorithms are modified to
  // take advantage of it.

  // We update the InteriorTracker as follows.  For every S2Cell in the index
  // we construct two edges: one edge from entry vertex of the cell to its
  // center, and one from the cell center to its exit vertex.  Here "entry"
  // and "exit" refer the S2CellId ordering, i.e. the order in which points
  // are encountered along the S2 space-filling curve.  The exit vertex then
  // becomes the entry vertex for the next cell in the index, unless there are
  // one or more empty intervening cells, in which case the InteriorTracker
  // state is unchanged because the intervening cells have no edges.

  // Shift the InteriorTracker focus point to the center of the current cell.
  if (tracker->is_active() && !edges.empty()) {
    if (!tracker->at_cellid(pcell.id())) {
      tracker->MoveTo(pcell.GetEntryVertex());
    }
    tracker->DrawTo(pcell.GetCenter());
    TestAllEdges(edges, tracker);
  }
  // Allocate and fill a new index cell.  To get the total number of shapes we
  // need to merge the shapes associated with the intersecting edges together
  // with the shapes that happen to contain the cell center.
  ShapeIdSet const& cshape_ids = tracker->shape_ids();
  int num_shapes = CountShapes(edges, cshape_ids);
  S2ShapeIndexCell* cell = new S2ShapeIndexCell;
  S2ClippedShape* base = cell->add_shapes(num_shapes);

  // To fill the index cell we merge the two sources of shapes: "edge shapes"
  // (those that have at least one edge that intersects this cell), and
  // "containing shapes" (those that contain the cell center).  We keep track
  // of the index of the next intersecting edge and the next containing shape
  // as we go along.  Both sets of shape ids are already sorted.
  int enext = 0;
  ShapeIdSet::const_iterator cnext = cshape_ids.begin();
  for (int i = 0; i < num_shapes; ++i) {
    S2ClippedShape* clipped = base + i;
    int eshape_id = num_shape_ids(), cshape_id = eshape_id;  // Sentinels
    if (enext != edges.size()) {
      eshape_id = edges[enext]->face_edge->shape_id;
    }
    if (cnext != cshape_ids.end()) {
      cshape_id = *cnext;
    }
    int ebegin = enext;
    if (cshape_id < eshape_id) {
      // The entire cell is in the shape interior.
      clipped->Init(cshape_id, 0);
      clipped->set_contains_center(true);
      ++cnext;
    } else {
      // Count the number of edges for this shape and allocate space for them.
      while (enext < edges.size() &&
             edges[enext]->face_edge->shape_id == eshape_id) {
        ++enext;
      }
      clipped->Init(eshape_id, enext - ebegin);
      for (int e = ebegin; e < enext; ++e) {
        clipped->set_edge(e - ebegin, edges[e]->face_edge->edge_id);
      }
      if (cshape_id == eshape_id) {
        clipped->set_contains_center(true);
        ++cnext;
      }
    }
  }
  // UpdateEdges() visits cells in increasing order of S2CellId, so during
  // initial construction of the index all insertions happen at the end.  It
  // is much faster to give an insertion hint in this case.  Otherwise the
  // hint doesn't do much harm.  With more effort we could provide a hint even
  // during incremental updates, but this is probably not worth the effort.
  cell_map_.insert(cell_map_.end(), std::make_pair(pcell.id(), cell));

  // Shift the InteriorTracker focus point to the exit vertex of this cell.
  if (tracker->is_active() && !edges.empty()) {
    tracker->DrawTo(pcell.GetExitVertex());
    TestAllEdges(edges, tracker);
    tracker->set_next_cellid(pcell.id().next());
  }
  return true;
}

// Call tracker->TestEdge() on all edges from shapes that have interiors.
/* static */
void S2ShapeIndex::TestAllEdges(vector<ClippedEdge const*> const& edges,
                                InteriorTracker* tracker) {
  for (int e = 0; e < edges.size(); ++e) {
    FaceEdge const* face_edge = edges[e]->face_edge;
    if (face_edge->has_interior) {
      tracker->TestEdge(face_edge->shape_id, face_edge->va, face_edge->vb);
    }
  }
}

// Return the number of distinct shapes that are either associated with the
// given edges, or that are currently stored in the InteriorTracker.
/* static */
int S2ShapeIndex::CountShapes(vector<ClippedEdge const*> const& edges,
                              ShapeIdSet const& cshape_ids) {
  int count = 0;
  int last_shape_id = -1;
  ShapeIdSet::const_iterator cnext = cshape_ids.begin();  // Next shape
  for (int e = 0; e < edges.size(); ++e) {
    if (edges[e]->face_edge->shape_id != last_shape_id) {
      ++count;
      last_shape_id = edges[e]->face_edge->shape_id;
      // Skip over any containing shapes up to and including this one,
      // updating "count" appropriately.
      for (; cnext != cshape_ids.end(); ++cnext) {
        if (*cnext > last_shape_id) break;
        if (*cnext < last_shape_id) ++count;
      }
    }
  }
  // Count any remaining containing shapes.
  count += (cshape_ids.end() - cnext);
  return count;
}

int S2ShapeIndex::GetNumEdges() const {
  // There is no need to apply updates before counting edges, since shapes are
  // added or removed from the "shapes_" vector immediately.
  int num_edges = 0;
  for (int id = 0; id < shapes_.size(); ++id) {
    if (shapes_[id] == nullptr) continue;
    num_edges += shapes_[id]->num_edges();
  }
  return num_edges;
}

size_t S2ShapeIndex::BytesUsed() const {
  MaybeApplyUpdates();
  size_t size = sizeof(*this);
  size += shapes_.capacity() * sizeof(shapes_[0]);
  // cell_map_ itself is already included in sizeof(*this).
  size += cell_map_.bytes_used() - sizeof(cell_map_);
  size += cell_map_.size() * sizeof(S2ShapeIndexCell);
  Iterator it;
  for (it.InitStale(*this); !it.Done(); it.Next()) {
    S2ShapeIndexCell const* cell = it.cell();
    size += cell->shapes_.capacity() * sizeof(S2ClippedShape);
    for (int s = 0; s < cell->num_shapes(); ++s) {
      S2ClippedShape const& clipped = cell->clipped(s);
      if (!clipped.is_inline()) {
        size += clipped.num_edges() * sizeof(int32);
      }
    }
  }
  // There are no pending removals because all updates have been applied.
  return size;
}
