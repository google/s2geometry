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

#include <math.h>
#include <algorithm>

#include "base/atomicops.h"
#include <gflags/gflags.h>
#include "base/spinlock.h"
#include "r1interval.h"
#include "r2.h"
#include "r2rect.h"
#include "s2cellid.h"
#include "s2edgeutil.h"
#include "s2paddedcell.h"

using std::max;
using std::vector;

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

// FLAGS_s2shapeindex_default_max_edges_per_cell:
//
// The default maximum number of edges per cell (not counting "long" edges).
// If a cell has more than this many edges, and it is not a leaf cell, then it
// is subdivided.  This flag can be overridden via S2ShapeIndexOptions.
DEFINE_int32(
    s2shapeindex_default_max_edges_per_cell, 10,  // TODO(ericv): Adjust
    "Default maximum number of edges (not counting 'long' edges) per cell");

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
  return NULL;
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
      index_status_(FRESH),
      update_state_(NULL) {
}

S2ShapeIndex::S2ShapeIndex(S2ShapeIndexOptions const& options)
    : options_(options),
      pending_additions_begin_(0),
      index_status_(FRESH),
      update_state_(NULL) {
}

void S2ShapeIndex::Init(S2ShapeIndexOptions const& options) {
  DCHECK(shapes_.empty());
  options_ = options;
}

S2ShapeIndex::~S2ShapeIndex() {
  Reset();
}

void S2ShapeIndex::Add(S2Shape const* shape) {
  // Additions are processed lazily by ApplyUpdates().
  shape->id_ = shapes_.size();
  shapes_.push_back(shape);
  base::subtle::NoBarrier_Store(&index_status_, STALE);
}

void S2ShapeIndex::Remove(S2Shape const* shape) {
  LOG(FATAL) << "Not implemented yet";
  base::subtle::NoBarrier_Store(&index_status_, STALE);
}

void S2ShapeIndex::Reset() {
  Iterator it;
  for (it.InitStale(*this); !it.Done(); it.Next()) {
    delete it.cell();
  }
  for (int id = 0; id < shapes_.size(); ++id) {
    S2Shape const* shape = shapes_[id];
    if (shape) shape->Release();
  }
  cell_map_.clear();
  // Note that vector::clear() does not actually free storage.
  vector<S2Shape const*> empty_shapes;
  shapes_.swap(empty_shapes);
  pending_additions_begin_ = 0;
  pending_removals_.reset();
  DCHECK(update_state_ == NULL);
  base::subtle::NoBarrier_Store(&index_status_, FRESH);
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
  int shape_id;       // The shape that this edge belongs to
  int edge_id;        // Edge id within that shape
  int max_level;      // Not desirable to subdivide this edge beyond this level
  R2Point a, b;       // The edge endpoints, clipped to a given face
  S2Point const* va;  // The original S2Loop vertices
  S2Point const* vb;
};

struct S2ShapeIndex::ClippedEdge {
  FaceEdge const* orig;   // The original unclipped edge
  R2Rect bound;           // Bounding box for the clipped portion
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

  // Return true if AddShape() has been called at least once.
  bool is_active() const { return is_active_; }

  // Add a shape whose interior should be tracked.  This should be followed by
  // calling TestEdge() with every edge of the given shape.
  void AddShape(S2Shape const* shape);

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
  inline void TestEdge(S2Shape const* shape,
                       S2Point const* c, S2Point const* d);

  // The set of shape ids that contain the current focus.
  ShapeIdSet const& shape_ids() const { return shape_ids_; }

  // Indicate that the caller has finished processing the given S2CellId.  By
  // using this method together with AtCellId(), the caller can avoid calling
  // MoveTo() in cases where the exit vertex of the previous cell is the same
  // as the entry vertex of the current cell.
  void DoneCellId(S2CellId cellid) {
    next_cellid_ = cellid.range_max().next();
  }

  // Return true if the focus is already at the entry vertex of the given
  // S2CellId (provided that the caller calls DoneCellId() as each cell is
  // processed).
  bool AtCellId(S2CellId cellid) const {
    return cellid.range_min() == next_cellid_;
  }

 private:
  // Remove "shape_id" from shape_ids_ if it exists, otherwise insert it.
  void ToggleShape(int shape_id);

  bool is_active_;
  S2Point a_, b_;
  S2Point const* c_;
  S2CellId next_cellid_;
  S2EdgeUtil::EdgeCrosser crosser_;
  ShapeIdSet shape_ids_;
};

// As shapes are added, we compute which ones contain the start of the
// S2CellId space-filling curve by drawing an edge from S2::Origin() to this
// point and counting how many shape edges cross this edge.
S2ShapeIndex::InteriorTracker::InteriorTracker()
    : is_active_(false), b_(S2::Origin()),
      next_cellid_(S2CellId::Begin(S2CellId::kMaxLevel)) {
  DrawTo(S2::FaceUVtoXYZ(0, -1, -1).Normalize());  // S2CellId curve start
}

void S2ShapeIndex::InteriorTracker::AddShape(S2Shape const* shape) {
  is_active_ = true;
  if (shape->contains_origin()) {
    ToggleShape(shape->id());
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
    vector<int>::iterator pos = shape_ids_.begin();
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
  c_ = NULL;
}

inline void S2ShapeIndex::InteriorTracker::TestEdge(S2Shape const* shape,
                                                    S2Point const* c,
                                                    S2Point const* d) {
  if (c != c_) { crosser_.RestartAt(c); }
  if (crosser_.EdgeOrVertexCrossing(d)) {
    ToggleShape(shape->id());
  }
  c_ = d;
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
    update_state_ = new UpdateState;
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
    delete update_state_;
    update_state_ = NULL;
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

// This method updates the index by applying all pending additions and
// removals.  It does *not* update index_status_ (see ApplyUpdatesThreadSafe).
void S2ShapeIndex::ApplyUpdatesInternal() {
  CHECK(cell_map_.empty()) << "Incremental updates not supported yet";
  vector<FaceEdge> all_edges[6];
  ReserveSpace(all_edges);

  InteriorTracker tracker;
  for (int id = pending_additions_begin_; id < shapes_.size(); ++id) {
    AddShapeEdges(id, all_edges, &tracker);
  }
  // Whether to insert or remove a given shape is represented implicitly:
  //  - if shape(id) is not NULL, its edges are inserted.
  //  - if shape(id) is NULL, its edges are removed.
  for (int face = 0; face < 6; ++face) {
    UpdateFaceEdges(face, all_edges[face], &tracker);
    // Save memory by clearing vectors after we are done with them.
    vector<FaceEdge> empty;
    all_edges[face].swap(empty);
  }
  pending_additions_begin_ = shapes_.size();
  // It is the caller's responsibility to update index_status_.
}

// Reserve an appropriate amount of space for the top-level face edges.  These
// vectors are responsible for most of the temporary memory usage during index
// construction.  Furthermore, if the arrays are grown via push_back() then up
// to 10% of the total run time consists of copying data as these arrays grow.
// So it is worthwhile to preallocate space via reserve().
void S2ShapeIndex::ReserveSpace(vector<FaceEdge> all_edges[6]) const {
  // If the number of edges is relatively small, then the fastest approach is
  // to simply reserve space on every face for the maximum possible number of
  // edges.
  int num_edges = 0;
  for (int id = pending_additions_begin_; id < shapes_.size(); ++id) {
    num_edges += shapes_[id]->num_edges();
  }
  int const kMaxCheapMemoryBytes = 10 << 20;  // 10MB
  int const kMaxCheapNumEdges = kMaxCheapMemoryBytes / (6 * sizeof(FaceEdge));
  if (num_edges <= kMaxCheapNumEdges) {
    for (int face = 0; face < 6; ++face) {
      all_edges[face].reserve(num_edges);
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
  int const sample_interval = max(1, num_edges / kDesiredSampleSize);

  // Initialize "edge_id" to be midway through the first sample interval.
  // Because samples are equally spaced the actual sample size may differ
  // slightly from the desired sample size.
  int edge_id = sample_interval / 2;
  int const actual_sample_size = (num_edges + edge_id) / sample_interval;
  int face_count[6] = { 0, 0, 0, 0, 0, 0 };
  for (int id = pending_additions_begin_; id < shapes_.size(); ++id) {
    S2Shape const* shape = shapes_[id];
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
    all_edges[face].reserve(fraction * num_edges);
  }
}

// Clip all edges of the given shape to the six cube faces, and add the
// clipped edges to "all_edges".
void S2ShapeIndex::AddShapeEdges(int id, vector<FaceEdge> all_edges[6],
                                 InteriorTracker* tracker) const {
  FaceEdge edge;
  edge.shape_id = id;
  S2Shape const* shape = shapes_[id];
  bool has_interior = shape->has_interior();
  if (has_interior) tracker->AddShape(shape);
  int num_edges = shape->num_edges();
  for (int e = 0; e < num_edges; ++e) {
    edge.edge_id = e;
    shape->GetEdge(e, &edge.va, &edge.vb);
    edge.max_level = GetEdgeMaxLevel(*edge.va, *edge.vb);
    if (has_interior) tracker->TestEdge(shape, edge.va, edge.vb);
    // Fast path: both endpoints are on the same face, and are far enough from
    // the edge of the face that don't intersect any (padded) adjacent face.
    int a_face = S2::GetFace(*edge.va);
    if (a_face == S2::GetFace(*edge.vb)) {
      S2::ValidFaceXYZtoUV(a_face, *edge.va, &edge.a);
      S2::ValidFaceXYZtoUV(a_face, *edge.vb, &edge.b);
      double const kMaxUV = 1 - kCellPadding;
      if (fabs(edge.a[0]) <= kMaxUV && fabs(edge.a[1]) <= kMaxUV &&
          fabs(edge.b[0]) <= kMaxUV && fabs(edge.b[1]) <= kMaxUV) {
        all_edges[a_face].push_back(edge);
        continue;
      }
    }
    // Otherwise we simply clip the edge to all six faces.
    for (int face = 0; face < 6; ++face) {
      if (S2EdgeUtil::ClipToPaddedFace(*edge.va, *edge.vb, face, kCellPadding,
                                       &edge.a, &edge.b)) {
        all_edges[face].push_back(edge);
      }
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
class S2ShapeIndex::EdgeAllocator {
 public:
  EdgeAllocator() : size_(0) {}

  // Return a pointer to a newly allocated edge.
  ClippedEdge* New() {
    if (size_ == edges_.size()) edges_.push_back(new ClippedEdge);
    return edges_[size_++];
  }
  // Return the number of allocated edges.
  size_t size() const { return size_; }

  // Reset the allocator to only contain the first "size" allocated edges.
  void Reset(size_t size) { size_ = size; }

  ~EdgeAllocator() {
    for (int i = 0; i < edges_.size(); ++i) {
      delete edges_[i];
    }
  }

 private:
  // We can't use vector<ClippedEdge> because edges are not allowed to move
  // once they have been allocated.  Instead we keep a pool of allocated edges
  // that are all deleted together at the end.
  size_t size_;
  vector<ClippedEdge*> edges_;

  DISALLOW_COPY_AND_ASSIGN(EdgeAllocator);
};

// Given a face and a vector of edges that intersect that face, add or
// remove all the edges from the index.  (An edge is added if shape(id) is
// not NULL, and removed otherwise.)
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
    clipped.orig = &face_edges[e];
    clipped.bound = R2Rect::FromPointPair(face_edges[e].a, face_edges[e].b);
    clipped_edge_storage.push_back(clipped);
    clipped_edges.push_back(&clipped_edge_storage.back());
    bound.AddRect(clipped.bound);
  }
  // Construct the initial face cell containing all the edges, and then update
  // all the edges in the index recursively.
  S2PaddedCell pcell(S2CellId::FromFace(face), kCellPadding);
  if (tracker->shape_ids().empty()) {
    S2CellId cellid = pcell.ShrinkToFit(bound);
    if (cellid != pcell.id()) pcell = S2PaddedCell(cellid, kCellPadding);
  }
  EdgeAllocator alloc;
  UpdateEdges(pcell, clipped_edges, tracker, &alloc);
}

// Given a cell and a set of ClippedEdges whose bounding boxes intersect that
// cell, add or remove all the edges from the index.  Temporary space for
// edges that need to be subdivided is allocated from the given EdgeAllocator.
void S2ShapeIndex::UpdateEdges(S2PaddedCell const& pcell,
                               vector<ClippedEdge const*> const& edges,
                               InteriorTracker* tracker,
                               EdgeAllocator* alloc) {
  // This function is recursive with a maximum recursion depth of 30
  // (S2CellId::kMaxLevel).  Note that using an explicit stack does not seem
  // to be any faster based on profiling.

  // Check whether we have few enough edges to make a leaf cell.
  if (MakeLeafCell(pcell, edges, tracker)) return;

  // Reserve space for the edges that will be passed to each child.  This is
  // important since otherwise the running time is dominated by the time
  // required to grow the vectors.  The amount of memory involved is
  // relatively small, so we simply reserve the maximum space for every child.
  vector<ClippedEdge const*> child_edges[2][2];  // [i][j]
  int num_edges = edges.size();
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
    ClippedEdge const* edge = edges[e];
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
  // Now recursively update the edges in each child.  We call the children in
  // increasing order of S2CellId so that when the index is first constructed,
  // all insertions into cell_map_ are at the end (which is much faster).
  for (int pos = 0; pos < 4; ++pos) {
    int i, j;
    pcell.GetChildIJ(pos, &i, &j);
    if (!child_edges[i][j].empty() || !tracker->shape_ids().empty()) {
      UpdateEdges(S2PaddedCell(pcell, i, j), child_edges[i][j], tracker, alloc);
    }
  }
  // Free any temporary edges that were allocated during clipping.
  alloc->Reset(alloc_size);
}

// Given an edge and an interval "middle" along the v-axis, clip the edge
// against the boundaries of "middle" and add the edge to the corresponding
// children.
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
  FaceEdge const& e = *edge->orig;
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
S2ShapeIndex::ClippedEdge const*
S2ShapeIndex::ClipVBound(ClippedEdge const* edge, int v_end, double v,
                         EdgeAllocator* alloc) {
  // See comments in ClipUBound.
  if (v_end == 0) {
    if (edge->bound[1].lo() >= v) return edge;
  } else {
    if (edge->bound[1].hi() <= v) return edge;
  }
  FaceEdge const& e = *edge->orig;
  double u = edge->bound[0].ClampPoint(
      S2EdgeUtil::InterpolateDouble(v, e.a[1], e.b[1], e.a[0], e.b[0]));
  int u_end = v_end ^ ((e.a[0] > e.b[0]) != (e.a[1] > e.b[1]));
  return UpdateBound(edge, u_end, u, v_end, v, alloc);
}

// Given an edge and two bound endpoints that need to be updated, allocate and
// return a new edge with the updated bound.
inline S2ShapeIndex::ClippedEdge const*
S2ShapeIndex::UpdateBound(ClippedEdge const* edge, int u_end, double u,
                          int v_end, double v, EdgeAllocator* alloc) {
  ClippedEdge* clipped = alloc->New();
  clipped->orig = edge->orig;
  clipped->bound[0][u_end] = u;
  clipped->bound[1][v_end] = v;
  clipped->bound[0][1-u_end] = edge->bound[0][1-u_end];
  clipped->bound[1][1-v_end] = edge->bound[1][1-v_end];
  DCHECK(!clipped->bound.is_empty());
  DCHECK(edge->bound.Contains(clipped->bound));
  return clipped;
}

// Attempt to build a leaf cell containing the given edges, and return true if
// successful.  (Otherwise the edges should be subdivided further.)
bool S2ShapeIndex::MakeLeafCell(S2PaddedCell const& pcell,
                                vector<ClippedEdge const*> const& edges,
                                InteriorTracker* tracker) {
  // Count the number of edges that have not reached their maximum level yet.
  // Return false if there are too many such edges.
  int count = 0;
  for (int e = 0; e < edges.size(); ++e) {
    count += (pcell.level() < edges[e]->orig->max_level);  // branchless
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
    if (!tracker->AtCellId(pcell.id())) {
      tracker->MoveTo(pcell.GetEntryVertex());
    }
    tracker->DrawTo(pcell.GetCenter());
    for (int e = 0; e < edges.size(); ++e) {
      FaceEdge const* orig = edges[e]->orig;
      tracker->TestEdge(shape(orig->shape_id), orig->va, orig->vb);
    }
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
  // as we go along.
  int enext = 0;
  ShapeIdSet::const_iterator cnext = cshape_ids.begin();
  for (int i = 0; i < num_shapes; ++i) {
    S2ClippedShape* clipped = base + i;
    int eshape_id = num_shape_ids(), cshape_id = eshape_id;  // Sentinels
    if (enext != edges.size()) eshape_id = edges[enext]->orig->shape_id;
    if (cnext != cshape_ids.end()) cshape_id = *cnext;
    int ebegin = enext;
    if (cshape_id < eshape_id) {
      // The entire cell is in the shape interior.
      clipped->Init(cshape_id, 0);
      clipped->set_contains_center(true);
      ++cnext;
    } else {
      // Count the number of edges for this shape and allocate space for them.
      while (enext < edges.size() && edges[enext]->orig->shape_id == eshape_id)
        ++enext;
      clipped->Init(eshape_id, enext - ebegin);
      for (int e = ebegin; e < enext; ++e) {
        clipped->set_edge(e - ebegin, edges[e]->orig->edge_id);
      }
      if (cshape_id == eshape_id) {
        clipped->set_contains_center(true);
        ++cnext;
      }
    }
  }
  // UpdateEdges() visits cells in increasing order of S2CellId, so during
  // initial construction of the index all insertions happen at the end.  It
  // is much faster to give an insertion hint in this case.
  cell_map_.insert(cell_map_.end(), std::make_pair(pcell.id(), cell));

  // Shift the InteriorTracker focus point to the exit vertex of this cell.
  if (tracker->is_active() && !edges.empty()) {
    tracker->DrawTo(pcell.GetExitVertex());
    for (int e = 0; e < edges.size(); ++e) {
      FaceEdge const* orig = edges[e]->orig;
      tracker->TestEdge(shape(orig->shape_id), orig->va, orig->vb);
    }
    tracker->DoneCellId(pcell.id());
  }
  return true;
}

// Return the number of distinct shapes that are either associated with the
// given edges, or that are currently stored in the InteriorTracker.
int S2ShapeIndex::CountShapes(vector<ClippedEdge const*> const& edges,
                              ShapeIdSet const& cshape_ids) {
  int count = 0;
  int last_shape_id = -1;
  ShapeIdSet::const_iterator cnext = cshape_ids.begin();  // Next shape
  for (int e = 0; e < edges.size(); ++e) {
    if (edges[e]->orig->shape_id != last_shape_id) {
      ++count;
      last_shape_id = edges[e]->orig->shape_id;
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

