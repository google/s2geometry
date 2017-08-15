// Copyright 2017 Google Inc. All Rights Reserved.
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

#ifndef S2_S2CLOSESTEDGEQUERY_BASE_H_
#define S2_S2CLOSESTEDGEQUERY_BASE_H_

#include <memory>
#include <vector>

#include <glog/logging.h>
#include "s2/third_party/absl/container/inlined_vector.h"
#include "s2/util/btree/btree_set.h"
#include "s2/_fpcontractoff.h"
#include "s2/s1angle.h"
#include "s2/s1chordangle.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cellid.h"
#include "s2/s2cellunion.h"
#include "s2/s2regioncoverer.h"
#include "s2/s2shapeindex.h"
#include "s2/s2shapeutil.h"
#include "s2/util/gtl/dense_hash_set.h"

// S2ClosestEdgeQueryBase is a templatized class for finding the closest
// edge(s) between two geometries.  It is not intended to be used directly,
// but rather to serve as the implementation of various specialized classes
// with more convenient APIs (such as S2ClosestEdgeQuery).  It is flexible
// enough so that it can be adapted to compute maximum distances and even
// potentially Hausdorff distances.
//
// By using the appropriate options, this class can answer questions such as:
//
//  - Find the minimum distance between two geometries A and B.
//  - Find all edges of geometry A that are within a distance D of geometry B.
//  - Find the k edges of geometry A that are closest to a given point P.
//
// Results can be filtered so that at most one edge per distinct S2Shape is
// returned, e.g. in order to find the closest 3 polygons to a given point
// rather than the closest 3 edges.  You can also specify whether polygons
// should include their interiors (i.e., if a point is contained by a polygon,
// should the distance be zero or should it be measured to the polygon
// boundary?)
//
// The input geometries may consist of any number of points, polylines, and
// polygons (collectively referred to as "shapes").  Shapes do not need to be
// disjoint; they may overlap or intersect arbitrarily.  The implementation is
// designed to be fast for both simple and complex geometries.
//
// The Distance template argument is used to represent distances.  Usually
// this type is a thin wrapper around S1ChordAngle, but another distance type
// may be substituted as long as it implements the API below.  This can be
// used to change the comparison function (e.g., to find the furthest edges
// from the target), or to get more accuracy if desired.
//
// The Distance concept is as follows:
//
// class Distance {
//  public:
//   // Default and copy constructors, assignment operator:
//   Distance();
//   Distance(Distance const&);
//   Distance& operator=(Distance const&);
//
//   // Factory methods:
//   static Distance Zero();      // Returns a zero distance.
//   static Distance Infinity();  // Larger than any valid distance.
//   static Distance Negative();  // Smaller than any valid distance.
//
//   // Comparison operators:
//   friend bool operator==(Distance x, Distance y);
//   friend bool operator<(Distance x, Distance y);
//
//   // Subtraction operator (needed to implement Options::max_error):
//   friend Distance operator-(Distance x, Distance y);
//
//   // Method that returns an upper bound on the S1Angle corresponding to the
//   // given Distance (needed to efficiently implement Options::max_distance).
//   // For example, if Distance measures WGS84 ellipsoid distance then the
//   // corresponding S1Angle would need to be 0.56% larger.
//   static S1Angle GetAngleBound(Distance x);
// };
template <class Distance>
class S2ClosestEdgeQueryBase {
 public:
  // Options that control the set of edges returned.  Note that by default
  // *all* edges are returned, so you will always want to set either the
  // max_edges() option or the max_distance() option (or both).
  class Options {
   public:
    Options();

    // Specifies that at most "max_edges" edges should be returned.
    //
    // REQUIRES: max_edges >= 1
    // DEFAULT: numeric_limits<int>::max()
    int max_edges() const;
    void set_max_edges(int max_edges);
    static constexpr int kMaxMaxEdges = std::numeric_limits<int>::max();

    // Specifies that only edges whose distance to the target is less than
    // "max_distance" should be returned.
    //
    // DEFAULT: Distance::Infinity()
    Distance max_distance() const;
    void set_max_distance(Distance max_distance);

    // Specifies that edges up to max_error() further away than the true
    // closest edges may be substituted in the result set, as long as such
    // edges satisfy all the remaining search criteria (such as max_distance).
    // This option only has an effect if max_edges() is also specified;
    // otherwise all edges closer than max_distance() will always be returned.
    //
    // Note that this does not affect how the distance between edges is
    // computed; it simply gives the algorithm permission to stop the search
    // early as soon as the best possible improvement drops below max_error().
    //
    // This can be used to implement distance predicates efficiently.  For
    // example, to determine whether the minimum distance is less than D, set
    // max_edges() == 1 and max_distance() == max_error() == D.  This causes
    // the algorithm to terminate as soon as it finds any edge whose distance
    // is less than D, rather than continuing to search for an edge that is
    // even closer.
    //
    // DEFAULT: Distance::Zero()
    Distance max_error() const;
    void set_max_error(Distance max_error);

    // Specifies that if the query target intersects the interior of a
    // polygonal shape, then the distance to that shape is considered to be
    // zero.  Such shapes are returned as a (shape_id, edge_id) pair with
    // (edge_id == -1).
    //
    // Polygon interiors are assumed to be semi-open sets.  (Note that since
    // all boundary edges are considered anyway, it doesn't really matter
    // whether they are considered to belong to the interior or not.)
    //
    // DEFAULT: false
    // TODO(ericv): Implement this, and possibly make true by default.
    bool include_interiors() const;
    void set_include_interiors(bool include_interiors);

    // Specifies that distances should be computed by examining every edge
    // rather than using the S2ShapeIndex.  This is useful for testing,
    // benchmarking, and debugging.
    //
    // DEFAULT: false
    bool use_brute_force() const;
    void set_use_brute_force(bool use_brute_force);

   private:
    int max_edges_ = kMaxMaxEdges;
    Distance max_distance_ = Distance::Infinity();
    Distance max_error_ = Distance::Zero();
    bool include_interiors_ = false;
    bool use_brute_force_ = false;
  };

  // The Target class represents the geometry to which the distance is
  // measured.  For example, there can be subtypes for measuring the distance
  // to a point, an edge, or to an S2ShapeIndex (an arbitrary collection of
  // geometry).
  //
  // Implementations do *not* need to be thread-safe.  They may cache data or
  // allocate temporary data structures in order to improve performance.
  class Target {
   public:
    virtual ~Target() {}

    // Specifies the maximum error allowed when computing distances in the
    // UpdateMinDistance() methods.  This method must return "true" if the
    // Target subtype takes advantage of this parameter.  (Most Target types
    // can use the default implementation which simply returns false.)
    virtual bool set_max_error(Distance const& max_error) { return false; }

    // Specifies the maximum number of edges for which distances should be
    // computed by examining every edge rather than using the S2ShapeIndex.
    // This can be estimated for a particular Target type using benchmarks.
    virtual int max_brute_force_edges() const = 0;

    // Returns an S2Cap that bounds the set of points whose distance to the
    // target is Distance::Zero().
    virtual S2Cap GetCapBound() const = 0;

    // If the distance to the edge (v0, v1) is less than "min_dist", then
    // updates "min_dist" and returns true.  Otherwise returns false.
    virtual bool UpdateMinDistance(S2Point const& v0, S2Point const& v1,
                                   Distance* min_dist) const = 0;

    // If the distance to the given S2Cell (including its interior) is less
    // than "min_dist", then updates "min_dist" and returns true.  Otherwise
    // returns false.
    virtual bool UpdateMinDistance(S2Cell const& cell,
                                   Distance* min_dist) const = 0;

    // TODO(ericv): Needed for Options::include_interiors().
    // Returns the set of polygonal shapes (but not more than "max_shapes")
    // whose interior intersects the target.  (Generally polygon interiors are
    // modeled as semi-open sets, but this is not a requirement.)
    // virtual std::vector<int> GetIntersectingShapes(S2ShapeIndex const& index,
    //                                                int max_shapes) const = 0;
  };

  struct Result {
    Distance distance;  // The distance from the target to this edge.
    int32 shape_id;     // Identifies a shape.
    int32 edge_id;      // Identifies an edge within the shape.

    // The default constructor yields an invalid result.
    Result() : distance(Distance::Infinity()), shape_id(-1), edge_id(-1) {}

    Result(Distance _distance, int32 _shape_id, int32 _edge_id)
        : distance(_distance), shape_id(_shape_id), edge_id(_edge_id) {}

    friend bool operator==(Result const& x, Result const& y) {
      return (x.distance == y.distance &&
              x.shape_id == y.shape_id &&
              x.edge_id == y.edge_id);
    }

    // Compares edges first by distance, then by (shape_id, edge_id).
    friend bool operator<(Result const& x, Result const& y) {
      if (x.distance < y.distance) return true;
      if (y.distance < x.distance) return false;
      if (x.shape_id < y.shape_id) return true;
      if (y.shape_id < x.shape_id) return false;
      return x.edge_id < y.edge_id;
    }

    // Indicates that linear rather than binary search should be used when this
    // type is used as the key in util::btree data structures.
    using goog_btree_prefer_linear_node_search = std::true_type;
  };

  // Convenience constructor that calls Init().
  explicit S2ClosestEdgeQueryBase(S2ShapeIndexBase const* index);

  // Default constructor; requires Init() to be called.
  S2ClosestEdgeQueryBase();
  ~S2ClosestEdgeQueryBase();

  // Initialize the query.
  // REQUIRES: Reset() must be called if "index" is modified.
  void Init(S2ShapeIndexBase const* index);

  // Reset the query state.  This method must be called whenever the
  // underlying index is modified.
  void Reset();

  // Returns a reference to the underlying S2ShapeIndex.
  S2ShapeIndexBase const& index() const;

  // Returns the closest edges to the given target that satisfy the given
  // options.  This method may be called multiple times.
  std::vector<Result> FindClosestEdges(Target* target, Options const& options);

  // This version can be more efficient when this method is called many times,
  // since it does not require allocating a new vector on each call.
  void FindClosestEdges(Target* target, Options const& options,
                        std::vector<Result>* results);

  // Convenience method that returns exactly one edge.  If no edges satisfy
  // the given search criteria, then a Result with distance == Infinity() and
  // shape_id == edge_id == -1 is returned.
  //
  // REQUIRES: options.max_edges() == 1
  Result FindClosestEdge(Target* target, Options const& options);

  ABSL_DEPRECATED("Use (Target *) version")
  std::vector<Result> FindClosestEdges(Target const& target,
                                       Options const& options) {
    return FindClosestEdges(const_cast<Target*>(&target), options);
  }

  ABSL_DEPRECATED("Use (Target *) version")
  void FindClosestEdges(Target const& target, Options const& options,
                        std::vector<Result>* results) {
    return FindClosestEdges(const_cast<Target*>(&target), options, results);
  }

  ABSL_DEPRECATED("Use (Target *) version")
  Result FindClosestEdge(Target const& target, Options const& options) {
    return FindClosestEdge(const_cast<Target*>(&target), options);
  }

 private:
  class QueueEntry;

  void FindClosestEdgesInternal(Target* target, Options const& options);
  void FindClosestEdgesBruteForce();
  void FindClosestEdgesOptimized();
  void InitQueue();
  void InitCovering();
  void AddInitialRange(S2ShapeIndexBase::Iterator const& first,
                       S2ShapeIndexBase::Iterator const& last);
  void MaybeAddResult(S2Shape const& shape, int edge_id);
  void ProcessEdges(QueueEntry const& entry);
  void EnqueueCell(S2CellId id, S2ShapeIndexCell const* index_cell);
  void EnqueueCurrentCell(S2CellId id);
  Options const& options() const { return *options_; }

  S2ShapeIndexBase const* index_;
  Options const* options_;
  Target* target_;

  // For the optimized algorihm we precompute the top-level S2CellIds that
  // will be added to the priority queue.  There can be at most 6 of these
  // cells.  Essentially this is just a covering of the indexed edges, except
  // that we also store pointers to the corresponding S2ShapeIndexCells to
  // reduce the number of index seeks required.
  //
  // The covering needs to be stored in a std::vector so that we can use
  // S2CellUnion::GetIntersection().
  int index_num_edges_;
  std::vector<S2CellId> index_covering_;
  absl::InlinedVector<S2ShapeIndexCell const*, 6> index_cells_;

  // The distance beyond which we can safely ignore further candidate edges.
  // (Candidates that are exactly at the limit are ignored; this is more
  // efficient for UpdateMinDistance() and should not affect clients since
  // distance measurements have a small amount of error anyway.)
  //
  // Initially this is the same as the maximum distance specified by the user,
  // but it can also be updated by the algorithm (see MaybeAddResult).
  Distance distance_limit_;

  // The current result set is stored in one of three ways:
  //
  //  - If max_edges() == 1, the best result is kept in result_singleton_.
  //
  //  - If max_edges() == "infinity", results are appended to result_vector_
  //    and sorted/uniqued at the end.
  //
  //  - Otherwise results are kept in a btree_set so that we can progressively
  //    reduce the distance limit once max_edges() results have been found.
  //    (A priority queue is not sufficient because we need to be able to
  //    check whether a candidate edge is already in the result set.)
  Result result_singleton_;
  std::vector<Result> result_vector_;
  util::btree::btree_set<Result> result_set_;

  // When the result edges are stored in a btree_set (see above), usually
  // duplicates can be removed simply by inserting candidate edges in the
  // current set.  However this is not true if Options::max_error() > 0 and
  // the Target subtype takes advantage of this by returning suboptimal
  // distances.  This is because when UpdateMinDistance() is called with
  // different "min_dist" parameters (i.e., the distance to beat), the
  // implementation may return a different distance for the same edge.  Since
  // the btree_set is keyed by (distance, shape_id, edge_id) this can create
  // duplicate edges in the results.
  //
  // The flag below is true when duplicates must be avoided explicitly.  This
  // is achieved by maintaining a separate set keyed by (shape_id, edge_id)
  // only, and checking whether each edge is in that set before computing the
  // distance to it.
  bool avoid_duplicates_;
  using ShapeEdgeId = s2shapeutil::ShapeEdgeId;
  google::dense_hash_set<
    ShapeEdgeId, s2shapeutil::ShapeEdgeIdHash> tested_edges_;

  // The algorithm maintains a priority queue of S2CellIds that contain at
  // least one S2ShapeIndexCell, sorted in increasing order of distance from
  // the target point.
  struct QueueEntry {
    // A lower bound on the distance from the target point to any edge point
    // within "id".  This is the key of the priority queue.
    Distance distance;

    // The cell being queued.
    S2CellId id;

    // If "id" belongs to the index, this field stores the corresponding
    // S2ShapeIndexCell.  Otherwise "id" is a proper ancestor of one or more
    // S2ShapeIndexCells and this field stores nullptr.  The purpose of this
    // field is to avoid an extra Seek() when the queue entry is processed.
    S2ShapeIndexCell const* index_cell;

    QueueEntry(Distance _distance, S2CellId _id,
               S2ShapeIndexCell const* _index_cell)
        : distance(_distance), id(_id), index_cell(_index_cell) {
    }
    bool operator<(QueueEntry const& other) const {
      // The priority queue returns the largest elements first, so we want the
      // "largest" entry to have the smallest distance.
      return other.distance < distance;
    }
  };
  using CellQueue =
      std::priority_queue<QueueEntry, absl::InlinedVector<QueueEntry, 16>>;
  CellQueue queue_;

  // Temporaries, defined here to avoid multiple allocations / initializations.

  S2ShapeIndexBase::Iterator iter_;
  std::vector<S2CellId> max_distance_covering_;
  std::vector<S2CellId> initial_cells_;

  S2ClosestEdgeQueryBase(S2ClosestEdgeQueryBase const&) = delete;
  void operator=(S2ClosestEdgeQueryBase const&) = delete;
};


//////////////////   Implementation details follow   ////////////////////


template <class Distance>
inline S2ClosestEdgeQueryBase<Distance>::Options::Options() {
}

template <class Distance>
inline int S2ClosestEdgeQueryBase<Distance>::Options::max_edges() const {
  return max_edges_;
}

template <class Distance>
inline void S2ClosestEdgeQueryBase<Distance>::Options::set_max_edges(
    int max_edges) {
  DCHECK_GE(max_edges, 1);
  max_edges_ = max_edges;
}

template <class Distance>
inline Distance S2ClosestEdgeQueryBase<Distance>::Options::max_distance()
    const {
  return max_distance_;
}

template <class Distance>
inline void S2ClosestEdgeQueryBase<Distance>::Options::set_max_distance(
    Distance max_distance) {
  max_distance_ = max_distance;
}

template <class Distance>
inline Distance S2ClosestEdgeQueryBase<Distance>::Options::max_error() const {
  return max_error_;
}

template <class Distance>
inline void S2ClosestEdgeQueryBase<Distance>::Options::set_max_error(
    Distance max_error) {
  max_error_ = max_error;
}

template <class Distance>
inline bool S2ClosestEdgeQueryBase<Distance>::Options::include_interiors()
    const {
  return include_interiors_;
}

template <class Distance>
inline void S2ClosestEdgeQueryBase<Distance>::Options::set_include_interiors(
    bool include_interiors) {
  include_interiors_ = include_interiors;
}

template <class Distance>
inline bool S2ClosestEdgeQueryBase<Distance>::Options::use_brute_force() const {
  return use_brute_force_;
}

template <class Distance>
inline void S2ClosestEdgeQueryBase<Distance>::Options::set_use_brute_force(
    bool use_brute_force) {
  use_brute_force_ = use_brute_force;
}

template <class Distance>
S2ClosestEdgeQueryBase<Distance>::S2ClosestEdgeQueryBase() {
  // Prevent inline constructor bloat by providing a definition.
}

template <class Distance>
S2ClosestEdgeQueryBase<Distance>::~S2ClosestEdgeQueryBase() {
  // Prevent inline destructor bloat by providing a definition.
}

template <class Distance>
inline S2ClosestEdgeQueryBase<Distance>::S2ClosestEdgeQueryBase(
    S2ShapeIndexBase const* index) {
  Init(index);
}

template <class Distance>
void S2ClosestEdgeQueryBase<Distance>::Init(S2ShapeIndexBase const* index) {
  index_ = index;
  Reset();
}

template <class Distance>
void S2ClosestEdgeQueryBase<Distance>::Reset() {
  index_num_edges_ = s2shapeutil::GetNumEdges(*index_);
  index_covering_.clear();
  index_cells_.clear();
}

template <class Distance>
inline S2ShapeIndexBase const& S2ClosestEdgeQueryBase<Distance>::index() const {
  return *index_;
}

template <class Distance>
inline std::vector<typename S2ClosestEdgeQueryBase<Distance>::Result>
S2ClosestEdgeQueryBase<Distance>::FindClosestEdges(Target* target,
                                                   Options const& options) {
  std::vector<Result> results;
  FindClosestEdges(target, options, &results);
  return results;
}

template <class Distance>
typename S2ClosestEdgeQueryBase<Distance>::Result
S2ClosestEdgeQueryBase<Distance>::FindClosestEdge(Target* target,
                                                  Options const& options) {
  DCHECK_EQ(options.max_edges(), 1);
  FindClosestEdgesInternal(target, options);
  return result_singleton_;
}

template <class Distance>
void S2ClosestEdgeQueryBase<Distance>::FindClosestEdges(
    Target* target, Options const& options,
    std::vector<Result>* results) {
  FindClosestEdgesInternal(target, options);
  results->clear();
  if (options.max_edges() == 1) {
    if (result_singleton_.shape_id >= 0) {
      results->push_back(result_singleton_);
    }
  } else if (options.max_edges() == Options::kMaxMaxEdges) {
    std::sort(result_vector_.begin(), result_vector_.end());
    std::unique_copy(result_vector_.begin(), result_vector_.end(),
                     std::back_inserter(*results));
    result_vector_.clear();
  } else {
    results->assign(result_set_.begin(), result_set_.end());
    result_set_.clear();
  }
}

template <class Distance>
void S2ClosestEdgeQueryBase<Distance>::FindClosestEdgesInternal(
    Target* target, Options const& options) {
  target_ = target;
  options_ = &options;
  distance_limit_ = options.max_distance();
  result_singleton_ = Result();
  DCHECK(result_vector_.empty());
  DCHECK(result_set_.empty());

  // If max_error() was specified and the target takes advantage of this
  // in its UpdateMinDistance() methods, then we need to avoid duplicate edges
  // in the results explicitly.  (Otherwise it happens automatically.)
  avoid_duplicates_ = (
      Distance::Zero() < options.max_error() &&
      target_->set_max_error(options.max_error()) &&
      options.max_edges() > 1 && !options.use_brute_force());

  if (options.use_brute_force() ||
      index_num_edges_ <= target_->max_brute_force_edges()) {
    FindClosestEdgesBruteForce();
  } else {
    FindClosestEdgesOptimized();
  }
}

template <class Distance>
void S2ClosestEdgeQueryBase<Distance>::FindClosestEdgesBruteForce() {
  int num_shape_ids = index_->num_shape_ids();
  for (int id = 0; id < num_shape_ids; ++id) {
    S2Shape const* shape = index_->shape(id);
    if (shape == nullptr) continue;
    int num_edges = shape->num_edges();
    for (int e = 0; e < num_edges; ++e) {
      MaybeAddResult(*shape, e);
    }
  }
}

template <class Distance>
void S2ClosestEdgeQueryBase<Distance>::FindClosestEdgesOptimized() {
  InitQueue();
  // Repeatedly find the closest S2Cell to "target" and either split it into
  // its four children or process all of its edges.
  while (!queue_.empty()) {
    // We need to copy the top entry before removing it, and we need to
    // remove it before adding any new entries to the queue.
    QueueEntry entry = queue_.top();
    queue_.pop();
    if (!(entry.distance < distance_limit_)) {
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
    if (!iter_.done() && iter_.id() <= id.child(1).range_max()) {
      EnqueueCurrentCell(id.child(1));
    }
    if (iter_.Prev() && iter_.id() >= id.range_min()) {
      EnqueueCurrentCell(id.child(0));
    }
    iter_.Seek(id.child(3).range_min());
    if (!iter_.done() && iter_.id() <= id.range_max()) {
      EnqueueCurrentCell(id.child(3));
    }
    if (iter_.Prev() && iter_.id() >= id.child(2).range_min()) {
      EnqueueCurrentCell(id.child(2));
    }
  }
}

template <class Distance>
void S2ClosestEdgeQueryBase<Distance>::InitQueue() {
  DCHECK(queue_.empty());
  iter_.Init(index_, S2ShapeIndex::UNPOSITIONED);

  // Optimization: if the user is searching for just the closest edge, and the
  // target happens to intersect an index cell, then we try to limit the search
  // region to a small disc by first processing the edges in that cell.  This
  // sets distance_limit_ based on the closest edge in that cell, which we can
  // then use to limit the search area.  This means that the cell containing
  // "target" will be processed twice, but in general this is still faster.
  S2Cap cap = target_->GetCapBound();
  if (options().max_edges() == 1 && iter_.Locate(cap.center())) {
    ProcessEdges(QueueEntry(Distance::Zero(), iter_.id(), &iter_.cell()));
    // Skip the rest of the algorithm if we found an intersecting edge.
    if (distance_limit_ == Distance::Zero()) return;
  }
  if (index_covering_.empty()) InitCovering();
  if (distance_limit_ == Distance::Infinity()) {
    // Start with the precomputed index covering.
    for (int i = 0; i < index_covering_.size(); ++i) {
      EnqueueCell(index_covering_[i], index_cells_[i]);
    }
  } else {
    // Compute a covering of the search disc and intersect it with the
    // precomputed index covering.
    S2RegionCoverer coverer;
    coverer.set_max_cells(4);
    S1Angle radius = cap.GetRadius() + Distance::GetAngleBound(distance_limit_);
    S2Cap search_cap(cap.center(), radius);
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
          EnqueueCell(iter_.id(), &iter_.cell());
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

template <class Distance>
void S2ClosestEdgeQueryBase<Distance>::InitCovering() {
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

  index_covering_.reserve(6);
  // Don't need to reserve index_cells_ since it is an InlinedVector.
  S2ShapeIndexBase::Iterator next(index_, S2ShapeIndex::BEGIN);
  S2ShapeIndexBase::Iterator last(index_, S2ShapeIndex::END);
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
      S2ShapeIndexBase::Iterator cell_first = next;
      next.Seek(id.range_max().next());
      S2ShapeIndexBase::Iterator cell_last = next;
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
template <class Distance>
void S2ClosestEdgeQueryBase<Distance>::AddInitialRange(
    S2ShapeIndexBase::Iterator const& first,
    S2ShapeIndexBase::Iterator const& last) {
  if (first.id() == last.id()) {
    // The range consists of a single index cell.
    index_covering_.push_back(first.id());
    index_cells_.push_back(&first.cell());
  } else {
    // Add the lowest common ancestor of the given range.
    int level = first.id().GetCommonAncestorLevel(last.id());
    DCHECK_GE(level, 0);
    index_covering_.push_back(first.id().parent(level));
    index_cells_.push_back(nullptr);
  }
}

template <class Distance>
void S2ClosestEdgeQueryBase<Distance>::MaybeAddResult(
    S2Shape const& shape, int edge_id) {
  if (avoid_duplicates_ &&
      !tested_edges_.insert(ShapeEdgeId(shape.id(), edge_id)).second) {
    return;
  }
  auto edge = shape.edge(edge_id);
  Distance distance = distance_limit_;
  if (!target_->UpdateMinDistance(edge.v0, edge.v1, &distance)) return;

  Result result(distance, shape.id(), edge_id);
  if (options().max_edges() == 1) {
    // Optimization for the common case where only the closest edge is wanted.
    result_singleton_ = result;
    distance_limit_ = distance - options().max_error();
  } else if (options().max_edges() == Options::kMaxMaxEdges) {
    result_vector_.push_back(result);  // Sort/unique at end.
  } else {
    // Add this edge to result_set_.  Note that even if we already have enough
    // edges, we can't erase an element before insertion because the "new"
    // edge might in fact be a duplicate.
    result_set_.insert(result);
    int size = result_set_.size();
    if (size >= options().max_edges()) {
      if (size > options().max_edges()) {
        result_set_.erase(--result_set_.end());
      }
      distance_limit_ = (--result_set_.end())->distance - options().max_error();
    }
  }
}

// Return the number of edges in the given index cell.
inline static int CountEdges(S2ShapeIndexCell const* cell) {
  int count = 0;
  for (int s = 0; s < cell->num_clipped(); ++s) {
    count += cell->clipped(s).num_edges();
  }
  return count;
}

// Process all the edges of the given index cell.
template <class Distance>
void S2ClosestEdgeQueryBase<Distance>::ProcessEdges(QueueEntry const& entry) {
  S2ShapeIndexCell const* index_cell = entry.index_cell;
  for (int s = 0; s < index_cell->num_clipped(); ++s) {
    S2ClippedShape const& clipped = index_cell->clipped(s);
    S2Shape const* shape = index_->shape(clipped.shape_id());
    for (int j = 0; j < clipped.num_edges(); ++j) {
      MaybeAddResult(*shape, clipped.edge(j));
    }
  }
}

// Enqueue the given cell id.
// REQUIRES: iter_ is positioned at a cell contained by "id".
template <class Distance>
inline void S2ClosestEdgeQueryBase<Distance>::EnqueueCurrentCell(S2CellId id) {
  DCHECK(id.contains(iter_.id()));
  if (iter_.id() == id) {
    EnqueueCell(id, &iter_.cell());
  } else {
    EnqueueCell(id, nullptr);
  }
}

// Add the given cell id to the queue.  "index_cell" is the corresponding
// S2ShapeIndexCell, or nullptr if "id" is not an index cell.
template <class Distance>
void S2ClosestEdgeQueryBase<Distance>::EnqueueCell(
    S2CellId id, S2ShapeIndexCell const* index_cell) {
  if (index_cell) {
    // If this index cell has only a few edges, then it is faster to check
    // them directly rather than computing the minimum distance to the S2Cell
    // and inserting it into the queue.
    static int const kMinEdgesToEnqueue = 10;
    int num_edges = CountEdges(index_cell);
    if (num_edges == 0) return;
    if (num_edges < kMinEdgesToEnqueue) {
      // Set "distance" to zero to avoid the expense of computing it.
      ProcessEdges(QueueEntry(Distance::Zero(), id, index_cell));
      return;
    }
  }
  // Otherwise compute the minimum distance to any point in the cell and add
  // it to the priority queue.
  S2Cell cell(id);
  Distance distance = distance_limit_;
  if (!target_->UpdateMinDistance(cell, &distance)) return;
  queue_.push(QueueEntry(distance, id, index_cell));
}

#endif  // S2_S2CLOSESTEDGEQUERY_BASE_H_
