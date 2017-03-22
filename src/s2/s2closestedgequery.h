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

#ifndef S2_S2CLOSESTEDGEQUERY_H_
#define S2_S2CLOSESTEDGEQUERY_H_

#include <memory>
#include <type_traits>
#include <vector>

#include <glog/logging.h>
#include "s2/third_party/absl/container/inlined_vector.h"
#include "s2/util/btree/btree_set.h"  // Like std::set, but faster and smaller.
#include "s2/fpcontractoff.h"
#include "s2/priority_queue_sequence.h"
#include "s2/s1angle.h"
#include "s2/s1chordangle.h"
#include "s2/s2cell.h"
#include "s2/s2cellid.h"
#include "s2/s2edgeutil.h"
#include "s2/s2shapeindex.h"

// S2ClosestEdgeQuery is a helper class for finding the closest edge(s) to a
// given query point or query edge.  For example, given a set of polylines,
// the following code efficiently finds the closest 5 edges to a query point:
//
// void Test(vector<S2Polyline*> const& polylines, S2Point const& target) {
//   S2ShapeIndex index;
//   for (S2Polyline* polyline : polylines) {
//     index.Add(new S2Polyline::Shape(polyline));
//   }
//   S2ClosestEdgeQuery query(index);
//   query.set_max_edges(5);
//   query.FindClosestEdges(target);
//   for (int i = 0; i < query.num_edges(); ++i) {
//     // query.shape_id(i) and query.edge_id(i) identify the edge.
//     // Convenience methods:
//     //   query.distance(i) is the distance to the target point.
//     //   query.GetEdge(i, &v0, &v1) retrieves the edge endpoints.
//     //   query.GetClosestPointOnEdge(i) returns the point on the edge
//     //     that is closest to the target point.
//     int polyline_index = query.shape_id(i);
//     int edge_index = query.edge_id(i);
//     S1Angle distance = query.distance(i);
//     S2Point closest_point = query.GetClosestPointOnEdge(i);
//   }
// }
//
// You can find either the k closest edges, or all edges within a given
// radius, or both (i.e., the k closest edges up to a given maximum radius).
// E.g. to find all the edges within 5 kilometers, call
//
//   query.set_max_distance(S2Earth::ToAngle(util::units::Kilometers(5)));
//
// To find the closest points to a query edge rather than a point, use:
//
//   query.FindClosestEdgesToEdge(v0, v1);
//
// The implementation is designed to be fast for both small and large sets of
// edges.
class S2ClosestEdgeQuery {
 public:
  // Convenience constructor that calls Init().
  explicit S2ClosestEdgeQuery(S2ShapeIndex const& index);

  // Default constructor; requires Init() to be called.
  S2ClosestEdgeQuery() {}
  ~S2ClosestEdgeQuery();

  // Initialize the query.
  // REQUIRES: Reset() must be called if "index" is modified.
  void Init(S2ShapeIndex const& index);

  // Reset the query state.  This method must be called whenever the
  // underlying index is modified.
  void Reset();

  // Return a reference to the underlying S2ShapeIndex.
  S2ShapeIndex const& index() const;

  // Only find the "max_edges" closest edges.
  // This value may be changed between calls to FindClosestEdges().
  //
  // Default value: numeric_limits<int>::max()
  // REQUIRES: max_edges >= 1
  int max_edges() const;
  void set_max_edges(int max_edges);

  // Find edges whose distance to the target is less than "max_distance".
  // (All edges whose true distance is less than the given value will be
  // returned, along with some edges whose true distance may be slightly
  // greater than the limit.)  This value may be changed between calls to
  // FindClosestEdges().
  //
  // Default value: S1Angle::Infinity().
  S1Angle max_distance() const;
  void set_max_distance(S1Angle max_distance);

  // Allow distance errors of up to "max_error" when determining which edge(s)
  // are closest.  Note that this does not affect how distances are measured,
  // it just gives the algorithm permission to not work so hard to find the
  // exact closest edge (or set of edges) when many edges are nearly equally
  // distant.  Specifically, if E1 is some edge in the "exact" result set,
  // then the algorithm is allowed to substitute an edge E2 that is up to
  // "max_error" further away than E1 (and that satisfies all the remaining
  // search criteria, including max_distance).
  //
  // Default value: S1Angle::Zero().
  S1Angle max_error() const;
  void set_max_error(S1Angle max_error);

  // Find the closest edges to "target" that satisfy the given max_distance()
  // and/or max_edges() criteria.  If neither of these is set, then all edges
  // are returned.  This method may be called multiple times.
  //
  // REQUIRES: "target" is unit length.
  void FindClosestEdges(S2Point const& target);

  // Find the closest edges to the given edge AB.  Otherwise similar to
  // FindClosestEdgess().
  void FindClosestEdgesToEdge(S2Point const& a, S2Point const& b);

  // Manually specify whether distances are computed using "brute force"
  // (i.e., by examining every edge) rather than using the S2ShapeIndex.
  // This is useful for testing, benchmarking, and debugging.
  //
  // REQUIRES: Init() has been called.
  void UseBruteForce(bool use_brute_force);

  ////////////////////// Result Interface ////////////////////////
  //
  // The result of a query consists of a set of edges which may be obtained
  // using the interface below.  Edges are sorted in order of increasing
  // distance to the target point.

  // The number of result edges.
  int num_edges() const;

  // The shape id of the given result edge.
  int shape_id(int i) const;

  // The edge id of the given result edge.
  int edge_id(int i) const;

  // The distance to the given result edge.
  S1Angle distance(int i) const;

  // Like distance(i), but expressed as an S1ChordAngle.
  S1ChordAngle distance_ca(int i) const;

  // Returns pointers to the vertices of the given result edge.  Example usage:
  //   S2Point const *v0, *v1;
  //   query.GetEdge(i, &v0, &v1);
  void GetEdge(int i, S2Point const** v0, S2Point const** v1) const;

  // Returns the point on given result edge that is closest to the target.
  S2Point GetClosestPointOnEdge(int i) const;

  //////////////////////// Convenience Methods ////////////////////////

  // Finds the closest edge to the target.  Equivalent to:
  //
  //   query.set_max_edges(1);
  //   query.FindClosestEdges();
  //
  // All other parameters are unchanged; e.g. max_distance() will still limit
  // the search radius.
  //
  // REQUIRES: "target" is unit length.
  void FindClosestEdge(S2Point const& target);

  // Return the minimum distance from the given point to any edge of the
  // S2ShapeIndex.  If there are no edges, return S1Angle::Infinity().
  //
  // SIDE EFFECTS: Calls set_max_edges(1); all other parameters are unchanged.
  S1Angle GetDistance(S2Point const& target);

  // Return the closest point on any edge of the S2ShapeIndex to the given
  // point.  If the index has no edges, return the input argument.
  //
  // SIDE EFFECTS: Calls set_max_edges(1); all other parameters are unchanged.
  S2Point Project(S2Point const& target);

 private:
  class QueueEntry;

  void AddInitialRange(S2ShapeIndex::Iterator const& first,
                       S2ShapeIndex::Iterator const& last);

  // TODO(ericv): Should the Target classes be factored out somewhere so that
  // they can be shared with the similar classes in S2ClosestPointQuery?
  class Target {
   public:
    virtual ~Target() {}
    virtual S2Point center() const = 0;
    virtual S1Angle radius() const = 0;
    virtual bool UpdateMinDistance(S2Point const& v0, S2Point const& v1,
                                   S1ChordAngle* min_dist) const = 0;
    virtual S1ChordAngle GetDistance(S2Cell const& cell) const = 0;
    virtual S2Point GetClosestPointOnEdge(S2Point const& v0,
                                          S2Point const& v1) const = 0;
  };

  class PointTarget : public Target {
   public:
    explicit PointTarget(S2Point const& point) : point_(point) {}
    S2Point center() const override { return point_; }
    S1Angle radius() const override { return S1Angle::Zero(); }
    bool UpdateMinDistance(S2Point const& v0, S2Point const& v1,
                           S1ChordAngle* min_dist) const override {
      return S2EdgeUtil::UpdateMinDistance(point_, v0, v1, min_dist);
    }
    S1ChordAngle GetDistance(S2Cell const& cell) const override {
      return cell.GetDistance(point_);
    }
    S2Point GetClosestPointOnEdge(S2Point const& v0,
                                  S2Point const& v1) const override {
      return S2EdgeUtil::GetClosestPoint(point_, v0, v1);
    }
   private:
    S2Point point_;
  };

  class EdgeTarget : public Target {
   public:
    EdgeTarget(S2Point const& a, S2Point const& b) : a_(a), b_(b) {}
    S2Point center() const override { return (a_ + b_).Normalize(); }
    S1Angle radius() const override { return 0.5 * S1Angle(a_, b_); }
    bool UpdateMinDistance(S2Point const& v0, S2Point const& v1,
                           S1ChordAngle* min_dist) const override {
      return S2EdgeUtil::UpdateEdgePairMinDistance(a_, b_, v0, v1, min_dist);
    }
    S1ChordAngle GetDistance(S2Cell const& cell) const override {
      return cell.GetDistanceToEdge(a_, b_);
    }
    S2Point GetClosestPointOnEdge(S2Point const& v0,
                                  S2Point const& v1) const override {
      return S2EdgeUtil::GetEdgePairClosestPoints(a_, b_, v0, v1).second;
    }
   private:
    S2Point a_, b_;
  };

  // Using templates rather than virtual functions speeds up the benchmarks
  // by 15-25% when the brute force algorithm is used (< 150 points).

  void FindClosestEdgesToTarget();
  void FindClosestEdgesBruteForce();
  void FindClosestEdgesOptimized();
  void InitQueue();
  void MaybeAddResult(S2Shape const& shape, int edge_id);
  void ProcessEdges(QueueEntry const& entry);
  void EnqueueCell(S2CellId id, S2ShapeIndexCell const* index_cell);
  void EnqueueCurrentCell(S2CellId id);

  //////////// Parameters /////////////

  int max_edges_;
  S1Angle max_distance_;
  S1Angle max_error_arg_;
  std::unique_ptr<const Target> target_;

  ////////// Fields that are constant after Init() is called /////////////

  S2ShapeIndex const* index_;

  // If the index has few edges, it is cheaper to use a brute force algorithm.
  bool use_brute_force_;

  // During Init() we precompute the top-level S2CellIds that will be added to
  // the priority queue.  There can be at most 6 of these cells.  Essentially
  // this is just a covering of the indexed edges, except that we also store
  // pointers to the corresponding S2ShapeIndexCells to reduce the number of
  // index seeks required.
  //
  // The covering needs to be stored in a std::vector so that we can use
  // S2CellUnion::GetIntersection().
  std::vector<S2CellId> index_covering_;
  absl::InlinedVector<S2ShapeIndexCell const*, 6> index_cells_;

  ////////// Fields that are updated during a query /////////////

  // The edges gathered so far are kept in a sorted container so that
  // duplicate edges can easily be removed.
  struct Result {
    // Default constructor needed by standard containers.
    Result() : distance(S1ChordAngle::Negative()), shape_id(-1), edge_id(-1) {}
    Result(S1ChordAngle _distance, int _shape_id, int _edge_id)
        : distance(_distance), shape_id(_shape_id), edge_id(_edge_id) {
    }
    bool operator<(Result const& other) const {
      // Edges are sorted first by distance and then by unique id.
      if (distance < other.distance) return true;
      if (distance > other.distance) return false;
      if (shape_id < other.shape_id) return true;
      if (shape_id > other.shape_id) return false;
      return edge_id < other.edge_id;
    }
    // Indicate that linear rather than binary search should be used within
    // btree nodes.  This is faster even for the comparison function above
    // since the result is determined by "distance" most of the time.
    using goog_btree_prefer_linear_node_search = std::true_type;

    S1ChordAngle distance;
    int shape_id, edge_id;
  };
  using ResultSet = util::btree::btree_set<Result>;
  ResultSet tmp_results_;

  // For efficiency, when max_edges() == 1 we keep the current best result in
  // this field rather than using "tmp_results_" above.
  Result tmp_result_singleton_;

  // Once all result edges have been found, they are copied into a vector.
  absl::InlinedVector<Result, 8> results_;

  Result const& result(int i) const;  // Internal accessor method.

  // The distance beyond which we can safely ignore further candidate edges.
  // (Candidates that are exactly at the limit are ignored; this is more
  // efficient for UpdateMinDistance() and should not affect clients since
  // distance measurements have a small amount of error anyway.)
  //
  // Initially this is the same as the maximum distance specified by the user,
  // but it can also be updated by the algorithm (see MaybeAddResult).
  S1ChordAngle max_distance_limit_;

  // Internally we convert max_error_arg_ to an S1ChordAngle for efficiency.
  S1ChordAngle max_error_;

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
    // S2ShapeIndexCells and this field stores nullptr.  The purpose of this
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
  using CellQueue =
      std::priority_queue<QueueEntry, absl::InlinedVector<QueueEntry, 16>>;
  CellQueue queue_;

  // Temporaries, defined here to avoid multiple allocations / initializations.

  S2ShapeIndex::Iterator iter_;
  std::vector<S2CellId> max_distance_covering_;
  std::vector<S2CellId> initial_cells_;

  S2ClosestEdgeQuery(S2ClosestEdgeQuery const&) = delete;
  void operator=(S2ClosestEdgeQuery const&) = delete;
};


//////////////////   Implementation details follow   ////////////////////


inline S2ClosestEdgeQuery::S2ClosestEdgeQuery(S2ShapeIndex const& index) {
  Init(index);
}

inline S2ShapeIndex const& S2ClosestEdgeQuery::index() const {
  return *index_;
}

inline int S2ClosestEdgeQuery::max_edges() const {
  return max_edges_;
}

inline void S2ClosestEdgeQuery::set_max_edges(int max_edges) {
  DCHECK_GE(max_edges, 1);
  max_edges_ = max_edges;
}

inline S1Angle S2ClosestEdgeQuery::max_distance() const {
  return max_distance_;
}

inline void S2ClosestEdgeQuery::set_max_distance(S1Angle max_distance) {
  max_distance_ = max_distance;
}

inline S1Angle S2ClosestEdgeQuery::max_error() const {
  return max_error_arg_;
}

inline void S2ClosestEdgeQuery::set_max_error(S1Angle max_error) {
  max_error_arg_ = max_error;
}

inline int S2ClosestEdgeQuery::num_edges() const {
  return results_.size();
}

inline int S2ClosestEdgeQuery::shape_id(int i) const {
  return results_[i].shape_id;
}

inline int S2ClosestEdgeQuery::edge_id(int i) const {
  return results_[i].edge_id;
}

inline S1ChordAngle S2ClosestEdgeQuery::distance_ca(int i) const {
  return results_[i].distance;
}

inline S1Angle S2ClosestEdgeQuery::distance(int i) const {
  return distance_ca(i).ToAngle();
}

inline void S2ClosestEdgeQuery::GetEdge(int i, S2Point const** v0,
                                        S2Point const** v1) const {
  index_->shape(shape_id(i))->GetEdge(edge_id(i), v0, v1);
}

inline void S2ClosestEdgeQuery::FindClosestEdge(S2Point const& target) {
  set_max_edges(1);
  FindClosestEdges(target);
}

#endif  // S2_S2CLOSESTEDGEQUERY_H_
