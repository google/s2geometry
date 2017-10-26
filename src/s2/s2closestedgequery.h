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
#include "s2/third_party/absl/base/macros.h"
#include "s2/third_party/absl/container/inlined_vector.h"
#include "s2/util/btree/btree_set.h"
#include "s2/_fpcontractoff.h"
#include "s2/priority_queue_sequence.h"
#include "s2/s1angle.h"
#include "s2/s1chordangle.h"
#include "s2/s2cell.h"
#include "s2/s2cellid.h"
#include "s2/s2closestedgequery_base.h"
#include "s2/s2edge_distances.h"
#include "s2/s2shapeindex.h"

// S2ClosestEdgeQuery is a helper class for finding the closest edge(s) to a
// given point, edge, S2Cell, or geometry collection.  For example, given a
// set of polylines, the following code efficiently finds the closest 5 edges
// to a query point:
//
// void Test(const vector<S2Polyline*>& polylines, const S2Point& point) {
//   MutableS2ShapeIndex index;
//   for (S2Polyline* polyline : polylines) {
//     index.Add(new S2Polyline::Shape(polyline));
//   }
//   S2ClosestEdgeQuery query(&index);
//   query.mutable_options()->set_max_edges(5);
//   S2ClosestEdgeQuery::PointTarget target(point);
//   for (const auto& result : query.FindClosestEdges(&target)) {
//     // The Result struct contains the following fields:
//     //   "distance" is the distance to the edge.
//     //   "shape_id" identifies the S2Shape containing the edge.
//     //   "edge_id" identifies the edge with the given shape.
//     // The following convenience methods may also be useful:
//     //   query.GetEdge(result) returns the endpoints of the edge.
//     //   query.Project(point, result) computes the closest point on the
//     //       result edge to the given target point.
//     int polyline_index = result.shape_id;
//     int edge_index = result.edge_id;
//     S1ChordAngle distance = result.distance;  // Use ToAngle() for S1Angle.
//     S2Shape::Edge edge = query.GetEdge(result);
//     S2Point closest_point = query.Project(point, result);
//   }
// }
//
// You can find either the k closest edges, or all edges within a given
// radius, or both (i.e., the k closest edges up to a given maximum radius).
// E.g. to find all the edges within 5 kilometers, call
//
//   query.mutable_options()->set_max_distance(
//       S2Earth::ToAngle(util::units::Kilometers(5)));
//
// By default *all* edges are returned, so you should always specify either
// max_edges() or max_distance() or both.  There is also a FindClosestEdge()
// convenience method that automatically sets max_edges() == 1 and returns
// only the closest edge.
//
// Note that by default, distances are measured to the boundary and interior
// of polygons.  For example, if a point is inside a polygon then its distance
// is zero.  To change this behavior, call set_include_interiors(false).
//
// If you only need to test whether the distance is above or below a given
// threshold (e.g., 10 km), you can use the IsDistanceLess() method.  This is
// much faster than actually calculating the distance with FindClosestEdge(),
// since the implementation can stop as soon as it can prove that the minimum
// distance is either above or below the threshold.
//
// To find the closest edges to a query edge rather than a point, use:
//
//   S2ClosestEdgeQuery::EdgeTarget target(v0, v1);
//   query.FindClosestEdges(&target);
//
// Similarly you can find the closest edges to an S2Cell by using an
// S2ClosestEdgeQuery::CellTarget, and you can find the closest edges to an
// arbitrary collection of points, polylines, and polygons by using an
// S2ClosestEdgeQuery::ShapeIndexTarget.
//
// The implementation is designed to be fast for both simple and complex
// geometric objects.
class S2ClosestEdgeQuery {
 public:
  // See S2ClosestEdgeQueryBase for full documentation.

  // A thin wrapper around S1ChordAngle that implements the Distance concept
  // required by S2ClosestEdgeQueryBase.
  class Distance : public S1ChordAngle {
   public:
    Distance() : S1ChordAngle() {}
    explicit Distance(S1Angle x) : S1ChordAngle(x) {}
    explicit Distance(S1ChordAngle x) : S1ChordAngle(x) {}
    static Distance Zero() { return Distance(S1ChordAngle::Zero()); }
    static Distance Infinity() { return Distance(S1ChordAngle::Infinity()); }
    static Distance Negative() { return Distance(S1ChordAngle::Negative()); }
    friend Distance operator-(Distance x, Distance y) {
      return Distance(S1ChordAngle(x) - y);
    }
    static S1Angle GetAngleBound(Distance x) {
      return x.PlusError(x.GetS1AngleConstructorMaxError()).ToAngle();
    }
    // If (dist < *this), updates *this and returns true.
    bool UpdateMin(const Distance& dist);
  };
  using Base = S2ClosestEdgeQueryBase<Distance>;
  using Result = Base::Result;

  // See S2ClosestEdgeQueryBase for full documentation of the available options.
  class Options : public Base::Options {
   public:
    // Versions of set_max_distance() that accept S1ChordAngle / S1Angle.
    //
    // Note that only edges whose distance is *less than* "max_distance" are
    // returned.  Normally this doesn't matter, because distances are not
    // computed exactly in the first place, but if such edges are needed then
    // you can retrieve them by specifying max_distance.Successor() instead.
    void set_max_distance(S1ChordAngle max_distance);
    void set_max_distance(S1Angle max_distance);

    // Versions of set_max_error() that accept S1ChordAngle / S1Angle.
    void set_max_error(S1ChordAngle max_error);
    void set_max_error(S1Angle max_error);

    // Like set_max_distance(), except that "max_distance" is increased by the
    // maximum error in the distance calculation.  This ensures that all edges
    // whose true distance is less than "max_distance" will be returned (along
    // with some edges whose true distance is slightly greater).
    //
    // Algorithms that need to do exact distance comparisons can use this
    // option to find a set of candidate edges that can then be filtered
    // further (e.g., using s2pred::CompareEdgeDistance).
    void set_conservative_max_distance(S1ChordAngle max_distance);
  };

  // "Target" represents the geometry that the distance is measured to.  There
  // are subtypes for measuring the distance to a point, an edge, an S2Cell,
  // or an S2ShapeIndex representing an arbitrary collection of geometry.
  using Target = Base::Target;

  // Target subtype that computes the closest distance to a point.
  class PointTarget final : public Target {
   public:
    explicit PointTarget(const S2Point& point);
    int max_brute_force_edges() const override;
    S2Cap GetCapBound() const override;
    bool UpdateMinDistance(const S2Point& v0, const S2Point& v1,
                           Distance* min_dist) const override;
    bool UpdateMinDistance(const S2Cell& cell,
                           Distance* min_dist) const override;
    std::vector<int> GetContainingShapes(const S2ShapeIndex& index,
                                         int max_shapes) const override;

   private:
    S2Point point_;
  };

  // Target subtype that computes the closest distance to an edge.
  class EdgeTarget final : public Target {
   public:
    EdgeTarget(const S2Point& a, const S2Point& b);
    int max_brute_force_edges() const override;
    S2Cap GetCapBound() const override;
    bool UpdateMinDistance(const S2Point& v0, const S2Point& v1,
                           Distance* min_dist) const override;
    bool UpdateMinDistance(const S2Cell& cell,
                           Distance* min_dist) const override;
    std::vector<int> GetContainingShapes(const S2ShapeIndex& index,
                                         int max_shapes) const override;

   private:
    S2Point a_, b_;
  };

  // Target subtype that computes the closest distance to an S2Cell
  // (including the interior of the cell).
  class CellTarget final : public Target {
   public:
    explicit CellTarget(const S2Cell& cell);
    int max_brute_force_edges() const override;
    S2Cap GetCapBound() const override;
    bool UpdateMinDistance(const S2Point& v0, const S2Point& v1,
                           Distance* min_dist) const override;
    bool UpdateMinDistance(const S2Cell& cell,
                           Distance* min_dist) const override;
    std::vector<int> GetContainingShapes(const S2ShapeIndex& index,
                                         int max_shapes) const override;

   private:
    S2Cell cell_;
  };

  // Target subtype that computes the closest distance to an S2ShapeIndex
  // (an arbitrary collection of points, polylines, and/or polygons).
  //
  // Note that ShapeIndexTarget has its own options:
  //
  //   include_interiors()
  //     - specifies that distances are measured to the boundary and interior
  //       of polygons in the S2ShapeIndex.  (If set to false, distance is
  //       measured to the polygon boundary only.)
  //       DEFAULT: true.
  //
  //   brute_force()
  //     - specifies that the distances should be computed by examining every
  //       edge in the S2ShapeIndex (for testing and debugging purposes).
  //       DEFAULT: false.
  //
  // These options are specified independently of the corresponding
  // S2ClosestEdgeQuery options.  For example, if include_interiors is true
  // for a ShapeIndexTarget but false for the S2ClosestEdgeQuery where the
  // target is used, then distances will be measured from the boundary of one
  // S2ShapeIndex to the boundary and interior of the other.
  //
  // The remaining S2ClosestEdgeQuery::Options are instead handled as follows:
  //
  //  - max_error() is copied from the current S2ClosestEdgeQuery, i.e. if you
  //    set query.options().max_error() then this value is automatically
  //    propagated to the ShapeIndexTarget.
  //
  //    Note that unlike the other Target subtypes, this option can affect the
  //    "distance" field of the results.  Suppose that max_edges() == 1 and
  //    max_error() == 0.01, and let the result edge be E with "distance"
  //    field d.  Then the implementation guarantees that the true distance
  //    from E to the target S2ShapeIndex is at least (d - 0.01), and
  //    furthermore no other edge E' of the query S2ShapeIndex is closer to
  //    the target S2ShapeIndex than (d - 0.01).
  //
  //    As always, this option does not affect max_distance().  Continuing the
  //    example above, if max_distance() == M then the "distance" field of the
  //    result edge satisfies (d < M) no matter how max_error() is set.
  //
  //  - max_edges() and max_distance() are set internally on every method call
  //    in order to implement the Target API.
  class ShapeIndexTarget final : public Target {
   public:
    explicit ShapeIndexTarget(const S2ShapeIndex* index);

    // Specifies that distance will be measured to the boundary and interior
    // of polygons in the S2ShapeIndex rather than to polygon boundaries only.
    //
    // DEFAULT: true
    bool include_interiors() const;
    void set_include_interiors(bool include_interiors);

    // Specifies that the distances should be computed by examining every edge
    // in the S2ShapeIndex (for testing and debugging purposes).
    //
    // DEFAULT: false
    bool use_brute_force() const;
    void set_use_brute_force(bool use_brute_force);

    bool set_max_error(const Distance& max_error) override;
    int max_brute_force_edges() const override;
    S2Cap GetCapBound() const override;
    bool UpdateMinDistance(const S2Point& v0, const S2Point& v1,
                           Distance* min_dist) const override;
    bool UpdateMinDistance(const S2Cell& cell,
                           Distance* min_dist) const override;
    std::vector<int> GetContainingShapes(const S2ShapeIndex& query_index,
                                         int max_shapes) const override;

   private:
    const S2ShapeIndex* index_;
    std::unique_ptr<S2ClosestEdgeQuery> query_;
  };

  // Convenience constructor that calls Init().  Options may be specified here
  // or changed at any time using the mutable_options() accessor method.
  explicit S2ClosestEdgeQuery(const S2ShapeIndex* index,
                              const Options& options = Options());

  // Default constructor; requires Init() to be called.
  S2ClosestEdgeQuery();
  ~S2ClosestEdgeQuery();

  // Initializes the query.  Options may be specified here or changed at any
  // time using the mutable_options() accessor method.
  //
  // REQUIRES: ReInit() must be called if "index" is modified.
  // REQUIRES: "index" must persist for the lifetime of this object.
  void Init(const S2ShapeIndex* index, const Options& options = Options());

  // Reinitialize the query.  This method must be called whenever the
  // underlying S2ShapeIndex is modified.
  void ReInit();

  // Return a reference to the underlying S2ShapeIndex.
  const S2ShapeIndex& index() const;

  // Returns the query options.  Options can be modifed between queries.
  const Options& options() const;
  Options* mutable_options();

  // Returns the closest edges to the given target that satisfy the given
  // options.  This method may be called multiple times.
  //
  // Note that if options().include_interiors() is true, the result vector may
  // include some entries with edge_id == -1.  This indicates that the target
  // intersects the indexed polygon with the given shape_id.
  std::vector<Result> FindClosestEdges(Target* target);

  // This version can be more efficient when this method is called many times,
  // since it does not require allocating a new vector on each call.
  void FindClosestEdges(Target* target, std::vector<Result>* results);

  //////////////////////// Convenience Methods ////////////////////////

  // Returns the closest edge to the target.  If no edge satisfies the search
  // criteria, then the Result object will have distance == Infinity() and
  // shape_id == edge_id == -1.
  //
  // Note that if options.include_interiors() is true, edge_id == -1 is also
  // used to indicate that the target intersects an indexed polygon (but in
  // that case distance == Zero() and shape_id >= 0).
  Result FindClosestEdge(Target* target);

  // Returns the minimum distance to the target.  If the index or target is
  // empty, returns S1ChordAngle::Infinity().
  S1ChordAngle GetDistance(Target* target);

  // Returns true if the distance to the target is less than "limit".
  //
  // This method is usually *much* faster than calling GetDistance(), since it
  // is much less work to determine whether the minimum distance is above or
  // below a threshold than it is to calculate the actual minimum distance.
  //
  // You can test whether the distance is less than or equal to "limit" by
  // testing whether IsDistanceLess(target, limit.Successor()).
  bool IsDistanceLess(Target* target, S1ChordAngle limit);

  // Like IsDistanceLess(), except that "limit" is increased by the maximum
  // error in the distance calculation.  This ensures that this function
  // returns true whenever the true, exact distance is less than "limit".
  //
  // For example, suppose that we want to test whether two geometries might
  // intersect each other after they are snapped together using S2Builder
  // (using the IdentitySnapFunction with a given "snap_radius").  Since
  // S2Builder uses exact distance predicates (s2predicates.h), we need to
  // measure the distance between the two geometries conservatively.  If the
  // distance is definitely greater than "snap_radius", then the geometries
  // are guaranteed to not intersect after snapping.
  bool IsConservativeDistanceLess(Target* target, S1ChordAngle limit);

  // Returns the endpoints of the given result edge.
  //
  // CAVEAT: If options().include_interiors() is true, then clients must not
  // pass this method any Result objects that correspond to shape interiors,
  // i.e. those where result.edge_id < 0.
  //
  // REQUIRES: result.edge_id >= 0
  S2Shape::Edge GetEdge(const Result& result) const;

  // Returns the point on given result edge that is closest to "point".
  S2Point Project(const S2Point& point, const Result& result) const;

 private:
  Options options_;
  Base base_;

  S2ClosestEdgeQuery(const S2ClosestEdgeQuery&) = delete;
  void operator=(const S2ClosestEdgeQuery&) = delete;
};


//////////////////   Implementation details follow   ////////////////////


inline bool S2ClosestEdgeQuery::Distance::UpdateMin(const Distance& dist) {
  if (dist < *this) {
    *this = dist;
    return true;
  }
  return false;
}

inline void S2ClosestEdgeQuery::Options::set_max_distance(
    S1ChordAngle max_distance) {
  Base::Options::set_max_distance(Distance(max_distance));
}

inline void S2ClosestEdgeQuery::Options::set_max_distance(
    S1Angle max_distance) {
  Base::Options::set_max_distance(Distance(max_distance));
}

inline void S2ClosestEdgeQuery::Options::set_max_error(S1ChordAngle max_error) {
  Base::Options::set_max_error(Distance(max_error));
}

inline void S2ClosestEdgeQuery::Options::set_max_error(S1Angle max_error) {
  Base::Options::set_max_error(Distance(max_error));
}

inline bool S2ClosestEdgeQuery::ShapeIndexTarget::include_interiors() const {
  return query_->options().include_interiors();
}

inline void S2ClosestEdgeQuery::ShapeIndexTarget::set_include_interiors(
    bool include_interiors) {
  query_->mutable_options()->set_include_interiors(include_interiors);
}

inline bool S2ClosestEdgeQuery::ShapeIndexTarget::use_brute_force() const {
  return query_->options().use_brute_force();
}

inline void S2ClosestEdgeQuery::ShapeIndexTarget::set_use_brute_force(
    bool use_brute_force) {
  query_->mutable_options()->set_use_brute_force(use_brute_force);
}

inline S2ClosestEdgeQuery::S2ClosestEdgeQuery(const S2ShapeIndex* index,
                                              const Options& options) {
  Init(index, options);
}

inline void S2ClosestEdgeQuery::Init(const S2ShapeIndex* index,
                                     const Options& options) {
  options_ = options;
  base_.Init(index);
}

inline void S2ClosestEdgeQuery::ReInit() {
  base_.ReInit();
}

inline const S2ShapeIndex& S2ClosestEdgeQuery::index() const {
  return base_.index();
}

inline const S2ClosestEdgeQuery::Options& S2ClosestEdgeQuery::options() const {
  return options_;
}

inline S2ClosestEdgeQuery::Options* S2ClosestEdgeQuery::mutable_options() {
  return &options_;
}

inline std::vector<S2ClosestEdgeQuery::Result>
S2ClosestEdgeQuery::FindClosestEdges(Target* target) {
  return base_.FindClosestEdges(target, options_);
}

inline void S2ClosestEdgeQuery::FindClosestEdges(Target* target,
                                                 std::vector<Result>* results) {
  base_.FindClosestEdges(target, options_, results);
}

inline S2ClosestEdgeQuery::Result S2ClosestEdgeQuery::FindClosestEdge(
    Target* target) {
  static_assert(sizeof(Options) <= 24, "Consider not copying Options here");
  Options tmp_options = options_;
  tmp_options.set_max_edges(1);
  return base_.FindClosestEdge(target, tmp_options);
}

inline S1ChordAngle S2ClosestEdgeQuery::GetDistance(Target* target) {
  return FindClosestEdge(target).distance;
}

inline S2Shape::Edge S2ClosestEdgeQuery::GetEdge(const Result& result) const {
  return index().shape(result.shape_id)->edge(result.edge_id);
}

inline S2Point S2ClosestEdgeQuery::Project(const S2Point& point,
                                           const Result& result) const {
  if (result.edge_id < 0) return point;
  auto edge = GetEdge(result);
  return S2::Project(point, edge.v0, edge.v1);
}

#endif  // S2_S2CLOSESTEDGEQUERY_H_
