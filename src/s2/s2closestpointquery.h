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

#ifndef S2_S2CLOSESTPOINTQUERY_H_
#define S2_S2CLOSESTPOINTQUERY_H_

#include <vector>

#include <glog/logging.h>
#include "s2/third_party/absl/container/inlined_vector.h"
#include "s2/priority_queue_sequence.h"
#include "s2/s1angle.h"
#include "s2/s1chordangle.h"
#include "s2/s2cap.h"
#include "s2/s2cellid.h"
#include "s2/s2cellunion.h"
#include "s2/s2edgeutil.h"
#include "s2/s2pointindex.h"
#include "s2/s2regioncoverer.h"

// Given a set of points stored in an S2PointIndex, S2ClosestPointQuery
// provides methods that find the closest point(s) to a given query point
// or query edge.  Example usage:
//
// void Test(vector<S2Point> const& points, vector<S2Point> const& targets) {
//   // The template argument allows auxiliary data to be attached to each
//   // point (in this case, the array index).
//   S2PointIndex<int> index;
//   for (S2Point const& point : points) {
//     index.Add(point, i);
//   }
//   S2ClosestPointQuery<int> query(index);
//   query.set_max_points(15);
//   for (S2Point const& target : targets) {
//     query.FindClosestPoints(target);
//     for (int j = 0; j < query.num_points(); ++j) {
//       // query.data(j) is the auxiliary data (the "points" array index).
//       // query.distance(j) is the distance to the target point.
//       DoSomething(target, query.point(j), query.data(j), query.distance(j));
//     }
//   }
// }
//
// You can find either the k closest points, or all points within a given
// radius, or both (i.e., the k closest points up to a given maximum radius).
// E.g. to find all the points within 5 kilometers, call
//
//   query.set_max_distance(S2Earth::ToAngle(util::units::Kilometers(5)));
//
// You can also restrict the results to an arbitrary S2Region, e.g.
//
//   S2LatLngRect rect(...);
//   query.set_region(&rect);  // Does *not* take ownership.
//
// To find the closest points to a query edge rather than a point, use:
//
//   query.FindClosestPointsToEdge(v0, v1);
//
// The implementation is designed to be very fast for both small and large
// point sets.

template <class Data>
class S2ClosestPointQuery {
 public:
  using Index = S2PointIndex<Data>;

  // Convenience constructor that calls Init().
  explicit S2ClosestPointQuery(Index const& index);

  // Default constructor; requires Init() to be called.
  S2ClosestPointQuery() {}
  ~S2ClosestPointQuery();

  // Initialize the query.
  // REQUIRES: Reset() must be called if "index" is modified.
  void Init(Index const& index);

  // Reset the query state.  This method must be called whenever the
  // underlying index is modified.
  void Reset();

  // Return a reference to the underlying S2PointIndex.
  Index const& index() const;

  // Only find the "max_points" closest points.
  // This value may be changed between calls to FindClosestPoints().
  //
  // Default value: numeric_limits<int>::max()
  // REQUIRES: max_points >= 1
  int max_points() const;
  void set_max_points(int max_points);

  // Only find points whose distance to the target is less than "max_distance".
  // This value may be changed between calls to FindClosestPoints().
  //
  // Default value: S1Angle::Infinity().
  S1Angle max_distance() const;
  void set_max_distance(S1Angle max_distance);

  // Only find points in the given S2Region.  "region" is owned by the caller
  // and must persist during any subsequent call(s) to FindClosestPoints.
  // This value may be changed between calls to FindClosestPoints(), or reset
  // by calling set_region(nullptr).
  //
  // Note that if you want to set the region to a disc around the target
  // point, it is faster to use set_max_distance() instead.  You can also call
  // both methods, e.g. if you want to limit the maximum distance to the
  // target and also require that points lie within a given rectangle.
  S2Region const* region() const;
  void set_region(S2Region const* region);

  // Find the closest points to "target" that satisfy the given max_distance()
  // and/or max_points() and/or region() criteria.  If none of these is set,
  // then all points are returned.
  void FindClosestPoints(S2Point const& target);

  // Convenience method that finds the closest point to the given target
  // point.  The result consists of either zero or one points, accessible
  // using the iterator interface below.  This method is equivalent to:
  //
  //   query.set_max_points(1);
  //   query.FindClosestPoints();
  //
  // Note that set_max_distance() will still limit the search radius.
  void FindClosestPoint(S2Point const& target);

  // Find the closest points to the given edge AB.  Otherwise similar to
  // FindClosestPoints().
  void FindClosestPointsToEdge(S2Point const& a, S2Point const& b);

  // Manually specify whether distances are computed using "brute force"
  // (i.e., by examining every point) rather than using the S2PointIndex.
  // This is useful for testing, benchmarking, and debugging.
  //
  // REQUIRES: Init() has been called.
  void UseBruteForce(bool use_brute_force);

  ////////////////////// Result Interface ////////////////////////
  //
  // The result of a query consists of a set of points which may be obtained
  // using the interface below.  Points are sorted in order of increasing
  // distance.

  // The number of result points.
  int num_points() const;

  // The result point at the given index.
  S2Point const& point(int i) const;

  // The distance to the result point at the given index.
  S1Angle distance(int i) const;

  // Like distance(i), but expressed as an S1ChordAngle.
  S1ChordAngle distance_ca(int i) const;

  // The client data associated with the result point at the given index.
  Data const& data(int i) const;

 private:
  using PointData = typename Index::PointData;
  using Iterator = typename Index::Iterator;

  class PointTarget {
   public:
    explicit PointTarget(S2Point const& point) : point_(point) {}
    S2Point center() const { return point_; }
    S1Angle radius() const { return S1Angle::Zero(); }
    bool UpdateMinDistance(S2Point const& x, S1ChordAngle* min_dist) const {
      S1ChordAngle distance(x, point_);
      // Only return true if the new distance is smaller.
      if (distance >= *min_dist) return false;
      *min_dist = distance;
      return true;
    }
    S1ChordAngle GetDistance(S2Cell const& cell) const {
      return cell.GetDistance(point_);
    }
   private:
    S2Point point_;
  };

  class EdgeTarget {
   public:
    EdgeTarget(S2Point const& a, S2Point const& b) : a_(a), b_(b) {}
    S2Point center() const { return (a_ + b_).Normalize(); }
    S1Angle radius() const { return 0.5 * S1Angle(a_, b_); }
    bool UpdateMinDistance(S2Point const& x, S1ChordAngle* min_dist) const {
      return S2EdgeUtil::UpdateMinDistance(x, a_, b_, min_dist);
    }
    S1ChordAngle GetDistance(S2Cell const& cell) const {
      return cell.GetDistanceToEdge(a_, b_);
    }
   private:
    S2Point a_, b_;
  };

  void InitIndexCovering();
  void CoverRange(Iterator const& first, Iterator const& last);

  // Using templates rather than virtual functions speeds up the benchmarks
  // by 15-25% when the brute force algorithm is used (< 150 points).

  template <class Target>
  void FindClosestPointsToTarget(Target const& target);

  template <class Target>
  void FindClosestPointsBruteForce(Target const& target);

  template <class Target>
  void FindClosestPointsOptimized(Target const& target);

  template <class Target>
  void MaybeAddResult(PointData const& point_data, Target const& target);

  template <class Target>
  void InitQueue(Target const& target);

  template <class Target>
  bool AddCell(S2CellId id, Iterator* iter, bool seek, Target const& target);

  //////////// Constants (tuned using the benchmarks) ////////////

  // The maximum number of points to process by brute force.
  static int const kMaxBruteForcePoints = 150;

  // The maximum number of points to process without subdividing further.
  static int const kMaxLeafPoints = 12;

  //////////// Parameters ////////////

  int max_points_;
  S1Angle max_distance_;
  S2Region const* region_;

  //////////// Fields that are constant during a query ////////////

  Index const* index_;

  // If the index has few edges, it is cheaper to use a brute force algorithm.
  bool use_brute_force_;

  // A small precomputed S2CellId covering of the indexed points.
  std::vector<S2CellId> index_covering_;

  //////////// Fields that are updated during a query ////////////

  // The result points gathered so far are kept in a priority queue.  Once
  // all points have been found, the underlying vector is sorted by distance
  // so that we can iterate through them sequentially.
  struct Result {
    S1ChordAngle distance;
    PointData const* point_data;
    Result(S1ChordAngle _distance, PointData const* _point_data)
        : distance(_distance), point_data(_point_data) {
    }
    bool operator<(Result const& other) const {
      // The algorithm works by replacing the result whose distance is largest
      // when a better candidate is found, so we keep the entries sorted such
      // that the largest distance is at the top of the heap.
      return distance < other.distance;
    }
  };
  using ResultHeap = priority_queue_sequence<Result>;
  ResultHeap results_;

  // The distance beyond which we can safely ignore further candidate points.
  // (Candidates that are exactly at the limit are ignored; this makes things
  // easier in the case of S2ClosestEdgeQuery and should not affect clients
  // since distance measurements have a small amount of error anyway.)
  //
  // Initially this is the same as the maximum distance specified by the user,
  // but it can also be updated by the algorithm (see MaybeAddResult).
  S1ChordAngle max_distance_limit_;

  // We also keep a priority queue of unprocessed S2Cells.
  struct QueueEntry {
    S1ChordAngle distance;  // Distance from target to any point in the cell.
    S2CellId id;
    QueueEntry(S1ChordAngle _distance, S2CellId _id)
      : distance(_distance), id(_id) {
    }
    bool operator<(QueueEntry const& other) const {
      // Sort the queue entries so that smaller distances are returned first.
      return distance > other.distance;
    }
  };
  using CellQueue =
      std::priority_queue<QueueEntry, gtl::InlinedVector<QueueEntry, 16>>;
  CellQueue queue_;

  // Temporaries, defined here to avoid multiple allocations / initializations.
  Iterator iter_;
  std::vector<S2CellId> region_covering_;
  std::vector<S2CellId> max_distance_covering_;
  std::vector<S2CellId> intersection_with_region_;
  std::vector<S2CellId> intersection_with_max_distance_;

  PointData const* tmp_point_data_[kMaxLeafPoints];

  S2ClosestPointQuery(S2ClosestPointQuery const&) = delete;
  void operator=(S2ClosestPointQuery const&) = delete;
};


//////////////////   Implementation details follow   ////////////////////


template <class Data>
S2ClosestPointQuery<Data>::~S2ClosestPointQuery() {
}

template <class Data>
S2ClosestPointQuery<Data>::S2ClosestPointQuery(Index const& index) {
  Init(index);
}

template <class Data>
void S2ClosestPointQuery<Data>::Init(Index const& index) {
  index_ = &index;
  max_points_ = std::numeric_limits<int>::max();
  max_distance_ = S1Angle::Infinity();
  region_ = nullptr;
  Reset();
}

template <class Data>
void S2ClosestPointQuery<Data>::Reset() {
  results_.mutable_rep()->clear();
  iter_.Init(*index_);
  UseBruteForce(index_->num_points() <= kMaxBruteForcePoints);
}

template <class Data>
S2PointIndex<Data> const& S2ClosestPointQuery<Data>::index() const {
  return *index_;
}

template <class Data>
int S2ClosestPointQuery<Data>::max_points() const {
  return max_points_;
}

template <class Data>
void S2ClosestPointQuery<Data>::set_max_points(int max_points) {
  DCHECK_GE(max_points, 1);
  max_points_ = max_points;
}

template <class Data>
S1Angle S2ClosestPointQuery<Data>:: max_distance() const {
  return max_distance_;
}

template <class Data>
void S2ClosestPointQuery<Data>::set_max_distance(S1Angle max_distance) {
  max_distance_ = max_distance;
}

template <class Data>
S2Region const* S2ClosestPointQuery<Data>::region() const {
  return region_;
}

template <class Data>
void S2ClosestPointQuery<Data>::set_region(S2Region const* region) {
  region_ = region;
}

template <class Data>
int S2ClosestPointQuery<Data>::num_points() const {
  return results_.size();
}

template <class Data>
S2Point const& S2ClosestPointQuery<Data>::point(int i) const {
  return results_.rep()[i].point_data->point();
}

template <class Data>
S1ChordAngle S2ClosestPointQuery<Data>::distance_ca(int i) const {
  return results_.rep()[i].distance;
}

template <class Data>
S1Angle S2ClosestPointQuery<Data>::distance(int i) const {
  return distance_ca(i).ToAngle();
}

template <class Data>
Data const& S2ClosestPointQuery<Data>::data(int i) const {
  return results_.rep()[i].point_data->data();
}

template <class Data>
void S2ClosestPointQuery<Data>::UseBruteForce(bool use_brute_force) {
  use_brute_force_ = use_brute_force;
  if (!use_brute_force) InitIndexCovering();
}

template <class Data>
void S2ClosestPointQuery<Data>::InitIndexCovering() {
  // Compute the "index covering", which is a small number of S2CellIds that
  // cover the indexed points.  There are two cases:
  //
  //  - If the index spans more than one face, then there is one covering cell
  // per spanned face, just big enough to cover the index cells on that face.
  //
  //  - If the index spans only one face, then we find the smallest cell "C"
  // that covers the index cells on that face (just like the case above).
  // Then for each of the 4 children of "C", if the child contains any index
  // cells then we create a covering cell that is big enough to just fit
  // those index cells (i.e., shrinking the child as much as possible to fit
  // its contents).  This essentially replicates what would happen if we
  // started with "C" as the covering cell, since "C" would immediately be
  // split, except that we take the time to prune the children further since
  // this will save work on every subsequent query.
  index_covering_.clear();
  iter_.Reset();
  if (iter_.Done()) return;  // Empty index.

  Iterator next = iter_, last = iter_;
  last.Finish();
  last.Prev();
  if (next.id() != last.id()) {
    // The index has at least two cells.  Choose a level such that the entire
    // index can be spanned with at most 6 cells (if the index spans multiple
    // faces) or 4 cells (it the index spans a single face).
    int level = next.id().GetCommonAncestorLevel(last.id()) + 1;

    // Visit each potential covering cell except the last (handled below).
    S2CellId last_id = last.id().parent(level);
    for (S2CellId id = next.id().parent(level); id != last_id; id = id.next()) {
      // Skip any covering cells that don't contain any index cells.
      if (id.range_max() < next.id()) continue;

      // Find the range of index cells contained by this covering cell and
      // then shrink the cell if necessary so that it just covers them.
      Iterator cell_first = next;
      next.Seek(id.range_max().next());
      Iterator cell_last = next;
      cell_last.Prev();
      CoverRange(cell_first, cell_last);
    }
  }
  CoverRange(next, last);
}

// Adds a cell to index_covering_ that covers the given inclusive range.
// REQUIRES: "first" and "last" have a common ancestor.
template <class Data>
void S2ClosestPointQuery<Data>::CoverRange(Iterator const& first,
                                           Iterator const& last) {
  // Add the lowest common ancestor of the given range.
  int level = first.id().GetCommonAncestorLevel(last.id());
  DCHECK_GE(level, 0);
  index_covering_.push_back(first.id().parent(level));
}

template <class Data>
void S2ClosestPointQuery<Data>::FindClosestPoint(S2Point const& target) {
  set_max_points(1);
  FindClosestPoints(target);
}

template <class Data>
void S2ClosestPointQuery<Data>::FindClosestPoints(S2Point const& point) {
  PointTarget target(point);
  FindClosestPointsToTarget(target);
}

template <class Data>
void S2ClosestPointQuery<Data>::FindClosestPointsToEdge(
    S2Point const& a, S2Point const& b) {
  EdgeTarget target(a, b);
  FindClosestPointsToTarget(target);
}

template <class Data> template <class Target>
void S2ClosestPointQuery<Data>::FindClosestPointsToTarget(
    Target const& target) {
  max_distance_limit_ = S1ChordAngle(max_distance_);
  max_distance_limit_ = max_distance_limit_.PlusError(
      max_distance_limit_.GetS1AngleConstructorMaxError() +
      max_distance_limit_.GetS2PointConstructorMaxError());
  results_.mutable_rep()->clear();
  if (use_brute_force_) {
    FindClosestPointsBruteForce(target);
  } else {
    FindClosestPointsOptimized(target);
  }
  std::sort(results_.mutable_rep()->begin(), results_.mutable_rep()->end());
}

template <class Data> template <class Target>
void S2ClosestPointQuery<Data>::FindClosestPointsBruteForce(
    Target const& target) {
  for (iter_.Reset(); !iter_.Done(); iter_.Next()) {
    MaybeAddResult(iter_.point_data(), target);
  }
}

template <class Data> template <class Target>
void S2ClosestPointQuery<Data>::FindClosestPointsOptimized(
    Target const& target) {
  InitQueue(target);
  while (!queue_.empty()) {
    // We need to copy the top entry before removing it, and we need to remove
    // it before adding any new entries to the queue.
    QueueEntry entry = queue_.top();
    queue_.pop();
    if (entry.distance >= max_distance_limit_) {
      queue_ = CellQueue();  // Clear any remaining entries.
      break;
    }
    S2CellId child = entry.id.child_begin();
    // We already know that it has too many points, so process its children.
    // Each child may either be processed directly or enqueued again.  The
    // loop is optimized so that we don't seek unnecessarily.
    bool seek = true;
    for (int i = 0; i < 4; ++i, child = child.next()) {
      seek = AddCell(child, &iter_, seek, target);
    }
  }
}

template <class Data> template <class Target>
void S2ClosestPointQuery<Data>::MaybeAddResult(
    PointData const& point_data, Target const& target) {
  S1ChordAngle distance = max_distance_limit_;
  if (!target.UpdateMinDistance(point_data.point(), &distance)) return;
  if (region_ && !region_->VirtualContainsPoint(point_data.point())) return;

  // Add this point to results_.
  if (results_.size() >= max_points_) {
    results_.pop();  // Replace the furthest result point.
  }
  results_.push(Result(distance, &point_data));
  if (results_.size() >= max_points_) {
    max_distance_limit_ = results_.top().distance;
  }
}

template <class Data> template <class Target>
void S2ClosestPointQuery<Data>::InitQueue(Target const& target) {
  DCHECK(queue_.empty());

  // Optimization: rather than starting with the entire index, see if we can
  // limit the search region to a small disc.  Then we can find a covering for
  // that disc and intersect it with the covering for the index.  This can
  // save a lot of work when the search region is small.

  if (max_points_ == 1) {
    // If the user is searching for just the closest point, we can compute an
    // upper bound on search radius by seeking to the target point in the
    // index and looking at the adjacent index points (in S2CellId order).
    // The minimum distance to either of these points is an upper bound on the
    // search radius.
    //
    // TODO(ericv): The same strategy would also work for small values of
    // max_points() > 1, e.g. max_points() == 20, except that we would need to
    // examine more neighbors (at least 20, and preferably 20 in each
    // direction).  It's not clear whether this is a common case, though, and
    // also this would require extending MaybeAddResult() so that it can
    // remove duplicate entries.  (The points added here may be re-added by
    // AddCell(), but this is okay when max_points() == 1.)
    iter_.Seek(S2CellId(target.center()));
    if (!iter_.Done()) {
      MaybeAddResult(iter_.point_data(), target);
    }
    if (!iter_.AtBegin()) {
      iter_.Prev();
      MaybeAddResult(iter_.point_data(), target);
    }
  }
  // We start with a covering of the set of indexed points, then intersect it
  // with the given region (if any) and maximum search radius disc (if any).
  std::vector<S2CellId> const* initial_cells = &index_covering_;
  if (region_) {
    S2RegionCoverer coverer;
    coverer.set_max_cells(4);
    coverer.GetCovering(*region_, &region_covering_);
    S2CellUnion::GetIntersection(index_covering_, region_covering_,
                                 &intersection_with_region_);
    initial_cells = &intersection_with_region_;
  }
  if (max_distance_limit_ < S1ChordAngle::Infinity()) {
    S2RegionCoverer coverer;
    coverer.set_max_cells(4);
    S2Cap search_cap(target.center(),
                     target.radius() + max_distance_limit_.ToAngle());
    coverer.GetFastCovering(search_cap, &max_distance_covering_);
    S2CellUnion::GetIntersection(*initial_cells, max_distance_covering_,
                                 &intersection_with_max_distance_);
    initial_cells = &intersection_with_max_distance_;
  }
  iter_.Reset();
  for (int i = 0; i < initial_cells->size() && !iter_.Done(); ++i) {
    S2CellId id = (*initial_cells)[i];
    AddCell(id, &iter_, id.range_min() > iter_.id() /*seek*/, target);
  }
}

// Either process the contents of the given cell immediately, or add it to the
// queue to be subdivided.  If "seek" is false, then "iter" must already be
// positioned at the first indexed point within this cell.
//
// Returns "true" if the cell was added to the queue, and "false" if it was
// processed immediately, in which case "iter" is left positioned at the next
// cell in S2CellId order.
template <class Data> template <class Target>
bool S2ClosestPointQuery<Data>::AddCell(S2CellId id, Iterator* iter,
                                        bool seek, Target const& target) {
  if (seek) iter->Seek(id.range_min());
  if (id.is_leaf()) {
    // Leaf cells can't be subdivided.
    for (; !iter->Done() && iter->id() == id; iter->Next()) {
      MaybeAddResult(iter->point_data(), target);
    }
    return false;  // No need to seek to next child.
  }
  S2CellId last = id.range_max();
  int num_points = 0;
  for (; !iter->Done() && iter->id() <= last; iter->Next()) {
    if (num_points == kMaxLeafPoints) {
      // This child cell has too many points, so enqueue it.
      S2Cell cell(id);
      S1ChordAngle distance = target.GetDistance(cell);
      if (distance < max_distance_limit_) {
        // We delay checking "region_" as long as possible because it may be
        // relatively expensive.
        if (region_ == nullptr || region_->MayIntersect(cell)) {
          queue_.push(QueueEntry(distance, id));
        }
      }
      return true;  // Seek to next child.
    }
    tmp_point_data_[num_points++] = &iter->point_data();
  }
  // There were few enough points that we might as well process them now.
  for (int i = 0; i < num_points; ++i) {
    MaybeAddResult(*tmp_point_data_[i], target);
  }
  return false;  // No need to seek to next child.
}

#endif  // S2_S2CLOSESTPOINTQUERY_H_
