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

#include "s2/s2closestedgequery.h"

#include <memory>
#include "s2/third_party/absl/memory/memory.h"
#include "s2/util/btree/btree_set.h"
#include "s2/s1angle.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cellid.h"
#include "s2/s2cellunion.h"
#include "s2/s2edge_distances.h"
#include "s2/s2regioncoverer.h"
#include "s2/s2shapeindex_region.h"

using std::vector;

S2ClosestEdgeQuery::PointTarget::PointTarget(S2Point const& point)
    : point_(point) {
}

int S2ClosestEdgeQuery::PointTarget::max_brute_force_edges() const {
  // Using BM_FindClosest (which finds the single closest edge), the
  // break-even points are approximately 80, 100, and 250 edges for point
  // cloud, fractal, and regular loop geometry respectively.
  return 120;
}

S2Cap S2ClosestEdgeQuery::PointTarget::GetCapBound() const {
  return S2Cap(point_, S1ChordAngle::Zero());
}

bool S2ClosestEdgeQuery::PointTarget::UpdateMinDistance(
    S2Point const& v0, S2Point const& v1, Distance* min_dist) const {
  return S2::UpdateMinDistance(point_, v0, v1, min_dist);
}

bool S2ClosestEdgeQuery::PointTarget::UpdateMinDistance(
    S2Cell const& cell, Distance* min_dist) const {
  Distance dist(cell.GetDistance(point_));
  if (dist > *min_dist) return false;
  *min_dist = dist;
  return true;
}

vector<int> S2ClosestEdgeQuery::PointTarget::GetContainingShapes(
    S2ShapeIndexBase const& index, int max_shapes) const {
  vector<int> results;
  MakeS2ContainsPointQuery(&index).VisitContainingShapes(
      point_, [&results, max_shapes](S2Shape* shape) {
        results.push_back(shape->id());
        return results.size() < max_shapes;
      });
  return results;
}

S2ClosestEdgeQuery::EdgeTarget::EdgeTarget(S2Point const& a, S2Point const& b)
    : a_(a), b_(b) {
}

int S2ClosestEdgeQuery::EdgeTarget::max_brute_force_edges() const {
  // Using BM_FindClosestToEdge (which finds the single closest edge), the
  // break-even points are approximately 40, 50, and 100 edges for point
  // cloud, fractal, and regular loop geometry respectively.
  return 60;
}

S2Cap S2ClosestEdgeQuery::EdgeTarget::GetCapBound() const {
  return S2Cap((a_ + b_).Normalize(), 0.5 * S1Angle(a_, b_));
}

bool S2ClosestEdgeQuery::EdgeTarget::UpdateMinDistance(
    S2Point const& v0, S2Point const& v1, Distance* min_dist) const {
  return S2::UpdateEdgePairMinDistance(a_, b_, v0, v1, min_dist);
}

bool S2ClosestEdgeQuery::EdgeTarget::UpdateMinDistance(
    S2Cell const& cell, Distance* min_dist) const {
  Distance dist(cell.GetDistance(a_, b_));
  if (dist > *min_dist) return false;
  *min_dist = dist;
  return true;
}

vector<int> S2ClosestEdgeQuery::EdgeTarget::GetContainingShapes(
    S2ShapeIndexBase const& index, int max_shapes) const {
  // We test the center of the edge in order to ensure that edge targets AB
  // and BA yield identical results (which is not guaranteed by the API but
  // users might expect).  Other options would be to test both endpoints, or
  // return different results for AB and BA in some cases.
  PointTarget target((a_ + b_).Normalize());
  return target.GetContainingShapes(index, max_shapes);
}

S2ClosestEdgeQuery::CellTarget::CellTarget(S2Cell const& cell)
    : cell_(cell) {
}

int S2ClosestEdgeQuery::CellTarget::max_brute_force_edges() const {
  // Using BM_FindClosestToCell (which finds the single closest edge), the
  // break-even points are approximately 20, 25, and 40 edges for point cloud,
  // fractal, and regular loop geometry respectively.
  return 30;
}

S2Cap S2ClosestEdgeQuery::CellTarget::GetCapBound() const {
  return cell_.GetCapBound();
}

bool S2ClosestEdgeQuery::CellTarget::UpdateMinDistance(
    S2Point const& v0, S2Point const& v1, Distance* min_dist) const {
  Distance dist(cell_.GetDistance(v0, v1));
  if (dist > *min_dist) return false;
  *min_dist = dist;
  return true;
}

bool S2ClosestEdgeQuery::CellTarget::UpdateMinDistance(
    S2Cell const& cell, Distance* min_dist) const {
  Distance dist(cell_.GetDistance(cell));
  if (dist > *min_dist) return false;
  *min_dist = dist;
  return true;
}

vector<int> S2ClosestEdgeQuery::CellTarget::GetContainingShapes(
    S2ShapeIndexBase const& index, int max_shapes) const {
  // We are required to return all polygons that contain the cell, and
  // optionally we can return any polygon that intersects the cell (to within
  // a small error margin).  There are many possible algorithms to do this;
  // the simplest is just to return the polygons that contain the cell center.
  // However the following is less work and returns more shapes in the case
  // where the target cell is large.
  S2ShapeIndexBase::Iterator it(&index);
  S2Point center = cell_.GetCenter();
  if (!it.Locate(center)) return vector<int>();

  if (it.id().level() < cell_.level()) {
    // The target cell is properly contained by the index cell.  The simplest
    // approach is just to use the cell center.
    return PointTarget(center).GetContainingShapes(index, max_shapes);
  }
  // Otherwise, the target cell contains the index cell.  In that case all
  // shapes of dimension 2 present in that cell may be returned.
  vector<int> results;
  S2ShapeIndexCell const& index_cell = it.cell();
  for (int s = 0; s < index_cell.num_clipped(); ++s) {
    S2ClippedShape const& clipped = index_cell.clipped(s);
    if (index.shape(clipped.shape_id())->has_interior()) {
      results.push_back(clipped.shape_id());
      if (results.size() >= max_shapes) break;
    }
  }
  return results;
}

S2ClosestEdgeQuery::ShapeIndexTarget::ShapeIndexTarget(
    S2ShapeIndex const* index)
    : index_(index), query_(absl::make_unique<S2ClosestEdgeQuery>(index)) {
}

bool S2ClosestEdgeQuery::ShapeIndexTarget::set_max_error(
    Distance const& max_error) {
  query_->mutable_options()->set_max_error(max_error);
  return true;  // Indicates that we may return suboptimal results.
}

int S2ClosestEdgeQuery::ShapeIndexTarget::max_brute_force_edges() const {
  // For BM_FindClosestToSameSizeAbuttingIndex (which uses two nearby indexes
  // with similar edge counts), the break-even points are approximately 20,
  // 30, and 40 edges for point cloud, fractal, and regular loop geometry
  // respectively.
  return 25;
}

S2Cap S2ClosestEdgeQuery::ShapeIndexTarget::GetCapBound() const {
  return MakeS2ShapeIndexRegion(index_).GetCapBound();
}

bool S2ClosestEdgeQuery::ShapeIndexTarget::UpdateMinDistance(
    S2Point const& v0, S2Point const& v1, Distance* min_dist) const {
  query_->mutable_options()->set_max_distance(*min_dist);
  EdgeTarget target(v0, v1);
  Result r = query_->FindClosestEdge(&target);
  if (r.shape_id < 0) return false;
  *min_dist = r.distance;
  return true;
}

bool S2ClosestEdgeQuery::ShapeIndexTarget::UpdateMinDistance(
    S2Cell const& cell, Distance* min_dist) const {
  query_->mutable_options()->set_max_distance(*min_dist);
  CellTarget target(cell);
  Result r = query_->FindClosestEdge(&target);
  if (r.shape_id < 0) return false;
  *min_dist = r.distance;
  return true;
}

vector<int> S2ClosestEdgeQuery::ShapeIndexTarget::GetContainingShapes(
    S2ShapeIndexBase const& query_index, int max_shapes) const {
  // It is sufficient to find the set of chain starts in the target index
  // (i.e., one vertex per connected component of edges) that are contained by
  // the query index.
  //
  // TODO(ericv): Do this by merge-joining the two S2ShapeIndexes, and share
  // the code with S2BoundaryOperation.
  util::btree::btree_set<int> results;
  auto query = MakeS2ContainsPointQuery(&query_index);
  int num_shape_ids = index_->num_shape_ids();
  for (int s = 0; s < num_shape_ids && max_shapes > 0; ++s) {
    S2Shape* shape = index_->shape(s);
    if (shape == nullptr) continue;
    int num_chains = shape->num_chains();
    for (int chain = 0; chain < num_chains && max_shapes > 0; ++chain) {
      query.VisitContainingShapes(
          shape->chain_edge(chain, 0).v0,
          [&results, &max_shapes](S2Shape* shape) {
            if (results.insert(shape->id()).second) --max_shapes;
            return max_shapes > 0;
          });
    }
  }
  return vector<int>(results.begin(), results.end());
}

void S2ClosestEdgeQuery::Options::set_conservative_max_distance(
    S1ChordAngle max_distance) {
  set_max_distance(Distance(max_distance.PlusError(
      S2::GetUpdateMinDistanceMaxError(max_distance))));
}

S2ClosestEdgeQuery::S2ClosestEdgeQuery() {
  // Prevent inline constructor bloat by defining here.
}

S2ClosestEdgeQuery::~S2ClosestEdgeQuery() {
  // Prevent inline destructor bloat by defining here.
}

bool S2ClosestEdgeQuery::IsDistanceLess(Target* target, S1ChordAngle limit) {
  options_.set_max_distance(limit);
  options_.set_max_error(limit);
  return FindClosestEdge(target).shape_id >= 0;
}
