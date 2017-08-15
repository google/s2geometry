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
#include "s2/s1angle.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cellid.h"
#include "s2/s2cellunion.h"
#include "s2/s2edge_distances.h"
#include "s2/s2regioncoverer.h"
#include "s2/s2shapeindex_region.h"

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

S2ClosestEdgeQuery::ShapeIndexTarget::ShapeIndexTarget(
    S2ShapeIndex const* index)
    : index_(index), query_(absl::MakeUnique<S2ClosestEdgeQuery>(index)) {
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

void S2ClosestEdgeQuery::FindClosestEdges(S2Point const& point) {
  target_.reset(new PointTarget(point));
  FindClosestEdges(target_.get(), &results_);
}

void S2ClosestEdgeQuery::FindClosestEdgesToEdge(S2Point const& a,
                                                S2Point const& b) {
  target_.reset(new EdgeTarget(a, b));
  FindClosestEdges(target_.get(), &results_);
}

void S2ClosestEdgeQuery::FindClosestEdge(S2Point const& point) {
  set_max_edges(1);
  FindClosestEdges(point);
}

S1Angle S2ClosestEdgeQuery::GetDistance(S2Point const& point) {
  PointTarget target(point);
  return GetDistance(&target).ToAngle();
}
