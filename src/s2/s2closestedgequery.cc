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

#include "s2/s1angle.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cellid.h"
#include "s2/s2cellunion.h"
#include "s2/s2edge_distances.h"
#include "s2/s2edgeutil.h"
#include "s2/s2regioncoverer.h"

S2ClosestEdgeQuery::PointTarget::PointTarget(S2Point const& point)
    : point_(point) {
}

int S2ClosestEdgeQuery::PointTarget::max_brute_force_edges() const {
  return 180;  // TODO(ericv): Re-tune using benchmarks.
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

S2Point S2ClosestEdgeQuery::PointTarget::GetClosestPointOnEdge(
    S2Shape::Edge const& edge) const {
  return S2::Project(point_, edge.v0, edge.v1);
}

S2ClosestEdgeQuery::EdgeTarget::EdgeTarget(S2Point const& a, S2Point const& b)
    : a_(a), b_(b) {
}

int S2ClosestEdgeQuery::EdgeTarget::max_brute_force_edges() const {
  return 180;  // XXX Re-tune using benchmarks.
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
  Distance dist(cell.GetDistanceToEdge(a_, b_));
  if (dist > *min_dist) return false;
  *min_dist = dist;
  return true;
}

S2Point S2ClosestEdgeQuery::EdgeTarget::GetClosestPointOnEdge(
    S2Shape::Edge const& edge) const {
  return S2::GetEdgePairClosestPoints(a_, b_, edge.v0, edge.v1).second;
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

S2Point S2ClosestEdgeQuery::GetClosestPointOnEdge(int i) const {
  return target_->GetClosestPointOnEdge(edge(i));
}

void S2ClosestEdgeQuery::FindClosestEdges(S2Point const& point) {
  target_.reset(new PointTarget(point));
  FindClosestEdges(*target_, &results_);
}

void S2ClosestEdgeQuery::FindClosestEdgesToEdge(S2Point const& a,
                                                S2Point const& b) {
  target_.reset(new EdgeTarget(a, b));
  FindClosestEdges(*target_, &results_);
}

void S2ClosestEdgeQuery::FindClosestEdge(S2Point const& point) {
  set_max_edges(1);
  FindClosestEdges(point);
}

S1Angle S2ClosestEdgeQuery::GetDistance(S2Point const& point) {
  return GetDistance(PointTarget(point)).ToAngle();
}

S2Point S2ClosestEdgeQuery::Project(S2Point const& point) {
  FindClosestEdge(point);
  if (num_edges() == 0) return point;
  return GetClosestPointOnEdge(0);
}
