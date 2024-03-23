// Copyright Google Inc. All Rights Reserved.
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

#include "s2/s2chain_interpolation_query.h"

#include <algorithm>
#include <vector>

#include "absl/log/absl_check.h"
#include "s2/s1angle.h"
#include "s2/s2edge_distances.h"
#include "s2/s2point.h"
#include "s2/s2shape.h"

S2ChainInterpolationQuery::S2ChainInterpolationQuery(const S2Shape* shape,
                                                     int chain_id) {
  Init(shape, chain_id);
}

void S2ChainInterpolationQuery::Init(const S2Shape* shape, int chain_id) {
  ABSL_CHECK(chain_id < shape->num_chains());

  int first_edge_id, last_edge_id;
  if (chain_id >= 0) {
    // If a valid chain id was provided, then the range of edge ids is defined
    // by the start and the length of the chain.
    const S2Shape::Chain chain = shape->chain(chain_id);
    first_edge_id = chain.start;
    last_edge_id = first_edge_id + chain.length - 1;
  } else {
    // If no valid chain id was provided then we use the whole range of shape's
    // edge ids.
    first_edge_id = 0;
    last_edge_id = shape->num_edges() - 1;
  }

  cumulative_values_.clear();
  cumulative_values_.reserve(last_edge_id + 2 - first_edge_id);

  S1Angle cumulative_angle = S1Angle::Zero();

  for (int i = first_edge_id; i <= last_edge_id; ++i) {
    cumulative_values_.push_back(cumulative_angle);
    const S2Shape::Edge edge = shape->edge(i);
    cumulative_angle += S1Angle(edge.v0, edge.v1);
  }
  if (!cumulative_values_.empty()) {
    cumulative_values_.push_back(cumulative_angle);
  }

  shape_ = shape;
  first_edge_id_ = first_edge_id;
  last_edge_id_ = last_edge_id;
}

S1Angle S2ChainInterpolationQuery::GetLength() const {
  // The total length equals the cumulative value at the end of the last
  // edge, iff there is at least one edge in the shape.
  return cumulative_values_.empty() ? S1Angle::Zero()
                                    : cumulative_values_.back();
}

S1Angle S2ChainInterpolationQuery::GetLengthAtEdgeEnd(const int edge_id) const {
  if (cumulative_values_.empty()) {
    return S1Angle::Zero();
  }

  if (edge_id < first_edge_id_ || edge_id > last_edge_id_) {
    return S1Angle::Infinity();
  }

  return cumulative_values_[edge_id - first_edge_id_ + 1];
}

S2ChainInterpolationQuery::Result S2ChainInterpolationQuery::AtDistance(
    const S1Angle& distance) const {
  // Return an invalid result for uninitialized queries or for shapes containing
  // no edges.
  if (cumulative_values_.empty()) {
    return Result();
  }

  // Binary search in the list of cumulative values, which by construction are
  // sorted in ascending order.
  const auto it = std::lower_bound(cumulative_values_.begin(),
                                   cumulative_values_.end(), distance);

  if (it == cumulative_values_.begin()) {
    // Corner case: the first vertex of the shape at distance = 0.
    return Result(shape_->edge(first_edge_id_).v0, first_edge_id_,
                  cumulative_values_.front());
  } else if (it == cumulative_values_.end()) {
    // Corner case: the input distance is greater than the total length, hence
    // we snap the result to the last vertex of the shape at distance = total
    // length.
    return Result(shape_->edge(last_edge_id_).v1, last_edge_id_,
                  cumulative_values_.back());
  } else {
    // Obtain the edge index and compute the interpolated result from the edge
    // vertices.
    int edge_id = it - cumulative_values_.begin() + first_edge_id_ - 1;
    const S2Shape::Edge edge = shape_->edge(edge_id);
    return Result(S2::GetPointOnLine(edge.v0, edge.v1, distance - *(it - 1)),
                  edge_id, distance);
  }
}

S2ChainInterpolationQuery::Result S2ChainInterpolationQuery::AtFraction(
    double fraction) const {
  return AtDistance(fraction * GetLength());
}

std::vector<S2Point> S2ChainInterpolationQuery::Slice(
    double begin_fraction, double end_fraction) const {
  std::vector<S2Point> slice;
  AddSlice(begin_fraction, end_fraction, slice);
  return slice;
}

void S2ChainInterpolationQuery::AddSlice(double begin_fraction,
                                         double end_fraction,
                                         std::vector<S2Point>& slice) const {
  if (cumulative_values_.empty()) {
    return;
  }

  const int original_size = slice.size();
  const bool reverse = begin_fraction > end_fraction;
  if (reverse) {
    std::swap(begin_fraction, end_fraction);
  }

  Result result = AtFraction(begin_fraction);
  const int begin_edge = result.edge_id();
  S2Point last_point = result.point();
  slice.push_back(last_point);

  result = AtFraction(end_fraction);
  const int end_edge = result.edge_id();
  for (int edge_id = begin_edge; edge_id < end_edge; ++edge_id) {
    const auto edge = shape_->edge(edge_id);
    if (last_point != edge.v1) {
      last_point = edge.v1;
      slice.push_back(last_point);
    }
  }
  slice.push_back(result.point());

  if (reverse) {
    std::reverse(slice.begin() + original_size, slice.end());
  }
}
