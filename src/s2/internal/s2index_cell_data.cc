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

#include "s2/internal/s2index_cell_data.h"

#include <algorithm>
#include <atomic>
#include <utility>

#include "absl/log/absl_check.h"
#include "s2/mutable_s2shape_index.h"
#include "s2/s2cell.h"
#include "s2/s2padded_cell.h"
#include "s2/s2shape.h"

namespace internal {

void S2IndexCellData::LoadCell(const S2ShapeIndex* index, S2CellId id,
                               const S2ShapeIndexCell* cell) {
  ABSL_DCHECK_NE(index, nullptr);

  if (index_ == index && cell_id_ == id) {
    return;
  }

  index_ = index;

  // Cache cell information.
  cell_ = cell;
  cell_id_ = id;

  // Reset atomic flags so we'll recompute cached values.  These form a write
  // barrier with the write to cell_id_ above and so should stay below it.
  s2cell_set_.store(false, std::memory_order_release);
  center_set_.store(false, std::memory_order_release);

  // Clear previous edges.  Reserve one slot in case we don't decode _any_ edges
  // (e.g. an interior cell or all the dimensions are unwanted, so we don't
  // try to create spans with a nullptr).
  edges_.clear();
  shape_regions_.clear();

  // Reset per-dimension region information.
  for (auto& region : dim_regions_) {
    region = Region();
  }

  int min_dim = 0;
  while (min_dim <= 2 && !dim_wanted(min_dim)) {
    ++min_dim;
  }

  int max_dim = 2;
  while (max_dim >= 0 && !dim_wanted(max_dim)) {
    --max_dim;
  }

  // No dimensions wanted, we're done.
  if (min_dim > 2 || max_dim < 0) {
    return;
  }

  for (int dim = min_dim; dim <= max_dim; ++dim) {
    int dim_start = edges_.size();

    for (const S2ClippedShape& clipped : cell_->clipped_shapes()) {
      int shape_id = clipped.shape_id();
      const S2Shape& shape = *index->shape(shape_id);

      // Only process the current dimension.
      if (shape.dimension() != dim) {
        continue;
      }

      // In the event we wanted dimensions 0 and 2, but not 1.
      if (!dim_wanted(dim)) {
        continue;
      }

      // Materialize clipped shape edges into the edges vector. Track where we
      // start so we can add information about the region for this shape.
      int shape_start = edges_.size();
      for (int i = 0; i < clipped.num_edges(); ++i) {
        int edge_id = clipped.edge(i);

        // Looking up an edge requires looking up which chain it's in, which is
        // often a binary search.  So let's manually lookup the chain
        // information and use that to find the edge, so we only have to do that
        // search once.
        S2Shape::ChainPosition position = shape.chain_position(edge_id);
        S2Shape::Edge edge =
            shape.chain_edge(position.chain_id, position.offset);
        edges_.emplace_back(edge, edge_id, position);
      }

      // Note which block of edges belongs to the shape.
      shape_regions_.emplace_back(std::make_pair(
          shape_id, Region{shape_start, edges_.size() - shape_start}));
    }

    // Save region information for the current dimension.
    dim_regions_[dim] = {dim_start, edges_.size() - dim_start};
  }
}

absl::Span<const S2IndexCellData::EdgeAndIdChain> S2IndexCellData::shape_edges(
    int id) const {
  for (const auto& iter : shape_regions_) {
    if (iter.first == id) {
      Region region = iter.second;
      if (static_cast<size_t>(region.start) < edges_.size()) {
        return {&edges_[region.start], region.size};
      }
      return {};
    }
  }
  return {};
}

absl::Span<const S2IndexCellData::EdgeAndIdChain> S2IndexCellData::dim_edges(
    int dim) const {
  ABSL_DCHECK(0 <= dim && dim <= 2);
  const Region& region = dim_regions_[dim];
  if (static_cast<size_t>(region.start) < edges_.size()) {
    return {&edges_[region.start], region.size};
  }
  return {};
}

absl::Span<const S2IndexCellData::EdgeAndIdChain>
S2IndexCellData::dim_range_edges(int dim0, int dim1) const {
  ABSL_DCHECK_LE(dim0, dim1);
  ABSL_DCHECK(0 <= dim0 && dim0 <= 2);
  ABSL_DCHECK(0 <= dim1 && dim1 <= 2);

  size_t start = dim_regions_[dim0].start;
  size_t size = 0;
  for (int dim = dim0; dim <= dim1; ++dim) {
    start = std::min(start, static_cast<size_t>(dim_regions_[dim].start));
    size += dim_regions_[dim].size;
  }

  if (start < edges_.size()) {
    return {&edges_[start], size};
  }
  return {};
}

bool S2IndexCellData::ShapeContains(const S2ClippedShape& clipped,
                                    const S2Point& point,
                                    S2VertexModel model) const {
#ifndef NDEBUG
  R2Rect bounds = S2PaddedCell(id(), MutableS2ShapeIndex::kCellPadding).bound();

  R2Point uv;
  S2::ValidFaceXYZtoUV(S2Cell(id()).face(), point, &uv);
  ABSL_DCHECK(bounds.Contains(uv));
#endif

  // Points and polylines don't contain anything except under the CLOSED model.
  const S2Shape& shape = *index_->shape(clipped.shape_id());
  if (shape.dimension() < 2) {
    if (model != S2VertexModel::CLOSED) {
      return false;
    }

    // Point/polyline only contains point if it's a vertex.
    for (const EdgeAndIdChain& edge : shape_edges(clipped.shape_id())) {
      if (edge.v0 == point || edge.v1 == point) {
        return true;
      }
    }
    return false;
  }

  // Test containment by drawing a line segment from the cell center to the
  // given point and counting edge crossings.
  S2Point center = this->center();
  S2EdgeCrosser crosser(&center, &point);

  bool inside = clipped.contains_center();
  for (const EdgeAndIdChain& edge : shape_edges(clipped.shape_id())) {
    int sign = crosser.CrossingSign(&edge.v0, &edge.v1);
    if (sign < 0) continue;
    if (sign == 0) {
      // For the OPEN and CLOSED models, check whether point is a vertex.
      if (model != S2VertexModel::SEMI_OPEN &&
          (edge.v0 == point || edge.v1 == point)) {
        return (model == S2VertexModel::CLOSED);
      }
      sign = S2::VertexCrossing(*crosser.a(), *crosser.b(), edge.v0, edge.v1);
    }
    inside ^= sign;
  }
  return inside;
}

}  // namespace internal
