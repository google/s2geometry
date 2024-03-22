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

#ifndef S2_S2INDEX_CELL_DATA_H_
#define S2_S2INDEX_CELL_DATA_H_

#include <atomic>
#include <utility>
#include <vector>

#include "absl/base/thread_annotations.h"
#include "absl/container/flat_hash_set.h"
#include "absl/log/absl_check.h"
#include "absl/synchronization/mutex.h"
#include "absl/types/span.h"
#include "s2/s2cell.h"
#include "s2/s2contains_point_query.h"
#include "s2/s2shape.h"
#include "s2/s2shape_index.h"

// A class for working with S2ShapeIndexCell data.  For larger queries like
// validation, we often end up looking up edges multiple times, and sometimes
// need to work with the edges themselves, their edge ids, or their chain and
// offset.
//
// S2ShapeIndexCell and the S2ClippedShape APIs fundamentally work with edge ids
// and can't be re-worked without significant effort and loss of efficiency.
// This class provides an alternative API than repeatedly querying through the
// shapes in the index.
//
// This is meant to support larger querying and validation operations such as
// S2ValidationQuery that have to proceed cell-by cell through an index.
//
// To use, simply call LoadCell() to decode the contents of a cell.
//
// The class promises that the edges will be looked up once when LoadCell() is
// called, and the edges, edge_ids, chain, and chain offsets are loaded into a
// contiguous chunk of memory that we can serve requests from via absl::Span.
// Since the chain and offset are computed anyways when looking up an edge via
// the shape.edge() API, we simply cache those values so the cost is minimal.
//
// The memory layout looks like this:
//
//   |     0D Shapes     |     1D Shapes     |     2D Shapes     |  Dimensions
//   |  5  |   1   |  3  |  2  |   7   |  0  |  6  |   4   |  8  |  Shapes
//   [ ......................... Edges ..........................]  Edges
//
// This allows us to look up individual shapes very quickly, as well as all
// shapes in a given dimension or contiguous range of dimensions:
//
//   Edges()        - Return Span over all edges.
//   ShapeEdges()   - Return Span over edges of a given shape.
//   DimEdges()     - Return Span over all edges of a given dimension.
//   DimRangeEges() - Return Span over all edges of a range of dimensions.
//
// We use a stable sort, so similarly to S2ShapeIndexCell, we promise that
// shapes _within a dimension_ are in the same order they are in the index
// itself, and the edges _within a shape_ are similarly in the same order.
//
// The clipped shapes in a cell are exposed through the Shapes() method.
//
class S2IndexCellData {
 public:
  // Extension of Edge with fields for the edge id, chain id, and offset.  It's
  // useful to bundle these together when decoding S2ShapeIndex cells because it
  // allows us to avoid repetitive edge and chain lookups in many cases.
  struct EdgeAndIdChain : S2Shape::Edge {
    EdgeAndIdChain() = default;

    EdgeAndIdChain(const Edge& edge, int edge_id, int chain, int offset)
        : Edge(edge), id(edge_id), chain(chain), offset(offset) {}

    EdgeAndIdChain(const Edge& edge, int edge_id, S2Shape::ChainPosition pos)
        : Edge(edge), id(edge_id), chain(pos.chain_id), offset(pos.offset) {}

    friend bool operator==(const EdgeAndIdChain& x, const EdgeAndIdChain& y) {
      return x.v0 == y.v0 && x.v1 == y.v1;
    }

    friend bool operator<(const EdgeAndIdChain& x, const EdgeAndIdChain& y) {
      return x.v0 < y.v0 || (x.v0 == y.v0 && x.v1 < y.v1);
    }

    int32 id;      // Id of the edge within its shape.
    int32 chain;   // Id of the chain the edge belongs to.
    int32 offset;  // Offset of the edge within the chain.
  };

  S2IndexCellData() = default;

  S2IndexCellData(const S2ShapeIndex* index, S2CellId id,
                  const S2ShapeIndexCell* cell) {
    LoadCell(index, id, cell);
  }

  // Resets internal state to defaults without de-allocating memory.  The next
  // call to LoadCell() will process the cell regardless of whether it was
  // already loaded.  Should be called when processing a new index.
  void Reset() {
    index_ = nullptr;
    cell_ = nullptr;
    edges_.clear();
    shape_regions_.clear();
  }

  // Returns true if a dimension (0, 1, or 2) is set to be decoded.
  bool dim_wanted(int dim) const {
    ABSL_DCHECK(0 <= dim && dim <= 2);
    return dim_wanted_[dim];
  }

  // Configures whether a particular dimension of shape should be decoded.
  void set_dim_wanted(int dim, bool wanted) {
    ABSL_DCHECK(0 <= dim && dim <= 2);
    dim_wanted_[dim] = wanted;
  }

  // Returns the id of the current cell.
  S2CellId id() const { return cell_id_; }

  // Returns a const pointer to the index of the current cell.
  const S2ShapeIndex* index() const { return index_; }

  // Returns an S2Cell instance for the current cell.
  const S2Cell& cell() const {
    if (!s2cell_set_.load(std::memory_order_acquire)) {
      absl::MutexLock lock(&lock_);
      if (!s2cell_set_.load(std::memory_order_relaxed)) {
        s2cell_ = S2Cell(cell_id_);
        s2cell_set_.store(true, std::memory_order_release);
      }
    }
    // `s2cell_` is set once an for all, it won't change after this function
    // returns.
    return ABSL_TS_UNCHECKED_READ(s2cell_);
  }

  // Returns the center point of the current cell.
  const S2Point& center() const {
    if (!center_set_.load(std::memory_order_acquire)) {
      absl::MutexLock lock(&lock_);
      if (!center_set_.load(std::memory_order_relaxed)) {
        cell_center_ = cell_id_.ToPoint();
        center_set_.store(true, std::memory_order_release);
      }
    }
    // `cell_center_` is set once an for all, it won't change after this
    // function returns.
    return ABSL_TS_UNCHECKED_READ(cell_center_);
  }

  // Loads the data from the given cell, previous cell data is cleared.  Both
  // the index and cell lifetimes must span the lifetime until the class is
  // destroyed or LoadCell() is called again, so we take by const pointer to
  // avoid binding temporaries.
  //
  // If the index, id and cell pointer are the same as in the previous call to
  // LoadCell, loading is not performed since we already have the data decoded.
  void LoadCell(const S2ShapeIndex* index, S2CellId id,
                const S2ShapeIndexCell* cell);

  // Returns the number of clipped shapes in the cell.
  int num_clipped() const { return cell_->num_clipped(); }

  // Returns a reference to the S2Shape for the given clipped shape.  The
  // returned reference is to a shape owned by the index passed to LoadCell().
  const S2Shape& shape(const S2ClippedShape& clipped) const {
    const S2Shape* shape = index_->shape(clipped.shape_id());
    ABSL_DCHECK_NE(shape, nullptr);
    return *shape;
  }

  // Same as above but takes the shape id directly.
  const S2Shape& shape(int shape_id) const {
    const S2Shape* shape = index_->shape(shape_id);
    ABSL_DCHECK_NE(shape, nullptr);
    return *shape;
  }

  // Return a constant view over the clipped shapes in the cell.
  absl::Span<const S2ClippedShape> clipped_shapes() const {
    return cell_->clipped_shapes();
  }

  // Returns a view into all the edges in the current cell.
  absl::Span<const EdgeAndIdChain> edges() const { return edges_; }

  // Returns a view into the edges in the current cell for a given shape.
  // REQUIRES: id is an id for one of the shapes in Shapes()
  absl::Span<const EdgeAndIdChain> shape_edges(int shape_id) const;

  // Returns a view into the edges in the current cell for a given dimension.
  // REQUIRES: dim is a valid dimension (0, 1, or 2)
  absl::Span<const EdgeAndIdChain> dim_edges(int dim) const;

  // Returns a view into the edges in the current cellfor a range of dimensions.
  // REQUIRES: dimensions are valid (0, 1, or 2)
  absl::Span<const EdgeAndIdChain> dim_range_edges(int dim0, int dim1) const;

  // Tests whether a shape in the current cell contains a point.  The logic here
  // is is the same as S2ContainsPointQuery::ShapeContains but only applies to
  // the current cell and doesn't have to lookup edges or cell centers again.
  //
  // Since this _only_ operates on edges within the current cell, it must not
  // be used to test points that are outside of the current cell, as there may
  // be intervening edges between the point and the cell center that we can't
  // see.
  //
  // By default operates with an OPEN vertex model.
  //
  // REQUIRES: The cell must contain point.
  bool ShapeContains(const S2ClippedShape& clipped, const S2Point& point,
                     S2VertexModel model = S2VertexModel::OPEN) const;

 private:
  const S2ShapeIndex* index_ = nullptr;
  const S2ShapeIndexCell* cell_ = nullptr;
  S2CellId cell_id_;

  // Computing the cell center and S2Cell can cost as much as looking up the
  // edge themselves, so defer doing it until needed.  We want to be able to
  // access the values through a const reference though, so we have to
  // synchronize access ourselves.
  //
  // TODO: We can use std::atomic_flag instead when C++20 is allowed.
  mutable std::atomic<bool> s2cell_set_ = false;
  mutable std::atomic<bool> center_set_ = false;
  mutable absl::Mutex lock_;
  mutable S2Cell s2cell_ ABSL_GUARDED_BY(lock_);
  mutable S2Point cell_center_ ABSL_GUARDED_BY(lock_);

  // Dimensions that we wish to decode, the default is all of them.
  bool dim_wanted_[3] = {true, true, true};

  // Storage space for edges of the current cell.
  std::vector<EdgeAndIdChain> edges_;

  // Simple pair for defining an integer valued region.
  struct Region {
    int start = 0;
    size_t size = 0;
  };

  // Map from shape id to the region of the edges_ array it's stored in.
  std::vector<std::pair<int, Region>> shape_regions_;

  // Region for each dimension we might encounter.
  Region dim_regions_[3];
};

#endif  // S2_S2INDEX_CELL_DATA_H_
