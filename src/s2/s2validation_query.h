// Copyright 2022 Google Inc. All Rights Reserved.
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


#ifndef S2_S2VALIDATION_QUERY_H_
#define S2_S2VALIDATION_QUERY_H_

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <utility>
#include <vector>

#include "absl/algorithm/container.h"
#include "absl/container/flat_hash_set.h"
#include "absl/container/inlined_vector.h"
#include "absl/log/absl_check.h"
#include "absl/strings/str_format.h"
#include "absl/types/span.h"
#include "s2/internal/s2disjoint_set.h"
#include "s2/internal/s2incident_edge_tracker.h"
#include "s2/internal/s2index_cell_data.h"
#include "s2/s2cell_id.h"
#include "s2/s2contains_point_query.h"
#include "s2/s2contains_vertex_query.h"
#include "s2/s2edge_crosser.h"
#include "s2/s2error.h"
#include "s2/s2point.h"
#include "s2/s2pointutil.h"
#include "s2/s2predicates.h"
#include "s2/s2shape.h"
#include "s2/s2shape_index.h"
#include "s2/s2shapeutil_edge_wrap.h"

// This header defines the following queries which each have slightly different
// correctness semantics for their particular domain.

// clang-format off
template <typename IndexType> class S2ValidQuery;
template <typename IndexType> class S2LegacyValidQuery;
// clang-format on

// Base class for validating geometry contained in an index according to a
// configurable model of correctness.  There are several different notions of
// "valid" geometry that we have to address, including the basic requirements
// for S2BooleanOperation, S2Polygon specific checks, OGC simple geometry
// standards, and STLib's requirements.
//
// Classes of geometry can be thought of as sets of all shapes that have certain
// invariants (such as interior-on-left or no-degenerate-edges).  Classes then
// naturally build on other classes by adding more rules to further restrict the
// allowed shapes.
//
// This extends-upon structure naturally lends itself to an inheritance based
// implementation.  Beginning with the most permissive class of geometry which
// just meets the criteria for S2BooleanOperation (S2ValidQuery), we can then
// build subsequent classes of geometry by inheriting from and extending the
// checks that are performed.
//
// Validation queries should be templated over the exact index type and can
// overload the virtual methods in the subclass API below.
//
// The Validate() driver function will call these functions in this order:
//
//   - Start()
//     - CheckShape()
//     - StartCell()
//       - StartShape()
//         - CheckEdge()
//       - FinishShape()
//   - Finish()
//
// Hoisting the loops into the base class like this allows us to fuse all the
// inner loops of the various queries so that we only have to iterate over
// the index and its geometry once.
//
// Several protected member functions are defined to give subclasses access to
// the data being validated:
//
//   const IndexType& Index();
//     -- Returns a reference to the current index being operated on.
//
//   const S2IndexCellData& CurrentCell();
//     -- Returns a reference to the decoded data for the current cell.
//
//   const IncidentEdgeSet& IncidentEdges();
//     -- Returns a reference to the current incident edge set.  We promise that
//     this set is updated with the current cell's edges before StartCell() is
//     called.
//
// A query then has a `bool Validate(const IndexType& index, S2Error* error)`
// method which validates the index and returns true if it's valid, otherwise
// false, with the validation failure details provided through the error
// parameter.
//
// This example validates an index as containing valid geometry for
// use with S2Polygon/S2Polyline:
//
//   S2LegacyValidQuery<MutableS2ShapeIndex> query;
//   if (!query.Validate(index, error)) {
//     ...
//   }
//
template <typename IndexType>
class S2ValidationQueryBase {
 public:
  using IncidentEdgeSet = internal::IncidentEdgeSet;
  using EdgeAndIdChain = typename internal::S2IndexCellData::EdgeAndIdChain;
  using S2IndexCellData = internal::S2IndexCellData;

  S2ValidationQueryBase() = default;
  virtual ~S2ValidationQueryBase() = default;

  // Validate the index by calling the hooks in the derived class.
  bool Validate(const IndexType& index, S2Error* error);

 protected:
  using Iterator = typename IndexType::Iterator;

  // Subclass API
  // Starts the validation process; called once per query.
  virtual bool Start(S2Error*) { return true; }

  // Validates each individual shape in the index; called once per shape.
  //
  // A reference to the current iterator state is passed in.  The function may
  // reposition the iterator in order to do shape checking.
  virtual bool CheckShape(Iterator& iter, const S2Shape& shape, int shape_id,
                          S2Error*) {
    return true;
  }

  // Starts processing of a cell in the index; called once per cell.
  virtual bool StartCell(S2Error*) { return true; }

  // Marks start of a clipped shape; called once per clipped shape in a cell.
  virtual bool StartShape(const S2Shape&, const S2ClippedShape&, S2Error*) {
    return true;
  }

  // Validates a single edge of a given shape; called at least once per edge of
  // each shape.
  virtual bool CheckEdge(const S2Shape&, const S2ClippedShape&,
                         const EdgeAndIdChain&, S2Error*) {
    return true;
  }

  // Marks end of a clipped shape; called once per clipped shape in a cell.
  virtual bool FinishShape(const S2Shape&, const S2ClippedShape&, S2Error*) {
    return true;
  }

  // Marks end of validation; called once per query.
  virtual bool Finish(S2Error* error) { return true; }

  // Returns a reference to the index we're validating.
  const IndexType& Index() const {
    ABSL_DCHECK(index_ != nullptr);
    return *index_;
  }

  // Returns current incident edge map.
  const IncidentEdgeSet& IncidentEdges() const {
    return incident_edge_tracker_.IncidentEdges();
  }

  // Returns a reference to the data for the current cell.
  const S2IndexCellData& CurrentCell() const { return cell_buffer_; }

 private:
  // Disallow copying and assignment through an abstract reference.
  S2ValidationQueryBase(const S2ValidationQueryBase&) = delete;
  S2ValidationQueryBase& operator=(const S2ValidationQueryBase&) = delete;

  // Sets the current index we're operating on, accessible through Index();
  void SetIndex(const IndexType* index) { index_ = index; }

  // Sets the current cell that we're operating on and loads its data.
  void SetCurrentCell(const Iterator& iter) {
    cell_buffer_.LoadCell(index_, iter.id(), &iter.cell());
  }

  S2IndexCellData cell_buffer_;
  internal::S2IncidentEdgeTracker incident_edge_tracker_;
  const IndexType* index_ = nullptr;
};

// This class represents the least strict class of geometry and corresponds to
// the requirements for compatibility with S2BooleanOperation, with these
// specific constraints:
//
// == General ==
//   * Points must be unit magnitude (according to S2::IsUnitLength()).
//   * Points must be finite (neither inf or nan).
//   * Edges must not have antipodal vertices.
//   * Degenerate edges of the form {A,A} are allowed by default.
//   * Reverse-duplicate edges of the form {A,B},{B,A} are allowed by default.
//   * Shape chains must have at least one edge, excluding the full and empty
//     polygon shapes, which may have at most one empty chain.
//
// == Polygons ==
//   * Polygon interiors must be disjoint from all other geometry.
//     * Polygon edges thus cannot cross any other edges, including at vertices.
//     * No geometry may test as contained in another polygon.
//   * Polygons can never have duplicate edges even among different polygons.
//   * Polygon edges must be connected and form a closed loop per chain.
//   * Polygon chains must oriented so that the interior is always on the left.
//
// == Polylines ==
//   * Polyline edges can cross by default.
//   * Polylines can have duplicate edges by default.
//
template <typename IndexType>
class S2ValidQuery : public S2ValidationQueryBase<IndexType> {
 protected:
  // Protected because options are only for subclasses to configure behavior.
  class Options {
   public:
    // Types of single-vertex touches allowed between shapes.  This is
    // configurable for each combination of dimension.  E.g. we can configure
    // polyline-polyline (1D-1D) and polyline-polygon (1D-2D) touches
    // separately.
    enum TouchType {
      kTouchNone = 0b00,      // Allow no touches between shapes.
      kTouchInterior = 0b01,  // Interior point may touch the other shape.
      kTouchBoundary = 0b10,  // Boundary point may touch the other shape.
      kTouchAny = 0b11        // Allow any touches between shapes.
    };

    Options() {
      // Default is to allow any touches between geometry.
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          set_allowed_touches(i, j, std::make_pair(kTouchAny, kTouchAny));
        }
      }
    }

    // Returns a TouchType with fields set indicating what types of vertex from
    // dimension A are allowed to touch vertices from dimension B.
    std::pair<TouchType, TouchType> allowed_touches(int dima, int dimb) const {
      ABSL_DCHECK_GE(dima, 0);
      ABSL_DCHECK_LE(dima, 2);
      ABSL_DCHECK_GE(dimb, 0);
      ABSL_DCHECK_LE(dimb, 2);

      // Just get the top half of the matrix.
      if (dima > dimb) {
        std::swap(dima, dimb);
      }

      return allowed_touches_[dima][dimb];
    }

    // Set the allowed touches based on geometry dimension.  We require that the
    // matrix of combinations be symmetric, so set_allowed_touches(1, 2, X, Y)
    // and set_allowed_touches(2, 1, Y, X) are equivalent.
    Options& set_allowed_touches(  //
        int dima, int dimb, std::pair<TouchType, TouchType> types) {
      ABSL_DCHECK(0 <= dima && dima <= 2);
      ABSL_DCHECK(0 <= dimb && dimb <= 2);

      // Just set the top half of the matrix.
      if (dima > dimb) {
        std::swap(dima, dimb);
      }

      allowed_touches_[dima][dimb] = types;
      return *this;
    }

    // Configures whether polylines can have duplicate edges.
    bool allow_duplicate_polyline_edges() const {
      return allow_duplicate_polyline_edges_;
    }

    Options& set_allow_duplicate_polyline_edges(bool flag) {
      allow_duplicate_polyline_edges_ = flag;
      return *this;
    }

    // Configures whether polyline edges can cross.
    bool allow_polyline_interior_crossings() const {
      return allow_polyline_interior_crossings_;
    }

    Options& set_allow_polyline_interior_crossings(bool flag) {
      allow_polyline_interior_crossings_ = flag;
      return *this;
    }

    // Configures whether reverse duplicate edges are allowed.
    bool allow_reverse_duplicates() const { return allow_reverse_duplicates_; }

    Options& set_allow_reverse_duplicates(bool flag) {
      allow_reverse_duplicates_ = flag;
      return *this;
    }

    // Configures whether degenerate (zero-length) edges are allowed.
    bool allow_degenerate_edges() const { return allow_degenerate_edges_; }

    Options& set_allow_degenerate_edges(bool flag) {
      allow_degenerate_edges_ = flag;
      return *this;
    }

   private:
    bool allow_degenerate_edges_ = true;
    bool allow_duplicate_polyline_edges_ = true;
    bool allow_reverse_duplicates_ = true;
    bool allow_polyline_interior_crossings_ = true;

    std::pair<TouchType, TouchType> allowed_touches_[3][3];
  };

  const Options& options() const { return options_; }
  Options& mutable_options() { return options_; }

  // We have to include unqualified names manually due to template inheritance.
  using Base = S2ValidationQueryBase<IndexType>;
  using Base::CurrentCell;
  using Base::IncidentEdges;
  using Base::Index;
  using typename Base::EdgeAndIdChain;
  using typename Base::Iterator;
  using typename Base::S2IndexCellData;

  bool CheckShape(Iterator& iter, const S2Shape& shape, int shape_id,
                  S2Error*) override;
  bool StartCell(S2Error*) override;
  bool CheckEdge(const S2Shape& shape, const S2ClippedShape& clipped,
                 const EdgeAndIdChain& edge, S2Error*) override;
  bool Finish(S2Error* error) override;

 private:
  // Returns true if the given S2Point is valid, meaning none of its components
  // are inf or NaN.
  static inline bool ValidPoint(const S2Point& p) {
    return std::isfinite(p.x()) && std::isfinite(p.y()) && std::isfinite(p.z());
  }

  // Checks that a point is fully outside any polygons.
  //
  // Returns false if a point is not interior to any polygons, true otherwise
  // with error populated.
  bool PointContained(S2CellId, int shape_id, const S2Point&, S2Error*);

  // Checks if any edge in the current cell is a duplicate (or reverse
  // duplicate).
  bool CheckForDuplicateEdges(S2Error* error) const;

  // Checks if any edges in the current cell have an interior crossing.
  bool CheckForInteriorCrossings(S2Error* error) const;

  // Checks that a given chain of a polygon is oriented properly relative to one
  // cell center in the index.  We only need to check one center because we
  // assume that the transition between cell indices is consistent, thus we just
  // need to make sure that the interior state isn't flipped.
  //
  // Returns false when the chain isn't properly oriented with S2Error set with
  // the details, true otherwise.
  //
  // REQUIRES: Shape must be two dimensional.
  // REQUIRES: Chain must be closed.
  // REQUIRES: Chain edges must be connected (no gaps).
  // REQUIRES: Chain must have at least three distinct points (cover some area).
  bool CheckChainOrientation(  //
      Iterator&, const S2Shape&, int shape_id, int chain_id, S2Error*);

  // Information on one vertex of an edge, including whether it's on the
  // boundary of its shape.  This is used by CheckTouchesAreValid() when we need
  // to check that vertex touches are valid.
  struct TestVertex {
    S2Point vertex;
    int32_t edge_id;
    int32_t shape_id;
    int32_t dim;
    bool on_boundary;
  };
  std::vector<TestVertex> test_vertices_;

  // Checks the shapes in the current cell to ensure that any touch points are
  // allowed under the configured semantics.
  //
  // Returns false if any vertices touch at an invalid point with the error
  // value set, true otherwise.
  bool CheckTouchesAreValid(S2Error*);

  Options options_;
};

//////////////////   Utility Implementation   ////////////////////

// Sorts a container of S2Shape::Edges in counter-clockwise order around an
// origin point.  By itself, ordering this way is ambiguous in terms of the
// starting point, so we take an edge to form a reference point to form a total
// ordering.
//
// Reverse duplicate edges are ordered so that the one with origin as v0 comes
// before the other.
//
// The edges and first edge must all contain origin as one of their vertices.
//
// Can sort any container with elements that are convertible to an S2Shape::Edge
// reference.  Both the origin point and first edge are taken by value so that
// it's safe to call using elements of the container for those values.
//
// Takes a Container to the start and end of the container.  std::sort requires
// random iterators and thus so do we.
template <typename Container>
void SortEdgesCcw(S2Point origin, S2Shape::Edge first, Container& data) {
  ABSL_DCHECK(first.v0 == origin || first.v1 == origin);
  const S2Point& first_vertex = (first.v0 == origin) ? first.v1 : first.v0;
  ABSL_DCHECK(first_vertex != origin);

  absl::c_sort(  //
      data,      //
      [&](const S2Shape::Edge& a, const S2Shape::Edge& b) {
        ABSL_DCHECK(a.v0 == origin || a.v1 == origin);
        ABSL_DCHECK(b.v0 == origin || b.v1 == origin);

        // `OrderedCCW` will return `true` if `a == b`, which will violate the
        // irreflexivity requirement of a strict weak ordering.
        if (a == b) {
          return false;
        }

        // Order reverse duplicates so that the one with edge.v0 ==
        // origin is first.
        if (a == b.Reversed()) {
          return a.v0 == origin;
        }

        if (a == first || b == first) {
          // If either edge is the first edge, compare so that the
          // first edge always comes first in the sorted output.
          return a == first;
        }

        // Otherwise check orientation of vertices.
        S2Point apnt = (a.v0 == origin) ? a.v1 : a.v0;
        S2Point bpnt = (b.v0 == origin) ? b.v1 : b.v0;
        return s2pred::OrderedCCW(first_vertex, apnt, bpnt, origin);
      });
}

template <typename IndexType>
bool S2ValidationQueryBase<IndexType>::Validate(const IndexType& index,
                                                S2Error* error) {
  SetIndex(&index);
  incident_edge_tracker_.Reset();
  cell_buffer_.Reset();

  if (!Start(error)) {
    return false;
  }

  // Run basic checks on individual shapes in the index.
  Iterator iter(&index, S2ShapeIndex::BEGIN);
  for (int shape_id = 0; shape_id < index.num_shape_ids(); ++shape_id) {
    const S2Shape* shape = index.shape(shape_id);
    if (shape != nullptr) {
      if (!CheckShape(iter, *shape, shape_id, error)) {
        return false;
      }
    }
  }

  for (iter.Begin(); !iter.done(); iter.Next()) {
    SetCurrentCell(iter);

    // Add two dimensional shape edges to the incident edge tracker to support
    // checks for things like crossing polygon chains and split interiors.
    for (const S2ClippedShape& clipped : CurrentCell().clipped_shapes()) {
      const int shape_id = clipped.shape_id();
      const S2Shape& shape = CurrentCell().shape(clipped);
      if (shape.dimension() < 2) {
        continue;
      }

      incident_edge_tracker_.StartShape(shape_id);
      for (const auto& edge : CurrentCell().shape_edges(shape_id)) {
        incident_edge_tracker_.AddEdge(edge.id, edge);
      }
      incident_edge_tracker_.FinishShape();
    }

    // Now notify that we're starting this cell.
    if (!StartCell(error)) {
      return false;
    }

    // Iterate the shapes and edges of the cell.
    for (const S2ClippedShape& clipped : CurrentCell().clipped_shapes()) {
      const int shape_id = clipped.shape_id();
      const S2Shape& shape = CurrentCell().shape(clipped);
      if (!StartShape(shape, clipped, error)) {
        return false;
      }

      for (const auto& edge : CurrentCell().shape_edges(shape_id)) {
        if (!CheckEdge(shape, clipped, edge, error)) {
          return false;
        }
      }

      if (!FinishShape(shape, clipped, error)) {
        return false;
      }
    }
  }

  // Run any final checks and finish validation
  return Finish(error);
}

// This class represents a semantic class of geometry that's compatible with the
// requirements of the S2Polygon and S2Polyline IsValid() methods.  It extends
// the S2ValidQuery requirements, specifically:
//
// == General ==
//   * Degenerate edges of the form {A,A} are not allowed.
//   * All the shapes in an S2ShapeIndex must be the same dimensionality.
//   * Duplicate vertices in a chain are not allowed.
//     I.e. a chain cannot touch itself even at one point.
//   * Different chains may touch, but only in such a way they don't create
//     duplicate edges.
//
template <typename IndexType>
class S2LegacyValidQuery final : public S2ValidQuery<IndexType> {
 public:
  S2LegacyValidQuery();

 protected:
  // We have to include unqualified names manually due to template inheritance.
  using Base = S2ValidQuery<IndexType>;
  using Base::CurrentCell;
  using Base::IncidentEdges;
  using Base::Index;
  using typename Base::EdgeAndIdChain;
  using typename Base::Iterator;
  using typename Base::S2IndexCellData;

  bool Start(S2Error*) final;
  bool CheckShape(Iterator&, const S2Shape& shape, int shape_id,
                  S2Error*) final;
  bool StartCell(S2Error*) final;
  bool CheckEdge(const S2Shape& shape, const S2ClippedShape& clipped,
                 const EdgeAndIdChain& edge, S2Error*) override;

 private:
  // Tuple of (shape, chain, vertex) for detecting duplicate vertices in the
  // same chain.
  struct ShapeChainVertex {
    // Constructor just to twiddle order of fields.
    ShapeChainVertex(int shape, int chain, S2Point vertex)
        : vertex(vertex), chain(chain), shape(shape) {}

    S2Point vertex;
    int chain;
    int shape;

    // Makes type hashable for use in Abseil containers.
    template <typename H>
    friend H AbslHashValue(H h, const ShapeChainVertex& a) {
      return H::combine(std::move(h), a.vertex, a.shape, a.chain);
    }
  };
};

//////////////////   S2ValidQuery Implementation   ////////////////////

template <typename IndexType>
bool S2ValidQuery<IndexType>::CheckShape(Iterator& iter, const S2Shape& shape,
                                         int shape_id, S2Error* error) {
  // Verify that shape isn't outright lying to us about its dimension.
  const int dim = shape.dimension();
  if (dim < 0 || dim > 2) {
    *error = S2Error(
        S2Error::INVALID_DIMENSION,
        absl::StrFormat("Shape %d has invalid dimension: %d", shape_id, dim));
    return false;
  }

  absl::InlinedVector<int, 4> chains_to_check;
  for (int chain_id = 0; chain_id < shape.num_chains(); ++chain_id) {
    const S2Shape::Chain& chain = shape.chain(chain_id);

    // Check that the first and last edges in a polygon chain connect to close
    // the chain.  This is true even for degenerate chains with one edge that's
    // a point, or two edges that are reverse duplicates.  Both are considered
    // closed.
    if (dim == 2 && chain.length > 0) {
      int edge_id = chain.start;
      int prev_id = s2shapeutil::PrevEdgeWrap(shape, edge_id);

      if (shape.edge(prev_id).v1 != shape.edge(edge_id).v0) {
        *error = S2Error(S2Error::LOOP_NOT_ENOUGH_VERTICES,
                         absl::StrFormat("Chain %d of shape %d isn't closed",
                                         chain_id, shape_id));
        return false;
      }
    }

    for (int offset = 0; offset < chain.length; ++offset) {
      const S2Shape::Edge& edge = shape.chain_edge(chain_id, offset);

      // Check that coordinates aren't inf/nan.
      if (!ValidPoint(edge.v0) || !ValidPoint(edge.v1)) {
        *error = S2Error(
            S2Error::INVALID_VERTEX,
            absl::StrFormat("Shape %d has invalid coordinates", shape_id));
        return false;
      }

      // Check that vertices are unit length.
      if (!S2::IsUnitLength(edge.v0) || !S2::IsUnitLength(edge.v1)) {
        *error = S2Error(
            S2Error::NOT_UNIT_LENGTH,
            absl::StrFormat("Shape %d has non-unit length vertices", shape_id));
        return false;
      }

      // (Optional) check polyline and polygon edges for degeneracy.
      if (dim > 0 && !options().allow_degenerate_edges()) {
        if (edge.IsDegenerate()) {
          *error = S2Error(
              S2Error::DUPLICATE_VERTICES,
              absl::StrFormat("Shape %d: chain %d, edge %d is degenerate",
                              shape_id, chain_id, chain.start + offset));
          return false;
        }
      }

      // Check that edge doesn't have antipodal vertices.
      if (edge.v0 == -edge.v1) {
        *error =
            S2Error(S2Error::ANTIPODAL_VERTICES,
                    absl::StrFormat("Shape %d has adjacent antipodal vertices",
                                    shape_id));
        return false;
      }

      // Check that chain edges are connected for polylines and polygons.
      if (dim > 0 && chain.length >= 2 && offset > 0) {
        const S2Shape::Edge last = shape.chain_edge(chain_id, offset - 1);
        if (last.v1 != edge.v0) {
          *error =
              S2Error(S2Error::NOT_CONTINUOUS,
                      absl::StrFormat("Chain %d of shape %d has neighboring "
                                      "edges that don't connect.",
                                      chain_id, shape_id));
          return false;
        }
      }
    }

    // The rest of the checks are for non-empty polygon chains only.
    if (dim != 2 || chain.length == 0) {
      continue;
    }

    // We need at least two distinct points in a chain before we can check its
    // orientation vs the cell center.  Scan until we find a vertex different
    // than the first.
    int unique_count = 1;
    const S2Point first = shape.chain_edge(chain_id, 0).v0;
    for (const S2Point& vertex : shape.vertices(chain_id)) {
      if (vertex != first) {
        ++unique_count;
        break;
      }
    }

    // Only a single unique point. A degenerate edge will never test as a vertex
    // crossing (because 3 out of 4 vertices to S2::VertexCrossing() would be
    // equal making it false), so they can't toggle interior state and we can
    // ignore them.
    if (unique_count == 1) {
      continue;
    }

    chains_to_check.emplace_back(chain_id);
  }

  // Check the selected chain to verify chain orientation.
  for (int chain_id : chains_to_check) {
    if (!CheckChainOrientation(iter, shape, shape_id, chain_id, error)) {
      return false;
    }
  }
  return true;
}

template <typename IndexType>
bool S2ValidQuery<IndexType>::CheckForDuplicateEdges(S2Error* error) const {
  int dim0 = options().allow_duplicate_polyline_edges() ? 2 : 1;
  int dim1 = 2;

  // This is O(N^2) but cells don't have many edges in them and benchmarks show
  // this to be faster than trying to sort the whole array ahead of time.
  absl::Span<const EdgeAndIdChain> edges =
      CurrentCell().dim_range_edges(dim0, dim1);
  int num_edges = edges.size();
  for (int i = 0; i < num_edges; ++i) {
    for (int j = i + 1; j < num_edges; ++j) {
      bool duplicate = (edges[i] == edges[j]);

      if (!options().allow_reverse_duplicates()) {
        duplicate |= (edges[i].Reversed() == edges[j]);
      }

      if (duplicate) {
        *error = S2Error(S2Error::OVERLAPPING_GEOMETRY,
                         "One or more duplicate polygon edges detected");
        return false;
      }
    }
  }

  return true;
}

template <typename IndexType>
bool S2ValidQuery<IndexType>::CheckForInteriorCrossings(S2Error* error) const {
  // Get all the polyline and polygon edges.
  absl::Span<const EdgeAndIdChain> edges = CurrentCell().dim_range_edges(1, 2);

  // If we're allowing polyline edges to cross polyline edges, then we only have
  // to check against polygon edges.
  size_t check_start = 0;
  if (options().allow_polyline_interior_crossings()) {
    check_start = CurrentCell().dim_edges(1).size();
  }

  // This can happen when we're allowing polyline crossings and only have
  // polylines, there's nothing to do.
  if (check_start >= edges.size()) {
    return true;
  }

  const size_t num_edges = edges.size();
  for (size_t i = 0; i + 1 < num_edges; ++i) {
    // We never have to check against edges at a lower index, because if we
    // intersect them, we'll have already checked that at this point.
    size_t j = std::max(check_start, i + 1);

    // We can skip adjacent edges.
    if (edges[i].v1 == edges[j].v0) {
      if (++j >= num_edges) {
        break;
      }
    }

    S2EdgeCrosser crosser(&edges[i].v0, &edges[i].v1);
    for (; j < num_edges; ++j) {
      if (crosser.c() == nullptr || *crosser.c() != edges[j].v0) {
        crosser.RestartAt(&edges[j].v0);
      }

      if (crosser.CrossingSign(&edges[j].v1) > 0) {
        *error =
            S2Error(S2Error::OVERLAPPING_GEOMETRY,
                    absl::StrFormat("Chain %d edge %d crosses chain %d edge %d",
                                    edges[i].chain, edges[i].offset,
                                    edges[j].chain, edges[j].offset));
        return false;
      }
    }
  }
  return true;
}

// Returns true if a vertex (0 or 1) of an edge of a polyline is a boundary
// point or not.  This is only true if the vertex is either the start or end
// point of a chain and the chain is open.  Returns false otherwise.
//
// REQUIRES: vertex == 0 or vertex == 1
inline bool PolylineVertexIsBoundaryPoint(
    const S2Shape& shape, const internal::S2IndexCellData::EdgeAndIdChain& edge,
    int vertex) {
  ABSL_DCHECK(vertex == 0 || vertex == 1);

  if (edge.offset == 0) {
    return s2shapeutil::PrevEdgeWrap(shape, edge.id) == -1 && vertex == 0;
  } else if (edge.offset == shape.chain(edge.chain).length - 1) {
    return s2shapeutil::NextEdgeWrap(shape, edge.id) == -1 && vertex == 1;
  }

  return false;
}

template <typename IndexType>
bool S2ValidQuery<IndexType>::CheckTouchesAreValid(S2Error* error) {
  using TouchType = typename Options::TouchType;

  const auto kAny = Options::kTouchAny;

  bool need_dim[3] = {true, true, true};
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j) {
      bool any_allowed =
          (options().allowed_touches(i, j) == std::make_pair(kAny, kAny));
      need_dim[i] &= !any_allowed;
    }
  }

  // If all touches are allowed, then we don't have to check, easy.
  if (!need_dim[0] && !need_dim[1] && !need_dim[2]) {
    return true;
  }

  // Gather unique vertices from each dimension that needs to be checked.
  test_vertices_.clear();
  for (const S2ClippedShape& clipped : CurrentCell().clipped_shapes()) {
    const int shape_id = clipped.shape_id();
    const S2Shape& shape = *Index().shape(shape_id);
    const int dim = shape.dimension();

    // Skip if we're not checking this dimension.
    if (!need_dim[dim]) {
      continue;
    }

    for (const EdgeAndIdChain& edge : CurrentCell().shape_edges(shape_id)) {
      // For polylines, we have to handle the start and ending edges specially,
      // since we can have open chains that have boundary points.
      if (dim == 1) {
        bool on_boundary = PolylineVertexIsBoundaryPoint(shape, edge, 0);
        test_vertices_.push_back(
            {edge.v0, edge.id, shape_id, dim, on_boundary});

        // Check vertex 1 too only if it's a boundary point, otherwise we would
        // test v1 twice after we grab v0 of the next edge.
        on_boundary = PolylineVertexIsBoundaryPoint(shape, edge, 1);
        if (on_boundary) {
          test_vertices_.push_back({edge.v1, edge.id, shape_id, dim, true});
        }
      } else {
        // All polygon vertices are on the boundary, all point vertices are not.
        test_vertices_.push_back({edge.v0, edge.id, shape_id, dim, dim == 2});
      }
    }
  }

  // For each test vertex, scan over the other shapes and verify that any vertex
  // touches are allowed under the current semantics.
  for (const TestVertex& testpnt : test_vertices_) {
    for (const S2ClippedShape& clipped : CurrentCell().clipped_shapes()) {
      const int shape_id = clipped.shape_id();
      const S2Shape& shape = *Index().shape(shape_id);
      const int dim = shape.dimension();

      for (const EdgeAndIdChain& edge : CurrentCell().shape_edges(shape_id)) {
        // Don't compare an edge against itself.
        if (testpnt.shape_id == shape_id && testpnt.edge_id == edge.id) {
          continue;
        }

        // Figure out which (if any) vertex of this edge we touched.
        int vertidx = -1;
        if (testpnt.vertex == edge.v0) vertidx = 0;
        if (testpnt.vertex == edge.v1) vertidx = 1;
        if (vertidx < 0) {
          continue;
        }

        // Closed polylines are always allowed, so if we get a hit on the same
        // shape and it's a polyline, check that it's not from the closure.
        if (testpnt.shape_id == shape_id && dim == 1) {
          if (vertidx == 0 &&
              s2shapeutil::PrevEdgeWrap(shape, edge.id) == testpnt.edge_id) {
            continue;
          }

          if (vertidx == 1 &&
              s2shapeutil::NextEdgeWrap(shape, edge.id) == testpnt.edge_id) {
            continue;
          }
        }

        // Points vertices are always interior, Polygon vertices are always
        // boundary, and it may be either for polylines.
        bool on_boundary = (dim == 2);
        if (dim == 1) {
          on_boundary = PolylineVertexIsBoundaryPoint(shape, edge, vertidx);
        }

        TouchType typea = Options::kTouchInterior;
        if (testpnt.on_boundary) {
          typea = Options::kTouchBoundary;
        }

        TouchType typeb = Options::kTouchInterior;
        if (on_boundary) {
          typeb = Options::kTouchBoundary;
        }

        using TouchPair = std::pair<TouchType, TouchType>;
        TouchPair allowed = options().allowed_touches(testpnt.dim, dim);

        // Returns true if both touches match the given allowed touches.
        const auto PermittedTouches = [](TouchPair allowed, TouchType typea,
                                         TouchType typeb) {
          return (allowed.first & typea) && (allowed.second & typeb);
        };

        // We require that touches be symmetric, so that the touch is invalid
        // only if it's invalid from A->B and B->A.  This lets us support
        // behavior such as requiring that one or the other point be an interior
        // point.
        if (!PermittedTouches(allowed, typea, typeb) &&
            !PermittedTouches(allowed, typeb, typea)) {
          *error = S2Error(S2Error::OVERLAPPING_GEOMETRY,
                           "Index has geometry with invalid vertex touches.");
          return false;
        }
      }
    }
  }
  return true;
}

template <typename IndexType>
bool S2ValidQuery<IndexType>::StartCell(S2Error* error) {
  if (!CheckForDuplicateEdges(error) || !CheckForInteriorCrossings(error)) {
    return false;
  }

  if (!CheckTouchesAreValid(error)) {
    return false;
  }

  return true;
}

template <typename IndexType>
bool S2ValidQuery<IndexType>::PointContained(S2CellId cell_id, int shape_id,
                                             const S2Point& point,
                                             S2Error* error) {
  const S2IndexCellData& cell = CurrentCell();
  for (const S2ClippedShape& clipped : cell.clipped_shapes()) {
    if (clipped.shape_id() == shape_id) {
      continue;
    }

    const S2Shape& shape = *Index().shape(clipped.shape_id());

    // Only check for containment in polygons.
    if (shape.dimension() != 2) {
      continue;
    }

    if (cell.ShapeContains(clipped, point)) {
      *error = S2Error(
          S2Error::OVERLAPPING_GEOMETRY,
          absl::StrFormat(
              "Shape %d has one or more edges contained in another shape.",
              shape_id));
      return true;
    }
  }

  return false;
}

template <typename IndexType>
bool S2ValidQuery<IndexType>::CheckChainOrientation(Iterator& iter,
                                                    const S2Shape& shape,
                                                    int shape_id, int chain_id,
                                                    S2Error* error) {
  const S2Shape::Chain& chain = shape.chain(chain_id);

  // Given that:
  //  1. Edges in the chain are connected continuously.
  //  2. The chain is closed.
  //  3. The chain at least two distinct points.
  //
  // Then we can test whether the chain is oriented properly relative to the
  // cell center by testing one edge of the chain for proper orientation.

  S2ContainsVertexQuery query;
  for (int offset = 0; offset < chain.length; ++offset) {
    const S2Point& vertex = shape.chain_edge(chain_id, offset).v0;
    query.Init(vertex);

    // Seek to the cell containing vertex and get the clipped shape.
    if (!iter.Locate(vertex)) {
      *error = S2Error::DataLoss("Shape vertex was not indexed");
      return false;
    }
    const S2Point& center = iter.id().ToPoint();

    const S2ClippedShape* clipped = iter.cell().find_clipped(shape_id);
    ABSL_DCHECK_NE(clipped, nullptr);

    // Compute winding number and vertex sign at the same time.
    int winding = clipped->contains_center();
    S2CopyingEdgeCrosser crosser(center, vertex);

    for (int i = 0; i < clipped->num_edges(); ++i) {
      const S2Shape::Edge& edge = shape.edge(clipped->edge(i));

      // Tally up the total change in winding number from center to vertex.
      winding += crosser.SignedEdgeOrVertexCrossing(edge.v0, edge.v1);

      // Include any edges incident on vertex in the contains vertex query.
      if (edge.IncidentOn(vertex)) {
        if (vertex == edge.v0) {
          query.AddEdge(edge.v1, +1);
        } else {
          query.AddEdge(edge.v0, -1);
        }
      }
    }

    bool duplicates = query.DuplicateEdges();
    int sign = 0;

    // If we have a sign of zero on the vertex, all the edges incident on it
    // were reverse duplicates and we can't use it to test orientation, continue
    // trying to find another vertex.
    if (!duplicates) {
      sign = query.ContainsSign();
      if (sign == 0) {
        continue;
      }
    }

    // The sign bit obtained by crossing edges should be consistent with the
    // sign produced by the S2ContainsVertexQuery.
    if (duplicates || winding != (sign < 0 ? 0 : 1)) {
      *error = S2Error(
          S2Error::POLYGON_INCONSISTENT_LOOP_ORIENTATIONS,
          absl::StrFormat(
              "Shape %d has one or more edges with interior on the right.",
              shape_id));
      return false;
    }
    return true;
  }
  return true;
}

template <typename IndexType>
bool S2ValidQuery<IndexType>::CheckEdge(const S2Shape& shape,
                                        const S2ClippedShape& clipped,
                                        const EdgeAndIdChain& edge,
                                        S2Error* error) {
  const S2IndexCellData& cell = CurrentCell();
  const int dim = shape.dimension();

  // For points, we can check that they're not contained in any other polygons
  // locally within the cell by crossing edges to the cell center.
  if (dim == 0 &&
      PointContained(cell.id(), clipped.shape_id(), edge.v0, error)) {
    return false;
  }

  // Edge is OK
  return true;
}

// Checks that edges of the given shape incident on vertex are ordered such that
// the incident chains do not cross.
//
// We can check this by looking at all the incident edges and making sure that,
// for each incoming edge, as we move counter-clockwise around the vertex, we
// encounter matching pairs of incoming/outgoing edges for each chain.
//
// Returns true if chains do not cross at the vertex, false otherwise.
inline bool CheckVertexCrossings(const S2Point& vertex, const S2Shape& shape,
                                 int shape_id,
                                 const absl::flat_hash_set<int32_t>& edge_ids,
                                 S2Error* error) {
  // Extend S2Shape::Edge to wrap an edge with its chain, id and previous id.
  struct EdgeWithInfo : public S2Shape::Edge {
    EdgeWithInfo(S2Shape::Edge edge, int id, int chain, int prev, int sign)
        : S2Shape::Edge(edge), id(id), chain(chain), prev(prev), sign(sign) {}

    int id;
    int chain;
    int prev;
    int sign;  // +1 for incoming and -1 for outgoing edges.
  };

  // Aggregate edges of the current shape and sort them CCW around the vertex.
  absl::InlinedVector<EdgeWithInfo, 6> edges;
  for (int32_t edge_id : edge_ids) {
    S2Shape::ChainPosition pos = shape.chain_position(edge_id);
    const S2Shape::Edge& edge = shape.edge(edge_id);
    edges.push_back({
        edge,                                       //
        edge_id,                                    //
        pos.chain_id,                               //
        s2shapeutil::PrevEdgeWrap(shape, edge_id),  //
        edge.v0 == vertex ? -1 : +1,                //
    });
  }
  SortEdgesCcw(vertex, edges[0], edges);

  // Mapping from chain_id => count.  We'll scan through and sum the signs on
  // the edges for each chain.  When we reach the incoming edge the sums should
  // be zero for every chain.
  absl::flat_hash_map<int, int> chain_sums;
  chain_sums.reserve(16);

  for (size_t i = 0; i < edges.size(); ++i) {
    const EdgeWithInfo& curr = edges[i];

    // Skip forward to next outgoing edge.
    if (curr.sign > 0) {
      continue;
    }

    // Scan until we find our incoming edge and tally chain counts.
    chain_sums.clear();
    size_t j;
    for (j = 1; j < edges.size(); ++j) {
      const EdgeWithInfo& edge = edges[(i + j) % edges.size()];
      if (curr.chain == edge.chain && curr.prev == edge.id) {
        for (const auto& sum : chain_sums) {
          if (sum.second != 0) {
            *error = S2Error(
                S2Error::OVERLAPPING_GEOMETRY,
                absl::StrFormat(
                    "Shape %d has one or more chains that cross at a vertex",
                    shape_id));
            return false;
          }
        }
        break;
      }

      chain_sums[edge.chain] += edge.sign;
    }

    // If we went all the way around and didn't find an incoming edge, then
    // the geometry must be malformed, and thus isn't valid.
    if (j == edges.size()) {
      *error = S2Error(S2Error::INVALID_VERTEX,
                       "Outgoing edge with no incoming edge");
      return false;
    }
  }

  return true;
}

template <typename IndexType>
bool S2ValidQuery<IndexType>::Finish(S2Error* error) {
  // We've checked edges having interiors on the right, and for crossings at
  // interior points.  The only case left is to check for chains that cross at a
  // vertex.
  for (const auto& item : IncidentEdges()) {
    const S2Shape& shape = *Index().shape(item.first.shape_id);
    if (shape.dimension() == 2) {
      if (!CheckVertexCrossings(item.first.vertex, shape, item.first.shape_id,
                                item.second, error)) {
        return false;
      }
    }
  }

  // If we get to this point we know that polygon edges don't cross any other
  // edges and that edges are properly oriented with the interior on the left.
  //
  // Since edges don't cross, any given chain must be entirely inside or outside
  // any other polygons.  Thus, to determine that polygon interiors are
  // disjoint, we only have to check one vertex of each chain in each shape for
  // containment.
  //
  // We use the OPEN containment model because we check elsewhere if the vertex
  // lands on another vertex.
  S2ContainsPointQueryOptions options;
  options.set_vertex_model(S2VertexModel::OPEN);

  S2ContainsPointQuery<IndexType> query(&Index(), options);
  for (int shape_id = 0; shape_id < Index().num_shape_ids(); ++shape_id) {
    const S2Shape* shape = Index().shape(shape_id);
    if (shape == nullptr || shape->dimension() == 0) {
      continue;
    }

    for (int chain = 0; chain < shape->num_chains(); ++chain) {
      if (shape->chain(chain).length < 1) {
        continue;
      }

      S2Point vertex = shape->chain_edge(chain, 0).v0;
      if (query.Contains(vertex)) {
        *error = S2Error(
            S2Error::OVERLAPPING_GEOMETRY,
            absl::StrFormat(
                "Shape %d has one or more edges contained in another shape.",
                shape_id));
        return false;
      }
    }
  }
  return true;
}

//////////////////   S2LegacyValidQuery Implementation   ////////////////////

template <typename IndexType>
S2LegacyValidQuery<IndexType>::S2LegacyValidQuery() : Base() {
  // Don't allow degenerate or reverse duplicates edges for legacy semantics.
  Base::mutable_options().set_allow_degenerate_edges(false);
  Base::mutable_options().set_allow_reverse_duplicates(false);
}

template <typename IndexType>
bool S2LegacyValidQuery<IndexType>::Start(S2Error* error) {
  if (!Base::Start(error)) {
    return false;
  }

  // We can't mix dimensions under legacy semantics.
  int dim = -1;
  for (const S2Shape* shape : Index()) {
    if (dim < 0) {
      dim = shape->dimension();
    }

    if (dim != shape->dimension()) {
      *error = S2Error(
          S2Error::INVALID_DIMENSION,
          "Mixed dimensional geometry is invalid for legacy semantics.");
      return false;
    }
  }

  return true;
}

template <typename IndexType>
bool S2LegacyValidQuery<IndexType>::CheckShape(Iterator& iter,
                                               const S2Shape& shape,
                                               int shape_id, S2Error* error) {
  // Count the number of empty chains.  Non-empty chains must have at least
  // three vertices.
  if (shape.dimension() == 2) {
    bool has_empty_loops = false;
    for (const S2Shape::Chain& chain : shape.chains()) {
      if (chain.length == 0) {
        has_empty_loops = true;
      } else if (chain.length < 3) {
        *error = S2Error(
            S2Error::LOOP_NOT_ENOUGH_VERTICES,
            absl::StrFormat(
                "Shape %d has a non-empty chain with less than three edges.",
                shape_id));
        return false;
      }
    }

    if (has_empty_loops && shape.num_chains() > 1) {
      *error = S2Error(
          S2Error::POLYGON_EMPTY_LOOP,
          absl::StrFormat("Shape %d has too many empty chains", shape_id));
      return false;
    }
  }

  if (!Base::CheckShape(iter, shape, shape_id, error)) {
    return false;
  }

  return true;
}

template <typename IndexType>
bool S2LegacyValidQuery<IndexType>::StartCell(S2Error* error) {
  // Check for duplicate vertices within a chain.
  const S2IndexCellData& cell = CurrentCell();
  for (const S2ClippedShape& clipped : cell.clipped_shapes()) {
    absl::Span<const EdgeAndIdChain> edges =
        cell.shape_edges(clipped.shape_id());

    for (size_t i = 0; i < edges.size(); ++i) {
      for (size_t j = i + 1; j < edges.size(); ++j) {
        if (edges[j].chain != edges[i].chain) {
          continue;
        }

        if (edges[j].v0 == edges[i].v0) {
          *error = S2Error(
              S2Error::DUPLICATE_VERTICES,
              absl::StrFormat("Chain %d of shape %d has duplicate vertices",
                              edges[i].chain, clipped.shape_id()));
          return false;
        }
      }
    }
  }

  return Base::StartCell(error);
}

template <typename IndexType>
bool S2LegacyValidQuery<IndexType>::CheckEdge(const S2Shape& shape,
                                              const S2ClippedShape& clipped,
                                              const EdgeAndIdChain& edge,
                                              S2Error* error) {
  if (!Base::CheckEdge(shape, clipped, edge, error)) {
    return false;
  }
  return true;
}

#endif  // S2_S2VALIDATION_QUERY_H_
