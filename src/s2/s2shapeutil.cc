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

#include "s2/s2shapeutil.h"

#include "s2/base/stringprintf.h"
#include "s2/third_party/absl/types/span.h"
#include "s2/s2crossingedgequery.h"
#include "s2/s2error.h"
#include "s2/s2loop.h"
#include "s2/s2paddedcell.h"
#include "s2/s2pointutil.h"
#include "s2/s2predicates.h"
#include "s2/s2shapeindex.h"

using absl::Span;
using std::pair;
using std::unique_ptr;
using std::vector;
using ChainPosition = S2Shape::ChainPosition;

namespace s2shapeutil {

LaxLoop::LaxLoop(vector<S2Point> const& vertices) {
  Init(vertices);
}

LaxLoop::LaxLoop(S2Loop const& loop) {
  Init(loop);
}

void LaxLoop::Init(vector<S2Point> const& vertices) {
  num_vertices_ = vertices.size();
  vertices_.reset(new S2Point[num_vertices_]);
  std::copy(vertices.begin(), vertices.end(), vertices_.get());
}

void LaxLoop::Init(S2Loop const& loop) {
  DCHECK(!loop.is_full()) << "Full loops not currently supported";
  if (loop.is_empty()) {
    num_vertices_ = 0;
    vertices_ = nullptr;
  } else {
    num_vertices_ = loop.num_vertices();
    vertices_.reset(new S2Point[num_vertices_]);
    std::copy(&loop.vertex(0), &loop.vertex(0) + num_vertices_,
              vertices_.get());
  }
}

S2Shape::Edge LaxLoop::edge(int e0) const {
  DCHECK_LT(e0, num_edges());
  int e1 = e0 + 1;
  if (e1 == num_vertices()) e1 = 0;
  return Edge(vertices_[e0], vertices_[e1]);
}

S2Shape::Edge LaxLoop::chain_edge(int i, int j) const {
  DCHECK_EQ(i, 0);
  DCHECK_LT(j, num_edges());
  int k = (j + 1 == num_vertices()) ? 0 : j + 1;
  return Edge(vertices_[j], vertices_[k]);
}

bool LaxLoop::contains_origin() const {
  // IsOriginOnLeft interprets a loop with no vertices as "full".
  if (num_vertices() == 0) return false;
  return IsOriginOnLeft(*this);
}

ClosedLaxPolyline::ClosedLaxPolyline(std::vector<S2Point> const& vertices) {
  Init(vertices);
}

ClosedLaxPolyline::ClosedLaxPolyline(S2Loop const& loop) {
  Init(loop);
}

LaxPolygon::LaxPolygon(vector<LaxPolygon::Loop> const& loops) {
  Init(loops);
}

LaxPolygon::LaxPolygon(S2Polygon const& polygon) {
  Init(polygon);
}

void LaxPolygon::Init(vector<LaxPolygon::Loop> const& loops) {
  vector<Span<S2Point const>> spans;
  for (LaxPolygon::Loop const& loop : loops) {
    spans.emplace_back(loop);
  }
  Init(spans);
}

void LaxPolygon::Init(S2Polygon const& polygon) {
  vector<Span<S2Point const>> spans;
  for (int i = 0; i < polygon.num_loops(); ++i) {
    S2Loop const* loop = polygon.loop(i);
    if (loop->is_full()) {
      spans.emplace_back();  // Empty span.
    } else {
      spans.emplace_back(&loop->vertex(0), loop->num_vertices());
    }
  }
  Init(spans);
}

void LaxPolygon::Init(vector<Span<S2Point const>> const& loops) {
  num_loops_ = loops.size();
  if (num_loops_ == 0) {
    num_vertices_ = 0;
    vertices_ = nullptr;
  } else if (num_loops_ == 1) {
    num_vertices_ = loops[0].size();
    vertices_.reset(new S2Point[num_vertices_]);
    std::copy(loops[0].begin(), loops[0].end(), vertices_.get());
  } else {
    cumulative_vertices_ = new int32[num_loops_ + 1];
    int32 num_vertices = 0;
    for (int i = 0; i < num_loops_; ++i) {
      cumulative_vertices_[i] = num_vertices;
      num_vertices += loops[i].size();
    }
    cumulative_vertices_[num_loops_] = num_vertices;
    vertices_.reset(new S2Point[num_vertices]);
    for (int i = 0; i < num_loops_; ++i) {
      std::copy(loops[i].begin(), loops[i].end(),
                vertices_.get() + cumulative_vertices_[i]);
    }
  }
}

LaxPolygon::~LaxPolygon() {
  if (num_loops() > 1) {
    delete[] cumulative_vertices_;
  }
}

int LaxPolygon::num_vertices() const {
  if (num_loops() <= 1) {
    return num_vertices_;
  } else {
    return cumulative_vertices_[num_loops()];
  }
}

int LaxPolygon::num_loop_vertices(int i) const {
  DCHECK_LT(i, num_loops());
  if (num_loops() == 1) {
    return num_vertices_;
  } else {
    return cumulative_vertices_[i + 1] - cumulative_vertices_[i];
  }
}

S2Point const& LaxPolygon::loop_vertex(int i, int j) const {
  DCHECK_LT(i, num_loops());
  DCHECK_LT(j, num_loop_vertices(i));
  if (num_loops() == 1) {
    return vertices_[j];
  } else {
    return vertices_[cumulative_vertices_[i] + j];
  }
}

S2Shape::Edge LaxPolygon::edge(int e0) const {
  DCHECK_LT(e0, num_edges());
  int e1 = e0 + 1;
  if (num_loops() == 1) {
    if (e1 == num_vertices_) { e1 = 0; }
  } else {
    // Find the index of the first vertex of the loop following this one.
    int const kMaxLinearSearchLoops = 12;  // From benchmarks.
    int* next = cumulative_vertices_ + 1;
    if (num_loops() <= kMaxLinearSearchLoops) {
      while (*next <= e0) ++next;
    } else {
      next = std::upper_bound(next, next + num_loops(), e0);
    }
    // Wrap around to the first vertex of the loop if necessary.
    if (e1 == *next) { e1 = next[-1]; }
  }
  return Edge(vertices_[e0], vertices_[e1]);
}

bool LaxPolygon::contains_origin() const {
  return IsOriginOnLeft(*this);
}

S2Shape::Chain LaxPolygon::chain(int i) const {
  DCHECK_LT(i, num_loops());
  if (num_loops() == 1) {
    return Chain(0, num_vertices_);
  } else {
    int start = cumulative_vertices_[i];
    return Chain(start, cumulative_vertices_[i + 1] - start);
  }
}

S2Shape::Edge LaxPolygon::chain_edge(int i, int j) const {
  DCHECK_LT(i, num_loops());
  DCHECK_LT(j, num_loop_vertices(i));
  int n = num_loop_vertices(i);
  int k = (j + 1 == n) ? 0 : j + 1;
  if (num_loops() == 1) {
    return Edge(vertices_[j], vertices_[k]);
  } else {
    int base = cumulative_vertices_[i];
    return Edge(vertices_[base + j], vertices_[base + k]);
  }
}

S2Shape::ChainPosition LaxPolygon::chain_position(int e) const {
  DCHECK_LT(e, num_edges());
  int const kMaxLinearSearchLoops = 12;  // From benchmarks.
  if (num_loops() == 1) {
    return ChainPosition(0, e);
  } else {
    // Find the index of the first vertex of the loop following this one.
    int* next = cumulative_vertices_ + 1;
    if (num_loops() <= kMaxLinearSearchLoops) {
      while (*next <= e) ++next;
    } else {
      next = std::upper_bound(next, next + num_loops(), e);
    }
    return ChainPosition(next - (cumulative_vertices_ + 1), e - next[-1]);
  }
}

LaxPolyline::LaxPolyline(vector<S2Point> const& vertices) {
  Init(vertices);
}

LaxPolyline::LaxPolyline(S2Polyline const& polyline) {
  Init(polyline);
}

void LaxPolyline::Init(vector<S2Point> const& vertices) {
  num_vertices_ = vertices.size();
  LOG_IF(WARNING, num_vertices_ == 1)
      << "s2shapeutil::LaxPolyline with one vertex has no edges";
  vertices_.reset(new S2Point[num_vertices_]);
  std::copy(vertices.begin(), vertices.end(), vertices_.get());
}

void LaxPolyline::Init(S2Polyline const& polyline) {
  num_vertices_ = polyline.num_vertices();
  LOG_IF(WARNING, num_vertices_ == 1)
      << "s2shapeutil::LaxPolyline with one vertex has no edges";
  vertices_.reset(new S2Point[num_vertices_]);
  std::copy(&polyline.vertex(0), &polyline.vertex(0) + num_vertices_,
            vertices_.get());
}

S2Shape::Edge LaxPolyline::edge(int e) const {
  DCHECK_LT(e, num_edges());
  return Edge(vertices_[e], vertices_[e + 1]);
}

int LaxPolyline::num_chains() const {
  return std::min(1, LaxPolyline::num_edges());  // Avoid virtual call.
}

S2Shape::Chain LaxPolyline::chain(int i) const {
  return Chain(0, LaxPolyline::num_edges());  // Avoid virtual call.
}

S2Shape::Edge LaxPolyline::chain_edge(int i, int j) const {
  DCHECK_EQ(i, 0);
  DCHECK_LT(j, num_edges());
  return Edge(vertices_[j], vertices_[j + 1]);
}

S2Shape::ChainPosition LaxPolyline::chain_position(int e) const {
  return ChainPosition(0, e);
}

VertexIdLaxLoop::VertexIdLaxLoop(std::vector<int32> const& vertex_ids,
                                 S2Point const* vertex_array) {
  Init(vertex_ids, vertex_array);
}

void VertexIdLaxLoop::Init(std::vector<int32> const& vertex_ids,
                           S2Point const* vertex_array) {
  num_vertices_ = vertex_ids.size();
  vertex_ids_.reset(new int32[num_vertices_]);
  std::copy(vertex_ids.begin(), vertex_ids.end(), vertex_ids_.get());
  vertex_array_ = vertex_array;
}

S2Shape::Edge VertexIdLaxLoop::edge(int e0) const {
  DCHECK_LT(e0, num_edges());
  int e1 = e0 + 1;
  if (e1 == num_vertices()) e1 = 0;
  return Edge(vertex(e0), vertex(e1));
}

S2Shape::Edge VertexIdLaxLoop::chain_edge(int i, int j) const {
  DCHECK_EQ(i, 0);
  DCHECK_LT(j, num_edges());
  int k = (j + 1 == num_vertices()) ? 0 : j + 1;
  return Edge(vertex(j), vertex(k));
}

bool VertexIdLaxLoop::contains_origin() const {
  // IsOriginOnLeft interprets a loop with no vertices as "full".
  if (num_vertices() == 0) return false;
  return IsOriginOnLeft(*this);
}

// This is a helper function for IsOriginOnLeft(), defined below.
//
// If the given vertex "vtest" is unbalanced (see definition below), sets
// "result" to indicate whether "shape" contains S2::Origin() and returns
// true.  Otherwise returns false.
static bool IsOriginOnLeftAtVertex(S2Shape const& shape,
                                   S2Point const& vtest, bool* result) {
  // Let P be an unbalanced vertex.  Vertex P is defined to be inside the
  // region if the region contains a particular direction vector starting from
  // P, namely the direction S2::Ortho(P).  Since the interior is defined as
  // the region to the left of all loops, this means we need to find the
  // unmatched edge incident to P that is immediately clockwise from
  // S2::Ortho(P).  P is contained by the region if and only if this edge is
  // outgoing.
  //
  // To convert this into a contains_origin() value, we count the number of
  // edges crossed between P and S2::Origin(), and invert the result for
  // every crossing.
  S2CopyingEdgeCrosser crosser(S2::Origin(), vtest);
  bool crossing_parity = false;
  util::btree::btree_map<S2Point, int> edge_map;
  int n = shape.num_edges();
  for (int e = 0; e < n; ++e) {
    auto edge = shape.edge(e);
    if (edge.v0 == edge.v1) continue;

    // Check whether this edge crosses the edge between P and S2::Origin().
    crossing_parity ^= crosser.EdgeOrVertexCrossing(edge.v0, edge.v1);

    // Keep track of (outgoing edges) - (incoming edges) for each vertex that
    // is adjacent to "vtest".
    if (edge.v0 == vtest) ++edge_map[edge.v1];
    if (edge.v1 == vtest) --edge_map[edge.v0];
  }
  // Find the unmatched edge that is immediately clockwise from S2::Ortho(P).
  S2Point reference_dir = S2::Ortho(vtest);
  pair<S2Point, int> best(reference_dir, 0);
  for (auto const& e : edge_map) {
    if (e.second == 0) continue;  // This is a "matched" edge.
    if (s2pred::OrderedCCW(reference_dir, best.first, e.first, vtest)) {
      best = e;
    }
  }
  if (best.second == 0) {
    return false;  // There are no unmatched edges incident to this vertex.
  }
  // Point P is contained by the shape if the edge immediately clockwise from
  // S2::Ortho(P) is an outgoing edge.  We then invert this result if the
  // number of crossings between P and S2::Origin() was odd.
  *result = crossing_parity != (best.second > 0);
  return true;
}

// See documentation in header file.
bool IsOriginOnLeft(S2Shape const& shape) {
  if (shape.num_edges() == 0) {
    // The shape is defined to be "full" if it contains an empty loop.
    return shape.num_chains() > 0;
  }
  // Define a "matched" edge as one that can be paired with a corresponding
  // reversed edge.  Define a vertex as "balanced" if all of its edges are
  // matched. In order to determine containment, we must find an unbalanced
  // vertex.  Often every vertex is unbalanced, so we start by trying an
  // arbitrary vertex.
  auto edge = shape.edge(0);
  bool result = false;
  if (IsOriginOnLeftAtVertex(shape, edge.v0, &result)) {
    return result;
  }
  // That didn't work, so now we do some extra work to find an unbalanced
  // vertex (if any).  Essentially we gather a list of edges and a list of
  // reversed edges, and then sort them.  The first edge that appears in one
  // list but not the other is guaranteed to be unmatched.
  int n = shape.num_edges();
  vector<S2Shape::Edge> edges, rev_edges;
  edges.reserve(n);
  rev_edges.reserve(n);
  for (int i = 0; i < n; ++i) {
    auto edge = shape.edge(i);
    edges.push_back(edge);
    rev_edges.push_back(S2Shape::Edge(edge.v1, edge.v0));
  }
  std::sort(edges.begin(), edges.end());
  std::sort(rev_edges.begin(), rev_edges.end());
  for (int i = 0; i < n; ++i) {
    if (edges[i] < rev_edges[i]) {  // edges[i] is unmatched
      CHECK(IsOriginOnLeftAtVertex(shape, edges[i].v0, &result));
      return result;
    }
    if (rev_edges[i] < edges[i]) {  // rev_edges[i] is unmatched
      CHECK(IsOriginOnLeftAtVertex(shape, rev_edges[i].v0, &result));
      return result;
    }
  }
  // All vertices are balanced, so this polygon is either empty or full.  By
  // convention it is defined to be "full" if it contains any empty loop.
  for (int i = 0; i < shape.num_chains(); ++i) {
    if (shape.chain(i).length == 0) return true;
  }
  return false;
}

std::ostream& operator<<(std::ostream& os, ShapeEdgeId id) {
  return os << id.shape_id << ":" << id.edge_id;
}

// Ensure that we don't usually need to allocate memory when collecting the
// edges in an S2ShapeIndex cell (which by default have about 10 edges).
using ShapeEdgeVector = absl::InlinedVector<ShapeEdge, 16>;

// Returns a vector containing all edges in the given S2ShapeIndexCell.
// (The result is returned as an output parameter so that the same storage can
// be reused, rather than allocating a new temporary vector each time.)
static void GetShapeEdges(S2ShapeIndex const& index,
                          S2ShapeIndexCell const& cell,
                          ShapeEdgeVector* shape_edges) {
  shape_edges->clear();
  for (int s = 0; s < cell.num_clipped(); ++s) {
    S2ClippedShape const& clipped = cell.clipped(s);
    S2Shape const& shape = *index.shape(clipped.shape_id());
    int num_clipped = clipped.num_edges();
    for (int i = 0; i < num_clipped; ++i) {
      shape_edges->push_back(ShapeEdge(shape, clipped.edge(i)));
    }
  }
}

// Given a vector of edges within an S2ShapeIndexCell, visit all pairs of
// crossing edges (of the given CrossingType).
static bool VisitCrossings(ShapeEdgeVector const& shape_edges,
                           CrossingType type, EdgePairVisitor const& visitor) {
  int min_crossing_sign = (type == CrossingType::INTERIOR) ? 1 : 0;
  int num_edges = shape_edges.size();
  for (int i = 0; i + 1 < num_edges; ++i) {
    ShapeEdge const& a = shape_edges[i];
    int j = i + 1;
    // A common situation is that an edge AB is followed by an edge BC.  We
    // only need to visit such crossings if CrossingType::ALL is specified
    // (even if AB and BC belong to different edge chains).
    if (type != CrossingType::ALL && a.v1() == shape_edges[j].v0()) {
      if (++j >= num_edges) break;
    }
    S2EdgeCrosser crosser(&a.v0(), &a.v1());
    for (; j < num_edges; ++j) {
      ShapeEdge const& b = shape_edges[j];
      int sign = crosser.CrossingSign(&b.v0(), &b.v1());
      if (sign >= min_crossing_sign) {
        if (!visitor(a, b, sign == 1)) return false;
      }
    }
  }
  return true;
}

bool VisitCrossings(S2ShapeIndex const& index, CrossingType type,
                    EdgePairVisitor const& visitor) {
  // TODO(ericv): Use brute force if the total number of edges is small enough
  // (using a larger threshold if the S2ShapeIndex is not constructed yet).
  ShapeEdgeVector shape_edges;
  for (S2ShapeIndex::Iterator it(&index); !it.Done(); it.Next()) {
    GetShapeEdges(index, it.cell(), &shape_edges);
    if (!VisitCrossings(shape_edges, type, visitor)) return false;
  }
  return true;
}

//////////////////////////////////////////////////////////////////////

// TODO(user): Change arg to pointer.
RangeIterator::RangeIterator(S2ShapeIndex const& index)
    : it_(&index), end_(S2CellId::End(0)) {
  Refresh();
}

void RangeIterator::Next() {
  it_.Next();
  Refresh();
}

void RangeIterator::SeekTo(RangeIterator const& target) {
  it_.Seek(target.range_min());
  // If the current cell does not overlap "target", it is possible that the
  // previous cell is the one we are looking for.  This can only happen when
  // the previous cell contains "target" but has a smaller S2CellId.
  if (it_.Done() || it_.id().range_min() > target.range_max()) {
    it_.Prev();
    if (it_.id().range_max() < target.id()) it_.Next();
  }
  Refresh();
}

void RangeIterator::SeekBeyond(RangeIterator const& target) {
  it_.Seek(target.range_max().next());
  if (!it_.Done() && it_.id().range_min() <= target.range_max()) {
    it_.Next();
  }
  Refresh();
}

// This method is inline, but is only called by non-inline methods defined in
// this file.  Putting the definition here enforces this requirement.
inline void RangeIterator::Refresh() {
  if (it_.Done()) {
    id_ = end_;
    cell_ = nullptr;
  } else {
    id_ = it_.id();
    cell_ = &it_.cell();
  }
  range_min_ = id_.range_min();
  range_max_ = id_.range_max();
}

// IndexCrosser is a helper class for determining whether two loops cross.
// It is instantiated twice for each pair of loops to be tested, once for the
// pair (A,B) and once for the pair (B,A), in order to be able to test edge
// crossings in the most efficient order.
class IndexCrosser {
 public:
  // If "swapped" is true, the loops A and B have been swapped.  This affects
  // how arguments are passed to the given loop relation, since for example
  // A.Contains(B) is not the same as B.Contains(A).
  IndexCrosser(S2ShapeIndex const& a_index, S2ShapeIndex const& b_index,
               CrossingType type, EdgePairVisitor const& visitor, bool swapped)
      : a_index_(a_index), b_index_(b_index), visitor_(visitor),
        min_crossing_sign_(type == CrossingType::INTERIOR ? 1 : 0),
        swapped_(swapped), b_query_(&b_index_) {
  }

  // Given two iterators positioned such that ai->id().Contains(bi->id()),
  // visit all crossings between edges of A and B that intersect a->id().
  // Terminates early and returns false if visitor_ returns false.
  // Advances both iterators past ai->id().
  bool VisitCrossings(RangeIterator* ai, RangeIterator* bi);

  // Given two index cells, visit all crossings between edges of those cells.
  // Terminates early and returns false if visitor_ returns false.
  bool VisitCellCellCrossings(S2ShapeIndexCell const& a_cell,
                              S2ShapeIndexCell const& b_cell);

 private:
  bool VisitEdgePair(ShapeEdge const& a, ShapeEdge const& b, bool is_interior);

  // Visit all crossings of the current edge with all edges of the given index
  // cell of B.  Terminates early and returns false if visitor_ returns false.
  bool VisitEdgeCellCrossings(ShapeEdge const& a,
                              S2ShapeIndexCell const& b_cell);

  // Visit all crossings of any edge in "a_cell" with any index cell of B that
  // is a descendant of "b_id".  Terminates early and returns false if
  // visitor_ returns false.
  bool VisitSubcellCrossings(S2ShapeIndexCell const& a_cell, S2CellId b_id);

  S2ShapeIndex const& a_index_;
  S2ShapeIndex const& b_index_;
  EdgePairVisitor const& visitor_;
  int const min_crossing_sign_;
  bool const swapped_;

  // Temporary data declared here to avoid repeated memory allocations.
  S2CrossingEdgeQuery b_query_;
  S2EdgeCrosser crosser_;
  vector<S2ShapeIndexCell const*> b_cells_;
  ShapeEdgeVector a_shape_edges_;
  ShapeEdgeVector b_shape_edges_;
};

inline bool IndexCrosser::VisitEdgePair(ShapeEdge const& a, ShapeEdge const& b,
                                        bool is_interior) {
  if (swapped_) {
    return visitor_(b, a, is_interior);
  } else {
    return visitor_(a, b, is_interior);
  }
}

bool IndexCrosser::VisitEdgeCellCrossings(ShapeEdge const& a,
                                          S2ShapeIndexCell const& b_cell) {
  // Test the current edge of A against all edges of "b_cell".
  GetShapeEdges(b_index_, b_cell, &b_shape_edges_);
  for (int j = 0; j < b_shape_edges_.size(); ++j) {
    ShapeEdge const& b = b_shape_edges_[j];
    int sign = crosser_.CrossingSign(&b.v0(), &b.v1());
    if (sign >= min_crossing_sign_) {
      if (!VisitEdgePair(a, b, sign == 1)) return false;
    }
  }
  return true;
}

bool IndexCrosser::VisitSubcellCrossings(S2ShapeIndexCell const& a_cell,
                                         S2CellId b_id) {
  // Test all edges of "a_cell" against the edges contained in B index cells
  // that are descendants of "b_id".
  GetShapeEdges(a_index_, a_cell, &a_shape_edges_);
  S2PaddedCell b_root(b_id, 0);
  for (int i = 0; i < a_shape_edges_.size(); ++i) {
    // Use an S2CrossingEdgeQuery starting at "b_root" to find the index cells
    // of B that might contain crossing edges.
    ShapeEdge const& a = a_shape_edges_[i];
    if (b_query_.GetCells(a.v0(), a.v1(), b_root, &b_cells_)) {
      crosser_.Init(&a.v0(), &a.v1());
      for (int c = 0; c < b_cells_.size(); ++c) {
        if (!VisitEdgeCellCrossings(a, *b_cells_[c])) return false;
      }
    }
  }
  return true;
}

bool IndexCrosser::VisitCellCellCrossings(S2ShapeIndexCell const& a_cell,
                                          S2ShapeIndexCell const& b_cell) {
  // Test all edges of "a_shape_edges_" against all edges of "b_cell".
  // TODO(ericv): Refactor this so that we don't need to gather the edges of
  // "a_cell" every time when testing against multiple B cells.
  // TODO(ericv): Process the cell with fewer edges in the outer loop.
  GetShapeEdges(a_index_, a_cell, &a_shape_edges_);
  GetShapeEdges(b_index_, b_cell, &b_shape_edges_);
  for (int i = 0; i < a_shape_edges_.size(); ++i) {
    ShapeEdge const& a = a_shape_edges_[i];
    crosser_.Init(&a.v0(), &a.v1());
    for (int j = 0; j < b_shape_edges_.size(); ++j) {
      ShapeEdge const& b = b_shape_edges_[j];
      int sign = crosser_.CrossingSign(&b.v0(), &b.v1());
      if (sign >= min_crossing_sign_) {
        if (!VisitEdgePair(a, b, sign == 1)) return false;
      }
    }
  }
  return true;
}

bool IndexCrosser::VisitCrossings(RangeIterator* ai, RangeIterator* bi) {
  DCHECK(ai->id().contains(bi->id()));
  if (ai->cell().num_edges() == 0) {
    // Skip over the cells of B using binary search.
    bi->SeekBeyond(*ai);
  } else {
    // If ai->id() intersects many edges of B, then it is faster to use
    // S2CrossingEdgeQuery to narrow down the candidates.  But if it
    // intersects only a few edges, it is faster to check all the crossings
    // directly.  We handle this by advancing "bi" and keeping track of how
    // many edges we would need to test.
    static int const kEdgeQueryMinEdges = 20;  // TODO: Tune using benchmarks.
    int b_edges = 0;
    b_cells_.clear();
    do {
      int cell_edges = bi->cell().num_edges();
      if (cell_edges > 0) {
        b_edges += cell_edges;
        if (b_edges >= kEdgeQueryMinEdges) {
          // There are too many edges, so use an S2CrossingEdgeQuery.
          if (!VisitSubcellCrossings(ai->cell(), ai->id())) return false;
          bi->SeekBeyond(*ai);
          return true;
        }
        b_cells_.push_back(&bi->cell());
      }
      bi->Next();
    } while (bi->id() <= ai->range_max());

    // Test all the edge crossings directly.
    for (int c = 0; c < b_cells_.size(); ++c) {
      if (!VisitCellCellCrossings(ai->cell(), *b_cells_[c])) {
        return false;
      }
    }
  }
  ai->Next();
  return true;
}

//////////////////////////////////////////////////////////////////////

bool VisitCrossings(S2ShapeIndex const& a_index, S2ShapeIndex const& b_index,
                    CrossingType type, EdgePairVisitor const& visitor) {
  DCHECK(type != CrossingType::NON_ADJACENT);
  // We look for S2CellId ranges where the indexes of A and B overlap, and
  // then test those edges for crossings.
  //
  // TODO(ericv): Use brute force if the total number of edges is small enough
  // (using a larger threshold if the S2ShapeIndex is not constructed yet).
  //
  // TODO(ericv): Consider using one crosser and passing "a_index", "b_index",
  // and "swapped" as parameters.  This would allow cell-cell intersections to
  // be optimized by processing the index with fewer edges in the outer loop.
  RangeIterator ai(a_index), bi(b_index);
  IndexCrosser ab(a_index, b_index, type, visitor, false); // Tests A against B
  IndexCrosser ba(b_index, a_index, type, visitor, true);  // Tests B against A
  while (!ai.Done() || !bi.Done()) {
    if (ai.range_max() < bi.range_min()) {
      // The A and B cells don't overlap, and A precedes B.
      ai.SeekTo(bi);
    } else if (bi.range_max() < ai.range_min()) {
      // The A and B cells don't overlap, and B precedes A.
      bi.SeekTo(ai);
    } else {
      // One cell contains the other.  Determine which cell is larger.
      int64 ab_relation = ai.id().lsb() - bi.id().lsb();
      if (ab_relation > 0) {
        // A's index cell is larger.
        if (!ab.VisitCrossings(&ai, &bi)) return false;
      } else if (ab_relation < 0) {
        // B's index cell is larger.
        if (!ba.VisitCrossings(&bi, &ai)) return false;
      } else {
        // The A and B cells are the same.
        if (ai.cell().num_edges() > 0 && bi.cell().num_edges() > 0) {
          if (!ab.VisitCellCellCrossings(ai.cell(), bi.cell())) return false;
        }
        ai.Next();
        bi.Next();
      }
    }
  }
  return true;
}

void ResolveComponents(vector<vector<S2Shape*>> const& components,
                       vector<vector<S2Shape*>>* faces) {
  faces->clear();
  if (components.empty()) return;

  // Since the loops do not overlap, any point on the sphere determines a loop
  // nesting hierarchy.  We choose to build a hierarchy around S2::Origin().
  // The hierarchy then determines the face structure.  Here are the details:
  //
  // 1. Build an S2ShapeIndex of all loops that do not contain S2::Origin().
  //    This leaves at most one unindexed loop per connected component
  //    (the "outer loop").
  // 2. For each component, choose a representative vertex and determine
  //    which indexed loops contain it.  The "depth" of this component is
  //    defined as the number of such loops.
  // 3. Assign the outer loop of each component to the containing loop whose
  //    depth is one less.  This generates a set of multi-loop faces.
  // 4. The output loops of all components at depth 0 become a single face.

  S2ShapeIndex index;
  // A map from shape.id() to the corresponding component number.
  vector<int> component_ids;
  vector<S2Shape*> outer_loops;
  for (int i = 0; i < components.size(); ++i) {
    auto const& component = components[i];
    for (S2Shape* loop : component) {
      if (component.size() > 1 && !loop->contains_origin()) {
        // Ownership is transferred back at the end of this function.
        index.Add(std::unique_ptr<S2Shape>(loop));
        component_ids.push_back(i);
      } else {
        outer_loops.push_back(loop);
      }
    }
    // Check that there is exactly one outer loop in each component.
    DCHECK_EQ(i + 1, outer_loops.size()) << "Component is not a subdivision";
  }
  // Find the loops containing each component.
  vector<vector<S2Shape*>> ancestors(components.size());
  for (int i = 0; i < outer_loops.size(); ++i) {
    auto loop = outer_loops[i];
    DCHECK_GT(loop->num_edges(), 0);
    index.GetContainingShapes(loop->edge(0).v0, &ancestors[i]);
  }
  // Assign each outer loop to the component whose depth is one less.
  // Components at depth 0 become a single face.
  util::btree::btree_map<S2Shape*, vector<S2Shape*>> children;
  for (int i = 0; i < outer_loops.size(); ++i) {
    S2Shape* ancestor = nullptr;
    int depth = ancestors[i].size();
    if (depth > 0) {
      for (auto candidate : ancestors[i]) {
        if (ancestors[component_ids[candidate->id()]].size() == depth - 1) {
          DCHECK(ancestor == nullptr);
          ancestor = candidate;
        }
      }
      DCHECK(ancestor != nullptr);
    }
    children[ancestor].push_back(outer_loops[i]);
  }
  // There is one face per loop that is not an outer loop, plus one for the
  // outer loops of components at depth 0.
  faces->resize(index.num_shape_ids() + 1);
  for (int i = 0; i < index.num_shape_ids(); ++i) {
    auto face = &(*faces)[i];
    auto loop = index.shape(i);
    auto itr = children.find(loop);
    if (itr != children.end()) {
      *face = itr->second;
    }
    face->push_back(loop);
  }
  faces->back() = children[nullptr];

  // Explicitly release the shapes from the index so they are not deleted.
  for (auto& ptr : index.ReleaseAll()) ptr.release();
}

// Helper function that formats a loop error message.  If the loop belongs to
// a multi-loop polygon, adds a prefix indicating which loop is affected.
static void InitLoopError(S2Error::Code code, char const* format,
                          ChainPosition ap, ChainPosition bp,
                          bool is_polygon, S2Error* error) {
  error->Init(code, format, ap.offset, bp.offset);
  if (is_polygon) {
    error->Init(code, "Loop %d: %s", ap.chain_id, error->text().c_str());
  }
}

// Given two loop edges that cross (including at a shared vertex), return true
// if there is a crossing error and set "error" to a human-readable message.
static bool FindCrossingError(S2Shape const& shape,
                              ShapeEdge const& a, ShapeEdge const& b,
                              bool is_interior, S2Error* error) {
  bool is_polygon = shape.num_chains() > 1;
  S2Shape::ChainPosition ap = shape.chain_position(a.id().edge_id);
  S2Shape::ChainPosition bp = shape.chain_position(b.id().edge_id);
  if (is_interior) {
    if (ap.chain_id != bp.chain_id) {
      error->Init(S2Error::POLYGON_LOOPS_CROSS,
                  "Loop %d edge %d crosses loop %d edge %d",
                  ap.chain_id, ap.offset, bp.chain_id, bp.offset);
    } else {
      InitLoopError(S2Error::LOOP_SELF_INTERSECTION,
                    "Edge %d crosses edge %d", ap, bp, is_polygon, error);
    }
    return true;
  }
  // Loops are not allowed to have duplicate vertices, and separate loops
  // are not allowed to share edges or cross at vertices.  We only need to
  // check a given vertex once, so we also require that the two edges have
  // the same end vertex.
  if (a.v1() != b.v1()) return false;
  if (ap.chain_id == bp.chain_id) {
    InitLoopError(S2Error::DUPLICATE_VERTICES,
                  "Edge %d has duplicate vertex with edge %d",
                  ap, bp, is_polygon, error);
    return true;
  }
  int a_len = shape.chain(ap.chain_id).length;
  int b_len = shape.chain(bp.chain_id).length;
  int a_next = (ap.offset + 1 == a_len) ? 0 : ap.offset + 1;
  int b_next = (bp.offset + 1 == b_len) ? 0 : bp.offset + 1;
  S2Point a2 = shape.chain_edge(ap.chain_id, a_next).v1;
  S2Point b2 = shape.chain_edge(bp.chain_id, b_next).v1;
  if (a.v0() == b.v0() || a.v0() == b2) {
    // The second edge index is sometimes off by one, hence "near".
    error->Init(S2Error::POLYGON_LOOPS_SHARE_EDGE,
                "Loop %d edge %d has duplicate near loop %d edge %d",
                ap.chain_id, ap.offset, bp.chain_id, bp.offset);
    return true;
  }
  // Since S2ShapeIndex loops are oriented such that the polygon interior is
  // always on the left, we need to handle the case where one wedge contains
  // the complement of the other wedge.  This is not specifically detected by
  // GetWedgeRelation, so there are two cases to check for.
  //
  // Note that we don't need to maintain any state regarding loop crossings
  // because duplicate edges are detected and rejected above.
  if (S2::GetWedgeRelation(a.v0(), a.v1(), a2, b.v0(), b2) ==
      S2::WEDGE_PROPERLY_OVERLAPS &&
      S2::GetWedgeRelation(a.v0(), a.v1(), a2, b2, b.v0()) ==
      S2::WEDGE_PROPERLY_OVERLAPS) {
    error->Init(S2Error::POLYGON_LOOPS_CROSS,
                "Loop %d edge %d crosses loop %d edge %d",
                ap.chain_id, ap.offset, bp.chain_id, bp.offset);
    return true;
  }
  return false;
}

bool FindAnyCrossing(S2ShapeIndex const& index, S2Error* error) {
  if (index.num_shape_ids() == 0) return false;
  DCHECK_EQ(1, index.num_shape_ids());
  S2Shape const& shape = *index.shape(0);
  return !VisitCrossings(
      index, CrossingType::NON_ADJACENT,
      [&](ShapeEdge const& a, ShapeEdge const& b, bool is_interior) {
        return !FindCrossingError(shape, a, b, is_interior, error);
      });
}

}  // namespace s2shapeutil
