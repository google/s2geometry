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
#include "s2/s2error.h"
#include "s2/s2loop.h"
#include "s2/s2predicates.h"
#include "s2/s2shapeindex.h"

using std::pair;
using std::unique_ptr;
using std::vector;

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

void LaxLoop::GetEdge(int e, S2Point const** a, S2Point const** b) const {
  DCHECK_LT(e, num_edges());
  *a = &vertices_[e];
  if (++e == num_vertices()) e = 0;
  *b = &vertices_[e];
}

bool LaxLoop::contains_origin() const {
  return IsOriginOnLeft(*this);
}

ClosedLaxPolyline::ClosedLaxPolyline(std::vector<S2Point> const& vertices) {
  Init(vertices);
}

ClosedLaxPolyline::ClosedLaxPolyline(S2Loop const& loop) {
  Init(loop);
}

// VertexArray points to an existing array of vertices.  It allows code
// sharing between the S2Polygon and vector<vector<S2Point>> constructors.
class LaxPolygon::VertexArray {
 public:
  VertexArray(S2Point const* begin, size_t size)
      : begin_(begin), size_(size) {
  }
  S2Point const* begin() const { return begin_; }
  S2Point const* end() const { return begin_ + size_; }
  size_t size() const { return size_; }
 private:
  S2Point const* begin_;
  size_t size_;
};

LaxPolygon::LaxPolygon(vector<LaxPolygon::Loop> const& loops) {
  Init(loops);
}

LaxPolygon::LaxPolygon(S2Polygon const& polygon) {
  Init(polygon);
}

void LaxPolygon::Init(vector<LaxPolygon::Loop> const& loops) {
  vector<VertexArray> v_arrays;
  for (LaxPolygon::Loop const& loop : loops) {
    v_arrays.push_back(VertexArray(loop.data(), loop.size()));
  }
  Init(v_arrays);
}

void LaxPolygon::Init(S2Polygon const& polygon) {
  vector<VertexArray> v_arrays;
  for (int i = 0; i < polygon.num_loops(); ++i) {
    S2Loop const* loop = polygon.loop(i);
    v_arrays.push_back(VertexArray(&loop->vertex(0), loop->num_vertices()));
  }
  Init(v_arrays);
}

void LaxPolygon::Init(vector<VertexArray> const& loops) {
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

void LaxPolygon::GetEdge(int e, S2Point const** a, S2Point const** b) const {
  DCHECK_LT(e, num_edges());
  int const kMaxLinearSearchLoops = 12;  // From benchmarks.

  // Since all loop vertices were concatenated, we don't need to figure out
  // which loop the edge belongs to in order to retrieve its first vertex.
  *a = &vertices_[e];
  if (num_loops() == 1) {
    if (++e == num_vertices_) { e = 0; }
  } else {
    // Find the index of the first vertex of the loop following this one.
    int* next = cumulative_vertices_ + 1;
    if (num_loops() <= kMaxLinearSearchLoops) {
      while (*next <= e) ++next;
    } else {
      next = std::upper_bound(next, next + num_loops(), e);
    }
    // Wrap around to the first vertex of the loop if necessary.
    if (++e == *next) { e = next[-1]; }
  }
  *b = &vertices_[e];
}

bool LaxPolygon::contains_origin() const {
  return IsOriginOnLeft(*this);
}

int LaxPolygon::chain_start(int i) const {
  DCHECK_LE(i, num_loops());
  if (num_loops() == 1) {
    return i == 0 ? 0 : num_vertices_;
  } else {
    return cumulative_vertices_[i];
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
  vertices_.reset(new S2Point[num_vertices_]);
  std::copy(vertices.begin(), vertices.end(), vertices_.get());
}

void LaxPolyline::Init(S2Polyline const& polyline) {
  num_vertices_ = polyline.num_vertices();
  vertices_.reset(new S2Point[num_vertices_]);
  std::copy(&polyline.vertex(0), &polyline.vertex(0) + num_vertices_,
            vertices_.get());
}

void LaxPolyline::GetEdge(int e, S2Point const** a,
                            S2Point const** b) const {
  DCHECK_LT(e, num_edges());
  *a = &vertices_[e];
  *b = &vertices_[e + 1];
}

int LaxPolyline::num_chains() const {
  return std::min(1, LaxPolyline::num_edges());  // Avoid virtual call.
}

int LaxPolyline::chain_start(int i) const {
  return i == 0 ? 0 : LaxPolyline::num_edges();  // Avoid virtual call.
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

void VertexIdLaxLoop::GetEdge(int e, S2Point const** a, S2Point const** b)
    const {
  DCHECK_LT(e, num_edges());
  *a = &vertex(e);
  if (++e == num_vertices()) e = 0;
  *b = &vertex(e);
}

bool VertexIdLaxLoop::contains_origin() const {
  return IsOriginOnLeft(*this);
}

// This is a helper function for IsOriginOnLeft(), defined below.
//
// If the given vertex "vtest" is degenerate (see definition in header file),
// returns false.  Otherwise sets "result" to indicate whether "shape"
// contains S2::Origin() and returns true.
static bool IsOriginOnLeftAtVertex(S2Shape const& shape,
                                   S2Point const& vtest, bool* result) {
  // Let P be a non-degenerate vertex.  Vertex P is defined to be inside the
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
  S2Point origin = S2::Origin();
  S2EdgeUtil::EdgeCrosser crosser(&origin, &vtest);
  bool crossing_parity = false;
  util::btree::btree_map<S2Point, int> edge_map;
  int n = shape.num_edges();
  for (int e = 0; e < n; ++e) {
    S2Point const *v0, *v1;
    shape.GetEdge(e, &v0, &v1);
    if (*v0 == *v1) continue;

    // Check whether this edge crosses the edge between P and S2::Origin().
    crossing_parity ^= crosser.EdgeOrVertexCrossing(v0, v1);

    // Keep track of (outgoing edges) - (incoming edges) for each vertex that
    // is adjacent to "vtest".
    if (*v0 == vtest) ++edge_map[*v1];
    if (*v1 == vtest) --edge_map[*v0];
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
    return false;
  }
  // Define a "matched" edge as one that can be paired with a corresponding
  // reversed edge.  Define a vertex as "degenerate" if all of its edges are
  // matched. In order to determine containment, we must find a non-degenerate
  // vertex.  Often every vertex is non-degenerate, so we start by trying an
  // arbitrary vertex.
  S2Point const *v0, *v1;
  shape.GetEdge(0, &v0, &v1);
  bool result = false;
  if (IsOriginOnLeftAtVertex(shape, *v0, &result)) {
    return result;
  }
  // That didn't work, so now we do some extra work to find a non-degenerate
  // vertex (if any).  Essentially we gather a list of edges and a list of
  // reversed edges, and then sort them.  The first edge that appears in one
  // list but not the other is guaranteed to be unmatched.
  class Edge {
   public:
    Edge(S2Point const* v0, S2Point const* v1) : v0_(v0), v1_(v1) {}
    // Define accessors to ensure that we don't accidentally compare pointers.
    S2Point const& v0() const { return *v0_; }
    S2Point const& v1() const { return *v1_; }
    bool operator<(Edge const& other) const {
      return v0() < other.v0() || (v0() == other.v0() && v1() < other.v1());
    }
   private:
    S2Point const *v0_, *v1_;
  };
  int n = shape.num_edges();
  vector<Edge> edges, rev_edges;
  edges.reserve(n);
  rev_edges.reserve(n);
  for (int i = 0; i < n; ++i) {
    shape.GetEdge(i, &v0, &v1);
    edges.push_back(Edge(v0, v1));
    rev_edges.push_back(Edge(v1, v0));
  }
  std::sort(edges.begin(), edges.end());
  std::sort(rev_edges.begin(), rev_edges.end());
  for (int i = 0; i < n; ++i) {
    if (edges[i] < rev_edges[i]) {  // edges[i] is unmatched
      CHECK(IsOriginOnLeftAtVertex(shape, edges[i].v0(), &result));
      return result;
    }
    if (rev_edges[i] < edges[i]) {  // rev_edges[i] is unmatched
      CHECK(IsOriginOnLeftAtVertex(shape, rev_edges[i].v0(), &result));
      return result;
    }
  }
  return false;  // All vertices are degenerate.
}

std::ostream& operator<<(std::ostream& os, ShapeEdgeId id) {
  return os << id.shape_id() << ":" << id.edge_id();
}

// A structure representing one edge within a S2ShapeIndexCell.
struct ShapeEdge {
  ShapeEdge() : edge_id(-1, -1), v0(nullptr), v1(nullptr) {}
  ShapeEdgeId edge_id;
  S2Point const *v0, *v1;
};

// Returns a vector representing all edges in the given S2ShapeIndexCell.
static void GetShapeEdges(S2ShapeIndex const& index,
                          S2ShapeIndexCell const& cell,
                          vector<ShapeEdge>* shape_edges) {
  int num_edges = 0;
  for (int s = 0; s < cell.num_clipped(); ++s) {
    num_edges += cell.clipped(s).num_edges();
  }
  shape_edges->resize(num_edges);
  int out = 0;
  for (int s = 0; s < cell.num_clipped(); ++s) {
    S2ClippedShape const& clipped = cell.clipped(s);
    int32 shape_id = clipped.shape_id();
    S2Shape* shape = index.shape(shape_id);
    int num_clipped = clipped.num_edges();
    for (int i = 0; i < num_clipped; ++i) {
      ShapeEdge* shape_edge = &(*shape_edges)[out++];
      int32 edge_id = clipped.edge(i);
      shape_edge->edge_id = ShapeEdgeId(shape_id, edge_id);
      shape->GetEdge(edge_id, &shape_edge->v0, &shape_edge->v1);
    }
  }
  DCHECK_EQ(num_edges, out);
}

// Given a vector of edges within an S2ShapeIndexCell, append all pairs of
// crossing edges (of the given CrossingType) to "edge_pairs".
static void AppendCrossingEdgePairs(vector<ShapeEdge> const& shape_edges,
                                    CrossingType type,
                                    EdgePairVector* edge_pairs) {
  int min_crossing_value = (type == CrossingType::ALL) ? 0 : 1;
  int num_edges = shape_edges.size();
  for (int i = 0; i + 1 < num_edges; ++i) {
    ShapeEdge const& a = shape_edges[i];
    S2EdgeUtil::EdgeCrosser crosser(a.v0, a.v1);
    for (int j = i + 1; j < num_edges; ++j) {
      ShapeEdge const& b = shape_edges[j];
      int sign = crosser.CrossingSign(b.v0, b.v1);
      if (sign >= min_crossing_value) {
        edge_pairs->push_back(std::make_pair(a.edge_id, b.edge_id));
      }
    }
  }
}

EdgePairVector GetCrossingEdgePairs(S2ShapeIndex const& index,
                                    CrossingType type) {
  // TODO(ericv): Use brute force if the total number of edges is small enough
  // (using a larger threshold if the S2ShapeIndex is not constructed yet).
  EdgePairVector edge_pairs;
  vector<ShapeEdge> shape_edges;
  for (S2ShapeIndex::Iterator it(index); !it.Done(); it.Next()) {
    GetShapeEdges(index, it.cell(), &shape_edges);
    AppendCrossingEdgePairs(shape_edges, type, &edge_pairs);
  }
  if (edge_pairs.size() > 1) {
    std::sort(edge_pairs.begin(), edge_pairs.end());
    edge_pairs.erase(std::unique(edge_pairs.begin(), edge_pairs.end()),
                     edge_pairs.end());
  }
  return edge_pairs;
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
        index.Add(loop);
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
    S2Point const *v0, *v1;
    DCHECK_GT(loop->num_edges(), 0);
    loop->GetEdge(0, &v0, &v1);
    index.GetContainingShapes(*v0, &ancestors[i]);
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

  // Explicitly remove the shapes from the index so they are not deleted.
  index.RemoveAll();
}

static bool FindSelfIntersection(S2ClippedShape const& a_clipped,
                                 S2Loop const& a_loop, S2Error* error) {
  // Test for crossings between all edge pairs that do not share a vertex.
  // This means that (a) the loop edge indices must differ by 2 or more, and
  // (b) the pair cannot consist of the first and last loop edges.  Part (b)
  // is worthwhile in the case of very small loops; e.g. it reduces the
  // number of crossing tests in a loop with four edges from two to one.

  int a_num_clipped = a_clipped.num_edges();
  for (int i = 0; i < a_num_clipped - 1; ++i) {
    int ai = a_clipped.edge(i);
    int j = i + 1;
    if (a_clipped.edge(j) == ai+1) {  // Adjacent edges
      if (++j >= a_num_clipped) continue;
    }
    S2EdgeUtil::EdgeCrosser crosser(&a_loop.vertex(ai), &a_loop.vertex(ai+1));
    for (int aj_prev = -2; j < a_num_clipped; ++j) {
      int aj = a_clipped.edge(j);
      if (aj - ai == a_loop.num_vertices() - 1) break;  // First and last edges
      if (aj != aj_prev + 1) crosser.RestartAt(&a_loop.vertex(aj));
      aj_prev = aj;
      // This test also catches duplicate vertices.
      int crossing = crosser.CrossingSign(&a_loop.vertex(aj + 1));
      if (crossing < 0) continue;
      if (crossing == 0) {
        error->Init(S2Error::DUPLICATE_VERTICES,
                    "Edge %d has duplicate vertex with edge %d", ai, aj);
      } else {
        error->Init(S2Error::LOOP_SELF_INTERSECTION,
                    "Edge %d crosses edge %d", ai, aj);
      }
      return true;
    }
  }
  return false;
}

// Return true if any of the given loops has a self-intersection (including a
// duplicate vertex), and set "error" to a human-readable error message.
// Otherwise return false and leave "error" unchanged.  All tests are limited
// to edges that intersect the given cell.
inline static bool FindSelfIntersection(
    vector<unique_ptr<S2Loop>> const& loops,
    S2ShapeIndexCell const& cell,
    S2Error* error) {
  for (int a = 0; a < cell.num_clipped(); ++a) {
    S2ClippedShape const& a_clipped = cell.clipped(a);
    if (FindSelfIntersection(a_clipped, *loops[a_clipped.shape_id()], error)) {
      error->Init(error->code(),
                  "Loop %d: %s", a_clipped.shape_id(), error->text().c_str());
      return true;
    }
  }
  return false;
}

// Given two loop edges for which CrossingSign returned a non-negative
// result "crossing", return true if there is a crossing and set "error" to a
// human-readable error message.
static bool GetCrossingError(S2Loop const& a_loop, int a_shape_id, int ai,
                             S2Loop const& b_loop, int b_shape_id, int bj,
                             int crossing, S2Error* error) {
  if (crossing > 0) {
    error->Init(S2Error::POLYGON_LOOPS_CROSS,
                "Loop %d edge %d crosses loop %d edge %d",
                a_shape_id, ai, b_shape_id, bj);
    return true;
  }
  // Loops are not allowed to share edges or cross at vertices.  We
  // only need to check this once per edge pair, so we also require
  // that the two edges have the same end vertex.  (This is only valid
  // because we are iterating over all the cells in the index.)
  if (a_loop.vertex(ai+1) == b_loop.vertex(bj+1)) {
    if (a_loop.vertex(ai) == b_loop.vertex(bj) ||
        a_loop.vertex(ai) == b_loop.vertex(bj+2)) {
      // The second edge index is sometimes off by one, hence "near".
      error->Init(S2Error::POLYGON_LOOPS_SHARE_EDGE,
                  "Loop %d edge %d has duplicate near loop %d edge %d",
                  a_shape_id, ai, b_shape_id, bj);
      return true;
    }
    // Note that we don't need to maintain any state regarding loop
    // crossings because duplicate edges are not allowed.
    if (S2EdgeUtil::GetWedgeRelation(
            a_loop.vertex(ai), a_loop.vertex(ai+1), a_loop.vertex(ai+2),
            b_loop.vertex(bj), b_loop.vertex(bj+2)) ==
        S2EdgeUtil::WEDGE_PROPERLY_OVERLAPS) {
      error->Init(S2Error::POLYGON_LOOPS_CROSS,
                  "Loop %d edge %d crosses loop %d edge %d",
                  a_shape_id, ai, b_shape_id, bj);
      return true;
    }
  }
  return false;
}

// Return true if any of the given loops crosses a different loop (including
// vertex crossings) or two loops share a common edge, and set "error" to a
// human-readable error message.  Otherwise return false and leave "error"
// unchanged.  All tests are limited to edges that intersect the given cell.
static bool FindLoopCrossing(vector<unique_ptr<S2Loop>> const& loops,
                             S2ShapeIndexCell const& cell, S2Error* error) {
  // Possible optimization:
  // Sort the ClippedShapes by edge count to reduce the number of calls to
  // s2pred::Sign.  If n is the total number of shapes in the cell, n_i is
  // the number of edges in shape i, and c_i is the number of continuous
  // chains formed by these edges, the total number of calls is
  //
  //   sum(n_i * (1 + c_j + n_j), i=0..n-2, j=i+1..n-1)
  //
  // So for example if n=2, shape 0 has one chain of 1 edge, and shape 1 has
  // one chain of 8 edges, the number of calls to Sign is 1*10=10 if the
  // shapes are sorted by edge count, and 8*3=24 otherwise.
  //
  // using SortedShape = pair<int, S2ShapeIndex::ClippedShape const*>;
  // vector<SortedShape> sorted_shapes;

  for (int a = 0; a < cell.num_clipped() - 1; ++a) {
    S2ClippedShape const& a_clipped = cell.clipped(a);
    S2Loop const& a_loop = *loops[a_clipped.shape_id()];
    int a_num_clipped = a_clipped.num_edges();
    for (int i = 0; i < a_num_clipped; ++i) {
      int ai = a_clipped.edge(i);
      S2EdgeUtil::EdgeCrosser crosser(&a_loop.vertex(ai), &a_loop.vertex(ai+1));
      for (int b = a + 1; b < cell.num_clipped(); ++b) {
        S2ClippedShape const& b_clipped = cell.clipped(b);
        S2Loop const& b_loop = *loops[b_clipped.shape_id()];
        int bj_prev = -2;
        int b_num_clipped = b_clipped.num_edges();
        for (int j = 0; j < b_num_clipped; ++j) {
          int bj = b_clipped.edge(j);
          if (bj != bj_prev + 1) crosser.RestartAt(&b_loop.vertex(bj));
          bj_prev = bj;
          int crossing = crosser.CrossingSign(&b_loop.vertex(bj + 1));
          if (crossing < 0) continue;  // No crossing
          if (GetCrossingError(a_loop, a_clipped.shape_id(), ai,
                               b_loop, b_clipped.shape_id(), bj,
                               crossing, error)) {
            return true;
          }
        }
      }
    }
  }
  return false;
}

bool FindSelfIntersection(S2ShapeIndex const& index, S2Loop const& loop,
                          S2Error* error) {
  DCHECK_EQ(1, index.num_shape_ids());
  for (S2ShapeIndex::Iterator it(index); !it.Done(); it.Next()) {
    if (FindSelfIntersection(it.cell().clipped(0), loop, error)) {
      return true;
    }
  }
  return false;
}

bool FindAnyCrossing(S2ShapeIndex const& index,
                     vector<unique_ptr<S2Loop>> const& loops,
                     S2Error* error) {
  for (S2ShapeIndex::Iterator it(index); !it.Done(); it.Next()) {
    if (FindSelfIntersection(loops, it.cell(), error)) {
      return true;
    }
    if (it.cell().num_clipped() >= 2 &&
        FindLoopCrossing(loops, it.cell(), error)) {
      return true;
    }
  }
  return false;
}

}  // namespace s2shapeutil
