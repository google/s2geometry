// Copyright 2017 Google Inc. All Rights Reserved.
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
//
// Boundary operations are implemented by constructing the boundary of the
// result and then using S2Builder to assemble the edges.  The boundary is
// obtained by clipping each of the two input regions to the interior or
// exterior of the other region.  For example, to compute the union of A and
// B, we clip the boundary of A to the exterior of B and the boundary of B to
// the exterior of A; the resulting set of edges defines the union of the two
// regions.
//
// We use exact predicates, but inexact constructions (e.g. computing the
// intersection point of two edges).  Nevertheless, the following algorithm is
// guaranteed to be 100% robust, in that the computed boundary stays within a
// small tolerance (snap_radius + S2::kIntersectionError) of the exact
// result, and also preserves the correct topology (i.e., no crossing edges).
//
// Unfortunately this robustness cannot quite be achieved using the strategy
// outlined above (clipping the two input regions and assembling the
// resulting edges).  Since computed intersection points are not exact, the
// input geometry passed to S2Builder might contain self-intersections, and
// these self-intersections cannot be eliminated reliably by snap rounding.
//
// So instead, we pass S2Builder the entire set of input edges where at least
// some portion of each edge belongs to the output boundary.  We allow
// S2Builder to compute the intersection points and snap round the edges
// (which it does in a way that is guaranteed to preserve the input topology).
// Then once this is finished, we remove the portions of each edge that would
// have been clipped if we had done the clipping first.  This step only
// involves deciding whether to keep or discard each edge in the output, since
// all intersection points have already been resolved, and therefore there is
// no risk of creating new self-intersections.
//
// This is implemented using the following classes:
//
//  - S2BoundaryOperation::Impl: the top-level class that clips each of
//                               the two regions to the other region.
//
//  - CrossingProcessor: a class that processes edge crossings and maintains
//                       the necessary state in order to clip the boundary
//                       of one region to the interior or exterior of the
//                       other region.
//
//  - EdgeClippingLayer: an S2Builder::Layer that removes graph edges that
//                       correspond to clipped portions of input edges, and
//                       passes the result to another layer for assembly.
//
//  - GraphEdgeClipper: a helper class that does the actual work of the
//                      EdgeClippingLayer.

#include "s2/s2boundary_operation.h"

#include <limits>
#include <memory>
#include <utility>

#include "s2/third_party/absl/memory/memory.h"
#include "s2/s2builder.h"
#include "s2/s2builder_layer.h"
#include "s2/s2builderutil_snap_functions.h"
#include "s2/s2crossingedgequery.h"
#include "s2/s2measures.h"
#include "s2/s2predicates.h"

// TODO(ericv): Remove this debugging output at some point.
extern bool s2builder_verbose;

namespace {  // Anonymous namespace for helper classes.

using absl::MakeUnique;
using std::make_pair;
using std::max;
using std::min;
using std::pair;
using std::swap;
using std::unique_ptr;
using std::vector;

using EdgeType = S2Builder::EdgeType;
using SnapFunction = S2Builder::SnapFunction;
using GraphOptions = S2Builder::GraphOptions;
using DegenerateEdges = GraphOptions::DegenerateEdges;
using DuplicateEdges = GraphOptions::DuplicateEdges;
using SiblingPairs = GraphOptions::SiblingPairs;

using Graph = S2Builder::Graph;
using EdgeId = Graph::EdgeId;
using VertexId = Graph::VertexId;
using InputEdgeId = Graph::InputEdgeId;
using InputEdgeIdSetId = Graph::InputEdgeIdSetId;

using PolygonModel = S2BoundaryOperation::PolygonModel;
using PolylineModel = S2BoundaryOperation::PolylineModel;

// A collection of special InputEdgeIds that allow the GraphEdgeClipper state
// modifications to be inserted into the list of edge crossings.
static InputEdgeId const kSetInside = -1;
static InputEdgeId const kSetInvertB = -2;
static InputEdgeId const kSetReverseA = -3;

// CrossingInputEdge represents an input edge B that crosses some other input
// edge A.  It stores the input edge id of edge B and also whether it crosses
// edge A from left to right (or vice versa).
class CrossingInputEdge {
 public:
  // Indicates that input edge "input_id" crosses another edge (from left to
  // right if "left_to_right" is true).
  CrossingInputEdge(InputEdgeId input_id, bool left_to_right)
      : left_to_right_(left_to_right), input_id_(input_id) {
  }

  InputEdgeId input_id() const { return input_id_; }
  bool left_to_right() const { return left_to_right_; }

  bool operator<(CrossingInputEdge const& other) const {
    return input_id_ < other.input_id_;
  }
  bool operator<(InputEdgeId const& other) const {
    return input_id_ < other;
  }

 private:
  bool left_to_right_ : 1;
  InputEdgeId input_id_ : 31;
};

// InputEdgeCrossings represents all pairs of intersecting input edges.
// It is sorted in lexicographic order.
using InputEdgeCrossings = vector<pair<InputEdgeId, CrossingInputEdge>>;

// Given two input edges A and B that intersect, suppose that A maps to a
// chain of snapped edges A_0, A_1, ..., A_m and B maps to a chain of snapped
// edges B_0, B_1, ..., B_n.  CrossingGraphEdge represents an edge from chain
// B that shares a vertex with chain A.  It is used as a temporary data
// representation while processing chain A.  The arguments are:
//
//   "id" - the Graph::EdgeId of an edge from chain B.
//   "a_index" - the index of the vertex (A_i) that is shared with chain A.
//   "outgoing" - true if the shared vertex is the first vertex of the B edge.
//   "dst" - the Graph::VertexId of the vertex that is not shared with chain A.
//
// Note that if an edge from the B chain shares both vertices with the A
// chain, there will be two entries: an outgoing edge that treats its first
// vertex as being shared, and an incoming edge that treats its second vertex
// as being shared.
struct CrossingGraphEdge {
  CrossingGraphEdge(EdgeId _id, int _a_index, bool _outgoing, VertexId _dst)
      : id(_id), a_index(_a_index), outgoing(_outgoing), dst(_dst) {
  }
  EdgeId id;
  int a_index;
  bool outgoing;
  VertexId dst;
};
using CrossingGraphEdgeVector = absl::InlinedVector<CrossingGraphEdge, 2>;

// Returns a vector of EdgeIds sorted by input edge id.  When more than one
// output edge has the same input edge id (i.e., the input edge snapped to a
// chain of edges), the edges are sorted so that they form a directed edge
// chain.
//
// This function could possibily be moved to S2Builder::Graph, but note that
// it has special requirements.  Namely, duplicate edges and sibling pairs
// must be kept in order to ensure that every output edge corresponds to
// exactly one input edge.  (See also S2Builder::Graph::GetInputEdgeOrder.)
static vector<EdgeId> GetInputEdgeChainOrder(
    Graph const& g, vector<InputEdgeId> const& input_ids) {

  DCHECK(g.options().edge_type() == EdgeType::DIRECTED);
  DCHECK(g.options().duplicate_edges() == DuplicateEdges::KEEP);
  DCHECK(g.options().sibling_pairs() == SiblingPairs::KEEP);

  // First, sort the edges so that the edges corresponding to each input edge
  // are consecutive.  (Each input edge was snapped to a chain of output
  // edges, or two chains in the case of undirected input edges.)
  vector<EdgeId> order = g.GetInputEdgeOrder(input_ids);

  // Now sort the group of edges corresponding to each input edge in edge
  // chain order (e.g.  AB, BC, CD).
  vector<pair<VertexId, EdgeId>> vmap;     // Map from source vertex to edge id.
  vector<int> indegree(g.num_vertices());  // Restricted to current input edge.
  for (int end, begin = 0; begin < order.size(); begin = end) {
    // Gather the edges that came from a single input edge.
    InputEdgeId input_id = input_ids[order[begin]];
    for (end = begin; end < order.size(); ++end) {
      if (input_ids[order[end]] != input_id) break;
    }
    if (end - begin == 1) continue;

    // Build a map from the source vertex of each edge to its edge id,
    // and also compute the indegree at each vertex considering only the edges
    // that came from the current input edge.
    for (int i = begin; i < end; ++i) {
      EdgeId e = order[i];
      vmap.push_back(make_pair(g.edge(e).first, e));
      indegree[g.edge(e).second] += 1;
    }
    std::sort(vmap.begin(), vmap.end());

    // Find the starting edge for building the edge chain.
    EdgeId next = g.num_edges();
    for (int i = begin; i < end; ++i) {
      EdgeId e = order[i];
      if (indegree[g.edge(e).first] == 0) next = e;
    }
    // Build the edge chain.
    for (int i = begin; ;) {
      order[i] = next;
      VertexId v = g.edge(next).second;
      indegree[v] = 0;  // Clear as we go along.
      if (++i == end) break;
      auto out = lower_bound(vmap.begin(), vmap.end(), make_pair(v, 0));
      DCHECK_EQ(v, out->first);
      next = out->second;
    }
    vmap.clear();
  }
  return order;
}

// Given a set of clipping instructions encoded as a set of InputEdgeCrossings,
// GraphEdgeClipper determines which graph edges correspond to clipped
// portions of input edges and removes them.
//
// The clipping model is as follows.  The input consists of edge chains.  The
// clipper maintains an "inside" boolean state as it clips each chain, and
// toggles this state whenever an input edge is crossed.  Any edges that are
// deemed to be "outside" after clipping are removed.
//
// The "inside" state can be reset when necessary (e.g., when jumping to the
// start of a new chain) by adding a special crossing marked kSetInside.
// There are also two other special "crossings" that modify the clipping
// parameters: kSetInvertB specifies that edges should be clipped to the
// exterior of the other region, and kSetReverseA specifies that edges should
// be reversed before emitting them (which is needed to implement difference
// operations).
class GraphEdgeClipper {
 public:
  // "input_dimensions" is a vector specifying the dimension of each input
  // edge (0, 1, or 2).  "input_crossings" is the set of all crossings to be
  // used when clipping the edges of "g", sorted in lexicographic order.
  //
  // The clipped set of edges and their corresponding set of input edge ids
  // are returned in "new_edges" and "new_input_edge_ids".  (These can be used
  // to construct a new S2Builder::Graph.)
  GraphEdgeClipper(Graph const& g, vector<int8> const& input_dimensions,
                   InputEdgeCrossings const& input_crossings,
                   vector<Graph::Edge>* new_edges,
                   vector<InputEdgeIdSetId>* new_input_edge_ids);
  void Run();

 private:
  void AddEdge(Graph::Edge edge, InputEdgeId input_edge_id);
  void GatherIncidentEdges(
      vector<VertexId> const& a, int ai,
      vector<CrossingInputEdge> const& b_input_edges,
      vector<CrossingGraphEdgeVector>* b_edges) const;
  int GetCrossedVertexIndex(
      vector<VertexId> const& a, CrossingGraphEdgeVector const& b,
      bool left_to_right) const;
  int GetVertexRank(CrossingGraphEdge const& e) const;
  bool EdgeChainOnLeft(vector<VertexId> const& a,
                       EdgeId b_first, EdgeId b_last) const;

  Graph const& g_;
  Graph::VertexInMap in_;
  Graph::VertexOutMap out_;
  vector<int8> const& input_dimensions_;
  InputEdgeCrossings const& input_crossings_;
  vector<Graph::Edge>* new_edges_;
  vector<InputEdgeIdSetId>* new_input_edge_ids_;

  // Every graph edge is associated with exactly one input edge in our case,
  // which means that we can declare g_.input_edge_id_set_ids() as a vector of
  // InputEdgeIds rather than a vector of InputEdgeIdSetIds.  (This also takes
  // advantage of the fact that IdSetLexicon represents a singleton set as the
  // value of its single element.)
  vector<InputEdgeId> const& input_ids_;

  vector<EdgeId> order_;  // Graph edges sorted in input edge id order.
  vector<int> rank_;      // The rank of each graph edge within order_.
};

GraphEdgeClipper::GraphEdgeClipper(
    Graph const& g, vector<int8> const& input_dimensions,
    InputEdgeCrossings const& input_crossings,
    vector<Graph::Edge>* new_edges,
    vector<InputEdgeIdSetId>* new_input_edge_ids)
    : g_(g), in_(g), out_(g),
      input_dimensions_(input_dimensions),
      input_crossings_(input_crossings),
      new_edges_(new_edges),
      new_input_edge_ids_(new_input_edge_ids),
      input_ids_(g.input_edge_id_set_ids()),
      order_(GetInputEdgeChainOrder(g_, input_ids_)),
      rank_(order_.size()) {
  for (int i = 0; i < order_.size(); ++i) {
    rank_[order_[i]] = i;
  }
}

inline void GraphEdgeClipper::AddEdge(Graph::Edge edge,
                                      InputEdgeId input_edge_id) {
  new_edges_->push_back(edge);
  new_input_edge_ids_->push_back(input_edge_id);
}

void GraphEdgeClipper::Run() {
  // Declare vectors here and reuse them to avoid reallocation.
  vector<VertexId> a_vertices;
  vector<int> a_num_crossings;
  vector<bool> a_isolated;
  vector<CrossingInputEdge> b_input_edges;
  vector<CrossingGraphEdgeVector> b_edges;

  bool inside = false;
  bool invert_b = false;
  bool reverse_a = false;
  auto next = input_crossings_.begin();
  for (int i = 0; i < order_.size(); ++i) {
    // For each input edge (the "A" input edge), gather all the input edges
    // that cross it (the "B" input edges).
    InputEdgeId a_input_id = input_ids_[order_[i]];
    Graph::Edge const& edge0 = g_.edge(order_[i]);
    b_input_edges.clear();
    for (; next != input_crossings_.end(); ++next) {
      if (next->first != a_input_id) break;
      if (next->second.input_id() >= 0) {
        b_input_edges.push_back(next->second);
      } else if (next->second.input_id() == kSetInside) {
        inside = next->second.left_to_right();
      } else if (next->second.input_id() == kSetInvertB) {
        invert_b = next->second.left_to_right();
      } else {
        DCHECK_EQ(next->second.input_id(), kSetReverseA);
        reverse_a = next->second.left_to_right();
      }
    }
    // Optimization for degenerate edges.
    // TODO(ericv): If the output layer for this edge dimension specifies
    // DegenerateEdges::DISCARD, then remove the edge here.
    if (edge0.first == edge0.second) {
      inside ^= (b_input_edges.size() & 1);
      AddEdge(edge0, a_input_id);
      continue;
    }
    // Optimization for the case where there are no crossings.
    if (b_input_edges.empty()) {
      // In general the caller only passes edges that are part of the output
      // (i.e., we could DCHECK(inside) here).  The one exception is for
      // polyline/polygon operations, where the polygon edges are needed to
      // compute the polyline output but are not emitted themselves.
      if (inside) {
        AddEdge(reverse_a ? Graph::reverse(edge0) : edge0, a_input_id);
      }
      continue;
    }
    // Walk along the chain of snapped edges for input edge A, and at each
    // vertex collect all the incident edges that belong to one of the
    // crossing edge chains (the "B" input edges).
    a_vertices.clear();
    a_vertices.push_back(edge0.first);
    b_edges.clear();
    b_edges.resize(b_input_edges.size());
    GatherIncidentEdges(a_vertices, 0, b_input_edges, &b_edges);
    for (; i < order_.size() && input_ids_[order_[i]] == a_input_id; ++i) {
      a_vertices.push_back(g_.edge(order_[i]).second);
      GatherIncidentEdges(a_vertices, a_vertices.size() - 1, b_input_edges,
                          &b_edges);
    }
    --i;
    if (s2builder_verbose) {
      std::cout << "input edge " << a_input_id << " (inside=" << inside << "):";
      for (VertexId id : a_vertices) std::cout << " " << id;
    }
    // Now for each B edge chain, decide which vertex of the A chain it
    // crosses, and keep track of the number of signed crossings at each A
    // vertex.  The sign of a crossing depends on whether the other edge
    // crosses from left to right or right to left.
    //
    // This would not be necessary if all calculations were done in exact
    // arithmetic, because crossings would have strictly alternating signs.
    // But because we have already snapped the result, some crossing locations
    // are ambiguous, and GetCrossedVertexIndex() handles this by choosing a
    // candidate vertex arbitrarily.  The end result is that rarely, we may
    // see two crossings in a row with the same sign.  We correct for this by
    // adding extra output edges that essentially link up the crossings in the
    // correct (alternating sign) order.  Compared to the "correct" behavior,
    // the only difference is that we have added some extra sibling pairs
    // (consisting of an edge and its corresponding reverse edge) which do not
    // affect the result.
    a_num_crossings.clear();
    a_num_crossings.resize(a_vertices.size());
    a_isolated.clear();
    a_isolated.resize(a_vertices.size());
    for (int bi = 0; bi < b_input_edges.size(); ++bi) {
      bool left_to_right = b_input_edges[bi].left_to_right();
      int a_index = GetCrossedVertexIndex(a_vertices, b_edges[bi],
                                          left_to_right);
      if (s2builder_verbose) {
        std::cout << std::endl << "  " << "b input edge "
                  << b_input_edges[bi].input_id() << " (l2r=" << left_to_right
                  << ", crossing=" << a_vertices[a_index] << ")";
        for (auto const& x : b_edges[bi]) {
          Graph::Edge const& e = g_.edge(x.id);
          std::cout << " (" << e.first << ", " << e.second << ")";
        }
      }
      // Keep track of the number of signed crossings (see above).
      bool is_line = input_dimensions_[b_input_edges[bi].input_id()] == 1;
      int sign = is_line ? 0 : (left_to_right == invert_b) ? -1 : 1;
      a_num_crossings[a_index] += sign;

      // Any polyline or polygon vertex that has at least one crossing but no
      // adjacent emitted edge may be emitted as an isolated vertex.
      a_isolated[a_index] = true;
    }
    if (s2builder_verbose) std::cout << std::endl;

    // Finally, we iterate through the A edge chain, keeping track of the
    // number of signed crossings as we go along.  The "multiplicity" is
    // defined as the cumulative number of signed crossings, and indicates how
    // many edges should be output (and in which direction) in order to link
    // up the edge crossings in the correct order.  (The multiplicity is
    // almost always either 0 or 1 except in very rare cases.)
    int multiplicity = inside + a_num_crossings[0];
    for (int ai = 1; ai < a_vertices.size(); ++ai) {
      if (multiplicity != 0) {
        a_isolated[ai - 1] = a_isolated[ai] = false;
      }
      int edge_count = reverse_a ? -multiplicity : multiplicity;
      // Output any forward edges required.
      for (int i = 0; i < edge_count; ++i) {
        AddEdge(Graph::Edge(a_vertices[ai - 1], a_vertices[ai]), a_input_id);
      }
      // Output any reverse edges required.
      for (int i = edge_count; i < 0; ++i) {
        AddEdge(Graph::Edge(a_vertices[ai], a_vertices[ai - 1]), a_input_id);
      }
      multiplicity += a_num_crossings[ai];
    }
    // Multiplicities other than 0 or 1 can only occur in the edge interior.
    DCHECK(multiplicity == 0 || multiplicity == 1);
    inside = (multiplicity != 0);

    // Output any isolated polyline vertices.
    // TODO(ericv): Only do this if an output layer wants degenerate edges.
    if (input_dimensions_[a_input_id] != 0) {
      for (int ai = 0; ai < a_vertices.size(); ++ai) {
        if (a_isolated[ai]) {
          AddEdge(Graph::Edge(a_vertices[ai], a_vertices[ai]), a_input_id);
        }
      }
    }
  }
}

// Given the vertices of the snapped edge chain for an input edge A and the
// set of input edges B that cross input edge A, this method gathers all of
// the snapped edges of B that are incident to a given snapped vertex of A.
// The incident edges for each input edge of B are appended to a separate
// output vector.  (A and B can refer to either the input edge or the
// corresponding snapped edge chain.)
void GraphEdgeClipper::GatherIncidentEdges(
    vector<VertexId> const& a, int ai,
    vector<CrossingInputEdge> const& b_input_edges,
    vector<CrossingGraphEdgeVector>* b_edges) const {
  // Examine all of the edges incident to the given vertex of A.  If any edge
  // comes from a B input edge, append it to the appropriate vector.
  DCHECK_EQ(b_input_edges.size(), b_edges->size());
  for (EdgeId e : in_.edge_ids(a[ai])) {
    InputEdgeId id = input_ids_[e];
    auto it = lower_bound(b_input_edges.begin(), b_input_edges.end(), id);
    if (it != b_input_edges.end() && it->input_id() == id) {
      auto& edges = (*b_edges)[it - b_input_edges.begin()];
      edges.push_back(CrossingGraphEdge(e, ai, false, g_.edge(e).first));
    }
  }
  for (EdgeId e : out_.edge_ids(a[ai])) {
    InputEdgeId id = input_ids_[e];
    auto it = lower_bound(b_input_edges.begin(), b_input_edges.end(), id);
    if (it != b_input_edges.end() && it->input_id() == id) {
      auto& edges = (*b_edges)[it - b_input_edges.begin()];
      edges.push_back(CrossingGraphEdge(e, ai, true, g_.edge(e).second));
    }
  }
}

// Returns the "vertex rank" of the shared vertex associated with the given
// CrossingGraphEdge.  Recall that graph edges are sorted in input edge order,
// and that the rank of an edge is its position in this order (rank_[e]).
// VertexRank(e) is defined such that VertexRank(e.src) == rank_[e] and
// VertexRank(e.dst) == rank_[e] + 1.  Note that the concept of "vertex rank"
// is only defined within a single edge chain (since different edge chains can
// have overlapping vertex ranks).
int GraphEdgeClipper::GetVertexRank(CrossingGraphEdge const& e) const {
  return rank_[e.id] + !e.outgoing;
}

// Given an edge chain A that is crossed by another edge chain B (where
// "left_to_right" indicates whether B crosses A from left to right), this
// method decides which vertex of A the crossing takes place at.  The
// parameters are the vertices of the A chain ("a") and the set of edges in
// the B chain ("b") that are incident to vertices of A.  The B chain edges
// are sorted in increasing order of (a_index, outgoing) tuple.
int GraphEdgeClipper::GetCrossedVertexIndex(
    vector<VertexId> const& a, CrossingGraphEdgeVector const& b,
    bool left_to_right) const {
  DCHECK(!a.empty());
  DCHECK(!b.empty());

  // The reason this calculation is tricky is that after snapping, the A and B
  // chains may meet and separate several times.  For example, if B crosses A
  // from left to right, then B may touch A, make an excursion to the left of
  // A, come back to A, then make an excursion to the right of A and come back
  // to A again, like this:
  //
  //  *--B--*-\             /-*-\
  //           B-\       /-B     B-\      6     7     8     9
  //  *--A--*--A--*-A,B-*--A--*--A--*-A,B-*--A--*--A--*-A,B-*
  //  0     1     2     3     4     5      \-B     B-/
  //                                          \-*-/
  //
  // (where "*" is a vertex, and "A" and "B" are edge labels).  Note that B
  // may also follow A for one or more edges whenever they touch (e.g. between
  // vertices 2 and 3 ).  In this case the only vertices of A where the
  // crossing could take place are 5 and 6, i.e. after all excursions of B to
  // the left of A, and before all excursions of B to the right of A.
  //
  // Other factors to consider are that the portion of B before and/or after
  // the crossing may be degenerate, and some or all of the B edges may be
  // reversed relative to the A edges.

  // First, check whether edge A is degenerate.
  int n = a.size();
  if (n == 1) return 0;

  // If edge chain B is incident to only one vertex of A, we're done.
  if (b[0].a_index == b.back().a_index) return b[0].a_index;

  // Determine whether the B chain visits the first and last vertices that it
  // shares with the A chain in the same order or the reverse order.  This is
  // only needed to implement one special case (see below).
  bool b_reversed = GetVertexRank(b[0]) > GetVertexRank(b.back());

  // Examine each incident B edge and use it to narrow the range of positions
  // where the crossing could occur in the B chain.  Vertex positions are
  // represented as a range [lo, hi] of vertex ranks in the B chain (see
  // GetVertexRank).
  //
  // Note that if an edge of B is incident to the first or last vertex of A,
  // we can't test which side of the A chain it is on.  There can be up to 4
  // such edges (one incoming and one outgoing edge at each vertex).  Two of
  // these edges logically extend past the end of the A chain and place no
  // restrictions on the crossing vertex.  The other two edges define the ends
  // of the subchain where B shares vertices with A.  We save these edges in
  // order to handle a special case (see below).
  int lo = -1, hi = order_.size();   // Vertex ranks of acceptable crossings
  EdgeId b_first = -1, b_last = -1;  // "b" subchain connecting "a" endpoints
  for (auto const& e : b) {
    int ai = e.a_index;
    if (ai == 0) {
      if (e.outgoing != b_reversed && e.dst != a[1]) b_first = e.id;
    } else if (ai == n - 1) {
      if (e.outgoing == b_reversed && e.dst != a[n - 2]) b_last = e.id;
    } else {
      // This B edge is incident to an interior vertex of the A chain.  First
      // check whether this edge is identical (or reversed) to an edge in the
      // A chain, in which case it does not create any restrictions.
      if (e.dst == a[ai - 1] || e.dst == a[ai + 1]) continue;

      // Otherwise we can test which side of the A chain the edge lies on.
      bool on_left = s2pred::OrderedCCW(g_.vertex(a[ai + 1]), g_.vertex(e.dst),
                                        g_.vertex(a[ai - 1]), g_.vertex(a[ai]));

      // Every B edge that is incident to an interior vertex of the A chain
      // places some restriction on where the crossing vertex could be.
      if (left_to_right == on_left) {
        // This is a pre-crossing edge, so the crossing cannot be before the
        // destination vertex of this edge.  (For example, the input B edge
        // crosses the input A edge from left to right and this edge of the B
        // chain is to the left of the A chain.)
        lo = max(lo, rank_[e.id] + 1);
      } else {
        // This is a post-crossing edge, so the crossing cannot be after the
        // source vertex of this edge.
        hi = min(hi, rank_[e.id]);
      }
    }
  }
  // There is one special case.  If there were no B edges incident to interior
  // vertices of A, then we can't reliably test which side of A the B edges
  // are on.  (An s2pred::Sign test doesn't work, since an edge of B can snap to
  // the "wrong" side of A while maintaining topological guarantees.)  So
  // instead we construct a loop consisting of the A edge chain plus the
  // portion of the B chain that connects the endpoints of A.  We can then
  // test the orientation of this loop.
  //
  // Note that it would be possible to avoid this in some situations by
  // testing whether either endpoint of the A chain has two incident B edges,
  // in which case we could check which side of the B chain the A edge is on
  // and use this to limit the possible crossing locations.
  if (lo < 0 && hi >= order_.size() && b_first >= 0 && b_last >= 0) {
    // Swap the edges if necessary so that they are in B chain order.
    if (b_reversed) swap(b_first, b_last);
    bool on_left = EdgeChainOnLeft(a, b_first, b_last);
    if (left_to_right == on_left) {
      lo = max(lo, rank_[b_last] + 1);
    } else {
      hi = min(hi, rank_[b_first]);
    }
  }

  // Otherwise we choose the smallest shared VertexId in the acceptable range,
  // in order to ensure that both chains choose the same crossing vertex.
  int best = -1;
  DCHECK_LE(lo, hi);
  for (auto const& e : b) {
    int ai = e.a_index;
    int vrank = GetVertexRank(e);
    if (vrank >= lo && vrank <= hi && (best < 0 || a[ai] < a[best])) {
      best = ai;
    }
  }
  return best;
}

// Given edge chains A and B that form a loop (after possibly reversing the
// direction of chain B), returns true if chain B is to the left of chain A.
// Chain A is given as a sequence of vertices, while chain B is specified as
// the first and last edges of the chain.
bool GraphEdgeClipper::EdgeChainOnLeft(
    vector<VertexId> const& a, EdgeId b_first, EdgeId b_last) const {
  // Gather all the interior vertices of the B subchain.
  vector<VertexId> loop;
  for (int i = rank_[b_first]; i < rank_[b_last]; ++i) {
    loop.push_back(g_.edge(order_[i]).second);
  }
  // Possibly reverse the chain so that it forms a loop when "a" is appended.
  if (g_.edge(b_last).second != a[0]) std::reverse(loop.begin(), loop.end());
  loop.insert(loop.end(), a.begin(), a.end());
  // Duplicate the first two vertices to simplify vertex indexing.
  loop.insert(loop.end(), loop.begin(), loop.begin() + 2);
  // Now B is to the left of A if and only if the loop is counterclockwise.
  double sum = 0;
  for (int i = 2; i < loop.size(); ++i) {
    sum += S2::TurnAngle(g_.vertex(loop[i - 2]), g_.vertex(loop[i - 1]),
                         g_.vertex(loop[i]));
  }
  return sum > 0;
}

// Given a set of clipping instructions encoded as a set of intersections
// between input edges, EdgeClippingLayer determines which graph edges
// correspond to clipped portions of input edges and removes them.  It
// assembles the remaining edges into a new S2Builder::Graph and passes the
// result to the given output layer for assembly.
class EdgeClippingLayer : public S2Builder::Layer {
 public:
  EdgeClippingLayer(vector<unique_ptr<S2Builder::Layer>> const* layers,
                    vector<int8> const* input_dimensions,
                    InputEdgeCrossings const* input_crossings)
      : layers_(*layers),
        input_dimensions_(*input_dimensions),
        input_crossings_(*input_crossings) {
  }

  // Layer interface:
  GraphOptions graph_options() const override;
  void Build(Graph const& g, S2Error* error) override;

 private:
  vector<unique_ptr<S2Builder::Layer>> const& layers_;
  vector<int8> const& input_dimensions_;
  InputEdgeCrossings const& input_crossings_;
};

GraphOptions EdgeClippingLayer::graph_options() const {
  // We keep all edges, including degenerate ones, so that we can figure out
  // the correspondence between input edge crossings and output edge
  // crossings.
  return GraphOptions(EdgeType::DIRECTED, DegenerateEdges::KEEP,
                      DuplicateEdges::KEEP, SiblingPairs::KEEP);
}

// Helper function (in anonymous namespace) to create an S2Builder::Graph from
// a vector of edges.
Graph MakeGraph(
    Graph const& g, GraphOptions* options, vector<Graph::Edge>* new_edges,
    vector<InputEdgeIdSetId>* new_input_edge_ids,
    IdSetLexicon* new_input_edge_id_set_lexicon, S2Error* error) {
  if (options->edge_type() == EdgeType::UNDIRECTED) {
    // Create a reversed edge for every edge.
    int n = new_edges->size();
    new_edges->reserve(2 * n);
    new_input_edge_ids->reserve(2 * n);
    for (int i = 0; i < n; ++i) {
      new_edges->push_back(Graph::reverse((*new_edges)[i]));
      new_input_edge_ids->push_back(IdSetLexicon::EmptySetId());
    }
  }
  Graph::ProcessEdges(options, new_edges, new_input_edge_ids,
                      new_input_edge_id_set_lexicon, error);
  return Graph(*options, &g.vertices(), new_edges, new_input_edge_ids,
               new_input_edge_id_set_lexicon, &g.label_set_ids(),
               &g.label_set_lexicon(), g.is_full_polygon_predicate());
}

void EdgeClippingLayer::Build(Graph const& g, S2Error* error) {
  // The bulk of the work is handled by GraphEdgeClipper.
  vector<Graph::Edge> new_edges;
  vector<InputEdgeIdSetId> new_input_edge_ids;
  // Destroy the GraphEdgeClipper immediately to save memory.
  GraphEdgeClipper(g, input_dimensions_, input_crossings_,
                   &new_edges, &new_input_edge_ids).Run();
  if (s2builder_verbose) {
    std::cout << "Edges after clipping: " << std::endl;
    for (int i = 0; i < new_edges.size(); ++i) {
      std::cout << "  " << new_input_edge_ids[i] << " (" << new_edges[i].first
                << ", " << new_edges[i].second << ")" << std::endl;
    }
  }
  // Construct one or more graphs from the clipped edges and pass them to the
  // given output layer(s).
  IdSetLexicon new_input_edge_id_set_lexicon;
  if (layers_.size() == 1) {
    GraphOptions options = layers_[0]->graph_options();
    Graph new_graph = MakeGraph(g, &options, &new_edges, &new_input_edge_ids,
                                &new_input_edge_id_set_lexicon, error);
    layers_[0]->Build(new_graph, error);
  } else {
    // The Graph objects must be valid until the last Build() call completes,
    // so we store all of the graph data in arrays with 3 elements.
    DCHECK_EQ(3, layers_.size());
    vector<Graph::Edge> layer_edges[3];
    vector<InputEdgeIdSetId> layer_input_edge_ids[3];
    S2Builder::GraphOptions layer_options[3];
    vector<S2Builder::Graph> layer_graphs;  // No default constructor.
    layer_graphs.reserve(3);
    // Separate the edges according to their dimension.
    for (int i = 0; i < new_edges.size(); ++i) {
      int d = input_dimensions_[new_input_edge_ids[i]];
      layer_edges[d].push_back(new_edges[i]);
      layer_input_edge_ids[d].push_back(new_input_edge_ids[i]);
    }
    // Clear variables to save space.
    vector<Graph::Edge>().swap(new_edges);
    vector<InputEdgeIdSetId>().swap(new_input_edge_ids);
    for (int d = 0; d < 3; ++d) {
      layer_options[d] = layers_[d]->graph_options();
      layer_graphs.push_back(MakeGraph(
          g, &layer_options[d], &layer_edges[d], &layer_input_edge_ids[d],
          &new_input_edge_id_set_lexicon, error));
      layers_[d]->Build(layer_graphs[d], error);
    }
  }
}

}  // namespace

// TODO(ericv): Expand this into a real public class, and remove the
// ShapeContains and GetContainingShapes methods from S2ShapeIndex.
class S2ContainsPointQuery {
 public:
  explicit S2ContainsPointQuery(S2ShapeIndex const* index) : index_(*index) {}
  bool Contains(S2Point const& p) const {
    return index_.GetContainingShapes(p, &dummy_);
  }
 private:
  S2ShapeIndex const& index_;
  mutable vector<S2Shape*> dummy_;
};

class S2BoundaryOperation::Impl {
 public:
  explicit Impl(S2BoundaryOperation* op)
      : op_(op), index_crossings_first_region_id_(-1) {
  }
  bool Build(S2Error* error);

 private:
  class CrossingIterator;
  class CrossingProcessor;
  using ShapeEdgeId = s2shapeutil::ShapeEdgeId;

  bool is_boolean_output() const { return op_->result_non_empty_ != nullptr; }
  bool AddBoundary(int a_region_id, bool invert_a, bool invert_b,
                   bool invert_result,
                   vector<ShapeEdgeId> const& a_chain_starts,
                   CrossingProcessor* cp);
  bool GetChainStarts(int a_region_id, bool invert_a, bool invert_b,
                      bool invert_result,
                      vector<ShapeEdgeId>* chain_starts) const;
  static bool HasInterior(S2ShapeIndex const& index);
  bool GetIndexCrossings(int region_id);
  bool AddBoundaryPair(bool invert_a, bool invert_b, bool invert_result,
                       CrossingProcessor* cp);
  bool BuildOpType(OpType op_type);

  S2BoundaryOperation* op_;

  // The S2Builder used to construct the output.
  unique_ptr<S2Builder> builder_;

  // A vector specifying the dimension of each edge added to S2Builder.
  vector<int8> input_dimensions_;

  // The set of all input edge crossings, which is used by EdgeClippingLayer
  // to construct the clipped output polygon.
  InputEdgeCrossings input_crossings_;

  // kSentinel is a sentinel value used to mark the end of vectors.
  static ShapeEdgeId const kSentinel;

  // A vector containing all pairs of crossing edges from the two input
  // regions (including edge pairs that share a common vertex).  The first
  // element of each pair is an edge from "index_crossings_first_region_id_",
  // while the second element of each pair is an edge from the other region.
  using EdgePairVector = vector<pair<ShapeEdgeId, ShapeEdgeId>>;
  EdgePairVector index_crossings_;

  // Indicates that the first element of each crossing edge pair in
  // "index_crossings_" corresponds to an edge from the given region.
  // This field is negative if index_crossings_ has not been computed yet.
  int index_crossings_first_region_id_;
};

s2shapeutil::ShapeEdgeId const S2BoundaryOperation::Impl::kSentinel(
    std::numeric_limits<int32>::max(), 0);

// A helper class for iterating through the edges from region B that cross a
// particular edge from region A.  It caches information from the current
// shape, chain, and edge so that it doesn't need to be looked up repeatedly.
// Typical usage:
//
//  void SomeFunction(ShapeEdgeId a_id, CrossingIterator *it) {
//    // Iterate through the edges that cross edge "a_id".
//    for (; !it->Done(a_id); it->Next()) {
//      ... use it->b_shape(), it->b_edge(), etc ...
//    }
class S2BoundaryOperation::Impl::CrossingIterator {
 public:
  CrossingIterator(S2ShapeIndex const* b_index,
                   EdgePairVector const& crossings)
      : b_index_(*b_index), it_(crossings.begin()), b_shape_id_(-1) {
    Update();
  }
  void Next() {
    ++it_;
    Update();
  }
  bool Done(ShapeEdgeId id) const { return a_id() != id; }

  ShapeEdgeId a_id() const { return it_->first; }
  ShapeEdgeId b_id() const { return it_->second; }
  S2ShapeIndex const& b_index() const { return b_index_; }
  S2Shape const& b_shape() const { return *b_shape_; }
  int b_dimension() const { return b_dimension_; }
  int b_shape_id() const { return b_shape_id_; }
  int b_edge_id() const { return b_id().edge_id; }

  S2Shape::Edge b_edge() const {
    return b_shape_->edge(b_edge_id());  // Opportunity to cache this.
  }

  // Information about the chain to which an edge belongs.
  struct ChainInfo {
    int chain_id;  // chain id
    int start;     // starting edge id
    int limit;     // limit edge id
  };
  // Returns a description of the chain to which the current B edge belongs.
  ChainInfo const& b_chain_info() const {
    if (b_info_.chain_id < 0) {
      b_info_.chain_id = b_shape().chain_position(b_edge_id()).chain_id;
      auto chain = b_shape().chain(b_info_.chain_id);
      b_info_.start = chain.start;
      b_info_.limit = chain.start + chain.length;
    }
    return b_info_;
  }

 private:
  // Updates information about the B shape whenever it changes.
  void Update() {
    if (a_id() != kSentinel && b_id().shape_id != b_shape_id_) {
      b_shape_id_ = b_id().shape_id;
      b_shape_ = b_index_.shape(b_shape_id_);
      b_dimension_ = b_shape_->dimension();
      b_info_.chain_id = -1;  // Computed on demand.
    }
  }

  S2ShapeIndex const& b_index_;
  EdgePairVector::const_iterator it_;
  int b_shape_id_;
  S2Shape const* b_shape_;
  int b_dimension_;
  mutable ChainInfo b_info_;  // Computed on demand.
};

// CrossingProcessor is a helper class that processes all the edges from one
// region that cross a specific edge of the other region.  It outputs the
// appropriate edges to an S2Builder, and outputs other information required
// by GraphEdgeClipper to the given vectors.
class S2BoundaryOperation::Impl::CrossingProcessor {
 public:
  // Prepares to build output for the given polygon and polyline boundary
  // models.  Edges are emitted to "builder", while other auxiliary data is
  // appended to the given vectors.
  //
  // If a predicate is being evaluated (i.e., we do not need to construct the
  // actual result), then "builder" and the various output vectors should all
  // be nullptr.
  CrossingProcessor(PolygonModel const& polygon_model,
                    PolylineModel const& polyline_model,
                    S2Builder* builder,
                    vector<int8>* input_dimensions,
                    InputEdgeCrossings *input_crossings)
      : polygon_model_(polygon_model), polyline_model_(polyline_model),
        builder_(builder), input_dimensions_(input_dimensions),
        input_crossings_(input_crossings), prev_inside_(false) {
  }

  // Starts processing edges from the given region.  "invert_a", "invert_b",
  // and "invert_result" indicate whether region A, region B, and/or the
  // result should be inverted, which allows operations such as union and
  // difference to be implemented.  For example, union is ~(~A & ~B).
  //
  // This method should be called in pairs, once to process the edges from
  // region A and once to process the edges from region B.
  void StartBoundary(int a_region_id, bool invert_a, bool invert_b,
                     bool invert_result);

  // Starts processing edges from the given shape.
  void StartShape(S2Shape const* a_shape);

  // Starts processing edges from the given chain.
  void StartChain(int chain_id, int chain_start, int chain_limit, bool inside);

  // Processes the given edge "a_id".  "it" should be positioned to the set of
  // edges from the other region that cross "a_id" (if any).
  bool ProcessEdge(ShapeEdgeId a_id, CrossingIterator* it);

  // This method should be called after each pair of calls to StartBoundary.
  // (The only operation that processes more than one pair of boundaries is
  // SYMMETRIC_DIFFERENCE, which computes the union of A-B and B-A.)
  void DoneBoundaryPair();

  // Indicates whether the point being processed along the current edge chain
  // is in the polygonal interior of the opposite region, using semi-open
  // boundaries.  If "invert_b_" is true then this field is inverted.
  //
  // This value along with the set of incident edges can be used to compute
  // whether the opposite region contains this point under any of the
  // supported boundary models (PolylineModel::CLOSED, etc).
  bool inside() const { return inside_; }

 private:
  // SourceEdgeCrossing represents an input edge that crosses some other
  // edge; it crosses the edge from left to right iff the second parameter
  // is "true".
  using SourceEdgeCrossing = pair<SourceId, bool>;
  class PointCrossingResult;
  class EdgeCrossingResult;

  InputEdgeId input_edge_id() const { return input_dimensions_->size(); }

  // Returns true if the edges on either side of the first vertex of the
  // current edge have not been emitted.
  //
  // REQUIRES: This method is called just after updating "inside_" for "v0".
  bool is_v0_isolated(ShapeEdgeId a_id) const {
    return !inside_ && v0_emitted_max_edge_id_ < a_id.edge_id;
  }

  // Returns true if "a_id" is the last edge of the current chain, and the
  // edges on either side of the last vertex have not been emitted (including
  // the possibility that the chain forms a loop).
  bool is_chain_last_vertex_isolated(ShapeEdgeId a_id) const {
    return (a_id.edge_id == chain_limit_ - 1 && !chain_v0_emitted_ &&
            v0_emitted_max_edge_id_ <= a_id.edge_id);
  }

  // Returns true if the given polyline edge contains "v0", taking into
  // account the specified PolylineModel.
  bool polyline_contains_v0(int edge_id, int chain_start) const {
    return (polyline_model_ != PolylineModel::OPEN || edge_id > chain_start);
  }

  void AddCrossing(SourceEdgeCrossing const& crossing) {
    source_edge_crossings_.push_back(make_pair(input_edge_id(), crossing));
  }

  void SetClippingState(InputEdgeId parameter, bool state) {
    AddCrossing(SourceEdgeCrossing(SourceId(parameter), state));
  }

  bool AddEdge(ShapeEdgeId a_id, S2Shape::Edge const& a,
               int dimension, int interior_crossings) {
    if (builder_ == nullptr) return false;  // Boolean output.
    if (interior_crossings > 0) {
      // Build a map that translates temporary edge ids (SourceId) to
      // the representation used by EdgeClippingLayer (InputEdgeId).
      SourceId src_id(a_region_id_, a_id.shape_id, a_id.edge_id);
      source_id_map_[src_id] = input_edge_id();
    }
    // Set the GraphEdgeClipper's "inside" state to match ours.
    if (inside_ != prev_inside_) SetClippingState(kSetInside, inside_);
    input_dimensions_->push_back(dimension);
    builder_->AddEdge(a.v0, a.v1);
    inside_ ^= (interior_crossings & 1);
    prev_inside_ = inside_;
    return true;
  }

  bool AddPointEdge(S2Point const& p, int dimension) {
    if (builder_ == nullptr) return false;  // Boolean output.
    if (!prev_inside_) SetClippingState(kSetInside, true);
    input_dimensions_->push_back(dimension);
    builder_->AddEdge(p, p);
    prev_inside_ = true;
    return true;
  }

  bool ProcessEdge0(ShapeEdgeId a_id, S2Shape::Edge const& a,
                    CrossingIterator* it);
  bool ProcessEdge1(ShapeEdgeId a_id, S2Shape::Edge const& a,
                    CrossingIterator* it);
  bool ProcessEdge2(ShapeEdgeId a_id, S2Shape::Edge const& a,
                    CrossingIterator* it);

  void SkipCrossings(ShapeEdgeId a_id, CrossingIterator* it);
  PointCrossingResult ProcessPointCrossings(
      ShapeEdgeId a_id, S2Point const& a0, CrossingIterator* it) const;
  EdgeCrossingResult ProcessEdgeCrossings(
      ShapeEdgeId a_id, S2Shape::Edge const& a, bool check_polygon_vertices,
      CrossingIterator* it);

  bool PolylineEdgeContainsVertex(S2Point const& v,
                                  CrossingIterator const& it) const;

  // Constructor parameters:

  PolygonModel polygon_model_;
  PolylineModel polyline_model_;

  // The output of the CrossingProcessor consists of a subset of the input
  // edges that are emitted to "builder_", and some auxiliary information
  // that allows GraphEdgeClipper to determine which segments of those input
  // edges belong to the output.  The auxiliary information consists of the
  // dimension of each input edge, and set of input edges from the other
  // region that cross each input input edge.
  S2Builder* builder_;
  vector<int8>* input_dimensions_;
  InputEdgeCrossings* input_crossings_;

  // Fields set by StartBoundary:

  int a_region_id_, b_region_id_;
  bool invert_a_, invert_b_, invert_result_;
  bool is_union_;  // True if this is a UNION operation.

  // Fields set by StartShape:

  S2Shape const* a_shape_;
  int a_dimension_;

  // Fields set by StartChain:

  int chain_id_;
  int chain_start_;
  int chain_limit_;

  // Fields updated by ProcessEdge:

  // A temporary representation of input_crossings_ that is used internally
  // until all necessary edges from *both* polygons have been emitted to the
  // S2Builder.  This field is then converted by DoneBoundaryPair() into
  // the InputEdgCrossings format expected by GraphEdgeClipper.
  //
  // The reason that we can't construct input_crossings_ directly is that it
  // uses InputEdgeIds to identify the edges from both polygons, and when we
  // are processing edges from the first polygon, InputEdgeIds have not yet
  // been assigned to the second polygon.  So instead this field identifies
  // edges from the first polygon using an InputEdgeId, and edges from the
  // second polygon using a (region_id, shape_id, edge_id) tuple (i.e., a
  // SourceId).
  //
  // All crossings are represented twice, once to indicate that an edge from
  // polygon 0 is crossed by an edge from polygon 1, and once to indicate that
  // an edge from polygon 1 is crossed by an edge from polygon 0.
  using SourceEdgeCrossings = vector<pair<InputEdgeId, SourceEdgeCrossing>>;
  SourceEdgeCrossings source_edge_crossings_;

  // A map that translates from SourceId (the (region_id, shape_id,
  // edge_id) triple that identifies an S2ShapeIndex edge) to InputEdgeId (the
  // sequentially increasing numbers assigned to input edges by S2Builder).
  using SourceIdMap = util::btree::btree_map<SourceId, InputEdgeId>;
  SourceIdMap source_id_map_;

  // Indicates whether the point being processed along the current edge chain
  // is in the polygonal interior of the opposite region, using semi-open
  // boundaries.  If "invert_b_" is true then this field is inverted.
  //
  // Equal to: b_index_.Contains(current point) ^ invert_b_
  bool inside_;

  // The value of that "inside_" would have just before the end of the
  // previous edge added to S2Builder.  This value is used to determine
  // whether the GraphEdgeClipper state needs to be updated when jumping from
  // one edge chain to another.
  bool prev_inside_;

  // The maximum edge id of any edge in the current chain whose v0 vertex has
  // already been emitted.  This is used to determine when an isolated vertex
  // needs to be emitted, e.g. when two closed polygons share only a vertex.
  int v0_emitted_max_edge_id_;

  // True if the first vertex of the current chain has been emitted.  This is
  // used when processing loops in order to determine whether the first/last
  // vertex of the loop should be emitted as an isolated vertex.
  bool chain_v0_emitted_;
};

// See documentation above.
void S2BoundaryOperation::Impl::CrossingProcessor::StartBoundary(
    int a_region_id, bool invert_a, bool invert_b, bool invert_result) {
  a_region_id_ = a_region_id;
  b_region_id_ = 1 - a_region_id;
  invert_a_ = invert_a;
  invert_b_ = invert_b;
  invert_result_ = invert_result;
  is_union_ = invert_b && invert_result;

  // Specify to GraphEdgeClipper how these edges should be clipped.
  SetClippingState(kSetReverseA, invert_a != invert_result);
  SetClippingState(kSetInvertB, invert_b);
}

// See documentation above.
inline void S2BoundaryOperation::Impl::CrossingProcessor::StartShape(
    S2Shape const* a_shape) {
  a_shape_ = a_shape;
  a_dimension_ = a_shape->dimension();
}

// See documentation above.
inline void S2BoundaryOperation::Impl::CrossingProcessor::StartChain(
    int chain_id, int chain_start, int chain_limit, bool inside) {
  chain_id_ = chain_id;
  chain_start_ = chain_start;
  chain_limit_ = chain_limit;
  inside_ = inside;
  v0_emitted_max_edge_id_ = chain_start - 1;  // No edges emitted yet.
  chain_v0_emitted_ = false;
}

// See documentation above.
bool S2BoundaryOperation::Impl::CrossingProcessor::ProcessEdge(
    ShapeEdgeId a_id, CrossingIterator* it) {
  // chain_edge() is faster than edge() when there are multiple chains.
  auto a = a_shape_->chain_edge(chain_id_, a_id.edge_id - chain_start_);
  if (a_dimension_ == 0) {
    return ProcessEdge0(a_id, a, it);
  } else if (a_dimension_ == 1) {
    return ProcessEdge1(a_id, a, it);
  } else {
    DCHECK_EQ(2, a_dimension_);
    return ProcessEdge2(a_id, a, it);
  }
}

// PointCrossingResult summarizes the relationship between a point from
// region A (the "test point") and a set of crossing edges from region B.
// For example, "has_polygon_match" indicates whether a polygon vertex from
// region B matches the test point.
struct S2BoundaryOperation::Impl::CrossingProcessor::PointCrossingResult {
  PointCrossingResult()
      : has_point_match(false), has_polyline_match(false),
        has_polygon_match(false) {
  }
  bool has_point_match;     // Matches point.
  bool has_polyline_match;  // Matches polyline vertex.
  bool has_polygon_match;   // Matches polygon vertex.
};

// Processes an edge of dimension 0 (i.e., a point) from region A.
bool S2BoundaryOperation::Impl::CrossingProcessor::ProcessEdge0(
    ShapeEdgeId a_id, S2Shape::Edge const& a, CrossingIterator* it) {
  DCHECK_EQ(a.v0, a.v1);
  // When a region is inverted, all points and polylines are discarded.
  if (invert_a_ != invert_result_) {
    SkipCrossings(a_id, it);
    return true;
  }
  PointCrossingResult r = ProcessPointCrossings(a_id, a.v0, it);

  // If force_emit > 0, the point will be emitted.  If force_emit < 0, the
  // point will be discarded.  Otherwise the current inside_ state is used.
  int force_emit = 0;
  if (polygon_model_ != PolygonModel::SEMI_OPEN && r.has_polygon_match) {
    force_emit = (polygon_model_ == PolygonModel::OPEN) ? -1 : 1;
  }
  if (r.has_polyline_match) force_emit = 1;

  // Note that the point/polyline output for UNION includes duplicates.
  if (!is_union_ && r.has_point_match) force_emit = 1;
  if (invert_b_) force_emit = -force_emit;
  if (force_emit > 0 || (force_emit == 0 && inside_)) {
    return AddPointEdge(a.v0, 0);
  }
  return true;
}

// Skip any crossings that were not needed to determine the result.
inline void S2BoundaryOperation::Impl::CrossingProcessor::SkipCrossings(
    ShapeEdgeId a_id, CrossingIterator* it) {
  while (!it->Done(a_id)) it->Next();
}

// Returns a summary of the relationship between a test point from region A and
// a set of crossing edges from region B (see PointCrossingResult).
S2BoundaryOperation::Impl::CrossingProcessor::PointCrossingResult
S2BoundaryOperation::Impl::CrossingProcessor::ProcessPointCrossings(
    ShapeEdgeId a_id, S2Point const& a0, CrossingIterator* it) const {
  PointCrossingResult r;
  for (; !it->Done(a_id); it->Next()) {
    if (it->b_dimension() == 0) {
      r.has_point_match = true;
    } else if (it->b_dimension() == 1) {
      if (PolylineEdgeContainsVertex(a0, *it)) {
        r.has_polyline_match = true;
      }
    } else {
      r.has_polygon_match = true;
    }
  }
  return r;
}

// EdgeCrossingResult summarizes the relationship between a edge from region A
// (the "test edge") and a set of crossing edges from region B.  For example,
// "has_polygon_match" indicates whether a polygon edge from region B matches
// the test edge.
struct S2BoundaryOperation::Impl::CrossingProcessor::EdgeCrossingResult {
  EdgeCrossingResult()
      : has_polyline_match(false), has_polygon_match(false),
        has_polygon_sibling(false), contains_a0(false), contains_a1(false),
        a0_crossings(0), a1_crossings(0), interior_crossings(0) {
  }
  bool has_polyline_match;   // Matches polyline edge (either direction).
  bool has_polygon_match;    // Matches polygon edge (same direction).
  bool has_polygon_sibling;  // Matches polygon edge (reverse direction).
  bool contains_a0;          // Start vertex is contained by matching vertex.
  bool contains_a1;          // End vertex is contained by matching vertex.
  int a0_crossings;          // Count of polygon crossings at start vertex.
  int a1_crossings;          // Count of polygon crossings at end vertex.
  int interior_crossings;    // Count of polygon crossings in edge interior.
};

// Processes an edge of dimension 1 (i.e., a polyline edge) from region A.
bool S2BoundaryOperation::Impl::CrossingProcessor::ProcessEdge1(
    ShapeEdgeId a_id, S2Shape::Edge const& a, CrossingIterator* it) {
  // When a region is inverted, all points and polylines are discarded.
  if (invert_a_ != invert_result_) {
    SkipCrossings(a_id, it);
    return true;
  }
  bool check_polygon_vertices = (polygon_model_ == PolygonModel::CLOSED);
  EdgeCrossingResult r = ProcessEdgeCrossings(a_id, a,
                                              check_polygon_vertices, it);
  if (polygon_model_ == PolygonModel::SEMI_OPEN) {
    r.contains_a0 |= inside_ ^ invert_b_;
  }
  inside_ ^= (r.a0_crossings & 1);

  // If force_emit > 0, the edge will be emitted.  If force_emit < 0, the
  // edge will be discarded.  Otherwise the current inside_ state is used.
  // (Note that in the SEMI_OPEN model, polygon sibling pairs cancel each
  // other and have no effect on point or edge containment.)
  int force_emit = 0;
  if (r.has_polygon_match) {
    force_emit += (polygon_model_ == PolygonModel::OPEN) ? -1 : 1;
  }
  if (r.has_polygon_sibling) {
    force_emit += (polygon_model_ == PolygonModel::CLOSED) ? 1 : -1;
  }
  // Note that the point/polyline output for UNION includes duplicates.
  if (!is_union_ && r.has_polyline_match) force_emit = 1;
  if (invert_b_) force_emit = -force_emit;

  if (force_emit != 0 && inside_ != (force_emit > 0)) {
    inside_ = (force_emit > 0);
    ++r.a1_crossings;  // Restores the correct (semi-open) inside_ state.
  }
  // If neither edge adjacent to v0 was emitted, and this polyline contains
  // v0, and the other region contains v0, then emit an isolated vertex.
  if (is_v0_isolated(a_id) &&
      polyline_contains_v0(a_id.edge_id, chain_start_) &&
      r.contains_a0 ^ invert_b_) {
    if (!AddPointEdge(a.v0, 1)) return false;
  }
  if (inside_ || r.interior_crossings > 0) {
    if (!AddEdge(a_id, a, 1 /*dimension*/, r.interior_crossings)) {
      return false;
    }
  }
  if (inside_) v0_emitted_max_edge_id_ = a_id.edge_id + 1;
  inside_ ^= (r.a1_crossings & 1);
  DCHECK_EQ(S2ContainsPointQuery(&it->b_index()).Contains(a.v1) ^ invert_b_,
            inside_);
  if (polygon_model_ == PolygonModel::SEMI_OPEN) {
    r.contains_a1 |= inside_ ^ invert_b_;
  }
  if (polyline_model_ == PolylineModel::CLOSED &&
      is_chain_last_vertex_isolated(a_id) && r.contains_a1 ^ invert_b_) {
    if (!AddPointEdge(a.v1, 1)) return false;
  }
  return true;
}

// Processes an edge of dimension 2 (i.e., a polygon edge) from region A.
bool S2BoundaryOperation::Impl::CrossingProcessor::ProcessEdge2(
    ShapeEdgeId a_id, S2Shape::Edge const& a, CrossingIterator* it) {
  // In order to keep only one copy of any shared polygon edges, we only
  // output shared edges when processing the second region.
  bool emit_shared = (a_region_id_ == 1);

  // Degeneracies such as isolated vertices and sibling pairs can only be
  // created by intersecting CLOSED polygons or unioning OPEN polygons.
  bool emit_degenerate =
      (polygon_model_ == PolygonModel::CLOSED && !invert_a_ && !invert_b_) ||
      (polygon_model_ == PolygonModel::OPEN && invert_a_ && invert_b_);

  EdgeCrossingResult r = ProcessEdgeCrossings(a_id, a, emit_degenerate, it);
  DCHECK(!r.has_polyline_match);
  inside_ ^= (r.a0_crossings & 1);

  // In order to keep only one copy of any shared polygon edges, we only
  // output shared edges when processing the second region.
  // TODO(ericv): Update this code to handle degenerate loops.
  DCHECK(!r.has_polygon_match || !r.has_polygon_sibling);
  int force_emit = 0;
  if (invert_a_ != invert_b_) swap(r.has_polygon_match, r.has_polygon_sibling);
  if (r.has_polygon_match) force_emit = emit_shared ? 1 : -1;
  if (r.has_polygon_sibling) force_emit = emit_degenerate ? 1 : -1;

  if (force_emit != 0 && inside_ != (force_emit > 0)) {
    inside_ = (force_emit > 0);
    ++r.a1_crossings;  // Restores the correct (semi-open) inside_ state.
  }
  if (a_id.edge_id == chain_start_) {
    chain_v0_emitted_ = inside_;
  } else if (emit_shared && r.contains_a0 && is_v0_isolated(a_id)) {
    if (!AddPointEdge(a.v0, 2)) return false;
  }
  if (inside_ || r.interior_crossings > 0) {
    if (!AddEdge(a_id, a, 2 /*dimension*/, r.interior_crossings)) {
      return false;
    }
  }
  if (inside_) v0_emitted_max_edge_id_ = a_id.edge_id + 1;
  inside_ ^= (r.a1_crossings & 1);
  DCHECK_EQ(S2ContainsPointQuery(&it->b_index()).Contains(a.v1) ^ invert_b_,
            inside_);
  if (emit_shared && r.contains_a1 && is_chain_last_vertex_isolated(a_id)) {
    if (!AddPointEdge(a.v1, 2)) return false;
  }
  return true;
}

// Returns a summary of the relationship between a test edge from region A and
// a set of crossing edges from region B (see EdgeCrossingResult).
// "check_polygon_vertices" specifies whether polygon vertices should be
// considered to contain the endpoints of the test edge.
S2BoundaryOperation::Impl::CrossingProcessor::EdgeCrossingResult
S2BoundaryOperation::Impl::CrossingProcessor::ProcessEdgeCrossings(
    ShapeEdgeId a_id, S2Shape::Edge const& a, bool check_polygon_vertices,
    CrossingIterator* it) {
  EdgeCrossingResult r;
  if (it->Done(a_id)) return r;

  // TODO(ericv): bool a_degenerate = (a.v0 == a.v1);
  S2CopyingEdgeCrosser crosser(a.v0, a.v1);
  for (; !it->Done(a_id); it->Next()) {
    if (it->b_dimension() == 0) continue;
    S2Shape::Edge b = it->b_edge();
    int crossing = crosser.CrossingSign(b.v0, b.v1);
    if (crossing < 0) continue;
    if (crossing > 0) {
      // The crossing occurs in the edge interior.  The condition below says
      // that (1) polyline crossings don't affect polygon output, and (2)
      // subtracting a crossing polyline from a polyline has no effect.
      if (a_dimension_ <= it->b_dimension() &&
          !(invert_b_ != invert_result_ && it->b_dimension() == 1)) {
        SourceId src_id(b_region_id_, it->b_shape_id(), it->b_edge_id());
        bool left_to_right = (s2pred::Sign(a.v0, a.v1, b.v0) > 0);
        AddCrossing(make_pair(src_id, left_to_right));
      }
      r.interior_crossings += (it->b_dimension() == 1) ? 2 : 1;
    } else if (it->b_dimension() == 1) {
      if (a_dimension_ == 2) continue;
      if ((a.v0 == b.v0 && a.v1 == b.v1) || (a.v0 == b.v1 && a.v1 == b.v0)) {
        r.has_polyline_match = true;
      }
      if (!is_union_) {
        if ((a.v0 == b.v0 || a.v0 == b.v1) &&
            PolylineEdgeContainsVertex(a.v0, *it)) {
          r.contains_a0 = true;
        }
        if ((a.v1 == b.v0 || a.v1 == b.v1) &&
            PolylineEdgeContainsVertex(a.v1, *it)) {
          r.contains_a1 = true;
        }
      }
    } else {
      if (S2::VertexCrossing(a.v0, a.v1, b.v0, b.v1)) {
        if (a.v0 == b.v0 || a.v0 == b.v1) {
          ++r.a0_crossings;
          if (a.v0 == b.v0 && a.v1 == b.v1) {
            r.has_polygon_match = true;
          } else if (a.v0 == b.v1 && a.v1 == b.v0) {
            r.has_polygon_sibling = true;
          }
        } else {
          ++r.a1_crossings;
        }
      }
      if (check_polygon_vertices) {
        if (a.v0 == b.v0 || a.v0 == b.v1) {
          r.contains_a0 = true;
        } else {
          r.contains_a1 = true;
        }
      }
    }
  }
  return r;
}

// Returns true if the vertex "v" is contained by the polyline edge referred
// to by the CrossingIterator "it", taking into account the PolylineModel.
//
// REQUIRES: it.b_dimension() == 1
// REQUIRES: "v" is an endpoint of it.b_edge()
bool S2BoundaryOperation::Impl::CrossingProcessor::PolylineEdgeContainsVertex(
    S2Point const& v, CrossingIterator const& it) const {
  DCHECK_EQ(1, it.b_dimension());
  DCHECK(it.b_edge().v0 == v || it.b_edge().v1 == v);

  // Closed polylines contain all their vertices.
  if (polyline_model_ == PolylineModel::CLOSED) return true;

  // All interior polyline vertices are contained.
  // The last polyline vertex is contained iff the model is CLOSED.
  // The first polyline vertex is contained iff the model is not OPEN.
  // The test below is written such that usually b_edge() is not needed.
  auto const& b_chain = it.b_chain_info();
  int b_edge_id = it.b_edge_id();
  return (b_edge_id < b_chain.limit - 1 || it.b_edge().v1 != v) &&
      (polyline_contains_v0(b_edge_id, b_chain.start) || it.b_edge().v0 != v);
}

// Translates the temporary representation of crossing edges (SourceId) into
// the format expected by EdgeClippingLayer (InputEdgeId).
void S2BoundaryOperation::Impl::CrossingProcessor::DoneBoundaryPair() {
  // Add entries that translate the "special" crossings.
  source_id_map_[SourceId(kSetInside)] = kSetInside;
  source_id_map_[SourceId(kSetInvertB)] = kSetInvertB;
  source_id_map_[SourceId(kSetReverseA)] = kSetReverseA;
  input_crossings_->reserve(input_crossings_->size() +
                            source_edge_crossings_.size());
  for (auto const& tmp : source_edge_crossings_) {
    auto it = source_id_map_.find(tmp.second.first);
    DCHECK(it != source_id_map_.end());
    input_crossings_->push_back(make_pair(
        tmp.first, CrossingInputEdge(it->second, tmp.second.second)));
  }
  source_edge_crossings_.clear();
  source_id_map_.clear();
}

// Clips the boundary of A to the interior of the opposite region B and adds
// the resulting edges to the output.  Optionally, any combination of region
// A, region B, and the result may be inverted, which allows operations such
// as union and difference to be implemented.
//
// Note that when an input region is inverted with respect to the output
// (e.g., invert_a != invert_result), all polygon edges are reversed and all
// points and polylines are discarded, since the complement of such objects
// cannot be represented.  (If you want to compute the complement of points
// or polylines, you can use LaxPolygon to represent your geometry as
// degenerate polygons instead.)
//
// This method must be called an even number of times (first to clip A to B
// and then to clip B to A), calling DoneBoundaryPair() after each pair.
bool S2BoundaryOperation::Impl::AddBoundary(
    int a_region_id, bool invert_a, bool invert_b, bool invert_result,
    vector<ShapeEdgeId> const& a_chain_starts, CrossingProcessor* cp) {
  S2ShapeIndex const& a_index = *op_->regions_[a_region_id];
  S2ShapeIndex const& b_index = *op_->regions_[1 - a_region_id];
  if (!GetIndexCrossings(a_region_id)) return false;
  cp->StartBoundary(a_region_id, invert_a, invert_b, invert_result);

  // Walk the boundary of region A and build a list of all edge crossings.
  // We also keep track of whether the current vertex is inside region B.
  auto next_start = a_chain_starts.begin();
  CrossingIterator next_crossing(&b_index, index_crossings_);
  ShapeEdgeId next_id = min(*next_start, next_crossing.a_id());
  while (next_id != kSentinel) {
    int a_shape_id = next_id.shape_id;
    S2Shape const& a_shape = *a_index.shape(a_shape_id);
    cp->StartShape(&a_shape);
    while (next_id.shape_id == a_shape_id) {
      // TODO(ericv): Special handling of dimension 0?  Can omit most of this
      // code, including the loop, since all chains are of length 1.
      int edge_id = next_id.edge_id;
      S2Shape::ChainPosition chain_position = a_shape.chain_position(edge_id);
      int chain_id = chain_position.chain_id;
      int chain_start = edge_id - chain_position.offset;
      int chain_limit = chain_start + a_shape.chain(chain_id).length;
      bool start_inside = (next_id == *next_start);
      if (start_inside) ++next_start;
      cp->StartChain(chain_id, chain_start, chain_limit, start_inside);
      while (edge_id < chain_limit) {
        ShapeEdgeId a_id(a_shape_id, edge_id);
        DCHECK(cp->inside() || next_crossing.a_id() == a_id);
        if (!cp->ProcessEdge(a_id, &next_crossing)) {
          return false;
        }
        if (cp->inside()) {
          ++edge_id;
        } else if (next_crossing.a_id().shape_id == a_shape_id &&
                   next_crossing.a_id().edge_id < chain_limit) {
          edge_id = next_crossing.a_id().edge_id;
        } else {
          break;
        }
      }
      next_id = min(*next_start, next_crossing.a_id());
    }
  }
  return true;
}

bool S2BoundaryOperation::Impl::GetChainStarts(
    int a_region_id, bool invert_a, bool invert_b, bool invert_result,
    vector<ShapeEdgeId>* chain_starts) const {
  S2ShapeIndex const& a_index = *op_->regions_[a_region_id];
  S2ShapeIndex const& b_index = *op_->regions_[1 - a_region_id];

  // Fast path for the case where region B has no two-dimensional shapes.
  bool b_has_interior = HasInterior(b_index);
  if (b_has_interior || invert_b) {
    S2ContainsPointQuery query(&b_index);
    int num_shape_ids = a_index.num_shape_ids();
    for (int s = 0; s < num_shape_ids; ++s) {
      S2Shape* a_shape = a_index.shape(s);
      if (a_shape == nullptr) continue;

      // When points or polylines are subtracted from another region, they
      // never belong to the output (they can only remove edges).  Therefore
      // we don't need to process these edge chains.
      if (invert_a != invert_result && a_shape->dimension() < 2) continue;
      int num_chains = a_shape->num_chains();
      for (int chain = 0; chain < num_chains; ++chain) {
        auto a = a_shape->chain_edge(chain, 0);
        // TODO(ericv): Optimization for boolean output: if a chain start is
        // contained by an empty cell, or more generally if it does not match
        // any B vertex, then we can exit early and return "false".
        if (!b_has_interior || query.Contains(a.v0) != invert_b) {
          int edge_id = a_shape->chain(chain).start;
          chain_starts->push_back(ShapeEdgeId(s, edge_id));
        }
      }
    }
  }
  chain_starts->push_back(kSentinel);
  return true;
}

#if 0
// Compute whether the start of each edge chain is "inside" with respect to
// the *polygons* of the other region (only).  If a region has no polygons,
// then all chain starts are inside/outside depending on the operation (e.g.,
// for union all chain starts are "inside").
//
// If we only want to know whether the result is empty, there are other
// optimizations.  For intersection, union, and difference we can define a
// boolean value for each region such that if there is a point in each region
// that matches the given containment status then the result is non-empty.
// Using the center of each region is always fine.  However it is not enough
// if one cell is empty and the other is non-empty, because cells are expanded
// slightly.  Therefore the non-empty part might be outside the cell.  However
// it is also OK if a point (such as a chain start) is inside the cell
// according to S2Cell::Contains() and the other cell is empty.
//
// For symmetric difference we can either call "difference" twice, or look for
// points that are known to have different containment statuses.
//
// Returns false if the client only wants to determine whether the result is
// empty, and the result is known to be non-empty.
bool S2BoundaryOperation::Impl::FindChainStarts() {
  S2ShapeIndex const& a_index = *op_->regions_[0];
  S2ShapeIndex const& b_index = *op_->regions_[1];
  s2shapeutil::RangeIterator ai(a_index), bi(b_index);
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

#endif

bool S2BoundaryOperation::Impl::HasInterior(S2ShapeIndex const& index) {
  for (int s = index.num_shape_ids(); --s >= 0; ) {
    S2Shape* shape = index.shape(s);
    if (shape && shape->has_interior()) return true;
  }
  return false;
}

// Initialize index_crossings_ to the set of crossing edge pairs such that the
// first element of each pair is an edge from "region_id".
bool S2BoundaryOperation::Impl::GetIndexCrossings(int region_id) {
  if (region_id == index_crossings_first_region_id_) return true;
  if (index_crossings_first_region_id_ < 0) {
    DCHECK_EQ(region_id, 0);  // For efficiency, not correctness.
    if (!s2shapeutil::VisitCrossings(
            *op_->regions_[0], *op_->regions_[1],
            s2shapeutil::CrossingType::ALL,
            [this](s2shapeutil::ShapeEdge const& a,
                   s2shapeutil::ShapeEdge const& b, bool is_interior) {
              if (is_interior && is_boolean_output()) return false;
              index_crossings_.push_back(make_pair(a.id(), b.id()));
              return true;  // Continue visiting.
            })) {
      return false;
    }
    if (index_crossings_.size() > 1) {
      std::sort(index_crossings_.begin(), index_crossings_.end());
      index_crossings_.erase(
          std::unique(index_crossings_.begin(), index_crossings_.end()),
          index_crossings_.end());
    }
    // Add a sentinel value to simplify the loop logic.
    index_crossings_.push_back(make_pair(kSentinel, kSentinel));
    index_crossings_first_region_id_ = 0;
  }
  if (region_id != index_crossings_first_region_id_) {
    for (auto& crossing : index_crossings_) {
      swap(crossing.first, crossing.second);
    }
    std::sort(index_crossings_.begin(), index_crossings_.end());
    index_crossings_first_region_id_ = region_id;
  }
  return true;
}

bool S2BoundaryOperation::Impl::AddBoundaryPair(
    bool invert_a, bool invert_b, bool invert_result, CrossingProcessor* cp) {
  // TODO(ericv): When computing a boolean result, there are other quick
  // checks that could be done here.
  vector<ShapeEdgeId> a_starts, b_starts;
  if (!GetChainStarts(0, invert_a, invert_b, invert_result, &a_starts) ||
      !GetChainStarts(1, invert_b, invert_a, invert_result, &b_starts) ||
      !AddBoundary(0, invert_a, invert_b, invert_result, a_starts, cp) ||
      !AddBoundary(1, invert_b, invert_a, invert_result, b_starts, cp)) {
    return false;
  }
  cp->DoneBoundaryPair();
  return true;
}

bool S2BoundaryOperation::Impl::BuildOpType(OpType op_type) {
  // CrossingProcessor does the real work of emitting the output edges.
  CrossingProcessor cp(op_->options_.polygon_model(),
                       op_->options_.polyline_model(),
                       builder_.get(), &input_dimensions_, &input_crossings_);
  switch (op_type) {
    case OpType::UNION:
      // A | B == ~(~A & ~B)
      return AddBoundaryPair(true, true, true, &cp);

    case OpType::INTERSECTION:
      // A & B
      return AddBoundaryPair(false, false, false, &cp);

    case OpType::DIFFERENCE:
      // A - B = A & ~B
      return AddBoundaryPair(false, true, false, &cp);

    case OpType::SYMMETRIC_DIFFERENCE:
      // Compute the union of (A - B) and (B - A).
      return (AddBoundaryPair(false, true, false, &cp) &&
              AddBoundaryPair(true, false, false, &cp));
  }
  LOG(FATAL) << "Invalid S2BoundaryOperation::OpType";
  return false;
}

bool S2BoundaryOperation::Impl::Build(S2Error* error) {
  error->Clear();
  if (is_boolean_output()) {
    *op_->result_non_empty_ = BuildOpType(op_->op_type());
    return true;
  }
  // TODO(ericv): Rather than having S2Builder split the edges, it would be
  // faster to call AddVertex() in this class and have a new S2Builder
  // option that increases the edge_snap_radius_ to account for errors in
  // the intersection point (the way that split_crossing_edges does).
  S2Builder::Options options(op_->options_.snap_function());
  options.set_split_crossing_edges(true);

  // TODO(ericv): Ideally idempotent() should be true, but existing clients
  // expect vertices closer than the full "snap_radius" to be snapped.
  options.set_idempotent(false);
  builder_ = MakeUnique<S2Builder>(options);
  builder_->StartLayer(MakeUnique<EdgeClippingLayer>(
      &op_->layers_, &input_dimensions_, &input_crossings_));
  (void) BuildOpType(op_->op_type());
  return builder_->Build(error);
}

S2BoundaryOperation::Options::Options()
    : Options(s2builderutil::IdentitySnapFunction(S1Angle::Zero())) {
}

S2BoundaryOperation::Options::Options(SnapFunction const& snap_function)
    : snap_function_(snap_function.Clone()),
      polygon_model_(PolygonModel::SEMI_OPEN),
      polyline_model_(PolylineModel::CLOSED) {
}

S2BoundaryOperation::Options::Options(Options const& options)
    :  snap_function_(options.snap_function_->Clone()),
       polygon_model_(options.polygon_model_),
       polyline_model_(options.polyline_model_) {
}

S2BoundaryOperation::Options& S2BoundaryOperation::Options::operator=(
    Options const& options) {
  snap_function_ = options.snap_function_->Clone();
  polygon_model_ = options.polygon_model_;
  polyline_model_ = options.polyline_model_;
  return *this;
}

SnapFunction const& S2BoundaryOperation::Options::snap_function() const {
  return *snap_function_;
}

void S2BoundaryOperation::Options::set_snap_function(
    SnapFunction const& snap_function) {
  snap_function_ = snap_function.Clone();
}

PolygonModel S2BoundaryOperation::Options::polygon_model() const {
  return polygon_model_;
}

void S2BoundaryOperation::Options::set_polygon_model(PolygonModel model) {
  polygon_model_ = model;
}

PolylineModel S2BoundaryOperation::Options::polyline_model() const {
  return polyline_model_;
}

void S2BoundaryOperation::Options::set_polyline_model(PolylineModel model) {
  polyline_model_ = model;
}

char const* S2BoundaryOperation::OpTypeToString(OpType op_type) {
  switch (op_type) {
    case OpType::UNION:                return "UNION";
    case OpType::INTERSECTION:         return "INTERSECTION";
    case OpType::DIFFERENCE:           return "DIFFERENCE";
    case OpType::SYMMETRIC_DIFFERENCE: return "SYMMETRIC DIFFERENCE";
    default:                           return "Unknown OpType";
  }
}

S2BoundaryOperation::S2BoundaryOperation(OpType op_type,
                                         Options const& options)
    : op_type_(op_type), options_(options), result_non_empty_(nullptr) {
}

S2BoundaryOperation::S2BoundaryOperation(OpType op_type, bool* result_non_empty,
                                         Options const& options)
    : op_type_(op_type), options_(options),
      result_non_empty_(result_non_empty) {
}

S2BoundaryOperation::S2BoundaryOperation(
    OpType op_type, unique_ptr<S2Builder::Layer> layer, Options const& options)
    : op_type_(op_type), options_(options), result_non_empty_(nullptr) {
  layers_.push_back(std::move(layer));
}

S2BoundaryOperation::S2BoundaryOperation(
    OpType op_type, vector<unique_ptr<S2Builder::Layer>> layers,
    Options const& options)
    : op_type_(op_type), options_(options), layers_(std::move(layers)),
      result_non_empty_(nullptr) {
}

bool S2BoundaryOperation::Build(S2ShapeIndex const& a, S2ShapeIndex const& b,
                                S2Error* error) {
  regions_[0] = &a;
  regions_[1] = &b;
  return Impl(this).Build(error);
}
