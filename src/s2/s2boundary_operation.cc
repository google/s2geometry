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
// small tolerance (snap_radius + S2EdgeUtil::kIntersectionError) of the exact
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
//  - CrossingQuery: a wrapper around S2CrossingEdgeQuery that finds all the
//                   input edges that cross a given query edge.
//
//  - EdgeClippingLayer: an S2Builder::Layer that removes graph edges that
//                       correspond to clipped portions of input edges, and
//                       passes the result to another layer for assembly.
//
//  - GraphEdgeClipper: a helper class that does the actual work of the
//                      EdgeClippingLayer.

#include "s2/s2boundary_operation.h"

#include <memory>

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
using SourceId = S2BoundaryOperation::SourceId;

// SourceEdgeCrossing represents an input edge that crosses some other edge;
// it crosses the edge from left to right iff the second parameter is "true".
using SourceEdgeCrossing = pair<SourceId, bool>;

// CrossingInputEdge represents an input edge B that crosses some other input
// edge A.  It stores the input edge id of edge B and also whether it crosses
// edge A from left to right (or vice versa).
//
// It also defines a special constant "kAtStartVertex" that represents an
// artificial extra crossing at the starting vertex of edge A.  This is used
// to toggle the inside/outside state of GraphEdgeClipper as necessary when
// making a discontinuous jump from one edge chain to another.
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

  // A special value that indicates there is a crossing at the starting vertex
  // of the edge being crossed.
  static CrossingInputEdge kAtStartVertex() {
    return CrossingInputEdge(-1, false);
  }
  bool at_start_vertex() const { return input_id_ < 0; }

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
using CrossingGraphEdgeVector = gtl::InlinedVector<CrossingGraphEdge, 2>;

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

// Given a set of clipping instructions encoded as a set of intersections
// between input edges, GraphEdgeClipper determines which graph edges
// correspond to clipped portions of input edges and removes them.
//
// The clipping model is as follows.  The input consists of edge chains.  The
// clipper maintains an "inside" boolean state as it clips each chain, and
// toggles this state whenever an input edge is crossed.  Any edges that are
// deemed to be "outside" after clipping are removed.  The "inside" state can
// be altered between chains when necessary by inserting an artificial extra
// crossing at the start of any edge (see CrossingEdge::kAtStartVertex).
//
// "clip_to_exterior" indicates whether the input edges of one region were
// originally clipped to the interior or exterior of the other region.
class GraphEdgeClipper {
 public:
  // "input_crossings" is the set of all crossings to be used when clipping
  // the edges of "g", sorted in lexicographic order.
  GraphEdgeClipper(Graph const& g, InputEdgeCrossings const* input_crossings,
                   bool clip_to_exterior);

  // Returns the clipped set of edges and their corresponding set of input
  // edge ids.  (This is used to construct a new S2Builder::Graph.)
  void Run(vector<Graph::Edge>* new_edges,
           vector<InputEdgeIdSetId>* new_input_edge_ids) const;

 private:
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
  InputEdgeCrossings const* input_crossings_;
  bool clip_to_exterior_;

  // Every graph edge is associated with exactly one input edge in our case,
  // which means that we can declare g_.input_edge_id_set_ids() as a vector of
  // InputEdgeIds rather than a vector of InputEdgeIdSetIds.  (This also takes
  // advantage of the fact that IdSetLexicon represents a singleton set as the
  // value of its single element.)
  vector<InputEdgeId> const& input_ids_;

  vector<EdgeId> order_;  // Graph edges sorted in input edge id order.
  vector<int> rank_;      // The rank of each graph edge within order_.
};

GraphEdgeClipper::GraphEdgeClipper(Graph const& g,
                                   InputEdgeCrossings const* input_crossings,
                                   bool clip_to_exterior)
    : g_(g), in_(g), out_(g), input_crossings_(input_crossings),
      clip_to_exterior_(clip_to_exterior),
      input_ids_(g.input_edge_id_set_ids()),
      order_(GetInputEdgeChainOrder(g_, input_ids_)),
      rank_(order_.size()) {
  for (int i = 0; i < order_.size(); ++i) {
    rank_[order_[i]] = i;
  }
}

void GraphEdgeClipper::Run(vector<Graph::Edge>* new_edges,
                           vector<InputEdgeIdSetId>* new_input_edge_ids) const {
  // Declare vectors here and reuse them to avoid reallocation.
  vector<VertexId> a_vertices;
  vector<int> a_num_crossings;
  vector<CrossingInputEdge> b_input_edges;
  vector<CrossingGraphEdgeVector> b_edges;

  bool inside = false;
  auto next = input_crossings_->begin();
  for (int i = 0; i < order_.size(); ++i) {
    // For each input edge (the "A" input edge), gather all the input edges
    // that cross it (the "B" input edges).
    InputEdgeId a_input_id = input_ids_[order_[i]];
    Graph::Edge const& edge0 = g_.edge(order_[i]);
    b_input_edges.clear();
    for (; next != input_crossings_->end(); ++next) {
      if (next->first != a_input_id) break;
      if (next->second.at_start_vertex()) {
        inside = !inside;
      } else {
        b_input_edges.push_back(next->second);
      }
    }
    // Optimization for degenerate edges.
    if (edge0.first == edge0.second) {
      inside ^= (b_input_edges.size() & 1);
      continue;
    }
    // Optimization for the case where there are no crossings.
    if (b_input_edges.empty()) {
      // In general the caller only passes edges that are part of the output
      // (i.e., we could DCHECK(inside) here).  The one exception is for
      // polyline/polygon operations, where the polygon edges are needed to
      // compute the polyline output but are not emitted themselves.
      if (inside) {
        new_edges->push_back(edge0);
        new_input_edge_ids->push_back(a_input_id);
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
      a_num_crossings[a_index] += (left_to_right == clip_to_exterior_) ? -1 : 1;
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
      // Output any forward edges required.
      for (int i = 0; i < multiplicity; ++i) {
        new_edges->push_back(Graph::Edge(a_vertices[ai - 1], a_vertices[ai]));
        new_input_edge_ids->push_back(a_input_id);
      }
      // Output any reverse edges required.
      for (int i = multiplicity; i < 0; ++i) {
        new_edges->push_back(Graph::Edge(a_vertices[ai], a_vertices[ai - 1]));
        new_input_edge_ids->push_back(a_input_id);
      }
      multiplicity += a_num_crossings[ai];
    }
    // Multiplicities other than 0 or 1 can only occur in the edge interior.
    DCHECK(multiplicity == 0 || multiplicity == 1);
    inside = (multiplicity != 0);
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
  // "clip_to_exterior" indicates whether the input edges of one region were
  // originally clipped to the interior or exterior of the other region.
  EdgeClippingLayer(unique_ptr<S2Builder::Layer> output_layer,
                    InputEdgeCrossings const* input_crossings,
                    bool clip_to_exterior)
      : output_layer_(std::move(output_layer)),
        input_crossings_(input_crossings),
        clip_to_exterior_(clip_to_exterior) {
  }

  // Layer interface:
  GraphOptions graph_options() const override;
  void Build(Graph const& g, S2Error* error) override;

 private:
  unique_ptr<S2Builder::Layer> output_layer_;
  InputEdgeCrossings const* input_crossings_;
  bool clip_to_exterior_;
};

GraphOptions EdgeClippingLayer::graph_options() const {
  // We keep all edges, including degenerate ones, so that we can figure out
  // the correspondence between input edge crossings and output edge
  // crossings.
  GraphOptions graph_options;
  graph_options.set_edge_type(EdgeType::DIRECTED);
  graph_options.set_degenerate_edges(DegenerateEdges::KEEP);
  graph_options.set_duplicate_edges(DuplicateEdges::KEEP);
  graph_options.set_sibling_pairs(SiblingPairs::KEEP);
  return graph_options;
}

void EdgeClippingLayer::Build(Graph const& g, S2Error* error) {
  // The bulk of the work is handled by GraphEdgeClipper.
  vector<Graph::Edge> new_edges;
  vector<InputEdgeIdSetId> new_input_edge_ids;
  GraphEdgeClipper clipper(g, input_crossings_, clip_to_exterior_);
  clipper.Run(&new_edges, &new_input_edge_ids);
  if (s2builder_verbose) {
    std::cout << "Edges after clipping: " << std::endl;
    for (int i = 0; i < new_edges.size(); ++i) {
      std::cout << "  " << new_input_edge_ids[i] << " (" << new_edges[i].first
                << ", " << new_edges[i].second << ")" << std::endl;
    }
  }
  // Construct a new graph from the clipped edges and pass it to the given
  // output layer.
  IdSetLexicon new_input_edge_id_set_lexicon;
  GraphOptions options = output_layer_->graph_options();
  Graph::ProcessEdges(&options, &new_edges, &new_input_edge_ids,
                      &new_input_edge_id_set_lexicon, error);
  Graph new_graph(options, g.vertices(), new_edges, new_input_edge_ids,
                  new_input_edge_id_set_lexicon, g.label_set_ids(),
                  g.label_set_lexicon());
  output_layer_->Build(new_graph, error);
}

}  // namespace

// CrossingQuery is a helper class that finds all the places where an edge
// from one region crosses the edges of the other region.
//
// TODO(ericv): Restructure this code so that it first finds all the pairs of
// intersecting edges, and then sorts the intersections to generate the
// crossings for each shape without visiting every edge.  This would allow
// boundary operations to be done in sublinear time in many cases (e.g.,
// can be found by traversing the region with fewer edges and querying the
// other region, or ideally by traversing the two S2ShapeIndexes in parallel
// (similar to S2Loop::HasCrossingRelation).
class S2BoundaryOperation::CrossingQuery {
 public:
  // Prepares to query the given region ("b").  "reverse_b" specifies whether
  // edges should be reversed before testing whether an edge is shared.
  // "keep_shared_edges" specifies whether a query edge that exactly matches a
  // polygon edge should be kept or discarded.  "clip_to_exterior" specifies
  // whether the query edges are clipped to the interior or exterior of "b".
  CrossingQuery(Region const& b, bool reverse_b,
                bool keep_shared_edges, bool clip_to_exterior)
      : b_(b), reverse_b_(reverse_b), b_query_(b.index),
        keep_shared_edges_(keep_shared_edges),
        clip_to_exterior_(clip_to_exterior) {
  }

  // Find all intersections of the edge (a0, a1) with edges in the index.
  // "a0_crossings" and "a1_crossings" are set to the number of crossings at
  // a0 and a1 respectively, while crossings that occur in the edge interior
  // are returned in "edge_crossings".
  void GetCrossings(S2Point const& a0, S2Point const& a1,
                    int* a0_crossings, int* a1_crossings,
                    vector<SourceEdgeCrossing>* edge_crossings,
                    bool* has_polyline_crossings);

 private:
  Region b_;
  bool reverse_b_;
  S2CrossingEdgeQuery b_query_;
  bool keep_shared_edges_;
  bool clip_to_exterior_;
  S2CrossingEdgeQuery::EdgeMap edge_map_;  // Temporary
};

void S2BoundaryOperation::CrossingQuery::GetCrossings(
    S2Point const& a0, S2Point const& a1, int* a0_crossings, int* a1_crossings,
    vector<SourceEdgeCrossing>* edge_crossings, bool* has_polyline_crossings) {
  *a0_crossings = 0;
  *a1_crossings = 0;
  edge_crossings->clear();
  *has_polyline_crossings = false;
  if (!b_query_.GetCandidates(a0, a1, &edge_map_)) return;

  // Iterate through the candidate loops, and then the candidate edges within
  // each loop.
  S2EdgeUtil::EdgeCrosser crosser(&a0, &a1);
  for (const auto& p : edge_map_) {
    S2Shape const& b_shape = *p.first;
    int b_dim = b_shape.dimension();
    if (b_dim == 0) continue;

    vector<int> const& b_candidates = p.second;
    int n = b_candidates.size();
    for (int j = 0; j < n; ++j) {
      S2Point const *b0, *b1;
      b_shape.GetEdge(b_candidates[j], &b0, &b1);
      int crossing = crosser.CrossingSign(b0, b1);
      if (b_dim == 1) {
        *has_polyline_crossings |= crossing > 0;
        continue;
      }
      if (crossing < 0) continue;

      // Check whether this is an edge or vertex crossing.
      if (b_.reverse_edges(b_shape) ^ reverse_b_) swap(b0, b1);
      if (crossing > 0) {
        // The crossing occurs in the edge interior.
        SourceId src_id(b_.id, b_shape.id(), b_candidates[j]);
        bool left_to_right = (s2pred::Sign(a0, a1, *b0) > 0);
        edge_crossings->push_back(make_pair(src_id, left_to_right));
      } else if (S2EdgeUtil::VertexCrossing(a0, a1, *b0, *b1)) {
        // There is a crossing at one of the vertices.  If only one vertex is
        // shared between the two edges, then the crossing obviously occurs at
        // that vertex, but if both vertices are shared (i.e., the edges are
        // identical or reversed) then we need to decide at which vertex the
        // crossing occurs (a0 or a1).
        //
        // If "clip_to_exterior" is false, it turns out that choosing a0 for
        // the crossing will keep shared edges and discard reversed edges.
        // The opposite is true if we choose a1 for the crossing.
        //
        // If "clip_to_exterior" is true, the cases above are reversed (i.e.,
        // choosing a1 will keep shared edges and discard reversed edges).
        if (a0 == *b0 && a1 == *b1) {
          // This is a shared edge.  To keep such edges, we choose the
          // crossing at a0 (or a1 if "clip_to_exterior" is true).
          if (keep_shared_edges_ ^ clip_to_exterior_) {
            ++*a0_crossings;
          } else {
            ++*a1_crossings;
          }
        } else if (a0 == *b1 && a1 == *b0) {
          // This is a reversed edge.  Such edges are discarded by choosing
          // the crossing at a0 (or a1 if "clip_to_exterior" is true).
          if (!clip_to_exterior_) {
            ++*a0_crossings;
          } else {
            ++*a1_crossings;
          }
        } else if (a0 == *b0 || a0 == *b1) {
          ++*a0_crossings;  // The two edges share only one vertex (a0).
        } else {
          ++*a1_crossings;  // The two edges share only one vertex (a1).
        }
      }
    }
  }
}

class S2BoundaryOperation::Impl {
 public:
  Impl(S2BoundaryOperation* op) : op_(op) {}

  bool Build(S2Error* error);

 private:
  // Returns a special constant indicating that there is an artificial
  // "crossing" at the starting vertex of a snapped edge chain.
  SourceEdgeCrossing kAtStartVertex() const {
    return SourceEdgeCrossing(SourceId(0, 0, -1), false);
  }
  void Start(bool clip_to_exterior);
  void AddBoundary(int a_region_id, bool reverse_a,
                   int b_region_id, bool reverse_b);
  void TranslateCrossings();

  S2BoundaryOperation* op_;

  // "clip_to_exterior" specifies whether the edges of one region are clipped
  // to the interior or exterior of the other region.
  bool clip_to_exterior_;

  // The S2Builder used to construct the output.
  S2Builder builder_;

  // Stores the set of all edge intersections between the two input polygons.
  // Each intersection is stored twice, once to indicate that an edge from
  // polygon 0 is crossed by an edge from polygon 1, and once to indicate that
  // an edge from polygon 1 is crossed by an edge from polygon 0.
  //
  // Because input edge ids for the second polygon are not assigned until it
  // is processed, edges from the two polygons are identified in different
  // ways: the first polygon uses the input edge id, while the second polygon
  // uses a (loop index, loop edge index) pair (SourceId).
  using SourceEdgeCrossings = vector<pair<InputEdgeId, SourceEdgeCrossing>>;
  SourceEdgeCrossings source_edge_crossings_;

  // A map that translates from SourceId (the (region_id, shape_id,
  // edge_id) triple that identifies an S2ShapeIndex edge) to InputEdgeId (the
  // sequentially increasing numbers assigned to input edges by S2Builder).
  using SourceEdgeMap = util::btree::btree_map<SourceId, InputEdgeId>;
  SourceEdgeMap source_edge_map_;

  // The number of edges added to S2Builder so far.
  InputEdgeId num_input_edges_;

  // The set of all input edge crossings, which is used by EdgeClippingLayer
  // to construct the clipped output polygon.
  InputEdgeCrossings input_crossings_;

  // Indicates whether the destination vertex of the previous edge added to
  // S2Builder was considered to be inside the polygon being clipped against.
  bool prev_inside_;
};

// Starts an operation.  "clip_to_exterior" specifies whether the edges of one
// region are clipped to the interior or exterior of the other region.
void S2BoundaryOperation::Impl::Start(bool clip_to_exterior) {
  clip_to_exterior_ = clip_to_exterior;
  num_input_edges_ = 0;
  prev_inside_ = false;
  if (!op_->layers_.empty()) {
    // TODO(ericv): Rather than having S2Builder split the edges, it would be
    // faster to call AddVertex() in this class and have a new S2Builder
    // option that increases the edge_snap_radius_ to account for errors in
    // the intersection point (the way that split_crossing_edges does).
    S2Builder::Options options(op_->options_.snap_function());
    options.set_split_crossing_edges(true);

    // The "idempotent" option doesn't do any snapping if the input already
    // meets the output guarantees; e.g. with the IdentitySnapFunction vertices
    // are not snapped if they are least snap_radius / 2 apart.  In our case we
    // want to snap all vertices that are close than the full "snap_radius".
    options.set_idempotent(false);

    builder_.Init(options);
    CHECK_EQ(op_->layers_.size(), 1);
    builder_.StartLayer(gtl::MakeUnique<EdgeClippingLayer>(
        std::move(op_->layers_[0]), &input_crossings_, clip_to_exterior_));
  }
}

// Clips the boundary of A to the interior or exterior of B (as specified by
// "clip_to_exterior_"), and adds the resulting edges to the output.  The region
// edges of A are reversed if "reverse_a" is true, and similarly for B.  This
// method must be called an even number of times (first to clip A to B and
// then to clip B to A), where each pair is followed by a call to either
// TranslateCrossings() or Build().
//
// If the same polygon edge appears in A and B, it is kept only if
// a_region_id == 0.  (This avoids having duplicate edges in the output.)
//
// If a polygon edge from A appears as a reversed polygon edge in B, it is
// discarded from both polygons unless polygon_model_ is PolygonModel::CLOSED.
void S2BoundaryOperation::Impl::AddBoundary(
    int a_region_id, bool reverse_a, int b_region_id, bool reverse_b) {
  DCHECK_EQ(op_->regions_.size(), 2);
  Region const& a = op_->regions_[a_region_id];
  Region const& b = op_->regions_[b_region_id];
  bool invert_b = reverse_b ^ clip_to_exterior_;

  // To keep only one copy of any polygon edges shared between A and B, we make
  // "keep_shared_edges" true for the first region and false for the second.
  bool keep_shared_edges = (a_region_id == 0);
  CrossingQuery query(b, reverse_b, keep_shared_edges, clip_to_exterior_);

  // Walk the boundary of region A and build a list of all edge crossings.
  // We also keep track of whether the current vertex is inside region B.
  vector<SourceEdgeCrossing> crossings;
  for (int i = 0; i < a.index.num_shape_ids(); ++i) {
    S2Shape* a_shape = a.index.shape(i);
    if (a_shape == nullptr) continue;
    int n = a_shape->num_edges();

    // When computing the difference or symmetric difference with a polyline,
    // subtracting a polyline from a polygon has no effect.
    if (reverse_a && a_shape->dimension() < 2) continue;

    // To preserve the cyclic order of the input vertices, for reversed loops
    // we generate a loop with vertex order (n-1, n-2, ..., 0) so that when it
    // is eventually inverted by S2Loop::Invert() we end up with the original
    // order (0, 1, ..., n-1).  The corresponding *edge* order of the reversed
    // loop is (n-2, n-3, ... 0, n-1).
    //
    // Example: Loop ABCD has edges AB, BC, CD, DA.
    // Reversed loop DCBA has edges DC, CB, BA, AD.
    bool reverse = a.reverse_edges(*a_shape) ^ reverse_a;
    bool inside;  // Does the current vertex belong to the result boundary?
    for (int j = 0; j < n; ++j) {
      int edge_id = j;
      S2Point const *a0, *a1;
      if (reverse) {  // See comments above.
        if (j <= n - 2) edge_id = (n - 2) - j;
        a_shape->GetEdge(edge_id, &a1, &a0);
      } else {
        a_shape->GetEdge(edge_id, &a0, &a1);
      }
      // Initialize "inside" on the first loop iteration.
      if (j == 0) inside = b.contains(*a0) ^ invert_b;
      int a0_crossings, a1_crossings;
      bool has_polyline_crossings;
      query.GetCrossings(*a0, *a1, &a0_crossings, &a1_crossings, &crossings,
                         &has_polyline_crossings);

      // Update "inside" to reflect any edge crossings at a0.
      inside ^= (a0_crossings & 1);
      if (inside || !crossings.empty() || has_polyline_crossings) {
        // "inside" may get out of sync when switching from one loop to the
        // next.  We fix this by adding an artificial extra crossing at the
        // start vertex of the current edge.
        if (inside != prev_inside_) {
          source_edge_crossings_.push_back(
              make_pair(num_input_edges_, kAtStartVertex()));
        }
        for (SourceEdgeCrossing const& crossing : crossings) {
          source_edge_crossings_.push_back(
              make_pair(num_input_edges_, crossing));
        }
        if (!crossings.empty() || has_polyline_crossings) {
          // Build a map that translates our temporary edge id (SourceId)
          // into the representation used by EdgeClippingLayer (InputEdgeId).
          source_edge_map_[SourceId(a.id, i, edge_id)] = num_input_edges_;
        }
        // Update "inside" to reflect any interior edge crossings.
        inside ^= (crossings.size() & 1);
        builder_.AddEdge(*a0, *a1);
        ++num_input_edges_;
        prev_inside_ = inside;
      }
      // Update "inside" to reflect any edge crossings at a1.
      inside ^= (a1_crossings & 1);
      DCHECK_EQ((b.contains(*a1) ^ invert_b), inside);
    }
  }
}

// Translates the temporary representation of crossing edges (SourceId) into
// the format expected by EdgeClippingLayer (InputEdgeId).
void S2BoundaryOperation::Impl::TranslateCrossings() {
  input_crossings_.reserve(source_edge_crossings_.size());
  for (auto const& tmp : source_edge_crossings_) {
    if (tmp.second == kAtStartVertex()) {
      input_crossings_.push_back(make_pair(
          tmp.first, CrossingInputEdge::kAtStartVertex()));
    } else {
      auto it = source_edge_map_.find(tmp.second.first);
      DCHECK(it != source_edge_map_.end());
      input_crossings_.push_back(make_pair(
          tmp.first, CrossingInputEdge(it->second, tmp.second.second)));
    }
  }
  source_edge_crossings_.clear();
  source_edge_map_.clear();
}

bool S2BoundaryOperation::Impl::Build(S2Error* error) {
  error->Clear();
  switch (op_->op_type_) {
    case OpType::UNION:
      // We want the boundary of A clipped to the exterior of B,
      // plus the boundary of B clipped to the exterior of A,
      // plus one copy of any directed edges that are in both boundaries.
      Start(true /*clip_to_exterior*/);
      AddBoundary(0, false, 1, false);
      AddBoundary(1, false, 0, false);
      break;

    case OpType::INTERSECTION:
      // We want the boundary of A clipped to the interior of B,
      // plus the boundary of B clipped to the interior of A,
      // plus one copy of any directed edges that are in both boundaries.
      Start(false /*clip_to_exterior*/);
      AddBoundary(0, false, 1, false);
      AddBoundary(1, false, 0, false);
      break;

    case OpType::DIFFERENCE:
      // Compute the intersection of A with the complement of B.
      Start(false /*clip_to_exterior*/);
      AddBoundary(0, false, 1, true);
      AddBoundary(1, true, 0, false);
      break;

    case OpType::SYMMETRIC_DIFFERENCE:
      // Compute the union of (A - B) and (B - A).
      Start(false /*clip_to_exterior*/);
      AddBoundary(0, false, 1, true);
      AddBoundary(1, true, 0, false);
      TranslateCrossings();
      AddBoundary(0, true, 1, false);
      AddBoundary(1, false, 0, true);
      break;
  }
  TranslateCrossings();
  std::sort(input_crossings_.begin(), input_crossings_.end());
  return builder_.Build(error);
}

S2BoundaryOperation::Options::Options()
    : Options(s2builderutil::IdentitySnapFunction(S1Angle::Zero())) {
}

S2BoundaryOperation::Options::Options(SnapFunction const& snap_function)
    : snap_function_(snap_function.Clone()),
      polygon_model_(PolygonModel::SEMI_OPEN),
      polyline_model_(PolylineModel::OPEN) {
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

PolylineModel S2BoundaryOperation::Options::polyline_model() const {
  return polyline_model_;
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
    : op_type_(op_type), options_(options), result_nonempty_(nullptr) {
}

S2BoundaryOperation::S2BoundaryOperation(
    OpType op_type, unique_ptr<S2Builder::Layer> layer, Options const& options)
    : S2BoundaryOperation(op_type, options) {
  layers_.push_back(std::move(layer));
}

static std::function<bool (S2Point const&)> ContainsUniquely(
    S2ShapeIndex const& index) {
  return [&index](S2Point const& p) {
    vector<S2Shape*> shapes;
    if (!index.GetContainingShapes(p, &shapes)) return false;
    DCHECK_EQ(shapes.size(), 1);
    return true;
  };
}

static bool ReverseNone(S2Shape const& shape) {
  return false;
}

void S2BoundaryOperation::AddRegion(S2ShapeIndex const& index) {
  AddRegion(index, ContainsUniquely(index), ReverseNone);
}

void S2BoundaryOperation::AddRegion(
    S2ShapeIndex const& index,
    std::function<bool (S2Point const&)> contains,
    std::function<bool (S2Shape const&)> reverse_edges) {
  DCHECK_LT(regions_.size(), 2);
  regions_.push_back(
      Region(regions_.size(), index, contains, reverse_edges));
}

bool S2BoundaryOperation::Build(S2Error* error) {
  return Impl(this).Build(error);
}
