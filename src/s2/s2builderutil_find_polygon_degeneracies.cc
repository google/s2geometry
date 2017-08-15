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

#include "s2/s2builderutil_find_polygon_degeneracies.h"

#include <cstdlib>
#include <utility>
#include <vector>

#include "s2/third_party/absl/memory/memory.h"
#include "s2/util/btree/btree_map.h"
#include "s2/s2builder_graph.h"
#include "s2/s2crossingedgequery.h"
#include "s2/s2edge_crosser.h"
#include "s2/s2pointutil.h"
#include "s2/s2predicates.h"

using absl::MakeUnique;
using std::make_pair;
using std::pair;
using std::vector;

using Graph = S2Builder::Graph;
using GraphOptions = S2Builder::GraphOptions;

using Edge = Graph::Edge;
using EdgeId = Graph::EdgeId;
using VertexId = Graph::VertexId;

namespace s2builderutil {

// This class determines whether a polygon contains one of its vertices given
// the edges incident to that vertex.  The result is +1 if the vertex is
// contained, -1 if it is not contained, and 0 if the incident edges consist
// of matched sibling pairs (in which case the result cannot be determined
// locally).
//
// Point containment is defined according to the "semi-open" boundary model
// (see s2boundary_operation.h), which means that if several polygons tile the
// region around a vertex, exactly one of those polygons contains that vertex.
//
// TODO(ericv): Make this class public and probably rename.
class S2ContainsVertexQuery {
 public:
  // "target" is the vertex whose containment will be determined.
  explicit S2ContainsVertexQuery(S2Point const& target) : target_(target) {}

  // Indicates that the polygon has an edge between "target" and "v" in the
  // given direction (+1 = outgoing, -1 = incoming).
  void AddEdge(S2Point const& v, int direction) {
    edge_map_[v] += direction;
  }

  // Returns +1 if the vertex is contained, -1 if it is not contained, and 0
  // if the incident edges consisted of matched sibling pairs.
  int ContainsSigned();

 private:
  S2Point target_;
  util::btree::btree_map<S2Point, int> edge_map_;
};

int S2ContainsVertexQuery::ContainsSigned() {
  // Find the unmatched edge that is immediately clockwise from S2::Ortho(P).
  S2Point reference_dir = S2::Ortho(target_);
  pair<S2Point, int> best(reference_dir, 0);
  for (auto const& e : edge_map_) {
    DCHECK_LE(std::abs(e.second), 1);
    if (e.second == 0) continue;  // This is a "matched" edge.
    if (s2pred::OrderedCCW(reference_dir, best.first, e.first, target_)) {
      best = e;
    }
  }
  return best.second;
}

namespace {

// The algorithm builds a set of connected components containing all edges
// that form degeneracies.  The shell/hole status of each degeneracy is
// initially unknown, and is expressed relative to the root vertex: "is_hole"
// means that the degeneracy is a hole if and only if the root vertex turns
// out to be inside the polygon.
struct Component {
  // The root vertex from which this component was built.
  VertexId root;

  // +1 if "root" inside the polygon, -1 if outside, and 0 if unknown.
  int root_sign = 0;

  // The degeneracies found in this component.  "is_hole" is expressed
  // relative to the root vertex: the degeneracy is a hole iff the root vertex
  // turns out to be inside the polygon (i.e., root_sign > 0).
  vector<PolygonDegeneracy> degeneracies;
};

// The actual implementation of FindPolygonDegeneracies.
class DegeneracyFinder {
 public:
  explicit DegeneracyFinder(S2Builder::Graph const* g)
      : g_(*g), in_(g_), out_(g_) {
  }
  vector<PolygonDegeneracy> Run(S2Error* error);

 private:
  // Methods are documented below.
  int ComputeDegeneracies();
  Component BuildComponent(VertexId root);
  bool CrossingParity(VertexId v0, VertexId v1, bool include_same) const;
  VertexId FindUnbalancedVertex() const;
  int ContainsVertexSigned(VertexId v0) const;
  void ComputeUnknownSignsBruteForce(VertexId known_vertex,
                                     int known_vertex_sign,
                                     vector<Component>* components) const;
  void ComputeUnknownSignsIndexed(VertexId known_vertex, int known_vertex_sign,
                                  vector<Component>* components) const;
  vector<PolygonDegeneracy> MergeDegeneracies(
      vector<Component> const& components) const;

  Graph const& g_;
  Graph::VertexInMap in_;
  Graph::VertexOutMap out_;
  vector<bool> is_vertex_used_;        // Has vertex been visited?
  vector<bool> is_edge_degeneracy_;    // Belongs to a degeneracy?
  vector<bool> is_vertex_unbalanced_;  // Has unbalanced sibling pairs?
};

vector<PolygonDegeneracy> DegeneracyFinder::Run(S2Error* error) {
  // Mark all degenerate edges and sibling pairs in the "is_edge_degeneracy_"
  // vector, and mark any vertices with unbalanced edges in the
  // "is_vertex_unbalanced_" vector.
  int num_degeneracies = ComputeDegeneracies();
  if (num_degeneracies == 0) {
    return vector<PolygonDegeneracy>();
  }

  // If all edges are degenerate, then use IsFullPolygon() to classify the
  // degeneracies (they are necessarily all the same type).
  if (num_degeneracies == g_.num_edges()) {
    bool is_hole = g_.IsFullPolygon(error);
    vector<PolygonDegeneracy> result(g_.num_edges());
    for (int e = 0; e < g_.num_edges(); ++e) {
      result[e] = PolygonDegeneracy(e, is_hole);
    }
    return result;
  }

  // Otherwise repeatedly build components starting from an unvisited
  // degeneracy.  (This avoids building components that don't contain any
  // degeneracies.)  Each component records the "is_hole" status of each
  // degeneracy relative to the root vertex of that component.  If the
  // component contains any non-degenerate portions, then we also determine
  // whether the root vertex is contained by the component (root_sign).
  // In addition we keep track of the number of components that were
  // completely degenerate (to help us decide whether to build an index).
  vector<Component> components;
  VertexId known_vertex = -1;
  int known_vertex_sign = 0;
  int num_unknown_signs = 0;
  is_vertex_used_.resize(g_.num_vertices());
  for (int e = 0; e < g_.num_edges(); ++e) {
    if (is_edge_degeneracy_[e]) {
      VertexId root = g_.edge(e).first;
      if (is_vertex_used_[root]) continue;
      Component component = BuildComponent(root);
      if (component.root_sign == 0) {
        ++num_unknown_signs;
      } else {
        known_vertex = root;
        known_vertex_sign = component.root_sign;
      }
      components.push_back(component);
    }
  }

  // If some components have an unknown root_sign (i.e., it is unknown whether
  // the root vertex is contained by the polygon or not), we determine the
  // sign of those root vertices by counting crossings starting from a vertex
  // whose sign is known.  Depending on how many components we need to do this
  // for, it may be worthwhile to build an index first.
  if (num_unknown_signs > 0) {
    if (known_vertex_sign == 0) {
      known_vertex = FindUnbalancedVertex();
      known_vertex_sign = ContainsVertexSigned(known_vertex);
    }
    int const kMaxUnindexedContainsCalls = 20;  // Tuned using benchmarks.
    if (num_unknown_signs <= kMaxUnindexedContainsCalls) {
      ComputeUnknownSignsBruteForce(known_vertex, known_vertex_sign,
                                    &components);
    } else {
      ComputeUnknownSignsIndexed(known_vertex, known_vertex_sign,
                                 &components);
    }
  }
  // Finally we convert the "is_hole" status of each degeneracy from a
  // relative value (compared to the component's root vertex) to an absolute
  // one, and sort all the degeneracies by EdgeId.
  return MergeDegeneracies(components);
}

int DegeneracyFinder::ComputeDegeneracies() {
  is_edge_degeneracy_.resize(g_.num_edges());
  is_vertex_unbalanced_.resize(g_.num_vertices());
  int num_degeneracies = 0;
  vector<EdgeId> const& in_edge_ids = in_.in_edge_ids();
  int n = g_.num_edges();
  for (int in = 0, out = 0; out < n; ++out) {
    Edge out_edge = g_.edge(out);
    if (out_edge.first == out_edge.second) {
      is_edge_degeneracy_[out] = true;
      ++num_degeneracies;
    } else {
      while (in < n && Graph::reverse(g_.edge(in_edge_ids[in])) < out_edge) {
        ++in;
      }
      if (in < n && Graph::reverse(g_.edge(in_edge_ids[in])) == out_edge) {
        is_edge_degeneracy_[out] = true;
        ++num_degeneracies;
      } else {
        // This edge does not have a sibling, which mean that we can determine
        // whether either vertex is contained by the polygon (using semi-open
        // boundaries) by examining only the edges incident to that vertex.
        // We only mark the first vertex since there is no advantage to
        // finding more than one unbalanced vertex per connected component.
        is_vertex_unbalanced_[out_edge.first] = true;
      }
    }
  }
  return num_degeneracies;
}

// Build a connected component starting at the given root vertex.  The
// information returned includes: the root vertex, whether the containment
// status of the root vertex could be determined using only the edges in this
// component, and a vector of the edges that belong to degeneracies along with
// the shell/hole status of each such edge relative to the root vertex.
Component DegeneracyFinder::BuildComponent(VertexId root) {
  Component result;
  result.root = root;
  // We keep track of the frontier of unexplored vertices, and whether each
  // vertex is on the same side of the polygon boundary as the root vertex.
  vector<pair<VertexId, bool>> frontier;
  frontier.push_back(make_pair(root, true));
  is_vertex_used_[root] = true;
  while (!frontier.empty()) {
    VertexId v0 = frontier.back().first;
    bool v0_same_inside = frontier.back().second;  // Same as root vertex?
    frontier.pop_back();
    if (result.root_sign == 0 && is_vertex_unbalanced_[v0]) {
      int v0_sign = ContainsVertexSigned(v0);
      DCHECK_NE(v0_sign, 0);
      result.root_sign = v0_same_inside ? v0_sign : -v0_sign;
    }
    for (EdgeId e : out_.edge_ids(v0)) {
      VertexId v1 = g_.edge(e).second;
      bool same_inside = v0_same_inside ^ CrossingParity(v0, v1, false);
      if (is_edge_degeneracy_[e]) {
        result.degeneracies.push_back(PolygonDegeneracy(e, same_inside));
      }
      if (is_vertex_used_[v1]) continue;
      same_inside ^= CrossingParity(v1, v0, true);
      frontier.push_back(make_pair(v1, same_inside));
      is_vertex_used_[v1] = true;
    }
  }
  return result;
}

// Counts the number of times that (v0, v1) crosses the edges incident to v0,
// and returns the result modulo 2.  This is equivalent to calling
// S2::VertexCrossing for the edges incident to v0, except that this
// implementation is more efficient (since it doesn't need to determine which
// two edge vertices are the same).
//
// If "include_same" is false, then the edge (v0, v1) and its sibling (v1, v0)
// (if any) are excluded from the parity calculation.
bool DegeneracyFinder::CrossingParity(VertexId v0, VertexId v1,
                                      bool include_same) const {
  int crossings = 0;
  S2Point p0 = g_.vertex(v0);
  S2Point p1 = g_.vertex(v1);
  S2Point p0_ref = S2::Ortho(p0);
  for (Edge const& edge : out_.edges(v0)) {
    if (edge.second == v1) {
      if (include_same) ++crossings;
    } else if (s2pred::OrderedCCW(p0_ref, g_.vertex(edge.second), p1, p0)) {
      ++crossings;
    }
  }
  for (EdgeId e : in_.edge_ids(v0)) {
    Edge edge = g_.edge(e);
    if (edge.first == v1) {
      if (include_same) ++crossings;
    } else if (s2pred::OrderedCCW(p0_ref, g_.vertex(edge.first), p1, p0)) {
      ++crossings;
    }
  }
  return crossings & 1;
}

VertexId DegeneracyFinder::FindUnbalancedVertex() const {
  for (VertexId v = 0; v < g_.num_vertices(); ++v) {
    if (is_vertex_unbalanced_[v]) return v;
  }
  LOG(DFATAL) << "Could not find previously marked unbalanced vertex";
  return -1;
}

int DegeneracyFinder::ContainsVertexSigned(VertexId v0) const {
  S2ContainsVertexQuery query(g_.vertex(v0));
  for (Edge const& edge : out_.edges(v0)) {
    query.AddEdge(g_.vertex(edge.second), 1);
  }
  for (EdgeId e : in_.edge_ids(v0)) {
    query.AddEdge(g_.vertex(g_.edge(e).first), -1);
  }
  return query.ContainsSigned();
}

// Determines any unknown signs of component root vertices by counting
// crossings starting from a vertex whose sign is known.  This version simply
// tests all edges for crossings.
void DegeneracyFinder::ComputeUnknownSignsBruteForce(
    VertexId known_vertex, int known_vertex_sign,
    vector<Component>* components) const {
  S2EdgeCrosser crosser;
  for (Component& component : *components) {
    if (component.root_sign != 0) continue;
    bool inside = known_vertex_sign > 0;
    crosser.Init(&g_.vertex(known_vertex), &g_.vertex(component.root));
    for (EdgeId e = 0; e < g_.num_edges(); ++e) {
      if (is_edge_degeneracy_[e]) continue;
      Edge const& edge = g_.edge(e);
      inside ^= crosser.EdgeOrVertexCrossing(&g_.vertex(edge.first),
                                             &g_.vertex(edge.second));
    }
    component.root_sign = inside ? 1 : -1;
  }
}

// An S2Shape representing the edges in an S2Builder::Graph.
class GraphShape : public S2Shape {
 public:
  explicit GraphShape(Graph const* g) : g_(*g) {}
  int num_edges() const final { return g_.num_edges(); }
  Edge edge(int e) const final {
    Graph::Edge g_edge = g_.edge(e);
    return Edge(g_.vertex(g_edge.first), g_.vertex(g_edge.second));
  }
  int dimension() const final { return 1; }
  bool contains_origin() const final { return false; }
  int num_chains() const final { return g_.num_edges(); }
  Chain chain(int i) const final { return Chain(i, 1); }
  Edge chain_edge(int i, int j) const final { return edge(i); }
  ChainPosition chain_position(int e) const final {
    return ChainPosition(e, 0);
  }

 private:
  Graph const& g_;
};

// Like ComputeUnknownSignsBruteForce, except that this method uses an index
// to find the set of edges that cross a given edge.
void DegeneracyFinder::ComputeUnknownSignsIndexed(
    VertexId known_vertex, int known_vertex_sign,
    vector<Component>* components) const {
  S2ShapeIndex index;
  index.Add(MakeUnique<GraphShape>(&g_));
  S2CrossingEdgeQuery query(&index);
  vector<EdgeId> crossing_edges;
  S2EdgeCrosser crosser;
  for (Component& component : *components) {
    if (component.root_sign != 0) continue;
    query.GetCandidates(g_.vertex(known_vertex), g_.vertex(component.root),
                       index.shape(0), &crossing_edges);
    bool inside = known_vertex_sign > 0;
    crosser.Init(&g_.vertex(known_vertex), &g_.vertex(component.root));
    for (EdgeId e : crossing_edges) {
      if (is_edge_degeneracy_[e]) continue;
      inside ^= crosser.EdgeOrVertexCrossing(&g_.vertex(g_.edge(e).first),
                                             &g_.vertex(g_.edge(e).second));
    }
    component.root_sign = inside ? 1 : -1;
  }
}

// Merges the degeneracies from all components together, and computes the
// final "is_hole" status of each edge (since up to this point, the "is_hole"
// value has been expressed relative to the root vertex of each component).
vector<PolygonDegeneracy> DegeneracyFinder::MergeDegeneracies(
    vector<Component> const& components) const {
  vector<PolygonDegeneracy> result;
  for (Component const& component : components) {
    DCHECK_NE(component.root_sign, 0);
    bool invert = component.root_sign < 0;
    for (auto const& d : component.degeneracies) {
      result.push_back(PolygonDegeneracy(d.edge_id, d.is_hole ^ invert));
    }
  }
  std::sort(result.begin(), result.end());
  return result;
}

}  // namespace

vector<PolygonDegeneracy> FindPolygonDegeneracies(Graph const& graph,
                                                  S2Error* error) {
  using DegenerateEdges = GraphOptions::DegenerateEdges;
  using SiblingPairs = GraphOptions::SiblingPairs;
  DCHECK(graph.options().degenerate_edges() == DegenerateEdges::DISCARD ||
         graph.options().degenerate_edges() == DegenerateEdges::DISCARD_EXCESS);
  DCHECK(graph.options().sibling_pairs() == SiblingPairs::DISCARD ||
         graph.options().sibling_pairs() == SiblingPairs::DISCARD_EXCESS);

  return DegeneracyFinder(&graph).Run(error);
}

}  // namespace s2builderutil
