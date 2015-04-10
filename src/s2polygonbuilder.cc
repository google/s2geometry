// Copyright 2006 Google Inc. All Rights Reserved.
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

#include "s2polygonbuilder.h"

#include <math.h>
#include <algorithm>
#include <ext/hash_map>
using __gnu_cxx::hash;
using __gnu_cxx::hash_map;
#include <ext/hash_set>
using __gnu_cxx::hash;
using __gnu_cxx::hash_set;
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <vector>

#include <glog/logging.h>
#include "base/macros.h"
#include "base/scoped_ptr.h"
#include "s2.h"
#include "s2cellid.h"
#include "s2edgeutil.h"
#include "s2latlng.h"
#include "s2loop.h"
#include "s2polygon.h"
#include "util/math/matrix3x3.h"

using std::min;
using std::pair;
using std::vector;

void S2PolygonBuilderOptions::set_undirected_edges(bool undirected_edges) {
  undirected_edges_ = undirected_edges;
}

void S2PolygonBuilderOptions::set_xor_edges(bool xor_edges) {
  xor_edges_ = xor_edges;
}

void S2PolygonBuilderOptions::set_validate(bool validate) {
  validate_ = validate;
  set_s2debug_override(DISABLE_S2DEBUG);
}

void S2PolygonBuilderOptions::set_s2debug_override(S2debugOverride override) {
  s2debug_override_ = override;
}

void S2PolygonBuilderOptions::set_snap_to_cell_centers(
    bool snap_to_cell_centers) {
  snap_to_cell_centers_ = snap_to_cell_centers;
}

void S2PolygonBuilderOptions::set_vertex_merge_radius(S1Angle angle) {
  vertex_merge_radius_ = angle;
}

void S2PolygonBuilderOptions::set_edge_splice_fraction(double fraction) {
  CHECK(fraction < sqrt(3) / 2);
  edge_splice_fraction_ = fraction;
}

void S2PolygonBuilderOptions::SetRobustnessRadius(S1Angle robustness_radius) {
  vertex_merge_radius_ = 2 * robustness_radius / edge_splice_fraction_;
}

S1Angle S2PolygonBuilderOptions::GetRobustnessRadius() const {
  return vertex_merge_radius_ * edge_splice_fraction_ / 2;
}

int S2PolygonBuilderOptions::GetSnapLevel() const {
  if (!snap_to_cell_centers())
    return -1;

  S1Angle const tolerance = GetRobustnessRadius();
  // The distance from a point in the cell to the center is at most
  // kMaxDiag / 2.  See the comment before S2::kMaxDiag.
  // TODO(user): Verify and understand this.
  int const level = S2::kMaxDiag.GetMinLevel(tolerance.radians() * 2.0);
  // GetMinLevel will return kMaxLevel even if the max level does not satisfy
  // the condition.
  if (level == S2CellId::kMaxLevel &&
      tolerance < S1Angle::Radians(S2::kMaxDiag.GetValue(level) / 2.0)) {
    return -1;
  }

  return level;
}

S2PolygonBuilder::S2PolygonBuilder(S2PolygonBuilderOptions const& options)
    : options_(options), edges_(new EdgeSet) {
}

S2PolygonBuilder::~S2PolygonBuilder() {
}

bool S2PolygonBuilder::HasEdge(S2Point const& v0, S2Point const& v1) {
  EdgeSet::const_iterator candidates = edges_->find(v0);
  return (candidates != edges_->end() &&
          candidates->second.find(v1) != candidates->second.end());
}

bool S2PolygonBuilder::AddEdge(S2Point const& v0, S2Point const& v1) {
  // If xor_edges is true, we look for an existing edge in the opposite
  // direction.  We either delete that edge or insert a new one.

  if (v0 == v1) return false;
  if (options_.xor_edges() && HasEdge(v1, v0)) {
    EraseEdge(v1, v0);
    return false;
  }
  if (edges_->find(v0) == edges_->end()) starting_vertices_.push_back(v0);
  (*edges_)[v0].insert(v1);
  if (options_.undirected_edges()) {
    if (edges_->find(v1) == edges_->end()) starting_vertices_.push_back(v1);
    (*edges_)[v1].insert(v0);
  }
  return true;
}

void S2PolygonBuilder::AddLoop(S2Loop const* loop) {
  // Ignore loops that do not have a boundary.
  if (loop->is_empty_or_full()) return;

  int sign = loop->sign();
  for (int i = loop->num_vertices(); i > 0; --i) {
    // Vertex indices need to be in the range [0, 2*num_vertices()-1].
    AddEdge(loop->vertex(i), loop->vertex(i + sign));
  }
}

void S2PolygonBuilder::AddPolygon(S2Polygon const* polygon) {
  for (int i = 0; i < polygon->num_loops(); ++i) {
    AddLoop(polygon->loop(i));
  }
}

void S2PolygonBuilder::EraseEdge(S2Point const& v0, S2Point const& v1) {
  // Note that there may be more than one copy of an edge if we are not XORing
  // them, so a VertexSet is a multiset.

  VertexSet* vset = &(*edges_)[v0];
  DCHECK(vset->find(v1) != vset->end());
  vset->erase(vset->find(v1));
  if (vset->empty()) edges_->erase(v0);

  if (options_.undirected_edges()) {
    vset = &(*edges_)[v1];
    DCHECK(vset->find(v0) != vset->end());
    vset->erase(vset->find(v0));
    if (vset->empty()) edges_->erase(v1);
  }
}

void S2PolygonBuilder::set_debug_matrix(Matrix3x3_d const& m) {
  debug_matrix_.reset(new Matrix3x3_d(m));
}

void S2PolygonBuilder::DumpVertex(S2Point const& v) const {
  if (debug_matrix_.get()) {
    // For orthonormal matrices, Inverse() == Transpose().
    std::cout << S2LatLng(debug_matrix_->Transpose() * v);
  } else {
    std::cout << std::setprecision(17) << v << std::setprecision(6);
  }
}

void S2PolygonBuilder::DumpEdges(S2Point const& v0) const {
  DumpVertex(v0);
  std::cout << ":\n";
  EdgeSet::const_iterator candidates = edges_->find(v0);
  if (candidates != edges_->end()) {
    VertexSet const& vset = candidates->second;
    for (VertexSet::const_iterator i = vset.begin(); i != vset.end(); ++i) {
      std::cout << "    ";
      DumpVertex(*i);
      std::cout << "\n";
    }
  }
}

void S2PolygonBuilder::Dump() const {
  for (EdgeSet::const_iterator i = edges_->begin(); i != edges_->end(); ++i) {
    DumpEdges(i->first);
  }
}

void S2PolygonBuilder::EraseLoop(S2Point const* v, int n) {
  for (int i = n - 1, j = 0; j < n; i = j++) {
    EraseEdge(v[i], v[j]);
  }
}

void S2PolygonBuilder::RejectLoop(S2Point const* v, int n,
                                  EdgeList* unused_edges) {
  for (int i = n - 1, j = 0; j < n; i = j++) {
    unused_edges->push_back(std::make_pair(v[i], v[j]));
  }
}

S2Loop* S2PolygonBuilder::AssembleLoop(S2Point const& v0, S2Point const& v1,
                                       EdgeList* unused_edges) {
  // We start at the given edge and assemble a loop taking left turns
  // whenever possible.  We stop the loop as soon as we encounter any
  // vertex that we have seen before *except* for the first vertex (v0).
  // This ensures that only CCW loops are constructed when possible.

  vector<S2Point> path;          // The path so far.
  // Maps a vertex to its index in "path".
  hash_map<S2Point, int, HashS2Point> index;

  path.push_back(v0);
  path.push_back(v1);
  index[v1] = 1;
  while (path.size() >= 2) {
    // Note that "v0" and "v1" become invalid if "path" is modified.
    S2Point const& v0 = path.end()[-2];
    S2Point const& v1 = path.end()[-1];
    S2Point v2;
    bool v2_found = false;
    EdgeSet::const_iterator candidates = edges_->find(v1);
    if (candidates != edges_->end()) {
      VertexSet const& vset = candidates->second;
      for (VertexSet::const_iterator i = vset.begin(); i != vset.end(); ++i) {
        // We prefer the leftmost outgoing edge, ignoring any reverse edges.
        if (*i == v0) continue;
        if (!v2_found || S2::OrderedCCW(v0, v2, *i, v1)) { v2 = *i; }
        v2_found = true;
      }
    }
    if (!v2_found) {
      // We've hit a dead end.  Remove this edge and backtrack.
      unused_edges->push_back(std::make_pair(v0, v1));
      EraseEdge(v0, v1);
      index.erase(v1);
      path.pop_back();
    } else if (index.insert(std::make_pair(v2, path.size())).second) {
      // This is the first time we've visited this vertex.
      path.push_back(v2);
    } else {
      // We've completed a loop. In a simple case last edge is the same as the
      // first one (since we did not add the very first vertex to the index we
      // would not know that the loop is completed till the second vertex is
      // examined). In this case we just remove the last edge to preserve the
      // original vertex order. In a more complicated case the edge that closed
      // the loop is different and we should remove initial vertices that are
      // not part of the loop.
      if (index[v2] == 1 && path[0]  == path[path.size() - 1]) {
        path.pop_back();
      } else {
        path.erase(path.begin(), path.begin() + index[v2]);
      }

      // In the case of undirected edges, we may have assembled a clockwise
      // loop while trying to assemble a CCW loop.  To fix this, we assemble
      // a new loop starting with an arbitrary edge in the reverse direction.
      // This is guaranteed to assemble a loop that is interior to the previous
      // one and will therefore eventually terminate.

      S2Loop* loop = new S2Loop(path, options_.s2debug_override());
      if (options_.validate() && !loop->IsValid()) {
        // We've constructed a loop that crosses itself, which can only
        // happen if there is bad input data.  Throw away the whole loop.
        RejectLoop(&path[0], path.size(), unused_edges);
        EraseLoop(&path[0], path.size());
        delete loop;
        return NULL;
      }

      if (options_.undirected_edges() && !loop->IsNormalized()) {
        scoped_ptr<S2Loop> deleter(loop);  // XXX for debugging
        return AssembleLoop(path[1], path[0], unused_edges);
      }
      return loop;
    }
  }
  return NULL;
}

class S2PolygonBuilder::PointIndex {
  // A PointIndex is a cheap spatial index to help us find mergeable
  // vertices.  Given a set of points, it can efficiently find all of the
  // points within a given search radius of an arbitrary query location.
  // It is essentially just a hash map from cell ids at a given fixed level to
  // the set of points contained by that cell id.
  //
  // This class is not suitable for general use because it only supports
  // fixed-radius queries and has various special-purpose operations to avoid
  // the need for additional data structures.
  //
  // Thread-unsafe as FindNearbyPoint modifies ids_.  This is of no consequence,
  // since the only use of PointIndex is as a local variable in the body of
  // AssembleLoops which is single-threaded.

 private:
  typedef std::multimap<S2CellId, S2Point> Map;
  Map map_;

  double vertex_radius_;
  double edge_fraction_;
  int level_;
  // Scratch space allocated here for efficiency.  Used by Insert, Erase,
  // and FindNearbyPoint.  Should be empty otherwise.
  mutable vector<S2CellId> ids_;

 public:
  PointIndex(double vertex_radius, double edge_fraction)
    : vertex_radius_(vertex_radius),
      edge_fraction_(edge_fraction),
      // We compute an S2CellId level such that the vertex neighbors at that
      // level of any point A are a covering for spherical cap (i.e. "disc")
      // of the given search radius centered at A.  This requires that the
      // minimum cell width at that level must be twice the search radius.
      level_(min(S2::kMinWidth.GetMaxLevel(2 * vertex_radius),
                 S2CellId::kMaxLevel - 1)) {
    // We insert a sentinel so that we don't need to test for map_.end().
    map_.insert(std::make_pair(S2CellId::Sentinel(), S2Point()));
  }

  void Insert(S2Point const& p) {
    DCHECK(ids_.empty());
    S2CellId::FromPoint(p).AppendVertexNeighbors(level_, &ids_);
    for (int i = ids_.size(); --i >= 0; ) {
      map_.insert(std::make_pair(ids_[i], p));
    }
    ids_.clear();
  }

  void Erase(S2Point const& p) {
    DCHECK(ids_.empty());
    S2CellId::FromPoint(p).AppendVertexNeighbors(level_, &ids_);
    for (int i = ids_.size(); --i >= 0; ) {
      Map::iterator j = map_.lower_bound(ids_[i]);
      for (; j->second != p; ++j) {
        DCHECK_EQ(ids_[i], j->first);
      }
      map_.erase(j);
    }
    ids_.clear();
  }

  void QueryCap(S2Point const& axis, vector<S2Point>* output) {
    // Return the set the points whose distance to "axis" is less than
    // vertex_radius_.

    output->clear();
    S2CellId id = S2CellId::FromPoint(axis).parent(level_);
    for (Map::const_iterator i = map_.lower_bound(id); i->first == id; ++i) {
      S2Point const& p = i->second;
      if (axis.Angle(p) < vertex_radius_) {
        output->push_back(p);
      }
    }
  }

  bool FindNearbyPoint(S2Point const& v0, S2Point const& v1,
                       S2Point* nearby) const {
    // Return a point whose distance from the edge (v0,v1) is less than
    // vertex_radius_, and which is not equal to v0 or v1.  The current
    // implementation returns the closest such point.
    //
    // Strategy: we compute a very cheap covering by approximating the edge as
    // two spherical caps, one around each endpoint, and then computing a
    // 4-cell covering of each one.  We could improve the quality of the
    // covering by using some intermediate points along the edge as well.
    //
    // Logically const but modifies ids_
    // which is required to be empty at start and end of this function.

    DCHECK(ids_.empty());
    double length = v0.Angle(v1);
    int level = min(level_, S2::kMinWidth.GetMaxLevel(length));
    S2CellId::FromPoint(v0).AppendVertexNeighbors(level, &ids_);
    S2CellId::FromPoint(v1).AppendVertexNeighbors(level, &ids_);

    // Sort the cell ids so that we can skip duplicates in the loop below.
    std::sort(ids_.begin(), ids_.end());

    double best_dist = 2 * vertex_radius_;
    for (int i = ids_.size(); --i >= 0; ) {
      if (i > 0 && ids_[i-1] == ids_[i]) continue;  // Skip duplicates.

      S2CellId max_id = ids_[i].range_max();
      for (Map::const_iterator j = map_.lower_bound(ids_[i].range_min());
           j->first <= max_id; ++j) {
        S2Point const& p = j->second;
        if (p == v0 || p == v1) continue;
        double dist = S2EdgeUtil::GetDistance(p, v0, v1).radians();
        if (dist < best_dist) {
          best_dist = dist;
          *nearby = p;
        }
      }
    }
    ids_.clear();
    return (best_dist < edge_fraction_ * vertex_radius_);
  }

 private:
  DISALLOW_COPY_AND_ASSIGN(PointIndex);
};

void S2PolygonBuilder::BuildMergeMap(PointIndex* index, MergeMap* merge_map) {
  // The overall strategy is to start from each vertex and grow a maximal
  // cluster of mergeable vertices.  In graph theoretic terms, we find the
  // connected components of the undirected graph whose edges connect pairs of
  // vertices that are separated by at most vertex_merge_radius().
  //
  // We then choose a single representative vertex for each cluster, and
  // update "merge_map" appropriately.  We choose an arbitrary existing
  // vertex rather than computing the centroid of all the vertices to avoid
  // creating new vertex pairs that need to be merged.  (We guarantee that all
  // vertex pairs are separated by at least the merge radius in the output.)

  // First, we build the set of all the distinct vertices in the input.
  // We need to include the source and destination of every edge.
  hash_set<S2Point, HashS2Point> vertices;
  for (EdgeSet::const_iterator i = edges_->begin(); i != edges_->end(); ++i) {
    vertices.insert(i->first);
    VertexSet const& vset = i->second;
    for (VertexSet::const_iterator j = vset.begin(); j != vset.end(); ++j)
      vertices.insert(*j);
  }

  // Build a spatial index containing all the distinct vertices.
  for (hash_set<S2Point, HashS2Point>::const_iterator i = vertices.begin();
       i != vertices.end(); ++i) {
    index->Insert(*i);
  }

  // Next, we loop through all the vertices and attempt to grow a maximial
  // mergeable group starting from each vertex.
  vector<S2Point> frontier, mergeable;
  for (hash_set<S2Point, HashS2Point>::const_iterator vstart = vertices.begin();
       vstart != vertices.end(); ++vstart) {
    // Skip any vertices that have already been merged with another vertex.
    if (merge_map->find(*vstart) != merge_map->end()) continue;

    // Grow a maximal mergeable component starting from "vstart", the
    // canonical representative of the mergeable group.
    frontier.push_back(*vstart);
    while (!frontier.empty()) {
      index->QueryCap(frontier.back(), &mergeable);
      frontier.pop_back();  // Do this before entering the loop below.
      for (int j = mergeable.size(); --j >= 0; ) {
        S2Point const& v1 = mergeable[j];
        if (v1 != *vstart) {
          // Erase from the index any vertices that will be merged.  This
          // ensures that we won't try to merge the same vertex twice.
          index->Erase(v1);
          frontier.push_back(v1);
          (*merge_map)[v1] = *vstart;
        }
      }
    }
  }
}

void S2PolygonBuilder::MoveVertices(MergeMap const& merge_map) {
  if (merge_map.empty()) return;

  // We need to copy the set of edges affected by the move, since
  // edges_ could be reallocated when we start modifying it.
  vector<pair<S2Point, S2Point> > edges;
  for (EdgeSet::const_iterator i = edges_->begin(); i != edges_->end(); ++i) {
    S2Point const& v0 = i->first;
    VertexSet const& vset = i->second;
    for (VertexSet::const_iterator j = vset.begin(); j != vset.end(); ++j) {
      S2Point const& v1 = *j;
      if (merge_map.find(v0) != merge_map.end() ||
          merge_map.find(v1) != merge_map.end()) {
        // We only need to modify one copy of each undirected edge.
        if (!options_.undirected_edges() || v0 < v1) {
          edges.push_back(std::make_pair(v0, v1));
        }
      }
    }
  }

  // Now erase all the old edges, and add all the new edges.  This will
  // automatically take care of any XORing that needs to be done, because
  // EraseEdge also erases the sibling of undirected edges.
  for (int i = 0; i < edges.size(); ++i) {
    S2Point v0 = edges[i].first;
    S2Point v1 = edges[i].second;
    EraseEdge(v0, v1);
    MergeMap::const_iterator new0 = merge_map.find(v0);
    if (new0 != merge_map.end()) v0 = new0->second;
    MergeMap::const_iterator new1 = merge_map.find(v1);
    if (new1 != merge_map.end()) v1 = new1->second;
    AddEdge(v0, v1);
  }
}

void S2PolygonBuilder::SpliceEdges(PointIndex const& index) {
  // We keep a stack of unprocessed edges.  Initially all edges are
  // pushed onto the stack.
  vector<pair<S2Point, S2Point> > edges;
  for (EdgeSet::const_iterator i = edges_->begin(); i != edges_->end(); ++i) {
    S2Point const& v0 = i->first;
    VertexSet const& vset = i->second;
    for (VertexSet::const_iterator j = vset.begin(); j != vset.end(); ++j) {
      S2Point const& v1 = *j;
      // We only need to modify one copy of each undirected edge.
      if (!options_.undirected_edges() || v0 < v1) {
        edges.push_back(std::make_pair(v0, v1));
      }
    }
  }

  // For each edge, we check whether there are any nearby vertices that should
  // be spliced into it.  If there are, we choose one such vertex, split
  // the edge into two pieces, and iterate on each piece.
  while (!edges.empty()) {
    S2Point v0 = edges.back().first;
    S2Point v1 = edges.back().second;
    edges.pop_back();  // Do this before pushing new edges.

    // If we are xoring, edges may be erased before we get to them.
    if (options_.xor_edges() && !HasEdge(v0, v1)) continue;

    S2Point vmid;
    if (!index.FindNearbyPoint(v0, v1, &vmid)) continue;

    EraseEdge(v0, v1);
    if (AddEdge(v0, vmid)) edges.push_back(std::make_pair(v0, vmid));
    if (AddEdge(vmid, v1)) edges.push_back(std::make_pair(vmid, v1));
  }
}

namespace {
// Snaps the point to the center of the cell containing it at the specified
// level.
S2Point SnapPointToLevel(S2Point const& p, int level) {
  return S2CellId::FromPoint(p).parent(level).ToPoint();
}

// Returns a newly allocated loop, whose vertices have been snapped to
// the centers of cells at the specified level.  The caller owns the
// returned loop and is responsible for deleting it.
S2Loop* SnapLoopToLevel(S2Loop const& loop, int level) {
  vector<S2Point> snapped_vertices(loop.num_vertices());
  for (int i = 0; i < loop.num_vertices(); ++i) {
    snapped_vertices[i] = SnapPointToLevel(loop.vertex(i), level);
  }
  return new S2Loop(snapped_vertices);
}
}

bool S2PolygonBuilder::AssembleLoops(vector<S2Loop*>* loops,
                                     EdgeList* unused_edges) {
  if (options_.vertex_merge_radius().radians() > 0) {
    PointIndex index(options_.vertex_merge_radius().radians(),
                     options_.edge_splice_fraction());
    MergeMap merge_map;
    BuildMergeMap(&index, &merge_map);
    MoveVertices(merge_map);
    if (options_.edge_splice_fraction() > 0) {
      SpliceEdges(index);
    }
  }

  int const snap_level = options_.GetSnapLevel();

  EdgeList dummy_unused_edges;
  if (unused_edges == NULL) unused_edges = &dummy_unused_edges;

  // We repeatedly choose an edge and attempt to assemble a loop
  // starting from that edge.  (This is always possible unless the
  // input includes extra edges that are not part of any loop.)  To
  // maintain a consistent scanning order over edges_ between
  // different machine architectures (e.g. 'clovertown' vs. 'opteron'),
  // we follow the order they were added to the builder.
  unused_edges->clear();
  for (int i = 0; i < starting_vertices_.size(); ) {
    S2Point const& v0 = starting_vertices_[i];
    EdgeSet::const_iterator candidates = edges_->find(v0);
    if (candidates == edges_->end()) {
      ++i;
      continue;
    }
    // NOTE(user): If we have such two S2Points a, b that:
    //
    //   a.x = b.x, a.y = b.y and
    //   -- a.z = b.z if CPU is Intel
    //   -- a.z <> b.z if CPU is AMD
    //
    // then the order of points picked up as v1 on the following line
    // can be inconsistent between different machine architectures.
    //
    // As of b/3088321 and of cl/17847332, it's not clear if such
    // input really exists in our input and probably it's O.K. not to
    // address it in favor of the speed.
    S2Point const& v1 = *(candidates->second.begin());
    S2Loop* loop = AssembleLoop(v0, v1, unused_edges);
    if (loop != NULL) {
      EraseLoop(&loop->vertex(0), loop->num_vertices());

      if (snap_level >= 0) {
        // TODO(user): Change AssembleLoop to return a vector<S2Point>,
        //   then optionally snap that before constructing the loop.
        //   This would prevent us from constructing two loops.
        S2Loop* snapped_loop = SnapLoopToLevel(*loop, snap_level);
        delete loop;
        loop = snapped_loop;
      }

      loops->push_back(loop);
    }
  }
  return unused_edges->empty();
}

bool S2PolygonBuilder::AssemblePolygon(S2Polygon* polygon,
                                       EdgeList* unused_edges) {
  vector<S2Loop*> loops;
  bool success = AssembleLoops(&loops, unused_edges);
  polygon->set_s2debug_override(options_.s2debug_override());
  if (options_.undirected_edges()) {
    // If edges are undirected, then all loops are already normalized (i.e.,
    // contain at most half the sphere).  This implies that no loop contains
    // the complement of any other loop, and therefore we can call the normal
    // Init() method.
    polygon->Init(&loops);
  } else {
    // If edges are directed, then shells and holes have opposite orientations
    // such that the polygon interior is always on their left-hand side.
    polygon->InitOriented(&loops);
  }
  if (options_.validate() && !polygon->IsValid()) {
    polygon->Release(&loops);
    if (unused_edges != NULL) {
      for (int i = 0; i < loops.size(); ++i) {
        RejectLoop(&loops[i]->vertex(0), loops[i]->num_vertices(),
                   unused_edges);
      }
    }
    for (int i = 0; i < loops.size(); ++i) {
      delete loops[i];
    }
    return false;
  }
  return success;
}