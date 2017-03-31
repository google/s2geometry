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

#ifndef S2_S2CROSSINGEDGEQUERY_H_
#define S2_S2CROSSINGEDGEQUERY_H_

#include <type_traits>
#include <vector>

#include "s2/third_party/absl/base/macros.h"
#include "s2/third_party/absl/container/inlined_vector.h"
#include "s2/util/btree/btree_map.h"  // Like std::map, but faster and smaller.
#include "s2/fpcontractoff.h"
#include "s2/r2.h"
#include "s2/r2rect.h"
#include "s2/s2shapeindex.h"
#include "s2/s2shapeutil.h"

class R2Rect;
class S2PaddedCell;
class S2Shape;

// S2CrossingEdgeQuery is used to find edges or shapes that are crossed by
// an edge.  Here is an example showing how to index a set of polylines,
// and then find the polylines that are crossed by a given edge AB:
//
// void Test(vector<S2Polyline*> const& polylines,
//           S2Point const& a0, S2Point const &a1) {
//   S2ShapeIndex index;
//   for (S2Polyline* polyline : polylines) {
//     index.Add(new S2Polyline::Shape(polyline));
//   }
//   S2CrossingEdgeQuery query(index);
//   S2CrossingEdgeQuery::EdgeMap edge_map;
//   if (!query.GetCrossings(a, b, s2shapeutil::CrossingType::ALL, &edge_map)) {
//     return;
//   }
//   for (auto const& it : edge_map) {
//     S2Shape const* shape = it.first;
//     S2Polyline* polyline = polylines[shape->id()];
//     vector<int> const& edges = it.second;
//     for (int edge : edges) {
//       S2Point const& b0 = polyline->vertex(edge);
//       S2Point const& b1 = polyline->vertex(edge) + 1);
//       CHECK_GE(S2EdgeUtil::CrossingSign(a0, a1, b0, b1), 0);
//     }
//   }
// }
//
// Note that if you need to query many edges, it is more efficient to declare
// a single S2CrossingEdgeQuery object and reuse it so that temporary storage
// does not need to be reallocated each time.
//
// If you want to find *all* pairs of crossing edges, it is more efficient to
// call s2shapeutil::GetCrossingEdgePairs().
class S2CrossingEdgeQuery {
 public:
  // An EdgeMap stores a sorted set of edge ids for each shape.  Its type
  // may change, but can be treated as "map<S2Shape const*, vector<int>>"
  // except that the shapes are sorted by shape id rather than pointer (to
  // ensure repeatable behavior).
  //
  // For simple comparisons (pointers in this case), it is much faster to use
  // linear rather than binary search within btree nodes.
  struct CompareBtreeLinearSearch {
    using goog_btree_prefer_linear_node_search = std::true_type;
    bool operator()(S2Shape const* x, S2Shape const* y) const {
      return x->id() < y->id();
    }
  };
  using EdgeMap = util::btree::btree_map<S2Shape const*, std::vector<int>,
                                         CompareBtreeLinearSearch>;

  // Convenience constructor that calls Init().
  explicit S2CrossingEdgeQuery(S2ShapeIndex const& index);

  // Default constructor; requires Init() to be called.
  S2CrossingEdgeQuery();

  // REQUIRES: "index" is not modified after this method is called.
  void Init(S2ShapeIndex const& index);

  // Given a query edge AB and a shape S, return all the edges of S that
  // intersect AB.  If "type" is CrossingType::INTERIOR, then only
  // intersections at a point interior to both edges are reported, while if
  // "type" is CrossingType::ALL then edges that share a vertex are also
  // reported.  Returns false if no edges cross AB.
  using CrossingType = s2shapeutil::CrossingType;
  bool GetCrossings(S2Point const& a, S2Point const& b,
                    S2Shape const* shape, CrossingType type,
                    std::vector<int>* edges);

  // Given a query edge AB, return all the edges that intersect AB.  If "type"
  // is CrossingType::INTERIOR, then only intersections at a point interior to
  // both edges are reported, while if "type" is CrossingType::ALL then edges
  // that share a vertex are also reported.
  //
  // The edges are returned as a map from shapes to the edges of that shape
  // that intersect AB.  Every returned shape has at least one crossing edge.
  // Returns false if no edges intersect AB.
  bool GetCrossings(S2Point const& a, S2Point const& b, CrossingType type,
                    EdgeMap* edge_map);

  /////////////////////////// Low-Level Methods ////////////////////////////
  //
  // Most clients will not need the following methods.  They can be slightly
  // more efficient but are harder to use, since they require the client to do
  // all the actual crossing tests.

  // Given a query edge AB and a shape S, return a superset of the edges of
  // S that intersect AB.  Returns false if "edges" is empty.
  bool GetCandidates(S2Point const& a, S2Point const& b, S2Shape const* shape,
                     std::vector<int>* edges);

  // Given a query edge AB, return a map from shapes to a superset of the
  // edges for that shape that intersect AB.  Returns false if no shapes
  // intersect AB.
  //
  // CAVEAT: This method may return shapes that have an empty set of candidate
  // edges, i.e. (*edge_map)[shape].size() == 0.  However the return value is
  // true only if at least one shape has a candidate edge.
  bool GetCandidates(S2Point const& a, S2Point const& b, EdgeMap* edge_map);

  // Given a query edge AB and a cell "root", return all S2ShapeIndex cells
  // within "root" that might contain edges intersecting AB.
  bool GetCells(S2Point const& a, S2Point const& b, S2PaddedCell const& root,
                std::vector<S2ShapeIndexCell const*>* cells);

 private:
  // Internal methods are documented with their definitions.
  void GetCells(S2Point const& a, S2Point const& b);
  void GetCells(S2PaddedCell const& pcell, R2Rect const& bound);
  void ClipVAxis(R2Rect const& edge_bound, double center, int i,
                 S2PaddedCell const& pcell);
  void SplitUBound(R2Rect const& edge_bound, double u,
                   R2Rect child_bounds[2]) const;
  void SplitVBound(R2Rect const& edge_bound, double v,
                   R2Rect child_bounds[2]) const;
  static void SplitBound(R2Rect const& edge_bound, int u_end, double u,
                         int v_end, double v, R2Rect child_bounds[2]);

  S2ShapeIndex const* index_;

  // Temporary storage used while processing a query.
  R2Point a_;
  R2Point b_;
  S2ShapeIndex::Iterator iter_;

  // This is a private field rather than a local variable to reduce memory
  // allocation when a single S2CrossingEdgeQuery object is queried many times.
  absl::InlinedVector<S2ShapeIndexCell const*, 8> cells_;

  S2CrossingEdgeQuery(S2CrossingEdgeQuery const&) = delete;
  void operator=(S2CrossingEdgeQuery const&) = delete;
};

//////////////////   Implementation details follow   ////////////////////

inline S2CrossingEdgeQuery::S2CrossingEdgeQuery() : index_(nullptr) {}
inline S2CrossingEdgeQuery::S2CrossingEdgeQuery(S2ShapeIndex const& index) {
  Init(index);
}
inline void S2CrossingEdgeQuery::Init(S2ShapeIndex const& index) {
  index_ = &index;
  iter_.Init(index);
}

#endif  // S2_S2CROSSINGEDGEQUERY_H_
