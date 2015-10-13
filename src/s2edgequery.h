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

#ifndef S2_GEOMETRY_S2EDGEQUERY_H_
#define S2_GEOMETRY_S2EDGEQUERY_H_

#include <vector>

#include "base/macros.h"
#include <map>
#include "fpcontractoff.h"
#include "r2.h"
#include "r2rect.h"
#include "s2.h"
#include "s2shapeindex.h"
#include "util/gtl/inlined_vector.h"

class R2Rect;
class S2PaddedCell;
class S2Shape;

// S2EdgeQuery is used to find edges or shapes that are crossed by an edge.
// If you need to query many edges, it is more efficient to declare a single
// S2EdgeQuery object and reuse it so that temporary storage does not need to
// be reallocated each time.
class S2EdgeQuery {
 public:
  // An EdgeMap stores a sorted set of edge ids for each shape.
  typedef std::map<S2Shape const*, std::vector<int> > EdgeMap;

  // Convenience constructor that calls Init().
  explicit S2EdgeQuery(S2ShapeIndex const& index);

  // Default constructor; requires Init() to be called.
  S2EdgeQuery();

  // REQUIRES: "index" is not modified after this method is called.
  void Init(S2ShapeIndex const& index);

  // Given a query edge AB and a shape S, return a superset of the edges of
  // S that intersect AB.  Returns false if "edges" is empty.
  bool GetCandidates(S2Point const& a, S2Point const& b,
                     S2Shape const* shape, std::vector<int>* edges);

  // Given a query edge AB, return a map from shapes to a superset of the
  // edges for that shape that intersect AB.  Returns false if no shapes
  // intersect AB.
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
  // allocation when a single S2EdgeQuery object is queried many times.
  util::gtl::InlinedVector<S2ShapeIndexCell const*, 8> cells_;

  DISALLOW_COPY_AND_ASSIGN(S2EdgeQuery);
};


//////////////////   Implementation details follow   ////////////////////


inline S2EdgeQuery::S2EdgeQuery()
    : index_(NULL) {
}
inline S2EdgeQuery::S2EdgeQuery(S2ShapeIndex const& index) {
  Init(index);
}
inline void S2EdgeQuery::Init(S2ShapeIndex const& index) {
  index_ = &index;
  iter_.Init(index);
}

#endif  // S2_GEOMETRY_S2EDGEQUERY_H_
