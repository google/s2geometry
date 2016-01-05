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

#include "s2shapeutil.h"

#include "base/stringprintf.h"
#include "s2error.h"
#include "s2loop.h"
#include "s2shapeindex.h"

using std::pair;
using std::vector;

namespace s2shapeutil {

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
inline static bool FindSelfIntersection(vector<S2Loop*> const& loops,
                                        S2ShapeIndexCell const& cell,
                                        S2Error* error) {
  for (int a = 0; a < cell.num_shapes(); ++a) {
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
static bool FindLoopCrossing(vector<S2Loop*> const& loops,
                             S2ShapeIndexCell const& cell, S2Error* error) {
  // Possible optimization:
  // Sort the ClippedShapes by edge count to reduce the number of calls to
  // S2::Sign.  If n is the total number of shapes in the cell, n_i is
  // the number of edges in shape i, and c_i is the number of continuous
  // chains formed by these edges, the total number of calls is
  //
  //   sum(n_i * (1 + c_j + n_j), i=0..n-2, j=i+1..n-1)
  //
  // So for example if n=2, shape 0 has one chain of 1 edge, and shape 1 has
  // one chain of 8 edges, the number of calls to Sign is 1*10=10 if the
  // shapes are sorted by edge count, and 8*3=24 otherwise.
  //
  // typedef pair<int, S2ShapeIndex::ClippedShape const*> SortedShape;
  // vector<SortedShape> sorted_shapes;

  for (int a = 0; a < cell.num_shapes() - 1; ++a) {
    S2ClippedShape const& a_clipped = cell.clipped(a);
    S2Loop const& a_loop = *loops[a_clipped.shape_id()];
    int a_num_clipped = a_clipped.num_edges();
    for (int i = 0; i < a_num_clipped; ++i) {
      int ai = a_clipped.edge(i);
      S2EdgeUtil::EdgeCrosser crosser(&a_loop.vertex(ai), &a_loop.vertex(ai+1));
      for (int b = a + 1; b < cell.num_shapes(); ++b) {
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
    if (FindSelfIntersection(it.cell()->clipped(0), loop, error)) {
      return true;
    }
  }
  return false;
}

bool FindAnyCrossing(S2ShapeIndex const& index, vector<S2Loop*> const& loops,
                     S2Error* error) {
  for (S2ShapeIndex::Iterator it(index); !it.Done(); it.Next()) {
    if (FindSelfIntersection(loops, *it.cell(), error)) {
      return true;
    }
    if (it.cell()->num_shapes() >= 2 &&
        FindLoopCrossing(loops, *it.cell(), error)) {
      return true;
    }
  }
  return false;
}

}  // namespace s2shapeutil
