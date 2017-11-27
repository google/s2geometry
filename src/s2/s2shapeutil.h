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
//
// This file contains miscellaneous S2ShapeIndex utility functions and classes:
//
//  - VisitCrossings: finds all edge intersections between S2ShapeIndexes.
//
// And functions that are mainly useful internally in the S2 library:
//
//  - FindAnyCrossing: helper function for S2Loop and S2Polygon validation.

#ifndef S2_S2SHAPEUTIL_H_
#define S2_S2SHAPEUTIL_H_

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "s2/third_party/absl/types/span.h"
#include "s2/_fpcontractoff.h"
#include "s2/s2edge_vector_shape.h"
#include "s2/s2lax_loop_shape.h"
#include "s2/s2lax_polygon_shape.h"
#include "s2/s2lax_polyline_shape.h"
#include "s2/s2loop.h"
#include "s2/s2point_vector_shape.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"
#include "s2/s2shapeindex.h"
#include "s2/s2shapeutil_shape_edge.h"

class S2Error;

namespace s2shapeutil {

/////////////////////////////////////////////////////////////////////////////
////////////////////// Utility Functions and Classes ////////////////////////

// A parameter that controls the reporting of edge intersections.
//
//  - CrossingType::INTERIOR reports intersections that occur at a point
//    interior to both edges (i.e., not at a vertex).
//
//  - CrossingType::ALL reports all intersections, even those where two edges
//    intersect only because they share a common vertex.
//
//  - CrossingType::NON_ADJACENT reports all intersections except for pairs of
//    the form (AB, BC) where both edges are from the same S2ShapeIndex.
//    (This is an optimization for checking polygon validity.)
enum class CrossingType { INTERIOR, ALL, NON_ADJACENT };

// A function that is called with pairs of crossing edges.  The function may
// return false in order to request that the algorithm should be terminated,
// i.e. no further crossings are needed.
//
// "is_interior" indicates that the crossing is at a point interior to both
// edges (i.e., not at a vertex).  (The calling function already has this
// information and it is moderately expensive to recompute.)
using EdgePairVisitor = std::function<
  bool (const ShapeEdge& a, const ShapeEdge& b, bool is_interior)>;

// Visits all pairs of crossing edges in the given S2ShapeIndex, terminating
// early if the given EdgePairVisitor function returns false (in which case
// VisitCrossings returns false as well).  "type" indicates whether all
// crossings should be visited, or only interior crossings.
//
// CAVEAT: Crossings may be visited more than once.
bool VisitCrossings(const S2ShapeIndex& index, CrossingType type,
                    const EdgePairVisitor& visitor);

// Like the above, but visits all pairs of crossing edges where one edge comes
// from each S2ShapeIndex.
//
// REQUIRES: type != CrossingType::NON_ADJACENT (not supported)
bool VisitCrossings(const S2ShapeIndex& a, const S2ShapeIndex& b,
                    CrossingType type, const EdgePairVisitor& visitor);


/////////////////////////////////////////////////////////////////////////////
///////////////// Methods used internally by the S2 library /////////////////


// Given an S2ShapeIndex containing a single polygonal shape (e.g., an
// S2Polygon or S2Loop), return true if any loop has a self-intersection
// (including duplicate vertices) or crosses any other loop (including vertex
// crossings and duplicate edges) and set "error" to a human-readable error
// message.  Otherwise return false and leave "error" unchanged.
//
// This method is used to implement the FindValidationError methods of S2Loop
// and S2Polygon.
//
// TODO(ericv): Add an option to support S2LaxPolygonShape rules (i.e.,
// duplicate vertices and edges are allowed, but loop crossings are not).
bool FindAnyCrossing(const S2ShapeIndex& index, S2Error* error);


}  // namespace s2shapeutil

#endif  // S2_S2SHAPEUTIL_H_
