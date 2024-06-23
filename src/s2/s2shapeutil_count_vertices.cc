// Copyright Google Inc. All Rights Reserved.
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

#include "s2/s2shapeutil_count_vertices.h"

#include <cstdint>

#include "absl/log/absl_check.h"
#include "s2/s2shape.h"
#include "s2/s2shape_index.h"

namespace s2shapeutil {

int64_t CountVertices(const S2Shape& shape) {
  int dimension = shape.dimension();
  ABSL_DCHECK_GE(dimension, 0);
  ABSL_DCHECK_LE(dimension, 2);
  switch (dimension) {
    case 0:
      // Dimension 0 shapes contain points, where each point is a degenerate
      // edge and each edge is in a separate chain.  Because of this,
      // num_edges() == num_chains() and we can use either for the number
      // of points.
      return shape.num_chains();
    case 1:
      // Dimension 1 shapes contain polylines, where each polyline is a
      // chain containing one or more edges.  Polylines aren't closed,
      // so they always have a form similar to: [AB, BC, CD] where the start
      // and end points are different.
      //
      // We can see the number of vertices is num_edges() + 1.  Summing over
      // all chains, we see the number of vertices is just the number of
      // edges + the number of chains (+1 for each).
      return shape.num_edges() + shape.num_chains();
    default:
      // Dimension 2 shapes contain polygons, which are closed, so they have
      // the form [AB, BC, CD, DA], which is even simpler than the
      // one-dimensional case above, the number of unique vertices simply
      // being the number of edges.
      return shape.num_edges();
  }
}

int64_t CountVertices(const S2ShapeIndex& index) {
  int64_t vertices = 0;
  for (const S2Shape* shape : index) {
    vertices += CountVertices(*shape);
  }
  return vertices;
}

}  // namespace s2shapeutil
