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

#include "s2/s2shapeutil_edge_wrap.h"

#include "absl/log/absl_check.h"
#include "s2/s2point.h"
#include "s2/s2shape.h"

namespace s2shapeutil {

int NextEdgeWrap(const S2Shape& shape, int edge_id) {
  ABSL_DCHECK_GE(edge_id, 0);
  ABSL_DCHECK_LT(edge_id, shape.num_edges());

  S2Shape::ChainPosition chainpos = shape.chain_position(edge_id);
  S2Shape::Chain chaininfo = shape.chain(chainpos.chain_id);

  // Polygon chains wrap around, point and polylines don't.
  int offset = chainpos.offset;
  switch (shape.dimension()) {
    case 2:
      // Polygon chains always wrap around.
      offset = (offset + 1) % chaininfo.length;
      break;

    case 1:
      // If we're at the end of a polyline, wrap around if it's closed.
      if (offset == chaininfo.length - 1) {
        const S2Shape::Edge& curr = shape.chain_edge(chainpos.chain_id, offset);
        const S2Shape::Edge& next = shape.chain_edge(chainpos.chain_id, 0);

        if (curr.v1 == next.v0) {
          offset = 0;
        } else {
          return -1;
        }
      } else {
        ++offset;
      }
      break;

    default:
      // Points are one per chain.
      return -1;
  }

  return chaininfo.start + offset;
}

int PrevEdgeWrap(const S2Shape& shape, int edge_id) {
  ABSL_DCHECK_GE(edge_id, 0);
  ABSL_DCHECK_LT(edge_id, shape.num_edges());

  S2Shape::ChainPosition chainpos = shape.chain_position(edge_id);
  S2Shape::Chain chaininfo = shape.chain(chainpos.chain_id);

  int offset = chainpos.offset;
  switch (shape.dimension()) {
    case 2:
      // Polygons always wrap around.
      --offset;
      if (offset < 0) {
        offset += chaininfo.length;
      }
      break;

    case 1:
      if (offset == 0) {
        int end = chaininfo.length - 1;
        const S2Shape::Edge& curr = shape.chain_edge(chainpos.chain_id, 0);
        const S2Shape::Edge& prev = shape.chain_edge(chainpos.chain_id, end);

        if (prev.v1 == curr.v0) {
          offset = end;
        } else {
          return -1;
        }
      } else {
        --offset;
      }
      break;

    default:
      // Points are one per chain.
      return -1;
  }

  return chaininfo.start + offset;
}

}  //  namespace s2shapeutil
