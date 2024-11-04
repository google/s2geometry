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

#include "s2/s2shapeutil_build_polygon_boundaries.h"

#include <memory>
#include <utility>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/log/absl_check.h"
#include "absl/types/span.h"
#include "s2/mutable_s2shape_index.h"
#include "s2/s2contains_point_query.h"
#include "s2/s2pointutil.h"
#include "s2/s2shape.h"
#include "s2/s2shapeutil_contains_brute_force.h"

using absl::WrapUnique;
using std::vector;

namespace s2shapeutil {

void BuildPolygonBoundaries(absl::Span<const vector<S2Shape*>> components,
                            vector<vector<S2Shape*>>* polygons) {
  polygons->clear();
  if (components.empty()) return;

  // Since the loop boundaries do not cross, a loop nesting hierarchy can be
  // defined by choosing any point on the sphere as the "point at infinity".
  // Loop A then contains loop B if (1) A contains the boundary of B and (2)
  // loop A does not contain the point at infinity.
  //
  // We choose S2::Origin() for this purpose.  The loop nesting hierarchy then
  // determines the face structure.  Here are the details:
  //
  // 1. Build an S2ShapeIndex of all loops that do not contain S2::Origin().
  //    This leaves at most one unindexed loop per connected component
  //    (the "outer loop").
  //
  // 2. For each component, choose a representative vertex and determine
  //    which indexed loops contain it.  The "depth" of this component is
  //    defined as the number of such loops.
  //
  // 3. Assign the outer loop of each component to the containing loop whose
  //    depth is one less.  This generates a set of multi-loop polygons.
  //
  // 4. The outer loops of all components at depth 0 become a single face.

  MutableS2ShapeIndex index;
  // A map from shape id to the corresponding component number.
  vector<int> component_ids;
  vector<S2Shape*> outer_loops;
  for (size_t i = 0; i < components.size(); ++i) {
    const auto& component = components[i];
    for (S2Shape* loop : component) {
      if (component.size() > 1 &&
          !s2shapeutil::ContainsBruteForce(*loop, S2::Origin())) {
        // Ownership is transferred back at the end of this function.
        index.Add(WrapUnique(loop));
        component_ids.push_back(i);
      } else {
        outer_loops.push_back(loop);
      }
    }
    // Check that there is exactly one outer loop in each component.
    ABSL_DCHECK_EQ(i + 1, outer_loops.size())
        << "Component is not a subdivision";
  }
  // Find the loops containing each component.
  vector<vector<int>> ancestors(components.size());
  auto contains_query = MakeS2ContainsPointQuery(&index);
  for (size_t i = 0; i < outer_loops.size(); ++i) {
    auto loop = outer_loops[i];
    ABSL_DCHECK_GT(loop->num_edges(), 0);
    ancestors[i] = contains_query.GetContainingShapeIds(loop->edge(0).v0);
  }
  // Assign each outer loop to the component whose depth is one less.
  // Components at depth 0 become a single face.
  absl::flat_hash_map<const S2Shape*, vector<S2Shape*>> children;
  for (size_t i = 0; i < outer_loops.size(); ++i) {
    int ancestor_id = -1;
    size_t depth = ancestors[i].size();
    if (depth > 0) {
      for (auto candidate_id : ancestors[i]) {
        if (ancestors[component_ids[candidate_id]].size() == depth - 1) {
          ABSL_DCHECK_EQ(ancestor_id, -1);
          ancestor_id = candidate_id;
        }
      }
      ABSL_DCHECK_GE(ancestor_id, 0);
    }
    const S2Shape* ancestor =
        (ancestor_id < 0) ? nullptr : index.shape(ancestor_id);
    children[ancestor].push_back(outer_loops[i]);
  }
  // There is one face per loop that is not an outer loop, plus one for the
  // outer loops of components at depth 0.
  polygons->resize(index.num_shape_ids() + 1);
  for (int i = 0; i < index.num_shape_ids(); ++i) {
    auto polygon = &(*polygons)[i];
    auto loop = index.shape(i);
    auto itr = children.find(loop);
    if (itr != children.end()) {
      *polygon = itr->second;
    }
    polygon->push_back(const_cast<S2Shape*>(loop));
  }
  polygons->back() = children[nullptr];

  // Explicitly release the shapes from the index so they are not deleted.
  for (auto& ptr : index.ReleaseAll()) ptr.release();
}

}  // namespace s2shapeutil
