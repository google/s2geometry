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

#ifndef S2_S2SHAPEUTIL_COUNT_VERTICES_H_
#define S2_S2SHAPEUTIL_COUNT_VERTICES_H_

#include <cstdint>

#include "s2/s2shape.h"
#include "s2/s2shape_index.h"

namespace s2shapeutil {

// Returns the total number of vertices in a single shape.
int64_t CountVertices(const S2Shape& shape);

// Returns the total number of vertices in all indexed shapes. This method takes
// time linear in the number of shapes.
int64_t CountVertices(const S2ShapeIndex& index);

}  // namespace s2shapeutil

#endif  // S2_S2SHAPEUTIL_COUNT_VERTICES_H_
