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

#ifndef S2_S2SHAPEUTIL_EDGE_WRAP_H_
#define S2_S2SHAPEUTIL_EDGE_WRAP_H_

#include "s2/s2shape.h"

namespace s2shapeutil {

// Convenience function that returns the edge id of next edge in a chain.  Wraps
// around at the start/end of any closed chains.
//
// This is intended for one-off lookups, as it has to look up the chain for
// the edge every time.  If you want many lookups or to iterate the edges of a
// chain, then it's better to do that directly.
//
// Returns -1 when the end of an open chain is reached. Polygon and closed
// polyline chains wrap around to the beginning and thus never return -1,
// while points always do.
int NextEdgeWrap(const S2Shape& shape, int edge_id);

// Convenience function that returns the edge id of previous edge in a chain.
// Wraps around at the start/end of any closed chains.
//
// This is intended for one-off lookups, as it has to look up the chain for
// the edge every time.  If you want many lookups or to iterate the edges of a
// chain, then it's better to do that directly.
//
// Returns -1 when the start of an open chain is reached. Polygon and closed
// polyline chains wrap around to the end and thus never return -1, while
// points always do.
int PrevEdgeWrap(const S2Shape& shape, int edge_id);

}  // namespace s2shapeutil

#endif  // S2_S2SHAPEUTIL_EDGE_WRAP_H_
