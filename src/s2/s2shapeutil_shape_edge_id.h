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

#ifndef S2_S2SHAPEUTIL_SHAPE_EDGE_ID_H_
#define S2_S2SHAPEUTIL_SHAPE_EDGE_ID_H_

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <ostream>

namespace s2shapeutil {

// ShapeEdgeId is a unique identifier for an edge within an S2ShapeIndex,
// consisting of a (shape_id, edge_id) pair.  It is similar to
// std::pair<int32_t, int32_t> except that it has named fields.
// It should be passed and returned by value.
struct ShapeEdgeId {
 public:
  int32_t shape_id, edge_id;
  ShapeEdgeId() : shape_id(-1), edge_id(-1) {}
  ShapeEdgeId(int32_t _shape_id, int32_t _edge_id);

  bool operator==(ShapeEdgeId other) const;
  bool operator!=(ShapeEdgeId other) const;
  bool operator<(ShapeEdgeId other) const;
  bool operator>(ShapeEdgeId other) const;
  bool operator<=(ShapeEdgeId other) const;
  bool operator>=(ShapeEdgeId other) const;

  template <typename H>
  friend H AbslHashValue(H h, ShapeEdgeId id) {
    return H::combine(std::move(h), id.shape_id, id.edge_id);
  }
};
std::ostream& operator<<(std::ostream& os, ShapeEdgeId id);


//////////////////   Implementation details follow   ////////////////////

inline ShapeEdgeId::ShapeEdgeId(int32_t _shape_id, int32_t _edge_id)
    : shape_id(_shape_id), edge_id(_edge_id) {}

inline bool ShapeEdgeId::operator==(ShapeEdgeId other) const {
  return shape_id == other.shape_id && edge_id == other.edge_id;
}

inline bool ShapeEdgeId::operator!=(ShapeEdgeId other) const {
  return !(*this == other);
}

inline bool ShapeEdgeId::operator<(ShapeEdgeId other) const {
  if (shape_id < other.shape_id) return true;
  if (other.shape_id < shape_id) return false;
  return edge_id < other.edge_id;
}

inline bool ShapeEdgeId::operator>(ShapeEdgeId other) const {
  return other < *this;
}

inline bool ShapeEdgeId::operator<=(ShapeEdgeId other) const {
  return !(other < *this);
}

inline bool ShapeEdgeId::operator>=(ShapeEdgeId other) const {
  return !(*this < other);
}

inline std::ostream& operator<<(std::ostream& os, ShapeEdgeId id) {
  return os << id.shape_id << ":" << id.edge_id;
}

}  // namespace s2shapeutil

#endif  // S2_S2SHAPEUTIL_SHAPE_EDGE_ID_H_
