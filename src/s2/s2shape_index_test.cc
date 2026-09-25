// Copyright 2012 Google Inc. All Rights Reserved.
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

#include "s2/s2shape_index.h"

#include <cstdint>
#include <initializer_list>
#include <string>

#include <gtest/gtest.h>
#include "s2/s2shape.h"
#include "s2/s2text_format.h"
#include "s2/util/coding/coder.h"

// TODO(ericv): Add tests for S2ShapeIndexCell and S2ClippedShape.
// (Currently these are tested indirectly by MutableS2ShapeIndex.)
// Also test the base Iterator type (which wraps another iterator).

namespace {

// Returns the body of an encoded S2ShapeIndexCell holding the given varint
// values.  Encoding them with Encoder keeps the test bytes in the same form
// that S2ShapeIndexCell::Encode() produces, without hand-encoding varints.
std::string CellBytes(std::initializer_list<uint64_t> values) {
  Encoder encoder;
  encoder.Ensure(values.size() * Encoder::kVarintMax64);
  for (uint64_t value : values) encoder.put_varint64(value);
  return std::string(encoder.base(), encoder.length());
}

// Decodes a cell whose body is exactly "values", and reports whether Decode()
// accepted it.
bool DecodeCell(S2ShapeIndexCell* cell, int num_shape_ids,
                std::initializer_list<uint64_t> values) {
  std::string bytes = CellBytes(values);
  Decoder decoder(bytes.data(), bytes.size());
  return cell->Decode(num_shape_ids, &decoder);
}

}  // namespace

TEST(S2ShapeIndexCell, DecodeRejectsEdgeCountOverflow) {
  // A single-shape cell in the "general case" encoding (bits 0-1 == 3) whose
  // num_edges field is 2**31, one more than the largest value an int32 can
  // hold.  Truncating it to int and then widening it to uint32_t for
  // S2ClippedShape::Init() used to understate the edge count, so the cell was
  // accepted with a bogus number of edges (and, for other out-of-range values,
  // requested a multi-gigabyte allocation).
  constexpr uint64_t kNumEdges = uint64_t{1} << 31;
  S2ShapeIndexCell cell;
  EXPECT_FALSE(DecodeCell(&cell, 1, {(kNumEdges << 3) | 3}));
}

TEST(S2ShapeIndexCell, DecodeRejectsTooManyClippedShapes) {
  // num_clipped == 3 in an index that only has two shapes.  Each clipped shape
  // is encoded with no edges and a shape delta of 0, so every shape id is in
  // range and only the clipped shape count is invalid; before the count was
  // checked, the cell decoded successfully with shape_id 0 repeated 3 times.
  constexpr uint64_t kNumClipped = 3;
  S2ShapeIndexCell cell;
  EXPECT_FALSE(DecodeCell(&cell, 2, {(kNumClipped << 3) | 3, 7, 7, 7}));
  EXPECT_EQ(0, cell.num_clipped());
}

TEST(S2ShapeIndexCell, DecodeRejectsShapeIdPastNumShapeIds) {
  // One malformed cell per branch that accumulates a shape id.  Each decodes
  // to shape id 5 in an index that only has two shapes (ids 0 and 1), which
  // used to be accepted and later indexed the shapes out of bounds.

  // Contiguous edge range: tag bit 0 == 0, and the shape delta is in the high
  // bits of the value that follows.
  S2ShapeIndexCell contiguous;
  EXPECT_FALSE(DecodeCell(&contiguous, 2, {0, 5 << 4}));

  // A clipped shape with no edges: tag 7, shape delta in the high bits.
  S2ShapeIndexCell no_edges;
  EXPECT_FALSE(DecodeCell(&no_edges, 2, {(uint64_t{5} << 4) | 7}));

  // A general edge list: tag 1, then the shape delta, then a single edge.
  S2ShapeIndexCell edges;
  EXPECT_FALSE(DecodeCell(&edges, 2, {1, 5, 0}));
}

TEST(S2ShapeIndexIterator, PrefixIncrement) {
  const auto index =
      s2textformat::MakeIndexOrDie("1:1 # 1:1, 2:2 # 2:2, 1:3, 2:4, 3:2");

  int count = 0;
  for (auto iter = index->begin(); iter != index->end(); ++iter) {
    ++count;
  }
  EXPECT_EQ(count, 3);
}

TEST(S2ShapeIndexIterator, PostfixIncrement) {
  const auto index =
      s2textformat::MakeIndexOrDie("1:1 # 1:1, 2:2 # 2:2, 1:3, 2:4, 3:2");

  int count = 0;
  for (auto iter = index->begin(); iter != index->end(); iter++) {
    ++count;
  }
  EXPECT_EQ(count, 3);
}
