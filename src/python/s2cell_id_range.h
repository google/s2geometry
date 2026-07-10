// Copyright 2026 Google LLC
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

// Helpers for S2CellId sequence iteration used by Python bindings.
//
// S2CellIdRange represents a contiguous range of S2CellIds at a fixed level,
// with forward and reverse iterators.  These types are exposed to pybind11 to
// support Python's sequence protocol (len, indexing, slicing, iteration,
// reversed) on S2CellId ranges returned by methods such as children() and
// cells_at_level().

#ifndef PYTHON_S2CELL_ID_RANGE_H_
#define PYTHON_S2CELL_ID_RANGE_H_

#include <cstdint>
#include <optional>

#include "s2/s2cell_id.h"

// A contiguous half-open range [begin, end) of S2CellIds at the same level
// along the Hilbert curve.  All operations are O(1).  The range may be empty
// (begin == end).
struct S2CellIdRange {
  S2CellId begin;
  S2CellId end;

  // Return the number of cells in the range.
  int64_t size() const;

  // Return the cell at 0-based offset i.  i must be in [0, size()); behavior
  // is undefined otherwise.
  S2CellId at(int64_t i) const;

  // Return the sub-range [start, stop).  start and stop must satisfy
  // 0 <= start <= stop <= size(); behavior is undefined otherwise.
  S2CellIdRange slice(int64_t start, int64_t stop) const;

  // Return true if cell is in the range and at the same level as begin.
  bool contains(S2CellId cell) const;
};

// Forward iterator over an S2CellIdRange.
struct S2CellIdForwardIterator {
  S2CellId cur;
  S2CellId end;

  // Return the next cell and advance, or std::nullopt if exhausted.
  std::optional<S2CellId> next();
};

// Reverse iterator over an S2CellIdRange.  Construct with cur = range.end;
// next() decrements cur before yielding, so the empty-range case (begin ==
// end) is handled without a separate flag.
struct S2CellIdReverseIterator {
  S2CellId cur;    // one-past-current; starts at range.end
  S2CellId begin;  // stop when cur reaches here

  // Return the next cell (in reverse order) and advance, or std::nullopt if
  // exhausted.
  std::optional<S2CellId> next();
};

#endif  // PYTHON_S2CELL_ID_RANGE_H_
