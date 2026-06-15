// Copyright 2024 Google Inc. All Rights Reserved.
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

#include "s2/s2cellid_range.h"

#include "absl/log/absl_check.h"

int64_t S2CellIdRange::size() const {
  return end.distance_from_begin() - begin.distance_from_begin();
}

S2CellId S2CellIdRange::at(int64_t i) const {
  ABSL_DCHECK_GE(i, 0);
  ABSL_DCHECK_LT(i, size());
  return begin.advance(i);
}

S2CellIdRange S2CellIdRange::slice(int64_t start, int64_t stop) const {
  ABSL_DCHECK_GE(start, 0);
  ABSL_DCHECK_LE(start, stop);
  ABSL_DCHECK_LE(stop, size());
  return S2CellIdRange{begin.advance(start), begin.advance(stop)};
}

bool S2CellIdRange::contains(S2CellId cell) const {
  return cell >= begin && cell < end && cell.level() == begin.level();
}

std::optional<S2CellId> S2CellIdForwardIterator::next() {
  if (cur == end) return std::nullopt;
  S2CellId result = cur;
  cur = cur.next();
  return result;
}

std::optional<S2CellId> S2CellIdReverseIterator::next() {
  if (done) return std::nullopt;
  S2CellId result = cur;
  if (cur == begin) {
    done = true;
  } else {
    cur = cur.prev();
  }
  return result;
}
