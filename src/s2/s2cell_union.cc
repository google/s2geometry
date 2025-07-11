// Copyright 2005 Google Inc. All Rights Reserved.
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

#include "s2/s2cell_union.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <string>
#include <utility>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/log/absl_check.h"
#include "absl/strings/str_cat.h"

#include "s2/base/commandlineflags.h"
#include "s2/s1angle.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2latlng_rect.h"
#include "s2/s2metrics.h"
#include "s2/s2point.h"
#include "s2/util/coding/coder.h"

using std::is_sorted;
using std::max;
using std::min;
using std::string;
using std::vector;

S2_DEFINE_int32(s2cell_union_decode_max_num_cells, 1000000,
                "The maximum number of cells allowed by S2CellUnion::Decode");

static const unsigned char kCurrentLosslessEncodingVersionNumber = 1;

vector<S2CellId> S2CellUnion::ToS2CellIds(const vector<uint64_t>& ids) {
  vector<S2CellId> cell_ids;
  cell_ids.reserve(ids.size());
  for (auto id : ids) cell_ids.push_back(S2CellId(id));
  return cell_ids;
}

S2CellUnion::S2CellUnion(const vector<uint64_t>& cell_ids)
    : cell_ids_(ToS2CellIds(cell_ids)) {
  Normalize();
}

S2CellUnion S2CellUnion::WholeSphere() {
  return S2CellUnion({S2CellId::FromFace(0), S2CellId::FromFace(1),
                      S2CellId::FromFace(2), S2CellId::FromFace(3),
                      S2CellId::FromFace(4), S2CellId::FromFace(5)});
}

S2CellUnion S2CellUnion::FromMinMax(S2CellId min_id, S2CellId max_id) {
  S2CellUnion result;
  result.InitFromMinMax(min_id, max_id);
  return result;
}

S2CellUnion S2CellUnion::FromBeginEnd(S2CellId begin, S2CellId end) {
  S2CellUnion result;
  result.InitFromBeginEnd(begin, end);
  return result;
}

void S2CellUnion::Init(const vector<uint64_t>& cell_ids) {
  cell_ids_ = ToS2CellIds(cell_ids);
  Normalize();
}

void S2CellUnion::InitFromMinMax(S2CellId min_id, S2CellId max_id) {
  ABSL_DCHECK(max_id.is_valid()) << max_id;
  InitFromBeginEnd(min_id, max_id.next());
}

void S2CellUnion::InitFromBeginEnd(S2CellId begin, S2CellId end) {
  ABSL_DCHECK(begin.is_leaf()) << begin;
  ABSL_DCHECK(end.is_leaf()) << end;
  const S2CellId kLeafEnd = S2CellId::End(S2CellId::kMaxLevel);
  ABSL_DCHECK(begin.is_valid() || begin == kLeafEnd) << begin;
  ABSL_DCHECK(end.is_valid() || end == kLeafEnd) << end;
  ABSL_DCHECK_LE(begin, end);

  // We repeatedly add the largest cell we can.
  cell_ids_.clear();
  for (S2CellId id = begin.maximum_tile(end);
       id != end; id = id.next().maximum_tile(end)) {
    cell_ids_.push_back(id);
  }
  // The output is already normalized.
  ABSL_DCHECK(IsNormalized());
}

void S2CellUnion::Pack(int excess) {
  if (cell_ids_.capacity() - cell_ids_.size() > static_cast<size_t>(excess)) {
    cell_ids_.shrink_to_fit();
  }
}

S2CellUnion* S2CellUnion::Clone() const {
  return new S2CellUnion(cell_ids_, VERBATIM);
}

// Returns true if the given four cells have a common parent.
// REQUIRES: The four cells are distinct.
inline static bool AreSiblings(S2CellId a, S2CellId b, S2CellId c, S2CellId d) {
  // A necessary (but not sufficient) condition is that the XOR of the
  // four cells must be zero.  This is also very fast to test.
  if ((a.id() ^ b.id() ^ c.id()) != d.id()) return false;

  // Now we do a slightly more expensive but exact test.  First, compute a
  // mask that blocks out the two bits that encode the child position of
  // "id" with respect to its parent, then check that the other three
  // children all agree with "mask".
  uint64_t mask = d.lsb() << 1;
  mask = ~(mask + (mask << 1));
  uint64_t id_masked = (d.id() & mask);
  return ((a.id() & mask) == id_masked &&
          (b.id() & mask) == id_masked &&
          (c.id() & mask) == id_masked &&
          !d.is_face());
}

bool S2CellUnion::IsValid() const {
  if (num_cells() > 0 && !cell_id(0).is_valid()) return false;
  for (int i = 1; i < num_cells(); ++i) {
    if (!cell_id(i).is_valid()) return false;
    if (cell_id(i - 1).range_max() >= cell_id(i).range_min()) return false;
  }
  return true;
}

bool S2CellUnion::IsNormalized() const {
  if (num_cells() > 0 && !cell_id(0).is_valid()) return false;
  for (int i = 1; i < num_cells(); ++i) {
    if (!cell_id(i).is_valid()) return false;
    if (cell_id(i - 1).range_max() >= cell_id(i).range_min()) return false;
    if (i >= 3 && AreSiblings(cell_id(i - 3), cell_id(i - 2),
                              cell_id(i - 1), cell_id(i))) {
      return false;
    }
  }
  return true;
}

void S2CellUnion::Normalize() {
  Normalize(&cell_ids_);
}

/*static*/ void S2CellUnion::Normalize(vector<S2CellId>* ids) {
  // Optimize the representation by discarding cells contained by other cells,
  // and looking for cases where all subcells of a parent cell are present.

  std::sort(ids->begin(), ids->end());
  size_t out = 0;
  for (S2CellId id : *ids) {
    ABSL_DCHECK(id.is_valid()) << id;
    // Check whether this cell is contained by the previous cell.
    if (out > 0 && (*ids)[out-1].contains(id)) continue;

    // Discard any previous cells contained by this cell.
    while (out > 0 && id.contains((*ids)[out-1])) --out;

    // Check whether the last 3 elements plus "id" can be collapsed into a
    // single parent cell.
    while (out >= 3 && AreSiblings((*ids)[out - 3], (*ids)[out - 2],
                                   (*ids)[out - 1], id)) {
      // Replace four children by their parent cell.
      id = id.parent();
      out -= 3;
    }
    (*ids)[out++] = id;
  }
  if (ids->size() != out)
    ids->resize(out);
}

void S2CellUnion::Denormalize(int min_level, int level_mod,
                              vector<S2CellId>* out) const {
  Denormalize(cell_ids_, min_level, level_mod, out);
}

void S2CellUnion::Denormalize(const vector<S2CellId>& in,
                              int min_level, int level_mod,
                              vector<S2CellId>* out) {
  ABSL_DCHECK_GE(min_level, 0);
  ABSL_DCHECK_LE(min_level, S2CellId::kMaxLevel);
  ABSL_DCHECK_GE(level_mod, 1);
  ABSL_DCHECK_LE(level_mod, 3);
  ABSL_DCHECK_NE(out, &in);

  out->clear();
  out->reserve(in.size());
  for (S2CellId id : in) {
    int level = id.level();
    int new_level = max(min_level, level);
    if (level_mod > 1) {
      // Round up so that (new_level - min_level) is a multiple of level_mod.
      // (Note that S2CellId::kMaxLevel is a multiple of 1, 2, and 3.)
      new_level += (S2CellId::kMaxLevel - (new_level - min_level)) % level_mod;
      new_level = min(S2CellId::kMaxLevel, new_level);
    }
    if (new_level == level) {
      out->push_back(id);
    } else {
      S2CellId end = id.child_end(new_level);
      for (id = id.child_begin(new_level); id != end; id = id.next()) {
        out->push_back(id);
      }
    }
  }
}

S2Cap S2CellUnion::GetCapBound() const {
  // Compute the approximate centroid of the region.  This won't produce the
  // bounding cap of minimal area, but it should be close enough.
  if (cell_ids_.empty()) return S2Cap::Empty();
  S2Point centroid(0, 0, 0);
  for (S2CellId id : *this) {
    double area = S2Cell::AverageArea(id.level());
    centroid += area * id.ToPoint();
  }
  if (centroid == S2Point(0, 0, 0)) {
    centroid = S2Point(1, 0, 0);
  } else {
    centroid = centroid.Normalize();
  }

  // Use the centroid as the cap axis, and expand the cap angle so that it
  // contains the bounding caps of all the individual cells.  Note that it is
  // *not* sufficient to just bound all the cell vertices because the bounding
  // cap may be concave (i.e. cover more than one hemisphere).
  S2Cap cap = S2Cap::FromPoint(centroid);
  for (S2CellId id : *this) {
    cap.AddCap(S2Cell(id).GetCapBound());
  }
  return cap;
}

S2LatLngRect S2CellUnion::GetRectBound() const {
  S2LatLngRect bound = S2LatLngRect::Empty();
  for (S2CellId id : *this) {
    bound = bound.Union(S2Cell(id).GetRectBound());
  }
  return bound;
}

void S2CellUnion::GetCellUnionBound(vector<S2CellId>* cell_ids) const {
  // TODO(user): Try returning the cells in the covering.  The comments
  // in s2region.h say `*cell_ids` should be small, however.
  GetCapBound().GetCellUnionBound(cell_ids);
}

// Returns true if "a" lies entirely before "b" on the Hilbert curve.  Note that
// this is not even a weak ordering, since incomparability is not transitive.
// Nevertheless, given a sorted vector of disjoint S2CellIds it can be used be
// used to find the first element that might intersect a given target S2CellId:
//
//   auto it = std::lower_bound(v.begin(), v.end(), target, EntirelyPrecedes);
//
// This works because std::lower_bound() only requires that the elements are
// partitioned with respect to the given predicate (which is true as long as the
// S2CellIds are sorted and disjoint).
static inline bool EntirelyPrecedes(S2CellId a, S2CellId b) {
  return a.range_max() < b.range_min();
}

bool S2CellUnion::Contains(S2CellId id) const {
  // This is an exact test.  Each cell occupies a linear span of the S2
  // space-filling curve, and the cell id is simply the position at the center
  // of this span.  The cell union ids are sorted in increasing order along
  // the space-filling curve.  So we simply find the first cell that might
  // intersect (is not entirely before) the target (using binary search).
  // There is containment if and only if this cell id contains the target id.
  ABSL_DCHECK(id.is_valid()) << id;

  const auto i = std::lower_bound(begin(), end(), id, EntirelyPrecedes);
  return i != end() && i->contains(id);
}

bool S2CellUnion::Intersects(S2CellId id) const {
  // This is an exact test; see the comments for Contains() above.
  ABSL_DCHECK(id.is_valid()) << id;

  const auto i = std::lower_bound(begin(), end(), id, EntirelyPrecedes);
  return i != end() && i->intersects(id);
}

bool S2CellUnion::Contains(const S2CellUnion& y) const {
  if (y.empty()) return true;
  if (empty()) return false;
  auto i = begin();
  for (S2CellId y_id : y) {
    // If our first cell ends before the one we need to contain, advance
    // where we start searching.
    if (EntirelyPrecedes(*i, y_id)) {
      i = std::lower_bound(i + 1, end(), y_id, EntirelyPrecedes);
      // If we're at the end, we don't contain the current y_id.
      if (i == end()) return false;
    }
    if (!i->contains(y_id)) return false;
  }
  return true;
}

bool S2CellUnion::Intersects(const S2CellUnion& y) const {
  // Walk along the two sorted vectors, looking for overlap.
  for (auto i = begin(), j = y.begin(); i != end() && j != y.end(); ) {
    if (EntirelyPrecedes(*i, *j)) {
      // Advance "i" to the first cell that might overlap *j.
      i = std::lower_bound(i + 1, end(), *j, EntirelyPrecedes);
      continue;
    }
    if (EntirelyPrecedes(*j, *i)) {
      // Advance "j" to the first cell that might overlap *i.
      j = std::lower_bound(j + 1, y.end(), *i, EntirelyPrecedes);
      continue;
    }
    // Neither cell is to the left of the other, so they must intersect.
    ABSL_DCHECK(i->intersects(*j));
    return true;
  }
  return false;
}

S2CellUnion S2CellUnion::Union(const S2CellUnion& y) const {
  vector<S2CellId> cell_ids;
  cell_ids.reserve(num_cells() + y.num_cells());
  cell_ids = cell_ids_;
  cell_ids.insert(cell_ids.end(), y.cell_ids_.begin(), y.cell_ids_.end());
  return S2CellUnion(std::move(cell_ids));
}

S2CellUnion S2CellUnion::Intersection(S2CellId id) const {
  ABSL_DCHECK(id.is_valid()) << id;
  S2CellUnion result;
  if (Contains(id)) {
    result.cell_ids_.push_back(id);
  } else {
    vector<S2CellId>::const_iterator i = std::lower_bound(
        cell_ids_.begin(), cell_ids_.end(), id.range_min());
    S2CellId id_max = id.range_max();
    while (i != cell_ids_.end() && *i <= id_max)
      result.cell_ids_.push_back(*i++);
  }
  ABSL_DCHECK(result.IsNormalized() || !IsNormalized());
  return result;
}

S2CellUnion S2CellUnion::Intersection(const S2CellUnion& y) const {
  S2CellUnion result;
  GetIntersection(cell_ids_, y.cell_ids_, &result.cell_ids_);
  // The output is normalized as long as both inputs are normalized.
  ABSL_DCHECK(result.IsNormalized() || !IsNormalized() || !y.IsNormalized());
  return result;
}

/*static*/ void S2CellUnion::GetIntersection(const vector<S2CellId>& x,
                                             const vector<S2CellId>& y,
                                             vector<S2CellId>* out) {
  ABSL_DCHECK_NE(out, &x);
  ABSL_DCHECK_NE(out, &y);
  ABSL_DCHECK(is_sorted(x.begin(), x.end()));
  ABSL_DCHECK(is_sorted(y.begin(), y.end()));

  // This is a fairly efficient calculation that uses binary search to skip
  // over sections of both input vectors.  It takes logarithmic time if all the
  // cells of "x" come before or after all the cells of "y" in S2CellId order.

  out->clear();
  vector<S2CellId>::const_iterator i = x.begin();
  vector<S2CellId>::const_iterator j = y.begin();
  while (i != x.end() && j != y.end()) {
    S2CellId imin = i->range_min();
    S2CellId jmin = j->range_min();
    if (imin > jmin) {
      // Either j->contains(*i) or the two cells are disjoint.
      if (*i <= j->range_max()) {
        out->push_back(*i++);
      } else {
        // Advance "j" to the first cell that might overlap *i.
        j = std::lower_bound(j + 1, y.end(), *i, EntirelyPrecedes);
      }
    } else if (jmin > imin) {
      // Identical to the code above with "i" and "j" reversed.
      if (*j <= i->range_max()) {
        out->push_back(*j++);
      } else {
        i = std::lower_bound(i + 1, x.end(), *j, EntirelyPrecedes);
      }
    } else {
      // "i" and "j" have the same range_min(), so one contains the other.
      if (*i < *j)
        out->push_back(*i++);
      else
        out->push_back(*j++);
    }
  }
  // The output is generated in sorted order.
  ABSL_DCHECK(is_sorted(out->begin(), out->end()));
}

static void GetDifferenceInternal(S2CellId cell,
                                  const S2CellUnion& y,
                                  vector<S2CellId>* cell_ids) {
  // Add the difference between cell and y to cell_ids.
  // If they intersect but the difference is non-empty, divide and conquer.
  if (!y.Intersects(cell)) {
    cell_ids->push_back(cell);
  } else if (!y.Contains(cell)) {
    S2CellId child = cell.child_begin();
    for (int i = 0; ; ++i) {
      GetDifferenceInternal(child, y, cell_ids);
      if (i == 3) break;  // Avoid unnecessary next() computation.
      child = child.next();
    }
  }
}

S2CellUnion S2CellUnion::Difference(const S2CellUnion& y) const {
  // TODO(ericv): this is approximately O(N*log(N)), but could probably
  // use similar techniques as GetIntersection() to be more efficient.

  S2CellUnion result;
  for (S2CellId id : *this) {
    GetDifferenceInternal(id, y, &result.cell_ids_);
  }
  // The output is normalized as long as the first argument is normalized.
  ABSL_DCHECK(result.IsNormalized() || !IsNormalized());
  return result;
}

void S2CellUnion::Expand(int expand_level) {
  vector<S2CellId> output;
  uint64_t level_lsb = S2CellId::lsb_for_level(expand_level);
  for (int i = num_cells(); --i >= 0; ) {
    S2CellId id = cell_id(i);
    if (id.lsb() < level_lsb) {
      id = id.parent(expand_level);
      // Optimization: skip over any cells contained by this one.  This is
      // especially important when very small regions are being expanded.
      while (i > 0 && id.contains(cell_id(i - 1))) --i;
    }
    output.push_back(id);
    id.AppendAllNeighbors(expand_level, &output);
  }
  Init(std::move(output));
}

void S2CellUnion::Expand(S1Angle min_radius, int max_level_diff) {
  int min_level = S2CellId::kMaxLevel;
  for (S2CellId id : *this) {
    min_level = min(min_level, id.level());
  }
  // Find the maximum level such that all cells are at least "min_radius" wide.
  int radius_level = S2::kMinWidth.GetLevelForMinValue(min_radius.radians());
  if (radius_level == 0 && min_radius.radians() > S2::kMinWidth.GetValue(0)) {
    // The requested expansion is greater than the width of a face cell.
    // The easiest way to handle this is to expand twice.
    Expand(0);
  }
  Expand(min(min_level + max_level_diff, radius_level));
}

uint64_t S2CellUnion::LeafCellsCovered() const {
  uint64_t num_leaves = 0;
  for (S2CellId id : *this) {
    const int inverted_level = S2CellId::kMaxLevel - id.level();
    num_leaves += (1ULL << (inverted_level << 1));
  }
  return num_leaves;
}

double S2CellUnion::AverageBasedArea() const {
  return S2Cell::AverageArea(S2CellId::kMaxLevel) * LeafCellsCovered();
}

double S2CellUnion::ApproxArea() const {
  double area = 0;
  for (S2CellId id : *this) {
    area += S2Cell(id).ApproxArea();
  }
  return area;
}

double S2CellUnion::ExactArea() const {
  double area = 0;
  for (S2CellId id : *this) {
    area += S2Cell(id).ExactArea();
  }
  return area;
}

bool operator==(const S2CellUnion& x, const S2CellUnion& y) {
  return x.cell_ids() == y.cell_ids();
}

bool operator!=(const S2CellUnion& x, const S2CellUnion& y) {
  return x.cell_ids() != y.cell_ids();
}

bool S2CellUnion::Contains(const S2Cell& cell) const {
  return Contains(cell.id());
}

bool S2CellUnion::MayIntersect(const S2Cell& cell) const {
  return Intersects(cell.id());
}

void S2CellUnion::Encode(Encoder* const encoder) const {
  // Unsigned char for version number, and N+1 uint64_t's for N cell_ids
  // (1 for vector length, N for the ids).
  encoder->Ensure(sizeof(unsigned char) +
                  sizeof(uint64_t) * (1 + cell_ids_.size()));

  encoder->put8(kCurrentLosslessEncodingVersionNumber);
  encoder->put64(uint64_t{cell_ids_.size()});
  for (const S2CellId cell_id : cell_ids_) {
    cell_id.Encode(encoder);
  }
}

bool S2CellUnion::Decode(Decoder* const decoder) {
  // Should contain at least version and vector length.
  if (decoder->avail() < sizeof(unsigned char) + sizeof(uint64_t)) return false;
  unsigned char version = decoder->get8();
  if (version > kCurrentLosslessEncodingVersionNumber) return false;

  uint64_t num_cells = decoder->get64();
  if (num_cells > static_cast<size_t>(
                      absl::GetFlag(FLAGS_s2cell_union_decode_max_num_cells))) {
    return false;
  }

  vector<S2CellId> temp_cell_ids(num_cells);
  for (size_t i = 0; i < num_cells; ++i) {
    if (!temp_cell_ids[i].Decode(decoder)) return false;
  }
  cell_ids_.swap(temp_cell_ids);
  return true;
}

bool S2CellUnion::Contains(const S2Point& p) const {
  return Contains(S2CellId(p));
}

std::string S2CellUnion::ToString() const {
  static const size_t kMaxCount = 500;
  string output = absl::StrCat("Size:", size(), " S2CellIds:");
  for (int i = 0, limit = min(kMaxCount, size()); i < limit; ++i) {
    if (i > 0) output += ",";
    output += cell_id(i).ToToken();
  }
  if (size() > kMaxCount) output += ",...";
  return output;
}
