// Copyright 2022 Google Inc. All Rights Reserved.
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

#ifndef S2_S2CELL_ITERATOR_JOIN_H_
#define S2_S2CELL_ITERATOR_JOIN_H_

#include <cstdint>
#include <functional>
#include <type_traits>
#include <vector>

#include "absl/base/optimization.h"
#include "absl/container/inlined_vector.h"
#include "absl/types/span.h"
#include "s2/internal/s2meta.h"
#include "s2/s1chord_angle.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2cell_iterator.h"
#include "s2/s2cell_range_iterator.h"
#include "s2/s2cell_union.h"

// Defines a class which can be used to perform an inner join operation on any
// two S2CellIterator iterators.  S2CellIteratorJoin takes an optional distance
// value which specifies a "buffer" around cells in the iterator.  If cells are
// within that distance, then they're considered to overlap.  This allows us to
// support "tolerant" versions of queries.  Currently the tolerance must not be
// negative.
//
// The iterator will find each pair of cells that are less than or equal to the
// tolerance distance from each other and call a visitor with a reference to
// each of the properly positioned iterators.  The visitor can then process the
// overlapping iterators however it wishes, returning true if it wishes to
// continue iterating, and false otherwise.
//
// Example usage:
//
//   // Process cell pairs until we find two cells with different edges.
//   bool ProcessRow(
//       const MutableS2ShapeIndex::Iterator& itera,
//       const MutableS2ShapeIndex::Iterator& iterb) {
//     return AllEdgesAreEqual(itera.cell(), iterb.cell());
//   }
//
//   MutableS2ShapeIndex index_a, index_b;
//   MakeS2CellIteratorJoin(&index_a, &index_b, distance).Join(ProcessRow);
//
// But we're not limited to joining iterators of the same type.  Any iterator
// implementing the S2CellIterator API will work.  For example, we can join an
// S2ShapeIndex and S2PointIndex:
//
//  bool ProcessCellPoints(
//     const MutableS2ShapeIndex::Iterator& itera,
//     const S2PointIndex<std::string>::Iterator& iterb) {
//       // Process shapes and points from overlapping cells.
//       return true;
//     }
//
//  MutableS2ShapeIndex index_a;
//  S2PointIndex<std::string> index_b;
//  MakeS2CellIteratorJoin(&index_a, &index_b).Join(ProcessCellPoints);
//
// We use a visitor pattern because we can inline the call to the visitor using
// a template parameter, and, as opposed to something like the Next/Value/Done
// pattern that S2CellIterator itself uses, we can pass the visitor to a
// different class or function to actually implement the join, invisibly to the
// user.
//
// In addition to allowing us to delegate iteration nicely, the visitor pattern
// lets us separate concerns when writing spatial algorithms.  The actual
// processing of indexed data is separated from the logic of positioning the
// iterators, which allows for much cleaner implementations.
//
// The visitor is a template parameter, which must ultimately be callable with
// two const references to each input iterator type, returning true if the
// iteration should continue, and false otherwise.
//
// Since the type of the visitor is inferred and we don't have access to C++
// concepts, we can't constrain that it's callable in the way we want directly.
// Instead, we check that whatever type is given would be convertible to a
// std::function with the appropriate signature:
//
//     std::function<bool(const IteratorA&, const IteratorB&)>
//
// Where IteratorA and IteratorB are the concrete types of our two iterators.
//
template <typename IterA, typename IterB>
class S2CellIteratorJoin {
 public:
  // Expose iterator types.
  using IteratorA = IterA;
  using IteratorB = IterB;

  S2CellIteratorJoin() = default;
  S2CellIteratorJoin(const IteratorA& iter_a, const IteratorB& iter_b,
                     S1ChordAngle tolerance = {})
      : iter_a_(MakeS2CellRangeIterator(iter_a)),
        iter_b_(MakeS2CellRangeIterator(iter_b)),
        tolerance_(tolerance) {}

  // Executes the join.  Explicitly supports type inference for the visitor.
  //
  // Returns false if the visitor ever does, true otherwise.
  template <typename Visitor>
  bool Join(Visitor visitor) {
    // We want to take the visitor as a template parameter for inlining, but
    // also be type safe, so check that we can call it the way we want and print
    // a reasonable error message if we can't.
    static_assert(
        std::is_convertible<
            Visitor, std::function<bool(const IteratorA&, const IteratorB&)>>{},
        "Visitor must return bool and be callable with two const "
        "references to the iterators");

    if (tolerance_ == S1ChordAngle::Zero()) {
      return ExactJoin(visitor);
    } else {
      return TolerantJoin(visitor);
    }
  }

 private:
  S2CellRangeIterator<IteratorA> iter_a_;
  S2CellRangeIterator<IteratorB> iter_b_;
  S1ChordAngle tolerance_;

  // Reusable storage for S2Cells on one side of the join.
  std::vector<S2Cell> matched_cells_;

  // Performs an exact inner join (when the tolerance is zero).
  template <typename Visitor>
  bool ExactJoin(Visitor& visitor);

  // ---- Tolerant join related code.

  // Maximum number of cross-terms before we recurse.
  static constexpr int kMaxCrossProduct = 25;
  static constexpr int kCoverLimit = kMaxCrossProduct / 2;

  // Does a tolerant join (when the tolerance is non-zero).
  template <typename Visitor>
  bool TolerantJoin(Visitor& visitor);

  // Positions an iterator and visits every position that intersects a given
  // S2CellId.  The iterator is passed to the visitor at each position.  Returns
  // false if the visitor ever does, otherwise true.
  template <typename Iterator, typename Visitor>
  static bool ScanCellRange(Iterator& iter, S2CellId id, Visitor&& visitor);

  // Given two lists of cells, filters the cells from B that are within distance
  // of each cell of A and passes them to ProcessCellPairs().  Returns false if
  // the visitor ever does, true otherwise.
  template <typename Visitor>
  bool ProcessNearby(absl::Span<const S2Cell> cells_a,
                     absl::Span<const S2Cell> cells_b, Visitor& visitor);

  // Processes a pair of cells known to be within tolerance of each other.
  //
  // If the portion of each index that's covered by the cells is small enough,
  // then we report pairs to the visitor, otherwise the cells are subdivided and
  // we recurse.
  //
  // Since there are only thirty levels to the cell hierarchy, this recursion is
  // safe as we'll never go more than 30 stack frames deep.
  template <typename Visitor>
  bool ProcessCellPairs(  //
      const S2Cell& cell_a, absl::Span<const S2Cell> cells_b, Visitor& visitor);

  // Estimates how many values are covered by a given S2CellId.
  template <typename Iterator>
  int EstimateCoveredCells(Iterator& iter, S2CellId cell);
};

// Clang format just can't handle templated methods like this well.
// clang-format off

// Factory function to build an S2CellIteratorJoin from two iterators that
// implement the S2CellIterator API.  Type inference is explicitly supported for
// the iterator types.
//
// The iterators are copied and a new S2CellIteratorJoin instance returned.
// Note that the underlying container the iterators were created from must exist
// for the lifetime of the join operation.
template <
    typename IterA, typename IterB,
    typename std::enable_if<
        s2meta::derived_from_v<IterA, S2CellIterator> &&
        s2meta::derived_from_v<IterB, S2CellIterator>, bool>::type = true>
S2CellIteratorJoin<IterA, IterB> MakeS2CellIteratorJoin(
    const IterA& iter_a, const IterB& iter_b, S1ChordAngle tolerance = {}) {
  return {iter_a, iter_b, tolerance};
}

// Factory function to build an S2CellIteratorJoin from any two types that have
// a nested ::Iterator class.  Two new iterators are instantiated and used to
// create a new S2CellIteratorJoin instance.  The underlying containers must
// live for the lifetime of the join operation, so we take them by const pointer
// to avoid binding temporaries.
template <typename IterableA, typename IterableB>
S2CellIteratorJoin<typename IterableA::Iterator, typename IterableB::Iterator>
MakeS2CellIteratorJoin(
  const IterableA* a, const IterableB* b, S1ChordAngle tolerance = {}) {
  return MakeS2CellIteratorJoin(typename IterableA::Iterator(a),
                                typename IterableB::Iterator(b), tolerance);
}
// clang-format on

//////////////////   Implementation details follow   ////////////////////

template <typename A, typename B>
template <typename Visitor>
bool S2CellIteratorJoin<A, B>::ExactJoin(Visitor& visitor) {
  iter_a_.Begin();
  iter_b_.Begin();

  // Iterate until we hit the end of an iterator or visitor tells us to stop.
  while (!iter_a_.done() && !iter_b_.done()) {
    int order = iter_a_.Relation(iter_b_);
    switch (order) {
      case -1:
        // A precedes B, seek A.
        iter_a_.SeekTo(iter_b_);
        break;

      case +1:
        // B precedes A, seek B.
        iter_b_.SeekTo(iter_a_);
        break;

      case 0: {
        // Iterators overlap.
        // If visitor rejects pair, then we're done.
        if (!visitor(iter_a_.iterator(), iter_b_.iterator())) {
          return false;
        }

        // Move the smaller of the cells forward.
        const uint64_t lsb_a = iter_a_.id().lsb();
        const uint64_t lsb_b = iter_b_.id().lsb();
        if (lsb_a < lsb_b) {
          iter_a_.Next();
        } else if (lsb_a > lsb_b) {
          iter_b_.Next();
        } else {
          // Cells are the same size and overlap, they must be the same cell.
          // Move both forward.
          iter_a_.Next();
          iter_b_.Next();
        }
        break;
      }
    }
  }
  return true;
}

template <typename A, typename B>
template <typename Iterator, typename Visitor>
inline bool S2CellIteratorJoin<A, B>::ScanCellRange(Iterator& iter, S2CellId id,
                                                    Visitor&& visitor) {
  iter.Locate(id);
  for (; !iter.done() && iter.id().intersects(id); iter.Next()) {
    if (!visitor(iter)) {
      return false;
    }
  }
  return true;
}

template <typename A, typename B>
template <typename Visitor>
bool S2CellIteratorJoin<A, B>::TolerantJoin(Visitor& visitor) {
  const auto Covering = [&](auto& iter) -> std::vector<S2Cell> {
    iter.Begin();
    if (iter.done()) {
      return {};
    }

    S2CellId min = iter.range_min();
    iter.Finish();
    if (!iter.Prev()) {
      return {};  // Empty iterator.
    }
    S2CellId max = iter.range_max();

    std::vector<S2Cell> cells;
    for (S2CellId id : S2CellUnion::FromMinMax(min, max)) {
      cells.emplace_back(id);
    }
    return cells;
  };

  // Seed the recursion with a coarse covering of each input iterator's current
  // position.
  return ProcessNearby(Covering(iter_a_), Covering(iter_b_), visitor);
}

template <typename A, typename B>
template <typename Visitor>
bool S2CellIteratorJoin<A, B>::ProcessNearby(absl::Span<const S2Cell> cells_a,
                                             absl::Span<const S2Cell> cells_b,
                                             Visitor& visitor) {
  std::vector<S2Cell> nearby_cells;
  for (const S2Cell& cell_a : cells_a) {
    nearby_cells.clear();
    for (const S2Cell& cell_b : cells_b) {
      if (cell_a.GetDistance(cell_b) <= tolerance_) {
        nearby_cells.emplace_back(cell_b);
      }
    }

    if (!nearby_cells.empty()) {
      if (!ProcessCellPairs(cell_a, nearby_cells, visitor)) {
        return false;
      }
    }
  }
  return true;
}

template <typename A, typename B>
template <typename Visitor>
bool S2CellIteratorJoin<A, B>::ProcessCellPairs(
    const S2Cell& cell_a, absl::Span<const S2Cell> cells_b, Visitor& visitor) {
  // Estimate how many index cells the A cell covers.
  int num_covered_a = EstimateCoveredCells(iter_a_, cell_a.id());
  if (num_covered_a == 0) {
    return true;
  }

  // If the A cell covers too many index cells, then subdivide it.
  absl::InlinedVector<S2Cell, 1> cells_a = {cell_a};
  if (num_covered_a >= kCoverLimit) {
    cells_a.resize(4);
    cell_a.Subdivide(&cells_a[0]);
  }

  // Scan the cells of the B union.  Prune any cells that don't cover any of the
  // index, and subdivide any that cover too much.
  bool subdivided = false;
  absl::InlinedVector<S2Cell, 1> subdivided_b;
  for (const S2Cell& cell_b : cells_b) {
    int num_covered_b = EstimateCoveredCells(iter_b_, cell_b.id());
    if (num_covered_b == 0) {
      continue;
    } else if (num_covered_b < kCoverLimit) {
      subdivided_b.emplace_back(cell_b);
    } else {
      S2Cell children[4];
      cell_b.Subdivide(children);
      subdivided_b.insert(subdivided_b.end(), children, children + 4);
      subdivided = true;
    }
  }

  // If A covers too many index cells or we had to subdivide B, then continue
  // the recursion by pairing the A cells B cells that are within the tolerance.
  if (num_covered_a >= kCoverLimit || subdivided) {
    return ProcessNearby(cells_a, subdivided_b, visitor);
  }

  // Otherwise A and the B union are small enough we can pair them up and report
  // nearby pairs to the visitor now.

  // Pre-compute the B S2Cells to avoid replicating work in the inner loop.
  matched_cells_.clear();
  for (const S2Cell& cell_b : cells_b) {
    ScanCellRange(iter_b_, cell_b.id(), [&](const auto& iter) {
      matched_cells_.emplace_back(iter.id());
      return true;
    });
  }

  return ScanCellRange(iter_a_, cell_a.id(), [&](const auto& iter_a) {
    // If the index cell doesn't have it's endpoint in this cell then ignore it.
    // This makes it so that we only see each index cell in a single A cell.
    if (!cell_a.id().intersects(iter_a.id().range_min())) {
      return true;
    }

    const S2Cell sub_cell_a(iter_a_.id());

    // For each sub cell in the cell_a range, scan each B cell in the cell union
    // and send the resulting (A,B) pairs to the visitor.  If it returns false
    // at any point then stop the join and return false.
    int idx = 0;
    bool success = true;
    for (int i = 0; i < cells_b.size() && success; ++i) {
      S2CellId id = cells_b[i].id();
      success &= ScanCellRange(iter_b_, id, [&](const auto& iter_b) {
        if (sub_cell_a.GetDistance(matched_cells_[idx++]) <= tolerance_) {
          if (!visitor(iter_a.iterator(), iter_b.iterator())) {
            return false;
          }
        }
        return true;
      });
    }
    return success;
  });
}

template <typename A, typename B>
template <typename Iterator>
int S2CellIteratorJoin<A, B>::EstimateCoveredCells(Iterator& iter,
                                                   S2CellId cell) {
  switch (iter.Locate(cell)) {
    case S2CellRelation::DISJOINT:
      return 0;

    case S2CellRelation::INDEXED:
      return 1;

    case S2CellRelation::SUBDIVIDED: {
      // Note: If we develop a Distance method for Iterators that can tell us
      // how many cells within a range efficiently, then we can avoid the linear
      // scan here that just has us repeatedly aborting.
      const S2CellId end = cell.range_max();
      int matches = 0;
      for (; !iter.done() && iter.id() <= end; iter.Next()) {
        ++matches;

        // There's too many matches, so give up. This is a heuristic that will
        // help us decide whether to recurse on the left or right.
        if (matches > kCoverLimit) {
          return kCoverLimit;
        }
      }
      return matches;
    }
  }
  ABSL_UNREACHABLE();
}

#endif  // S2_S2CELL_ITERATOR_JOIN_H_
