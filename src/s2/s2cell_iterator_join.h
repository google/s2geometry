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

#include <algorithm>
#include <array>
#include <bitset>
#include <functional>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include "s2/base/logging.h"
#include "absl/base/optimization.h"
#include "s2/s1chord_angle.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2cell_iterator.h"
#include "s2/s2cell_range_iterator.h"

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

  // Does a tolerant join (when the tolerance is non-zero).
  template <typename Visitor>
  bool TolerantJoin(Visitor& visitor);

  // Returns a minimal set of cells that cover the range of the iterator to seed
  // the tolerant join set.  Returns empty list if the iterator has no cells.
  template <typename Iterator>
  static absl::InlinedVector<S2CellId, 4> GetIteratorCovering(Iterator& iter);

  // Forms all possible pairs of cells from two lists.  Visits any pairs that
  // are closer to each other than the tolerance.
  template <typename Visitor>
  bool ProcessNearbyCellPairs(absl::InlinedVector<S2CellId, 4> cells_a,
                              absl::InlinedVector<S2CellId, 4> cells_b,
                              Visitor& visitor);

  // Processes a pair of cells known to be within tolerance of each other.
  //
  // We find the cells in the indexes that are covered by cell_a and cell_b.  If
  // the number of resulting (a,b) cells isn't too large, then we report them to
  // the visitor for processing. Otherwise, we push whichever cell has _fewer_
  // matches onto the cell stack along with the children of the other cell, and
  // process them recursively.
  //
  // Since there are only thirty levels to the cell hierarchy, this recursion is
  // safe as we'll never go more than 30 stack frames deep.
  template <typename Visitor>
  bool ProcessCellPair(const S2Cell& cell_a, const S2Cell& cell_b,
                       Visitor& visitor);

  // Find cells from a given iterator that are covered by a given cell.
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
    typename std::enable_if<S2CellIterator::ImplementedBy<IterA>{} &&
                            S2CellIterator::ImplementedBy<IterB>{},
                            bool>::type = true>
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
        const uint64 lsb_a = iter_a_.id().lsb();
        const uint64 lsb_b = iter_b_.id().lsb();
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
template <typename Iterator>
absl::InlinedVector<S2CellId, 4> S2CellIteratorJoin<A, B>::GetIteratorCovering(
    Iterator& iter) {
  absl::InlinedVector<S2CellId, 4> cells;

  // Position at the front, if the iterator's empty, we're done.
  iter.Begin();
  if (iter.done()) {
    return cells;
  }

  // Find the range spanned by the iterator.
  // TODO(b/248535702): Once GetCellUnionBound is migrated we can replace this.
  S2CellId min_cell = iter.range_min();
  iter.Finish();
  if (!iter.Prev()) {
    return cells;  // Empty index, we're done.
  }
  S2CellId max_cell = iter.range_max();

  // Find cell level that encompasses all of the iterator range.
  int level = min_cell.GetCommonAncestorLevel(max_cell);

  if (level < 0) {
    // The ends of the cell range don't have a common ancestor, meaning
    // they're on different faces, so iterate the faces and add the ones that
    // overlap.
    const S2CellId end_id = max_cell.parent(0);
    for (S2CellId id = min_cell.parent(0); id <= end_id; id = id.next()) {
      cells.emplace_back(id);
    }
  } else {
    cells.emplace_back(min_cell.parent(level));
  }

  return cells;
}

template <typename A, typename B>
template <typename Visitor>
bool S2CellIteratorJoin<A, B>::TolerantJoin(Visitor& visitor) {
  absl::InlinedVector<S2CellId, 4> cells_a, cells_b;

  // Seed cells with a covering for iterator A, quit if empty.
  cells_a = GetIteratorCovering(iter_a_);
  if (cells_a.empty()) {
    return true;
  }

  // Seed cells with a covering for iterator B, quit if empty.
  cells_b = GetIteratorCovering(iter_b_);
  if (cells_b.empty()) {
    return true;
  }

  return ProcessNearbyCellPairs(cells_a, cells_b, visitor);
}

template <typename A, typename B>
template <typename Visitor>
bool S2CellIteratorJoin<A, B>::ProcessNearbyCellPairs(
    absl::InlinedVector<S2CellId, 4> cells_a,
    absl::InlinedVector<S2CellId, 4> cells_b, Visitor& visitor) {
  for (const auto& id_a : cells_a) {
    const S2Cell cell_a(id_a);

    for (const auto& id_b : cells_b) {
      const S2Cell cell_b(id_b);

      if (cell_a.GetDistance(cell_b) <= tolerance_) {
        if (!ProcessCellPair(cell_a, cell_b, visitor)) {
          return false;
        }
      }
    }
  }
  return true;
}

template <typename A, typename B>
template <typename Visitor>
bool S2CellIteratorJoin<A, B>::ProcessCellPair(const S2Cell& cell_a,
                                               const S2Cell& cell_b,
                                               Visitor& visitor) {
  // Find matches for the A cell.
  int num_covered_a = EstimateCoveredCells(iter_a_, cell_a.id());
  if (num_covered_a == 0) {
    return true;
  }

  // Find matches for the B cell.
  int num_covered_b = EstimateCoveredCells(iter_b_, cell_b.id());
  if (num_covered_b == 0) {
    return true;
  }

  // If there's not too many matches, we can just process them directly
  if (num_covered_a * num_covered_b < kMaxCrossProduct) {
    // Pre-compute the S2Cells for the B side to avoid having to do it on every
    // iteration of the outer loop.
    matched_cells_.clear();
    iter_b_.Locate(cell_b.id());
    for (int i = 0; i < num_covered_b; ++i, iter_b_.Next()) {
      matched_cells_.emplace_back(iter_b_.id());
    }

    iter_a_.Locate(cell_a.id());
    for (int i = 0; i < num_covered_a; ++i, iter_a_.Next()) {
      iter_b_.Locate(cell_b.id());
      for (int j = 0; j < num_covered_b; ++j, iter_b_.Next()) {
        if (S2Cell(iter_a_.id()).GetDistance(matched_cells_[j]) <= tolerance_) {
          if (!visitor(iter_a_.iterator(), iter_b_.iterator())) {
            return false;
          }
        }
      }
    }
    return true;
  }

  // Too many matches, recurse on A or B, whichever's larger.
  absl::InlinedVector<S2CellId, 4> cells_a, cells_b;

  if (num_covered_a > num_covered_b) {
    cells_b.push_back(cell_b.id());
    for (int i = 0; i < 4; ++i) {
      cells_a.push_back(cell_a.id().child(i));
    }
  } else {
    cells_a.push_back(cell_a.id());
    for (int i = 0; i < 4; ++i) {
      cells_b.push_back(cell_b.id().child(i));
    }
  }
  return ProcessNearbyCellPairs(cells_a, cells_b, visitor);
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
      for (; iter.id() <= end; iter.Next()) {
        ++matches;

        // Give up after max/2 matches.  There's too many matches anyways, so
        // this is a heuristic that will help us decide whether to recurse on
        // the left or right.
        if (matches >= kMaxCrossProduct / 2) {
          return matches;
        }
      }
      return matches;
    }
  }
  ABSL_UNREACHABLE();
}

#endif  // S2_S2CELL_ITERATOR_JOIN_H_
