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

#ifndef S2_S2CELL_RANGE_ITERATOR_H_
#define S2_S2CELL_RANGE_ITERATOR_H_

#include "absl/meta/type_traits.h"
#include "s2/s2cell_id.h"
#include "s2/s2cell_iterator.h"
#include "s2/s2shape_index.h"

// S2CellRangeIterator is a wrapper around an S2CellIterator that caches the
// range_min() and range_max() that each cell covers as it iterates.  This lets
// us define range based methods such as SeekTo() and Locate() efficiently.
//
// Computing the range_max() and range_min() of a cell isn't expensive but it's
// not free either, so we extend the S2CellIterator interface instead of
// integrating this functionality there, allowing the user to pay only for what
// they use.
//
// An S2CellRangeIterator wraps an S2CellIterator, but is also itself an
// S2CellIterator and thus can be used anywhere one is required.
template <typename Iterator>
class S2CellRangeIterator final : public S2CellIterator {
 public:
  static_assert(S2CellIterator::ImplementedBy<Iterator>{},
                "S2CellRangeIterator requires an S2CellIterator");

  // Construct a new S2CellRangeIterator positioned at the beginning.
  explicit S2CellRangeIterator(Iterator iter) : it_(std::move(iter)) {
    Begin();
  }

  // Returns a const reference to the underlying iterator.
  const Iterator& iterator() const { return it_; }

  // The current S2CellId and cell contents.
  S2CellId id() const override { return it_.id(); }

  // The min and max leaf cell ids covered by the current cell.  If done() is
  // true, these methods return a value larger than any valid cell id.
  S2CellId range_min() const { return range_min_; }
  S2CellId range_max() const { return range_max_; }

  // Queries the relationship between two range iterators.  Returns -1 if this
  // iterator's current position entirely precedes the other iterator's current
  // position, +1 if it entirely follows, and 0 if they overlap.
  template <typename T>
  int Relation(const S2CellRangeIterator<T>& b) {
    if (range_max() < b.range_min()) return -1;
    if (range_min() > b.range_max()) return +1;
    return 0;
  }

  // S2CellIterator API.
  void Begin() override;
  void Next() override;
  bool Prev() override;
  void Seek(S2CellId target) override;
  void Finish() override;
  bool done() const override { return it_.done(); }
  bool Locate(const S2Point& target) override;
  S2CellRelation Locate(S2CellId target) override;

  // The same as above, but uses another S2CellRangeIterator as the target.
  template <typename T>
  S2CellRelation Locate(const S2CellRangeIterator<T>& target);

  // Position the iterator at the first cell that overlaps or follows
  // "target", i.e. such that range_max() >= target.range_min().
  template <typename T>
  void SeekTo(const S2CellRangeIterator<T>& target);

  // Position the iterator at the first cell that follows "target", i.e. the
  // first cell such that range_min() > target.range_max().
  template <typename T>
  void SeekBeyond(const S2CellRangeIterator<T>& target);

 private:
  // Updates internal state after the iterator has been repositioned.
  void Refresh();

  Iterator it_;
  S2CellId range_min_, range_max_;
};

// Builds a new S2CellRangeIterator from an index, supporting type inference.
//
// We may wish to provide overloads for other types in the future, so we
// disqualify this function from overload resolution using std::enable_if when
// the type isn't an S2ShapeIndex.
//
// The index must live for the duration of the iterator, so we take it by const
// pointer instead of reference to avoid binding to temporaries.
template <typename IndexType,
          typename std::enable_if<S2ShapeIndex::ImplementedBy<IndexType>{},
                                  bool>::type = true>
S2CellRangeIterator<typename IndexType::Iterator> MakeS2CellRangeIterator(
    const IndexType* index) {
  using Iterator = typename IndexType::Iterator;
  return S2CellRangeIterator<Iterator>(Iterator(index, S2ShapeIndex::BEGIN));
}

// Builds a new S2CellRangeIterator from an S2CellIterator, explicitly supports
// type inference for the Iterator parameter.
template <typename Iterator,
          typename std::enable_if<S2CellIterator::ImplementedBy<Iterator>{},
                                  bool>::type = true>
auto MakeS2CellRangeIterator(Iterator&& iter) {
  return S2CellRangeIterator<absl::decay_t<Iterator>>(
      std::forward<Iterator>(iter));
}

//////////////////   Implementation details follow   ////////////////////

template <typename Iterator>
void S2CellRangeIterator<Iterator>::S2CellRangeIterator::Begin() {
  it_.Begin();
  Refresh();
}

template <typename Iterator>
void S2CellRangeIterator<Iterator>::S2CellRangeIterator::Next() {
  it_.Next();
  Refresh();
}

template <typename Iterator>
bool S2CellRangeIterator<Iterator>::S2CellRangeIterator::Prev() {
  bool status = it_.Prev();
  Refresh();
  return status;
}

template <typename Iterator>
void S2CellRangeIterator<Iterator>::S2CellRangeIterator::Seek(S2CellId target) {
  it_.Seek(target);
  Refresh();
}

template <typename Iterator>
void S2CellRangeIterator<Iterator>::S2CellRangeIterator::Finish() {
  it_.Finish();
  Refresh();
}

template <typename Iterator>
bool S2CellRangeIterator<Iterator>::S2CellRangeIterator::Locate(
    const S2Point& target) {
  bool status = it_.Locate(target);
  Refresh();
  return status;
}

template <typename Iterator>
S2CellRelation S2CellRangeIterator<Iterator>::Locate(S2CellId target) {
  // Let T be the target cell id, let I = Seek(T.range_min()) and let Prev(I) be
  // the predecessor of I.  If T contains any index cells, then T contains I.
  // Similarly, if T is contained by an index cell, then the containing cell is
  // either I or Prev(I).  We test for containment by comparing the ranges of
  // leaf cells spanned by T, I, and Prev(I).
  Seek(target.range_min());
  if (!done()) {
    // The target is contained by the cell we landed on, so it's indexed.
    if (id() >= target && range_min() <= target) {
      return S2CellRelation::INDEXED;
    }

    // The cell we landed on is contained by the target, so it's subdivided.
    if (id() <= target.range_max()) {
      return S2CellRelation::SUBDIVIDED;
    }
  }

  // Otherwise check the previous cell (if it exists).  If it contains the
  // target then it's indexed, otherwise the target cell is disjoint.
  if (Prev() && range_max() >= target) {
    return S2CellRelation::INDEXED;
  }
  return S2CellRelation::DISJOINT;
}

// Convenience re-implementation of the above function, see it for details.
template <typename Iterator>
template <typename T>
S2CellRelation S2CellRangeIterator<Iterator>::S2CellRangeIterator::Locate(
    const S2CellRangeIterator<T>& target) {
  Seek(target.range_min());
  if (!done()) {
    // The target is contained by the cell we landed on, so it's indexed.
    if (id() >= target.id() && range_min() <= target.id()) {
      return S2CellRelation::INDEXED;
    }

    // The cell we landed on is contained by the target, so it's subdivided.
    if (id() <= target.range_max()) {
      return S2CellRelation::SUBDIVIDED;
    }
  }

  // Otherwise check the previous cell (if it exists).  If it contains the
  // target then it's indexed, otherwise the target cell is disjoint.
  if (Prev() && range_max() >= target.id()) {
    return S2CellRelation::INDEXED;
  }
  return S2CellRelation::DISJOINT;
}

template <typename Iterator>
template <typename T>
void S2CellRangeIterator<Iterator>::SeekTo(
    const S2CellRangeIterator<T>& target) {
  Seek(target.range_min());

  // If the current cell does not overlap "target", it is possible that the
  // previous cell is the one we are looking for.  This can only happen when
  // the previous cell contains "target" but has a smaller S2CellId.
  if (done() || range_min() > target.range_max()) {
    if (Prev() && range_max() < target.id()) {
      Next();
    }
  }
  Refresh();
}

template <typename Iterator>
template <typename T>
void S2CellRangeIterator<Iterator>::SeekBeyond(
    const S2CellRangeIterator<T>& target) {
  Seek(target.range_max().next());
  if (!done() && range_min() <= target.range_max()) {
    Next();
  }
  Refresh();
}

// This method is inline, but is only called by non-inline methods defined in
// this file.  Putting the definition here enforces this requirement.
template <typename Iterator>
inline void S2CellRangeIterator<Iterator>::Refresh() {
  if (done()) {
    range_min_ = S2CellId::Sentinel().range_min();
    range_max_ = S2CellId::Sentinel().range_max();
  } else {
    range_min_ = id().range_min();
    range_max_ = id().range_max();
  }
}

#endif  // S2_S2CELL_RANGE_ITERATOR_H_
