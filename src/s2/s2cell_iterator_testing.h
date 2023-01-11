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

#ifndef S2_S2CELL_ITERATOR_TESTING_H_
#define S2_S2CELL_ITERATOR_TESTING_H_

#include "absl/container/btree_map.h"
#include "s2/s2cell_iterator.h"

// A mock iterator for testing.  Iterates an absl::btree_map mapping S2CellId to
// another type.  Rather than instantiating directly, consider using
// MakeMockS2CellIterator for convenience.

template <typename T>
class MockS2CellIterator final : public S2CellIterator {
 public:
  // We use a btree_map because we require that S2CellIterator containers be
  // ordered for fast seek and locate.
  using CellMap = absl::btree_map<S2CellId, T>;

  // Build an iterator from a btree map.  The map must persist for the lifetime
  // of the iterator, so we take it by const pointer to avoid binding
  // temporaries.
  explicit MockS2CellIterator(const CellMap* map)
      : map_(map), begin_(map->begin()), end_(map->end()) {
    Begin();
  }

  // Returns the current value the iterator is positioned at.
  //
  // REQUIRES: !done()
  const T& value() const {
    S2_DCHECK(!done());
    return iter_->second;
  }

  // Returns the current S2CellId at which the iterator is positioned.
  S2CellId id() const override {
    if (done()) {
      return S2CellId::Sentinel();
    }
    return iter_->first;
  }

  // Returns true if the iterator has reached the end of the input.
  bool done() const override { return iter_ == end_; }

  // Positions the iterator at the first position.
  void Begin() override { iter_ = begin_; }

  // Positions the iterator past the last value.
  void Finish() override { iter_ = end_; }

  // Positions the iterator at the next value.
  //
  // REQUIRES: !done()
  void Next() override {
    S2_DCHECK(!done());
    ++iter_;
  }

  // Positions the iterator at the previous value, returning false if it's
  // already positioned at the beginning, true otherwise.
  bool Prev() override {
    if (iter_ == begin_) {
      return false;
    }
    --iter_;
    return true;
  }

  // Seeks the iterator to the first cell with id() >= target or the end
  // of the iterator if no such cell exists.
  void Seek(S2CellId target) override { iter_ = map_->lower_bound(target); }

  // Positions the iterator at the cell containing target and returns true. If
  // no such cell exists, then returns false and leaves the iterator in an
  // undefined (but valid) state.
  bool Locate(const S2Point& target) override {
    return LocateImpl(*this, target);
  }

  // Let T be the target S2CellId.  If T is contained by some index cell I
  // (including equality), this method positions the iterator at I and returns
  // INDEXED.  Otherwise if T contains one or more (smaller) index cells, it
  // positions the iterator at the first such cell I and returns SUBDIVIDED.
  // Otherwise it returns DISJOINT and leaves the iterator in an undefined
  // (but valid) state.
  S2CellRelation Locate(S2CellId target) override {
    return LocateImpl(*this, target);
  }

 private:
  const CellMap* map_;
  typename CellMap::const_iterator iter_, begin_, end_;
};

// Construct a MockS2CellIterator from a btree_map of S2CellIDs to any type.
//
// This function explicitly supports inferring the template parameter T.
template <typename T>
MockS2CellIterator<T> MakeMockS2CellIterator(
    const absl::btree_map<S2CellId, T>* map) {
  return MockS2CellIterator<T>(map);
}

#endif  // S2_S2CELL_ITERATOR_TESTING_H_
