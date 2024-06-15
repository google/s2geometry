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

#ifndef S2_S2CELL_ITERATOR_H_
#define S2_S2CELL_ITERATOR_H_

#include <ostream>

#include "absl/meta/type_traits.h"
#include "s2/s2cell_id.h"

// Possible relationships between two S2CellIds in an index.
enum class S2CellRelation : int {
  INDEXED,     // Target is contained by an index cell.
  SUBDIVIDED,  // Target is subdivided into one or more index cells.
  DISJOINT     // Target does not intersect any index cells.
};

inline std::ostream& operator<<(std::ostream& os, S2CellRelation r) {
  switch (r) {
    case S2CellRelation::INDEXED:
      os << "INDEXED";
      break;

    case S2CellRelation::SUBDIVIDED:
      os << "SUBDIVIDED";
      break;

    case S2CellRelation::DISJOINT:
      os << "DISJOINT";
      break;
  }
  return os;
}

// An abstract base class for iterators over any sorted collection keyed by
// S2CellId.
//
// This is intentionally opaque to the type of the value we might be mapping to,
// which can be anything (index cells, points, integers, other S2CellIds, etc),
// and instead only defines the common subset of functionality for positioning
// the iterator and querying the current cell id.
//
// This class is generally not used directly via pointer but instead to type
// check that a given template parameter implements the interface.
//
// A default implementation for Locate is given as a static method.  Inheritors
// should call it directly when defining their overloads of the API and mark
// their methods final to ensure de-virtualization when using sub-classes
// directly.
//
// A canonical implementation for Locate thus might look like:
//
//   class MyIterator : public S2CellIterator {
//    public:
//     void Locate(const S2Point& point) final {
//       return LocateImpl(*this, point);
//     }
//   };
//
// This will ensure code that uses MyIterator directly and calls Locate will
// directly use the methods of MyIterator instead of calling through the vtable.
class S2CellIterator {
 public:
  // A type function to check if a type is derived from S2CellIterator.  This is
  // useful for writing static checks on template parameters when we want to
  // inline a particular iterator call, but we need to make sure it implements
  // the interface that we want.  We don't have access to c++ concepts, so this
  // is the next best thing:
  //
  //   template <typename Iterator>
  //   void Frobnicate(Iterator& iter) {
  //     static_assert(S2CellIterator::ImplementedBy<Iterator>{},
  //       "We require an object implementing the S2CellIterator API.");
  //   }
  template <typename T>
  using ImplementedBy = std::is_convertible<absl::decay_t<T>*, S2CellIterator*>;

  S2CellIterator() = default;
  virtual ~S2CellIterator() = default;

  // Returns the current S2CellId that the iterator is positioned at.  This
  // function should be cheap to call (ideally directly returning the current
  // value).  When the iterator is done, this should return S2CellId::Sentinel.
  virtual S2CellId id() const = 0;

  // Return true if the iterator has reached the end of the input. This function
  // should be cheap to call to check if iteration has ended.
  virtual bool done() const = 0;

  // Positions the iterator at the first position.
  virtual void Begin() = 0;

  // Positions the iterator past the last value.  After calling this function,
  // the done() method should return true.
  virtual void Finish() = 0;

  // Positions the iterator at the next value.  Must not be called when done()
  // is true.
  virtual void Next() = 0;

  // Positions the iterator at the previous value.  Returns false if the
  // iterator is already at the start.
  virtual bool Prev() = 0;

  // Seeks the iterator to the first cell with id() >= target or the end
  // of the iterator if no such cell exists.
  virtual void Seek(S2CellId target) = 0;

  // Positions the iterator at the cell containing target and returns true. If
  // no such cell exists, return false and leave the iterator in an undefined
  // (but valid) state.
  virtual bool Locate(const S2Point& target) = 0;

  // Let T be the target S2CellId.  If T is contained by some index cell I
  // (including equality), this method positions the iterator at I and returns
  // INDEXED.  Otherwise if T contains one or more (smaller) index cells, it
  // positions the iterator at the first such cell I and returns SUBDIVIDED.
  // Otherwise it returns DISJOINT and leaves the iterator in an undefined
  // (but valid) state.
  virtual S2CellRelation Locate(S2CellId target) = 0;

 protected:
  template <typename Iterator>
  static inline bool LocateImpl(Iterator& iter, const S2Point& point);

  template <typename Iterator>
  static inline S2CellRelation LocateImpl(Iterator& iter, S2CellId target);

  // Disable public copying and assigning via abstract base class pointer.
  S2CellIterator(const S2CellIterator&) = default;
  S2CellIterator& operator=(const S2CellIterator&) = default;
};

//////////////////   Implementation details follow   ////////////////////

template <typename Iterator>
inline bool S2CellIterator::LocateImpl(Iterator& iter, const S2Point& point) {
  static_assert(S2CellIterator::ImplementedBy<Iterator>{},
                "Iterator must implement the S2CellIterator API.");

  // Let I = Seek(T), where T is the leaf cell containing the target point, and
  // let Prev(I) be the predecessor of I.  If T is contained by an index cell,
  // then the containing cell is either I or Prev(I).  We test for containment
  // by comparing the ranges of leaf cells spanned by T, I, and Prev(I).
  S2CellId target(point);

  iter.Seek(target);
  if (!iter.done() && iter.id().range_min() <= target) {
    return true;
  }

  if (iter.Prev() && iter.id().range_max() >= target) {
    return true;
  }
  return false;
}

template <typename Iterator>
inline S2CellRelation S2CellIterator::LocateImpl(Iterator& iter,
                                                 S2CellId target) {
  static_assert(S2CellIterator::ImplementedBy<Iterator>{},
                "Iterator must implement the S2CellIterator API.");

  // Let T be the target cell id, let I = Seek(T.range_min()) and let Prev(I) be
  // the predecessor of I.  If T contains any index cells, then T contains I.
  // Similarly, if T is contained by an index cell, then the containing cell is
  // either I or Prev(I).  We test for containment by comparing the ranges of
  // leaf cells spanned by T, I, and Prev(I).
  iter.Seek(target.range_min());
  if (!iter.done()) {
    // The target is contained by the cell we landed on, so it's indexed.
    if (iter.id() >= target && iter.id().range_min() <= target) {
      return S2CellRelation::INDEXED;
    }

    // The cell we landed on is contained by the target, so it's subdivided.
    if (iter.id() <= target.range_max()) {
      return S2CellRelation::SUBDIVIDED;
    }
  }

  // Otherwise check the previous cell (if it exists).  If it contains the
  // target then it's indexed, otherwise the target cell is disjoint.
  if (iter.Prev() && iter.id().range_max() >= target) {
    return S2CellRelation::INDEXED;
  }
  return S2CellRelation::DISJOINT;
}

#endif  // S2_S2CELL_ITERATOR_H_
