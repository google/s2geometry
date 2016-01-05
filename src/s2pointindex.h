// Copyright 2015 Google Inc. All Rights Reserved.
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

#ifndef S2GEOMETRY_S2POINTINDEX_H_
#define S2GEOMETRY_S2POINTINDEX_H_

#include "util/btree/btree_map.h"  // Like std::map, but faster and smaller.
#include "s2.h"
#include "s2cellid.h"

// S2PointIndex maintains an index of points sorted by leaf S2CellId.  Each
// point has some associated client-supplied data, such as an integer or
// pointer.  This can be used to map results back to client data structures.
//
// The class supports adding or removing points dynamically, and provides a
// seekable iterator interface for navigating the index.
//
// You can use this class in conjuction with S2ClosestPointQuery to find the
// closest index points to a given query point.  For example:
//
// void Test(vector<S2Point> const& points, S2Point const& target) {
//   // The template argument allows auxiliary data to be attached to each
//   // point (in this case, the array index).
//   S2PointIndex<int> index;
//   for (int i = 0; i < points.size(); ++i) {
//     index.Add(points[i], i);
//   }
//   S2ClosestPointQuery<int> query(index);
//   query.FindClosestPoint(target);
//   if (query.num_points() > 0) {
//     // query.point(0) is the closest point (result 0).
//     // query.distance(0) is the distance to the target.
//     // query.data(0) is the auxiliary data (the array index set above).
//     DoSomething(query.point(0), query.data(0), query.distance(0));
//   }
// }
//
// Alternatively, you can access the index directly using the iterator
// interface.  For example, here is how to iterate through all the points in a
// given S2CellId "target_id":
//
//   S2PointIndex<int>::Iterator it(index);
//   it.Seek(target_id.range_min());
//   for (; !it.Done() && it.id() <= target_id.range_max(); it.Next()) {
//     DoSomething(it.id(), it.point(), it.data());
//   }
//
// Points can be added or removed from the index at any time by calling Add()
// or Remove().  However when the index is modified, you must call Reset() on
// each iterator before using it again (or simply create a new iterator).
//
//   index.Add(new_point, 123456);
//   it.Reset();
//   it.Seek(target.range_min());
//
// TODO(ericv): Make this a subtype of S2Region, so that it can also be used
// to efficiently compute coverings of a collection of S2Points.
//
// REQUIRES: "Data" has default and copy constructors.
// REQUIRES: "Data" has operator== and operator<.
template <class Data>
class S2PointIndex {
 public:
  // PointData is essentially std::pair with named fields.  It stores an
  // S2Point and its associated client data.
  class PointData {
   public:
    PointData() {}  // Needed by STL
    PointData(S2Point const& point, Data const& data)
        : point_(point), data_(data) {
    }
    S2Point const& point() const { return point_; }
    Data const& data() const { return data_; }
    bool operator==(PointData const& other) const {
      return point_ == other.point_ && data_ == other.data_;
    }
    // Not required by S2PointIndex but useful for tests.
    bool operator<(PointData const& other) const {
      return point_ < other.point_ || (point_ == other.point_ &&
                                       data_ < other.data_);
    }
   private:
    S2Point point_;
    Data data_;
  };

  // Default constructor.
  S2PointIndex();

  // Returns the number of points in the index.
  int num_points() const;

  // Adds the given point to the index.  Invalidates all iterators; each
  // iterator's Reset() method must be called before it is used again.
  void Add(S2Point const& point, Data const& data);
  void Add(PointData const& point_data);

  // Removes the given point from the index.  Both the "point" and "data"
  // fields must match the point to be removed.  Returns false if the given
  // point was not present.  Invalidates all iterators; each iterator's
  // Reset() method must be called before it is used again.
  bool Remove(S2Point const& point, Data const& data);
  bool Remove(PointData const& point_data);

  // Resets the index to its original empty state.  Invalidates all iterators;
  // each iterator's Reset() method must be called before it is used again.
  void Reset();

 private:
  // Defined here because the Iterator class below uses it.
  typedef util::btree::btree_multimap<S2CellId, PointData> Map;

 public:
  class Iterator {
   public:
    // Default constructor; must be followed by a call to Init().
    Iterator();

    // Convenience constructor that calls Init().
    explicit Iterator(S2PointIndex const& index);

    // Initialize an iterator for the given S2ShapeIndex.  If the index is
    // non-empty, the iterator is positioned at the first cell.
    void Init(S2PointIndex const& index);

    // Reset the iterator to its original state (positioned at the first point
    // in the index).  This method must be called after the index is modified
    // to restore the iterator to a valid state.
    void Reset();

    // The S2CellId for the current index entry.
    // REQUIRES: !Done()
    S2CellId id() const;

    // The point associated with the current index entry.
    // REQUIRES: !Done()
    S2Point const& point() const;

    // The client-supplied data associated with the current index entry.
    // REQUIRES: !Done()
    Data const& data() const;

    // The (S2Point, data) pair associated with the current index entry.
    PointData const& point_data() const;

    // Advance the iterator to the next index entry.
    // REQUIRES: !Done()
    void Next();

    // Position the iterator at the previous index entry.
    // REQUIRES: !AtBegin()
    void Prev();

    // Return true if the iterator is positioned past the last index entry.
    bool Done() const;

    // Return true if the iterator is positioned at the first index entry.
    bool AtBegin() const;

    // Position the iterator at the first entry with id() >= target, or at the
    // end of the index if no such entry exists.
    void Seek(S2CellId target);

    // Advance the iterator to the next entry with id() >= target.  If the
    // iterator is Done() or already satisfies id() >= target, do nothing.
    void SeekForward(S2CellId target);

    // Position the iterator so that Done() is true.
    void Finish();

   private:
    Map const* map_;
    typename Map::const_iterator iter_, end_;
  };

 private:
  friend class Iterator;
  Map map_;

  S2PointIndex(S2PointIndex const&) = delete;
  void operator=(S2PointIndex const&) = delete;
};


//////////////////   Implementation details follow   ////////////////////


template <class Data>
S2PointIndex<Data>::S2PointIndex() {
}

template <class Data>
inline int S2PointIndex<Data>::num_points() const {
  return map_.size();
}

template <class Data>
void S2PointIndex<Data>::Add(PointData const& point_data) {
  S2CellId id = S2CellId::FromPoint(point_data.point());
  map_.insert(std::make_pair(id, point_data));
}

template <class Data>
void S2PointIndex<Data>::Add(S2Point const& point, Data const& data) {
  Add(PointData(point, data));
}

template <class Data>
bool S2PointIndex<Data>::Remove(PointData const& point_data) {
  S2CellId id = S2CellId::FromPoint(point_data.point());
  for (typename Map::iterator it = map_.lower_bound(id), end = map_.end();
       it != end && it->first == id; ++it) {
    if (it->second == point_data) {
      map_.erase(it);
      return true;
    }
  }
  return false;
}

template <class Data>
bool S2PointIndex<Data>::Remove(S2Point const& point, Data const& data) {
  return Remove(PointData(point, data));
}

template <class Data>
void S2PointIndex<Data>::Reset() {
  map_.clear();
}

template <class Data>
inline S2PointIndex<Data>::Iterator::Iterator() : map_(nullptr) {
}

template <class Data>
inline S2PointIndex<Data>::Iterator::Iterator(
    S2PointIndex<Data> const& index) {
  Init(index);
}

template <class Data>
inline void S2PointIndex<Data>::Iterator::Init(
    S2PointIndex<Data> const& index) {
  map_ = &index.map_;
  Reset();
}

template <class Data>
inline void S2PointIndex<Data>::Iterator::Reset() {
  iter_ = map_->begin();
  end_ = map_->end();
}

template <class Data>
inline S2CellId S2PointIndex<Data>::Iterator::id() const {
  DCHECK(!Done());
  return iter_->first;
}

template <class Data>
inline S2Point const& S2PointIndex<Data>::Iterator::point() const {
  DCHECK(!Done());
  return iter_->second.point();
}

template <class Data>
inline Data const& S2PointIndex<Data>::Iterator::data() const {
  DCHECK(!Done());
  return iter_->second.data();
}

template <class Data>
inline typename S2PointIndex<Data>::PointData const&
S2PointIndex<Data>::Iterator::point_data() const {
  DCHECK(!Done());
  return iter_->second;
}

template <class Data>
inline void S2PointIndex<Data>::Iterator::Next() {
  DCHECK(!Done());
  ++iter_;
}

template <class Data>
inline void S2PointIndex<Data>::Iterator::Prev() {
  DCHECK(!AtBegin());
  --iter_;
}

template <class Data>
inline bool S2PointIndex<Data>::Iterator::Done() const {
  return iter_ == end_;
}

template <class Data>
inline bool S2PointIndex<Data>::Iterator::AtBegin() const {
  return iter_ == map_->begin();
}

template <class Data>
inline void S2PointIndex<Data>::Iterator::Seek(S2CellId target) {
  iter_ = map_->lower_bound(target);
}

template <class Data>
inline void S2PointIndex<Data>::Iterator::SeekForward(
    S2CellId target) {
  if (!Done() && id() < target) Seek(target);
}

template <class Data>
inline void S2PointIndex<Data>::Iterator::Finish() {
  iter_ = end_;
}

#endif  // S2GEOMETRY_S2POINTINDEX_H_
