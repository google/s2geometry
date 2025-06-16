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
//
// S2ShapeIndex is an abstract base class for indexing polygonal geometry in
// memory.  The main documentation is with the class definition below.
// (Some helper classes are defined first.)
//
// S2ShapeIndex has two major subtypes:
//
//  - MutableS2ShapeIndex is for building new S2ShapeIndexes.  Indexes may
//    also be updated dynamically by inserting or deleting shapes.  Once an
//    index has been built it can be encoded compactly and later decoded as
//    either a MutableS2ShapeIndex or an EncodedS2ShapeIndex.
//
//  - EncodedS2ShapeIndex is an S2ShapeIndex implementation that works
//    directly with encoded data.  Rather than decoding everything in advance,
//    geometry is decoded incrementally (down to individual edges) as needed.
//    It can be initialized from a single block of data in nearly constant
//    time.  This saves large amounts of memory and is also much faster in the
//    common situation where geometric data is loaded from somewhere, decoded,
//    and then only a single operation is performed on it.  (Speedups for
//    large geometric objects can be over 1000x.)  It supports all
//    S2ShapeIndex operations including boolean operations, measuring
//    distances, etc.

#ifndef S2_S2SHAPE_INDEX_H_
#define S2_S2SHAPE_INDEX_H_

#include <array>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <type_traits>
#include <utility>

#include "absl/log/absl_check.h"
#include "absl/types/span.h"

#include "s2/_fp_contract_off.h"  // IWYU pragma: keep
#include "s2/s2cell_id.h"
#include "s2/s2cell_iterator.h"
#include "s2/s2point.h"
#include "s2/s2shape.h"
#include "s2/util/coding/coder.h"
#include "s2/util/gtl/compact_array.h"

class R1Interval;
class S2PaddedCell;

// S2ClippedShape represents the part of a shape that intersects an S2Cell.
// It consists of the set of edge ids that intersect that cell, and a boolean
// indicating whether the center of the cell is inside the shape (for shapes
// that have an interior).
//
// Note that the edges themselves are not clipped; we always use the original
// edges for intersection tests so that the results will be the same as the
// original shape.
class S2ClippedShape {
 public:
  S2ClippedShape() : num_edges_(0) {}

  // The shape id of the clipped shape.
  int shape_id() const;

  // Returns true if the center of the S2CellId is inside the shape.  Returns
  // false for shapes that do not have an interior.
  bool contains_center() const;

  // The number of edges that intersect the S2CellId.
  int num_edges() const;

  // Returns the edge id of the given edge in this clipped shape.  Edges are
  // sorted in increasing order of edge id.
  //
  // REQUIRES: 0 <= i < num_edges()
  int edge(int i) const;

  // Returns true if the clipped shape contains the given edge id.
  bool ContainsEdge(int id) const;

 private:
  // This class may be copied by value, but note that it does *not* own its
  // underlying data.  (It is owned by the containing S2ShapeIndexCell.)

  friend class MutableS2ShapeIndex;
  friend class S2ShapeIndexCell;
  friend class S2Stats;

  // Internal methods are documented with their definition.
  void Init(int32_t shape_id, uint32_t num_edges);
  void Destruct();
  bool is_inline() const;
  void set_contains_center(bool contains_center);
  void set_edge(int i, int edge);

  // All fields are packed into 16 bytes (assuming 64-bit pointers).  Up to
  // two edge ids are stored inline; this is an important optimization for
  // clients that use S2Shapes consisting of a single edge.
  int32_t shape_id_;
  uint32_t contains_center_ : 1;  // shape contains the cell center
  // TODO(user): Use in-class initializer when C++20 is allowed in
  // opensource version.
  uint32_t num_edges_ : 31;

  // The maximum number of edges that we can store inline in the union.
  static constexpr int kMaxInlineEdges = 2;

  // If there are more than two edges, this field holds a pointer.
  // Otherwise it holds an array of edge ids.
  union {
    int32_t* edges_;  // Owned by the containing S2ShapeIndexCell.
    std::array<int32_t, kMaxInlineEdges> inline_edges_;
  };
};

// S2ShapeIndexCell stores the index contents for a particular S2CellId.
// It consists of a set of clipped shapes.
class S2ShapeIndexCell {
 public:
  S2ShapeIndexCell() = default;
  ~S2ShapeIndexCell();

  // Returns the number of clipped shapes in this cell.
  int num_clipped() const { return shapes_.size(); }

  // Returns the clipped shape at the given index.  Shapes are kept sorted in
  // increasing order of shape id.
  //
  // REQUIRES: 0 <= i < num_clipped()
  const S2ClippedShape& clipped(int i) const { return shapes_[i]; }

  // Returns a pointer to the clipped shape corresponding to the given shape,
  // or nullptr if the shape does not intersect this cell.
  const S2ClippedShape* find_clipped(int shape_id) const;

  // Returns a read-only span over the clipped shapes in the cell.
  absl::Span<const S2ClippedShape> clipped_shapes() const {
    if (!shapes_.empty()) {
      return {shapes_.begin(), shapes_.size()};
    }
    return {};
  }

  // Convenience method that returns the total number of edges in all clipped
  // shapes.
  int num_edges() const;

  // Appends an encoded representation of the S2ShapeIndexCell to "encoder".
  // "num_shape_ids" should be set to index.num_shape_ids(); this information
  // allows the encoding to be more compact in some cases.
  //
  // REQUIRES: "encoder" uses the default constructor, so that its buffer
  //           can be enlarged as necessary by calling Ensure(int).
  void Encode(int num_shape_ids, Encoder* encoder) const;

  // Decodes an S2ShapeIndexCell, returning true on success.
  // "num_shape_ids" should be set to index.num_shape_ids().
  bool Decode(int num_shape_ids, Decoder* decoder);

 private:
  friend class MutableS2ShapeIndex;
  friend class EncodedS2ShapeIndex;
  friend class S2Stats;

  // Internal methods are documented with their definitions.
  S2ClippedShape* add_shapes(int n);
  static void EncodeEdges(const S2ClippedShape& clipped, Encoder* encoder);
  static bool DecodeEdges(int num_edges, S2ClippedShape* clipped,
                          Decoder* decoder);

  using S2ClippedShapeSet = gtl::compact_array<S2ClippedShape>;
  S2ClippedShapeSet shapes_;

  S2ShapeIndexCell(const S2ShapeIndexCell&) = delete;
  void operator=(const S2ShapeIndexCell&) = delete;
};

// S2ShapeIndex is an abstract base class for indexing polygonal geometry in
// memory.  The objects in the index are known as "shapes", and may consist of
// points, polylines, and/or polygons, possibly overlapping.  The index makes
// it very fast to answer queries such as finding nearby shapes, measuring
// distances, testing for intersection and containment, etc.
//
// Each object in the index implements the S2Shape interface.  An S2Shape is a
// collection of edges that optionally defines an interior.  The edges do not
// need to be connected, so for example an S2Shape can represent a polygon
// with multiple shells and/or holes, or a set of polylines, or a set of
// points.  All geometry within a single S2Shape must have the same dimension,
// so for example if you want to create an S2ShapeIndex containing a polyline
// and 10 points, then you will need at least two different S2Shape objects.
//
// There are two important types of S2ShapeIndex.  MutableS2ShapeIndex allows
// you to build an index incrementally by adding or removing shapes, whereas
// EncodedS2ShapeIndex works very efficiently with existing indexes by keeping
// the index data in encoded form (see introduction at the top of this file).
//
// Code that only needs read-only ("const") access to an index should use the
// S2ShapeIndex base class as the parameter type, so that it will work with
// any S2ShapeIndex subtype.  For example:
//
//   void DoSomething(const S2ShapeIndex& index) {
//     ... works with MutableS2ShapeIndex or EncodedS2ShapeIndex ...
//   }
//
// There are a number of built-in classes that work with S2ShapeIndex objects.
// Generally these classes accept any collection of geometry that can be
// represented by an S2ShapeIndex, i.e. any combination of points, polylines,
// and polygons.  Such classes include:
//
// - S2ContainsPointQuery: returns the shape(s) that contain a given point.
//
// - S2ClosestEdgeQuery: returns the closest edge(s) to a given point, edge,
//                       S2CellId, or S2ShapeIndex.
//
// - S2CrossingEdgeQuery: returns the edge(s) that cross a given edge.
//
// - S2BooleanOperation: computes boolean operations such as union,
//                       and boolean predicates such as containment.
//
// - S2ShapeIndexRegion: can be used together with S2RegionCoverer to
//                       approximate geometry as a set of S2CellIds.
//
// - S2ShapeIndexBufferedRegion: computes approximations that have been
//                               expanded by a given radius.
//
// Here is an example showing how to index a set of polygons and then
// determine which polygon(s) contain each of a set of query points:
//
//   void TestContainment(const vector<S2Point>& points,
//                        const vector<S2Polygon*>& polygons) {
//     MutableS2ShapeIndex index;
//     for (auto polygon : polygons) {
//       index.Add(std::make_unique<S2Polygon::Shape>(polygon));
//     }

//     auto query = MakeS2ContainsPointQuery(&index);
//     for (const auto& point : points) {
//       for (int shape_id : query.GetContainingShapeIds(point)) {
//         S2Polygon* polygon = polygons[shape_id];
//         ... do something with (point, polygon) ...
//       }
//     }
//   }
//
// This example uses S2Polygon::Shape, which is one example of an S2Shape
// object.  S2Polyline and S2Loop also have nested Shape classes, and there are
// additional S2Shape types defined in *_shape.h.
//
// Internally, an S2ShapeIndex is essentially a map from S2CellIds to the set
// of shapes that intersect each S2CellId.  It is adaptively refined to ensure
// that no cell contains more than a small number of edges.
//
// In addition to implementing a shared set of virtual methods, all
// S2ShapeIndex subtypes define an Iterator type with the same API.  This
// makes it easy to convert code that uses a particular S2ShapeIndex subtype
// to instead use the abstract base class (or vice versa).  You can also
// choose to avoid the overhead of virtual method calls by making the
// S2ShapeIndex type a template argument, like this:
//
//   template <class IndexType>
//   void DoSomething(const IndexType& index) {
//     for (typename IndexType::Iterator it(&index, S2ShapeIndex::BEGIN);
//          !it.done(); it.Next()) {
//       ...
//     }
//   }
//
// The S2ShapeIndex subtypes provided by the S2 library are thread-compatible,
// meaning that const methods may be called concurrently from multiple threads,
// while non-const methods require exclusive access to the S2ShapeIndex.
class S2ShapeIndex {
 public:
  // Each subtype of S2ShapeIndex should define an Iterator type derived
  // from the following base class.
  class IteratorBase : public S2CellIterator {
   public:
    // Returns a reference to the contents of the current index cell.
    // REQUIRES: !done()
    virtual const S2ShapeIndexCell& cell() const = 0;

    // Returns a newly allocated copy of this iterator.
    virtual std::unique_ptr<IteratorBase> Clone() const = 0;

   protected:
    // Protect the default constructor and move/copy constructors and operators.
    // This allows sub-classes to still be default/copy/move-constructible but
    // prevents accidental slicing through a base class pointer.
    IteratorBase() = default;
    IteratorBase(const IteratorBase&) = default;
    IteratorBase& operator=(const IteratorBase&) = default;
    IteratorBase(IteratorBase&&) = default;
    IteratorBase& operator=(IteratorBase&&) = default;
  };

  // A type function to check if a type is derived from S2ShapeIndex.  This is
  // useful for writing static checks on template parameters when we want to
  // inline a particular iterator call, but we need to make sure it implements
  // the interface that we want.  We don't have access to c++ concepts, so this
  // is the next best thing:
  //
  //   template <typename Index>
  //   void Process(Index& iter) {
  //     static_assert(S2ShapeIndex::ImplementedBy<Index>{},
  //       "We require an S2ShapeIndex.");
  //   }
  template <typename T>
  using ImplementedBy = std::is_convertible<std::decay_t<T>*, S2ShapeIndex*>;

  virtual ~S2ShapeIndex() = default;

  // Returns the number of distinct shape ids in the index.  This is the same
  // as the number of shapes provided that no shapes have ever been removed.
  // (Shape ids are never reused.)
  virtual int num_shape_ids() const = 0;

  // Returns a pointer to the shape with the given id, or nullptr if the shape
  // has been removed from the index.
  virtual const S2Shape* shape(int id) const = 0;

  // Stores an encoded representation of the index into the given encoder.
  virtual void Encode(Encoder* encoder) const = 0;

  // Allows iterating over the indexed shapes using range-based for loops:
  //
  //   for (S2Shape* shape : index) { ... }
  //
  // CAVEAT: Returns nullptr for shapes that have been removed from the index.
  class ShapeIterator {
   public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = const S2Shape*;
    using difference_type = std::ptrdiff_t;
    using pointer = const S2Shape**;
    using reference = const S2Shape*&;

    ShapeIterator() = default;
    const S2Shape* operator*() const;
    ShapeIterator& operator++();
    ShapeIterator operator++(int);

    // REQUIRES: "it" and *this must reference the same S2ShapeIndex.
    bool operator==(ShapeIterator it) const;

    // REQUIRES: "it" and *this must reference the same S2ShapeIndex.
    bool operator!=(ShapeIterator it) const;

   private:
    friend class S2ShapeIndex;
    ShapeIterator(const S2ShapeIndex* index, int shape_id)
        : index_(index), shape_id_(shape_id) {}

    const S2ShapeIndex* index_ = nullptr;
    int shape_id_ = 0;
  };
  ShapeIterator begin() const;
  ShapeIterator end() const;

  // Returns the number of bytes currently occupied by the index (including any
  // unused space at the end of vectors, etc).
  virtual size_t SpaceUsed() const = 0;

  // Minimizes memory usage by requesting that any data structures that can be
  // rebuilt should be discarded.  This method invalidates all iterators.
  //
  // Like all non-const methods, this method is not thread-safe.
  virtual void Minimize() = 0;

  // When passed to an Iterator constructor, specifies whether the iterator
  // should be positioned at the beginning of the index (BEGIN), the end of
  // the index (END), or arbitrarily (UNPOSITIONED).  By default iterators are
  // unpositioned, since this avoids an extra seek in this situation where one
  // of the seek methods (such as Locate) is immediately called.
  enum InitialPosition { BEGIN, END, UNPOSITIONED };

  // A random access iterator that provides low-level access to the cells of
  // the index.  Cells are sorted in increasing order of S2CellId.
  class Iterator final : public IteratorBase {
   public:
    // Default constructor; must be followed by a call to Init().
    Iterator() = default;

    // Constructs an iterator positioned as specified.  By default iterators
    // are unpositioned, since this avoids an extra seek in this situation
    // where one of the seek methods (such as Locate) is immediately called.
    //
    // If you want to position the iterator at the beginning, e.g. in order to
    // loop through the entire index, do this instead:
    //
    //   for (S2ShapeIndex::Iterator it(&index, S2ShapeIndex::BEGIN);
    //        !it.done(); it.Next()) { ... }
    explicit Iterator(const S2ShapeIndex* index,
                      InitialPosition pos = UNPOSITIONED) {
      Init(index, pos);
    }

    // Initializes an iterator for the given S2ShapeIndex.  This method may
    // also be called in order to restore an iterator to a valid state after
    // the underlying index has been updated (although it is usually easier
    // just to declare a new iterator whenever required, since iterator
    // construction is cheap).
    void Init(const S2ShapeIndex* index,
              InitialPosition pos = UNPOSITIONED) {
      iter_ = index->NewIterator(pos);
    }

    // Iterators are copyable and movable.
    Iterator(const Iterator&);
    Iterator& operator=(const Iterator&);
    Iterator(Iterator&&) = default;
    Iterator& operator=(Iterator&&) = default;

    // Returns the S2CellId of the current index cell.  If done() is true,
    // returns a value larger than any valid S2CellId (S2CellId::Sentinel()).
    S2CellId id() const override { return iter_->id(); }

    // Returns the center point of the cell.
    // REQUIRES: !done()
    S2Point center() const { return id().ToPoint(); }

    // Returns a reference to the contents of the current index cell.
    // REQUIRES: !done()
    const S2ShapeIndexCell& cell() const override { return iter_->cell(); }

    // Returns true if the iterator is positioned past the last index cell.
    bool done() const override { return iter_->done(); }

    // Positions the iterator at the first index cell (if any).
    void Begin() override { iter_->Begin(); }

    // Positions the iterator past the last index cell.
    void Finish() override { iter_->Finish(); }

    // Positions the iterator at the next index cell.
    // REQUIRES: !done()
    void Next() override { iter_->Next(); }

    // If the iterator is already positioned at the beginning, returns false.
    // Otherwise positions the iterator at the previous entry and returns true.
    bool Prev() override { return iter_->Prev(); }

    // Positions the iterator at the first cell with id() >= target, or at the
    // end of the index if no such cell exists.
    void Seek(S2CellId target) override { iter_->Seek(target); }

    // Positions the iterator at the cell containing target and returns true. If
    // no such cell exists, return false and leave the iterator in an undefined
    // (but valid) state.
    bool Locate(const S2Point& target) override {
      return LocateImpl(*this, target);
    }

    // Let T be the target S2CellId.  If T is contained by some index cell I
    // (including equality), this method positions the iterator at I and
    // returns INDEXED.  Otherwise if T contains one or more (smaller) index
    // cells, it positions the iterator at the first such cell I and returns
    // SUBDIVIDED.  Otherwise it returns DISJOINT and leaves the iterator
    // positioned arbitrarily.
    S2CellRelation Locate(S2CellId target) override {
      return LocateImpl(*this, target);
    }

    // Clone our underlying iterator to clone this iterator.
    std::unique_ptr<IteratorBase> Clone() const override {
      return std::make_unique<Iterator>(Iterator(iter_->Clone()));
    };

   private:
    // Although S2ShapeIndex::Iterator can be used to iterate over any
    // index subtype, it is more efficient to use the subtype's iterator when
    // the subtype is known at compile time.  For example, MutableS2ShapeIndex
    // should use a MutableS2ShapeIndex::Iterator.
    //
    // The following declarations prevent accidental use of
    // S2ShapeIndex::Iterator when the actual subtype is known.  (If you
    // really want to do this, you can down_cast the index argument to
    // S2ShapeIndex.)
    template <class T>
    explicit Iterator(const T* index, InitialPosition pos = UNPOSITIONED) {}

    explicit Iterator(std::unique_ptr<IteratorBase> iter)
        : iter_(std::move(iter)) {}

    std::unique_ptr<IteratorBase> iter_;
  };

  // ShapeFactory is an interface for decoding vectors of S2Shapes.  It allows
  // random access to the shapes in order to support lazy decoding.  See
  // s2shapeutil_coding.h for useful subtypes.
  class ShapeFactory {
   public:
    virtual ~ShapeFactory() = default;

    // Returns the number of S2Shapes in the vector.
    virtual int size() const = 0;

    // Returns the S2Shape object corresponding to the given "shape_id".
    // Returns nullptr if a shape cannot be decoded or a shape is missing
    // (e.g., because MutableS2ShapeIndex::Release() was called).
    virtual std::unique_ptr<S2Shape> operator[](int shape_id) const = 0;

    // Returns a deep copy of this ShapeFactory.
    virtual std::unique_ptr<ShapeFactory> Clone() const = 0;
  };

 protected:
  // Returns a new iterator positioned as specified.
  virtual std::unique_ptr<IteratorBase> NewIterator(InitialPosition pos)
      const = 0;

  S2ShapeIndex() = default;
  S2ShapeIndex(const S2ShapeIndex&) = delete;
  S2ShapeIndex(S2ShapeIndex&&) = delete;
  S2ShapeIndex& operator=(const S2ShapeIndex&) = delete;
  S2ShapeIndex& operator=(S2ShapeIndex&&) = delete;
};

//////////////////   Implementation details follow   ////////////////////


inline int S2ClippedShape::shape_id() const {
  return shape_id_;
}

inline bool S2ClippedShape::contains_center() const {
  return contains_center_;
}

inline int S2ClippedShape::num_edges() const {
  return num_edges_;
}

inline int S2ClippedShape::edge(int i) const {
  return is_inline() ? inline_edges_[i] : edges_[i];
}

// Initialize an S2ClippedShape to hold the given number of edges.
inline void S2ClippedShape::Init(int32_t shape_id, uint32_t num_edges) {
  shape_id_ = shape_id;
  num_edges_ = num_edges;
  contains_center_ = false;
  if (!is_inline()) {
    edges_ = new int32_t[num_edges];
  }
}

// Free any memory allocated by this S2ClippedShape.  We don't do this in
// the destructor because S2ClippedShapes are copied by STL code, and we
// don't want to repeatedly copy and free the edge data.  Instead the data
// is owned by the containing S2ShapeIndexCell.
inline void S2ClippedShape::Destruct() {
  if (!is_inline()) delete[] edges_;
}

inline bool S2ClippedShape::is_inline() const {
  return num_edges_ <= kMaxInlineEdges;
}

// Set "contains_center_" to indicate whether this clipped shape contains the
// center of the cell to which it belongs.
inline void S2ClippedShape::set_contains_center(bool contains_center) {
  contains_center_ = contains_center;
}

// Set the i-th edge of this clipped shape to be the given edge of the
// original shape.
inline void S2ClippedShape::set_edge(int i, int edge) {
  if (is_inline()) {
    inline_edges_[i] = edge;
  } else {
    edges_[i] = edge;
  }
}

// Inline because an index cell frequently contains just one shape.
inline int S2ShapeIndexCell::num_edges() const {
  int n = 0;
  for (int i = 0; i < num_clipped(); ++i) n += clipped(i).num_edges();
  return n;
}

inline const S2Shape* S2ShapeIndex::ShapeIterator::operator*() const {
  return index_->shape(shape_id_);
}

inline S2ShapeIndex::ShapeIterator& S2ShapeIndex::ShapeIterator::operator++() {
  ++shape_id_;
  return *this;
}

inline S2ShapeIndex::ShapeIterator S2ShapeIndex::ShapeIterator::operator++(
    int) {
  return ShapeIterator(index_, shape_id_++);
}

inline bool S2ShapeIndex::ShapeIterator::operator==(ShapeIterator it) const {
  ABSL_DCHECK_EQ(index_, it.index_);
  return shape_id_ == it.shape_id_;
}

inline bool S2ShapeIndex::ShapeIterator::operator!=(ShapeIterator it) const {
  ABSL_DCHECK_EQ(index_, it.index_);
  return shape_id_ != it.shape_id_;
}

inline S2ShapeIndex::ShapeIterator S2ShapeIndex::begin() const {
  return ShapeIterator(this, 0);
}

inline S2ShapeIndex::ShapeIterator S2ShapeIndex::end() const {
  return ShapeIterator(this, num_shape_ids());
}

inline S2ShapeIndex::Iterator::Iterator(const Iterator& other)
    : iter_(other.iter_->Clone()) {
}

inline S2ShapeIndex::Iterator& S2ShapeIndex::Iterator::operator=(
    const Iterator& other) {
  iter_ = other.iter_->Clone();
  return *this;
}

#endif  // S2_S2SHAPE_INDEX_H_
