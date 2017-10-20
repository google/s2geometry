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
// S2ShapeIndex indexes a set of shapes.  A shape is a collection of edges
// that optionally defines an interior.  It can be used to represent a set of
// points, a set of polylines, or a set of polygons.  The shapes in the index
// are allowed to intersect arbitrarily.  The index makes it very fast to
// answer queries such as finding the closest shape(s) to a given point,
// determing which shape(s) contain a given point, determining whether two
// S2ShapeIndexes intersect, etc.
//
// For example, to index a set of polygons and then determine which polygons
// contain various query points:
//
// void Test(const vector<S2Polygon*>& polygons,
//           const vector<S2Point>& points) {
//   S2ShapeIndex index;
//   for (auto polygon : polygons) {
//     index.Add(absl::make_unique<S2Polygon::Shape>(polygon));
//   }
//   for (const auto& point: points) {
//     MakeS2ContainsPointQuery(&index).VisitContainingShapes(
//         point, [](S2Shape* shape) {
//           Output(point, down_cast<S2Polygon::Shape*>(shape)->polygon());
//         });
//   }
// }
//
// TODO(ericv): There are currently two S2ShapeIndex implementations, both
// derived from a common base class (S2ShapeIndexBase).  S2ShapeIndex itself
// represents a mutable index that can be updated incrementally by adding or
// removing shapes.  It also has an Encode method that allows the index to be
// serialized.  An encoded S2ShapeIndex can be decoded either into its
// original form (S2ShapeIndex) or into an EncodedS2ShapeIndex (the other
// subtype of S2ShapeIndexBase).  The main advantage of EncodedS2ShapeIndex is
// that it can be constructed instantaneously, since the index is kept in its
// original encoded form.  Data is decoded only when an operation needs it.
// For example, to determine which shapes(s) contain a given query point only
// requires decoding the data in the S2ShapeIndexCell that contains that point.
//
// S2ShapeIndex is designed to handle up to hundreds of millions of edges.
// All data structures are designed to be small, so the index is compact;
// generally it is smaller than the underlying data being indexed.  The index
// is also fast to construct.
//
// All "const" methods are thread-safe provided that they do not overlap with
// calls to non-const methods.  Non-const methods are not thread-safe.  This
// means that if you update the index, you need to ensure that no other thread
// is reading or updating the index (including through Iterator objects).
//
// S2Polygon, S2Loop, and S2Polyline define S2Shape classes that allow these
// objects to be indexed easily.  Additional S2Shape classes are defined in
// s2shapeutil.  You can find useful query methods in S2CrossingEdgeQuery and
// S2ClosestEdgeQuery.
//
// Example showing how to build an index of S2Polylines:
//
// void Test(const vector<S2Polyline*>& polylines) {
//   S2ShapeIndex index;
//   for (S2Polyline* polyline : polylines) {
//     index.Add(absl::make_unique<S2Polyline::Shape>(polyline));
//   }
//   // Now use an S2CrossingEdgeQuery or S2ClosestEdgeQuery here ...
// }

#ifndef S2_S2SHAPEINDEX_H_
#define S2_S2SHAPEINDEX_H_

#include <array>
#include <atomic>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "s2/third_party/absl/base/integral_types.h"
#include <glog/logging.h>
#include "s2/third_party/absl/base/macros.h"
#include "s2/base/mutex.h"
#include "s2/base/spinlock.h"
#include "s2/third_party/absl/base/thread_annotations.h"
#include "s2/third_party/absl/memory/memory.h"
#include "s2/util/btree/btree_map.h"
#include "s2/_fpcontractoff.h"
#include "s2/s2cellid.h"
#include "s2/s2pointutil.h"
#include "s2/s2shape.h"
#include "s2/util/gtl/compact_array.h"

class R1Interval;
class S2PaddedCell;
class S2ShapeIndex;

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
  // The shape id of the clipped shape.
  int shape_id() const;

  // Return true if the center of the S2CellId is inside the shape.  Returns
  // false for shapes that do not have an interior.
  bool contains_center() const;

  // The number of edges that intersect the S2CellId.
  int num_edges() const;

  // Return the edge id of the given edge in this clipped shape.  Edges are
  // sorted in increasing order of edge id.
  //
  // REQUIRES: 0 <= i < num_edges()
  int edge(int i) const;

  // Return true if the clipped shape contains the given edge id.
  bool ContainsEdge(int id) const;

 private:
  // This class may be copied by value, but note that it does *not* own its
  // underlying data.  (It is owned by the containing S2ShapeIndexCell.)

  friend class S2ShapeIndex;      // Init(), etc.
  friend class S2ShapeIndexCell;  // Destruct()
  friend class S2Stats;

  // Internal methods are documented with their definition.
  void Init(int32 shape_id, int32 num_edges);
  void Destruct();
  bool is_inline() const;
  void set_contains_center(bool contains_center);
  void set_edge(int i, int edge);

  // All fields are packed into 16 bytes (assuming 64-bit pointers).  Up to
  // two edge ids are stored inline; this is an important optimization for
  // clients that use S2Shapes consisting of a single edge.
  int32 shape_id_;
  uint32 contains_center_ : 1;  // shape contains the cell center
  uint32 num_edges_ : 31;

  // If there are more than two edges, this field holds a pointer.
  // Otherwise it holds an array of edge ids.
  union {
    int32* edges_;           // This pointer is owned by the containing Cell.
    std::array<int32, 2> inline_edges_;
  };
};

// S2ShapeIndexCell stores the index contents for a particular S2CellId.
// It consists of a set of clipped shapes.
class S2ShapeIndexCell {
 public:
  // Return the number of clipped shapes in this cell.
  int num_clipped() const { return shapes_.size(); }

  // Return the clipped shape at the given index.  Shapes are kept sorted in
  // increasing order of shape id.
  //
  // REQUIRES: 0 <= i < num_clipped()
  const S2ClippedShape& clipped(int i) const { return shapes_[i]; }

  // Return a pointer to the clipped shape corresponding to the given shape,
  // or nullptr if the shape does not intersect this cell.
  const S2ClippedShape* find_clipped(const S2Shape* shape) const;
  const S2ClippedShape* find_clipped(int shape_id) const;

  // Convenience method that returns the total number of edges in all clipped
  // shapes.
  int num_edges() const;

 private:
  friend class S2ShapeIndex;  // shapes_ write access
  friend class S2Stats;

  // Internal methods are documented with their definitions.
  S2ShapeIndexCell() {}
  ~S2ShapeIndexCell();
  S2ClippedShape* add_shapes(int n);

  using S2ClippedShapeSet = compact_array<S2ClippedShape>;
  S2ClippedShapeSet shapes_;

  S2ShapeIndexCell(const S2ShapeIndexCell&) = delete;
  void operator=(const S2ShapeIndexCell&) = delete;
};

// S2ShapeIndexBase is an abstract base class for all S2ShapeIndex
// implementations.  An S2ShapeIndex is essentially a map from S2CellIds to
// the set of clipped shapes that intersect each S2CellId.  It is adaptively
// refined to ensure that no cell contains more than a small number of edges.
//
// In addition to implementing a shared set of virtual methods, all
// S2ShapeIndex subtypes define an Iterator type with the same API.  This
// makes it easy to convert code that uses a particular S2ShapeIndex subtype
// to instead use the abstract base class (or vice versa).  It also makes it
// possible to avoid the overhead of virtual method calls if desired by
// making the S2ShapeIndex subtype a template argument.
//
// TODO(ericv/jrosenstock): Rename this class to S2ShapeIndex.
class S2ShapeIndexBase {
 protected:
  class IteratorBase;

 public:
  virtual ~S2ShapeIndexBase() {}

  // Returns the number of distinct shape ids in the index.  This is the same
  // as the number of shapes provided that no shapes have ever been removed.
  // (Shape ids are never reused.)
  virtual int num_shape_ids() const = 0;

  // Returns a pointer to the shape with the given id, or nullptr if the shape
  // has been removed from the index.
  virtual S2Shape* shape(int id) const = 0;

  // Minimizes memory usage by requesting that any data structures that can be
  // rebuilt should be discarded.  This method invalidates all iterators.
  //
  // Like all non-const methods, this method is not thread-safe.
  virtual void Minimize() = 0;

  // The possible relationships between a "target" cell and the cells of the
  // S2ShapeIndex.  If the target is an index cell or is contained by an index
  // cell, it is "INDEXED".  If the target is subdivided into one or more
  // index cells, it is "SUBDIVIDED".  Otherwise it is "DISJOINT".
  enum CellRelation {
    INDEXED,       // Target is contained by an index cell
    SUBDIVIDED,    // Target is subdivided into one or more index cells
    DISJOINT       // Target does not intersect any index cells
  };

  // When passed to an Iterator constructor, specifies whether the iterator
  // should be positioned at the beginning of the index (BEGIN), the end of
  // the index (END), or arbitrarily (UNPOSITIONED).  By default iterators are
  // unpositioned, since this avoids an extra seek in this situation where one
  // of the seek methods (such as Locate) is immediately called.
  enum InitialPosition { BEGIN, END, UNPOSITIONED };

  // A random access iterator that provides low-level access to the cells of
  // the index.  Cells are sorted in increasing order of S2CellId.
  class Iterator {
   public:
    // Default constructor; must be followed by a call to Init().
    Iterator() : iter_(nullptr) {}

    // Constructs an iterator positioned as specified.  By default iterators
    // are unpositioned, since this avoids an extra seek in this situation
    // where one of the seek methods (such as Locate) is immediately called.
    //
    // If you want to position the iterator at the beginning, e.g. in order to
    // loop through the entire index, do this instead:
    //
    //   for (S2ShapeIndexBase::Iterator it(&index, S2ShapeIndex::BEGIN);
    //        !it.done(); it.Next()) { ... }
    explicit Iterator(const S2ShapeIndexBase* index,
                      InitialPosition pos = UNPOSITIONED)
        : iter_(index->NewIterator(pos)) {}

    // Initializes an iterator for the given S2ShapeIndex.  This method may
    // also be called in order to restore an iterator to a valid state after
    // the underlying index has been updated (although it is usually easier
    // just to declare a new iterator whenever required, since iterator
    // construction is cheap).
    void Init(const S2ShapeIndexBase* index,
              InitialPosition pos = UNPOSITIONED) {
      iter_ = index->NewIterator(pos);
    }

    // Iterators are copyable and moveable.
    Iterator(const Iterator&);
    Iterator& operator=(const Iterator&);
    Iterator(Iterator&&);
    Iterator& operator=(Iterator&&);

    // Returns the S2CellId of the current index cell.  If done() is true,
    // returns a value larger than any valid S2CellId (S2CellId::Sentinel()).
    S2CellId id() const { return iter_->id(); }

    // Returns the center point of the cell.
    // REQUIRES: !done()
    S2Point center() const { return id().ToPoint(); }

    // Returns a reference to the contents of the current index cell.
    // REQUIRES: !done()
    const S2ShapeIndexCell& cell() const { return iter_->cell(); }

    // Returns true if the iterator is positioned past the last index cell.
    bool done() const { return iter_->done(); }

    // Positions the iterator at the first index cell (if any).
    void Begin() { iter_->Begin(); }

    // Positions the iterator past the last index cell.
    void Finish() { iter_->Finish(); }

    // Positions the iterator at the next index cell.
    // REQUIRES: !done()
    void Next() { iter_->Next(); }

    // If the iterator is already positioned at the beginning, returns false.
    // Otherwise positions the iterator at the previous entry and returns true.
    bool Prev() { return iter_->Prev(); }

    // Positions the iterator at the first cell with id() >= target, or at the
    // end of the index if no such cell exists.
    void Seek(S2CellId target) { iter_->Seek(target); }

    // Positions the iterator at the cell containing "target".  If no such cell
    // exists, returns false and leaves the iterator positioned arbitrarily.
    // The returned index cell is guaranteed to contain all edges that might
    // intersect the line segment between "target" and the cell center.
    bool Locate(const S2Point& target) {
      return IteratorBase::LocateImpl(target, this);
    }

    // Let T be the target S2CellId.  If T is contained by some index cell I
    // (including equality), this method positions the iterator at I and
    // returns INDEXED.  Otherwise if T contains one or more (smaller) index
    // cells, it positions the iterator at the first such cell I and returns
    // SUBDIVIDED.  Otherwise it returns DISJOINT and leaves the iterator
    // positioned arbitrarily.
    CellRelation Locate(S2CellId target) {
      return IteratorBase::LocateImpl(target, this);
    }

   private:
    // Although S2ShapeIndexBase::Iterator can be used to iterate over any
    // index subtype, it is more efficient to use the subtype's iterator when
    // the subtype is known at compile time.  For example, MutableS2ShapeIndex
    // should use a MutableS2ShapeIndex::Iterator.
    //
    // The following declarations prevent accidental use of
    // S2ShapeIndexBase::Iterator with known subtypes.  (If you really want to
    // do this, you can down_cast the index argument to S2ShapeIndexBase.)
    template <class T>
    explicit Iterator(const T* index, InitialPosition pos = UNPOSITIONED) {}

    template <class T>
    void Init(const T* index, InitialPosition pos = UNPOSITIONED) {}

    std::unique_ptr<IteratorBase> iter_;
  };

 protected:
  // Each subtype of S2ShapeIndexBase should define an Iterator type derived
  // from the following base class.
  class IteratorBase {
   public:
    virtual ~IteratorBase() {}

    IteratorBase(const IteratorBase&);
    IteratorBase& operator=(const IteratorBase&);

    // Returns the S2CellId of the current index cell.  If done() is true,
    // returns a value larger than any valid S2CellId (S2CellId::Sentinel()).
    S2CellId id() const;

    // Returns the center point of the cell.
    // REQUIRES: !done()
    S2Point center() const;

    // Returns a reference to the contents of the current index cell.
    // REQUIRES: !done()
    const S2ShapeIndexCell& cell() const;

    // Returns true if the iterator is positioned past the last index cell.
    bool done() const;

    // Positions the iterator at the first index cell (if any).
    virtual void Begin() = 0;

    // Positions the iterator past the last index cell.
    virtual void Finish() = 0;

    // Positions the iterator at the next index cell.
    // REQUIRES: !done()
    virtual void Next() = 0;

    // If the iterator is already positioned at the beginning, returns false.
    // Otherwise positions the iterator at the previous entry and returns true.
    virtual bool Prev() = 0;

    // Positions the iterator at the first cell with id() >= target, or at the
    // end of the index if no such cell exists.
    virtual void Seek(S2CellId target) = 0;

    // Positions the iterator at the cell containing "target".  If no such cell
    // exists, returns false and leaves the iterator positioned arbitrarily.
    // The returned index cell is guaranteed to contain all edges that might
    // intersect the line segment between "target" and the cell center.
    virtual bool Locate(const S2Point& target) = 0;

    // Let T be the target S2CellId.  If T is contained by some index cell I
    // (including equality), this method positions the iterator at I and
    // returns INDEXED.  Otherwise if T contains one or more (smaller) index
    // cells, it positions the iterator at the first such cell I and returns
    // SUBDIVIDED.  Otherwise it returns DISJOINT and leaves the iterator
    // positioned arbitrarily.
    virtual CellRelation Locate(S2CellId target) = 0;

   protected:
    IteratorBase() : id_(S2CellId::Sentinel()), cell_(nullptr) {}

    // Sets the iterator state.  "cell" typically points to the cell contents,
    // but may also be given as "nullptr" in order to implement decoding on
    // demand.  In that situation, the first that the client attempts to
    // access the cell contents, the GetCell() method is called and "cell_" is
    // updated in a thread-safe way.
    void set_state(S2CellId id, const S2ShapeIndexCell* cell);

    // Sets the iterator state so that done() is true.
    void set_finished();

    // Returns the current contents of the "cell_" field, which may be nullptr
    // if the cell contents have not been decoded yet.
    const S2ShapeIndexCell* raw_cell() const;

    // This method is called to decode the contents of the current cell, if
    // set_state() was previously called with a nullptr "cell" argument.  This
    // allows decoding on demand for subtypes that keep the cell contents in
    // an encoded state.  It does not need to be implemented at all if
    // set_state() is always called with (cell != nullptr).
    //
    // REQUIRES: This method is thread-safe.
    // REQUIRES: Multiple calls to this method return the same value.
    virtual const S2ShapeIndexCell* GetCell() const = 0;

    // Returns an exact copy of this iterator.
    virtual std::unique_ptr<IteratorBase> Clone() const = 0;

    // Makes a copy of the given source iterator.
    // REQUIRES: "other" has the same concrete type as "this".
    virtual void Copy(const IteratorBase& other) = 0;

    // The default implementation of Locate(S2Point).  It is instantiated by
    // each subtype in order to (1) minimize the number of virtual method
    // calls (since subtypes are typically "final") and (2) ensure that the
    // correct versions of non-virtual methods such as cell() are called.
    template <class Iter>
    static bool LocateImpl(const S2Point& target, Iter* it);

    // The default implementation of Locate(S2CellId) (see comments above).
    template <class Iter>
    static CellRelation LocateImpl(S2CellId target, Iter* it);

   private:
    friend class Iterator;

    // This method is "const" because it is used internally by "const" methods
    // in order to implement decoding on demand.
    void set_cell(const S2ShapeIndexCell* cell) const;

    S2CellId id_;
    mutable std::atomic<const S2ShapeIndexCell*> cell_;
  };

  // Returns a new iterator positioned as specified.
  virtual std::unique_ptr<IteratorBase> NewIterator(InitialPosition pos)
      const = 0;
};

// S2ShapeIndex is essentially a map from S2CellIds to the set of clipped
// shaped that intersect each S2CellId.  It is adaptively refined to ensure
// that no cell contains more than a small number of edges.
//
// TODO(ericv/jrosenstock): Rename this class to MutableS2ShapeIndex, and
// probably split into a different file.
class S2ShapeIndex final : public S2ShapeIndexBase {
 private:
  using CellMap = util::btree::btree_map<S2CellId, S2ShapeIndexCell*>;

 public:
  // Options that affect construction of the S2ShapeIndex.
  class Options {
   public:
    Options();

    // The maximum number of edges per cell.  If a cell has more than this
    // many edges that are not considered "long" relative to the cell size,
    // then it is subdivided.  (Whether an edge is considered "long" is
    // controlled by --s2shapeindex_cell_size_to_long_edge_ratio flag.)
    //
    // Values between 10 and 50 represent a reasonable balance between memory
    // usage, construction time, and query time.  Small values make queries
    // faster, while large values make construction faster and use less memory.
    // Values higher than 50 do not save significant additional memory, and
    // query times can increase substantially, especially for algorithms that
    // visit all pairs of potentially intersecting edges (such as polygon
    // validation), since this is quadratic in the number of edges per cell.
    //
    // Note that the *average* number of edges per cell is generally slightly
    // less than half of the maximum value defined here.
    //
    // Defaults to value given by --s2shapeindex_default_max_edges_per_cell.
    int max_edges_per_cell() const { return max_edges_per_cell_; }
    void set_max_edges_per_cell(int max_edges_per_cell);

   private:
    int max_edges_per_cell_;
  };

  // Create an S2ShapeIndex that uses the default option settings.  Option
  // values may be changed by calling Init().
  S2ShapeIndex();

  // Create an S2ShapeIndex with the given options.
  explicit S2ShapeIndex(const Options& options);

  ~S2ShapeIndex();

  // Initialize an S2ShapeIndex with the given options.  This method may only
  // be called when the index is empty (i.e. newly created or Reset() has
  // just been called).
  void Init(const Options& options);

  const Options& options() const { return options_; }

  // The number of distinct shape ids that have been assigned.  This equals
  // the number of shapes in the index provided that no shapes have ever been
  // removed.  (Shape ids are not reused.)
  int num_shape_ids() const override { return shapes_.size(); }

  // Return a pointer to the shape with the given id, or nullptr if the shape
  // has been removed from the index.
  S2Shape* shape(int id) const override { return shapes_[id].get(); }

  // Minimizes memory usage by requesting that any data structures that can be
  // rebuilt should be discarded.  This method invalidates all iterators.
  //
  // Like all non-const methods, this method is not thread-safe.
  void Minimize() override;

  class Iterator final : public IteratorBase {
   public:
    // Default constructor; must be followed by a call to Init().
    Iterator();

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
                      InitialPosition pos = UNPOSITIONED);

    // Initializes an iterator for the given S2ShapeIndex.  This method may
    // also be called in order to restore an iterator to a valid state after
    // the underlying index has been updated (although it is usually easier
    // just to declare a new iterator whenever required, since iterator
    // construction is cheap).
    void Init(const S2ShapeIndex* index, InitialPosition pos = UNPOSITIONED);

    // Initialize an iterator for the given S2ShapeIndex without applying any
    // pending updates.  This can be used to observe the actual current state
    // of the index without modifying it in any way.
    void InitStale(const S2ShapeIndex* index,
                   InitialPosition pos = UNPOSITIONED);

    // Inherited non-virtual methods:
    //   S2CellId id() const;
    //   bool done() const;
    //   S2Point center() const;
    const S2ShapeIndexCell& cell() const;

    // IteratorBase API:
    void Begin() override;
    void Finish() override;
    void Next() override;
    bool Prev() override;
    void Seek(S2CellId target) override;
    bool Locate(const S2Point& target) override;
    CellRelation Locate(S2CellId target) override;

   protected:
    const S2ShapeIndexCell* GetCell() const override;
    std::unique_ptr<IteratorBase> Clone() const override;
    void Copy(const IteratorBase& other) override;

   private:
    void Refresh();  // Updates the IteratorBase fields.
    const S2ShapeIndex* index_;
    CellMap::const_iterator iter_, end_;
  };

  // Takes ownership of the given shape and adds it to the index.  Also assigns
  // a unique id to the shape (shape->id()) and returns that id.  Shape ids
  // are assigned sequentially starting from 0 in the order shapes are added.
  // Invalidates all iterators and their associated data.
  int Add(std::unique_ptr<S2Shape> shape);

  // Remove the given shape from the index and return ownership to the caller.
  // Invalidates all iterators and their associated data.
  std::unique_ptr<S2Shape> Release(int shape_id);

  // Resets the index to its original state and returns ownership of all
  // shapes to the caller.  This method is much more efficient than removing
  // all shapes one at a time.
  std::vector<std::unique_ptr<S2Shape>> ReleaseAll();

  // Resets the index to its original state and deletes all shapes.  Any
  // options specified via Init() are preserved.
  void Clear();

  // Returns the number of bytes currently occupied by the index (including any
  // unused space at the end of vectors, etc). It has the same thread safety
  // as the other "const" methods (see introduction).
  size_t SpaceUsed() const;

  // Calls to Add() and Release() are normally queued and processed on the
  // first subsequent query (in a thread-safe way).  This has many advantages,
  // the most important of which is that sometimes there *is* no subsequent
  // query, which lets us avoid building the index completely.
  //
  // This method forces any pending updates to be applied immediately.
  // Calling this method is rarely a good idea.  (One valid reason is to
  // exclude the cost of building the index from benchmark results.)
  void ForceApplyUpdates();

  // Return true if there are no pending updates that need to be applied.
  // This can be useful to avoid building the index unnecessarily, or for
  // choosing between two different algorithms depending on whether the index
  // is available.
  //
  // The returned index status may be slightly out of date if the index was
  // built in a different thread.  This is fine for the intended use (as an
  // efficiency hint), but it should not be used by internal methods  (see
  // MaybeApplyUpdates).
  bool is_fresh() const;

  ABSL_DEPRECATED("Use Clear()")
  void Reset() { Clear(); }

 protected:
  std::unique_ptr<IteratorBase> NewIterator(InitialPosition pos) const override;

 private:
  friend class Iterator;           // cell_map_
  friend class S2ShapeIndexTest;
  friend class S2Stats;

  class BatchDescriptor;
  class ClippedEdge;
  class EdgeAllocator;
  class FaceEdge;
  class InteriorTracker;
  struct RemovedShape;

  using ShapeIdSet = std::vector<int>;

  // Internal methods are documented with their definitions.
  bool is_first_update() const;
  bool is_shape_being_removed(int shape_id) const;
  void MaybeApplyUpdates() const;
  void ApplyUpdatesThreadSafe();
  void ApplyUpdatesInternal();
  void GetUpdateBatches(std::vector<BatchDescriptor>* batches) const;
  static void GetBatchSizes(int num_items, int max_batches,
                            double final_bytes_per_item,
                            double high_water_bytes_per_item,
                            double preferred_max_bytes_per_batch,
                            std::vector<int>* batch_sizes);
  void ReserveSpace(const BatchDescriptor& batch,
                    std::vector<FaceEdge> all_edges[6]) const;
  void AddShape(int id, std::vector<FaceEdge> all_edges[6],
                InteriorTracker* tracker) const;
  void RemoveShape(const RemovedShape& removed,
                   std::vector<FaceEdge> all_edges[6],
                   InteriorTracker* tracker) const;
  void AddFaceEdge(FaceEdge* edge, std::vector<FaceEdge> all_edges[6]) const;
  void UpdateFaceEdges(int face, const std::vector<FaceEdge>& face_edges,
                       InteriorTracker* tracker);
  S2CellId ShrinkToFit(const S2PaddedCell& pcell, const R2Rect& bound) const;
  void SkipCellRange(S2CellId begin, S2CellId end, InteriorTracker* tracker,
                     EdgeAllocator* alloc, bool disjoint_from_index);
  void UpdateEdges(const S2PaddedCell& pcell,
                   std::vector<const ClippedEdge*>* edges,
                   InteriorTracker* tracker, EdgeAllocator* alloc,
                   bool disjoint_from_index);
  void AbsorbIndexCell(const S2PaddedCell& pcell,
                       const Iterator& iter,
                       std::vector<const ClippedEdge*>* edges,
                       InteriorTracker* tracker,
                       EdgeAllocator* alloc);
  int GetEdgeMaxLevel(const S2Shape::Edge& edge) const;
  static int CountShapes(const std::vector<const ClippedEdge*>& edges,
                         const ShapeIdSet& cshape_ids);
  bool MakeIndexCell(const S2PaddedCell& pcell,
                     const std::vector<const ClippedEdge*>& edges,
                     InteriorTracker* tracker);
  static void TestAllEdges(const std::vector<const ClippedEdge*>& edges,
                           InteriorTracker* tracker);
  inline static const ClippedEdge* UpdateBound(const ClippedEdge* edge,
                                               int u_end, double u,
                                               int v_end, double v,
                                               EdgeAllocator* alloc);
  static const ClippedEdge* ClipUBound(const ClippedEdge* edge,
                                       int u_end, double u,
                                       EdgeAllocator* alloc);
  static const ClippedEdge* ClipVBound(const ClippedEdge* edge,
                                       int v_end, double v,
                                       EdgeAllocator* alloc);
  static void ClipVAxis(const ClippedEdge* edge, const R1Interval& middle,
                        std::vector<const ClippedEdge*> child_edges[2],
                        EdgeAllocator* alloc);

  // The amount by which cells are "padded" to compensate for numerical errors
  // when clipping line segments to cell boundaries.
  static const double kCellPadding;

  // The shapes in the index, accessed by their shape id.  Removed shapes are
  // replaced by nullptr pointers.
  std::vector<std::unique_ptr<S2Shape>> shapes_;

  // A map from S2CellId to the set of clipped shapes that intersect that
  // cell.  The cell ids cover a set of non-overlapping regions on the
  // sphere.  Note that this field is updated lazily (see below).  Const
  // methods *must* call MaybeApplyUpdates() before accessing this field.
  // (The easiest way to achieve this is simply to use an Iterator.)
  CellMap cell_map_;

  // The options supplied for this index.
  Options options_;

  // The id of the first shape that has been queued for addition but not
  // processed yet.
  int pending_additions_begin_ = 0;

  // The representation of an edge that has been queued for removal.
  struct RemovedShape {
    int32 shape_id;
    bool has_interior;
    bool contains_tracker_origin;
    std::vector<S2Shape::Edge> edges;
  };

  // The set of shapes that have been queued for removal but not processed
  // yet.  Note that we need to copy the edge data since the caller is free to
  // destroy the shape once Release() has been called.  This field is present
  // only when there are removed shapes to process (to save memory).
  std::unique_ptr<std::vector<RemovedShape>> pending_removals_;

  // Additions and removals are queued and processed on the first subsequent
  // query.  There are several reasons to do this:
  //
  //  - It is significantly more efficient to process updates in batches.
  //  - Often the index will never be queried, in which case we can save both
  //    the time and memory required to build it.  Examples:
  //     + S2Loops that are created simply to pass to an S2Polygon.  (We don't
  //       need the S2Loop index, because S2Polygon builds its own index.)
  //     + Applications that load a database of geometry and then query only
  //       a small fraction of it.
  //     + Applications that only read and write geometry (Decode/Encode).
  //
  // The main drawback is that we need to go to some extra work to ensure that
  // "const" methods are still thread-safe.  Note that the goal is *not* to
  // make this class thread-safe in general, but simply to hide the fact that
  // we defer some of the indexing work until query time.
  //
  // The textbook approach to this problem would be to use a mutex and a
  // condition variable.  Unfortunately pthread mutexes are huge (40 bytes).
  // Instead we use spinlock (which is only 4 bytes) to guard a few small
  // fields representing the current update status, and only create additional
  // state while the update is actually occurring.
  mutable SpinLock lock_;

  enum IndexStatus {
    STALE,     // There are pending updates.
    UPDATING,  // Updates are currently being applied.
    FRESH,     // There are no pending updates.
  };
  // Reads and writes to this field are guarded by "lock_".
  std::atomic<IndexStatus> index_status_;

  // UpdateState holds temporary data related to thread synchronization.  It
  // is only allocated while updates are being applied.
  struct UpdateState {
    // This mutex is used as a condition variable.  It is locked by the
    // updating thread for the entire duration of the update; other threads
    // lock it in order to wait until the update is finished.
    Mutex wait_mutex;

    // The number of threads currently waiting on "wait_mutex_".  The
    // UpdateState can only be freed when this number reaches zero.
    //
    // Reads and writes to this field are guarded by "lock_".
    int num_waiting;

    UpdateState() : num_waiting(0) {
    }

    ~UpdateState() {
      DCHECK_EQ(0, num_waiting);
    }
  };
  std::unique_ptr<UpdateState> update_state_;

  // Documented in the .cc file.
  void UnlockAndSignal()
      UNLOCK_FUNCTION(lock_)
      UNLOCK_FUNCTION(update_state_->wait_mutex);

  S2ShapeIndex(const S2ShapeIndex&) = delete;
  void operator=(const S2ShapeIndex&) = delete;
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

inline bool S2ClippedShape::is_inline() const {
  return num_edges_ <= inline_edges_.size();
}

inline const S2ClippedShape* S2ShapeIndexCell::find_clipped(
    const S2Shape* shape) const {
  return find_clipped(shape->id());
}

// Inline because an index cell frequently contains just one shape.
inline int S2ShapeIndexCell::num_edges() const {
  int n = 0;
  for (int i = 0; i < num_clipped(); ++i) n += clipped(i).num_edges();
  return n;
}

inline S2ShapeIndexBase::IteratorBase::IteratorBase(const IteratorBase& other)
    : id_(other.id_), cell_(other.raw_cell()) {
}

inline S2ShapeIndexBase::IteratorBase&
S2ShapeIndexBase::IteratorBase::operator=(const IteratorBase& other) {
  id_ = other.id_;
  set_cell(other.raw_cell());
  return *this;
}

inline S2CellId S2ShapeIndexBase::IteratorBase::id() const {
  return id_;
}

inline const S2ShapeIndexCell& S2ShapeIndexBase::IteratorBase::cell() const {
  DCHECK(!done());
  auto cell = raw_cell();
  if (cell == nullptr) {
    cell = GetCell();
    set_cell(cell);
  }
  return *cell;
}

inline bool S2ShapeIndexBase::IteratorBase::done() const {
  return id_ == S2CellId::Sentinel();
}

inline S2Point S2ShapeIndexBase::IteratorBase::center() const {
  DCHECK(!done());
  return id().ToPoint();
}

inline void S2ShapeIndexBase::IteratorBase::set_state(
    S2CellId id, const S2ShapeIndexCell* cell) {
  id_ = id;
  set_cell(cell);
}

inline void S2ShapeIndexBase::IteratorBase::set_finished() {
  id_ = S2CellId::Sentinel();
  set_cell(nullptr);
}

inline const S2ShapeIndexCell* S2ShapeIndexBase::IteratorBase::raw_cell()
    const {
  return cell_.load(std::memory_order_relaxed);
}

inline void S2ShapeIndexBase::IteratorBase::set_cell(
    const S2ShapeIndexCell* cell) const {
  cell_.store(cell, std::memory_order_relaxed);
}

template <class Iter>
inline bool S2ShapeIndexBase::IteratorBase::LocateImpl(
    const S2Point& target_point, Iter* it) {
  // Let I = cell_map_->lower_bound(T), where T is the leaf cell containing
  // "target_point".  Then if T is contained by an index cell, then the
  // containing cell is either I or I'.  We test for containment by comparing
  // the ranges of leaf cells spanned by T, I, and I'.

  S2CellId target(target_point);
  it->Seek(target);
  if (!it->done() && it->id().range_min() <= target) return true;
  if (it->Prev() && it->id().range_max() >= target) return true;
  return false;
}

template <class Iter>
inline S2ShapeIndexBase::CellRelation
S2ShapeIndexBase::IteratorBase::LocateImpl(S2CellId target, Iter* it) {
  // Let T be the target, let I = cell_map_->lower_bound(T.range_min()), and
  // let I' be the predecessor of I.  If T contains any index cells, then T
  // contains I.  Similarly, if T is contained by an index cell, then the
  // containing cell is either I or I'.  We test for containment by comparing
  // the ranges of leaf cells spanned by T, I, and I'.

  it->Seek(target.range_min());
  if (!it->done()) {
    if (it->id() >= target && it->id().range_min() <= target) return INDEXED;
    if (it->id() <= target.range_max()) return SUBDIVIDED;
  }
  if (it->Prev() && it->id().range_max() >= target) return INDEXED;
  return DISJOINT;
}

inline S2ShapeIndexBase::Iterator::Iterator(const Iterator& other)
    : iter_(other.iter_->Clone()) {
}

inline S2ShapeIndexBase::Iterator& S2ShapeIndexBase::Iterator::operator=(
    const Iterator& other) {
  iter_->Copy(*other.iter_);
  return *this;
}

inline S2ShapeIndexBase::Iterator::Iterator(Iterator&& other)
    : iter_(std::move(other.iter_)) {
}

inline S2ShapeIndexBase::Iterator& S2ShapeIndexBase::Iterator::operator=(
    Iterator&& other) {
  iter_ = std::move(other.iter_);
  return *this;
}

inline S2ShapeIndex::Iterator::Iterator() : index_(nullptr) {
}

inline S2ShapeIndex::Iterator::Iterator(const S2ShapeIndex* index,
                                        InitialPosition pos) {
  Init(index, pos);
}

inline void S2ShapeIndex::Iterator::Init(const S2ShapeIndex* index,
                                         InitialPosition pos) {
  index->MaybeApplyUpdates();
  InitStale(index, pos);
}

inline void S2ShapeIndex::Iterator::InitStale(const S2ShapeIndex* index,
                                              InitialPosition pos) {
  index_ = index;
  end_ = index_->cell_map_.end();
  if (pos == BEGIN) {
    iter_ = index_->cell_map_.begin();
  } else {
    iter_ = end_;
  }
  Refresh();
}

inline const S2ShapeIndexCell& S2ShapeIndex::Iterator::cell() const {
  // Since S2ShapeIndex always sets the "cell_" field, we can skip the logic
  // in the base class that conditionally calls GetCell().
  return *raw_cell();
}

inline void S2ShapeIndex::Iterator::Refresh() {
  if (iter_ == end_) {
    set_finished();
  } else {
    set_state(iter_->first, iter_->second);
  }
}

inline void S2ShapeIndex::Iterator::Begin() {
  // Make sure that the index has not been modified since Init() was called.
  DCHECK(index_->is_fresh());
  iter_ = index_->cell_map_.begin();
  end_ = index_->cell_map_.end();
  Refresh();
}

inline void S2ShapeIndex::Iterator::Finish() {
  iter_ = end_;
  Refresh();
}

inline void S2ShapeIndex::Iterator::Next() {
  DCHECK(!done());
  ++iter_;
  Refresh();
}

inline bool S2ShapeIndex::Iterator::Prev() {
  if (iter_ == index_->cell_map_.begin()) return false;
  --iter_;
  Refresh();
  return true;
}

inline void S2ShapeIndex::Iterator::Seek(S2CellId target) {
  iter_ = index_->cell_map_.lower_bound(target);
  Refresh();
}

inline std::unique_ptr<S2ShapeIndex::IteratorBase> S2ShapeIndex::NewIterator(
    InitialPosition pos) const {
  return absl::make_unique<Iterator>(this, pos);
}

inline bool S2ShapeIndex::is_fresh() const {
  return index_status_.load(std::memory_order_relaxed) == FRESH;
}

// Return true if this is the first update to the index.
inline bool S2ShapeIndex::is_first_update() const {
  // Note that it is not sufficient to check whether cell_map_ is empty, since
  // entries are added during the update process.
  return pending_additions_begin_ == 0;
}

// Given that the given shape is being updated, return true if it is being
// removed (as opposed to being added).
inline bool S2ShapeIndex::is_shape_being_removed(int shape_id) const {
  // All shape ids being removed are less than all shape ids being added.
  return shape_id < pending_additions_begin_;
}

// Ensure that any pending updates have been applied.  This method must be
// called before accessing the cell_map_ field, even if the index_status_
// appears to be FRESH, because a memory barrier is required in order to
// ensure that all the index updates are visible if the updates were done in
// another thread.
inline void S2ShapeIndex::MaybeApplyUpdates() const {
  // To avoid acquiring and releasing the spinlock on every query, we use
  // atomic operations when testing whether the status is FRESH and when
  // updating the status to be FRESH.  This guarantees that any thread that
  // sees a status of FRESH will also see the corresponding index updates.
  if (index_status_.load(std::memory_order_acquire) != FRESH) {
    const_cast<S2ShapeIndex*>(this)->ApplyUpdatesThreadSafe();
  }
}

#endif  // S2_S2SHAPEINDEX_H_
