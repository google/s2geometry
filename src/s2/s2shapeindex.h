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
// S2ShapeIndex indexes a set of "shapes", where a shape is a collection of
// edges that optionally defines an interior.  A shape can be as simple as a
// single edge, or as complex as a collection of loops.  For shapes that have
// interiors, the index makes it very fast to determine the shape(s) that
// contain a given point or region.
//
// For example, to index a set of polygons and then determine which polygons
// contain various query points:
//
// void Test(vector<S2Polygon*> const& polygons,
//           vector<S2Point> const& points) {
//   S2ShapeIndex index;
//   for (auto polygon : polygons) {
//     index.Add(new S2Polygon::Shape(polygon));
//   }
//   for (auto const& point: points) {
//     vector<S2Shape*> shapes;
//     index.GetContainingShapes(point, &shapes);
//     for (auto shape: shapes) {
//       Output(point, down_cast<S2Polygon::Shape const*>(shape)->polygon());
//     }
//   }
// }
//
// The index can be updated incrementally by adding or removing shapes
// (although each individual shape is immutable).  It is designed to handle up
// to hundreds of millions of edges.  All data structures are designed to be
// small, so the index is compact; generally it is smaller than the underlying
// data being indexed.  The index is also fast to construct.
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
// void Test(vector<S2Polyline*> const& polylines) {
//   S2ShapeIndex index;
//   for (S2Polyline* polyline : polylines) {
//     index.Add(new S2Polyline::Shape(polyline));
//   }
//   // Now use an S2CrossingEdgeQuery or S2ClosestEdgeQuery here ...
// }

#ifndef S2_S2SHAPEINDEX_H_
#define S2_S2SHAPEINDEX_H_

#include <array>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include "s2/base/atomicops.h"
#include "s2/base/integral_types.h"
#include <glog/logging.h>
#include "s2/base/macros.h"
#include "s2/base/mutex.h"
#include "s2/base/spinlock.h"
#include "s2/util/btree/btree_map.h"  // Like std::map, but faster and smaller.
#include "s2/fpcontractoff.h"
#include "s2/s2.h"
#include "s2/s2cellid.h"
#include "s2/util/gtl/compact_array.h"

class R1Interval;
class S2PaddedCell;
class S2ShapeIndex;

// S2Shape is an abstract base class that defines a shape.  Typically it
// wraps some other geometric object in order to provide access to its edges
// without duplicating the edge data.  Shapes are immutable once they have
// been indexed; to modify a shape it must be removed and then added again.
// A shape can be removed from one S2ShapeIndex and then added to another,
// but it can belong to only one index at a time.
//
// There are various useful subtypes defined in s2shapeutil.h, and also in the
// polygonal geometry classes such as S2Polygon and S2Polyline.
//
// Sometimes you will want to add your own data to shapes, so that when an
// S2Shape is returned by a method (such as GetContainingShapes) then you can
// map it back to your source data.  Here are some reasonable ways to do this:
//
//  - [Easy] Every shape has a unique id() assigned by S2ShapeIndex.  Ids are
//    assigned sequentially starting from 0 in the order the shapes are added
//    to the index.  You can use this id to look up arbitrary data stored in
//    your own vector.
//
//  - [Easy] If all of your shapes are the same type, then you can create your
//    own subclass of some existing S2Shape type (such as S2Polyline::Shape)
//    and add your own methods and fields.  You can access this data by
//    downcasting the S2Shape pointers returned by S2ShapeIndex methods.
//
//  - [More general] If you need shapes of different types, then you can still
//    attach extra data in a uniform way by overriding the user_data() method.
//    This is a virtual method declared in the S2Shape base class that returns
//    an arbitrary pointer.  Because it exists in all S2Shapes, you can
//    override it in each type of shape you need to use, and given a S2Shape
//    pointer, you can call this method without knowing the actual S2Shape
//    subtype.  For example, if you need to use polyline and polygon shapes,
//    you can do this:
//
//      class MyPolyline : public S2Polyline::Shape {
//       public:
//        virtual void* mutable_user_data() { return &my_data_; }
//       private:
//        MyData my_data_;
//      };
//      class MyPolygon : public S2Polygon::Shape { ... };
//      ...
//      S2Shape* shape = index.shape(id);
//      MyData* my_data = static_cast<MyData*>(shape->mutable_user_data());
class S2Shape {
 public:
  S2Shape() : id_(-1) {}
  virtual ~S2Shape() {}

  // Return the number of edges in this shape.  Edges have ids ranging from 0
  // to num_edges() - 1.
  virtual int num_edges() const = 0;

  // Return pointers to the edge endpoints for the given edge id.
  // REQUIRES: 0 <= id < num_edges()
  // PROVIDES: (**a) != (**b), i.e. zero-length edges are not allowed.
  virtual void GetEdge(int id, S2Point const** a, S2Point const** b) const = 0;

  // Return true if this shape has an interior, i.e. the shape consists of one
  // or more closed non-intersecting loops.
  virtual bool has_interior() const = 0;

  // Returns true if this shape contains S2::Origin().  Should return false
  // for shapes that do not have an interior.
  //
  // TODO(ericv): Consider allowing shapes to also return their own choice of
  // origin(), to make this method easier to implement for arbitrary shapes.
  virtual bool contains_origin() const = 0;

  // A unique id assigned to this shape by S2ShapeIndex.  Shape ids are
  // assigned sequentially starting from 0 in the order shapes are added.
  int id() const { return id_; }

  // Virtual methods that return pointers of your choice.  These methods are
  // intended to help with the problem of attaching additional data to S2Shape
  // objects (see overview comments).  For example, you could return a pointer
  // to a source object, or a pointer to a bundle of additional data allocated
  // directly with the S2Shape.  Example usage:
  //
  //  MyData const* data = static_cast<MyData const*>(shape->user_data());
  virtual void const* user_data() const { return nullptr; }
  virtual void* mutable_user_data() { return nullptr; }

 private:
  friend class S2ShapeIndex;
  int id_;  // Assigned by S2ShapeIndex when the shape is added.

  S2Shape(S2Shape const&) = delete;
  void operator=(S2Shape const&) = delete;
};

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
// Currently it consists of a set of clipped shapes.
class S2ShapeIndexCell {
 public:
  // Return the number of clipped shapes in this cell.
  int num_shapes() const { return shapes_.size(); }

  // Return the clipped shape at the given index.  Shapes are kept sorted in
  // increasing order of shape id.
  //
  // REQUIRES: 0 <= i < num_shapes()
  S2ClippedShape const& clipped(int i) const { return shapes_[i]; }

  // Return a pointer to the clipped shape corresponding to the given shape,
  // or nullptr if the shape does not intersect this cell.
  S2ClippedShape const* find_clipped(S2Shape const* shape) const;
  S2ClippedShape const* find_clipped(int shape_id) const;

 private:
  friend class S2ShapeIndex;  // shapes_ write access
  friend class S2Stats;

  // Internal methods are documented with their definitions.
  S2ShapeIndexCell() {}
  ~S2ShapeIndexCell();
  S2ClippedShape* add_shapes(int n);

  using S2ClippedShapeSet = compact_array<S2ClippedShape>;
  S2ClippedShapeSet shapes_;

  S2ShapeIndexCell(S2ShapeIndexCell const&) = delete;
  void operator=(S2ShapeIndexCell const&) = delete;
};

// Options that affect construction of the S2ShapeIndex.
// This class is intended to be copied by value as desired.
class S2ShapeIndexOptions {
 public:
  S2ShapeIndexOptions();

  // The maximum number of edges per cell.  If a cell has more than this many
  // edges that are "long" relative to the cell size, and it is not a leaf
  // cell, then it is subdivided.  (Whether an edge is considered "long" is
  // controlled by the --s2shapeindex_min_cell_size_for_edge flag.)
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

// The shape index is essentially a map from S2CellId to a set of clipped
// shaped that intersect that cell id.  It is adaptively refined to ensure
// that no cell contains more than a small number of edges.
class S2ShapeIndex {
 private:
  using CellMap = util::btree::btree_map<S2CellId, S2ShapeIndexCell*>;

 public:
  // Create an S2ShapeIndex that uses the default option settings.  Option
  // values may be changed by calling Init().
  S2ShapeIndex();

  // Create an S2ShapeIndex with the given options.
  explicit S2ShapeIndex(S2ShapeIndexOptions const& options);

  ~S2ShapeIndex();

  // Initialize an S2ShapeIndex with the given options.  This method may only
  // be called when the index is empty (i.e. newly created or Reset() has
  // just been called).
  void Init(S2ShapeIndexOptions const& options);

  S2ShapeIndexOptions const& options() const { return options_; }

  // The number of distinct shape ids that have been assigned.  This equals
  // the number of shapes in the index provided that no shapes have ever been
  // removed.  (Shape ids are not reused.)
  int num_shape_ids() const { return shapes_.size(); }

  // Return a pointer to the shape with the given id, or nullptr if the shape
  // has been removed from the index.
  S2Shape* shape(int id) const { return shapes_[id]; }

  // Take ownership of the given shape and add it to the index.  Also assigns
  // a unique id to the shape (shape->id()).  Shape ids are assigned
  // sequentially starting from 0 in the order shapes are added.  Invalidates
  // all iterators and their associated data.
  //
  // REQUIRES: "shape" is not currently in any other S2ShapeIndex.
  // REQUIRES: "shape" persists for the lifetime of the index or until
  //           Remove(shape) is called.
  void Add(S2Shape* shape);

  // Remove the given shape from the index and return ownership to the caller.
  // Invalidates all iterators and their associated data.
  void Remove(S2Shape* shape);

  // Resets the index to its original state and returns ownership of all
  // shapes to the caller.  This method is much more efficient than removing
  // all shapes one at a time.
  void RemoveAll();

  // Resets the index to its original state and deletes all shapes.  Any
  // options specified via Init() are preserved.
  void Reset();

  // Return true if "shape" contains the given point P.
  bool ShapeContains(S2Shape const* shape, S2Point const& p) const;

  // Return true if the given point P is contained by at least one shape,
  // and return a list of the containing shapes.
  bool GetContainingShapes(S2Point const& p,
                           std::vector<S2Shape*>* shapes) const;

  // The possible relationships between a "target" cell and the cells of the
  // S2ShapeIndex.  If the target is an index cell or is contained by an index
  // cell, it is "INDEXED".  If the target is subdivided into one or more
  // index cells, it is "SUBDIVIDED".  Otherwise it is "DISJOINT".
  enum CellRelation {
    INDEXED,       // Target is contained by an index cell
    SUBDIVIDED,    // Target is subdivided into one or more index cells
    DISJOINT       // Target does not intersect any index cells
  };

  // A random access iterator that provides low-level access to the cells of
  // the index.  Cells are sorted in increasing order of S2CellId.
  //
  // This class is intended to be copied by value as desired.
  class Iterator {
   public:
    // Default constructor; must be followed by a call to Init().
    Iterator();

    // Convenience constructor that calls Init().
    explicit Iterator(S2ShapeIndex const& index);

    // Initialize an iterator for the given S2ShapeIndex.  If the index is
    // non-empty, the iterator is positioned at the first cell.  This method
    // may also be called in order to restore an iterator to a valid state
    // after the underlying index has been updated (although it is usually
    // easier just to declare a new iterator whenever required, since iterator
    // construction is cheap).
    void Init(S2ShapeIndex const& index);

    // Reset the iterator to its original state (positioned at the first cell
    // in the index).  Note that this method does *not* restore the iterator
    // to a valid state after the underlying index has been updated; to do
    // that you must call Init().
    void Reset();

    // The cell id for this cell.
    S2CellId id() const;

    // Pointer to the cell contents.
    S2ShapeIndexCell const* cell() const;

    // The center of the cell (used as a reference point for shape interiors).
    S2Point center() const;

    // Advance the iterator to the next cell in the index.
    // REQUIRES: !Done()
    void Next();

    // Position the iterator at the previous cell in the index.
    // REQUIRES: !AtBegin()
    void Prev();

    // Return true if the iterator is positioned past the last index cell.
    bool Done() const;

    // Return true if the iterator is positioned at the first index cell.
    bool AtBegin() const;

    // Position the iterator at the first cell with id() >= target, or at the
    // end of the index if no such cell exists.
    void Seek(S2CellId target);

    // Advance the iterator to the next cell with id() >= target.  If the
    // iterator is Done() or already satisfies id() >= target, do nothing.
    void SeekForward(S2CellId target);

    // Position the iterator so that Done() is true.
    void Finish();

    // Position the iterator at the cell containing "target".  If no such cell
    // exists, returns false and leaves the iterator positioned arbitrarily.
    // The returned index cell is guaranteed to contain all edges that might
    // intersect the line segment between "target" and the cell center.
    bool Locate(S2Point const& target);

    // Let T be the target S2CellId.  If T is contained by some index cell I
    // (including equality), then position the iterator at I and return
    // INDEXED.  Otherwise if T contains one or more (smaller) index cells,
    // then position the iterator at the first such cell I and return
    // SUBDIVIDED.  Otherwise return DISJOINT and leave the iterator
    // positioned arbitrarily.
    CellRelation Locate(S2CellId target);

    // Initialize an iterator for the given S2ShapeIndex without applying any
    // pending updates.  This can be used to observe the actual current state
    // of the index without modifying it in any way.
    void InitStale(S2ShapeIndex const& index);

   private:
    S2ShapeIndex const* index_;
    CellMap::const_iterator iter_, end_;
  };

  // Return the total number of edges in all indexed shapes.  This method
  // takes time linear in the number of shapes.
  int GetNumEdges() const;

  // Return the number of bytes occupied by the index (including any unused
  // space at the end of vectors, etc).  This method applies any pending
  // updates before measuring memory usage.  It has the same thread safety as
  // the other "const" methods (see introduction).
  size_t BytesUsed() const;

  // Calls to Add() and Remove() are normally queued and processed on the
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

 private:
  friend class Iterator;           // cell_map_
  friend class S2ShapeIndexTest;
  friend class S2Stats;

  class BatchDescriptor;
  class ClippedEdge;
  class EdgeAllocator;
  class FaceEdge;
  class InteriorTracker;
  class RemovedShape;

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
  void ReserveSpace(BatchDescriptor const& batch,
                    std::vector<FaceEdge> all_edges[6]) const;
  void AddShape(int id, std::vector<FaceEdge> all_edges[6],
                InteriorTracker* tracker) const;
  void RemoveShape(RemovedShape const& removed,
                   std::vector<FaceEdge> all_edges[6],
                   InteriorTracker* tracker) const;
  void AddFaceEdge(FaceEdge* edge, std::vector<FaceEdge> all_edges[6]) const;
  void UpdateFaceEdges(int face, std::vector<FaceEdge> const& face_edges,
                       InteriorTracker* tracker);
  S2CellId ShrinkToFit(S2PaddedCell const& pcell, R2Rect const& bound) const;
  void SkipCellRange(S2CellId begin, S2CellId end, InteriorTracker* tracker,
                     EdgeAllocator* alloc, bool disjoint_from_index);
  void UpdateEdges(S2PaddedCell const& pcell,
                   std::vector<ClippedEdge const*>* edges,
                   InteriorTracker* tracker, EdgeAllocator* alloc,
                   bool disjoint_from_index);
  void AbsorbIndexCell(S2PaddedCell const& pcell,
                       Iterator const& iter,
                       std::vector<ClippedEdge const*>* edges,
                       InteriorTracker* tracker,
                       EdgeAllocator* alloc);
  int GetEdgeMaxLevel(S2Point const& a, S2Point const& b) const;
  static int CountShapes(std::vector<ClippedEdge const*> const& edges,
                         ShapeIdSet const& cshape_ids);
  bool MakeIndexCell(S2PaddedCell const& pcell,
                     std::vector<ClippedEdge const*> const& edges,
                     InteriorTracker* tracker);
  static void TestAllEdges(std::vector<ClippedEdge const*> const& edges,
                           InteriorTracker* tracker);
  inline static ClippedEdge const* UpdateBound(ClippedEdge const* edge,
                                               int u_end, double u,
                                               int v_end, double v,
                                               EdgeAllocator* alloc);
  static ClippedEdge const* ClipUBound(ClippedEdge const* edge,
                                       int u_end, double u,
                                       EdgeAllocator* alloc);
  static ClippedEdge const* ClipVBound(ClippedEdge const* edge,
                                       int v_end, double v,
                                       EdgeAllocator* alloc);
  static void ClipVAxis(ClippedEdge const* edge, R1Interval const& middle,
                        std::vector<ClippedEdge const*> child_edges[2],
                        EdgeAllocator* alloc);

  // The amount by which cells are "padded" to compensate for numerical errors
  // when clipping line segments to cell boundaries.
  static double const kCellPadding;

  // The shapes in the index, accessed by their shape id.  Removed shapes are
  // replaced by nullptr pointers.
  std::vector<S2Shape*> shapes_;

  // A map from S2CellId to the set of clipped shapes that intersect that
  // cell.  The cell ids cover a set of non-overlapping regions on the
  // sphere.  Note that this field is updated lazily (see below).  Const
  // methods *must* call MaybeApplyUpdates() before accessing this field.
  // (The easiest way to achieve this is simply to use an Iterator.)
  CellMap cell_map_;

  // The options supplied for this index.
  S2ShapeIndexOptions options_;

  // The id of the first shape that has been queued for addition but not
  // processed yet.
  int pending_additions_begin_;

  // The representation of an edge that has been queued for removal.
  struct RemovedShape {
    int32 shape_id;
    bool has_interior;
    bool contains_origin;
    std::vector<std::pair<S2Point, S2Point>> edges;
  };

  // The set of shapes that have been queued for removal but not processed
  // yet.  Note that we need to copy the edge data since the caller is free to
  // destroy the shape once Remove() has been called.  This field is present
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
  Atomic32 index_status_;  // Stores an IndexStatus

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
  void UnlockAndSignal();

  S2ShapeIndex(S2ShapeIndex const&) = delete;
  void operator=(S2ShapeIndex const&) = delete;
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

inline S2ClippedShape const* S2ShapeIndexCell::find_clipped(
    S2Shape const* shape) const {
  return find_clipped(shape->id());
}

inline S2ShapeIndex::Iterator::Iterator() : index_(nullptr) {
}
inline S2ShapeIndex::Iterator::Iterator(S2ShapeIndex const& index) {
  Init(index);
}
inline void S2ShapeIndex::Iterator::Init(S2ShapeIndex const& index) {
  index.MaybeApplyUpdates();
  InitStale(index);
}
inline void S2ShapeIndex::Iterator::InitStale(S2ShapeIndex const& index) {
  index_ = &index;
  iter_ = index_->cell_map_.begin();
  end_ = index_->cell_map_.end();
}
inline S2CellId S2ShapeIndex::Iterator::id() const {
  DCHECK(!Done());
  return iter_->first;
}
inline S2ShapeIndexCell const* S2ShapeIndex::Iterator::cell() const {
  DCHECK(!Done());
  return iter_->second;
}
inline void S2ShapeIndex::Iterator::Next() {
  DCHECK(!Done());
  ++iter_;
}
inline void S2ShapeIndex::Iterator::Prev() {
  DCHECK(!AtBegin());
  --iter_;
}
inline bool S2ShapeIndex::Iterator::Done() const {
  return iter_ == end_;
}
inline bool S2ShapeIndex::Iterator::AtBegin() const {
  return iter_ == index_->cell_map_.begin();
}
inline void S2ShapeIndex::Iterator::Reset() {
  // Make sure that the index has not been modified since Init() was called.
  DCHECK(index_->is_fresh());
  iter_ = index_->cell_map_.begin();
  end_ = index_->cell_map_.end();
}
inline void S2ShapeIndex::Iterator::Seek(S2CellId target) {
  iter_ = index_->cell_map_.lower_bound(target);
}
inline void S2ShapeIndex::Iterator::SeekForward(S2CellId target) {
  if (!Done() && id() < target) Seek(target);
}
inline void S2ShapeIndex::Iterator::Finish() {
  iter_ = end_;
}

inline bool S2ShapeIndex::is_fresh() const {
  return base::subtle::NoBarrier_Load(&index_status_) == FRESH;
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
  if (base::subtle::Acquire_Load(&index_status_) != FRESH) {
    const_cast<S2ShapeIndex*>(this)->ApplyUpdatesThreadSafe();
  }
}

#endif  // S2_S2SHAPEINDEX_H_
