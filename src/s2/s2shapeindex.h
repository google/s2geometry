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
// points, a set of polylines, or a set of polygons.  For shapes that have
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
#include <atomic>
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>

#include <glog/logging.h>

#include "s2/base/mutex.h"
#include "s2/base/spinlock.h"
#include "s2/_fpcontractoff.h"
#include "s2/s2cellid.h"
#include "s2/third_party/absl/base/integral_types.h"
#include "s2/third_party/absl/base/macros.h"
#include "s2/util/btree/btree_map.h"
#include "s2/util/gtl/compact_array.h"

class R1Interval;
class S2PaddedCell;
class S2ShapeIndex;

// An S2Shape is an abstract base class that defines a shape.  A shape is a
// collection of edges that optionally defines an interior.  It can be used to
// represent a set of points, a set of polylines, or a set of polygons.
//
// The edges of an S2Shape are indexed by a contiguous range of "edge ids"
// starting at 0.  The edges are further subdivided into "chains", where each
// chain consists of a sequence of edges connected end-to-end (a polyline).
// S2Shape has methods that allow edges to be accessed either using the global
// numbering (edge id) or within a particular chain.  The global numbering is
// sufficient for most purposes, but the chain representation is useful for
// certain algorithms such as intersection (see S2BoundaryOperation).
//
// Typically an S2Shape wraps some other geometric object in order to provide
// access to its edges without duplicating the edge data.  Shapes are
// immutable once they have been indexed; to modify a shape it must be
// removed from the index and then added again.  A shape can be removed from
// one S2ShapeIndex and then added to another, but it can belong to only one
// index at a time.
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
//  - [More general] If you need shapes of different types in the same index,
//    you can still attach extra data in a uniform way by overriding the
//    user_data() method.  This is a virtual method declared in the S2Shape
//    base class that returns an arbitrary pointer.  Because it exists in all
//    S2Shapes, you can override it in each type of shape you have and call it
//    without knowing the concrete S2Shape subtype.  For example, if you have
//    polyline and polygon shapes, you can do this:
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
  // An edge, consisting of two vertices "v0" and "v1".  Zero-length edges are
  // allowed, and can be used to represent points.
  struct Edge {
    S2Point v0, v1;
    Edge() = default;
    Edge(S2Point const& _v0, S2Point const& _v1) : v0(_v0), v1(_v1) {}

    // TODO(ericv): Convert to functions, define all 6 comparisons.
    bool operator==(Edge const& other) const {
      return v0 == other.v0 && v1 == other.v1;
    }
    bool operator<(Edge const& other) const {
      return v0 < other.v0 || (v0 == other.v0 && v1 < other.v1);
    }
  };

  // A range of edge ids corresponding to a chain of connected edges,
  // specified as a (start, length) pair.  The chain is defined to consist of
  // edge ids {start, start + 1, ..., start + length - 1}.
  struct Chain {
    int32 start, length;
    Chain() = default;
    Chain(int32 _start, int32 _length) : start(_start), length(_length) {}
  };

  // The position of an edge within a given edge chain, specified as a
  // (chain_id, offset) pair.  Chains are numbered sequentially starting from
  // zero, and offsets are measured from the start of each chain.
  struct ChainPosition {
    int32 chain_id, offset;
    ChainPosition() = default;
    ChainPosition(int32 _chain_id, int32 _offset)
        : chain_id(_chain_id), offset(_offset) {}
  };

  S2Shape() : id_(-1) {}
  virtual ~S2Shape() {}

  // Return the number of edges in this shape.  Edges have ids ranging from 0
  // to num_edges() - 1.
  virtual int num_edges() const = 0;

  // Return the endpoints of the given edge id.
  //
  // REQUIRES: 0 <= id < num_edges()
  virtual Edge edge(int edge_id) const = 0;

  // Returns the dimension of the geometry represented by this shape.
  //
  //  0 - Point geometry.  Each point is represented as a degenerate edge.
  //
  //  1 - Polyline geometry.  Polyline edges may be degenerate.  A shape may
  //      represent any number of polylines.  Polylines edges may intersect.
  //
  //  2 - Polygon geometry.  The edges may be returned in any order, but it
  //      must be possible to assemble them into a collection of non-crossing
  //      loops such that the polygon interior is always on the left.
  //      Polygons may have degeneracies (e.g., degenerate edges or sibling
  //      pairs consisting of an edge and its corresponding reversed edge).
  //
  // Note that this method allows degenerate geometry of different dimensions
  // to be distinguished, e.g. it allows a point to be distinguished from a
  // polyline or polygon that has been simplified to a single point.
  virtual int dimension() const = 0;

  // Convenience function that returns true if this shape has an interior.
  inline bool has_interior() const { return dimension() == 2; }

  // Returns true if this shape contains S2::Origin().  Should return false
  // for shapes that do not have an interior.
  //
  // TODO(ericv): Consider allowing shapes to also return their own choice of
  // origin(), to make this method easier to implement for arbitrary shapes.
  virtual bool contains_origin() const = 0;

  // Returns the number of contiguous edge chains in the shape.  For example,
  // a shape whose edges are [AB, BC, CD, AE, EF] would consist of two chains
  // (AB,BC,CD and AE,EF).  Every chain is assigned a "chain id" numbered
  // sequentially starting from zero.
  //
  // Note that it is always acceptable to implement this method by returning
  // num_edges() (i.e. every chain consists of a single edge), but this may
  // reduce the efficiency of some algorithms.
  virtual int num_chains() const = 0;

  // Returns the range of edge ids corresponding to the given edge chain.
  // Edge chains must consist of contiguous, non-overlapping ranges that cover
  // the entire range of edge ids.  This is spelled out more formally below:
  //
  // REQUIRES: 0 <= i < num_chains()
  // REQUIRES: chain(i).length > 0, for all i
  // REQUIRES: chain(0).start == 0
  // REQUIRES: chain(i).start + chain(i).length == chain(i+1).start,
  //           for i < num_chains() - 1
  // REQUIRES: chain(i).start + chain(i).length == num_edges(),
  //           for i == num_chains() - 1
  virtual Chain chain(int chain_id) const = 0;

  // Returns the edge at offset "offset" within edge chain "chain_id".
  // Equivalent to "shape.edge(shape.chain(chain_id).start + offset)"
  // but may be more efficient.
  virtual Edge chain_edge(int chain_id, int offset) const = 0;

  // Finds the chain containing the given edge, and returns the position of
  // that edge as a (chain_id, offset) pair.
  //
  // REQUIRES: shape.chain(pos.chain_id).start + pos.offset == edge_id
  // REQUIRES: shape.chain(pos.chain_id + 1).start > edge_id
  //
  // where     pos == shape.chain_position(edge_id).
  virtual ChainPosition chain_position(int edge_id) const = 0;

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
  int num_clipped() const { return shapes_.size(); }

  // Return the clipped shape at the given index.  Shapes are kept sorted in
  // increasing order of shape id.
  //
  // REQUIRES: 0 <= i < num_clipped()
  S2ClippedShape const& clipped(int i) const { return shapes_[i]; }

  // Return a pointer to the clipped shape corresponding to the given shape,
  // or nullptr if the shape does not intersect this cell.
  S2ClippedShape const* find_clipped(S2Shape const* shape) const;
  S2ClippedShape const* find_clipped(int shape_id) const;

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

  S2ShapeIndexCell(S2ShapeIndexCell const&) = delete;
  void operator=(S2ShapeIndexCell const&) = delete;
};

// Options that affect construction of the S2ShapeIndex.
// This class is intended to be copied by value as desired.
class S2ShapeIndexOptions {
 public:
  S2ShapeIndexOptions();

  // The maximum number of edges per cell.  If a cell has more than this many
  // edges that are not considered "long" relative to the cell size, and it is
  // not a leaf cell, then it is subdivided.  (Whether an edge is considered
  // "long" is controlled by the --s2shapeindex_cell_size_to_long_edge_ratio
  // flag.)
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
  S2Shape* shape(int id) const { return shapes_[id].get(); }

  // Take ownership of the given shape and add it to the index.  Also assigns
  // a unique id to the shape (shape->id()).  Shape ids are assigned
  // sequentially starting from 0 in the order shapes are added.  Invalidates
  // all iterators and their associated data.
  void Add(std::unique_ptr<S2Shape> shape);

  // Remove the given shape from the index and return ownership to the caller.
  // Invalidates all iterators and their associated data.
  std::unique_ptr<S2Shape> Release(int shape_id);

  // Resets the index to its original state and returns ownership of all
  // shapes to the caller.  This method is much more efficient than removing
  // all shapes one at a time.
  std::vector<std::unique_ptr<S2Shape>> ReleaseAll();

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

    // Reference to the cell contents.
    S2ShapeIndexCell const& cell() const;

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

  ABSL_DEPRECATED("Use Add(unique_ptr<S2Shape>) instead.")
  void Add(S2Shape* shape) { Add(std::unique_ptr<S2Shape>(shape)); }

  ABSL_DEPRECATED("Use Release instead.")
  void Remove(S2Shape* shape) { Release(shape->id()).release(); }

  ABSL_DEPRECATED("Use ReleaseAll instead.")
  void RemoveAll() { for (auto& ptr : ReleaseAll()) ptr.release(); }

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
  int GetEdgeMaxLevel(S2Shape::Edge const& edge) const;
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
  std::vector<std::unique_ptr<S2Shape>> shapes_;

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

// Inline because an index cell frequently contains just one shape.
inline int S2ShapeIndexCell::num_edges() const {
  int n = 0;
  for (int i = 0; i < num_clipped(); ++i) n += clipped(i).num_edges();
  return n;
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
inline S2ShapeIndexCell const& S2ShapeIndex::Iterator::cell() const {
  DCHECK(!Done());
  return *iter_->second;
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
