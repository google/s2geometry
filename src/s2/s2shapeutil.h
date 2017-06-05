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

// Author: ericv@google.com (Eric Veach)
//
// This file contains various subclasses of S2Shape:
//
//  - LaxLoop: like S2Loop::Shape but allows duplicate vertices & edges,
//             more compact representation, faster to initialize.
//  - LaxPolygon: like S2Polygon::Shape but allows duplicate vertices & edges,
//                more compact representation, faster to initialize.
//  - LaxPolyline: like S2Polyline::Shape but slightly more compact.
//  - ClosedLaxPolyline: like LaxPolyline but joins last vertex to first.
//  - VertexIdLaxLoop: like LaxLoop, but vertices are specified as indices
//                     into a vertex array.
//  - EdgeVectorShape: represents an arbitrary collection of edges.
//  - S2LoopOwningShape: like S2SLoop::Shape but owns the underlying S2Loop.
//  - S2PolygonOwningShape: like S2SPolygon::Shape but owns the S2Polygon.
//  - S2PolylineOwningShape: like S2SPolyline::Shape but owns the S2Polyline.
//
// It also contains miscellaneous S2ShapeIndex utility functions and classes:
//
//  - RangeIterator: useful for merging two or more S2ShapeIndexes.
//  - VisitCrossings: finds all edge intersections between S2ShapeIndexes.
//  - IsOriginOnLeft: helper for implementing S2Shape::contains_origin.
//
// And functions that are mainly useful internally in the S2 library:
//
//  - ResolveComponents: computes containment relationships among loops.
//  - FindAnyCrossing: helper function for S2Loop and S2Polygon validation.
//
// TODO(ericv): The *OwningShapes could probably be removed by allowing the
// underlying shapes (e.g., S2Polygon::Shape) to take ownership.

#ifndef S2_S2SHAPEUTIL_H_
#define S2_S2SHAPEUTIL_H_

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "s2/fpcontractoff.h"
#include "s2/s2loop.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"
#include "s2/s2shapeindex.h"

class S2Loop;
class S2Error;

namespace s2shapeutil {

/////////////////////////////////////////////////////////////////////////////
///////////////////////////  Utility Shape Types  ///////////////////////////

// LaxLoop represents a closed loop of edges surrounding an interior
// region.  It is similar to S2Loop::Shape except that this class allows
// duplicate vertices and edges.  Loops may have any number of vertices,
// including 0, 1, or 2.  (A one-vertex loop defines a degenerate edge
// consisting of a single point.)
//
// Note that LaxLoop is faster to initialize and more compact than
// S2Loop::Shape, but does not support the same operations as S2Loop.
class LaxLoop : public S2Shape {
 public:
  LaxLoop() {}  // Must call Init().
  explicit LaxLoop(std::vector<S2Point> const& vertices);
  explicit LaxLoop(S2Loop const& loop);  // Copies the loop data.

  // Initialize the shape.
  void Init(std::vector<S2Point> const& vertices);

  // Initialize from an S2Loop (by copying its data).
  // REQUIRES: !loop->is_full()
  void Init(S2Loop const& loop);

  int num_vertices() const { return num_vertices_; }
  S2Point const& vertex(int i) const { return vertices_[i]; }

  // S2Shape interface:
  int num_edges() const final { return num_vertices(); }
  Edge edge(int e) const final;
  void GetEdge(int e, S2Point const** a, S2Point const** b) const final;
  // Not final; overridden by ClosedLaxPolyline.
  int dimension() const override { return 2; }
  // Not final; overridden by ClosedLaxPolyline.
  bool contains_origin() const override;
  int num_chains() const final { return std::min(1, num_vertices_); }
  Chain chain(int i) const final { return Chain(0, num_vertices_); }
  Edge chain_edge(int i, int j) const final;
  ChainPosition chain_position(int e) const final {
    return ChainPosition(0, e);
  }

 private:
  // For clients that have many small loops, we save some memory by
  // representing the vertices as an array rather than using std::vector.
  int32 num_vertices_;
  std::unique_ptr<S2Point[]> vertices_;
};

// LaxPolygon represents a region defined by a collection of zero or more
// closed loops.  The interior is the region to the left of all loops.  This
// is similar to S2Polygon::Shape except that this class supports polygons
// with degeneracies.  Degeneracies are of two types: degenerate edges (from a
// vertex to itself) and sibling edge pairs (consisting of two oppositely
// oriented edges).  Degeneracies can represent either "shells" or "holes"
// depending on the loop they are contained by.  For example, a degenerate
// edge or sibling pair contained by a "shell" would be interpreted as a
// degenerate hole.  Such edges form part of the boundary of the polygon.
//
// Loops with fewer than three vertices are interpreted as follows:
//  - A loop with two vertices defines two edges (in opposite directions).
//  - A loop with one vertex defines a single degenerate edge.
//  - A loop with no vertices is interpreted as the "full loop" containing
//    all points on the sphere.  If this loop is present, then all other loops
//    must form degeneracies (i.e., degenerate edges or sibling pairs).  For
//    example, two loops {} and {X} would be interpreted as the full polygon
//    with a degenerate single-point hole at X.
//
// LaxPolygon does not have any error checking, and it is perfectly fine to
// create LaxPolygon objects that do not meet the requirements below (e.g., in
// order to analyze or fix those problems).  However, LaxPolygons must satisfy
// some additional conditions in order to perform certain operations:
//
//  - In order to be valid for point containment tests, the polygon must
//    satisfy the "interior is on the left" rule.  This means that there must
//    not be any crossing edges, and if there are duplicate edges then all but
//    at most one of thm must belong to a sibling pair (i.e., the number of
//    edges in opposite directions must differ by at most one).
//
//  - To be valid for polygon operations (S2BoundaryOperation), degenerate
//    edges and sibling pairs cannot coincide with any other edges.  For
//    example, the following situations are not allowed:
//
//      {AA, AA}      // degenerate edge coincides with another edge
//      {AA, AB}      // degenerate edge coincides with another edge
//      {AB, BA, AB}  // sibling pair coincides with another edge
//
// Note that LaxPolygon is must faster to initialize and is more compact than
// S2Polygon, but unlike S2Polygon it does not have any built-in operations.
// Instead you should use S2ShapeIndex operations (S2BoundaryOperation,
// S2ClosestEdgeQuery, etc).
class LaxPolygon : public S2Shape {
 public:
  LaxPolygon() {}  // Must call Init().
  using Loop = std::vector<S2Point>;
  explicit LaxPolygon(std::vector<Loop> const& loops);
  explicit LaxPolygon(S2Polygon const& polygon);  // Copies the polygon data.
  ~LaxPolygon() override;

  // Initialize the shape.
  void Init(std::vector<Loop> const& loops);

  // Initialize from an S2Polygon (by copying its data).
  // REQUIRES: !polygon->is_full()
  void Init(S2Polygon const& polygon);

  // Return the number of loops.
  int num_loops() const { return num_loops_; }

  // Return the total number of vertices in all loops.
  int num_vertices() const;

  // Return the number of vertices in the given loop.
  int num_loop_vertices(int i) const;

  // Return the vertex from loop "i" at index "j".
  // REQUIRES: 0 <= i < num_loops()
  // REQUIRES: 0 <= j < num_loop_vertices(i)
  S2Point const& loop_vertex(int i, int j) const;

  // S2Shape interface:
  int num_edges() const final { return num_vertices(); }
  Edge edge(int e) const final;
  void GetEdge(int e, S2Point const** a, S2Point const** b) const final;
  int dimension() const final { return 2; }
  bool contains_origin() const final;
  int num_chains() const final { return num_loops(); }
  Chain chain(int i) const final;
  Edge chain_edge(int i, int j) const final;
  ChainPosition chain_position(int e) const final;

 private:
  class VertexArray;
  void Init(std::vector<VertexArray> const& loops);

  int32 num_loops_;
  std::unique_ptr<S2Point[]> vertices_;
  // If num_loops_ <= 1, this union stores the number of vertices.
  // Otherwise it points to an array of size (num_loops + 1) where element "i"
  // is the total number of vertices in loops 0..i-1.
  union {
    int32 num_vertices_;
    int32* cumulative_vertices_;  // Don't use unique_ptr in unions.
  };
};

// LaxPolyline represents a polyline.  It is similar to S2Polyline::Shape
// except that duplicate vertices are allowed, and the representation is
// slightly more compact.  Polylines may have any number of vertices, but note
// that polylines with fewer than 2 vertices do not define any edges.
class LaxPolyline : public S2Shape {
 public:
  LaxPolyline() {}  // Must call Init().
  explicit LaxPolyline(std::vector<S2Point> const& vertices);
  explicit LaxPolyline(S2Polyline const& polyline);  // Copies the data.

  // Initialize the shape.
  void Init(std::vector<S2Point> const& vertices);

  // Initialize from an S2Polyline (by copying its data).
  void Init(S2Polyline const& polyline);

  int num_vertices() const { return num_vertices_; }
  S2Point const& vertex(int i) const { return vertices_[i]; }

  // S2Shape interface:
  int num_edges() const final { return std::max(0, num_vertices() - 1); }
  Edge edge(int e) const final;
  void GetEdge(int e, S2Point const** a, S2Point const** b) const final;
  int dimension() const final { return 1; }
  bool contains_origin() const final { return false; }
  int num_chains() const final;
  Chain chain(int i) const final;
  Edge chain_edge(int i, int j) const final;
  ChainPosition chain_position(int e) const final;

 private:
  // For clients that have many small polylines, we save some memory by
  // representing the vertices as an array rather than using std::vector.
  int32 num_vertices_;
  std::unique_ptr<S2Point[]> vertices_;
};

// ClosedLaxPolyline is like LaxPolyline except that the last vertex is
// implicitly joined to the first.  It is also like LaxLoop except that it
// does not have an interior (which makes it more efficient to index).
class ClosedLaxPolyline : public LaxLoop {
 public:
  ClosedLaxPolyline() {}  // Must call Init().
  explicit ClosedLaxPolyline(std::vector<S2Point> const& vertices);
  explicit ClosedLaxPolyline(S2Loop const& loop);  // Copies the loop data.
  int dimension() const final { return 1; }
  bool contains_origin() const final { return false; }
};

// EdgeVectorShape is an S2Shape representing an arbitrary set of edges.  It
// is mainly used for testing, but it can also be useful if you have, say, a
// collection of polylines and don't care about memory efficiency (since this
// class would store most of the vertices twice).
//
// Note that if you already have data stored in an S2Loop, S2Polyline, or
// S2Polygon, then you would be better off using the "Shape" class defined
// within those classes (e.g., S2Loop::Shape).  Similarly, if the vertex data
// is stored in your own data structures, you can easily write your own
// subclass of S2Shape that points to the existing vertex data rather than
// copying it.
class EdgeVectorShape : public S2Shape {
 public:
  EdgeVectorShape() {}
  // Convenience constructor for creating a vector of length 1.
  EdgeVectorShape(S2Point const& a, S2Point const& b) {
    edges_.push_back(std::make_pair(a, b));
  }
  // Add an edge to the vector.  IMPORTANT: This method should only be called
  // *before* adding the EdgeVectorShape to an S2ShapeIndex.  S2Shapes can
  // only be modified by removing them from the index, making changes, and
  // then adding them back again.
  void Add(S2Point const& a, S2Point const& b) {
    edges_.push_back(std::make_pair(a, b));
  }

  // S2Shape interface:
  int num_edges() const final { return edges_.size(); }
  Edge edge(int e) const final {
    return Edge(edges_[e].first, edges_[e].second);
  }
  void GetEdge(int e, S2Point const** a, S2Point const** b) const final {
    *a = &edges_[e].first;
    *b = &edges_[e].second;
  }
  int dimension() const final { return 1; }
  bool contains_origin() const final { return false; }
  int num_chains() const final { return edges_.size(); }
  Chain chain(int i) const final { return Chain(i, 1); }
  Edge chain_edge(int i, int j) const final {
    DCHECK_EQ(j, 0);
    return Edge(edges_[i].first, edges_[i].second);
  }
  ChainPosition chain_position(int e) const final {
    return ChainPosition(e, 0);
  }

 private:
  std::vector<std::pair<S2Point, S2Point>> edges_;
};

// PointVectorShape is an S2Shape representing a set of S2Points. Each point
// is reprsented as a degenerate edge with the same starting and ending
// vertices.
//
// This class is useful for adding a collection of points to an S2ShapeIndex.
class PointVectorShape : public S2Shape {
 public:
  ~PointVectorShape() override = default;

  // TODO(user): Change the signature to
  // PointVectorShape(std::vector<S2Point> &&points)
  // once we get an exception to the style rules.

  // Transfers the contents of "points" to this class and leaves "points" empty.
  explicit PointVectorShape(std::vector<S2Point>* points) {
    points_.swap(*points);
  }

  // S2Shape interface:
  int num_edges() const final { return points_.size(); }
  Edge edge(int e) const final {
    return Edge(points_[e], points_[e]);
  }
  void GetEdge(int e, S2Point const** a, S2Point const** b) const final {
    *a = &points_[e];
    *b = &points_[e];
  }
  int dimension() const final { return 0; }
  bool contains_origin() const final { return false; }
  int num_chains() const final { return points_.size(); }
  Chain chain(int i) const final { return Chain(i, 1); }
  Edge chain_edge(int i, int j) const final {
    DCHECK_EQ(j, 0);
    return Edge(points_[i], points_[i]);
  }
  ChainPosition chain_position(int e) const final {
    return ChainPosition(e, 0);
  }

 private:
  std::vector<S2Point> points_;

  PointVectorShape(const PointVectorShape&) = delete;
  PointVectorShape& operator=(const PointVectorShape&) = delete;
};

// VertexIdLaxLoop is just like LaxLoop, except that vertices are
// specified as indices into a vertex array.  This representation can be more
// compact when many loops are arranged in a mesh structure.
class VertexIdLaxLoop : public S2Shape {
 public:
  VertexIdLaxLoop() {}  // Must call Init().
  explicit VertexIdLaxLoop(std::vector<int32> const& vertex_ids,
                           S2Point const* vertex_array);

  // Initialize the shape.  "vertex_ids" is a vector of indices into
  // "vertex_array".
  void Init(std::vector<int32> const& vertex_ids,
            S2Point const* vertex_array);

  // Returns the number of vertices in the loop.
  int num_vertices() const { return num_vertices_; }
  int32 vertex_id(int i) const { return vertex_ids_[i]; }
  S2Point const& vertex(int i) const { return vertex_array_[vertex_id(i)]; }

  // S2Shape interface:
  int num_edges() const final { return num_vertices(); }
  Edge edge(int e) const final;
  void GetEdge(int e, S2Point const** a, S2Point const** b) const final;
  int dimension() const final { return 2; }
  bool contains_origin() const final;
  int num_chains() const final { return 1; }
  Chain chain(int i) const final { return Chain(0, num_vertices_); }
  Edge chain_edge(int i, int j) const final;
  ChainPosition chain_position(int e) const final {
    return ChainPosition(0, e);
  }

 private:
  int32 num_vertices_;
  std::unique_ptr<int32[]> vertex_ids_;
  S2Point const* vertex_array_;
};

// Like S2Loop::Shape, except that the referenced S2Loop is automatically
// deleted when this object is released by the S2ShapeIndex.  This is useful
// when an S2Loop is constructed solely for the purpose of indexing it.
class S2LoopOwningShape : public S2Loop::Shape {
 public:
  S2LoopOwningShape() {}  // Must call Init().
  explicit S2LoopOwningShape(std::unique_ptr<S2Loop const> loop)
      : S2Loop::Shape(loop.release()) {
  }
  ~S2LoopOwningShape() override { delete loop(); }
};

// Like S2Polygon::Shape, except that the referenced S2Polygon is
// automatically deleted when this object is released by the S2ShapeIndex.
// This is useful when an S2Polygon is constructed solely for the purpose of
// indexing it.
class S2PolygonOwningShape : public S2Polygon::Shape {
 public:
  S2PolygonOwningShape() {}  // Must call Init().
  explicit S2PolygonOwningShape(std::unique_ptr<S2Polygon const> polygon)
      : S2Polygon::Shape(polygon.release()) {
  }
  ~S2PolygonOwningShape() override { delete polygon(); }
};

// Like S2Polyline::Shape, except that the referenced S2Polyline is
// automatically deleted when this object is released by the S2ShapeIndex.
// This is useful when an S2Polyline is constructed solely for the purpose of
// indexing it.
class S2PolylineOwningShape : public S2Polyline::Shape {
 public:
  S2PolylineOwningShape() {}  // Must call Init().
  explicit S2PolylineOwningShape(std::unique_ptr<S2Polyline const> polyline)
      : S2Polyline::Shape(polyline.release()) {
  }
  ~S2PolylineOwningShape() override { delete polyline(); }
};


/////////////////////////////////////////////////////////////////////////////
////////////////////// Utility Functions and Classes ////////////////////////

// RangeIterator is a wrapper over S2ShapeIndex::Iterator with extra methods
// that are useful for merging the contents of two or more S2ShapeIndexes.
class RangeIterator {
 public:
  // Construct a new RangeIterator positioned at the first cell of the index.
  explicit RangeIterator(S2ShapeIndex const& index);

  // The current S2CellId and cell contents.
  S2CellId id() const { return id_; }
  S2ShapeIndexCell const& cell() const { return it_.cell(); }

  // The min and max leaf cell ids covered by the current cell.  If Done() is
  // true, these methods return a value larger than any valid cell id.
  S2CellId range_min() const { return range_min_; }
  S2CellId range_max() const { return range_max_; }

  void Next();
  bool Done() { return id_ == end_; }

  // Position the iterator at the first cell that overlaps or follows
  // "target", i.e. such that range_max() >= target.range_min().
  void SeekTo(RangeIterator const& target);

  // Position the iterator at the first cell that follows "target", i.e. the
  // first cell such that range_min() > target.range_max().
  void SeekBeyond(RangeIterator const& target);

 private:
  // Updates internal state after the iterator has been repositioned.
  void Refresh();
  S2ShapeIndex::Iterator it_;
  S2CellId const end_;
  S2CellId id_, range_min_, range_max_;
  S2ShapeIndexCell const* cell_;
};

// A parameter that controls the reporting of edge intersections.
//
//  - CrossingType::INTERIOR reports intersections that occur at a point
//    interior to both edges (i.e., not at a vertex).
//
//  - CrossingType::ALL reports all intersections, even those where two edges
//    intersect only because they share a common vertex.
//
// - CrossingType::NON_ADJACENT reports all intersections except for pairs of
//    the form (AB, BC) where both edges are from the same S2ShapeIndex.
//    (This is an optimization for checking polygon validity.)
enum class CrossingType { INTERIOR, ALL, NON_ADJACENT };

// ShapeEdgeId is a unique identifier for an edge within an S2ShapeIndex,
// consisting of a (shape_id, edge_id) pair.  It is similar to
// std::pair<int32, int32> except that it has named fields.
// It should be passed and returned by value.
struct ShapeEdgeId {
 public:
  int32 shape_id, edge_id;
  ShapeEdgeId(int32 _shape_id, int32 _edge_id);
  bool operator==(ShapeEdgeId other) const;
  bool operator!=(ShapeEdgeId other) const;
  bool operator<(ShapeEdgeId other) const;
  bool operator>(ShapeEdgeId other) const;
  bool operator<=(ShapeEdgeId other) const;
  bool operator>=(ShapeEdgeId other) const;
};
std::ostream& operator<<(std::ostream& os, ShapeEdgeId id);

// A class representing a ShapeEdgeId together with the two endpoints of that
// edge.  It should be passed by reference.
struct ShapeEdge {
 public:
  ShapeEdge() : id_(-1, -1) {}
  ShapeEdge(S2Shape const& shape, int32 edge_id);
  ShapeEdgeId id() const { return id_; }
  S2Point const& v0() const { return edge_.v0; }
  S2Point const& v1() const { return edge_.v1; }

 private:
  ShapeEdgeId id_;
  S2Shape::Edge edge_;
};

// A function that is called with pairs of crossing edges.  The function may
// return false in order to request that the algorithm should be terminated,
// i.e. no further crossings are needed.
//
// "is_interior" indicates that the crossing is at a point interior to both
// edges (i.e., not at a vertex).  (The calling function already has this
// information and it is moderately expensive to recompute.)
using EdgePairVisitor = std::function<
  bool (ShapeEdge const& a, ShapeEdge const& b, bool is_interior)>;

// Visits all pairs of crossing edges in the given S2ShapeIndex, terminating
// early if the given EdgePairVisitor function returns false (in which case
// VisitCrossings returns false as well).  "type" indicates whether all
// crossings should be visited, or only interior crossings.
//
// CAVEAT: Crossings may be visited more than once.
bool VisitCrossings(S2ShapeIndex const& index, CrossingType type,
                    EdgePairVisitor const& visitor);

// Like the above, but visits all pairs of crossing edges where one edge comes
// from each S2ShapeIndex.
//
// REQUIRES: type != CrossingType::NON_ADJACENT (not supported)
bool VisitCrossings(S2ShapeIndex const& a, S2ShapeIndex const& b,
                    CrossingType type, EdgePairVisitor const& visitor);

// This is a helper function for implementing S2Shape::contains_origin().
//
// Given a shape consisting of closed polygonal loops, define the interior of
// the shape to be the region to the left of all edges (which must be oriented
// consistently).  This function then returns true if S2::Origin() is
// contained by the shape.
//
// Unlike S2Loop and S2Polygon, this method allows duplicate vertices and
// edges, which requires some extra care with definitions.  The rule that we
// apply is that an edge and its reverse edge "cancel" each other: the result
// is the same as if that edge pair were not present.  Therefore shapes that
// consist only of degenerate loop(s) are either empty or full; by convention,
// the shape is considered full if and only if it contains an empty loop (see
// LaxPolygon for details).
//
// Determining whether a loop on the sphere contains a point is harder than
// the corresponding problem in 2D plane geometry.  It cannot be implemented
// just by counting edge crossings because there is no such thing as a "point
// at infinity" that is guaranteed to be outside the loop.
bool IsOriginOnLeft(S2Shape const& shape);


/////////////////////////////////////////////////////////////////////////////
///////////////// Methods used internally by the S2 library /////////////////


// The purpose of this function is to construct polygons consisting of
// multiple loops.  It takes as input a collection of loops whose boundaries
// do not cross, and groups them into polygons whose interiors do not
// intersect (where the boundary of each polygon may consist of multiple
// loops).
//
// some of those islands have lakes, then the input to this function would
// islands, and their lakes.  Each loop would actually be present twice, once
// in each direction (see below).  The output would consist of one polygon
// representing each lake, one polygon representing each island not including
// islands or their lakes, and one polygon representing the rest of the world
//
// This method is intended for internal use; external clients should use
// S2Builder, which has more convenient interface.
//
// The input consists of a set of connected components, where each component
// consists of one or more loops that satisfy the following properties:
//
//  - The loops in each component must form a subdivision of the sphere (i.e.,
//    they must cover the entire sphere without overlap), except that a
//    component may consist of a single loop if and only if that loop is
//    degenerate (i.e., its interior is empty).
//
//  - The boundaries of different connected components must be disjoint
//    (i.e. no crossing edges or shared vertices).
//
//  - No component should be empty, and no loop should have zero edges.
//
// The output consists of a set of polygons, where each polygon is defined by
// the collection of loops that form its boundary.  This function does not
// actually construct any S2Shapes; it simply identifies the loops that belong
// to each polygon.
void ResolveComponents(
    std::vector<std::vector<S2Shape*>> const& components,
    std::vector<std::vector<S2Shape*>>* polygons);

// Given an S2ShapeIndex containing a single polygonal shape (e.g., an
// S2Polygon or S2Loop), return true if any loop has a self-intersection
// (including duplicate vertices) or crosses any other loop (including vertex
// crossings and duplicate edges) and set "error" to a human-readable error
// message.  Otherwise return false and leave "error" unchanged.
//
// This method is used to implement the FindValidationError methods of S2Loop
// and S2Polygon.
//
// TODO(ericv): Add an option to support LaxPolygon rules (i.e., duplicate
// vertices and edges are allowed, but loop crossings are not).
bool FindAnyCrossing(S2ShapeIndex const& index, S2Error* error);


//////////////////   Implementation details follow   ////////////////////


inline ShapeEdgeId::ShapeEdgeId(int32 _shape_id, int32 _edge_id)
    : shape_id(_shape_id), edge_id(_edge_id) {
}

inline bool ShapeEdgeId::operator==(ShapeEdgeId other) const {
  return shape_id == other.shape_id && edge_id == other.edge_id;
}

inline bool ShapeEdgeId::operator!=(ShapeEdgeId other) const {
  return !(*this == other);
}

inline bool ShapeEdgeId::operator<(ShapeEdgeId other) const {
  if (shape_id < other.shape_id) return true;
  if (other.shape_id < shape_id) return false;
  return edge_id < other.edge_id;
}

inline bool ShapeEdgeId::operator>(ShapeEdgeId other) const {
  return other < *this;
}

inline bool ShapeEdgeId::operator<=(ShapeEdgeId other) const {
  return !(other < *this);
}

inline bool ShapeEdgeId::operator>=(ShapeEdgeId other) const {
  return !(*this < other);
}

inline ShapeEdge::ShapeEdge(S2Shape const& shape, int32 edge_id)
    : id_(shape.id(), edge_id), edge_(shape.edge(edge_id)) {
}

}  // namespace s2shapeutil

#endif  // S2_S2SHAPEUTIL_H_
