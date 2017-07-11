// Copyright 2017 Google Inc. All Rights Reserved.
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

#ifndef S2_S2BOUNDARY_OPERATION_H_
#define S2_S2BOUNDARY_OPERATION_H_

#include <memory>
#include <utility>
#include <vector>
#include "s2/util/btree/btree_map.h"
#include "s2/s2builder.h"
#include "s2/s2builder_graph.h"
#include "s2/value_lexicon.h"

// This class implements boolean operations (intersection, union, difference,
// and symmetric difference) for regions whose boundaries are defined by
// geodesic edges.
//
// S2BoundaryOperation operates on exactly two input regions at a time.  Each
// region is represented as an S2ShapeIndex and may contain any number of
// points, polylines, and polygons.  The region is essentially the union of
// these objects, except that polygon interiors must be disjoint from all
// other geometry (including other polygon interiors).  If the input geometry
// for a region does not meet this condition, it can be normalized by
// computing its union first.  Note that points or polylines are allowed to
// coincide with the boundaries of polygons.
//
// Degeneracies are supported.  A polygon loop or polyline may consist of a
// single edge from a vertex to itself, and polygons may contain "sibling
// pairs" consisting of an edge and its corresponding reverse edge.  Polygons
// must not have any duplicate edges (due to the requirement that polygon
// interiors are disjoint), but polylines may have duplicate edges or can even
// be self-intersecting.
//
// Points and polyline edges are treated as multisets: if the same point or
// polyline edge appears multiple times in the input, it will appear multiple
// times in the output.  This is useful for modeling large sets of points or
// polylines as a single region while maintaining their distinct identities,
// even when the points or polylines intersect each other.  It is also useful
// for reconstructing polylines that loop back on themselves.  (Duplicate
// edges can be merged using GraphOptions::DuplicateEdges::MERGE in the
// S2Builder output layer, if desired.)
//
// Polylines are always considered to be directed.  Polyline edges between the
// same pair of vertices are defined to intersect even if the two edges are in
// opposite directions.  (Undirected polylines can be modeled by specifying
// GraphOptions::EdgeType::UNDIRECTED in the S2Builder output layer.)
//
// The output of each operation is sent to an S2Builder::Layer provided by the
// client.  This allows clients to build any representation of the geometry
// they choose.  It also allows the client to do additional postprocessing of
// the output before building data structures; for example, the client can
// easily discard degeneracies or convert them to another data type.
//
// The boundaries of polygons and polylines can be modeled as open, semi-open,
// or closed.  For polylines these options are defined as follows:
//
//  - In the OPEN model, polylines do not contain their first or last vertex.
//
//  - In the SEMI_OPEN model, polylines contain vertices except the last.
//    Therefore if one polyline starts where another polyline stops, the two
//    polylines do not intersect.
//
//  - In the CLOSED model, polylines contain all of their vertices.
//
// When multiple polylines are present, they are processed independently and
// have no effect on each other.  For example, in the OPEN boundary model the
// polyline ABC contains the vertex B, while set of polylines {AB, BC} does
// not.  (If you want to treat the polylines as a union instead, with
// boundaries merged according to the "mod 2" rule, this can be achieved by
// reassembling the edges into maximal polylines using S2PolylineVectorLayer
// with EdgeType::UNDIRECTED, DuplicateEdges::MERGE, and PolylineType::WALK.)
//
// For polygons:
//
//  - In the OPEN model, polygons do not contain their vertices or edges.
//    This implies that a polyline that follows the boundary of a polygon will
//    not intersect it.
//
//  - In the SEMI_OPEN model, polygon point containment is defined such that
//    if several polygons tile the region around a vertex, then exactly one of
//    those polygons contains that vertex.  Similarly polygons contain all of
//    their edges, but none of their reversed edges.  This implies that a
//    polyline and polygon edge with the same endpoints intersect if and only
//    if they are in the same direction.  (This rule ensures that if a
//    polyline is intersected with a polygon and its complement, the two
//    resulting polylines do not have any edges in common.)
//
//  - In the CLOSED model, polygons contain all their vertices, edges, and
//    reversed edges.  This implies that a polyline that shares an edge (in
//    either direction) with a polygon is defined to intersect it.  Similarly,
//    this is the only model where polygons that touch at a vertex or along an
//    edge intersect.
//
// Operations between geometry of different dimensions are defined as follows:
//
//  - For UNION, the higher-dimensional primitive always wins.  For example
//    the union of a closed polygon A with a polyline B that coincides with
//    the boundary of A consists only of the polygon A.
//
//  - For INTERSECTION, the lower-dimensional primitive always wins.  For
//    example, the intersection of a closed polygon A with a point B that
//    coincides with a vertex of A consists only of the point B.
//
//  - For DIFFERENCE, higher-dimensional primitives are not affected by
//    subtracting lower-dimensional subsets.  For example, subtracting a point
//    or polyline from a polygon A yields the original polygon A.  This rule
//    exists because in general, it is impossible to represent the output
//    using the specified boundary model(s).  (Consider subtracting one vertex
//    from a PolylineModel::CLOSED polyline, or subtracting one edge from a
//    PolygonModel::CLOSED polygon.)  If you want to perform operations like
//    this, consider representing all boundaries explicitly (topological
//    boundaries) using OPEN boundary models.  Another option for polygons is
//    to subtract a degenerate loop, which yields a polygon with a degenerate
//    hole (see s2shapeutil::LaxPolygon).
//
// Note that in the case of Precision::EXACT operations, the above remarks
// only apply to the output before snapping.  Snapping may cause nearby
// distinct edges to become coincident, e.g. a polyline may become coincident
// with a polygon boundary.  However also note that S2BoundaryOperation is
// perfectly happy to accept such geometry as input.
//
// Note the following differences between S2BoundaryOperation and the similar
// S2MultiBoundaryOperation class:
//
//  - S2BoundaryOperation operates on exactly two regions at a time, whereas
//    S2MultiBoundaryOperation operates on any number of regions.
//
//  - S2BoundaryOperation is potentially much faster when the input is already
//    represented as S2ShapeIndexes.  The algorithm is output sensitive and is
//    often sublinear in the input size.  This can be a big advantage if, say,
//
//  - S2BoundaryOperation supports exact predicates and the corresponding
//    exact operations (i.e., operations that are equivalent to computing the
//    exact result and then snap rounding it).
//
//  - S2MultiBoundaryOperation has better error guarantees when there are many
//    regions, since it requires only one snapping operation for any number of
//    input regions.
//
// Example usage:
//   S2ShapeIndex a, b;  // Input geometry, e.g. containing polygons.
//   S2Polygon polygon;  // Output geometry.
//   S2BoundaryOperation::Options options;
//   options.set_snap_function(snap_function);
//   S2BoundaryOperation op(S2BoundaryOperation::OpType::INTERSECTION,
//                          absl::MakeUnique<S2PolygonLayer>(&polygon),
//                          options);
//   S2Error error;
//   if (!op.Build(a, b, &error)) {
//     LOG(ERROR) << error.text();
//     ...
//   }
//
// If the output includes objects of different dimensions, they can be
// assembled into different layers with code like this:
//
//   vector<S2Point> points;
//   vector<unique_ptr<S2Polyline>> polylines;
//   S2Polygon polygon;
//   S2BoundaryOperation op(
//       S2BoundaryOperation::OpType::UNION,
//       absl::MakeUnique<s2builderutil::PointVectorLayer>(&points),
//       absl::MakeUnique<s2builderutil::S2PolylineVectorLayer>(&polylines),
//       absl::MakeUnique<S2PolygonLayer>(&polygon));

class S2BoundaryOperation {
 public:
  // The supported operation types.
  enum class OpType {
    UNION,                // Contained by either region.
    INTERSECTION,         // Contained by both regions.
    DIFFERENCE,           // Contained by the first region but not the second.
    SYMMETRIC_DIFFERENCE  // Contained by one region but not the other.
  };
  // Translates OpType to one of the strings above.
  static char const* OpTypeToString(OpType op_type);

  // Defines whether polygons are considered to contain their vertices and/or
  // edges (see definitions above).
  enum class PolygonModel { OPEN, SEMI_OPEN, CLOSED };

  // Defines whether polylines are considered to contain their endpoints
  // (see definitions above).
  enum class PolylineModel { OPEN, SEMI_OPEN, CLOSED };

  // With Precision::EXACT, the boundary operation is evaluated using the
  // exact input geometry.  Predicates that use this option will produce exact
  // results; for example, they can distinguish between a polyline that barely
  // intersects a polygon from one that barely misses it.  Constructive
  // operations (ones that yield new geometry, as opposed to predicates) are
  // implemented by computing the exact result and then snap rounding it
  // according to the given snap_function() (see below).  This is as close as
  // it is possible to get to the exact result while requiring that vertex
  // coordinates have type "double".
  //
  // With Precision::SNAPPED, the input regions are snapped together *before*
  // the boundary operation is evaluated.  So for example, two polygons that
  // overlap slightly will be treated as though they share a common boundary,
  // and similarly two polygons that are slightly separated from each other
  // will be treated as though they share a common boundary.  Snapped results
  // are useful for dealing with points, since in S2 the only points that lie
  // exactly on a polyline or polygon edge are the endpoints of that edge.
  //
  // Conceptually, the difference between these two options is that with
  // Precision::SNAPPED, the inputs are snap rounded (together), whereas with
  // Precision::EXACT only the result is snap rounded.
  enum class Precision { EXACT, SNAPPED };

  // Forward declaration for Options::set_source_id_lexicon().  (This option
  // lets you determine the mapping from input edges to output edges.)
  class SourceId;

  class Options {
   public:
    Options();

    // Convenience constructor that calls set_snap_function().
    explicit Options(S2Builder::SnapFunction const& snap_function);

    // Specifies the function to be used for snap rounding.
    //
    // DEFAULT: s2builderutil::IdentitySnapFunction(S1Angle::Zero())
    // [This does no snapping and preserves all input vertices exactly unless
    //  there are crossing edges, in which case the snap radius is increased
    //  to the maximum intersection point error (S2::kIntersectionError.]
    S2Builder::SnapFunction const& snap_function() const;
    void set_snap_function(S2Builder::SnapFunction const& snap_function);

    // Defines whether polygons are considered to contain their vertices
    // and/or edges.
    //
    // DEFAULT: PolygonModel::SEMI_OPEN
    PolygonModel polygon_model() const;
    void set_polygon_model(PolygonModel model);

    // Defines whether polylines are considered to contain their vertices.
    //
    // DEFAULT: PolylineModel::CLOSED
    PolylineModel polyline_model() const;
    void set_polyline_model(PolylineModel model);

    // Specifies whether the operation should use the exact input geometry
    // (Precision::EXACT), or whether the two input regions should be snapped
    // together first (Precision::SNAPPED).
    //
    // DEFAULT: Precision::EXACT
    Precision precision() const;
    // void set_precision(Precision precision);

    // If true, the input geometry is interpreted as representing nearby
    // geometry that has been snapped or simplified.  It then outputs a
    // conservative result based on the value of polygon_model() and
    // polyline_model().  For the most part, this only affects the handling of
    // degeneracies.
    //
    // - If the model is OPEN, the result is as open as possible.  For
    //   example, the intersection of two identical degenerate shells is empty
    //   under PolygonModel::OPEN because they could have been disjoint before
    //   snapping.  Similarly, two identical degenerate polylines have an
    //   empty intersection under PolylineModel::OPEN.
    //
    // - If the model is CLOSED, the result is as closed as possible.  In the
    //   case of the DIFFERENCE operation, this is equivalent to evaluating
    //   A - B as Closure(A) - Interior(B).  For other operations, it affects
    //   only the handling of degeneracies.  For example, the union of two
    //   identical degenerate holes is empty under PolygonModel::CLOSED
    //   (i.e., the hole disappears) because the holes could have been
    //   disjoint before snapping.
    //
    // - If the model is SEMI_OPEN, the result is as degenerate as possible.
    //   New degeneracies will not be created, but all degeneracies that
    //   coincide with the opposite region's boundary are retained unless this
    //   would cause a duplicate polygon edge to be created.  This model is
    //   is very useful for working with input data that has both positive and
    //   negative degeneracies (i.e., degenerate shells and holes).
    //
    // DEFAULT: false
    bool conservative_output() const;
    // void set_conservative_output(bool conservative);

    // If specified, then each output edge will be labelled with one or more
    // SourceIds indicating which input edge(s) it corresponds to.  This
    // can be useful if your input geometry has additional data that needs to
    // be propagated from the input to the output (e.g., elevations).
    //
    // You can access the labels by using an S2Builder::Layer type that
    // supports labels, such as S2PolygonLayer.  The layer outputs a
    // "label_set_lexicon" and an "label_set_id" for each edge.  You can then
    // look up the source information for each edge like this:
    //
    // for (int32 label : label_set_lexicon.id_set(label_set_id)) {
    //   SourceId const& src = source_id_lexicon.value(label);
    //   // region_id() specifies which S2ShapeIndex the edge is from (0 or 1).
    //   DoSomething(src.region_id(), src.shape_id(), src.edge_id());
    // }
    //
    // DEFAULT: nullptr
    ValueLexicon<SourceId>* source_id_lexicon() const;
    // void set_source_id_lexicon(ValueLexicon<SourceId>* source_id_lexicon);

    // Options may be assigned and copied.
    Options(Options const& options);
    Options& operator=(Options const& options);

   private:
    std::unique_ptr<S2Builder::SnapFunction> snap_function_;
    PolygonModel polygon_model_;
    PolylineModel polyline_model_;
    Precision precision_;
    bool conservative_;
    ValueLexicon<SourceId>* source_id_lexicon_;
  };

  S2BoundaryOperation(OpType op_type,
                      std::unique_ptr<S2Builder::Layer> layer,
                      Options const& options = Options());

  // Specifies that "result_non_empty" should be set to indicate whether the
  // exact result of the operation is non-empty (contains at least one edge).
  // This can be used to efficiently implement boolean relationship
  // predicates.  For example, you can check whether a polygon contains a
  // polyline by subtracting the polygon from the polyline (BuildDifference)
  // and testing whether the result is empty.
  S2BoundaryOperation(OpType op_type, bool* result_non_empty,
                      Options const& options = Options());

  // Specifies that the output boundary edges should be sent to three
  // different layers according to their dimension.  Points (represented by
  // degenerate edges) are sent to layer 0, polyline edges are sent to
  // layer 1, and polygon edges are sent to layer 2.
  //
  // The dimension of an edge is defined as the minimum dimension of the two
  // input edges that produced it.  For example, the intersection of two
  // crossing polyline edges is a considered to be a degenerate polyline
  // rather than a point, so it is sent to layer 1.  Clients can easily
  // reclassify such polylines as points if desired, but this rule makes it
  // easier for clients that want to process point, polyline, and polygon
  // inputs differently.
  //
  // The layers are always built in the order 0, 1, 2, and all arguments to
  // the Build() calls are guaranteed to be valid until the last call returns.
  // All Graph objects have the same set of vertices and the same lexicon
  // objects, in order to make it easier to write classes that process all the
  // edges in parallel.
  S2BoundaryOperation(OpType op_type,
                      std::vector<std::unique_ptr<S2Builder::Layer>> layers,
                      Options const& options = Options());

  OpType op_type() const { return op_type_; }

  // Executes the given operation.  Returns true on success, and otherwise
  // sets "error" appropriately.  (This class does not generate any errors
  // itself, but the S2Builder::Layer might.)
  bool Build(S2ShapeIndex const& a, S2ShapeIndex const& b, S2Error* error);

  // SourceId identifies an edge from one of the two input S2ShapeIndexes.
  // It consists of a region id (0 or 1), a shape id within that region's
  // S2ShapeIndex, and an edge id within that shape.
  class SourceId {
   public:
    SourceId();
    SourceId(int region_id, int32 shape_id, int32 edge_id);
    explicit SourceId(int32 special_edge_id);
    int region_id() const { return region_id_; }
    int32 shape_id() const { return shape_id_; }
    int32 edge_id() const { return edge_id_; }
    // TODO(ericv): Convert to functions, define all 6 comparisons.
    bool operator==(SourceId other) const;
    bool operator<(SourceId other) const;

   private:
    uint32 region_id_ : 1;
    uint32 shape_id_ : 31;
    int32 edge_id_;
  };

 private:
  class Impl;  // The actual implementation.

  // Internal constructor to reduce code duplication.
  S2BoundaryOperation(OpType op_type, Options const& options);

  OpType op_type_;
  Options options_;

  // The input regions.
  S2ShapeIndex const* regions_[2];

  // The output consists either of zero layers, one layer, or three layers.
  std::vector<std::unique_ptr<S2Builder::Layer>> layers_;

  // The following field is set if and only if there are no output layers.
  bool* result_non_empty_;
};


//////////////////   Implementation details follow   ////////////////////


inline S2BoundaryOperation::SourceId::SourceId()
    : region_id_(0), shape_id_(0), edge_id_(-1) {
}

inline S2BoundaryOperation::SourceId::SourceId(
    int region_id, int32 shape_id, int32 edge_id)
    : region_id_(region_id), shape_id_(shape_id), edge_id_(edge_id) {
}

inline S2BoundaryOperation::SourceId::SourceId(int special_edge_id)
    : region_id_(0), shape_id_(0), edge_id_(special_edge_id) {
}

inline bool S2BoundaryOperation::SourceId::operator==(SourceId other) const {
  return (region_id_ == other.region_id_ &&
          shape_id_ == other.shape_id_ &&
          edge_id_ == other.edge_id_);
}

inline bool S2BoundaryOperation::SourceId::operator<(SourceId other) const {
  if (region_id_ < other.region_id_) return true;
  if (region_id_ > other.region_id_) return false;
  if (shape_id_ < other.shape_id_) return true;
  if (shape_id_ > other.shape_id_) return false;
  return edge_id_ < other.edge_id_;
}

#endif  // S2_S2BOUNDARY_OPERATION_H_
