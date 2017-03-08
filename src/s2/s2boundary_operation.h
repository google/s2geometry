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
#include "s2/util/btree/btree_map.h"  // Like std::map, but faster and smaller.
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
// computing its union first.
//
// Degeneracies are supported.  A polygon loop or polyline may consist of a
// single edge from a vertex to itself, and polygons may contain "sibling
// pairs" consisting of an edge and its corresponding reverse edge.  Polylines
// may also be self-intersecting or contain duplicate edges, however each
// polygon edge may occur only once in a given input region.
//
// Points and polyline edges are treated as multisets: if the same point or
// polyline edge appears multiple times in the input, it will appear multiple
// times in the output.  This is useful for reconstructing polylines that loop
// back on themselves, or distinct polylines that have been snapped together.
// (Duplicate edges can be merged using GraphOptions::DuplicateEdges::MERGE in
// the S2Builder output layer, if desired.)
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
//  - In the SEMI_OPEN model, polylines contain vertices except the last.
//    Therefore if one polyline starts where another polyline stops, the two
//    polylines do not intersect.
//  - In the CLOSED model, polylines contain all of their vertices.
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
// CAVEATS:
//
//  - Currently this class requires that every S2Shape must consist of a
//    single chain of edges, i.e. it must be a closed loop, a single polyline,
//    or a single point (represented as a degenerate edge).
//    TODO(ericv): Consider relaxing this requirement.
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
  // intersects a polygon from one that barely misses it.  The corresponding
  // operations (yielding a new geometric object) are implemented by computing
  // the exact result and then snap rounding it.  (This is as close as it is
  // possible to get to the exact result while requiring that vertex
  // coordinates have type "double".)
  //
  // With Precision::SNAPPED, the input regions are snapped together before
  // the boundary operation is evaluated.  So for example, two polygons that
  // overlap slightly will be treated as though they share a common boundary,
  // and similarly two polygons that are slightly separated from each other
  // will be treated as though they share a common boundary.  Snapped results
  // are useful for dealing with points, since in S2 the only points that lie
  // exactly on a polyline or polygon edge are the endpoints of that edge.
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
    // Default value: s2builderutil::IdentitySnapFunction(S1Angle::Zero()).
    // This does no snapping and preserves all input vertices exactly unless
    // there are crossing edges, in which case the snap radius is increased to
    // the maximum intersection point error (S2EdgeUtil::kIntersectionError).
    S2Builder::SnapFunction const& snap_function() const;
    void set_snap_function(S2Builder::SnapFunction const& snap_function);

    // Defines whether polygons are considered to contain their vertices
    // and/or edges.
    //
    // Default value: SEMI_OPEN.
    PolygonModel polygon_model() const;
    // void set_polygon_model(PolygonModel model);

    // Defines whether polylines are considered to contain their vertices.
    //
    // Default value: CLOSED.
    PolylineModel polyline_model() const;
    // void set_polyline_model(PolylineModel model);

    // Specifies whether the operation should use the exact input geometry
    // (Precision::EXACT), or whether the two input regions should be snapped
    // together first (Precision::SNAPPED).
    //
    // Default value: EXACT.
    Precision precision() const;
    // void set_precision(Precision precision);

    // If true, the input geometry is interpreted as representing nearby
    // geometry that has been snapped or simplified.  It then outputs a
    // conservative result based on the value of polygon_model():
    //
    // - If polygon_model() is OPEN, the result is as open as possible.  For
    //   example, points never intersect polylines or other inputs points
    //   under this model since if the points are moved slightly then the
    //   intersection vanishes.  Similarly, two polylines intersect only if
    //   one polyline properly crosses the other, since if the two polylines
    //   touch without crossing then the intersection can be removed by
    //   perturbing the input geometry slightly.
    //
    // - In the CLOSED model, the result is as closed as possible.  For
    //   example subtracting point or polylines from other points or polylines
    //   will not have any effect, since the intersection can be removed by
    //   perturbing the geometry slightly.  A more surprising result is that
    //   subtracting a polygon from itself yields the boundary of that
    //   polygon, because the first polygon may have properly contained the
    //   second polygon before snapping.
    //
    // - In the SEMI_OPEN model, the result is as degenerate as possible.  In
    //   particular, degeneracies that coincide with the opposite region's
    //   boundary are retained unless this would cause a duplicate edge to be
    //   emitted.
    //
    // Default value: false.
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
    // Default value: nullptr.
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

#if 0
  // Specifies that "result_empty" should be set to indicate whether the exact
  // result of the operation is empty (contains no edges).  This can be used
  // to efficiently implement boolean relationship predicates.  For example,
  // you can check whether a polygon contains a polyline by subtracting the
  // polygon from the polyline (BuildDifference) and testing whether the
  // result is empty.
  S2BoundaryOperation(OpType op_type, bool* result_empty,
                      Options const& options = Options());

  // Specifies that the output boundary edges should be sent to three
  // different layers according to their dimension.  Points (represented by
  // degenerate edges) are sent to "layer0", polyline edges are sent to
  // "layer1", and polygon edges are sent to "layer2".
  //
  // The dimension of an edge is defined as the minimum dimension of the two
  // input edges that produced it.  For example, the intersection of two
  // crossing polyline edges is a considered to be a degenerate polyline
  // rather than a point, so it is sent to "layer1".  Clients can easily
  // reclassify such polylines as points if desired, but this rule makes it
  // easier for clients that want to process point, polyline, and polygon
  // inputs differently.
  S2BoundaryOperation(OpType op_type,
                      std::unique_ptr<S2Builder::Layer> layer0,
                      std::unique_ptr<S2Builder::Layer> layer1,
                      std::unique_ptr<S2Builder::Layer> layer2,
                      Options const& options = Options());
#endif

  // Adds an input region consisting of the geometry in the given
  // S2ShapeIndex.  There must be exactly two calls to this method, since this
  // class only supports operations on pairs of regions.
  //
  // The index may contain any mixture of 0D, 1D, and 2D geometry (i.e.,
  // points, polyline, and polygons).  The rules for combining geometric
  // primitives in each dimension are documented above.
  void AddRegion(S2ShapeIndex const& index);

  // Executes the given operation.  Returns true on success, and otherwise
  // sets "error" appropriately.  (This class does not generate any errors
  // itself, but the S2Builder::Layer might.)
  bool Build(S2Error* error);

  // SourceId identifies an edge from one of the two input S2ShapeIndexes.
  // It consists of a region id (0 or 1), a shape id within that region's
  // S2ShapeIndex, and an edge id within that shape.
  class SourceId {
   public:
    SourceId();
    SourceId(int region_id, int32 shape_id, int32 edge_id);
    int region_id() const { return region_id_; }
    int32 shape_id() const { return shape_id_; }
    int32 edge_id() const { return edge_id_; }
    bool operator==(SourceId other) const;
    bool operator<(SourceId other) const;

   private:
    uint32 region_id_ : 1;
    uint32 shape_id_ : 31;
    int32 edge_id_;
  };

  // This version of AddRegion() is used internally by S2Polygon to adapt its
  // internal S2ShapeIndex to the requirements of this class.
  //
  // The function "contains" indicates whether the given point is contained by
  // the region.  The result must depend only on the shapes in the index where
  // has_interior() is true (i.e., polygons).  Furthermore, the result must
  // change whenever an edge from one of these shapes is crossed.  This
  // implies that for a given set of indexed shapes, only two (complementary)
  // definitions of the interior are possible.
  //
  // The function "reverse_edges" indicates whether the edges of the given
  // shape should be reversed before using them.  After any such reversal, it
  // must be true that the interior of the region is to the left of all edges
  // that define the interior (i.e., shapes where has_interior() is true).
  // (Shapes that define polylines or point sets have no such restrictions.)
  void AddRegion(S2ShapeIndex const& index,
                 std::function<bool (S2Point const&)> contains,
                 std::function<bool (S2Shape const&)> reverse_edges);

 private:
  class Impl;  // The actual implementation.
  class CrossingQuery;  // Helper class.

  // A Region consists of an S2ShapeIndex together with the functions that
  // define point containment and whether any shape edges are reversed.
  struct Region {
    int id;  // The input region id (0 or 1).
    S2ShapeIndex const& index;
    std::function<bool (S2Point const&)> contains;
    std::function<bool (S2Shape const&)> reverse_edges;

    Region(int _id, S2ShapeIndex const& _index,
           std::function<bool (S2Point const&)> _contains,
           std::function<bool (S2Shape const&)> _reverse_edges)
        : id(_id), index(_index),
          contains(std::move(_contains)),
          reverse_edges(std::move(_reverse_edges)) {
    }
  };

  // Internal constructor to reduce code duplication.
  S2BoundaryOperation(OpType op_type, Options const& options);

  OpType op_type_;
  Options options_;

  // The input regions (there are always exactly two when Build is called).
  std::vector<Region> regions_;

  // The output consists either of zero layers, one layer, or three layers.
  std::vector<std::unique_ptr<S2Builder::Layer>> layers_;

  // The following field is set if and only if there are no output layers.
  bool* result_nonempty_;
};


//////////////////   Implementation details follow   ////////////////////


inline S2BoundaryOperation::SourceId::SourceId()
    : region_id_(0), shape_id_(0), edge_id_(-1) {
}

inline S2BoundaryOperation::SourceId::SourceId(
    int region_id, int32 shape_id, int32 edge_id)
    : region_id_(region_id), shape_id_(shape_id), edge_id_(edge_id) {
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
