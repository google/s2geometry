// Copyright 2016 Google Inc. All Rights Reserved.
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
// This class is a replacement for S2PolygonBuilder.  Once all clients have
// been updated to use this class, S2PolygonBuilder will be removed.

#ifndef S2_S2BUILDER_H_
#define S2_S2BUILDER_H_

#include <algorithm>
#include <cstdint>
#include <functional>
#include <memory>
#include <utility>
#include <vector>

#include "absl/base/macros.h"
#include "absl/container/flat_hash_set.h"
#include "absl/log/absl_check.h"
#include "absl/log/absl_log.h"
#include "absl/types/span.h"

#include "s2/_fp_contract_off.h"  // IWYU pragma: keep
#include "s2/id_set_lexicon.h"
#include "s2/mutable_s2shape_index.h"
#include "s2/s1angle.h"
#include "s2/s1chord_angle.h"
#include "s2/s2cell_id.h"
#include "s2/s2edge_crossings.h"
#include "s2/s2edge_distances.h"
#include "s2/s2error.h"
#include "s2/s2memory_tracker.h"
#include "s2/s2point.h"
#include "s2/s2point_index.h"
#include "s2/s2point_span.h"
#include "s2/s2shape.h"
#include "s2/s2shape_index.h"
#include "s2/util/gtl/compact_array.h"

class S2Loop;
class S2Polygon;
class S2Polyline;

// S2Builder is a tool for assembling polygonal geometry from edges.  Here are
// some of the things it is designed for:
//
// 1. Building polygons, polylines, and polygon meshes from unsorted
//    collections of edges.
//
// 2. Snapping geometry to discrete representations (such as S2CellId centers
//    or E7 lat/lng coordinates) while preserving the input topology and with
//    guaranteed error bounds.
//
// 3. Simplifying geometry (e.g. for indexing, display, or storage).
//
// 4. Importing geometry from other formats, including repairing geometry
//    that has errors.
//
// 5. As a tool for implementing more complex operations such as polygon
//    intersections and unions.
//
// The implementation is based on the framework of "snap rounding".  Unlike
// most snap rounding implementations, S2Builder defines edges as geodesics on
// the sphere (straight lines) and uses the topology of the sphere (i.e.,
// there are no "seams" at the poles or 180th meridian).  The algorithm is
// designed to be 100% robust for arbitrary input geometry.  It offers the
// following properties:
//
//   - Guaranteed bounds on how far input vertices and edges can move during
//     the snapping process (i.e., at most the given "snap_radius").
//
//   - Guaranteed minimum separation between edges and vertices other than
//     their endpoints (similar to the goals of Iterated Snap Rounding).  In
//     other words, edges that do not intersect in the output are guaranteed
//     to have a minimum separation between them.
//
//   - Idempotency (similar to the goals of Stable Snap Rounding), i.e. if the
//     input already meets the output criteria then it will not be modified.
//
//   - Preservation of the input topology (up to the creation of
//     degeneracies).  This means that there exists a continuous deformation
//     from the input to the output such that no vertex crosses an edge.  In
//     other words, self-intersections won't be created, loops won't change
//     orientation, etc.
//
//   - The ability to snap to arbitrary discrete point sets (such as S2CellId
//     centers, E7 lat/lng points on the sphere, or simply a subset of the
//     input vertices), rather than being limited to an integer grid.
//
// Here are some of its other features:
//
//  - It can handle both directed and undirected edges.  Undirected edges can
//    be useful for importing data from other formats, e.g. where loops have
//    unspecified orientations.
//
//  - It can eliminate self-intersections by finding all edge pairs that cross
//    and adding a new vertex at each intersection point.
//
//  - It can simplify polygons to within a specified tolerance.  For example,
//    if two vertices are close enough they will be merged, and if an edge
//    passes nearby a vertex then it will be rerouted through that vertex.
//    Optionally, it can also detect nearly straight chains of short edges and
//    replace them with a single long edge, while maintaining the same
//    accuracy, separation, and topology guarantees ("simplify_edge_chains").
//
//  - It supports many different output types through the concept of "layers"
//    (polylines, polygons, polygon meshes, etc).  You can build multiple
//    layers at once in order to ensure that snapping does not create
//    intersections between different objects (for example, you can simplify a
//    set of contour lines without the risk of having them cross each other).
//
//  - It supports edge labels, which allow you to attach arbitrary information
//    to edges and have it preserved during the snapping process.  (This can
//    also be achieved using layers, at a coarser level of granularity.)
//
// Caveats:
//
//  - Because S2Builder only works with edges, it cannot distinguish between
//    the empty and full polygons.  If your application can generate both the
//    empty and full polygons, you must implement logic outside of this class.
//
// Example showing how to snap a polygon to E7 coordinates:
//
//  using s2builderutil::IntLatLngSnapFunction;
//  S2Builder builder(S2Builder::Options(IntLatLngSnapFunction(7)));
//  S2Polygon output;
//  builder.StartLayer(std::make_unique<s2builderutil::S2PolygonLayer>(&output));
//  builder.AddPolygon(input);
//  S2Error error;
//  if (!builder.Build(&error)) {
//    ABSL_LOG(ERROR) << error;
//    ...
//  }
class S2Builder {
 public:
  // Indicates whether the input edges are undirected.  Typically this is
  // specified for each output layer (e.g., s2builderutil::S2PolygonLayer).
  //
  // Directed edges are preferred, since otherwise the output is ambiguous.
  // For example, output polygons may be the *inverse* of the intended result
  // (e.g., a polygon intended to represent the world's oceans may instead
  // represent the world's land masses).  Directed edges are also somewhat
  // more efficient.
  //
  // However even with undirected edges, most S2Builder layer types try to
  // preserve the input edge direction whenever possible.  Generally, edges
  // are reversed only when it would yield a simpler output.  For example,
  // S2PolygonLayer assumes that polygons created from undirected edges should
  // cover at most half of the sphere.  Similarly, S2PolylineVectorLayer
  // assembles edges into as few polylines as possible, even if this means
  // reversing some of the "undirected" input edges.
  //
  // For shapes with interiors, directed edges should be oriented so that the
  // interior is to the left of all edges.  This means that for a polygon with
  // holes, the outer loops ("shells") should be directed counter-clockwise
  // while the inner loops ("holes") should be directed clockwise.  Note that
  // S2Builder::AddPolygon() follows this convention automatically.
  enum class EdgeType : uint8_t { DIRECTED, UNDIRECTED };

  // A SnapFunction restricts the locations of the output vertices.  For
  // example, there are predefined snap functions that require vertices to be
  // located at S2CellId centers or at E5/E6/E7 coordinates.  The SnapFunction
  // can also specify a minimum spacing between vertices (the "snap radius").
  //
  // A SnapFunction defines the following methods:
  //
  // 1. The SnapPoint() method, which snaps a point P to a nearby point (the
  //    "candidate snap site").  Any point may be returned, including P
  //    itself (this is the "identity snap function").
  //
  // 2. "snap_radius", the maximum distance that vertices can move when
  //    snapped.  The snap_radius must be at least as large as the maximum
  //    distance between P and SnapPoint(P) for any point P.
  //
  //    Note that the maximum distance that edge interiors can move when
  //    snapped is slightly larger than "snap_radius", and is returned by the
  //    function S2Builder::Options::max_edge_deviation() (see there for
  //    details).
  //
  // 3. "min_vertex_separation", the guaranteed minimum distance between
  //    vertices in the output.  This is generally a fraction of
  //    "snap_radius" where the fraction depends on the snap function.
  //
  // 4. "min_edge_vertex_separation", the guaranteed minimum distance between
  //    edges and non-incident vertices in the output.  This is generally a
  //    fraction of "snap_radius" where the fraction depends on the snap
  //    function.
  //
  // It is important to note that SnapPoint() does not define the actual
  // mapping from input vertices to output vertices, since the points it
  // returns (the candidate snap sites) are further filtered to ensure that
  // they are separated by at least the snap radius.  For example, if you
  // specify E7 coordinates (2cm resolution) and a snap radius of 10m, then a
  // subset of points returned by SnapPoint will be chosen (the "snap sites"),
  // and each input vertex will be mapped to the closest site.  Therefore you
  // cannot assume that P is necessarily snapped to SnapPoint(P).
  //
  // S2Builder makes the following guarantees:
  //
  // 1. Every vertex is at a location returned by SnapPoint().
  //
  // 2. Vertices are within "snap_radius" of the corresponding input vertex.
  //
  // 3. Edges are within "max_edge_deviation" of the corresponding input edge
  //    (a distance slightly larger than "snap_radius").
  //
  // 4. Vertices are separated by at least "min_vertex_separation"
  //    (a fraction of "snap_radius" that depends on the snap function).
  //
  // 5. Edges and non-incident vertices are separated by at least
  //    "min_edge_vertex_separation" (a fraction of "snap_radius").
  //
  // 6. Vertex and edge locations do not change unless one of the conditions
  //    above is not already met (idempotency / stability).
  //
  // 7. The topology of the input geometry is preserved (up to the creation
  //    of degeneracies).  This means that there exists a continuous
  //    deformation from the input to the output such that no vertex
  //    crosses an edge.
  class SnapFunction {
   public:
    virtual ~SnapFunction() = default;

    // The maximum distance that vertices can move when snapped.  The snap
    // radius can be any value between zero and SnapFunction::kMaxSnapRadius().
    //
    // If the snap radius is zero, then vertices are snapped together only if
    // they are identical.  Edges will not be snapped to any vertices other
    // than their endpoints, even if there are vertices whose distance to the
    // edge is zero, unless split_crossing_edges() is true (see below).
    //
    // REQUIRES: snap_radius() <= kMaxSnapRadius
    virtual S1Angle snap_radius() const = 0;

    // The maximum supported snap radius (equivalent to about 7800km).
    static S1Angle kMaxSnapRadius();

    // The guaranteed minimum distance between vertices in the output.
    // This is generally some fraction of "snap_radius".
    virtual S1Angle min_vertex_separation() const = 0;

    // The guaranteed minimum spacing between edges and non-incident vertices
    // in the output.  This is generally some fraction of "snap_radius".
    virtual S1Angle min_edge_vertex_separation() const = 0;

    // Returns a candidate snap site for the given point.  The final vertex
    // locations are a subset of the snap sites returned by this function
    // (spaced at least "min_vertex_separation" apart).
    //
    // The only requirement is that SnapPoint(x) must return a point whose
    // distance from "x" is no greater than "snap_radius".
    virtual S2Point SnapPoint(const S2Point& point) const = 0;

    // Returns a deep copy of this SnapFunction.
    virtual std::unique_ptr<SnapFunction> Clone() const = 0;
  };

  class Options {
   public:
    Options();

    // Convenience constructor that calls set_snap_function().
    explicit Options(const SnapFunction& snap_function);

    // Sets the desired snap function.  The snap function is copied
    // internally, so you can safely pass a temporary object.
    //
    // Note that if your input data includes vertices that were created using
    // S2::GetIntersection(), then you should use a "snap_radius" of
    // at least S2::kIntersectionMergeRadius, e.g. by calling
    //
    //  options.set_snap_function(s2builderutil::IdentitySnapFunction(
    //      S2::kIntersectionMergeRadius));
    //
    // DEFAULT: s2builderutil::IdentitySnapFunction(S1Angle::Zero())
    // [This does no snapping and preserves all input vertices exactly.]
    const SnapFunction& snap_function() const;
    void set_snap_function(const SnapFunction& snap_function);

    // The maximum distance from snapped edge vertices to the original edge.
    // This is the same as snap_function().snap_radius() except when
    // split_crossing_edges() is true (see below), in which case the edge snap
    // radius is increased by S2::kIntersectionError.
    S1Angle edge_snap_radius() const;

    // The maximum distance that any point along an edge can move when snapped.
    // It is slightly larger than edge_snap_radius() because when a geodesic
    // edge is snapped, the edge center moves further than its endpoints.
    // S2Builder ensures that this distance is at most 10% larger than
    // edge_snap_radius().
    S1Angle max_edge_deviation() const;

    // If true, then detect all pairs of crossing edges and eliminate them by
    // adding a new vertex at their intersection point.  See also the
    // AddIntersection() method which allows intersection points to be added
    // selectively.
    //
    // When this option is true, intersection_tolerance() is automatically set
    // to a minimum of S2::kIntersectionError (see intersection_tolerance()
    // for why this is necessary).  Note that this means that edges can move
    // by up to S2::kIntersectionError even when the specified snap radius is
    // zero.  The exact distance that edges can move is always given by
    // max_edge_deviation() defined above.
    //
    // Undirected edges should always be used when the output is a polygon,
    // since splitting a directed loop at a self-intersection converts it into
    // two loops that don't define a consistent interior according to the
    // "interior is on the left" rule.  (On the other hand, it is fine to use
    // directed edges when defining a polygon *mesh* because in that case the
    // input consists of sibling edge pairs.)
    //
    // Self-intersections can also arise when importing data from a 2D
    // projection.  You can minimize this problem by subdividing the input
    // edges so that the S2 edges (which are geodesics) stay close to the
    // original projected edges (which are curves on the sphere).  This can
    // be done using S2EdgeTessellator, for example.
    //
    // DEFAULT: false
    bool split_crossing_edges() const;
    void set_split_crossing_edges(bool split_crossing_edges);

    // Specifies the maximum allowable distance between a vertex added by
    // AddIntersection() and the edge(s) that it is intended to snap to.  This
    // method must be called before AddIntersection() can be used.  It has the
    // effect of increasing the snap radius for edges (but not vertices) by
    // the given distance.
    //
    // The intersection tolerance should be set to the maximum error in the
    // intersection calculation used.  For example, if S2::GetIntersection()
    // is used then the error should be set to S2::kIntersectionError.  If
    // S2::GetPointOnLine() is used then the error should be set to
    // S2::kGetPointOnLineError.  If S2::Project() is used then the error
    // should be set to S2::kProjectPerpendicularError.  If more than one
    // method is used then the intersection tolerance should be set to the
    // maximum such error.
    //
    // The reason this option is necessary is that computed intersection
    // points are not exact.  For example, S2::GetIntersection(a, b, c, d)
    // returns a point up to S2::kIntersectionError away from the true
    // mathematical intersection of the edges AB and CD.  Furthermore such
    // intersection points are subject to further snapping in order to ensure
    // that no pair of vertices is closer than the specified snap radius.  For
    // example, suppose the computed intersection point X of edges AB and CD
    // is 1 nanonmeter away from both edges, and the snap radius is 1 meter.
    // In that case X might snap to another vertex Y exactly 1 meter away,
    // which would leave us with a vertex Y that could be up to 1.000000001
    // meters from the edges AB and/or CD.  This means that AB and/or CD might
    // not snap to Y leaving us with two edges that still cross each other.
    //
    // However if the intersection tolerance is set to 1 nanometer then the
    // snap radius for edges is increased to 1.000000001 meters ensuring that
    // both edges snap to a common vertex even in this worst case.  (This
    // technique does not work if the vertex snap radius is increased as well;
    // it requires edges and vertices to be handled differently.)
    //
    // Note that this option allows edges to move by up to the given
    // intersection tolerance even when the snap radius is zero.  The exact
    // distance that edges can move is always given by max_edge_deviation()
    // defined above.
    //
    // When split_crossing_edges() is true, the intersection tolerance is
    // automatically set to a minimum of S2::kIntersectionError.  A larger
    // value can be specified by calling this method explicitly.
    //
    // DEFAULT: S1Angle::Zero()
    S1Angle intersection_tolerance() const;
    void set_intersection_tolerance(S1Angle intersection_tolerance);

    // If true, then simplify the output geometry by replacing nearly straight
    // chains of short edges with a single long edge.
    //
    // The combined effect of snapping and simplifying will not change the
    // input by more than the guaranteed tolerances (see the list documented
    // with the SnapFunction class).  For example, simplified edges are
    // guaranteed to pass within snap_radius() of the *original* positions of
    // all vertices that were removed from that edge.  This is a much tighter
    // guarantee than can be achieved by snapping and simplifying separately.
    //
    // However, note that this option does not guarantee idempotency.  In
    // other words, simplifying geometry that has already been simplified once
    // may simplify it further.  (This is unavoidable, since tolerances are
    // measured with respect to the original geometry, which is no longer
    // available when the geometry is simplified a second time.)
    //
    // When the output consists of multiple layers, simplification is
    // guaranteed to be consistent: for example, edge chains are simplified in
    // the same way across layers, and simplification preserves topological
    // relationships between layers (e.g., no crossing edges will be created).
    // Note that edge chains in different layers do not need to be identical
    // (or even have the same number of vertices, etc) in order to be
    // simplified together.  All that is required is that they are close
    // enough together so that the same simplified edge can meet all of their
    // individual snapping guarantees.
    //
    // Note that edge chains are approximated as parametric curves rather than
    // point sets.  This means that if an edge chain backtracks on itself (for
    // example, ABCDEFEDCDEFGH) then such backtracking will be preserved to
    // within snap_radius() (for example, if the preceding point were all in a
    // straight line then the edge chain would be simplified to ACFCFH, noting
    // that C and F have degree > 2 and therefore can't be simplified away).
    //
    // Simplified edges are assigned all labels associated with the edges of
    // the simplified chain.
    //
    // For this option to have any effect, a SnapFunction with a non-zero
    // snap_radius() must be specified.  Also note that vertices specified
    // using ForceVertex are never simplified away.
    //
    // DEFAULT: false
    bool simplify_edge_chains() const;
    void set_simplify_edge_chains(bool simplify_edge_chains);

    // If true, then snapping occurs only when the input geometry does not
    // already meet the S2Builder output guarantees (see the SnapFunction
    // class description for details).  This means that if all input vertices
    // are at snapped locations, all vertex pairs are separated by at least
    // min_vertex_separation(), and all edge-vertex pairs are separated by at
    // least min_edge_vertex_separation(), then no snapping is done.
    //
    // If false, then all vertex pairs and edge-vertex pairs closer than
    // "snap_radius" will be considered for snapping.  This can be useful, for
    // example, if you know that your geometry contains errors and you want to
    // make sure that features closer together than "snap_radius" are merged.
    //
    // This option is automatically turned off by simplify_edge_chains(),
    // since simplifying edge chains is never guaranteed to be idempotent.
    //
    // DEFAULT: true
    bool idempotent() const;
    void set_idempotent(bool idempotent);

    // Specifies that internal memory usage should be tracked using the given
    // S2MemoryTracker.  If a memory limit is specified and more more memory
    // than this is required then an error will be returned.  Example usage:
    //
    //   S2MemoryTracker tracker;
    //   tracker.set_limit(500 << 20);  // 500 MB
    //   S2Builder::Options options;
    //   options.set_memory_tracker(&tracker);
    //   S2Builder builder{options};
    //   ...
    //   S2Error error;
    //   if (!builder.Build(&error)) {
    //     if (error.code() == S2Error::RESOURCE_EXHAUSTED) {
    //       ABSL_LOG(ERROR) << error;  // Memory limit exceeded
    //     }
    //   }
    //
    // CAVEATS:
    //
    //  - Memory allocated by the output S2Builder layers is not tracked.
    //
    //  - While memory tracking is reasonably complete and accurate, it does
    //    not account for every last byte.  It is intended only for the
    //    purpose of preventing clients from running out of memory.
    //
    // DEFAULT: nullptr (memory tracking disabled)
    S2MemoryTracker* memory_tracker() const;
    void set_memory_tracker(S2MemoryTracker* tracker);

    // Options may be assigned and copied.
    Options(const Options& options);
    Options& operator=(const Options& options);

   private:
    std::unique_ptr<SnapFunction> snap_function_;
    bool split_crossing_edges_ = false;
    S1Angle intersection_tolerance_ = S1Angle::Zero();
    bool simplify_edge_chains_ = false;
    bool idempotent_ = true;
    S2MemoryTracker* memory_tracker_ = nullptr;
  };

  class Graph;
  // The following classes are only needed by Layer implementations.
  class GraphOptions;

  // For output layers that represent polygons, there is an ambiguity inherent
  // in spherical geometry that does not exist in planar geometry.  Namely, if
  // a polygon has no edges, does it represent the empty polygon (containing
  // no points) or the full polygon (containing all points)?  This ambiguity
  // also occurs for polygons that consist only of degeneracies, e.g. a
  // degenerate loop with only two edges could be either a degenerate shell in
  // the empty polygon or a degenerate hole in the full polygon.
  //
  // To resolve this ambiguity, an IsFullPolygonPredicate may be specified for
  // each output layer (see AddIsFullPolygonPredicate below).  If the output
  // after snapping consists only of degenerate edges and/or sibling pairs
  // (including the case where there are no edges at all), then the layer
  // implementation calls the given predicate to determine whether the polygon
  // is empty or full except for those degeneracies.  The predicate is given
  // an S2Builder::Graph containing the output edges, but note that in general
  // the predicate must also have knowledge of the input geometry in order to
  // determine the correct result.
  //
  // This predicate is only needed by layers that are assembled into polygons.
  // It is not used by other layer types.
  using IsFullPolygonPredicate =
      std::function<bool (const Graph& g, S2Error* error)>;

  // Default constructor; requires Init() to be called.
  S2Builder();

  // Convenience constructor that calls Init().  Note that to use the default
  // options, C++ syntax requires an extra layer of parentheses:
  //
  //   S2Builder builder{S2Builder::Options()};
  explicit S2Builder(const Options& options);

  // Initializes an S2Builder with the given options.
  void Init(const Options& options);
  const Options& options() const { return options_; }

  // Starts a new output layer.  This method must be called before adding any
  // edges to the S2Builder.  You may call this method multiple times to build
  // multiple geometric objects that are snapped to the same set of sites.
  //
  // For example, if you have a set of contour lines, then you could put each
  // contour line in a separate layer.  This keeps the contour lines separate
  // from each other, while also ensuring that no crossing edges are created
  // when they are snapped and/or simplified.  (This is not true if the
  // contour lines are snapped or simplified independently.)
  //
  // Similarly, if you have a set of polygons that share common boundaries
  // (e.g., countries), you can snap and/or simplify them at the same time by
  // putting them in different layers, while ensuring that their boundaries
  // remain consistent (i.e., no crossing edges or T-vertices are introduced).
  //
  // Ownership of the layer is transferred to the S2Builder.  Example usage:
  //
  // S2Polyline line1, line2;
  // builder.StartLayer(make_unique<s2builderutil::S2PolylineLayer>(&line1)));
  // ... Add edges using builder.AddEdge(), etc ...
  // builder.StartLayer(make_unique<s2builderutil::S2PolylineLayer>(&line2)));
  // ... Add edges using builder.AddEdge(), etc ...
  // S2Error error;
  // ABSL_CHECK(builder.Build(&error)) << error;  // Builds "line1" & "line2"
  class Layer;

  void StartLayer(std::unique_ptr<Layer> layer);

  // Adds a degenerate edge (representing a point) to the current layer.
  void AddPoint(const S2Point& v);

  // Adds the given edge to the current layer.
  void AddEdge(const S2Point& v0, const S2Point& v1);

  // Adds the edges in the given polyline to the current layer.  Note that
  // polylines with 0 or 1 vertices are defined to have no edges.
  void AddPolyline(S2PointSpan polyline);
  void AddPolyline(const S2Polyline& polyline);

  // Adds the edges in the given loop to the current layer.  Note that a loop
  // consisting of one vertex adds a single degenerate edge.
  void AddLoop(S2PointLoopSpan loop);

  // Adds the edges in the given loop to the current layer. Loops with a single
  // vertex are ignored.
  //
  // If the sign() of an S2Loop is negative (i.e. the loop represents a hole
  // within a polygon), the edge directions are automatically reversed to
  // ensure that the polygon interior is always to the left of every edge.
  void AddLoop(const S2Loop& loop);

  // Adds the loops in the given polygon to the current layer.  Loops
  // representing holes have their edge directions automatically reversed as
  // described for AddLoop().  Note that this method does not distinguish
  // between the empty and full polygons, i.e. adding a full polygon has the
  // same effect as adding an empty one.
  void AddPolygon(const S2Polygon& polygon);

  // Adds the edges of the given shape to the current layer.
  void AddShape(const S2Shape& shape);

  // If "vertex" is the intersection point of two edges AB and CD (as computed
  // by S2::GetIntersection()), this method ensures that AB and CD snap to a
  // common vertex.  (Note that the common vertex may be different than
  // "vertex" in order to ensure that no pair of vertices is closer than the
  // given snap radius.)  Unlike Options::split_crossing_edges(), this method
  // may be used to split crossing edge pairs selectively.
  //
  // This method can also be used to tessellate edges using S2::GetPointOnLine()
  // or S2::Project() provided that a suitable intersection tolerance is
  // specified (see intersection_tolerance() for details).
  //
  // This method implicitly overrides the idempotent() option, since adding an
  // intersection point implies a desire to have nearby edges snapped to it
  // even if these edges already satisfy the S2Builder output guarantees.
  // (Otherwise for example edges would never be snapped to nearby
  // intersection points when the snap radius is zero.)
  //
  // Note that unlike ForceVertex(), this method maintains all S2Builder
  // guarantees regarding minimum vertex-vertex separation, minimum
  // edge-vertex separation, and edge chain simplification.
  //
  // REQUIRES: options().intersection_tolerance() > S1Angle::Zero()
  // REQUIRES: "vertex" was computed by S2::GetIntersection() (in order to
  //           guarantee that both edges snap to a common vertex)
  void AddIntersection(const S2Point& vertex);

  // For layers that are assembled into polygons, this method specifies a
  // predicate that is called when the output consists entirely of degenerate
  // edges and/or sibling pairs.  The predicate is given an S2Builder::Graph
  // containing the output edges (if any) and is responsible for deciding
  // whether this graph represents the empty polygon (possibly with degenerate
  // shells) or the full polygon (possibly with degenerate holes).  Note that
  // this cannot be determined from the output edges alone; it also requires
  // knowledge of the input geometry.  (Also see IsFullPolygonPredicate above.)
  //
  // This method should be called at most once per layer; additional calls
  // simply overwrite the previous value for the current layer.
  //
  // The default predicate simply returns false (i.e., degenerate polygons are
  // assumed to be empty).  Arguably it would better to return an error in
  // this case, but the fact is that relatively few clients need to be able to
  // construct full polygons, and it is unreasonable to expect all such
  // clients to supply an appropriate predicate.
  //
  // The reason for having a predicate rather than a boolean value is that the
  // predicate is responsible for determining whether the output polygon is
  // empty or full.  In general the input geometry is not degenerate, but
  // rather collapses into a degenerate configuration due to snapping and/or
  // simplification.
  //
  // TODO(ericv): Provide standard predicates to handle common cases,
  // e.g. valid input geometry that becomes degenerate due to snapping.
  void AddIsFullPolygonPredicate(IsFullPolygonPredicate predicate);

  // A predicate that returns an error indicating that no polygon predicate
  // has been specified.
  static bool IsFullPolygonUnspecified(const S2Builder::Graph& g,
                                       S2Error* error);

  // Returns a predicate that returns a constant value (true or false);
  static IsFullPolygonPredicate IsFullPolygon(bool is_full);

  // Forces a vertex to be located at the given position.  This can be used to
  // prevent certain input vertices from moving.  However if you are trying to
  // preserve input edges, be aware that this option does not prevent edges from
  // being split by new vertices.
  //
  // Forced vertices are subject to the following limitations:
  //
  //  - Forced vertices are never snapped.  This is true even when the given
  //    position is not allowed by the given snap function (e.g. you can force
  //    a vertex at a non-S2CellId center when using S2CellIdSnapFunction).
  //    If you want to ensure that forced vertices obey the snap function
  //    restrictions, you must call snap_function().SnapPoint() explicitly.
  //
  //  - There is no guaranteed minimum separation between pairs of forced
  //    vertices, i.e. snap_function().min_vertex_separation() does not apply.
  //    (This must be true because forced vertices can be placed arbitrarily.)
  //
  //  - There is no guaranteed minimum separation between forced vertices and
  //    non-incident edges, i.e. snap_function().min_edge_vertex_separation()
  //    does not apply.
  //
  //  - Forced vertices are never simplified away (i.e. when simplification is
  //    requested using options().simplify_edge_chains()).
  //
  // All other guarantees continue to hold, e.g. the input topology will always
  // be preserved.
  void ForceVertex(const S2Point& vertex);

  // Every edge can have a set of non-negative integer labels attached to it.
  // When used with an appropriate layer type, you can then retrieve the
  // labels associated with each output edge.  This can be useful when merging
  // or combining data from several sources.  (Note that in many cases it is
  // easier to use separate output layers rather than labels.)
  //
  // Labels are 32-bit non-negative integers.  To support other label types,
  // you can use ValueLexicon to store the set of unique labels seen so far:
  //
  //   ValueLexicon<MyLabel> my_label_lexicon;
  //   builder.set_label(my_label_lexicon.Add(label));
  //
  // The current set of labels is represented as a stack.  This makes it easy
  // to add and remove labels hierarchically (e.g., polygon 5, loop 2).  Use
  // set_label() and clear_labels() if you need at most one label per edge.
  //
  using Label = int32_t;

  // Clear the stack of labels.
  void clear_labels();

  // Add a label to the stack.
  // REQUIRES: label >= 0.
  void push_label(Label label);

  // Remove a label from the stack.
  void pop_label();

  // Convenience function that clears the stack and adds a single label.
  // REQUIRES: label >= 0.
  void set_label(Label label);

  // Performs the requested edge splitting, snapping, simplification, etc, and
  // then assembles the resulting edges into the requested output layers.
  //
  // Returns true if all edges were assembled; otherwise sets "error"
  // appropriately.  Depending on the error, some or all output layers may
  // have been created.  Automatically resets the S2Builder state so that it
  // can be reused.
  //
  // REQUIRES: error != nullptr.
  bool Build(S2Error* error);

  // Clears all input data and resets the builder state.  Any options
  // specified are preserved.
  void Reset();

  ///////////////////////////////////////////////////////////////////////////
  // The following methods may be called at any time, including from
  // S2Builder::Layer implementations.

  // Returns the number of input edges.
  int num_input_edges() const;

  // Returns the endpoints of the given input edge.
  //
  // REQUIRES: 0 <= input_edge_id < num_input_edges()
  S2Shape::Edge input_edge(int input_edge_id) const;

 private:
  //////////////////////  Input Types  /////////////////////////
  // All types associated with the S2Builder inputs are prefixed with "Input".

  // Identifies an input vertex.
  using InputVertexId = int32_t;

  // Defines an input edge.
  using InputEdge = std::pair<InputVertexId, InputVertexId>;

  // Identifies an input edge.
  using InputEdgeId = int32_t;

  // Identifies the set of input edge ids that were snapped to a given edge.
  using InputEdgeIdSetId = int32_t;

  // Sort key for prioritizing input vertices.  (Note that keys are *not*
  // compared using std::less; see SortInputVertices for details.)
  using InputVertexKey = std::pair<S2CellId, InputVertexId>;

  //////////////////////  Output Types  /////////////////////////
  // These types define the output vertices and edges.

  // Identifies a snapped vertex ("snap site").  If there is only one layer,
  // than SiteId is the same as Graph::VertexId, but if there are many layers
  // then each Graph may contain only a subset of the sites.  Also see
  // GraphOptions::allow_vertex_filtering().
  using SiteId = int32_t;

  // Defines an output edge.
  using Edge = std::pair<SiteId, SiteId>;

  // Identifies an output edge.
  using EdgeId = int32_t;

  // Identifies an output edge in a particular layer.
  using LayerEdgeId = std::pair<int, EdgeId>;

  //////////////////////  Internal Types  /////////////////////////
  class EdgeChainSimplifier;

  // MemoryTracker is a helper class to measure S2Builder memory usage.  It is
  // based on a detailed analysis of the data structures used.  This approach
  // is fragile because the memory tracking code needs to be updated whenever
  // S2Builder is modified, however S2Builder has been quite stable and this
  // approach allows the memory usage to be measured quite accurately.
  //
  // CAVEATS:
  //
  //  - Does not track memory used by edge labels.  (It is tricky to do this
  //    accurately because they are stored in an IdSetLexicon, and labels
  //    are typically a tiny fraction of the total space used.)
  //
  //  - Does not track memory used to represent layers internally.  (The
  //    number of layers is typically small compared to the numbers of
  //    vertices and edges, and the amount of memory used by the Layer and
  //    IsFullPolygonPredicate objects is difficult to measure.)
  //
  //  - Does not track memory used by the output layer Build() methods.  (This
  //    includes both temporary space, e.g. due to calling S2Builder::Graph
  //    methods, and also any geometric objects created by these layers.)
  class MemoryTracker : public S2MemoryTracker::Client {
   public:
    bool TallyEdgeSites(const gtl::compact_array<SiteId>& sites);
    bool ReserveEdgeSite(gtl::compact_array<SiteId>* sites);
    bool ClearEdgeSites(std::vector<gtl::compact_array<SiteId>>* edge_sites);

    bool TallyIndexedSite();
    bool FixSiteIndexTally(const S2PointIndex<SiteId>& index);
    bool DoneSiteIndex(const S2PointIndex<SiteId>& index);

    bool TallySimplifyEdgeChains(
        absl::Span<const gtl::compact_array<InputVertexId>> site_vertices,
        absl::Span<const std::vector<Edge>> layer_edges);

    bool TallyFilterVertices(int num_sites,
                             absl::Span<const std::vector<Edge>> layer_edges);
    bool DoneFilterVertices();

   private:
    // The amount of non-inline memory used to store edge sites.
    int64_t edge_sites_bytes_ = 0;

    // The amount of memory used by the S2PointIndex for sites.
    int64_t site_index_bytes_ = 0;

    // The amount of temporary memory used by Graph::FilterVertices().
    int64_t filter_vertices_bytes_ = 0;
  };

  InputVertexId AddVertex(const S2Point& v);
  void ChooseSites();
  void ChooseAllVerticesAsSites();
  std::vector<InputVertexKey> SortInputVertices();
  void AddEdgeCrossings(const MutableS2ShapeIndex& input_edge_index);
  void AddForcedSites(S2PointIndex<SiteId>* site_index);
  bool is_forced(SiteId v) const;
  void ChooseInitialSites(S2PointIndex<SiteId>* site_index);
  S2Point SnapSite(const S2Point& point) const;
  void CollectSiteEdges(const S2PointIndex<SiteId>& site_index);
  void SortSitesByDistance(const S2Point& x,
                           gtl::compact_array<SiteId>* sites) const;
  void InsertSiteByDistance(SiteId new_site_id, const S2Point& x,
                            gtl::compact_array<SiteId>* sites);
  void AddExtraSites(const MutableS2ShapeIndex& input_edge_index);
  void MaybeAddExtraSites(InputEdgeId edge_id, absl::Span<const SiteId> chain,
                          const MutableS2ShapeIndex& input_edge_index,
                          absl::flat_hash_set<InputEdgeId>* edges_to_resnap);
  void AddExtraSite(const S2Point& new_site,
                    const MutableS2ShapeIndex& input_edge_index,
                    absl::flat_hash_set<InputEdgeId>* edges_to_resnap);
  S2Point GetSeparationSite(const S2Point& site_to_avoid,
                            const S2Point& v0, const S2Point& v1,
                            InputEdgeId input_edge_id) const;
  S2Point GetCoverageEndpoint(const S2Point& p, const S2Point& n) const;
  void SnapEdge(InputEdgeId e, std::vector<SiteId>* chain) const;

  void BuildLayers();
  void BuildLayerEdges(
      std::vector<std::vector<Edge>>* layer_edges,
      std::vector<std::vector<InputEdgeIdSetId>>* layer_input_edge_ids,
      IdSetLexicon* input_edge_id_set_lexicon);
  void AddSnappedEdges(
      InputEdgeId begin, InputEdgeId end, const GraphOptions& options,
      std::vector<Edge>* edges, std::vector<InputEdgeIdSetId>* input_edge_ids,
      IdSetLexicon* input_edge_id_set_lexicon,
      std::vector<gtl::compact_array<InputVertexId>>* site_vertices);
  void MaybeAddInputVertex(
      InputVertexId v, SiteId id,
      std::vector<gtl::compact_array<InputVertexId>>* site_vertices) const;
  void AddSnappedEdge(SiteId src, SiteId dst, InputEdgeIdSetId id,
                      EdgeType edge_type, std::vector<Edge>* edges,
                      std::vector<InputEdgeIdSetId>* input_edge_ids) const;
  void SimplifyEdgeChains(
      const std::vector<gtl::compact_array<InputVertexId>>& site_vertices,
      std::vector<std::vector<Edge>>* layer_edges,
      std::vector<std::vector<InputEdgeIdSetId>>* layer_input_edge_ids,
      IdSetLexicon* input_edge_id_set_lexicon);
  void MergeLayerEdges(
      absl::Span<const std::vector<Edge>> layer_edges,
      absl::Span<const std::vector<InputEdgeIdSetId>> layer_input_edge_ids,
      std::vector<Edge>* edges, std::vector<InputEdgeIdSetId>* input_edge_ids,
      std::vector<int>* edge_layers) const;
  static bool StableLessThan(const Edge& a, const Edge& b,
                             const LayerEdgeId& ai, const LayerEdgeId& bi);

  //////////// Parameters /////////////

  // S2Builder options.
  Options options_;

  // The maximum distance (inclusive) that a vertex can move when snapped,
  // equal to S1ChordAngle(options_.snap_function().snap_radius()).
  S1ChordAngle site_snap_radius_ca_;

  // The maximum distance (inclusive) that an edge can move when snapping to a
  // snap site.  It can be slightly larger than the site snap radius when
  // edges are being split at crossings.
  S1ChordAngle edge_snap_radius_ca_;

  // True if we need to check that snapping has not changed the input topology
  // around any vertex (i.e. Voronoi site).  Normally this is only necessary for
  // forced vertices, but if the snap radius is very small (e.g., zero) and
  // split_crossing_edges() is true then we need to do this for all vertices.
  // In all other situations, any snapped edge that crosses a vertex will also
  // be closer than min_edge_vertex_separation() to that vertex, which will
  // cause us to add a separation site anyway.
  bool check_all_site_crossings_;

  S1Angle max_edge_deviation_;
  S1ChordAngle edge_site_query_radius_ca_;
  S1ChordAngle min_edge_length_to_split_ca_;

  S1Angle min_site_separation_;
  S1ChordAngle min_site_separation_ca_;
  S1ChordAngle min_edge_site_separation_ca_;
  S1ChordAngle min_edge_site_separation_ca_limit_;

  S1ChordAngle max_adjacent_site_separation_ca_;

  // The squared sine of the edge snap radius.  This is equivalent to the snap
  // radius (squared) for distances measured through the interior of the
  // sphere to the plane containing an edge.  This value is used only when
  // interpolating new points along edges (see GetSeparationSite).
  double edge_snap_radius_sin2_;

  // A copy of the argument to Build().
  S2Error* error_;

  // True if snapping was requested.  This is true if either snap_radius() is
  // positive, or split_crossing_edges() is true (which implicitly requests
  // snapping to ensure that both crossing edges are snapped to the
  // intersection point).
  bool snapping_requested_;

  // Initially false, and set to true when it is discovered that at least one
  // input vertex or edge does not meet the output guarantees (e.g., that
  // vertices are separated by at least snap_function.min_vertex_separation).
  bool snapping_needed_;

  //////////// Input Data /////////////

  // A flag indicating whether label_set_ has been modified since the last
  // time label_set_id_ was computed.
  bool label_set_modified_;

  std::vector<S2Point> input_vertices_;
  std::vector<InputEdge> input_edges_;

  std::vector<std::unique_ptr<Layer>> layers_;
  std::vector<GraphOptions> layer_options_;
  std::vector<InputEdgeId> layer_begins_;
  std::vector<IsFullPolygonPredicate> layer_is_full_polygon_predicates_;

  // Each input edge has "label set id" (an int32_t) representing the set of
  // labels attached to that edge.  This vector is populated only if at least
  // one label is used.
  using LabelSetId = int32_t;
  std::vector<LabelSetId> label_set_ids_;
  IdSetLexicon label_set_lexicon_;

  // The current set of labels (represented as a stack).
  std::vector<Label> label_set_;

  // The LabelSetId corresponding to the current label set, computed on demand
  // (by adding it to label_set_lexicon()).
  LabelSetId label_set_id_;

  ////////////// Data for Snapping and Simplifying //////////////

  // The number of sites specified using ForceVertex().  These sites are
  // always at the beginning of the sites_ vector.
  SiteId num_forced_sites_;

  // The set of snapped vertex locations ("sites").
  std::vector<S2Point> sites_;

  // A map from each input edge to the set of sites "nearby" that edge,
  // defined as the set of sites that are candidates for snapping and/or
  // avoidance.  Note that compact_array will inline up to two sites, which
  // usually takes care of the vast majority of edges.  Sites are kept sorted
  // by increasing distance from the origin of the input edge.
  //
  // Once snapping is finished, this field is discarded unless edge chain
  // simplification was requested, in which case instead the sites are
  // filtered by removing the ones that each edge was snapped to, leaving only
  // the "sites to avoid" (needed for simplification).
  std::vector<gtl::compact_array<SiteId>> edge_sites_;

  // An object to track the memory usage of this class.
  MemoryTracker tracker_;

  S2Builder(const S2Builder&) = delete;
  S2Builder& operator=(const S2Builder&) = delete;
};

// This class is only needed by S2Builder::Layer implementations.  A layer is
// responsible for assembling an S2Builder::Graph of snapped edges into the
// desired output format (e.g., an S2Polygon).  The GraphOptions class allows
// each Layer type to specify requirements on its input graph: for example, if
// DegenerateEdges::DISCARD is specified, then S2Builder will ensure that all
// degenerate edges are removed before passing the graph to the S2Layer::Build
// method.
class S2Builder::GraphOptions {
 public:
  using EdgeType = S2Builder::EdgeType;
  enum class DegenerateEdges : uint8_t;
  enum class DuplicateEdges : uint8_t;
  enum class SiblingPairs : uint8_t;

  // All S2Builder::Layer subtypes should specify GraphOptions explicitly
  // using this constructor, rather than relying on default values.
  GraphOptions(EdgeType edge_type, DegenerateEdges degenerate_edges,
               DuplicateEdges duplicate_edges, SiblingPairs sibling_pairs)
      : edge_type_(edge_type), degenerate_edges_(degenerate_edges),
        duplicate_edges_(duplicate_edges), sibling_pairs_(sibling_pairs),
        allow_vertex_filtering_(true) {
  }

  // The default options specify that all edges should be kept, since this
  // produces the least surprising output and makes it easier to diagnose the
  // problem when an option is left unspecified.
  GraphOptions() : edge_type_(EdgeType::DIRECTED),
                   degenerate_edges_(DegenerateEdges::KEEP),
                   duplicate_edges_(DuplicateEdges::KEEP),
                   sibling_pairs_(SiblingPairs::KEEP),
                   allow_vertex_filtering_(true) {
  }

  // Specifies whether the S2Builder input edges should be treated as
  // undirected.  If true, then all input edges are duplicated into pairs
  // consisting of an edge and a sibling (reverse) edge.  Note that the
  // automatically created sibling edge has an empty set of labels and does
  // not have an associated InputEdgeId.
  //
  // The layer implementation is responsible for ensuring that exactly one
  // edge from each pair is used in the output, i.e. *only half* of the graph
  // edges will be used.  (Note that some values of the sibling_pairs() option
  // automatically take care of this issue by removing half of the edges and
  // changing edge_type() to DIRECTED.)
  //
  // DEFAULT: EdgeType::DIRECTED
  EdgeType edge_type() const;
  void set_edge_type(EdgeType edge_type);

  // Controls how degenerate edges (i.e., an edge from a vertex to itself) are
  // handled.  Such edges may be present in the input, or they may be created
  // when both endpoints of an edge are snapped to the same output vertex.
  // The options available are:
  //
  // DISCARD: Discards all degenerate edges.  This is useful for layers that
  //          do not support degeneracies, such as S2PolygonLayer.
  //
  // DISCARD_EXCESS: Discards all degenerate edges that are connected to
  //                 non-degenerate edges and merges any remaining duplicate
  //                 degenerate edges.  This is useful for simplifying
  //                 polygons while ensuring that loops that collapse to a
  //                 single point do not disappear.
  //
  // KEEP: Keeps all degenerate edges.  Be aware that this may create many
  //       redundant edges when simplifying geometry (e.g., a polyline of the
  //       form AABBBBBCCCCCCDDDD).  DegenerateEdges::KEEP is mainly useful
  //       for algorithms that require an output edge for every input edge.
  //
  // DEFAULT: DegenerateEdges::KEEP
  enum class DegenerateEdges : uint8_t { DISCARD, DISCARD_EXCESS, KEEP };
  DegenerateEdges degenerate_edges() const;
  void set_degenerate_edges(DegenerateEdges degenerate_edges);

  // Controls how duplicate edges (i.e., edges that are present multiple
  // times) are handled.  Such edges may be present in the input, or they can
  // be created when vertices are snapped together.  When several edges are
  // merged, the result is a single edge labelled with all of the original
  // input edge ids.
  //
  // DEFAULT: DuplicateEdges::KEEP
  enum class DuplicateEdges : uint8_t { MERGE, KEEP };
  DuplicateEdges duplicate_edges() const;
  void set_duplicate_edges(DuplicateEdges duplicate_edges);

  // Controls how sibling edge pairs (i.e., pairs consisting of an edge and
  // its reverse edge) are handled.  Layer types that define an interior
  // (e.g., polygons) normally discard such edge pairs since they do not
  // affect the result (i.e., they define a "loop" with no interior).  The
  // various options include:
  //
  // DISCARD: Discards all sibling edge pairs.
  //
  // DISCARD_EXCESS: Like DISCARD, except that a single sibling pair is kept
  //                 if the result would otherwise be empty.  This is useful
  //                 for polygons with degeneracies (S2LaxPolygonShape), and
  //                 for simplifying polylines while ensuring that they are
  //                 not split into multiple disconnected pieces.
  //
  // KEEP: Keeps sibling pairs.  This can be used to create polylines that
  //       double back on themselves, or degenerate loops (with a layer type
  //       such as S2LaxPolygonShape).
  //
  // REQUIRE: Requires that all edges have a sibling (and returns an error
  //          otherwise).  This is useful with layer types that create a
  //          collection of adjacent polygons (a polygon mesh).
  //
  // CREATE: Ensures that all edges have a sibling edge by creating them if
  //         necessary.  This is useful with polygon meshes where the input
  //         polygons do not cover the entire sphere.  Such edges always have
  //         an empty set of labels and do not have an associated InputEdgeId.
  //
  // If edge_type() is EdgeType::UNDIRECTED, a sibling edge pair is considered
  // to consist of four edges (two duplicate edges and their siblings), since
  // only two of these four edges will be used in the final output.
  //
  // Furthermore, since the options REQUIRE and CREATE guarantee that all
  // edges will have siblings, S2Builder implements these options for
  // undirected edges by discarding half of the edges in each direction and
  // changing the edge_type() to EdgeType::DIRECTED.  For example, two
  // undirected input edges between vertices A and B would first be converted
  // into two directed edges in each direction, and then one edge of each pair
  // would be discarded leaving only one edge in each direction.
  //
  // Degenerate edges are considered not to have siblings.  If such edges are
  // present, they are passed through unchanged by SiblingPairs::DISCARD.  For
  // SiblingPairs::REQUIRE or SiblingPairs::CREATE with undirected edges, the
  // number of copies of each degenerate edge is reduced by a factor of two.
  //
  // Any of the options that discard edges (DISCARD, DISCARD_EXCESS, and
  // REQUIRE/CREATE in the case of undirected edges) have the side effect that
  // when duplicate edges are present, all of the corresponding edge labels
  // are merged together and assigned to the remaining edges.  (This avoids
  // the problem of having to decide which edges are discarded.)  Note that
  // this merging takes place even when all copies of an edge are kept.  For
  // example, consider the graph {AB1, AB2, AB3, BA4, CD5, CD6} (where XYn
  // denotes an edge from X to Y with label "n").  With SiblingPairs::DISCARD,
  // we need to discard one of the copies of AB.  But which one?  Rather than
  // choosing arbitrarily, instead we merge the labels of all duplicate edges
  // (even ones where no sibling pairs were discarded), yielding {AB123,
  // AB123, CD45, CD45} (assuming that duplicate edges are being kept).
  // Notice that the labels of duplicate edges are merged even if no siblings
  // were discarded (such as CD5, CD6 in this example), and that this would
  // happen even with duplicate degenerate edges (e.g. the edges EE7, EE8).
  //
  // DEFAULT: SiblingPairs::KEEP
  enum class SiblingPairs : uint8_t {
    DISCARD,
    DISCARD_EXCESS,
    KEEP,
    REQUIRE,
    CREATE
  };
  SiblingPairs sibling_pairs() const;
  void set_sibling_pairs(SiblingPairs sibling_pairs);

  // This is a specialized option that is only needed by clients that want to
  // work with the graphs for multiple layers at the same time (e.g., in order
  // to check whether the same edge is present in two different graphs).  [Note
  // that if you need to do this, usually it is easier just to build a single
  // graph with suitable edge labels.]
  //
  // When there are a large number of layers, then by default S2Builder builds
  // a minimal subgraph for each layer containing only the vertices needed by
  // the edges in that layer.  This ensures that layer types that iterate over
  // the vertices run in time proportional to the size of that layer rather
  // than the size of all layers combined.  (For example, if there are a
  // million layers with one edge each, then each layer would be passed a
  // graph with 2 vertices rather than 2 million vertices.)
  //
  // If this option is set to false, this optimization is disabled.  Instead
  // the graph passed to this layer will contain the full set of vertices.
  // (This is not recommended when the number of layers could be large.)
  //
  // DEFAULT: true
  bool allow_vertex_filtering() const;
  void set_allow_vertex_filtering(bool allow_vertex_filtering);

 private:
  EdgeType edge_type_;
  DegenerateEdges degenerate_edges_;
  DuplicateEdges duplicate_edges_;
  SiblingPairs sibling_pairs_;
  bool allow_vertex_filtering_;
};

bool operator==(const S2Builder::GraphOptions& x,
                const S2Builder::GraphOptions& y);


//////////////////   Implementation details follow   ////////////////////


// The maximum snap radius is just large enough to support snapping to
// S2CellId level 0.  It is equivalent to 7800km on the Earth's surface.
inline S1Angle S2Builder::SnapFunction::kMaxSnapRadius() {
  // This value can't be larger than 85.7 degrees without changing the code
  // related to min_edge_length_to_split_ca_, and increasing it to 90 degrees
  // or more would most likely require significant changes to the algorithm.
  return S1Angle::Degrees(70);
}

inline const S2Builder::SnapFunction& S2Builder::Options::snap_function()
    const {
  return *snap_function_;
}

inline void S2Builder::Options::set_snap_function(
    const SnapFunction& snap_function) {
  snap_function_ = snap_function.Clone();
}

inline bool S2Builder::Options::split_crossing_edges() const {
  return split_crossing_edges_;
}

inline void S2Builder::Options::set_split_crossing_edges(
    bool split_crossing_edges) {
  split_crossing_edges_ = split_crossing_edges;
}

inline S1Angle S2Builder::Options::intersection_tolerance() const {
  if (!split_crossing_edges()) return intersection_tolerance_;
  return std::max(intersection_tolerance_, S2::kIntersectionError);
}

inline void S2Builder::Options::set_intersection_tolerance(
    S1Angle intersection_tolerance) {
  ABSL_DCHECK_GE(intersection_tolerance, S1Angle::Zero());
  intersection_tolerance_ = intersection_tolerance;
}

inline bool S2Builder::Options::simplify_edge_chains() const {
  return simplify_edge_chains_;
}

inline void S2Builder::Options::set_simplify_edge_chains(
    bool simplify_edge_chains) {
  simplify_edge_chains_ = simplify_edge_chains;

  // Simplification requires a non-zero snap radius, and while it might be
  // possible to do some simplifying without snapping, it is much simpler to
  // always snap (even if the input geometry already meets the other output
  // requirements).  We need to compute edge_sites_ in order to avoid
  // approaching non-incident vertices too closely, for example.
  set_idempotent(false);
}

inline bool S2Builder::Options::idempotent() const {
  return idempotent_;
}

inline void S2Builder::Options::set_idempotent(bool idempotent) {
  idempotent_ = idempotent;
}

inline S2MemoryTracker* S2Builder::Options::memory_tracker() const {
  return memory_tracker_;
}

inline void S2Builder::Options::set_memory_tracker(
    S2MemoryTracker* tracker) {
  memory_tracker_ = tracker;
}

inline S2Builder::GraphOptions::EdgeType
S2Builder::GraphOptions::edge_type() const {
  return edge_type_;
}

inline void S2Builder::GraphOptions::set_edge_type(EdgeType edge_type) {
  edge_type_ = edge_type;
}

inline S2Builder::GraphOptions::DegenerateEdges
S2Builder::GraphOptions::degenerate_edges() const {
  return degenerate_edges_;
}

inline void S2Builder::GraphOptions::set_degenerate_edges(
    DegenerateEdges degenerate_edges) {
  degenerate_edges_ = degenerate_edges;
}

inline S2Builder::GraphOptions::DuplicateEdges
S2Builder::GraphOptions::duplicate_edges() const {
  return duplicate_edges_;
}

inline void S2Builder::GraphOptions::set_duplicate_edges(
    DuplicateEdges duplicate_edges) {
  duplicate_edges_ = duplicate_edges;
}

inline S2Builder::GraphOptions::SiblingPairs
S2Builder::GraphOptions::sibling_pairs() const {
  return sibling_pairs_;
}

inline void S2Builder::GraphOptions::set_sibling_pairs(
    SiblingPairs sibling_pairs) {
  sibling_pairs_ = sibling_pairs;
}

inline bool S2Builder::GraphOptions::allow_vertex_filtering() const {
  return allow_vertex_filtering_;
}

inline void S2Builder::GraphOptions::set_allow_vertex_filtering(
    bool allow_vertex_filtering) {
  allow_vertex_filtering_ = allow_vertex_filtering;
}

inline bool S2Builder::is_forced(SiteId v) const {
  return v < num_forced_sites_;
}

inline void S2Builder::AddPoint(const S2Point& v) {
  AddEdge(v, v);
}

inline int S2Builder::num_input_edges() const {
  return input_edges_.size();
}

inline S2Shape::Edge S2Builder::input_edge(int input_edge_id) const {
  const InputEdge& edge = input_edges_[input_edge_id];
  return S2Shape::Edge(input_vertices_[edge.first],
                       input_vertices_[edge.second]);
}

#endif  // S2_S2BUILDER_H_
