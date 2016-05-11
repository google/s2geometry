// Copyright 2006 Google Inc. All Rights Reserved.
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

#ifndef S2_S2POLYGONBUILDER_H_
#define S2_S2POLYGONBUILDER_H_

#include <unordered_map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "s2/base/macros.h"
#include "s2/fpcontractoff.h"
#include "s2/s1angle.h"
#include "s2/s2.h"
#include "s2/util/math/matrix3x3.h"

class S2Loop;
class S2Polygon;

// This is a simple class for assembling polygons out of edges.  It requires
// that no two edges cross.  It can handle both directed and undirected edges,
// and optionally it can also remove duplicate edge pairs (consisting of two
// identical edges or an edge and its reverse edge).  This is useful for
// computing seamless unions of polygons that have been cut into pieces.
//
// Here are some of the situations this class was designed to handle:
//
// 1. Computing the union of disjoint polygons that may share part of their
//    boundaries.  For example, reassembling a lake that has been split into
//    two loops by a state boundary.
//
// 2. Constructing polygons from input data that does not follow S2
//    conventions, i.e. where loops may have repeated vertices, or distinct
//    loops may share edges, or shells and holes have opposite or unspecified
//    orientations.
//
// 3. Computing the symmetric difference of a set of polygons whose edges
//    intersect only at vertices.  This can be used to implement a limited
//    form of polygon intersection or subtraction as well as unions.
//
// 4. As a tool for implementing other polygon operations by generating a
//    collection of directed edges and then assembling them into loops.
//
// 5. Making a polygon robust to serialization operations
//    that change precision.  See set_robustness_radius() below.
//
// Caveat: Because S2PolygonBuilder only works with polygon boundaries, it
// cannot distinguish between the empty and full polygon.  If your application
// can generate both the empty and full polygons, you must implement logic
// outside of this class (see AssemblePolygon below).

// S2PolygonBuilderOptions is intended to be copied by value as desired.
class S2PolygonBuilderOptions {
 public:
  S2PolygonBuilderOptions() :
      undirected_edges_(false), xor_edges_(true), validate_(false),
      s2debug_override_(S2Debug::ALLOW), snap_to_cell_centers_(false),
      vertex_merge_radius_(S1Angle::Radians(0)), edge_splice_fraction_(0.866) {
  }

  // These are the options that should be used for assembling well-behaved
  // input data into polygons.  All edges should be directed such that
  // "shells" and "holes" have opposite orientations (typically CCW shells and
  // clockwise holes), unless it is known that shells and holes do not share
  // any edges.
  inline static S2PolygonBuilderOptions DIRECTED_XOR();

  // These are the options that should be used for assembling polygons that do
  // not follow the conventions above, e.g. where edge directions may vary
  // within a single loop, or shells and holes are not oppositely oriented and
  // may share edges.
  inline static S2PolygonBuilderOptions UNDIRECTED_XOR();

  // These are the options that should be used for assembling edges where the
  // desired output is a collection of loops rather than a polygon, and edges
  // may occur more than once.  Edges are treated as undirected and are not
  // XORed together.
  inline static S2PolygonBuilderOptions UNDIRECTED_UNION();

  // Default value: false.
  //
  // If "undirected_edges" is false, then the input is assumed to consist of
  // edges that can be assembled into oriented loops without reversing any of
  // the edges.  Otherwise, "undirected_edges" should be set to true.
  bool undirected_edges() const { return undirected_edges_; }
  void set_undirected_edges(bool undirected_edges);

  // Default value: true.
  //
  // If "xor_edges" is true, then any duplicate edge pairs are removed.  This
  // is useful for computing the union of a collection of polygons whose
  // interiors are disjoint but whose boundaries may share some common edges
  // (e.g. computing the union of South Africa, Lesotho, and Swaziland).
  //
  // Note that for directed edges, a "duplicate edge pair" consists of an edge
  // and its corresponding reverse edge.  This means that either (a) "shells"
  // and "holes" must have opposite orientations, or (b) shells and holes do
  // not share edges.  Otherwise undirected_edges() should be specified.
  //
  // There are only two reasons to turn off xor_edges():
  //
  //  (1) AssemblePolygon() will be called, and you want to assert that there
  //      are no duplicate edge pairs in the input.
  //
  //  (2) AssembleLoops() will be called, and you want to keep abutting loops
  //      separate in the output rather than merging their regions together
  //      (e.g. assembling loops for Kansas City, KS and Kansas City, MO
  //      simultaneously).
  bool xor_edges() const { return xor_edges_; }
  void set_xor_edges(bool xor_edges);

  // Default value: false.
  //
  // If true, IsValid() is called on all loops and polygons before
  // constructing them.  If any loop is invalid (e.g. self-intersecting), it
  // is rejected and returned as a set of "unused edges".  Any remaining valid
  // loops are kept.  If the entire polygon is invalid (e.g. two loops
  // intersect), then all loops are rejected and returned as unused edges.
  //
  // When this option is turned on, the --s2debug flag is automatically
  // disabled (see set_s2debug_override) to avoid fatal assertions in debug
  // mode when invalid geometry is encountered.
  bool validate() const { return validate_; }
  void set_validate(bool validate);

  // Default value: S2Debug::ALLOW
  //
  // This method can be used to turn off the automatic validity checks
  // controlled by the --s2debug flag (which is on by default in debug mode).
  // The only reason to turn off these checks is if you intend to call
  // IsValid() or FindValidationError() explicitly, e.g. because you would
  // like to report a specific error.
  //
  // Note that AssemblePolygon() always sets polygon->s2debug_override() to
  // the value specified here (any previous value is ignored).
  S2Debug s2debug_override() const { return s2debug_override_; }
  void set_s2debug_override(S2Debug override);

  // Default value: false.
  //
  // If true, the built polygon will have its vertices snapped to the centers
  // of s2 cells at the smallest level number such that no vertex will move
  // by more than the robustness radius.  If the robustness radius is smaller
  // than half the leaf cell diameter, no snapping will occur.  This is useful
  // because snapped polygons can be Encode()d using less space.
  // See S2EncodePointsCompressed.
  bool snap_to_cell_centers() const { return snap_to_cell_centers_; }
  void set_snap_to_cell_centers(bool snap_to_cell_centers);

  // Default value: 0.
  //
  // If set to a positive value, all vertex pairs that are separated by less
  // than this distance will be merged together.  Note that vertices can move
  // arbitrarily far if there is a long chain of vertices separated by less
  // than this distance.
  //
  // This method is useful for assembling polygons out of input data where
  // vertices and/or edges may not be perfectly aligned.
  S1Angle vertex_merge_radius() const { return vertex_merge_radius_; }
  void set_vertex_merge_radius(S1Angle vertex_merge_radius);

  // Default value: 0.866 (approximately sqrt(3)/2).
  //
  // The edge splice radius is automatically set to this fraction of the vertex
  // merge radius.  If the edge splice radius is positive, then all vertices
  // that are closer than this distance to an edge are spliced into that edge.
  // Note that edges can move arbitrarily far if there is a long chain of
  // vertices in just the right places.
  //
  // You can turn off edge splicing by setting this value to zero.  This will
  // save some time if you don't need this feature, or you don't want vertices
  // to be spliced into nearby edges for some reason.
  //
  // Note that the edge splice fraction must be less than sqrt(3)/2 in order to
  // avoid infinite loops in the merge algorithm.  The default value is very
  // close to this maximum and therefore results in the maximum amount of edge
  // splicing for a given vertex merge radius.
  //
  // The only reason to reduce the edge splice fraction is if you want to limit
  // changes in edge direction due to splicing.  The direction of an edge can
  // change by up to asin(edge_splice_fraction) due to each splice.  Thus by
  // default, edges are allowed to change direction by up to 60 degrees per
  // splice.  However, note that most direction changes are much smaller than
  // this: the worst case occurs only if the vertex being spliced is just
  // outside the vertex merge radius from one of the edge endpoints.
  double edge_splice_fraction() const { return edge_splice_fraction_; }
  void set_edge_splice_fraction(double edge_splice_fraction);

  // The lossless way to serialize a polygon takes 24 bytes per vertex
  // (3 doubles).  If you want a representation with fewer bits, you
  // need to snap your vertices to a grid.  If a vertex is very close
  // to an edge, the snapping operation can lead to self-intersection.
  // To guarantee that a polygon remains valid when its vertices are
  // moved by an angle of up to epsilon, you need:
  //       vertex_merge_radius * edge_splice_fraction >= 2 * epsilon,
  // as the edge and the vertex can each move by epsilon upon snapping.
  // The edge_splice_fraction cannot be zero to have a
  // robustness guarantee.
  //
  // If your grid has maximum diameter d, call set_robustness_radius(d/2).
  //
  // This getter/setter computes vertex_merge_radius based on the
  // robustness requirement.
  S1Angle GetRobustnessRadius() const;
  void SetRobustnessRadius(S1Angle robustness_radius);

  // If snap_to_cell_centers is true, returns the minimum level at which
  // snapping a point to the center of a cell at that level will move the
  // point by no more than the robustness radius.  Returns -1 if there is
  // no such level, or if snap_to_cell_centers is false.
  int GetSnapLevel() const;

 private:
  bool undirected_edges_;
  bool xor_edges_;
  bool validate_;
  S2Debug s2debug_override_;
  bool snap_to_cell_centers_;
  S1Angle vertex_merge_radius_;
  double edge_splice_fraction_;
};

class S2PolygonBuilder {
 public:
  explicit S2PolygonBuilder(S2PolygonBuilderOptions const& options);
  ~S2PolygonBuilder();

  S2PolygonBuilderOptions const& options() const { return options_; }

  // Add the given edge to the polygon builder.  This method should be used
  // for input data that may not follow S2 polygon conventions.  Note that
  // edges are not allowed to cross each other.  Also note that as a
  // convenience, edges where v0 == v1 are ignored.
  //
  // Returns true if an edge was added, and false if an edge was erased
  // (due to XORing) or not added (if both endpoints were the same).
  bool AddEdge(S2Point const& v0, S2Point const& v1);

  // Add all edges in the given loop.  If the sign() of the loop is negative
  // (i.e. this loop represents a hole), the reverse edges are added instead.
  // This implies that "shells" are CCW and "holes" are CW, as required for
  // the directed edges convention described above.
  //
  // This method does not take ownership of the loop.
  void AddLoop(S2Loop const* loop);

  // Add all loops in the given polygon.  Shells and holes are added with
  // opposite orientations as described for AddLoop().  Note that this method
  // does not distinguish between the empty and full polygons, i.e. adding a
  // full polygon has the same effect as adding an empty one.
  //
  // This method does not take ownership of the polygon.
  void AddPolygon(S2Polygon const* polygon);

  // This type is used to return any edges that could not be assembled.
  using EdgeList = std::vector<std::pair<S2Point, S2Point>>;

  // Assembles the given edges into as many non-crossing loops as possible.
  // When there is a choice about how to assemble the loops, then CCW loops
  // are preferred.  Returns true if all edges were assembled.  If
  // "unused_edges" is not nullptr, it is initialized to the set of edges that
  // could not be assembled into loops.
  //
  // Note that if xor_edges() is false and duplicate edge pairs may be
  // present, then undirected_edges() should be specified unless all loops can
  // be assembled in a counter-clockwise direction.  Otherwise this method may
  // not be able to assemble all loops due to its preference for CCW loops.
  //
  // This method resets the S2PolygonBuilder state so that it can be reused.
  bool AssembleLoops(std::vector<S2Loop*>* loops, EdgeList* unused_edges);

  // Like AssembleLoops, but then assembles the loops into a polygon.  If the
  // edges were directed, then it is expected that holes and shells will have
  // opposite orientations, and the polygon interior is to the left of all
  // edges.  If the edges were undirected, then all loops are first normalized
  // so that they enclose at most half of the sphere, and the polygon interior
  // consists of points enclosed by an odd number of loops.
  //
  // Returns true if all of the edges were assembled into loops.  Note that a
  // partial polygon is returned even when the return value is false.
  //
  // If the validate() option is set and the assembled polygon is not valid,
  // "polygon" is reset to empty and false is returned.  Note that for the
  // polygon to be valid, there must not be any duplicate edges in the input.
  // Duplicate edges can be eliminated using the xor_edges() option.
  //
  // Note that because the polygon is constructed from its boundary, this
  // method cannot distinguish between the empty and full polygons.  An empty
  // boundary always yields an empty polygon.  If the result should sometimes
  // be the full polygon, such logic must be implemented outside of this class
  // (and will need to consider factors other than the polygon's boundary).
  // For example, it is often possible to estimate the polygon area.
  bool AssemblePolygon(S2Polygon* polygon, EdgeList* unused_edges);

  // This function is only for debugging.  If it is called, then points will
  // be transformed by the inverse of the given matrix before being printed as
  // lat-lng coordinates in degrees.  "m" should be orthonormal.
  void set_debug_matrix(Matrix3x3_d const& m);

 protected:
  // These functions print either a single vertex, all edges from a single
  // vertex, or all edges in the builder.
  void DumpVertex(S2Point const& v) const;
  void DumpEdges(S2Point const& v0) const;
  void Dump() const;

 private:
  // Return true if the given edge exists.
  bool HasEdge(S2Point const& v0, S2Point const& v1);

  // Erase an edge or an entire loop.  The edge/loop must exist.
  void EraseEdge(S2Point const& v0, S2Point const& v1);
  void EraseLoop(S2Point const* v, int n);

  // Assembles and returns a single loop starting with the given edge.
  // If a loop cannot be assembled starting from this edge, returns nullptr
  // and updates "unused_edges".
  std::unique_ptr<S2Loop> AssembleLoop(S2Point const& v0, S2Point const& v1,
                                       EdgeList* unused_edges);

  // Adds all the given edges to "unused_edges".
  void RejectLoop(S2Point const* v, int n, EdgeList* unused_edges);

  // Builds a map indicating which vertices need to be moved from their
  // current position to a new position, and also returns a spatial index
  // containing all of the vertices that do not need to be moved.
  class PointIndex;

  using MergeMap = std::unordered_map<S2Point, S2Point, S2PointHash>;
  void BuildMergeMap(PointIndex* index, MergeMap* merge_map);

  // Moves a set of vertices from old to new positions.
  void MoveVertices(MergeMap const& merge_map);

  // Modifies each edge by splicing in any vertices whose distance to the edge
  // is at most (edge_splice_fraction() * vertex_merge_radius()).
  void SpliceEdges(PointIndex const& index);

  S2PolygonBuilderOptions options_;

  // This is only used for debugging purposes.
  std::unique_ptr<Matrix3x3_d> debug_matrix_;

  // The current set of edges, grouped by origin.  The set of destination
  // vertices is a multiset so that the same edge can be present more than
  // once.  We could have also used a multiset<pair<S2Point, S2Point>>,
  // but this representation is a bit more convenient.
  using VertexSet = std::multiset<S2Point>;
  using EdgeSet = std::unordered_map<S2Point, VertexSet, S2PointHash>;
  EdgeSet edges_;

  // Unique collection of the starting (first) vertex of all edges,
  // in the order they are added to edges_.
  std::vector<S2Point> starting_vertices_;

  S2PolygonBuilder(S2PolygonBuilder const&) = delete;
  void operator=(S2PolygonBuilder const&) = delete;
};

inline S2PolygonBuilderOptions S2PolygonBuilderOptions::DIRECTED_XOR() {
  return S2PolygonBuilderOptions();
}

inline S2PolygonBuilderOptions S2PolygonBuilderOptions::UNDIRECTED_XOR() {
  S2PolygonBuilderOptions options;
  options.set_undirected_edges(true);
  return options;
}

inline S2PolygonBuilderOptions S2PolygonBuilderOptions::UNDIRECTED_UNION() {
  S2PolygonBuilderOptions options;
  options.set_undirected_edges(true);
  options.set_xor_edges(false);
  return options;
}

#endif  // S2_S2POLYGONBUILDER_H_
