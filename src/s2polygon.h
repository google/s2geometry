// Copyright 2005 Google Inc. All Rights Reserved.
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

#ifndef S2GEOMETRY_S2POLYGON_H_
#define S2GEOMETRY_S2POLYGON_H_

#include <cstddef>
#include <map>
#include <vector>

#include "base/integral_types.h"
#include "base/macros.h"
#include "fpcontractoff.h"
#include "s2.h"
#include "s2cellid.h"
#include "s2latlngrect.h"
#include "s2loop.h"
#include "s2polyline.h"
#include "s2region.h"
#include "s2shapeindex.h"

class Decoder;
class Encoder;
class S1Angle;
class S2Cap;
class S2Cell;
class S2CellUnion;
class S2Error;
class S2Loop;
class S2PolygonBuilder;
class S2Polyline;
class S2VertexFilter;
class S2XYZFaceSiTi;

// An S2Polygon is an S2Region object that represents a polygon.  A polygon is
// defined by zero or more loops; recall that the interior of a loop is
// defined to be its left-hand side (see S2Loop).  There are two different
// conventions for creating an S2Polygon:
//
//   - InitNested() expects the input loops to be nested hierarchically.  The
//     polygon interior then consists of the set of points contained by an odd
//     number of loops.  So for example, a circular region with a hole in it
//     would be defined as two CCW loops, with one loop containing the other.
//     The loops can be provided in any order.
//
//     When the orientation of the input loops is unknown, the nesting
//     requirement is typically met by calling S2Loop::Normalize() on each
//     loop (which inverts the loop if necessary so that it encloses at most
//     half the sphere).  But in fact any set of loops can be used as long as
//     (1) there is no pair of loops that cross, and (2) there is no pair of
//     loops whose union is the entire sphere.
//
//   - InitOriented() expects the input loops to be oriented such that the
//     polygon interior is on the left-hand side of every loop.  So for
//     example, a circular region with a hole in it would be defined using a
//     CCW outer loop and a CW inner loop.  The loop orientations must all be
//     consistent; for example, it is not valid to have one CCW loop nested
//     inside another CCW loop, because the region between the two loops is on
//     the left-hand side of one loop and the right-hand side of the other.
//
// Most clients will not call these methods directly; instead they should use
// S2PolygonBuilder, which has better support for dealing with imperfect data.
//
// When the polygon is initialized, the given loops are automatically
// converted into a canonical form consisting of "shells" and "holes".  Shells
// and holes are both oriented CCW, and are nested hierarchically.  The loops
// are reordered to correspond to a preorder traversal of the nesting
// hierarchy; InitOriented may also invert some loops.
//
// Polygons may represent any region of the sphere with a polygonal boundary,
// including the entire sphere (known as the "full" polygon).  The full
// polygon consists of a single full loop (see S2Loop), whereas the empty
// polygon has no loops at all.
//
// Polygons have the following restrictions:
//
//  - Loops may not cross, i.e. the boundary of a loop may not intersect
//    both the interior and exterior of any other loop.
//
//  - Loops may not share edges, i.e. if a loop contains an edge AB, then
//    no other loop may contain AB or BA.
//
//  - Loops may share vertices, however no vertex may appear twice in a
//    single loop (see S2Loop).
//
//  - No loop may be empty.  The full loop may appear only in the full polygon.

class S2Polygon : public S2Region {
 public:
  // The default constructor creates an empty polygon.  It can be made
  // non-empty by calling Init(), Decode(), etc.
  S2Polygon();

  // Convenience constructor that calls InitNested() with the given loops.
  // Takes ownership of the loops and clears the vector.
  explicit S2Polygon(std::vector<S2Loop*>* loops);

  // Convenience constructor that creates a polygon with a single loop
  // corresponding to the given cell.
  explicit S2Polygon(S2Cell const& cell);

  // Convenience constructor that calls Init(S2Loop*).  Note that this method
  // automatically converts the special empty loop (see S2Loop) into an empty
  // polygon, unlike the vector-of-loops constructor which does not allow
  // empty loops at all.  Takes ownership of the loop.
  explicit S2Polygon(S2Loop* loop);

  // Convenience constructor to disable the automatic validity checking
  // controlled by the --s2debug flag.  Example:
  //
  //   S2Polygon* polygon = new S2Polygon(loops, DISABLE_S2DEBUG);
  //
  // This is equivalent to:
  //
  //   S2Polygon* polygon = new S2Polygon;
  //   polygon->set_s2debug_override(DISABLE_S2DEBUG);
  //   polygon->Init(loops);
  //
  // The main reason to use this constructor is if you intend to call
  // IsValid() explicitly.  See set_s2debug_override() for details.
  S2Polygon(std::vector<S2Loop*>* loops, S2debugOverride override);

  // Create a polygon from a set of hierarchically nested loops.  The polygon
  // interior consists of the points contained by an odd number of loops.
  // (Recall that a loop contains the set of points on its left-hand side.)
  //
  // This method takes ownership of the given loops and clears the given
  // vector.  It then figures out the loop nesting hierarchy and assigns every
  // loop a depth.  Shells have even depths, and holes have odd depths.  Note
  // that the loops are reordered so the hierarchy can be traversed more
  // easily (see GetParent(), GetLastDescendant(), and S2Loop::depth()).
  //
  // This method may be called more than once, in which case any existing
  // loops are deleted before being replaced by the input loops.
  void InitNested(std::vector<S2Loop*>* loops);

  // Like InitNested(), but expects loops to be oriented such that the polygon
  // interior is on the left-hand side of all loops.  This implies that shells
  // and holes should have opposite orientations in the input to this method.
  // (During initialization, loops representing holes will automatically be
  // inverted.)
  void InitOriented(std::vector<S2Loop*>* loops);

  // Initialize a polygon from a single loop.  Note that this method
  // automatically converts the special empty loop (see S2Loop) into an empty
  // polygon, unlike the vector-of-loops Init() method which does not allow
  // empty loops at all.  Takes ownership of the loop.
  void Init(S2Loop* loop);

  // Releases ownership of the loops of this polygon, appends them to "loops" if
  // non-nullptr, and resets the polygon to be empty.  Note that the caller is
  // responsible for deleting the loops whether they pass nullptr or a valid
  // pointer as the loops parameter.  In the former case, they should copy the
  // loop pointers beforehand.
  void Release(std::vector<S2Loop*>* loops);

  // Makes a deep copy of the given source polygon.  The destination polygon
  // will be cleared if necessary.
  void Copy(S2Polygon const* src);

  // Destroys the polygon and frees its loops.
  ~S2Polygon();

  // Allows overriding the automatic validity checks controlled by the
  // --s2debug flag.  If this flag is true, then polygons are automatically
  // checked for validity as they are initialized.  The main reason to disable
  // this flag is if you intend to call IsValid() explicitly, like this:
  //
  //   S2Polygon polygon;
  //   polygon.set_s2debug_override(DISABLE_S2DEBUG);
  //   polygon.Init(...);
  //   if (!polygon.IsValid()) { ... }
  //
  // Without the call to set_s2debug_override(), invalid data would cause a
  // fatal error in Init() whenever the --s2debug flag is enabled.
  //
  // This setting is preserved across calls to Init() and Decode().
  void set_s2debug_override(S2debugOverride override);
  S2debugOverride s2debug_override() const;

  // Returns true if this is a valid polygon (including checking whether all
  // the loops are themselves valid).  Note that validity is checked
  // automatically during initialization when --s2debug is enabled (true by
  // default in debug binaries).
  bool IsValid() const;

  // Returns true if this is *not* a valid polygon and sets "error"
  // appropriately.  Otherwise returns false and leaves "error" unchanged.
  //
  // REQUIRES: error != nullptr
  bool FindValidationError(S2Error* error) const;

  // Return true if this is the empty polygon (consisting of no loops).
  bool is_empty() const { return loops_.empty(); }

  // Return true if this is the full polygon (consisting of a single loop that
  // encompasses the entire sphere).
  bool is_full() const { return num_loops() == 1 && loop(0)->is_full(); }

  // Return the number of loops in this polygon.
  int num_loops() const { return loops_.size(); }

  // Total number of vertices in all loops.
  int num_vertices() const { return num_vertices_; }

  // Return the loop at the given index.  Note that during initialization, the
  // given loops are reordered according to a preorder traversal of the loop
  // nesting hierarchy.  This implies that every loop is immediately followed
  // by its descendants.  This hierarchy can be traversed using the methods
  // GetParent(), GetLastDescendant(), and S2Loop::depth().
  S2Loop const* loop(int k) const { return loops_[k]; }
  S2Loop* loop(int k) { return loops_[k]; }

  // Return the index of the parent of loop k, or -1 if it has no parent.
  int GetParent(int k) const;

  // Return the index of the last loop that is contained within loop k.
  // Returns num_loops() - 1 if k < 0.  Note that loops are indexed according
  // to a preorder traversal of the nesting hierarchy, so the immediate
  // children of loop k can be found by iterating over loops
  // (k+1)..GetLastDescendant(k) and selecting those whose depth is equal to
  // (loop(k)->depth() + 1).
  int GetLastDescendant(int k) const;

  // Return the area of the polygon interior, i.e. the region on the left side
  // of an odd number of loops.  The return value is between 0 and 4*Pi.
  double GetArea() const;

  // Return the true centroid of the polygon multiplied by the area of the
  // polygon (see s2.h for details on centroids).  The result is not unit
  // length, so you may want to normalize it.  Also note that in general, the
  // centroid may not be contained by the polygon.
  //
  // We prescale by the polygon area for two reasons: (1) it is cheaper to
  // compute this way, and (2) it makes it easier to compute the centroid of
  // more complicated shapes (by splitting them into disjoint regions and
  // adding their centroids).
  S2Point GetCentroid() const;

  // If all of the polygon's vertices happen to be the centers of S2Cells at
  // some level, then return that level, otherwise return -1.  See also
  // InitToSnapped() and S2PolygonBuilderOptions::snap_to_cell_centers().
  // Returns -1 if the polygon has no vertices.
  int GetSnapLevel() const;

  // Return the distance from the given point to the polygon interior.  If the
  // polygon is empty, return S1Angle::Infinity().  "x" should be unit length.
  S1Angle GetDistance(S2Point const& x) const;

  // Return the distance from the given point to the polygon boundary.  If the
  // polygon is empty or full, return S1Angle::Infinity() (since the polygon
  // has no boundary).  "x" should be unit length.
  S1Angle GetDistanceToBoundary(S2Point const& x) const;

  // Return the overlap fractions between two polygons, i.e. the ratios of the
  // area of intersection to the area of each polygon.
  static std::pair<double, double> GetOverlapFractions(S2Polygon const* a,
                                                       S2Polygon const* b);

  // If the given point is contained by the polygon, return it.  Otherwise
  // return the closest point on the polygon boundary.  If the polygon is
  // empty, return the input argument.  Note that the result may or may not be
  // contained by the polygon.  "x" should be unit length.
  S2Point Project(S2Point const& x) const;

  // Return the closest point on the polygon boundary to the given point.  If
  // the polygon is empty or full, return the input argument (since the
  // polygon has no boundary).  "x" should be unit length.
  S2Point ProjectToBoundary(S2Point const& x) const;

  // Return true if this polygon contains the given other polygon, i.e.
  // if polygon A contains all points contained by polygon B.
  bool Contains(S2Polygon const* b) const;

  // Returns true if this polgyon (A) approximately contains the given other
  // polygon (B). This is true if it is possible to move the vertices of B
  // no further than "tolerance" such that A contains the modified B.
  //
  // For example, the empty polygon will contain any polygon whose maximum
  // width is no more than "tolerance".
  bool ApproxContains(S2Polygon const* b, S1Angle tolerance) const;

  // Return true if this polygon intersects the given other polygon, i.e.
  // if there is a point that is contained by both polygons.
  bool Intersects(S2Polygon const* b) const;

  // Returns true if this polgyon (A) and the given polygon (B) are
  // approximately disjoint.  This is true if it is possible to ensure that A
  // and B do not intersect by moving their vertices no further than
  // "tolerance".
  //
  // For example, any polygon is approximately disjoint from a polygon whose
  // maximum width is no more than "tolerance".
  bool ApproxDisjoint(S2Polygon const* b, S1Angle tolerance) const;

  // Initialize this polygon to the intersection, union, or difference
  // (A - B) of the given two polygons.  The "vertex_merge_radius" determines
  // how close two vertices must be to be merged together and how close a
  // vertex must be to an edge in order to be spliced into it (see
  // S2PolygonBuilder for details).  By default, the merge radius is just
  // large enough to compensate for errors that occur when computing
  // intersection points between edges (S2EdgeUtil::kIntersectionMergeRadius).
  //
  // If you are going to convert the resulting polygon to a lower-precision
  // format, it is necessary to increase the merge radius in order to get a
  // valid result after rounding (i.e. no duplicate vertices, etc).  For
  // example, if you are going to convert them to geostore::PolygonProto
  // format, then S1Angle::E7(1) is a good value for "vertex_merge_radius".
  void InitToIntersection(S2Polygon const* a, S2Polygon const* b);
  void InitToApproxIntersection(S2Polygon const* a, S2Polygon const* b,
                                S1Angle vertex_merge_radius);
  void InitToUnion(S2Polygon const* a, S2Polygon const* b);
  void InitToApproxUnion(S2Polygon const* a, S2Polygon const* b,
                         S1Angle vertex_merge_radius);
  void InitToDifference(S2Polygon const* a, S2Polygon const* b);
  void InitToApproxDifference(S2Polygon const* a, S2Polygon const* b,
                              S1Angle vertex_merge_radius);

  // Initializes this polygon to a polygon that contains fewer vertices and is
  // within tolerance of the polygon a, with some caveats.
  // If snap_to_cell_centers, the vertices of this polygon will be snapped to
  // the centers of leaf cells at the smallest level that is guaranteed to
  // the polygon valid given the specified tolerance.
  //
  // - If there is a very small island in the original polygon, it may
  //   disappear completely.  Thus some parts of the original polygon
  //   may not be close to the simplified polygon.  Those parts are small,
  //   though, and arguably don't need to be kept.
  // - However, if there are dense islands, they may all disappear, instead
  //   of replacing them by a big simplified island.
  // - For small tolerances (compared to the polygon size), it may happen that
  //   the simplified polygon has more vertices than the original, if the
  //   first step of the simplification creates too many self-intersections.
  //   One can construct irrealistic cases where that happens to an extreme
  //   degree.
  void InitToSimplified(S2Polygon const* a, S1Angle tolerance,
                        bool snap_to_cell_centers);

  // Initializes this polygon to a polygon that contains fewer vertices and is
  // within tolerance of the polygon a, while ensuring that the vertices on the
  // cell boundary are preserved.
  // Precondition: Polygon a is contained in the cell.
  void InitToSimplifiedInCell(S2Polygon const* a, S2Cell const& cell,
                              S1Angle tolerance);

  // Use S2PolygonBuilder to build this polygon by assembling the edges of a
  // given polygon after snapping its vertices to the center of leaf cells.
  // This will simplify the polygon with a tolerance of
  // S2::kMaxDiag.GetValue(S2CellId::kMaxLevel), or approximately
  // 0.13 microdegrees, or 1.5cm on the surface of the Earth.
  // Such a polygon can be efficiently compressed when serialized.
  // The snap level can be changed to a non-leaf level if needed.
  void InitToSnapped(S2Polygon const* polygon,
                     int snap_level = S2CellId::kMaxLevel);

  // Initialize this polygon to the complement of the given polygon.
  void InitToComplement(S2Polygon const* a);

  // Invert the polygon (replace it by its complement).
  void Invert();

  // Intersect this polygon with the polyline "in" and append the resulting
  // zero or more polylines to "out" (which must be empty).  The polylines
  // are appended in the order they would be encountered by traversing "in"
  // from beginning to end.  Note that the output may include polylines with
  // only one vertex, but there will not be any zero-vertex polylines.
  //
  // This is equivalent to calling ApproxIntersectWithPolyline() with the
  // "vertex_merge_radius" set to S2EdgeUtil::kIntersectionMergeRadius.
  void IntersectWithPolyline(S2Polyline const* in,
                             std::vector<S2Polyline*> *out) const;

  // Similar to IntersectWithPolyline(), except that vertices will be
  // dropped as necessary to ensure that all adjacent vertices in the
  // sequence obtained by concatenating the output polylines will be
  // farther than "vertex_merge_radius" apart.  Note that this can change
  // the number of output polylines and/or yield single-vertex polylines.
  void ApproxIntersectWithPolyline(S2Polyline const* in,
                                   std::vector<S2Polyline*> *out,
                                   S1Angle vertex_merge_radius) const;

  // Same as IntersectWithPolyline, but subtracts this polygon from
  // the given polyline.
  void SubtractFromPolyline(S2Polyline const* in,
                            std::vector<S2Polyline*> *out) const;

  // Same as IntersectWithPolylineSloppy, but subtracts this polygon
  // from the given polyline.
  void ApproxSubtractFromPolyline(S2Polyline const* in,
                                  std::vector<S2Polyline*> *out,
                                  S1Angle vertex_merge_radius) const;

  // Return a polygon which is the union of the given polygons.
  // Clears the vector and deletes the polygons!
  static S2Polygon* DestructiveUnion(std::vector<S2Polygon*>* polygons);
  static S2Polygon* DestructiveApproxUnion(std::vector<S2Polygon*>* polygons,
                                           S1Angle vertex_merge_radius);

  // Initialize this polygon to the outline of the given cell union.
  // In principle this polygon should exactly contain the cell union and
  // this polygon's inverse should not intersect the cell union, but rounding
  // issues may cause this not to be the case.
  void InitToCellUnionBorder(S2CellUnion const& cells);

  // Return true if every loop of this polygon shares at most one vertex with
  // its parent loop.  Every polygon has a unique normalized form.  Normalized
  // polygons are useful for testing since it is easy to compare whether two
  // polygons represent the same region.
  bool IsNormalized() const;

  // Return true if two polygons have exactly the same loops.  The loops must
  // appear in the same order, and corresponding loops must have the same
  // linear vertex ordering (i.e., cyclic rotations are not allowed).
  bool Equals(S2Polygon const* b) const;

  // Return true if two polygons have the same boundary.  More precisely, this
  // method requires that both polygons have loops with the same cyclic vertex
  // order and the same nesting hierarchy.  (This implies that vertices may be
  // cyclically rotated between corresponding loops, and the loop ordering may
  // be different between the two polygons as long as the nesting hierarchy is
  // the same.)
  bool BoundaryEquals(S2Polygon const* b) const;

  // Return true if two polygons have the same boundary except for vertex
  // perturbations.  Both polygons must have loops with the same cyclic vertex
  // order and the same nesting hierarchy, but the vertex locations are
  // allowed to differ by up to "max_error".
  bool BoundaryApproxEquals(S2Polygon const* b, double max_error = 1e-15) const;

  // Return true if two polygons have boundaries that are within "max_error"
  // of each other along their entire lengths.  More precisely, there must be
  // a bijection between the two sets of loops such that for each pair of
  // loops, "a_loop->BoundaryNear(b_loop)" is true.
  bool BoundaryNear(S2Polygon const* b, double max_error = 1e-15) const;

  // Return the total number of bytes used by the polygon.
  size_t BytesUsed() const;

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  // GetRectBound() returns essentially tight results, while GetCapBound()
  // might have a lot of extra padding.  Both bounds are conservative in that
  // if the loop contains a point P, then the bound contains P also.
  virtual S2Polygon* Clone() const;
  virtual S2Cap GetCapBound() const;  // Cap surrounding rect bound.
  virtual S2LatLngRect GetRectBound() const { return bound_; }

  virtual bool Contains(S2Cell const& cell) const;
  virtual bool MayIntersect(S2Cell const& cell) const;
  virtual bool VirtualContainsPoint(S2Point const& p) const;

  // The point 'p' does not need to be normalized.
  bool Contains(S2Point const& p) const;

  //  Encode the polygon with about 4 bytes per vertex on Google's geographic
  //  repository, assuming the vertices have all been snapped to the centers of
  //  S2Cells at a given level (typically with InitToSnapped). The other
  //  vertices are stored using 24 bytes. Decoding a polygon encoded this way
  //  always returns the original polygon, without any loss of precision. The
  //  snap level is chosen to be the one that has the most vertices snapped to
  //  S2Cells at that level. If most vertices need 24 bytes, then all vertices
  //  are encoded this way (this method automatically chooses the encoding that
  //  has the best chance of giving the smaller output size).
  virtual void Encode(Encoder* const encoder) const;

  // Decode a polygon encoded with Encode().
  virtual bool Decode(Decoder* const decoder);

  // This function decodes a polygon by pointing the S2Loop vertices directly to
  // the memory of decoder (that needs to be persistent). It is much faster than
  // Decode(). However, it does this only if all the polygon vertices were
  // encoded exactly: this happens only if the polygon wasn't snapped beforehand
  // to a given level. Otherwise it falls back to Decode().
  virtual bool DecodeWithinScope(Decoder* const decoder);

  // Wrapper class for indexing a polygon (see S2ShapeIndex).  Once this
  // object is inserted into an S2ShapeIndex it is owned by that index, and
  // will be automatically deleted when no longer needed by the index.  Note
  // that this class does not take ownership of the polygon; if you want this
  // behavior, see s2shapeutil::S2PolygonOwningShape.  You can also subtype
  // this class to store additional data (see S2Shape for details).

#ifndef SWIG
  class Shape : public S2Shape {
   public:
    Shape() {}  // Must call Init().
    ~Shape();

    // Initialization.  Does not take ownership of "loop".
    explicit Shape(S2Polygon const* polygon);
    void Init(S2Polygon const* polygon);

    S2Polygon const* polygon() const { return polygon_; }

    // S2Shape interface:
    int num_edges() const { return num_edges_; }
    void GetEdge(int e, S2Point const** a, S2Point const** b) const;
    bool has_interior() const { return true; }
    bool contains_origin() const;

   private:
    S2Polygon const* polygon_;
    int num_edges_;
    // An array where element "i" is the total number of edges in loops 0..i-1.
    // This field is only used for polygons that have a large number of loops.
    int* cumulative_edges_;
  };
#endif  // SWIG

 private:
  friend class S2Stats;

  // Given that loops_ contains a single loop, initialize all other fields.
  void InitOneLoop();

  // Compute has_holes_, num_vertices_, bound_, subregion_bound_.
  void InitLoopProperties();

  // Deletes the contents of the loops_ vector and resets the polygon state.
  void ClearLoops();

  // Return true if there is an error in the loop nesting hierarchy.
  bool FindLoopNestingError(S2Error* error) const;

  // A map from each loop to its immediate children with respect to nesting.
  // This map is built during initialization of multi-loop polygons to
  // determine which are shells and which are holes, and then discarded.
  typedef std::map<S2Loop*, std::vector<S2Loop*> > LoopMap;

  void InsertLoop(S2Loop* new_loop, S2Loop* parent, LoopMap* loop_map);
  void InitLoop(S2Loop* loop, int depth, LoopMap* loop_map);

  // Add the polygon's loops to the S2ShapeIndex.  (The actual work of
  // building the index only happens when the index is first used.)
  void InitIndex();

  // Simplifies the polygon either normally or with a vertex filter specifying
  // some vertices to keep in the simplified polygon. In the second case
  // snap_to_cell_centers is ignored. Does not take ownership of a nor of
  // should_keep.
  void InitToSimplifiedInternal(S2Polygon const* a,
                                S1Angle tolerance,
                                bool snap_to_cell_centers,
                                S2VertexFilter const* should_keep);

  // Given an iterator that is already positioned at the S2ShapeIndexCell
  // containing "p", returns Contains(p).
  bool Contains(S2ShapeIndex::Iterator const& it, S2Point const& p) const;

  // Return true if the polygon boundary intersects "target".  It may also
  // return true when the polygon boundary does not intersect "target" but
  // some edge comes within the worst-case error tolerance.
  //
  // REQUIRES: it.id().contains(target.id())
  // [This condition is true whenever it.Locate(target) returns INDEXED.]
  bool BoundaryApproxIntersects(S2ShapeIndex::Iterator const& it,
                                S2Cell const& target) const;

  // Encode the polygon's S2Points directly as three doubles using
  // (40 + 43 * num_loops + 24 * num_vertices) bytes.
  void EncodeLossless(Encoder* encoder) const;

  // Decode a polygon encoded with EncodeLossless().  Used by the Decode and
  // DecodeWithinScope methods above.  The within_scope parameter specifies
  // whether to call DecodeWithinScope on the loops.
  bool DecodeLossless(Decoder* const decoder, bool within_scope);

  // Encode the polygon's vertices using about 4 bytes / vertex plus 24 bytes /
  // unsnapped vertex. All the loop vertices must be converted first to the
  // S2XYZFaceSiTi format using S2Loop::GetXYZFaceSiTiVertices, and concatenated
  // in the all_vertices array.
  //
  // REQUIRES: snap_level >= 0.
  void EncodeCompressed(Encoder* encoder, S2XYZFaceSiTi const* all_vertices,
                        int snap_level) const;

  // Decode a polygon encoded with EncodeCompressed().
  bool DecodeCompressed(Decoder* decoder);

  // Internal implementation of intersect/subtract polyline functions above.
  void InternalClipPolyline(bool invert,
                            S2Polyline const* a,
                            std::vector<S2Polyline*> *out,
                            S1Angle vertex_merge_radius) const;

  int CompareBoundary(S2Loop const* b) const;
  bool ContainsBoundary(S2Polygon const* b) const;
  bool ExcludesBoundary(S2Polygon const* b) const;
  bool ContainsNonCrossingBoundary(S2Loop const* b, bool reverse_b) const;
  bool ExcludesNonCrossingShells(S2Polygon const* b) const;
  bool ExcludesNonCrossingComplementShells(S2Polygon const* b) const;
  bool AnyLoopContains(S2Loop const* b) const;
  bool AnyLoopIntersects(S2Loop const* b) const;

  static void ClipBoundary(S2Polygon const* a, bool reverse_a,
                           S2Polygon const* b, bool invert_b,
                           bool add_shared_edges, S2PolygonBuilder* builder);

  std::vector<S2Loop*> loops_;
  bool has_holes_;          // Polygon has at least one hole.

  // Allows overriding the automatic validity checking controlled by the
  // --s2debug flag.
  uint8 s2debug_override_;  // Store enum in 1 byte rather than 4.

  // True if InitOriented() was called and the given loops had inconsistent
  // orientations (i.e., it is not possible to construct a polygon such that
  // the interior is on the left-hand side of all loops).  We need to remember
  // this error so that it can be returned later by FindValidationError(),
  // since it is not possible to detect this error once the polygon has been
  // initialized.  This field is not preserved by Encode/Decode.
  uint8 error_inconsistent_loop_orientations_;

  // Cache for num_vertices().
  int num_vertices_;

  // In general we build the index the first time it is needed, but we make an
  // exception for Contains(S2Point) because this method has a simple brute
  // force implementation that is also relatively cheap.  For this one method
  // we keep track of the number of calls made and only build the index once
  // enough calls have been made that we think an index would be worthwhile.
  mutable Atomic32 unindexed_contains_calls_;

  // "bound_" is a conservative bound on all points contained by this polygon:
  // if A.Contains(P), then A.bound_.Contains(S2LatLng(P)).
  S2LatLngRect bound_;

  // Since "bound_" is not exact, it is possible that a polygon A contains
  // another polygon B whose bounds are slightly larger.  "subregion_bound_"
  // has been expanded sufficiently to account for this error, i.e.
  // if A.Contains(B), then A.subregion_bound_.Contains(B.bound_).
  S2LatLngRect subregion_bound_;

  // Spatial index of all the polygon loops.
  S2ShapeIndex index_;

#ifndef SWIG
  S2Polygon(S2Polygon const&) = delete;
  void operator=(S2Polygon const&) = delete;
#endif
};

#endif  // S2GEOMETRY_S2POLYGON_H_
