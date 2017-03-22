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

#ifndef S2_S2EDGEUTIL_H_
#define S2_S2EDGEUTIL_H_

#include <cmath>

#include <glog/logging.h>

#include "s2/base/macros.h"
#include "s2/third_party/absl/container/inlined_vector.h"
#include "s2/fpcontractoff.h"
#include "s2/r2.h"
#include "s2/r2rect.h"
#include "s2/s1angle.h"
#include "s2/s1chordangle.h"
#include "s2/s1interval.h"
#include "s2/s2latlng.h"
#include "s2/s2latlngrect.h"
#include "s2/s2pointutil.h"
#include "s2/s2predicates.h"
#include "s2/util/math/vector.h"

class ExactFloat;

// This class contains various utility functions related to edges.  It
// collects together common code that is needed to implement polygonal
// geometry such as polylines, loops, and general polygons.
class S2EdgeUtil {
 public:
  class CopyingEdgeCrosser;  // Forward declaration

  // This class allows edges to be efficiently tested for intersection with a
  // given fixed edge AB.  It is especially efficient when testing for
  // intersection with an edge chain connecting vertices v0, v1, v2, ...
  //
  // Example usage:
  //
  //   void CountIntersections(S2Point const& a, S2Point const& b,
  //                           vector<pair<S2Point, S2Point>> const& edges) {
  //     int count = 0;
  //     EdgeCrosser crosser(&a, &b);
  //     for (auto const& edge : edges) {
  //       if (crosser.CrossingSign(&edge.first, &edge.second) >= 0) {
  //         ++count;
  //       }
  //     }
  //     return count;
  //   }
  //
  // This class expects that the client already has all the necessary vertices
  // stored in memory, so that this class can refer to them with pointers and
  // does not need to make its own copies.  If this is not the case (e.g., you
  // want to pass temporary objects as vertices), see CopyingEdgeCrosser.
  class EdgeCrosser {
   public:
    // Default constructor; must be followed by a call to Init().
    EdgeCrosser() {}

    // Convenience constructor that calls Init() with the given fixed edge AB.
    // The arguments "a" and "b" must point to values that persist for the
    // lifetime of the EdgeCrosser object (or until the next Init() call).
    EdgeCrosser(S2Point const* a, S2Point const* b);

    // Initialize the EdgeCrosser with the given fixed edge AB.  The arguments
    // "a" and "b" must point to values that persist for the lifetime of the
    // EdgeCrosser object (or until the next Init() call).
    void Init(S2Point const* a, S2Point const* b);

    // This function determines whether the edge AB intersects the edge CD.
    // Returns +1 if AB crosses CD at a point that is interior to both edges.
    // Returns  0 if any two vertices from different edges are the same.
    // Returns -1 if there is no crossing.
    // Returns -1 or 0 if either edge is degenerate (A == B or C == D).
    //
    // Properties of CrossingSign:
    //
    //  (1) CrossingSign(b,a,c,d) == CrossingSign(a,b,c,d)
    //  (2) CrossingSign(c,d,a,b) == CrossingSign(a,b,c,d)
    //  (3) CrossingSign(a,b,c,d) == 0 if a==c, a==d, b==c, b==d
    //  (3) CrossingSign(a,b,c,d) <= 0 if a==b or c==d
    //
    // This function implements an exact, consistent perturbation model such
    // that no three points are ever considered to be collinear.  This means
    // that even if you have 4 points A, B, C, D that lie exactly in a line
    // (say, around the equator), C and D will be treated as being slightly to
    // one side or the other of AB.  This is done in a way such that the
    // results are always consistent (see s2pred::Sign).
    //
    // Note that if you want to check an edge against a chain of other edges,
    // it is slightly more efficient to use the single-argument version of
    // CrossingSign below.
    //
    // The arguments must point to values that persist until the next call.
    int CrossingSign(S2Point const* c, S2Point const* d);

    // This method extends the concept of a "crossing" to the case where AB
    // and CD have a vertex in common.  The two edges may or may not cross,
    // according to the rules defined in VertexCrossing() below.  The rules
    // are designed so that point containment tests can be implemented simply
    // by counting edge crossings.  Similarly, determining whether one edge
    // chain crosses another edge chain can be implemented by counting.
    //
    // Returns true if CrossingSign(c, d) > 0, or AB and CD share a vertex
    // and VertexCrossing(a, b, c, d) returns true.
    //
    // The arguments must point to values that persist until the next call.
    bool EdgeOrVertexCrossing(S2Point const* c, S2Point const* d);

    ///////////////////////// Edge Chain Methods ///////////////////////////
    //
    // You don't need to use these unless you're trying to squeeze out every
    // last drop of performance.  Essentially all you are saving is a test
    // whether the first vertex of the current edge is the same as the second
    // vertex of the previous edge.  Example usage:
    //
    //   vector<S2Point> chain;
    //   crosser.RestartAt(&chain[0]);
    //   for (int i = 1; i < chain.size(); ++i) {
    //     if (crosser.EdgeOrVertexCrossing(&chain[i])) { ++count; }
    //   }

    // Convenience constructor that uses AB as the fixed edge, and C as the
    // first vertex of the vertex chain (equivalent to calling RestartAt(c)).
    //
    // The arguments must point to values that persist until the next call.
    EdgeCrosser(S2Point const* a, S2Point const* b, S2Point const* c);

    // Call this method when your chain 'jumps' to a new place.
    // The argument must point to a value that persists until the next call.
    void RestartAt(S2Point const* c);

    // Like CrossingSign above, but uses the last vertex passed to one of
    // the crossing methods (or RestartAt) as the first vertex of the current
    // edge.
    //
    // The argument must point to a value that persists until the next call.
    int CrossingSign(S2Point const* d);

    // Like EdgeOrVertexCrossing above, but uses the last vertex passed to one
    // of the crossing methods (or RestartAt) as the first vertex of the
    // current edge.
    //
    // The argument must point to a value that persists until the next call.
    bool EdgeOrVertexCrossing(S2Point const* d);

   private:
    friend class CopyingEdgeCrosser;

    // These functions handle the "slow path" of CrossingSign().
    int CrossingSignInternal(S2Point const* d);
    int CrossingSignInternal2(S2Point const& d);

    // Used internally by CopyingEdgeCrosser.  Updates "c_" only.
    void set_c(S2Point const* c) { c_ = c; }

    // The fields below are constant after the call to Init().
    S2Point const* a_;
    S2Point const* b_;
    Vector3_d a_cross_b_;

    // To reduce the number of calls to s2pred::ExpensiveSign(), we compute an
    // outward-facing tangent at A and B if necessary.  If the plane
    // perpendicular to one of these tangents separates AB from CD (i.e., one
    // edge on each side) then there is no intersection.
    bool have_tangents_;  // True if the tangents have been computed.
    S2Point a_tangent_;   // Outward-facing tangent at A.
    S2Point b_tangent_;   // Outward-facing tangent at B.

    // The fields below are updated for each vertex in the chain.
    S2Point const* c_;       // Previous vertex in the vertex chain.
    int acb_;                // The orientation of triangle ACB.

    // The field below is a temporary used by CrossingSignInternal().
    int bda_;                // The orientation of triangle BDA.

    EdgeCrosser(EdgeCrosser const&) = delete;
    void operator=(EdgeCrosser const&) = delete;
  };

  // CopyingEdgeCrosser is exactly like EdgeCrosser, except that it makes its
  // own copy of all arguments so that they do not need to persist between
  // calls.  This is less efficient, but makes it possible to use points that
  // are generated on demand and cannot conveniently be stored by the client.
  class CopyingEdgeCrosser {
   public:
    // These methods are all exactly like EdgeCrosser, except that the
    // arguments can be temporaries.
    CopyingEdgeCrosser() {}
    CopyingEdgeCrosser(S2Point const& a, S2Point const& b);
    void Init(S2Point const& a, S2Point const& b);
    int CrossingSign(S2Point const& c, S2Point const& d);
    bool EdgeOrVertexCrossing(S2Point const& c, S2Point const& d);
    CopyingEdgeCrosser(S2Point const& a, S2Point const& b, S2Point const& c);
    void RestartAt(S2Point const& c);
    int CrossingSign(S2Point const& d);
    bool EdgeOrVertexCrossing(S2Point const& d);

   private:
    S2Point a_, b_, c_;
    EdgeCrosser crosser_;

    CopyingEdgeCrosser(CopyingEdgeCrosser const&) = delete;
    void operator=(CopyingEdgeCrosser const&) = delete;
  };

  // This class computes a bounding rectangle that contains all edges defined
  // by a vertex chain v0, v1, v2, ...  All vertices must be unit length.
  // Note that the bounding rectangle of an edge can be larger than the
  // bounding rectangle of its endpoints, e.g. consider an edge that passes
  // through the north pole.
  //
  // The bounds are calculated conservatively to account for numerical errors
  // when S2Points are converted to S2LatLngs.  More precisely, this class
  // guarantees the following.  Let L be a closed edge chain (loop) such that
  // the interior of the loop does not contain either pole.  Now if P is any
  // point such that L.Contains(P), then RectBound(L).Contains(S2LatLng(P)).
  class RectBounder {
   public:
    RectBounder() : bound_(S2LatLngRect::Empty()) {}

    // This method is called to add each vertex to the chain.  Requires that 'b'
    // has unit length.
    void AddPoint(S2Point const& b);

    // Return the bounding rectangle of the edge chain that connects the
    // vertices defined so far.  This bound satisfies the guarantee made
    // above, i.e. if the edge chain defines a loop, then the bound contains
    // the S2LatLng coordinates of all S2Points contained by the loop.
    S2LatLngRect GetBound() const;

    // Expand a bound returned by GetBound() so that it is guaranteed to
    // contain the bounds of any subregion whose bounds are computed using
    // this class.  For example, consider a loop L that defines a square.
    // GetBound() ensures that if a point P is contained by this square, then
    // S2LatLng(P) is contained by the bound.  But now consider a diamond
    // shaped loop S contained by L.  It is possible that GetBound() returns a
    // *larger* bound for S than it does for L, due to rounding errors.  This
    // method expands the bound for L so that it is guaranteed to contain the
    // bounds of any subregion S.
    //
    // More precisely, if L is a loop that does not contain either pole, and S
    // is a loop such that L.Contains(S), then
    //
    //   ExpandForSubregions(RectBound(L)).Contains(RectBound(S)).
    static S2LatLngRect ExpandForSubregions(S2LatLngRect const& bound);

    // Return the maximum error in GetBound() provided that the result does
    // not include either pole.  It is only to be used for testing purposes
    // (e.g., by passing it to S2LatLngRect::ApproxEquals).
    static S2LatLng MaxErrorForTests();

   private:
    S2Point a_;             // The previous vertex in the chain.
    S2LatLng a_latlng_;     // The corresponding latitude-longitude.
    S2LatLngRect bound_;    // The current bounding rectangle.

    RectBounder(RectBounder const&) = delete;
    void operator=(RectBounder const&) = delete;
  };

  // The purpose of this class is to find edges that intersect a given
  // longitude interval.  It can be used as an efficient rejection test when
  // attempting to find edges that intersect a given region.  It accepts a
  // vertex chain v0, v1, v2, ...  and returns a boolean value indicating
  // whether each edge intersects the specified longitude interval.
  class LongitudePruner {
   public:
    // 'interval' is the longitude interval to be tested against, and
    // 'v0' is the first vertex of edge chain.
    LongitudePruner(S1Interval const& interval, S2Point const& v0);

    // Returns true if the edge (v0, v1) intersects the given longitude
    // interval, and then saves 'v1' to be used as the next 'v0'.
    inline bool Intersects(S2Point const& v1);

   private:
    S1Interval interval_;    // The interval to be tested against.
    double lng0_;            // The longitude of the next v0.

    LongitudePruner(LongitudePruner const&) = delete;
    void operator=(LongitudePruner const&) = delete;
  };

  // Return true if edge AB crosses CD at a point that is interior
  // to both edges.  Properties:
  //
  //  (1) SimpleCrossing(b,a,c,d) == SimpleCrossing(a,b,c,d)
  //  (2) SimpleCrossing(c,d,a,b) == SimpleCrossing(a,b,c,d)
  static bool SimpleCrossing(S2Point const& a, S2Point const& b,
                             S2Point const& c, S2Point const& d);

  // Like SimpleCrossing, except that points that lie exactly on a line are
  // arbitrarily classified as being on one side or the other (according to
  // the rules of s2pred::Sign).  It returns +1 if there is a crossing, -1
  // if there is no crossing, and 0 if any two vertices from different edges
  // are the same.  Returns 0 or -1 if either edge is degenerate.
  // Properties of CrossingSign:
  //
  //  (1) CrossingSign(b,a,c,d) == CrossingSign(a,b,c,d)
  //  (2) CrossingSign(c,d,a,b) == CrossingSign(a,b,c,d)
  //  (3) CrossingSign(a,b,c,d) == 0 if a==c, a==d, b==c, b==d
  //  (3) CrossingSign(a,b,c,d) <= 0 if a==b or c==d
  //
  // Note that if you want to check an edge against a *chain* of other
  // edges, it is much more efficient to use an EdgeCrosser (above).
  static int CrossingSign(S2Point const& a, S2Point const& b, S2Point const& c,
                          S2Point const& d);

  // Given two edges AB and CD where at least two vertices are identical
  // (i.e. CrossingSign(a,b,c,d) == 0), this function defines whether the
  // two edges "cross" in a such a way that point-in-polygon containment tests
  // can be implemented by counting the number of edge crossings.  The basic
  // rule is that a "crossing" occurs if AB is encountered after CD during a
  // CCW sweep around the shared vertex starting from a fixed reference point.
  //
  // Note that according to this rule, if AB crosses CD then in general CD
  // does not cross AB.  However, this leads to the correct result when
  // counting polygon edge crossings.  For example, suppose that A,B,C are
  // three consecutive vertices of a CCW polygon.  If we now consider the edge
  // crossings of a segment BP as P sweeps around B, the crossing number
  // changes parity exactly when BP crosses BA or BC.
  //
  // Useful properties of VertexCrossing (VC):
  //
  //  (1) VC(a,a,c,d) == VC(a,b,c,c) == false
  //  (2) VC(a,b,a,b) == VC(a,b,b,a) == true
  //  (3) VC(a,b,c,d) == VC(a,b,d,c) == VC(b,a,c,d) == VC(b,a,d,c)
  //  (3) If exactly one of a,b equals one of c,d, then exactly one of
  //      VC(a,b,c,d) and VC(c,d,a,b) is true
  //
  // It is an error to call this method with 4 distinct vertices.
  static bool VertexCrossing(S2Point const& a, S2Point const& b,
                             S2Point const& c, S2Point const& d);

  // A convenience function that calls CrossingSign() to handle cases
  // where all four vertices are distinct, and VertexCrossing() to handle
  // cases where two or more vertices are the same.  This defines a crossing
  // function such that point-in-polygon containment tests can be implemented
  // by simply counting edge crossings.
  static bool EdgeOrVertexCrossing(S2Point const& a, S2Point const& b,
                                   S2Point const& c, S2Point const& d);

  // Given two edges AB and CD such that CrossingSign(A, B, C, D) > 0, return
  // their intersection point.  Useful properties of GetIntersection (GI):
  //
  //  (1) GI(b,a,c,d) == GI(a,b,d,c) == GI(a,b,c,d)
  //  (2) GI(c,d,a,b) == GI(a,b,c,d)
  //
  // The returned intersection point X is guaranteed to be very close to the
  // true intersection point of AB and CD, even if the edges intersect at a
  // very small angle.  See "kIntersectionError" below for details.
  static S2Point GetIntersection(S2Point const& a, S2Point const& b,
                                 S2Point const& c, S2Point const& d);

  // kIntersectionError is an upper bound on the distance from the intersection
  // point returned by GetIntersection() to the true intersection point.
  static S1Angle const kIntersectionError;

  // When using S2PolygonBuilder with computed intersection points, the
  // vertex_merge_radius() should be at least this large in order to avoid
  // incorrect output.
  static S1Angle const kIntersectionMergeRadius;  // 2 * kIntersectionError

  // Given a point X and an edge AB, return the distance ratio AX / (AX + BX).
  // If X happens to be on the line segment AB, this is the fraction "t" such
  // that X == Interpolate(t, A, B).  Requires that A and B are distinct.
  static double GetDistanceFraction(S2Point const& x,
                                    S2Point const& a, S2Point const& b);

  // Return the point X along the line segment AB whose distance from A is the
  // given fraction "t" of the distance AB.  Does NOT require that "t" be
  // between 0 and 1.  Note that all distances are measured on the surface of
  // the sphere, so this is more complicated than just computing (1-t)*a + t*b
  // and normalizing the result.
  static S2Point Interpolate(double t, S2Point const& a, S2Point const& b);

  // Like Interpolate(), except that the parameter "ax" represents the desired
  // distance from A to the result X rather than a fraction between 0 and 1.
  static S2Point InterpolateAtDistance(S1Angle ax,
                                       S2Point const& a, S2Point const& b);

  // Return the minimum distance from X to any point on the edge AB.  All
  // arguments should be unit length.  The result is very accurate for small
  // distances but may have some numerical error if the distance is large
  // (approximately Pi/2 or greater).  The case A == B is handled correctly.
  //
  // If you want to compare a distance against a fixed threshold, e.g.
  //    if (S2EdgeUtil::GetDistance(x, a, b) < limit)
  // then it is significantly faster to use UpdateMinDistance() below.
  static S1Angle GetDistance(S2Point const& x,
                             S2Point const& a, S2Point const& b);

  // Returns true if the distance from X to the edge AB is less than "limit".
  // This method is significantly faster than GetDistance().  If you want to
  // compare against a fixed S1Angle, you should convert it to an S1ChordAngle
  // once and save the value, since this step is relatively expensive.
  static bool IsDistanceLess(S2Point const& x,
                             S2Point const& a, S2Point const& b,
                             S1ChordAngle limit);

  // If the distance from X to the edge AB is less then "min_dist", this
  // method updates "min_dist" and returns true.  Otherwise it returns false.
  // The case A == B is handled correctly.
  //
  // Use this method when you want to compute many distances and keep track of
  // the minimum.  It is significantly faster than using GetDistance(),
  // because (1) using S1ChordAngle is much faster than S1Angle, and (2) it
  // can save a lot of work by not actually computing the distance when it is
  // obviously larger than the current minimum.
  static bool UpdateMinDistance(S2Point const& x,
                                S2Point const& a, S2Point const& b,
                                S1ChordAngle* min_dist);

  // Returns the maximum error in the result of UpdateMinDistance (and
  // associated functions such as UpdateMinInteriorDistance, IsDistanceLess,
  // etc), assuming that all input points are normalized to within the bounds
  // guaranteed by S2Point::Normalize().  The error can be added or subtracted
  // from an S1ChordAngle "x" using x.PlusError(error).
  static double GetUpdateMinDistanceMaxError(S1ChordAngle dist);

  // Returns true if the minimum distance from X to the edge AB is attained at
  // an interior point of AB (i.e., not an endpoint), and that distance is
  // less than "limit".
  static bool IsInteriorDistanceLess(S2Point const& x,
                                     S2Point const& a, S2Point const& b,
                                     S1ChordAngle limit);

  // If the minimum distance from X to AB is attained at an interior point of
  // AB (i.e., not an endpoint), and that distance is less than "min_dist",
  // then update "min_dist" and return true.  Otherwise return false.
  static bool UpdateMinInteriorDistance(S2Point const& x,
                                        S2Point const& a, S2Point const& b,
                                        S1ChordAngle* min_dist);

  // Return the point along the edge AB that is closest to the point X.
  // The fractional distance of this point along the edge AB can be obtained
  // using GetDistanceFraction() above.  Requires that all vectors have
  // unit length.
  static S2Point GetClosestPoint(S2Point const& x,
                                 S2Point const& a, S2Point const& b);

  // A slightly more efficient version of GetClosestPoint() where the cross
  // product of the two endpoints has been precomputed.  The cross product
  // does not need to be normalized, but should be computed using
  // S2::RobustCrossProd() for the most accurate results.  Requires that
  // x, a, and b have unit length.
  static S2Point GetClosestPoint(S2Point const& x,
                                 S2Point const& a, S2Point const& b,
                                 Vector3_d const& a_cross_b);

  /////////////////////////////////////////////////////////////////////
  ///////////////     Methods for pairs of edges      /////////////////

  // Like UpdateMinDistance(), but computes the minimum distance between the
  // given pair of edges.  (If the two edges cross, the distance is zero.)
  // The cases a0 == a1 and b0 == b1 are handled correctly.
  static bool UpdateEdgePairMinDistance(S2Point const& a0, S2Point const& a1,
                                        S2Point const& b0, S2Point const& b1,
                                        S1ChordAngle* min_dist);

  // Return the pair of points (a, b) that achieves the minimum distance
  // between edges a0a1 and b0b1, where "a" is a point on a0a1 and "b" is a
  // point on b0b1.  If the two edges intersect, "a" and "b" are both equal to
  // the intersection point.  Handles a0 == a1 and b0 == b1 correctly.
  static std::pair<S2Point, S2Point> GetEdgePairClosestPoints(
      S2Point const& a0, S2Point const& a1,
      S2Point const& b0, S2Point const& b1);

  // Return true if every point on edge B=b0b1 is no further than "tolerance"
  // from some point on edge A=a0a1.
  // Requires that tolerance is less than 90 degrees.
  static bool IsEdgeBNearEdgeA(S2Point const& a0, S2Point const& a1,
                               S2Point const& b0, S2Point const& b1,
                               S1Angle tolerance);

  // Given an edge chain (x0, x1, x2), the wedge at x1 is the region to the
  // left of the edges.  More precisely, it is the set of all rays from x1x0
  // (inclusive) to x1x2 (exclusive) in the *clockwise* direction.
  //
  // The following functions compare two *non-empty* wedges that share the
  // same middle vertex: A=(a0, ab1, a2) and B=(b0, ab1, b2).

  // Detailed relation from one wedge A to another wedge B.
  enum WedgeRelation {
    WEDGE_EQUALS,                 // A and B are equal.
    WEDGE_PROPERLY_CONTAINS,      // A is a strict superset of B.
    WEDGE_IS_PROPERLY_CONTAINED,  // A is a strict subset of B.
    WEDGE_PROPERLY_OVERLAPS,      // A-B, B-A, and A intersect B are non-empty.
    WEDGE_IS_DISJOINT,            // A and B are disjoint.
  };

  // Return the relation from wedge A to B.
  // REQUIRES: A and B are non-empty.
  static WedgeRelation GetWedgeRelation(
      S2Point const& a0, S2Point const& ab1, S2Point const& a2,
      S2Point const& b0, S2Point const& b2);

  // Returns true if wedge A contains wedge B.  Equivalent to but faster than
  // GetWedgeRelation() == WEDGE_PROPERLY_CONTAINS || WEDGE_EQUALS.
  // REQUIRES: A and B are non-empty.
  static bool WedgeContains(S2Point const& a0, S2Point const& ab1,
                            S2Point const& a2, S2Point const& b0,
                            S2Point const& b2);

  // Returns true if wedge A intersects wedge B.  Equivalent to but faster
  // than GetWedgeRelation() != WEDGE_IS_DISJOINT.
  // REQUIRES: A and B are non-empty.
  static bool WedgeIntersects(S2Point const& a0, S2Point const& ab1,
                              S2Point const& a2, S2Point const& b0,
                              S2Point const& b2);

  // FaceSegment represents an edge AB clipped to an S2 cube face.  It is
  // represented by a face index and a pair of (u,v) coordinates.
  struct FaceSegment {
    int face;
    R2Point a, b;
  };
  using FaceSegmentVector = absl::InlinedVector<FaceSegment, 6>;

  // Subdivide the given edge AB at every point where it crosses the boundary
  // between two S2 cube faces and return the corresponding FaceSegments.  The
  // segments are returned in order from A toward B.  The input points must be
  // unit length.
  //
  // This method guarantees that the returned segments form a continuous path
  // from A to B, and that all vertices are within kFaceClipErrorUVDist of the
  // line AB.  All vertices lie within the [-1,1]x[-1,1] cube face rectangles.
  // The results are consistent with s2pred::Sign(), i.e. the edge is
  // well-defined even its endpoints are antipodal.  [TODO(ericv): Extend the
  // implementation of S2::RobustCrossProd so that this statement is true.]
  static void GetFaceSegments(S2Point const& a, S2Point const& b,
                              FaceSegmentVector* segments);

  // Given an edge AB and a face, return the (u,v) coordinates for the portion
  // of AB that intersects that face.  This method guarantees that the clipped
  // vertices lie within the [-1,1]x[-1,1] cube face rectangle and are within
  // kFaceClipErrorUVDist of the line AB, but the results may differ from
  // those produced by GetFaceSegments.  Returns false if AB does not
  // intersect the given face.
  static bool ClipToFace(S2Point const& a, S2Point const& b, int face,
                         R2Point* a_uv, R2Point* b_uv);

  // Like ClipToFace, but rather than clipping to the square [-1,1]x[-1,1]
  // in (u,v) space, this method clips to [-R,R]x[-R,R] where R=(1+padding).
  static bool ClipToPaddedFace(S2Point const& a, S2Point const& b, int face,
                               double padding, R2Point* a_uv, R2Point* b_uv);

  // The maximum error in the vertices returned by GetFaceSegments and
  // ClipToFace (compared to an exact calculation):
  //
  //  - kFaceClipErrorRadians is the maximum angle between a returned vertex
  //    and the nearest point on the exact edge AB.  It is equal to the
  //    maximum directional error in S2::RobustCrossProd, plus the error when
  //    projecting points onto a cube face.
  //
  //  - kFaceClipErrorDist is the same angle expressed as a maximum distance
  //    in (u,v)-space.  In other words, a returned vertex is at most this far
  //    from the exact edge AB projected into (u,v)-space.

  //  - kFaceClipErrorUVCoord is the same angle expressed as the maximum error
  //    in an individual u- or v-coordinate.  In other words, for each
  //    returned vertex there is a point on the exact edge AB whose u- and
  //    v-coordinates differ from the vertex by at most this amount.

  static double const kFaceClipErrorRadians;
  static double const kFaceClipErrorUVDist;
  static double const kFaceClipErrorUVCoord;

  // Return true if the edge AB intersects the given (closed) rectangle to
  // within the error bound below.
  static bool IntersectsRect(R2Point const& a, R2Point const& b,
                             R2Rect const& rect);

  // The maximum error in IntersectRect.  If some point of AB is inside the
  // rectangle by at least this distance, the result is guaranteed to be true;
  // if all points of AB are outside the rectangle by at least this distance,
  // the result is guaranteed to be false.  This bound assumes that "rect" is
  // a subset of the rectangle [-1,1]x[-1,1] or extends slightly outside it
  // (e.g., by 1e-10 or less).
  static double const kIntersectsRectErrorUVDist;

  // Given an edge AB, return the portion of AB that is contained by the given
  // rectangle "clip".  Returns false if there is no intersection.
  static bool ClipEdge(R2Point const& a, R2Point const& b, R2Rect const& clip,
                       R2Point* a_clipped, R2Point* b_clipped);

  // Given an edge AB and a rectangle "clip", return the bounding rectangle of
  // the portion of AB intersected by "clip".  The resulting bound may be
  // empty.  This is a convenience function built on top of ClipEdgeBound.
  static R2Rect GetClippedEdgeBound(R2Point const& a, R2Point const& b,
                                    R2Rect const& clip);

  // This function can be used to clip an edge AB to sequence of rectangles
  // efficiently.  It represents the clipped edges by their bounding boxes
  // rather than as a pair of endpoints.  Specifically, let A'B' be some
  // portion of an edge AB, and let "bound" be a tight bound of A'B'.  This
  // function updates "bound" (in place) to be a tight bound of A'B'
  // intersected with a given rectangle "clip".  If A'B' does not intersect
  // "clip", returns false and does not necessarily update "bound".
  //
  // REQUIRES: "bound" is a tight bounding rectangle for some portion of AB.
  // (This condition is automatically satisfied if you start with the bounding
  // box of AB and clip to a sequence of rectangles, stopping when the method
  // returns false.)
  static bool ClipEdgeBound(R2Point const& a, R2Point const& b,
                            R2Rect const& clip, R2Rect* bound);

  // The maximum error in the vertices generated by ClipEdge and the bounds
  // generated by ClipEdgeBound (compared to an exact calculation):
  //
  //  - kEdgeClipErrorUVCoord is the maximum error in a u- or v-coordinate
  //    compared to the exact result, assuming that the points A and B are in
  //    the rectangle [-1,1]x[1,1] or slightly outside it (by 1e-10 or less).
  //
  //  - kEdgeClipErrorUVDist is the maximum distance from a clipped point to
  //    the corresponding exact result.  It is equal to the error in a single
  //    coordinate because at most one coordinate is subject to error.

  static double const kEdgeClipErrorUVCoord;
  static double const kEdgeClipErrorUVDist;

  // Given a value x that is some linear combination of a and b, return the
  // value x1 that is the same linear combination of a1 and b1.  This function
  // makes the following guarantees:
  //  - If x == a, then x1 = a1 (exactly).
  //  - If x == b, then x1 = b1 (exactly).
  //  - If a <= x <= b, then a1 <= x1 <= b1 (even if a1 == b1).
  // REQUIRES: a != b
  inline static double InterpolateDouble(double x, double a, double b,
                                         double a1, double b1);

  // Convert a point from the ExactFloat representation to the closest
  // double-precision value.  The result is normalized to be unit length.
  static S2Point S2PointFromExact(Vector3<ExactFloat> const& x);

 private:
  friend class GetIntersectionStats;
  friend class S2EdgeUtilTesting;

  // This is an *internal only* method declared here for testing purposes.
  static S2Point GetIntersectionExact(S2Point const& a0, S2Point const& a1,
                                      S2Point const& b0, S2Point const& b1);

  // The maximum error in the method above.
  static S1Angle const kIntersectionExactError;

  // The following field is used exclusively by S2EdgeUtilTesting in order to
  // measure how often each intersection method is used by GetIntersection().
  // If non-nullptr, then it points to an array of integers indexed by an
  // IntersectionMethod enum value.  Each call to GetIntersection() increments
  // the array entry corresponding to the intersection method that was used.
  static int* intersection_method_tally_;

  // The list of intersection methods implemented by GetIntersection().
  enum class IntersectionMethod {
    SIMPLE,
    SIMPLE_LD,
    STABLE,
    STABLE_LD,
    EXACT,
    NUM_METHODS
  };
  static char const* GetIntersectionMethodName(IntersectionMethod method);

  // Contains only static members.
  S2EdgeUtil() = delete;
  S2EdgeUtil(S2EdgeUtil const&) = delete;
  void operator=(S2EdgeUtil const&) = delete;
};


//////////////////   Implementation details follow   ////////////////////


inline S2EdgeUtil::EdgeCrosser::EdgeCrosser(S2Point const* a, S2Point const* b)
    : a_(a), b_(b), a_cross_b_(a_->CrossProd(*b_)), have_tangents_(false),
      c_(nullptr) {
  DCHECK(S2::IsUnitLength(*a));
  DCHECK(S2::IsUnitLength(*b));
}

inline void S2EdgeUtil::EdgeCrosser::Init(S2Point const* a, S2Point const* b) {
  a_ = a;
  b_ = b;
  a_cross_b_ = a->CrossProd(*b_);
  have_tangents_ = false;
  c_ = nullptr;
}

inline int S2EdgeUtil::EdgeCrosser::CrossingSign(S2Point const* c,
                                                 S2Point const* d) {
  if (c != c_) RestartAt(c);
  return CrossingSign(d);
}

inline bool S2EdgeUtil::EdgeCrosser::EdgeOrVertexCrossing(S2Point const* c,
                                                          S2Point const* d) {
  if (c != c_) RestartAt(c);
  return EdgeOrVertexCrossing(d);
}

inline S2EdgeUtil::EdgeCrosser::EdgeCrosser(
    S2Point const* a, S2Point const* b, S2Point const* c)
    : a_(a), b_(b), a_cross_b_(a_->CrossProd(*b_)), have_tangents_(false) {
  DCHECK(S2::IsUnitLength(*a));
  DCHECK(S2::IsUnitLength(*b));
  RestartAt(c);
}

inline void S2EdgeUtil::EdgeCrosser::RestartAt(S2Point const* c) {
  DCHECK(S2::IsUnitLength(*c));
  c_ = c;
  acb_ = -s2pred::TriageSign(*a_, *b_, *c_, a_cross_b_);
}

inline int S2EdgeUtil::EdgeCrosser::CrossingSign(S2Point const* d) {
  DCHECK(S2::IsUnitLength(*d));
  // For there to be an edge crossing, the triangles ACB, CBD, BDA, DAC must
  // all be oriented the same way (CW or CCW).  We keep the orientation of ACB
  // as part of our state.  When each new point D arrives, we compute the
  // orientation of BDA and check whether it matches ACB.  This checks whether
  // the points C and D are on opposite sides of the great circle through AB.

  // Recall that TriageSign is invariant with respect to rotating its
  // arguments, i.e. ABD has the same orientation as BDA.
  int bda = s2pred::TriageSign(*a_, *b_, *d, a_cross_b_);
  if (acb_ == -bda && bda != 0) {
    // The most common case -- triangles have opposite orientations.  Save the
    // current vertex D as the next vertex C, and also save the orientation of
    // the new triangle ACB (which is opposite to the current triangle BDA).
    c_ = d;
    acb_ = -bda;
    return -1;
  }
  bda_ = bda;
  return CrossingSignInternal(d);
}

inline bool S2EdgeUtil::EdgeCrosser::EdgeOrVertexCrossing(S2Point const* d) {
  // We need to copy c_ since it is clobbered by CrossingSign().
  S2Point const* c = c_;
  int crossing = CrossingSign(d);
  if (crossing < 0) return false;
  if (crossing > 0) return true;
  return VertexCrossing(*a_, *b_, *c, *d);
}

inline S2EdgeUtil::CopyingEdgeCrosser::CopyingEdgeCrosser(S2Point const& a,
                                                          S2Point const& b)
    : a_(a), b_(b), c_(S2Point()), crosser_(&a_, &b_) {
}

inline void S2EdgeUtil::CopyingEdgeCrosser::Init(S2Point const& a,
                                                 S2Point const& b) {
  a_ = a;
  b_ = b;
  c_ = S2Point();
  crosser_.Init(&a_, &b_);
}

inline int S2EdgeUtil::CopyingEdgeCrosser::CrossingSign(S2Point const& c,
                                                        S2Point const& d) {
  if (c != c_) RestartAt(c);
  return CrossingSign(d);
}

inline bool S2EdgeUtil::CopyingEdgeCrosser::EdgeOrVertexCrossing(
    S2Point const& c, S2Point const& d) {
  if (c != c_) RestartAt(c);
  return EdgeOrVertexCrossing(d);
}

inline S2EdgeUtil::CopyingEdgeCrosser::CopyingEdgeCrosser(
    S2Point const& a, S2Point const& b, S2Point const& c)
    : a_(a), b_(b), c_(c), crosser_(&a_, &b_, &c) {
}

inline void S2EdgeUtil::CopyingEdgeCrosser::RestartAt(S2Point const& c) {
  c_ = c;
  crosser_.RestartAt(&c_);
}

inline int S2EdgeUtil::CopyingEdgeCrosser::CrossingSign(S2Point const& d) {
  int result = crosser_.CrossingSign(&d);
  c_ = d;
  crosser_.set_c(&c_);
  return result;
}

inline bool S2EdgeUtil::CopyingEdgeCrosser::EdgeOrVertexCrossing(
    S2Point const& d) {
  bool result = crosser_.EdgeOrVertexCrossing(&d);
  c_ = d;
  crosser_.set_c(&c_);
  return result;
}

inline bool S2EdgeUtil::LongitudePruner::Intersects(S2Point const& v1) {
  double lng1 = S2LatLng::Longitude(v1).radians();
  bool result = interval_.Intersects(S1Interval::FromPointPair(lng0_, lng1));
  lng0_ = lng1;
  return result;
}

inline bool S2EdgeUtil::IsDistanceLess(S2Point const& x,
                                       S2Point const& a, S2Point const& b,
                                       S1ChordAngle limit) {
  return UpdateMinDistance(x, a, b, &limit);
}

inline bool S2EdgeUtil::IsInteriorDistanceLess(
    S2Point const& x, S2Point const& a, S2Point const& b,
    S1ChordAngle limit) {
  return UpdateMinInteriorDistance(x, a, b, &limit);
}

inline bool S2EdgeUtil::ClipToFace(S2Point const& a, S2Point const& b, int face,
                                   R2Point* a_uv, R2Point* b_uv) {
  return ClipToPaddedFace(a, b, face, 0.0, a_uv, b_uv);
}

inline double S2EdgeUtil::InterpolateDouble(double x, double a, double b,
                                            double a1, double b1) {
  DCHECK_NE(a, b);
  // To get results that are accurate near both A and B, we interpolate
  // starting from the closer of the two points.
  if (std::fabs(a - x) <= std::fabs(b - x)) {
    return a1 + (b1 - a1) * (x - a) / (b - a);
  } else {
    return b1 + (a1 - b1) * (x - b) / (a - b);
  }
}

#endif  // S2_S2EDGEUTIL_H_
