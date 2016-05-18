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

#ifndef S2_S2_H_
#define S2_S2_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "s2/base/integral_types.h"
#include "s2/base/macros.h"
#include "s2/base/port.h"
#include "s2/fpcontractoff.h"
#include "s2/r2.h"
#include "s2/util/math/mathutil.h"  // IWYU pragma: export
#include "s2/util/math/matrix3x3.h"
#include "s2/util/math/vector.h"  // IWYU pragma: export
#include "s2/util/math/vector3_hash.h"  // IWYU pragma: export

class S1Angle;

// This flag enables extra validity checking throughout the S2 code.  It is
// turned on by default in debug-mode builds.  You should turn it off in tests
// that want to create invalid geometry:
//
// TEST(MyClass, InvalidGeometry) {
//   FLAGS_s2debug = false;  // Automatically restored between tests
//   ...
// }
DECLARE_bool(s2debug);

// Alternatively, some classes have methods that allow the --s2debug flag to
// be disabled for specific objects (e.g., see S2Polygon).
enum class S2Debug : uint8 {
  ALLOW,    // Validity checks are controlled by --s2debug
  DISABLE   // No validity checks even when --s2debug is true
};

// An S2Point represents a point on the unit sphere as a 3D vector.  Usually
// points are normalized to be unit length, but some methods do not require
// this.  See util/math/vector.h for the methods available.  Among other
// things, there are overloaded operators that make it convenient to write
// arithmetic expressions (e.g. (1-x)*p1 + x*p2).
using S2Point = Vector3_d;

using S2PointHash = GoodFastHash<S2Point>;

// The S2 class is simply a namespace for constants and static utility
// functions related to spherical geometry, such as area calculations and edge
// intersection tests.  The name "S2" is derived from the mathematical symbol
// for the two-dimensional unit sphere (note that the "2" refers to the
// dimension of the surface, not the space it is embedded in).
//
// This class also defines a framework for decomposing the unit sphere into a
// hierarchy of "cells".  Each cell is a quadrilateral bounded by four
// geodesics.  The top level of the hierarchy is obtained by projecting the
// six faces of a cube onto the unit sphere, and lower levels are obtained by
// subdividing each cell into four children recursively.
//
// This class specifies the details of how the cube faces are projected onto
// the unit sphere.  This includes getting the face ordering and orientation
// correct so that sequentially increasing cell ids follow a continuous
// space-filling curve over the entire sphere, and defining the
// transformation from cell-space to cube-space in order to make the cells
// more uniform in size.
//
// This file also contains documentation of the various coordinate systems
// and conventions used.
//
class S2 {
 public:
  // Return a unique "origin" on the sphere for operations that need a fixed
  // reference point.  In particular, this is the "point at infinity" used for
  // point-in-polygon testing (by counting the number of edge crossings).
  //
  // It should *not* be a point that is commonly used in edge tests in order
  // to avoid triggering code to handle degenerate cases.  (This rules out the
  // north and south poles.)  It should also not be on the boundary of any
  // low-level S2Cell for the same reason.
  inline static S2Point Origin();

  // Return true if the given point is approximately unit length
  // (this is mainly useful for assertions).
  static bool IsUnitLength(S2Point const& p);

  // Return a unit-length vector that is orthogonal to "a".  Satisfies
  // Ortho(-a) = -Ortho(a) for all a.
  static S2Point Ortho(S2Point const& a);

  // Given a point "z" on the unit sphere, extend this into a right-handed
  // coordinate frame of unit-length column vectors m = (x,y,z).  Note that
  // the vectors (x,y) are an orthonormal frame for the tangent space at "z",
  // while "z" itself is an orthonormal frame for the normal space at "z".
  static Matrix3x3_d GetFrame(S2Point const& z);
  static void GetFrame(S2Point const& z, Matrix3x3_d* m);

  // Given an orthonormal basis "m" of column vectors and a point "p", return
  // the coordinates of "p" with respect to the basis "m".  The resulting
  // point "q" satisfies the identity (m * q == p).
  static S2Point ToFrame(Matrix3x3_d const& m, S2Point const& p);

  // Given an orthonormal basis "m" of column vectors and a point "q" with
  // respect to that basis, return the equivalent point "p" with respect to
  // the standard axis-aligned basis.  The result satisfies (p == m * q).
  static S2Point FromFrame(Matrix3x3_d const& m, S2Point const& q);

  // Rotate the given point about the given axis by the given angle.  "p" and
  // "axis" must be unit length; "angle" has no restrictions (e.g., it can be
  // positive, negative, greater than 360 degrees, etc).
  static S2Point Rotate(S2Point const& p, S2Point const& axis, S1Angle angle);

  // Return true if two points are within the given distance of each other
  // (this is mainly useful for testing).
  static bool ApproxEquals(S2Point const& a, S2Point const& b,
                           double max_error = 1e-15);

  // Return true if two arrays of points have the same sequence of points
  // except for perturbations.  More precisely, corresponding pairs of points
  // must be separated by no more than "max_error".  (Mostly for testing
  // purposes.)
  static bool PointsApproxEqual(S2Point const* a, int num_a,
                                S2Point const* b, int num_b,
                                double max_error = 1e-15);

  // Return a vector "c" that is orthogonal to the given unit-length vectors
  // "a" and "b".  This function is similar to a.CrossProd(b) except that it
  // does a better job of ensuring orthogonality when "a" is nearly parallel
  // to "b", and it returns a non-zero result even when a == b or a == -b.
  //
  // It satisfies the following properties (RCP == RobustCrossProd):
  //
  //   (1) RCP(a,b) != 0 for all a, b
  //   (2) RCP(b,a) == -RCP(a,b) unless a == b or a == -b
  //   (3) RCP(-a,b) == -RCP(a,b) unless a == b or a == -b
  //   (4) RCP(a,-b) == -RCP(a,b) unless a == b or a == -b
  static Vector3_d RobustCrossProd(S2Point const& a, S2Point const& b);

  // Return true if the points A, B, C are strictly counterclockwise.  Return
  // false if the points are clockwise or collinear (i.e. if they are all
  // contained on some great circle).
  //
  // Due to numerical errors, situations may arise that are mathematically
  // impossible, e.g. ABC may be considered strictly CCW while BCA is not.
  // However, the implementation guarantees the following:
  //
  //   If SimpleCCW(a,b,c), then !SimpleCCW(c,b,a) for all a,b,c.
  static bool SimpleCCW(S2Point const& a, S2Point const& b, S2Point const& c);

  // Returns +1 if the points A, B, C are counterclockwise, -1 if the points
  // are clockwise, and 0 if any two points are the same.  This function is
  // essentially like taking the sign of the determinant of ABC, except that
  // it has additional logic to make sure that the above properties hold even
  // when the three points are coplanar, and to deal with the limitations of
  // floating-point arithmetic.
  //
  // Sign satisfies the following conditions:
  //
  //  (1) Sign(a,b,c) == 0 if and only if a == b, b == c, or c == a
  //  (2) Sign(b,c,a) == Sign(a,b,c) for all a,b,c
  //  (3) Sign(c,b,a) == -Sign(a,b,c) for all a,b,c
  //
  // In other words:
  //
  //  (1) The result is zero if and only if two points are the same.
  //  (2) Rotating the order of the arguments does not affect the result.
  //  (3) Exchanging any two arguments inverts the result.
  //
  // On the other hand, note that it is not true in general that
  // Sign(-a,b,c) == -Sign(a,b,c), or any similar identities
  // involving antipodal points.
  static int Sign(S2Point const& a, S2Point const& b, S2Point const& c);

  // A more efficient version of Sign that allows the precomputed
  // cross-product of A and B to be specified.  (Unlike the 3 argument
  // version this method is also inlined.)
  inline static int Sign(S2Point const& a, S2Point const& b, S2Point const& c,
                         Vector3_d const& a_cross_b);

  // This version of Sign returns +1 if the points are definitely CCW,
  // -1 if they are definitely CW, and 0 if two points are identical or the
  // result is uncertain.  Uncertain certain cases can be resolved, if
  // desired, by calling ExpensiveSign.
  //
  // The purpose of this method is to allow additional cheap tests to be done,
  // where possible, in order to avoid calling ExpensiveSign unnecessarily.
  inline static int TriageSign(S2Point const& a, S2Point const& b,
                               S2Point const& c, Vector3_d const& a_cross_b);

  // This function is invoked by Sign() if the sign of the determinant is
  // uncertain.  It always returns a non-zero result unless two of the input
  // points are the same.  It uses a combination of multiple-precision
  // arithmetic and symbolic perturbations to ensure that its results are
  // always self-consistent (cf. Simulation of Simplicity, Edelsbrunner and
  // Muecke).  The basic idea is to assign an infinitesimal symbolic
  // perturbation to every possible S2Point such that no three S2Points are
  // collinear and no four S2Points are coplanar.  These perturbations are so
  // small that they do not affect the sign of any determinant that was
  // non-zero before the perturbations.
  //
  // Unlike Sign(), this method does not require the input points to be
  // normalized.
  static int ExpensiveSign(S2Point const& a, S2Point const& b,
                           S2Point const& c);

  // Given 4 points on the unit sphere, return true if the edges OA, OB, and
  // OC are encountered in that order while sweeping CCW around the point O.
  // You can think of this as testing whether A <= B <= C with respect to the
  // CCW ordering around O that starts at A, or equivalently, whether B is
  // contained in the range of angles (inclusive) that starts at A and extends
  // CCW to C.  Properties:
  //
  //  (1) If OrderedCCW(a,b,c,o) && OrderedCCW(b,a,c,o), then a == b
  //  (2) If OrderedCCW(a,b,c,o) && OrderedCCW(a,c,b,o), then b == c
  //  (3) If OrderedCCW(a,b,c,o) && OrderedCCW(c,b,a,o), then a == b == c
  //  (4) If a == b or b == c, then OrderedCCW(a,b,c,o) is true
  //  (5) Otherwise if a == c, then OrderedCCW(a,b,c,o) is false
  static bool OrderedCCW(S2Point const& a, S2Point const& b, S2Point const& c,
                         S2Point const& o);

  // Return the interior angle at the vertex B in the triangle ABC.  The
  // return value is always in the range [0, Pi].  All points should be
  // normalized.  Ensures that Angle(a,b,c) == Angle(c,b,a) for all a,b,c.
  //
  // The angle is undefined if A or C is diametrically opposite from B, and
  // becomes numerically unstable as the length of edge AB or BC approaches
  // 180 degrees.
  static double Angle(S2Point const& a, S2Point const& b, S2Point const& c);

  // Return the exterior angle at vertex B in the triangle ABC.  The return
  // value is positive if ABC is counterclockwise and negative otherwise.  If
  // you imagine an ant walking from A to B to C, this is the angle that the
  // ant turns at vertex B (positive = left = CCW, negative = right = CW).
  // This quantity is also known as the "geodesic curvature" at B.
  //
  // Ensures that TurnAngle(a,b,c) == -TurnAngle(c,b,a) for all distinct
  // a,b,c. The result is undefined if (a == b || b == c), but is either
  // -Pi or Pi if (a == c).  All points should be normalized.
  static double TurnAngle(S2Point const& a, S2Point const& b,
                          S2Point const& c);

  // Return the area of triangle ABC.  This method combines two different
  // algorithms to get accurate results for both large and small triangles.
  // The maximum error is about 5e-15 (about 0.25 square meters on the Earth's
  // surface), the same as GirardArea() below, but unlike that method it is
  // also accurate for small triangles.  Example: when the true area is 100
  // square meters, Area() yields an error about 1 trillion times smaller than
  // GirardArea().
  //
  // All points should be unit length, and no two points should be antipodal.
  // The area is always positive.
  static double Area(S2Point const& a, S2Point const& b, S2Point const& c);

  // Return the area of the triangle computed using Girard's formula.  All
  // points should be unit length, and no two points should be antipodal.
  //
  // This method is about twice as fast as Area() but has poor relative
  // accuracy for small triangles.  The maximum error is about 5e-15 (about
  // 0.25 square meters on the Earth's surface) and the average error is about
  // 1e-15.  These bounds apply to triangles of any size, even as the maximum
  // edge length of the triangle approaches 180 degrees.  But note that for
  // such triangles, tiny perturbations of the input points can change the
  // true mathematical area dramatically.
  static double GirardArea(S2Point const& a, S2Point const& b,
                           S2Point const& c);

  // Like Area(), but returns a positive value for counterclockwise triangles
  // and a negative value otherwise.
  static double SignedArea(S2Point const& a, S2Point const& b,
                           S2Point const& c);

  // About centroids:
  // ----------------
  //
  // There are several notions of the "centroid" of a triangle.  First, there
  // is the planar centroid, which is simply the centroid of the ordinary
  // (non-spherical) triangle defined by the three vertices.  Second, there is
  // the surface centroid, which is defined as the intersection of the three
  // medians of the spherical triangle.  It is possible to show that this
  // point is simply the planar centroid projected to the surface of the
  // sphere.  Finally, there is the true centroid (mass centroid), which is
  // defined as the surface integral over the spherical triangle of (x,y,z)
  // divided by the triangle area.  This is the point that the triangle would
  // rotate around if it was spinning in empty space.
  //
  // The best centroid for most purposes is the true centroid.  Unlike the
  // planar and surface centroids, the true centroid behaves linearly as
  // regions are added or subtracted.  That is, if you split a triangle into
  // pieces and compute the average of their centroids (weighted by triangle
  // area), the result equals the centroid of the original triangle.  This is
  // not true of the other centroids.
  //
  // Also note that the surface centroid may be nowhere near the intuitive
  // "center" of a spherical triangle.  For example, consider the triangle
  // with vertices A=(1,eps,0), B=(0,0,1), C=(-1,eps,0) (a quarter-sphere).
  // The surface centroid of this triangle is at S=(0, 2*eps, 1), which is
  // within a distance of 2*eps of the vertex B.  Note that the median from A
  // (the segment connecting A to the midpoint of BC) passes through S, since
  // this is the shortest path connecting the two endpoints.  On the other
  // hand, the true centroid is at M=(0, 0.5, 0.5), which when projected onto
  // the surface is a much more reasonable interpretation of the "center" of
  // this triangle.

  // Return the centroid of the planar triangle ABC.  This can be normalized
  // to unit length to obtain the "surface centroid" of the corresponding
  // spherical triangle, i.e. the intersection of the three medians.  However,
  // note that for large spherical triangles the surface centroid may be
  // nowhere near the intuitive "center" (see example above).
  static S2Point PlanarCentroid(S2Point const& a, S2Point const& b,
                                S2Point const& c);

  // Returns the true centroid of the spherical triangle ABC multiplied by the
  // signed area of spherical triangle ABC.  The reasons for multiplying by
  // the signed area are (1) this is the quantity that needs to be summed to
  // compute the centroid of a union or difference of triangles, and (2) it's
  // actually easier to calculate this way.  All points must have unit length.
  static S2Point TrueCentroid(S2Point const& a, S2Point const& b,
                              S2Point const& c);

  ////////////////////////// S2Cell Decomposition /////////////////////////
  //
  // The following methods define the cube-to-sphere projection used by
  // the S2Cell decomposition.
  //
  // In the process of converting a latitude-longitude pair to a 64-bit cell
  // id, the following coordinate systems are used:
  //
  //  (id)
  //    An S2CellId is a 64-bit encoding of a face and a Hilbert curve position
  //    on that face.  The Hilbert curve position implicitly encodes both the
  //    position of a cell and its subdivision level (see s2cellid.h).
  //
  //  (face, i, j)
  //    Leaf-cell coordinates.  "i" and "j" are integers in the range
  //    [0,(2**30)-1] that identify a particular leaf cell on the given face.
  //    The (i, j) coordinate system is right-handed on each face, and the
  //    faces are oriented such that Hilbert curves connect continuously from
  //    one face to the next.
  //
  //  (face, s, t)
  //    Cell-space coordinates.  "s" and "t" are real numbers in the range
  //    [0,1] that identify a point on the given face.  For example, the point
  //    (s, t) = (0.5, 0.5) corresponds to the center of the top-level face
  //    cell.  This point is also a vertex of exactly four cells at each
  //    subdivision level greater than zero.
  //
  //  (face, si, ti)
  //    Discrete cell-space coordinates.  These are obtained by multiplying
  //    "s" and "t" by 2**31 and rounding to the nearest unsigned integer.
  //    Discrete coordinates lie in the range [0,2**31].  This coordinate
  //    system can represent the edge and center positions of all cells with
  //    no loss of precision (including non-leaf cells).  In binary, each
  //    coordinate of a level-k cell center ends with a 1 followed by
  //    (30 - k) 0s.  The coordinates of its edges end with (at least)
  //    (31 - k) 0s.
  //
  //  (face, u, v)
  //    Cube-space coordinates in the range [-1,1].  To make the cells at each
  //    level more uniform in size after they are projected onto the sphere,
  //    we apply a nonlinear transformation of the form u=f(s), v=f(t).
  //    The (u, v) coordinates after this transformation give the actual
  //    coordinates on the cube face (modulo some 90 degree rotations) before
  //    it is projected onto the unit sphere.
  //
  //  (face, u, v, w)
  //    Per-face coordinate frame.  This is an extension of the (face, u, v)
  //    cube-space coordinates that adds a third axis "w" in the direction of
  //    the face normal.  It is always a right-handed 3D coordinate system.
  //    Cube-space coordinates can be converted to this frame by setting w=1,
  //    while (u,v,w) coordinates can be projected onto the cube face by
  //    dividing by w, i.e. (face, u/w, v/w).
  //
  //  (x, y, z)
  //    Direction vector (S2Point).  Direction vectors are not necessarily unit
  //    length, and are often chosen to be points on the biunit cube
  //    [-1,+1]x[-1,+1]x[-1,+1].  They can be be normalized to obtain the
  //    corresponding point on the unit sphere.
  //
  //  (lat, lng)
  //    Latitude and longitude (S2LatLng).  Latitudes must be between -90 and
  //    90 degrees inclusive, and longitudes must be between -180 and 180
  //    degrees inclusive.
  //
  // Note that the (i, j), (s, t), (si, ti), and (u, v) coordinate systems are
  // right-handed on all six faces.

  // This is the number of levels needed to specify a leaf cell.  This
  // constant is defined here so that the S2::Metric class and the conversion
  // functions below can be implemented without including s2cellid.h.  Please
  // see s2cellid.h for other useful constants and conversion functions.
  static int const kMaxCellLevel = 30;

  // The maximum index of a valid leaf cell plus one.  The range of valid leaf
  // cell indices is [0..kLimitIJ-1].
  static int const kLimitIJ = 1 << kMaxCellLevel;  // == S2CellId::kMaxSize

  // The maximum value of an si- or ti-coordinate.  The range of valid (si,ti)
  // values is [0..kMaxSiTi].
  static unsigned int const kMaxSiTi = 1U << (kMaxCellLevel + 1);

  // Convert an s- or t-value to the corresponding u- or v-value.  This is
  // a non-linear transformation from [0,1] to [-1,1] that attempts to
  // make the cell sizes more uniform.
  inline static double STtoUV(double s);

  // The inverse of the STtoUV transformation.  Note that it is not always
  // true that UVtoST(STtoUV(x)) == x due to numerical errors.
  inline static double UVtoST(double u);

  // Convert the i- or j-index of a leaf cell to the minimum corresponding s-
  // or t-value contained by that cell.  The argument must be in the range
  // [0..2**30], i.e. up to one position beyond the normal range of valid leaf
  // cell indices.
  inline static double IJtoSTMin(int i);

  // Return the i- or j-index of the leaf cell containing the given
  // s- or t-value.  If the argument is outside the range spanned by valid
  // leaf cell indices, return the index of the closest valid leaf cell (i.e.,
  // return values are clamped to the range of valid leaf cell indices).
  inline static int STtoIJ(double s);

  // Convert an si- or ti-value to the corresponding s- or t-value.
  inline static double SiTitoST(unsigned int si);

  // Return the si- or ti-coordinate that is nearest to the given s- or
  // t-value.  The result may be outside the range of valid (si,ti)-values.
  inline static unsigned int STtoSiTi(double s);

  // Convert (face, u, v) coordinates to a direction vector (not
  // necessarily unit length).
  inline static S2Point FaceUVtoXYZ(int face, double u, double v);
  inline static S2Point FaceUVtoXYZ(int face, R2Point const& uv);

  // If the dot product of p with the given face normal is positive,
  // set the corresponding u and v values (which may lie outside the range
  // [-1,1]) and return true.  Otherwise return false.
  inline static bool FaceXYZtoUV(int face, S2Point const& p,
                                 double* pu, double* pv);
  inline static bool FaceXYZtoUV(int face, S2Point const& p, R2Point* puv);

  // Given a *valid* face for the given point p (meaning that dot product
  // of p with the face normal is positive), return the corresponding
  // u and v values (which may lie outside the range [-1,1]).
  inline static void ValidFaceXYZtoUV(int face, S2Point const& p,
                                      double* pu, double* pv);
  inline static void ValidFaceXYZtoUV(int face, S2Point const& p, R2Point* puv);

  // Transform the given point P to the (u,v,w) coordinate frame of the given
  // face (where the w-axis represents the face normal).
  static S2Point FaceXYZtoUVW(int face, S2Point const& p);

  // Return the face containing the given direction vector.  (For points on
  // the boundary between faces, the result is arbitrary but repeatable.)
  inline static int GetFace(S2Point const& p);

  // Convert a direction vector (not necessarily unit length) to
  // (face, u, v) coordinates.
  inline static int XYZtoFaceUV(S2Point const& p, double* pu, double* pv);
  inline static int XYZtoFaceUV(S2Point const& p, R2Point* puv);

  // Convert a direction vector (not necessarily unit length) to
  // (face, si, ti) coordinates and, if p is exactly equal to the center of a
  // cell, return the level of this cell (-1 otherwise).
  static int XYZtoFaceSiTi(S2Point const& p, int* face,
                           unsigned int* si, unsigned int* ti);

  // Convert (face, si, ti) coordinates to a direction vector (not necessarily
  // unit length).
  static S2Point FaceSiTitoXYZ(int face, unsigned int si, unsigned int ti);

  // Return the right-handed normal (not necessarily unit length) for an
  // edge in the direction of the positive v-axis at the given u-value on
  // the given face.  (This vector is perpendicular to the plane through
  // the sphere origin that contains the given edge.)
  inline static S2Point GetUNorm(int face, double u);

  // Return the right-handed normal (not necessarily unit length) for an
  // edge in the direction of the positive u-axis at the given v-value on
  // the given face.
  inline static S2Point GetVNorm(int face, double v);

  // Return the unit-length normal, u-axis, or v-axis for the given face.
  inline static S2Point GetNorm(int face);
  inline static S2Point GetUAxis(int face);
  inline static S2Point GetVAxis(int face);

  // Return the given axis of the given face (u=0, v=1, w=2).
  inline static S2Point GetUVWAxis(int face, int axis);

  // With respect to the (u,v,w) coordinate system of a given face, return the
  // face that lies in the given direction (negative=0, positive=1) of the
  // given axis (u=0, v=1, w=2).  For example, GetUVWFace(4, 0, 1) returns the
  // face that is adjacent to face 4 in the positive u-axis direction.
  inline static int GetUVWFace(int face, int axis, int direction);

  ////////////////////////////////////////////////////////////////////////
  // The canonical Hilbert traversal order looks like an inverted 'U':
  // the subcells are visited in the order (0,0), (0,1), (1,1), (1,0).
  // The following tables encode the traversal order for various
  // orientations of the Hilbert curve (axes swapped and/or directions
  // of the axes reversed).

  // Together these flags define a cell orientation.  If 'kSwapMask'
  // is true, then canonical traversal order is flipped around the
  // diagonal (i.e. i and j are swapped with each other).  If
  // 'kInvertMask' is true, then the traversal order is rotated by 180
  // degrees (i.e. the bits of i and j are inverted, or equivalently,
  // the axis directions are reversed).
  static int const kSwapMask = 0x01;
  static int const kInvertMask = 0x02;

  // kIJtoPos[orientation][ij] -> pos
  //
  // Given a cell orientation and the (i,j)-index of a subcell (0=(0,0),
  // 1=(0,1), 2=(1,0), 3=(1,1)), return the order in which this subcell is
  // visited by the Hilbert curve (a position in the range [0..3]).
  static int const kIJtoPos[4][4];

  // kPosToIJ[orientation][pos] -> ij
  //
  // Return the (i,j) index of the subcell at the given position 'pos' in the
  // Hilbert curve traversal order with the given orientation.  This is the
  // inverse of the previous table:
  //
  //   kPosToIJ[r][kIJtoPos[r][ij]] == ij
  static int const kPosToIJ[4][4];

  // kPosToOrientation[pos] -> orientation_modifier
  //
  // Return a modifier indicating how the orientation of the child subcell
  // with the given traversal position [0..3] is related to the orientation
  // of the parent cell.  The modifier should be XOR-ed with the parent
  // orientation to obtain the curve orientation in the child.
  static int const kPosToOrientation[4];

  ////////////////////////// S2Cell Metrics //////////////////////////////
  //
  // The following are various constants that describe the shapes and sizes of
  // cells.  They are useful for deciding which cell level to use in order to
  // satisfy a given condition (e.g. that cell vertices must be no further
  // than "x" apart).  All of the raw constants are differential quantities;
  // you can use the GetValue(level) method to compute the corresponding length
  // or area on the unit sphere for cells at a given level.  The minimum and
  // maximum bounds are valid for cells at all levels, but they may be
  // somewhat conservative for very large cells (e.g. face cells).

  // Defines a cell metric of the given dimension (1 == length, 2 == area).
  template <int dim> class Metric {
   public:
    explicit constexpr Metric(double deriv) : deriv_(deriv) {}

    // The "deriv" value of a metric is a derivative, and must be multiplied by
    // a length or area in (s,t)-space to get a useful value.
    double deriv() const { return deriv_; }

    // Return the value of a metric for cells at the given level. The value is
    // either a length or an area on the unit sphere, depending on the
    // particular metric.
    double GetValue(int level) const { return ldexp(deriv_, - dim * level); }

    // Return the level at which the metric has approximately the given
    // value.  For example, S2::kAvgEdge.GetClosestLevel(0.1) returns the
    // level at which the average cell edge length is approximately 0.1.
    // The return value is always a valid level.
    int GetClosestLevel(double value) const;

    // Return the minimum level such that the metric is at most the given
    // value, or S2CellId::kMaxLevel if there is no such level.  For example,
    // S2::kMaxDiag.GetMinLevel(0.1) returns the minimum level such that all
    // cell diagonal lengths are 0.1 or smaller.  The return value is always a
    // valid level.
    // TODO(ericv): Rename to GetCeilingLevel(), or GetLevelForMaxValue()
    int GetMinLevel(double value) const;

    // Return the maximum level such that the metric is at least the given
    // value, or zero if there is no such level.  For example,
    // S2::kMinWidth.GetMaxLevel(0.1) returns the maximum level such that all
    // cells have a minimum width of 0.1 or larger.  The return value is
    // always a valid level.
    // TODO(ericv): Rename to GetFloorLevel(), or GetLevelForMinValue()
    int GetMaxLevel(double value) const;

   private:
    double const deriv_;

    Metric(Metric const&) = delete;
    void operator=(Metric const&) = delete;
  };
  using LengthMetric = Metric<1>;
  using AreaMetric = Metric<2>;

  // Each cell is bounded by four planes passing through its four edges and
  // the center of the sphere.  These metrics relate to the angle between each
  // pair of opposite bounding planes, or equivalently, between the planes
  // corresponding to two different s-values or two different t-values.  For
  // example, the maximum angle between opposite bounding planes for a cell at
  // level k is kMaxAngleSpan.GetValue(k), and the average angle span for all
  // cells at level k is approximately kAvgAngleSpan.GetValue(k).
  static LengthMetric const kMinAngleSpan;
  static LengthMetric const kMaxAngleSpan;
  static LengthMetric const kAvgAngleSpan;

  // The width of geometric figure is defined as the distance between two
  // parallel bounding lines in a given direction.  For cells, the minimum
  // width is always attained between two opposite edges, and the maximum
  // width is attained between two opposite vertices.  However, for our
  // purposes we redefine the width of a cell as the perpendicular distance
  // between a pair of opposite edges.  A cell therefore has two widths, one
  // in each direction.  The minimum width according to this definition agrees
  // with the classic geometric one, but the maximum width is different.  (The
  // maximum geometric width corresponds to kMaxDiag defined below.)
  //
  // For a cell at level k, the distance between opposite edges is at least
  // kMinWidth.GetValue(k) and at most kMaxWidth.GetValue(k).  The average
  // width in both directions for all cells at level k is approximately
  // kAvgWidth.GetValue(k).
  //
  // The width is useful for bounding the minimum or maximum distance from a
  // point on one edge of a cell to the closest point on the opposite edge.
  // For example, this is useful when "growing" regions by a fixed distance.
  static LengthMetric const kMinWidth;
  static LengthMetric const kMaxWidth;
  static LengthMetric const kAvgWidth;

  // The minimum edge length of any cell at level k is at least
  // kMinEdge.GetValue(k), and the maximum is at most kMaxEdge.GetValue(k).
  // The average edge length is approximately kAvgEdge.GetValue(k).
  //
  // The edge length metrics can also be used to bound the minimum, maximum,
  // or average distance from the center of one cell to the center of one of
  // its edge neighbors.  In particular, it can be used to bound the distance
  // between adjacent cell centers along the space-filling Hilbert curve for
  // cells at any given level.
  static LengthMetric const kMinEdge;
  static LengthMetric const kMaxEdge;
  static LengthMetric const kAvgEdge;

  // The minimum diagonal length of any cell at level k is at least
  // kMinDiag.GetValue(k), and the maximum is at most kMaxDiag.GetValue(k).
  // The average diagonal length is approximately kAvgDiag.GetValue(k).
  //
  // The maximum diagonal also happens to be the maximum diameter of any cell,
  // and also the maximum geometric width (see the discussion above).  So for
  // example, the distance from an arbitrary point to the closest cell center
  // at a given level is at most half the maximum diagonal length.
  static LengthMetric const kMinDiag;
  static LengthMetric const kMaxDiag;
  static LengthMetric const kAvgDiag;

  // The minimum area of any cell at level k is at least kMinArea.GetValue(k),
  // and the maximum is at most kMaxArea.GetValue(k).  The average area of all
  // cells at level k is exactly kAvgArea.GetValue(k).
  static AreaMetric const kMinArea;
  static AreaMetric const kMaxArea;
  static AreaMetric const kAvgArea;

  // This is the maximum edge aspect ratio over all cells at any level, where
  // the edge aspect ratio of a cell is defined as the ratio of its longest
  // edge length to its shortest edge length.
  static double const kMaxEdgeAspect;

  // This is the maximum diagonal aspect ratio over all cells at any level,
  // where the diagonal aspect ratio of a cell is defined as the ratio of its
  // longest diagonal length to its shortest diagonal length.
  static double const kMaxDiagAspect;

 private:
  // Documented in the .cc file; exposed here for testing.
  static int StableSign(S2Point const& a, S2Point const& b, S2Point const& c);
  static int ExactSign(S2Point const& a, S2Point const& b, S2Point const& c);
  friend class StableSignTest;

  // The value below is the maximum error in computing the determinant
  // a.CrossProd(b).DotProd(c).
  static double const kMaxDetError;

  // The U,V,W axes for each face.
  static double const kFaceUVWAxes[6][3][3];

  // The precomputed neighbors of each face (see GetUVWFace).
  static int const kFaceUVWFaces[6][3][2];

  // Contains only static methods.
  S2() = delete;
  S2(S2 const&) = delete;
  void operator=(S2 const&) = delete;
};


//////////////////   Implementation details follow   ////////////////////


// Uncomment the following line for testing purposes only.  It greatly
// increases the number of degenerate cases that need to be handled using
// ExpensiveSign().
// #define S2_TEST_DEGENERACIES

inline S2Point S2::Origin() {
#ifdef S2_TEST_DEGENERACIES
  return S2Point(0, 0, 1);  // This makes polygon operations much slower.
#else
  // This point is about 66km from the north pole towards the East Siberian
  // Sea.  See the unittest for more details.  It is written here using
  // floating-point literals because the optimizer doesn't seem willing to
  // evaluate Normalize() at compile time.
  return S2Point(-0.0099994664350250197, 0.0025924542609324121,
                 0.99994664350250195);
#endif
}

inline int S2::TriageSign(S2Point const& a, S2Point const& b,
                          S2Point const& c, Vector3_d const& a_cross_b) {
  DCHECK(IsUnitLength(a));
  DCHECK(IsUnitLength(b));
  DCHECK(IsUnitLength(c));
  double det = a_cross_b.DotProd(c);

  // Double-check borderline cases in debug mode.
  DCHECK(!FLAGS_s2debug ||
         fabs(det) <= kMaxDetError ||
         fabs(det) >= 100 * kMaxDetError ||
         det * ExpensiveSign(a, b, c) > 0);

  if (det > kMaxDetError) return 1;
  if (det < -kMaxDetError) return -1;
  return 0;
}

inline int S2::Sign(S2Point const& a, S2Point const& b, S2Point const& c,
                    Vector3_d const& a_cross_b) {
  int sign = TriageSign(a, b, c, a_cross_b);
  if (sign == 0) sign = ExpensiveSign(a, b, c);
  return sign;
}

// We have implemented three different projections from cell-space (s,t) to
// cube-space (u,v): linear, quadratic, and tangent.  They have the following
// tradeoffs:
//
//   Linear - This is the fastest transformation, but also produces the least
//   uniform cell sizes.  Cell areas vary by a factor of about 5.2, with the
//   largest cells at the center of each face and the smallest cells in
//   the corners.
//
//   Tangent - Transforming the coordinates via atan() makes the cell sizes
//   more uniform.  The areas vary by a maximum ratio of 1.4 as opposed to a
//   maximum ratio of 5.2.  However, each call to atan() is about as expensive
//   as all of the other calculations combined when converting from points to
//   cell ids, i.e. it reduces performance by a factor of 3.
//
//   Quadratic - This is an approximation of the tangent projection that
//   is much faster and produces cells that are almost as uniform in size.
//   It is about 3 times faster than the tangent projection for converting
//   cell ids to points or vice versa.  Cell areas vary by a maximum ratio of
//   about 2.1.
//
// Here is a table comparing the cell uniformity using each projection.  "Area
// ratio" is the maximum ratio over all subdivision levels of the largest cell
// area to the smallest cell area at that level, "edge ratio" is the maximum
// ratio of the longest edge of any cell to the shortest edge of any cell at
// the same level, and "diag ratio" is the ratio of the longest diagonal of
// any cell to the shortest diagonal of any cell at the same level.  "ToPoint"
// and "FromPoint" are the times in microseconds required to convert cell ids
// to and from points (unit vectors) respectively.  "ToPointRaw" is the time
// to convert to a non-unit-length vector, which is all that is needed for
// some purposes.
//
//               Area    Edge    Diag   ToPointRaw  ToPoint  FromPoint
//              Ratio   Ratio   Ratio             (microseconds)
// -------------------------------------------------------------------
// Linear:      5.200   2.117   2.959      0.020     0.087     0.085
// Tangent:     1.414   1.414   1.704      0.237     0.299     0.258
// Quadratic:   2.082   1.802   1.932      0.033     0.096     0.108
//
// The worst-case cell aspect ratios are about the same with all three
// projections.  The maximum ratio of the longest edge to the shortest edge
// within the same cell is about 1.4 and the maximum ratio of the diagonals
// within the same cell is about 1.7.
//
// This data was produced using s2cell_unittest and s2cellid_unittest.

#define S2_LINEAR_PROJECTION    0
#define S2_TAN_PROJECTION       1
#define S2_QUADRATIC_PROJECTION 2

#define S2_PROJECTION S2_QUADRATIC_PROJECTION

#if S2_PROJECTION == S2_LINEAR_PROJECTION

inline double S2::STtoUV(double s) {
  return 2 * s - 1;
}

inline double S2::UVtoST(double u) {
  return 0.5 * (u + 1);
}

#elif S2_PROJECTION == S2_TAN_PROJECTION

inline double S2::STtoUV(double s) {
  // Unfortunately, tan(M_PI_4) is slightly less than 1.0.  This isn't due to
  // a flaw in the implementation of tan(), it's because the derivative of
  // tan(x) at x=pi/4 is 2, and it happens that the two adjacent floating
  // point numbers on either side of the infinite-precision value of pi/4 have
  // tangents that are slightly below and slightly above 1.0 when rounded to
  // the nearest double-precision result.

  s = tan(M_PI_2 * s - M_PI_4);
  return s + (1.0 / (GG_LONGLONG(1) << 53)) * s;
}

inline double S2::UVtoST(double u) {
  volatile double a = atan(u);
  return (2 * M_1_PI) * (a + M_PI_4);
}

#elif S2_PROJECTION == S2_QUADRATIC_PROJECTION

inline double S2::STtoUV(double s) {
  if (s >= 0.5) return (1/3.) * (4*s*s - 1);
  else          return (1/3.) * (1 - 4*(1-s)*(1-s));
}

inline double S2::UVtoST(double u) {
  if (u >= 0) return 0.5 * sqrt(1 + 3*u);
  else        return 1 - 0.5 * sqrt(1 - 3*u);
}

#else

#error Unknown value for S2_PROJECTION

#endif

inline double S2::IJtoSTMin(int i) {
  DCHECK(i >= 0 && i <= kLimitIJ);
  return (1.0 / kLimitIJ) * i;
}

inline int S2::STtoIJ(double s) {
  return std::max(0, std::min(kLimitIJ - 1,
                              MathUtil::FastIntRound(kLimitIJ * s - 0.5)));
}

inline double S2::SiTitoST(unsigned int si) {
  DCHECK(si >= 0 && si <= kMaxSiTi);
  return (1.0 / kMaxSiTi) * si;
}

inline unsigned int S2::STtoSiTi(double s) {
  // kMaxSiTi == 2^31, so the result doesn't fit in an int32 when s == 1.
  return static_cast<unsigned int>(MathUtil::FastInt64Round(s * kMaxSiTi));
}

inline S2Point S2::FaceUVtoXYZ(int face, double u, double v) {
  switch (face) {
    case 0:  return S2Point( 1,  u,  v);
    case 1:  return S2Point(-u,  1,  v);
    case 2:  return S2Point(-u, -v,  1);
    case 3:  return S2Point(-1, -v, -u);
    case 4:  return S2Point( v, -1, -u);
    default: return S2Point( v,  u, -1);
  }
}

inline S2Point S2::FaceUVtoXYZ(int face, R2Point const& uv) {
  return FaceUVtoXYZ(face, uv[0], uv[1]);
}

inline void S2::ValidFaceXYZtoUV(int face, S2Point const& p,
                                 double* pu, double* pv) {
  DCHECK_GT(p.DotProd(GetNorm(face)), 0);
  switch (face) {
    case 0:  *pu =  p[1] / p[0]; *pv =  p[2] / p[0]; break;
    case 1:  *pu = -p[0] / p[1]; *pv =  p[2] / p[1]; break;
    case 2:  *pu = -p[0] / p[2]; *pv = -p[1] / p[2]; break;
    case 3:  *pu =  p[2] / p[0]; *pv =  p[1] / p[0]; break;
    case 4:  *pu =  p[2] / p[1]; *pv = -p[0] / p[1]; break;
    default: *pu = -p[1] / p[2]; *pv = -p[0] / p[2]; break;
  }
}

inline void S2::ValidFaceXYZtoUV(int face, S2Point const& p, R2Point* puv) {
  ValidFaceXYZtoUV(face, p, &(*puv)[0], &(*puv)[1]);
}

inline int S2::GetFace(S2Point const& p) {
  int face = p.LargestAbsComponent();
  if (p[face] < 0) face += 3;
  return face;
}

inline int S2::XYZtoFaceUV(S2Point const& p, double* pu, double* pv) {
  int face = GetFace(p);
  ValidFaceXYZtoUV(face, p, pu, pv);
  return face;
}

inline int S2::XYZtoFaceUV(S2Point const& p, R2Point* puv) {
  return XYZtoFaceUV(p, &(*puv)[0], &(*puv)[1]);
}

inline bool S2::FaceXYZtoUV(int face, S2Point const& p,
                            double* pu, double* pv) {
  if (face < 3) {
    if (p[face] <= 0) return false;
  } else {
    if (p[face-3] >= 0) return false;
  }
  ValidFaceXYZtoUV(face, p, pu, pv);
  return true;
}

inline bool S2::FaceXYZtoUV(int face, S2Point const& p, R2Point* puv) {
  return FaceXYZtoUV(face, p, &(*puv)[0], &(*puv)[1]);
}

inline S2Point S2::GetUNorm(int face, double u) {
  switch (face) {
    case 0:  return S2Point( u, -1,  0);
    case 1:  return S2Point( 1,  u,  0);
    case 2:  return S2Point( 1,  0,  u);
    case 3:  return S2Point(-u,  0,  1);
    case 4:  return S2Point( 0, -u,  1);
    default: return S2Point( 0, -1, -u);
  }
}

inline S2Point S2::GetVNorm(int face, double v) {
  switch (face) {
    case 0:  return S2Point(-v,  0,  1);
    case 1:  return S2Point( 0, -v,  1);
    case 2:  return S2Point( 0, -1, -v);
    case 3:  return S2Point( v, -1,  0);
    case 4:  return S2Point( 1,  v,  0);
    default: return S2Point( 1,  0,  v);
  }
}

inline S2Point S2::GetNorm(int face) {
  return GetUVWAxis(face, 2);
}

inline S2Point S2::GetUAxis(int face) {
  return GetUVWAxis(face, 0);
}

inline S2Point S2::GetVAxis(int face) {
  return GetUVWAxis(face, 1);
}

inline S2Point S2::GetUVWAxis(int face, int axis) {
  double const* p = kFaceUVWAxes[face][axis];
  return S2Point(p[0], p[1], p[2]);
}

inline int S2::GetUVWFace(int face, int axis, int direction) {
  DCHECK(face >= 0 && face <= 5);
  DCHECK(axis >= 0 && axis <= 2);
  DCHECK(direction >= 0 && direction <= 1);
  return kFaceUVWFaces[face][axis][direction];
}

template <int dim>
int S2::Metric<dim>::GetMinLevel(double value) const {
  if (value <= 0) return S2::kMaxCellLevel;

  // This code is equivalent to computing a floating-point "level"
  // value and rounding up.  ilogb() returns the exponent corresponding to a
  // fraction in the range [1,2).
  int level = ilogb(value / deriv_);
  level = std::max(0, std::min(S2::kMaxCellLevel, -(level >> (dim - 1))));
  DCHECK(level == S2::kMaxCellLevel || GetValue(level) <= value);
  DCHECK(level == 0 || GetValue(level - 1) > value);
  return level;
}

template <int dim>
int S2::Metric<dim>::GetMaxLevel(double value) const {
  if (value <= 0) return S2::kMaxCellLevel;

  // This code is equivalent to computing a floating-point "level"
  // value and rounding down.
  int level = ilogb(deriv_ / value);
  level = std::max(0, std::min(S2::kMaxCellLevel, level >> (dim - 1)));
  DCHECK(level == 0 || GetValue(level) >= value);
  DCHECK(level == S2::kMaxCellLevel || GetValue(level + 1) < value);
  return level;
}

template <int dim>
int S2::Metric<dim>::GetClosestLevel(double value) const {
  return GetMinLevel((dim == 1 ? M_SQRT2 : 2) * value);
}

#endif  // S2_S2_H_
