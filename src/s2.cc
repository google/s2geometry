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

#include "s2.h"

#include <float.h>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "util/bits/bits.h"
#include "s1angle.h"
#include "util/math/exactfloat/exactfloat.h"
#include "util/math/matrix3x3.h"

using std::max;
using std::min;

// Define storage for header file constants (the values are not needed
// here for integral constants).
int const S2::kMaxCellLevel;
int const S2::kSwapMask;
int const S2::kInvertMask;
int const S2::kLimitIJ;
unsigned int const S2::kMaxSiTi;

// kMaxDetError is the maximum error in computing (AxB).C where all vectors
// are unit length.  Using standard inequalities, it can be shown that
//
//  fl(AxB) = AxB + D where |D| <= (|AxB| + (2/sqrt(3))*|A|*|B|) * e
//
// where "fl()" denotes a calculation done in floating-point arithmetic,
// |x| denotes either absolute value or the L2-norm as appropriate, and
// e = 0.5*DBL_EPSILON.  Similarly,
//
//  fl(B.C) = B.C + d where |d| <= (1.5*|B.C| + 1.5*|B|*|C|) * e .
//
// Applying these bounds to the unit-length vectors A,B,C and neglecting
// relative error (which does not affect the sign of the result), we get
//
//  fl((AxB).C) = (AxB).C + d where |d| <= (2.5 + 2/sqrt(3)) * e
//
// which is about 3.6548 * e, or 1.8274 * DBL_EPSILON.
double const S2::kMaxDetError = 1.8274 * DBL_EPSILON;

int const S2::kFaceUVWFaces[6][3][2] = {
  { { 4, 1 }, { 5, 2 }, { 3, 0 } },
  { { 0, 3 }, { 5, 2 }, { 4, 1 } },
  { { 0, 3 }, { 1, 4 }, { 5, 2 } },
  { { 2, 5 }, { 1, 4 }, { 0, 3 } },
  { { 2, 5 }, { 3, 0 }, { 1, 4 } },
  { { 4, 1 }, { 3, 0 }, { 2, 5 } }
};

double const S2::kFaceUVWAxes[6][3][3] = {
  {
    { 0,  1,  0 },
    { 0,  0,  1 },
    { 1,  0,  0 }
  },
  {
    {-1,  0,  0 },
    { 0,  0,  1 },
    { 0,  1,  0 }
  },
  {
    {-1,  0,  0 },
    { 0, -1,  0 },
    { 0,  0,  1 }
  },
  {
    { 0,  0, -1 },
    { 0, -1,  0 },
    {-1,  0,  0 }
  },
  {
    { 0,  0, -1 },
    { 1,  0,  0 },
    { 0, -1,  0 }
  },
  {
    { 0,  1,  0 },
    { 1,  0,  0 },
    { 0,  0, -1 }
  }
};

COMPILE_ASSERT(S2::kSwapMask == 0x01 && S2::kInvertMask == 0x02,
               masks_changed);

DEFINE_bool(s2debug, !!google::DEBUG_MODE,
            "Enable debugging checks in s2 code");

S2Point S2::FaceXYZtoUVW(int face, S2Point const& p) {
  // The result coordinates are simply the dot products of P with the (u,v,w)
  // axes for the given face (see kFaceUVWAxes).
  switch (face) {
    case 0:  return S2Point( p.y(),  p.z(),  p.x());
    case 1:  return S2Point(-p.x(),  p.z(),  p.y());
    case 2:  return S2Point(-p.x(), -p.y(),  p.z());
    case 3:  return S2Point(-p.z(), -p.y(), -p.x());
    case 4:  return S2Point(-p.z(),  p.x(), -p.y());
    default: return S2Point( p.y(),  p.x(), -p.z());
  }
}

bool S2::IsUnitLength(S2Point const& p) {
  // Normalize() is guaranteed to return a vector whose L2-norm differs from 1
  // by less than 2 * DBL_EPSILON.  Thus the squared L2-norm differs by less
  // than 4 * DBL_EPSILON.  The actual calculated Norm2() can have up to 1.5 *
  // DBL_EPSILON of additional error.  The total error of 5.5 * DBL_EPSILON
  // can then be rounded down since the result must be a representable
  // double-precision value.
  return fabs(p.Norm2() - 1) <= 5 * DBL_EPSILON;  // About 1.11e-15
}

S2Point S2::Ortho(S2Point const& a) {
#ifdef S2_TEST_DEGENERACIES
  // Vector3::Ortho() always returns a point on the X-Y, Y-Z, or X-Z planes.
  // This leads to many more degenerate cases in polygon operations.
  return a.Ortho();
#else
  int k = a.LargestAbsComponent() - 1;
  if (k < 0) k = 2;
  S2Point temp(0.012, 0.0053, 0.00457);
  temp[k] = 1;
  return a.CrossProd(temp).Normalize();
#endif
}

void S2::GetFrame(S2Point const& z, Matrix3x3_d* m) {
  DCHECK(IsUnitLength(z));
  m->SetCol(2, z);
  m->SetCol(1, Ortho(z));
  m->SetCol(0, m->Col(1).CrossProd(z));  // Already unit-length.
}

S2Point S2::ToFrame(Matrix3x3_d const& m, S2Point const& p) {
  // The inverse of an orthonormal matrix is its transpose.
  return m.Transpose() * p;
}

S2Point S2::FromFrame(Matrix3x3_d const& m, S2Point const& q) {
  return m * q;
}

S2Point S2::Rotate(S2Point const& p, S2Point const& axis, S1Angle angle) {
  DCHECK(IsUnitLength(p));
  DCHECK(IsUnitLength(axis));
  // Let M be the plane through P that is perpendicular to "axis", and let
  // "center" be the point where M intersects "axis".  We construct a
  // right-handed orthogonal frame (dx, dy, center) such that "dx" is the
  // vector from "center" to P, and "dy" has the same length as "dx".  The
  // result can then be expressed as (cos(angle)*dx + sin(angle)*dy + center).
  S2Point center = p.DotProd(axis) * axis;
  S2Point dx = p - center;
  S2Point dy = axis.CrossProd(p);
  // Mathematically the result is unit length, but normalization is necessary
  // to ensure that numerical errors don't accumulate.
  return (cos(angle) * dx + sin(angle) * dy + center).Normalize();
}

bool S2::ApproxEquals(S2Point const& a, S2Point const& b, double max_error) {
  return a.Angle(b) <= max_error;
}

bool S2::PointsApproxEqual(S2Point const* a, int num_a,
                           S2Point const* b, int num_b,
                           double max_error) {
  if (num_a != num_b) return false;
  for (int i = 0; i < num_a; ++i) {
    if (!S2::ApproxEquals(a[i], b[i], max_error)) {
      return false;
    }
  }
  return true;
}

Vector3_d S2::RobustCrossProd(S2Point const& a, S2Point const& b) {
  // The direction of a.CrossProd(b) becomes unstable as (a + b) or (a - b)
  // approaches zero.  This leads to situations where a.CrossProd(b) is not
  // very orthogonal to "a" and/or "b".  We could fix this using Gram-Schmidt,
  // but we also want b.RobustCrossProd(a) == -a.RobustCrossProd(b).
  //
  // The easiest fix is to just compute the cross product of (b+a) and (b-a).
  // Mathematically, this cross product is exactly twice the cross product of
  // "a" and "b", but it has the numerical advantage that (b+a) and (b-a)
  // are always perpendicular (since "a" and "b" are unit length).  This
  // yields a result that is nearly orthogonal to both "a" and "b" even if
  // these two values differ only in the lowest bit of one component.

  DCHECK(IsUnitLength(a));
  DCHECK(IsUnitLength(b));
  Vector3_d x = (b + a).CrossProd(b - a);
  if (x != S2Point(0, 0, 0)) return x;

  // The only result that makes sense mathematically is to return zero, but
  // we find it more convenient to return an arbitrary orthogonal vector.
  return Ortho(a);
}

bool S2::SimpleCCW(S2Point const& a, S2Point const& b, S2Point const& c) {
  // We compute the signed volume of the parallelepiped ABC.  The usual
  // formula for this is (AxB).C, but we compute it here using (CxA).B
  // in order to ensure that ABC and CBA are not both CCW.  This follows
  // from the following identities (which are true numerically, not just
  // mathematically):
  //
  //     (1) x.CrossProd(y) == -(y.CrossProd(x))
  //     (2) (-x).DotProd(y) == -(x.DotProd(y))

  return c.CrossProd(a).DotProd(b) > 0;
}

int S2::RobustCCW(S2Point const& a, S2Point const& b, S2Point const& c) {
  // We don't need RobustCrossProd() here because RobustCCW() does its own
  // error estimation and calls ExpensiveCCW() if there is any uncertainty
  // about the result.
  return RobustCCW(a, b, c, a.CrossProd(b));
}

// ExpensiveCCW() uses arbitrary-precision arithmetic and the "simulation of
// simplicity" technique in order to be completely robust (i.e., to return
// consistent results for all possible inputs).
//
// Below we define a floating-point type with enough precision so that it can
// represent the exact determinant of any 3x3 matrix of floating-point
// numbers.  It uses ExactFloat, which is based on the OpenSSL Bignum library
// and therefore has a permissive BSD-style license.  (At one time we also
// supported an option based on MPFR, but that has an LGPL license and is
// therefore not suited for some applications.)

typedef Vector3<ExactFloat> Vector3_xf;

int S2::ExpensiveCCW(S2Point const& a, S2Point const& b, S2Point const& c) {
  // Return zero if and only if two points are the same.  This ensures (1).
  if (a == b || b == c || c == a) return 0;

  // Next we try recomputing the determinant still using floating-point
  // arithmetic but in a more precise way.  This is more expensive than the
  // simple calculation done by TriageCCW(), but it is still *much* cheaper
  // than using arbitrary-precision arithmetic.  This optimization is able to
  // compute the correct determinant sign in virtually all cases except when
  // the three points are truly collinear (e.g., three points on the equator).
  int det_sign = StableCCW(a, b, c);
  if (det_sign != 0) return det_sign;

  // Otherwise fall back to exact arithmetic and symbolic permutations.
  return ExactCCW(a, b, c);
}

// Compute the determinant in a numerically stable way.  Unlike TriageCCW(),
// this method can usually compute the correct determinant sign even when all
// three points are as collinear as possible.  For example if three points are
// spaced 1km apart along a random line on the Earth's surface using the
// nearest representable points, there is only a 0.4% chance that this method
// will not be able to find the determinant sign.  The probability of failure
// decreases as the points get closer together; if the collinear points are
// 1 meter apart, the failure rate drops to 0.0004%.
//
// This method could be extended to also handle nearly-antipodal points (and
// in fact an earlier version of this code did exactly that), but antipodal
// points are rare in practice so it seems better to simply fall back to
// exact arithmetic in that case.
int S2::StableCCW(S2Point const& a, S2Point const& b, S2Point const& c) {
  Vector3_d ab = b - a;
  Vector3_d bc = c - b;
  Vector3_d ca = a - c;
  double ab2 = ab.Norm2();
  double bc2 = bc.Norm2();
  double ca2 = ca.Norm2();

  // Now compute the determinant ((A-C)x(B-C)).C, where the vertices have been
  // cyclically permuted if necessary so that AB is the longest edge.  (This
  // minimizes the magnitude of cross product.)  At the same time we also
  // compute the maximum error in the determinant.  Using a similar technique
  // to the one used for kMaxDetError, the error is at most
  //
  //   |d| <= (3 + 6/sqrt(3)) * |A-C| * |B-C| * e
  //
  // where e = 0.5 * DBL_EPSILON.  If the determinant magnitude is larger than
  // this value then we know its sign with certainty.
  double const kDetErrorMultiplier = 3.2321 * DBL_EPSILON;  // see above
  double det, max_error;
  if (ab2 >= bc2 && ab2 >= ca2) {
    // AB is the longest edge, so compute (A-C)x(B-C).C.
    det = -(ca.CrossProd(bc).DotProd(c));
    max_error = kDetErrorMultiplier * sqrt(ca2 * bc2);
  } else if (bc2 >= ca2) {
    // BC is the longest edge, so compute (B-A)x(C-A).A.
    det = -(ab.CrossProd(ca).DotProd(a));
    max_error = kDetErrorMultiplier * sqrt(ab2 * ca2);
  } else {
    // CA is the longest edge, so compute (C-B)x(A-B).B.
    det = -(bc.CrossProd(ab).DotProd(b));
    max_error = kDetErrorMultiplier * sqrt(bc2 * ab2);
  }
  return (fabs(det) <= max_error) ? 0 : (det > 0) ? 1 : -1;
}

// Forward declaration.
static int SymbolicallyPerturbedCCW(
    Vector3_xf const& a, Vector3_xf const& b,
    Vector3_xf const& c, Vector3_xf const& b_cross_c);

// Compute the determinant using exact arithmetic and/or symbolic
// permutations.  Requires that the three points are distinct.
int S2::ExactCCW(S2Point const& a, S2Point const& b, S2Point const& c) {
  DCHECK(a != b && b != c && c != a);

  // Sort the three points in lexicographic order, keeping track of the sign
  // of the permutation.  (Each exchange inverts the sign of the determinant.)
  int perm_sign = 1;
  const S2Point *pa = &a, *pb = &b, *pc = &c;
  using std::swap;
  if (*pa > *pb) { swap(pa, pb); perm_sign = -perm_sign; }
  if (*pb > *pc) { swap(pb, pc); perm_sign = -perm_sign; }
  if (*pa > *pb) { swap(pa, pb); perm_sign = -perm_sign; }
  DCHECK(*pa < *pb && *pb < *pc);

  // Construct multiple-precision versions of the sorted points and compute
  // their exact 3x3 determinant.
  Vector3_xf xa = Vector3_xf::Cast(*pa);
  Vector3_xf xb = Vector3_xf::Cast(*pb);
  Vector3_xf xc = Vector3_xf::Cast(*pc);
  Vector3_xf xb_cross_xc = xb.CrossProd(xc);
  ExactFloat det = xa.DotProd(xb_cross_xc);

  // The precision of ExactFloat is high enough that the result should always
  // be exact (no rounding was performed).
  DCHECK(!det.is_nan());
  DCHECK_LT(det.prec(), det.max_prec());

  // If the exact determinant is non-zero, we're done.
  int det_sign = det.sgn();
  if (det_sign == 0) {
    // Otherwise, we need to resort to symbolic perturbations to resolve the
    // sign of the determinant.
    det_sign = SymbolicallyPerturbedCCW(xa, xb, xc, xb_cross_xc);
  }
  DCHECK_NE(0, det_sign);
  return perm_sign * det_sign;
}

// The following function returns the sign of the determinant of three points
// A, B, C under a model where every possible S2Point is slightly perturbed by
// a unique infinitesmal amount such that no three perturbed points are
// collinear and no four points are coplanar.  The perturbations are so small
// that they do not change the sign of any determinant that was non-zero
// before the perturbations, and therefore can be safely ignored unless the
// determinant of three points is exactly zero (using multiple-precision
// arithmetic).
//
// Since the symbolic perturbation of a given point is fixed (i.e., the
// perturbation is the same for all calls to this method and does not depend
// on the other two arguments), the results of this method are always
// self-consistent.  It will never return results that would correspond to an
// "impossible" configuration of non-degenerate points.
//
// Requirements:
//   The 3x3 determinant of A, B, C must be exactly zero.
//   The points must be distinct, with A < B < C in lexicographic order.
//
// Returns:
//   +1 or -1 according to the sign of the determinant after the symbolic
// perturbations are taken into account.
//
// Reference:
//   "Simulation of Simplicity" (Edelsbrunner and Muecke, ACM Transactions on
//   Graphics, 1990).
//
static int SymbolicallyPerturbedCCW(
    Vector3_xf const& a, Vector3_xf const& b,
    Vector3_xf const& c, Vector3_xf const& b_cross_c) {
  // This method requires that the points are sorted in lexicographically
  // increasing order.  This is because every possible S2Point has its own
  // symbolic perturbation such that if A < B then the symbolic perturbation
  // for A is much larger than the perturbation for B.
  //
  // Alternatively, we could sort the points in this method and keep track of
  // the sign of the permutation, but it is more efficient to do this before
  // converting the inputs to the multi-precision representation, and this
  // also lets us re-use the result of the cross product B x C.
  DCHECK(a < b && b < c);

  // Every input coordinate x[i] is assigned a symbolic perturbation dx[i].
  // We then compute the sign of the determinant of the perturbed points,
  // i.e.
  //               | a[0]+da[0]  a[1]+da[1]  a[2]+da[2] |
  //               | b[0]+db[0]  b[1]+db[1]  b[2]+db[2] |
  //               | c[0]+dc[0]  c[1]+dc[1]  c[2]+dc[2] |
  //
  // The perturbations are chosen such that
  //
  //   da[2] > da[1] > da[0] > db[2] > db[1] > db[0] > dc[2] > dc[1] > dc[0]
  //
  // where each perturbation is so much smaller than the previous one that we
  // don't even need to consider it unless the coefficients of all previous
  // perturbations are zero.  In fact, it is so small that we don't need to
  // consider it unless the coefficient of all products of the previous
  // perturbations are zero.  For example, we don't need to consider the
  // coefficient of db[1] unless the coefficient of db[2]*da[0] is zero.
  //
  // The follow code simply enumerates the coefficients of the perturbations
  // (and products of perturbations) that appear in the determinant above, in
  // order of decreasing perturbation magnitude.  The first non-zero
  // coefficient determines the sign of the result.  The easiest way to
  // enumerate the coefficients in the correct order is to pretend that each
  // perturbation is some tiny value "eps" raised to a power of two:
  //
  // eps**    1      2      4      8     16     32     64     128    256
  //        da[2]  da[1]  da[0]  db[2]  db[1]  db[0]  dc[2]  dc[1]  dc[0]
  //
  // Essentially we can then just count in binary and test the corresponding
  // subset of perturbations at each step.  So for example, we must test the
  // coefficient of db[2]*da[0] before db[1] because eps**12 > eps**16.
  //
  // Of course, not all products of these perturbations appear in the
  // determinant above, since the determinant only contains the products of
  // elements in distinct rows and columns.  Thus we don't need to consider
  // da[2]*da[1], db[1]*da[1], etc.  Furthermore, sometimes different pairs of
  // perturbations have the same coefficient in the determinant; for example,
  // da[1]*db[0] and db[1]*da[0] have the same coefficient (c[2]).  Therefore
  // we only need to test this coefficient the first time we encounter it in
  // the binary order above (which will be db[1]*da[0]).
  //
  // The sequence of tests below also appears in Table 4-ii of the paper
  // referenced above, if you just want to look it up, with the following
  // translations: [a,b,c] -> [i,j,k] and [0,1,2] -> [1,2,3].  Also note that
  // some of the signs are different because the opposite cross product is
  // used (e.g., B x C rather than C x B).

  int det_sign = b_cross_c[2].sgn();            // da[2]
  if (det_sign != 0) return det_sign;
  det_sign = b_cross_c[1].sgn();                // da[1]
  if (det_sign != 0) return det_sign;
  det_sign = b_cross_c[0].sgn();                // da[0]
  if (det_sign != 0) return det_sign;

  det_sign = (c[0]*a[1] - c[1]*a[0]).sgn();     // db[2]
  if (det_sign != 0) return det_sign;
  det_sign = c[0].sgn();                        // db[2] * da[1]
  if (det_sign != 0) return det_sign;
  det_sign = -(c[1].sgn());                     // db[2] * da[0]
  if (det_sign != 0) return det_sign;
  det_sign = (c[2]*a[0] - c[0]*a[2]).sgn();     // db[1]
  if (det_sign != 0) return det_sign;
  det_sign = c[2].sgn();                        // db[1] * da[0]
  if (det_sign != 0) return det_sign;
  // The following test is listed in the paper, but it is redundant because
  // the previous tests guarantee that C == (0, 0, 0).
  DCHECK_EQ(0, (c[1]*a[2] - c[2]*a[1]).sgn());  // db[0]

  det_sign = (a[0]*b[1] - a[1]*b[0]).sgn();     // dc[2]
  if (det_sign != 0) return det_sign;
  det_sign = -(b[0].sgn());                     // dc[2] * da[1]
  if (det_sign != 0) return det_sign;
  det_sign = b[1].sgn();                        // dc[2] * da[0]
  if (det_sign != 0) return det_sign;
  det_sign = a[0].sgn();                        // dc[2] * db[1]
  if (det_sign != 0) return det_sign;
  return 1;                                     // dc[2] * db[1] * da[0]
}

double S2::Angle(S2Point const& a, S2Point const& b, S2Point const& c) {
  // RobustCrossProd() is necessary to get good accuracy when two of the input
  // points are very close together.
  return RobustCrossProd(a, b).Angle(RobustCrossProd(c, b));
}

double S2::TurnAngle(S2Point const& a, S2Point const& b, S2Point const& c) {
  // We use RobustCrossProd() to get good accuracy when two points are very
  // close together, and RobustCCW() to ensure that the sign is correct for
  // turns that are close to 180 degrees.
  //
  // Unfortunately we can't save RobustCrossProd(a, b) and pass it as the
  // optional 4th argument to RobustCCW(), because RobustCCW() requires
  // a.CrossProd(b) exactly (the robust version differs in magnitude).
  double angle = RobustCrossProd(a, b).Angle(RobustCrossProd(b, c));

  // Don't return RobustCCW() * angle because it is legal to have (a == c).
  return (RobustCCW(a, b, c) > 0) ? angle : -angle;
}

double S2::Area(S2Point const& a, S2Point const& b, S2Point const& c) {
  DCHECK(IsUnitLength(a));
  DCHECK(IsUnitLength(b));
  DCHECK(IsUnitLength(c));
  // This method is based on l'Huilier's theorem,
  //
  //   tan(E/4) = sqrt(tan(s/2) tan((s-a)/2) tan((s-b)/2) tan((s-c)/2))
  //
  // where E is the spherical excess of the triangle (i.e. its area),
  //       a, b, c, are the side lengths, and
  //       s is the semiperimeter (a + b + c) / 2 .
  //
  // The only significant source of error using l'Huilier's method is the
  // cancellation error of the terms (s-a), (s-b), (s-c).  This leads to a
  // *relative* error of about 1e-16 * s / min(s-a, s-b, s-c).  This compares
  // to a relative error of about 1e-15 / E using Girard's formula, where E is
  // the true area of the triangle.  Girard's formula can be even worse than
  // this for very small triangles, e.g. a triangle with a true area of 1e-30
  // might evaluate to 1e-5.
  //
  // So, we prefer l'Huilier's formula unless dmin < s * (0.1 * E), where
  // dmin = min(s-a, s-b, s-c).  This basically includes all triangles
  // except for extremely long and skinny ones.
  //
  // Since we don't know E, we would like a conservative upper bound on
  // the triangle area in terms of s and dmin.  It's possible to show that
  // E <= k1 * s * sqrt(s * dmin), where k1 = 2*sqrt(3)/Pi (about 1).
  // Using this, it's easy to show that we should always use l'Huilier's
  // method if dmin >= k2 * s^5, where k2 is about 1e-2.  Furthermore,
  // if dmin < k2 * s^5, the triangle area is at most k3 * s^4, where
  // k3 is about 0.1.  Since the best case error using Girard's formula
  // is about 1e-15, this means that we shouldn't even consider it unless
  // s >= 3e-4 or so.

  // We use volatile doubles to force the compiler to truncate all of these
  // quantities to 64 bits.  Otherwise it may compute a value of dmin > 0
  // simply because it chose to spill one of the intermediate values to
  // memory but not one of the others.
  volatile double sa = b.Angle(c);
  volatile double sb = c.Angle(a);
  volatile double sc = a.Angle(b);
  volatile double s = 0.5 * (sa + sb + sc);
  if (s >= 3e-4) {
    // Consider whether Girard's formula might be more accurate.
    double s2 = s * s;
    double dmin = s - max(sa, max(sb, sc));
    if (dmin < 1e-2 * s * s2 * s2) {
      // This triangle is skinny enough to consider Girard's formula.
      double area = GirardArea(a, b, c);
      if (dmin < s * (0.1 * area)) return area;
    }
  }
  // Use l'Huilier's formula.
  return 4 * atan(sqrt(max(0.0, tan(0.5 * s) * tan(0.5 * (s - sa)) *
                           tan(0.5 * (s - sb)) * tan(0.5 * (s - sc)))));
}

double S2::GirardArea(S2Point const& a, S2Point const& b, S2Point const& c) {
  // This is equivalent to the usual Girard's formula but is slightly more
  // accurate, faster to compute, and handles a == b == c without a special
  // case.  RobustCrossProd() is necessary to get good accuracy when two of
  // the input points are very close together.

  Vector3_d ab = RobustCrossProd(a, b);
  Vector3_d bc = RobustCrossProd(b, c);
  Vector3_d ac = RobustCrossProd(a, c);
  return max(0.0, ab.Angle(ac) - ab.Angle(bc) + bc.Angle(ac));
}

double S2::SignedArea(S2Point const& a, S2Point const& b, S2Point const& c) {
  return RobustCCW(a, b, c) * Area(a, b, c);
}

S2Point S2::PlanarCentroid(S2Point const& a, S2Point const& b,
                           S2Point const& c) {
  return (1./3) * (a + b + c);
}

S2Point S2::TrueCentroid(S2Point const& a, S2Point const& b,
                         S2Point const& c) {
  DCHECK(IsUnitLength(a));
  DCHECK(IsUnitLength(b));
  DCHECK(IsUnitLength(c));

  // I couldn't find any references for computing the true centroid of a
  // spherical triangle...  I have a truly marvellous demonstration of this
  // formula which this margin is too narrow to contain :)

  // Use Angle() in order to get accurate results for small triangles.
  double angle_a = b.Angle(c);
  double angle_b = c.Angle(a);
  double angle_c = a.Angle(b);
  double ra = (angle_a == 0) ? 1 : (angle_a / sin(angle_a));
  double rb = (angle_b == 0) ? 1 : (angle_b / sin(angle_b));
  double rc = (angle_c == 0) ? 1 : (angle_c / sin(angle_c));

  // Now compute a point M such that:
  //
  //  [Ax Ay Az] [Mx]                       [ra]
  //  [Bx By Bz] [My]  = 0.5 * det(A,B,C) * [rb]
  //  [Cx Cy Cz] [Mz]                       [rc]
  //
  // To improve the numerical stability we subtract the first row (A) from the
  // other two rows; this reduces the cancellation error when A, B, and C are
  // very close together.  Then we solve it using Cramer's rule.
  //
  // TODO(ericv): This code still isn't as numerically stable as it could be.
  // The biggest potential improvement is to compute B-A and C-A more
  // accurately so that (B-A)x(C-A) is always inside triangle ABC.
  S2Point x(a.x(), b.x() - a.x(), c.x() - a.x());
  S2Point y(a.y(), b.y() - a.y(), c.y() - a.y());
  S2Point z(a.z(), b.z() - a.z(), c.z() - a.z());
  S2Point r(ra, rb - ra, rc - ra);
  return 0.5 * S2Point(y.CrossProd(z).DotProd(r),
                       z.CrossProd(x).DotProd(r),
                       x.CrossProd(y).DotProd(r));
}

bool S2::OrderedCCW(S2Point const& a, S2Point const& b, S2Point const& c,
                    S2Point const& o) {
  // The last inequality below is ">" rather than ">=" so that we return true
  // if A == B or B == C, and otherwise false if A == C.  Recall that
  // RobustCCW(x,y,z) == -RobustCCW(z,y,x) for all x,y,z.

  int sum = 0;
  if (RobustCCW(b, o, a) >= 0) ++sum;
  if (RobustCCW(c, o, b) >= 0) ++sum;
  if (RobustCCW(a, o, c) > 0) ++sum;
  return sum >= 2;
}

int S2::XYZtoFaceSiTi(S2Point const& p, int* face, unsigned int* si,
                      unsigned int* ti) {
  double u, v;
  *face = XYZtoFaceUV(p, &u, &v);
  *si = STtoSiTi(UVtoST(u));
  *ti = STtoSiTi(UVtoST(v));
  // If the levels corresponding to si,ti are not equal, then p is not a cell
  // center.  The si,ti values 0 and kMaxSiTi need to be handled specially
  // because they do not correspond to cell centers at any valid level; they
  // are mapped to level -1 by the code below.
  int level = kMaxCellLevel - Bits::FindLSBSetNonZero(*si | kMaxSiTi);
  if (level < 0 ||
      level != kMaxCellLevel - Bits::FindLSBSetNonZero(*ti | kMaxSiTi)) {
    return -1;
  }
  DCHECK_LE(level, kMaxCellLevel);
  // In infinite precision, this test could be changed to ST == SiTi. However,
  // due to rounding errors, UVtoST(XYZtoFaceUV(FaceUVtoXYZ(STtoUV(...)))) is
  // not idempotent. On the other hand, center_raw is computed exactly the same
  // way p was originally computed (if it is indeed the center of an S2Cell):
  // the comparison can be exact.
  S2Point center = FaceSiTitoXYZ(*face, *si, *ti).Normalize();
  return p == center ? level : -1;
}

S2Point S2::FaceSiTitoXYZ(int face, unsigned int si, unsigned int ti) {
  double u = STtoUV(S2::SiTitoST(si));
  double v = STtoUV(S2::SiTitoST(ti));
  return FaceUVtoXYZ(face, u, v);
}

// kIJtoPos[orientation][ij] -> pos
int const S2::kIJtoPos[4][4] = {
  // (0,0) (0,1) (1,0) (1,1)
  {     0,    1,    3,    2  },  // canonical order
  {     0,    3,    1,    2  },  // axes swapped
  {     2,    3,    1,    0  },  // bits inverted
  {     2,    1,    3,    0  },  // swapped & inverted
};

// kPosToIJ[orientation][pos] -> ij
int const S2::kPosToIJ[4][4] = {
  // 0  1  2  3
  {  0, 1, 3, 2 },    // canonical order:    (0,0), (0,1), (1,1), (1,0)
  {  0, 2, 3, 1 },    // axes swapped:       (0,0), (1,0), (1,1), (0,1)
  {  3, 2, 0, 1 },    // bits inverted:      (1,1), (1,0), (0,0), (0,1)
  {  3, 1, 0, 2 },    // swapped & inverted: (1,1), (0,1), (0,0), (1,0)
};

// kPosToOrientation[pos] -> orientation_modifier
int const S2::kPosToOrientation[4] = {
  kSwapMask,
  0,
  0,
  kInvertMask + kSwapMask,
};

// All of the values below were obtained by a combination of hand analysis and
// Mathematica.  In general, S2_TAN_PROJECTION produces the most uniform
// shapes and sizes of cells, S2_LINEAR_PROJECTION is considerably worse, and
// S2_QUADRATIC_PROJECTION is somewhere in between (but generally closer to
// the tangent projection than the linear one).

S2::LengthMetric const S2::kMinAngleSpan(
    S2_PROJECTION == S2_LINEAR_PROJECTION ? 1.0 :                      // 1.000
    S2_PROJECTION == S2_TAN_PROJECTION ? M_PI / 2 :                    // 1.571
    S2_PROJECTION == S2_QUADRATIC_PROJECTION ? 4. / 3 :                // 1.333
    0);

S2::LengthMetric const S2::kMaxAngleSpan(
    S2_PROJECTION == S2_LINEAR_PROJECTION ? 2 :                        // 2.000
    S2_PROJECTION == S2_TAN_PROJECTION ? M_PI / 2 :                    // 1.571
    S2_PROJECTION == S2_QUADRATIC_PROJECTION ? 1.704897179199218452 :  // 1.705
    0);

S2::LengthMetric const S2::kAvgAngleSpan(M_PI / 2);                    // 1.571
// This is true for all projections.

S2::LengthMetric const S2::kMinWidth(
    S2_PROJECTION == S2_LINEAR_PROJECTION ? sqrt(2. / 3) :             // 0.816
    S2_PROJECTION == S2_TAN_PROJECTION ? M_PI / (2 * sqrt(2)) :        // 1.111
    S2_PROJECTION == S2_QUADRATIC_PROJECTION ? 2 * sqrt(2) / 3 :       // 0.943
    0);

S2::LengthMetric const S2::kMaxWidth(S2::kMaxAngleSpan.deriv());
// This is true for all projections.

S2::LengthMetric const S2::kAvgWidth(
    S2_PROJECTION == S2_LINEAR_PROJECTION ? 1.411459345844456965 :     // 1.411
    S2_PROJECTION == S2_TAN_PROJECTION ? 1.437318638925160885 :        // 1.437
    S2_PROJECTION == S2_QUADRATIC_PROJECTION ? 1.434523672886099389 :  // 1.435
    0);

S2::LengthMetric const S2::kMinEdge(
    S2_PROJECTION == S2_LINEAR_PROJECTION ? 2 * sqrt(2) / 3 :          // 0.943
    S2_PROJECTION == S2_TAN_PROJECTION ? M_PI / (2 * sqrt(2)) :        // 1.111
    S2_PROJECTION == S2_QUADRATIC_PROJECTION ? 2 * sqrt(2) / 3 :       // 0.943
    0);

S2::LengthMetric const S2::kMaxEdge(S2::kMaxAngleSpan.deriv());
// This is true for all projections.

S2::LengthMetric const S2::kAvgEdge(
    S2_PROJECTION == S2_LINEAR_PROJECTION ? 1.440034192955603643 :     // 1.440
    S2_PROJECTION == S2_TAN_PROJECTION ? 1.461667032546739266 :        // 1.462
    S2_PROJECTION == S2_QUADRATIC_PROJECTION ? 1.459213746386106062 :  // 1.459
    0);

S2::LengthMetric const S2::kMinDiag(
    S2_PROJECTION == S2_LINEAR_PROJECTION ? 2 * sqrt(2) / 3 :          // 0.943
    S2_PROJECTION == S2_TAN_PROJECTION ? M_PI * sqrt(2) / 3 :          // 1.481
    S2_PROJECTION == S2_QUADRATIC_PROJECTION ? 8 * sqrt(2) / 9 :       // 1.257
    0);

S2::LengthMetric const S2::kMaxDiag(
    S2_PROJECTION == S2_LINEAR_PROJECTION ? 2 * sqrt(2) :              // 2.828
    S2_PROJECTION == S2_TAN_PROJECTION ? M_PI * sqrt(2. / 3) :         // 2.565
    S2_PROJECTION == S2_QUADRATIC_PROJECTION ? 2.438654594434021032 :  // 2.439
    0);

S2::LengthMetric const S2::kAvgDiag(
    S2_PROJECTION == S2_LINEAR_PROJECTION ? 2.031817866418812674 :     // 2.032
    S2_PROJECTION == S2_TAN_PROJECTION ? 2.063623197195635753 :        // 2.064
    S2_PROJECTION == S2_QUADRATIC_PROJECTION ? 2.060422738998471683 :  // 2.060
    0);

S2::AreaMetric const S2::kMinArea(
    S2_PROJECTION == S2_LINEAR_PROJECTION ? 4 / (3 * sqrt(3)) :        // 0.770
    S2_PROJECTION == S2_TAN_PROJECTION ? (M_PI*M_PI) / (4*sqrt(2)) :   // 1.745
    S2_PROJECTION == S2_QUADRATIC_PROJECTION ? 8 * sqrt(2) / 9 :       // 1.257
    0);

S2::AreaMetric const S2::kMaxArea(
    S2_PROJECTION == S2_LINEAR_PROJECTION ? 4 :                        // 4.000
    S2_PROJECTION == S2_TAN_PROJECTION ? M_PI * M_PI / 4 :             // 2.467
    S2_PROJECTION == S2_QUADRATIC_PROJECTION ? 2.635799256963161491 :  // 2.636
    0);

S2::AreaMetric const S2::kAvgArea(4 * M_PI / 6);                       // 2.094
// This is true for all projections.

double const S2::kMaxEdgeAspect = (
    S2_PROJECTION == S2_LINEAR_PROJECTION ? sqrt(2) :                  // 1.414
    S2_PROJECTION == S2_TAN_PROJECTION ?  sqrt(2) :                    // 1.414
    S2_PROJECTION == S2_QUADRATIC_PROJECTION ? 1.442615274452682920 :  // 1.443
    0);

double const S2::kMaxDiagAspect = sqrt(3);                             // 1.732
// This is true for all projections.
