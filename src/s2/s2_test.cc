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

#include "s2/s2.h"

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iterator>
#include <unordered_set>
#include <vector>

#include "s2/base/integral_types.h"
#include <glog/logging.h>
#include <gtest/gtest.h>
#include "s2/s2cell.h"
#include "s2/s2cellid.h"
#include "s2/s2edgeutil.h"
#include "s2/s2latlng.h"
#include "s2/s2testing.h"
#include "s2/util/math/matrix3x3.h"

using std::max;
using std::min;
using std::unordered_set;
using std::vector;

inline static int SwapAxes(int ij) {
  return ((ij >> 1) & 1) + ((ij & 1) << 1);
}

inline static int InvertBits(int ij) {
  return ij ^ 3;
}

TEST(S2, TraversalOrder) {
  for (int r = 0; r < 4; ++r) {
    for (int i = 0; i < 4; ++i) {
      // Check consistency with respect to swapping axes.
      EXPECT_EQ(S2::kIJtoPos[r][i],
                S2::kIJtoPos[r ^ S2::kSwapMask][SwapAxes(i)]);
      EXPECT_EQ(S2::kPosToIJ[r][i],
                SwapAxes(S2::kPosToIJ[r ^ S2::kSwapMask][i]));

      // Check consistency with respect to reversing axis directions.
      EXPECT_EQ(S2::kIJtoPos[r][i],
                S2::kIJtoPos[r ^ S2::kInvertMask][InvertBits(i)]);
      EXPECT_EQ(S2::kPosToIJ[r][i],
                InvertBits(S2::kPosToIJ[r ^ S2::kInvertMask][i]));

      // Check that the two tables are inverses of each other.
      EXPECT_EQ(S2::kIJtoPos[r][S2::kPosToIJ[r][i]], i);
      EXPECT_EQ(S2::kPosToIJ[r][S2::kIJtoPos[r][i]], i);
    }
  }
}

TEST(S2, ST_UV_Conversions) {
  // Check boundary conditions.
  for (double s = 0; s <= 1; s += 0.5) {
    volatile double u = S2::STtoUV(s);
    EXPECT_EQ(u, 2 * s - 1);
  }
  for (double u = -1; u <= 1; ++u) {
    volatile double s = S2::UVtoST(u);
    EXPECT_EQ(s, 0.5 * (u + 1));
  }
  // Check that UVtoST and STtoUV are inverses.
  for (double x = 0; x <= 1; x += 0.0001) {
    EXPECT_NEAR(S2::UVtoST(S2::STtoUV(x)), x, 1e-15);
    EXPECT_NEAR(S2::STtoUV(S2::UVtoST(2*x-1)), 2*x-1, 1e-15);
  }
}

TEST(S2, FaceUVtoXYZ) {
  // Check that each face appears exactly once.
  S2Point sum;
  for (int face = 0; face < 6; ++face) {
    S2Point center = S2::FaceUVtoXYZ(face, 0, 0);
    EXPECT_EQ(S2::GetNorm(face), center);
    EXPECT_EQ(fabs(center[center.LargestAbsComponent()]), 1);
    sum += center.Fabs();
  }
  EXPECT_EQ(sum, S2Point(2, 2, 2));

  // Check that each face has a right-handed coordinate system.
  for (int face = 0; face < 6; ++face) {
    EXPECT_EQ(S2::GetUAxis(face).CrossProd(S2::GetVAxis(face))
              .DotProd(S2::FaceUVtoXYZ(face, 0, 0)), 1);
  }

  // Check that the Hilbert curves on each face combine to form a
  // continuous curve over the entire cube.
  for (int face = 0; face < 6; ++face) {
    // The Hilbert curve on each face starts at (-1,-1) and terminates
    // at either (1,-1) (if axes not swapped) or (-1,1) (if swapped).
    int sign = (face & S2::kSwapMask) ? -1 : 1;
    EXPECT_EQ(S2::FaceUVtoXYZ(face, sign, -sign),
              S2::FaceUVtoXYZ((face + 1) % 6, -1, -1));
  }
}

TEST(S2, FaceXYZtoUVW) {
  for (int face = 0; face < 6; ++face) {
    EXPECT_EQ(S2Point( 0,  0,  0), S2::FaceXYZtoUVW(face, S2Point(0, 0, 0)));
    EXPECT_EQ(S2Point( 1,  0,  0), S2::FaceXYZtoUVW(face,  S2::GetUAxis(face)));
    EXPECT_EQ(S2Point(-1,  0,  0), S2::FaceXYZtoUVW(face, -S2::GetUAxis(face)));
    EXPECT_EQ(S2Point( 0,  1,  0), S2::FaceXYZtoUVW(face,  S2::GetVAxis(face)));
    EXPECT_EQ(S2Point( 0, -1,  0), S2::FaceXYZtoUVW(face, -S2::GetVAxis(face)));
    EXPECT_EQ(S2Point( 0,  0,  1), S2::FaceXYZtoUVW(face,  S2::GetNorm(face)));
    EXPECT_EQ(S2Point( 0,  0, -1), S2::FaceXYZtoUVW(face, -S2::GetNorm(face)));
  }
}

TEST(S2, XYZToFaceSiTi) {
  // Check the conversion of random cells to center points and back.
  for (int level = 0; level <= S2CellId::kMaxLevel; ++level) {
    for (int i = 0; i < 1000; ++i) {
      S2CellId id = S2Testing::GetRandomCellId(level);
      int face;
      unsigned int si, ti;
      int actual_level = S2::XYZtoFaceSiTi(id.ToPoint(), &face, &si, &ti);
      EXPECT_EQ(level, actual_level);
      S2CellId actual_id =
          S2CellId::FromFaceIJ(face, si / 2, ti / 2).parent(level);
      EXPECT_EQ(id, actual_id);

      // Now test a point near the cell center but not equal to it.
      S2Point p_moved = id.ToPoint() + S2Point(1e-13, 1e-13, 1e-13);
      int face_moved;
      unsigned int si_moved, ti_moved;
      actual_level = S2::XYZtoFaceSiTi(p_moved, &face_moved, &si_moved,
                                       &ti_moved);
      EXPECT_EQ(-1, actual_level);
      EXPECT_EQ(face, face_moved);
      EXPECT_EQ(si, si_moved);
      EXPECT_EQ(ti, ti_moved);

      // Finally, test some random (si,ti) values that may be at different
      // levels, or not at a valid level at all (for example, si == 0).
      int face_random = S2Testing::rnd.Uniform(S2CellId::kNumFaces);
      unsigned int si_random, ti_random;
      unsigned int mask = -1 << (S2CellId::kMaxLevel - level);
      do {
        si_random = S2Testing::rnd.Rand32() & mask;
        ti_random = S2Testing::rnd.Rand32() & mask;
      } while (si_random > S2::kMaxSiTi || ti_random > S2::kMaxSiTi);
      S2Point p_random = S2::FaceSiTitoXYZ(face_random, si_random, ti_random);
      actual_level = S2::XYZtoFaceSiTi(p_random, &face, &si, &ti);
      if (face != face_random) {
        // The chosen point is on the edge of a top-level face cell.
        EXPECT_EQ(-1, actual_level);
        EXPECT_TRUE(si == 0 || si == S2::kMaxSiTi ||
                    ti == 0 || ti == S2::kMaxSiTi);
      } else {
        EXPECT_EQ(si_random, si);
        EXPECT_EQ(ti_random, ti);
        if (actual_level >= 0) {
          EXPECT_EQ(p_random,
                    S2CellId::FromFaceIJ(face, si / 2, ti / 2).
                    parent(actual_level).ToPoint());
        }
      }
    }
  }
}

TEST(S2, UVNorms) {
  // Check that GetUNorm and GetVNorm compute right-handed normals for
  // an edge in the increasing U or V direction.
  for (int face = 0; face < 6; ++face) {
    for (double x = -1; x <= 1; x += 1/1024.) {
      EXPECT_DOUBLE_EQ(S2::FaceUVtoXYZ(face, x, -1)
                       .CrossProd(S2::FaceUVtoXYZ(face, x, 1))
                       .Angle(S2::GetUNorm(face, x)), 0);
      EXPECT_DOUBLE_EQ(S2::FaceUVtoXYZ(face, -1, x)
                       .CrossProd(S2::FaceUVtoXYZ(face, 1, x))
                       .Angle(S2::GetVNorm(face, x)), 0);
    }
  }
}

TEST(S2, UVWAxis) {
  for (int face = 0; face < 6; ++face) {
    // Check that axes are consistent with FaceUVtoXYZ.
    EXPECT_EQ(S2::FaceUVtoXYZ(face, 1, 0) - S2::FaceUVtoXYZ(face, 0, 0),
              S2::GetUAxis(face));
    EXPECT_EQ(S2::FaceUVtoXYZ(face, 0, 1) - S2::FaceUVtoXYZ(face, 0, 0),
              S2::GetVAxis(face));
    EXPECT_EQ(S2::FaceUVtoXYZ(face, 0, 0), S2::GetNorm(face));

    // Check that every face coordinate frame is right-handed.
    EXPECT_EQ(1, S2::GetUAxis(face).CrossProd(S2::GetVAxis(face))
              .DotProd(S2::GetNorm(face)));

    // Check that GetUVWAxis is consistent with GetUAxis, GetVAxis, GetNorm.
    EXPECT_EQ(S2::GetUAxis(face), S2::GetUVWAxis(face, 0));
    EXPECT_EQ(S2::GetVAxis(face), S2::GetUVWAxis(face, 1));
    EXPECT_EQ(S2::GetNorm(face), S2::GetUVWAxis(face, 2));
  }
}

TEST(S2, UVWFace) {
  // Check that GetUVWFace is consistent with GetUVWAxis.
  for (int face = 0; face < 6; ++face) {
    for (int axis = 0; axis < 3; ++axis) {
      EXPECT_EQ(S2::GetFace(-S2::GetUVWAxis(face, axis)),
                S2::GetUVWFace(face, axis, 0));
      EXPECT_EQ(S2::GetFace(S2::GetUVWAxis(face, axis)),
                S2::GetUVWFace(face, axis, 1));
    }
  }
}

TEST(S2, AngleMethods) {
  S2Point pz(0, 0, 1);
  S2Point p000(1, 0, 0);
  S2Point p045 = S2Point(1, 1, 0).Normalize();
  S2Point p090(0, 1, 0);
  S2Point p180(-1, 0, 0);

  EXPECT_DOUBLE_EQ(S2::Angle(p000, pz, p045), M_PI_4);
  EXPECT_DOUBLE_EQ(S2::TurnAngle(p000, pz, p045), -3 * M_PI_4);

  EXPECT_DOUBLE_EQ(S2::Angle(p045, pz, p180), 3 * M_PI_4);
  EXPECT_DOUBLE_EQ(S2::TurnAngle(p045, pz, p180), -M_PI_4);

  EXPECT_DOUBLE_EQ(S2::Angle(p000, pz, p180), M_PI);
  EXPECT_DOUBLE_EQ(S2::TurnAngle(p000, pz, p180), 0);

  EXPECT_DOUBLE_EQ(S2::Angle(pz, p000, p045), M_PI_2);
  EXPECT_DOUBLE_EQ(S2::TurnAngle(pz, p000, p045), M_PI_2);

  EXPECT_DOUBLE_EQ(S2::Angle(pz, p000, pz), 0);
  EXPECT_DOUBLE_EQ(fabs(S2::TurnAngle(pz, p000, pz)), M_PI);
}

TEST(S2, AreaMethods) {
  S2Point pz(0, 0, 1);
  S2Point p000(1, 0, 0);
  S2Point p045 = S2Point(1, 1, 0).Normalize();
  S2Point p090(0, 1, 0);
  S2Point p180(-1, 0, 0);

  EXPECT_DOUBLE_EQ(S2::Area(p000, p090, pz), M_PI_2);
  EXPECT_DOUBLE_EQ(S2::Area(p045, pz, p180), 3 * M_PI_4);

  // Make sure that Area() has good *relative* accuracy even for
  // very small areas.
  static double const eps = 1e-10;
  S2Point pepsx = S2Point(eps, 0, 1).Normalize();
  S2Point pepsy = S2Point(0, eps, 1).Normalize();
  double expected1 = 0.5 * eps * eps;
  EXPECT_NEAR(S2::Area(pepsx, pepsy, pz), expected1, 1e-14 * expected1);

  // Make sure that it can handle degenerate triangles.
  S2Point pr = S2Point(0.257, -0.5723, 0.112).Normalize();
  S2Point pq = S2Point(-0.747, 0.401, 0.2235).Normalize();
  EXPECT_EQ(S2::Area(pr, pr, pr), 0);
  // The following test is not exact due to rounding error.
  EXPECT_NEAR(S2::Area(pr, pq, pr), 0, 1e-15);
  EXPECT_EQ(S2::Area(p000, p045, p090), 0);

  double max_girard = 0;
  for (int i = 0; i < 10000; ++i) {
    S2Point p0 = S2Testing::RandomPoint();
    S2Point d1 = S2Testing::RandomPoint();
    S2Point d2 = S2Testing::RandomPoint();
    S2Point p1 = (p0 + 1e-15 * d1).Normalize();
    S2Point p2 = (p0 + 1e-15 * d2).Normalize();
    // The actual displacement can be as much as 1.2e-15 due to roundoff.
    // This yields a maximum triangle area of about 0.7e-30.
    EXPECT_LE(S2::Area(p0, p1, p2), 0.7e-30);
    max_girard = max(max_girard, S2::GirardArea(p0, p1, p2));
  }
  // This check only passes if GirardArea() uses RobustCrossProd().
  LOG(INFO) << "Worst case Girard for triangle area 1e-30: " << max_girard;
  EXPECT_LE(max_girard, 1e-14);

  // Try a very long and skinny triangle.
  S2Point p045eps = S2Point(1, 1, eps).Normalize();
  double expected2 = 5.8578643762690495119753e-11;  // Mathematica.
  EXPECT_NEAR(S2::Area(p000, p045eps, p090), expected2, 1e-9 * expected2);

  // Triangles with near-180 degree edges that sum to a quarter-sphere.
  static double const eps2 = 1e-14;
  S2Point p000eps2 = S2Point(1, 0.1*eps2, eps2).Normalize();
  double quarter_area1 = S2::Area(p000eps2, p000, p045) +
                         S2::Area(p000eps2, p045, p180) +
                         S2::Area(p000eps2, p180, pz) +
                         S2::Area(p000eps2, pz, p000);
  EXPECT_DOUBLE_EQ(quarter_area1, M_PI);

  // Four other triangles that sum to a quarter-sphere.
  S2Point p045eps2 = S2Point(1, 1, eps2).Normalize();
  double quarter_area2 = S2::Area(p045eps2, p000, p045) +
                         S2::Area(p045eps2, p045, p180) +
                         S2::Area(p045eps2, p180, pz) +
                         S2::Area(p045eps2, pz, p000);
  EXPECT_DOUBLE_EQ(quarter_area2, M_PI);

  // Compute the area of a hemisphere using four triangles with one near-180
  // degree edge and one near-degenerate edge.
  for (int i = 0; i < 100; ++i) {
    double lng = 2 * M_PI * S2Testing::rnd.RandDouble();
    S2Point p0 = S2LatLng::FromRadians(1e-20, lng).Normalized().ToPoint();
    S2Point p1 = S2LatLng::FromRadians(0, lng).Normalized().ToPoint();
    double p2_lng = lng + S2Testing::rnd.RandDouble();
    S2Point p2 = S2LatLng::FromRadians(0, p2_lng).Normalized().ToPoint();
    S2Point p3 = S2LatLng::FromRadians(0, lng + M_PI).Normalized().ToPoint();
    S2Point p4 = S2LatLng::FromRadians(0, lng + 5.0).Normalized().ToPoint();
    double area = (S2::Area(p0, p1, p2) + S2::Area(p0, p2, p3) +
                   S2::Area(p0, p3, p4) + S2::Area(p0, p4, p1));
    EXPECT_NEAR(area, 2 * M_PI, 2e-15);
  }
}

TEST(S2, TrueCentroid) {
  // Test TrueCentroid() with very small triangles.  This test assumes that
  // the triangle is small enough so that it is nearly planar.
  for (int i = 0; i < 100; ++i) {
    Vector3_d p, x, y;
    S2Testing::GetRandomFrame(&p, &x, &y);
    double d = 1e-4 * pow(1e-4, S2Testing::rnd.RandDouble());
    S2Point p0 = (p - d * x).Normalize();
    S2Point p1 = (p + d * x).Normalize();
    S2Point p2 = (p + 3 * d * y).Normalize();
    S2Point centroid = S2::TrueCentroid(p0, p1, p2).Normalize();

    // The centroid of a planar triangle is at the intersection of its
    // medians, which is two-thirds of the way along each median.
    S2Point expected_centroid = (p + d * y).Normalize();
    EXPECT_LE(centroid.Angle(expected_centroid), 2e-8);
  }
}

TEST(Sign, CollinearPoints) {
  // The following points happen to be *exactly collinear* along a line that it
  // approximate tangent to the surface of the unit sphere.  In fact, C is the
  // exact midpoint of the line segment AB.  All of these points are close
  // enough to unit length to satisfy S2::IsUnitLength().
  S2Point a(0.72571927877036835, 0.46058825605889098, 0.51106749730504852);
  S2Point b(0.7257192746638208, 0.46058826573818168, 0.51106749441312738);
  S2Point c(0.72571927671709457, 0.46058826089853633, 0.51106749585908795);
  EXPECT_EQ(c - a, b - c);
  EXPECT_NE(0, S2::Sign(a, b, c));
  EXPECT_EQ(S2::Sign(a, b, c), S2::Sign(b, c, a));
  EXPECT_EQ(S2::Sign(a, b, c), -S2::Sign(c, b, a));

  // The points "x1" and "x2" are exactly proportional, i.e. they both lie
  // on a common line through the origin.  Both points are considered to be
  // normalized, and in fact they both satisfy (x == x.Normalize()).
  // Therefore the triangle (x1, x2, -x1) consists of three distinct points
  // that all lie on a common line through the origin.
  S2Point x1(0.99999999999999989, 1.4901161193847655e-08, 0);
  S2Point x2(1, 1.4901161193847656e-08, 0);
  EXPECT_EQ(x1, x1.Normalize());
  EXPECT_EQ(x2, x2.Normalize());
  EXPECT_NE(0, S2::Sign(x1, x2, -x1));
  EXPECT_EQ(S2::Sign(x1, x2, -x1), S2::Sign(x2, -x1, x1));
  EXPECT_EQ(S2::Sign(x1, x2, -x1), -S2::Sign(-x1, x2, x1));

  // Here are two more points that are distinct, exactly proportional, and
  // that satisfy (x == x.Normalize()).
  S2Point x3 = S2Point(1, 1, 1).Normalize();
  S2Point x4 = 0.99999999999999989 * x3;
  EXPECT_EQ(x3, x3.Normalize());
  EXPECT_EQ(x4, x4.Normalize());
  EXPECT_NE(x3, x4);
  EXPECT_NE(0, S2::Sign(x3, x4, -x3));

  // The following two points demonstrate that Normalize() is not idempotent,
  // i.e. y0.Normalize() != y0.Normalize().Normalize().  Both points satisfy
  // S2::IsNormalized(), though, and the two points are exactly proportional.
  S2Point y0 = S2Point(1, 1, 0);
  S2Point y1 = y0.Normalize();
  S2Point y2 = y1.Normalize();
  EXPECT_NE(y1, y2);
  EXPECT_EQ(y2, y2.Normalize());
  EXPECT_NE(0, S2::Sign(y1, y2, -y1));
  EXPECT_EQ(S2::Sign(y1, y2, -y1), S2::Sign(y2, -y1, y1));
  EXPECT_EQ(S2::Sign(y1, y2, -y1), -S2::Sign(-y1, y2, y1));
}

// Given 3 points A, B, C that are exactly coplanar with the origin and where
// A < B < C in lexicographic order, verify that ABC is counterclockwise (if
// expected == 1) or clockwise (if expected == -1) using S2::ExpensiveSign().
//
// This method is intended specifically for checking the cases where
// symbolic perturbations are needed to break ties.
static void CheckSymbolicSign(int expected, S2Point const& a,
                              S2Point const& b, S2Point const& c) {
  CHECK_LT(a, b);
  CHECK_LT(b, c);
  CHECK_EQ(0, a.DotProd(b.CrossProd(c)));

  // Use ASSERT rather than EXPECT to suppress spurious error messages.
  ASSERT_EQ(expected, S2::ExpensiveSign(a, b, c));
  ASSERT_EQ(expected, S2::ExpensiveSign(b, c, a));
  ASSERT_EQ(expected, S2::ExpensiveSign(c, a, b));
  ASSERT_EQ(-expected, S2::ExpensiveSign(c, b, a));
  ASSERT_EQ(-expected, S2::ExpensiveSign(b, a, c));
  ASSERT_EQ(-expected, S2::ExpensiveSign(a, c, b));
}

TEST(Sign, SymbolicPerturbationCodeCoverage) {
  // The purpose of this test is simply to get code coverage of
  // SymbolicallyPerturbedCCW().  Let M_1, M_2, ... be the sequence of
  // submatrices whose determinant sign is tested by that function.  Then the
  // i-th test below is a 3x3 matrix M (with rows A, B, C) such that:
  //
  //    det(M) = 0
  //    det(M_j) = 0 for j < i
  //    det(M_i) != 0
  //    A < B < C in lexicographic order.
  //
  // I checked that reversing the sign of any of the "return" statements in
  // SymbolicallyPerturbedCCW() will cause this test to fail.

  // det(M_1) = b0*c1 - b1*c0
  CheckSymbolicSign(1,
                    S2Point(-3, -1, 0), S2Point(-2, 1, 0), S2Point(1, -2, 0));

  // det(M_2) = b2*c0 - b0*c2
  CheckSymbolicSign(1,
                    S2Point(-6, 3, 3), S2Point(-4, 2, -1), S2Point(-2, 1, 4));

  // det(M_3) = b1*c2 - b2*c1
  CheckSymbolicSign(1, S2Point(0, -1, -1), S2Point(0, 1, -2), S2Point(0, 2, 1));
  // From this point onward, B or C must be zero, or B is proportional to C.

  // det(M_4) = c0*a1 - c1*a0
  CheckSymbolicSign(1, S2Point(-1, 2, 7), S2Point(2, 1, -4), S2Point(4, 2, -8));

  // det(M_5) = c0
  CheckSymbolicSign(1,
                    S2Point(-4, -2, 7), S2Point(2, 1, -4), S2Point(4, 2, -8));

  // det(M_6) = -c1
  CheckSymbolicSign(1, S2Point(0, -5, 7), S2Point(0, -4, 8), S2Point(0, -2, 4));

  // det(M_7) = c2*a0 - c0*a2
  CheckSymbolicSign(1,
                    S2Point(-5, -2, 7), S2Point(0, 0, -2), S2Point(0, 0, -1));

  // det(M_8) = c2
  CheckSymbolicSign(1, S2Point(0, -2, 7), S2Point(0, 0, 1), S2Point(0, 0, 2));
  // From this point onward, C must be zero.

  // det(M_9) = a0*b1 - a1*b0
  CheckSymbolicSign(1, S2Point(-3, 1, 7), S2Point(-1, -4, 1), S2Point(0, 0, 0));

  // det(M_10) = -b0
  CheckSymbolicSign(1,
                    S2Point(-6, -4, 7), S2Point(-3, -2, 1), S2Point(0, 0, 0));

  // det(M_11) = b1
  CheckSymbolicSign(-1, S2Point(0, -4, 7), S2Point(0, -2, 1), S2Point(0, 0, 0));

  // det(M_12) = a0
  CheckSymbolicSign(-1,
                    S2Point(-1, -4, 5), S2Point(0, 0, -3), S2Point(0, 0, 0));

  // det(M_13) = 1
  CheckSymbolicSign(1, S2Point(0, -4, 5), S2Point(0, 0, -5), S2Point(0, 0, 0));
}

// This test repeatedly constructs some number of points that are on or nearly
// on a given great circle.  Then it chooses one of these points as the
// "origin" and sorts the other points in CCW order around it.  Of course,
// since the origin is on the same great circle as the points being sorted,
// nearly all of these tests are degenerate.  It then does various consistency
// checks to verify that the points are indeed sorted in CCW order.
//
// It is easier to think about what this test is doing if you imagine that the
// points are in general position rather than on a great circle.
class SignTest : public testing::Test {
 protected:
  // The following method is used to sort a collection of points in CCW order
  // around a given origin.  It returns true if A comes before B in the CCW
  // ordering (starting at an arbitrary fixed direction).
  class LessCCW {
   public:
    LessCCW(S2Point const& origin, S2Point const& start)
        : origin_(origin), start_(start) {
    }
    bool operator()(S2Point const& a, S2Point const& b) const {
      // OrderedCCW() acts like "<=", so we need to invert the comparison.
      return !S2::OrderedCCW(start_, b, a, origin_);
    }
   private:
    S2Point const origin_;
    S2Point const start_;
  };

  // Given a set of points with no duplicates, first remove "origin" from
  // "points" (if it exists) and then sort the remaining points in CCW order
  // around "origin" putting the result in "sorted".
  static void SortCCW(vector<S2Point> const& points, S2Point const& origin,
                      vector<S2Point>* sorted) {
    // Make a copy of the points with "origin" removed.
    sorted->clear();
    std::remove_copy(points.begin(), points.end(), back_inserter(*sorted),
                     origin);

    // Sort the points CCW around the origin starting at (*sorted)[0].
    LessCCW less(origin, (*sorted)[0]);
    std::sort(sorted->begin(), sorted->end(), less);
  }

  // Given a set of points sorted circularly CCW around "origin", and the
  // index "start" of a point A, count the number of CCW triangles OAB over
  // all sorted points B not equal to A.  Also check that the results of the
  // CCW tests are consistent with the hypothesis that the points are sorted.
  static int CountCCW(vector<S2Point> const& sorted, S2Point const& origin,
                      int start) {
    int num_ccw = 0;
    int last_sign = 1;
    int const n = sorted.size();
    for (int j = 1; j < n; ++j) {
      int sign = S2::Sign(origin, sorted[start], sorted[(start + j) % n]);
      EXPECT_NE(0, sign);
      if (sign > 0) ++num_ccw;

      // Since the points are sorted around the origin, we expect to see a
      // (possibly empty) sequence of CCW triangles followed by a (possibly
      // empty) sequence of CW triangles.
      EXPECT_FALSE(sign > 0 && last_sign < 0);
      last_sign = sign;
    }
    return num_ccw;
  }

  // Test exhaustively whether the points in "sorted" are sorted circularly
  // CCW around "origin".
  static void TestCCW(vector<S2Point> const& sorted, S2Point const& origin) {
    int const n = sorted.size();
    int total_num_ccw = 0;
    int last_num_ccw = CountCCW(sorted, origin, n - 1);
    for (int start = 0; start < n; ++start) {
      int num_ccw = CountCCW(sorted, origin, start);
      // Each iteration we increase the start index by 1, therefore the number
      // of CCW triangles should decrease by at most 1.
      EXPECT_GE(num_ccw, last_num_ccw - 1);
      total_num_ccw += num_ccw;
      last_num_ccw = num_ccw;
    }
    // We have tested all triangles of the form OAB.  Exactly half of these
    // should be CCW.
    EXPECT_EQ(n * (n-1) / 2, total_num_ccw);
  }

  static void AddNormalized(S2Point const& a, vector<S2Point>* points) {
    points->push_back(a.Normalize());
  }

  // Add two points A1 and A2 that are slightly offset from A along the
  // tangent toward B, and such that A, A1, and A2 are exactly collinear
  // (i.e. even with infinite-precision arithmetic).
  static void AddTangentPoints(S2Point const& a, S2Point const& b,
                               vector<S2Point>* points) {
    Vector3_d dir = S2::RobustCrossProd(a, b).CrossProd(a).Normalize();
    if (dir == S2Point(0, 0, 0)) return;
    for (;;) {
      S2Point delta = 1e-15 * S2Testing::rnd.RandDouble() * dir;
      if ((a + delta) != a && (a + delta) - a == a - (a - delta) &&
          S2::IsUnitLength(a + delta) && S2::IsUnitLength(a - delta)) {
        points->push_back(a + delta);
        points->push_back(a - delta);
        return;
      }
    }
  }

  // Add zero or more (but usually one) point that is likely to trigger
  // Sign() degeneracies among the given points.
  static void AddDegeneracy(vector<S2Point>* points) {
    S2Testing::Random* rnd = &S2Testing::rnd;
    S2Point a = (*points)[rnd->Uniform(points->size())];
    S2Point b = (*points)[rnd->Uniform(points->size())];
    int coord = rnd->Uniform(3);
    switch (rnd->Uniform(8)) {
      case 0:
        // Add a random point (not uniformly distributed) along the great
        // circle AB.
        AddNormalized(rnd->UniformDouble(-1, 1) * a +
                      rnd->UniformDouble(-1, 1) * b, points);
        break;
      case 1:
        // Perturb one coordinate by the minimum amount possible.
        a[coord] = nextafter(a[coord], rnd->OneIn(2) ? 2 : -2);
        AddNormalized(a, points);
        break;
      case 2:
        // Perturb one coordinate by up to 1e-15.
        a[coord] += 1e-15 * rnd->UniformDouble(-1, 1);
        AddNormalized(a, points);
        break;
      case 3:
        // Scale a point just enough so that it is different while still being
        // considered normalized.
        a *= rnd->OneIn(2) ? (1 + 2e-16) : (1 - 1e-16);
        if (S2::IsUnitLength(a)) points->push_back(a);
        break;
      case 4: {
        // Add the intersection point of AB with X=0, Y=0, or Z=0.
        S2Point dir(0, 0, 0);
        dir[coord] = rnd->OneIn(2) ? 1 : -1;
        Vector3_d norm = S2::RobustCrossProd(a, b).Normalize();
        if (norm.Norm2() > 0) {
          AddNormalized(S2::RobustCrossProd(dir, norm), points);
        }
        break;
      }
      case 5:
        // Add two closely spaced points along the tangent at A to the great
        // circle through AB.
        AddTangentPoints(a, b, points);
        break;
      case 6:
        // Add two closely spaced points along the tangent at A to the great
        // circle through A and the X-axis.
        AddTangentPoints(a, S2Point(1, 0, 0), points);
        break;
      case 7:
        // Add the negative of a point.
        points->push_back(-a);
        break;
    }
  }

  // Sort the points around the given origin, and then do some consistency
  // checks to verify that they are actually sorted.
  static void SortAndTest(vector<S2Point> const& points,
                          S2Point const& origin) {
    vector<S2Point> sorted;
    SortCCW(points, origin, &sorted);
    TestCCW(sorted, origin);
  }

  // Construct approximately "n" points near the great circle through A and B,
  // then sort them and test whether they are sorted.
  static void TestGreatCircle(S2Point a, S2Point b, int n) {
    a = a.Normalize();
    b = b.Normalize();
    vector<S2Point> points;
    points.push_back(a);
    points.push_back(b);
    while (points.size() < n) {
      AddDegeneracy(&points);
    }
    // Remove any (0, 0, 0) points that were accidentically created, then sort
    // the points and remove duplicates.
    points.erase(std::remove(points.begin(), points.end(), S2Point(0, 0, 0)),
                 points.end());
    std::sort(points.begin(), points.end());
    points.erase(std::unique(points.begin(), points.end()), points.end());
    EXPECT_GE(points.size(), n / 2);

    SortAndTest(points, a);
    SortAndTest(points, b);
    for (int k = 0; k < points.size(); ++k) {
      SortAndTest(points, points[k]);
    }
  }
};

TEST_F(SignTest, StressTest) {
  // The run time of this test is *cubic* in the parameter below.
  static int const kNumPointsPerCircle = 20;

  // This test is randomized, so it is beneficial to run it several times.
  for (int iter = 0; iter < 3; ++iter) {
    // The most difficult great circles are the ones in the X-Y, Y-Z, and X-Z
    // planes, for two reasons.  First, when one or more coordinates are close
    // to zero then the perturbations can be much smaller, since floating
    // point numbers are spaced much more closely together near zero.  (This
    // tests the handling of things like underflow.)  The second reason is
    // that most of the cases of SymbolicallyPerturbedCCW() can only be
    // reached when one or more input point coordinates are zero.
    TestGreatCircle(S2Point(1, 0, 0), S2Point(0, 1, 0), kNumPointsPerCircle);
    TestGreatCircle(S2Point(1, 0, 0), S2Point(0, 0, 1), kNumPointsPerCircle);
    TestGreatCircle(S2Point(0, -1, 0), S2Point(0, 0, 1), kNumPointsPerCircle);

    // This tests a great circle where at least some points have X, Y, and Z
    // coordinates with exactly the same mantissa.  One useful property of
    // such points is that when they are scaled (e.g. multiplying by 1+eps),
    // all such points are exactly collinear with the origin.
    TestGreatCircle(S2Point(1 << 25, 1, -8), S2Point(-4, -(1 << 20), 1),
                    kNumPointsPerCircle);
  }
}

class StableSignTest : public testing::Test {
 protected:
  // Estimate the probability that S2::StableSign() will not be able to compute
  // the determinant sign of a triangle A, B, C consisting of three points
  // that are as collinear as possible and spaced the given distance apart.
  double GetFailureRate(double km) {
    int const kIters = 1000;
    int failure_count = 0;
    double m = tan(S2Testing::KmToAngle(km).radians());
    for (int iter = 0; iter < kIters; ++iter) {
      S2Point a, x, y;
      S2Testing::GetRandomFrame(&a, &x, &y);
      S2Point b = (a - m * x).Normalize();
      S2Point c = (a + m * x).Normalize();
      int sign = S2::StableSign(a, b, c);
      if (sign != 0) {
        EXPECT_EQ(S2::ExactSign(a, b, c), sign);
      } else {
        ++failure_count;
      }
    }
    double rate = static_cast<double>(failure_count) / kIters;
    LOG(INFO) << "StableSign failure rate for " << km << " km = " << rate;
    return rate;
  }
};

TEST_F(StableSignTest, FailureRate) {
  // Verify that StableSign() is able to handle most cases where the three
  // points are as collinear as possible.  (For reference, TriageSign() fails
  // virtually 100% of the time on this test.)
  //
  // Note that the failure rate *decreases* as the points get closer together,
  // and the decrease is approximately linear.  For example, the failure rate
  // is 0.4% for collinear points spaced 1km apart, but only 0.0004% for
  // collinear points spaced 1 meter apart.

  EXPECT_LT(GetFailureRate(1.0), 0.01);  //  1km spacing: <  1% (actual 0.4%)
  EXPECT_LT(GetFailureRate(10.0), 0.1);  // 10km spacing: < 10% (actual 4%)
}

// Note: obviously, I could have defined a bundle of metrics like this in the
// S2 class itself rather than just for testing.  However, it's not clear that
// this is useful other than for testing purposes, and I find
// S2::kMinWidth.GetMaxLevel(width) to be slightly more readable than
// than S2::kWidth.min().GetMaxLevel(width).  Also, there is no fundamental
// reason that we need to analyze the minimum, maximum, and average values of
// every metric; it would be perfectly reasonable to just define one of these.

template<int dim>
class MetricBundle {
 public:
  using Metric = S2::Metric<dim>;
  MetricBundle(Metric const& min, Metric const& max, Metric const& avg) :
    min_(min), max_(max), avg_(avg) {}
  Metric const& min_;
  Metric const& max_;
  Metric const& avg_;

 private:
  MetricBundle(MetricBundle const&) = delete;
  void operator=(MetricBundle const&) = delete;
};

template<int dim>
static void CheckMinMaxAvg(MetricBundle<dim> const& bundle) {
  EXPECT_LE(bundle.min_.deriv(), bundle.avg_.deriv());
  EXPECT_LE(bundle.avg_.deriv(), bundle.max_.deriv());
}

template<int dim>
static void CheckLessOrEqual(MetricBundle<dim> const& a,
                             MetricBundle<dim> const& b) {
  EXPECT_LE(a.min_.deriv(), b.min_.deriv());
  EXPECT_LE(a.max_.deriv(), b.max_.deriv());
  EXPECT_LE(a.avg_.deriv(), b.avg_.deriv());
}

TEST(S2, Metrics) {
  MetricBundle<1> angle_span(S2::kMinAngleSpan, S2::kMaxAngleSpan,
                             S2::kAvgAngleSpan);
  MetricBundle<1> width(S2::kMinWidth, S2::kMaxWidth, S2::kAvgWidth);
  MetricBundle<1> edge(S2::kMinEdge, S2::kMaxEdge, S2::kAvgEdge);
  MetricBundle<1> diag(S2::kMinDiag, S2::kMaxDiag, S2::kAvgDiag);
  MetricBundle<2> area(S2::kMinArea, S2::kMaxArea, S2::kAvgArea);

  // First, check that min <= avg <= max for each metric.
  CheckMinMaxAvg(angle_span);
  CheckMinMaxAvg(width);
  CheckMinMaxAvg(edge);
  CheckMinMaxAvg(diag);
  CheckMinMaxAvg(area);

  // Check that the maximum aspect ratio of an individual cell is consistent
  // with the global minimums and maximums.
  EXPECT_GE(S2::kMaxEdgeAspect, 1);
  EXPECT_LE(S2::kMaxEdgeAspect, S2::kMaxEdge.deriv() / S2::kMinEdge.deriv());
  EXPECT_GE(S2::kMaxDiagAspect, 1);
  EXPECT_LE(S2::kMaxDiagAspect, S2::kMaxDiag.deriv() / S2::kMinDiag.deriv());

  // Check various conditions that are provable mathematically.
  CheckLessOrEqual(width, angle_span);
  CheckLessOrEqual(width, edge);
  CheckLessOrEqual(edge, diag);

  EXPECT_GE(S2::kMinArea.deriv(),
            S2::kMinWidth.deriv() * S2::kMinEdge.deriv() - 1e-15);
  EXPECT_LE(S2::kMaxArea.deriv(),
            S2::kMaxWidth.deriv() * S2::kMaxEdge.deriv() + 1e-15);

  // GetMinLevelForLength() and friends have built-in assertions, we just need
  // to call these functions to test them.
  //
  // We don't actually check that the metrics are correct here, e.g. that
  // GetMinWidth(10) is a lower bound on the width of cells at level 10.
  // It is easier to check these properties in s2cell_unittest, since
  // S2Cell has methods to compute the cell vertices, etc.

  for (int level = -2; level <= S2CellId::kMaxLevel + 3; ++level) {
    double width = S2::kMinWidth.deriv() * pow(2, -level);
    if (level >= S2CellId::kMaxLevel + 3) width = 0;

    // Check boundary cases (exactly equal to a threshold value).
    int expected_level = max(0, min(S2CellId::kMaxLevel, level));
    EXPECT_EQ(S2::kMinWidth.GetMinLevel(width), expected_level);
    EXPECT_EQ(S2::kMinWidth.GetMaxLevel(width), expected_level);
    EXPECT_EQ(S2::kMinWidth.GetClosestLevel(width), expected_level);

    // Also check non-boundary cases.
    EXPECT_EQ(S2::kMinWidth.GetMinLevel(1.2 * width), expected_level);
    EXPECT_EQ(S2::kMinWidth.GetMaxLevel(0.8 * width), expected_level);
    EXPECT_EQ(S2::kMinWidth.GetClosestLevel(1.2 * width), expected_level);
    EXPECT_EQ(S2::kMinWidth.GetClosestLevel(0.8 * width), expected_level);

    // Same thing for area.
    double area = S2::kMinArea.deriv() * pow(4, -level);
    if (level <= -3) area = 0;
    EXPECT_EQ(S2::kMinArea.GetMinLevel(area), expected_level);
    EXPECT_EQ(S2::kMinArea.GetMaxLevel(area), expected_level);
    EXPECT_EQ(S2::kMinArea.GetClosestLevel(area), expected_level);
    EXPECT_EQ(S2::kMinArea.GetMinLevel(1.2 * area), expected_level);
    EXPECT_EQ(S2::kMinArea.GetMaxLevel(0.8 * area), expected_level);
    EXPECT_EQ(S2::kMinArea.GetClosestLevel(1.2 * area), expected_level);
    EXPECT_EQ(S2::kMinArea.GetClosestLevel(0.8 * area), expected_level);
  }
}

TEST(S2, Frames) {
  Matrix3x3_d m;
  S2Point z = S2Point(0.2, 0.5, -3.3).Normalize();
  S2::GetFrame(z, &m);
  EXPECT_TRUE(S2::ApproxEquals(m.Col(2), z));
  EXPECT_TRUE(S2::IsUnitLength(m.Col(0)));
  EXPECT_TRUE(S2::IsUnitLength(m.Col(1)));
  EXPECT_DOUBLE_EQ(m.Det(), 1);

  EXPECT_TRUE(S2::ApproxEquals(S2::ToFrame(m, m.Col(0)), S2Point(1, 0, 0)));
  EXPECT_TRUE(S2::ApproxEquals(S2::ToFrame(m, m.Col(1)), S2Point(0, 1, 0)));
  EXPECT_TRUE(S2::ApproxEquals(S2::ToFrame(m, m.Col(2)), S2Point(0, 0, 1)));

  EXPECT_TRUE(S2::ApproxEquals(S2::FromFrame(m, S2Point(1, 0, 0)), m.Col(0)));
  EXPECT_TRUE(S2::ApproxEquals(S2::FromFrame(m, S2Point(0, 1, 0)), m.Col(1)));
  EXPECT_TRUE(S2::ApproxEquals(S2::FromFrame(m, S2Point(0, 0, 1)), m.Col(2)));
}

static void TestRotate(S2Point const& p, S2Point const& axis, S1Angle angle) {
  S2Point result = S2::Rotate(p, axis, angle);

  // "result" should be unit length.
  EXPECT_TRUE(S2::IsUnitLength(result));

  // "result" and "p" should be the same distance from "axis".
  double kMaxPositionError = 1e-15;
  EXPECT_LE((S1Angle(result, axis) - S1Angle(p, axis)).abs().radians(),
            kMaxPositionError);

  // Check that the rotation angle is correct.  We allow a fixed error in the
  // *position* of the result, so we need to convert this into a rotation
  // angle.  The allowable error can be very large as "p" approaches "axis".
  double axis_distance = p.CrossProd(axis).Norm();
  double max_rotation_error;
  if (axis_distance < kMaxPositionError) {
    max_rotation_error = 2 * M_PI;
  } else {
    max_rotation_error = asin(kMaxPositionError / axis_distance);
  }
  double actual_rotation = S2::TurnAngle(p, axis, result) + M_PI;
  double rotation_error = remainder(angle.radians() - actual_rotation,
                                    2 * M_PI);
  EXPECT_LE(rotation_error, max_rotation_error);
}

TEST(S2, Rotate) {
  for (int iter = 0; iter < 1000; ++iter) {
    S2Point axis = S2Testing::RandomPoint();
    S2Point target = S2Testing::RandomPoint();
    // Choose a distance whose logarithm is uniformly distributed.
    double distance = M_PI * pow(1e-15, S2Testing::rnd.RandDouble());
    // Sometimes choose points near the far side of the axis.
    if (S2Testing::rnd.OneIn(5)) distance = M_PI - distance;
    S2Point p = S2EdgeUtil::InterpolateAtDistance(S1Angle::Radians(distance),
                                                  axis, target);
    // Choose the rotation angle.
    double angle = 2 * M_PI * pow(1e-15, S2Testing::rnd.RandDouble());
    if (S2Testing::rnd.OneIn(3)) angle = -angle;
    if (S2Testing::rnd.OneIn(10)) angle = 0;
    TestRotate(p, axis, S1Angle::Radians(angle));
  }
}

TEST(S2, S2PointHashSpreads) {
  int kTestPoints = 1 << 16;
  unordered_set<size_t> set;
  unordered_set<S2Point, S2PointHash> points;
  S2PointHash hasher;
  S2Point base = S2Point(1, 1, 1);
  for (int i = 0; i < kTestPoints; ++i) {
    // All points in a tiny cap to test avalanche property of hash
    // function (the cap would be of radius 1mm on Earth (4*10^9/2^35).
    S2Point perturbed = base + S2Testing::RandomPoint() / (1ULL << 35);
    perturbed = perturbed.Normalize();
    set.insert(hasher(perturbed));
    points.insert(perturbed);
  }
  // A real collision is extremely unlikely.
  EXPECT_EQ(0, kTestPoints - points.size());
  // Allow a few for the hash.
  EXPECT_GE(10, kTestPoints - set.size());
}

// Given a point P, return the minimum level at which an edge of some S2Cell
// parent of P is nearly collinear with S2::Origin().  This is the minimum
// level for which Sign() may need to resort to expensive calculations in
// order to determine which side of an edge the origin lies on.
static int GetMinExpensiveLevel(S2Point const& p) {
  S2CellId id = S2CellId::FromPoint(p);
  for (int level = 0; level <= S2CellId::kMaxLevel; ++level) {
    S2Cell cell(id.parent(level));
    for (int k = 0; k < 4; ++k) {
      S2Point a = cell.GetVertex(k);
      S2Point b = cell.GetVertex((k + 1) & 3);
      if (S2::TriageSign(a, b, S2::Origin(), a.CrossProd(b)) == 0) {
        return level;
      }
    }
  }
  return S2CellId::kMaxLevel + 1;
}

TEST(S2, OriginTest) {
  // To minimize the number of expensive Sign() calculations,
  // S2::Origin() should not be nearly collinear with any commonly used edges.
  // Two important categories of such edges are:
  //
  //  - edges along a line of longitude (reasonably common geographically)
  //  - S2Cell edges (used extensively when computing S2Cell coverings)
  //
  // This implies that the origin:
  //
  //  - should not be too close to either pole (since all lines of longitude
  //    converge at the poles)
  //  - should not be colinear with edges of any S2Cell except for very small
  //    ones (which are used less frequently)
  //
  // The point chosen below is about 66km from the north pole towards the East
  // Siberian Sea.  The purpose of the STtoUV(2/3) calculation is to keep the
  // origin as far away as possible from the longitudinal edges of large
  // S2Cells.  (The line of longitude through the chosen point is always 1/3
  // or 2/3 of the way across any S2Cell with longitudinal edges that it
  // passes through.)

  EXPECT_EQ(S2Point(-0.01, 0.01 * S2::STtoUV(2./3), 1).Normalize(),
            S2::Origin());

  // Check that the origin is not too close to either pole.  (We don't use
  // S2Earth because we don't want to depend on that package.)
  double distance_km = acos(S2::Origin().z()) * S2Testing::kEarthRadiusKm;
  EXPECT_GE(distance_km, 50.0);
  LOG(INFO) << "\nS2::Origin() coordinates: " << S2LatLng(S2::Origin())
            << ", distance from pole: " << distance_km << " km";

  // Check that S2::Origin() is not collinear with the edges of any large
  // S2Cell.  We do this is two parts.  For S2Cells that belong to either
  // polar face, we simply need to check that S2::Origin() is not nearly
  // collinear with any edge of any cell that contains it (except for small
  // cells < 3 meters across).
  EXPECT_GE(GetMinExpensiveLevel(S2::Origin()), 22);

  // For S2Cells that belong to the four non-polar faces, only longitudinal
  // edges can possibly be colinear with S2::Origin().  We check these edges
  // by projecting S2::Origin() onto the equator, and then testing all S2Cells
  // that contain this point to make sure that none of their edges are nearly
  // colinear with S2::Origin() (except for small cells < 3 meters across).
  S2Point equator_point(S2::Origin().x(), S2::Origin().y(), 0);
  EXPECT_GE(GetMinExpensiveLevel(equator_point), 22);
}
