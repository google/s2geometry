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
//

#include "s2polygon.h"

#include <math.h>
#include <stdlib.h>
#include <algorithm>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "base/macros.h"
#include "base/stringprintf.h"
#include "strings/serialize.h"
#include <gtest/gtest.h>
#include "util/coding/coder.h"
#include "r1interval.h"
#include "s1angle.h"
#include "s2.h"
#include "s2cap.h"
#include "s2cell.h"
#include "s2cellid.h"
#include "s2cellunion.h"
#include "s2edgeutil.h"
#include "s2error.h"
#include "s2latlng.h"
#include "s2loop.h"
#include "s2polygonbuilder.h"
#include "s2polyline.h"
#include "s2regioncoverer.h"
#include "s2testing.h"
#include "s2textformat.h"
#include "util/gtl/fixedarray.h"
#include "util/gtl/stl_util.h"
#include "util/math/matrix3x3.h"

using std::max;
using std::min;
using std::numeric_limits;
using std::pair;
using std::swap;
using std::unique_ptr;
using std::vector;


// A set of nested loops around the point 0:0 (lat:lng).
// Every vertex of kNear0 is a vertex of kNear1.
char const kNearPoint[] = "0:0";
string const kNear0 = "-1:0, 0:1, 1:0, 0:-1;";
string const kNear1 = "-1:-1, -1:0, -1:1, 0:1, 1:1, 1:0, 1:-1, 0:-1;";
string const kNear2 = "-1:-2, -2:5, 5:-2;";
string const kNear3 = "-2:-2, -3:6, 6:-3;";
string const kNearHemi = "0:-90, -90:0, 0:90, 90:0;";

// A set of nested loops around the point 0:180 (lat:lng).
// Every vertex of kFar0 and kFar2 belongs to kFar1, and all
// the loops except kFar2 are non-convex.
string const kFar0 = "0:179, 1:180, 0:-179, 2:-180;";
string const kFar1 =
  "0:179, -1:179, 1:180, -1:-179, 0:-179, 3:-178, 2:-180, 3:178;";
string const kFar2 = "3:-178, 3:178, -1:179, -1:-179;";
string const kFar3 = "-3:-178, 4:-177, 4:177, -3:178, -2:179;";
string const kFarHemi = "0:-90, 60:90, -60:90;";

// A set of nested loops around the point -90:0 (lat:lng).
string const kSouthPoint = "-89.9999:0.001";
string const kSouth0a = "-90:0, -89.99:0.01, -89.99:0;";
string const kSouth0b = "-90:0, -89.99:0.03, -89.99:0.02;";
string const kSouth0c = "-90:0, -89.99:0.05, -89.99:0.04;";
string const kSouth1 = "-90:0, -89.9:0.1, -89.9:-0.1;";
string const kSouth2 = "-90:0, -89.8:0.2, -89.8:-0.2;";
string const kSouthHemi = "0:-180, 0:60, 0:-60;";

// Two different loops that surround all the Near and Far loops except
// for the hemispheres.
string const kNearFar1 = "-1:-9, -9:-9, -9:9, 9:9, 9:-9, 1:-9, "
                         "1:-175, 9:-175, 9:175, -9:175, -9:-175, -1:-175;";
string const kNearFar2 = "-2:15, -2:170, -8:-175, 8:-175, "
                         "2:170, 2:15, 8:-4, -8:-4;";

// Loops that result from intersection of other loops.
string const kFarHSouthH = "0:-180, 0:90, -60:90, 0:-90;";

// Rectangles that form a cross, with only shared vertices, no crossing edges.
// Optional holes outside the intersecting region.
string const kCross1 = "-2:1, -1:1, 1:1, 2:1, 2:-1, 1:-1, -1:-1, -2:-1;";
string const kCross1SideHole = "-1.5:0.5, -1.2:0.5, -1.2:-0.5, -1.5:-0.5;";
string const kCross2 = "1:-2, 1:-1, 1:1, 1:2, -1:2, -1:1, -1:-1, -1:-2;";
string const kCross2SideHole = "0.5:-1.5, 0.5:-1.2, -0.5:-1.2, -0.5:-1.5;";
string const kCrossCenterHole = "-0.5:0.5, 0.5:0.5, 0.5:-0.5, -0.5:-0.5;";

// Two rectangles that intersect, but no edges cross and there's always
// local containment (rather than crossing) at each shared vertex.
// In this ugly ASCII art, 1 is A+B, 2 is B+C:
//      +---+---+---+
//      | A | B | C |
//      +---+---+---+
string const kOverlap1 = "0:1, 1:1, 2:1, 2:0, 1:0, 0:0;";
string const kOverlap1SideHole = "0.2:0.8, 0.8:0.8, 0.8:0.2, 0.2:0.2;";
string const kOverlap2 = "1:1, 2:1, 3:1, 3:0, 2:0, 1:0;";
string const kOverlap2SideHole = "2.2:0.8, 2.8:0.8, 2.8:0.2, 2.2:0.2;";
string const kOverlapCenterHole = "1.2:0.8, 1.8:0.8, 1.8:0.2, 1.2:0.2;";

// An empty polygon.
string const kEmpty = "";
// By symmetry, the intersection of the two polygons has almost half the area
// of either polygon.
string const kOverlap3 = "-10:10, 0:10, 0:-10, -10:-10, -10:0";
string const kOverlap4 = "-10:0, 10:0, 10:-10, -10:-10";

class S2PolygonTestBase : public testing::Test {
 public:
  S2PolygonTestBase();
  ~S2PolygonTestBase();

 protected:
  // Some standard polygons to use in the tests.
  S2Polygon const* const empty;
  S2Polygon const* const full;
  S2Polygon const* const near_0;
  S2Polygon const* const near_10;
  S2Polygon const* const near_30;
  S2Polygon const* const near_32;
  S2Polygon const* const near_3210;
  S2Polygon const* const near_H3210;

  S2Polygon const* const far_10;
  S2Polygon const* const far_21;
  S2Polygon const* const far_321;
  S2Polygon const* const far_H20;
  S2Polygon const* const far_H3210;

  S2Polygon const* const south_0ab;
  S2Polygon const* const south_2;
  S2Polygon const* const south_210b;
  S2Polygon const* const south_H21;
  S2Polygon const* const south_H20abc;

  S2Polygon const* const nf1_n10_f2_s10abc;

  S2Polygon const* const nf2_n2_f210_s210ab;

  S2Polygon const* const f32_n0;
  S2Polygon const* const n32_s0b;

  S2Polygon const* const cross1;
  S2Polygon const* const cross1_side_hole;
  S2Polygon const* const cross1_center_hole;
  S2Polygon const* const cross2;
  S2Polygon const* const cross2_side_hole;
  S2Polygon const* const cross2_center_hole;

  S2Polygon const* const overlap1;
  S2Polygon const* const overlap1_side_hole;
  S2Polygon const* const overlap1_center_hole;
  S2Polygon const* const overlap2;
  S2Polygon const* const overlap2_side_hole;
  S2Polygon const* const overlap2_center_hole;

  S2Polygon const* const far_H;
  S2Polygon const* const south_H;
  S2Polygon const* const far_H_south_H;
};

static bool TestEncodeDecode(const S2Polygon* src) {
  Encoder encoder;
  src->Encode(&encoder);
  Decoder decoder(encoder.base(), encoder.length());
  S2Polygon dst;
  dst.Decode(&decoder);
  return src->Equals(&dst);
}

static S2Polygon* MakePolygon(string const& str) {
  unique_ptr<S2Polygon> polygon(s2textformat::MakeVerbatimPolygon(str));

#if 0
  // TODO(ericv): Ensure that S2PolygonBuilder2 passes this test.
  // Check that InitToSnapped() is idempotent.
  S2Polygon snapped1, snapped2;
  snapped1.InitToSnapped(polygon.get());
  snapped2.InitToSnapped(&snapped1);
  EXPECT_TRUE(snapped1.Equals(&snapped2));
#endif

  // Check that Decode(Encode(x)) is the identity function.
  EXPECT_TRUE(TestEncodeDecode(polygon.get()));
  return polygon.release();
}

static void CheckContains(string const& a_str, string const& b_str) {
  S2Polygon* a = MakePolygon(a_str);
  S2Polygon* b = MakePolygon(b_str);
  unique_ptr<S2Polygon> delete_a(a);
  unique_ptr<S2Polygon> delete_b(b);
  EXPECT_TRUE(a->Contains(b));
  EXPECT_TRUE(a->ApproxContains(b, S1Angle::Radians(1e-15)));
  EXPECT_FALSE(a->ApproxDisjoint(b, S1Angle::Radians(1e-15)));
}

static void CheckContainsPoint(string const& a_str, string const& b_str) {
  unique_ptr<S2Polygon> a(s2textformat::MakePolygon(a_str));
  EXPECT_TRUE(a->VirtualContainsPoint(s2textformat::MakePoint(b_str)))
    << " " << a_str << " did not contain " << b_str;
}

TEST(S2Polygon, Init) {
  CheckContains(kNear1, kNear0);
  CheckContains(kNear2, kNear1);
  CheckContains(kNear3, kNear2);
  CheckContains(kNearHemi, kNear3);
  CheckContains(kFar1, kFar0);
  CheckContains(kFar2, kFar1);
  CheckContains(kFar3, kFar2);
  CheckContains(kFarHemi, kFar3);
  CheckContains(kSouth1, kSouth0a);
  CheckContains(kSouth1, kSouth0b);
  CheckContains(kSouth1, kSouth0c);
  CheckContains(kSouthHemi, kSouth2);
  CheckContains(kNearFar1, kNear3);
  CheckContains(kNearFar1, kFar3);
  CheckContains(kNearFar2, kNear3);
  CheckContains(kNearFar2, kFar3);

  CheckContainsPoint(kNear0, kNearPoint);
  CheckContainsPoint(kNear1, kNearPoint);
  CheckContainsPoint(kNear2, kNearPoint);
  CheckContainsPoint(kNear3, kNearPoint);
  CheckContainsPoint(kNearHemi, kNearPoint);
  CheckContainsPoint(kSouth0a, kSouthPoint);
  CheckContainsPoint(kSouth1, kSouthPoint);
  CheckContainsPoint(kSouth2, kSouthPoint);
  CheckContainsPoint(kSouthHemi, kSouthPoint);
}

TEST(S2Polygon, OverlapFractions) {
  unique_ptr<S2Polygon> a(MakePolygon(kEmpty));
  unique_ptr<S2Polygon> b(MakePolygon(kEmpty));
  auto result = S2Polygon::GetOverlapFractions(a.get(), b.get());
  EXPECT_DOUBLE_EQ(1.0, result.first);
  EXPECT_DOUBLE_EQ(1.0, result.second);

  b.reset(MakePolygon(kOverlap3));
  result = S2Polygon::GetOverlapFractions(a.get(), b.get());
  EXPECT_DOUBLE_EQ(1.0, result.first);
  EXPECT_DOUBLE_EQ(0.0, result.second);

  a.reset(MakePolygon(kOverlap4));
  result = S2Polygon::GetOverlapFractions(a.get(), b.get());
  EXPECT_NEAR(0.5, result.first, 1e-14);
  EXPECT_NEAR(0.5, result.second, 1e-14);
}

TEST(S2Polygon, OriginNearPole) {
  // S2Polygon operations are more efficient if S2::Origin() is near a pole.
  // (Loops that contain a pole tend to have very loose bounding boxes because
  // they span the full longitude range.  S2Polygon canonicalizes all loops so
  // that they don't contain S2::Origin(), thus by placing S2::Origin() near a
  // pole we minimize the number of canonical loops which contain that pole.)
  EXPECT_GE(S2LatLng::Latitude(S2::Origin()).degrees(), 80);
}

S2PolygonTestBase::S2PolygonTestBase() :
    empty(new S2Polygon()),
    full(MakePolygon("full")),
    near_0(MakePolygon(kNear0)),
    near_10(MakePolygon(kNear0 + kNear1)),
    near_30(MakePolygon(kNear3 + kNear0)),
    near_32(MakePolygon(kNear2 + kNear3)),
    near_3210(MakePolygon(kNear0 + kNear2 + kNear3 + kNear1)),
    near_H3210(MakePolygon(kNear0 + kNear2 + kNear3 + kNearHemi + kNear1)),

    far_10(MakePolygon(kFar0 + kFar1)),
    far_21(MakePolygon(kFar2 + kFar1)),
    far_321(MakePolygon(kFar2 + kFar3 + kFar1)),
    far_H20(MakePolygon(kFar2 + kFarHemi + kFar0)),
    far_H3210(MakePolygon(kFar2 + kFarHemi + kFar0 + kFar1 + kFar3)),

    south_0ab(MakePolygon(kSouth0a + kSouth0b)),
    south_2(MakePolygon(kSouth2)),
    south_210b(MakePolygon(kSouth2 + kSouth0b + kSouth1)),
    south_H21(MakePolygon(kSouth2 + kSouthHemi + kSouth1)),
    south_H20abc(MakePolygon(kSouth2 + kSouth0b + kSouthHemi +
                             kSouth0a + kSouth0c)),

    nf1_n10_f2_s10abc(MakePolygon(kSouth0c + kFar2 + kNear1 + kNearFar1 +
                                  kNear0 + kSouth1 + kSouth0b + kSouth0a)),

    nf2_n2_f210_s210ab(MakePolygon(kFar2 + kSouth0a + kFar1 + kSouth1 + kFar0 +
                                   kSouth0b + kNearFar2 + kSouth2 + kNear2)),

    f32_n0(MakePolygon(kFar2 + kNear0 + kFar3)),
    n32_s0b(MakePolygon(kNear3 + kSouth0b + kNear2)),

    cross1(MakePolygon(kCross1)),
    cross1_side_hole(MakePolygon(kCross1 + kCross1SideHole)),
    cross1_center_hole(MakePolygon(kCross1 + kCrossCenterHole)),
    cross2(MakePolygon(kCross2)),
    cross2_side_hole(MakePolygon(kCross2 + kCross2SideHole)),
    cross2_center_hole(MakePolygon(kCross2 + kCrossCenterHole)),

    overlap1(MakePolygon(kOverlap1)),
    overlap1_side_hole(MakePolygon(kOverlap1 + kOverlap1SideHole)),
    overlap1_center_hole(MakePolygon(kOverlap1 + kOverlapCenterHole)),
    overlap2(MakePolygon(kOverlap2)),
    overlap2_side_hole(MakePolygon(kOverlap2 + kOverlap2SideHole)),
    overlap2_center_hole(MakePolygon(kOverlap2 + kOverlapCenterHole)),

    far_H(MakePolygon(kFarHemi)),
    south_H(MakePolygon(kSouthHemi)),
    far_H_south_H(MakePolygon(kFarHSouthH)) {
}

S2PolygonTestBase::~S2PolygonTestBase() {
  delete empty;
  delete full;
  delete near_0;
  delete near_10;
  delete near_30;
  delete near_32;
  delete near_3210;
  delete near_H3210;

  delete far_10;
  delete far_21;
  delete far_321;
  delete far_H20;
  delete far_H3210;

  delete south_0ab;
  delete south_2;
  delete south_210b;
  delete south_H21;
  delete south_H20abc;

  delete nf1_n10_f2_s10abc;

  delete nf2_n2_f210_s210ab;

  delete f32_n0;
  delete n32_s0b;

  delete cross1;
  delete cross1_side_hole;
  delete cross1_center_hole;
  delete cross2;
  delete cross2_side_hole;
  delete cross2_center_hole;

  delete overlap1;
  delete overlap1_side_hole;
  delete overlap1_center_hole;
  delete overlap2;
  delete overlap2_side_hole;
  delete overlap2_center_hole;

  delete far_H;
  delete south_H;
  delete far_H_south_H;
}

static void CheckEqual(S2Polygon const* a, S2Polygon const* b,
                       double max_error = 0) {
  if (a->BoundaryApproxEquals(b, max_error)) return;
  S2PolygonBuilder builder(S2PolygonBuilderOptions::DIRECTED_XOR());
  S2Polygon a2, b2;
  builder.AddPolygon(a);
  ASSERT_TRUE(builder.AssemblePolygon(&a2, NULL));
  builder.AddPolygon(b);
  ASSERT_TRUE(builder.AssemblePolygon(&b2, NULL));
  EXPECT_TRUE(a2.BoundaryApproxEquals(&b2, max_error));
}

static void CheckComplementary(S2Polygon const* a, S2Polygon const* b) {
  S2Polygon b1;
  b1.InitToComplement(b);
  CheckEqual(a, &b1);
}

TEST(S2Polygon, TestApproxContainsAndDisjoint) {
  // We repeatedly choose a random cell id and intersect its bounding polygon
  // "A" with the bounding polygon "B" of one its child cells.  The result may
  // not be contained by either A or B, because the vertices of B near the
  // edge midpoints of A may be slightly outside A, and even when the crossing
  // edges are intersected, the intersection point may also be slightly
  // outside A and/or B.
  //
  // We repeat the test many times and expect that some fraction of the exact
  // tests should fail, while all of the approximate test should succeed.
  int const kIters = 1000;
  int exact_contains = 0, exact_disjoint = 0;
  for (int iter = 0; iter < kIters; ++iter) {
    S2CellId id = S2Testing::GetRandomCellId(10);
    S2Polygon parent_polygon((S2Cell(id)));
    S2Polygon child_polygon(S2Cell(id.child(0)));

    // Get the intersection.  There is no guarantee that the intersection will
    // be contained by A or B.  Similarly, the intersection may slightly
    // overlap an adjacent disjoint polygon C.
    S2Polygon intersection;
    intersection.InitToIntersection(&parent_polygon, &child_polygon);
    if (parent_polygon.Contains(&intersection)) {
      ++exact_contains;
    }
    EXPECT_TRUE(parent_polygon.ApproxContains(
        &intersection, S2EdgeUtil::kIntersectionTolerance));

    S2Polygon adjacent_polygon(S2Cell(id.child(1)));
    if (!adjacent_polygon.Intersects(&intersection)) {
      ++exact_disjoint;
    }
    EXPECT_TRUE(adjacent_polygon.ApproxDisjoint(
        &intersection, S2EdgeUtil::kIntersectionTolerance));
  }
  // Back-of-the-envelope calculations show that about 56% of the exact
  // containment tests should succeed, and about 75% of the exact disjoint
  // tests should succeed.  These are close to the measured values.
  EXPECT_LT(exact_contains, 0.60 * kIters);  // 51.8% succeed
  EXPECT_LT(exact_disjoint, 0.86 * kIters);  // 78.6% succeed
}

// Given a pair of polygons where A contains B, check that various identities
// involving union, intersection, and difference operations hold true.
static void TestOneNestedPair(S2Polygon const* a, S2Polygon const* b) {
  EXPECT_TRUE(a->Contains(b));
  EXPECT_EQ(!b->is_empty(), a->Intersects(b));
  EXPECT_EQ(!b->is_empty(), b->Intersects(a));

  S2Polygon c, d, e;
  c.InitToUnion(a, b);
  CheckEqual(&c, a);

  d.InitToIntersection(a, b);
  CheckEqual(&d, b);

  e.InitToDifference(b, a);
  EXPECT_TRUE(e.is_empty());
}

// Given a pair of disjoint polygons A and B, check that various identities
// involving union, intersection, and difference operations hold true.
static void TestOneDisjointPair(S2Polygon const* a, S2Polygon const* b) {
  EXPECT_FALSE(a->Intersects(b));
  EXPECT_FALSE(b->Intersects(a));
  EXPECT_EQ(b->is_empty(), a->Contains(b));
  EXPECT_EQ(a->is_empty(), b->Contains(a));

  S2Polygon ab, c, d, e, f;
  S2PolygonBuilder builder(S2PolygonBuilderOptions::DIRECTED_XOR());
  builder.AddPolygon(a);
  builder.AddPolygon(b);
  ASSERT_TRUE(builder.AssemblePolygon(&ab, NULL));

  c.InitToUnion(a, b);
  CheckEqual(&c, &ab);

  d.InitToIntersection(a, b);
  EXPECT_TRUE(d.is_empty());

  e.InitToDifference(a, b);
  CheckEqual(&e, a);

  f.InitToDifference(b, a);
  CheckEqual(&f, b);
}

// Given polygons A and B whose union covers the sphere, check that various
// identities involving union, intersection, and difference hold true.
static void TestOneCoveringPair(S2Polygon const* a, S2Polygon const* b) {
  EXPECT_EQ(a->is_full(), a->Contains(b));
  EXPECT_EQ(b->is_full(), b->Contains(a));

  S2Polygon c, d, e, f;
  c.InitToUnion(a, b);
  EXPECT_TRUE(c.is_full());
}

// Given polygons A and B such that both A and its complement intersect both B
// and its complement, check that various identities involving union,
// intersection, and difference hold true.
static void TestOneOverlappingPair(S2Polygon const* a, S2Polygon const* b) {
  EXPECT_FALSE(a->Contains(b));
  EXPECT_FALSE(b->Contains(a));
  EXPECT_TRUE(a->Intersects(b));

  S2Polygon c, d, e;
  c.InitToUnion(a, b);
  EXPECT_FALSE(c.is_full());

  d.InitToIntersection(a, b);
  EXPECT_FALSE(d.is_empty());

  e.InitToDifference(b, a);
  EXPECT_FALSE(e.is_empty());
}

// Given a pair of polygons where A contains B, test various identities
// involving A, B, and their complements.
static void TestNestedPair(S2Polygon const* a, S2Polygon const* b) {
  S2Polygon a1, b1;
  a1.InitToComplement(a);
  b1.InitToComplement(b);

  TestOneNestedPair(a, b);
  TestOneNestedPair(&b1, &a1);
  TestOneDisjointPair(&a1, b);
  TestOneCoveringPair(a, &b1);
}

// Given a pair of disjoint polygons A and B, test various identities
// involving A, B, and their complements.
static void TestDisjointPair(S2Polygon const* a, S2Polygon const* b) {
  S2Polygon a1, b1;
  a1.InitToComplement(a);
  b1.InitToComplement(b);

  TestOneDisjointPair(a, b);
  TestOneCoveringPair(&a1, &b1);
  TestOneNestedPair(&a1, b);
  TestOneNestedPair(&b1, a);
}

// Given polygons A and B such that both A and its complement intersect both B
// and its complement, test various identities involving these four polygons.
static void TestOverlappingPair(S2Polygon const* a, S2Polygon const* b) {
  S2Polygon a1, b1;
  a1.InitToComplement(a);
  b1.InitToComplement(b);

  TestOneOverlappingPair(a, b);
  TestOneOverlappingPair(&a1, &b1);
  TestOneOverlappingPair(&a1, b);
  TestOneOverlappingPair(a, &b1);
}

// "a1" is the complement of "a", and "b1" is the complement of "b".
static void TestOneComplementPair(S2Polygon const* a, S2Polygon const* a1,
                                  S2Polygon const* b, S2Polygon const* b1) {
  // Check DeMorgan's Law and that subtraction is the same as intersection
  // with the complement.  This function is called multiple times in order to
  // test the various combinations of complements.

  S2Polygon a1_or_b, a_and_b1, a_minus_b;
  a_and_b1.InitToIntersection(a, b1);
  a1_or_b.InitToUnion(a1, b);
  a_minus_b.InitToDifference(a, b);

  CheckComplementary(&a1_or_b, &a_and_b1);
  CheckEqual(&a_minus_b, &a_and_b1);
}

// Test identities that should hold for any pair of polygons A, B and their
// complements.
static void TestComplements(S2Polygon const* a, S2Polygon const* b) {
  S2Polygon a1, b1;
  a1.InitToComplement(a);
  b1.InitToComplement(b);

  TestOneComplementPair(a, &a1, b, &b1);
  TestOneComplementPair(&a1, a, b, &b1);
  TestOneComplementPair(a, &a1, &b1, b);
  TestOneComplementPair(&a1, a, &b1, b);
}

static void TestDestructiveUnion(S2Polygon const* a, S2Polygon const* b) {
  S2Polygon c;
  c.InitToUnion(a, b);
  vector<S2Polygon*> polygons;
  polygons.push_back(a->Clone());
  polygons.push_back(b->Clone());
  unique_ptr<S2Polygon> c_destructive(S2Polygon::DestructiveUnion(&polygons));
  CheckEqual(&c, c_destructive.get());
}

static void TestRelationWithDesc(S2Polygon const* a, S2Polygon const* b,
                                 bool contains, bool contained,
                                 bool intersects, const char* description) {
  SCOPED_TRACE(description);
  EXPECT_EQ(contains, a->Contains(b));
  EXPECT_EQ(contained, b->Contains(a));
  EXPECT_EQ(intersects, a->Intersects(b));
  if (contains) TestNestedPair(a, b);
  if (contained) TestNestedPair(b, a);
  if (!intersects) TestDisjointPair(a, b);
  if (intersects && !(contains | contained)) {
    TestOverlappingPair(a, b);  // See TestOverlappingPair for definition
  }
  TestDestructiveUnion(a, b);
  TestComplements(a, b);
}

TEST_F(S2PolygonTestBase, Relations) {
#define TestRelation(a, b, contains, contained, intersects)   \
  TestRelationWithDesc(a, b, contains, contained, intersects, \
                       "args " #a ", " #b)

  TestRelation(near_10, empty, true, false, false);
  TestRelation(near_10, near_10, true, true, true);
  TestRelation(full, near_10, true, false, true);
  TestRelation(near_10, near_30, false, true, true);
  TestRelation(near_10, near_32, false, false, false);
  TestRelation(near_10, near_3210, false, true, true);
  TestRelation(near_10, near_H3210, false, false, false);
  TestRelation(near_30, near_32, true, false, true);
  TestRelation(near_30, near_3210, true, false, true);
  TestRelation(near_30, near_H3210, false, false, true);
  TestRelation(near_32, near_3210, false, true, true);
  TestRelation(near_32, near_H3210, false, false, false);
  TestRelation(near_3210, near_H3210, false, false, false);

  TestRelation(far_10, far_21, false, false, false);
  TestRelation(far_10, far_321, false, true, true);
  TestRelation(far_10, far_H20, false, false, false);
  TestRelation(far_10, far_H3210, false, false, false);
  TestRelation(far_21, far_321, false, false, false);
  TestRelation(far_21, far_H20, false, false, false);
  TestRelation(far_21, far_H3210, false, true, true);
  TestRelation(far_321, far_H20, false, false, true);
  TestRelation(far_321, far_H3210, false, false, true);
  TestRelation(far_H20, far_H3210, false, false, true);

  TestRelation(south_0ab, south_2, false, true, true);
  TestRelation(south_0ab, south_210b, false, false, true);
  TestRelation(south_0ab, south_H21, false, true, true);
  TestRelation(south_0ab, south_H20abc, false, true, true);
  TestRelation(south_2, south_210b, true, false, true);
  TestRelation(south_2, south_H21, false, false, true);
  TestRelation(south_2, south_H20abc, false, false, true);
  TestRelation(south_210b, south_H21, false, false, true);
  TestRelation(south_210b, south_H20abc, false, false, true);
  TestRelation(south_H21, south_H20abc, true, false, true);

  TestRelation(nf1_n10_f2_s10abc, nf2_n2_f210_s210ab, false, false, true);
  TestRelation(nf1_n10_f2_s10abc, near_32, true, false, true);
  TestRelation(nf1_n10_f2_s10abc, far_21, false, false, false);
  TestRelation(nf1_n10_f2_s10abc, south_0ab, false, false, false);
  TestRelation(nf1_n10_f2_s10abc, f32_n0, true, false, true);

  TestRelation(nf2_n2_f210_s210ab, near_10, false, false, false);
  TestRelation(nf2_n2_f210_s210ab, far_10, true, false, true);
  TestRelation(nf2_n2_f210_s210ab, south_210b, true, false, true);
  TestRelation(nf2_n2_f210_s210ab, south_0ab, true, false, true);
  TestRelation(nf2_n2_f210_s210ab, n32_s0b, true, false, true);

  TestRelation(cross1, cross2, false, false, true);
  TestRelation(cross1_side_hole, cross2, false, false, true);
  TestRelation(cross1_center_hole, cross2, false, false, true);
  TestRelation(cross1, cross2_side_hole, false, false, true);
  TestRelation(cross1, cross2_center_hole, false, false, true);
  TestRelation(cross1_side_hole, cross2_side_hole, false, false, true);
  TestRelation(cross1_center_hole, cross2_side_hole, false, false, true);
  TestRelation(cross1_side_hole, cross2_center_hole, false, false, true);
  TestRelation(cross1_center_hole, cross2_center_hole, false, false, true);

  // These cases, when either polygon has a hole, test a different code path
  // from the other cases.
  TestRelation(overlap1, overlap2, false, false, true);
  TestRelation(overlap1_side_hole, overlap2, false, false, true);
  TestRelation(overlap1_center_hole, overlap2, false, false, true);
  TestRelation(overlap1, overlap2_side_hole, false, false, true);
  TestRelation(overlap1, overlap2_center_hole, false, false, true);
  TestRelation(overlap1_side_hole, overlap2_side_hole, false, false, true);
  TestRelation(overlap1_center_hole, overlap2_side_hole, false, false, true);
  TestRelation(overlap1_side_hole, overlap2_center_hole, false, false, true);
  TestRelation(overlap1_center_hole, overlap2_center_hole, false, false, true);
#undef TestRelation
}

TEST_F(S2PolygonTestBase, EmptyAndFull) {
  EXPECT_TRUE(empty->is_empty());
  EXPECT_FALSE(full->is_empty());
  EXPECT_FALSE(empty->is_full());
  EXPECT_TRUE(full->is_full());

  TestNestedPair(empty, empty);
  TestNestedPair(full, empty);
  TestNestedPair(full, full);
}

struct TestCase {
  char const* a;
  char const* b;
  char const* a_and_b;
  char const* a_or_b;
  char const* a_minus_b;
};

TestCase test_cases[] = {
  // Two triangles that share an edge.
  { "4:2, 3:1, 3:3;",

    "3:1, 2:2, 3:3;",

    "",

    "4:2, 3:1, 2:2, 3:3;",

    "4:2, 3:1, 3:3;"
  },

  // Two vertical bars and a horizontal bar connecting them.
  { "0:0, 0:2, 3:2, 3:0;   0:3, 0:5, 3:5, 3:3;",

    "1:1, 1:4, 2:4, 2:1;",

    "1:1, 1:2, 2:2, 2:1;   1:3, 1:4, 2:4, 2:3;",

    "0:0, 0:2, 1:2, 1:3, 0:3, 0:5, 3:5, 3:3, 2:3, 2:2, 3:2, 3:0;",

    "0:0, 0:2, 1:2, 1:1, 2:1, 2:2, 3:2, 3:0;   "
    "0:3, 0:5, 3:5, 3:3, 2:3, 2:4, 1:4, 1:3;"
  },

  // Two vertical bars and two horizontal bars centered around S2::Origin().
  { "1:88, 1:93, 2:93, 2:88;   -1:88, -1:93, 0:93, 0:88;",

    "-2:89, -2:90, 3:90, 3:89;   -2:91, -2:92, 3:92, 3:91;",

    "1:89, 1:90, 2:90, 2:89;   1:91, 1:92, 2:92, 2:91;   "
    "-1:89, -1:90, 0:90, 0:89;   -1:91, -1:92, 0:92, 0:91;",

    "-1:88, -1:89, -2:89, -2:90, -1:90, -1:91, -2:91, -2:92, -1:92, -1:93,"
    "0:93, 0:92, 1:92, 1:93, 2:93, 2:92, 3:92, 3:91, 2:91, 2:90, 3:90,"
    "3:89, 2:89, 2:88, 1:88, 1:89, 0:89, 0:88;   "
    "0:90, 0:91, 1:91, 1:90;",

    "1:88, 1:89, 2:89, 2:88;   1:90, 1:91, 2:91, 2:90;   "
    "1:92, 1:93, 2:93, 2:92;   -1:88, -1:89, 0:89, 0:88;   "
    "-1:90, -1:91, 0:91, 0:90;   -1:92, -1:93, 0:93, 0:92;"
  },

  // Two interlocking square doughnuts centered around -S2::Origin().
  { "-1:-93, -1:-89, 3:-89, 3:-93;   0:-92, 0:-90, 2:-90, 2:-92;",

    "-3:-91, -3:-87, 1:-87, 1:-91;   -2:-90, -2:-88, 0:-88, 0:-90;",

    "-1:-91, -1:-90, 0:-90, 0:-91;   0:-90, 0:-89, 1:-89, 1:-90;",

    "-1:-93, -1:-91, -3:-91, -3:-87, 1:-87, 1:-89, 3:-89, 3:-93;   "
    "0:-92, 0:-91, 1:-91, 1:-90, 2:-90, 2:-92;   "
    "-2:-90, -2:-88, 0:-88, 0:-89, -1:-89, -1:-90;",

    "-1:-93, -1:-91, 0:-91, 0:-92, 2:-92, 2:-90, 1:-90, 1:-89, 3:-89, 3:-93;   "
    "-1:-90, -1:-89, 0:-89, 0:-90;"
  },

  // An incredibly thin triangle intersecting a square, such that the two
  // intersection points of the triangle with the square are identical.
  // This results in a degenerate loop that needs to be handled correctly.
  { "10:44, 10:46, 12:46, 12:44;",

    "11:45, 89:45.00000000000001, 90:45;",

    "",  // Empty intersection!

    // Original square with extra vertex, and triangle disappears (due to
    // default vertex_merge_radius of S2EdgeUtil::kIntersectionTolerance).
    "10:44, 10:46, 12:46, 12:45, 12:44;",

    "10:44, 10:46, 12:46, 12:45, 12:44;"
  },
};

TEST_F(S2PolygonTestBase, Operations) {
  S2Polygon far_south;
  far_south.InitToIntersection(far_H, south_H);
  CheckEqual(&far_south, far_H_south_H, 1e-15);

  for (int i = 0; i < arraysize(test_cases); ++i) {
    SCOPED_TRACE(StringPrintf("Polygon operation test case %d", i));
    TestCase* test = test_cases + i;
    unique_ptr<S2Polygon> a(MakePolygon(test->a));
    unique_ptr<S2Polygon> b(MakePolygon(test->b));
    unique_ptr<S2Polygon> expected_a_and_b(MakePolygon(test->a_and_b));
    unique_ptr<S2Polygon> expected_a_or_b(MakePolygon(test->a_or_b));
    unique_ptr<S2Polygon> expected_a_minus_b(MakePolygon(test->a_minus_b));

    // The intersections in the "expected" data were computed in lat-lng
    // space, while the actual intersections are computed using geodesics.
    // The error due to this depends on the length and direction of the line
    // segment being intersected, and how close the intersection is to the
    // endpoints of the segment.  The worst case is for a line segment between
    // two points at the same latitude, where the intersection point is in the
    // middle of the segment.  In this case the error is approximately
    // (p * t^2) / 8, where "p" is the absolute latitude in radians, "t" is
    // the longitude difference in radians, and both "p" and "t" are small.
    // The test cases all have small latitude and longitude differences.
    // If "p" and "t" are converted to degrees, the following error bound is
    // valid as long as (p * t^2 < 150).

    static double const kMaxError = 1e-4;

    S2Polygon a_and_b, a_or_b, a_minus_b;
    a_and_b.InitToIntersection(a.get(), b.get());
    CheckEqual(&a_and_b, expected_a_and_b.get(), kMaxError);
    a_or_b.InitToUnion(a.get(), b.get());
    CheckEqual(&a_or_b, expected_a_or_b.get(), kMaxError);
    TestDestructiveUnion(a.get(), b.get());
    a_minus_b.InitToDifference(a.get(), b.get());
    CheckEqual(&a_minus_b, expected_a_minus_b.get(), kMaxError);
  }
}

void ClearPolylineVector(vector<S2Polyline*>* polylines) {
  for (vector<S2Polyline*>::const_iterator it = polylines->begin();
       it != polylines->end(); ++it) {
    delete *it;
  }
  polylines->clear();
}

static void PolylineIntersectionSharedEdgeTest(const S2Polygon *p,
                                               int start_vertex,
                                               int direction) {
  SCOPED_TRACE(StringPrintf("Polyline intersection shared edge test "
                            " start=%d direction=%d",
                            start_vertex, direction));
  vector<S2Point> points;
  points.push_back(p->loop(0)->vertex(start_vertex));
  points.push_back(p->loop(0)->vertex(start_vertex + direction));
  S2Polyline polyline(points);
  vector<S2Polyline*> polylines;
  if (direction < 0) {
    p->IntersectWithPolyline(&polyline, &polylines);
    EXPECT_EQ(0, polylines.size());
    ClearPolylineVector(&polylines);
    p->SubtractFromPolyline(&polyline, &polylines);
    ASSERT_EQ(1, polylines.size());
    ASSERT_EQ(2, polylines[0]->num_vertices());
    EXPECT_EQ(points[0], polylines[0]->vertex(0));
    EXPECT_EQ(points[1], polylines[0]->vertex(1));
  } else {
    p->IntersectWithPolyline(&polyline, &polylines);
    ASSERT_EQ(1, polylines.size());
    ASSERT_EQ(2, polylines[0]->num_vertices());
    EXPECT_EQ(points[0], polylines[0]->vertex(0));
    EXPECT_EQ(points[1], polylines[0]->vertex(1));
    ClearPolylineVector(&polylines);
    p->SubtractFromPolyline(&polyline, &polylines);
    EXPECT_EQ(0, polylines.size());
  }
  ClearPolylineVector(&polylines);
}

// This tests polyline-polyline intersections.
// It covers the same edge cases as TestOperations and also adds some
// extra tests for shared edges.
TEST_F(S2PolygonTestBase, PolylineIntersection) {
  for (int v = 0; v < 3; ++v) {
    PolylineIntersectionSharedEdgeTest(cross1, v, 1);
    PolylineIntersectionSharedEdgeTest(cross1, v + 1, -1);
    PolylineIntersectionSharedEdgeTest(cross1_side_hole, v, 1);
    PolylineIntersectionSharedEdgeTest(cross1_side_hole, v + 1, -1);
  }

  // See comments in TestOperations about the vlue of this constant.
  static double const kMaxError = 1e-4;

  // This duplicates some of the tests in TestOperations by
  // converting the outline of polygon A to a polyline then intersecting
  // it with the polygon B. It then converts B to a polyline and intersects
  // it with A. It then feeds all of the results into a polygon builder and
  // tests that the output is equal to doing an intersection between A and B.
  for (int i = 0; i < arraysize(test_cases); ++i) {
    SCOPED_TRACE(StringPrintf("Polyline intersection test case %d", i));
    TestCase* test = test_cases + i;
    unique_ptr<S2Polygon> a(MakePolygon(test->a));
    unique_ptr<S2Polygon> b(MakePolygon(test->b));
    unique_ptr<S2Polygon> expected_a_and_b(MakePolygon(test->a_and_b));

    vector<S2Point> points;
    vector<S2Polyline *> polylines;
    for (int ab = 0; ab < 2; ab++) {
      S2Polygon *tmp = ab ? a.get() : b.get();
      S2Polygon *tmp2 = ab ? b.get() : a.get();
      for (int l = 0; l < tmp->num_loops(); l++) {
        points.clear();
        if (tmp->loop(l)->is_hole()) {
          for (int v = tmp->loop(l)->num_vertices(); v >=0 ; v--) {
            points.push_back(tmp->loop(l)->vertex(v));
          }
        } else {
          for (int v = 0; v <= tmp->loop(l)->num_vertices(); v++) {
            points.push_back(tmp->loop(l)->vertex(v));
          }
        }
        S2Polyline polyline(points);
        vector<S2Polyline *> tmp;
        tmp2->IntersectWithPolyline(&polyline, &tmp);
        polylines.insert(polylines.end(), tmp.begin(), tmp.end());
      }
    }

    S2PolygonBuilder builder(S2PolygonBuilderOptions::DIRECTED_XOR());
    for (int i = 0; i < polylines.size(); i++) {
      for (int j = 0; j < polylines[i]->num_vertices() - 1; j++) {
        builder.AddEdge(polylines[i]->vertex(j), polylines[i]->vertex(j + 1));
        VLOG(3) << " ... Adding edge: " << polylines[i]->vertex(j) << " - " <<
            polylines[i]->vertex(j + 1);
      }
    }
    ClearPolylineVector(&polylines);

    S2Polygon a_and_b;
    ASSERT_TRUE(builder.AssemblePolygon(&a_and_b, NULL));
    CheckEqual(&a_and_b, expected_a_and_b.get(), kMaxError);
  }
}

static void CheckCoveringIsConservative(S2Polygon const& polygon,
                                        vector<S2CellId> const& cells) {
  // Check that Contains(S2Cell) and MayIntersect(S2Cell) are implemented
  // conservatively, by comparing against the Contains/Intersect result with
  // the "cell polygon" defined by the four cell vertices.  Please note that
  // the cell polygon is *not* an exact representation of the S2Cell: cell
  // vertices are rounded from their true mathematical positions, which leads
  // to tiny cracks and overlaps between the cell polygons at different cell
  // levels.  That is why Contains(S2Cell) and MayIntersect(S2Cell) cannot be
  // implemented by simply converting the cell to an S2Polygon.  But it is
  // still useful to do this as a sanity check.  In particular:
  //
  //  - If Contains(cell) is true, the polygon must contain the cell polygon.
  //  - If the polygon intersects the cell polygon, then MayIntersect(cell)
  //    must return true.
  //
  for (int i = 0; i < cells.size(); ++i) {
    S2Cell cell(cells[i]);
    S2Polygon cell_poly(cell);
    if (polygon.Contains(cell)) {
      EXPECT_TRUE(polygon.Contains(&cell_poly));
    }
    if (polygon.Intersects(&cell_poly)) {
      EXPECT_TRUE(polygon.MayIntersect(cell));
    }
  }
}

// Remove a random polygon from "pieces" and return it.
static S2Polygon* ChoosePiece(vector<S2Polygon*> *pieces) {
  int i = S2Testing::rnd.Uniform(pieces->size());
  S2Polygon* result = (*pieces)[i];
  pieces->erase(pieces->begin() + i);
  return result;
}

static void SplitAndAssemble(S2Polygon const* polygon) {
  S2PolygonBuilder builder(S2PolygonBuilderOptions::DIRECTED_XOR());
  S2Polygon expected;
  builder.AddPolygon(polygon);
  ASSERT_TRUE(builder.AssemblePolygon(&expected, NULL));

  for (int iter = 0; iter < 10; ++iter) {
    S2RegionCoverer coverer;
    // Compute the minimum level such that the polygon's bounding
    // cap is guaranteed to be cut.
    double diameter = 2 * polygon->GetCapBound().GetRadius().radians();
    int min_level = S2::kMaxWidth.GetMinLevel(diameter);

    // TODO(ericv): Choose a level that will have up to 256 cells in the
    // covering.
    int level = min_level + S2Testing::rnd.Uniform(4);
    coverer.set_min_level(min_level);
    coverer.set_max_level(level);
    coverer.set_max_cells(500);

    vector<S2CellId> cells;
    coverer.GetCovering(*polygon, &cells);
    S2CellUnion covering;
    covering.Init(cells);
    S2Testing::CheckCovering(*polygon, covering, false);
    CheckCoveringIsConservative(*polygon, cells);
    VLOG(2) << cells.size() << " cells in covering";
    vector<S2Polygon*> pieces;
    for (int i = 0; i < cells.size(); ++i) {
      S2Cell cell(cells[i]);
      S2Polygon window(cell);
      S2Polygon* piece = new S2Polygon;
      piece->InitToIntersection(polygon, &window);
      pieces.push_back(piece);
      VLOG(4) << "\nPiece " << i << ":\n  Window: "
              << s2textformat::ToString(&window)
              << "\n  Piece: " << s2textformat::ToString(piece);
    }

    // Now we repeatedly remove two random pieces, compute their union, and
    // insert the result as a new piece until only one piece is left.
    //
    // We don't use S2Polygon::DestructiveUnion() because it joins the pieces
    // in a mostly deterministic order.  We don't just call random_shuffle()
    // on the pieces and repeatedly join the last two pieces in the vector
    // because this always joins a single original piece to the current union
    // rather than doing the unions according to a random tree structure.
    while (pieces.size() > 1) {
      unique_ptr<S2Polygon> a(ChoosePiece(&pieces));
      unique_ptr<S2Polygon> b(ChoosePiece(&pieces));
      S2Polygon* c = new S2Polygon;
      c->InitToUnion(a.get(), b.get());
      pieces.push_back(c);
      VLOG(4) << "\nJoining piece a: " << s2textformat::ToString(a.get())
              << "\n  With piece b: " << s2textformat::ToString(b.get())
              << "\n  To get piece c: " << s2textformat::ToString(c);
    }
    unique_ptr<S2Polygon> result(pieces[0]);
    pieces.pop_back();

    // The moment of truth!
    EXPECT_TRUE(expected.BoundaryNear(result.get()))
        << "\nActual:\n" << s2textformat::ToString(result.get())
        << "\nExpected:\n" << s2textformat::ToString(&expected);
  }
}

TEST_F(S2PolygonTestBase, Splitting) {
  // It takes too long to test all the polygons in debug mode, so we just pick
  // out some of the more interesting ones.

  SplitAndAssemble(near_10);
  SplitAndAssemble(near_H3210);
  SplitAndAssemble(far_H3210);
  SplitAndAssemble(south_0ab);
  SplitAndAssemble(south_210b);
  SplitAndAssemble(south_H20abc);
  SplitAndAssemble(nf1_n10_f2_s10abc);
  SplitAndAssemble(nf2_n2_f210_s210ab);
  SplitAndAssemble(far_H);
  SplitAndAssemble(south_H);
  SplitAndAssemble(far_H_south_H);
}

TEST(S2Polygon, InitToCellUnionBorder) {
  // Test S2Polygon::InitToCellUnionBorder().
  // The main thing to check is that adjacent cells of different sizes get
  // merged correctly.  To do this we generate two random adjacent cells,
  // convert to polygon, and make sure the polygon only has a single loop.
  for (int iter = 0; iter < 500; ++iter) {
    SCOPED_TRACE(StringPrintf("Iteration %d", iter));

    // Choose a random non-leaf cell.
    S2CellId big_cell =
        S2Testing::GetRandomCellId(S2Testing::rnd.Uniform(S2CellId::kMaxLevel));
    // Get all neighbors at some smaller level.
    int small_level = big_cell.level() +
        S2Testing::rnd.Uniform(min(16, S2CellId::kMaxLevel - big_cell.level()));
    vector<S2CellId> neighbors;
    big_cell.AppendAllNeighbors(small_level, &neighbors);
    // Pick one at random.
    S2CellId small_cell = neighbors[S2Testing::rnd.Uniform(neighbors.size())];
    // If it's diagonally adjacent, bail out.
    S2CellId edge_neighbors[4];
    big_cell.GetEdgeNeighbors(edge_neighbors);
    bool diagonal = true;
    for (int i = 0; i < 4; ++i) {
      if (edge_neighbors[i].contains(small_cell)) {
        diagonal = false;
      }
    }
    VLOG(3) << iter << ": big_cell " << big_cell <<
        " small_cell " << small_cell;
    if (diagonal) {
      VLOG(3) << "  diagonal - bailing out!";
      continue;
    }

    vector<S2CellId> cells;
    cells.push_back(big_cell);
    cells.push_back(small_cell);
    S2CellUnion cell_union;
    cell_union.Init(cells);
    EXPECT_EQ(2, cell_union.num_cells());
    S2Polygon poly;
    poly.InitToCellUnionBorder(cell_union);
    EXPECT_EQ(1, poly.num_loops());
    // If the conversion were perfect we could test containment, but due to
    // rounding the polygon won't always exactly contain both cells.  We can
    // at least test intersection.
    EXPECT_TRUE(poly.MayIntersect(S2Cell(big_cell)));
    EXPECT_TRUE(poly.MayIntersect(S2Cell(small_cell)));
  }
}

TEST(S2Polygon, UnionWithAmbgiuousCrossings) {
  S2Point const a_data[] = {
    S2Point(0.044856812877680216, -0.80679210859571904, 0.5891301722422051),
    S2Point(0.044851868273159699, -0.80679240802900054, 0.5891301386444033),
    S2Point(0.044854246527738666, -0.80679240292188514, 0.58912996457145106)
  };
  S2Point const b_data[3] = {
    S2Point(0.044849715793028468, -0.80679253837178111, 0.58913012401412856),
    S2Point(0.044855344598821352, -0.80679219751320641, 0.589130162266992),
    S2Point(0.044854017712818696, -0.80679210327223405, 0.58913039235179754)
  };
  vector<S2Point> a_vertices(a_data, a_data + arraysize(a_data));
  vector<S2Point> b_vertices(b_data, b_data + arraysize(b_data));
  vector<S2Loop*> loops;
  loops.push_back(new S2Loop(a_vertices));
  S2Polygon a(&loops);
  loops.push_back(new S2Loop(b_vertices));
  S2Polygon b(&loops);
  S2Polygon c;
  c.InitToUnion(&a, &b);
  EXPECT_FALSE(c.is_empty());
}

TEST(S2Polygon, InitToSloppySupportsEmptyPolygons) {
  S2Polygon empty_polygon;
  S2Polygon polygon;
  polygon.InitToSnapped(&empty_polygon);
  // InitToSloppy is further tested by SnapSplitsPolygon.
}

TEST(S2Polygon, InitToSnappedDoesNotRotateVertices) {
  // This particular example came from MapFacts, but in fact InitToSnapped
  // used to cyclically rotate the vertices of all "hole" loops.
  unique_ptr<S2Polygon> polygon(s2textformat::MakePolygon(
      "49.9305505:-124.8345463, 49.9307448:-124.8299657, "
      "49.9332101:-124.8301996, 49.9331224:-124.8341368; "
      "49.9311087:-124.8327042, 49.9318176:-124.8312621, "
      "49.9318866:-124.8334451"));
  S2Polygon polygon2, polygon3;
  polygon2.InitToSnapped(polygon.get());

  // Check that the first vertex is the same when converted to E7.
  EXPECT_EQ(S2LatLng::Latitude(polygon->loop(0)->vertex(0)).e7(),
            S2LatLng::Latitude(polygon2.loop(0)->vertex(0)).e7());
  EXPECT_EQ(S2LatLng::Longitude(polygon->loop(0)->vertex(0)).e7(),
            S2LatLng::Longitude(polygon2.loop(0)->vertex(0)).e7());

  // Check that snapping twice doesn't rotate the vertices.
  polygon3.InitToSnapped(&polygon2);
  EXPECT_TRUE(polygon2.Equals(&polygon3));
}

TEST(S2Polygon, InitToSnappedWithSnapLevel) {
  const unique_ptr<const S2Polygon> polygon(
      s2textformat::MakePolygon("0:0, 0:2, 2:0; 0:0, 0:-2, -2:-2, -2:0"));
  for (int level = 0; level <= S2CellId::kMaxLevel; ++level) {
    S2Polygon snapped_polygon;
    snapped_polygon.InitToSnapped(polygon.get(), level);
    EXPECT_TRUE(snapped_polygon.IsValid());
    double cell_angle = S2::kMaxDiag.GetValue(level);
    S1Angle merge_radius = S1Angle::Radians(cell_angle);
    EXPECT_TRUE(snapped_polygon.ApproxContains(polygon.get(), merge_radius));
  }
}

TEST(S2Polygon, MultipleInit) {
  unique_ptr<S2Polygon> polygon(s2textformat::MakePolygon("0:0, 0:2, 2:0"));
  EXPECT_EQ(1, polygon->num_loops());
  EXPECT_EQ(3, polygon->num_vertices());
  S2LatLngRect bound1 = polygon->GetRectBound();

  vector<S2Loop*> loops;
  loops.push_back(s2textformat::MakeLoop("10:0, -10:-20, -10:20"));
  loops.push_back(s2textformat::MakeLoop("40:30, 20:10, 20:50"));
  polygon->InitNested(&loops);
  EXPECT_TRUE(polygon->IsValid());
  EXPECT_TRUE(loops.empty());
  EXPECT_EQ(2, polygon->num_loops());
  EXPECT_EQ(6, polygon->num_vertices());
  EXPECT_TRUE(bound1 != polygon->GetRectBound());
}

TEST_F(S2PolygonTestBase, TestSimpleEncodeDecode) {
  Encoder encoder;
  cross1->Encode(&encoder);
  Decoder decoder(encoder.base(), encoder.length());
  S2Polygon decoded_polygon;
  ASSERT_TRUE(decoded_polygon.Decode(&decoder));
  EXPECT_TRUE(cross1->BoundaryEquals(&decoded_polygon));
  EXPECT_EQ(cross1->GetRectBound(), decoded_polygon.GetRectBound());
}

TEST_F(S2PolygonTestBase, TestEncodeDecodeDefaultPolygon) {
  S2Polygon polygon;
  EXPECT_TRUE(TestEncodeDecode(&polygon));
}

TEST_F(S2PolygonTestBase, CompressedEmptyPolygonRequires3Bytes) {
  S2Polygon empty_polygon;
  Encoder encoder;

  S2Polygon snapped_empty_polygon;
  snapped_empty_polygon.InitToSnapped(&empty_polygon);

  snapped_empty_polygon.Encode(&encoder);
  // 1 byte for version, 1 for the level, 1 for the length.
  EXPECT_EQ(1 + 1 + 1, encoder.length());

  EXPECT_TRUE(snapped_empty_polygon.is_empty());
  EXPECT_EQ(S2LatLngRect::Empty(), snapped_empty_polygon.GetRectBound());
}

TEST_F(S2PolygonTestBase, CompressedEncodedPolygonRequires69Bytes) {
  const unique_ptr<const S2Polygon> polygon(
      s2textformat::MakePolygon("0:0, 0:2, 2:0; 0:0, 0:-2, -2:-2, -2:0"));

  S2Polygon snapped_polygon;
  snapped_polygon.InitToSnapped(polygon.get());

  Encoder encoder;
  snapped_polygon.Encode(&encoder);

  // 2 loops, one with 3 vertices, one with 4.
  // Polygon:
  //   1 byte for version
  //   1 byte for level
  //   1 byte for num_loops
  // Loops:
  //   5 bytes overhead
  //   8 bytes per vertex
  EXPECT_EQ(1 + 1 + 1 + 2 * 5 + 7 * 8, encoder.length());
}

TEST_F(S2PolygonTestBase, CompressedEncodedPolygonDecodesApproxEqual) {
  // To compare the boundaries, etc we want to snap first.
  S2Polygon snapped;
  snapped.InitToSnapped(near_30);
  ASSERT_EQ(2, snapped.num_loops());
  EXPECT_EQ(0, snapped.loop(0)->depth());
  EXPECT_EQ(1, snapped.loop(1)->depth());

  Encoder encoder;
  snapped.Encode(&encoder);

  Decoder decoder(encoder.base(), encoder.length());

  S2Polygon decoded_polygon;
  ASSERT_TRUE(decoded_polygon.Decode(&decoder));
  ASSERT_TRUE(decoded_polygon.IsValid());
  EXPECT_TRUE(snapped.BoundaryEquals(&decoded_polygon));
  EXPECT_EQ(snapped.GetRectBound(), decoded_polygon.GetRectBound());
  EXPECT_EQ(snapped.num_vertices(), decoded_polygon.num_vertices());
  EXPECT_EQ(2, decoded_polygon.num_loops());
  EXPECT_EQ(0, decoded_polygon.loop(0)->depth());
  EXPECT_EQ(1, decoded_polygon.loop(1)->depth());
}

// This test checks that S2Polygons created directly from S2Cells behave
// identically to S2Polygons created from the vertices of those cells; this
// previously was not the case, because S2Cells calculate their bounding
// rectangles slightly differently, and S2Polygons created from them just
// copied the S2Cell bounds.
TEST(S2Polygon, TestS2CellConstructorAndContains) {
  S2LatLng latlng(S1Angle::E6(40565459), S1Angle::E6(-74645276));
  S2Cell cell(S2CellId::FromLatLng(latlng));
  S2Polygon cell_as_polygon(cell);
  S2Polygon empty;
  S2Polygon polygon_copy;
  polygon_copy.InitToUnion(&cell_as_polygon, &empty);
  EXPECT_TRUE(polygon_copy.Contains(&cell_as_polygon));
  EXPECT_TRUE(cell_as_polygon.Contains(&polygon_copy));
}

TEST(S2PolygonTest, Project) {
  unique_ptr<S2Polygon> polygon(MakePolygon(kNear0 + kNear2));
  S2Point point;
  S2Point projected;

  // The point inside the polygon should be projected into itself.
  point = s2textformat::MakePoint("1.1:0");
  projected = polygon->Project(point);
  EXPECT_TRUE(S2::ApproxEquals(point, projected));

  // The point is on the outside of the polygon.
  point = s2textformat::MakePoint("5.1:-2");
  projected = polygon->Project(point);
  EXPECT_TRUE(S2::ApproxEquals(s2textformat::MakePoint("5:-2"), projected));

  // The point is inside the hole in the polygon.
  point = s2textformat::MakePoint("-0.49:-0.49");
  projected = polygon->Project(point);
  EXPECT_TRUE(S2::ApproxEquals(s2textformat::MakePoint("-0.5:-0.5"),
                               projected, 1e-6));

  point = s2textformat::MakePoint("0:-3");
  projected = polygon->Project(point);
  EXPECT_TRUE(S2::ApproxEquals(s2textformat::MakePoint("0:-2"), projected));
}

// Helper function for testing the distance methods.  "boundary_x" is the
// expected result of projecting "x" onto the polygon boundary.  For
// convenience it can be set to S2Point() to indicate that (boundary_x == x).
static void TestDistanceMethods(S2Polygon const& polygon, S2Point const& x,
                                S2Point boundary_x) {
  // This error is not guaranteed by the implementation but is okay for tests.
  S1Angle const kMaxError = S1Angle::Radians(1e-15);

  if (boundary_x == S2Point()) boundary_x = x;
  EXPECT_LE(S1Angle(boundary_x, polygon.ProjectToBoundary(x)), kMaxError);

  if (polygon.is_empty() || polygon.is_full()) {
    EXPECT_EQ(S1Angle::Infinity(), polygon.GetDistanceToBoundary(x));
  } else {
    // EXPECT_NEAR only works with doubles.
    EXPECT_NEAR(S1Angle(x, boundary_x).degrees(),
                polygon.GetDistanceToBoundary(x).degrees(),
                kMaxError.degrees());
  }
  if (polygon.Contains(x)) {
    EXPECT_EQ(S1Angle::Zero(), polygon.GetDistance(x));
    EXPECT_EQ(x, polygon.Project(x));
  } else {
    EXPECT_EQ(polygon.GetDistanceToBoundary(x), polygon.GetDistance(x));
    EXPECT_EQ(polygon.ProjectToBoundary(x), polygon.Project(x));
  }
}

TEST_F(S2PolygonTestBase, GetDistance) {
  // The empty and full loops don't have boundaries.
  TestDistanceMethods(*empty, S2Point(0, 1, 0), S2Point());
  TestDistanceMethods(*full, S2Point(0, 1, 0), S2Point());

  // A polygon consisting of two nested rectangles centered around
  // S2LatLng(0,0).  Note that because lines of latitude are curved on the
  // sphere, it is not straightforward to project points onto any edge except
  // along the equator.  (The equator is the only line of latitude that is
  // also a geodesic.)
  unique_ptr<S2Polygon> nested(s2textformat::MakePolygon(
      "3:1, 3:-1, -3:-1, -3:1; 4:2, 4:-2, -4:-2, -4:2;"));

  // All points on the boundary of the polygon should be at distance zero.
  for (int i = 0; i < nested->num_loops(); i++) {
    S2Loop const* loop = nested->loop(i);
    for (int j = 0; j < loop->num_vertices(); j++) {
      // A vertex.
      TestDistanceMethods(*nested, loop->vertex(j), S2Point());
      // A point along an edge.
      TestDistanceMethods(*nested, S2EdgeUtil::Interpolate(
          S2Testing::rnd.RandDouble(), loop->vertex(j), loop->vertex(j+1)),
                          S2Point());
    }
  }
  // A point outside the outer shell that projects to an edge.
  TestDistanceMethods(*nested, S2LatLng::FromDegrees(0, -4.7).ToPoint(),
                      S2LatLng::FromDegrees(0, -2).ToPoint());
  // A point outside the outer shell that projects to a vertex.
  TestDistanceMethods(*nested, S2LatLng::FromDegrees(6, -3).ToPoint(),
                      S2LatLng::FromDegrees(4, -2).ToPoint());
  // A point inside the polygon that projects to an outer edge.
  TestDistanceMethods(*nested, S2LatLng::FromDegrees(0, 1.7).ToPoint(),
                      S2LatLng::FromDegrees(0, 2).ToPoint());
  // A point inside the polygon that projects to an inner vertex.
  TestDistanceMethods(*nested, S2LatLng::FromDegrees(-3.3, -1.3).ToPoint(),
                      S2LatLng::FromDegrees(-3, -1).ToPoint());
  // A point inside the inner hole.
  TestDistanceMethods(*nested, S2LatLng::FromDegrees(0, 0.1).ToPoint(),
                      S2LatLng::FromDegrees(0, 1).ToPoint());
}

class IsValidTest : public testing::Test {
 public:
  IsValidTest() {
    init_oriented_ = false;
    modify_polygon_hook_ = NULL;
    rnd_ = &S2Testing::rnd;
    rnd_->Reset(FLAGS_s2_random_seed);
  }

  ~IsValidTest() {
    Reset();
  }

  vector<S2Point>* AddLoop() {
    vector<S2Point>* vloop = new vector<S2Point>;
    vloops_.push_back(vloop);
    return vloop;
  }

  // Create "num_loops" nested regular loops around a common center point.
  // All loops have the same number of vertices (at least "min_vertices").
  // Furthermore, the vertices at the same index position are collinear with
  // the common center point of all the loops.  The loop radii decrease
  // exponentially in order to prevent accidental loop crossings when one of
  // the loops is modified.
  void AddConcentricLoops(int num_loops, int min_vertices) {
    DCHECK_LE(num_loops, 10);  // Because radii decrease exponentially.
    S2Point center = S2Testing::RandomPoint();
    int num_vertices = min_vertices + rnd_->Uniform(10);
    for (int i = 0; i < num_loops; ++i) {
      S1Angle radius = S1Angle::Degrees(80 * pow(0.1, i));
      *AddLoop() = S2Testing::MakeRegularPoints(center, radius, num_vertices);
    }
  }

  void Reset() {
    for (int i = 0; i < vloops_.size(); ++i) {
      delete vloops_[i];
    }
    vloops_.clear();
  }

  void CheckInvalid(string const& snippet) {
    vector<S2Loop*> loops;
    for (int i = 0; i < vloops_.size(); ++i) {
      loops.push_back(new S2Loop(*vloops_[i], DISABLE_S2DEBUG));
    }
    std::random_shuffle(loops.begin(), loops.end(), *rnd_);
    S2Polygon polygon;
    polygon.set_s2debug_override(DISABLE_S2DEBUG);
    if (init_oriented_) {
      polygon.InitOriented(&loops);
    } else {
      polygon.InitNested(&loops);
    }
    if (modify_polygon_hook_) (*modify_polygon_hook_)(&polygon);
    S2Error error;
    EXPECT_TRUE(polygon.FindValidationError(&error));
    EXPECT_TRUE(error.text().find(snippet) != string::npos)
        << "\nActual error: " << error << "\nExpected substring: " << snippet;
    Reset();
  }

 protected:
  static int const kIters = 100;

  bool init_oriented_;
  void (*modify_polygon_hook_)(S2Polygon*);
  S2Testing::Random* rnd_;
  vector<vector<S2Point>*> vloops_;
};

TEST_F(IsValidTest, UnitLength) {
  // This test can only be run in optimized builds because there are
  // DCHECK(IsUnitLength()) calls scattered throughout the S2 code.
  if (google::DEBUG_MODE) return;
  for (int iter = 0; iter < kIters; ++iter) {
    AddConcentricLoops(1 + rnd_->Uniform(6), 3 /*min_vertices*/);
    vector<S2Point>* vloop = vloops_[rnd_->Uniform(vloops_.size())];
    S2Point* p = &(*vloop)[rnd_->Uniform(vloop->size())];
    switch (rnd_->Uniform(3)) {
      case 0: *p = S2Point(0, 0, 0); break;
      case 1: *p *= 1e-30 * pow(1e60, rnd_->RandDouble()); break;
      case 2: *p = numeric_limits<double>::quiet_NaN() * S2Point(); break;
    }
    CheckInvalid("unit length");
  }
}

TEST_F(IsValidTest, VertexCount) {
  for (int iter = 0; iter < kIters; ++iter) {
    vector<S2Point>* vloop = AddLoop();
    if (rnd_->OneIn(2)) {
      vloop->push_back(S2Testing::RandomPoint());
      vloop->push_back(S2Testing::RandomPoint());
    }
    CheckInvalid("at least 3 vertices");
  }
}

TEST_F(IsValidTest, DuplicateVertex) {
  for (int iter = 0; iter < kIters; ++iter) {
    AddConcentricLoops(1, 3 /*min_vertices*/);
    vector<S2Point>* vloop = vloops_[0];
    int n = vloop->size();
    int i = rnd_->Uniform(n);
    int j = rnd_->Uniform(n - 1);
    (*vloop)[i] = (*vloop)[j + (j >= i)];
    CheckInvalid("duplicate vertex");
  }
}

TEST_F(IsValidTest, SelfIntersection) {
  for (int iter = 0; iter < kIters; ++iter) {
    // Use multiple loops so that we can test both holes and shells.  We need
    // at least 5 vertices so that the modified edges don't intersect any
    // nested loops.
    AddConcentricLoops(1 + rnd_->Uniform(6), 5 /*min_vertices*/);
    vector<S2Point>* vloop = vloops_[rnd_->Uniform(vloops_.size())];
    int n = vloop->size();
    int i = rnd_->Uniform(n);
    swap((*vloop)[i], (*vloop)[(i+1) % n]);
    CheckInvalid("crosses edge");
  }
}

TEST_F(IsValidTest, EmptyLoop) {
  for (int iter = 0; iter < kIters; ++iter) {
    AddConcentricLoops(rnd_->Uniform(5), 3 /*min_vertices*/);
    *AddLoop() = S2Loop::kEmpty();
    CheckInvalid("empty loop");
  }
}

TEST_F(IsValidTest, FullLoop) {
  for (int iter = 0; iter < kIters; ++iter) {
    // This is only an error if there is at least one other loop.
    AddConcentricLoops(1 + rnd_->Uniform(5), 3 /*min_vertices*/);
    *AddLoop() = S2Loop::kFull();
    CheckInvalid("full loop");
  }
}

TEST_F(IsValidTest, LoopsCrossing) {
  for (int iter = 0; iter < kIters; ++iter) {
    AddConcentricLoops(2, 4 /*min_vertices*/);
    // Both loops have the same number of vertices, and vertices at the same
    // index position are collinear with the center point, so we can create a
    // crossing by simply exchanging two vertices at the same index position.
    int n = vloops_[0]->size();
    int i = rnd_->Uniform(n);
    swap((*vloops_[0])[i], (*vloops_[1])[i]);
    if (rnd_->OneIn(2)) {
      // By copy the two adjacent vertices from one loop to the other, we can
      // ensure that the crossings happen at vertices rather than edges.
      (*vloops_[0])[(i+1) % n] = (*vloops_[1])[(i+1) % n];
      (*vloops_[0])[(i+n-1) % n] = (*vloops_[1])[(i+n-1) % n];
    }
    CheckInvalid("crosses loop");
  }
}

TEST_F(IsValidTest, DuplicateEdge) {
  for (int iter = 0; iter < kIters; ++iter) {
    AddConcentricLoops(2, 4 /*min_vertices*/);
    int n = vloops_[0]->size();
    if (rnd_->OneIn(2)) {
      // Create a shared edge (same direction in both loops).
      int i = rnd_->Uniform(n);
      (*vloops_[0])[i] = (*vloops_[1])[i];
      (*vloops_[0])[(i+1) % n] = (*vloops_[1])[(i+1) % n];
    } else {
      // Create a reversed edge (opposite direction in each loop) by cutting
      // loop 0 into two halves along one of its diagonals and replacing both
      // loops with the result.
      int split = 2 + rnd_->Uniform(n - 3);
      vloops_[1]->clear();
      vloops_[1]->push_back((*vloops_[0])[0]);
      for (int i = split; i < n; ++i) {
        vloops_[1]->push_back((*vloops_[0])[i]);
      }
      vloops_[0]->resize(split + 1);
    }
    CheckInvalid("has duplicate");
  }
}

TEST_F(IsValidTest, InconsistentOrientations) {
  for (int iter = 0; iter < kIters; ++iter) {
    AddConcentricLoops(2 + rnd_->Uniform(5), 3 /*min_vertices*/);
    init_oriented_ = true;
    CheckInvalid("Inconsistent loop orientations");
  }
}

static void SetInvalidLoopDepth(S2Polygon* polygon) {
  int i = S2Testing::rnd.Uniform(polygon->num_loops());
  if (i == 0 || S2Testing::rnd.OneIn(3)) {
    polygon->loop(i)->set_depth(-1);
  } else {
    polygon->loop(i)->set_depth(polygon->loop(i-1)->depth() + 2);
  }
}

TEST_F(IsValidTest, LoopDepthNegative) {
  modify_polygon_hook_ = SetInvalidLoopDepth;
  for (int iter = 0; iter < kIters; ++iter) {
    AddConcentricLoops(1 + rnd_->Uniform(4), 3 /*min_vertices*/);
    CheckInvalid("invalid loop depth");
  }
}

static void SetInvalidLoopNesting(S2Polygon* polygon) {
  int i = S2Testing::rnd.Uniform(polygon->num_loops());
  polygon->loop(i)->Invert();
}

TEST_F(IsValidTest, LoopNestingInvalid) {
  modify_polygon_hook_ = SetInvalidLoopNesting;
  for (int iter = 0; iter < kIters; ++iter) {
    AddConcentricLoops(2 + rnd_->Uniform(4), 3 /*min_vertices*/);
    CheckInvalid("Invalid nesting");
  }
}

TEST_F(IsValidTest, FuzzTest) {
  // Check that the S2Loop/S2Polygon constructors and IsValid() don't crash
  // when they receive arbitrary invalid input.  (We don't test large inputs;
  // it is assumed that the client enforces their own size limits before even
  // attempting to construct geometric objects.)
  if (google::DEBUG_MODE)
    return;  // Requires unit length vertices.
  for (int iter = 0; iter < kIters; ++iter) {
    int num_loops = 1 + rnd_->Uniform(10);
    for (int i = 0; i < num_loops; ++i) {
      int num_vertices = rnd_->Uniform(10);
      vector<S2Point>* vloop = AddLoop();
      while (vloop->size() < num_vertices) {
        // Since the number of vertices is random, we automatically test empty
        // loops, full loops, and invalid vertex counts.  Also since most
        // vertices are random, we automatically get self-intersections and
        // loop crossings.  That leaves zero and NaN vertices, duplicate
        // vertices, and duplicate edges to be created explicitly.
        if (rnd_->OneIn(10)) {
          // Zero vertex.
          vloop->push_back(S2Point(0, 0, 0));
        } else if (rnd_->OneIn(10)) {
          // NaN vertex.
          vloop->push_back(numeric_limits<double>::quiet_NaN() * S2Point());
        } else if (rnd_->OneIn(10) && !vloop->empty()) {
          // Duplicate vertex.
          vloop->push_back((*vloop)[rnd_->Uniform(vloop->size())]);
        } else if (rnd_->OneIn(10) && vloop->size() + 2 <= num_vertices) {
          // Try to copy an edge from a random loop.
          vector<S2Point>* other = vloops_[rnd_->Uniform(vloops_.size())];
          int n = other->size();
          if (n >= 2) {
            int k0 = rnd_->Uniform(n);
            int k1 = (k0 + 1) % n;
            if (rnd_->OneIn(2)) swap(k0, k1);  // Copy reversed edge.
            vloop->push_back((*other)[k0]);
            vloop->push_back((*other)[k1]);
          }
        } else {
          // Random non-unit-length point.
          S2Point p = S2Testing::RandomPoint();
          vloop->push_back(1e-30 * pow(1e60, rnd_->RandDouble()) * p);
        }
      }
    }
    CheckInvalid("");  // We could get any error message.
  }
}

// Returns the diameter of a loop (maximum distance between any two
// points in the loop).
S1Angle LoopDiameter(S2Loop const& loop) {
  S1Angle diameter;
  for (int i = 0; i < loop.num_vertices(); ++i) {
    S2Point test_point = loop.vertex(i);
    for (int j = i + 1; j < loop.num_vertices(); ++j) {
      diameter = max(diameter,
                     S2EdgeUtil::GetDistance(test_point, loop.vertex(j),
                                             loop.vertex(j+1)));
    }
  }
  return diameter;
}

// Returns the maximum distance from any vertex of poly_a to poly_b, that is,
// the directed Haussdorf distance of the set of vertices of poly_a to the
// boundary of poly_b.
//
// Doesn't consider loops from poly_a that have diameter less than min_diameter
// in degrees.
double MaximumDistanceInDegrees(S2Polygon const& poly_a,
                                S2Polygon const& poly_b,
                                double min_diameter_in_degrees) {
  double min_distance = 360;
  bool has_big_loops = false;
  for (int l = 0; l < poly_a.num_loops(); ++l) {
    S2Loop const* a_loop = poly_a.loop(l);
    if (LoopDiameter(*a_loop).degrees() <= min_diameter_in_degrees) {
      continue;
    }
    has_big_loops = true;
    for (int v = 0; v < a_loop->num_vertices(); ++v) {
      double distance = poly_b.GetDistance(a_loop->vertex(v)).degrees();
      if (distance < min_distance) {
        min_distance = distance;
      }
    }
  }
  if (has_big_loops) {
    return min_distance;
  } else {
    return 0.;  // As if the first polygon were empty.
  }
}

class S2PolygonSimplifierTest : public ::testing::Test {
 protected:
  S2PolygonSimplifierTest() {
    simplified = NULL;
    original = NULL;
  }

  ~S2PolygonSimplifierTest() {
    delete simplified;
    delete original;
  }

  // Owns poly.
  void SetInput(S2Polygon* poly, double tolerance_in_degrees) {
    delete original;
    delete simplified;
    original = poly;

    simplified = new S2Polygon();
    simplified->InitToSimplified(original,
                                 S1Angle::Degrees(tolerance_in_degrees),
                                 false);  // snap_to_cell_centers
  }

  void SetInput(string const& poly, double tolerance_in_degrees) {
    SetInput(s2textformat::MakePolygon(poly), tolerance_in_degrees);
  }

  S2Polygon* simplified;
  S2Polygon* original;
};

TEST_F(S2PolygonSimplifierTest, NoSimplification) {
  SetInput("0:0, 0:20, 20:20, 20:0", 1.0);
  EXPECT_EQ(4, simplified->num_vertices());

  EXPECT_EQ(0, MaximumDistanceInDegrees(*simplified, *original, 0));
  EXPECT_EQ(0, MaximumDistanceInDegrees(*original, *simplified, 0));
}

// Here, 10:-2 will be removed and  0:0-20:0 will intersect two edges.
// (The resulting polygon will in fact probably have more edges.)
TEST_F(S2PolygonSimplifierTest, SimplifiedLoopSelfIntersects) {
  SetInput("0:0, 0:20, 10:-0.1, 20:20, 20:0, 10:-0.2", 0.22);

  // To make sure that this test does something, we check
  // that the vertex 10:-0.2 is not in the simplification anymore.
  S2Point test_point = s2textformat::MakePoint("10:-0.2");
  EXPECT_LT(0.05, simplified->GetDistance(test_point).degrees());

  EXPECT_GE(0.22, MaximumDistanceInDegrees(*simplified, *original, 0));
  EXPECT_GE(0.22, MaximumDistanceInDegrees(*original, *simplified, 0.22));
}

TEST_F(S2PolygonSimplifierTest, NoSimplificationManyLoops) {
  SetInput("0:0,    0:1,   1:0;   0:20, 0:21, 1:20; "
           "20:20, 20:21, 21:20; 20:0, 20:1, 21:0", 0.01);
  EXPECT_EQ(0, MaximumDistanceInDegrees(*simplified, *original, 0));
  EXPECT_EQ(0, MaximumDistanceInDegrees(*original, *simplified, 0));
}

TEST_F(S2PolygonSimplifierTest, TinyLoopDisappears) {
  SetInput("0:0, 0:1, 1:1, 1:0", 1.1);
  EXPECT_TRUE(simplified->is_empty());
}

TEST_F(S2PolygonSimplifierTest, StraightLinesAreSimplified) {
  SetInput("0:0, 1:0, 2:0, 3:0, 4:0, 5:0, 6:0,"
           "6:1, 5:1, 4:1, 3:1, 2:1, 1:1, 0:1", 0.01);
  EXPECT_EQ(4, simplified->num_vertices());
}

TEST_F(S2PolygonSimplifierTest, EdgeSplitInManyPieces) {
  // near_square's right four-point side will be simplified to a vertical
  // line at lng=7.9, that will cut the 9 teeth of the saw (the edge will
  // therefore be broken into 19 pieces).
  const string saw =
      "1:1, 1:8, 2:2, 2:8, 3:2, 3:8, 4:2, 4:8, 5:2, 5:8,"
      "6:2, 6:8, 7:2, 7:8, 8:2, 8:8, 9:2, 9:8, 10:1";
  const string near_square =
      "0:0, 0:7.9, 1:8.1, 10:8.1, 11:7.9, 11:0";
  SetInput(saw + ";" + near_square, 0.21);

  EXPECT_TRUE(simplified->IsValid());
  EXPECT_GE(0.11, MaximumDistanceInDegrees(*simplified, *original, 0));
  EXPECT_GE(0.11, MaximumDistanceInDegrees(*original, *simplified, 0));
  // The resulting polygon's 9 little teeth are very small and disappear
  // due to the vertex_merge_radius of the polygon builder.  There remains
  // nine loops.
  EXPECT_EQ(9, simplified->num_loops());
}

TEST_F(S2PolygonSimplifierTest, EdgesOverlap) {
  // Two loops, One edge of the second one ([0:1 - 0:2]) is part of an
  // edge of the first one..
  SetInput("0:0, 0:3, 1:0; 0:1, -1:1, 0:2", 0.01);
  unique_ptr<S2Polygon> true_poly(
      s2textformat::MakePolygon("0:3, 1:0, 0:0, 0:1, -1:1, 0:2"));
  EXPECT_TRUE(simplified->BoundaryApproxEquals(true_poly.get()));
}

// Creates a loop from a comma separated list of u:v coordinates relative to a
// cell. The loop "0:0, 1:0, 1:1, 0:1" is counter-clockwise.
S2Polygon* MakeCellLoop(const S2Cell& cell, string const& str) {
  vector<S2Point> vertices;
  for (int i = 0; i < 4; ++i) {
    vertices.push_back(cell.GetVertex(i));
  }
  vector<pair<string, string> > p;
  CHECK(DictionaryParse(str, &p)) << ": str == \"" << str << "\"";
  vector<S2Point> loop_vertices;
  for (int i = 0; i < p.size(); ++i) {
    double u = strtod(p[i].first.c_str(), NULL);
    double v = strtod(p[i].second.c_str(), NULL);
    Vector3_d p = (vertices[0] * (1.0 - u) + vertices[1] * u) * (1.0 - v) +
        (vertices[2] * (1.0 - u) + vertices[3] * u) * v;
    loop_vertices.push_back(p.Normalize());
  }
  vector<S2Loop*> loops;
  loops.push_back(new S2Loop(loop_vertices));
  return new S2Polygon(&loops);
}

TEST(InitToSimplifiedInCell, PointsOnCellBoundaryKept) {
  S2CellId cell_id = S2CellId::FromToken("89c25c");
  S2Cell cell(cell_id);
  unique_ptr<S2Polygon> loop(MakeCellLoop(cell, "0.1:0, 0.2:0, 0.2:0.5"));
  S1Angle tolerance =
      S1Angle(loop->loop(0)->vertex(0), loop->loop(0)->vertex(1)) * 1.1;
  S2Polygon simplified_loop;
  simplified_loop.InitToSimplified(loop.get(), tolerance,
                                   false);  // snap_to_cell_centers
  EXPECT_TRUE(simplified_loop.is_empty());
  S2Polygon simplified_loop_in_cell;
  simplified_loop_in_cell.InitToSimplifiedInCell(loop.get(), cell, tolerance);
  EXPECT_TRUE(loop->BoundaryEquals(&simplified_loop_in_cell));
  EXPECT_EQ(3, simplified_loop_in_cell.num_vertices());
  EXPECT_EQ(-1, simplified_loop.GetSnapLevel());
}

TEST(InitToSimplifiedInCell, PointsInsideCellSimplified) {
  S2CellId cell_id = S2CellId::FromToken("89c25c");
  S2Cell cell(cell_id);
  unique_ptr<S2Polygon> loop(MakeCellLoop(cell,
                                          "0.1:0, 0.2:0, 0.2:0.5, 0.2:0.8"));
  S1Angle tolerance =
      S1Angle(loop->loop(0)->vertex(0), loop->loop(0)->vertex(1)) * 1.1;
  S2Polygon simplified_loop;
  simplified_loop.InitToSimplifiedInCell(loop.get(), cell, tolerance);
  EXPECT_TRUE(loop->BoundaryNear(&simplified_loop));
  EXPECT_EQ(3, simplified_loop.num_vertices());
  EXPECT_EQ(-1, simplified_loop.GetSnapLevel());
}

S2Polygon* MakeRegularPolygon(const string& center,
                              int num_points, double radius_in_degrees) {
  S1Angle radius = S1Angle::Degrees(radius_in_degrees);
  S2Loop* l = S2Testing::MakeRegularLoop(s2textformat::MakePoint(center),
                                         radius, num_points);
  vector<S2Loop*> loops;
  loops.push_back(l);
  return new S2Polygon(&loops);
}

// Tests that a regular polygon with many points gets simplified
// enough.
TEST_F(S2PolygonSimplifierTest, LargeRegularPolygon) {
  const double kRadius = 2.;  // in degrees
  const int num_initial_points = 1000;
  const int num_desired_points = 250;
  double tolerance = 1.05 * kRadius * (1 - cos(M_PI / num_desired_points));

  S2Polygon* p = MakeRegularPolygon("0:0", num_initial_points, kRadius);
  SetInput(p, tolerance);

  EXPECT_GE(tolerance, MaximumDistanceInDegrees(*simplified, *original, 0));
  EXPECT_GE(tolerance, MaximumDistanceInDegrees(*original, *simplified, 0));
  EXPECT_GE(250, simplified->num_vertices());
  EXPECT_LE(200, simplified->num_vertices());
}

class S2PolygonDecodeTest : public ::testing::Test {
 protected:
  S2PolygonDecodeTest() : data_array_(kMaxBytes) {
    encoder_.reset(data_array_.data(), kMaxBytes);
  }

  ~S2PolygonDecodeTest() override {}

  void AppendByte(int value) {
    encoder_.put8(value);
  }

  void AppendInt32(int value) {
    encoder_.put32(value);
  }

  void AppendRandomData(int size) {
    for (int i = 0; i < size && encoder_.avail() > 0; ++i) {
      AppendByte(random_.Uniform(256));
    }
  }

  void AppendRandomData() {
    AppendRandomData(random_.Uniform(kMaxBytes));
  }

  void AppendFakeLosslessEncodingData() {
    AppendByte(1);                      // polygon number
    AppendByte(0);                      // unused
    AppendByte(0);                      // "has holes" flag
    AppendInt32(PickRandomCount());     // num loops
    AppendByte(1);                      // loop version
    AppendInt32(PickRandomCount());     // num vertices
    AppendRandomData();                 // junk to fill out the buffer
  }

  void AppendFakeCompressedEncodingData() {
    AppendByte(4);                      // polygon number
    AppendByte(random_.Uniform(50));    // snap level
    AppendInt32(PickRandomCount());     // num loops
    AppendInt32(PickRandomCount());     // num vertices
    AppendRandomData();                 // junk to fill out the buffer
  }

  int32 PickRandomCount() {
    if (random_.OneIn(10)) {
      return -1;
    }
    if (random_.OneIn(10)) {
      return 0;
    }
    if (random_.OneIn(10)) {
      return 1000000000;
    }
    if (random_.OneIn(2)) {
      return random_.Uniform(1000000000);
    }
    return random_.Uniform(1000);
  }

  bool Test() {
    decoder_.reset(data_array_.data(), encoder_.length());
    encoder_.clear();
    S2Polygon polygon;
    polygon.set_s2debug_override(DISABLE_S2DEBUG);
    return polygon.Decode(&decoder_);
  }

  // Random number generator.
  S2Testing::Random random_;

  // Maximum size of the data array.
  const int kMaxBytes = 256;

  // The data array.
  FixedArray<int8> data_array_;

  // Encoder that is used to put data into the array.
  Encoder encoder_;

  // Decoder used to extract data from the array.
  Decoder decoder_;
};

TEST_F(S2PolygonDecodeTest, FuzzLosslessEncoding) {
  // Some parts of the S2 library DCHECK on invalid data, even if we set
  // FLAGS_s2debug to false or use S2Polygon::set_s2debug_override. So we
  // only run this test in opt mode.
#ifdef NDEBUG
  for (int i = 0; i < 100000; ++i) {
    AppendFakeLosslessEncodingData();
    Test();
  }
#endif
}

TEST_F(S2PolygonDecodeTest, FuzzCompressedEncoding) {
  // Some parts of the S2 library DCHECK on invalid data, even if we set
  // FLAGS_s2debug to false or use S2Polygon::set_s2debug_override. So we
  // only run this test in opt mode.
#ifdef NDEBUG
  for (int i = 0; i < 100000; ++i) {
    AppendFakeCompressedEncodingData();
    Test();
  }
#endif
}

TEST_F(S2PolygonDecodeTest, FuzzEverything) {
  // Some parts of the S2 library DCHECK on invalid data, even if we set
  // FLAGS_s2debug to false or use S2Polygon::set_s2debug_override. So we
  // only run this test in opt mode.
#ifdef NDEBUG
  for (int i = 0; i < 100000; ++i) {
    AppendRandomData();
    Test();
  }
#endif
}

TEST_F(S2PolygonTestBase, FullPolygonShape) {
  S2Polygon::Shape shape(full);
  EXPECT_EQ(0, shape.num_edges());
  EXPECT_TRUE(shape.has_interior());
  EXPECT_TRUE(shape.contains_origin());
}

TEST_F(S2PolygonTestBase, EmptyPolygonShape) {
  S2Polygon::Shape shape(empty);
  EXPECT_EQ(0, shape.num_edges());
  EXPECT_TRUE(shape.has_interior());
  EXPECT_FALSE(shape.contains_origin());
}

void TestPolygonShape(S2Polygon const* polygon) {
  DCHECK(!polygon->is_full());
  // Allocate the shape so that we can Release() it below.
  S2Polygon::Shape* shape = new S2Polygon::Shape(polygon);
  EXPECT_EQ(polygon, shape->polygon());
  EXPECT_EQ(polygon->num_vertices(), shape->num_edges());
  for (int e = 0, i = 0; i < polygon->num_loops(); ++i) {
    S2Loop const* loop_i = polygon->loop(i);
    for (int j = 0; j < loop_i->num_vertices(); ++j, ++e) {
      S2Point const *v0, *v1;
      shape->GetEdge(e, &v0, &v1);
      EXPECT_EQ(&loop_i->vertex(j), v0);
      EXPECT_EQ(&loop_i->vertex(j+1), v1);
    }
  }
  EXPECT_TRUE(shape->has_interior());
  EXPECT_EQ(polygon->Contains(S2::Origin()), shape->contains_origin());
  // In debug mode this tests that "shape" is deleted but "polygon" is not.
  shape->Release();
}

TEST_F(S2PolygonTestBase, OneLoopPolygonShape) {
  TestPolygonShape(near_0);
}

TEST_F(S2PolygonTestBase, SeveralLoopPolygonShape) {
  TestPolygonShape(near_3210);
}

TEST_F(S2PolygonTestBase, ManyLoopPolygonShape) {
  int const kNumLoops = 1000;
  int const kNumVerticesPerLoop = 6;
  S2Polygon polygon;
  S2Testing::ConcentricLoopsPolygon(S2Point(1, 0, 0), kNumLoops,
                                    kNumVerticesPerLoop, &polygon);
  TestPolygonShape(&polygon);
}

