// Copyright 2018 Google Inc. All Rights Reserved.
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

#include "s2/s2loop_measures.h"

#include <cfloat>
#include <cstdint>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include "absl/log/absl_check.h"
#include "absl/log/absl_log.h"
#include "absl/log/absl_vlog_is_on.h"
#include "absl/log/log_streamer.h"
#include "absl/random/random.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "s2/s1angle.h"
#include "s2/s2debug.h"
#include "s2/s2latlng.h"
#include "s2/s2loop.h"
#include "s2/s2measures.h"
#include "s2/s2point.h"
#include "s2/s2point_span.h"
#include "s2/s2random.h"
#include "s2/s2testing.h"
#include "s2/s2text_format.h"
#include "s2/util/math/mathutil.h"

using absl::string_view;
using s2textformat::ParsePointsOrDie;
using std::fabs;
using std::min;
using std::string;
using std::vector;

namespace {

// Simple inefficient reference implementation of PruneDegeneracies() on a
// character string.  The result will be some cyclic permutation of what
// PruneDegeneracies() produces.
std::string BruteForceQuadraticPrune(const string_view s) {
  std::string answer(s);
  // Make repeated passes over answer, reducing AAs and ABAs to As
  // (and, incidentally, if the whole string has length 1 or 2, turn it
  // into the empty string).
  // Return answer as soon as we see that we've reached a fixed point.
  while (true) {
    bool something_changed = false;
    for (int i = 0; i < answer.size(); ++i) {
      if (answer[i] == answer[(i + 1) % answer.size()]) {
        // AA starting at i (or A is the whole string).
        // Remove answer[i].
        answer.erase(i, /*count=*/1);
        something_changed = true;
        break;
      } else if (answer[i] == answer[(i + 2) % answer.size()]) {
        // ABA starting at i (or AB is the whole string).
        if (i + 1 < answer.size()) {
          // Remove the AB at answer[i] and answer[i+1].
          answer.erase(i, /*count=*/2);
        } else {
          // Remove the BA at answer[0] and answer[1].
          answer.erase(0, /*count=*/2);
        }
        something_changed = true;
        break;
      }
    }
    if (!something_changed) return answer;
  }
}

// Return the lexicographically least cyclic permutation of s.
std::string BruteForceQuadraticCyclicallyCanonicalize(const string_view s) {
  std::string answer;
  for (int i = 0; i < s.size(); ++i) {
    const std::string candidate = absl::StrCat(s.substr(i), s.substr(0, i));
    if (i == 0 || candidate < answer) answer = candidate;
  }
  return answer;
}

// Given a string where each character "ch" represents a vertex (such as
// "abac"), returns a vector of S2Points of the form (ch, 0, 0).  Note that
// these points are not unit length and therefore are not suitable for general
// use; however, they are useful for testing certain functions below.
vector<S2Point> MakeTestLoop(string_view loop_str) {
  vector<S2Point> loop;
  for (const char ch : loop_str) {
    loop.push_back(S2Point(ch, 0, 0));
  }
  return loop;
}

// Given a loop whose vertices are represented as characters (such as "abcd" or
// "abccb"), verify that S2::PruneDegeneracies() yields a loop cyclically
// equivalent to "expected_str".
std::string TestPruneDegeneracies(string_view input_str,
                                  string_view expected_str) {
  const vector<S2Point> input = MakeTestLoop(input_str);
  vector<S2Point> new_vertices;
  string actual_str;
  for (const S2Point& p : S2::PruneDegeneracies(input, &new_vertices)) {
    actual_str.push_back(static_cast<char>(p[0]));
  }
  EXPECT_EQ(BruteForceQuadraticCyclicallyCanonicalize(actual_str),
            BruteForceQuadraticCyclicallyCanonicalize(expected_str))
      << "input_str = \"" << input_str << "\"";
  return actual_str;  // in case caller wants to show it
}

TEST(PruneDegeneracies, CompletelyDegenerate) {
  TestPruneDegeneracies("", "");
  TestPruneDegeneracies("a", "");
  TestPruneDegeneracies("aaaaa", "");
  TestPruneDegeneracies("ab", "");
  TestPruneDegeneracies("abb", "");
  TestPruneDegeneracies("aab", "");
  TestPruneDegeneracies("aba", "");
  TestPruneDegeneracies("abba", "");
  TestPruneDegeneracies("abcb", "");
  TestPruneDegeneracies("abcba", "");
  TestPruneDegeneracies("abcdcdedefedcbcdcb", "");
}

TEST(PruneDegeneracies, PartiallyDegenerate) {
  TestPruneDegeneracies("abc", "abc");
  TestPruneDegeneracies("abca", "abc");
  TestPruneDegeneracies("abcc", "abc");
  TestPruneDegeneracies("abccaa", "abc");
  TestPruneDegeneracies("aabbcc", "abc");
  TestPruneDegeneracies("abcdedca", "abc");
  TestPruneDegeneracies("abcbabcbcdc", "abc");
  TestPruneDegeneracies("xyzabcazy", "abc");
  TestPruneDegeneracies("xxyyzzaabbccaazzyyxx", "abc");
  TestPruneDegeneracies("abcdb", "bcd");
  TestPruneDegeneracies("abcdecb", "cde");
  TestPruneDegeneracies("abcdefdcb", "def");
  TestPruneDegeneracies("abcad", "bca");
  TestPruneDegeneracies("abcdbae", "cdb");
  TestPruneDegeneracies("abcdecbaf", "dec");
}

TEST(PruneDegeneracies, AllSmallCases) {
  for (int base = 0; base <= 10; ++base) {
    for (int exponent = 0; exponent <= 12; ++exponent) {
      // Do all base^exponent strings of length `exponent` using a `base`-letter
      // alphabet (as long as there are 5000 or fewer of them).

      const int64_t num_strings = MathUtil::IPow<int64_t>(base, exponent);
      if (num_strings > 5000) break;   // getting too many, for this base
      if (num_strings == 0) continue;  // base=0 and exponent>0 is not useful
      if (base > exponent) continue;  // more chars than positions is not useful

      ABSL_LOG(INFO) << "      pruning " << base << "^" << exponent << "="
                     << num_strings << " string"
                     << (num_strings == 1 ? "" : "s") << " of length "
                     << exponent << " from " << base << " character"
                     << (base == 1 ? "" : "s");
      for (int64_t i_string = 0; i_string < num_strings; ++i_string) {
        // Construct the i_string'th string to be the `exponent`-digit numeral
        // representing the integer i_string in base `base`, little-endian.
        std::string string;
        {
          int64_t scratch = i_string;
          for (int position = 0; position < exponent; ++position) {
            string.push_back(static_cast<char>('a' + (scratch % base)));
            scratch /= base;
          }
          ABSL_CHECK_EQ(scratch, 0);
        }
        const std::string result =
            TestPruneDegeneracies(string, BruteForceQuadraticPrune(string));

        // If log level >=3, show all strings of this alphabet and length;
        // if log level ==2, then just show the first 12 and last 2;
        // if log level <=1, don't show them at all.
        const int max_to_print_at_beginning =
            ABSL_VLOG_IS_ON(3) ? num_strings : 12;
        if (i_string < max_to_print_at_beginning ||
            (num_strings - i_string - 1) < 2) {
          // Output at log level >=2 looks like, e.g.:
          //     ...
          //     pruning 5^5=3125 strings
          //         ...
          //         10/3125: "acaaa" -> ""
          //         11/3125: "bcaaa" -> "bca"
          //         ...    (<-- actual ellipsis emitted here if log level is 2)
          //     ...
          ABSL_VLOG(2) << "          " << i_string << "/" << num_strings
                       << ": \"" << string << "\" -> \"" << result << "\"";
        } else if (i_string == max_to_print_at_beginning) {
          ABSL_VLOG(2) << "          ...";
        }
      }
    }
  }
}

// Given a loop whose vertices are represented as characters (such as "abcd" or
// "abccb"), verify that S2::GetCanonicalLoopOrder returns the given result.
void TestCanonicalLoopOrder(string_view input_str,
                            S2::LoopOrder expected_order) {
  EXPECT_EQ(expected_order, S2::GetCanonicalLoopOrder(MakeTestLoop(input_str)));
}

TEST(GetCanonicalLoopOrder, AllDegeneracies) {
  TestCanonicalLoopOrder("", S2::LoopOrder(0, 1));
  TestCanonicalLoopOrder("a", S2::LoopOrder(0, 1));
  TestCanonicalLoopOrder("aaaaa", S2::LoopOrder(0, 1));
  TestCanonicalLoopOrder("ba", S2::LoopOrder(1, 1));
  TestCanonicalLoopOrder("bab", S2::LoopOrder(1, 1));
  TestCanonicalLoopOrder("cbab", S2::LoopOrder(2, 1));
  TestCanonicalLoopOrder("bacbcab", S2::LoopOrder(8, -1));
}

TEST(GetPerimeter, Empty) {
  EXPECT_EQ(S1Angle::Zero(), S2::GetPerimeter(vector<S2Point>{}));
}

TEST(GetPerimeter, Octant) {
  auto loop = ParsePointsOrDie("0:0, 0:90, 90:0");
  EXPECT_DOUBLE_EQ(3 * M_PI_2, S2::GetPerimeter(loop).radians());
}

TEST(GetPerimeter, MoreThanTwoPi) {
  // Make sure that GetPerimeter doesn't use S1ChordAngle, which can only
  // represent distances up to 2*Pi.
  auto loop = ParsePointsOrDie("0:0, 0:90, 0:180, 90:0, 0:-90");
  EXPECT_DOUBLE_EQ(5 * M_PI_2, S2::GetPerimeter(loop).radians());
}

TEST(GetSignedArea, Underflow) {
  auto loop = ParsePointsOrDie("0:0, 0:1e-88, 1e-88:1e-88, 1e-88:0");
  EXPECT_GT(S2::GetSignedArea(loop), 0);
}

TEST(GetSignedArea, ErrorAccumulation) {
#if defined(__APPLE__) && defined(__aarch64__)
  GTEST_SKIP() << "https://github.com/google/s2geometry/issues/395";
#endif
  // Loop encompassing half an octant of the sphere.
  vector<S2Point> loop{
      {1.0, 0.0, 0.0},
      {M_SQRT1_2, M_SQRT1_2, 0.0},
      {0.0, 0.0, 1.0},
  };

  // Area of just one loop.
  const double expected_area = S2::GetSignedArea(loop);

  // Repeat the loop kIters times.  We shouldn't accumulate significant error.
  constexpr int kIters = 100001;
  ABSL_CHECK_EQ(kIters % 16, 1);  // Area wraps every sixteen loops.

  loop.reserve(3 * kIters);
  for (int i = 0; i < kIters - 1; ++i) {
    loop.push_back(loop[0]);
    loop.push_back(loop[1]);
    loop.push_back(loop[2]);
  }
  double actual_area = S2::GetSignedArea(loop);

  // For well conditioned sums, Kahan has ~constant error of 2 epsilon.  To get
  // absolute error we multiply by the sum of absolute values of the terms, so
  // kIters*expected_area.
  //
  // see: https://en.wikipedia.org/wiki/Kahan_summation_algorithm#Accuracy
  double allowed_error = 2 * DBL_EPSILON * (kIters * std::fabs(expected_area));
  EXPECT_NEAR(expected_area, actual_area, allowed_error);
}

class LoopTestBase : public testing::Test {
 protected:
  // Some standard loops to use in the tests (see descriptions below).
  vector<S2Point> full_;
  vector<S2Point> v_loop_;
  vector<S2Point> north_hemi_;
  vector<S2Point> north_hemi3_;
  vector<S2Point> west_hemi_;
  vector<S2Point> east_hemi_;
  vector<S2Point> candy_cane_;
  vector<S2Point> line_triangle_;
  vector<S2Point> skinny_chevron_;
  vector<S2Point> three_leaf_clover_;
  vector<S2Point> tessellated_loop_;

 public:
  LoopTestBase() :
      // The full loop is represented as a loop with no vertices.
      full_(),

      // A degenerate loop in the shape of a "V".
      v_loop_(ParsePointsOrDie("5:1, 0:2, 5:3, 0:2")),

      // The northern hemisphere, defined using two pairs of antipodal points.
      north_hemi_(ParsePointsOrDie("0:-180, 0:-90, 0:0, 0:90")),

      // The northern hemisphere, defined using three points 120 degrees apart.
      north_hemi3_(ParsePointsOrDie("0:-180, 0:-60, 0:60")),

      // The western hemisphere, defined using two pairs of antipodal points.
      west_hemi_(ParsePointsOrDie("0:-180, -90:0, 0:0, 90:0")),

      // The eastern hemisphere, defined using two pairs of antipodal points.
      east_hemi_(ParsePointsOrDie("90:0, 0:0, -90:0, 0:-180")),

      // A spiral stripe that slightly over-wraps the equator.
      candy_cane_(ParsePointsOrDie(
          "-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70")),

      // A completely degenerate triangle along the equator that Sign()
      // considers to be CCW.
      line_triangle_(ParsePointsOrDie("0:1, 0:2, 0:3")),

      // A nearly-degenerate CCW chevron near the equator with very long sides
      // (about 80 degrees).  Its area is less than 1e-640, which is too small
      // to represent in double precision.
      skinny_chevron_(ParsePointsOrDie("0:0, -1e-320:80, 0:1e-320, 1e-320:80")),

      // A loop where the same vertex appears three times.
      three_leaf_clover_(ParsePointsOrDie(
          "0:0, -3:3, 3:3, 0:0, 3:0, 3:-3, 0:0, -3:-3, -3:0")),

      // A loop with groups of 3 or more vertices in a straight line.
      tessellated_loop_(ParsePointsOrDie(
          "10:34, 5:34, 0:34, -10:34, -10:36, -5:36, 0:36, 10:36")) {
  }
};

static void TestAreaConsistentWithCurvature(const vector<S2Point>& loop) {
  // Check that the area computed using GetArea() is consistent with the loop
  // curvature.  According to the Gauss-Bonnet theorem, the area of the loop
  // equals 2*Pi minus its curvature.
  double area = S2::GetArea(loop);
  double gauss_area = 2 * M_PI - S2::GetCurvature(loop);
  // The error bound below is sufficient for current tests but not guaranteed.
  EXPECT_LE(fabs(area - gauss_area), 1e-14)
      << "Failed loop: " << s2textformat::ToString(loop)
      << "\nArea = " << area << ", Gauss Area = " << gauss_area;
}

TEST_F(LoopTestBase, GetAreaConsistentWithCurvature) {
  TestAreaConsistentWithCurvature(full_);
  TestAreaConsistentWithCurvature(north_hemi_);
  TestAreaConsistentWithCurvature(north_hemi3_);
  TestAreaConsistentWithCurvature(west_hemi_);
  TestAreaConsistentWithCurvature(east_hemi_);
  TestAreaConsistentWithCurvature(candy_cane_);
  TestAreaConsistentWithCurvature(line_triangle_);
  TestAreaConsistentWithCurvature(skinny_chevron_);
  TestAreaConsistentWithCurvature(three_leaf_clover_);
  TestAreaConsistentWithCurvature(tessellated_loop_);
}

TEST_F(LoopTestBase, GetSurfaceIntegralGreaterThan4Pi) {
  // This test demonstrates that even when GetSurfaceIntegral() returns a an
  // area greater than 4*Pi, GetSignedArea() still returns an accurate result.

  // GetSurfaceIntegral() returns an area > 4 * Pi for this loop.  (Note that
  // the result of GetSurfaceIntegral is only correct modulo 4 * Pi, and that
  // S2::GetSignedArea() automatically corrects for this.)
  const vector<S2Point> loop1 = {
    {1, 0, 0}, {0, 1, 1e-150}, S2Point{-1, -2, 0}.Normalize(),
    {-1, 0, 1e-50}, {0, 0, 1}
  };
  ASSERT_TRUE(S2Loop(loop1).IsValid());
  EXPECT_GT(S2::GetSurfaceIntegral(loop1, S2::SignedArea), 4 * M_PI + 0.1);
  TestAreaConsistentWithCurvature(loop1);
}

TEST_F(LoopTestBase, GetAreaConsistentWithOrientation) {
  // Test that GetArea() returns an area near 0 for degenerate loops that
  // contain almost no points, and an area near 4*Pi for degenerate loops that
  // contain almost all points.

  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "GET_AREA_CONSISTENT_WITH_ORIENTATION",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  static constexpr int kMaxVertices = 6;
  for (int i = 0; i < 50; ++i) {
    int num_vertices = absl::Uniform(bitgen, 3, kMaxVertices + 1);
    // Repeatedly choose N vertices that are exactly on the equator until we
    // find some that form a valid loop.
    vector<S2Point> loop;
    do {
      loop.clear();
      for (int i = 0; i < num_vertices; ++i) {
        // We limit longitude to the range [0, 90] to ensure that the loop is
        // degenerate (as opposed to following the entire equator).
        loop.push_back(
            S2LatLng::FromRadians(0, absl::Uniform(bitgen, 0.0, M_PI_2))
                .ToPoint());
      }
    } while (!S2Loop(loop, S2Debug::DISABLE).IsValid());
    bool ccw = S2::IsNormalized(loop);
    // The error bound is sufficient for current tests but not guaranteed.
    EXPECT_NEAR(ccw ? 0 : 4 * M_PI, S2::GetArea(loop), 1e-14)
        << "Failed loop " << i << ": " << s2textformat::ToString(loop);
    EXPECT_EQ(!ccw, S2Loop(loop).Contains(S2Point(0, 0, 1)));
  }
}

TEST_F(LoopTestBase, GetAreaAccuracy) {
  // TODO(ericv): Test that GetArea() has an accuracy significantly better
  // than 1e-15 on loops whose area is small.
}

TEST_F(LoopTestBase, GetAreaAndCentroid) {
  EXPECT_EQ(4 * M_PI, S2::GetArea(full_));
  EXPECT_EQ(S2Point(0, 0, 0), S2::GetCentroid(full_));

  EXPECT_DOUBLE_EQ(S2::GetArea(north_hemi_), 2 * M_PI);
  EXPECT_NEAR(2 * M_PI, S2::GetArea(east_hemi_), 1e-12);

  // Construct spherical caps of random height, and approximate their boundary
  // with closely spaces vertices.  Then check that the area and centroid are
  // correct.
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "GET_AREA_AND_CENTROID",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int iter = 0; iter < 50; ++iter) {
    // Choose a coordinate frame for the spherical cap.
    S2Point x, y, z;
    s2random::Frame(bitgen, x, y, z);

    // Given two points at latitude phi and whose longitudes differ by dtheta,
    // the geodesic between the two points has a maximum latitude of
    // atan(tan(phi) / cos(dtheta/2)).  This can be derived by positioning
    // the two points at (-dtheta/2, phi) and (dtheta/2, phi).
    //
    // We want to position the vertices close enough together so that their
    // maximum distance from the boundary of the spherical cap is kMaxDist.
    // Thus we want fabs(atan(tan(phi) / cos(dtheta/2)) - phi) <= kMaxDist.
    static const double kMaxDist = 1e-6;
    double height = absl::Uniform(bitgen, 0.0, 2.0);
    double phi = asin(1 - height);
    double max_dtheta = 2 * acos(tan(fabs(phi)) / tan(fabs(phi) + kMaxDist));
    max_dtheta = min(M_PI, max_dtheta);  // At least 3 vertices.

    vector<S2Point> loop;
    for (double theta = 0; theta < 2 * M_PI;
         theta += absl::Uniform(bitgen, 0.0, max_dtheta)) {
      loop.push_back(cos(theta) * cos(phi) * x +
                     sin(theta) * cos(phi) * y +
                     sin(phi) * z);
    }
    double area = S2::GetArea(loop);
    S2Point centroid = S2::GetCentroid(loop);
    double expected_area = 2 * M_PI * height;
    EXPECT_LE(fabs(area - expected_area), 2 * M_PI * kMaxDist);
    S2Point expected_centroid = expected_area * (1 - 0.5 * height) * z;
    EXPECT_LE((centroid - expected_centroid).Norm(), 2 * kMaxDist);
  }
}

static void ExpectSameOrder(S2PointLoopSpan loop1, S2::LoopOrder order1,
                            S2PointLoopSpan loop2, S2::LoopOrder order2) {
  ABSL_DCHECK_EQ(loop1.size(), loop2.size());
  int i1 = order1.first, i2 = order2.first;
  int dir1 = order1.dir, dir2 = order2.dir;
  for (int n = loop1.size(); --n >= 0; ) {
    ASSERT_EQ(loop1[i1], loop2[i2]) << ": " << order1 << " vs. " << order2;
    i1 += dir1;
    i2 += dir2;
  }
}

// Check that the curvature is *identical* when the vertex order is
// rotated, and that the sign is inverted when the vertices are reversed.
static void CheckCurvatureInvariants(const vector<S2Point>& loop_in) {
  S2::LoopOrder order_in = S2::GetCanonicalLoopOrder(loop_in);
  auto loop = loop_in;
  double expected = S2::GetCurvature(loop);
  for (int i = 0; i < loop.size(); ++i) {
    std::reverse(loop.begin(), loop.end());
    EXPECT_EQ((expected == 2 * M_PI) ? expected : -expected,
              S2::GetCurvature(loop));
    ExpectSameOrder(loop_in, order_in, loop, S2::GetCanonicalLoopOrder(loop));
    std::reverse(loop.begin(), loop.end());
    std::rotate(loop.begin(), loop.begin() + 1, loop.end());
    EXPECT_EQ(expected, S2::GetCurvature(loop));
    ExpectSameOrder(loop_in, order_in, loop, S2::GetCanonicalLoopOrder(loop));
  }
}

TEST_F(LoopTestBase, GetCurvature) {
  EXPECT_EQ(-2 * M_PI, S2::GetCurvature(full_));

  EXPECT_EQ(2 * M_PI, S2::GetCurvature(v_loop_));
  CheckCurvatureInvariants(v_loop_);

  // This curvature should be computed exactly.
  EXPECT_EQ(0, S2::GetCurvature(north_hemi3_));
  CheckCurvatureInvariants(north_hemi3_);

  EXPECT_NEAR(0, S2::GetCurvature(west_hemi_), 1e-15);
  CheckCurvatureInvariants(west_hemi_);

  // We don't have an easy way to estimate the curvature of these loops, but
  // we can still check that the expected invariants hold.
  CheckCurvatureInvariants(candy_cane_);
  CheckCurvatureInvariants(three_leaf_clover_);

  EXPECT_DOUBLE_EQ(2 * M_PI, S2::GetCurvature(line_triangle_));
  CheckCurvatureInvariants(line_triangle_);

  EXPECT_DOUBLE_EQ(2 * M_PI, S2::GetCurvature(skinny_chevron_));
  CheckCurvatureInvariants(skinny_chevron_);

  // Build a narrow spiral loop starting at the north pole.  This is designed
  // to test that the error in GetCurvature is linear in the number of
  // vertices even when the partial sum of the curvatures gets very large.
  // The spiral consists of two "arms" defining opposite sides of the loop.
  // This is a pathological loop that contains many long parallel edges.
  constexpr int kArmPoints = 10000;  // Number of vertices in each "arm"
  const double kArmRadius = 0.01;  // Radius of spiral.
  vector<S2Point> spiral(2 * kArmPoints);
  spiral[kArmPoints] = S2Point(0, 0, 1);
  for (int i = 0; i < kArmPoints; ++i) {
    double angle = (2 * M_PI / 3) * i;
    double x = cos(angle);
    double y = sin(angle);
    double r1 = i * kArmRadius / kArmPoints;
    double r2 = (i + 1.5) * kArmRadius / kArmPoints;
    spiral[kArmPoints - i - 1] = S2Point(r1 * x, r1 * y, 1).Normalize();
    spiral[kArmPoints + i] = S2Point(r2 * x, r2 * y, 1).Normalize();
  }

  // Check that GetCurvature() is consistent with GetArea() to within the
  // error bound of the former.  We actually use a tiny fraction of the
  // worst-case error bound, since the worst case only happens when all the
  // roundoff errors happen in the same direction and this test is not
  // designed to achieve that.  The error in GetArea() can be ignored for the
  // purposes of this test since it is generally much smaller.
  EXPECT_NEAR(2 * M_PI - S2::GetArea(spiral), S2::GetCurvature(spiral),
              0.01 * S2::GetCurvatureMaxError(spiral));
}

TEST(KahanSum, DefaultValue) {
  S2::internal::KahanSum<double> sum;
  EXPECT_EQ(0.0, (double)sum);
  EXPECT_EQ(0.0, sum.Compensation());
}

TEST(KahanSum, SingleValue) {
  S2::internal::KahanSum<double> sum;
  sum += -3;
  EXPECT_EQ(-3.0, (double)sum);
  EXPECT_EQ(0.0, sum.Compensation());
}

// Summing the squares has particularly bad behavior when doing the natural sum,
// and it has a closed formula for verifying the correctness.
TEST(KahanSum, SumOfSquares) {
  for (int direction = 0; direction < 2; ++direction) {
    S2::internal::KahanSum<double> safe_sum;
    double sum = 0;
    const int64_t n = 1000000;
    for (int64_t i = 0; i <= n; ++i) {
      const int64_t v = direction == 0 ? i : n - i;
      safe_sum += (v * v);
      sum += v * v;
    }
    const int64_t expected_sum = (2 * n + 1) * n * (n + 1) / 6;
    // Yes we *do* want a *strict* equality check here.
    EXPECT_EQ(static_cast<double>(expected_sum), (double)safe_sum);
    // To show that the trivial summation really gets it wrong!
    EXPECT_GE(fabs(sum - expected_sum), expected_sum / (n * n));
  }
}

}  // namespace
