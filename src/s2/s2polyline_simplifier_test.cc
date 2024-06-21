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

#include "s2/s2polyline_simplifier.h"

#include <cfloat>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include "absl/log/log_streamer.h"
#include "absl/random/random.h"
#include "absl/strings/string_view.h"
#include "s2/s1angle.h"
#include "s2/s1chord_angle.h"
#include "s2/s2edge_crossings.h"
#include "s2/s2edge_distances.h"
#include "s2/s2point.h"
#include "s2/s2random.h"
#include "s2/s2testing.h"
#include "s2/s2text_format.h"

using absl::string_view;
using std::vector;

void CheckSimplify(string_view src, string_view dst, string_view target,
                   string_view avoid, const vector<bool>& disc_on_left,
                   double radius_degrees, bool expected_result) {
  S1ChordAngle radius(S1Angle::Degrees(radius_degrees));
  S2PolylineSimplifier s;
  s.Init(s2textformat::MakePointOrDie(src));
  for (const S2Point& p : s2textformat::ParsePointsOrDie(target)) {
    s.TargetDisc(p, radius);
  }
  int i = 0;
  for (const S2Point& p : s2textformat::ParsePointsOrDie(avoid)) {
    s.AvoidDisc(p, radius, disc_on_left[i++]);
  }
  EXPECT_EQ(expected_result, s.Extend(s2textformat::MakePointOrDie(dst)))
      << "\nsrc = " << src << "\ndst = " << dst << "\ntarget = " << target
      << "\navoid = " << avoid;
}

TEST(S2PolylineSimplifier, Src) {
  S2PolylineSimplifier s;
  s.Init(S2Point(1, 0, 0));
  EXPECT_EQ(s.src(), S2Point(1, 0, 0));
}

TEST(S2PolylineSimplifier, Reuse) {
  // Check that Init() can be called more than once.
  S2PolylineSimplifier s;
  S1ChordAngle radius(S1Angle::Degrees(10));
  s.Init(S2Point(1, 0, 0));
  EXPECT_TRUE(s.TargetDisc(S2Point(1, 1, 0).Normalize(), radius));
  EXPECT_TRUE(s.TargetDisc(S2Point(1, 1, 0.1).Normalize(), radius));
  EXPECT_FALSE(s.Extend(S2Point(1, 1, 0.4).Normalize()));

  s.Init(S2Point(0, 1, 0));
  EXPECT_TRUE(s.TargetDisc(S2Point(1, 1, 0.3).Normalize(), radius));
  EXPECT_TRUE(s.TargetDisc(S2Point(1, 1, 0.2).Normalize(), radius));
  EXPECT_FALSE(s.Extend(S2Point(1, 1, 0).Normalize()));
}

TEST(S2PolylineSimplifier, NoConstraints) {
  // No constraints, dst == src.
  CheckSimplify("0:1", "0:1", "", "", {}, 0, true);

  // No constraints, dst != src.
  CheckSimplify("0:1", "1:0", "", "", {}, 0, true);

  // No constraints, (src, dst) longer than 90 degrees (not supported).
  CheckSimplify("0:0", "0:91", "", "", {}, 0, false);
}

TEST(S2PolylineSimplifier, TargetOnePoint) {
  // Three points on a straight line.  In theory zero tolerance should work,
  // but in practice there are floating point errors.
  CheckSimplify("0:0", "0:2", "0:1", "", {}, 1e-10, true);

  // Three points where the middle point is too far away.
  CheckSimplify("0:0", "0:2", "1:1", "", {}, 0.9, false);

  // A target disc that contains the source vertex.
  CheckSimplify("0:0", "0:2", "0:0.1", "", {}, 1.0, true);

  // A target disc that contains the destination vertex.
  CheckSimplify("0:0", "0:2", "0:2.1", "", {}, 1.0, true);
}

TEST(S2PolylineSimplifier, AvoidOnePoint) {
  // Three points on a straight line, attempting to avoid the middle point.
  CheckSimplify("0:0", "0:2", "", "0:1", {true}, 1e-10, false);

  // Three points where the middle point can be successfully avoided.
  CheckSimplify("0:0", "0:2", "", "1:1", {true}, 0.9, true);

  // Three points where the middle point is on the left, but where the client
  // requires the point to be on the right of the edge.
  CheckSimplify("0:0", "0:2", "", "1:1", {false}, 1e-10, false);

  // Check cases where the point to be avoided is behind the source vertex.
  // In this situation "disc_on_left" should not affect the result.
  CheckSimplify("0:0", "0:2", "", "1:-1", {false}, 1.4, true);
  CheckSimplify("0:0", "0:2", "", "1:-1", {true}, 1.4, true);
  CheckSimplify("0:0", "0:2", "", "-1:-1", {false}, 1.4, true);
  CheckSimplify("0:0", "0:2", "", "-1:-1", {true}, 1.4, true);
}

TEST(S2PolylineSimplifier, AvoidSeveralPoints) {
  // Tests that when several discs are avoided but none are targeted, the
  // range of acceptable edge directions is not represented a single interval.
  // (The output edge can proceed through any gap between the discs as long as
  // the "disc_on_left" criteria are satisfied.)
  //
  // This test involves 3 very small discs spaced 120 degrees apart around the
  // source vertex, where "disc_on_left" is true for all discs.  This means
  // that each disc blocks the 90 degrees of potential output directions just
  // to its left, leave 3 gaps measuring about 30 degrees each.
  for (const auto& dst : {"0:2", "1.732:-1", "-1.732:-1"}) {
    CheckSimplify("0:0", dst, "", "0.01:2, 1.732:-1.01, -1.732:-0.99",
                  {true, true, true}, 0.00001, true);

    // Also test that directions prohibited by "disc_on_left" are avoided.
    CheckSimplify("0:0", dst, "", "0.01:2, 1.732:-1.01, -1.732:-0.99",
                  {false, false, false}, 0.00001, false);
  }
}

TEST(S2PolylineSimplifier, TargetAndAvoid) {
  // Target several points that are separated from the proposed edge by about
  // 0.7 degrees, and avoid several points that are separated from the
  // proposed edge by about 1.4 degrees.
  CheckSimplify("0:0", "10:10", "2:3, 4:3, 7:8",
                "4:2, 7:5, 7:9", {true, true, false}, 1.0, true);

  // The same example, but one point to be targeted is 1.4 degrees away.
  CheckSimplify("0:0", "10:10", "2:3, 4:6, 7:8",
                "4:2, 7:5, 7:9", {true, true, false}, 1.0, false);

  // The same example, but one point to be avoided is 0.7 degrees away.
  CheckSimplify("0:0", "10:10", "2:3, 4:3, 7:8",
                "4:2, 6:5, 7:9", {true, true, false}, 1.0, false);
}

TEST(S2PolylineSimplifier, Precision) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "PRECISION",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  // This is a rough upper bound on both the error in constructing the disc
  // locations (i.e., S2::GetPointOnLine, etc.), and also on the
  // padding that S2PolylineSimplifier uses to ensure that its results are
  // conservative (i.e., the error calculated by GetSemiwidth).
  const S1Angle kMaxError = S1Angle::Radians(25 * DBL_EPSILON);

  // We repeatedly generate a random edge.  We then target several discs that
  // barely overlap the edge, and avoid several discs that barely miss the
  // edge.  About half the time, we choose one disc and make it slightly too
  // large or too small so that targeting fails.
  constexpr int kIters = 1000;  // Passes with 1 million iterations.
  S2PolylineSimplifier simplifier;
  for (int iter = 0; iter < kIters; ++iter) {
    S2Point src = s2random::Point(bitgen);
    simplifier.Init(src);
    S2Point dst =
        S2::GetPointOnLine(src, s2random::Point(bitgen),
                           S1Angle::Radians(absl::Uniform(bitgen, 0.0, 1.0)));
    S2Point n = S2::RobustCrossProd(src, dst).Normalize();

    // If bad_disc >= 0, then we make targeting fail for that disc.
    constexpr int kNumDiscs = 5;
    int bad_disc = absl::Uniform(bitgen, -kNumDiscs, kNumDiscs);
    for (int i = 0; i < kNumDiscs; ++i) {
      // The center of the disc projects to a point that is the given fraction
      // "f" along the edge (src, dst).  If f < 0, the center is located
      // behind "src" (in order to test this case).
      double f = absl::Uniform(bitgen, -0.5, 1.0);
      S2Point a = ((1 - f) * src + f * dst).Normalize();
      S1Angle r = S1Angle::Radians(absl::Uniform(bitgen, 0.0, 1.0));
      bool on_left = absl::Bernoulli(bitgen, 0.5);
      S2Point x = S2::GetPointOnLine(a, on_left ? n : -n, r);
      // If the disc is behind "src", adjust its radius so that it just
      // touches "src" rather than just touching the line through (src, dst).
      if (f < 0) r = S1Angle(src, x);
      // We grow the radius slightly if we want to target the disc and shrink
      // it otherwise, *unless* we want targeting to fail for this disc, in
      // which case these actions are reversed.
      bool avoid = absl::Bernoulli(bitgen, 0.5);
      bool grow_radius = (avoid == (i == bad_disc));
      S1ChordAngle radius(grow_radius ? r + kMaxError : r - kMaxError);
      if (avoid) {
        simplifier.AvoidDisc(x, radius, on_left);
      } else {
        simplifier.TargetDisc(x, radius);
      }
    }
    // The result is true iff all the disc constraints were satisfiable.
    EXPECT_EQ(bad_disc < 0, simplifier.Extend(dst));
  }
}
