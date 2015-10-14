// Copyright 2015 Google Inc. All Rights Reserved.
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

#include "s1chordangle.h"

#include <limits>

#include <gtest/gtest.h>
#include "s1angle.h"
#include "s2testing.h"

using std::numeric_limits;

TEST(S1ChordAngle, DefaultConstructor) {
  // Check that the default constructor returns an angle of 0.
  S1ChordAngle a;
  EXPECT_EQ(S1ChordAngle::Zero(), a);
}

TEST(S1ChordAngle, TwoPointConstructor) {
  for (int iter = 0; iter < 100; ++iter) {
    S2Point x, y, z;
    S2Testing::GetRandomFrame(&x, &y, &z);
    EXPECT_EQ(S1Angle::Zero(), S1ChordAngle(z, z).ToAngle());
    EXPECT_NEAR(M_PI, S1ChordAngle(-z, z).ToAngle().radians(), 1e-7);
    EXPECT_DOUBLE_EQ(M_PI_2, S1ChordAngle(x, z).ToAngle().radians());
    S2Point w = (y + z).Normalize();
    EXPECT_DOUBLE_EQ(M_PI_4, S1ChordAngle(w, z).ToAngle().radians());
  }
}

TEST(S1ChordAngle, FromLength2) {
  EXPECT_EQ(0, S1ChordAngle::FromLength2(0).ToAngle().degrees());
  EXPECT_DOUBLE_EQ(60, S1ChordAngle::FromLength2(1).ToAngle().degrees());
  EXPECT_DOUBLE_EQ(90, S1ChordAngle::FromLength2(2).ToAngle().degrees());
  EXPECT_EQ(180, S1ChordAngle::FromLength2(4).ToAngle().degrees());
  EXPECT_EQ(180, S1ChordAngle::FromLength2(5).ToAngle().degrees());
}

TEST(S1ChordAngle, Zero) {
  EXPECT_EQ(S1Angle::Zero(), S1ChordAngle::Zero().ToAngle());
}

TEST(S1ChordAngle, Straight) {
  EXPECT_EQ(S1Angle::Degrees(180), S1ChordAngle::Straight().ToAngle());
}

TEST(S1ChordAngle, Infinity) {
  EXPECT_LT(S1ChordAngle::Straight(), S1ChordAngle::Infinity());
  EXPECT_EQ(S1ChordAngle::Infinity(), S1ChordAngle::Infinity());
  EXPECT_EQ(S1Angle::Infinity(), S1ChordAngle::Infinity().ToAngle());
}

TEST(S1ChordAngle, Negative) {
  EXPECT_LT(S1ChordAngle::Negative(), S1ChordAngle::Zero());
  EXPECT_EQ(S1ChordAngle::Negative(), S1ChordAngle::Negative());
  EXPECT_LT(S1ChordAngle::Negative().ToAngle(), S1Angle::Zero());
}

TEST(S1ChordAngle, Predicates) {
  EXPECT_TRUE(S1ChordAngle::Zero().is_zero());
  EXPECT_FALSE(S1ChordAngle::Zero().is_negative());
  EXPECT_FALSE(S1ChordAngle::Zero().is_special());
  EXPECT_FALSE(S1ChordAngle::Straight().is_special());
  EXPECT_TRUE(S1ChordAngle::Negative().is_negative());
  EXPECT_TRUE(S1ChordAngle::Negative().is_special());
  EXPECT_TRUE(S1ChordAngle::Infinity().is_infinity());
  EXPECT_TRUE(S1ChordAngle::Infinity().is_special());
}

TEST(S1ChordAngle, ToFromS1Angle) {
  EXPECT_EQ(0, S1ChordAngle(S1Angle::Zero()).ToAngle().radians());
  EXPECT_EQ(4, S1ChordAngle(S1Angle::Radians(M_PI)).length2());
  EXPECT_EQ(M_PI, S1ChordAngle(S1Angle::Radians(M_PI)).ToAngle().radians());
  EXPECT_EQ(S1Angle::Infinity(), S1ChordAngle(S1Angle::Infinity()).ToAngle());
  EXPECT_LT(S1ChordAngle(S1Angle::Radians(-1)).ToAngle().radians(), 0);
  EXPECT_DOUBLE_EQ(1.0,
                   S1ChordAngle(S1Angle::Radians(1.0)).ToAngle().radians());
}

TEST(S1ChordAngle, Arithmetic) {
  S1ChordAngle zero = S1ChordAngle::Zero();
  S1ChordAngle degree30 = S1ChordAngle(S1Angle::Degrees(30));
  S1ChordAngle degree60 = S1ChordAngle(S1Angle::Degrees(60));
  S1ChordAngle degree90 = S1ChordAngle(S1Angle::Degrees(90));
  S1ChordAngle degree120 = S1ChordAngle(S1Angle::Degrees(120));
  S1ChordAngle degree180 = S1ChordAngle::Straight();
  EXPECT_EQ(0, (zero + zero).ToAngle().degrees());
  EXPECT_EQ(0, (zero - zero).ToAngle().degrees());
  EXPECT_EQ(0, (degree60 - degree60).ToAngle().degrees());
  EXPECT_EQ(0, (degree180 - degree180).ToAngle().degrees());
  EXPECT_EQ(0, (zero - degree60).ToAngle().degrees());
  EXPECT_EQ(0, (degree30 - degree90).ToAngle().degrees());
  EXPECT_DOUBLE_EQ(60, (degree60 + zero).ToAngle().degrees());
  EXPECT_DOUBLE_EQ(60, (degree60 - zero).ToAngle().degrees());
  EXPECT_DOUBLE_EQ(60, (zero + degree60).ToAngle().degrees());
  EXPECT_DOUBLE_EQ(90, (degree30 + degree60).ToAngle().degrees());
  EXPECT_DOUBLE_EQ(90, (degree60 + degree30).ToAngle().degrees());
  EXPECT_DOUBLE_EQ(60, (degree90 - degree30).ToAngle().degrees());
  EXPECT_DOUBLE_EQ(30, (degree90 - degree60).ToAngle().degrees());
  EXPECT_EQ(180, (degree180 + zero).ToAngle().degrees());
  EXPECT_EQ(180, (degree180 - zero).ToAngle().degrees());
  EXPECT_EQ(180, (degree90 + degree90).ToAngle().degrees());
  EXPECT_EQ(180, (degree120 + degree90).ToAngle().degrees());
  EXPECT_EQ(180, (degree120 + degree120).ToAngle().degrees());
  EXPECT_EQ(180, (degree30 + degree180).ToAngle().degrees());
  EXPECT_EQ(180, (degree180 + degree180).ToAngle().degrees());
}

TEST(S1ChordAngle, Trigonometry) {
  static int const kIters = 20;
  for (int iter = 0; iter <= kIters; ++iter) {
    double radians = M_PI * iter / kIters;
    S1ChordAngle angle(S1Angle::Radians(radians));
    EXPECT_NEAR(sin(radians), sin(angle), 1e-15);
    EXPECT_NEAR(cos(radians), cos(angle), 1e-15);
    // Since the tan(x) is unbounded near Pi/4, we map the result back to an
    // angle before comparing.  (The assertion is that the result is equal to
    // the tangent of a nearby angle.)
    EXPECT_NEAR(atan(tan(radians)), atan(tan(angle)), 1e-15);
  }

  // Unlike S1Angle, S1ChordAngle can represent 90 and 180 degrees exactly.
  S1ChordAngle angle90 = S1ChordAngle::FromLength2(2);
  S1ChordAngle angle180 = S1ChordAngle::FromLength2(4);
  EXPECT_EQ(1, sin(angle90));
  EXPECT_EQ(0, cos(angle90));
  EXPECT_EQ(numeric_limits<double>::infinity(), tan(angle90));
  EXPECT_EQ(0, sin(angle180));
  EXPECT_EQ(-1, cos(angle180));
  EXPECT_EQ(0, tan(angle180));
}
