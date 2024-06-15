// Copyright Google Inc. All Rights Reserved.
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

#include "s2/gmock_matchers.h"

#include <limits>
#include <memory>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "s2/s1angle.h"
#include "s2/s2latlng.h"
#include "s2/s2latlng_rect.h"
#include "s2/s2loop.h"
#include "s2/s2point.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"
#include "s2/s2text_format.h"

namespace S2 {

// Tests that the correct matchers are being invoked. Does not test for
// the correctness of the s2 comparison methods, those are assumed to be
// correct.

using S2::LatLngEquals;
using S2::PointApproxEquals;
using S2::PointEquals;
using S2::PolygonApproxEquals;
using S2::PolygonBoundaryApproxEquals;
using S2::PolygonBoundaryEquals;
using S2::PolygonBoundaryNear;
using S2::PolygonEquals;
using S2::PolylineApproxEquals;
using S2::PolylineEquals;
using s2textformat::MakeLatLngOrDie;
using s2textformat::MakeLatLngRectOrDie;
using s2textformat::MakeLoopOrDie;
using s2textformat::MakePointOrDie;
using s2textformat::MakePolygonOrDie;
using s2textformat::MakePolylineOrDie;
using std::string;
using ::testing::Not;

TEST(S2testmatchersTest, PointEquals) {
  auto expected = MakePointOrDie("10:10");

  // As object.
  EXPECT_THAT(s2textformat::MakePointOrDie("10:10"), PointEquals(expected));

  // As string.
  EXPECT_THAT("10:10", PointEquals("10:10"));

  // Negation.
  EXPECT_THAT("10.01:10", Not(PointEquals("10:10")));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2PointMatcher m("10:10", S1Angle::Zero());
  m.MatchAndExplain("10.01:10", &listener);
  EXPECT_EQ("angle in degrees: 0.0100000", listener.str());
}

TEST(S2testmatchersTest, PointApproxEquals) {
  auto expected = MakePointOrDie("10:10");

  // As object.
  EXPECT_THAT(MakePointOrDie("10:10"),
              PointApproxEquals(expected, S1Angle::Radians(1e-10)));

  // As string.
  EXPECT_THAT("10:10", PointApproxEquals("10:10", S1Angle::Radians(1e-10)));

  // Negation.
  EXPECT_THAT("11:10",
              Not(PointApproxEquals("10:10", S1Angle::Radians(1e-10))));

  // Explain.
  {
    ::testing::StringMatchResultListener listener;
    S2PointMatcher m("11:10", S1Angle::Radians(0.01));
    m.MatchAndExplain("10:10", &listener);
    EXPECT_EQ("angle in degrees: 1.0000000", listener.str());
  }

  // Default (uninitialized) point is not approximately equal to anything.
  S2Point default_point;
  EXPECT_THAT(default_point,
              Not(PointApproxEquals("10:10", S1Angle::Radians(1e-10))));
  {
    ::testing::StringMatchResultListener listener;
    S2PointMatcher m("10:10", S1Angle::Radians(0.01));
    m.MatchAndExplain(default_point, &listener);
    EXPECT_EQ("point is uninitialized or almost zero", listener.str());
  }

  // A point represented by a vector that is almost the zero-vector should also
  // not equal anything.
  S2Point nearly_zero_vector(std::numeric_limits<double>::min(),
                             std::numeric_limits<double>::min(),
                             std::numeric_limits<double>::min());
  EXPECT_THAT(nearly_zero_vector,
              Not(PointApproxEquals("10:10", S1Angle::Radians(1e-10))));
}

TEST(S2testmatchersTest, PolylineEquals) {
  auto expected = MakePolylineOrDie("10:10,20:20");
  auto value = MakePolylineOrDie("10:10,20:20");

  // As object.
  EXPECT_THAT(*value, PolylineEquals(*expected));

  // As string.
  EXPECT_THAT("10:10,20:20", PolylineEquals("10:10,20:20"));

  // Negation.
  EXPECT_THAT("10:10,5:5", Not(PolylineEquals("10:10,20:20")));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2PolylineMatcher m("10:10,20:20",
                      S2PolylineMatcher::Options(
                          S2PolylineMatcher::Options::EQUALS, S1Angle::Zero()));
  m.MatchAndExplain("10:10,5:5", &listener);
  EXPECT_EQ("S2Polylines don't match", listener.str());
}

TEST(S2testmatchersTest, PolylineApproxEquals) {
  auto expected = MakePolylineOrDie("10:10,20:20");
  auto value = MakePolylineOrDie("10:10,20:20");

  // As object.
  EXPECT_THAT(*value, PolylineApproxEquals(*expected, S1Angle::Radians(1e-10)));
  // As string.
  EXPECT_THAT("10:10,20:20",
              PolylineApproxEquals("10:10,20:20", S1Angle::Radians(1e-10)));

  // Negation.
  EXPECT_THAT("10:10,5:5", Not(PolylineApproxEquals("10:10,20:20",
                                                    S1Angle::Radians(1e-10))));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2PolylineMatcher m(
      "10:10,20:20",
      S2PolylineMatcher::Options(S2PolylineMatcher::Options::APPROX_EQUALS,
                                 S1Angle::Degrees(0.1)));
  m.MatchAndExplain("10:10,5:5", &listener);
  EXPECT_EQ("S2Polylines don't match with tolerance in degrees 0.1000000",
            listener.str());
}

TEST(S2testmatchersTest, PolygonEquals) {
  auto expected = MakePolygonOrDie("10:10,10:0,0:0");
  auto value = MakePolygonOrDie("10:10,10:0,0:0");

  // As object.
  EXPECT_THAT(*value, PolygonEquals(*expected));
  // As string.
  EXPECT_THAT("10:10,10:0,0:0", PolygonEquals("10:10,10:0,0:0"));

  // Negation.
  EXPECT_THAT("10:10,10:0,0:0", Not(PolygonEquals("10:10,10:0,0:1")));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2PolygonMatcher m("10:10,10:0,0:0",
                     S2PolygonMatcher::Options(
                         S2PolygonMatcher::Options::EQUALS, S1Angle::Zero()));
  m.MatchAndExplain("10:10,10:0,0:1", &listener);
  EXPECT_EQ("S2Polygons don't match", listener.str());
}

TEST(S2testmatchersTest, PolygonApproxEquals) {
  auto expected = MakePolygonOrDie("10:10,10:0,0:0");
  auto value = MakePolygonOrDie("10:10,10:0,0:0");

  // As object.
  EXPECT_THAT(*value, PolygonApproxEquals(*expected, S1Angle::Radians(1e-10)));
  // As string.
  EXPECT_THAT("10:10,10:0,0:0",
              PolygonApproxEquals("10:10,10:0,0:0", S1Angle::Degrees(0.1)));

  // Negation.
  EXPECT_THAT("10:10,10:0,0:0", Not(PolygonApproxEquals(
                                    "10:10,10:0,0:1", S1Angle::Degrees(0.1))));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2PolygonMatcher m(
      "10:10,10:0,0:0",
      S2PolygonMatcher::Options(S2PolygonMatcher::Options::APPROX_EQUALS,
                                S1Angle::Degrees(0.1)));
  m.MatchAndExplain("10:10,10:0,0:1", &listener);
  EXPECT_EQ("S2Polygons don't match with tolerance in degrees 0.1000000",
            listener.str());
}

TEST(S2testmatchersTest, PolygonBoundaryEquals) {
  auto expected = MakePolygonOrDie("10:10,10:0,0:0");
  auto value = MakePolygonOrDie("10:10,10:0,0:0");

  // As object.
  EXPECT_THAT(*value, PolygonBoundaryEquals(*expected));
  // As string.
  EXPECT_THAT("10:10,10:0,0:0", PolygonBoundaryEquals("10:10,10:0,0:0"));

  // Negation.
  EXPECT_THAT("10:10,10:0,0:0", Not(PolygonBoundaryEquals("10:10,10:0,0:1")));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2PolygonMatcher m(
      "10:10,10:0,0:0",
      S2PolygonMatcher::Options(S2PolygonMatcher::Options::BOUNDARY_EQUALS,
                                S1Angle::Zero()));
  m.MatchAndExplain("10:10,10:0,0:1", &listener);
  EXPECT_EQ("S2Polygons don't match", listener.str());
}

TEST(S2testmatchersTest, PolygonBoundaryApproxEquals) {
  auto expected = MakePolygonOrDie("10:10,10:0,0:0");
  auto value = MakePolygonOrDie("10:10,10:0,0:0");

  // As object.
  EXPECT_THAT(*value,
              PolygonBoundaryApproxEquals(*expected, S1Angle::Radians(1e-10)));
  // As string.
  EXPECT_THAT("10:10,10:0,0:0", PolygonBoundaryApproxEquals(
                                    "10:10,10:0,0:0", S1Angle::Radians(1e-10)));

  // Negation.
  EXPECT_THAT("10:10,10:0,0:0", Not(PolygonBoundaryApproxEquals(
                                    "10:10,10:0,0:1", S1Angle::Degrees(0.1))));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2PolygonMatcher m("10:10,10:0,0:0",
                     S2PolygonMatcher::Options(
                         S2PolygonMatcher::Options::BOUNDARY_APPROX_EQUALS,
                         S1Angle::Degrees(0.1)));
  m.MatchAndExplain("10:10,10:0,0:1", &listener);
  EXPECT_EQ("S2Polygons don't match with tolerance in degrees 0.1000000",
            listener.str());
}

TEST(S2testmatchersTest, PolygonBoundaryNear) {
  auto expected = MakePolygonOrDie("10:10,10:0,0:0");
  auto value = MakePolygonOrDie("10:10,10:0,0:0");

  // As object.
  EXPECT_THAT(*value, PolygonBoundaryNear(*expected, S1Angle::Radians(1e-10)));
  // As string.
  EXPECT_THAT("10:10,10:0,0:0",
              PolygonBoundaryNear("10:10,10:0,0:0", S1Angle::Radians(1e-10)));

  // Negation.
  EXPECT_THAT("10:10,10:0,0:0", Not(PolygonBoundaryNear(
                                    "10:10,10:0,0:1", S1Angle::Degrees(0.1))));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2PolygonMatcher m(
      "10:10,10:0,0:0",
      S2PolygonMatcher::Options(S2PolygonMatcher::Options::BOUNDARY_NEAR,
                                S1Angle::Degrees(0.1)));
  m.MatchAndExplain("10:10,10:0,0:1", &listener);
  EXPECT_EQ("S2Polygons don't match with tolerance in degrees 0.1000000",
            listener.str());
}

TEST(S2testmatchersTest, LoopEquals) {
  auto expected = MakeLoopOrDie("10:10,10:0,0:0");
  auto value = MakeLoopOrDie("10:10,10:0,0:0");

  // As object.
  EXPECT_THAT(*value, LoopEquals(*expected));
  // As string.
  EXPECT_THAT("10:10,10:0,0:0", LoopEquals("10:10,10:0,0:0"));

  // Negation.
  EXPECT_THAT("10:10,10:0,0:0", Not(LoopEquals("10:10,10:0,0:1")));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2LoopMatcher m(
      "10:10,10:0,0:0",
      S2LoopMatcher::Options(S2LoopMatcher::Options::EQUALS, S1Angle::Zero()));
  m.MatchAndExplain("10:10,10:0,0:1", &listener);
  EXPECT_EQ("S2Loops don't match", listener.str());
}

TEST(S2testmatchersTest, LoopBoundaryEquals) {
  auto expected = MakeLoopOrDie("10:10,10:0,0:0");
  auto value = MakeLoopOrDie("10:10,10:0,0:0");

  // As object.
  EXPECT_THAT(*value, LoopBoundaryEquals(*expected));
  // As string.
  EXPECT_THAT("10:10,10:0,0:0", LoopBoundaryEquals("10:10,10:0,0:0"));

  // Negation.
  EXPECT_THAT("10:10,10:0,0:0", Not(LoopBoundaryEquals("10:10,10:0,0:1")));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2LoopMatcher m("10:10,10:0,0:0", S2LoopMatcher::Options(
                                        S2LoopMatcher::Options::BOUNDARY_EQUALS,
                                        S1Angle::Zero()));
  m.MatchAndExplain("10:10,10:0,0:1", &listener);
  EXPECT_EQ("S2Loops don't match", listener.str());
}

TEST(S2testmatchersTest, LoopBoundaryEqualsRotatedVertices) {
  string expected_str = "10:10,10:0,0:0";
  auto expected = MakeLoopOrDie(expected_str);
  auto value = MakeLoopOrDie("10:0,0:0,10:10");

  // As object.
  EXPECT_THAT(*value, LoopBoundaryEquals(*expected));
  // As string.
  EXPECT_THAT("10:10,10:0,0:0", LoopBoundaryEquals(expected_str));

  // Negation.
  EXPECT_THAT("10:10,10:0,0:0", Not(LoopBoundaryEquals("10:10,10:0,0:1")));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2LoopMatcher m(
      "10:10,10:0,0:0",
      S2LoopMatcher::Options(S2LoopMatcher::Options::EQUALS, S1Angle::Zero()));
  m.MatchAndExplain("10:10,10:0,0:1", &listener);
  EXPECT_EQ("S2Loops don't match", listener.str());
}

TEST(S2testmatchersTest, LoopBoundaryApproxEquals) {
  auto expected = MakeLoopOrDie("10:10,10:0,0:0");
  auto value = MakeLoopOrDie("10:10,10:0,0:0");

  // As object.
  EXPECT_THAT(*value,
              LoopBoundaryApproxEquals(*expected, S1Angle::Radians(1e-10)));
  // As string.
  EXPECT_THAT("10:10,10:0,0:0", LoopBoundaryApproxEquals(
                                    "10:10,10:0,0:0", S1Angle::Radians(1e-10)));

  // Negation.
  EXPECT_THAT("10:10,10:0,0:0", Not(LoopBoundaryApproxEquals(
                                    "10:10,10:0,0:1", S1Angle::Degrees(0.1))));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2LoopMatcher m(
      "10:10,10:0,0:0",
      S2LoopMatcher::Options(S2LoopMatcher::Options::BOUNDARY_APPROX_EQUALS,
                             S1Angle::Degrees(0.1)));
  m.MatchAndExplain("10:10,10:0,0:1", &listener);
  EXPECT_EQ("S2Loops don't match with tolerance in degrees 0.1000000",
            listener.str());
}

TEST(S2testmatchersTest, LoopBoundaryNear) {
  auto expected = MakeLoopOrDie("10:10,10:0,0:0");
  auto value = MakeLoopOrDie("10:10,10:0,0:0");

  // As object.
  EXPECT_THAT(*value, LoopBoundaryNear(*expected, S1Angle::Radians(1e-10)));
  // As string.
  EXPECT_THAT("10:10,10:0,0:0",
              LoopBoundaryNear("10:10,10:0,0:0", S1Angle::Radians(1e-10)));

  // Negation.
  EXPECT_THAT("10:10,10:0,0:0",
              Not(LoopBoundaryNear("10:10,10:0,0:1", S1Angle::Degrees(0.1))));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2LoopMatcher m("10:10,10:0,0:0",
                  S2LoopMatcher::Options(S2LoopMatcher::Options::BOUNDARY_NEAR,
                                         S1Angle::Degrees(0.1)));
  m.MatchAndExplain("10:10,10:0,0:1", &listener);
  EXPECT_EQ("S2Loops don't match with tolerance in degrees 0.1000000",
            listener.str());
}

TEST(S2testmatchersTest, LatLngRectEquals) {
  auto expected = MakeLatLngRectOrDie("0:0,10:10");
  auto value = MakeLatLngRectOrDie("0:0,10:10");

  // As object.
  EXPECT_THAT(value, LatLngRectEquals(expected));

  // As string.
  EXPECT_THAT("0:0,10:10", LatLngRectEquals("0:0,10:10"));

  // Negation.
  EXPECT_THAT("0:0,10:10", Not(LatLngRectEquals("0:0,11:10")));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2LatLngRectMatcher m(
      "0:0,10:10", S2LatLngRectMatcher::Options(
                       S2LatLngRectMatcher::Options::EQUALS, S1Angle::Zero()));
  m.MatchAndExplain("0:0,11:10", &listener);
  EXPECT_EQ("S2LatLngRects don't match", listener.str());
}

TEST(S2testmatchersTest, LatLngRectApproxEquals) {
  auto expected = MakeLatLngRectOrDie("0:0,10:10");
  auto value = MakeLatLngRectOrDie("0:0,10:10");

  // As object.
  EXPECT_THAT(value, LatLngRectApproxEquals(expected, S1Angle::Degrees(0.1)));

  // As string.
  EXPECT_THAT("0:0,10:10",
              LatLngRectApproxEquals("0:0,10:10", S1Angle::Degrees(0.1)));

  // Negation.
  EXPECT_THAT("0:0,10:10",
              Not(LatLngRectApproxEquals("0:0,11:10", S1Angle::Degrees(0.1))));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2LatLngRectMatcher m(
      "0:0,10:10",
      S2LatLngRectMatcher::Options(S2LatLngRectMatcher::Options::APPROX_EQUALS,
                                   S1Angle::Degrees(0.1)));
  m.MatchAndExplain("0:0,11:10", &listener);
  EXPECT_EQ("S2LatLngRects don't match with tolerance in degrees 0.1000000",
            listener.str());
}

TEST(S2testmatchersTest, LatLngRectApproxEquals2) {
  auto expected = MakeLatLngRectOrDie("0:0,10:10");
  auto value = MakeLatLngRectOrDie("0:0,10:10");

  auto tolerance = S2LatLng::FromDegrees(0.1, 0.2);

  // As object.
  EXPECT_THAT(value, LatLngRectApproxEquals(expected, tolerance));

  // As string.
  EXPECT_THAT("0:0,10:10", LatLngRectApproxEquals("0:0,10:10", tolerance));

  // Negation.
  EXPECT_THAT("0:0,10:10", Not(LatLngRectApproxEquals("0:0,11:10", tolerance)));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2LatLngRectMatcher m(
      "0:0,10:10",
      S2LatLngRectMatcher::Options(
          S2LatLngRectMatcher::Options::APPROX_EQUALS_2, tolerance));
  m.MatchAndExplain("0:0,11:10", &listener);
  EXPECT_EQ(
      "S2LatLngRects don't match with tolerance in degrees [0.1000000, "
      "0.2000000]",
      listener.str());
}

TEST(S2testmatchersTest, LatLngEquals) {
  auto expected = MakeLatLngOrDie("10:10");
  auto value = MakeLatLngOrDie("10:10");

  // As object.
  EXPECT_THAT(value, LatLngEquals(expected));

  // As string.
  EXPECT_THAT("10:10", LatLngEquals("10:10"));

  // Negation.
  EXPECT_THAT("10:10", Not(LatLngEquals("11:10")));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2LatLngMatcher m("10:10",
                    S2LatLngMatcher::Options(S2LatLngMatcher::Options::EQUALS,
                                             S1Angle::Zero()));
  m.MatchAndExplain("11:10", &listener);
  EXPECT_EQ("S2LatLngs don't match", listener.str());
}

TEST(S2testmatchersTest, LatLngApproxEquals) {
  auto expected = MakeLatLngOrDie("10:10");
  auto value = MakeLatLngOrDie("10:10");

  // As object.
  EXPECT_THAT(value, LatLngApproxEquals(expected, S1Angle::Degrees(0.1)));

  // As string, testing tolerances in both directions
  EXPECT_THAT("10.095:9.95",
              LatLngApproxEquals("10:10", S1Angle::Degrees(0.1)));
  EXPECT_THAT("9.95:10.095",
              LatLngApproxEquals("10:10", S1Angle::Degrees(0.1)));

  // Negation.
  EXPECT_THAT("10:10", Not(LatLngApproxEquals("11:10", S1Angle::Degrees(0.1))));

  // Explain.
  ::testing::StringMatchResultListener listener;
  S2LatLngMatcher m(
      "10:10", S2LatLngMatcher::Options(S2LatLngMatcher::Options::APPROX_EQUALS,
                                        S1Angle::Degrees(0.1)));
  m.MatchAndExplain("11:10", &listener);
  EXPECT_EQ("S2LatLngs don't match with tolerance in degrees 0.1000000",
            listener.str());
}

}  // namespace S2
