// Copyright 2017 Google Inc. All Rights Reserved.
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

#include "s2/s2textformat.h"

#include <vector>
#include <gtest/gtest.h>
#include "s2/third_party/absl/strings/str_split.h"
#include "s2/s1angle.h"
#include "s2/s2latlng.h"
#include "s2/s2loop.h"
#include "s2/s2polyline.h"
#include "s2/s2shapeutil.h"
#include "s2/s2testing.h"

using std::unique_ptr;
using std::vector;

using LaxPolygon = s2shapeutil::LaxPolygon;

namespace {

static int const kIters = 10000;

// Verify that s2textformat::ToString() formats the given lat/lng with at most
// "max_digits" after the decimal point and has no trailing zeros.
void ExpectMaxDigits(S2LatLng const& ll, int max_digits) {
  string result = s2textformat::ToString(ll.ToPoint());
  vector<string> values = strings::Split(result, ':', strings::SkipEmpty());
  EXPECT_EQ(2, values.size()) << result;
  for (auto const& value : values) {
    int num_digits = 0;
    if (value.find('.') != string::npos) {
      num_digits = value.size() - value.find('.') - 1;
      EXPECT_NE('0', value.back());
    }
    EXPECT_LE(num_digits, max_digits) << value;
  }
}

void ExpectString(string const& expected, S2LatLng const& ll) {
  EXPECT_EQ(expected, s2textformat::ToString(ll.ToPoint()));
}

TEST(ToString, SpecialCases) {
  ExpectString("0:0", S2LatLng::FromDegrees(0, 0));
  ExpectString("1e-20:1e-30", S2LatLng::FromDegrees(1e-20, 1e-30));
}

TEST(ToString, MinimalDigitsE5) {
  for (int iter = 0; iter < kIters; ++iter) {
    S2LatLng ll(S2Testing::RandomPoint());
    S2LatLng ll_e5 = S2LatLng::FromE5(ll.lat().e5(), ll.lng().e5());
    ExpectMaxDigits(ll_e5, 5);
  }
}

TEST(ToString, MinimalDigitsE6) {
  for (int iter = 0; iter < kIters; ++iter) {
    S2LatLng ll(S2Testing::RandomPoint());
    S2LatLng ll_e6 = S2LatLng::FromE6(ll.lat().e6(), ll.lng().e6());
    ExpectMaxDigits(ll_e6, 6);
  }
}

TEST(ToString, MinimalDigitsE7) {
  ExpectMaxDigits(S2LatLng::FromDegrees(0, 0), 7);
  for (int iter = 0; iter < kIters; ++iter) {
    S2LatLng ll(S2Testing::RandomPoint());
    S2LatLng ll_e7 = S2LatLng::FromE7(ll.lat().e7(), ll.lng().e7());
    ExpectMaxDigits(ll_e7, 7);
  }
}

TEST(ToString, MinimalDigitsDoubleConstants) {
  // Verify that points specified as floating-point literals in degrees using
  // up to 10 digits after the decimal point are formatted with the minimal
  // number of digits.
  for (int iter = 0; iter < kIters; ++iter) {
    int max_digits = S2Testing::rnd.Uniform(11);
    int64 scale = MathUtil::FastInt64Round(pow(10, max_digits));
    int64 lat = MathUtil::FastInt64Round(
        S2Testing::rnd.UniformDouble(-90 * scale, 90 * scale));
    int64 lng = MathUtil::FastInt64Round(
        S2Testing::rnd.UniformDouble(-180 * scale, 180 * scale));
    S2LatLng ll = S2LatLng::FromDegrees(lat / static_cast<double>(scale),
                                        lng / static_cast<double>(scale));
    ExpectMaxDigits(ll, max_digits);
  }
}

TEST(ToString, EmptyLoop) {
  S2Loop loop;
  EXPECT_EQ("", s2textformat::ToString(loop));
}

TEST(ToString, EmptyPolyline) {
  S2Polyline polyline;
  EXPECT_EQ("", s2textformat::ToString(polyline));
}

TEST(ToString, EmptyPointVector) {
  vector<S2Point> points;
  EXPECT_EQ("", s2textformat::ToString(points));
}

TEST(MakeLaxPolygon, Empty) {
  unique_ptr<LaxPolygon> shape(s2textformat::MakeLaxPolygon(""));
  EXPECT_EQ(0, shape->num_loops());
}

TEST(MakeLaxPolygon, Full) {
  unique_ptr<LaxPolygon> shape(s2textformat::MakeLaxPolygon("full"));
  EXPECT_EQ(1, shape->num_loops());
  EXPECT_EQ(0, shape->num_loop_vertices(0));
}

TEST(MakeLaxPolygon, FullWithHole) {
  unique_ptr<LaxPolygon> shape(s2textformat::MakeLaxPolygon("full; 0:0"));
  EXPECT_EQ(2, shape->num_loops());
  EXPECT_EQ(0, shape->num_loop_vertices(0));
  EXPECT_EQ(1, shape->num_loop_vertices(1));
  EXPECT_EQ(1, shape->num_edges());
}

}  // namespace
