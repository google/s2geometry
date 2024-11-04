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

#include <cfloat>
#include <cmath>
#include <memory>
#include <ostream>
#include <string>

#include <gtest/gtest.h>
#include "absl/log/absl_check.h"
#include "absl/strings/string_view.h"
#include "s2/s1angle.h"
#include "s2/s2latlng.h"
#include "s2/s2latlng_rect.h"
#include "s2/s2loop.h"
#include "s2/s2point.h"
#include "s2/s2pointutil.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"
#include "s2/s2text_format.h"

using absl::string_view;
using std::string;
using std::unique_ptr;

// Printing methods need to be in the global namespace.

void PrintTo(const S2Point& p, std::ostream* os) {
  *os << s2textformat::ToString(p);
}

void PrintTo(const S2Polyline& l, std::ostream* os) {
  *os << s2textformat::ToString(l);
}

void PrintTo(const S2Polygon& p, std::ostream* os) {
  *os << s2textformat::ToString(p);
}

void PrintTo(const S2LatLngRect& r, std::ostream* os) {
  *os << s2textformat::ToString(r);
}

void PrintTo(const S2Loop& l, std::ostream* os) {
  *os << s2textformat::ToString(l);
}

namespace S2 {

namespace {

// Returns whether the point is approximately the zero-length vector.
bool IsAlmostZeroLength(const S2Point& p) {
  // Implementation is similar to S2::IsUnitLength. See comments there.
  // TODO(ericv): Verify if this is actually a reasonable value.
  return fabs(p.Norm2()) <= 5 * DBL_EPSILON;
}

}  // namespace

// S2PointMatcher methods.

S2PointMatcher::S2PointMatcher(const S2Point& expected, S1Angle tolerance)
    : expected_(expected), tolerance_(tolerance) {
  ABSL_CHECK(!IsAlmostZeroLength(expected));
}

bool S2PointMatcher::MatchAndExplain(
    string_view s, ::testing::MatchResultListener* listener) const {
  return MatchAndExplain(s2textformat::MakePointOrDie(s), listener);
}

bool S2PointMatcher::MatchAndExplain(
    const S2Point& g, ::testing::MatchResultListener* listener) const {
  // S2::ApproxEquals will either check fail (debug mode) to unconditionally
  // return true if either point is the zero-length vector (default S2Point), so
  // handle that explicitly.
  if (IsAlmostZeroLength(g)) {
    *listener << "point is uninitialized or almost zero";
    return false;
  }
  if ((tolerance_ == S1Angle::Zero() && expected_ != g) ||
      !S2::ApproxEquals(g, expected_, tolerance_)) {
    *listener << "angle in degrees: " << S1Angle(g, expected_);
    return false;
  }

  return true;
}

// S2PolylineMatcher methods.

bool S2PolylineMatcher::MatchAndExplain(
    string_view s, ::testing::MatchResultListener* listener) const {
  unique_ptr<S2Polyline> line = s2textformat::MakePolylineOrDie(s);
  return MatchAndExplain(*line, listener);
}

bool S2PolylineMatcher::MatchAndExplain(
    const S2Polyline& g, ::testing::MatchResultListener* listener) const {
  switch (options_.type) {
    case Options::EQUALS: {
      if (!expected_->Equals(g)) {
        *listener << "S2Polylines don't match";
        return false;
      }
      break;
    }
    case Options::APPROX_EQUALS:
      if (!expected_->ApproxEquals(g, options_.tolerance)) {
        *listener << "S2Polylines don't match with tolerance in degrees "
                  << options_.tolerance;
        return false;
      }
      break;
  }

  return true;
}

// S2PolygonMatcher methods.

bool S2PolygonMatcher::MatchAndExplain(
    string_view s, ::testing::MatchResultListener* listener) const {
  unique_ptr<S2Polygon> polygon =
      s2textformat::MakePolygonOrDie(s, S2Debug::ALLOW);
  return MatchAndExplain(*polygon, listener);
}

bool S2PolygonMatcher::MatchAndExplain(
    const S2Polygon& g, ::testing::MatchResultListener* listener) const {
  switch (options_.type) {
    case Options::EQUALS: {
      if (!expected_->Equals(g)) {
        *listener << "S2Polygons don't match";
        return false;
      }
      break;
    }
    case Options::APPROX_EQUALS: {
      if (!expected_->ApproxEquals(g, options_.tolerance)) {
        *listener << "S2Polygons don't match with tolerance in degrees "
                  << options_.tolerance;
        return false;
      }
      break;
    }
    case Options::BOUNDARY_EQUALS: {
      if (!expected_->BoundaryEquals(g)) {
        *listener << "S2Polygons don't match";
        return false;
      }
      break;
    }
    case Options::BOUNDARY_APPROX_EQUALS: {
      if (!expected_->BoundaryApproxEquals(g, options_.tolerance)) {
        *listener << "S2Polygons don't match with tolerance in degrees "
                  << options_.tolerance;
        return false;
      }
      break;
    }
    case Options::BOUNDARY_NEAR: {
      if (!expected_->BoundaryNear(g, options_.tolerance)) {
        *listener << "S2Polygons don't match with tolerance in degrees "
                  << options_.tolerance;
        return false;
      }
      break;
    }
  }

  return true;
}

// S2LoopMatcher methods.

bool S2LoopMatcher::MatchAndExplain(
    string_view s, ::testing::MatchResultListener* listener) const {
  unique_ptr<S2Loop> loop = s2textformat::MakeLoopOrDie(s);
  return MatchAndExplain(*loop, listener);
}

bool S2LoopMatcher::MatchAndExplain(
    const S2Loop& g, ::testing::MatchResultListener* listener) const {
  switch (options_.type) {
    case Options::EQUALS: {
      if (!expected_->Equals(g)) {
        *listener << "S2Loops don't match";
        return false;
      }
      break;
    }
    case Options::BOUNDARY_EQUALS: {
      if (!expected_->BoundaryEquals(g)) {
        *listener << "S2Loops don't match";
        return false;
      }
      break;
    }
    case Options::BOUNDARY_APPROX_EQUALS: {
      if (!expected_->BoundaryApproxEquals(g, options_.tolerance)) {
        *listener << "S2Loops don't match with tolerance in degrees "
                  << options_.tolerance;
        return false;
      }
      break;
    }
    case Options::BOUNDARY_NEAR: {
      if (!expected_->BoundaryNear(g, options_.tolerance)) {
        *listener << "S2Loops don't match with tolerance in degrees "
                  << options_.tolerance;
        return false;
      }
      break;
    }
  }

  return true;
}

// S2LatLngRectMatcher methods.

bool S2LatLngRectMatcher::MatchAndExplain(
    string_view s, ::testing::MatchResultListener* listener) const {
  S2LatLngRect rect = s2textformat::MakeLatLngRectOrDie(s);
  return MatchAndExplain(rect, listener);
}

bool S2LatLngRectMatcher::MatchAndExplain(
    const S2LatLngRect& g, ::testing::MatchResultListener* listener) const {
  switch (options_.type) {
    case Options::EQUALS: {
      if (!(expected_ == g)) {
        *listener << "S2LatLngRects don't match";
        return false;
      }
      break;
    }
    case Options::APPROX_EQUALS: {
      if (!expected_.ApproxEquals(g, options_.tolerance)) {
        *listener << "S2LatLngRects don't match with tolerance in degrees "
                  << options_.tolerance;
        return false;
      }
      break;
    }
    case Options::APPROX_EQUALS_2: {
      if (!expected_.ApproxEquals(g, options_.tolerance_lat_lng)) {
        *listener << "S2LatLngRects don't match with tolerance in degrees "
                  << options_.tolerance_lat_lng;
        return false;
      }
      break;
    }
  }
  return true;
}

// S2LatLngMatcher methods.

bool S2LatLngMatcher::MatchAndExplain(
    string_view s, ::testing::MatchResultListener* listener) const {
  S2LatLng latlng = s2textformat::MakeLatLngOrDie(s);
  return MatchAndExplain(latlng, listener);
}

bool S2LatLngMatcher::MatchAndExplain(
    const S2LatLng& latlng, ::testing::MatchResultListener* listener) const {
  switch (options_.type) {
    case Options::EQUALS: {
      if (!(expected_ == latlng)) {
        *listener << "S2LatLngs don't match";
        return false;
      }
      break;
    }
    case Options::APPROX_EQUALS: {
      if (!expected_.ApproxEquals(latlng, options_.tolerance)) {
        *listener << "S2LatLngs don't match with tolerance in degrees "
                  << options_.tolerance;
        return false;
      }
      break;
    }
  }
  return true;
}

// Predefined matchers.

// Matchers for S2Points provided as objects.

::testing::PolymorphicMatcher<S2PointMatcher> PointEquals(
    const S2Point& expected) {
  return ::testing::MakePolymorphicMatcher(
      S2PointMatcher(expected, S1Angle::Zero()));
}

::testing::PolymorphicMatcher<S2PointMatcher> PointApproxEquals(
    const S2Point& expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2PointMatcher(expected, tolerance));
}

// Matchers for S2Points provided as strings.

::testing::PolymorphicMatcher<S2PointMatcher> PointEquals(
    string_view expected) {
  return ::testing::MakePolymorphicMatcher(
      S2PointMatcher(expected, S1Angle::Zero()));
}

::testing::PolymorphicMatcher<S2PointMatcher> PointApproxEquals(
    string_view expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2PointMatcher(expected, tolerance));
}

// Matchers for S2Polylines provided as objects.

::testing::PolymorphicMatcher<S2PolylineMatcher> PolylineEquals(
    const S2Polyline& expected) {
  return ::testing::MakePolymorphicMatcher(S2PolylineMatcher(
      expected, S2PolylineMatcher::Options(S2PolylineMatcher::Options::EQUALS,
                                           S1Angle::Zero())));
}

::testing::PolymorphicMatcher<S2PolylineMatcher> PolylineApproxEquals(
    const S2Polyline& expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2PolylineMatcher(
      expected, S2PolylineMatcher::Options(
                    S2PolylineMatcher::Options::APPROX_EQUALS, tolerance)));
}

// Matchers for S2Polylines provided as strings.

::testing::PolymorphicMatcher<S2PolylineMatcher> PolylineEquals(
    string_view expected) {
  return ::testing::MakePolymorphicMatcher(S2PolylineMatcher(
      expected, S2PolylineMatcher::Options(S2PolylineMatcher::Options::EQUALS,
                                           S1Angle::Zero())));
}

::testing::PolymorphicMatcher<S2PolylineMatcher> PolylineApproxEquals(
    string_view expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2PolylineMatcher(
      expected, S2PolylineMatcher::Options(
                    S2PolylineMatcher::Options::APPROX_EQUALS, tolerance)));
}

// Matchers for S2Polygons provided as objects.

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonEquals(
    const S2Polygon& expected) {
  return ::testing::MakePolymorphicMatcher(S2PolygonMatcher(
      expected, S2PolygonMatcher::Options(S2PolygonMatcher::Options::EQUALS,
                                          S1Angle::Zero())));
}

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonApproxEquals(
    const S2Polygon& expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2PolygonMatcher(
      expected, S2PolygonMatcher::Options(
                    S2PolygonMatcher::Options::APPROX_EQUALS, tolerance)));
}

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonBoundaryEquals(
    const S2Polygon& expected) {
  return ::testing::MakePolymorphicMatcher(S2PolygonMatcher(
      expected,
      S2PolygonMatcher::Options(S2PolygonMatcher::Options::BOUNDARY_EQUALS,
                                S1Angle::Zero())));
}

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonBoundaryApproxEquals(
    const S2Polygon& expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2PolygonMatcher(
      expected,
      S2PolygonMatcher::Options(
          S2PolygonMatcher::Options::BOUNDARY_APPROX_EQUALS, tolerance)));
}

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonBoundaryNear(
    const S2Polygon& expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2PolygonMatcher(
      expected, S2PolygonMatcher::Options(
                    S2PolygonMatcher::Options::BOUNDARY_NEAR, tolerance)));
}

// Matchers for S2Polygons provided as string.

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonEquals(
    string_view expected) {
  return ::testing::MakePolymorphicMatcher(S2PolygonMatcher(
      expected, S2PolygonMatcher::Options(S2PolygonMatcher::Options::EQUALS,
                                          S1Angle::Zero())));
}

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonApproxEquals(
    string_view expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2PolygonMatcher(
      expected, S2PolygonMatcher::Options(
                    S2PolygonMatcher::Options::APPROX_EQUALS, tolerance)));
}

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonBoundaryEquals(
    string_view expected) {
  return ::testing::MakePolymorphicMatcher(S2PolygonMatcher(
      expected,
      S2PolygonMatcher::Options(S2PolygonMatcher::Options::BOUNDARY_EQUALS,
                                S1Angle::Zero())));
}

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonBoundaryApproxEquals(
    string_view expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2PolygonMatcher(
      expected,
      S2PolygonMatcher::Options(
          S2PolygonMatcher::Options::BOUNDARY_APPROX_EQUALS, tolerance)));
}

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonBoundaryNear(
    string_view expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2PolygonMatcher(
      expected, S2PolygonMatcher::Options(
                    S2PolygonMatcher::Options::BOUNDARY_NEAR, tolerance)));
}

// Matchers for S2Loops provided as objects.

::testing::PolymorphicMatcher<S2LoopMatcher> LoopEquals(
    const S2Loop& expected) {
  return ::testing::MakePolymorphicMatcher(S2LoopMatcher(
      expected,
      S2LoopMatcher::Options(S2LoopMatcher::Options::EQUALS, S1Angle::Zero())));
}

::testing::PolymorphicMatcher<S2LoopMatcher> LoopBoundaryEquals(
    const S2Loop& expected) {
  return ::testing::MakePolymorphicMatcher(S2LoopMatcher(
      expected, S2LoopMatcher::Options(S2LoopMatcher::Options::BOUNDARY_EQUALS,
                                       S1Angle::Zero())));
}

::testing::PolymorphicMatcher<S2LoopMatcher> LoopBoundaryApproxEquals(
    const S2Loop& expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2LoopMatcher(
      expected,
      S2LoopMatcher::Options(S2LoopMatcher::Options::BOUNDARY_APPROX_EQUALS,
                             tolerance)));
}

::testing::PolymorphicMatcher<S2LoopMatcher> LoopBoundaryNear(
    const S2Loop& expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2LoopMatcher(
      expected, S2LoopMatcher::Options(S2LoopMatcher::Options::BOUNDARY_NEAR,
                                       tolerance)));
}

// Matchers for S2Loops provided as string.

::testing::PolymorphicMatcher<S2LoopMatcher> LoopEquals(string_view expected) {
  return ::testing::MakePolymorphicMatcher(S2LoopMatcher(
      expected,
      S2LoopMatcher::Options(S2LoopMatcher::Options::EQUALS, S1Angle::Zero())));
}

::testing::PolymorphicMatcher<S2LoopMatcher> LoopBoundaryEquals(
    string_view expected) {
  return ::testing::MakePolymorphicMatcher(S2LoopMatcher(
      expected, S2LoopMatcher::Options(S2LoopMatcher::Options::BOUNDARY_EQUALS,
                                       S1Angle::Zero())));
}

::testing::PolymorphicMatcher<S2LoopMatcher> LoopBoundaryApproxEquals(
    string_view expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2LoopMatcher(
      expected,
      S2LoopMatcher::Options(S2LoopMatcher::Options::BOUNDARY_APPROX_EQUALS,
                             tolerance)));
}

::testing::PolymorphicMatcher<S2LoopMatcher> LoopBoundaryNear(
    string_view expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2LoopMatcher(
      expected, S2LoopMatcher::Options(S2LoopMatcher::Options::BOUNDARY_NEAR,
                                       tolerance)));
}

// Matchers for S2LatLngRect provided as objects.

::testing::PolymorphicMatcher<S2LatLngRectMatcher> LatLngRectEquals(
    const S2LatLngRect& expected) {
  return ::testing::MakePolymorphicMatcher(S2LatLngRectMatcher(
      expected, S2LatLngRectMatcher::Options(
                    S2LatLngRectMatcher::Options::EQUALS, S1Angle::Zero())));
}

::testing::PolymorphicMatcher<S2LatLngRectMatcher> LatLngRectApproxEquals(
    const S2LatLngRect& expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2LatLngRectMatcher(
      expected, S2LatLngRectMatcher::Options(
                    S2LatLngRectMatcher::Options::APPROX_EQUALS, tolerance)));
}

::testing::PolymorphicMatcher<S2LatLngRectMatcher> LatLngRectApproxEquals(
    const S2LatLngRect& expected, S2LatLng tolerance) {
  return ::testing::MakePolymorphicMatcher(S2LatLngRectMatcher(
      expected, S2LatLngRectMatcher::Options(
                    S2LatLngRectMatcher::Options::APPROX_EQUALS_2, tolerance)));
}

// Matchers for S2LatLngRect provided as string.

::testing::PolymorphicMatcher<S2LatLngRectMatcher> LatLngRectEquals(
    string_view expected) {
  return ::testing::MakePolymorphicMatcher(S2LatLngRectMatcher(
      expected, S2LatLngRectMatcher::Options(
                    S2LatLngRectMatcher::Options::EQUALS, S1Angle::Zero())));
}

::testing::PolymorphicMatcher<S2LatLngRectMatcher> LatLngRectApproxEquals(
    string_view expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2LatLngRectMatcher(
      expected, S2LatLngRectMatcher::Options(
                    S2LatLngRectMatcher::Options::APPROX_EQUALS, tolerance)));
}

::testing::PolymorphicMatcher<S2LatLngRectMatcher> LatLngRectApproxEquals(
    string_view expected, S2LatLng tolerance) {
  return ::testing::MakePolymorphicMatcher(S2LatLngRectMatcher(
      expected, S2LatLngRectMatcher::Options(
                    S2LatLngRectMatcher::Options::APPROX_EQUALS_2, tolerance)));
}

// Matchers for S2LatLng provided as objects.

::testing::PolymorphicMatcher<S2LatLngMatcher> LatLngEquals(
    const S2LatLng& expected) {
  return ::testing::MakePolymorphicMatcher(S2LatLngMatcher(
      expected, S2LatLngMatcher::Options(S2LatLngMatcher::Options::EQUALS,
                                         S1Angle::Zero())));
}

::testing::PolymorphicMatcher<S2LatLngMatcher> LatLngApproxEquals(
    const S2LatLng& expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2LatLngMatcher(
      expected, S2LatLngMatcher::Options(
                    S2LatLngMatcher::Options::APPROX_EQUALS, tolerance)));
}

// Matchers for S2LatLng provided as string.

::testing::PolymorphicMatcher<S2LatLngMatcher> LatLngEquals(
    string_view expected) {
  return ::testing::MakePolymorphicMatcher(S2LatLngMatcher(
      expected, S2LatLngMatcher::Options(S2LatLngMatcher::Options::EQUALS,
                                         S1Angle::Zero())));
}

::testing::PolymorphicMatcher<S2LatLngMatcher> LatLngApproxEquals(
    string_view expected, S1Angle tolerance) {
  return ::testing::MakePolymorphicMatcher(S2LatLngMatcher(
      expected, S2LatLngMatcher::Options(
                    S2LatLngMatcher::Options::APPROX_EQUALS, tolerance)));
}

}  // namespace S2
