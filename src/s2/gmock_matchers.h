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

#ifndef S2_GMOCK_MATCHERS_H_
#define S2_GMOCK_MATCHERS_H_

#include <memory>
#include <ostream>
#include <string>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
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

// GMock matchers for S2LatLng, S2LatLngRect, S2Point, S2Polygon, S2Polyline,
// S2Loop.
//
// Defines polymorphic matchers that can accept a string representation
// of each object or a const reference to the object. See gmock_matchers_test.cc
// for additional usage examples.
//
// ====================
// =     S2LatLng     =
// ====================
//
// Matchers:
//   * LatLngEquals: defined as S2LatLng::operator==
//   * LatLngApproxEquals: defined as S2LatLng::ApproxEquals
//
// Examples:
//   EXPECT_THAT(latlng, S2::LatLngEquals("10:10"));
//   EXPECT_THAT(latlng, S2::LatLngApproxEquals(expected,
//                                              S1Angle::Degrees(0.1)));
//
// ========================
// =     S2LatLngRect     =
// ========================
//
// Matchers:
//   * LatLngRectEquals: defined as S2LatLngRect::operator==
//   * LatLngRectApproxEquals: defined as S2LatLngRect::ApproxEquals. Accepts
//         a single tolerance or individual tolerances for lat and lng.
//
// Examples:
//   EXPECT_THAT(rect, S2::LatLngRectEquals(expected));
//   EXPECT_THAT(rect, S2::LatLngRectApproxEquals(
//                         "0:0,10:10", S1Angle::Radians(1e-10)));
//   EXPECT_THAT(rect, S2::LatLngRectApproxEquals(
//                         expected, S2LatLng::FromDegrees(0.1, 0.2)));
//
// ===================
// =     S2Point     =
// ===================
//
// Matchers:
//   * PointEquals: exact equality check.
//   * PointApproxEquals: defined as s2pointutil::ApproxEquals with tolerance
//     with the exception that the default S2Point is not approximately equal
//     to anything. Expected value is not permitted to be zero or near zero.
//
// Examples:
//   EXPECT_THAT(point, S2::PointEquals("10:10");
//   EXPECT_THAT(point, S2::PointApproxEquals(
//                          expected, S1Angle::Radians(1e-10));
//
// =====================
// =     S2Polygon     =
// =====================
//
// Matchers:
//   * PolygonEquals: defined as S2Polygon::Equals
//   * PolygonApproxEquals: defined as S2Polygon::ApproxEquals
//   * PolygonBoundaryEquals: defined as S2Polygon::BoundaryEquals
//   * PolygonBoundaryApproxEquals: defined as S2Polygon::BoundaryApproxEquals
//   * PolygonBoundaryNear: defined as S2Polygon::BoundaryNear
//
// Examples:
//   EXPECT_THAT(polygon, S2::PolygonEquals("10:10,10:0,0:0");
//   EXPECT_THAT(polygon, S2::PolygonApproxEquals(
//                            expected, S1Angle::Radians(1e-10));
//   EXPECT_THAT(polygon, S2::PolygonBoundaryEquals(expected);
//   EXPECT_THAT(polygon, S2::PolygonBoundaryApproxEquals(
//                            expected, S1Angle::Radians(1e-10));
//   EXPECT_THAT(polygon, S2::PolygonBoundaryNear(
//                            expected, S1Angle::Radians(1e-10));
//
// ======================
// =     S2Polyline     =
// ======================
//
// Matchers:
//   * PolylineEquals: defined as S2Polyline::Equals
//   * PolylineApproxEquals: defined as S2Polyline::ApproxEquals
//
// Examples:
//   EXPECT_THAT(line, S2::PolylineEquals(expected));
//   EXPECT_THAT(line, S2::PolylineApproxEquals(
//                         "10:10,20:20", S1Angle::Radians(1e-10)));
//
// ==================
// =     S2Loop     =
// ==================
//
// Matchers:
//   * LoopEquals: defined as S2Loop::Equals
//   * LoopBoundaryEquals: defined as S2Loop::BoundaryEquals
//   * LoopBoundaryApproxEquals: defined as S2Loop::BoundaryApproxEquals
//   * LoopBoundaryNear: defined as S2Loop::BoundaryNear
//
// Examples:
//   EXPECT_THAT(loop, S2::LoopEquals("10:10,10:0,0:0"));
//   EXPECT_THAT(loop, S2::LoopBoundaryEquals("10:10,10:0,0:0"));
//   EXPECT_THAT(loop, S2::LoopBoundaryApproxEquals(
//                         expected, S1Angle::Radians(1e-10)));
//   EXPECT_THAT(loop, S2::LoopBoundaryNear(
//                         expected, S1Angle::Radians(1e-10)));

// Print methods necessary for gmock to print values of S2Point, S2Polyline,
// S2Polygon, and S2LatLngRect.

void PrintTo(const S2Point& p, std::ostream* os);
void PrintTo(const S2Polyline& l, std::ostream* os);
void PrintTo(const S2Polygon& p, std::ostream* os);
void PrintTo(const S2LatLngRect& r, std::ostream* os);
void PrintTo(const S2Loop& l, std::ostream* os);

namespace S2 {

// Matcher for S2Point.
class S2PointMatcher {
 public:
  S2PointMatcher(const S2Point& expected, S1Angle tolerance);

  S2PointMatcher(absl::string_view expected, S1Angle tolerance)
      : expected_(s2textformat::MakePointOrDie(expected)),
        tolerance_(tolerance) {}

  bool MatchAndExplain(absl::string_view s,
                       ::testing::MatchResultListener* listener) const;

  // When tolerance_ == S1Angle::Zero() the points will be matched with the
  // S2Point::operator== method.
  bool MatchAndExplain(const S2Point& g,
                       ::testing::MatchResultListener* listener) const;

  void DescribeTo(::std::ostream* os) const {
    *os << "== " << s2textformat::ToString(expected_);
  }

  void DescribeNegationTo(std::ostream* os) const {
    *os << "!= " << s2textformat::ToString(expected_);
  }

 private:
  S2Point expected_;
  S1Angle tolerance_;
};

// Matcher for S2Polyline.
class S2PolylineMatcher {
 public:
  struct Options {
    // The types of comparison correspond to methods in S2Polyline. Each type
    // uses the same name as the corresponding method in S2Polyline.
    enum ComparisonType {
      EQUALS,        // S2Polyline::Equals
      APPROX_EQUALS  // S2Polyline::ApproxEquals
    };

    Options() = default;

    Options(ComparisonType a_type, S1Angle a_tolerance)
        : type(a_type), tolerance(a_tolerance) {}

    ComparisonType type;
    S1Angle tolerance;
  };

  S2PolylineMatcher(const S2PolylineMatcher& other) {
    options_ = other.options_;
    expected_.reset(other.expected_->Clone());
  }

  S2PolylineMatcher& operator=(const S2PolylineMatcher& other) {
    options_ = other.options_;
    expected_.reset(other.expected_->Clone());
    return *this;
  }

  S2PolylineMatcher(const S2Polyline& expected, Options options)
      : expected_(expected.Clone()), options_(options) {}

  S2PolylineMatcher(absl::string_view expected, Options options)
      : expected_(s2textformat::MakePolylineOrDie(expected)),
        options_(options) {}

  virtual ~S2PolylineMatcher() = default;

  bool MatchAndExplain(absl::string_view s,
                       ::testing::MatchResultListener* listener) const;

  bool MatchAndExplain(const S2Polyline& g,
                       ::testing::MatchResultListener* listener) const;

  void DescribeTo(::std::ostream* os) const {
    *os << "== " << s2textformat::ToString(*expected_);
  }

  void DescribeNegationTo(std::ostream* os) const {
    *os << "!= " << s2textformat::ToString(*expected_);
  }

 private:
  std::unique_ptr<S2Polyline> expected_;
  Options options_;
};

// Matcher for S2Polygon.
class S2PolygonMatcher {
 public:
  struct Options {
    // The types of comparison correspond to methods in S2Polygon. Each type
    // uses the same name as the corresponding method in S2Polygon.
    enum ComparisonType {
      EQUALS,                  // S2Polygon::Equals
      APPROX_EQUALS,           // S2Polygon::ApproxEquals
      BOUNDARY_EQUALS,         // S2Polygon::BoundaryEquals
      BOUNDARY_APPROX_EQUALS,  // S2Polygon::BoundaryApproxEquals
      BOUNDARY_NEAR            // S2Polygon::BoundaryNear
    };

    Options() = default;

    Options(ComparisonType a_type, S1Angle a_tolerance)
        : type(a_type), tolerance(a_tolerance) {}

    ComparisonType type;
    S1Angle tolerance;
  };

  S2PolygonMatcher(const S2Polygon& expected, Options options)
      : expected_(expected.Clone()), options_(options) {}

  S2PolygonMatcher(absl::string_view expected, Options options)
      : expected_(s2textformat::MakePolygonOrDie(expected, S2Debug::ALLOW)),
        options_(options) {}

  S2PolygonMatcher(const S2PolygonMatcher& other) {
    options_ = other.options_;
    expected_.reset(other.expected_->Clone());
  }

  S2PolygonMatcher& operator=(const S2PolygonMatcher& other) {
    options_ = other.options_;
    expected_.reset(other.expected_->Clone());
    return *this;
  }

  bool MatchAndExplain(absl::string_view s,
                       ::testing::MatchResultListener* listener) const;

  bool MatchAndExplain(const S2Polygon& g,
                       ::testing::MatchResultListener* listener) const;

  void DescribeTo(::std::ostream* os) const {
    *os << "== " << s2textformat::ToString(*expected_);
  }

  void DescribeNegationTo(std::ostream* os) const {
    *os << "!= " << s2textformat::ToString(*expected_);
  }

 private:
  std::unique_ptr<S2Polygon> expected_;
  Options options_;
};

// Matcher for S2Loop.
class S2LoopMatcher {
 public:
  struct Options {
    // The types of comparison correspond to methods in S2Loop. Each type
    // uses the same name as the corresponding method in S2Loop.
    enum ComparisonType {
      EQUALS,                  // S2Loop::Equals
      BOUNDARY_EQUALS,         // S2Loop::BoundaryEquals
      BOUNDARY_APPROX_EQUALS,  // S2Loop::BoundaryApproxEquals
      BOUNDARY_NEAR            // S2Loop::BoundaryNear
    };

    Options() = default;

    Options(ComparisonType a_type, S1Angle a_tolerance)
        : type(a_type), tolerance(a_tolerance) {}

    ComparisonType type;
    S1Angle tolerance;
  };

  S2LoopMatcher(const S2Loop& expected, Options options)
      : expected_(expected.Clone()), options_(options) {}

  S2LoopMatcher(absl::string_view expected, Options options)
      : expected_(s2textformat::MakeLoopOrDie(expected)), options_(options) {}

  S2LoopMatcher(const S2LoopMatcher& other) {
    options_ = other.options_;
    expected_.reset(other.expected_->Clone());
  }

  S2LoopMatcher& operator=(const S2LoopMatcher& other) {
    options_ = other.options_;
    expected_.reset(other.expected_->Clone());
    return *this;
  }

  bool MatchAndExplain(absl::string_view s,
                       ::testing::MatchResultListener* listener) const;

  bool MatchAndExplain(const S2Loop& g,
                       ::testing::MatchResultListener* listener) const;

  void DescribeTo(::std::ostream* os) const {
    *os << "== " << s2textformat::ToString(*expected_);
  }

  void DescribeNegationTo(std::ostream* os) const {
    *os << "!= " << s2textformat::ToString(*expected_);
  }

 private:
  std::unique_ptr<S2Loop> expected_;
  Options options_;
};

// Matcher for S2LatLngRect.
class S2LatLngRectMatcher {
 public:
  struct Options {
    // The types of comparison correspond to methods in S2LatLngRect. Each type
    // uses the same name as the corresponding method in S2LatLngRect.
    enum ComparisonType {
      EQUALS,           // S2LatLngRect::operator==
      APPROX_EQUALS,    // S2LatLngRect::ApproxEquals.
                        // Same tolerance for lat and lng.
      APPROX_EQUALS_2,  // s2LatLngRect::ApproxEquals.
                        // Separate tolerance for lat and lng.
    };

    Options() = default;

    Options(ComparisonType a_type, S1Angle a_tolerance)
        : type(a_type), tolerance(a_tolerance) {}

    Options(ComparisonType a_type, S2LatLng tolerance)
        : type(a_type), tolerance_lat_lng(tolerance) {}

    ComparisonType type;
    S1Angle tolerance;
    S2LatLng tolerance_lat_lng;
  };

  S2LatLngRectMatcher(const S2LatLngRect& expected, Options options)
      : expected_(expected), options_(options) {}

  S2LatLngRectMatcher(absl::string_view expected, Options options)
      : expected_(s2textformat::MakeLatLngRectOrDie(expected)),
        options_(options) {}

  S2LatLngRectMatcher(const S2LatLngRectMatcher& other) {
    options_ = other.options_;
    expected_ = other.expected_;
  }

  S2LatLngRectMatcher& operator=(const S2LatLngRectMatcher& other) = default;

  bool MatchAndExplain(absl::string_view s,
                       ::testing::MatchResultListener* listener) const;

  bool MatchAndExplain(const S2LatLngRect& g,
                       ::testing::MatchResultListener* listener) const;

  void DescribeTo(::std::ostream* os) const {
    *os << "== " << s2textformat::ToString(expected_);
  }

  void DescribeNegationTo(std::ostream* os) const {
    *os << "!= " << s2textformat::ToString(expected_);
  }

 private:
  S2LatLngRect expected_;
  Options options_;
};

// Matcher for S2LatLng.
class S2LatLngMatcher {
 public:
  struct Options {
    // The types of comparison correspond to methods in S2LatLng. Each type
    // uses the same name as the corresponding method in S2LatLng.
    enum ComparisonType {
      EQUALS,         // S2LatLng::operator==
      APPROX_EQUALS,  // S2LatLng::ApproxEquals.
                      // Same tolerance for lat and lng.
    };
    Options() = default;

    Options(ComparisonType a_type, S1Angle a_tolerance)
        : type(a_type), tolerance(a_tolerance) {}

    ComparisonType type;
    S1Angle tolerance;
  };

  S2LatLngMatcher(const S2LatLng& expected, Options options)
      : expected_(expected), options_(options) {}

  S2LatLngMatcher(absl::string_view expected, Options options)
      : expected_(s2textformat::MakeLatLngOrDie(expected)), options_(options) {}

  S2LatLngMatcher(const S2LatLngMatcher& other) {
    options_ = other.options_;
    expected_ = other.expected_;
  }

  S2LatLngMatcher& operator=(const S2LatLngMatcher& other) = default;

  bool MatchAndExplain(absl::string_view s,
                       ::testing::MatchResultListener* listener) const;

  bool MatchAndExplain(const S2LatLng& g,
                       ::testing::MatchResultListener* listener) const;

  void DescribeTo(::std::ostream* os) const {
    *os << "== " << s2textformat::ToString(expected_);
  }

  void DescribeNegationTo(std::ostream* os) const {
    *os << "!= " << s2textformat::ToString(expected_);
  }

 private:
  S2LatLng expected_;
  Options options_;
};

// Predefined matchers for S2Point, S2Polyline, S2Polygon, S2Loop, S2LatLngRect.
// These matchers can be used in EXPECT_THAT statements.

// Matchers for S2Points provided as objects.

::testing::PolymorphicMatcher<S2PointMatcher> PointEquals(
    const S2Point& expected);

::testing::PolymorphicMatcher<S2PointMatcher> PointApproxEquals(
    const S2Point& expected, S1Angle tolerance);

// Matchers for S2Points provided as strings.

::testing::PolymorphicMatcher<S2PointMatcher> PointEquals(
    absl::string_view expected);

::testing::PolymorphicMatcher<S2PointMatcher> PointApproxEquals(
    absl::string_view expected, S1Angle tolerance);

// Matchers for S2Polylines provided as objects.

::testing::PolymorphicMatcher<S2PolylineMatcher> PolylineEquals(
    const S2Polyline& expected);

::testing::PolymorphicMatcher<S2PolylineMatcher> PolylineApproxEquals(
    const S2Polyline& expected, S1Angle tolerance);

// Matchers for S2Polylines provided as strings.

::testing::PolymorphicMatcher<S2PolylineMatcher> PolylineEquals(
    absl::string_view expected);

::testing::PolymorphicMatcher<S2PolylineMatcher> PolylineApproxEquals(
    absl::string_view expected, S1Angle tolerance);

// Matchers for S2Polygons provided as objects.

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonEquals(
    const S2Polygon& expected);

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonApproxEquals(
    const S2Polygon& expected, S1Angle tolerance);

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonBoundaryEquals(
    const S2Polygon& expected);

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonBoundaryApproxEquals(
    const S2Polygon& expected, S1Angle tolerance);

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonBoundaryNear(
    const S2Polygon& expected, S1Angle tolerance);

// Matchers for S2Polygons provided as string.

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonEquals(
    absl::string_view expected);

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonApproxEquals(
    absl::string_view expected, S1Angle tolerance);

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonBoundaryEquals(
    absl::string_view expected);

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonBoundaryApproxEquals(
    absl::string_view expected, S1Angle tolerance);

::testing::PolymorphicMatcher<S2PolygonMatcher> PolygonBoundaryNear(
    absl::string_view expected, S1Angle tolerance);

// Matchers for S2Loops provided as objects.

::testing::PolymorphicMatcher<S2LoopMatcher> LoopEquals(const S2Loop& expected);

::testing::PolymorphicMatcher<S2LoopMatcher> LoopBoundaryEquals(
    const S2Loop& expected);

::testing::PolymorphicMatcher<S2LoopMatcher> LoopBoundaryApproxEquals(
    const S2Loop& expected, S1Angle tolerance);

::testing::PolymorphicMatcher<S2LoopMatcher> LoopBoundaryNear(
    const S2Loop& expected, S1Angle tolerance);

// Matchers for S2Loops provided as string.

::testing::PolymorphicMatcher<S2LoopMatcher> LoopEquals(
    absl::string_view expected);

::testing::PolymorphicMatcher<S2LoopMatcher> LoopBoundaryEquals(
    absl::string_view expected);

::testing::PolymorphicMatcher<S2LoopMatcher> LoopBoundaryApproxEquals(
    absl::string_view expected, S1Angle tolerance);

::testing::PolymorphicMatcher<S2LoopMatcher> LoopBoundaryNear(
    absl::string_view expected, S1Angle tolerance);

// Matchers for S2LatLngRect provided as objects.

::testing::PolymorphicMatcher<S2LatLngRectMatcher> LatLngRectEquals(
    const S2LatLngRect& expected);

::testing::PolymorphicMatcher<S2LatLngRectMatcher> LatLngRectApproxEquals(
    const S2LatLngRect& expected, S1Angle tolerance);

::testing::PolymorphicMatcher<S2LatLngRectMatcher> LatLngRectApproxEquals(
    const S2LatLngRect& expected, S2LatLng tolerance);

// Matchers for S2LatLngRect provided as string.

::testing::PolymorphicMatcher<S2LatLngRectMatcher> LatLngRectEquals(
    absl::string_view expected);

::testing::PolymorphicMatcher<S2LatLngRectMatcher> LatLngRectApproxEquals(
    absl::string_view expected, S1Angle tolerance);

::testing::PolymorphicMatcher<S2LatLngRectMatcher> LatLngRectApproxEquals(
    absl::string_view expected, S2LatLng tolerance);

// Matchers for S2LatLng provided as objects.

::testing::PolymorphicMatcher<S2LatLngMatcher> LatLngEquals(
    const S2LatLng& expected);

::testing::PolymorphicMatcher<S2LatLngMatcher> LatLngApproxEquals(
    const S2LatLng& expected, S1Angle tolerance = S1Angle::Radians(1e-15));

// Matchers for S2LatLng provided as string.

::testing::PolymorphicMatcher<S2LatLngMatcher> LatLngEquals(
    absl::string_view expected);

::testing::PolymorphicMatcher<S2LatLngMatcher> LatLngApproxEquals(
    absl::string_view expected, S1Angle tolerance = S1Angle::Radians(1e-15));

}  // namespace S2

#endif  // S2_GMOCK_MATCHERS_H_
