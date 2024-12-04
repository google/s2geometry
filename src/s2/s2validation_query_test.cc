// Copyright 2022 Google Inc. All Rights Reserved.
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


#include "s2/s2validation_query.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/algorithm/container.h"
#include "absl/flags/flag.h"
#include "absl/log/absl_check.h"
#include "absl/log/log_streamer.h"
#include "absl/random/bit_gen_ref.h"
#include "absl/random/random.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "s2/util/coding/coder.h"
#include "s2/mutable_s2shape_index.h"
#include "s2/s2debug.h"
#include "s2/s2error.h"
#include "s2/s2fractal.h"
#include "s2/s2latlng.h"
#include "s2/s2lax_polygon_shape.h"
#include "s2/s2point.h"
#include "s2/s2random.h"
#include "s2/s2shape.h"
#include "s2/s2shapeutil_coding.h"
#include "s2/s2testing.h"
#include "s2/s2text_format.h"

namespace {

using ::absl::StrCat;
using ::absl::string_view;
using ::std::make_unique;
using ::std::string;
using ::std::unique_ptr;
using ::std::vector;
using ::testing::Contains;
using ::testing::Eq;
using ::testing::IsFalse;
using ::testing::IsTrue;

// Converts latitude/longitude in radians to an S2Point
S2Point LatLngToPoint(double lat, double lng) {
  return S2LatLng::FromRadians(lat, lng).ToPoint();
}

// Diamond shaped pattern of points around (0,0).
const S2Point kDiamondPoints[] = {S2LatLng::FromDegrees(+1.0, +0.0).ToPoint(),
                                  S2LatLng::FromDegrees(-1.0, +0.0).ToPoint(),
                                  S2LatLng::FromDegrees(+0.0, -1.0).ToPoint(),
                                  S2LatLng::FromDegrees(+0.0, +1.0).ToPoint()};

// An S2Shape subclass with one chain that's open.
class OpenShape : public S2Shape {
 public:
  int dimension() const final { return 2; }

  int num_edges() const final { return 3; }
  Edge chain_edge(int, int id) const final { return edge(id); }
  Edge edge(int id) const final {
    ABSL_CHECK_GE(id, 0);
    ABSL_CHECK_LT(id, 3);

    switch (id) {
      case 0:
        return Edge(kDiamondPoints[0], kDiamondPoints[1]);
      case 1:
        return Edge(kDiamondPoints[1], kDiamondPoints[2]);
      case 2:
        return Edge(kDiamondPoints[2], kDiamondPoints[3]);
      default:
        ABSL_CHECK(false);
    }
  }

  int num_chains() const final { return 1; }
  Chain chain(int) const final { return {1, num_edges()}; }

  ReferencePoint GetReferencePoint() const final {
    return {S2Point(1, 0, 0), true};
  }
  ChainPosition chain_position(int id) const final { return {0, id}; }
  TypeTag type_tag() const final { return 1; }
};

// An S2Shape subclass returning an invalid dimension value.
class BadDimensionShape : public S2Shape {
 public:
  int dimension() const final { return 42; }
  int num_edges() const final { return 0; }
  Edge edge(int) const final { return {}; }
  Edge chain_edge(int, int) const final { return {}; }
  Chain chain(int) const final { return {}; }
  int num_chains() const final { return 1; }
  ChainPosition chain_position(int) const final { return {}; }
  ReferencePoint GetReferencePoint() const final {
    return {S2Point(1, 0, 0), true};
  }
  TypeTag type_tag() const final { return 1; }
};

// An S2Shape subclass returning a chain that's too short.
class BadChainLengthShape : public S2Shape {
 public:
  int dimension() const final { return 2; }
  int num_edges() const final { return 2; }

  Edge edge(int id) const final {
    ABSL_CHECK_GE(id, 0);
    ABSL_CHECK_LT(id, 2);

    switch (id) {
      case 0:
        return Edge(kDiamondPoints[0], kDiamondPoints[1]);
      case 1:
        return Edge(kDiamondPoints[1], kDiamondPoints[2]);
      default:
        ABSL_CHECK(false);
    }
  }

  Edge chain_edge(int, int) const final { return {}; }
  int num_chains() const final { return 1; }
  Chain chain(int) const final { return {0, num_edges()}; }
  ChainPosition chain_position(int) const final { return {}; }
  ReferencePoint GetReferencePoint() const final {
    return {S2Point(1, 0, 0), true};
  }
  TypeTag type_tag() const final { return 1; }
};

// Returns num evenly spaced edges all sharing a common center point.
vector<S2Shape::Edge> CcwEdgesAbout(S2Point center, int num = 10) {
  vector<S2Shape::Edge> edges;
  for (int i = 0; i < num; ++i) {
    double angle = 2 * M_PI / num * i;
    edges.emplace_back(center, LatLngToPoint(std::sin(angle), std::cos(angle)));
  }
  return edges;
}

// Builds a "quilt" test shape which is a series of rings that stretch from the
// south to north pole and form a grid where every vertex has at least two
// chains incident on it.
std::unique_ptr<S2LaxPolygonShape> MakeQuilt() {
  // Defines a grid where points are separated by 15 degrees in longitude and 30
  // degrees in latitude.  The integral indices for x thus range from [0, 24)
  // and y ranges from [0, 13].  X values are automatically wrapped into range.
  const auto GridPoint = [](int x, int y) {
    constexpr double kLatDelta = 15;  // degrees;
    constexpr double kLonDelta = 15;  // degrees;

    ABSL_DCHECK(0 <= y && y <= 12);
    ABSL_DCHECK_LE(0, x);
    x %= 24;

    // Return exact polar points to avoid numerical issues.
    if (y == 0) return S2Point(0, 0, -1);
    if (y == 12) return S2Point(0, 0, +1);

    return S2LatLng::FromDegrees(-90 + kLatDelta * y, -180 + kLonDelta * x)
        .ToPoint();
  };

  // Now build loops that touch at every-other point.  This will give us a
  // diamond quilt pattern where every vertex has two chains incident on it.
  std::vector<std::vector<S2Point>> loops;

  for (int x = 0; x < 24; x += 2) {
    for (int y = 0; y < 12; y += 2) {
      std::vector<S2Point> loop;
      loop.emplace_back(GridPoint(x + 0, y + 1));
      loop.emplace_back(GridPoint(x + 1, y + 2));
      loop.emplace_back(GridPoint(x + 2, y + 1));
      loop.emplace_back(GridPoint(x + 1, y + 0));
      loops.emplace_back(std::move(loop));
    }
  }
  return std::make_unique<S2LaxPolygonShape>(std::move(loops));
}

TEST(SortEdgesCcw, SortsEdges) {
  const S2Point kOrigin = LatLngToPoint(0.0, 0.0);
  constexpr int kNumEdges = 10;

  // Generate edges in order around the origin.
  vector<S2Shape::Edge> sorted = CcwEdgesAbout(kOrigin, kNumEdges);

  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "SORT_EDGES_CCW_SORTS_EDGES",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < kNumEdges; ++i) {
    // Shift sorted edges by one.
    absl::c_rotate(sorted, sorted.begin() + 1);

    // Randomize the edges.
    vector<S2Shape::Edge> shuffled = sorted;
    absl::c_shuffle(shuffled, bitgen);

    // Sorting about the first edge should always give back the sorted array.
    SortEdgesCcw(kOrigin, sorted[0], shuffled);
    EXPECT_THAT(shuffled, Eq(sorted));
  }
}

TEST(SortEdgesCcw, SortsEdgesFlipped) {
  const S2Point kOrigin = LatLngToPoint(0.0, 0.0);
  constexpr int kNumEdges = 10;

  // Generate edges in order around the origin.
  vector<S2Shape::Edge> sorted = CcwEdgesAbout(kOrigin, kNumEdges);

  // Flip the orientation of some of the edges.  This only changes their
  // direction, not their ordering around the vertex.
  sorted[3] = sorted[3].Reversed();
  sorted[8] = sorted[8].Reversed();

  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "SORT_EDGES_CCW_SORTS_EDGES_FLIPPED",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < kNumEdges; ++i) {
    // Shift sorted edges by one.
    absl::c_rotate(sorted, sorted.begin() + 1);

    // randomize the edges.
    vector<S2Shape::Edge> shuffled = sorted;
    absl::c_shuffle(shuffled, bitgen);

    // Sorting about the first edge should always give back the sorted array.
    SortEdgesCcw(kOrigin, sorted[0], shuffled);
    EXPECT_THAT(shuffled, Eq(sorted));
  }
}

TEST(SortEdgesCcw, StartEdgeAlwaysFirst) {
  const S2Point kOrigin = LatLngToPoint(0.0, 0.0);
  constexpr int kNumEdges = 10;

  // Generate edges in order around the origin.
  vector<S2Shape::Edge> sorted = CcwEdgesAbout(kOrigin, kNumEdges);

  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "SORT_EDGES_CCW_START_EDGE_ALWAYS_FIRST",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < kNumEdges; ++i) {
    // Randomize the edges.
    vector<S2Shape::Edge> shuffled = sorted;
    absl::c_shuffle(shuffled, bitgen);

    // The reference edge we sort around should always come out first.
    SortEdgesCcw(kOrigin, sorted[i], shuffled);
    EXPECT_THAT(shuffled[0], Eq(sorted[i]));
  }
}

TEST(SortEdgesCcw, ReverseDuplicatesOrdered) {
  const S2Point kOrigin = LatLngToPoint(0.0, 0.0);
  constexpr int kNumEdges = 10;

  // Generate edges in order around the origin.
  vector<S2Shape::Edge> sorted = CcwEdgesAbout(kOrigin, kNumEdges);

  // Insert two reverse duplicate edges.
  sorted.insert(sorted.begin() + 8, sorted[8].Reversed());
  sorted.insert(sorted.begin() + 3, sorted[3].Reversed());

  // Randomize the edges.
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "SORT_EDGES_CCW_REVERSE_DUPLICATES_ORDERED",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  vector<S2Shape::Edge> shuffled = sorted;
  absl::c_shuffle(shuffled, bitgen);

  // After sorting, the reverse duplicates should be right after their sibling.
  SortEdgesCcw(kOrigin, sorted[4], shuffled);

  S2Point common = sorted[4].v0;
  EXPECT_THAT(shuffled[0], Eq(shuffled[1].Reversed()));
  EXPECT_THAT(shuffled[0].v0, Eq(common));
  EXPECT_THAT(shuffled[6], Eq(shuffled[7].Reversed()));
  EXPECT_THAT(shuffled[6].v0, Eq(common));
}

// Common logic to test across various validation queries.
template <typename QueryType>
class ValidationQueryTest : public ::testing::Test {
  using IndexPtr = unique_ptr<MutableS2ShapeIndex>;

 protected:
  // Checks that a particular shape is valid.
  void ExpectShapeValid(unique_ptr<S2Shape> shape) {
    MutableS2ShapeIndex index;
    index.Add(std::move(shape));

    S2Error error;
    EXPECT_THAT(query_.Validate(index, &error), IsTrue());
    EXPECT_THAT(error.code(), Eq(S2Error::OK));
  }

  // Checks that a particular S2ShapeIndex fails to validate, and that the error
  // matches one of the given codes.
  void ExpectIndexInvalid(const MutableS2ShapeIndex& index,
                          absl::Span<const S2Error::Code> codes = {}) {
    S2Error error;
    EXPECT_THAT(query_.Validate(index, &error), IsFalse());
    if (!codes.empty()) {
      EXPECT_THAT(codes, Contains(error.code()));
    }
  }

  // Checks that a particular shape fails to validate, and that the error
  // matches the given code.
  void ExpectShapeInvalid(unique_ptr<S2Shape> shape,
                          absl::Span<const S2Error::Code> codes = {}) {
    MutableS2ShapeIndex index;
    index.Add(std::move(shape));
    ExpectIndexInvalid(index, codes);
  }

  void ExpectShapeInvalid(unique_ptr<S2Shape> shape, S2Error::Code code) {
    ExpectShapeInvalid(std::move(shape), std::vector<S2Error::Code>{code});
  }

  // Checks that the given geometry string validates.
  void ExpectGeometryValid(string_view geometry) {
    IndexPtr index = s2textformat::MakeIndexOrDie(geometry);
    S2Error error;
    EXPECT_THAT(query_.Validate(*index, &error), IsTrue()) << error;
    EXPECT_THAT(error.code(), Eq(S2Error::OK));
  }

  // Checks that the given geometry string does not validate, and that the given
  // error code is returned.
  void ExpectGeometryInvalid(string_view geometry, S2Error::Code code) {
    IndexPtr index = s2textformat::MakeIndexOrDie(geometry);
    S2Error error;
    EXPECT_THAT(query_.Validate(*index, &error), IsFalse());
    EXPECT_THAT(error.code(), Eq(code));
  }

  // Create "num_loops" nested regular loops around a common center point.  All
  // loops have the same number of vertices (at least "min_vertices").
  // Furthermore, the vertices at the same index position are collinear with the
  // common center point of all the loops.  The loop radii decrease
  // exponentially in order to prevent accidental loop crossings when one of the
  // loops is modified.
  void AddConcentricLoops(  //
      absl::BitGenRef bitgen, int num_loops, int min_vertices) {
    ABSL_DCHECK_LE(num_loops, 10);  // Because radii decrease exponentially.
    S2Point center = s2random::Point(bitgen);
    int num_vertices = min_vertices + absl::Uniform(bitgen, 0, 10);
    for (int i = 0; i < num_loops; ++i) {
      S1Angle radius = S1Angle::Degrees(80 * pow(0.1, i));
      loops_.emplace_back(  //
          S2Testing::MakeRegularPoints(center, radius, num_vertices));
    }
  }

  // Creates and returns a new polygon from the given loops.
  S2Polygon MakePolygon(absl::BitGenRef bitgen) {
    std::vector<std::unique_ptr<S2Loop>> loops;
    for (const auto& loop : loops_) {
      loops.push_back(make_unique<S2Loop>(loop, S2Debug::DISABLE));
    }

    absl::c_shuffle(loops, bitgen);
    S2Polygon polygon;
    polygon.set_s2debug_override(S2Debug::DISABLE);
    polygon.InitNested(std::move(loops));
    return polygon;
  }

  std::vector<std::vector<S2Point>> loops_;

 private:
  QueryType query_;
};

// Define test suite for all validation queries.
template <typename IndexType>
using AllValidationQueries = ValidationQueryTest<IndexType>;

using AllQueryTypes = ::testing::Types<  //
    S2ValidQuery<S2ShapeIndex>,          //
    S2LegacyValidQuery<S2ShapeIndex>     //
  >;
TYPED_TEST_SUITE(AllValidationQueries, AllQueryTypes);

TYPED_TEST(AllValidationQueries, BasicGeometryOk) {
  // Basic polygon should be valid.
  this->ExpectGeometryValid("## 1:0, 0:-1, -1:0, 0:1");

  // Basic polyline should be valid.
  this->ExpectGeometryValid("# 0:0, 1:0, 0:-1, -1:0, 0:1 #");

  // Basic multipoint should be valid.
  this->ExpectGeometryValid("0:0 | 1:0 | 0:-1 | -1:0 | 0:1 ##");

  // Basic polygon with a hole should be valid.
  this->ExpectGeometryValid(      //
      "## 2:0, 0:-2, -2:0, 0:2;"  //
      "   0:1, -1:0, 0:-1, 1:0;"  //
  );

  // Polygon with improperly oriented hole should fail.
  this->ExpectGeometryInvalid(                         //
      "## 2:0, 0:-2, -2:0, 0:2;"                       //
      "   1:0, 0:-1, -1:0, 0:1;",                      //
      S2Error::POLYGON_INCONSISTENT_LOOP_ORIENTATIONS  //
  );
}

TYPED_TEST(AllValidationQueries, EmptyGeometryOk) {
  this->ExpectGeometryValid("##");
}

TYPED_TEST(AllValidationQueries, FullGeometryOk) {
  this->ExpectGeometryValid("## full");
}

TYPED_TEST(AllValidationQueries, InteriorOnRightRegression) {
  // Found via fuzzing, this should not test as invalid.  The root cause was
  // that we weren't clearing the incident edge list before gathering edges,
  // causing a duplicate edge to appear and give us the wrong interior state.
  this->ExpectGeometryValid(  //
      "## 0:4, 3:128, 4:2, 0:0");
}

TYPED_TEST(AllValidationQueries, TangentPolygonsOk) {
  // Two polygons touching at one vertex should be valid.
  this->ExpectGeometryValid(     //
      "## 1:0, 0:-1, -1:0, 0:1"  //
      "|  0:1, -1:2,  0:3, 1:2"  //
  );
}

TYPED_TEST(AllValidationQueries, AntipodalEdgeFails) {
  // We need to specify the points instead of using MakeIndexOrDie to ensure
  // they're exactly opposite sign to be antipodal.
  const vector<vector<S2Point>> kPoints = {{S2Point(M_SQRT1_2, M_SQRT1_2, 0),
                                            S2Point(0, 1, 0), S2Point(-1, 0, 0),
                                            S2Point(1, 0, 0)}};

  this->ExpectShapeInvalid(make_unique<S2LaxPolygonShape>(kPoints),
                           S2Error::ANTIPODAL_VERTICES);
}

TYPED_TEST(AllValidationQueries, BadlyDimensionedFails) {
  this->ExpectShapeInvalid(make_unique<BadDimensionShape>(),
                           S2Error::INVALID_DIMENSION);
}

TYPED_TEST(AllValidationQueries, OpenChainFails) {
  this->ExpectShapeInvalid(make_unique<OpenShape>(),
                           S2Error::LOOP_NOT_ENOUGH_VERTICES);
}

TYPED_TEST(AllValidationQueries, DuplicatePolygonEdgesFail) {
  // Two polygons sharing an edge is invalid.
  this->ExpectGeometryInvalid(   //
      "## 2:0, 0:-2, -2:0, 0:2"  //
      " | 2:0, 0:-2,  0:0",      //
      S2Error::OVERLAPPING_GEOMETRY);
}

TYPED_TEST(AllValidationQueries, ChainsTouchingOk) {
  // Polygon with a hole with chains touching at a point should be valid.
  this->ExpectGeometryValid(      //
      "## 2:0, 0:-2, -2:0, 0:2;"  //
      "   0:2, -1:0, 0:-1, 1:0;"  //
  );

  this->ExpectGeometryValid(      //
      "## 2:0, 0:-2, -2:0, 0:2;"  //
      "   0:1, -2:0, 0:-1, 1:0;"  //
  );

  this->ExpectGeometryInvalid(                         //
      "## 2:0,  0:-2, -2:0, 0:2;"                      //
      "   1:0,  0:-2, -1:0, 0:2;",                     //
      S2Error::POLYGON_INCONSISTENT_LOOP_ORIENTATIONS  //
  );
}

TYPED_TEST(AllValidationQueries, NestedShellsFail) {
  // Polygon with improperly oriented hole should fail
  this->ExpectGeometryInvalid(                         //
      "## 2:0, 0:-2, -2:0, 0:2;"                       //
      "   1:0, 0:-1, -1:0, 0:1",                       //
      S2Error::POLYGON_INCONSISTENT_LOOP_ORIENTATIONS  //
  );

  // Even if they're touching.
  this->ExpectGeometryInvalid(                         //
      "## 2:0, 0:-2, -2:0, 0:2;"                       //
      "   2:0, 0:-1, -1:0, 0:1",                       //
      S2Error::POLYGON_INCONSISTENT_LOOP_ORIENTATIONS  //
  );

  this->ExpectGeometryInvalid(                         //
      "## 2:0, 0:-2, -2:0, 0:2;"                       //
      "   2:0, 0:-1, -2:0, 0:1",                       //
      S2Error::POLYGON_INCONSISTENT_LOOP_ORIENTATIONS  //
  );

  this->ExpectGeometryInvalid(                         //
      "## 2:0, 0:-2, -2:0, 0:2;"                       //
      "   1:0, 0:-2, -1:0, 0:1",                       //
      S2Error::POLYGON_INCONSISTENT_LOOP_ORIENTATIONS  //
  );

  this->ExpectGeometryInvalid(                         //
      "## 2:0, 0:-2, -2:0, 0:2;"                       //
      "   1:0, 0:-1, -2:0, 0:1",                       //
      S2Error::POLYGON_INCONSISTENT_LOOP_ORIENTATIONS  //
  );

  this->ExpectGeometryInvalid(                         //
      "## 2:0, 0:-2, -2:0, 0:2;"                       //
      "   1:0, 0:-1, -1:0, 0:2",                       //
      S2Error::POLYGON_INCONSISTENT_LOOP_ORIENTATIONS  //
  );
}

TYPED_TEST(AllValidationQueries, ChainsCannotCross) {
  // Polygon chains crossing is an error, regardless of orientation.
  this->ExpectGeometryInvalid(     //
      "## 3:0, 0:-3, -3:0, 0:+3;"  //
      "   3:2, 0:-1, -3:2, 0:+5",  //
      S2Error::POLYGON_INCONSISTENT_LOOP_ORIENTATIONS);

  this->ExpectGeometryInvalid(      //
      "## 0:3, 3:0,   0:-3, -3:0;"  //
      "   3:2, 0:+5, -3:2,  0:-1",  //
      S2Error::OVERLAPPING_GEOMETRY);

  // Crossing at a vertex isn't allowed either.  This test vector is carefully
  // constructed to avoid triggering the interior-on-the-right check.
  this->ExpectGeometryInvalid(                        //
      "## 0:-6, -6:0, 0:6, 6:0 ;"                     //
      "   0:0,   3:0, 6:0, 6:3, 6:6, 3:6, 0:6, 0:3",  //
      S2Error::OVERLAPPING_GEOMETRY);
}

TYPED_TEST(AllValidationQueries, ShellInHoleFails) {
  this->ExpectGeometryInvalid(  //
      "## 0:0, 10:10, 10:0; 5:21, 8:21, 6:23",
      S2Error::POLYGON_INCONSISTENT_LOOP_ORIENTATIONS);
}

TYPED_TEST(AllValidationQueries, LoopsCrossing) {
  // None of the validation semantics allow loops to cross.
  static constexpr int kIters = 100;

  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "LOOPS_CROSSING",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  for (int iter = 0; iter < kIters; ++iter) {
    this->AddConcentricLoops(bitgen, 2, 4 /*min_vertices*/);

    // Grab the two loops that were just added.
    std::vector<S2Point>& loop0 = this->loops_[this->loops_.size() - 2];
    std::vector<S2Point>& loop1 = this->loops_[this->loops_.size() - 1];

    // Both loops have the same number of vertices, and vertices at the same
    // index position are collinear with the center point, so we can create a
    // crossing by simply exchanging two vertices at the same index position.
    int n = loop0.size();
    int i = absl::Uniform(bitgen, 0, n);
    std::swap(loop0[i], loop1[i]);
    if (absl::Bernoulli(bitgen, 0.5)) {
      // By copy the two adjacent vertices from one loop to the other, we can
      // ensure that the crossings happen at vertices rather than edges.
      loop0[(i + 1) % n] = loop1[(i + 1) % n];
      loop0[(i + n - 1) % n] = loop1[(i + n - 1) % n];
    }

    S2Polygon polygon = this->MakePolygon(bitgen);

    // Don't check the error code, just that the polygon is invalid.  There's
    // several different types of invalidity that can result from crossing loops
    // this way.
    this->ExpectShapeInvalid(std::make_unique<S2Polygon::Shape>(&polygon));
  }
}

TYPED_TEST(AllValidationQueries, IndexWithUnindexVerticesFails) {
  // This was found by fuzz testing.  It contains a chain vertex that's not
  // actually indexed, so checking chain orientations will fail, which we should
  // catch.
  const std::string data(
      "\010\212\005\001\002Y\000\010\177\210[\226\001\020\002\003\004\215\026)"
      "!\373\362\362\341\257?\013\212t\250\343\304\000\210\363\213~\032:\306?"
      "\027\034\201\214\213\203\357?\000\000\000p\341\000\000\000\211s\013~"
      "\032:\306?[j\260}\330d\357?\320\246\303os$\306?"
      "\002\201\302\270\000O\266?\242\025\030\264\323\340\357\002\201\302\270("
      "\326O\266?\000\000\000\000\000\000$\000\000\027?\211s\013]]]~\032:\306?"
      "\000\000\000\000\000\000\000\000\014\000\003\006("
      "\340\010\020\010\001\020",
      147);

  S2Error error;
  Decoder decoder(data.data(), data.size());

  MutableS2ShapeIndex index;
  EXPECT_TRUE(index.Init(&decoder,
                         s2shapeutil::FullDecodeShapeFactory(&decoder, error)));
  EXPECT_TRUE(error.ok()) << error.message();
  this->ExpectIndexInvalid(index, {S2Error::DATA_LOSS});
}

TYPED_TEST(AllValidationQueries, OutgoingEdgeButNoIncomingEdge) {
  // This was found by fuzz testing.  It contains a vertex with an outgoing edge
  // but no corresponding incoming edge.
  const std::string data(
      "\010\212\005\001\002Y\000\010\177\020\000\001\020\000M\004\215\026)!"
      "\373\010\357?\364\013\212t\250\343\305?\211s\013~\032:\306?"
      "\027\034\201\214\213\203\357?\205\000\000\000\000\000\000\000\211s\013~"
      "\032:\306?[j\260}\330d\357?\320\246\303os$\306?"
      "\002\201\302\270\326O\266?\242\025\030\264\323\340\357?"
      "\002\201\302\270\326O\266?\000\000\000\000\003 "
      "\000\000\027\034\201\214\213\203\357?\211s\013~\032:\306?"
      "\377\000\363\000\336\000\000\000\014\000\003\006("
      "\340\010\020\010\001\014\020\264",
      149);

  S2Error error;
  Decoder decoder(data.data(), data.size());

  MutableS2ShapeIndex index;
  EXPECT_TRUE(index.Init(&decoder,
                         s2shapeutil::FullDecodeShapeFactory(&decoder, error)));
  EXPECT_TRUE(error.ok()) << error.message();
  this->ExpectIndexInvalid(index, {S2Error::INVALID_VERTEX});
}

TYPED_TEST(AllValidationQueries, InvalidChainNearChain) {
  // This was found by fuzz testing.  It contains a cell with a chain that has
  // valid points but also a chain with unnormalized points, which would fail
  // if we didn't defer chain orientation checks.
  const std::string data(
      "\010\212\005\001\002Y\000\010\177\020\000\001\020\002\003\004\215\026)!"
      "\373\010\357?\364\013\212t\250\343\305?\211s\013~\032:\306?"
      "\027\034\201\214\213\203\357?\000\000\000\000\000\000\000\000\211s\013~"
      "\032:\306?[j\260~\330d\357?\320\246\303os$\306?"
      "\002\201\302\270\326O\266?\242\025\030\264\323\340\357?"
      "\002\201\302\270\326O\266?"
      "\000\000\000\000\000\000\000\000\027\034\201\214\213\203\357?\211s\013~"
      "\032:\306?\000\000\000\000\000\000\000\000\014\000\003\006("
      "\340\010\020\010\001\020",
      147);

  S2Error error;
  Decoder decoder(data.data(), data.size());

  MutableS2ShapeIndex index;
  EXPECT_TRUE(index.Init(&decoder,
                         s2shapeutil::FullDecodeShapeFactory(&decoder, error)));
  EXPECT_TRUE(error.ok()) << error.message();
  this->ExpectIndexInvalid(index, {S2Error::NOT_UNIT_LENGTH});
}

//////////////////   Multidimensional Query Tests   ////////////////////

// Define test suite for validation queries that allow multidimensional geometry
template <typename IndexType>
using MultiDimensionalQueries = ValidationQueryTest<IndexType>;

using MultiDimensionalTypes = ::testing::Types<  //
    S2ValidQuery<MutableS2ShapeIndex>            //
    >;

TYPED_TEST_SUITE(MultiDimensionalQueries, MultiDimensionalTypes);

TYPED_TEST(MultiDimensionalQueries, BasicGeometryOk) {
  // Basic multi-dimensional geometry should be valid;
  this->ExpectGeometryValid(    //
      "  3:0| 0:-3| -3:0| 0:3"  //
      "# 2:0, 0:-2, -2:0, 0:2"  //
      "# 1:0, 0:-1, -1:0, 0:1"  //
  );
}

TYPED_TEST(MultiDimensionalQueries, ContainedGeometryFails) {
  // Point contained in polygon is invalid.
  this->ExpectGeometryInvalid("0:0 ## 2:0, 0:-2, -2:0, 0:2",
                              S2Error::OVERLAPPING_GEOMETRY);

  // Polyline contained in polygon is invalid.
  this->ExpectGeometryInvalid("# 0:-1, 0:1 # 2:0, 0:-2, -2:0, 0:2",
                              S2Error::OVERLAPPING_GEOMETRY);

  // Polygon contained in polygon is invalid.
  this->ExpectGeometryInvalid(       //
      "## 2:0, 0:-2, -2:0, 0:2"      //
      " | 1:0, 0:-1, -1:0, 0:1",     //
      S2Error::OVERLAPPING_GEOMETRY  //
  );

  // Either end of a polyline contained in a polygon is invalid.
  this->ExpectGeometryInvalid("# 0:-3, 0:1 # 2:0, 0:-2, -2:0, 0:2",
                              S2Error::OVERLAPPING_GEOMETRY);
  this->ExpectGeometryInvalid("# 0:-1, 0:3 # 2:0, 0:-2, -2:0, 0:2",
                              S2Error::OVERLAPPING_GEOMETRY);

  // Polyline edges crossing are fine though.
  this->ExpectGeometryValid("# 0:-1, 0:1 | 1:0, -1:0 #");
}

//////////////////   S2ValidQuery Tests   ////////////////////

// Define test suite for s2valid validation query only.
template <typename IndexType>
using S2ValidTest = ValidationQueryTest<IndexType>;

TYPED_TEST_SUITE(S2ValidTest,
                 ::testing::Types<S2ValidQuery<MutableS2ShapeIndex>>);

TYPED_TEST(S2ValidTest, QuiltIsValid) { this->ExpectShapeValid(MakeQuilt()); }

TYPED_TEST(S2ValidTest, DegenerateRingsAllowed) {
  // Point loops with a single edge should be allowed.
  this->ExpectGeometryValid("## 0:0");

  // Degenerate loops with a sibling edge pair should be allowed.
  this->ExpectGeometryValid("## 0:0, 1:1");
}

TYPED_TEST(S2ValidTest, SplitInteriorsOk) {
  // Chains touching at two points (splitting interior), should be fine.
  this->ExpectGeometryValid(       //
      "## 3:0, 0:-3, -3:0, 0:+3;"  //
      "   3:0, 0:+1, -3:0, 0:-1"   //
  );
}

TYPED_TEST(S2ValidTest, PolylineEdgesCrossSemanticsOk) {
  // Interior crossings between polylines are ok.
  this->ExpectGeometryValid(         //
      "# 0:0, 1:1, 0:2, 1:3, 0:4 "   //
      "| 1:0, 0:1, 1:2, 0:3, 1:4 #"  //
  );

  // Vertex crossings between polylines are ok.
  this->ExpectGeometryValid(                             //
      "# 0:0, 1:1, 2:2, 1:3, 0:4, 1:5, 2:6, 1:7, 0:8"    //
      "| 2:0, 1:1, 0:2, 1:3, 2:4, 1:5, 0:6, 1:7, 2:8 #"  //
  );

  // Interior crossings within a polyline are ok.
  this->ExpectGeometryValid(                                  //
      "# 0:0, 1:1, 0:2, 1:3, 0:4, 1:4, 0:3, 1:2, 0:1, 1:0 #"  //
  );

  // Vertex crossings within a polyline are ok.
  this->ExpectGeometryValid(                             //
      "# 0:0, 1:1, 2:2, 1:3, 0:4, 1:5, 2:6, 1:7, 0:8,"   //
      "  2:0, 1:1, 0:2, 1:3, 2:4, 1:5, 0:6, 1:7, 2:8 #"  //
  );

  // A closed loop that touches at the endpoints is ok.
  this->ExpectGeometryValid(         //
      "# 2:1, 1:0, 0:1, 1:2, 2:1 #"  //
  );

  // Two polylines touching at endpoints is also OK;
  this->ExpectGeometryValid(  //
      "# 0:0, 1:1, 0:2"       //
      "| 1:3, 0:4, 1:5 #"     //
  );
}

TYPED_TEST(S2ValidTest, ReverseDuplicateOnCenterWorks) {
  // A pair of reverse duplicate edges touching the cell center should work.
  this->ExpectGeometryValid(      //
      "## 2:0, 0:-2, -2:0, 0:2;"  //
      "   0:0, 1:1");
}

TYPED_TEST(S2ValidTest, PolygonOnCentersWorks) {
  // Draw a nested-diamond polygon using 8 cells straddling the equator/prime
  // meridian.  The chain orientation check will have to fall back to a brute
  // force check but these should still work.
  const vector<vector<S2Point>> kPoints = {
      {
          S2Cell(S2CellId::FromToken("0ec")).GetCenter(),
          S2Cell(S2CellId::FromToken("044")).GetCenter(),
          S2Cell(S2CellId::FromToken("1bc")).GetCenter(),
          S2Cell(S2CellId::FromToken("114")).GetCenter(),
      },
      {
          S2Cell(S2CellId::FromToken("104")).GetCenter(),
          S2Cell(S2CellId::FromToken("1ac")).GetCenter(),
          S2Cell(S2CellId::FromToken("054")).GetCenter(),
          S2Cell(S2CellId::FromToken("0fc")).GetCenter(),
      }};

  this->ExpectShapeValid(make_unique<S2LaxPolygonShape>(kPoints));
}

TYPED_TEST(S2ValidTest, DegeneratePolygonOnCentersworks) {
  // Now to be really rude and test a polygon consisting of nothing but
  // reverse-duplicate pairs between cell centers.
  const vector<vector<S2Point>> kPoints = {{
      S2Cell(S2CellId::FromToken("0ec")).GetCenter(),
      S2Cell(S2CellId::FromToken("044")).GetCenter(),
      S2Cell(S2CellId::FromToken("1bc")).GetCenter(),
      S2Cell(S2CellId::FromToken("114")).GetCenter(),
      S2Cell(S2CellId::FromToken("1bc")).GetCenter(),
      S2Cell(S2CellId::FromToken("044")).GetCenter(),
  }};

  this->ExpectShapeValid(make_unique<S2LaxPolygonShape>(kPoints));

  // Same thing but the polygon just goes out at a diagonal and back.
  const std::string tokens[] = {"1004", "1014", "1044", "1054", "1104", "1114"};

  vector<vector<S2Point>> loops;
  vector<S2Point> loop;
  for (int i = 0; i < 6; ++i) {
    loop.push_back(S2Cell(S2CellId::FromToken(tokens[i])).GetCenter());
  }

  for (int i = 4; i > 0; --i) {
    loop.push_back(S2Cell(S2CellId::FromToken(tokens[i])).GetCenter());
  }
  loops.push_back(std::move(loop));

  this->ExpectShapeValid(make_unique<S2LaxPolygonShape>(loops));
}

//////////////////   S2LegacyValid Tests   ////////////////////

// Define test suite for legacy validation query only.
template <typename IndexType>
using S2LegacyValidTest = ValidationQueryTest<IndexType>;

TYPED_TEST_SUITE(S2LegacyValidTest,                                         //
                 ::testing::Types<S2LegacyValidQuery<MutableS2ShapeIndex>>  //
);

TYPED_TEST(S2LegacyValidTest, QuiltIsNotValid) {
  // The quilt has reverse duplicate edges near the poles.
  this->ExpectShapeInvalid(MakeQuilt(), S2Error::OVERLAPPING_GEOMETRY);
}

TYPED_TEST(S2LegacyValidTest, MultiDimensionalFails) {
  // Multi-dimensional geometry shouldn't be allowed.
  this->ExpectGeometryInvalid(    //
      "  3:0| 0:-3| -3:0| 0:3"    //
      "# 2:0, 0:-2, -2:0, 0:2"    //
      "# 1:0, 0:-1, -1:0, 0:1",   //
      S2Error::INVALID_DIMENSION  //
  );
}

TYPED_TEST(S2LegacyValidTest, SplitInteriorsOk) {
  // Chains touching at two points (splitting interior), should be fine.
  this->ExpectGeometryValid(       //
      "## 3:0, 0:-3, -3:0, 0:+3;"  //
      "   3:0, 0:+1, -3:0, 0:-1"   //
  );
}

TYPED_TEST(S2LegacyValidTest, SelfTouchingLoopFails) {
  // A loop that touches itself should fail.
  this->ExpectGeometryInvalid("## 2:0, 0:-2, -2:0, -1:1, 0:-2, 1:1",
                              S2Error::DUPLICATE_VERTICES);
}

TYPED_TEST(S2LegacyValidTest, DegenerateEdgesFail) {
  // A degenerate polygon edge should fail.
  this->ExpectGeometryInvalid("## 2:0, 2:0, 0:-2, -2:0, 0:-2",
                              S2Error::DUPLICATE_VERTICES);

  // A degenerate polyline edge should fail.
  this->ExpectGeometryInvalid("# 0:0, 0:0, 1:1, 2:2 #",
                              S2Error::DUPLICATE_VERTICES);
}

TYPED_TEST(S2LegacyValidTest, ShortChainsFail) {
  const S2Error::Code kCode = S2Error::LOOP_NOT_ENOUGH_VERTICES;

  // A polygon chain with one or two edges should fail.
  this->ExpectGeometryInvalid("## 0:0", kCode);
  this->ExpectGeometryInvalid("## 0:0, 1:1", kCode);
}

}  // namespace
