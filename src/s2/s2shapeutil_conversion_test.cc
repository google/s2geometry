// Copyright 2021 Google Inc. All Rights Reserved.
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

#include "s2/s2shapeutil_conversion.h"

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include "absl/types/span.h"
#include "s2/s2lax_polygon_shape.h"
#include "s2/s2lax_polyline_shape.h"
#include "s2/s2loop.h"
#include "s2/s2point.h"
#include "s2/s2point_vector_shape.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"
#include "s2/s2text_format.h"

namespace s2shapeutil {
namespace {

using std::unique_ptr;
using std::vector;

// ShapeToS2Points tests

TEST(S2ShapeConversionUtilTest, PointVectorShapeToPoints) {
  vector<S2Point> points = s2textformat::ParsePointsOrDie("11:11, 10:0, 5:5");
  S2PointVectorShape point_vector(points);
  vector<S2Point> extract = ShapeToS2Points(point_vector);
  // TODO(user,b/205813109): Use gmock ASSERT_THAT.
  ASSERT_EQ(extract.size(), 3);
  for (int i = 0; i < extract.size(); i++) {
    EXPECT_EQ(extract[i], points[i]);
  }
}

// ShapeToS2Polyline tests

TEST(S2ShapeConversionUtilTest, LineToS2Polyline) {
  vector<S2Point> points = s2textformat::ParsePointsOrDie("11:11, 10:0, 5:5");
  S2LaxPolylineShape lax_polyline(points);
  unique_ptr<S2Polyline> polyline = ShapeToS2Polyline(lax_polyline);
  EXPECT_EQ(polyline->num_vertices(), 3);
  for (int i = 0; i < polyline->num_vertices(); i++) {
    EXPECT_EQ(polyline->vertex(i), points[i]);
  }
}

TEST(S2ShapeConversionUtilTest, ClosedLineToS2Polyline) {
  vector<S2Point> points =
      s2textformat::ParsePointsOrDie("0:0, 0:10, 10:10, 0:0");
  S2LaxPolylineShape lax_polyline(points);
  unique_ptr<S2Polyline> polyline = ShapeToS2Polyline(lax_polyline);
  EXPECT_EQ(polyline->num_vertices(), 4);
  for (int i = 0; i < polyline->num_vertices(); i++) {
    EXPECT_EQ(polyline->vertex(i), points[i]);
  }
}

// ShapeToS2Polygon tests

// Creates a (lax) polygon shape from the provided loops, and ensures that the
// S2Polygon produced by ShapeToS2Polygon represents the same polygon.
void VerifyShapeToS2Polygon(absl::Span<const S2LaxPolygonShape::Loop> loops,
                            int expected_num_loops, int expected_num_vertices) {
  S2LaxPolygonShape lax_polygon(loops);
  unique_ptr<S2Polygon> polygon = ShapeToS2Polygon(lax_polygon);

  EXPECT_EQ(polygon->num_loops(), expected_num_loops);
  EXPECT_EQ(polygon->num_vertices(), expected_num_vertices);
  for (int i = 0; i < polygon->num_loops(); i++) {
    const S2Loop& loop = *polygon->loop(i);
    for (int j = 0; j < loop.num_vertices(); j++) {
      EXPECT_EQ(loop.oriented_vertex(j), loops[i][j]);
    }
  }
}

TEST(S2ShapeConversionUtilTest, PolygonWithHoleToS2Polygon) {
  // a polygon with one shell and one hole
  S2LaxPolygonShape::Loop shell(
      s2textformat::ParsePointsOrDie("0:0, 0:10, 10:10, 10:0"));
  S2LaxPolygonShape::Loop hole(
      s2textformat::ParsePointsOrDie("4:4, 6:4, 6:6, 4:6"));
  vector<S2LaxPolygonShape::Loop> loops{shell, hole};

  VerifyShapeToS2Polygon(loops, 2, 8);
}

TEST(S2ShapeConversionUtilTest, MultiPolygonToS2Polygon) {
  // a polygon with multiple shells
  S2LaxPolygonShape::Loop shell1(
      s2textformat::ParsePointsOrDie("0:0, 0:2, 2:2, 2:0"));
  S2LaxPolygonShape::Loop shell2(
      s2textformat::ParsePointsOrDie("0:4, 0:6, 3:6"));
  vector<S2LaxPolygonShape::Loop> loops{shell1, shell2};

  VerifyShapeToS2Polygon(loops, 2, 7);
}

TEST(S2ShapeConversionUtilTest, TwoHolesToS2Polygon) {
  // a polygon shell with two holes
  S2LaxPolygonShape::Loop shell(s2textformat::ParsePointsOrDie(
      "0:0, 0:10, 10:10, 10:0"));
  S2LaxPolygonShape::Loop hole1(s2textformat::ParsePointsOrDie(
      "1:1, 3:3, 1:3"));
  S2LaxPolygonShape::Loop hole2(s2textformat::ParsePointsOrDie(
      "2:6, 4:7, 2:8"));
  vector<S2LaxPolygonShape::Loop> loops{shell, hole1, hole2};

  VerifyShapeToS2Polygon(loops, 3, 10);
}

TEST(S2ShapeConversionUtilTest, FullPolygonToS2Polygon) {
  // verify that a full polygon is converted correctly
  const absl::Span<const S2Point> full_loop = S2Loop::kFull();
  S2LaxPolygonShape::Loop loop1(full_loop.begin(), full_loop.end());
  vector<S2LaxPolygonShape::Loop> loops{loop1};

  auto lax_polygon = s2textformat::MakeLaxPolygonOrDie("full");
  unique_ptr<S2Polygon> polygon = ShapeToS2Polygon(*lax_polygon);
  EXPECT_EQ(polygon->num_loops(), 1);
  EXPECT_EQ(polygon->num_vertices(), 1);
  EXPECT_TRUE(polygon->is_full());
  for (int i = 0; i < polygon->num_loops(); i++) {
    const S2Loop& loop = *polygon->loop(i);
    for (int j = 0; j < loop.num_vertices(); j++) {
      EXPECT_EQ(loop.oriented_vertex(j), loops[i][j]);
    }
  }
}

}  // namespace
}  // namespace s2shapeutil
