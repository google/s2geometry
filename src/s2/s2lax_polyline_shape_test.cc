// Copyright 2013 Google Inc. All Rights Reserved.
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

#include "s2/s2lax_polyline_shape.h"

#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>
#include "s2/util/coding/coder.h"
#include "s2/s2coder.h"
#include "s2/s2coder_testing.h"
#include "s2/s2error.h"
#include "s2/s2point.h"
#include "s2/s2shape.h"
#include "s2/s2shapeutil_testing.h"
#include "s2/s2text_format.h"

using std::vector;

TEST(S2LaxPolylineShape, NoVertices) {
  vector<S2Point> vertices;
  S2LaxPolylineShape shape(vertices);
  EXPECT_EQ(0, shape.num_edges());
  EXPECT_EQ(0, shape.num_chains());
  EXPECT_EQ(1, shape.dimension());
  EXPECT_TRUE(shape.is_empty());
  EXPECT_FALSE(shape.is_full());
  EXPECT_FALSE(shape.GetReferencePoint().contained);
}

TEST(S2LaxPolylineShape, OneVertex) {
  vector<S2Point> vertices = {S2Point(1, 0, 0)};
  S2LaxPolylineShape shape(vertices);
  EXPECT_EQ(0, shape.num_edges());
  EXPECT_EQ(0, shape.num_chains());
  EXPECT_EQ(1, shape.dimension());
  EXPECT_TRUE(shape.is_empty());
  EXPECT_FALSE(shape.is_full());
}

TEST(S2LaxPolylineShape, Move) {
  // Construct a shape to use as the correct answer and a second identical shape
  // to be moved.
  const std::vector<S2Point> vertices =
      s2textformat::ParsePointsOrDie("1:1, 4:4, 2:2, 3:3");
  const S2LaxPolylineShape correct(vertices);
  S2LaxPolylineShape to_move(vertices);

  // Test the move constructor.
  S2LaxPolylineShape move1(std::move(to_move));
  s2testing::ExpectEqual(correct, move1);
  EXPECT_EQ(correct.id(), move1.id());
  ASSERT_EQ(vertices.size(), move1.num_vertices());
  for (int i = 0; i < move1.num_vertices(); ++i) {
    ASSERT_EQ(vertices[i], move1.vertex(i));
  }

  // Test the move-assignment operator.
  S2LaxPolylineShape move2;
  move2 = std::move(move1);
  s2testing::ExpectEqual(correct, move2);
  EXPECT_EQ(correct.id(), move2.id());
  ASSERT_EQ(vertices.size(), move2.num_vertices());
  for (int i = 0; i < move2.num_vertices(); ++i) {
    ASSERT_EQ(vertices[i], move2.vertex(i));
  }
}

TEST(S2LaxPolylineShape, EdgeAccess) {
  vector<S2Point> vertices = s2textformat::ParsePointsOrDie("0:0, 0:1, 1:1");
  S2LaxPolylineShape shape(vertices);
  EXPECT_EQ(2, shape.num_edges());
  EXPECT_EQ(1, shape.num_chains());
  EXPECT_EQ(0, shape.chain(0).start);
  EXPECT_EQ(2, shape.chain(0).length);
  EXPECT_EQ(1, shape.dimension());
  EXPECT_FALSE(shape.is_empty());
  EXPECT_FALSE(shape.is_full());
  auto edge0 = shape.edge(0);
  EXPECT_EQ(vertices[0], edge0.v0);
  EXPECT_EQ(vertices[1], edge0.v1);
  auto edge1 = shape.edge(1);
  EXPECT_EQ(vertices[1], edge1.v0);
  EXPECT_EQ(vertices[2], edge1.v1);
}

TEST(EncodedS2LaxPolylineShape, RoundtripEncoding) {
  vector<S2Point> vertices = s2textformat::ParsePointsOrDie("0:0, 0:1, 1:1");
  S2LaxPolylineShape shape(vertices);

  Encoder encoder;
  shape.Encode(&encoder, s2coding::CodingHint::COMPACT);
  Decoder a_decoder(encoder.base(), encoder.length());
  EncodedS2LaxPolylineShape a_shape;
  ASSERT_TRUE(a_shape.Init(&a_decoder));

  Encoder b_encoder;
  a_shape.Encode(&b_encoder, s2coding::CodingHint::COMPACT);
  Decoder b_decoder(b_encoder.base(), b_encoder.length());
  EncodedS2LaxPolylineShape b_shape;
  ASSERT_TRUE(b_shape.Init(&b_decoder));
  s2testing::ExpectEqual(shape, b_shape);
}

// TODO(b/222446546): Decoding EncodedS2PointVector on ARM isn't currently
// supported, so comment out S2Coder test on ARM for now.
#ifndef __arm__

TEST(S2LaxPolylineShape, S2CoderWorks) {
  vector<S2Point> vertices = s2textformat::ParsePointsOrDie("0:0, 0:1, 1:1");
  S2LaxPolylineShape shape(vertices);

  S2Error error;
  auto decoded = s2coding::RoundTrip(S2LaxPolylineShape::Coder(), shape, error);
  s2testing::ExpectEqual(decoded, shape);
}

TEST(S2LaxPolylineShape, ChainIteratorWorks) {
  S2LaxPolylineShape empty;
  std::vector<S2Point> points = s2textformat::ParsePointsOrDie("0:0, 0:1, 1:1");
  S2LaxPolylineShape shape(points);

  S2Shape::ChainIterator empty_begin = empty.chains().begin();
  S2Shape::ChainIterator empty_end = empty.chains().end();
  S2Shape::ChainIterator it1 = shape.chains().begin();
  S2Shape::ChainIterator it2 = shape.chains().begin();
  S2Shape::ChainIterator end = shape.chains().end();

  for (const auto& chain : shape.chains()) {
    EXPECT_EQ(chain.start, 0);
    EXPECT_EQ(chain.length, 2);
  }

  EXPECT_EQ(empty_begin, empty_end);
  EXPECT_EQ(it1, it2);
  EXPECT_NE(it1, end);
  EXPECT_EQ((*it1).start, 0);
  EXPECT_EQ((*it1).length, 2);
  EXPECT_EQ(++it1, end);
  EXPECT_NE(it1, it2);
  EXPECT_NE(it2++, end);
  EXPECT_EQ(it2, end);
}

TEST(S2LaxPolylineShape, ChainVertexIteratorWorks) {
  std::vector<std::vector<S2Point>> test_sets;
  test_sets.push_back(s2textformat::ParsePointsOrDie("0:0, 0:0"));
  test_sets.push_back(s2textformat::ParsePointsOrDie("0:0, 0:1"));
  test_sets.push_back(s2textformat::ParsePointsOrDie("0:0, 0:1, 1:1"));
  test_sets.push_back(s2textformat::ParsePointsOrDie("0:0, 0:1, 1:1, 2:2"));

  for (int i = 0; i < test_sets.size(); ++i) {
    const auto& points = test_sets[i];
    S2LaxPolylineShape shape(points);

    S2Shape::Chain chain = *shape.chains().begin();
    S2Shape::ChainVertexRange vertices(&shape, chain);
    EXPECT_EQ(vertices.num_vertices(), points.size());

    auto it1 = vertices.begin();
    auto it2 = it1;
    int vertex_index = 0;
    for (const S2Point& p : shape.vertices(0)) {
      EXPECT_EQ(p, points[vertex_index]);
      EXPECT_EQ(p, *S2Shape::ChainVertexIterator(&shape, chain, vertex_index));

      EXPECT_NE(it1, vertices.end());
      EXPECT_NE(it2, vertices.end());
      EXPECT_NE(it1++, vertices.end());
      ++it2;

      ++vertex_index;
    }
    EXPECT_EQ(it1, vertices.end());
    EXPECT_EQ(it2, vertices.end());
  }
}

#endif
