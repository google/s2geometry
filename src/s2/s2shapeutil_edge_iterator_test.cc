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

#include "s2/s2shapeutil_edge_iterator.h"

#include <memory>
#include <string>
#include <vector>

#include <gtest/gtest.h>
#include "s2/mutable_s2shape_index.h"
#include "s2/s2shape.h"
#include "s2/s2shape_index.h"
#include "s2/s2shapeutil_shape_edge_id.h"
#include "s2/s2text_format.h"

namespace s2shapeutil {

namespace {

using std::vector;

// Returns the full list of edges in the given S2ShapeIndex.
// The edges are collected from points, lines, and polygons in that order.
vector<S2Shape::Edge> GetEdges(const S2ShapeIndex* index) {
  vector<S2Shape::Edge> result;
  for (const S2Shape* shape : *index) {
    if (shape == nullptr) continue;
    for (int j = 0; j < shape->num_edges(); ++j) {
      result.push_back(shape->edge(j));
    }
  }
  return result;
}

// Verifies that the edges produced by an EdgeIterator matches GetEdges.
void Verify(const S2ShapeIndex* index) {
  vector<S2Shape::Edge> expected = GetEdges(index);

  int i = 0;
  int shape_id = -1;
  int edge_id = -1;
  for (EdgeIterator it(index); !it.Done(); it.Next(), ++edge_id, ++i) {
    // The iterator visits the edges of each shape in order.  When we see a new
    // shape id, reset the edge_id count.
    if (it.shape_id() != shape_id) {
      shape_id = it.shape_id();
      edge_id = 0;
    }

    ASSERT_TRUE(i < expected.size());
    EXPECT_EQ(expected[i], it.edge());
    EXPECT_EQ(it.edge_id(), edge_id);
    EXPECT_EQ(it.shape_edge_id(), ShapeEdgeId(shape_id, edge_id));
  }
}

}  // namespace

TEST(S2ShapeutilEdgeIteratorTest, Empty) {
  auto index = s2textformat::MakeIndexOrDie("##");
  Verify(index.get());
}

TEST(S2ShapeutilEdgeIteratorTest, Points) {
  auto index = s2textformat::MakeIndexOrDie("0:0|1:1##");
  Verify(index.get());
}

TEST(S2ShapeutilEdgeIteratorTest, Lines) {
  auto index = s2textformat::MakeIndexOrDie("#0:0,10:10|5:5,5:10|1:2,2:1#");
  Verify(index.get());
}

TEST(S2ShapeutilEdgeIteratorTest, Polygons) {
  auto index =
      s2textformat::MakeIndexOrDie("##10:10,10:0,0:0|-10:-10,-10:0,0:0,0:-10");
  Verify(index.get());
}

TEST(S2ShapeutilEdgeIteratorTest, Collection) {
  auto index = s2textformat::MakeIndexOrDie(
      "1:1|7:2#1:1,2:2,3:3|2:2,1:7#"
      "10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");
  Verify(index.get());
}

TEST(S2ShapeutilEdgeIteratorTest, Remove) {
  auto index = s2textformat::MakeIndexOrDie(
      "1:1|7:2#1:1,2:2,3:3|2:2,1:7#"
      "10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");
  index->Release(0);

  Verify(index.get());
}

TEST(S2ShapeutilEdgeIteratorTest, AssignmentAndEquality) {
  auto index1 = s2textformat::MakeIndexOrDie(
      "1:1|7:2#1:1,2:2,3:3|2:2,1:7#"
      "10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");

  auto index2 = s2textformat::MakeIndexOrDie(
      "1:1|7:2#1:1,2:2,3:3|2:2,1:7#"
      "10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");

  EdgeIterator it1(index1.get());
  EdgeIterator it2(index2.get());

  // Different indices.
  EXPECT_TRUE(it1 != it2);

  it1 = it2;
  EXPECT_EQ(it1, it2);

  it1.Next();
  EXPECT_TRUE(it1 != it2);

  it2.Next();
  EXPECT_EQ(it1, it2);
}

}  // namespace s2shapeutil
