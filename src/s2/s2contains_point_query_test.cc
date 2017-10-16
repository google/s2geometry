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

// Author: ericv@google.com (Eric Veach)

#include "s2/s2contains_point_query.h"

#include <memory>
#include <gtest/gtest.h>
#include "s2/third_party/absl/memory/memory.h"
#include "s2/s2cap.h"
#include "s2/s2loop.h"
#include "s2/s2shapeindex.h"
#include "s2/s2testing.h"
#include "s2/s2textformat.h"

using absl::make_unique;
using s2shapeutil::ShapeEdge;
using s2shapeutil::ShapeEdgeId;
using s2textformat::MakeIndex;
using s2textformat::MakePoint;
using std::vector;

TEST(S2ContainsPointQuery, VertexModelOpen) {
  auto index = MakeIndex("0:0 # 0:1, 0:2 # 0:5, 0:7, 2:6");
  S2ContainsPointQueryOptions options(S2VertexModel::OPEN);
  auto q = MakeS2ContainsPointQuery(index.get(), options);
  EXPECT_FALSE(q.Contains(MakePoint("0:0")));
  EXPECT_FALSE(q.Contains(MakePoint("0:1")));
  EXPECT_FALSE(q.Contains(MakePoint("0:2")));
  EXPECT_FALSE(q.Contains(MakePoint("0:5")));
  EXPECT_FALSE(q.Contains(MakePoint("0:7")));
  EXPECT_FALSE(q.Contains(MakePoint("2:6")));
  EXPECT_TRUE(q.Contains(MakePoint("1:6")));
  EXPECT_FALSE(q.Contains(MakePoint("10:10")));

  // Test the last few cases using the Init() method instead.
  S2ContainsPointQuery<S2ShapeIndex> q2;
  q2.Init(index.get(), options);
  EXPECT_FALSE(q2.ShapeContains(*index->shape(1), MakePoint("1:6")));
  EXPECT_TRUE(q2.ShapeContains(*index->shape(2), MakePoint("1:6")));
  EXPECT_FALSE(q2.ShapeContains(*index->shape(2), MakePoint("0:5")));
  EXPECT_FALSE(q2.ShapeContains(*index->shape(2), MakePoint("0:7")));
}

TEST(S2ContainsPointQuery, VertexModelSemiOpen) {
  auto index = MakeIndex("0:0 # 0:1, 0:2 # 0:5, 0:7, 2:6");
  S2ContainsPointQueryOptions options(S2VertexModel::SEMI_OPEN);
  auto q = MakeS2ContainsPointQuery(index.get(), options);
  EXPECT_FALSE(q.Contains(MakePoint("0:0")));
  EXPECT_FALSE(q.Contains(MakePoint("0:1")));
  EXPECT_FALSE(q.Contains(MakePoint("0:2")));
  EXPECT_FALSE(q.Contains(MakePoint("0:5")));
  EXPECT_TRUE(q.Contains(MakePoint("0:7")));  // Contained vertex.
  EXPECT_FALSE(q.Contains(MakePoint("2:6")));
  EXPECT_TRUE(q.Contains(MakePoint("1:6")));
  EXPECT_FALSE(q.Contains(MakePoint("10:10")));

  // Test the last few cases using the Init() method instead.
  S2ContainsPointQuery<S2ShapeIndex> q2;
  q2.Init(index.get(), options);
  EXPECT_FALSE(q2.ShapeContains(*index->shape(1), MakePoint("1:6")));
  EXPECT_TRUE(q2.ShapeContains(*index->shape(2), MakePoint("1:6")));
  EXPECT_FALSE(q2.ShapeContains(*index->shape(2), MakePoint("0:5")));
  EXPECT_TRUE(q2.ShapeContains(*index->shape(2), MakePoint("0:7")));
}

TEST(S2ContainsPointQuery, VertexModelClosed) {
  auto index = MakeIndex("0:0 # 0:1, 0:2 # 0:5, 0:7, 2:6");
  S2ContainsPointQueryOptions options(S2VertexModel::CLOSED);
  auto q = MakeS2ContainsPointQuery(index.get(), options);
  EXPECT_TRUE(q.Contains(MakePoint("0:0")));
  EXPECT_TRUE(q.Contains(MakePoint("0:1")));
  EXPECT_TRUE(q.Contains(MakePoint("0:2")));
  EXPECT_TRUE(q.Contains(MakePoint("0:5")));
  EXPECT_TRUE(q.Contains(MakePoint("0:7")));
  EXPECT_TRUE(q.Contains(MakePoint("2:6")));
  EXPECT_TRUE(q.Contains(MakePoint("1:6")));
  EXPECT_FALSE(q.Contains(MakePoint("10:10")));

  // Test the last few cases using the Init() method instead.
  S2ContainsPointQuery<S2ShapeIndex> q2;
  q2.Init(index.get(), options);
  EXPECT_FALSE(q2.ShapeContains(*index->shape(1), MakePoint("1:6")));
  EXPECT_TRUE(q2.ShapeContains(*index->shape(2), MakePoint("1:6")));
  EXPECT_TRUE(q2.ShapeContains(*index->shape(2), MakePoint("0:5")));
  EXPECT_TRUE(q2.ShapeContains(*index->shape(2), MakePoint("0:7")));
}

TEST(S2ContainsPointQuery, GetContainingShapes) {
  // Also tests ShapeContains().
  int const kNumVerticesPerLoop = 10;
  S1Angle const kMaxLoopRadius = S2Testing::KmToAngle(10);
  S2Cap const center_cap(S2Testing::RandomPoint(), kMaxLoopRadius);
  S2ShapeIndex index;
  for (int i = 0; i < 100; ++i) {
    std::unique_ptr<S2Loop> loop = S2Loop::MakeRegularLoop(
        S2Testing::SamplePoint(center_cap),
        S2Testing::rnd.RandDouble() * kMaxLoopRadius, kNumVerticesPerLoop);
    index.Add(make_unique<S2Loop::OwningShape>(std::move(loop)));
  }
  auto query = MakeS2ContainsPointQuery(&index);
  for (int i = 0; i < 100; ++i) {
    S2Point p = S2Testing::SamplePoint(center_cap);
    vector<S2Shape*> expected;
    for (int j = 0; j < index.num_shape_ids(); ++j) {
      S2Shape* shape = index.shape(j);
      S2Loop const* loop = down_cast<S2Loop::Shape const*>(shape)->loop();
      if (loop->Contains(p)) {
        EXPECT_TRUE(query.ShapeContains(*shape, p));
        expected.push_back(shape);
      } else {
        EXPECT_FALSE(query.ShapeContains(*shape, p));
      }
    }
    vector<S2Shape*> actual = query.GetContainingShapes(p);
    EXPECT_EQ(expected, actual);
  }
}

using EdgeIdVector = vector<ShapeEdgeId>;

void ExpectIncidentEdgeIds(EdgeIdVector const& expected,
                           S2ShapeIndex const& index, S2Point const& p) {
  EdgeIdVector actual;
  auto q = MakeS2ContainsPointQuery(&index);
  EXPECT_TRUE(
      q.VisitIncidentEdges(p, [&actual](s2shapeutil::ShapeEdge const& e) {
          actual.push_back(e.id());
          return true;
        }));
  EXPECT_EQ(expected, actual);
}

TEST(S2ContainsPointQuery, VisitIncidentEdges) {
  auto index = MakeIndex("0:0 | 1:1 # 1:1, 1:2 # 1:2, 1:3, 2:2");
  ExpectIncidentEdgeIds({{0, 0}}, *index, MakePoint("0:0"));
  ExpectIncidentEdgeIds({{0, 1}, {1, 0}}, *index, MakePoint("1:1"));
  ExpectIncidentEdgeIds({{1, 0}, {2, 0}, {2, 2}}, *index, MakePoint("1:2"));
  ExpectIncidentEdgeIds({{2, 0}, {2, 1}}, *index, MakePoint("1:3"));
  ExpectIncidentEdgeIds({{2, 1}, {2, 2}}, *index, MakePoint("2:2"));
}
