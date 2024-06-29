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
#include <string>
#include <utility>
#include <vector>

#include "s2/base/casts.h"
#include <gtest/gtest.h>
#include "absl/log/log_streamer.h"
#include "absl/random/random.h"
#include "s2/mutable_s2shape_index.h"
#include "s2/s1angle.h"
#include "s2/s2cap.h"
#include "s2/s2loop.h"
#include "s2/s2point.h"
#include "s2/s2random.h"
#include "s2/s2shape.h"
#include "s2/s2shapeutil_shape_edge.h"
#include "s2/s2shapeutil_shape_edge_id.h"
#include "s2/s2testing.h"
#include "s2/s2text_format.h"

using s2shapeutil::ShapeEdgeId;
using s2textformat::MakeIndexOrDie;
using s2textformat::MakePointOrDie;
using std::make_unique;
using std::unique_ptr;
using std::vector;

TEST(S2ContainsPointQuery, VertexModelOpen) {
  auto index = MakeIndexOrDie("0:0 # -1:1, 1:1 # 0:5, 0:7, 2:6");
  S2ContainsPointQueryOptions options(S2VertexModel::OPEN);
  auto q = MakeS2ContainsPointQuery(index.get(), options);
  EXPECT_FALSE(q.Contains(MakePointOrDie("0:0")));
  EXPECT_FALSE(q.Contains(MakePointOrDie("-1:1")));
  EXPECT_FALSE(q.Contains(MakePointOrDie("1:1")));
  EXPECT_FALSE(q.Contains(MakePointOrDie("0:2")));
  EXPECT_FALSE(q.Contains(MakePointOrDie("0:3")));
  EXPECT_FALSE(q.Contains(MakePointOrDie("0:5")));
  EXPECT_FALSE(q.Contains(MakePointOrDie("0:7")));
  EXPECT_FALSE(q.Contains(MakePointOrDie("2:6")));
  EXPECT_TRUE(q.Contains(MakePointOrDie("1:6")));
  EXPECT_FALSE(q.Contains(MakePointOrDie("10:10")));

  // Test the last few cases using the Init() method instead.
  S2ContainsPointQuery<MutableS2ShapeIndex> q2;
  q2.Init(index.get(), options);
  EXPECT_FALSE(q2.ShapeContains(1, MakePointOrDie("1:6")));
  EXPECT_TRUE(q2.ShapeContains(2, MakePointOrDie("1:6")));
  EXPECT_FALSE(q2.ShapeContains(2, MakePointOrDie("0:5")));
  EXPECT_FALSE(q2.ShapeContains(2, MakePointOrDie("0:7")));
}

TEST(S2ContainsPointQuery, VertexModelSemiOpen) {
  auto index = MakeIndexOrDie("0:0 # -1:1, 1:1 # 0:5, 0:7, 2:6");
  S2ContainsPointQueryOptions options(S2VertexModel::SEMI_OPEN);
  auto q = MakeS2ContainsPointQuery(index.get(), options);
  EXPECT_FALSE(q.Contains(MakePointOrDie("0:0")));
  EXPECT_FALSE(q.Contains(MakePointOrDie("-1:1")));
  EXPECT_FALSE(q.Contains(MakePointOrDie("1:1")));
  EXPECT_FALSE(q.Contains(MakePointOrDie("0:2")));
  EXPECT_FALSE(q.Contains(MakePointOrDie("0:5")));
  EXPECT_TRUE(q.Contains(MakePointOrDie("0:7")));  // Contained vertex.
  EXPECT_FALSE(q.Contains(MakePointOrDie("2:6")));
  EXPECT_TRUE(q.Contains(MakePointOrDie("1:6")));
  EXPECT_FALSE(q.Contains(MakePointOrDie("10:10")));

  // Test the last few cases using the Init() method instead.
  S2ContainsPointQuery<MutableS2ShapeIndex> q2;
  q2.Init(index.get(), options);
  EXPECT_FALSE(q2.ShapeContains(1, MakePointOrDie("1:6")));
  EXPECT_TRUE(q2.ShapeContains(2, MakePointOrDie("1:6")));
  EXPECT_FALSE(q2.ShapeContains(2, MakePointOrDie("0:5")));
  EXPECT_TRUE(q2.ShapeContains(2, MakePointOrDie("0:7")));
}

TEST(S2ContainsPointQuery, VertexModelClosed) {
  auto index = MakeIndexOrDie("0:0 # -1:1, 1:1 # 0:5, 0:7, 2:6");
  S2ContainsPointQueryOptions options(S2VertexModel::CLOSED);
  auto q = MakeS2ContainsPointQuery(index.get(), options);
  EXPECT_TRUE(q.Contains(MakePointOrDie("0:0")));
  EXPECT_TRUE(q.Contains(MakePointOrDie("-1:1")));
  EXPECT_TRUE(q.Contains(MakePointOrDie("1:1")));
  EXPECT_FALSE(q.Contains(MakePointOrDie("0:2")));
  EXPECT_TRUE(q.Contains(MakePointOrDie("0:5")));
  EXPECT_TRUE(q.Contains(MakePointOrDie("0:7")));
  EXPECT_TRUE(q.Contains(MakePointOrDie("2:6")));
  EXPECT_TRUE(q.Contains(MakePointOrDie("1:6")));
  EXPECT_FALSE(q.Contains(MakePointOrDie("10:10")));

  // Test the last few cases using the Init() method instead.
  S2ContainsPointQuery<MutableS2ShapeIndex> q2;
  q2.Init(index.get(), options);
  EXPECT_FALSE(q2.ShapeContains(1, MakePointOrDie("1:6")));
  EXPECT_TRUE(q2.ShapeContains(2, MakePointOrDie("1:6")));
  EXPECT_TRUE(q2.ShapeContains(2, MakePointOrDie("0:5")));
  EXPECT_TRUE(q2.ShapeContains(2, MakePointOrDie("0:7")));
}

TEST(S2ContainsPointQuery, VisitContainingShapesCanStopEarly) {
  auto index = MakeIndexOrDie("0:0 # 0:0, 1:1 # -1:0, 0:1, 1:0, 0:-1");
  const S2Point kPoint = S2LatLng::FromDegrees(0, 0).ToPoint();

  // Under a closed vertex model there should be 3 shapes that contain 0,0.
  S2ContainsPointQueryOptions options;
  options.set_vertex_model(S2VertexModel::CLOSED);
  auto query = MakeS2ContainsPointQuery(index.get(), options);

  // But if we return false, we should only see the first one.
  int count = 0;
  bool status = query.VisitContainingShapes(kPoint, [&](const S2Shape*) {
    ++count;
    return false;
  });
  EXPECT_FALSE(status);
  EXPECT_EQ(count, 1);
}

TEST(S2ContainsPointQuery, GetContainingShapes) {
  // Also tests ShapeContains().
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "GET_CONTAINING_SHAPES",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  constexpr int kNumVerticesPerLoop = 10;
  const S1Angle kMaxLoopRadius = S2Testing::KmToAngle(10);
  const S2Cap center_cap(s2random::Point(bitgen), kMaxLoopRadius);
  MutableS2ShapeIndex index;
  for (int i = 0; i < 100; ++i) {
    unique_ptr<S2Loop> loop = S2Loop::MakeRegularLoop(
        s2random::SamplePoint(bitgen, center_cap),
        absl::Uniform(bitgen, 0.0, 1.0) * kMaxLoopRadius, kNumVerticesPerLoop);
    index.Add(make_unique<S2Loop::OwningShape>(std::move(loop)));
  }
  auto query = MakeS2ContainsPointQuery(&index);
  for (int i = 0; i < 100; ++i) {
    S2Point p = s2random::SamplePoint(bitgen, center_cap);

    vector<int> expected_ids;
    vector<const S2Shape*> expected_shapes;

    int shape_id = 0;
    for (const S2Shape* shape : index) {
      const S2Loop* loop = down_cast<const S2Loop::Shape*>(shape)->loop();
      if (loop->Contains(p)) {
        EXPECT_TRUE(query.ShapeContains(shape_id, p));
        expected_shapes.push_back(shape);
        expected_ids.push_back(shape_id);
      } else {
        EXPECT_FALSE(query.ShapeContains(shape_id, p));
      }
      ++shape_id;
    }
    EXPECT_EQ(query.GetContainingShapes(p), expected_shapes);
    EXPECT_EQ(query.GetContainingShapeIds(p), expected_ids);
  }
}

using EdgeIdVector = vector<ShapeEdgeId>;

void ExpectIncidentEdgeIds(const EdgeIdVector& expected,
                           const MutableS2ShapeIndex& index, const S2Point& p) {
  EdgeIdVector actual;
  auto q = MakeS2ContainsPointQuery(&index);
  EXPECT_TRUE(
      q.VisitIncidentEdges(p, [&actual](const s2shapeutil::ShapeEdge& e) {
          actual.push_back(e.id());
          return true;
        }));
  EXPECT_EQ(expected, actual);
}

TEST(S2ContainsPointQuery, VisitIncidentEdges) {
  auto index = MakeIndexOrDie("0:0 | 1:1 # 1:1, 1:2 # 1:2, 1:3, 2:2");
  ExpectIncidentEdgeIds({{0, 0}}, *index, MakePointOrDie("0:0"));
  ExpectIncidentEdgeIds({{0, 1}, {1, 0}}, *index, MakePointOrDie("1:1"));
  ExpectIncidentEdgeIds({{1, 0}, {2, 0}, {2, 2}}, *index,
                        MakePointOrDie("1:2"));
  ExpectIncidentEdgeIds({{2, 0}, {2, 1}}, *index, MakePointOrDie("1:3"));
  ExpectIncidentEdgeIds({{2, 1}, {2, 2}}, *index, MakePointOrDie("2:2"));
}
