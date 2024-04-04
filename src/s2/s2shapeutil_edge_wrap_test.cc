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


#include "s2/s2shapeutil_edge_wrap.h"

#include <memory>

#include <gtest/gtest.h>
#include "s2/mutable_s2shape_index.h"
#include "s2/s2shape.h"
#include "s2/s2text_format.h"

namespace {

using s2shapeutil::NextEdgeWrap;
using s2shapeutil::PrevEdgeWrap;

TEST(S2Shape, NextPrevEdgePointDoesNotWrap) {
  auto index = s2textformat::MakeIndexOrDie("1:1 | 2:2 ##");
  const S2Shape& shape = *index->shape(0);

  // Points have one chain per point so we should always get -1.
  EXPECT_EQ(PrevEdgeWrap(shape, 0), -1);
  EXPECT_EQ(NextEdgeWrap(shape, 0), -1);

  EXPECT_EQ(PrevEdgeWrap(shape, 1), -1);
  EXPECT_EQ(NextEdgeWrap(shape, 1), -1);
}

TEST(S2Shape, NextPrevEdgeOpenPolylineDoesNotWrap) {
  auto index = s2textformat::MakeIndexOrDie("# 1:1, 2:2, 3:3 #");
  const S2Shape& shape = *index->shape(0);

  // Open polylines should not wrap around.
  EXPECT_EQ(PrevEdgeWrap(shape, 0), -1);
  EXPECT_EQ(NextEdgeWrap(shape, 0), 1);

  EXPECT_EQ(PrevEdgeWrap(shape, 1), 0);
  EXPECT_EQ(NextEdgeWrap(shape, 1), -1);
}

TEST(S2Shape, NextPrevEdgeClosedPolylineWraps) {
  auto index = s2textformat::MakeIndexOrDie("# 0:0, 1:1, 0:2, -1:1, 0:0 #");
  const S2Shape& shape = *index->shape(0);

  // Closed polylines should wrap around.
  EXPECT_EQ(PrevEdgeWrap(shape, 0), 3);
  EXPECT_EQ(NextEdgeWrap(shape, 0), 1);

  EXPECT_EQ(PrevEdgeWrap(shape, 3), 2);
  EXPECT_EQ(NextEdgeWrap(shape, 3), 0);
}

TEST(S2Shape, NextPrevEdgePolygonWraps) {
  auto index = s2textformat::MakeIndexOrDie("## 0:0, 1:1, 0:2, -1:1");
  const S2Shape& shape = *index->shape(0);

  // Polygons should always wrap.
  EXPECT_EQ(PrevEdgeWrap(shape, 0), 3);
  EXPECT_EQ(NextEdgeWrap(shape, 0), 1);

  EXPECT_EQ(PrevEdgeWrap(shape, 3), 2);
  EXPECT_EQ(NextEdgeWrap(shape, 3), 0);
}

}  // namespace
