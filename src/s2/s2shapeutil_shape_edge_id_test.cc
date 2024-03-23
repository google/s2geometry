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

#include "s2/s2shapeutil_shape_edge_id.h"

#include <gmock/gmock.h>
#include <gtest/gtest.h>

namespace {

using ::s2shapeutil::ShapeEdgeId;
using ::testing::Eq;
using ::testing::Ge;
using ::testing::Gt;
using ::testing::Le;
using ::testing::Lt;
using ::testing::Ne;
using ::testing::Not;

TEST(ShapeEdgeIdTest, BothFieldsEqualIsEqual) {
  EXPECT_THAT(ShapeEdgeId(10, 20), Eq(ShapeEdgeId(10, 20)));
  EXPECT_THAT(ShapeEdgeId(10, 20), Not(Ne(ShapeEdgeId(10, 20))));
}

TEST(ShapeEdgeIdTest, BothShapeIdUnequalIsUnequal) {
  EXPECT_THAT(ShapeEdgeId(11, 20), Ne(ShapeEdgeId(10, 20)));
  EXPECT_THAT(ShapeEdgeId(11, 20), Not(Eq(ShapeEdgeId(10, 20))));
}

TEST(ShapeEdgeIdTest, BothEdgeIdUnequalIsUnequal) {
  EXPECT_THAT(ShapeEdgeId(10, 21), Ne(ShapeEdgeId(10, 20)));
  EXPECT_THAT(ShapeEdgeId(10, 21), Not(Eq(ShapeEdgeId(10, 20))));
}

TEST(ShapeEdgeIdTest, LessThanIsLexicographicShapeIdFirst) {
  EXPECT_THAT(ShapeEdgeId(10, 20), Not(Lt(ShapeEdgeId(10, 20))));
  EXPECT_THAT(ShapeEdgeId(10, 20), Lt(ShapeEdgeId(11, 20)));

  EXPECT_THAT(ShapeEdgeId(10, 20), Not(Lt(ShapeEdgeId(10, 20))));
  EXPECT_THAT(ShapeEdgeId(10, 20), Lt(ShapeEdgeId(10, 21)));
}

TEST(ShapeEdgeIdTest, LessEqIsLexicographicShapeIdFirst) {
  EXPECT_THAT(ShapeEdgeId(10, 20), Not(Le(ShapeEdgeId(9, 20))));
  EXPECT_THAT(ShapeEdgeId(10, 20), Le(ShapeEdgeId(10, 20)));

  EXPECT_THAT(ShapeEdgeId(10, 20), Not(Le(ShapeEdgeId(10, 19))));
  EXPECT_THAT(ShapeEdgeId(10, 20), Le(ShapeEdgeId(10, 20)));
}

TEST(ShapeEdgeIdTest, GreaterThanIsLexicographicShapeIdFirst) {
  EXPECT_THAT(ShapeEdgeId(10, 20), Not(Gt(ShapeEdgeId(10, 20))));
  EXPECT_THAT(ShapeEdgeId(10, 20), Gt(ShapeEdgeId(9, 20)));

  EXPECT_THAT(ShapeEdgeId(10, 20), Not(Gt(ShapeEdgeId(10, 20))));
  EXPECT_THAT(ShapeEdgeId(10, 20), Gt(ShapeEdgeId(10, 19)));
}

TEST(ShapeEdgeIdTest, GreaterEqIsLexicographicShapeIdFirst) {
  EXPECT_THAT(ShapeEdgeId(10, 20), Not(Ge(ShapeEdgeId(11, 20))));
  EXPECT_THAT(ShapeEdgeId(10, 20), Ge(ShapeEdgeId(10, 20)));

  EXPECT_THAT(ShapeEdgeId(10, 20), Not(Ge(ShapeEdgeId(10, 21))));
  EXPECT_THAT(ShapeEdgeId(10, 20), Ge(ShapeEdgeId(10, 20)));
}

}  // namespace
