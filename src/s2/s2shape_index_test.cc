// Copyright 2012 Google Inc. All Rights Reserved.
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

#include "s2/s2shape_index.h"

#include <gtest/gtest.h>
#include "s2/s2shape.h"
#include "s2/s2text_format.h"

// TODO(ericv): Add tests for S2ShapeIndexCell and S2ClippedShape.
// (Currently these are tested indirectly by MutableS2ShapeIndex.)
// Also test the base Iterator type (which wraps another iterator).

TEST(S2ShapeIndexIterator, PrefixIncrement) {
  const auto index =
      s2textformat::MakeIndexOrDie("1:1 # 1:1, 2:2 # 2:2, 1:3, 2:4, 3:2");

  int count = 0;
  for (auto iter = index->begin(); iter != index->end(); ++iter) {
    ++count;
  }
  EXPECT_EQ(count, 3);
}

TEST(S2ShapeIndexIterator, PostfixIncrement) {
  const auto index =
      s2textformat::MakeIndexOrDie("1:1 # 1:1, 2:2 # 2:2, 1:3, 2:4, 3:2");

  int count = 0;
  for (auto iter = index->begin(); iter != index->end(); iter++) {
    ++count;
  }
  EXPECT_EQ(count, 3);
}
