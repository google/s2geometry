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

#include "s2/s2shapeutil_count_vertices.h"

#include <memory>

#include <gtest/gtest.h>
#include "s2/mutable_s2shape_index.h"
#include "s2/s2text_format.h"

namespace {

using ::s2shapeutil::CountVertices;
using ::std::unique_ptr;

TEST(CountVertices, CountsCorrectly) {
  unique_ptr<MutableS2ShapeIndex> index;

  // Test index built only out of three points.
  index = s2textformat::MakeIndexOrDie("1:1 | 2:2 | 3:3 # #");
  ASSERT_EQ(3, CountVertices(*index));

  // Test index built out of two points and a two edge polyline.
  index = s2textformat::MakeIndexOrDie("1:1 | 2:2 # 3:3, 4:4, 5:5 #");
  ASSERT_EQ(5, CountVertices(*index));

  // Test index built out of two points, one two edge polyline, and four edge
  // polygon.
  index = s2textformat::MakeIndexOrDie(
      "1:1 | 2:2 # 3:3, 4:4, 5:5 # 6:6, 7:7, 8:8, 9:9");
  ASSERT_EQ(9, CountVertices(*index));

  // Test that degenerate polylines count correctly.
  index = s2textformat::MakeIndexOrDie("# 3:3, 3:3, 3:3 #");
  ASSERT_EQ(3, CountVertices(*index));

  // Test that degenerate polygons count correctly.
  index = s2textformat::MakeIndexOrDie("# # 4:4, 4:4, 4:4, 4:4");
  ASSERT_EQ(4, CountVertices(*index));
}

}  // namespace
