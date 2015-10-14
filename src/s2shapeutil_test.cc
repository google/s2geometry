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

#include "s2shapeutil.h"

#include <gtest/gtest.h>
#include "s2testing.h"

namespace {

using s2shapeutil::S2EdgeVectorShape;

TEST(S2EdgeVectorShape, EdgeAccess) {
  S2EdgeVectorShape shape;
  S2Testing::rnd.Reset(FLAGS_s2_random_seed);
  int const kNumEdges = 100;
  for (int i = 0; i < kNumEdges; ++i) {
    S2Point a = S2Testing::RandomPoint();  // Control the evaluation order
    shape.Add(a, S2Testing::RandomPoint());
  }
  EXPECT_EQ(kNumEdges, shape.num_edges());
  S2Testing::rnd.Reset(FLAGS_s2_random_seed);
  for (int i = 0; i < kNumEdges; ++i) {
    S2Point const *a, *b;
    shape.GetEdge(i, &a, &b);
    EXPECT_EQ(S2Testing::RandomPoint(), *a);
    EXPECT_EQ(S2Testing::RandomPoint(), *b);
  }
}

TEST(S2EdgeVectorShape, SingletonConstructor) {
  S2Point a(1, 0, 0), b(0, 1, 0);
  S2EdgeVectorShape shape(a, b);
  EXPECT_EQ(1, shape.num_edges());
  S2Point const *pa, *pb;
  shape.GetEdge(0, &pa, &pb);
  EXPECT_EQ(a, *pa);
  EXPECT_EQ(b, *pb);
}

TEST(S2EdgeVectorShape, Ownership) {
  S2EdgeVectorShape* shape = new S2EdgeVectorShape;
  shape->Release();  // Verify there is no memory leak.
}

}  // namespace
