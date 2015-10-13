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

#include "s2closestedgequery.h"

#include <memory>

#include <gtest/gtest.h>
#include "s1angle.h"
#include "s2cap.h"
#include "s2loop.h"
#include "s2testing.h"

using std::unique_ptr;

class S2ClosestEdgeQueryTest : public ::testing::Test {
 protected:
  S2ShapeIndex index_;
  void TestDistance(S2Point const& target, S1Angle max_error) const;
};

void S2ClosestEdgeQueryTest::TestDistance(S2Point const& target,
                                          S1Angle max_error) const {
  S2ClosestEdgeQuery query(index_);
  query.UseBruteForce(true);
  S1Angle true_min_dist = query.GetDistance(target);
  query.UseBruteForce(false);
  S1Angle min_dist = query.ApproxGetDistance(target, max_error);
  EXPECT_LE(min_dist, true_min_dist + max_error + S1Angle::Radians(1e-15));
  S2Point closest = query.ApproxProject(target, max_error);
  EXPECT_NEAR(min_dist.radians(), S1Angle(target, closest).radians(), 1e-15);
}

TEST_F(S2ClosestEdgeQueryTest, RegularLoop) {
  // The running time of this test is proportional to
  //    kLoopIters * kQueryIters * kNumVertices.
  static int const kLoopIters = 5;
  static int const kQueryIters = 20;
  static int const kNumVertices = 1000;
  static S1Angle const kRadius = S2Testing::KmToAngle(10);
  for (int i_loop = 0; i_loop < kLoopIters; ++i_loop) {
    S2Point center = S2Testing::RandomPoint();
    unique_ptr<S2Loop> loop(S2Testing::MakeRegularLoop(center, kRadius,
                                                       kNumVertices));
    index_.Add(new S2Loop::Shape(loop.get()));
    // Choose query points from a region whose area is approximately 4 times
    // larger than the loop, i.e. such that about 1/4 of the query points will
    // be inside the loop.
    S2Cap query_cap(center, 2 * kRadius);
    for (int i_query = 0; i_query < kQueryIters; ++i_query) {
      S2Point target = S2Testing::SamplePoint(query_cap);
      // Choose a maximum error whose logarithm is uniformly distributed over
      // a reasonable range, except that it is sometimes zero.
      S1Angle max_error = kRadius * pow(1e-4, S2Testing::rnd.RandDouble());
      if (S2Testing::rnd.OneIn(5)) max_error = S1Angle::Zero();
      TestDistance(target, max_error);
    }
    index_.Reset();
  }
}

TEST_F(S2ClosestEdgeQueryTest, NoEdges) {
  S2ShapeIndex index;
  S2ClosestEdgeQuery query(index);
  EXPECT_EQ(S1Angle::Infinity(), query.GetDistance(S2Point(1, 0, 0)));
  EXPECT_EQ(-1, query.shape_id());
  EXPECT_EQ(-1, query.edge_id());
  EXPECT_EQ(S2Point(1, 0, 0), query.Project(S2Point(1, 0, 0)));
  EXPECT_EQ(-1, query.shape_id());
  EXPECT_EQ(-1, query.edge_id());
}

