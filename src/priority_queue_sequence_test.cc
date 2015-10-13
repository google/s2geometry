// Copyright 2015 Google Inc. All Rights Reserved.
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

#include "priority_queue_sequence.h"

#include <algorithm>
#include "s2testing.h"
#include <gtest/gtest.h>

using std::vector;

TEST(PriorityQueueSequence, RandomOrder) {
  for (int iter = 0; iter < 1000; ++iter) {
    // Insert the numbers from 0 to n-1 in a random order, using the priority
    // queue to keep track of the k smallest numbers.  The verify that at
    // the end, the priority queue contains the number from 0 to k-1.
    int n = 1 + S2Testing::rnd.Uniform(100);
    int k = 1 + S2Testing::rnd.Uniform(n);
    vector<int> values(n);
    for (int i = 0; i < n; ++i) {
      values[i] = i;
    }
    std::random_shuffle(values.begin(), values.end(), S2Testing::rnd);
    priority_queue_sequence<int> queue;
    for (int i = 0; i < k; ++i) {
      queue.mutable_rep()->push_back(values[i]);
    }
    queue.make_heap();
    for (int i = k; i < n; ++i) {
      if (values[i] >= queue.top()) continue;
      queue.pop();
      queue.push(values[i]);
    }
    // Access the sequence through rep().
    EXPECT_EQ(k, queue.rep().size());
    for (int i = k; i > 0; --i) {
      EXPECT_EQ(i, queue.size());
      EXPECT_EQ(i - 1, queue.top());
      queue.pop();
    }
    EXPECT_TRUE(queue.empty());
  }
}
