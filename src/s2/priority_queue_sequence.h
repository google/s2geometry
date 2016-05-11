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

#ifndef S2_PRIORITY_QUEUE_SEQUENCE_H_
#define S2_PRIORITY_QUEUE_SEQUENCE_H_

#include <functional>
#include <queue>

#include "s2/util/gtl/inlined_vector.h"

// Like std::priority_queue, except that:
//
//  (1) the underlying sequence type is accessible via rep() and mutable_rep()
//      methods, so that clients can manipulate it directly if desired; and
//  (2) the default sequence type is gtl::InlinedVector.
//
// Exposing the sequence type gives more flexibility in how it is used,
// e.g. the queue can be emptied by calling clear(), the elements in the
// queue can be examined without removing them, etc.  Note that if you modify
// the underlying sequence, then you must call make_heap() before calling any
// of the priority queue methods (push, pop, etc).
//
// Using InlinedVector as the default sequence type increases efficiency in
// cases where the maximum queue size is small.

template <class T, class Sequence = gtl::InlinedVector<T, 16>,
          class Compare = std::less<typename Sequence::value_type>>
class priority_queue_sequence
    : public std::priority_queue<T, Sequence, Compare> {
 public:
  explicit priority_queue_sequence(const Compare& cmp = Compare(),
                                   const Sequence& seq = Sequence())
      : std::priority_queue<T, Sequence, Compare>(cmp, seq) {
  }

  template <typename InputIterator>
  priority_queue_sequence(InputIterator first, InputIterator last,
                          const Compare& cmp = Compare(),
                          const Sequence& seq = Sequence())
      : std::priority_queue<T, Sequence, Compare>(first, last, cmp, seq) {
  }

  // Returns the underlying sequence.
  Sequence const& rep() const { return this->c; }

  // Allows the underlying sequence to be modified.  If you do this, you must
  // call make_heap() before calling any of the other priority queue methods
  // defined below.
  Sequence* mutable_rep() { return &this->c; }

  // Restores the priority queue invariants after the underlying sequence has
  // been modified using mutable_rep().
  void make_heap() {
    std::make_heap(this->c.begin(), this->c.end(), this->comp);
  }
};

#endif  // S2_PRIORITY_QUEUE_SEQUENCE_H_
