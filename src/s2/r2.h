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

#ifndef S2_R2_H_
#define S2_R2_H_

#include "s2/_fp_contract_off.h"  // IWYU pragma: keep
#include "s2/util/math/vector.h"  // IWYU pragma: export

using R2Point = Vector2_d;

// An edge in R2 space.
struct R2Edge {
  R2Edge() = default;
  R2Edge(const R2Point& v0, const R2Point& v1) : v0(v0), v1(v1) {}

  bool operator==(const R2Edge& b) const { return v0 == b.v0 && v1 == b.v1; }

  template <typename H>
  friend H AbslHashValue(H h, const R2Edge& e) {
    return H::combine(std::move(h), e.v0, e.v1);
  }

  R2Point v0;
  R2Point v1;
};

#endif  // S2_R2_H_
