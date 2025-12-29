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

#ifndef S2_S2POINT_ARRAY_H_
#define S2_S2POINT_ARRAY_H_

// `S2Point` has a default constructor that zero-initializes the memory,
// which is wasteful when it will be overwritten.  This file provides
// `MakeS2PointArrayForOverwrite()` to create an array of `S2Point`s that
// is uninitialized.

#include <cstdint>
#include <memory>
#include <new>

#include "absl/log/absl_check.h"
#include "s2/s2point.h"

namespace s2internal {

struct S2PointArrayDeleter {
  void operator()(S2Point* p) const {
    ::operator delete(p, std::align_val_t(alignof(*p)));
  }
};

// TODO(b/469042890): Make this a class with a `num_vertices_` field.
// Use `std::allocator` to allocate/deallocate.  This saves some space
// since the size does not have to be stored by the allocator.
using UniqueS2PointArray = std::unique_ptr<S2Point[], S2PointArrayDeleter>;

inline UniqueS2PointArray MakeS2PointArrayForOverwrite(int32_t size) {
  // Note that even `std::make_unique_for_overwrite<S2Point[]>(size)`
  // zero-initializes the memory, so we use `::operator new` directly.
  ABSL_DCHECK_GE(size, 0);
  S2Point* p = static_cast<S2Point*>(
      ::operator new(size * sizeof(*p), std::align_val_t(alignof(*p))));
  return UniqueS2PointArray(p);
}

}  // namespace s2internal

#endif  // S2_S2POINT_ARRAY_H_
