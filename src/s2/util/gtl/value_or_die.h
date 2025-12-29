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

#ifndef S2_UTIL_GTL_VALUE_OR_DIE_H_
#define S2_UTIL_GTL_VALUE_OR_DIE_H_

#include <utility>

namespace gtl {

// Equivalent to `std::move(v).value()`, but makes explicit that the behavior
// is to die if `v` does not contain a value.
template <typename T>
T ValueOrDie(absl::StatusOr<T>&& v) {
  return std::move(v).value();
}

}  // namespace gtl

#endif  // S2_UTIL_GTL_VALUE_OR_DIE_H_
