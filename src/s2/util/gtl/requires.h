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

#ifndef S2_UTIL_GTL_REQUIRES_H_
#define S2_UTIL_GTL_REQUIRES_H_

#include <type_traits>

namespace gtl {

// C++17 port of C++20's `requires` expressions.
// https://en.cppreference.com/w/cpp/language/constraints#Requires_expressions
template <typename... T, typename F>
constexpr bool Requires(F) {
  return std::is_invocable_v<F, T...>;
}

}  // namespace gtl

#endif  // S2_UTIL_GTL_REQUIRES_H_
