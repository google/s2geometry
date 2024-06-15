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

#ifndef S2_UTIL_GTL_UNALIGNED_INTERNAL_H_
#define S2_UTIL_GTL_UNALIGNED_INTERNAL_H_

#include <array>
#include <cstddef>
#include <type_traits>
#include <utility>

#include "absl/base/casts.h"
#include "s2/util/gtl/requires.h"

namespace gtl {

template <typename T>
class Unaligned;

namespace internal_unaligned {

template <typename T>
constexpr bool IsBitCastableTo = gtl::Requires<std::array<char, sizeof(T)>>(
    [](auto&& x) -> decltype(absl::bit_cast<T>(x)) {});

struct NoDefaultConstructor {
  explicit NoDefaultConstructor(std::in_place_t) {}
};
struct TrivialDefaultConstructor {
  TrivialDefaultConstructor() = default;
  explicit TrivialDefaultConstructor(std::in_place_t) {}
};
template <typename Derived>
struct UserDefaultConstructor {
  UserDefaultConstructor() { static_cast<Derived*>(this)->Store({}); }
  explicit UserDefaultConstructor(std::in_place_t) {}
};

template <typename T>
using UnalignedBase = std::conditional_t<
    std::is_default_constructible_v<T>,
    std::conditional_t<std::is_trivially_default_constructible_v<T>,
                       TrivialDefaultConstructor,
                       UserDefaultConstructor<Unaligned<T>>>,
    NoDefaultConstructor>;

}  // namespace internal_unaligned
}  // namespace gtl

#endif  // S2_UTIL_GTL_UNALIGNED_INTERNAL_H_
