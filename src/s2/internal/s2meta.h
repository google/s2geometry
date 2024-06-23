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

#ifndef S2_INTERNAL_S2META_H_
#define S2_INTERNAL_S2META_H_

#include <type_traits>

#include "absl/meta/type_traits.h"

// This file contains functions for meta programming, primarily related to
// checking that the type of an opaque template parameter is an instance of some
// S2 base type.  These are defined internally to avoid clients relying on them.

// A variable template to check if a type is derived from another type ignoring
// cv and reference qualifiers. This is similar to std::derived_from which was
// added in C++20.
//
// TODO(b/335677213): Remove this when we have access to C++20.

namespace s2meta {

template <typename Derived, typename Base>
using derived_from = absl::conjunction<
    std::is_base_of<std::decay_t<Base>, std::decay_t<Derived>>,
    std::is_convertible<const volatile std::decay_t<Derived>*,
                        const volatile std::decay_t<Base>*>>;

template <typename Derived, typename Base>
constexpr bool derived_from_v = derived_from<Derived, Base>::value;

}  // namespace s2meta

#endif  // S2_INTERNAL_S2META_H_
