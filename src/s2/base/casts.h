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

#ifndef S2_BASE_CASTS_H_
#define S2_BASE_CASTS_H_

#include <cassert>
#include <type_traits>

#include "absl/base/config.h"

template <typename To, typename From>
inline To down_cast(From* f) {
  static_assert((std::is_base_of<From, std::remove_pointer_t<To>>::value),
                "target type not derived from source type");

  // We depend on abseil-cpp internals, but if this doesn't get defined,
  // we will just not have the assert.
#if ABSL_INTERNAL_HAS_RTTI
  assert(f == nullptr || dynamic_cast<To>(f) != nullptr);
#endif

  return static_cast<To>(f);
}

// Overload of down_cast for references. Use like this: down_cast<T&>(foo).
// The code is slightly convoluted because we're still using the pointer
// form of dynamic cast. (The reference form throws an exception if it
// fails.)
//
// There's no need for a special const overload either for the pointer
// or the reference form. If you call down_cast with a const T&, the
// compiler will just bind From to const T.
template <typename To, typename From>
inline To down_cast(From& f) {
  static_assert(std::is_lvalue_reference<To>::value,
                "target type not a reference");
  static_assert((std::is_base_of<From, std::remove_reference_t<To>>::value),
                "target type not derived from source type");

  // We skip the assert and hence the dynamic_cast if RTTI is disabled.
#if ABSL_INTERNAL_HAS_RTTI
  assert(dynamic_cast<std::remove_reference_t<To>*>(&f) != nullptr);
#endif

  return static_cast<To>(f);
}

#endif  // S2_BASE_CASTS_H_
