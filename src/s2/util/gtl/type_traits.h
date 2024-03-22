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

#ifndef S2_UTIL_GTL_TYPE_TRAITS_H_
#define S2_UTIL_GTL_TYPE_TRAITS_H_

#include "absl/meta/type_traits.h"

#ifdef SWIG
%include "third_party/absl/meta/type_traits.h"
#endif

namespace gtl {

// Provides a member using-declaration for T, i.e. the identity transformation.
// Drop-in replacement for C++20's std::type_identity.
template <typename T>
struct type_identity {
  using type = T;
};

template <typename T>
using type_identity_t = typename type_identity<T>::type;

// Removes any references and cv qualifiers from the provided type. Drop-in
// replacement for C++20's std::remove_cvref.
template <typename T>
struct remove_cvref {
  using type = absl::remove_cv_t<absl::remove_reference_t<T>>;
};
template <typename T>
using remove_cvref_t = typename remove_cvref<T>::type;

}  // namespace gtl

#endif  // S2_UTIL_GTL_TYPE_TRAITS_H_
