// Copyright 2016 Google Inc. All Rights Reserved.
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

// Utility functions for pointers.

#ifndef S2_UTIL_GTL_PTR_UTIL_H_
#define S2_UTIL_GTL_PTR_UTIL_H_

#include <cstddef>

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

#include "s2/util/gtl/gtl_namespace.h"

GTL_NAMESPACE_BEGIN

// Transfers ownership of a raw pointer to a std::unique_ptr of deduced type.
// Example:
//   X* NewX(int, int);
//   auto x = WrapUnique(NewX(1, 2));  // 'x' is std::unique_ptr<X>.
//
// WrapUnique is useful for capturing the output of a raw pointer factory.
// However, prefer 'MakeUnique<T>(args...) over 'WrapUnique(new T(args...))'.
//   auto x = WrapUnique(new X(1, 2));  // works, but nonideal.
//   auto x = MakeUnique<X>(1, 2);  // safer, standard, avoids raw 'new'.
//
// Note: Cannot wrap pointers to array of unknown bound (i.e. U(*)[]).
template <typename T>
std::unique_ptr<T> WrapUnique(T* ptr) {
  static_assert(
      !std::is_array<T>::value || std::extent<T>::value != 0,
      "types T[0] or T[] are unsupported");
  return std::unique_ptr<T>(ptr);
}


namespace internal {

// Trait to select overloads and return types for MakeUnique.
template <typename T>
struct MakeUniqueResult {
  using scalar = std::unique_ptr<T>;
};
template <typename T>
struct MakeUniqueResult<T[]> {
  using array = std::unique_ptr<T[]>;
};
template <typename T, size_t N>
struct MakeUniqueResult<T[N]> {
  using invalid = void;
};

}  // namespace internal

// MakeUnique<T>(...) is an early implementation of C++14 std::make_unique.
// It is designed to be 100% compatible with std::make_unique so that the
// eventual switchover will be a simple renaming operation.
//
// For a std::make_unique API description, see
// [http://en.cppreference.com/w/cpp/memory/unique_ptr/make_unique].
//
// For some technical motivation to prefer 'MakeUnique<T>(a,b)' over
// 'std::unique_ptr<T>(new T(a,b))', see [http://herbsutter.com/gotw/_102/].
// There's reviewer value in treating 'new T(a,b)' with scrutiny.
//
// Example usage:
//
//    auto p = MakeUnique<X>(args...);  // 'p' is std::unique_ptr<X>
//    auto pa = MakeUnique<X[]>(5);  // 'pa' is std::unique_ptr<X[]>
//
// Three overloads are needed:
//
//   - For non-array T:
//
//      template <typename T, typename... Args>
//      std::unique_ptr<T> MakeUnique<T>(Args&&... args);
//
//        Allocates a T with 'new T(std::forward<Args> args...),
//        forwarding all 'args' to T's constructor.
//        Returns a std::unique_ptr<T> owning that object.
//
//   - For T that are array of unknown bound, say U[]:
//
//      template <typename T>
//      std::unique_ptr<T> MakeUnique(size_t n);
//
//     Suppose U such that T is U[]. Allocates an array with 'new U[n]()'.
//     Returns a std::unique_ptr<U[]> owning that array.
//     Note that 'U[n]()' is different from 'U[n]', and elements will be
//     value-initialized.
//
//   - Like std::make_unique, MakeUnique is deleted for T that are array
//     of known bound. Suppose U and N such that T is U[N].
//     std::unique_ptr<U[N]> is not useful, and MakeUnique<T> should return
//     std::unique_ptr<T> for any T.
//     See [http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2013/n3588.txt].
//     You can use MakeUnique<std::array<U, N>>() to approximate this if needed.
//
// (Rvalue references used by arbiter approval CL/62311378)
// Overload for non-array types.
template <typename T, typename... Args>
typename internal::MakeUniqueResult<T>::scalar
MakeUnique(Args&&... args) {  // NOLINT
  return std::unique_ptr<T>(
      new T(std::forward<Args>(args)...));  // NOLINT(build/c++11)
}

// Overload for array of unknown bound.
// The allocation of arrays needs to use the array form of new,
// and cannot take element constructor arguments.
template <typename T>
typename internal::MakeUniqueResult<T>::array
MakeUnique(size_t n) {
  return std::unique_ptr<T>(new typename std::remove_extent<T>::type[n]());
}

// Reject arrays of known bound.
template <typename T, typename... Args>
typename internal::MakeUniqueResult<T>::invalid
MakeUnique(Args&&... /* args */) = delete;  // NOLINT

// RawPtr(ptr) extracts the raw pointer from pointer-like 'ptr'.
template <typename T>
auto RawPtr(T&& ptr) -> decltype(&*ptr) {  // NOLINT
  // ptr is a universal reference to support Ts with non-const operators.
  return (ptr != nullptr) ? &*ptr : nullptr;
}
inline std::nullptr_t RawPtr(std::nullptr_t) { return nullptr; }

GTL_NAMESPACE_END

#endif  // S2_UTIL_GTL_PTR_UTIL_H_
