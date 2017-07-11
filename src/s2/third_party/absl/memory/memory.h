// Copyright 2017 Google Inc. All Rights Reserved.
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

// -----------------------------------------------------------------------------
// File: memory.h
// -----------------------------------------------------------------------------
//
// This header file contains utility functions for managing the creation and
// conversion of smart pointers. This file is an extension to the C++
// standard <memory> library header file.

#ifndef S2_THIRD_PARTY_ABSL_MEMORY_MEMORY_H_
#define S2_THIRD_PARTY_ABSL_MEMORY_MEMORY_H_

#include <cstddef>

#include <cstddef>
#include <memory>
#include <type_traits>
#include <utility>

#include "s2/third_party/absl/meta/type_traits.h"

namespace absl {

// -----------------------------------------------------------------------------
// Function Template: WrapUnique()
// -----------------------------------------------------------------------------
//
// Transfers ownership of a raw pointer to a `std::unique_ptr`. The returned
// value is a `std::unique_ptr` of deduced type.
//
// Example:
//   X* NewX(int, int);
//   auto x = WrapUnique(NewX(1, 2));  // 'x' is std::unique_ptr<X>.
//
// `absl::WrapUnique` is useful for capturing the output of a raw pointer
// factory. However, prefer 'absl::MakeUnique<T>(args...) over
// 'absl::WrapUnique(new T(args...))'.
//
//   auto x = WrapUnique(new X(1, 2));  // works, but nonideal.
//   auto x = MakeUnique<X>(1, 2);      // safer, standard, avoids raw 'new'.
//
// Note that `absl::WrapUnique(p)` is valid only if `delete p` is a valid
// expression. In particular, `absl::WrapUnique()` cannot wrap pointers to
// arrays, functions or void, and it must not be used to capture pointers
// obtained from array-new expressions (even though that would compile!).
template <typename T>
std::unique_ptr<T> WrapUnique(T* ptr) {
  static_assert(!std::is_array<T>::value, "array types are unsupported");
  static_assert(std::is_object<T>::value, "non-object types are unsupported");
  return std::unique_ptr<T>(ptr);
}

namespace internal {

// Traits to select proper overload and return type for `absl::MakeUnique<>`.
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

// -----------------------------------------------------------------------------
// Function Template: MakeUnique<T>() / make_unique<T>()
// -----------------------------------------------------------------------------
//
// Creates a `std::unique_ptr<>`, while avoiding issues creating temporaries
// during the construction process. `absl::MakeUnique<>` also avoids redundant
// type declarations, by avoiding the need to explicitly use the `new` operator.
//
// This implementation of `absl::MakeUnique<>` is designed for C++11 code and
// will be replaced in C++14 by the equivalent `std::make_unique<>` abstraction.
// `absl::MakeUnique<>` is designed  to be 100% compatible with
// `std::make_unique<>` so that the eventual migration will involve a simple
// rename operation.
//
// For more background on why `std::unique_ptr<T>(new T(a,b))` is problematic,
// see Herb Sutter's explanation on
// (Exception-Safe Function Calls)[http://herbsutter.com/gotw/_102/].
// (In general, reviewers should treat `new T(a,b)` with scrutiny.)
//
// Example usage:
//
//    auto p = MakeUnique<X>(args...);  // 'p'  is a std::unique_ptr<X>
//    auto pa = MakeUnique<X[]>(5);     // 'pa' is a std::unique_ptr<X[]>
//
// Three overloads of `absl::MakeUnique` are required:
//
//   - For non-array T:
//
//       Allocates a T with `new T(std::forward<Args> args...)`,
//       forwarding all `args` to T's constructor.
//       Returns a `std::unique_ptr<T>` owning that object.
//
//   - For an array of unknown bounds T[]:
//
//       `absl::MakeUnique<>` will allocate an array T of type U[] with
//       `new U[n]()` and return a `std::unique_ptr<U[]>` owning that array.
//
//       Note that 'U[n]()' is different from 'U[n]', and elements will be
//       value-initialized. Note as well that `std::unique_ptr` will perform its
//       own destruction of the array elements upon leaving scope, even though
//       the array [] does not have a default destructor.
//
//       NOTE: an array of unknown bounds T[] may still be (and often will be)
//       initialized to have a size, and will still use this overload. E.g:
//
//         auto my_array = absl::MakeUnique<int[]>(10);
//
//   - For an array of known bounds T[N]:
//
//       `absl::MakeUnique<>` is deleted (like with `std::make_unique<>`) as
//       this overload is not useful.
//
//       NOTE: an array of known bounds T[N] is not considered a useful
//       construction, and may cause undefined behavior in templates. E.g:
//
//         auto my_array = absl::MakeUnique<int[10]>();
//
//       In those cases, of course, you can still use the overload above and
//       simply intialize it to its desired size:
//
//         auto my_array = absl::MakeUnique<int[]>(10);

// `absl::MakeUnique` overload for non-array types.
template <typename T, typename... Args>
typename internal::MakeUniqueResult<T>::scalar
MakeUnique(Args&&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

// `absl::MakeUnique` overload for an array T[] of unknown bounds.
// The array allocation needs to use the `new T[size]` form and cannot take
// element constructor arguments. The `std::unique_ptr` will manage destructing
// these array elements.
template <typename T>
typename internal::MakeUniqueResult<T>::array
MakeUnique(size_t n) {
  return std::unique_ptr<T>(new typename absl::remove_extent_t<T>[n]());
}

// `absl::MakeUnique` overload for an array T[N] of known bounds.
// This construction will be rejected.
template <typename T, typename... Args>
typename internal::MakeUniqueResult<T>::invalid
MakeUnique(Args&&... /* args */) = delete;

// `absl::make_unique` overload for non-array types.
template <typename T, typename... Args>
typename internal::MakeUniqueResult<T>::scalar
make_unique(Args&&... args) {
  return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

// `absl::make_unique` overload for an array T[] of unknown bounds.
// The array allocation needs to use the `new T[size]` form and cannot take
// element constructor arguments. The `std::unique_ptr` will manage destructing
// these array elements.
template <typename T>
typename internal::MakeUniqueResult<T>::array
make_unique(size_t n) {
  return std::unique_ptr<T>(new typename absl::remove_extent_t<T>[n]());
}

// `absl::make_unique` overload for an array T[N] of known bounds.
// This construction will be rejected.
template <typename T, typename... Args>
typename internal::MakeUniqueResult<T>::invalid
make_unique(Args&&... /* args */) = delete;

// -----------------------------------------------------------------------------
// Function Template: RawPtr()
// -----------------------------------------------------------------------------
//
// Extracts the raw pointer from a pointer-like 'ptr'. `absl::RawPtr` is useful
// within templates that need to handle a complement of raw pointers,
// `std::nullptr_t`, and smart pointers.
template <typename T>
auto RawPtr(T&& ptr) -> decltype(&*ptr) {  // NOLINT
  // ptr is a forwarding reference to support Ts with non-const operators.
  return (ptr != nullptr) ? &*ptr : nullptr;
}
inline std::nullptr_t RawPtr(std::nullptr_t) { return nullptr; }

// -----------------------------------------------------------------------------
// Function Template: ShareUniquePtr()
// -----------------------------------------------------------------------------
//
// Transforms a `std::unique_ptr` rvalue into a `std::shared_ptr`. The returned
// value is a `std::shared_ptr` of deduced type and ownership is transferred to
// the shared pointer.
//
// Example:
//
//     auto up = absl::MakeUnique<int>(10);
//     auto sp = absl::ShareUniquePtr(std::move(up));  // shared_ptr<int>
//     CHECK_EQ(*sp, 10);
//     CHECK(up == nullptr);
//
// Note that this conversion is correct even when T is an array type, although
// the resulting shared pointer may not be very useful.
//
// Implements the resolution of [LWG 2415](http://wg21.link/lwg2415), by which a
// null shared pointer does not attempt to call the deleter.
template <typename T, typename D>
std::shared_ptr<T> ShareUniquePtr(std::unique_ptr<T, D>&& ptr) {  // NOLINT
  return ptr ? std::shared_ptr<T>(std::move(ptr)) : std::shared_ptr<T>();
}

// -----------------------------------------------------------------------------
// Function Template: WeakenPtr()
// -----------------------------------------------------------------------------
//
// Creates a weak pointer associated with a given shared pointer. The returned
// value is a `std::weak_ptr` of deduced type.
//
// Example:
//
//    auto sp = std::make_shared<int>(10);
//    auto wp = absl::WeakenPtr(sp);
//    CHECK_EQ(sp.get(), wp.lock().get());
//    sp.reset();
//    CHECK(wp.lock() == nullptr);
//
template <typename T>
std::weak_ptr<T> WeakenPtr(const std::shared_ptr<T>& ptr) {
  return std::weak_ptr<T>(ptr);
}

}  // namespace absl

// TODO(user): Delete temporary aliases after namespace update.
namespace gtl {

template <typename T>
std::unique_ptr<T> WrapUnique(T* ptr) {
  return absl::WrapUnique(ptr);
}

using absl::MakeUnique;

}  // namespace gtl

#endif  // S2_THIRD_PARTY_ABSL_MEMORY_MEMORY_H_
