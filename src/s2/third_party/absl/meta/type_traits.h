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

// This file contains local extensions to the standard <type_traits> header,
// including extensions which are in existing or upcoming standards, but aren't
//
// Utilities should be added to this header only if they are clearly on
// track for standardization (e.g. if they've been voted into the working
// paper), or otherwise well-baked and unlikely to change substantially.
//
// WARNING: use of many of the constructs in this header will count as "complex
// template metaprogramming", so before proceeding, please carefully consider

#ifndef S2_THIRD_PARTY_ABSL_META_TYPE_TRAITS_H_
#define S2_THIRD_PARTY_ABSL_META_TYPE_TRAITS_H_

#include <type_traits>

#include "s2/third_party/absl/base/config.h"

namespace absl {

namespace internal_type_traits {
template <typename... Ts>
struct VoidTImpl {
  using type = void;
};
}  // namespace internal_type_traits

// This is like C++17's std::void_t
// (http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2014/n3911.pdf).
// It's a dummy metafunction that ignores its arguments and returns `void`.
//
// We do not use the standard-specified implementation in order to preserve
// compatibility with gcc < 5.1. Unfortunately, our workaround provides
// slightly different behavior in some circumstances, such as when ordering
// partial specializations.
template <typename... Ts>
using void_t = typename internal_type_traits::VoidTImpl<Ts...>::type;

// conjunction, disjunction, and negation are equivalent to the constructs with
// the same names in C++17
// (http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2015/p0013r1.html).
// They implement logical AND, OR, and NOT operations as compile-time
// metafunctions. Their template arguments must have `value` members that
// are convertible to bool. However, conjunction and disjunction have
// short-circuiting behavior: argument evaluation stops once the result is
// known, so the following arguments need not have the `value` member.
template <typename... Ts>
struct conjunction;

template <typename T, typename... Ts>
struct conjunction<T, Ts...>
    : std::conditional<T::value, conjunction<Ts...>, T>::type {};

template <typename T>
struct conjunction<T> : T {};

template <>
struct conjunction<> : std::true_type {};

template <typename... Ts>
struct disjunction;

template <typename T, typename... Ts>
struct disjunction<T, Ts...> :
      std::conditional<T::value, T, disjunction<Ts...>>::type {};

template <typename T>
struct disjunction<T> : T {};

template <>
struct disjunction<> : std::false_type {};

template <typename T>
struct negation : std::integral_constant<bool, !T::value> {};

// The following is_trivially traits are defined because they are not fully
// supported in libstdc++ 4.x.
// The extensions (__has_trivial_xxx) are implemented in gcc (version >= 4.3)
// and clang. Since we are supporting libstdc++ > 4.7, they should always be
// present.
// The extensions are documented at
// https://gcc.gnu.org/onlinedocs/gcc/Type-Traits.html#Type-Traits.
// We add static assertion to detect any incompliance with std:: traits (if
// they are implemented by compiler and standard library).
template <typename T>
struct is_trivially_destructible
    : std::integral_constant<bool, __has_trivial_destructor(T) &&
                                       std::is_destructible<T>::value> {
#ifdef GOOGLE_HAVE_STD_IS_TRIVIALLY_DESTRUCTIBLE
  static_assert(std::is_trivially_destructible<T>::value ==
                    is_trivially_destructible::value,
                "Not compliant with std::is_trivially_destructible");
#endif  // GOOGLE_HAVE_STD_IS_TRIVIALLY_DESTRUCTIBLE
};

// According to C++ standard, Section: 20.15.4.3 [meta.unary.prop]
// The predicate condition for a template specialization is_constructible<T,
// Args...> shall be satisfied if and only if the following variable
// definition would be well-formed for some invented variable t:
//
// T t(declval<Args>()...);
//
// is_trivially_constructible<T, Args...> additionally requires that the
// variable definition does not call any operation that is not trivial.
// For the purposes of this check, the call to std::declval is considered
// trivial.
//
// Notes from http://en.cppreference.com/w/cpp/types/is_constructible:
// In many implementations, is_nothrow_constructible also checks if the
// destructor throws because it is effectively noexcept(T(arg)). Same
// applies to is_trivially_constructible, which, in these implementations, also
// requires that the destructor is trivial.
// GCC bug 51452: https://gcc.gnu.org/bugzilla/show_bug.cgi?id=51452
// LWG issue 2116: http://cplusplus.github.io/LWG/lwg-active.html#2116.

// "T obj();" need to be well-formed and not call any non-trivial operation.
// Nontrivally destructible types will cause the expression to be nontrivial.
template <typename T>
struct is_trivially_default_constructible
    : std::integral_constant<bool,
                             __has_trivial_constructor(T) &&
                                 std::is_default_constructible<T>::value &&
                                 is_trivially_destructible<T>::value> {
#ifdef GOOGLE_HAVE_STD_IS_TRIVIALLY_CONSTRUCTIBLE
  static_assert(std::is_trivially_default_constructible<T>::value ==
                    is_trivially_default_constructible::value,
                "Not compliant with std::is_trivially_default_constructible");
#endif  // GOOGLE_HAVE_STD_IS_TRIVIALLY_CONSTRUCTIBLE
};

// "T obj(declval<const T&>());" need to be well-formed and not call any
// nontrivial operation.
// Nontrivally destructible types will cause the expression to be nontrivial.
template <typename T>
struct is_trivially_copy_constructible
    : std::integral_constant<bool, __has_trivial_copy(T) &&
                                       std::is_copy_constructible<T>::value &&
                                       is_trivially_destructible<T>::value> {
#ifdef GOOGLE_HAVE_STD_IS_TRIVIALLY_CONSTRUCTIBLE
  static_assert(std::is_trivially_copy_constructible<T>::value ==
                    is_trivially_copy_constructible::value,
                "Not compliant with std::is_trivially_copy_constructible");
#endif  // GOOGLE_HAVE_STD_IS_TRIVIALLY_CONSTRUCTIBLE
};

// is_assignable<T, U>::value is true if the expression "declval<T>() =
// declval<U>()" is well-formed when treated as an unevaluated operand.
// is_trivially_assignable<T, U> requires the assignment to call no operation
// that is not trivial.
// is_trivially_copy_assignable<T> is just is_trivially_assignable<T, const T&>.
template <typename T>
struct is_trivially_copy_assignable
    : std::integral_constant<bool, __has_trivial_assign(T) &&
                                       std::is_copy_assignable<T>::value> {
#ifdef GOOGLE_HAVE_STD_IS_TRIVIALLY_ASSIGNABLE
  static_assert(std::is_trivially_copy_assignable<T>::value ==
                    is_trivially_copy_assignable::value,
                "Not compliant with std::is_trivially_copy_assignable");
#endif  // GOOGLE_HAVE_STD_IS_TRIVIALLY_ASSIGNABLE
};

}  // namespace absl

// Temporary aliases while moving into the absl namespace.
// TODO(user): Delete temporary aliases after namespace update.
namespace gtl {
template <typename... Ts>
using void_t = typename absl::void_t<Ts...>;

using absl::conjunction;
using absl::disjunction;
using absl::negation;
using absl::is_trivially_destructible;
using absl::is_trivially_default_constructible;
using absl::is_trivially_copy_constructible;
using absl::is_trivially_copy_assignable;
}  // namespace gtl

#endif  // S2_THIRD_PARTY_ABSL_META_TYPE_TRAITS_H_
