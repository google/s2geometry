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

//
// Copyright 2017 The Abseil Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// -----------------------------------------------------------------------------
// span.h
// -----------------------------------------------------------------------------
//
// This header file defines a `Span<T>` type for holding a view of an existing
// array of data. The `Span` object, much like the `absl::string_view` object,
// does not own such data itself. A span provides a lightweight way to pass
// around view of such data.
//
// Additionally, this header file defines `MakeSpan()` and `MakeConstSpan()`
// factory functions, for clearly creating spans of type `Span<T>` or read-only
// `Span<const T>` when such types may be difficult to identify due to issues
// with implicit conversion.
//
// The C++ standards committee currently has a proposal for a `std::span` type,
// (http://wg21.link/p0122), which is not yet part of the standard (though may
// become part of C++20). As of August 2017, the differences between
// `absl::Span` and this proposal are:
//    * `absl::Span` uses `std::size_t` for `size_type`
//    * `absl::Span` has no `operator()`
//    * `absl::Span` has no constructors for `std::unique_ptr` or
//      `std::shared_ptr`
//    * `absl::span` has the factory functions `MakeSpan()` and
//      `MakeConstSpan()`
//    * `absl::Span` has `front()` and `back()` methods
//    * bounds-checked access to `absl::Span` is accomplished with `at()`
//    * `absl::Span` has compiler-provided move and copy constructors and
//      assignment. This is due to them being specified as `constexpr`, but that
//      implies const in C++11.
//    * `absl::Span` has no `element_type` or `index_type` typedefs
//    * A read-only `absl::Span<const T>` can be implicitly constructed from an
//      initializer list.
//    * `absl::Span` has no `bytes()`, `size_bytes()`, `as_bytes()`, or
//      `as_mutable_bytes()` methods
//    * `absl::Span` has no static extent template parameter, nor constructors
//      which exist only because of the static extent parameter.
//    * `absl::Span` has an explicit mutable-reference constructor
//
// For more information, see the class comments below.
#ifndef S2_THIRD_PARTY_ABSL_TYPES_SPAN_H_
#define S2_THIRD_PARTY_ABSL_TYPES_SPAN_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <iterator>
#include <string>
#include <type_traits>
#include <utility>

#include "s2/third_party/absl/base/gdb_scripting.h"
#include "s2/third_party/absl/algorithm/algorithm.h"
#include "s2/third_party/absl/base/macros.h"
#include "s2/third_party/absl/base/port.h"
#include "s2/third_party/absl/meta/type_traits.h"

namespace absl {

template <typename T>
class Span;

namespace internal_span {
// A constexpr min function
constexpr std::size_t Min(std::size_t a, std::size_t b) noexcept {
  return a < b ? a : b;
}

// Wrappers for access to container data pointers.
template <typename C>
constexpr auto GetDataImpl(C& c, char) noexcept  // NOLINT(runtime/references)
    -> decltype(c.data()) {
  return c.data();
}
template <typename C>
constexpr auto GetDataImpl(C& c, int) noexcept  // NOLINT(runtime/references)
    -> decltype(c.mutable_data()) {
  return c.mutable_data();
}

// Before C++17, string::data returns a const char* in all cases.
inline char* GetDataImpl(string& s,
                         int) noexcept {  // NOLINT(runtime/references)
  return &s[0];
}

template <typename C>
constexpr auto GetData(C& c) noexcept
    -> decltype(GetDataImpl(c, 0)) {  // NOLINT(runtime/references)
  return GetDataImpl(c, 0);
}

// Detection idioms for size() and data().
template <typename C>
using HasSize =
    std::is_integral<absl::decay_t<decltype(std::declval<C&>().size())>>;

// We want to enable conversion from vector<T*> to Span<const T* const> but
// disable conversion from vector<Derived> to Span<Base>. Here we use
// the fact that U** is convertible to Q* const* if and only if Q is the same
// type or a more cv-qualified version of U.  We also decay the result type of
// data() to avoid problems with classes which have a member function data()
// which returns a reference.
template <typename T, typename C>
using HasData =
    std::is_convertible<absl::decay_t<decltype(GetData(std::declval<C&>()))>*,
                        T* const*>;

// Extracts value type from a Container
template <typename C>
struct ElementType {
  using type = typename absl::remove_reference_t<C>::value_type;
};

template <typename T, size_t N>
struct ElementType<T (&)[N]> {
  using type = T;
};

template <typename C>
using ElementT = typename ElementType<C>::type;

template <typename T>
using EnableIfMutable = absl::enable_if_t<!std::is_const<T>::value, int>;

// Used in defining relationals for Span.  Identity<T> in an argument list puts
// T in a non-deduced context.
template <typename SpanT>
using Identity = absl::decay_t<SpanT>;

template <typename T>
bool EqualImpl(Span<T> a, Span<T> b) {
  static_assert(std::is_const<T>::value, "");
  return absl::equal(a.begin(), a.end(), b.begin(), b.end());
}

template <typename T>
bool LessThanImpl(Span<T> a, Span<T> b) {
  static_assert(std::is_const<T>::value, "");
  return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}
}  // namespace internal_span

//------------------------------------------------------------------------------
// Span
//------------------------------------------------------------------------------
//
// A `Span` is an "array view" type for holding a view of a contiguous data
// array; the `Span` object does not and cannot own such data itself. A span
// provides an easy way to provide overloads for anything operating on
// contiguous sequences without needing to manage pointers and array lengths
// manually.

// A span is conceptually a pointer (ptr) and a length (size) into an already
// existing array of contiguous memory; the array it represents references the
// elements "ptr[0] .. ptr[size-1]". Passing a properly-constructed `Span`
// instead of raw pointers avoids many issues related to index out of bounds
// errors.
//
// Spans may also be constructed from containers holding contiguous sequences.
// Such containers must supply `data()` and `size() const` methods (e.g
// `std::vector<T>`, `absl::InlinedVector<T, N>`). All implicit conversions to
// `absl::Span` from such containers will create spans of type `const T`;
// spans which can mutate their values (of type `T`) must use explicit
// constructors.
//
// A `Span<T>` is somewhat analogous to an `absl::string_view`, but for an array
// of elements of type `T`. A user of `Span` must ensure that the data being
// pointed to outlives the `Span` itself.
//
// You can construct a `Span<T>` in several ways:
//
//   * Explicitly from a reference to a container type
//   * Explicitly from a pointer and size
//   * Implicitly from a container type (but only for spans of type `const T`)
//   * Using the `MakeSpan()` or `MakeConstSpan()` factory functions.
//
// Examples:
//
//   // Construct a Span explicitly from a container:
//   std::vector<int> v = {1, 2, 3, 4, 5};
//   auto span = absl::Span<const int>(v);
//
//   // Construct a Span explicitly from a C-style array:
//   int a[5] =  {1, 2, 3, 4, 5};
//   auto span = absl::Span<const int>(a);
//
//   // Construct a Span implicitly from a container
//   void MyRoutine(absl::Span<const int> a) {
//     ...
//   };
//   std::vector v = {1,2,3,4,5};
//   MyRoutine(v)                     // convert to Span<const T>
//
// Note that `Span` objects, in addition to requiring that the memory they
// point to remains alive, must also ensure that such memory does not get
// reallocated. Therefore, to avoid undefined behavior, containers with
// associated span views should not invoke operations that may reallocate memory
// (such as resizing) or invalidate iterarors into the container.
//
// One common use for a `Span` is when passing arguments to a routine that can
// accept a variety of array types (e.g. a `std::vector`, `absl::InlinedVector`,
// a C-style array, etc.). Instead of creating overloads for each case, you
// can simply specify a `Span` as the argument to such a routine.
//
// Example:
//
//   void MyRoutine(absl::Span<const int> a) {
//     ...
//   };
//
//   std::vector v = {1,2,3,4,5};
//   MyRoutine(v);
//
//   absl::InlinedVector<int, 4> my_inline_vector;
//   MyRoutine(my_inline_vector);
//
//   // Explicit constructor from pointer,size
//   int my_array = new int[10];
//   MyRoutine(absl::Span<const int>(my_array, 10));
template <typename T>
class Span {
 private:
  // Used to determine whether a Span can be constructed from a container of
  // type C.
  template <typename C>
  using EnableIfConvertibleFrom =
      absl::enable_if_t<internal_span::HasData<T, C>::value &&
                        internal_span::HasSize<C>::value>;

  // Used to SFINAE-enable a function when the slice elements are const.
  template <typename U>
  using EnableIfConstView = absl::enable_if_t<std::is_const<T>::value, U>;

  // Used to SFINAE-enable a function when the slice elements are mutable.
  template <typename U>
  using EnableIfMutableView = absl::enable_if_t<!std::is_const<T>::value, U>;

 public:
  using value_type = absl::remove_cv_t<T>;
  using pointer = T*;
  using const_pointer = const T*;
  using reference = T&;
  using const_reference = const T&;
  using iterator = pointer;
  using const_iterator = const_pointer;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using size_type = size_t;
  using difference_type = ptrdiff_t;

  static const size_type npos = -1;

  constexpr Span() noexcept : Span(nullptr, 0) {}
  constexpr Span(pointer array, size_type length) noexcept
      : ptr_(array), len_(length) {}

  // Implicit conversion constructors
  template <size_t N>
  constexpr Span(T (&a)[N]) noexcept  // NOLINT(runtime/explicit)
      : Span(a, N) {}

  // Substring of another Span.
  // pos must be non-negative and <= x.length().
  // len must be non-negative and will be pinned to at most x.length() - pos.
  // If len==npos, the substring continues till the end of x.
  ABSL_DEPRECATED("Prefer Span(...).subspan(pos, len)")
  constexpr Span(Span x, size_type pos, size_type len) noexcept
      : ptr_(x.ptr_ + pos), len_(internal_span::Min(x.len_ - pos, len)) {}

  // The constructor for any class supplying 'T* data()' or 'T* mutable_data()'
  // (the former is called if both exist), and 'some_integral_type size()
  // const'. proto2::RepeatedField is an example of this. Also supports string
  // arguments, when T==char. The appropriate ctor is selected using SFINAE.
  template <typename V, typename = EnableIfConvertibleFrom<V>,
            typename = EnableIfMutableView<V>>
  ABSL_DEPRECATED("Use the explicit mutable ref conversion ctor")
  Span(V* v)  // NOLINT(runtime/explicit)
      : Span(internal_span::GetData(*v), v->size()) {}

  // Explicit reference constructor for a mutable `Span<T>` type
  template <typename V, typename = EnableIfConvertibleFrom<V>,
            typename = EnableIfMutableView<V>>
  explicit Span(V& v) noexcept  // NOLINT(runtime/references)
      : Span(internal_span::GetData(v), v.size()) {}

  // Implicit reference constructor for a read-only `Span<const T>` type
  template <typename V, typename = EnableIfConvertibleFrom<V>,
            typename = EnableIfConstView<V>>
  constexpr Span(const V& v) noexcept  // NOLINT(runtime/explicit)
      : Span(internal_span::GetData(v), v.size()) {}

  // Implicit constructor from an initializer list, making it possible to pass a
  // brace-enclosed initializer list to a function expecting a `Span`. Such
  // spans constructed from an initializer list must be of type `Span<const T>`.
  //
  //   void Process(absl::Span<const int> x);
  //   Process({1, 2, 3});
  //
  // Note: the data referenced by the initializer list cannot be lvalues as such
  // data must outlive this `Span`.
  //
  // Example:
  //    Span<const int> s={1,2};       // Error. lvalues won't outlive Span
  //    return Span<const int>({3,4}); // Error. lvalues won't outlive Span
  template <typename LazyT = T,
            typename = EnableIfConstView<LazyT>>
  constexpr Span(
      std::initializer_list<value_type> v) noexcept  // NOLINT(runtime/explicit)
      : Span(v.begin(), v.size()) {}

  // Accessors

  // Span::data()
  //
  // Returns a pointer to the span's underlying array of data (which is held
  // outside the span).
  constexpr pointer data() const noexcept { return ptr_; }

  // Span::size()
  //
  // Returns the size of this span.
  constexpr size_type size() const noexcept { return len_; }

  // Span::length()
  //
  // Returns the length (size) of this span.
  constexpr size_type length() const noexcept { return size(); }

  // Span::empty()
  //
  // Returns a boolean indicating whether or not this span is considered empty.
  constexpr bool empty() const noexcept { return size() == 0; }

  // Span::operator[]
  //
  // Returns a reference to the i'th element of this span.
  constexpr reference operator[](size_type i) const noexcept { return ptr_[i]; }

  // Span::at()
  //
  // Returns a reference to the i'th element of this span.
  constexpr reference at(size_type i) const {
    return ABSL_ASSERT(i < size()), ptr_[i];
  }

  // Span::front()
  //
  // Returns a reference to the first element of this span.
  constexpr reference front() const noexcept {
    return ABSL_ASSERT(size() > 0), ptr_[0];
  }

  // Span::back()
  //
  // Returns a reference to the last element of this span.
  constexpr reference back() const noexcept {
    return ABSL_ASSERT(size() > 0), ptr_[size() - 1];
  }

  // Span::begin()
  //
  // Returns an iterator to the first element of this span.
  constexpr iterator begin() const noexcept { return ptr_; }

  // Span::cbegin()
  //
  // Returns a const iterator to the first element of this span.
  constexpr const_iterator cbegin() const noexcept { return ptr_; }

  // Span::end()
  //
  // Returns an iterator to the last element of this span.
  constexpr iterator end() const noexcept { return ptr_ + len_; }

  // Span::cend()
  //
  // Returns a const iterator to the last element of this span.
  constexpr const_iterator cend() const noexcept { return end(); }

  // Span::rbegin()
  //
  // Returns a reverse iterator starting at the last element of this span.
  reverse_iterator rbegin() const noexcept { return reverse_iterator(end()); }

  // Span::crbegin()
  //
  // Returns a reverse const iterator starting at the last element of this span.
  const_reverse_iterator crbegin() const noexcept { return rbegin(); }

  // Span::rend()
  //
  // Returns a reverse iterator starting at the first element of this span.
  reverse_iterator rend() const noexcept { return reverse_iterator(begin()); }

  // Span::crend()
  //
  // Returns a reverse iterator starting at the first element of this span.
  const_reverse_iterator crend() const noexcept { return rend(); }

  // Span mutations

  // Span::remove_prefix()
  //
  // Removes the first `n` elements from the span.
  void remove_prefix(size_type n) noexcept {
    assert(len_ >= n);
    ptr_ += n;
    len_ -= n;
  }

  // Span::remove_suffix()
  //
  // Removes the last `n` elements from the span.
  void remove_suffix(size_type n) noexcept {
    assert(len_ >= n);
    len_ -= n;
  }

  ABSL_DEPRECATED("Use remove_suffix(1) instead.")
  void pop_back() { remove_suffix(1); }
  ABSL_DEPRECATED("Use remove_prefix(1) instead.")
  void pop_front() { remove_prefix(1); }

  // Span::subspan()
  //
  // Returns a `Span` starting at element `pos` and of length `len`, with
  // proper bounds checking to ensure `len` does not exceed the ptr+size of the
  // original array. (Spans whose `len` would point past the end of the array
  // will throw a `std::out_of_range`.)
  constexpr Span subspan(size_type pos = 0, size_type len = npos) const
      noexcept {
    return Span(*this, pos, len);
  }

 private:
  pointer ptr_;
  size_type len_;
};

template <typename T>
const typename Span<T>::size_type Span<T>::npos;

// Span relationals

// Equality is compared element-by-element, while ordering is lexicographical.
// We provide three overloads for each operator to cover any combination on the
// left or right hand side of mutable Span<T>, read-only Span<const T>, and
// convertible-to-read-only Span<T>.

// operator==
template <typename T>
constexpr bool operator==(Span<T> a, Span<T> b) {
  return internal_span::EqualImpl<const T>(a, b);
}
template <typename T>
constexpr bool operator==(internal_span::Identity<Span<const T>> a, Span<T> b) {
  return internal_span::EqualImpl<const T>(a, b);
}
template <typename T>
constexpr bool operator==(Span<T> a, internal_span::Identity<Span<const T>> b) {
  return internal_span::EqualImpl<const T>(a, b);
}

// operator!=
template <typename T>
constexpr bool operator!=(Span<T> a, Span<T> b) {
  return !(a == b);
}
template <typename T>
constexpr bool operator!=(internal_span::Identity<Span<const T>> a, Span<T> b) {
  return !(a == b);
}
template <typename T>
constexpr bool operator!=(Span<T> a, internal_span::Identity<Span<const T>> b) {
  return !(a == b);
}

// operator<
template <typename T>
constexpr bool operator<(Span<T> a, Span<T> b) {
  return internal_span::LessThanImpl<const T>(a, b);
}
template <typename T>
constexpr bool operator<(internal_span::Identity<Span<const T>> a, Span<T> b) {
  return internal_span::LessThanImpl<const T>(a, b);
}
template <typename T>
constexpr bool operator<(Span<T> a, internal_span::Identity<Span<const T>> b) {
  return internal_span::LessThanImpl<const T>(a, b);
}

// operator>
template <typename T>
constexpr bool operator>(Span<T> a, Span<T> b) {
  return b < a;
}
template <typename T>
constexpr bool operator>(internal_span::Identity<Span<const T>> a, Span<T> b) {
  return b < a;
}
template <typename T>
constexpr bool operator>(Span<T> a, internal_span::Identity<Span<const T>> b) {
  return b < a;
}

// operator<=
template <typename T>
constexpr bool operator<=(Span<T> a, Span<T> b) {
  return !(b < a);
}
template <typename T>
constexpr bool operator<=(internal_span::Identity<Span<const T>> a, Span<T> b) {
  return !(b < a);
}
template <typename T>
constexpr bool operator<=(Span<T> a, internal_span::Identity<Span<const T>> b) {
  return !(b < a);
}

// operator>=
template <typename T>
constexpr bool operator>=(Span<T> a, Span<T> b) {
  return !(a < b);
}
template <typename T>
constexpr bool operator>=(internal_span::Identity<Span<const T>> a, Span<T> b) {
  return !(a < b);
}
template <typename T>
constexpr bool operator>=(Span<T> a, internal_span::Identity<Span<const T>> b) {
  return !(a < b);
}

// MakeSpan()
//
// Constructs a mutable `Span<T>`, deducing `T` automatically from either a
// container or pointer+size.
//
// Because a read-only `Span<const T>` is implicitly constructed from container
// types regardless of whether the container itself is a const container,
// constructing mutable spans of type `Span<T>` from containers requires
// explicit constructors. The container-accepting version of `MakeSpan()`
// deduces the type of `T` by the constness of the pointer received from the
// container's `data()` member. Similarly, the pointer-accepting version returns
// a `Span<const T>` if `T` is `const`, and a `Span<T>` otherwise.
//
// Examples:
//
//   void MyRoutine(absl::Span<MyComplicatedType> a) {
//     ...
//   };
//   // my_vector is a container of non-const types
//   std::vector<MyComplicatedType> my_vector;
//
//   // Constructing a Span implicitly attempts to create a Span of type
//   // `Span<const T>`
//   MyRoutine(my_vector);                // error, type mismatch
//
//   // Explicitly constructing the Span is verbose
//   MyRoutine(absl::Span<MyComplicatedType>(my_vector);
//
//   // Use MakeSpan() to make an absl::Span<T>
//   MyRoutine(absl::MakeSpan(my_vector));
//
//   // Construct a span from an array ptr+size
//   absl::Span<T> my_span() {
//     return absl::MakeSpan(&array[0], num_elements_);
//   }
template <int&... ExplicitArgumentBarrier, typename T>
constexpr Span<T> MakeSpan(T* ptr, size_t size) noexcept {
  return Span<T>(ptr, size);
}

template <int&... ExplicitArgumentBarrier, typename T>
constexpr Span<T> MakeSpan(T* begin, T* end) noexcept {
  return ABSL_ASSERT(begin <= end), Span<T>(begin, end - begin);
}

template <int&... ExplicitArgumentBarrier, typename C>
constexpr auto MakeSpan(C& c) noexcept  // NOLINT(runtime/references)
    -> decltype(absl::MakeSpan(internal_span::GetData(c), c.size())) {
  return MakeSpan(internal_span::GetData(c), c.size());
}

template <int&... ExplicitArgumentBarrier, typename T, size_t N>
constexpr Span<T> MakeSpan(T (&array)[N]) noexcept {
  return Span<T>(array, N);
}

// MakeConstSpan()
//
// Constructs a `Span<const T>` as with `MakeSpan`, deducing `T` automatically,
// but always returning a `Span<const T>`.
//
// Example:
//
//   absl::Span<const T> my_span() {
//     return absl::MakeConstSpan(&array_[0], num_elements_);
//   }
//
// MakeSpan differs from gtl::MakeSlice in the following ways:
//   * Automatically deduces the constness of T instead of always returning
//   Span<const T>. Use absl::MakeConstSpan to get this behavior.  In the future
//   migration away from gtl::MakeSlice, existing callers of MakeSlice will be
//   migrated to absl::MakeConstSpan, while existing callers of MakeMutableSlice
//   will be migrated to absl::MakeSpan.
//
//   * Accepts an lvalue reference instead of a forwarding reference, preventing
//   the use of MakeSpan on temporaries.
template <int&... ExplicitArgumentBarrier, typename T>
constexpr Span<const T> MakeConstSpan(T* ptr, size_t size) noexcept {
  return Span<const T>(ptr, size);
}

template <int&... ExplicitArgumentBarrier, typename T>
constexpr Span<const T> MakeConstSpan(T* begin, T* end) noexcept {
  return ABSL_ASSERT(begin <= end), Span<const T>(begin, end - begin);
}

template <int&... ExplicitArgumentBarrier, typename C>
constexpr auto MakeConstSpan(const C& c) noexcept -> decltype(MakeSpan(c)) {
  return MakeSpan(c);
}

template <int&... ExplicitArgumentBarrier, typename T, size_t N>
constexpr Span<const T> MakeConstSpan(const T (&array)[N]) noexcept {
  return Span<const T>(array, N);
}

}  // namespace absl
namespace gtl {
namespace internal_span = absl::internal_span;

// NOTE: gtl::MutableArraySlice is deprecated, absl::Span should be used
// instead.
//
// MutableArraySlice<T> represents a mutable array of elements, and, like
// Span, does not own the backing store. The implicit constructors it
// provides allow functions not to worry about whether their mutable arguments
// refer to vectors, arrays, proto2::RepeatedFields, etc.:
//
//   void MyMutatingRoutine(gtl::MutableArraySlice<int> a) {
//     for (int i = 0; i < a.size(); i++) { .. mutate a[i] .. }
//   }
//
//   std::vector<int> my_vector;
//   MyMutatingRoutine(&my_vector);
//
//   int my_array[10];
//   MyMutatingRoutine(my_array);
//
//   int* my_array = new int[10];
//   MyMutatingRoutine(gtl::MutableArraySlice<int>(my_array, 10));
//
//   MyProto my_proto;
//   for (int i = 0; i < 10; ++i) { my_proto.add_value(i); }
//   MyMutatingRoutine(my_proto.mutable_value());
//
// There are also factory functions which allow deduction of types of produced
// slices. For example:
//
//    std::vector<int> my_vector(0, 10);
//    auto my_mutable_slice = MakeMutableArraySlice(&my_vector);
//
//    static const int kMyArray[] = {1, 2, 3, 4, 5};
//    auto boxed_array = MakeArraySlice(kMyArray);
//    auto partial_array = MakeArraySlice(kMyArray, 3);
template <typename T>
using ArraySlice = absl::Span<const T>;

// Mutable version of ArraySlice, which allows the clients to mutate the
// underlying data. It is implicitly convertible to ArraySlice since it provides
// the data() and size() methods with correct signatures. When a
// MutableArraySlice is created from a pointer to a container (as opposed to raw
// memory pointer), the pointer must not be null.
//
// A note on const-ness: "mutable" here refers to the mutability of the
// underlying data, not of the slice itself. It is perfectly reasonable to have
// a variable of type "const MutableArraySlice<T>"; this means that the bounds
// of the view on the array cannot be changed, but the underlying data in the
// array still may be modified. This is akin to a "T* const" pointer, as opposed
// to a "const T*" pointer (corresponding to a non-const ArraySlice<T>).
template <typename T>
using MutableArraySlice = absl::Span<T>;

// Returns an ArraySlice referring to a container or array. Does not participate
// in overload resolution unless a Span can be constructed from 'c'.
template <int&... ExplicitArgumentBarrier, typename C>
auto MakeArraySlice(C&& c)
    -> decltype(absl::MakeConstSpan(std::forward<C>(c))) {
  return absl::MakeConstSpan(std::forward<C>(c));
}

template <int&... ExplicitArgumentBarrier, typename T>
ArraySlice<T> MakeArraySlice(T* ptr, size_t len) {
  return absl::MakeConstSpan(ptr, len);
}

template <int&... ExplicitArgumentBarrier, typename T>
ArraySlice<T> MakeArraySlice(T* begin, T* end) {
  assert(begin <= end);
  return absl::MakeConstSpan(begin, end - begin);
}

// Returns a MutableArraySlice referring to a container or array. Does not
// participate in overload resolution unless a Span can be constructed
// from 'c' or '&c'.
template <int&... ExplicitArgumentBarrier, typename C,
          internal_span::EnableIfMutable<internal_span::ElementT<C>> = 0>
auto MakeMutableArraySlice(C&& c)
    -> decltype(absl::MakeSpan(std::forward<C>(c))) {
  return absl::MakeSpan(std::forward<C>(c));
}

template <int&... ExplicitArgumentBarrier, typename T,
          internal_span::EnableIfMutable<T> = 0>
MutableArraySlice<T> MakeMutableArraySlice(T* ptr, size_t len) {
  return absl::MakeSpan(ptr, len);
}

template <int&... ExplicitArgumentBarrier, typename T,
          internal_span::EnableIfMutable<T> = 0>
MutableArraySlice<T> MakeMutableArraySlice(T* begin, T* end) {
  return absl::MakeSpan(begin, end - begin);
}

}  // namespace gtl

// Link the pretty printer for gdb.  Only affects debug builds.
DEFINE_GDB_AUTO_SCRIPT("devtools/gdb/component/core/array_slice.py")
#endif  // S2_THIRD_PARTY_ABSL_TYPES_SPAN_H_
