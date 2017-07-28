// Copyright 2009 Google Inc. All Rights Reserved.
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
// Based on ideas suggested by Sanjay Ghemawat (sanjay@google.com)
//
// A Span<T> represents an array of elements of type T.  It is conceptually a
// pointer (ptr) and a length (size), and the array it represents contains the
// elements "ptr[0] .. ptr[size-1]". The backing store for the array is *not*
// owned by the Span object, and clients must arrange for the backing store to
// remain live while the Span object is in use.
//
// A Span<T> is somewhat analogous to an absl::string_view, but for array
// elements of type T, and the elements can be mutated if T is mutable.
//
// Implicit conversion operations are provided from any type supplying data()
// and size() const methods (e.g. std::vector<T> and absl::InlinedVector<T, N>).
// Note that Span objects constructed from types in this way may be
// invalidated by any operations that mutate the underlying vector.
//
// One common use for Span is when passing arguments to a routine where you want
// to be able to accept a variety of array types (e.g. a vector, an
// absl::InlinedVector, a C-style array, etc.).  The usual approach here is to
// have the client explicitly pass in a pointer and a length, as in:
//
//   void MyRoutine(const int* elems, int N) {
//     for (int i = 0; i < N; i++) { .. do something with elems[i] .. }
//   }
//
// Unfortunately, this leads to ugly and error-prone code at the call site:
//
//   std::vector<int> my_vector;
//   MyRoutine(my_vector.data(), my_vector.size());
//
//   absl::InlinedVector<int, 4> my_inline_vector;
//   MyRoutine(my_inline_vector.array(), my_inline_vector.size());
//
//   absl::FixedArray<int> my_fixed_array(10);
//   MyRoutine(my_fixed_array.get(), my_fixed_array.size());
//
//   int my_array[10];
//   MyRoutine(my_array, 10);
//
// Instead, you can use a Span as the argument to the routine. If you only need
// to view the elements, take a Span<const T>:
//
//   void MyRoutine(absl::Span<const int> a) {
//     for (int i = 0; i < a.size(); i++) { .. do something with a[i] .. }
//   }
//
// This makes the call sites cleaner, for the most part:
//
//   std::vector<int> my_vector;
//   MyRoutine(my_vector);
//
//   absl::InlinedVector<int, 4> my_inline_vector;
//   MyRoutine(my_inline_vector);
//
//   absl::FixedArray<int> my_fixed_array(10);
//   MyRoutine(my_fixed_array);
//
//   int my_array[10];
//   MyRoutine(my_array);
//
//   int* my_array = new int[10];
//   MyRoutine(absl::Span<int>(my_array, 10));
//
//  If, instead, you want to mutate the values through the Span, callers must
//  use the by-reference constructor:
//
//    void MyRoutine(absl::Span<int> a) { a[0] = 1; }
//
// Note that this constructor is explicit.
//
// For convenience, you may also use any of the `MakeSpan` or `MakeConstSpan`
// factory functions. This is recommended to avoid issues of implicit conversion
// when using the explicit constructor above is unwieldy.
//
//    std::vector<int> my_vector;
//    MyRoutine(Span<int>(my_vector));
//    MyRoutine(absl::MakeSpan(my_vector));
//
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

#ifndef S2_THIRD_PARTY_ABSL_TYPES_SPAN_H_
#define S2_THIRD_PARTY_ABSL_TYPES_SPAN_H_

#include <cstddef>
#include <algorithm>
#include <cassert>
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
  return s.empty() ? nullptr : &*s.begin();
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

// We want to enable conversion from vector<T*> to Span<const T*> but
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

// This struct allows us to write Span's relational operators in such a way that
// we can compare any type convertible to a Span and a Span, regardless of left
// or right-hand positioning or the const-ness of the underlying elements.
template <typename T>
struct SpanRelationals {
  static_assert(
      std::is_const<T>::value,
      "SpanRelationals must be used with a const parameter type to work.");

  // These relational operators have the same semantics as the
  // vector<T> relational operators: they do deep (elementwise)
  // comparisons.  Spans are equal iff their size is the same and all their
  // elements are equal.  Spans are ordered by lexicographical comparison.
  friend bool operator==(Span<T> lhs, Span<T> rhs) noexcept {
    return absl::equal(lhs.begin(), lhs.end(), rhs.begin(), rhs.end());
  }

  friend bool operator!=(Span<T> lhs, Span<T> rhs) noexcept {
    return !(lhs == rhs);
  }

  friend bool operator<(Span<T> lhs, Span<T> rhs) noexcept {
    return std::lexicographical_compare(lhs.begin(), lhs.end(), rhs.begin(),
                                        rhs.end());
  }

  // We can't use std::lexicographical_compare with std::greater, because then
  // an empty range would be both less than and greater than a non-empty range.
  friend bool operator>(Span<T> lhs, Span<T> rhs) noexcept { return rhs < lhs; }

  friend bool operator<=(Span<T> lhs, Span<T> rhs) noexcept {
    return !(rhs < lhs);
  }

  friend bool operator>=(Span<T> lhs, Span<T> rhs) noexcept {
    return !(lhs < rhs);
  }
};
}  // namespace internal_span

template <typename T>
class Span : private internal_span::SpanRelationals<const T> {
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
  using value_type = absl::remove_const_t<T>;
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
  constexpr explicit Span(std::nullptr_t) noexcept : Span() {}
  constexpr Span(pointer array, size_type length) noexcept
      : ptr_(array), len_(length) {}

  // Implicit conversion constructors
  template <size_t N>
  constexpr Span(T (&a)[N]) noexcept  // NOLINT(runtime/explicit)
      : Span(a, N) {}

  // The constructor for any class supplying 'T* data()' or 'T* mutable_data()'
  // (the former is called if both exist), and 'some_integral_type size()
  // const'. proto2::RepeatedField is an example of this. Also supports string
  // arguments, when T==char. The appropriate ctor is selected using SFINAE.
  template <typename V, typename = EnableIfConvertibleFrom<V>,
            typename = EnableIfMutableView<V>>
  ABSL_DEPRECATED("Use the explicit mutable ref conversion ctor")
  Span(V* v)  // NOLINT(runtime/explicit)
      : Span(internal_span::GetData(*v), v->size()) {}

  template <typename V, typename = EnableIfConvertibleFrom<V>,
            typename = EnableIfMutableView<V>>
  explicit Span(V& v) noexcept  // NOLINT(runtime/references)
      : Span(internal_span::GetData(v), v.size()) {}

  template <typename V, typename = EnableIfConvertibleFrom<V>,
            typename = EnableIfConstView<V>>
  constexpr Span(const V& v) noexcept  // NOLINT(runtime/explicit)
      : Span(internal_span::GetData(v), v.size()) {}

  // Implicitly constructs a Span from an initializer list, making it
  // possible to pass a brace-enclosed initializer list to a function expecting
  // a Span:
  //
  //   void Process(Span<const int> x);
  //   Process({1, 2, 3});
  //
  // The data referenced by the initializer list must outlive this
  // Span. For example, "Span<int> s={1,2};" and "return
  // Span<int>({3,4});" are errors, as the resulting Span may
  // reference data that is no longer valid.
  //
  // Only read-only Spans may be constructed from initializer_lists.
  template <typename LazyT = T,
            typename = EnableIfConstView<LazyT>>
  constexpr Span(
      std::initializer_list<value_type> v) noexcept  // NOLINT(runtime/explicit)
      : Span(v.begin(), v.size()) {}

  // Substring of another Span.
  // pos must be non-negative and <= x.length().
  // len must be non-negative and will be pinned to at most x.length() - pos.
  // If len==npos, the substring continues till the end of x.
  constexpr Span(Span x, size_type pos, size_type len) noexcept
      : ptr_(x.ptr_ + pos), len_(internal_span::Min(x.len_ - pos, len)) {}

  // Accessors.
  constexpr pointer data() const noexcept { return ptr_; }
  constexpr size_type size() const noexcept { return len_; }
  constexpr size_type length() const noexcept { return size(); }
  constexpr bool empty() const noexcept { return size() == 0; }

  constexpr reference operator[](size_type i) const noexcept { return ptr_[i]; }
  constexpr reference at(size_type i) const {
    return ABSL_ASSERT(i < size()), ptr_[i];
  }
  constexpr reference front() const noexcept {
    return ABSL_ASSERT(size() > 0), ptr_[0];
  }
  constexpr reference back() const noexcept {
    return ABSL_ASSERT(size() > 0), ptr_[size() - 1];
  }

  constexpr iterator begin() const noexcept { return ptr_; }
  constexpr const_iterator cbegin() const noexcept { return ptr_; }

  constexpr iterator end() const noexcept { return ptr_ + len_; }
  constexpr const_iterator cend() const noexcept { return end(); }

  reverse_iterator rbegin() const noexcept { return reverse_iterator(end()); }
  const_reverse_iterator crbegin() const noexcept { return rbegin(); }

  reverse_iterator rend() const noexcept { return reverse_iterator(begin()); }
  const_reverse_iterator crend() const noexcept { return rend(); }

  void remove_prefix(size_type n) noexcept {
    assert(len_ >= n);
    ptr_ += n;
    len_ -= n;
  }
  void remove_suffix(size_type n) noexcept {
    assert(len_ >= n);
    len_ -= n;
  }

  ABSL_DEPRECATED("Use remove_suffix(1) instead.")
  void pop_back() { remove_suffix(1); }
  ABSL_DEPRECATED("Use remove_prefix(1) instead.")
  void pop_front() { remove_prefix(1); }

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

// MakeSpan constructs a Span<T>, deducing T automatically.  The
// pointer-accepting versions return a Span<const T> if T is const, and a
// Span<T> otherwise.  Similarly, the container-accepting version of MakeSpan
// deduces the type of T by the constness of the pointer received from .data().
// MakeConstSpan always returns a Span<const T>.
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
