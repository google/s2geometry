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

// A string_view points to part or all of a string, Cord, double-quoted string
// literal, or other string-like object.  A string_view does *not* own the
// string to which it points.  A string_view is not null-terminated.
//
// You can use string_view as a function or method parameter.  A string_view
// parameter can receive a double-quoted string literal argument, a
// "const char*" argument, a string argument, or a string_view argument with no
// data copying.  Systematic use of string_view for arguments reduces data
// copies and strlen() calls.
//
// Prefer passing string_view by value:
//   void MyFunction(string_view arg);
// If circumstances require, you may also pass by const reference:
//   void MyFunction(const string_view& arg);  // not preferred
// Both of these have the same lifetime semantics.  Passing by value generates
// slightly smaller code.
// For more discussion, see the thread go/stringpiecebyvalue on c-users.
//
// string_view is also suitable for local variables if you know that the
// lifetime of the underlying object is longer than the lifetime of your
// string_view variable.
//
// Beware of binding a string_view to a temporary:
//   string_view sv = obj.ReturnAString();  // BAD: lifetime problem
//
// This code is okay:
//   string str = obj.ReturnAString();  // str owns its contents
//   string_view sv(str);               // GOOD, because str outlives sv
//
// string_view is sometimes a poor choice for a return value and usually a poor
// choice for a data member.  If you do use a string_view this way, it is your
// responsibility to ensure that the object pointed to by the string_view
// outlives the string_view.
//
// A string_view may represent just part of a string.  For example, when
// splitting a string, std::vector<string_view> is a natural data type for the
// output.  For another example, a Cord is a non-contiguous, potentially very
// long string-like object.  The Cord class has an interface that iteratively
// provides string_view objects that point to the successive pieces of a Cord
// object.
//
// A string_view is not null-terminated.  If you write code that scans a
// string_view, you must check its length before reading any characters. Common
// idioms that work on null-terminated strings do not work on string_view
// objects.
//
// There are several ways to create a null string_view:
//   string_view()
//   string_view(nullptr, 0)
// For all of the above, sv.data() == nullptr, sv.length() == 0, and
// sv.empty() == true.  Also, if you create a string_view with a non-null
// pointer then sv.data() != nullptr.
//
// Thus, you can use string_view() to signal an out-of-band value that
// is different from other string_view values.  This is similar to the way that
// const char* p1 = nullptr; is different from const char* p2 = "";.
//
// There are many ways to create an empty string_view:
//   const char* nullcp = nullptr;
//   string_view()
//   string_view(nullcp, 0)
//   string_view("")
//   string_view("", 0)
//   string_view("abcdef", 0)
//   string_view("abcdef"+6, 0)
// For all of the above, size will be 0.
//
// Be careful not to confuse the string_view nullity with emptiness. The set of
// empty string_views properly includes the null string_view. That is, the null
// string_view is an empty string_view, but some non-null string_views are empty
// string_views too.
//
// All empty string_views whether null or not, are equal.
// string_view:
//   string_view() == string_view("", 0)
//   string_view(nullptr, 0) == string_view("abcdef"+6, 0)
//
// Look carefully at this example:
//   string_view("") == nullptr
// True or false?  TRUE, because string_view::operator== converts the right-hand
// side from nullptr to string_view(nullptr), and then compares two zero-length
// spans of characters. However, we are working to make this example produce a
// compile error.
//
// Suppose you want to write:
//   bool TestUnclear(string_view sv) { return sv == nullptr; }  // BAD
// Do not do that.  Write one of these instead:
//   bool TestNull(string_view sv) { return sv.data() == nullptr; }
//   bool TestEmpty(string_view sv) { return sv.empty(); }
// The intent of TestUnclear is unclear.  Did you mean TestNull or TestEmpty?
// Right now, TestUnclear behaves likes TestEmpty.
// We are working to make TestUnclear produce a compile error.
// TestNull is good to test for an out-of-band signal.
// TestEmpty is good to test for an empty string_view.
//
// Caveats (again):
// (1) The lifetime of the pointed-to string (or piece of a string)
//     must be longer than the lifetime of the string_view.
// (2) There may or may not be a '\0' character after the end of
//     string_view data.
// (3) A null string_view is empty.
//     An empty string_view may or may not be a null string_view.

#ifndef S2_THIRD_PARTY_ABSL_STRINGS_STRING_VIEW_H_
#define S2_THIRD_PARTY_ABSL_STRINGS_STRING_VIEW_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <iosfwd>
#include <iterator>
#include <limits>
#include <string>

#include "s2/third_party/absl/base/integral_types.h"
#include "s2/third_party/absl/base/internal/throw_delegate.h"
#include "s2/third_party/absl/base/macros.h"
#include "s2/third_party/absl/base/port.h"
#include "s2/third_party/absl/strings/internal/fastmem.h"

// string_view has *two* size types.
// string_view::size_type
//   is unsigned
//   is 32 bits in LP32, 64 bits in LP64, 64 bits in LLP64
//   no future changes intended
// stringpiece_ssize_type
//   is signed
//   is 32 bits in LP32, 64 bits in LP64, 64 bits in LLP64
//   future changes intended
//
typedef string::difference_type stringpiece_ssize_type;

namespace absl {

// Please use absl::string_view in place of StringPiece in new code.
class string_view {
 public:
  using traits_type = std::char_traits<char>;
  using value_type = char;
  using pointer = char*;
  using const_pointer = const char*;
  using reference = char&;
  using const_reference = const char&;
  using const_iterator = const char*;
  using iterator = const_iterator;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using reverse_iterator = const_reverse_iterator;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;

  static constexpr size_type npos = static_cast<size_type>(-1);

  // [string.view.cons]
  // We provide implicit conversion constructors so users can pass
  // in a "const char*" or a "string" wherever a "string_view" is
  // expected.
  // Style waiver: go/style-guide-exception-20978288
  constexpr string_view() noexcept : ptr_(nullptr), length_(0) {}

  template <class Allocator>
  string_view(  // NOLINT(runtime/explicit)
      const std::basic_string<char, std::char_traits<char>, Allocator>&
          str) noexcept
      : ptr_(str.data()), length_(str.size()) {}
#if defined(HAS_GLOBAL_STRING)
  template <class Allocator>
  string_view(  // NOLINT(runtime/explicit)
      const basic_string<char, std::char_traits<char>, Allocator>& str) noexcept
      : ptr_(str.data()), length_(str.size()) {}
#endif

  // Implicit conversion from null-terminated `str`. Argument can be
  // null-valued for now, but std::string_view will not accept null-valued
  // `str`. Please use `absl::NullSafeStringView(str)` (below) when accepting
  // a possibly-null pointer.
  constexpr string_view(const char* str)  // NOLINT(runtime/explicit)
      : ptr_(str), length_(StrLenInternal(str)) {}

  constexpr string_view(const char* data, size_type len)
      : ptr_(data), length_(CheckLengthInternal(len)) {}

  // TODO(b/36227513): harmlessly omitted to work around gdb bug.
  // constexpr string_view(const string_view&) noexcept = default;
  // string_view& operator=(const string_view&) noexcept = default;

  // [string.view.iterators]
  constexpr const_iterator begin() const noexcept { return ptr_; }
  constexpr const_iterator end() const noexcept { return ptr_ + length_; }
  constexpr const_iterator cbegin() const noexcept { return begin(); }
  constexpr const_iterator cend() const noexcept { return end(); }
  const_reverse_iterator rbegin() const noexcept {
    return const_reverse_iterator(end());
  }
  const_reverse_iterator rend() const noexcept {
    return const_reverse_iterator(begin());
  }
  const_reverse_iterator crbegin() const noexcept { return rbegin(); }
  const_reverse_iterator crend() const noexcept { return rend(); }

  // [string.view.capacity]
  constexpr size_type size() const noexcept {
    return length_;
  }
  constexpr size_type length() const noexcept { return size(); }
  constexpr size_type max_size() const noexcept { return kMaxSize; }
  constexpr bool empty() const noexcept { return length_ == 0; }

  // [string.view.access]
  constexpr const_reference operator[](size_type i) const { return ptr_[i]; }

  constexpr const_reference front() const { return ptr_[0]; }

  constexpr const_reference back() const { return ptr_[size() - 1]; }

  // data() may return a pointer to a buffer with embedded null characters, and
  // the returned buffer may or may not be null terminated.  Therefore it is
  // typically a mistake to pass data() to a routine that expects a null
  // terminated string.
  constexpr const_pointer data() const noexcept { return ptr_; }

  // [string.view.modifiers]

  void remove_prefix(size_type n) {
    assert(n <= length_);
    ptr_ += n;
    length_ -= n;
  }

  void remove_suffix(size_type n) {
    assert(n <= length_);
    length_ -= n;
  }

  void swap(string_view& s) noexcept {
    auto t = *this;
    *this = s;
    s = t;
  }

  // Checks whether string_view starts with x and if so advances the beginning
  // of it to past the match.  It's basically a shortcut for starts_with
  // followed by remove_prefix.
  ABSL_DEPRECATED("Use absl::ConsumePrefix")
  bool Consume(string_view x);
  // Like above but for the end of the string.
  ABSL_DEPRECATED("Use absl::ConsumeSuffix")
  bool ConsumeFromEnd(string_view x);

  // [string.view.ops]
  // Explicit string conversion operator. Supports conversion to both
  // std::basic_string and ::basic_string where available.
  template <class A>
  explicit operator std::basic_string<char, traits_type, A>() const {
    if (!data()) return {};
    return std::basic_string<char, traits_type, A>(data(), size());
  }
#ifdef HAS_GLOBAL_STRING
  template <class A>
  explicit operator ::basic_string<char, traits_type, A>() const {
    if (!data()) return {};
    return ::basic_string<char, traits_type, A>(data(), size());
  }
#endif  // HAS_GLOBAL_STRING

// Note that std::string_view::to_string returns std::basic_string.
// But string_view::to_string returns a ::basic_string whenever
// that template is available.
#ifdef HAS_GLOBAL_STRING
  template <class A = std::allocator<char>>
  ABSL_DEPRECATED("Use string(sp)")
  ::basic_string<char, traits_type, A> to_string(const A& a = A()) const {
    if (!data()) return ::basic_string<char, traits_type, A>(a);
    return ::basic_string<char, traits_type, A>(data(), size(), a);
  }
#else   // !HAS_GLOBAL_STRING

template <class A = std::allocator<char>>
ABSL_DEPRECATED("Use std::string(sp)")
std::basic_string<char, traits_type, A> to_string(const A& a = A()) const {
  if (!data()) return std::basic_string<char, traits_type, A>(a);
  return std::basic_string<char, traits_type, A>(data(), size(), a);
  }
#endif  // HAS_GLOBAL_STRING

  ABSL_DEPRECATED("Use string(sp)")
  string as_string() const {
    if (!data()) return {};
    return string(data(), size());
  }
  ABSL_DEPRECATED("Use string(sp)")
  string ToString() const {
    if (!data()) return {};
    return string(data(), size());
  }

  ABSL_DEPRECATED("use absl::CopyToString")
  void CopyToString(string* target) const;

  size_type copy(char* buf, size_type n, size_type pos = 0) const;

  // substr throws if `pos > size'
  // Use absl::ClippedSubstr if you need a truncating substr operation.
  string_view substr(size_type pos, size_type n = npos) const {
    if (PREDICT_FALSE(pos > length_))
      absl::ThrowStdOutOfRange("absl::string_view::substr");
    n = std::min(n, length_ - pos);
    return string_view(ptr_ + pos, n);
  }

  // returns {-1, 0, 1}
  // See http://en.cppreference.com/w/cpp/string/basic_string_view/compare
  int compare(string_view x) const noexcept {
    auto min_length = std::min(length_, x.length_);
    if (min_length > 0) {
      int r = memcmp(ptr_, x.ptr_, min_length);
      if (r < 0) return -1;
      if (r > 0) return 1;
    }
    if (length_ < x.length_) return -1;
    if (length_ > x.length_) return 1;
    return 0;
  }
  // Substring variants of `compare`.
  int compare(size_type pos1, size_type count1, string_view v) const {
    return substr(pos1, count1).compare(v);
  }
  int compare(size_type pos1, size_type count1, string_view v, size_type pos2,
              size_type count2) const {
    return substr(pos1, count1).compare(v.substr(pos2, count2));
  }
  int compare(const char* s) const { return compare(string_view(s)); }
  int compare(size_type pos1, size_type count1, const char* s) const {
    return substr(pos1, count1).compare(string_view(s));
  }
  int compare(size_type pos1, size_type count1, const char* s,
              size_type count2) const {
    return substr(pos1, count1).compare(string_view(s, count2));
  }

  // Returns whether this begins with x.
  ABSL_DEPRECATED("Use absl::StartsWith")
  bool starts_with(string_view x) const {
    return length_ >= x.length_ &&
           absl::strings_internal::memeq(ptr_, x.ptr_, x.length_);
  }

  // Returns whether this ends with x.
  ABSL_DEPRECATED("Use absl::EndsWith")
  bool ends_with(string_view x) const {
    return length_ >= x.length_ &&
           absl::strings_internal::memeq(ptr_ + (length_ - x.length_), x.ptr_,
                                         x.length_);
  }

  ABSL_DEPRECATED("Use absl::StrContains")
  bool contains(string_view s) const;

  // [string.view.find]

  size_type find(string_view s, size_type pos = 0) const noexcept;
  size_type find(char c, size_type pos = 0) const noexcept;

  size_type rfind(string_view s, size_type pos = npos) const
      noexcept;
  size_type rfind(char c, size_type pos = npos) const noexcept;

  size_type find_first_of(string_view s, size_type pos = 0) const
      noexcept;
  size_type find_first_of(char c, size_type pos = 0) const
      noexcept {
    return find(c, pos);
  }

  size_type find_last_of(string_view s, size_type pos = npos) const
      noexcept;
  size_type find_last_of(char c, size_type pos = npos) const
      noexcept {
    return rfind(c, pos);
  }

  size_type find_first_not_of(string_view s,
                                           size_type pos = 0) const noexcept;
  size_type find_first_not_of(char c, size_type pos = 0) const
      noexcept;

  size_type find_last_not_of(string_view s,
                                          size_type pos = npos) const noexcept;
  size_type find_last_not_of(char c, size_type pos = npos) const
      noexcept;

  // Legacy variants of the size and find family returning their old signed
  // type.
  // Will eventually be removed by go/lsc-stringpiece-size as call sites are
  // migrated to safely accept the size_type return type.
  constexpr stringpiece_ssize_type ssize() const noexcept { return size(); }
  constexpr stringpiece_ssize_type slength() const noexcept { return length(); }
  stringpiece_ssize_type sfind(string_view s, size_type pos = 0) const
      noexcept {
    return find(s, pos);
  }
  stringpiece_ssize_type sfind(char c, size_type pos = 0) const noexcept {
    return find(c, pos);
  }
  stringpiece_ssize_type srfind(string_view s, size_type pos = npos) const
      noexcept {
    return rfind(s, pos);
  }
  stringpiece_ssize_type srfind(char c, size_type pos = npos) const noexcept {
    return rfind(c, pos);
  }
  stringpiece_ssize_type sfind_first_of(string_view s, size_type pos = 0) const
      noexcept {
    return find_first_of(s, pos);
  }
  stringpiece_ssize_type sfind_first_of(char c, size_type pos = 0) const
      noexcept {
    return find_first_of(c, pos);
  }
  stringpiece_ssize_type sfind_last_of(string_view s,
                                       size_type pos = npos) const noexcept {
    return find_last_of(s, pos);
  }
  stringpiece_ssize_type sfind_last_of(char c, size_type pos = npos) const
      noexcept {
    return find_last_of(c, pos);
  }
  stringpiece_ssize_type sfind_first_not_of(string_view s,
                                            size_type pos = 0) const noexcept {
    return find_first_not_of(s, pos);
  }
  stringpiece_ssize_type sfind_first_not_of(char c, size_type pos = 0) const
      noexcept {
    return find_first_not_of(c, pos);
  }
  stringpiece_ssize_type sfind_last_not_of(string_view s,
                                           size_type pos = npos) const
      noexcept {
    return find_last_not_of(s, pos);
  }
  stringpiece_ssize_type sfind_last_not_of(char c, size_type pos = npos) const
      noexcept {
    return find_last_not_of(c, pos);
  }

 private:
  static constexpr size_type kMaxSize =
      std::numeric_limits<size_type>::max() / 2 + 1;

  static constexpr size_type StrLenInternal(const char* str) {
    return str ?
// check whether __builtin_strlen is provided by the compiler.
// GCC doesn't have __has_builtin()
// (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=66970),
// but has __builtin_strlen according to
// https://gcc.gnu.org/onlinedocs/gcc-4.7.0/gcc/Other-Builtins.html.
#if ABSL_HAVE_BUILTIN(__builtin_strlen) || \
    (defined(__GNUC__) && !defined(__clang__))
               __builtin_strlen(str)
#else
               strlen(str)
#endif
               : 0;
  }

  static constexpr size_type CheckLengthInternal(size_type len) {
    return ABSL_ASSERT(len <= kMaxSize), len;
  }

  const char* ptr_;
  size_type length_;
};

// This large function is defined inline so that in a fairly common case where
// one of the arguments is a literal, the compiler can elide a lot of the
// following comparisons.
inline bool operator==(string_view x, string_view y) noexcept {
  auto len = x.size();
  if (len != y.size()) {
    return false;
  }

  return x.data() == y.data() || len <= 0 ||
         absl::strings_internal::memeq(x.data(), y.data(), len);
  /* absl:oss-replace-with
  return x.data() == y.data() || len <= 0 ||
         memcmp(x.data(), y.data(), len) == 0;
  absl:oss-replace-end */
}

inline bool operator!=(string_view x, string_view y) noexcept {
  return !(x == y);
}

inline bool operator<(string_view x, string_view y) noexcept {
  auto min_size = std::min(x.size(), y.size());
  const int r = min_size == 0 ? 0 : memcmp(x.data(), y.data(), min_size);
  return (r < 0) || (r == 0 && x.size() < y.size());
}

inline bool operator>(string_view x, string_view y) noexcept { return y < x; }

inline bool operator<=(string_view x, string_view y) noexcept {
  return !(y < x);
}

inline bool operator>=(string_view x, string_view y) noexcept {
  return !(x < y);
}

// [string.view.io]
std::ostream& operator<<(std::ostream& o, string_view piece);

// Like `s.substr(pos, n)`, but clips `pos` to an upper bound of `s.size()`.
// Provided because std::string_view::substr throws if `pos > size()`,
// to support b/37991613.
inline string_view ClippedSubstr(string_view s, size_t pos,
                                 size_t n = string_view::npos) {
  pos = std::min(pos, static_cast<size_t>(s.size()));
  return s.substr(pos, n);
}

// Create an absl::string_view from a `p` even if it's null-valued. Should be
// used where a absl::string_view must be created from a possibly-null pointer.
// Our absl::string_view has historically been constructible from a null-valued
// pointer, but the same null value isn't valid for std::string_view.
inline string_view NullSafeStringView(const char* p) {
  return p ? string_view(p) : string_view();
}

}  // namespace absl


#endif  // S2_THIRD_PARTY_ABSL_STRINGS_STRING_VIEW_H_
