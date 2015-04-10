// Copyright 2012 Google Inc. All Rights Reserved.
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
// This file declares INTERNAL parts of the Split API that are inline/templated
// or otherwise need to be available at compile time. The main two abstractions
// defined in here are
//
//   - SplitIterator<>
//   - Splitter<>
//
// Everything else is plumbing for those two.
//
// DO NOT INCLUDE THIS FILE DIRECTLY. Use this file by including
// strings/split.h.
//
// IWYU pragma: private, include "strings/split.h"

#ifndef STRINGS_SPLIT_INTERNAL_H_
#define STRINGS_SPLIT_INTERNAL_H_

#ifdef _GLIBCXX_DEBUG
#include <glibcxx_debug_traits.h>
#endif  // _GLIBCXX_DEBUG

#include <algorithm>
#include <iterator>
#include <map>
#include <vector>

#include "base/port.h"  // for LANG_CXX11
#include "base/template_util.h"
#include "base/type_traits.h"  // for base::enable_if
#include "strings/stringpiece.h"

#ifdef LANG_CXX11
// This must be included after "base/port.h", which defines LANG_CXX11.
#include <initializer_list>  // NOLINT(build/include_order)
#endif  // LANG_CXX11

namespace strings {

namespace internal {

#ifdef _GLIBCXX_DEBUG
using ::glibcxx_debug_traits::IsStrictlyDebugWrapperBase;
#else  // _GLIBCXX_DEBUG
template <typename T> struct IsStrictlyDebugWrapperBase : base::false_ {};
#endif  // _GLIBCXX_DEBUG

// The default Predicate object, which doesn't filter out anything.
struct NoFilter {
  bool operator()(StringPiece /* ignored */) {
    return true;
  }
};

// This class is implicitly constructible from const char*,
// StringPiece, and string, but not from other types that are
// themselves convertible to StringPiece (because C++ won't implicitly
// use a sequence of 2 user-defined conversion).  If it's constructed
// from a temporary string, it moves that temporary into local storage
// so that the temporary isn't destroyed before Split()'s result.
class ConvertibleToStringPiece {
 private:
  // Used for temporary string arguments.  Must be declared before
  // 'value' in order to construct 'value' from it.
  string copy_;

 public:
  StringPiece value;
  ConvertibleToStringPiece(const char* s) : value(s) {}    // NOLINT
  ConvertibleToStringPiece(char* s) : value(s) {}          // NOLINT
  ConvertibleToStringPiece(StringPiece s) : value(s) {}    // NOLINT
  ConvertibleToStringPiece(const string& s) : value(s) {}  // NOLINT

#if LANG_CXX11
  // This overload captures temporary arguments.
  // gave approval for this use in cl/48636778.
  ConvertibleToStringPiece(string&& s)  // NOLINT
      : copy_(std::move(s)),            // NOLINT
        value(copy_) {}
#endif

  ConvertibleToStringPiece(const ConvertibleToStringPiece& other)
      : copy_(other.copy_) {
    if (other.copy_.empty()) {
      value = other.value;
    } else {
      value = copy_;
    }
  }

  // Effectively a move constructor, using C++98 features.
  explicit ConvertibleToStringPiece(ConvertibleToStringPiece* other) {
    if (other->copy_.empty()) {
      value = other->value;
    } else {
      using std::swap;
      swap(copy_, other->copy_);
      value = copy_;
      other->value = StringPiece();
    }
  }

  void swap(ConvertibleToStringPiece& other) {
    using std::swap;
    swap(copy_, other.copy_);
    swap(value, other.value);

    // Fix up StringPieces when the real data is stored in the string.
    if (!copy_.empty()) value = copy_;
    if (!other.copy_.empty()) other.value = other.copy_;
  }

  ConvertibleToStringPiece& operator=(ConvertibleToStringPiece other) {
    swap(other);
    return *this;
  }
};

// This class splits a string using the given delimiter, returning the split
// substrings via an iterator interface. An optional Predicate functor may be
// supplied, which will be used to filter the split strings: strings for which
// the predicate returns false will be skipped. A Predicate object is any
// functor that takes a StringPiece and returns bool. By default, the NoFilter
// Predicate is used, which does not filter out anything.
//
// This class is NOT part of the public splitting API.
//
// Usage:
//
//   using strings::delimiter::Literal;
//   Literal d(",");
//   for (SplitIterator<Literal> it("a,b,c", d), end(d); it != end; ++it) {
//     StringPiece substring = *it;
//     DoWork(substring);
//   }
//
// The explicit single-argument constructor is used to create an "end" iterator.
// The two-argument constructor is used to split the given text using the given
// delimiter.
template <typename Delimiter, typename Predicate = NoFilter>
class SplitIterator
    : public std::iterator<std::forward_iterator_tag, StringPiece> {
 public:
  // Two constructors for "end" iterators.
  explicit SplitIterator(Delimiter d)
      : delimiter_(d), predicate_(), is_end_(true) {}
  SplitIterator(Delimiter d, Predicate p)
      : delimiter_(d), predicate_(p), is_end_(true) {}

  // Two constructors taking the text to iterate.
  SplitIterator(StringPiece text, Delimiter d)
      : text_(text), position_(0), delimiter_(d), predicate_(), is_end_(false) {
    ++(*this);
  }
  SplitIterator(StringPiece text, Delimiter d, Predicate p)
      : text_(text), position_(0), delimiter_(d), predicate_(p),
        is_end_(false) {
    ++(*this);
  }

  // Constructor copying the delimiter and predicate out of an
  // existing SplitIterator, but applying them to new text.
  SplitIterator(StringPiece text, const SplitIterator& other)
      : text_(text),
        position_(0),
        delimiter_(other.delimiter_),
        predicate_(other.predicate_),
        is_end_(false) {
    ++(*this);
  }

  StringPiece operator*() const { return curr_piece_; }
  const StringPiece* operator->() const { return &curr_piece_; }

  SplitIterator& operator++() {
    do {
      if (text_.end() == curr_piece_.end()) {
        // Already consumed all of text_, so we're done.
        is_end_ = true;
        return *this;
      }
      StringPiece found_delimiter = delimiter_.Find(text_, position_);
      assert(found_delimiter.data() != NULL);
      assert(text_.begin() + position_ <= found_delimiter.begin());
      assert(found_delimiter.end() <= text_.end());
      // found_delimiter is allowed to be empty.
      // Sets curr_piece_ to all text up to but excluding the delimiter itself.
      // Increments position_ by the length of curr_piece_ and found_delimiter.
      size_t curr_size = found_delimiter.begin() - (text_.begin() + position_);
      curr_piece_.set(text_.begin() + position_, curr_size);
      position_ += curr_piece_.size() + found_delimiter.size();
    } while (!predicate_(curr_piece_));
    return *this;
  }

  SplitIterator operator++(int /* postincrement */) {
    SplitIterator old(*this);
    ++(*this);
    return old;
  }

  bool operator==(const SplitIterator& other) const {
    // Two "end" iterators are always equal. If the two iterators being compared
    // aren't both end iterators, then we fallback to comparing their fields.
    // Importantly, the text being split must be equal and the current piece
    // within the text being split must also be equal. The delimiter_ and
    // predicate_ fields need not be checked here because they're template
    // parameters that are already part of the SplitIterator's type.
    return (is_end_ && other.is_end_) ||
           (is_end_ == other.is_end_ &&
            text_ == other.text_ &&
            text_.data() == other.text_.data() &&
            position_ == other.position_ &&
            curr_piece_ == other.curr_piece_ &&
            curr_piece_.data() == other.curr_piece_.data());
  }

  bool operator!=(const SplitIterator& other) const {
    return !(*this == other);
  }

 private:
  // The text being split.
  StringPiece text_;
  size_t position_;
  Delimiter delimiter_;
  Predicate predicate_;
  bool is_end_;
  // Holds the currently split piece of text. Will always refer to string data
  // within text_. This value is returned when the iterator is dereferenced.
  StringPiece curr_piece_;
};

// Declares a functor that can convert a StringPiece to another type. This works
// for any type that has a constructor (explicit or not) taking a single
// StringPiece argument. A specialization exists for converting to string
// because the underlying data needs to be copied. In theory, these
// specializations could be extended to work with other types (e.g., int32), but
// then a solution for error reporting would need to be devised.
template <typename To>
struct StringPieceTo {
  To operator()(StringPiece from) const {
    return To(from);
  }
};

// Specialization for converting to string.
template <>
struct StringPieceTo<string> {
  string operator()(StringPiece from) const {
    return from.ToString();
  }
};

// Specialization for converting to *const* string.
template <>
struct StringPieceTo<const string> {
  string operator()(StringPiece from) const {
    return from.ToString();
  }
};

// HasMappedType<T>::value is true iff there exists a type T::mapped_type.
template <typename T>
struct HasMappedType {
  template <typename U> static base::small_ check(...);  // default: No
  template <typename U> static base::big_ check(typename U::mapped_type*);
  enum { value = sizeof(base::big_) == sizeof(check<T>(0)) };
};

// HasValueType<T>::value is true iff there exists a type T::value_type.
template <typename T>
struct HasValueType {
  template <typename U> static base::small_ check(...);  // default: No
  template <typename U> static base::big_ check(typename U::value_type*);
  enum { value = sizeof(base::big_) == sizeof(check<T>(0)) };
};

// HasConstIterator<T>::value is true iff there exists a type
// T::const_iterator.
template <typename T>
struct HasConstIterator {
  template <typename U> static base::small_ check(...);  // default: No
  template <typename U> static base::big_ check(typename U::const_iterator*);
  enum { value = sizeof(base::big_) == sizeof(check<T>(0)) };
};

// IsInitializerList<T>::value is true iff T is an initializer_list.
// More details below in Splitter<> where this is used.
template <typename T>
struct IsInitializerList {
  static base::small_ check(...);  // default: No
#ifdef LANG_CXX11
  template <typename U> static base::big_ check(std::initializer_list<U>*);
#endif  // LANG_CXX11
  enum { value = sizeof(base::big_) == sizeof(check(static_cast<T*>(0))) };
};

// A SplitterIsConvertibleTo<C>::type typedef exists iff the specified
// condition is true for type 'C'.
//
// Restricts conversion to container-like types (by testing for the presence of
// a const_iterator member type) and also to disable conversion to an
// initializer_list (which also has a const_iterator). Otherwise, code compiled
// in C++11 will get an error due to ambiguous conversion paths (in C++11
// vector<T>::operator= is overloaded to take either a vector<T> or an
// initializer_list<T>).
//
// This trick was taken from util/gtl/container_literal.h
template <typename C>
struct SplitterIsConvertibleTo
    : base::enable_if<
          !IsStrictlyDebugWrapperBase<C>::value &&
          !IsInitializerList<C>::value &&
          HasValueType<C>::value &&
          HasConstIterator<C>::value> {};  // NOLINT

// This class implements the behavior of the split API by giving callers access
// to the underlying split substrings in various convenient ways, such as
// through iterators or implicit conversion functions. Do not construct this
// class directly, rather use the Split() function instead.
//
// Output containers can be collections of either StringPiece or string objects.
// StringPiece is more efficient because the underlying data will not need to be
// copied; the returned StringPieces will all refer to the data within the
// original input string. If a collection of string objects is used, then each
// substring will be copied.
//
// An optional Predicate functor may be supplied. This predicate will be used to
// filter the split strings: only strings for which the predicate returns true
// will be kept. A Predicate object is any unary functor that takes a
// StringPiece and returns bool. By default, the NoFilter predicate is used,
// which does not filter out anything.
template <typename Delimiter, typename Predicate = NoFilter>
class Splitter {
 public:
  typedef std::forward_iterator_tag iterator_category;
  typedef internal::SplitIterator<Delimiter, Predicate> iterator;
  typedef internal::SplitIterator<Delimiter, Predicate> const_iterator;

  Splitter(ConvertibleToStringPiece* input_text, Delimiter d)
      : text_(input_text), begin_(text_.value, d), end_(d) {}

  Splitter(ConvertibleToStringPiece* input_text, Delimiter d, Predicate p)
      : text_(input_text), begin_(text_.value, d, p), end_(d, p) {}

  Splitter(const Splitter& other)
      : text_(other.text_),
        begin_(text_.value, other.begin_),
        end_(other.end_) {}

  // Range functions that iterate the split substrings as StringPiece objects.
  // These methods enable a Splitter to be used in a range-based for loop in
  // C++11, for example:
  //
  //   for (StringPiece sp : my_splitter) {
  //     DoWork(sp);
  //   }
  const_iterator begin() const { return begin_; }
  const_iterator end() const { return end_; }

#ifdef LANG_CXX11
  // An implicit conversion operator that uses a default template argument to
  // restrict its use to only those containers that the splitter is convertible
  // to. Default template arguments on function templates are a C++11 feature,
  // which is why this code is under an #ifdef LANG_CXX11.
  template <typename Container,
            typename OnlyIf =
                typename SplitterIsConvertibleTo<Container>::type>
  operator Container() const {
    return ContainerConverter<
        Container,
        typename Container::value_type,
        HasMappedType<Container>::value>()(*this);
  }
#else
  // Not under LANG_CXX11.
  // Same conversion operator as in the C++11 block above, except this one
  // doesn't use SplitterIsConvertibleTo<>.
  template <typename Container>
  operator Container() const {
    return ContainerConverter<
        Container,
        typename Container::value_type,
        HasMappedType<Container>::value>()(*this);
  }
#endif  // LANG_CXX11

  // Returns a pair with its .first and .second members set to the first two
  // strings returned by the begin() iterator. Either/both of .first and .second
  // will be constructed with empty strings if the iterator doesn't have a
  // corresponding value.
  template <typename First, typename Second>
  operator std::pair<First, Second>() const {
    StringPieceTo<First> first_converter;
    StringPieceTo<Second> second_converter;
    StringPiece first, second;
    const_iterator it = begin();
    if (it != end()) {
      first = *it;
      if (++it != end()) {
        second = *it;
      }
    }
    return std::make_pair(first_converter(first), second_converter(second));
  }

 private:
  // ContainerConverter is a functor converting a Splitter to the requested
  // Container and ValueType. It can be specialized to optimize splitting to
  // certain combinations of Container and ValueType.
  //
  // This base template handles the generic case of storing the split results in
  // the requested non-map-like container and converting the split substrings to
  // the requested type.
  template <typename Container, typename ValueType, bool is_map = false>
  struct ContainerConverter {
    Container operator()(const Splitter& splitter) const {
      Container c;
      std::transform(splitter.begin(),
                     splitter.end(),
                     std::inserter(c, c.end()),
                     StringPieceTo<typename Container::value_type>());
      return c;
    }
  };

  // Partial specialization for a vector<string>.
  //
  // Optimized for the common case of splitting to a vector<string>. In this
  // case we first split the results to a vector<StringPiece> so the returned
  // vector<string> can have space reserved to avoid string copies.
  template <typename A>
  struct ContainerConverter<std::vector<string, A>, string, false> {
    std::vector<string, A> operator()(const Splitter& splitter) const {
      std::vector<StringPiece> vsp(splitter.begin(), splitter.end());
      std::vector<string, A> v(vsp.size());
      for (size_t i = 0; i < vsp.size(); ++i) {
        vsp[i].CopyToString(&v[i]);
      }
      return v;
    }
  };

  // Partial specialization for containers of pairs (e.g., maps).
  //
  // The algorithm is to insert a new pair into the map for each even-numbered
  // item, with the even-numbered item as the key with a default-constructed
  // value. Each odd-numbered item will then be assigned to the last pair's
  // value.
  template <typename Container, typename First, typename Second>
  struct ContainerConverter<Container, std::pair<First, Second>, true> {
    Container operator()(const Splitter& splitter) const {
      Container m;
      StringPieceTo<First> key_converter;
      StringPieceTo<Second> val_converter;
      typename Container::iterator curr_pair;
      bool is_even = true;
      for (const_iterator it = splitter.begin(); it != splitter.end(); ++it) {
        if (is_even) {
          curr_pair = InsertInMap(std::make_pair(key_converter(*it), Second()),
                                  &m);
        } else {
          curr_pair->second = val_converter(*it);
        }
        is_even = !is_even;
      }
      return m;
    }

    // Overloaded InsertInMap() function. The first overload is the commonly
    // used one for most map-like objects. The second overload is a special case
    // for multimap, because multimap's insert() member function directly
    // returns an iterator, rather than a pair<iterator, bool> like map's.
    template <typename Map>
    typename Map::iterator InsertInMap(
        const typename Map::value_type& value, Map* map) const {
      return map->insert(value).first;
    }

    // InsertInMap overload for multimap.
    template <typename K, typename T, typename C, typename A>
    typename std::multimap<K, T, C, A>::iterator InsertInMap(
        const typename std::multimap<K, T, C, A>::value_type& value,
        typename std::multimap<K, T, C, A>* map) const {
      return map->insert(value);
    }
  };

  const ConvertibleToStringPiece text_;
  const const_iterator begin_;
  const const_iterator end_;
};

}  // namespace internal

}  // namespace strings

#endif  // STRINGS_SPLIT_INTERNAL_H_