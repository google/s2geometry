// Copyright 2015 Google Inc. All Rights Reserved.
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
// This file declares INTERNAL parts of the Join API that are inlined/templated
// or otherwise need to be available at compile time. The main abstractions
// defined in this file are:
//
//   - A handful of default Formatters
//   - JoinAlgorithm() overloads
//   - JoinRange() overloads
//   - JoinTuple()
//
// DO NOT INCLUDE THIS FILE DIRECTLY. Use this file by including
// strings/join.h.
//
// IWYU pragma: private, include "strings/join.h"

#ifndef STRINGS_JOIN_INTERNAL_H_
#define STRINGS_JOIN_INTERNAL_H_

#include <iterator>
#include <string>

#include "strings/strcat.h"

#ifdef LANG_CXX11
#include <memory>
#include <utility>
#endif

namespace strings {
namespace internal {

//
// Formatter objects
//
// The following are implementation classes for standard Formatter objects. The
// factory functions that users will call to create and use these formatters are
// defined and documented in strings/join.h.
//

// The default formatter. Converts alpha-numeric types to strings.
struct AlphaNumFormatterImpl {
  template <typename T>
  void operator()(string* out, const T& t) const {
    AlphaNum(t).Piece().AppendToString(out);
  }

  void operator()(string* out, const AlphaNum& t) const {
    t.Piece().AppendToString(out);
  }
};


// A type that's used to overload the JoinAlgorithm() function (defined below)
// for ranges that do not require additional formatting (e.g., a range of
// strings).

struct NoFormatter : public AlphaNumFormatterImpl {};

// Formats a std::pair<>. The 'first' member is formatted using f1_ and the
// 'second' member is formatted using f2_. sep_ is the separator.
template <typename F1, typename F2>
class PairFormatterImpl {
 public:
  PairFormatterImpl(F1 f1, StringPiece sep, F2 f2)
      : f1_(f1), sep_(sep.ToString()), f2_(f2) {}

  template <typename T>
  void operator()(string* out, const T& p) const {
    f1_(out, p.first);
    out->append(sep_);
    f2_(out, p.second);
  }

 private:
  F1 f1_;
  string sep_;
  F2 f2_;
};

// Wraps another formatter and dereferences the argument to operator() then
// passes the dereferenced argument to the wrapped formatter. This can be
// useful, for example, to join a vector<int*>.
template <typename Formatter>
class DereferenceFormatterImpl {
 public:
  DereferenceFormatterImpl() : f_() {}
  explicit DereferenceFormatterImpl(Formatter f) : f_(f) {}

  template <typename T>
  void operator()(string* out, const T& t) const {
    f_(out, *t);
  }

 private:
  Formatter f_;
};

// DefaultFormatter<T> is a traits class that selects a default Formatter to use
// for the given type T. The ::Type member names the Formatter to use. This is
// used by the strings::Join() functions that do NOT take a Formatter argument,
// in which case a default Formatter must be chosen.
//
// AlphaNumFormatterImpl is the default in the base template, followed by
// specializations for other types.
template <typename ValueType>
struct DefaultFormatter {
  typedef AlphaNumFormatterImpl Type;
};
template <>
struct DefaultFormatter<const char*> {
  typedef AlphaNumFormatterImpl Type;
};
template <>
struct DefaultFormatter<char*> {
  typedef AlphaNumFormatterImpl Type;
};
template <>
struct DefaultFormatter<string> {
  typedef NoFormatter Type;
};
template <>
struct DefaultFormatter<StringPiece> {
  typedef NoFormatter Type;
};
template <typename ValueType>
struct DefaultFormatter<ValueType*> {
  typedef DereferenceFormatterImpl<typename DefaultFormatter<ValueType>::Type>
      Type;
};

#ifdef LANG_CXX11
template <typename ValueType>
struct DefaultFormatter<std::unique_ptr<ValueType>>
    : public DefaultFormatter<ValueType*> {};
#endif

//
// JoinAlgorithm() functions
//

// The main joining algorithm. This simply joins the elements in the given
// iterator range, each separated by the given separator, into an output string,
// and formats each element using the provided Formatter object.
template <typename Iterator, typename Formatter>
string JoinAlgorithm(Iterator start, Iterator end, StringPiece s, Formatter f) {
  string result;
  StringPiece sep("");
  for (Iterator it = start; it != end; ++it) {
    result.append(sep.data(), sep.size());
    f(&result, *it);
    sep = s;
  }
  return result;
}

// No-op placeholder for input iterators which can not be iterated over.
template <typename Iterator>
size_t GetResultSize(Iterator, Iterator, size_t, std::input_iterator_tag) {
  return 0;
}

// Calculates space to reserve, if the iterator supports multiple passes.
template <typename Iterator>
size_t GetResultSize(Iterator start, Iterator end, size_t s_size,
                     std::forward_iterator_tag) {
  size_t length = 0, num_elements = 0;
  for (Iterator it = start; it != end; ++it) {
    length += it->size();
    ++num_elements;
  }
  // Adds the size of all the separators
  return length + s_size * (num_elements - 1);
}

// A joining algorithm that's optimized for an iterator range of string-like
// objects that do not need any additional formatting. This is to optimize the
// common case of joining, say, a vector<string> or a vector<StringPiece>.
//
// This is an overload of the previous JoinAlgorithm() function. Here the
// Formatter argument is of type NoFormatter. Since NoFormatter is an internal
// type, this overload is only invoked when strings::Join() is called with a
// range of string-like objects (e.g., string, StringPiece), and an explicit
// Formatter argument was NOT specified.
//
// The optimization is that the needed space will be reserved in the output
// string to avoid the need to resize while appending. To do this, the iterator
// range will be traversed twice: once to calculate the total needed size, and
// then again to copy the elements and delimiters to the output string.
template <typename Iterator>
string JoinAlgorithm(Iterator start, Iterator end, StringPiece s, NoFormatter) {
  string result;
  if (start != end) {
    typename std::iterator_traits<Iterator>::iterator_category iterator_tag;
    result.reserve(GetResultSize(start, end, s.size(), iterator_tag));

    // Joins strings
    StringPiece sep("", 0);
    for (Iterator it = start; it != end; ++it) {
      result.append(sep.data(), sep.size());
      result.append(it->data(), it->size());
      sep = s;
    }
  }

  return result;
}

#ifdef LANG_CXX11
// JoinTupleLoop implements a loop over the elements of a std::tuple, which
// are heterogeneous. The primary template matches the tuple interior case. It
// continues the iteration after appending a separator (for nonzero indices)
// and formatting an element of the tuple. The specialization for the I=N case
// matches the end-of-tuple, and terminates the iteration.
template <typename Tup, size_t I, size_t N, typename Formatter>
struct JoinTupleLoop {
  void operator()(string* out, const Tup& tup, StringPiece sep,
                  Formatter fmt) const {
    if (I > 0) out->append(sep.data(), sep.size());
    fmt(out, std::get<I>(tup));
    JoinTupleLoop<Tup, I + 1, N, Formatter>()(out, tup, sep, fmt);
  }
};
template <typename Tup, size_t N, typename Formatter>
struct JoinTupleLoop<Tup, N, N, Formatter> {
  void operator()(string* out, const Tup& tup, StringPiece sep,
                  Formatter fmt) const {}
};

template <typename... T, typename Formatter>
string JoinAlgorithm(const std::tuple<T...>& tup, StringPiece sep,
                     Formatter fmt) {
  typedef typename std::tuple<T...> Tup;
  const size_t kTupSize = std::tuple_size<Tup>::value;
  string result;
  JoinTupleLoop<Tup, 0, kTupSize, Formatter>()(&result, tup, sep, fmt);
  return result;
}
#endif  // LANG_CXX11

template <typename Iterator>
string JoinRange(Iterator first, Iterator last, StringPiece separator) {
  // No formatter was explicitly given, so a default must be chosen.
  typedef typename std::iterator_traits<Iterator>::value_type ValueType;
  typedef typename DefaultFormatter<ValueType>::Type Formatter;
  return JoinAlgorithm(first, last, separator, Formatter());
}

// In C++11, use std::begin(c) and std::end(c).
// Pre-11, fall back to member functions c.begin() and c.end().
#ifdef LANG_CXX11
#define STRINGS_JOIN_RANGE_BEGIN(c) std::begin(c)
#define STRINGS_JOIN_RANGE_END(c) std::end(c)
#else
#define STRINGS_JOIN_RANGE_BEGIN(c) c.begin()
#define STRINGS_JOIN_RANGE_END(c) c.end()
#endif

template <typename Range, typename Formatter>
string JoinRange(const Range& range, StringPiece separator, Formatter fmt) {
  return JoinAlgorithm(STRINGS_JOIN_RANGE_BEGIN(range),
                       STRINGS_JOIN_RANGE_END(range),
                       separator, fmt);
}

template <typename Range>
string JoinRange(const Range& range, StringPiece separator) {
  return JoinRange(STRINGS_JOIN_RANGE_BEGIN(range),
                   STRINGS_JOIN_RANGE_END(range),
                   separator);
}

#undef STRINGS_JOIN_RANGE_BEGIN
#undef STRINGS_JOIN_RANGE_END

}  // namespace internal
}  // namespace strings

#endif  // STRINGS_JOIN_INTERNAL_H_