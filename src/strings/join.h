// Copyright 2008 Google Inc. All Rights Reserved.
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
// #status: RECOMMENDED
// #category: operations on strings
// #summary: Functions for joining ranges of elements with an element separator.
//
#ifndef STRINGS_JOIN_H_
#define STRINGS_JOIN_H_

#include <stdio.h>
#include <string.h>

#include <ext/hash_map>
using __gnu_cxx::hash;
using __gnu_cxx::hash_map;  // Not used in this file.
#include <ext/hash_set>
using __gnu_cxx::hash;
using __gnu_cxx::hash_set;  // Not used in this file.
#include <iterator>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "base/integral_types.h"
#include "base/macros.h"
#include "base/port.h"
#include "base/template_util.h"
#include "strings/join_internal.h"
#include "strings/numbers.h"
#include "strings/strcat.h"    // For backward compatibility.
#include "strings/stringpiece.h"
#include "util/hash/hash.h"

#ifdef LANG_CXX11
#include <initializer_list>  // NOLINT(build/include_order)
#include <tuple>  // NOLINT(build/include_order)
#endif  // LANG_CXX11

namespace strings {

//                              strings::Join()
//
// The strings::Join() function joins the given range of elements, with each
// element separated by the given separator string, and returns the result as a
// string. Ranges may be specified by passing a container or array compatible
// with std::begin() and std::end(), a brace-initialized std::initializer_list,
// as individual begin and end iterators, or as a std::tuple of heterogeneous
// objects. The separator string is taken as a StringPiece, which means it can
// be specified as a string literal, C-string, C++ string, etc. By default,
// non-string elements are converted to strings using AlphaNum, which yields
// the same behavior as using StrCat(). This means that strings::Join() works
// out-of-the-box on collections of strings, ints, floats, doubles, etc.  An
// optional final argument of a "Formatter" (details below) function object
// may be given. This object will be responsible for converting each
// argument in the Range to a string.
//
// Example 1:
//   // Joins a collection of strings. This also works with a collection of
//   // StringPiece or even const char*.
//   vector<string> v = util::gtl::Container("foo", "bar", "baz");
//   string s = strings::Join(v, "-");
//   EXPECT_EQ("foo-bar-baz", s);
//
// Example 2:
//   // Joins the values in the given std::initializer_list<> specified using
//   // brace initialization. This also works with an initializer_list of ints
//   // or StringPiece--any AlphaNum-compatible type.
//   string s = strings::Join({"foo", "bar", "baz"}, "-");
//   EXPECT_EQ("foo-bar-baz", s);
//
// Example 3:
//   // Joins a collection of ints. This also works with floats, doubles,
//   // int64s; any StrCat-compatible type.
//   vector<int> v = util::gtl::Container(1, 2, 3, -4);
//   string s = strings::Join(v, "-");
//   EXPECT_EQ("1-2-3--4", s);
//
// Example 4:
//   // Joins a collection of pointer-to-int. By default, pointers are
//   // dereferenced and the pointee is formatted using the default format for
//   // that type. The upshot of this is that all levels of inderection are
//   // collapsed, so this works just as well for vector<int***> as
//   // vector<int*>.
//   int x = 1, y = 2, z = 3;
//   vector<int*> v = util::gtl::Container(&x, &y, &z);
//   string s = strings::Join(v, "-");
//   EXPECT_EQ("1-2-3", s);
//
// Example 5:
//   //   In C++11, dereferecing std::unique_ptr is also supported:
//   vector<std::unique_ptr<int>> v
//   v.emplace_back(new int(1));
//   v.emplace_back(new int(2));
//   v.emplace_back(new int(3));
//   string s = strings::Join(v, "-");
//   EXPECT_EQ("1-2-3", s);
//
// Example 6:
//   // Joins a map, with each key-value pair separated by an equals sign.
//   // This would also work with, say, a vector<pair<>>.
//    map<string, int> m = util::gtl::Container(
//        std::make_pair("a", 1),
//        std::make_pair("b", 2),
//        std::make_pair("c", 3));
//   string s = strings::Join(m, ",", strings::PairFormatter("="));
//   EXPECT_EQ("a=1,b=2,c=3", s);
//
// Example 7:
//   // These examples show how strings::Join() handles a few common edge cases.
//   vector<string> v_empty;
//   EXPECT_EQ("", strings::Join(v_empty, "-"));
//
//   vector<string> v_one_item = util::gtl::Container("foo");
//   EXPECT_EQ("foo", strings::Join(v_one_item, "-"));
//
//   vector<string> v_empty_string = util::gtl::Container("");
//   EXPECT_EQ("", strings::Join(v_empty_string, "-"));
//
//   vector<string> v_one_item_empty_string = util::gtl::Container("a", "");
//   EXPECT_EQ("a-", strings::Join(v_one_item_empty_string, "-"));
//
//   vector<string> v_two_empty_string = util::gtl::Container("", "");
//   EXPECT_EQ("-", strings::Join(v_two_empty_string, "-"));
//
// Example 8:
//   // Join a std::tuple<T...>.
//   string s = strings::Join(std::make_tuple(123, "abc", 0.456), "-");
//   EXPECT_EQ("123-abc-0.456", s);
//

//
// Formatters
//
// A Formatter is a function object that is responsible for formatting its
// argument as a string and appending it to the given output string. Formatters
// are an extensible part of the Join2 API: They allow callers to provide their
// own conversion function to enable strings::Join() work with arbitrary types.
//
// The following is an example Formatter that simply uses StrAppend to format an
// integer as a string.
//
//   struct MyFormatter {
//     void operator()(string* out, int i) const {
//       StrAppend(out, i);
//     }
//   };
//
// You would use the above formatter by passing an instance of it as the final
// argument to strings::Join():
//
//   vector<int> v = util::gtl::Container(1, 2, 3, 4);
//   string s = strings::Join(v, "-", MyFormatter());
//   EXPECT_EQ("1-2-3-4", s);
//
// The following standard formatters are provided with the Join2 API.
//
// - AlphaNumFormatter (the default)
// - PairFormatter
// - DereferenceFormatter
//

// AlphaNumFormatter()
//
// Default formatter used if none is specified. Uses AlphaNum to convert numeric
// arguments to strings.
inline internal::AlphaNumFormatterImpl AlphaNumFormatter() {
  return internal::AlphaNumFormatterImpl();
}

// PairFormatter()
//
// Formats a std::pair by putting the given separator between the pair's .first
// and .second members. The separator argument is required. By default, the
// first and second members are themselves formatted using AlphaNumFormatter(),
// but the caller may specify other formatters to use for the members.
template <typename FirstFormatter, typename SecondFormatter>
inline internal::PairFormatterImpl<FirstFormatter, SecondFormatter>
PairFormatter(FirstFormatter f1, StringPiece sep, SecondFormatter f2) {
  return internal::PairFormatterImpl<FirstFormatter, SecondFormatter>(
      f1, sep, f2);
}
inline internal::PairFormatterImpl<
  internal::AlphaNumFormatterImpl,
  internal::AlphaNumFormatterImpl>
PairFormatter(StringPiece sep) {
  return PairFormatter(AlphaNumFormatter(), sep, AlphaNumFormatter());
}

// DereferenceFormatter()
//
// Dereferences its argument then formats it using AlphaNumFormatter (by
// default), or the given formatter if one is explicitly given. This is useful
// for formatting a container of pointer-to-T. This pattern often shows up when
// joining repeated fields in protocol buffers.
template <typename Formatter>
internal::DereferenceFormatterImpl<Formatter>
DereferenceFormatter(Formatter f) {
  return internal::DereferenceFormatterImpl<Formatter>(f);
}
inline internal::DereferenceFormatterImpl<internal::AlphaNumFormatterImpl>
DereferenceFormatter() {
  return internal::DereferenceFormatterImpl<internal::AlphaNumFormatterImpl>(
      AlphaNumFormatter());
}

//
// strings::Join() overloads
//

template <typename Iterator, typename Formatter>
string Join(Iterator start, Iterator end, StringPiece sep, Formatter fmt) {
  return internal::JoinAlgorithm(start, end, sep, fmt);
}

template <typename Range, typename Formatter>
string Join(const Range& range, StringPiece separator, Formatter fmt) {
  return internal::JoinRange(range, separator, fmt);
}

#ifdef LANG_CXX11
template <typename T, typename Formatter>
string Join(std::initializer_list<T> il, StringPiece separator, Formatter fmt) {
  return internal::JoinRange(il, separator, fmt);
}

template <typename... T, typename Formatter>
string Join(const std::tuple<T...>& value, StringPiece separator,
            Formatter fmt) {
  return internal::JoinAlgorithm(value, separator, fmt);
}
#endif  // LANG_CXX11

template <typename Iterator>
string Join(Iterator start, Iterator end, StringPiece separator) {
  return internal::JoinRange(start, end, separator);
}

template <typename Range>
string Join(const Range& range, StringPiece separator) {
  return internal::JoinRange(range, separator);
}

#ifdef LANG_CXX11
template <typename T>
string Join(std::initializer_list<T> il, StringPiece separator) {
  return internal::JoinRange(il, separator);
}

template <typename... T>
string Join(const std::tuple<T...>& value, StringPiece separator) {
  return internal::JoinAlgorithm(value, separator, AlphaNumFormatter());
}
#endif  // LANG_CXX11

}  // namespace strings

// ----------------------------------------------------------------------
// LEGACY(jgm): Utilities provided in util/csv/writer.h are now preferred for
//
// Example for CSV formatting a single record (a sequence container of string,
// char*, or StringPiece values) using the util::csv::WriteRecordToString helper
// function:
//   std::vector<string> record = ...;
//   string line = util::csv::WriteRecordToString(record);
//
// NOTE: When writing many records, use the util::csv::Writer class directly.
//
// JoinCSVLineWithDelimiter()
//    This function is the inverse of SplitCSVLineWithDelimiter() in that the
//    string returned by JoinCSVLineWithDelimiter() can be passed to
//    SplitCSVLineWithDelimiter() to get the original string vector back.
//    Quotes and escapes the elements of original_cols according to CSV quoting
//    rules, and the joins the escaped quoted strings with commas using
//    JoinStrings().  Note that JoinCSVLineWithDelimiter() will not necessarily
//    return the same string originally passed in to
//    SplitCSVLineWithDelimiter(), since SplitCSVLineWithDelimiter() can handle
//    gratuitous spacing and quoting. 'output' must point to an empty string.
//
//    Example:
//     [Google], [x], [Buchheit, Paul], [string with " quote in it], [ space ]
//     --->  [Google,x,"Buchheit, Paul","string with "" quote in it"," space "]
//
// JoinCSVLine()
//    A convenience wrapper around JoinCSVLineWithDelimiter which uses
//    ',' as the delimiter.
// ----------------------------------------------------------------------
void JoinCSVLine(const std::vector<string>& original_cols, string* output);
string JoinCSVLine(const std::vector<string>& original_cols);
void JoinCSVLineWithDelimiter(const std::vector<string>& original_cols,
                              char delimiter,
                              string* output);

template <class CONTAINER>
GOOGLE_DEPRECATED("Use strings::Join()")
void JoinStrings(const CONTAINER& components,
                 StringPiece delim,
                 string* result) {
  *result = strings::Join(components, delim);
}

template <class ITERATOR>
GOOGLE_DEPRECATED("Use strings::Join()")
void JoinStringsIterator(const ITERATOR& start,
                         const ITERATOR& end,
                         StringPiece delim,
                         string* result) {
  *result = strings::Join(start, end, delim);
}

#endif  // STRINGS_JOIN_H_