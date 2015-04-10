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
// #summary: Merges strings or numbers with no delimiter.
//
#ifndef STRINGS_STRCAT_H_
#define STRINGS_STRCAT_H_

#include <string>

#include "base/integral_types.h"
#include "base/macros.h"
#include "base/port.h"
#include "strings/numbers.h"
#include "strings/stringpiece.h"

// The AlphaNum type was designed to be used as the parameter type for StrCat().
// Any routine accepting either a string or a number may accept it.
// The basic idea is that by accepting a "const AlphaNum &" as an argument
// to your function, your callers will automagically convert bools, integers,
// and floating point values to strings for you.
//
// NOTE: Use of AlphaNum outside of the //strings package is unsupported except
// for the specific case of function parameters of type "AlphaNum" or "const
// AlphaNum &". In particular, instantiating AlphaNum directly as a stack
// variable is not supported.
//
// Conversion from 8-bit values is not accepted because if it were, then an
// attempt to pass ':' instead of ":" might result in a 58 ending up in your
// result.
//
// Bools convert to "0" or "1".
//
// Floating point values are converted to a string which, if passed to strtod(),
// would produce the exact same original double (except in case of NaN; all NaNs
// are considered the same value). We try to keep the string short but it's not
// guaranteed to be as short as possible.
//
// You can convert to Hexadecimal output rather than Decimal output using Hex.
// To do this, pass strings::Hex(my_int) as a parameter to StrCat. You may
// specify a minimum field width using a separate parameter, so the equivalent
// of StringPrintf("%04x", my_int) is StrCat(Hex(my_int, Hex::ZERO_PAD_4))
//
// This class has implicit constructors.
// Style guide exception granted:
// http://goto/style-guide-exception-20978288
//
namespace strings {

struct Hex {
  uint64 value;
  enum PadSpec {
    NONE = 1,
    ZERO_PAD_2,
    ZERO_PAD_3,
    ZERO_PAD_4,
    ZERO_PAD_5,
    ZERO_PAD_6,
    ZERO_PAD_7,
    ZERO_PAD_8,
    ZERO_PAD_9,
    ZERO_PAD_10,
    ZERO_PAD_11,
    ZERO_PAD_12,
    ZERO_PAD_13,
    ZERO_PAD_14,
    ZERO_PAD_15,
    ZERO_PAD_16,
  } spec;
  template <class Int>
  explicit Hex(Int v, PadSpec s = NONE)
      : spec(s) {
    // Prevent sign-extension by casting integers to
    // their unsigned counterparts.
    static_assert(
        sizeof(v) == 1 || sizeof(v) == 2 || sizeof(v) == 4 || sizeof(v) == 8,
        "Unknown integer type");
    value = sizeof(v) == 1 ? static_cast<uint8>(v)
          : sizeof(v) == 2 ? static_cast<uint16>(v)
          : sizeof(v) == 4 ? static_cast<uint32>(v)
          : static_cast<uint64>(v);
  }
};

class AlphaNum {
 public:
  // No bool ctor -- bools convert to an integral type.
  // A bool ctor would also convert incoming pointers (bletch).

  AlphaNum(int32 i32)  // NOLINT(runtime/explicit)
      : piece_(digits_, FastInt32ToBufferLeft(i32, digits_) - &digits_[0]) {}
  AlphaNum(uint32 u32)  // NOLINT(runtime/explicit)
      : piece_(digits_, FastUInt32ToBufferLeft(u32, digits_) - &digits_[0]) {}
  AlphaNum(int64 i64)  // NOLINT(runtime/explicit)
      : piece_(digits_, FastInt64ToBufferLeft(i64, digits_) - &digits_[0]) {}
  AlphaNum(uint64 u64)  // NOLINT(runtime/explicit)
      : piece_(digits_, FastUInt64ToBufferLeft(u64, digits_) - &digits_[0]) {}

#ifdef _LP64
  AlphaNum(long x)  // NOLINT(runtime/explicit)
    : piece_(digits_, FastInt64ToBufferLeft(x, digits_) - &digits_[0]) {}
  AlphaNum(unsigned long x)  // NOLINT(runtime/explicit)
    : piece_(digits_, FastUInt64ToBufferLeft(x, digits_) - &digits_[0]) {}
#else
  AlphaNum(long x)  // NOLINT(runtime/explicit)
    : piece_(digits_, FastInt32ToBufferLeft(x, digits_) - &digits_[0]) {}
  AlphaNum(unsigned long x)  // NOLINT(runtime/explicit)
    : piece_(digits_, FastUInt32ToBufferLeft(x, digits_) - &digits_[0]) {}
#endif

  AlphaNum(float f)  // NOLINT(runtime/explicit)
    : piece_(digits_, strlen(FloatToBuffer(f, digits_))) {}
  AlphaNum(double f)  // NOLINT(runtime/explicit)
    : piece_(digits_, strlen(DoubleToBuffer(f, digits_))) {}

  AlphaNum(Hex hex);  // NOLINT(runtime/explicit)

  AlphaNum(const char *c_str) : piece_(c_str) {}   // NOLINT(runtime/explicit)
  AlphaNum(const StringPiece &pc) : piece_(pc) {}  // NOLINT(runtime/explicit)

#if defined(HAS_GLOBAL_STRING)
  template <class Allocator>
  AlphaNum(const basic_string<char, std::char_traits<char>,
                              Allocator> &str)
      : piece_(str) {}
#endif
  template <class Allocator>
  AlphaNum(const std::basic_string<char, std::char_traits<char>,
                                   Allocator> &str)  // NOLINT(runtime/explicit)
      : piece_(str) {}


  StringPiece::size_type size() const { return piece_.size(); }
  const char *data() const { return piece_.data(); }
  StringPiece Piece() const { return piece_; }

 private:
  StringPiece piece_;
  char digits_[kFastToBufferSize];

  // Use ":" not ':'
  AlphaNum(char c);  // NOLINT(runtime/explicit)

  DISALLOW_COPY_AND_ASSIGN(AlphaNum);
};

extern AlphaNum gEmptyAlphaNum;

}  // namespace strings

using strings::AlphaNum;
using strings::gEmptyAlphaNum;

// ----------------------------------------------------------------------
// StrCat()
//    This merges the given strings or numbers, with no delimiter.  This
//    is designed to be the fastest possible way to construct a string out
//    of a mix of raw C strings, StringPieces, strings, bool values,
//    and numeric values.
//
//    Don't use this for user-visible strings.  The localization process
//    works poorly on strings built up out of fragments.
//
//    For clarity and performance, don't use StrCat when appending to a
//    string.  In particular, avoid using any of these (anti-)patterns:
//      str.append(StrCat(...))
//      str += StrCat(...)
//      str = StrCat(str, ...)
//    where the last is the worse, with the potential to change a loop
//    from a linear time operation with O(1) dynamic allocations into a
//    quadratic time operation with O(n) dynamic allocations.  StrAppend
//    is a better choice than any of the above, subject to the restriction
//    of StrAppend(&str, a, b, c, ...) that none of the a, b, c, ... may
//    be a reference into str.
// ----------------------------------------------------------------------

string StrCat(const AlphaNum &a) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c)
    MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f)
    MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i)
    MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k, const AlphaNum &l)
    MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
              const AlphaNum &m) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
              const AlphaNum &m, const AlphaNum &n) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
              const AlphaNum &m, const AlphaNum &n, const AlphaNum &o)
    MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
              const AlphaNum &m, const AlphaNum &n, const AlphaNum &o,
              const AlphaNum &p) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
              const AlphaNum &m, const AlphaNum &n, const AlphaNum &o,
              const AlphaNum &p, const AlphaNum &q) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
              const AlphaNum &m, const AlphaNum &n, const AlphaNum &o,
              const AlphaNum &p, const AlphaNum &q, const AlphaNum &r)
    MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
              const AlphaNum &m, const AlphaNum &n, const AlphaNum &o,
              const AlphaNum &p, const AlphaNum &q, const AlphaNum &r,
              const AlphaNum &s) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
              const AlphaNum &m, const AlphaNum &n, const AlphaNum &o,
              const AlphaNum &p, const AlphaNum &q, const AlphaNum &r,
              const AlphaNum &s, const AlphaNum &t) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
              const AlphaNum &m, const AlphaNum &n, const AlphaNum &o,
              const AlphaNum &p, const AlphaNum &q, const AlphaNum &r,
              const AlphaNum &s, const AlphaNum &t, const AlphaNum &u)
    MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
              const AlphaNum &m, const AlphaNum &n, const AlphaNum &o,
              const AlphaNum &p, const AlphaNum &q, const AlphaNum &r,
              const AlphaNum &s, const AlphaNum &t, const AlphaNum &u,
              const AlphaNum &v) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
              const AlphaNum &m, const AlphaNum &n, const AlphaNum &o,
              const AlphaNum &p, const AlphaNum &q, const AlphaNum &r,
              const AlphaNum &s, const AlphaNum &t, const AlphaNum &u,
              const AlphaNum &v, const AlphaNum &w) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
              const AlphaNum &m, const AlphaNum &n, const AlphaNum &o,
              const AlphaNum &p, const AlphaNum &q, const AlphaNum &r,
              const AlphaNum &s, const AlphaNum &t, const AlphaNum &u,
              const AlphaNum &v, const AlphaNum &w, const AlphaNum &x)
    MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
              const AlphaNum &m, const AlphaNum &n, const AlphaNum &o,
              const AlphaNum &p, const AlphaNum &q, const AlphaNum &r,
              const AlphaNum &s, const AlphaNum &t, const AlphaNum &u,
              const AlphaNum &v, const AlphaNum &w, const AlphaNum &x,
              const AlphaNum &y) MUST_USE_RESULT;
string StrCat(const AlphaNum &a, const AlphaNum &b, const AlphaNum &c,
              const AlphaNum &d, const AlphaNum &e, const AlphaNum &f,
              const AlphaNum &g, const AlphaNum &h, const AlphaNum &i,
              const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
              const AlphaNum &m, const AlphaNum &n, const AlphaNum &o,
              const AlphaNum &p, const AlphaNum &q, const AlphaNum &r,
              const AlphaNum &s, const AlphaNum &t, const AlphaNum &u,
              const AlphaNum &v, const AlphaNum &w, const AlphaNum &x,
              const AlphaNum &y, const AlphaNum &z) MUST_USE_RESULT;

inline string StrCat(const AlphaNum &a) { return string(a.data(), a.size()); }

namespace strings {
namespace internal {

// Do not call directly - this is not part of the public API.
string StrCatNineOrMore(const AlphaNum *a1, ...);

}  // namespace internal
}  // namespace strings

// Support 9 or more arguments
inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k,
    const AlphaNum &l) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, &l, null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
    const AlphaNum &m) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, &l, &m, null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
    const AlphaNum &m, const AlphaNum &n) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, &l, &m, &n, null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
    const AlphaNum &m, const AlphaNum &n, const AlphaNum &o) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, &l, &m, &n, &o,
                                             null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
    const AlphaNum &m, const AlphaNum &n, const AlphaNum &o,
    const AlphaNum &p) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, &l, &m, &n, &o, &p,
                                             null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
    const AlphaNum &m, const AlphaNum &n, const AlphaNum &o, const AlphaNum &p,
    const AlphaNum &q) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, &l, &m, &n, &o, &p, &q,
                                             null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
    const AlphaNum &m, const AlphaNum &n, const AlphaNum &o, const AlphaNum &p,
    const AlphaNum &q, const AlphaNum &r) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, &l, &m, &n, &o, &p, &q, &r,
                                             null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
    const AlphaNum &m, const AlphaNum &n, const AlphaNum &o, const AlphaNum &p,
    const AlphaNum &q, const AlphaNum &r, const AlphaNum &s) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, &l, &m, &n, &o, &p, &q, &r,
                                             &s, null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
    const AlphaNum &m, const AlphaNum &n, const AlphaNum &o, const AlphaNum &p,
    const AlphaNum &q, const AlphaNum &r, const AlphaNum &s,
    const AlphaNum &t) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, &l, &m, &n, &o, &p, &q, &r,
                                             &s, &t, null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
    const AlphaNum &m, const AlphaNum &n, const AlphaNum &o, const AlphaNum &p,
    const AlphaNum &q, const AlphaNum &r, const AlphaNum &s, const AlphaNum &t,
    const AlphaNum &u) {
    const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, &l, &m, &n, &o, &p, &q, &r,
                                             &s, &t, &u, null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
    const AlphaNum &m, const AlphaNum &n, const AlphaNum &o, const AlphaNum &p,
    const AlphaNum &q, const AlphaNum &r, const AlphaNum &s, const AlphaNum &t,
    const AlphaNum &u, const AlphaNum &v) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, &l, &m, &n, &o, &p, &q, &r,
                                             &s, &t, &u, &v, null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
    const AlphaNum &m, const AlphaNum &n, const AlphaNum &o, const AlphaNum &p,
    const AlphaNum &q, const AlphaNum &r, const AlphaNum &s, const AlphaNum &t,
    const AlphaNum &u, const AlphaNum &v, const AlphaNum &w) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, &l, &m, &n, &o, &p, &q, &r,
                                             &s, &t, &u, &v, &w, null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
    const AlphaNum &m, const AlphaNum &n, const AlphaNum &o, const AlphaNum &p,
    const AlphaNum &q, const AlphaNum &r, const AlphaNum &s, const AlphaNum &t,
    const AlphaNum &u, const AlphaNum &v, const AlphaNum &w,
    const AlphaNum &x) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, &l, &m, &n, &o, &p, &q, &r,
                                             &s, &t, &u, &v, &w, &x,
                                             null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
    const AlphaNum &m, const AlphaNum &n, const AlphaNum &o, const AlphaNum &p,
    const AlphaNum &q, const AlphaNum &r, const AlphaNum &s, const AlphaNum &t,
    const AlphaNum &u, const AlphaNum &v, const AlphaNum &w, const AlphaNum &x,
    const AlphaNum &y) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, &l, &m, &n, &o, &p, &q, &r,
                                             &s, &t, &u, &v, &w, &x, &y,
                                             null_alphanum);
}

inline string StrCat(
    const AlphaNum &a, const AlphaNum &b, const AlphaNum &c, const AlphaNum &d,
    const AlphaNum &e, const AlphaNum &f, const AlphaNum &g, const AlphaNum &h,
    const AlphaNum &i, const AlphaNum &j, const AlphaNum &k, const AlphaNum &l,
    const AlphaNum &m, const AlphaNum &n, const AlphaNum &o, const AlphaNum &p,
    const AlphaNum &q, const AlphaNum &r, const AlphaNum &s, const AlphaNum &t,
    const AlphaNum &u, const AlphaNum &v, const AlphaNum &w, const AlphaNum &x,
    const AlphaNum &y, const AlphaNum &z) {
  const AlphaNum* null_alphanum = NULL;
  return strings::internal::StrCatNineOrMore(&a, &b, &c, &d, &e, &f, &g, &h, &i,
                                             &j, &k, &l, &m, &n, &o, &p, &q, &r,
                                             &s, &t, &u, &v, &w, &x, &y, &z,
                                             null_alphanum);
}

// ----------------------------------------------------------------------
// StrAppend()
//    Same as above, but adds the output to the given string.
//    WARNING: For speed, StrAppend does not try to check each of its input
//    arguments to be sure that they are not a subset of the string being
//    appended to.  That is, while this will work:
//
//    string s = "foo";
//    s += s;
//
//    This will not (necessarily) work:
//
//    string s = "foo";
//    StrAppend(&s, s);
//
//    Note: while StrCat supports appending up to 26 arguments, StrAppend
//    is currently limited to 9.  That's rarely an issue except when
//    automatically transforming StrCat to StrAppend, and can easily be
//    worked around as consecutive calls to StrAppend are quite efficient.
// ----------------------------------------------------------------------

void StrAppend(string *dest,      const AlphaNum &a);
void StrAppend(string *dest,      const AlphaNum &a, const AlphaNum &b);
void StrAppend(string *dest,      const AlphaNum &a, const AlphaNum &b,
               const AlphaNum &c);
void StrAppend(string *dest,      const AlphaNum &a, const AlphaNum &b,
               const AlphaNum &c, const AlphaNum &d);

// Support up to 9 params by using a default empty AlphaNum.
void StrAppend(string *dest,      const AlphaNum &a, const AlphaNum &b,
               const AlphaNum &c, const AlphaNum &d, const AlphaNum &e,
               const AlphaNum &f = gEmptyAlphaNum,
               const AlphaNum &g = gEmptyAlphaNum,
               const AlphaNum &h = gEmptyAlphaNum,
               const AlphaNum &i = gEmptyAlphaNum);

#endif  // STRINGS_STRCAT_H_