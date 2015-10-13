// Copyright 2010 Google Inc. All Rights Reserved.
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

// Author: jyrki@google.com (Jyrki Alakuijala)
// Refactored from contributions of various authors in strings/strutil.h
//
// This file contains string processing functions related to
// uppercase, lowercase, etc.
//
// These functions are for ASCII only. If you need to process UTF8 strings,
// take a look at files in i18n/utf8.

#ifndef STRINGS_CASE_H_
#define STRINGS_CASE_H_

#include <string.h>
#ifndef _MSC_VER
#include <strings.h>  // for strcasecmp, but msvc does not have this header
#endif
#include <functional>
#include <string>

#include "strings/memutil.h"
#include "strings/stringpiece.h"

// These are function objects (this is the kind of thing that STL
// uses).  This provides a comparison function appropriate for
// char *s.  A common use is in hashtables:
//
//   hash_map<const char*, int, CStringCaseHash, strcaseeq> ht;
//
// Note that your hash function must also be case-insensitive.  CStringCaseHash
// is suitable, and defined in util/hash/case_insensitive_hash.h
//
// Case-insensitive hashing is not recommended for new code. See
// util/hash/case_insensitive_hash.h
struct strcaseeq : public std::binary_function<const char*, const char*, bool> {
  bool operator()(const char* s1, const char* s2) const {
    return ((s1 == 0 && s2 == 0) ||
            (s1 && s2 && strcasecmp(s1, s2) == 0));
  }
};

// For strcaselt, sorting would put NULL string last.
struct strcaselt : public std::binary_function<const char*, const char*, bool> {
  bool operator()(const char* s1, const char* s2) const {
    return (s1 != s2) && (s2 == 0 || (s1 != 0 && strcasecmp(s1, s2) < 0));
  }
};

// ----------------------------------------------------------------------
// GetCapitalization()
//    Return a value indicating whether the string is entirely
//    lowercase, entirely uppercase, first letter uppercase, or
//    mixed case.  As returned by ascii_islower() and ascii_isupper().
// ----------------------------------------------------------------------

enum CapsType {
  CAPS_LOWER,
  CAPS_UPPER,
  CAPS_FIRST,
  CAPS_MIXED,
  CAPS_NOALPHA,
};

CapsType GetCapitalization(const char* s);

// Case-insensitive string comparison, uses C/POSIX locale.
// Returns:
//    less than 0:    if s1 < s2
//    equal to 0:     if s1 == s2
//    greater than 0: if s1 > s2
inline int StringCaseCompare(const string& s1, const string& s2) {
  return strcasecmp(s1.c_str(), s2.c_str());
}

// Returns true if the two strings are equal, case-insensitively speaking.
// Uses C/POSIX locale.
inline bool StringCaseEqual(const string& s1, const string& s2) {
  return strcasecmp(s1.c_str(), s2.c_str()) == 0;
}

// Case-insensitive less-than string comparison, uses C/POSIX locale.
// Useful as a template parameter for STL set/map of strings, if uniqueness of
// keys is case-insensitive.
struct StringCaseLess {
  bool operator()(const string& s1, const string& s2) const {
    return strcasecmp(s1.c_str(), s2.c_str()) < 0;
  }
};

// Case-insensitive StringPiece comparison.
// Returns:
//    less than 0:    if s1 < s2
//    equal to 0:     if s1 == s2
//    greater than 0: if s1 > s2
int CaseCompare(StringPiece s1, StringPiece s2);

// Returns true if the two StringPieces are equal, case-insensitively speaking.
inline bool CaseEqual(StringPiece s1, StringPiece s2) {
  if (s1.size() != s2.size()) return false;
  return memcasecmp(s1.data(), s2.data(), s1.size()) == 0;
}

// Case-insensitive less-than StringPiece comparison.
// Useful as a template parameter for STL set/map of StringPiece-compatible
// types, if uniqueness of keys is case-insensitive.
struct CaseLess {
  bool operator()(StringPiece s1, StringPiece s2) const {
    return CaseCompare(s1, s2) < 0;
  }
};


// ----------------------------------------------------------------------
// LowerString()
// LowerStringToBuf()
//    Convert the characters in "s" to lowercase.
//    Works only with ASCII strings; for UTF8, see ToLower in
//    util/utf8/public/unilib.h
//    Changes contents of "s".  LowerStringToBuf copies at most
//    "n" characters (including the terminating '\0')  from "s"
//    to another buffer.
// ----------------------------------------------------------------------

void LowerString(char* s);
void LowerString(string* s);
void LowerStringToBuf(const char* s, char* buf, int n);

namespace strings {
inline string ToLower(StringPiece s) {
  string out(s.data(), s.size());
  LowerString(&out);
  return out;
}
}  // namespace strings

// ----------------------------------------------------------------------
// UpperString()
// UpperStringToBuf()
//    Convert the characters in "s" to uppercase.
//    Works only with ASCII strings; for UTF8, see ToUpper in
//    util/utf8/public/unilib.h
//    UpperString changes "s". UpperStringToBuf copies at most
//    "n" characters (including the terminating '\0')  from "s"
//    to another buffer.
// ----------------------------------------------------------------------

void UpperString(char* s);
void UpperString(string* s);
void UpperStringToBuf(const char* s, char* buf, int n);

namespace strings {
inline string ToUpper(StringPiece s) {
  string out(s.data(), s.size());
  UpperString(&out);
  return out;
}
}  // namespace strings

// ---------------------------------------------------------------------
// TitlecaseString()
//     Capitalize first character of each word in a string.
//     delimiters is a set of characters that can be used as word
//     boundaries.
//     This function can be implemented using regular expression,
//     but this version is supposed to be more efficient.
// ---------------------------------------------------------------------
void TitlecaseString(string* s, StringPiece delimiters);

#endif  // STRINGS_CASE_H_
