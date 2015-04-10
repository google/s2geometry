// Copyright 2003 Google Inc. All Rights Reserved.
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
// Utility functions for operating on StringPieces
// Collected here for convenience
//

#ifndef STRINGS_STRINGPIECE_UTILS_H_
#define STRINGS_STRINGPIECE_UTILS_H_

#include <ctype.h>
#include <stddef.h>
#include <vector>

#include "base/integral_types.h"
#include "strings/split.h"
#include "strings/stringpiece.h"
#include "util/gtl/charmap.h"

namespace strings {

// Removes leading ascii_isspace() characters.
// Returns number of characters removed.
stringpiece_ssize_type RemoveLeadingWhitespace(StringPiece* text);

// Removes trailing ascii_isspace() characters.
// Returns number of characters removed.
stringpiece_ssize_type RemoveTrailingWhitespace(StringPiece* text);

// Removes leading and trailing ascii_isspace() chars.
// Returns number of chars removed.
stringpiece_ssize_type RemoveWhitespaceContext(StringPiece* text);

// Removes all characters up to and including specified char.
// Returns number of characters removed.
stringpiece_ssize_type RemoveUntil(StringPiece* text, char sentinel);

// Consume a leading positive integer value.  If any digits
// were found, store the value of the leading unsigned number in
// "*val", advance "*s" past the consumed number, and return true.
// If overflow occurred, returns false.
// Otherwise, returns false.
// Equivalent to RE::Consume(s, "(\\d+)", val) but significantly
// faster:
//
//   Run on panther (4 X 2200 MHz CPUs); 2008/10/21-09:20:57
//   CPU: AMD Opteron Engineering Sample (4 cores) dL1:64KB dL2:1024KB
//   Benchmark                     Time(ns)    CPU(ns) Iterations
//   ------------------------------------------------------------
//   BM_ConsumeDigits                    20         20   38049268
//   BM_ConsumeDigitsWithRE             420        416    1678915
//
// (i.e. 20 ns vs. 420 ns (even with the regexp pre-created)).
bool ConsumeLeadingDigits(StringPiece* s, uint64* val);

// If *s starts with 'expected', consume it and return true.
// Otherwise, return false.
inline bool ConsumeLeadingChar(StringPiece* s, char expected) {
  if (!s->empty() && (*s)[0] == expected) {
    s->remove_prefix(1);
    return true;
  } else {
    return false;
  }
}

// If "*s" starts with "expected", consume it and return true.
// Otherwise, return false.
inline bool ConsumePrefix(StringPiece* s, StringPiece expected) {
  if (s->starts_with(expected)) {
    s->remove_prefix(expected.size());
    return true;
  }
  return false;
}

// If "*s" ends with "expected", remove it and return true.
// Otherwise, return false.
inline bool ConsumeSuffix(StringPiece* s, StringPiece expected) {
  if (s->ends_with(expected)) {
    s->remove_suffix(expected.size());
    return true;
  }
  return false;
}

// Consume a leading component of "*s" whose first character
// is in "first_char_set" and is followed by zero or more occurrences of
// characters in "rest_char_set".  Returns true if we matched
// the pattern, sets "*result" to the consumed data (if "result" is not NULL),
// and removes the consumed characters from the front of "*s".  Otherwise,
// returns false (and "*s" is unmodified).
inline bool ConsumeComponent(const Charmap& first_char_set,
                             const Charmap& rest_char_set,
                             StringPiece* s,
                             StringPiece* result) {
  const char* p = s->data();
  const char* limit = p + s->size();
  if (p == limit || !first_char_set.contains(*p)) {
    return false;
  }
  p++;
  while (p < limit && rest_char_set.contains(*p)) {
    p++;
  }
  // Consume some digits
  const stringpiece_ssize_type N = p - s->data();
  if (result) {
    result->set(s->data(), N);
  }
  s->remove_prefix(N);
  return true;
}


// Checks if two stringpiece values are equal ignoring case.
bool EqualIgnoreCase(StringPiece piece1,
                     StringPiece piece2);

// Returns true if "text" starts with "starts_with". The comparison ignores
// case.
bool StartsWithIgnoreCase(StringPiece text,
                          StringPiece starts_with);

// Returns true if "text" ends with "ends_with". The comparison ignores
// case.
bool EndsWithIgnoreCase(StringPiece text, StringPiece ends_with);

// This is similar to gstrncasestr() in strutil.h, except that it works with
// StringPieces.  It acts the same as StringPiece::find(), except that it is
// case insensitive.
stringpiece_ssize_type FindIgnoreCase(StringPiece haystack, StringPiece needle);

// Like ConsumePrefix, but case insensitive.
inline bool ConsumeCasePrefix(StringPiece* s, StringPiece expected) {
  if (StartsWithIgnoreCase(*s, expected)) {
    s->remove_prefix(expected.size());
    return true;
  }
  return false;
}

// Yields the longest prefix in common between both input strings.
// Pointer-wise, the returned result is a subset of input "a".
StringPiece FindLongestCommonPrefix(StringPiece a, StringPiece b);

// strto32 and strto64
int32 ParseInt32Prefix(StringPiece str, size_t* len, int radix);
int64 ParseInt64Prefix(StringPiece str, size_t* len, int radix);

}  // namespace strings

// ----------------------------------------------------------------------
// StringPieceCaseHash
// StringPieceCaseEqual
//
// Function objects for case-insensitive hash_map from StringPiece.  E.g.,
// hash_map<StringPiece, int, StringPieceCaseHash, StringPieceCaseEqual> ht;
// ----------------------------------------------------------------------

struct StringPieceCaseHash {
  size_t operator()(StringPiece sp) const;
};

struct StringPieceCaseEqual {
  bool operator()(StringPiece piece1, StringPiece piece2) const {
    return strings::EqualIgnoreCase(piece1, piece2);
  }
};

#endif  // STRINGS_STRINGPIECE_UTILS_H_