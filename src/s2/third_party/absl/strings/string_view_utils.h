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
// Utility functions for operating on absl::string_views
// Collected here for convenience
//

#ifndef S2_THIRD_PARTY_ABSL_STRINGS_STRING_VIEW_UTILS_H_
#define S2_THIRD_PARTY_ABSL_STRINGS_STRING_VIEW_UTILS_H_

#include <ctype.h>
#include <cstddef>
#include <vector>

#include "s2/third_party/absl/base/integral_types.h"
#include "s2/third_party/absl/strings/internal/fastmem.h"
#include "s2/third_party/absl/strings/match.h"
#include "s2/third_party/absl/strings/str_split.h"
#include "s2/third_party/absl/strings/string_view.h"

namespace absl {

// If "*s" starts with "expected", consume it and return true.
// Otherwise, return false.
inline bool ConsumePrefix(absl::string_view* s, absl::string_view expected) {
  if (!StartsWith(*s, expected))
    return false;
  s->remove_prefix(expected.size());
  return true;
}

// If "*s" ends with "expected", remove it and return true.
// Otherwise, return false.
inline bool ConsumeSuffix(absl::string_view* s, absl::string_view expected) {
  if (!EndsWith(*s, expected))
    return false;
  s->remove_suffix(expected.size());
  return true;
}

// Assigns 'src' to '*dest'.
inline void CopyToString(absl::string_view src, string* dest) {
  dest->assign(src.begin(), src.end());
}

}  // namespace absl

namespace strings {

using absl::StartsWith;
using absl::EndsWith;
using absl::ConsumePrefix;
using absl::ConsumeSuffix;
using absl::CopyToString;

inline bool Contains(absl::string_view s, absl::string_view x) {
  return absl::StrContains(s, x);
}

// Removes leading ascii_isspace() characters.
// Returns number of characters removed.
stringpiece_ssize_type RemoveLeadingWhitespace(absl::string_view* text);

// Removes trailing ascii_isspace() characters.
// Returns number of characters removed.
stringpiece_ssize_type RemoveTrailingWhitespace(absl::string_view* text);

// Removes leading and trailing ascii_isspace() chars.
// Returns number of chars removed.
stringpiece_ssize_type RemoveWhitespaceContext(absl::string_view* text);

// Removes all characters up to and including specified char.
// Returns number of characters removed.
stringpiece_ssize_type RemoveUntil(absl::string_view* text, char sentinel);

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
bool ConsumeLeadingDigits(absl::string_view* s, uint64* val);

// If *s starts with 'expected', consume it and return true.
// Otherwise, return false.
inline bool ConsumeLeadingChar(absl::string_view* s, char expected) {
  if (!s->empty() && (*s)[0] == expected) {
    s->remove_prefix(1);
    return true;
  } else {
    return false;
  }
}

// Checks if two stringpiece values are equal ignoring case.
bool EqualIgnoreCase(absl::string_view piece1, absl::string_view piece2);

// This is similar to gstrncasestr() in strutil.h, except that it works with
// absl::string_views. It acts the same as absl::string_view::find(), except
// that it is case insensitive.
stringpiece_ssize_type FindIgnoreCase(absl::string_view haystack,
                                      absl::string_view needle);

// Like ConsumePrefix, but case insensitive.
inline bool ConsumeCasePrefix(absl::string_view* s,
                              absl::string_view expected) {
  if (StartsWithIgnoreCase(*s, expected)) {
    s->remove_prefix(expected.size());
    return true;
  }
  return false;
}

// Like ConsumeSuffix, but case insensitive.
inline bool ConsumeCaseSuffix(absl::string_view* s,
                              absl::string_view expected) {
  if (EndsWithIgnoreCase(*s, expected)) {
    s->remove_suffix(expected.size());
    return true;
  }
  return false;
}

// Yields the longest prefix in common between both input strings.
// Pointer-wise, the returned result is a subset of input "a".
absl::string_view FindLongestCommonPrefix(absl::string_view a,
                                          absl::string_view b);

using absl::EndsWithIgnoreCase;
using absl::StartsWithIgnoreCase;

}  // namespace strings

// ----------------------------------------------------------------------
// StringPieceCaseHash
// StringPieceCaseEqual
//
// Function objects for case-insensitive hash_map from absl::string_view.  E.g.,
// hash_map<absl::string_view, int, StringPieceCaseHash, StringPieceCaseEqual>
// ht;
// ----------------------------------------------------------------------

struct StringPieceCaseHash {
  size_t operator()(absl::string_view sp) const;
};

struct StringPieceCaseEqual {
  bool operator()(absl::string_view piece1, absl::string_view piece2) const {
    return strings::EqualIgnoreCase(piece1, piece2);
  }
};

#endif  // S2_THIRD_PARTY_ABSL_STRINGS_STRING_VIEW_UTILS_H_
