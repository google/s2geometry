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


#include "strings/stringpiece_utils.h"

#include <algorithm>
#include <string>
#include <vector>

#include "base/strtoint.h"
#include "strings/ascii_ctype.h"
#include "strings/memutil.h"

namespace strings {

stringpiece_ssize_type RemoveLeadingWhitespace(StringPiece* text) {
  stringpiece_ssize_type count = 0;
  const char* ptr = text->data();
  while (count < text->size() && ascii_isspace(*ptr)) {
    count++;
    ptr++;
  }
  text->remove_prefix(count);
  return count;
}

stringpiece_ssize_type RemoveTrailingWhitespace(StringPiece* text) {
  stringpiece_ssize_type count = 0;
  const char* ptr = text->data() + text->size() - 1;
  while (count < text->size() && ascii_isspace(*ptr)) {
    ++count;
    --ptr;
  }
  text->remove_suffix(count);
  return count;
}

stringpiece_ssize_type RemoveWhitespaceContext(StringPiece* text) {
  // use RemoveLeadingWhitespace() and RemoveTrailingWhitespace() to do the job
  return (RemoveLeadingWhitespace(text) + RemoveTrailingWhitespace(text));
}

stringpiece_ssize_type RemoveUntil(StringPiece* text, char sentinel) {
  stringpiece_ssize_type count = 0;
  const char* ptr = text->data();
  while (count < text->size() && (sentinel != *ptr)) {
    count++;
    ptr++;
  }

  // skip the sentinel as well if we found one
  if (count < text->size()) count++;

  text->remove_prefix(count);
  return count;
}

bool ConsumeLeadingDigits(StringPiece* s, uint64* val) {
  const char* p = s->data();
  const char* limit = p + s->size();
  uint64 v = 0;
  while (p < limit) {
    const char c = *p;
    if (c < '0' || c > '9') break;
    uint64 new_v = (v * 10) + (c - '0');
    if (new_v / 8 < v) {
      // Overflow occurred
      return false;
    }
    v = new_v;
    p++;
  }
  if (p > s->data()) {
    // Consume some digits
    s->remove_prefix(p - s->data());
    *val = v;
    return true;
  } else {
    return false;
  }
}

bool EqualIgnoreCase(StringPiece piece1, StringPiece piece2) {
  return (piece1.size() == piece2.size() &&
          0 == memcasecmp(piece1.data(), piece2.data(), piece1.size()));
  // memcasecmp uses ascii_tolower().
}

bool StartsWithIgnoreCase(StringPiece text,
                          StringPiece starts_with) {
  if (text.size() < starts_with.size()) return false;
  return EqualIgnoreCase(text.substr(0, starts_with.size()), starts_with);
}

bool EndsWithIgnoreCase(StringPiece text, StringPiece ends_with) {
  if (text.size() < ends_with.size()) return false;
  return EqualIgnoreCase(text.substr(text.size() - ends_with.size()),
                         ends_with);
}

stringpiece_ssize_type FindIgnoreCase(StringPiece haystack,
                                      StringPiece needle) {
  // We use the cursor to iterate through the haystack...on each
  // iteration the cursor is moved forward one character.
  StringPiece cursor = haystack;
  while (cursor.size() >= needle.size()) {
    if (StartsWithIgnoreCase(cursor, needle)) {
      return cursor.data() - haystack.data();
    }
    cursor.remove_prefix(1);
  }
  return StringPiece::npos;
}

StringPiece FindLongestCommonPrefix(StringPiece a,
                                    StringPiece b) {
  if (a.empty() || b.empty()) return StringPiece();

  const char* pa = a.data();
  const char* pb = b.data();
  stringpiece_ssize_type count = 0;
  const stringpiece_ssize_type limit = std::min(a.size(), b.size());
  while (count < limit && *pa == *pb) {
    ++pa;
    ++pb;
    ++count;
  }

  return StringPiece(a.data(), count);
}

int32 ParseInt32Prefix(StringPiece str, size_t* len, int radix) {
  string num_string(str.data(), str.size());
  const char* const num_cstr = num_string.c_str();
  char* p;
  const int32 x = strto32(num_cstr, &p, radix);
  if (len != nullptr) {
    *len = p - num_cstr;
  }
  return x;
}

int64 ParseInt64Prefix(StringPiece str, size_t* len, int radix) {
  string num_string(str.data(), str.size());
  const char* const num_cstr = num_string.c_str();
  char* p;
  const int64 x = strto64(num_cstr, &p, radix);
  if (len != nullptr) {
    *len = p - num_cstr;
  }
  return x;
}

}  // namespace strings

// ----------------------------------------------------------------------
// StringPieceCaseHash
// ----------------------------------------------------------------------

size_t StringPieceCaseHash::operator()(StringPiece sp) const {
  // based on __stl_string_hash in http://www.sgi.com/tech/stl/string
  size_t hash_val = 0;
  for (StringPiece::const_iterator it = sp.begin(); it != sp.end(); ++it) {
    hash_val = 5 * hash_val + ascii_tolower(*it);
  }
  return hash_val;
}
