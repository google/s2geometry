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



#include "s2/third_party/absl/strings/string_view_utils.h"

#include <algorithm>
#include <string>
#include <vector>

#include "s2/third_party/absl/strings/ascii.h"
#include "s2/third_party/absl/strings/internal/memutil.h"

namespace strings {

stringpiece_ssize_type RemoveLeadingWhitespace(absl::string_view* text) {
  stringpiece_ssize_type count = 0;
  const char* ptr = text->data();
  while (count < text->size() && absl::ascii_isspace(*ptr)) {
    count++;
    ptr++;
  }
  text->remove_prefix(count);
  return count;
}

stringpiece_ssize_type RemoveTrailingWhitespace(absl::string_view* text) {
  stringpiece_ssize_type count = 0;
  const char* ptr = text->data() + text->size() - 1;
  while (count < text->size() && absl::ascii_isspace(*ptr)) {
    ++count;
    --ptr;
  }
  text->remove_suffix(count);
  return count;
}

stringpiece_ssize_type RemoveWhitespaceContext(absl::string_view* text) {
  // use RemoveLeadingWhitespace() and RemoveTrailingWhitespace() to do the job
  return (RemoveLeadingWhitespace(text) + RemoveTrailingWhitespace(text));
}

stringpiece_ssize_type RemoveUntil(absl::string_view* text, char sentinel) {
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

bool ConsumeLeadingDigits(absl::string_view* s, uint64* val) {
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

bool EqualIgnoreCase(absl::string_view piece1, absl::string_view piece2) {
  return (piece1.size() == piece2.size() &&
          0 == absl::strings_internal::memcasecmp(piece1.data(), piece2.data(),
                                                  piece1.size()));
  // memcasecmp uses absl::ascii_tolower().
}

stringpiece_ssize_type FindIgnoreCase(absl::string_view haystack,
                                      absl::string_view needle) {
  // We use the cursor to iterate through the haystack...on each
  // iteration the cursor is moved forward one character.
  absl::string_view cursor = haystack;
  while (cursor.size() >= needle.size()) {
    if (StartsWithIgnoreCase(cursor, needle)) {
      return cursor.data() - haystack.data();
    }
    cursor.remove_prefix(1);
  }
  return absl::string_view::npos;
}

absl::string_view FindLongestCommonPrefix(absl::string_view a,
                                          absl::string_view b) {
  if (a.empty() || b.empty()) return absl::string_view();

  const char* pa = a.data();
  const char* pb = b.data();
  stringpiece_ssize_type count = 0;
  const stringpiece_ssize_type limit = std::min(a.size(), b.size());
  while (count < limit && *pa == *pb) {
    ++pa;
    ++pb;
    ++count;
  }

  return absl::string_view(a.data(), count);
}

}  // namespace strings

// ----------------------------------------------------------------------
// StringPieceCaseHash
// ----------------------------------------------------------------------

size_t StringPieceCaseHash::operator()(absl::string_view sp) const {
  // based on __stl_string_hash in http://www.sgi.com/tech/stl/string
  size_t hash_val = 0;
  for (const char c : sp) {
    hash_val = 5 * hash_val + absl::ascii_tolower(c);
  }
  return hash_val;
}
