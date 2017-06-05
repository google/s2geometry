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

#include "s2/third_party/absl/strings/strip.h"

#include <cstddef>
#include <string>

#include "s2/third_party/absl/base/port.h"
#include "s2/third_party/absl/strings/ascii_ctype.h"

void StripWhitespace(string* str) {
  size_t str_length = str->length();

  // Strip off leading whitespace.
  size_t first = 0;
  while (first < str_length && ascii_isspace(str->at(first))) {
    ++first;
  }
  // If entire string is white space.
  if (first == str_length) {
    str->clear();
    return;
  }
  if (first > 0) {
    str->erase(0, first);
    str_length -= first;
  }

  // Strip off trailing whitespace.
  // Here, the string is not empty and the first character is not whitespace.
  size_t last = str_length - 1;
  while (ascii_isspace(str->at(last))) {
    --last;
  }
  if (last != (str_length - 1)) {
    str->erase(last + 1, string::npos);
  }
}
