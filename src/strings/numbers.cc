// Copyright 2016 Google Inc. All Rights Reserved.
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

#include "strings/numbers.h"

#include "base/strtoint.h"
#include "strings/ascii_ctype.h"

int HexDigitsPrefix(const char* buf, int num_digits) {
  for (int i = 0; i < num_digits; i++)
    if (!ascii_isxdigit(buf[i]))
      return 0;  // This also detects end of string as '\0' is not xdigit.
  return 1;
}

uint64 ParseLeadingHex64Value(const char *str, uint64 deflt) {
  char *error = nullptr;
  const uint64 value = strtou64(str, &error, 16);
  return (error == str) ? deflt : value;
}
