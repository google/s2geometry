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
// Architecture-neutral plug compatible replacements for strtol() friends.
//
// Long's have different lengths on ILP-32 and LP-64 platforms, and so overflow
// behavior across the two varies when strtol() and similar are used to parse
// 32-bit integers.  Similar problems exist with atoi(), because although it
// has an all-integer interface, it uses strtol() internally, and so suffers
// from the same narrowing problems on assignments to int.
//
// Examples:
//   errno = 0;
//   i = strtol("3147483647", nullptr, 10);
//   printf("%d, errno %d\n", i, errno);
//   //   32-bit platform: 2147483647, errno 34
//   //   64-bit platform: -1147483649, errno 0
//
//   printf("%d\n", atoi("3147483647"));
//   //   32-bit platform: 2147483647
//   //   64-bit platform: -1147483649
//
// A way round this is to define local replacements for these, and use them
// instead of the standard libc functions.
//
// In most 32-bit cases the replacements can be inlined away to a call to the
// libc function.  In a couple of 64-bit cases, however, adapters are required,
// to provide the right overflow and errno behavior.
//

#ifndef S2_BASE_STRTOINT_H_
#define S2_BASE_STRTOINT_H_

#include <cstdlib> // For strtol* functions.
#include <string>
#include <cstdint>

#include "s2/base/port.h"
#include "s2/third_party/absl/base/macros.h"

// Adapter functions for handling overflow and errno.
std::int32_t strto32_adapter(const char *nptr, char **endptr, int base);
std::uint32_t strtou32_adapter(const char *nptr, char **endptr, int base);

// Conversions to a 32-bit integer can pass the call to strto[u]l on 32-bit
// platforms, but need a little extra work on 64-bit platforms.
inline std::int32_t strto32(const char *nptr, char **endptr, int base) {
  if (sizeof(std::int32_t) == sizeof(long))
    return static_cast<std::int32_t>(strtol(nptr, endptr, base));
  else
    return strto32_adapter(nptr, endptr, base);
}

inline std::uint32_t strtou32(const char *nptr, char **endptr, int base) {
  if (sizeof(std::uint32_t) == sizeof(unsigned long))
    return static_cast<std::uint32_t>(strtoul(nptr, endptr, base));
  else
    return strtou32_adapter(nptr, endptr, base);
}

// For now, long long is 64-bit on all the platforms we care about, so these
// functions can simply pass the call to strto[u]ll.
inline std::int64_t strto64(const char *nptr, char **endptr, int base) {
  static_assert(sizeof(std::int64_t) == sizeof(long long),
                "sizeof std::int64_t is not sizeof long long");
  return strtoll(nptr, endptr, base);
}

inline std::uint64_t strtou64(const char *nptr, char **endptr, int base) {
  static_assert(sizeof(std::uint64_t) == sizeof(unsigned long long),
                "sizeof std::uint64_t is not sizeof long long");
  return strtoull(nptr, endptr, base);
}

// Although it returns an int, atoi() is implemented in terms of strtol, and
// so has differing overflow and underflow behavior.  atol is the same.
inline std::int32_t atoi32(const char *nptr) {
  return strto32(nptr, nullptr, 10);
}

inline std::int64_t atoi64(const char *nptr) {
  return strto64(nptr, nullptr, 10);
}

// Convenience versions of the above that take a string argument.
inline std::int32_t atoi32(const string &s) {
  return atoi32(s.c_str());
}

inline std::int64_t atoi64(const string &s) {
  return atoi64(s.c_str());
}

#endif  // S2_BASE_STRTOINT_H_
