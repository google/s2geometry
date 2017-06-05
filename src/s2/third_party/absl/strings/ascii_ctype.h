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

// -----------------------------------------------------------------------------
// File: ascii_ctype.h
// -----------------------------------------------------------------------------
//
// This package contains character classification functions analogous to those
// found in the ANSI C Standard Library <ctype.h> header file.
//
// C++ implementations provide <ctype.h> functionality based on their
// C environment locale. In general, reliance on such a locale is not ideal, as
// the locale standard is problematic (and may not return invariant information
// for the same character set, for example). These `ascii_*()` functions are
// hard-wired for standard ASCII, much faster, and guaranteed to behave
// consistently.
//
// `ascii_isalnum()`, `ascii_isalpha()`, `ascii_isascii()`, `ascii_isblank()`,
// `ascii_iscntrl()`, `ascii_isdigit()`, `ascii_isgraph()`, `ascii_islower()`,
// `ascii_isprint()`, `ascii_ispunct()`, `ascii_isspace()`, `ascii_isupper()`,
// `ascii_isxdigit()`
//   Analagous to the <ctype.h> functions with similar names, these
//   functions take an unsigned char and return a bool, based on whether the
//   character matches the condition specified.
//
//   If the input character has a numerical value greater than 127, these
//   functions return `false`.
//
// `ascii_tolower()`, `ascii_toupper()`
//   Analagous to the <ctype.h> functions with similar names, these functions
//   take an unsigned char and return a char.
//
//   If the input character is not an ASCII {lower,upper}-case letter (including
//   numerical values greater than 127) then the functions return the same value
//   as the input character.

#ifndef S2_THIRD_PARTY_ABSL_STRINGS_ASCII_CTYPE_H_
#define S2_THIRD_PARTY_ABSL_STRINGS_ASCII_CTYPE_H_

namespace absl {
namespace internal {

// Declaration for an array of bitfields holding character information.
extern const unsigned char kAsciiPropertyBits[256];

// Declaration for the array of characters to upper-case characters.
extern const char kAsciiToUpper[256];

// Declaration for the array of characters to lower-case characters.
extern const char kAsciiToLower[256];

}  // namespace internal
}  // namespace absl

// Public functions.

// -----------------------------------------------------------------------------
// ascii_isalpha()
// -----------------------------------------------------------------------------
//
// Determines whether the given character is an alphabetic character.
static inline bool ascii_isalpha(unsigned char c) {
  return (absl::internal::kAsciiPropertyBits[c] & 0x01) != 0;
}

// -----------------------------------------------------------------------------
// ascii_isalnum()
// -----------------------------------------------------------------------------
//
// Determines whether the given character is an alphanumeric character.
static inline bool ascii_isalnum(unsigned char c) {
  return (absl::internal::kAsciiPropertyBits[c] & 0x04) != 0;
}

// -----------------------------------------------------------------------------
// ascii_isspace()
// -----------------------------------------------------------------------------
//
// Determines whether the given character is a whitespace character (space,
// tab, vertical tab, formfeed, linefeed, or carriage return).
static inline bool ascii_isspace(unsigned char c) {
  return (absl::internal::kAsciiPropertyBits[c] & 0x08) != 0;
}

// -----------------------------------------------------------------------------
// ascii_ispunct()
// -----------------------------------------------------------------------------
//
// Determines whether the given character is a punctuation character.
static inline bool ascii_ispunct(unsigned char c) {
  return (absl::internal::kAsciiPropertyBits[c] & 0x10) != 0;
}

// -----------------------------------------------------------------------------
// ascii_isblank()
// -----------------------------------------------------------------------------
//
// Determines whether the given character is a blank character (tab or space).
static inline bool ascii_isblank(unsigned char c) {
  return (absl::internal::kAsciiPropertyBits[c] & 0x20) != 0;
}

// -----------------------------------------------------------------------------
// ascii_iscntrl()
// -----------------------------------------------------------------------------
//
// Determines whether the given character is a control character.
static inline bool ascii_iscntrl(unsigned char c) {
  return (absl::internal::kAsciiPropertyBits[c] & 0x40) != 0;
}

// -----------------------------------------------------------------------------
// ascii_isxdigit()
// -----------------------------------------------------------------------------
//
// Determines whether the given character can be represented as a hexadecimal
// digit character (i.e. {0-9} or {A-F}).
static inline bool ascii_isxdigit(unsigned char c) {
  return (absl::internal::kAsciiPropertyBits[c] & 0x80) != 0;
}

// -----------------------------------------------------------------------------
// ascii_isdigit()
// -----------------------------------------------------------------------------
//
// Determines whether the given character can be represented as a decimal
// digit character (i.e. {0-9}).
static inline bool ascii_isdigit(unsigned char c) {
  return c >= '0' && c <= '9';
}

// -----------------------------------------------------------------------------
// ascii_isprint()
// -----------------------------------------------------------------------------
//
// Determines whether the given character is printable, including whitespace.
static inline bool ascii_isprint(unsigned char c) {
  return c >= 32 && c < 127;
}

// -----------------------------------------------------------------------------
// ascii_isgraph()
// -----------------------------------------------------------------------------
//
// Determines whether the given character has a graphical representation.
static inline bool ascii_isgraph(unsigned char c) {
  return c >  32 && c < 127;
}

// -----------------------------------------------------------------------------
// ascii_isupper()
// -----------------------------------------------------------------------------
//
// Determines whether the given character is uppercase.
static inline bool ascii_isupper(unsigned char c) {
  return c >= 'A' && c <= 'Z';
}

// -----------------------------------------------------------------------------
// ascii_islower()
// -----------------------------------------------------------------------------
//
// Determines whether the given character is lowercase.
static inline bool ascii_islower(unsigned char c) {
  return c >= 'a' && c <= 'z';
}

// -----------------------------------------------------------------------------
// ascii_isascii()
// -----------------------------------------------------------------------------
//
// Determines whether the given character is ASCII.
static inline bool ascii_isascii(unsigned char c) {
  return c < 128;
}

// -----------------------------------------------------------------------------
// ascii_tolower()
// -----------------------------------------------------------------------------
//
// Returns an ASCII character, converting to lowercase if uppercase is
// passed. Note that character values > 127 are simply returned.
static inline char ascii_tolower(unsigned char c) {
  return absl::internal::kAsciiToLower[c];
}

// -----------------------------------------------------------------------------
// ascii_toupper()
// -----------------------------------------------------------------------------
//
// Returns the ASCII character, converting to upper-case if lower-case is
// passed. Note that characters values > 127 are simply returned.
static inline char ascii_toupper(unsigned char c) {
  return absl::internal::kAsciiToUpper[c];
}

#endif  // S2_THIRD_PARTY_ABSL_STRINGS_ASCII_CTYPE_H_
