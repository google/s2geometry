//
// Copyright 2017 The Abseil Authors.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// -----------------------------------------------------------------------------
// File: strip.h
// -----------------------------------------------------------------------------
//
// This file contains various functions for stripping substrings from a string.
#ifndef S2_THIRD_PARTY_ABSL_STRINGS_STRIP_H_
#define S2_THIRD_PARTY_ABSL_STRINGS_STRIP_H_

#include <cstddef>
#include <string>

#include "s2/third_party/absl/base/macros.h"
#include "s2/third_party/absl/strings/ascii.h"
#include "s2/third_party/absl/strings/ascii_ctype.h"
#include "s2/third_party/absl/strings/match.h"
#include "s2/third_party/absl/strings/string_view.h"

namespace absl {

// ConsumePrefix()
//
// Strips the `expected` prefix from the start of the given string, returning
// `true` if the strip operation succeeded or false otherwise.
//
// Example:
//
//   absl::string_view input("abc");
//   EXPECT_TRUE(absl::ConsumePrefix(&input, "a"));
//   EXPECT_EQ(input, "bc");
inline bool ConsumePrefix(absl::string_view* str, absl::string_view expected) {
  if (!absl::StartsWith(*str, expected)) return false;
  str->remove_prefix(expected.size());
  return true;
}
// ConsumeSuffix()
//
// Strips the `expected` suffix from the end of the given string, returning
// `true` if the strip operation succeeded or false otherwise.
//
// Example:
//
//   absl::string_view input("abcdef");
//   EXPECT_TRUE(absl::ConsumeSuffix(&input, "def"));
//   EXPECT_EQ(input, "abc");
inline bool ConsumeSuffix(absl::string_view* str, absl::string_view expected) {
  if (!absl::EndsWith(*str, expected)) return false;
  str->remove_suffix(expected.size());
  return true;
}

// StripPrefix()
//
// Returns a view into the input string 'str' with the given 'prefix' removed,
// but leaving the original string intact. If the prefix does not match at the
// start of the string, returns the original string instead.
inline absl::string_view StripPrefix(absl::string_view str,
                                     absl::string_view prefix) {
  if (absl::StartsWith(str, prefix)) str.remove_prefix(prefix.size());
  return str;
}

// StripSuffix()
//
// Returns a view into the input string 'str' with the given 'suffix' removed,
// but leaving the original string intact. If the suffix does not match at the
// end of the string, returns the original string instead.
inline absl::string_view StripSuffix(absl::string_view str,
                                     absl::string_view suffix) {
  if (absl::EndsWith(str, suffix)) str.remove_suffix(suffix.size());
  return str;
}

}  // namespace absl


// Returns a copy of the input string 'str' with the given 'prefix' removed. If
// the prefix doesn't match, returns a copy of the original string.
//
// The "Try" version stores the stripped string in the 'result' out-param and
// returns true iff the prefix was found and removed. It is safe for 'result' to
// point back to the input string.
//
// See also absl::ConsumePrefix().
inline string StripPrefixString(absl::string_view str,
                                absl::string_view prefix) {
  return string(absl::StripPrefix(str, prefix));
}
inline bool TryStripPrefixString(absl::string_view str,
                                 absl::string_view prefix, string* result) {
  bool res = absl::ConsumePrefix(&str, prefix);
  result->assign(str.begin(), str.end());
  return res;
}

// Returns a copy of the input string 'str' with the given 'suffix' removed. If
// the suffix doesn't match, returns a copy of the original string.
//
// The "Try" version stores the stripped string in the 'result' out-param and
// returns true iff the suffix was found and removed. It is safe for 'result' to
// point back to the input string.
//
// See also absl::ConsumeSuffix().
inline string StripSuffixString(absl::string_view str,
                                absl::string_view suffix) {
  return string(absl::StripSuffix(str, suffix));
}
inline bool TryStripSuffixString(absl::string_view str,
                                 absl::string_view suffix, string* result) {
  bool res = absl::ConsumeSuffix(&str, suffix);
  result->assign(str.begin(), str.end());
  return res;
}

// Replaces any of the characters in 'remove' with the character 'replace_with'.
//
void ReplaceCharacters(char* str, size_t len, absl::string_view remove,
                       char replace_with);
void ReplaceCharacters(string* s, absl::string_view remove, char replace_with);

// Replaces the character 'remove' with the character 'replace_with'.
//
inline void ReplaceCharacter(char* str, size_t len, char remove,
                             char replace_with) {
  for (char* end = str + len; str != end; ++str) {
    if (*str == remove) *str = replace_with;
  }
}

// Replaces runs of one or more 'dup_char' with a single occurrence, and returns
// the number of characters that were removed.
//
// Example:
//       StripDupCharacters("a//b/c//d", '/', 0) => "a/b/c/d"
ptrdiff_t StripDupCharacters(string* s, char dup_char, ptrdiff_t start_pos);

// Removes whitespace from both ends of the given string. This function has
// various overloads for passing the input string as a C-string + length,
// string, or absl::string_view. If the caller is using NUL-terminated strings,
// it is the caller's responsibility to insert the NUL character at the end of
// the substring.
ABSL_DEPRECATED("Use absl::StripAsciiWhitespace() instead")
inline void StripWhitespace(const char** str, ptrdiff_t* len) {
  auto stripped = absl::StripAsciiWhitespace(absl::string_view(*str, *len));
  *len = stripped.size();
  if (*len) *str = &*stripped.begin();
}
ABSL_DEPRECATED("Use absl::StripAsciiWhitespace() instead")
inline void StripWhitespace(char** str, ptrdiff_t* len) {
  auto stripped = absl::StripAsciiWhitespace(absl::string_view(*str, *len));
  *len = stripped.size();
  if (*len) *str = const_cast<char*>(&*stripped.begin());
}
ABSL_DEPRECATED("Use absl::StripAsciiWhitespace() instead")
inline void StripWhitespace(string* str) { absl::StripAsciiWhitespace(str); }

ABSL_DEPRECATED("Use absl::StripAsciiWhitespace() instead")
inline void StripWhitespace(absl::string_view* str) {
  *str = absl::StripAsciiWhitespace(*str);
}

namespace strings {

// Calls StripWhitespace() on each element in the given collection.
//
// Note: this implementation is conceptually similar to
//
//   std::for_each(c.begin(), c.end(), StripWhitespace);
//
// except that StripWhitespace requires a *pointer* to the element, so the above
// std::for_each solution wouldn't work.
template <typename Collection>
inline void StripWhitespaceInCollection(Collection* collection) {
  for (typename Collection::iterator it = collection->begin();
       it != collection->end(); ++it)
    StripWhitespace(&(*it));
}

}  // namespace strings

// Removes whitespace from the beginning of the given string. The version that
// takes a string modifies the given input string.
//
// The versions that take C-strings return a pointer to the first non-whitespace
// character if one is present or nullptr otherwise. 'line' must be
// NUL-terminated.
inline void StripLeadingWhitespace(string* str) {
  absl::StripLeadingAsciiWhitespace(str);
}

inline const char* StripLeadingWhitespace(const char* line) {
  auto stripped = absl::StripLeadingAsciiWhitespace(line);

  if (stripped.empty()) return nullptr;

  return stripped.begin();
}
// StripLeadingWhitespace for non-const strings.
inline char* StripLeadingWhitespace(char* line) {
  auto stripped = absl::StripLeadingAsciiWhitespace(line);

  if (stripped.empty()) return nullptr;

  return const_cast<char*>(stripped.begin());
}

// Removes whitespace from the end of the given string.
inline void StripTrailingWhitespace(string* s) {
  absl::StripTrailingAsciiWhitespace(s);
}

// Removes the trailing '\n' or '\r\n' from 's', if one exists. Returns true if
// a newline was found and removed.
bool StripTrailingNewline(string* s);

// Removes leading, trailing, and duplicate internal whitespace.
inline void RemoveExtraWhitespace(string* s) {
  absl::RemoveExtraAsciiWhitespace(s);
}

// Returns a pointer to the first non-whitespace character in 'str'. Never
// returns nullptr. 'str' must be NUL-terminated.
inline const char* SkipLeadingWhitespace(const char* str) {
  while (absl::ascii_isspace(*str)) ++str;
  return str;
}
inline char* SkipLeadingWhitespace(char* str) {
  while (absl::ascii_isspace(*str)) ++str;
  return str;
}

// Strips everything enclosed in pairs of curly braces ('{' and '}') and the
// curly braces themselves. Doesn't touch open braces without a closing brace.
// Does not handle nesting.
void StripCurlyBraces(string* s);

// Performs the same operation as StripCurlyBraces, but allows the caller to
// specify different left and right bracket characters, such as '(' and ')'.
void StripBrackets(char left, char right, string* s);

// Strips everything between a right angle bracket ('<') and left angle bracket
// ('>') including the brackets themselves, e.g.
// "the quick <b>brown</b> fox" --> "the quick brown fox".
//
// This does not understand HTML nor does it know anything about HTML tags or
// comments. This is simply a text processing function that removes text between
// pairs of angle brackets. Note that in the example above the word "brown" is
// not removed because it is not between pairs of angle brackets.
//
// This is NOT safe for security and this will NOT prevent against XSS.
//
// For a more full-featured HTML parser, see //webutil/pageutil/pageutil.h.
void StripMarkupTags(string* s);
string OutputWithMarkupTagsStripped(const string& s);

// Removes any occurrences of the characters in 'remove' from the:
//
//   - start of the string "Left"
//   - end of the string "Right"
//   - both ends of the string
//
// Returns the number of chars removed.
ptrdiff_t TrimStringLeft(string* s, absl::string_view remove);
ptrdiff_t TrimStringRight(string* s, absl::string_view remove);
inline ptrdiff_t TrimString(string* s, absl::string_view remove) {
  return TrimStringRight(s, remove) + TrimStringLeft(s, remove);
}
ptrdiff_t TrimStringLeft(absl::string_view* s, absl::string_view remove);
ptrdiff_t TrimStringRight(absl::string_view* s, absl::string_view remove);
inline ptrdiff_t TrimString(absl::string_view* s, absl::string_view remove) {
  return TrimStringRight(s, remove) + TrimStringLeft(s, remove);
}

// Removes leading and trailing runs, and collapses middle runs of a set of
// characters into a single character (the first one specified in 'remove').
// E.g.: TrimRunsInString(&s, " :,()") removes leading and trailing delimiter
// chars and collapses and converts internal runs of delimiters to single ' '
// characters, so, for example, "  a:(b):c  " -> "a b c".
void TrimRunsInString(string* s, absl::string_view remove);

// Removes all internal '\0' characters from the string.
void RemoveNullsInString(string* s);

// Removes all occurrences of the given character from the given string. Returns
// the new length.
ptrdiff_t strrm(char* str, char c);
ptrdiff_t memrm(char* str, ptrdiff_t strlen, char c);

// Removes all occurrences of any character from 'chars' from the given string.
// Returns the new length.
ptrdiff_t strrmm(char* str, const char* chars);
ptrdiff_t strrmm(string* str, const string& chars);

template <typename T>
ABSL_DEPRECATED("stop passing int*")
inline typename std::enable_if<std::is_same<T, int>::value, void>::type
StripWhitespace(const char** str, T* len) {
  ptrdiff_t pdt = *len;
  StripWhitespace(str, &pdt);
  *len = static_cast<int>(pdt);
}
template <typename T>
ABSL_DEPRECATED("stop passing int*")
inline typename std::enable_if<std::is_same<T, int>::value, void>::type
StripWhitespace(char** str, T* len) {
  ptrdiff_t pdt = *len;
  StripWhitespace(str, &pdt);
  *len = static_cast<int>(pdt);
}


#endif  // S2_THIRD_PARTY_ABSL_STRINGS_STRIP_H_
