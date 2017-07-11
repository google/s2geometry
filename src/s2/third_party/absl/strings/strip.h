// Copyright 2011 Google Inc. All Rights Reserved.
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
//
// This file contains various functions for "stripping" aka removing various
// characters and substrings from a string. Much of the code in here is old and
// operates on C-style strings such as char* and const char*. Prefer the
// interfaces that take and return absl::string_view and C++ string objects. See
// //third_party/absl/strings/string_view_utils.h for similar functions that
// operate on absl::string_view.

#ifndef S2_THIRD_PARTY_ABSL_STRINGS_STRIP_H_
#define S2_THIRD_PARTY_ABSL_STRINGS_STRIP_H_

#include <cstddef>

#include <string>

#include "s2/third_party/absl/base/macros.h"
#include "s2/third_party/absl/strings/ascii_ctype.h"
#include "s2/third_party/absl/strings/string_view.h"

// Returns a copy of the input string 'str' with the given 'prefix' removed. If
// the prefix doesn't match, returns a copy of the original string.
//
// The "Try" version stores the stripped string in the 'result' out-param and
// returns true iff the prefix was found and removed. It is safe for 'result' to
// point back to the input string.
//
// See also absl::ConsumePrefix().
string StripPrefixString(absl::string_view str, absl::string_view prefix);
bool TryStripPrefixString(absl::string_view str, absl::string_view prefix,
                          string* result);

// Returns a copy of the input string 'str' with the given 'suffix' removed. If
// the suffix doesn't match, returns a copy of the original string.
//
// The "Try" version stores the stripped string in the 'result' out-param and
// returns true iff the suffix was found and removed. It is safe for 'result' to
// point back to the input string.
//
// See also absl::ConsumeSuffix().
string StripSuffixString(absl::string_view str, absl::string_view suffix);
bool TryStripSuffixString(absl::string_view str, absl::string_view suffix,
                          string* result);

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

ABSL_DEPRECATED("Use ReplaceCharacters()")
inline void StripString(char* str, absl::string_view remove,
                        char replace_with) {
  ReplaceCharacters(str, strlen(str), remove, replace_with);
}
ABSL_DEPRECATED("Use ReplaceCharacters()")
inline void StripString(char* str, size_t len, absl::string_view remove,
                        char replace_with) {
  ReplaceCharacters(str, len, remove, replace_with);
}
ABSL_DEPRECATED("Use ReplaceCharacters()")
inline void StripString(string* s, absl::string_view remove,
                        char replace_with) {
  ReplaceCharacters(s, remove, replace_with);
}
ABSL_DEPRECATED("Use ReplaceCharacter()")
inline void StripString(char* str, char remove, char replace_with) {
  ReplaceCharacter(str, strlen(str), remove, replace_with);
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
void StripWhitespace(const char** str, ptrdiff_t* len);
inline void StripWhitespace(char** str, ptrdiff_t* len) {
  // The "real" type for StripWhitespace is ForAll char types C, take
  // (C, int) as input and return (C, int) as output.  We're using the
  // cast here to assert that we can take a char*, even though the
  // function thinks it's assigning to const char*.
  StripWhitespace(const_cast<const char**>(str), len);
}
void StripWhitespace(string* str);
inline void StripWhitespace(absl::string_view* str) {
  const char* data = str->data();
  ptrdiff_t len = str->size();
  StripWhitespace(&data, &len);
  *str = absl::string_view(data, len);
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
void StripLeadingWhitespace(string* str);
inline const char* StripLeadingWhitespace(const char* line) {
  // skip leading whitespace
  while (ascii_isspace(*line))
    ++line;

  if ('\0' == *line)  // end of line, no non-whitespace
    return nullptr;

  return line;
}
// StripLeadingWhitespace for non-const strings.
inline char* StripLeadingWhitespace(char* line) {
  return const_cast<char*>(
      StripLeadingWhitespace(const_cast<const char*>(line)));
}

// Removes whitespace from the end of the given string.
void StripTrailingWhitespace(string* s);

// Removes the trailing '\n' or '\r\n' from 's', if one exists. Returns true if
// a newline was found and removed.
bool StripTrailingNewline(string* s);

// Removes leading, trailing, and duplicate internal whitespace.
void RemoveExtraWhitespace(string* s);

// Returns a pointer to the first non-whitespace character in 'str'. Never
// returns nullptr. 'str' must be NUL-terminated.
inline const char* SkipLeadingWhitespace(const char* str) {
  while (ascii_isspace(*str))
    ++str;
  return str;
}
inline char* SkipLeadingWhitespace(char* str) {
  while (ascii_isspace(*str))
    ++str;
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
