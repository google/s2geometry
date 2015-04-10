// Copyright 1999 Google Inc. All Rights Reserved.
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
//

// TODO(user): visit each const_cast.  Some of them are no longer necessary
// because last Single Unix Spec and grte v2 are more const-y.

#include "strings/util.h"

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <time.h>           // for FastTimeToBuffer()
#include <algorithm>
#include <string>
#include <vector>

#include <glog/logging.h>

#include "strings/ascii_ctype.h"
#include "strings/numbers.h"
#include "strings/stringpiece.h"
#include "util/gtl/stl_util.h"  // for STLAppendToString

#ifdef OS_WINDOWS
#ifdef min  // windows.h defines this to something silly
#undef min
#endif
#endif

using std::vector;

// Use this instead of gmtime_r if you want to build for Windows.
// Windows doesn't have a 'gmtime_r', but it has the similar 'gmtime_s'.
// TODO(user): Probably belongs in //base:time_support.{cc|h}.
static struct tm* PortableSafeGmtime(const time_t* timep, struct tm* result) {
#ifdef OS_WINDOWS
  return gmtime_s(result, timep) == 0 ? result : NULL;
#else
  return gmtime_r(timep, result);
#endif  // OS_WINDOWS
}

char* strnstr(const char* haystack, const char* needle,
                     size_t haystack_len) {
  if (*needle == '\0') {
    return const_cast<char*>(haystack);
  }
  size_t needle_len = strlen(needle);
  char* where;
  while ((where = strnchr(haystack, *needle, haystack_len)) != NULL) {
    if (where - haystack + needle_len > haystack_len) {
      return NULL;
    }
    if (strncmp(where, needle, needle_len) == 0) {
      return where;
    }
    haystack_len -= where + 1 - haystack;
    haystack = where + 1;
  }
  return NULL;
}

const char* strnprefix(const char* haystack, int haystack_size,
                              const char* needle, int needle_size) {
  if (needle_size > haystack_size) {
    return NULL;
  } else {
    if (strncmp(haystack, needle, needle_size) == 0) {
      return haystack + needle_size;
    } else {
      return NULL;
    }
  }
}

const char* strncaseprefix(const char* haystack, int haystack_size,
                                  const char* needle, int needle_size) {
  if (needle_size > haystack_size) {
    return NULL;
  } else {
    if (strncasecmp(haystack, needle, needle_size) == 0) {
      return haystack + needle_size;
    } else {
      return NULL;
    }
  }
}

char* strcasesuffix(char* str, const char* suffix) {
  const int lenstr = strlen(str);
  const int lensuffix = strlen(suffix);
  char* strbeginningoftheend = str + lenstr - lensuffix;

  if (lenstr >= lensuffix && 0 == strcasecmp(strbeginningoftheend, suffix)) {
    return (strbeginningoftheend);
  } else {
    return (NULL);
  }
}

const char* strnsuffix(const char* haystack, int haystack_size,
                              const char* needle, int needle_size) {
  if (needle_size > haystack_size) {
    return NULL;
  } else {
    const char* start = haystack + haystack_size - needle_size;
    if (strncmp(start, needle, needle_size) == 0) {
      return start;
    } else {
      return NULL;
    }
  }
}

const char* strncasesuffix(const char* haystack, int haystack_size,
                           const char* needle, int needle_size) {
  if (needle_size > haystack_size) {
    return NULL;
  } else {
    const char* start = haystack + haystack_size - needle_size;
    if (strncasecmp(start, needle, needle_size) == 0) {
      return start;
    } else {
      return NULL;
    }
  }
}

char* strchrnth(const char* str, const char& c, int n) {
  if (str == NULL)
    return NULL;
  if (n <= 0)
    return const_cast<char*>(str);
  const char* sp;
  int k = 0;
  for (sp = str; *sp != '\0'; sp ++) {
    if (*sp == c) {
      ++k;
      if (k >= n)
        break;
    }
  }
  return (k < n) ? NULL : const_cast<char*>(sp);
}

char* AdjustedLastPos(const char* str, char separator, int n) {
  if ( str == NULL )
    return NULL;
  const char* pos = NULL;
  if ( n > 0 )
    pos = strchrnth(str, separator, n);

  // if n <= 0 or separator appears fewer than n times, get the last occurrence
  if ( pos == NULL)
    pos = strrchr(str, separator);
  return const_cast<char*>(pos);
}


// ----------------------------------------------------------------------
// Misc. routines
// ----------------------------------------------------------------------

bool IsAscii(const char* str, int len) {
  const char* end = str + len;
  while (str < end) {
    if (!ascii_isascii(*str++)) {
      return false;
    }
  }
  return true;
}

namespace strings {

bool IsPrint(StringPiece str) {
  const char* strp = str.data();
  const char* end = strp + str.size();
  while (strp < end) {
    if (!ascii_isprint(*strp++)) {
      return false;
    }
  }
  return true;
}

}  // namespace strings

string StringReplace(StringPiece s, StringPiece oldsub,
                     StringPiece newsub, bool replace_all) {
  string ret;
  StringReplace(s, oldsub, newsub, replace_all, &ret);
  return ret;
}

void StringReplace(StringPiece s, StringPiece oldsub,
                   StringPiece newsub, bool replace_all,
                   string* res) {
  if (oldsub.empty()) {
    res->append(s.data(), s.length());  // If empty, append the given string.
    return;
  }

  StringPiece::size_type start_pos = 0;
  StringPiece::size_type pos;
  do {
    pos = s.find(oldsub, start_pos);
    if (pos == StringPiece::npos) {
      break;
    }
    res->append(s.data() + start_pos, pos - start_pos);
    res->append(newsub.data(), newsub.length());
    // Start searching again after the "old".
    start_pos = pos + oldsub.length();
  } while (replace_all);
  res->append(s.data() + start_pos, s.length() - start_pos);
}

int GlobalReplaceSubstring(StringPiece substring,
                           StringPiece replacement,
                           string* s) {
  CHECK(s != NULL);
  if (s->empty() || substring.empty())
    return 0;
  string tmp;
  int num_replacements = 0;
  int pos = 0;
  for (int match_pos = s->find(substring.data(), pos, substring.length());
       match_pos != string::npos;
       pos = match_pos + substring.length(),
           match_pos = s->find(substring.data(), pos, substring.length())) {
    ++num_replacements;
    // Append the original content before the match.
    tmp.append(*s, pos, match_pos - pos);
    // Append the replacement for the match.
    tmp.append(replacement.begin(), replacement.end());
  }
  // Append the content after the last match. If no replacements were made, the
  // original string is left untouched.
  if (num_replacements > 0) {
    tmp.append(*s, pos, s->length() - pos);
    s->swap(tmp);
  }
  return num_replacements;
}

char *gstrcasestr(const char* haystack, const char* needle) {
  char c, sc;
  size_t len;

  if ((c = *needle++) != 0) {
    c = ascii_tolower(c);
    len = strlen(needle);
    do {
      do {
        if ((sc = *haystack++) == 0)
          return NULL;
      } while (ascii_tolower(sc) != c);
    } while (strncasecmp(haystack, needle, len) != 0);
    haystack--;
  }
  // This is a const violation but strstr() also returns a char*.
  return const_cast<char*>(haystack);
}

const char *gstrncasestr(const char* haystack, const char* needle, size_t len) {
  char c, sc;

  if ((c = *needle++) != 0) {
    c = ascii_tolower(c);
    size_t needle_len = strlen(needle);
    do {
      do {
        if (len-- <= needle_len
            || 0 == (sc = *haystack++))
          return NULL;
      } while (ascii_tolower(sc) != c);
    } while (strncasecmp(haystack, needle, needle_len) != 0);
    haystack--;
  }
  return haystack;
}

char *gstrncasestr_split(const char* str,
                         const char* prefix, char non_alpha,
                         const char* suffix,
                         size_t n) {
  int prelen = prefix == NULL ? 0 : strlen(prefix);
  int suflen = suffix == NULL ? 0 : strlen(suffix);

  // adjust the string and its length to avoid unnessary searching.
  // an added benefit is to avoid unnecessary range checks in the if
  // statement in the inner loop.
  if (suflen + prelen >= n)  return NULL;
  str += prelen;
  n -= prelen;
  n -= suflen;

  const char* where = NULL;

  // for every occurance of non_alpha in the string ...
  while ((where = static_cast<const char*>(
            memchr(str, non_alpha, n))) != NULL) {
    // ... test whether it is followed by suffix and preceded by prefix
    if ((!suflen || strncasecmp(where + 1, suffix, suflen) == 0) &&
        (!prelen || strncasecmp(where - prelen, prefix, prelen) == 0)) {
      return const_cast<char*>(where - prelen);
    }
    // if not, advance the pointer, and adjust the length according
    n -= (where + 1) - str;
    str = where + 1;
  }

  return NULL;
}

char *strcasestr_alnum(const char *haystack, const char *needle) {
  const char *haystack_ptr;
  const char *needle_ptr;

  // Skip non-alnums at beginning
  while ( !ascii_isalnum(*needle) )
    if ( *needle++ == '\0' )
      return const_cast<char*>(haystack);
  needle_ptr = needle;

  // Skip non-alnums at beginning
  while ( !ascii_isalnum(*haystack) )
    if ( *haystack++ == '\0' )
      return NULL;
  haystack_ptr = haystack;

  while ( *needle_ptr != '\0' ) {
    // Non-alnums - advance
    while ( !ascii_isalnum(*needle_ptr) )
      if ( *needle_ptr++ == '\0' )
        return const_cast<char *>(haystack);

    while ( !ascii_isalnum(*haystack_ptr) )
      if ( *haystack_ptr++ == '\0' )
        return NULL;

    if ( ascii_tolower(*needle_ptr) == ascii_tolower(*haystack_ptr) ) {
      // Case-insensitive match - advance
      needle_ptr++;
      haystack_ptr++;
    } else {
      // No match - rollback to next start point in haystack
      haystack++;
      while ( !ascii_isalnum(*haystack) )
        if ( *haystack++ == '\0' )
          return NULL;
      haystack_ptr = haystack;
      needle_ptr = needle;
    }
  }
  return const_cast<char *>(haystack);
}

int CountSubstring(StringPiece text, StringPiece substring) {
  CHECK_GT(substring.length(), 0);

  int count = 0;
  StringPiece::size_type curr = 0;
  while (StringPiece::npos != (curr = text.find(substring, curr))) {
    ++count;
    ++curr;
  }
  return count;
}

const char* strstr_delimited(const char* haystack,
                             const char* needle,
                             char delim) {
  if (!needle || !haystack) return NULL;
  if (*needle == '\0') return haystack;

  int needle_len = strlen(needle);

  while (true) {
    // Skip any leading delimiters.
    while (*haystack == delim) ++haystack;

    // Walk down the haystack, matching every character in the needle.
    const char* this_match = haystack;
    int i = 0;
    for (; i < needle_len; i++) {
      if (*haystack != needle[i]) {
        // We ran out of haystack or found a non-matching character.
        break;
      }
      ++haystack;
    }

    // If we matched the whole needle, ensure that it's properly delimited.
    if (i == needle_len && (*haystack == '\0' || *haystack == delim)) {
      return this_match;
    }

    // No match. Consume non-delimiter characters until we run out of them.
    while (*haystack != delim) {
      if (*haystack == '\0') return NULL;
      ++haystack;
    }
  }
  LOG(FATAL) << "Unreachable statement";
  return NULL;
}

// ----------------------------------------------------------------------
// Older versions of libc have a buggy strsep.
// ----------------------------------------------------------------------

char* gstrsep(char** stringp, const char* delim) {
  char *s;
  const char *spanp;
  int c, sc;
  char *tok;

  if ((s = *stringp) == NULL)
    return NULL;

  tok = s;
  while (true) {
    c = *s++;
    spanp = delim;
    do {
      if ((sc = *spanp++) == c) {
        if (c == 0)
          s = NULL;
        else
          s[-1] = 0;
        *stringp = s;
        return tok;
      }
    } while (sc != 0);
  }

  return NULL; /* should not happen */
}

void FastStringAppend(string* s, const char* data, int len) {
  STLAppendToString(s, data, len);
}

// TODO(jyrki): add a microbenchmark and revisit
// the optimizations done here.
//
// Several converters use this table to reduce
// division and modulo operations.
extern const char two_ASCII_digits[100][2];

const char two_ASCII_digits[100][2] = {
  {'0', '0'}, {'0', '1'}, {'0', '2'}, {'0', '3'}, {'0', '4'},
  {'0', '5'}, {'0', '6'}, {'0', '7'}, {'0', '8'}, {'0', '9'},
  {'1', '0'}, {'1', '1'}, {'1', '2'}, {'1', '3'}, {'1', '4'},
  {'1', '5'}, {'1', '6'}, {'1', '7'}, {'1', '8'}, {'1', '9'},
  {'2', '0'}, {'2', '1'}, {'2', '2'}, {'2', '3'}, {'2', '4'},
  {'2', '5'}, {'2', '6'}, {'2', '7'}, {'2', '8'}, {'2', '9'},
  {'3', '0'}, {'3', '1'}, {'3', '2'}, {'3', '3'}, {'3', '4'},
  {'3', '5'}, {'3', '6'}, {'3', '7'}, {'3', '8'}, {'3', '9'},
  {'4', '0'}, {'4', '1'}, {'4', '2'}, {'4', '3'}, {'4', '4'},
  {'4', '5'}, {'4', '6'}, {'4', '7'}, {'4', '8'}, {'4', '9'},
  {'5', '0'}, {'5', '1'}, {'5', '2'}, {'5', '3'}, {'5', '4'},
  {'5', '5'}, {'5', '6'}, {'5', '7'}, {'5', '8'}, {'5', '9'},
  {'6', '0'}, {'6', '1'}, {'6', '2'}, {'6', '3'}, {'6', '4'},
  {'6', '5'}, {'6', '6'}, {'6', '7'}, {'6', '8'}, {'6', '9'},
  {'7', '0'}, {'7', '1'}, {'7', '2'}, {'7', '3'}, {'7', '4'},
  {'7', '5'}, {'7', '6'}, {'7', '7'}, {'7', '8'}, {'7', '9'},
  {'8', '0'}, {'8', '1'}, {'8', '2'}, {'8', '3'}, {'8', '4'},
  {'8', '5'}, {'8', '6'}, {'8', '7'}, {'8', '8'}, {'8', '9'},
  {'9', '0'}, {'9', '1'}, {'9', '2'}, {'9', '3'}, {'9', '4'},
  {'9', '5'}, {'9', '6'}, {'9', '7'}, {'9', '8'}, {'9', '9'}
};

static inline void PutTwoDigits(int i, char* p) {
  DCHECK_GE(i, 0);
  DCHECK_LT(i, 100);
  p[0] = two_ASCII_digits[i][0];
  p[1] = two_ASCII_digits[i][1];
}

char* FastTimeToBuffer(time_t s, char* buffer) {
  if (s == 0) {
    time(&s);
  }

  struct tm tm;
  if (PortableSafeGmtime(&s, &tm) == NULL) {
    // Error message must fit in 30-char buffer.
    memcpy(buffer, "Invalid:", sizeof("Invalid:"));
    FastInt64ToBufferLeft(s, buffer+strlen(buffer));
    return buffer;
  }

  // strftime format: "%a, %d %b %Y %H:%M:%S GMT",
  // but strftime does locale stuff which we do not want
  // plus strftime takes > 10x the time of hard code

  const char* weekday_name = "Xxx";
  switch (tm.tm_wday) {
    default: { DLOG(FATAL) << "tm.tm_wday: " << tm.tm_wday; } break;
    case 0:  weekday_name = "Sun"; break;
    case 1:  weekday_name = "Mon"; break;
    case 2:  weekday_name = "Tue"; break;
    case 3:  weekday_name = "Wed"; break;
    case 4:  weekday_name = "Thu"; break;
    case 5:  weekday_name = "Fri"; break;
    case 6:  weekday_name = "Sat"; break;
  }

  const char* month_name = "Xxx";
  switch (tm.tm_mon) {
    default:  { DLOG(FATAL) << "tm.tm_mon: " << tm.tm_mon; } break;
    case 0:   month_name = "Jan"; break;
    case 1:   month_name = "Feb"; break;
    case 2:   month_name = "Mar"; break;
    case 3:   month_name = "Apr"; break;
    case 4:   month_name = "May"; break;
    case 5:   month_name = "Jun"; break;
    case 6:   month_name = "Jul"; break;
    case 7:   month_name = "Aug"; break;
    case 8:   month_name = "Sep"; break;
    case 9:   month_name = "Oct"; break;
    case 10:  month_name = "Nov"; break;
    case 11:  month_name = "Dec"; break;
  }

  // Write out the buffer.

  memcpy(buffer+0, weekday_name, 3);
  buffer[3] = ',';
  buffer[4] = ' ';

  PutTwoDigits(tm.tm_mday, buffer+5);
  buffer[7] = ' ';

  memcpy(buffer+8, month_name, 3);
  buffer[11] = ' ';

  int32 year = tm.tm_year + 1900;
  if (year >= 0 && year <= 9999) {
    PutTwoDigits(year/100, buffer+12);
  } else {
    memcpy(buffer, "Invalid:", sizeof("Invalid:"));
    FastInt64ToBufferLeft(s, buffer+strlen(buffer));
    return buffer;
  }
  PutTwoDigits(year%100, buffer+14);
  buffer[16] = ' ';

  PutTwoDigits(tm.tm_hour, buffer+17);
  buffer[19] = ':';

  PutTwoDigits(tm.tm_min, buffer+20);
  buffer[22] = ':';

  PutTwoDigits(tm.tm_sec, buffer+23);

  // includes ending NUL
  memcpy(buffer+25, " GMT", 5);

  return buffer;
}

char* strdup_with_new(const char* the_string) {
  if (the_string == NULL)
    return NULL;
  else
    return strndup_with_new(the_string, strlen(the_string));
}

char* strndup_with_new(const char* the_string, int max_length) {
  if (the_string == NULL)
    return NULL;

  char* result = new char[max_length + 1];
  result[max_length] = '\0';  // terminate the string because strncpy might not
  return strncpy(result, the_string, max_length);
}

const char* ScanForFirstWord(const char* the_string, const char** end_ptr) {
  CHECK(end_ptr != NULL) << ": precondition violated";

  if (the_string == NULL)  // empty string
    return NULL;

  const char* curr = the_string;
  while ((*curr != '\0') && ascii_isspace(*curr))  // skip initial spaces
    ++curr;

  if (*curr == '\0')  // no valid word found
    return NULL;

  // else has a valid word
  const char* first_word = curr;

  // now locate the end of the word
  while ((*curr != '\0') && !ascii_isspace(*curr))
    ++curr;

  *end_ptr = curr;
  return first_word;
}

namespace strings {

StringPiece ScanForFirstWord(StringPiece input) {
  const char* curr = input.data();
  const char* const end = curr + input.size();

  // Skip initial spaces to locate the start of the word.
  while ((curr < end) && ascii_isspace(*curr))
    ++curr;
  const char* const word = curr;

  // Skip subsequent non-spaces to locate the end of the word.
  while ((curr < end) && !ascii_isspace(*curr))
    ++curr;
  return StringPiece(word, curr - word);
}

}  // namespace strings

const char *AdvanceIdentifier(const char *str) {
  // Not using isalpha and isalnum so as not to rely on the locale.
  // We could have used ascii_isalpha and ascii_isalnum.
  char ch = *str++;
  if (!((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') || ch == '_'))
    return NULL;
  while (true) {
    ch = *str;
    if (!((ch >= 'a' && ch <= 'z') || (ch >= 'A' && ch <= 'Z') ||
          (ch >= '0' && ch <= '9') || ch == '_'))
      return str;
    str++;
  }
}

bool IsIdentifier(const char *str) {
  const char *end = AdvanceIdentifier(str);
  return end && *end == '\0';
}

void UniformInsertString(string* s, int interval, const char* separator) {
  const size_t separator_len = strlen(separator);

  if (interval < 1 ||      // invalid interval
      s->empty() ||        // nothing to do
      separator_len == 0)  // invalid separator
    return;

  int num_inserts = (s->size() - 1) / interval;  // -1 to avoid appending at end
  if (num_inserts == 0)  // nothing to do
    return;

  string tmp;
  tmp.reserve(s->size() + num_inserts * separator_len + 1);

  for (int i = 0; i < num_inserts ; ++i) {
    // append this interval
    tmp.append(*s, i * interval, interval);
    // append a separator
    tmp.append(separator, separator_len);
  }

  // append the tail
  const size_t tail_pos = num_inserts * interval;
  tmp.append(*s, tail_pos, s->size() - tail_pos);

  s->swap(tmp);
}

void InsertString(string *const s,
                  const vector<uint32> &indices,
                  char const *const separator) {
  const unsigned num_indices(indices.size());
  if (num_indices == 0) {
    return;  // nothing to do...
  }

  const unsigned separator_len(strlen(separator));
  if (separator_len == 0) {
    return;  // still nothing to do...
  }

  string tmp;
  const unsigned s_len(s->size());
  tmp.reserve(s_len + separator_len * num_indices);

  vector<uint32>::const_iterator const ind_end(indices.end());
  vector<uint32>::const_iterator ind_pos(indices.begin());

  uint32 last_pos(0);
  while (ind_pos != ind_end) {
    const uint32 pos(*ind_pos);
    DCHECK_GE(pos, last_pos);
    DCHECK_LE(pos, s_len);

    tmp.append(s->substr(last_pos, pos - last_pos));
    tmp.append(separator);

    last_pos = pos;
    ++ind_pos;
  }
  tmp.append(s->substr(last_pos));

  s->swap(tmp);
}

int FindNth(StringPiece s, char c, int n) {
  int pos = -1;

  for ( int i = 0; i < n; ++i ) {
    pos = s.find_first_of(c, pos + 1);
    if ( pos == StringPiece::npos ) {
      break;
    }
  }
  return pos;
}

int ReverseFindNth(StringPiece s, char c, int n) {
  if ( n <= 0 ) {
    return static_cast<int>(StringPiece::npos);
  }

  int pos = s.size();

  for ( int i = 0; i < n; ++i ) {
    // If pos == 0, we return StringPiece::npos right away. Otherwise,
    // the following find_last_of call would take (pos - 1) as string::npos,
    // which means it would again search the entire input string.
    if (pos == 0) {
      return static_cast<int>(StringPiece::npos);
    }
    pos = s.find_last_of(c, pos - 1);
    if ( pos == string::npos ) {
      break;
    }
  }
  return pos;
}

namespace strings {

// FindEol()
// Returns the location of the next end-of-line sequence.

StringPiece FindEol(StringPiece s) {
  for (size_t i = 0; i < s.length(); ++i) {
    if (s[i] == '\n') {
      return StringPiece(s.data() + i, 1);
    }
    if (s[i] == '\r') {
      if (i+1 < s.length() && s[i+1] == '\n') {
        return StringPiece(s.data() + i, 2);
      } else {
        return StringPiece(s.data() + i, 1);
      }
    }
  }
  return StringPiece(s.data() + s.length(), 0);
}

}  // namespace strings

bool OnlyWhitespace(StringPiece s) {
  for ( int i = 0; i < s.size(); ++i ) {
    if ( !ascii_isspace(s[i]) ) return false;
  }
  return true;
}

string PrefixSuccessor(StringPiece prefix) {
  // We can increment the last character in the string and be done
  // unless that character is 255 (0xff), in which case we have to erase the
  // last character and increment the previous character, unless that
  // is 255, etc. If the string is empty or consists entirely of
  // 255's, we just return the empty string.
  bool done = false;
  string limit(prefix.data(), prefix.size());
  int index = limit.length() - 1;
  while (!done && index >= 0) {
    if (limit[index] == '\xff') {  // char literal avoids signed/unsigned.
      limit.erase(index);
      index--;
    } else {
      limit[index]++;
      done = true;
    }
  }
  if (!done) {
    return "";
  } else {
    return limit;
  }
}

string ImmediateSuccessor(StringPiece s) {
  // Return the input string, with an additional NUL byte appended.
  string out;
  out.reserve(s.size() + 1);
  out.append(s.data(), s.size());
  out.push_back('\0');
  return out;
}

void FindShortestSeparator(StringPiece start,
                           StringPiece limit,
                           string* separator) {
  // Find length of common prefix
  size_t min_length = std::min(start.size(), limit.size());
  size_t diff_index = 0;
  while ((diff_index < min_length) &&
         (start[diff_index] == limit[diff_index])) {
    diff_index++;
  }

  if (diff_index >= min_length) {
    // Handle the case where either string is a prefix of the other
    // string, or both strings are identical.
    start.CopyToString(separator);
    return;
  }

  if (diff_index+1 == start.size()) {
    // If the first difference is in the last character, do not bother
    // incrementing that character since the separator will be no
    // shorter than "start".
    start.CopyToString(separator);
    return;
  }

  if (start[diff_index] == '\xff') {  // char literal avoids signed/unsigned.
    // Avoid overflow when incrementing start[diff_index]
    start.CopyToString(separator);
    return;
  }

  separator->assign(start.data(), diff_index);
  separator->push_back(start[diff_index] + 1);
  if (*separator >= limit) {
    // Never pick a separator that causes confusion with "limit"
    start.CopyToString(separator);
  }
}

int SafeSnprintf(char *str, size_t size, const char *format, ...) {
  va_list printargs;
  va_start(printargs, format);
  int ncw = vsnprintf(str, size, format, printargs);
  va_end(printargs);
  return (ncw < size && ncw >= 0) ? ncw : 0;
}

bool GetlineFromStdioFile(FILE* file, string* str, char delim) {
  str->erase();
  while (true) {
    if (feof(file) || ferror(file)) {
      return false;
    }
    int c = getc(file);
    if (c == EOF) return false;
    if (c == delim) return true;
    str->push_back(c);
  }
}
