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
// File: numbers.h
// -----------------------------------------------------------------------------
//
// This package contains functions for converting strings to numbers. For
// converting numbers to strings, use `StrCat()` or `StrAppend()` in str_cat.h,
// which automatically detect and convert most number values appropriately.

#ifndef S2_THIRD_PARTY_ABSL_STRINGS_NUMBERS_H_
#define S2_THIRD_PARTY_ABSL_STRINGS_NUMBERS_H_

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <limits>
#include <string>
#include <type_traits>

#include "s2/third_party/absl/base/macros.h"
#include "s2/third_party/absl/base/port.h"
#include "s2/third_party/absl/numeric/int128.h"
#include "s2/third_party/absl/strings/ascii.h"
#include "s2/third_party/absl/strings/ascii_ctype.h"
#include "s2/third_party/absl/strings/string_view.h"

namespace absl {

// SimpleAtoi()
//
// Converts the given string into an integer value, returning `true` if
// successful. The string must reflect a base-10 integer (optionally followed or
// preceded by ASCII whitespace) whose value falls within the range of the
// integer type,
template <typename int_type>
ABSL_MUST_USE_RESULT bool SimpleAtoi(absl::string_view s, int_type* out);

// SimpleAtof()
//
// Converts the given string (optionally followed or preceded by ASCII
// whitespace) into a float, which may be rounded on overflow or underflow.
ABSL_MUST_USE_RESULT bool SimpleAtof(absl::string_view str, float* value);

// SimpleAtod()
//
// Converts the given string (optionally followed or preceded by ASCII
// whitespace) into a double, which may be rounded on overflow or underflow.
ABSL_MUST_USE_RESULT bool SimpleAtod(absl::string_view str, double* value);

// SimpleAtob()
//
// Converts the given string into a boolean, returning `true` if successful.
// The following case-insensitive strings are interpreted as boolean `true`:
// "true", "t", "yes", "y", "1". The following case-insensitive strings
// are interpreted as boolean `false`: "false", "f", "no", "n", "0".
ABSL_MUST_USE_RESULT bool SimpleAtob(absl::string_view str, bool* value);

}  // namespace absl

// End of public API.  Implementation details follow.

namespace absl {
namespace numbers_internal {

// safe_strto?() functions for implementing SimpleAtoi()
bool safe_strto32_base(absl::string_view text, int32_t* value, int base);
bool safe_strto64_base(absl::string_view text, int64_t* value, int base);
bool safe_strtou32_base(absl::string_view text, uint32_t* value, int base);
bool safe_strtou64_base(absl::string_view text, uint64_t* value, int base);

// These functions are intended for speed. All functions take an output buffer
// as an argument and return a pointer to the last byte they wrote, which is the
// terminating '\0'. At most `kFastToBufferSize` bytes are written.
char* FastInt32ToBuffer(int32_t i, char* buffer);
char* FastUInt32ToBuffer(uint32_t i, char* buffer);
char* FastInt64ToBuffer(int64_t i, char* buffer);
char* FastUInt64ToBuffer(uint64_t i, char* buffer);

static const int kFastToBufferSize = 32;
static const int kSixDigitsToBufferSize = 16;

char* RoundTripDoubleToBuffer(double d, char* buffer);
char* RoundTripFloatToBuffer(float f, char* buffer);

// Helper function for fast formatting of floating-point values.
// The result is the same as printf's "%g", a.k.a. "%.6g"; that is, six
// significant digits are returned, trailing zeros are removed, and numbers
// outside the range 0.0001-999999 are output using scientific notation
// (1.23456e+06). This routine is heavily optimized.
// Required buffer size is `kSixDigitsToBufferSize`.
size_t SixDigitsToBuffer(double d, char* buffer);

template <typename int_type>
char* FastIntToBuffer(int_type i, char* buffer) {
  static_assert(sizeof(i) <= 64 / 8,
                "FastIntToBuffer works only with 64-bit-or-less integers.");
  // TODO(user): This signed-ness check is used because it works correctly
  // with enums, and it also serves to check that int_type is not a pointer.
  // If one day something like std::is_signed<enum E> works, switch to it.
  if (static_cast<int_type>(1) - 2 < 0) {  // Signed
    if (sizeof(i) > 32 / 8) {           // 33-bit to 64-bit
      return numbers_internal::FastInt64ToBuffer(i, buffer);
    } else {  // 32-bit or less
      return numbers_internal::FastInt32ToBuffer(i, buffer);
    }
  } else {                     // Unsigned
    if (sizeof(i) > 32 / 8) {  // 33-bit to 64-bit
      return numbers_internal::FastUInt64ToBuffer(i, buffer);
    } else {  // 32-bit or less
      return numbers_internal::FastUInt32ToBuffer(i, buffer);
    }
  }
}

}  // namespace numbers_internal

// SimpleAtoi()
//
// Converts a string to an integer, using `safe_strto?()` functions for actual
// parsing, returning `true` if successful. The `safe_strto?()` functions apply
// strict checking; the string must be a base-10 integer, optionally followed or
// preceded by ASCII whitespace, with a value in the range of the corresponding
// integer type.
template <typename int_type>
ABSL_MUST_USE_RESULT bool SimpleAtoi(absl::string_view s, int_type* out) {
  static_assert(sizeof(*out) == 4 || sizeof(*out) == 8,
                "SimpleAtoi works only with 32-bit or 64-bit integers.");
  static_assert(!std::is_floating_point<int_type>::value,
                "Use SimpleAtof or SimpleAtod instead.");
  bool parsed;
  // TODO(user): This signed-ness check is used because it works correctly
  // with enums, and it also serves to check that int_type is not a pointer.
  // If one day something like std::is_signed<enum E> works, switch to it.
  if (static_cast<int_type>(1) - 2 < 0) {  // Signed
    if (sizeof(*out) == 64 / 8) {       // 64-bit
      int64_t val;
      parsed = numbers_internal::safe_strto64_base(s, &val, 10);
      *out = static_cast<int_type>(val);
    } else {  // 32-bit
      int32_t val;
      parsed = numbers_internal::safe_strto32_base(s, &val, 10);
      *out = static_cast<int_type>(val);
    }
  } else {                         // Unsigned
    if (sizeof(*out) == 64 / 8) {  // 64-bit
      uint64_t val;
      parsed = numbers_internal::safe_strtou64_base(s, &val, 10);
      *out = static_cast<int_type>(val);
    } else {  // 32-bit
      uint32_t val;
      parsed = numbers_internal::safe_strtou32_base(s, &val, 10);
      *out = static_cast<int_type>(val);
    }
  }
  return parsed;
}

}  // namespace absl

inline bool safe_strtof(absl::string_view str, float* value) {
  return absl::SimpleAtof(str, value);
}
inline bool safe_strtod(absl::string_view str, double* value) {
  return absl::SimpleAtod(str, value);
}
inline bool safe_strtob(absl::string_view str, bool* value) {
  return absl::SimpleAtob(str, value);
}

using absl::numbers_internal::  // NOLINT(readability/namespace)
    SixDigitsToBuffer;
using absl::numbers_internal::  // NOLINT(readability/namespace)
    kFastToBufferSize;
using absl::numbers_internal::  // NOLINT(readability/namespace)
    kSixDigitsToBufferSize;

inline char* FastInt32ToBufferLeft(int32_t i, char* buffer) {
  return absl::numbers_internal::FastInt32ToBuffer(i, buffer);
}
inline char* FastUInt32ToBufferLeft(uint32_t i, char* buffer) {
  return absl::numbers_internal::FastUInt32ToBuffer(i, buffer);
}
inline char* FastInt64ToBufferLeft(int64 i, char* buffer) {
  return absl::numbers_internal::FastInt64ToBuffer(i, buffer);
}
inline char* FastUInt64ToBufferLeft(uint64 i, char* buffer) {
  return absl::numbers_internal::FastUInt64ToBuffer(i, buffer);
}

inline char* FastInt32ToBuffer(int32_t i, char* buffer) {
  absl::numbers_internal::FastInt32ToBuffer(i, buffer);
  return buffer;
}
inline char* FastUInt32ToBuffer(uint32_t i, char* buffer) {
  absl::numbers_internal::FastUInt32ToBuffer(i, buffer);
  return buffer;
}
inline char* FastInt64ToBuffer(int64 i, char* buffer) {
  absl::numbers_internal::FastInt64ToBuffer(i, buffer);
  return buffer;
}
inline char* FastUInt64ToBuffer(uint64 i, char* buffer) {
  absl::numbers_internal::FastUInt64ToBuffer(i, buffer);
  return buffer;
}

template <typename int_type>
char* FastIntToBufferLeft(int_type i, char* buffer) {
  return absl::numbers_internal::FastIntToBuffer(i, buffer);
}

using absl::SimpleAtoi;  // NOLINT(readability/namespace)

// Converts a double or float into a string which, if passed to `strtod()` or
// `strtof()` respectively, will produce the exact same original double or
// float.
//
// Exception: for NaN values,` strtod(RoundTripDtoa(NaN))` or
// `strtof(RoundTripFtoa(NaN))` may produce any NaN value, not necessarily the
// exact same original NaN value.
//
// Note: Calls to `RoundTrip*toa()` should preferably be replaced with
// `absl::StrCat(absl::LegacyPrecision(d))`.
//
// This routine attempts to produce a short output string; however it is not
// guaranteed to be as short as possible.
string RoundTripDtoa(double value);
string RoundTripFtoa(float value);
inline string RoundTripFtoa(double value) {
  // Cast to float so we get single-precision output.
  return RoundTripFtoa(static_cast<float>(value));
}

// Overloads of RoundTrip(*toa() to prevent accidentally passing integers to the
// RoundTrip formatters. It is expensive to use these functions to convert
// integers to strings. Instead, please use the implicit formatting provided by
// `StrCat()` and `StrAppend()`.
string RoundTripDtoa(int value) = delete;
string RoundTripFtoa(int value) = delete;

// Converts an integer into a string, and is much faster than `printf("%d")`.
//
// Note: Calls to `SimpleItoa()` should preferably be replaced with
// `absl::StrCat(i)`.
template <typename int_type>
string SimpleItoa(int_type i);

ABSL_MUST_USE_RESULT inline string SimpleFtoa(float f) {
  return RoundTripFtoa(f);
}
ABSL_MUST_USE_RESULT inline string SimpleDtoa(double d) {
  return RoundTripDtoa(d);
}

// Converts a boolean into a string, which if passed to `safe_strtob()` will
// produce the exact same original boolean, returning `true` if the value ==
// true, and `false` otherwise.
string SimpleBtoa(bool value);

ABSL_DEPRECATED("Use absl::StrCat to convert numbers to strings")
char* DoubleToBuffer(double i, char* buffer);
ABSL_DEPRECATED("Use absl::StrCat to convert numbers to strings")
char* FloatToBuffer(float i, char* buffer);

using absl::numbers_internal::  // NOLINT(readability/namespace)
    safe_strto32_base;
using absl::numbers_internal::  // NOLINT(readability/namespace)
    safe_strtou32_base;

bool safe_strtosize_t_base(absl::string_view text, size_t* value, int base);

inline bool safe_strtou64_base(absl::string_view text, uint64* value,
                               int base) {
  uint64_t val;
  bool result = absl::numbers_internal::safe_strtou64_base(text, &val, base);
  *value = val;
  return result;
}
inline bool safe_strto64_base(absl::string_view text, int64* value, int base) {
  int64_t val;
  bool result = absl::numbers_internal::safe_strto64_base(text, &val, base);
  *value = val;
  return result;
}

// Convenience functions for a base value == 10.
inline bool safe_strto32(absl::string_view text, int32_t* value) {
  return safe_strto32_base(text, value, 10);
}

inline bool safe_strto64(absl::string_view text, int64* value) {
  return safe_strto64_base(text, value, 10);
}

inline bool safe_strtou32(absl::string_view text, uint32_t* value) {
  return safe_strtou32_base(text, value, 10);
}

inline bool safe_strtou64(absl::string_view text, uint64* value) {
  return safe_strtou64_base(text, value, 10);
}

inline bool safe_strtosize_t(absl::string_view text, size_t* value) {
  return safe_strtosize_t_base(text, value, 10);
}

// Converts a fingerprint to 16 hex digits.
string FpToString(Fprint fp);

// Converts a string of 16 hex digits to a fingerprint, returning `true` on
// successful conversion, or `false` on invalid input.
bool StringToFp(absl::string_view hex, Fprint* fp);

// Converts between uint128 and 32-digit hex string. Note: no leading 0x or
// 0X prefix is generated.
string Uint128ToHexString(uint128 ui128);

// Converts a hex string to an uint128 value, returning `true` on successful
// conversion, or `false` on invalid input. If the string contains a character
// which is not a valid hexadecimal digit, this is deemed an invalid input,
// including "x", so a leading 0x or 0X prefix is not allowed.
bool HexStringToUint128(absl::string_view hex, uint128* value);

// Converts the number argument to a string representation in base-36, returning
// the number of bytes written, not including terminating NUL. Conversion fails
// if the buffer is too small to to hold the string and terminating NUL. A
// return value of 0 indicates an error.
size_t u64tostr_base36(uint64_t number, size_t buf_size, char* buffer);

// Converts the given string representation into a 64-bit unsigned integer
// value similar to `atoi(s)`, except `s` may refer to metric size suffixes for
// kilo, mega, giga, and tera. (E.g. 16k", "32M", "2G", "4t").
uint64 atoi_kmgt(const char* s);
inline uint64 atoi_kmgt(const string& s) { return atoi_kmgt(s.c_str()); }

// These functions are intended for speed. `FastTimeToBuffer()` puts the output
// into RFC822 format.
//
// `FastIntToBuffer()` uses at most 22 bytes; `FastTimeToBuffer()` uses exactly
// 30 bytes. They all return a pointer to the beginning of the output, which is
// the same as the beginning of the input buffer.
//
// NOTE: In 64-bit land, `sizeof(time_t)` is 8, so it is possible to pass to
// `FastTimeToBuffer()` a time whose year cannot be represented in 4 digits. In
// this case, the output buffer will contain the string "Invalid:<value>"
//
// WARNING: `FastTimeToBuffer(0, ...)` returns the current time, not 1970.
// WARNING: This "0" behavior is deprecated.  Please pass `time(nullptr)`
//          if you want a string from the current time.
//
// Previously documented minimums -- the buffers provided must be at least this
// long, though these numbers are subject to change:
//
//     Int32, UInt32:                   12 bytes
//     Int64, UInt64, Int, Uint:        22 bytes
//     Time:                            30 bytes
//
// Use `kFastToBufferSize` rather than hardcoding constants.
char* FastUInt32ToBuffer(uint32_t i, char* buffer);
char* FastUInt64ToBuffer(uint64 i, char* buffer);

ABSL_DEPRECATED("Use FastFormatRFC1123GMT() from util/time/time.h")
char* FastTimeToBuffer(time_t t, char* buffer);

// Returns 1 if `buf` is prefixed by `num_digits` of hex digits; returns 0
// otherwise. The function checks for '\0' for string termination.
int HexDigitsPrefix(const char* buf, ptrdiff_t num_digits);

// Eliminates all leading zeroes (unless the string itself is composed
// of nothing but zeroes, in which case one is kept: 0...0 becomes 0).
void ConsumeStrayLeadingZeroes(string* str);

// A simple parser for int32 values. Returns the parsed value
// if a valid integer is found; else returns deflt. It does not
// check if str is entirely consumed.
// This cannot handle decimal numbers with leading 0s, since they will be
// treated as octal.  If you know it's decimal, use ParseLeadingDec32Value.
int32_t ParseLeadingInt32Value(absl::string_view str, int32_t deflt);

// A simple parser for uint32 values. Returns the parsed value
// if a valid integer is found; else returns deflt. It does not
// check if str is entirely consumed.
// This cannot handle decimal numbers with leading 0s, since they will be
// treated as octal.  If you know it's decimal, use ParseLeadingUDec32Value.
uint32_t ParseLeadingUInt32Value(absl::string_view str, uint32_t deflt);

// A simple parser for decimal int32 values. Returns the parsed value
// if a valid integer is found; else returns deflt. It does not
// check if str is entirely consumed.
// The string passed in is treated as *10 based*.
// This can handle strings with leading 0s.
// See also: ParseLeadingDec64Value
int32_t ParseLeadingDec32Value(absl::string_view str, int32_t deflt);

// A simple parser for decimal uint32 values. Returns the parsed value
// if a valid integer is found; else returns deflt. It does not
// check if str is entirely consumed.
// The string passed in is treated as *10 based*.
// This can handle strings with leading 0s.
// See also: ParseLeadingUDec64Value
uint32_t ParseLeadingUDec32Value(absl::string_view str, uint32_t deflt);

// A simple parser for long long values.
// Returns the parsed value if a
// valid integer is found; else returns deflt
uint64 ParseLeadingUInt64Value(absl::string_view str, uint64 deflt);
int64 ParseLeadingInt64Value(absl::string_view str, int64 deflt);
uint64 ParseLeadingHex64Value(absl::string_view str, uint64 deflt);
int64 ParseLeadingDec64Value(absl::string_view str, int64 deflt);
uint64 ParseLeadingUDec64Value(absl::string_view str, uint64 deflt);

// A simple parser for double values. Returns the parsed value
// if a valid double is found; else returns deflt. It does not
// check if str is entirely consumed.
double ParseLeadingDoubleValue(const char* str, double deflt);
inline double ParseLeadingDoubleValue(const string& str, double deflt) {
  return ParseLeadingDoubleValue(str.c_str(), deflt);
}

// A recognizer of boolean string values. Returns the parsed value
// if a valid value is found; else returns deflt.  This skips leading
// whitespace, is case insensitive, and recognizes these forms:
// 0/1, false/true, no/yes, n/y
bool ParseLeadingBoolValue(absl::string_view str, bool deflt);

// -----------------------------------------------------------------------------
// Natural Sort Order Utilities
// -----------------------------------------------------------------------------
//
// A Natural Sort Order sorts strings containing multi-digit characters ordered
// as if those digits were considered as one character, in numerical order.
// For example, "image9" is considered before "image10".  (This goes by the
// name `natsort()` in Go and PHP.)
//
// For more information, see https://en.wikipedia.org/wiki/Natural_sort_order.

// These are like std::less<string> and std::greater<string>, except when a
// run of digits is encountered at corresponding points in the two
// arguments.  Such digit strings are compared numerically instead
// of lexicographically.  Therefore if you sort by
// "autodigit_less", some machine names might get sorted as:
//    exaf1
//    exaf2
//    exaf10
// When using "strict" comparison (AutoDigitStrCmp with the strict flag
// set to true, or the strict version of the other functions),
// strings that represent equal numbers will not be considered equal if
// the string representations are not identical.  That is, "01" < "1" in
// strict mode, but "01" == "1" otherwise.

int AutoDigitStrCmp(absl::string_view a, absl::string_view b, bool strict);
inline bool AutoDigitLessThan(absl::string_view a, absl::string_view b) {
  return AutoDigitStrCmp(a, b, false) < 0;
}
inline bool StrictAutoDigitLessThan(absl::string_view a, absl::string_view b) {
  return AutoDigitStrCmp(a, b, true) < 0;
}

// For speed, when you know you have zero-terminated strings.
int AutoDigitStrCmpZ(const char* a, const char* b, bool strict);

struct autodigit_less {
  bool operator()(const string& a, const string& b) const {
    return AutoDigitStrCmp(a, b, /*strict=*/false) < 0;
  }
};

struct autodigit_greater {
  bool operator()(const string& a, const string& b) const {
    return AutoDigitStrCmp(a, b, /*strict=*/false) > 0;
  }
};

struct strict_autodigit_less {
  bool operator()(const string& a, const string& b) const {
    return AutoDigitStrCmp(a, b, /*strict=*/true) < 0;
  }
};

struct strict_autodigit_greater {
  bool operator()(const string& a, const string& b) const {
    return AutoDigitStrCmp(a, b, /*strict=*/true) > 0;
  }
};

// Converts an integer to a string. Puts commas every 3 spaces. Faster than
// printf("%d")?
template <typename IntType>
inline string SimpleItoaWithCommas(IntType ii) {
  string s1 = SimpleItoa(ii);
  absl::string_view sp1(s1);
  string output;
  // Copy leading non-digit characters unconditionally.
  // This picks up the leading sign.
  while (!sp1.empty() && !absl::ascii_isdigit(sp1[0])) {
    output.push_back(sp1[0]);
    sp1.remove_prefix(1);
  }
  // Copy rest of input characters.
  for (absl::string_view::size_type i = 0; i < sp1.size(); ++i) {
    if (i > 0 && (sp1.size() - i) % 3 == 0) {
      output.push_back(',');
    }
    output.push_back(sp1[i]);
  }
  return output;
}

// Converts an integer to a string. Truncates values to K, G, M or T as
// appropriate. Opposite of atoi_kmgt() E.g. 3000 -> 2K   57185920 -> 45M
string ItoaKMGT(int64 i);

// Parses an expression in 'text' of the form: <double><sep><double> where
// <double> may be a double-precision number and <sep> is a single char or "..",
// and must be one of the chars in parameter 'separators', which may contain '-'
// or '.' (which means "..") or any chars not allowed in a double. If
// allow_unbounded_markers, <double> may also be a '?' to indicate unboundedness
// (if on the left of <sep>, means unbounded below; if on the right, means
// unbounded above). Depending on num_required_bounds, which may be 0, 1, or 2,
// <double> may also be the empty string, indicating unboundedness. If
// require_separator is false, then a single <double> is acceptable and is
// parsed as a range bounded from below. We also check that the character
// following the range must be in acceptable_terminators. If null_terminator_ok,
// then it is also OK if the range ends in \0 or after len chars. If
// allow_currency is true, the first <double> may be optionally preceded by a
// '$', in which case *is_currency will be true, and the second <double> may
// similarly be preceded by a '$'. In these cases, the '$' will be ignored
// (otherwise it's an error). If allow_comparators is true, the expression in
// 'text' may also be of the form <comparator><double>, where <comparator> is
// '<' or '>' or '<=' or '>='. separators and require_separator are ignored in
// this format, but all other parameters function as for the first format.
// Returns true if the expression is parsed successfully; false otherwise. If
// successful, output params are: 'end', which points to the char just beyond
// the expression; 'from' and 'to' are set to the values of the <double>s, and
// are -inf and inf (or unchanged, depending on dont_modify_unbounded) if
// unbounded. Output params are undefined if false is returned. len is the input
// length, or -1 if text is '\0'-terminated, which is more efficient.
struct DoubleRangeOptions {
  const char* separators;
  bool require_separator;
  const char* acceptable_terminators;
  bool null_terminator_ok;
  bool allow_unbounded_markers;
  uint32_t num_required_bounds;
  bool dont_modify_unbounded;
  bool allow_currency;
  bool allow_comparators;
};

bool ParseDoubleRange(const char* text, ptrdiff_t len, const char** end,
                      double* from, double* to, bool* is_currency,
                      const DoubleRangeOptions& opts);

template <typename int_type>
string SimpleItoa(int_type i) {
  char buf[absl::numbers_internal::kFastToBufferSize];
  return string(buf, absl::numbers_internal::FastIntToBuffer(i, buf));
}


#endif  // S2_THIRD_PARTY_ABSL_STRINGS_NUMBERS_H_
