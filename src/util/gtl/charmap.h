// Copyright 2015 Google Inc. All Rights Reserved.
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

// Character Map Class
//
// A fast, bit-vector map for 8-bit unsigned characters.
// This class is useful for non-character purposes as well.
//

#ifndef UTIL_GTL_CHARMAP_H_
#define UTIL_GTL_CHARMAP_H_

#include <string.h>

#include "base/integral_types.h"
#include "base/macros.h"
#include "base/port.h"

#ifdef LANG_CXX11
#define UTIL_GTL_CHARMAP_CONSTEXPR_ constexpr
#else
#define UTIL_GTL_CHARMAP_CONSTEXPR_
#endif  // !LANG_CXX11

class Charmap {
 public:
  UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap() : m_() {}

  // Initializes with a given char*.  Note that NUL is not treated as
  // a terminator, but rather a char to be flicked.
  Charmap(const char* str, int len) : m_() {
    while (len--) SetChar(*str++);
  }

  // Initializes with a given char*.  NUL is treated as a terminator
  // and will not be in the charmap.
  explicit Charmap(const char* str) : m_() {
    while (*str) SetChar(*str++);
  }

  UTIL_GTL_CHARMAP_CONSTEXPR_
  bool contains(unsigned char c) const {
    return (m_[c / 64] >> (c % 64)) & 0x1;
  }

  // Returns true if and only if a character exists in both maps.
  bool IntersectsWith(const Charmap& c) const {
    for (int i = 0; i < arraysize(m_); ++i) {
      if ((m_[i] & c.m_[i]) != 0)
        return true;
    }
    return false;
  }

  bool IsZero() const {
    for (int i = 0; i < arraysize(m_); ++i) {
      if (m_[i] != 0)
        return false;
    }
    return true;
  }

  // Containing only a single specified char.
  static UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap
  Char(char x) {
    return Charmap(CharMaskForWord(x, 0), CharMaskForWord(x, 1),
                   CharMaskForWord(x, 2), CharMaskForWord(x, 3));
  }

  // Containing all the chars in the C-string 's'.
  // CAUTION: deeply recursive. Use only in constexpr initializers.
  static UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap
  FromString(const char* s) {
    return *s == 0 ? Charmap() : (Char(*s) | FromString(s + 1));
  }

  // Containing all the chars in the closed interval [lo,hi].
  static UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap
  Range(char lo, char hi) {
    return Charmap(RangeForWord(lo, hi, 0), RangeForWord(lo, hi, 1),
                   RangeForWord(lo, hi, 2), RangeForWord(lo, hi, 3));
  }

  friend UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap
  operator&(const Charmap& a, const Charmap& b) {
    return Charmap(a.m_[0] & b.m_[0], a.m_[1] & b.m_[1],
                   a.m_[2] & b.m_[2], a.m_[3] & b.m_[3]);
  }

  friend UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap
  operator|(const Charmap& a, const Charmap& b) {
    return Charmap(a.m_[0] | b.m_[0], a.m_[1] | b.m_[1],
                   a.m_[2] | b.m_[2], a.m_[3] | b.m_[3]);
  }

  friend UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap
  operator~(const Charmap& a) {
    return Charmap(~a.m_[0], ~a.m_[1], ~a.m_[2], ~a.m_[3]);
  }

 private:
  UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap(
      uint64 b0, uint64 b1, uint64 b2, uint64 b3)
#ifdef LANG_CXX11
      : m_{b0, b1, b2, b3} {}
#else  // !LANG_CXX11  (C++98 version)
  {
    m_[0] = b0;
    m_[1] = b1;
    m_[2] = b2;
    m_[3] = b3;
  }
#endif  // !LANG_CXX11

  static UTIL_GTL_CHARMAP_CONSTEXPR_ uint64
  RangeForWord(unsigned char lo, unsigned char hi, uint64 word) {
    return OpenRangeFromZeroForWord(hi + 1, word)
           & ~OpenRangeFromZeroForWord(lo, word);
  }

  // All the chars in the specified word of the range [0, upper).
  static UTIL_GTL_CHARMAP_CONSTEXPR_ uint64
  OpenRangeFromZeroForWord(uint64 upper, uint64 word) {
    return (upper <= 64 * word) ? 0 :
           (upper >= 64 * (word + 1)) ? ~static_cast<uint64>(0) :
           (~static_cast<uint64>(0) >> (64 - upper % 64));
  }

  static UTIL_GTL_CHARMAP_CONSTEXPR_ uint64
  CharMaskForWord(unsigned char x, uint64 word) {
    return (x / 64 == word) ? (static_cast<uint64>(1) << (x % 64)) : 0;
  }

 private:
  void SetChar(unsigned char c) {
    m_[c / 64] |= static_cast<uint64>(1) << (c % 64);
  }

  uint64 m_[4];
};

namespace util {
namespace gtl {
namespace charmaps {
// Mirror the char-classifying predicates in <cctype>
inline UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap Upper() {
  return Charmap::Range('A', 'Z');
}
inline UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap Lower() {
  return Charmap::Range('a', 'z');
}
inline UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap Digit() {
  return Charmap::Range('0', '9');
}
inline UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap Alpha() {
  return Lower() | Upper();
}
inline UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap Alnum() {
  return Digit() | Alpha();
}
inline UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap XDigit() {
  return Digit() | Charmap::Range('A', 'F') | Charmap::Range('a', 'f');
}
inline UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap Print() {
  return Charmap::Range(0x20, 0x7e);
}
inline UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap Space() {
  return Charmap::FromString("\t\n\v\f\r ");
}
inline UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap Cntrl() {
  return Charmap::Range(0, 0x7f) & ~Print();
}
inline UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap Blank() {
  return Charmap::FromString("\t ");
}
inline UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap Graph() {
  return Print() & ~Space();
}
inline UTIL_GTL_CHARMAP_CONSTEXPR_ Charmap Punct() {
  return Graph() & ~Alnum();
}
}  // namespace charmaps
}  // namespace gtl
}  // namespace util
#undef UTIL_GTL_CHARMAP_CONSTEXPR_
#endif  // UTIL_GTL_CHARMAP_H_
