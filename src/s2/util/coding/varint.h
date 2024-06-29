// Copyright 2001 Google Inc. All Rights Reserved.
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
// Raw support for varint encoding.  Higher level interfaces are
// provided by Encoder/Decoder/IOBuffer.  Clients should typically use
// those interfaces, unless speed is paramount.
//
// If decoding speed is very important, consider using PrefixVarint instead.
// It has the same compression ratio, but generally faster decoding.
//
// Provided routines:
//      vi_parse32_unchecked
//      vi_parse64_unchecked
//      vi_encode32_unchecked
//      vi_encode64_unchecked

#ifndef S2_UTIL_CODING_VARINT_H_
#define S2_UTIL_CODING_VARINT_H_

// Avoid adding expensive includes here.
#include <cassert>
#include <cstddef>

#include <cstdint>
#include <limits>
#include <string>
#include <utility>

#include "absl/base/macros.h"
#include "absl/numeric/bits.h"
#include "s2/util/bits/bits.h"

// Just a namespace, not a real class
class Varint {
 public:
  // Maximum lengths of varint encoding of uint32 and uint64
  static constexpr int kMax32 = 5;
  static constexpr int kMax64 = 10;

  // The decoder does not read past the end of the encoded data.
  static constexpr int kSlopBytes = 0;

  // REQUIRES   "ptr" points to a buffer of length at least kMaxXX
  // EFFECTS    Scan next varint from "ptr" and store in OUTPUT.
  //            Returns pointer just past last read byte.  Returns
  //            nullptr if a valid varint value was not found.
  static const char* Parse32(const char* ptr, uint32_t* OUTPUT);
  static const char* Parse64(const char* ptr, uint64_t* OUTPUT);

  // A fully inlined version of Parse32: useful in the most time critical
  // routines, but its code size is large
  static const char* Parse32Inline(const char* ptr, uint32_t* OUTPUT);

  // REQUIRES   "ptr" points just past the last byte of a varint-encoded value.
  // REQUIRES   A second varint must be encoded just before the one we parse,
  //            OR "base" must point to the first byte of the one we parse.
  // REQUIRES   Bytes [base, ptr-1] are readable
  //
  // EFFECTS    Scan backwards from "ptr" and store in OUTPUT. Stop at the last
  //            byte of the previous varint, OR at "base", whichever one comes
  //            first. Returns pointer to the first byte of the decoded varint
  //            nullptr if a valid varint value was not found.
  static const char* Parse32Backward(const char* ptr, const char* base,
                                     uint32_t* OUTPUT);
  static const char* Parse64Backward(const char* ptr, const char* base,
                                     uint64_t* OUTPUT);

  // Attempts to parse a varint32 from a prefix of the bytes in [ptr,limit-1].
  // Never reads a character at or beyond limit.  If a valid/terminated varint32
  // was found in the range, stores it in *OUTPUT and returns a pointer just
  // past the last byte of the varint32. Else returns nullptr.  On success,
  // "result <= limit".
  static const char* Parse32WithLimit(const char* ptr, const char* limit,
                                      uint32_t* OUTPUT);
  static const char* Parse64WithLimit(const char* ptr, const char* limit,
                                      uint64_t* OUTPUT);

  // REQUIRES   "ptr" points to the first byte of a varint-encoded value.
  // EFFECTS     Scans until the end of the varint and returns a pointer just
  //             past the last byte. Returns nullptr if "ptr" does not point to
  //             a valid varint value.
  static const char* Skip32(const char* ptr);
  static const char* Skip64(const char* ptr);

  // REQUIRES   "ptr" points just past the last byte of a varint-encoded value.
  // REQUIRES   A second varint must be encoded just before the one we parse,
  //            OR "base" must point to the first byte of the one we parse.
  // REQUIRES   Bytes [base, ptr-1] are readable
  //
  // EFFECTS    Scan backwards from "ptr" and stop at the last byte of the
  //            previous varint, OR at "base", whichever one comes first.
  //            Returns pointer to the first byte of the skipped varint or
  //            nullptr if a valid varint value was not found.
  static const char* Skip32Backward(const char* ptr, const char* base);
  static const char* Skip64Backward(const char* ptr, const char* base);

  // REQUIRES   "ptr" points to a buffer of length sufficient to hold "v".
  // EFFECTS    Encodes "v" into "ptr" and returns a pointer to the
  //            byte just past the last encoded byte.
  static char* Encode32(char* ptr, uint32_t v);
  static char* Encode64(char* ptr, uint64_t v);

  // A fully inlined version of Encode32: useful in the most time critical
  // routines, but its code size is large
  static char* Encode32Inline(char* ptr, uint32_t v);

  // EFFECTS    Returns the encoding length of the specified value.
  static int Length32(uint32_t v);
  static int Length64(uint64_t v);

  // EFFECTS    Appends the varint representation of "value" to "*s".
  static void Append32(std::string* s, uint32_t value);
  static void Append64(std::string* s, uint64_t value);

  // EFFECTS    Encodes a pair of values to "*s".  The encoding
  //            is done by weaving together 4 bit groups of
  //            each number into a single 64 bit value, and then
  //            encoding this value as a Varint64 value.  This means
  //            that if both a and b are small, both values can be
  //            encoded in a single byte.
  ABSL_DEPRECATED("Use TwoValuesVarint::Encode32.")
  static void EncodeTwo32Values(std::string* s, uint32_t a, uint32_t b);

  // Decode and sum up a sequence of deltas until the sum >= goal.
  // It is significantly faster than calling ParseXXInline in a loop.
  // NOTE(user): The code does NO error checking, it assumes all the
  // deltas are valid and the sum of deltas will never exceed
  // numeric_limits<int64>::max(). The code works for both 32bits and
  // 64bits varint, and on 64 bits machines, the 64 bits version is
  // almost always faster. Thus we only have a 64 bits interface here.
  // The interface is slightly different from the other functions in that
  // it requires *signed* integers.
  // REQUIRES   "ptr" points to the first byte of a varint-encoded delta.
  //            The sum of deltas >= goal (the code does NO boundary check).
  //            goal is positive and fit into a signed int64.
  // EFFECTS    Returns a pointer just past last read byte.
  //            "out" stores the actual sum.
  static const char* FastDecodeDeltas(const char* ptr, int64_t goal,
                                      int64_t* out);

 private:
  static const char* Parse32FallbackInline(const char* p, uint32_t* val);
  static const char* Parse32Fallback(const char* p, uint32_t* val);
#if defined(__x86_64__)
  static std::pair<const char*, uint64_t> Parse64FallbackPair(const char* p,
                                                              int64_t res1);
#endif
  template <typename T>
  static int LengthImpl(T);
  static const char* Parse64Fallback(const char* p, uint64_t* val);
  static char* Encode32Fallback(char* ptr, uint32_t v);

  static const char* Parse32BackwardSlow(const char* ptr, const char* base,
                                         uint32_t* OUTPUT);
  static const char* Parse64BackwardSlow(const char* ptr, const char* base,
                                         uint64_t* OUTPUT);
  static const char* Skip32BackwardSlow(const char* ptr, const char* base);
  static const char* Skip64BackwardSlow(const char* ptr, const char* base);

  static void Append32Slow(std::string* s, uint32_t value);
  static void Append64Slow(std::string* s, uint64_t value);
};

/***** Implementation details; clients should ignore *****/

inline const char* Varint::Parse32FallbackInline(const char* p,
                                                 uint32_t* OUTPUT) {
  // Fast path
  const unsigned char* ptr = reinterpret_cast<const unsigned char*>(p);
  uint32_t byte, result;
  byte = *(ptr++); result = byte & 127;
  assert(byte >= 128);   // Already checked in inlined prelude
  byte = *(ptr++); result |= (byte & 127) <<  7; if (byte < 128) goto done;
  byte = *(ptr++); result |= (byte & 127) << 14; if (byte < 128) goto done;
  byte = *(ptr++); result |= (byte & 127) << 21; if (byte < 128) goto done;
  byte = *(ptr++); result |= (byte & 127) << 28; if (byte < 16) goto done;
  return nullptr;       // Value is too long to be a varint32
 done:
  *OUTPUT = result;
  return reinterpret_cast<const char*>(ptr);
 }

 inline const char* Varint::Parse32(const char* p, uint32_t* OUTPUT) {
   // Fast path for inlining
   const unsigned char* ptr = reinterpret_cast<const unsigned char*>(p);
   uint32_t byte = *ptr;
   if (byte < 128) {
     *OUTPUT = byte;
     return reinterpret_cast<const char*>(ptr) + 1;
   } else {
     return Parse32Fallback(p, OUTPUT);
   }
 }

 inline const char* Varint::Parse32Inline(const char* p, uint32_t* OUTPUT) {
   // Fast path for inlining
   const unsigned char* ptr = reinterpret_cast<const unsigned char*>(p);
   uint32_t byte = *ptr;
   if (byte < 128) {
     *OUTPUT = byte;
     return reinterpret_cast<const char*>(ptr) + 1;
   } else {
     return Parse32FallbackInline(p, OUTPUT);
   }
 }

inline const char* Varint::Skip32(const char* p) {
  const unsigned char* ptr = reinterpret_cast<const unsigned char*>(p);
  if (*ptr++ < 128) return reinterpret_cast<const char*>(ptr);
  if (*ptr++ < 128) return reinterpret_cast<const char*>(ptr);
  if (*ptr++ < 128) return reinterpret_cast<const char*>(ptr);
  if (*ptr++ < 128) return reinterpret_cast<const char*>(ptr);
  if (*ptr++ < 16) return reinterpret_cast<const char*>(ptr);
  return nullptr; // value is too long to be a varint32
}

inline const char* Varint::Parse32Backward(const char* p, const char* base,
                                           uint32_t* OUTPUT) {
  if (p > base + kMax32) {
    // Fast path
    const unsigned char* ptr = reinterpret_cast<const unsigned char*>(p);
    uint32_t byte, result;
    byte = *(--ptr); if (byte > 127) return nullptr;
    result = byte;
    byte = *(--ptr); if (byte < 128) goto done;
    result <<= 7; result |= (byte & 127);
    byte = *(--ptr); if (byte < 128) goto done;
    result <<= 7; result |= (byte & 127);
    byte = *(--ptr); if (byte < 128) goto done;
    result <<= 7; result |= (byte & 127);
    byte = *(--ptr); if (byte < 128) goto done;
    result <<= 7; result |= (byte & 127);
    byte = *(--ptr); if (byte < 128) goto done;
    return nullptr; // Value is too long to be a varint32
 done:
    *OUTPUT = result;
    return reinterpret_cast<const char*>(ptr+1);
  } else {
    return Parse32BackwardSlow(p, base, OUTPUT);
  }
}

inline const char* Varint::Skip32Backward(const char* p, const char* base) {
  if (p > base + kMax32) {
    const unsigned char* ptr = reinterpret_cast<const unsigned char*>(p);
    if (*(--ptr) > 127) return nullptr;
    if (*(--ptr) < 128) return reinterpret_cast<const char*>(ptr+1);
    if (*(--ptr) < 128) return reinterpret_cast<const char*>(ptr+1);
    if (*(--ptr) < 128) return reinterpret_cast<const char*>(ptr+1);
    if (*(--ptr) < 128) return reinterpret_cast<const char*>(ptr+1);
    if (*(--ptr) < 128) return reinterpret_cast<const char*>(ptr+1);
    return nullptr; // value is too long to be a varint32
  } else {
    return Skip32BackwardSlow(p, base);
  }
}

inline const char* Varint::Parse32WithLimit(const char* p, const char* l,
                                            uint32_t* OUTPUT) {
  // Version with bounds checks.
  // This formerly had an optimization to inline the non-bounds checking Parse32
  // but it was found to be slower than the straightforward implementation.
  const unsigned char* ptr = reinterpret_cast<const unsigned char*>(p);
  const unsigned char* limit = reinterpret_cast<const unsigned char*>(l);
  uint32_t b, result;
  if (ptr >= limit) return nullptr;
  b = *(ptr++); result = b & 127;          if (b < 128) goto done;
  if (ptr >= limit) return nullptr;
  b = *(ptr++); result |= (b & 127) <<  7; if (b < 128) goto done;
  if (ptr >= limit) return nullptr;
  b = *(ptr++); result |= (b & 127) << 14; if (b < 128) goto done;
  if (ptr >= limit) return nullptr;
  b = *(ptr++); result |= (b & 127) << 21; if (b < 128) goto done;
  if (ptr >= limit) return nullptr;
  b = *(ptr++); result |= (b & 127) << 28; if (b < 16) goto done;
  return nullptr;       // Value is too long to be a varint32
 done:
  *OUTPUT = result;
  return reinterpret_cast<const char*>(ptr);
 }

 inline const char* Varint::Parse64(const char* p, uint64_t* OUTPUT) {
#if defined(__x86_64__)
  auto ptr = reinterpret_cast<const int8_t*>(p);
  int64_t byte = *ptr;
  if (byte >= 0) {
    *OUTPUT = static_cast<uint64_t>(byte);
    return reinterpret_cast<const char*>(ptr) + 1;
  } else {
    auto tmp = Parse64FallbackPair(p, byte);
    if (ABSL_PREDICT_TRUE(tmp.first)) *OUTPUT = tmp.second;
    return tmp.first;
  }
#else
  const unsigned char* ptr = reinterpret_cast<const unsigned char*>(p);
  uint32_t byte = *ptr;
  if (byte < 128) {
    *OUTPUT = byte;
    return reinterpret_cast<const char*>(ptr) + 1;
  } else {
    return Parse64Fallback(p, OUTPUT);
  }
#endif
 }

inline const char* Varint::Skip64(const char* p) {
  const unsigned char* ptr = reinterpret_cast<const unsigned char*>(p);
  if (*ptr++ < 128) return reinterpret_cast<const char*>(ptr);
  if (*ptr++ < 128) return reinterpret_cast<const char*>(ptr);
  if (*ptr++ < 128) return reinterpret_cast<const char*>(ptr);
  if (*ptr++ < 128) return reinterpret_cast<const char*>(ptr);
  if (*ptr++ < 128) return reinterpret_cast<const char*>(ptr);
  if (*ptr++ < 128) return reinterpret_cast<const char*>(ptr);
  if (*ptr++ < 128) return reinterpret_cast<const char*>(ptr);
  if (*ptr++ < 128) return reinterpret_cast<const char*>(ptr);
  if (*ptr++ < 128) return reinterpret_cast<const char*>(ptr);
  if (*ptr++ < 2) return reinterpret_cast<const char*>(ptr);
  return nullptr; // value is too long to be a varint64
}

inline const char* Varint::Parse64Backward(const char* p, const char* b,
                                           uint64_t* OUTPUT) {
  if (p > b + kMax64) {
    // Fast path
    const unsigned char* ptr = reinterpret_cast<const unsigned char*>(p);
    uint32_t byte;
    uint64_t res;

    byte = *(--ptr); if (byte > 127) return nullptr;

    res = byte;
    byte = *(--ptr); if (byte < 128) goto done;
    res <<= 7; res |= (byte & 127);
    byte = *(--ptr); if (byte < 128) goto done;
    res <<= 7; res |= (byte & 127);
    byte = *(--ptr); if (byte < 128) goto done;
    res <<= 7; res |= (byte & 127);
    byte = *(--ptr); if (byte < 128) goto done;
    res <<= 7; res |= (byte & 127);
    byte = *(--ptr); if (byte < 128) goto done;
    res <<= 7; res |= (byte & 127);
    byte = *(--ptr); if (byte < 128) goto done;
    res <<= 7; res |= (byte & 127);
    byte = *(--ptr); if (byte < 128) goto done;
    res <<= 7; res |= (byte & 127);
    byte = *(--ptr); if (byte < 128) goto done;
    res <<= 7; res |= (byte & 127);
    byte = *(--ptr); if (byte < 128) goto done;
    res <<= 7; res |= (byte & 127);
    byte = *(--ptr); if (byte < 128) goto done;

    return nullptr;       // Value is too long to be a varint64

 done:
    *OUTPUT = res;
    return reinterpret_cast<const char*>(ptr + 1);
  } else {
    return Parse64BackwardSlow(p, b, OUTPUT);
  }
}

inline const char* Varint::Skip64Backward(const char* p, const char* b) {
  if (p > b + kMax64) {
    // Fast path
    const unsigned char* ptr = reinterpret_cast<const unsigned char*>(p);
    if (*(--ptr) > 127) return nullptr;
    if (*(--ptr) < 128) return reinterpret_cast<const char*>(ptr+1);
    if (*(--ptr) < 128) return reinterpret_cast<const char*>(ptr+1);
    if (*(--ptr) < 128) return reinterpret_cast<const char*>(ptr+1);
    if (*(--ptr) < 128) return reinterpret_cast<const char*>(ptr+1);
    if (*(--ptr) < 128) return reinterpret_cast<const char*>(ptr+1);
    if (*(--ptr) < 128) return reinterpret_cast<const char*>(ptr+1);
    if (*(--ptr) < 128) return reinterpret_cast<const char*>(ptr+1);
    if (*(--ptr) < 128) return reinterpret_cast<const char*>(ptr+1);
    if (*(--ptr) < 128) return reinterpret_cast<const char*>(ptr+1);
    if (*(--ptr) < 128) return reinterpret_cast<const char*>(ptr+1);
    return nullptr; // value is too long to be a varint64
  } else {
    return Skip64BackwardSlow(p, b);
  }
}

template <typename T>
int Varint::LengthImpl(T v) {
  static_assert(std::numeric_limits<T>::digits == 32 ||
                    std::numeric_limits<T>::digits == 64,
                "Unexpected length. This implementation works up to 448 bits, "
                "but performance has not been tested other than 32 and 64.");
  // This computes value == 0 ? 1 : floor(log2(v)) / 7 + 1
  // Rather than using division by 7 to accomplish this, we use multiplication
  // by 9/64. This has a number of important properties:
  //   * It's roughly 1/7.111111. This makes the 0 bits set case have the same
  //     value as the 7 bits set case, so offsetting by 1 gives us the correct
  //     value for integers up to 448 bits.
  //   * Multiplying by 9 is special. x * 9 = x << 3 + x, and so this
  //     multiplication can be done by a single shifted add on arm:
  //       add w0, w0, w0, lsl #3
  //     or a single lea instruction on x86/64:
  //       leal (%rax,%rax,8), %eax)
  //   * Dividing by 64 is a 6 bit right shift.
  //
  // Summarized:
  //
  // value == 0 ? 1 : (floor(log2(v)) / 7 + 1)
  // floor((32 - clz) / 7.1111) + 1
  // ((32 - clz) * 9) / 64 + 1
  // (((32 - clz) * 9) >> 6) + 1
  // ((32 - clz) * 9 + (1 << 6)) >> 6
  // (32 * 9 + 64 - clz * 9) >> 6
  // (352 - clz * 9) >> 6
  // on arm:
  //  (352 - (clz + (clz << 3))) >> 6
  // on x86:
  //  (352 - lea(clz, clz, 8)) >> 6
  uint32_t clz = absl::countl_zero(v);
  return static_cast<int>(
      ((std::numeric_limits<T>::digits * 9 + 64) - (clz * 9)) >> 6);
}

// X86 CPUs lacking the lzcnt instruction are faster with the bsr-based
// implementation. MSVC does not define __LZCNT__, the nearest option that
// it interprets as lzcnt availability is __AVX2__.
#if (defined(__x86__) || defined(__x86_64__) || defined(_M_IX86) || \
     defined(_M_X64)) &&                                            \
    !(defined(__LZCNT__) || defined(__AVX2__))
#define UTIL_CODING_VARINT_H_PREFER_BSR 1
#else
#define UTIL_CODING_VARINT_H_PREFER_BSR 0
#endif
inline int Varint::Length32(uint32_t v) {
#if UTIL_CODING_VARINT_H_PREFER_BSR
  // Use bsr instruction
  uint32_t log2value = static_cast<uint32_t>(Bits::Log2FloorNonZero(v | 0x1));
  return static_cast<int>((log2value * 9 + (64 + 9)) / 64);
#else
  return LengthImpl(v);
#endif
}

inline int Varint::Length64(uint64_t v) {
#if UTIL_CODING_VARINT_H_PREFER_BSR
  // Use bsr instruction
  uint32_t log2value = static_cast<uint32_t>(Bits::Log2FloorNonZero64(v | 0x1));
  return static_cast<int>((log2value * 9 + (64 + 9)) / 64);
#else
  return LengthImpl(v);
#endif
}
#undef UTIL_CODING_VARINT_H_PREFER_BSR

inline void Varint::Append32(std::string* s, uint32_t value) {
  // Inline the fast-path for single-character output, but fall back to the .cc
  // file for the full version. The size<capacity check is so the compiler can
  // optimize out the string resize code.
  if (value < 128 && s->size() < s->capacity()) {
    s->push_back(static_cast<unsigned char>(value));
  } else {
    Append32Slow(s, value);
  }
}

inline void Varint::Append64(std::string* s, uint64_t value) {
  // Inline the fast-path for single-character output, but fall back to the .cc
  // file for the full version. The size<capacity check is so the compiler can
  // optimize out the string resize code.
  if (value < 128 && s->size() < s->capacity()) {
    s->push_back(static_cast<unsigned char>(value));
  } else {
    Append64Slow(s, value);
  }
}

inline char* Varint::Encode32Inline(char* sptr, uint32_t v) {
  // Operate on characters as unsigneds
  uint8_t* ptr = reinterpret_cast<uint8_t*>(sptr);
  static const uint32_t B = 128;
  if (v < (1<<7)) {
    *(ptr++) = static_cast<uint8_t>(v);
  } else if (v < (1<<14)) {
    *(ptr++) = static_cast<uint8_t>(v | B);
    *(ptr++) = static_cast<uint8_t>(v >> 7);
  } else if (v < (1<<21)) {
    *(ptr++) = static_cast<uint8_t>(v | B);
    *(ptr++) = static_cast<uint8_t>((v >> 7) | B);
    *(ptr++) = static_cast<uint8_t>(v >> 14);
  } else if (v < (1<<28)) {
    *(ptr++) = static_cast<uint8_t>(v | B);
    *(ptr++) = static_cast<uint8_t>((v >> 7) | B);
    *(ptr++) = static_cast<uint8_t>((v >> 14) | B);
    *(ptr++) = static_cast<uint8_t>(v >> 21);
  } else {
    *(ptr++) = static_cast<uint8_t>(v | B);
    *(ptr++) = static_cast<uint8_t>((v >> 7) | B);
    *(ptr++) = static_cast<uint8_t>((v >> 14) | B);
    *(ptr++) = static_cast<uint8_t>((v >> 21) | B);
    *(ptr++) = static_cast<uint8_t>(v >> 28);
  }
  return reinterpret_cast<char*>(ptr);
}

#if (-1 >> 1) != -1
#error FastDecodeDeltas() needs right-shift to sign-extend.
#endif
inline const char* Varint::FastDecodeDeltas(const char* ptr, int64_t goal,
                                            int64_t* out) {
  int64_t value;
  int64_t sum = -goal;
  int64_t shift = 0;
  // Make decoding faster by eliminating unpredictable branching.
  do {
    value = static_cast<int8_t>(*ptr++);  // sign extend one byte of data
    sum += (value & 0x7F) << shift;
    shift += 7;
    // (value >> 7) is either -1(continuation byte) or 0 (stop byte)
    shift &= value >> 7;
    // Loop if we haven't reached goal (sum < 0) or we haven't finished
    // parsing current delta (value < 0). We write it in the form of
    // (a | b) < 0 as opposed to (a < 0 || b < 0) as the former one is
    // usually as fast as a test for (a < 0).
  } while ((sum | value) < 0);

  *out = goal + sum;
  return ptr;
}

#endif  // S2_UTIL_CODING_VARINT_H_
