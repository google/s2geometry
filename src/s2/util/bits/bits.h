// Copyright 2002 Google Inc. All Rights Reserved.
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

#ifndef S2_UTIL_BITS_BITS_H_
#define S2_UTIL_BITS_BITS_H_

//
// Various bit-twiddling functions, all of which are static members of the Bits
// class (making it effectively a namespace). Operands are unsigned integers.
// Munging bits in _signed_ integers is fraught with peril! For example,
// -5 << n has undefined behavior (for some values of n).
//
// Bits provide the following:
//
//   * Count(Ones.*|LeadingZeros.*)? . In a similar vein, there's also the
//     Find[LM]SBSetNonZero.* family of functions. You can think of them as
//     (trailing|leading) zero bit count + 1. Also in a similar vein,
//     (Capped)?Difference, which count the number of one bits in foo ^ bar.
//
//   * ReverseBits${power_of_two}
//
//   * Log2(Floor|Ceiling)(NonZero)?.* - The NonZero variants have undefined
//     behavior if argument is 0.
//
//   * Bytes(ContainByte(LessThan)?|AllInRange) - These scan a sequence of bytes
//     looking for one with(out)? some property.
//
//   * (Get|Set|Copy)Bits
//
// The only other thing is BitPattern, which is a trait class template (not in
// Bits) containing a few bit patterns (which vary based on value of template
// parameter).

#if defined(__i386__) || defined(__x86_64__)
#include <x86intrin.h>
#endif

#include <glog/logging.h>

#include "s2/base/casts.h"
#include "s2/base/int128.h"
#include "s2/base/macros.h"
#include "s2/third_party/absl/base/integral_types.h"
#include "s2/third_party/absl/base/port.h"

class Bits {
 public:
  // A traits class template for unsigned integer type sizes. Primary
  // information contained herein is corresponding (unsigned) integer type.
  // E.g. UnsignedTypeBySize<32>::Type is uint32. Used by UnsignedType.
  template<int size /* in bytes */>
  struct UnsignedTypeBySize;

  // Auxilliary struct for figuring out an unsigned type for a given type.
  template<typename T> struct UnsignedType {
    typedef typename UnsignedTypeBySize<sizeof(T)>::Type Type;
  };

  // Return the number of one bits in the given integer.
  static int CountOnesInByte(unsigned char n);

  static int CountOnes(uint32 n) {
#if defined(__powerpc64__) && defined(__GNUC__)
    // Use popcount builtin if we know it is inlined and fast.
    return PopcountWithBuiltin<uint32>(n);
#elif (defined(__i386__) || defined(__x86_64__)) && defined(__POPCNT__) && \
    defined(__GNUC__)
    return PopcountWithBuiltin<uint32>(n);
#else
    n -= ((n >> 1) & 0x55555555);
    n = ((n >> 2) & 0x33333333) + (n & 0x33333333);
    return static_cast<int>((((n + (n >> 4)) & 0xF0F0F0F) * 0x1010101) >> 24);
#endif
  }

  // Count bits using sideways addition [WWG'57]. See Knuth TAOCP v4 7.1.3(59)
  static inline int CountOnes64(uint64 n) {
#if defined(__powerpc64__) && defined(__GNUC__)
    return PopcountWithBuiltin<uint64>(n);
#elif defined(__x86_64__) && defined(__POPCNT__) && defined(__GNUC__)
    return PopcountWithBuiltin<uint64>(n);
#elif defined(_LP64)
    n -= (n >> 1) & 0x5555555555555555ULL;
    n = ((n >> 2) & 0x3333333333333333ULL) + (n & 0x3333333333333333ULL);
    return static_cast<int>(
        (((n + (n >> 4)) & 0xF0F0F0F0F0F0F0FULL) * 0x101010101010101ULL) >> 56);
#else
    return CountOnes(n >> 32) + CountOnes(n & 0xffffffff);
#endif
  }

  // Count bits in uint128
  static inline int CountOnes128(uint128 n) {
    return Bits::CountOnes64(Uint128High64(n)) +
           Bits::CountOnes64(Uint128Low64(n));
  }

  // Count bits using popcnt instruction (available on argo machines).
  // Doesn't check if the instruction exists.
  // Please use TestCPUFeature(POPCNT) from base/cpuid/cpuid.h before using
  // this.
  static inline int CountOnes64withPopcount(uint64 n) {
#if defined(__x86_64__) && defined __GNUC__
    int64 count = 0;
    asm("popcnt %1,%0" : "=r"(count) : "rm"(n) : "cc");
    return static_cast<int>(count);
#else
    return CountOnes64(n);
#endif
  }

  // Count leading zeroes.  This is similar to wordsize - 1 - floor(log2(n)).
  // Returns number of bits if n is 0.
  static inline int CountLeadingZeros32(uint32 n) {
    // Instead of using __builtin_clz(), we explicitly write target specific
    // assembly because we want to handle n == 0.  If we used __builtin_clz(),
    // we would need to use something like "n ? __builtin_clz(n) : 32".  The
    // check is not necessary on POWER and aarch64 but we cannot depend on
    // that becasue __builtin_clz(0) is documented to be undefined.
#if defined(__aarch64__) && defined(__GNUC__)
    int32 count;
    asm("clz %w0,%w1" : "=r"(count) : "r"(n));
    return count;
#elif (defined(__i386__) || defined(__x86_64__)) && defined(__LZCNT__) && \
    defined(__GNUC__)
    return __lzcnt32(n);
#elif (defined(__i386__) || defined(__x86_64__)) && defined(__GNUC__)
    if (n == 0) return 32;
    int32 idx;
    asm("bsr %1, %0"
        : "=r"(idx)
        : "ro"(n)
        : "cc");              // bsr writes Z flag
    return 31 ^ idx;
#elif defined(__powerpc64__) && defined(__GNUC__)
    int32 count;
    asm("cntlzw %0,%1" : "=r"(count) : "r"(n));
    return count;
#elif defined(__GNUC__)
    return CountLeadingZerosWithBuiltin<uint32>(n);
#else
    return CountLeadingZeros32_Portable(n);
#endif
  }

  static inline int CountLeadingZeros64(uint64 n) {
#if defined(__aarch64__) && defined(__GNUC__)
    int64 count;
    asm("clz %0,%1" : "=r"(count) : "r"(n));
    return count;
#elif defined(__powerpc64__) && defined(__GNUC__)
    int64 count;
    asm("cntlzd %0,%1" : "=r"(count) : "r"(n));
    return count;
#elif (defined(__i386__) || defined(__x86_64__)) && defined(__LZCNT__) && \
    defined(__GNUC__)
    return __lzcnt64(n);
#elif defined(__x86_64__) && defined(__GNUC__)
    if (n == 0) return 64;
    int64 idx;
    asm ("bsr %1, %0"
         : "=r"(idx)
         : "ro"(n)
         : "cc");              // bsr writes Z flag
    return static_cast<int>(63 ^ idx);
#elif defined(__GNUC__)
    return CountLeadingZerosWithBuiltin<uint64>(n);
#else
    return CountLeadingZeros64_Portable(n);
#endif
  }

  static inline int CountLeadingZeros128(uint128 n) {
    if (uint64 hi = Uint128High64(n))
      return Bits::CountLeadingZeros64(hi);
    return Bits::CountLeadingZeros64(Uint128Low64(n)) + 64;
  }

  // Reverse the bits in the given integer.
  static uint8 ReverseBits8(uint8 n);
  static uint32 ReverseBits32(uint32 n);
  static uint64 ReverseBits64(uint64 n);
  static uint128 ReverseBits128(uint128 n);

  // Return the number of one bits in the byte sequence.
  static int Count(const void *m, int num_bytes);

  // Return the number of different bits in the given byte sequences.
  // (i.e., the Hamming distance)
  static int Difference(const void *m1, const void *m2, int num_bytes);

  // Return the number of different bits in the given byte sequences,
  // up to a maximum.  Values larger than the maximum may be returned
  // (because multiple bits are checked at a time), but the function
  // may exit early if the cap is exceeded.
  static int CappedDifference(const void *m1, const void *m2,
                              int num_bytes, int cap);

  // Return floor(log2(n)) for positive integer n.  Returns -1 iff n == 0.
  static int Log2Floor(uint32 n);
  static int Log2Floor64(uint64 n);
  static int Log2Floor128(uint128 n);

  // Potentially faster version of Log2Floor() that returns an
  // undefined value if n == 0
  static int Log2FloorNonZero(uint32 n);
  static int Log2FloorNonZero64(uint64 n);
  static int Log2FloorNonZero128(uint128 n);

  // Return ceiling(log2(n)) for positive integer n.  Returns -1 iff n == 0.
  static int Log2Ceiling(uint32 n);
  static int Log2Ceiling64(uint64 n);
  static int Log2Ceiling128(uint128 n);

  // Return the first set least / most significant bit, 0-indexed.  Returns an
  // undefined value if n == 0.  FindLSBSetNonZero() is similar to ffs() except
  // that it's 0-indexed, while FindMSBSetNonZero() is the same as
  // Log2FloorNonZero().
  static int FindLSBSetNonZero(uint32 n);
  static int FindLSBSetNonZero64(uint64 n);
  static int FindLSBSetNonZero128(uint128 n);
  static int FindMSBSetNonZero(uint32 n) { return Log2FloorNonZero(n); }
  static int FindMSBSetNonZero64(uint64 n) { return Log2FloorNonZero64(n); }
  static int FindMSBSetNonZero128(uint128 n) {
    return Log2FloorNonZero128(n);
  }

  // Viewing bytes as a stream of unsigned bytes, does that stream
  // contain any byte equal to c?
  template <class T> static bool BytesContainByte(T bytes, uint8 c);

  // Viewing bytes as a stream of unsigned bytes, does that stream
  // contain any byte b < c?
  template <class T> static bool BytesContainByteLessThan(T bytes, uint8 c);

  // Viewing bytes as a stream of unsigned bytes, are all elements of that
  // stream in [lo, hi]?
  template <class T> static bool BytesAllInRange(T bytes, uint8 lo, uint8 hi);

  // Extract 'nbits' consecutive bits from 'src'.  Position of bits are
  // specified by 'offset' from the LSB.  'T' is a scalar type (integral,
  // float or pointer) whose size is the same as one of the unsigned types.
  // The return type is an unsigned type having the same size as T.
  template<typename T>
  static typename UnsignedType<T>::Type GetBits(const T src,
                                                const int offset,
                                                const int nbits) {
    typedef typename UnsignedType<T>::Type UnsignedT;
    const UnsignedT unsigned_src = bit_cast<UnsignedT>(src);
    DCHECK_GT(sizeof(UnsignedT) * 8, offset);
    DCHECK_GE(sizeof(UnsignedT) * 8, offset + nbits);
    const UnsignedT result =
        (unsigned_src >> offset) & NBitsFromLSB<UnsignedT>(nbits);
    return result;
  }

  // Overwrite 'nbits' consecutive bits of 'dest.'.  Position of bits are
  // specified by an offset from the LSB.  'T' is a scalar type (integral,
  // float or pointer) whose size is the same as one of the unsigned types.
  template<typename T>
  static void SetBits(const typename UnsignedType<T>::Type value,
                      const int offset,
                      const int nbits,
                      T* const dest) {
    typedef typename UnsignedType<T>::Type UnsignedT;
    const UnsignedT unsigned_dest = bit_cast<UnsignedT>(*dest);
    DCHECK_GT(sizeof(UnsignedT) * 8, offset);
    DCHECK_GE(sizeof(UnsignedT) * 8, offset + nbits);
    const UnsignedT mask = NBitsFromLSB<UnsignedT>(nbits);
    const UnsignedT unsigned_result =
        (unsigned_dest & ~(mask << offset)) | ((value & mask) << offset);
    *dest = bit_cast<T>(unsigned_result);
  }

  // Combine SetBits and GetBits for convenience.  This is meant to be a
  // replacement for BitCopy() for some use cases.  Unlike BitCopy(),
  // Bits::CopyBits() operating on multibyte types has the same behavior on
  // big-endian and little-endian machines. Sample usage:
  //
  // uint32 a, b;
  // Bits::CopyBits(&a, 0, b, 12, 3);
  template<typename DestType, typename SrcType>
  static void CopyBits(DestType* const dest,
                       const int dest_offset,
                       const SrcType src,
                       const int src_offset,
                       const int nbits) {
    const typename UnsignedType<SrcType>::Type value =
        GetBits(src, src_offset, nbits);
    SetBits(value, dest_offset, nbits, dest);
  }

 private:
  // We only use this for unsigned types and for 0 <= n <= sizeof(UnsignedT).
  template<typename UnsignedT>
  static UnsignedT NBitsFromLSB(const int nbits) {
    const UnsignedT all_ones = ~static_cast<UnsignedT>(0);
    return nbits == 0 ? static_cast<UnsignedT>(0)
                      : all_ones >> (sizeof(UnsignedT) * 8 - nbits);
  }

#ifdef __GNUC__
  template<typename UnsignedT>
  static inline int CountLeadingZerosWithBuiltin(UnsignedT n);
  template<typename UnsignedT>
  static inline int PopcountWithBuiltin(UnsignedT n);
#endif

  // Portable implementations.
  static int Log2Floor_Portable(uint32 n);
  static int Log2Floor64_Portable(uint64 n);
  static int Log2FloorNonZero_Portable(uint32 n);
  static int Log2FloorNonZero64_Portable(uint64 n);
  static int CountLeadingZeros32_Portable(uint32 n);
  static int CountLeadingZeros64_Portable(uint64 n);
  static int FindLSBSetNonZero_Portable(uint32 n);
  static int FindLSBSetNonZero64_Portable(uint64 n);

  static const char num_bits[];
  Bits(Bits const&) = delete;
  void operator=(Bits const&) = delete;

  // Allow tests to call _Portable variants directly.
  // Originally, I wanted to depend on //testing/production_stub/public
  // so that I would be able to
  // #include "s2/testing/production_stub/public/gunit_prod.h", which provides
  // the previously mentioned header file says to instead
  // #include "s2/testing/base/gunit_prod.h". I cannot find a library that a)
  // provides the alternative header file b) is visisble. Thus, I have thrown
  // my hands up in frustration, and gone with raw friend instead of using the
  // usual FRIEND_TEST macro. Hate the game, not the player.
  friend class Bits_Port32_Test;
  friend class Bits_Port64_Test;
};

// A utility class for some handy bit patterns.  The names l and h
// were chosen to match Knuth Volume 4: l is 0x010101... and h is 0x808080...;
// half_ones is ones in the lower half only.  We assume sizeof(T) is 1 or even.
template <class T> struct BitPattern {
  static const T half_ones = (static_cast<T>(1) << (sizeof(T)*4)) - 1;
  static const T l = (sizeof(T) == 1) ? 1 :
                       (half_ones / 0xff * (half_ones + 2));
  static const T h = ~(l * 0x7f);
};

// ------------------------------------------------------------------------
// Implementation details follow
// ------------------------------------------------------------------------

#if defined(__GNUC__)

inline int Bits::Log2Floor(uint32 n) {
  return n == 0 ? -1 : 31 ^ __builtin_clz(n);
}

inline int Bits::Log2FloorNonZero(uint32 n) {
  return 31 ^ __builtin_clz(n);
}

inline int Bits::FindLSBSetNonZero(uint32 n) {
  return __builtin_ctz(n);
}

inline int Bits::Log2Floor64(uint64 n) {
  return n == 0 ? -1 : 63 ^ __builtin_clzll(n);
}

inline int Bits::Log2FloorNonZero64(uint64 n) {
  return 63 ^ __builtin_clzll(n);
}

inline int Bits::FindLSBSetNonZero64(uint64 n) {
  return __builtin_ctzll(n);
}

#elif defined(_MSC_VER)

inline int Bits::FindLSBSetNonZero(uint32 n) {
  return Bits::FindLSBSetNonZero_Portable(n);
}

inline int Bits::FindLSBSetNonZero64(uint64 n) {
  return Bits::FindLSBSetNonZero64_Portable(n);
}

inline int Bits::Log2FloorNonZero(uint32 n) {
#ifdef _M_IX86
  _asm {
    bsr ebx, n
    mov n, ebx
  }
  return n;
#else
  return Bits::Log2FloorNonZero_Portable(n);
#endif
}

inline int Bits::Log2Floor(uint32 n) {
#ifdef _M_IX86
  _asm {
    xor ebx, ebx
    mov eax, n
    and eax, eax
    jz return_ebx
    bsr ebx, eax
return_ebx:
    mov n, ebx
  }
  return n;
#else
  return Bits::Log2Floor_Portable(n);
#endif
}

inline int Bits::Log2Floor64(uint64 n) {
  return Bits::Log2Floor64_Portable(n);
}

inline int Bits::Log2FloorNonZero64(uint64 n) {
  return Bits::Log2FloorNonZero64_Portable(n);
}

#else  // !__GNUC__ && !_MSC_VER

inline int Bits::Log2Floor(uint32 n) {
  return Bits::Log2Floor_Portable(n);
}

inline int Bits::Log2FloorNonZero(uint32 n) {
  return Bits::Log2FloorNonZero_Portable(n);
}

inline int Bits::FindLSBSetNonZero(uint32 n) {
  return Bits::FindLSBSetNonZero_Portable(n);
}

inline int Bits::Log2Floor64(uint64 n) {
  return Bits::Log2Floor64_Portable(n);
}

inline int Bits::Log2FloorNonZero64(uint64 n) {
  return Bits::Log2FloorNonZero64_Portable(n);
}

inline int Bits::FindLSBSetNonZero64(uint64 n) {
  return Bits::FindLSBSetNonZero64_Portable(n);
}

#endif

inline int Bits::Log2Floor128(uint128 n) {
  if (uint64 hi = Uint128High64(n))
    return 64 + Log2FloorNonZero64(hi);
  return Log2Floor64(Uint128Low64(n));
}

inline int Bits::Log2FloorNonZero128(uint128 n) {
  if (uint64 hi = Uint128High64(n))
    return 64 + Log2FloorNonZero64(hi);
  return Log2FloorNonZero64(Uint128Low64(n));
}

inline int Bits::FindLSBSetNonZero128(uint128 n) {
  if (uint64 lo = Uint128Low64(n))
    return Bits::FindLSBSetNonZero64(lo);
  return 64 + Bits::FindLSBSetNonZero64(Uint128High64(n));
}

inline int Bits::CountOnesInByte(unsigned char n) {
  return num_bits[n];
}

inline uint8 Bits::ReverseBits8(unsigned char n) {
#if defined(__aarch64__) && defined(__GNUC__)
  // aarch64 has a reverse bits instruction but there is no gcc builtin.
  uint32 result;
  const uint32 n_shifted = static_cast<uint32>(n) << 24;
  asm("rbit %w0, %w1" : "=r"(result) : "r"(n_shifted));
  return result;
#else
  n = static_cast<unsigned char>(((n >> 1) & 0x55) | ((n & 0x55) << 1));
  n = static_cast<unsigned char>(((n >> 2) & 0x33) | ((n & 0x33) << 2));
  return static_cast<unsigned char>(((n >> 4) & 0x0f)  | ((n & 0x0f) << 4));
#endif
}

inline uint32 Bits::ReverseBits32(uint32 n) {
#if defined(__aarch64__) && defined(__GNUC__)
  uint32 result;
  asm("rbit %w0, %w1" : "=r"(result) : "r"(n));
  return result;
#else
  n = ((n >> 1) & 0x55555555) | ((n & 0x55555555) << 1);
  n = ((n >> 2) & 0x33333333) | ((n & 0x33333333) << 2);
  n = ((n >> 4) & 0x0F0F0F0F) | ((n & 0x0F0F0F0F) << 4);
  return bswap_32(n);
#endif
}

inline uint64 Bits::ReverseBits64(uint64 n) {
#if defined(__aarch64__) && defined(__GNUC__)
  uint64 result;
  asm("rbit %0, %1" : "=r"(result) : "r"(n));
  return result;
#elif defined(_LP64)
  n = ((n >> 1) & 0x5555555555555555ULL) | ((n & 0x5555555555555555ULL) << 1);
  n = ((n >> 2) & 0x3333333333333333ULL) | ((n & 0x3333333333333333ULL) << 2);
  n = ((n >> 4) & 0x0F0F0F0F0F0F0F0FULL) | ((n & 0x0F0F0F0F0F0F0F0FULL) << 4);
  return bswap_64(n);
#else
  return ReverseBits32( n >> 32 ) |
         (static_cast<uint64>(ReverseBits32(n &  0xffffffff)) << 32);
#endif
}

inline uint128 Bits::ReverseBits128(uint128 n) {
  return absl::MakeUint128(ReverseBits64(Uint128Low64(n)),
                           ReverseBits64(Uint128High64(n)));
}

inline int Bits::Log2FloorNonZero_Portable(uint32 n) {
  // Just use the common routine
  return Log2Floor(n);
}

// Log2Floor64() is defined in terms of Log2Floor32(), Log2FloorNonZero32()
inline int Bits::Log2Floor64_Portable(uint64 n) {
  const uint32 topbits = static_cast<uint32>(n >> 32);
  if (topbits == 0) {
    // Top bits are zero, so scan in bottom bits
    return Log2Floor(static_cast<uint32>(n));
  } else {
    return 32 + Log2FloorNonZero(topbits);
  }
}

// Log2FloorNonZero64() is defined in terms of Log2FloorNonZero32()
inline int Bits::Log2FloorNonZero64_Portable(uint64 n) {
  const uint32 topbits = static_cast<uint32>(n >> 32);
  if (topbits == 0) {
    // Top bits are zero, so scan in bottom bits
    return Log2FloorNonZero(static_cast<uint32>(n));
  } else {
    return 32 + Log2FloorNonZero(topbits);
  }
}

// FindLSBSetNonZero64() is defined in terms of FindLSBSetNonZero()
inline int Bits::FindLSBSetNonZero64_Portable(uint64 n) {
  const uint32 bottombits = static_cast<uint32>(n);
  if (bottombits == 0) {
    // Bottom bits are zero, so scan in top bits
    return 32 + FindLSBSetNonZero(static_cast<uint32>(n >> 32));
  } else {
    return FindLSBSetNonZero(bottombits);
  }
}

template <class T>
inline bool Bits::BytesContainByteLessThan(T bytes, uint8 c) {
  T l = BitPattern<T>::l;
  T h = BitPattern<T>::h;
  // The c <= 0x80 code is straight out of Knuth Volume 4.
  // Usually c will be manifestly constant.
  return c <= 0x80 ?
      ((h & (bytes - l * c) & ~bytes) != 0) :
      ((((bytes - l * c) | (bytes ^ h)) & h) != 0);
}

template <class T> inline bool Bits::BytesContainByte(T bytes, uint8 c) {
  // Usually c will be manifestly constant.
  return Bits::BytesContainByteLessThan<T>(bytes ^ (c * BitPattern<T>::l), 1);
}

template <class T>
inline bool Bits::BytesAllInRange(T bytes, uint8 lo, uint8 hi) {
  T l = BitPattern<T>::l;
  T h = BitPattern<T>::h;
  // In the common case, lo and hi are manifest constants.
  if (lo > hi) {
    return false;
  }
  if (hi - lo < 128) {
    T x = bytes - l * lo;
    T y = bytes + l * (127 - hi);
    return ((x | y) & h) == 0;
  }
  return !Bits::BytesContainByteLessThan(bytes + (255 - hi) * l,
                                         lo + (255 - hi));
}

// Specializations for Bits::UnsignedTypeBySize.  For unsupported type
// sizes, a compile-time error will be generated.
template<>
struct Bits::UnsignedTypeBySize<1> {
  typedef uint8 Type;
};

template<>
struct Bits::UnsignedTypeBySize<2> {
  typedef uint16 Type;
};

template<>
struct Bits::UnsignedTypeBySize<4> {
  typedef uint32 Type;
};

template<>
struct Bits::UnsignedTypeBySize<8> {
  typedef uint64 Type;
};

template<>
struct Bits::UnsignedTypeBySize<16> {
  typedef uint128 Type;
};

#ifdef __GNUC__
// Select a correct gcc/llvm builtin for clz based on size.
template<typename UnsignedT>
int Bits::CountLeadingZerosWithBuiltin(UnsignedT n) {
  if (n == 0)   // __builtin_clz(0) is officially undefined.
    return sizeof(n) * 8;
  constexpr size_t size = sizeof(n);
  static_assert(size == sizeof(unsigned)
                || size == sizeof(unsigned long)        // NOLINT(runtime/int)
                || size == sizeof(unsigned long long),  // NOLINT(runtime/int)
                "No __builtin_clz for this size");
  if (size == sizeof(unsigned)) {
    return __builtin_clz(n);
  } else if (size == sizeof(unsigned long)) {           // NOLINT(runtime/int)
    return __builtin_clzl(n);
  } else {
    return __builtin_clzll(n);
  }
}

// Select a correct gcc/llvm builtin for popcount based on size.
template<typename UnsignedT>
int Bits::PopcountWithBuiltin(UnsignedT n) {
  constexpr size_t size = sizeof(n);
  static_assert(size == sizeof(unsigned)
                || size == sizeof(unsigned long)        // NOLINT(runtime/int)
                || size == sizeof(unsigned long long),  // NOLINT(runtime/int)
                "No __builtin_popcount for this size");
  if (size == sizeof(unsigned)) {
    return __builtin_popcount(n);
  } else if (size == sizeof(unsigned long)) {           // NOLINT(runtime/int)
    return __builtin_popcountl(n);
  } else {
    return __builtin_popcountll(n);
  }
}
#endif  // __GNUC__

#endif  // S2_UTIL_BITS_BITS_H_
