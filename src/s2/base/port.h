// Copyright Google Inc. All Rights Reserved.
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

#ifndef S2_BASE_PORT_H_
#define S2_BASE_PORT_H_

// This file contains things that are not used in third_party/absl but needed by
// - Platform specific requirement
//   - MSVC
// - Utility macros
// - Endianness
// - Hash
// - Global variables
// - Type alias
// - Predefined system/language macros
// - Predefined system/language functions
// - Performance optimization (alignment)
// - Obsolete

#include <cassert>
#include <climits>
#include <cstdlib>
#include <cstring>

#include "s2/base/integral_types.h"
#include "absl/base/config.h"
#include "absl/base/port.h"

#ifdef SWIG
%include "third_party/absl/base/port.h"
#endif

// -----------------------------------------------------------------------------
// MSVC Specific Requirements
// -----------------------------------------------------------------------------

#ifdef _MSC_VER /* if Visual C++ */

#include <winsock2.h>  // Must come before <windows.h>
#include <intrin.h>
#include <process.h>  // _getpid()
#include <windows.h>
#undef ERROR
#undef DELETE
#undef DIFFERENCE
#define STDIN_FILENO 0
#define STDOUT_FILENO 1
#define STDERR_FILENO 2
#define S_IRUSR 00400
#define S_IWUSR 00200
#define S_IXUSR 00100
#define S_IRGRP 00040
#define S_IWGRP 00020
#define S_IXGRP 00010
#define S_IROTH 00004
#define S_IWOTH 00002
#define S_IXOTH 00001

// This compiler flag can be easily overlooked on MSVC.
// _CHAR_UNSIGNED gets set with the /J flag.
#ifndef _CHAR_UNSIGNED
#error chars must be unsigned!  Use the /J flag on the compiler command line.  // NOLINT
#endif

// Allow comparisons between signed and unsigned values.
//
// Lots of Google code uses this pattern:
//   for (int i = 0; i < container.size(); ++i)
// Since size() returns an unsigned value, this warning would trigger
// frequently.  Very few of these instances are actually bugs since containers
// rarely exceed MAX_INT items.  Unfortunately, there are bugs related to
// signed-unsigned comparisons that have been missed because we disable this
// warning.  For example:
//   const long stop_time = os::GetMilliseconds() + kWaitTimeoutMillis;
//   while (os::GetMilliseconds() <= stop_time) { ... }
#pragma warning(disable : 4018)  // level 3
#pragma warning(disable : 4267)  // level 3

// Don't warn about unused local variables.
//
// extension to silence particular instances of this warning.  There's no way
// to define ABSL_ATTRIBUTE_UNUSED to quiet particular instances of this warning
// in VC++, so we disable it globally.  Currently, there aren't many false
// positives, so perhaps we can address those in the future and re-enable these
// warnings, which sometimes catch real bugs.
#pragma warning(disable : 4101)  // level 3

// Allow initialization and assignment to a smaller type without warnings about
// possible loss of data.
//
// There is a distinct warning, 4267, that warns about size_t conversions to
// smaller types, but we don't currently disable that warning.
//
// Correct code can be written in such a way as to avoid false positives
// by making the conversion explicit, but Google code isn't usually that
// verbose.  There are too many false positives to address at this time.  Note
// that this warning triggers at levels 2, 3, and 4 depending on the specific
// type of conversion.  By disabling it, we not only silence minor narrowing
// conversions but also serious ones.
#pragma warning(disable : 4244)  // level 2, 3, and 4

// Allow silent truncation of double to float.
//
// Silencing this warning has caused us to miss some subtle bugs.
#pragma warning(disable : 4305)  // level 1

// Allow a constant to be assigned to a type that is too small.
//
// I don't know why we allow this at all.  I can't think of a case where this
// wouldn't be a bug, but enabling the warning breaks many builds today.
#pragma warning(disable : 4307)  // level 2

// Allow passing the this pointer to an initializer even though it refers
// to an uninitialized object.
//
// Some observer implementations rely on saving the this pointer.  Those are
// safe because the pointer is not dereferenced until after the object is fully
// constructed.  This could however, obscure other instances.  In the future, we
// should look into disabling this warning locally rather globally.
#pragma warning(disable : 4355)  // level 1 and 4

// Allow implicit coercion from an integral type to a bool.
//
// These could be avoided by making the code more explicit, but that's never
// been the style here, so there would be many false positives.  It's not
// obvious if a true positive would ever help to find an actual bug.
#pragma warning(disable : 4800)  // level 3

#endif  // _MSC_VER

// LANG_CXX11
// GXX_EXPERIMENTAL_CXX0X is defined by gcc and clang up to at least
// gcc-4.7 and clang-3.1 (2011-12-13).  __cplusplus was defined to 1
// in gcc before 4.7 (Crosstool 16) and clang before 3.1, but is
// defined according to the language version in effect thereafter.
// Microsoft Visual Studio 14 (2015) sets __cplusplus==199711 despite
// reasonably good C++11 support, so we set LANG_CXX for it and
// newer versions (_MSC_VER >= 1900).
#if (defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L || \
     (defined(_MSC_VER) && _MSC_VER >= 1900))
// DEPRECATED: Do not key off LANG_CXX11. Instead, write more accurate condition
// that checks whether the C++ feature you need is available or missing, and
// define a more specific feature macro (GOOGLE_HAVE_FEATURE_FOO). You can check
// http://en.cppreference.com/w/cpp/compiler_support for compiler support on C++
// features.
// Define this to 1 if the code is compiled in C++11 mode; leave it
// undefined otherwise.  Do NOT define it to 0 -- that causes
// '#ifdef LANG_CXX11' to behave differently from '#if LANG_CXX11'.
#define LANG_CXX11 1
#endif

// This sanity check can be removed when all references to
// LANG_CXX11 is removed from the code base.
#if defined(__cplusplus) && !defined(LANG_CXX11) && !defined(SWIG)
#error "LANG_CXX11 is required."
#endif

// GOOGLE_OBSCURE_SIGNAL
#if defined(__APPLE__)
// No SIGPWR on MacOSX.  SIGINFO seems suitably obscure.
#define GOOGLE_OBSCURE_SIGNAL SIGINFO
#else
/* We use SIGPWR since that seems unlikely to be used for other reasons. */
#define GOOGLE_OBSCURE_SIGNAL SIGPWR
#endif

// -----------------------------------------------------------------------------
// Endianness
// -----------------------------------------------------------------------------

// IS_LITTLE_ENDIAN, IS_BIG_ENDIAN
#if defined __linux__ || defined OS_ANDROID || defined(__ANDROID__)
// TODO(user): http://b/21460321; use one of OS_ANDROID or __ANDROID__.
// _BIG_ENDIAN
#include <endian.h>

#elif defined(__APPLE__)

// BIG_ENDIAN
#include <machine/endian.h>  // NOLINT(build/include)
/* Let's try and follow the Linux convention */
#define __BYTE_ORDER  BYTE_ORDER
#define __LITTLE_ENDIAN LITTLE_ENDIAN
#define __BIG_ENDIAN BIG_ENDIAN

#endif

// defines __BYTE_ORDER
#ifdef _WIN32
#define __BYTE_ORDER __LITTLE_ENDIAN
#define IS_LITTLE_ENDIAN
#else // _WIN32

// define the macros IS_LITTLE_ENDIAN or IS_BIG_ENDIAN
// using the above endian definitions from endian.h if
// endian.h was included
#ifdef __BYTE_ORDER
#if __BYTE_ORDER == __LITTLE_ENDIAN
#define IS_LITTLE_ENDIAN
#endif

#if __BYTE_ORDER == __BIG_ENDIAN
#define IS_BIG_ENDIAN
#endif

#else  // __BYTE_ORDER

#if defined(__LITTLE_ENDIAN__)
#define IS_LITTLE_ENDIAN
#elif defined(__BIG_ENDIAN__)
#define IS_BIG_ENDIAN
#endif

#endif  // __BYTE_ORDER
#endif  // _WIN32

// byte swap functions (bswap_16, bswap_32, bswap_64).

// The following guarantees declaration of the byte swap functions
#ifdef _MSC_VER
#include <cstdlib>  // NOLINT(build/include)
#define bswap_16(x) _byteswap_ushort(x)
#define bswap_32(x) _byteswap_ulong(x)
#define bswap_64(x) _byteswap_uint64(x)

#elif defined(__APPLE__)
// Mac OS X / Darwin features
#include <libkern/OSByteOrder.h>
#define bswap_16(x) OSSwapInt16(x)
#define bswap_32(x) OSSwapInt32(x)
#define bswap_64(x) OSSwapInt64(x)

#elif defined(__GLIBC__) || defined(__BIONIC__) || defined(__ASYLO__)
#include <byteswap.h>  // IWYU pragma: export

#else

static inline uint16 bswap_16(uint16 x) {
#ifdef __cplusplus
  return static_cast<uint16>(((x & 0xFF) << 8) | ((x & 0xFF00) >> 8));
#else
  return (uint16)(((x & 0xFF) << 8) | ((x & 0xFF00) >> 8));  // NOLINT
#endif  // __cplusplus
}
#define bswap_16(x) bswap_16(x)
static inline uint32 bswap_32(uint32 x) {
  return (((x & 0xFF) << 24) |
          ((x & 0xFF00) << 8) |
          ((x & 0xFF0000) >> 8) |
          ((x & 0xFF000000) >> 24));
}
#define bswap_32(x) bswap_32(x)
static inline uint64 bswap_64(uint64 x) {
  return (((x & 0xFFULL) << 56) |
          ((x & 0xFF00ULL) << 40) |
          ((x & 0xFF0000ULL) << 24) |
          ((x & 0xFF000000ULL) << 8) |
          ((x & 0xFF00000000ULL) >> 8) |
          ((x & 0xFF0000000000ULL) >> 24) |
          ((x & 0xFF000000000000ULL) >> 40) |
          ((x & 0xFF00000000000000ULL) >> 56));
}
#define bswap_64(x) bswap_64(x)

#endif

// printf macros
// __STDC_FORMAT_MACROS must be defined before inttypes.h inclusion */
#if defined(__APPLE__)
/* From MacOSX's inttypes.h:
 * "C++ implementations should define these macros only when
 *  __STDC_FORMAT_MACROS is defined before <inttypes.h> is included." */
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif /* __STDC_FORMAT_MACROS */
#endif /* __APPLE__ */

// printf macros for size_t, in the style of inttypes.h
#if defined(_LP64) || defined(__APPLE__)
#define __PRIS_PREFIX "z"
#else
#define __PRIS_PREFIX
#endif

// Use these macros after a % in a printf format string
// to get correct 32/64 bit behavior, like this:
// size_t size = records.size();
// printf("%" PRIuS "\n", size);
#define PRIdS __PRIS_PREFIX "d"
#define PRIxS __PRIS_PREFIX "x"
#define PRIuS __PRIS_PREFIX "u"
#define PRIXS __PRIS_PREFIX "X"
#define PRIoS __PRIS_PREFIX "o"

// -----------------------------------------------------------------------------
// Performance Optimization
// -----------------------------------------------------------------------------

// Alignment

// Unaligned APIs

// Portable handling of unaligned loads, stores, and copies.
// On some platforms, like ARM, the copy functions can be more efficient
// then a load and a store.
//
// It is possible to implement all of these these using constant-length memcpy
// calls, which is portable and will usually be inlined into simple loads and
// stores if the architecture supports it. However, such inlining usually
// happens in a pass that's quite late in compilation, which means the resulting
// loads and stores cannot participate in many other optimizations, leading to
// overall worse code.
// TODO(user): These APIs are forked in Abseil, see
// LLVM, we should reimplement these APIs with functions calling memcpy(), and
// maybe publish them in Abseil.

// The unaligned API is C++ only.  The declarations use C++ features
// (namespaces, inline) which are absent or incompatible in C.
#if defined(__cplusplus)

#if defined(ADDRESS_SANITIZER) || defined(THREAD_SANITIZER) || \
    defined(MEMORY_SANITIZER)
// Consider we have an unaligned load/store of 4 bytes from address 0x...05.
// AddressSanitizer will treat it as a 3-byte access to the range 05:07 and
// will miss a bug if 08 is the first unaddressable byte.
// ThreadSanitizer will also treat this as a 3-byte access to 05:07 and will
// miss a race between this access and some other accesses to 08.
// MemorySanitizer will correctly propagate the shadow on unaligned stores
// and correctly report bugs on unaligned loads, but it may not properly
// update and report the origin of the uninitialized memory.
// For all three tools, replacing an unaligned access with a tool-specific
// callback solves the problem.

// Make sure uint16_t/uint32_t/uint64_t are defined.
#include <cstdint>

extern "C" {
uint16_t __sanitizer_unaligned_load16(const void *p);
uint32_t __sanitizer_unaligned_load32(const void *p);
uint64_t __sanitizer_unaligned_load64(const void *p);
void __sanitizer_unaligned_store16(void *p, uint16_t v);
void __sanitizer_unaligned_store32(void *p, uint32_t v);
void __sanitizer_unaligned_store64(void *p, uint64_t v);
}  // extern "C"

inline uint16 UNALIGNED_LOAD16(const void *p) {
  return __sanitizer_unaligned_load16(p);
}

inline uint32 UNALIGNED_LOAD32(const void *p) {
  return __sanitizer_unaligned_load32(p);
}

inline uint64 UNALIGNED_LOAD64(const void *p) {
  return __sanitizer_unaligned_load64(p);
}

inline void UNALIGNED_STORE16(void *p, uint16 v) {
  __sanitizer_unaligned_store16(p, v);
}

inline void UNALIGNED_STORE32(void *p, uint32 v) {
  __sanitizer_unaligned_store32(p, v);
}

inline void UNALIGNED_STORE64(void *p, uint64 v) {
  __sanitizer_unaligned_store64(p, v);
}

#elif defined(UNDEFINED_BEHAVIOR_SANITIZER)

inline uint16 UNALIGNED_LOAD16(const void *p) {
  uint16 t;
  memcpy(&t, p, sizeof t);
  return t;
}

inline uint32 UNALIGNED_LOAD32(const void *p) {
  uint32 t;
  memcpy(&t, p, sizeof t);
  return t;
}

inline uint64 UNALIGNED_LOAD64(const void *p) {
  uint64 t;
  memcpy(&t, p, sizeof t);
  return t;
}

inline void UNALIGNED_STORE16(void *p, uint16 v) { memcpy(p, &v, sizeof v); }

inline void UNALIGNED_STORE32(void *p, uint32 v) { memcpy(p, &v, sizeof v); }

inline void UNALIGNED_STORE64(void *p, uint64 v) { memcpy(p, &v, sizeof v); }

#elif defined(__x86_64__) || defined(_M_X64) || defined(__i386) || \
    defined(_M_IX86) || defined(__ppc__) || defined(__PPC__) ||    \
    defined(__ppc64__) || defined(__PPC64__)

// x86 and x86-64 can perform unaligned loads/stores directly;
// modern PowerPC hardware can also do unaligned integer loads and stores;
// but note: the FPU still sends unaligned loads and stores to a trap handler!

#define UNALIGNED_LOAD16(_p) (*reinterpret_cast<const uint16 *>(_p))
#define UNALIGNED_LOAD32(_p) (*reinterpret_cast<const uint32 *>(_p))
#define UNALIGNED_LOAD64(_p) (*reinterpret_cast<const uint64 *>(_p))

#define UNALIGNED_STORE16(_p, _val) (*reinterpret_cast<uint16 *>(_p) = (_val))
#define UNALIGNED_STORE32(_p, _val) (*reinterpret_cast<uint32 *>(_p) = (_val))
#define UNALIGNED_STORE64(_p, _val) (*reinterpret_cast<uint64 *>(_p) = (_val))

#elif defined(__arm__) && !defined(__ARM_ARCH_5__) &&          \
    !defined(__ARM_ARCH_5T__) && !defined(__ARM_ARCH_5TE__) && \
    !defined(__ARM_ARCH_5TEJ__) && !defined(__ARM_ARCH_6__) && \
    !defined(__ARM_ARCH_6J__) && !defined(__ARM_ARCH_6K__) &&  \
    !defined(__ARM_ARCH_6Z__) && !defined(__ARM_ARCH_6ZK__) && \
    !defined(__ARM_ARCH_6T2__)

// ARMv7 and newer support native unaligned accesses, but only of 16-bit
// and 32-bit values (not 64-bit); older versions either raise a fatal signal,
// do an unaligned read and rotate the words around a bit, or do the reads very
// slowly (trip through kernel mode). There's no simple #define that says just
// “ARMv7 or higher”, so we have to filter away all ARMv5 and ARMv6
// sub-architectures. Newer gcc (>= 4.6) set an __ARM_FEATURE_ALIGNED #define,
// so in time, maybe we can move on to that.
//
// This is a mess, but there's not much we can do about it.
//
// To further complicate matters, only LDR instructions (single reads) are
// allowed to be unaligned, not LDRD (two reads) or LDM (many reads). Unless we
// explicitly tell the compiler that these accesses can be unaligned, it can and
// will combine accesses. On armcc, the way to signal this is done by accessing
// through the type (uint32 __packed *), but GCC has no such attribute
// (it ignores __attribute__((packed)) on individual variables). However,
// we can tell it that a _struct_ is unaligned, which has the same effect,
// so we do that.

namespace base {
namespace internal {

struct Unaligned16Struct {
  uint16 value;
  uint8 dummy;  // To make the size non-power-of-two.
} ABSL_ATTRIBUTE_PACKED;

struct Unaligned32Struct {
  uint32 value;
  uint8 dummy;  // To make the size non-power-of-two.
} ABSL_ATTRIBUTE_PACKED;

}  // namespace internal
}  // namespace base

#define UNALIGNED_LOAD16(_p) \
  ((reinterpret_cast<const ::base::internal::Unaligned16Struct *>(_p))->value)
#define UNALIGNED_LOAD32(_p) \
  ((reinterpret_cast<const ::base::internal::Unaligned32Struct *>(_p))->value)

#define UNALIGNED_STORE16(_p, _val)                                        \
  ((reinterpret_cast< ::base::internal::Unaligned16Struct *>(_p))->value = \
       (_val))
#define UNALIGNED_STORE32(_p, _val)                                        \
  ((reinterpret_cast< ::base::internal::Unaligned32Struct *>(_p))->value = \
       (_val))

// TODO(user): NEON supports unaligned 64-bit loads and stores.
// See if that would be more efficient on platforms supporting it,
// at least for copies.

inline uint64 UNALIGNED_LOAD64(const void *p) {
  uint64 t;
  memcpy(&t, p, sizeof t);
  return t;
}

inline void UNALIGNED_STORE64(void *p, uint64 v) { memcpy(p, &v, sizeof v); }

#else

#define NEED_ALIGNED_LOADS

// These functions are provided for architectures that don't support
// unaligned loads and stores.

inline uint16 UNALIGNED_LOAD16(const void *p) {
  uint16 t;
  memcpy(&t, p, sizeof t);
  return t;
}

inline uint32 UNALIGNED_LOAD32(const void *p) {
  uint32 t;
  memcpy(&t, p, sizeof t);
  return t;
}

inline uint64 UNALIGNED_LOAD64(const void *p) {
  uint64 t;
  memcpy(&t, p, sizeof t);
  return t;
}

inline void UNALIGNED_STORE16(void *p, uint16 v) { memcpy(p, &v, sizeof v); }

inline void UNALIGNED_STORE32(void *p, uint32 v) { memcpy(p, &v, sizeof v); }

inline void UNALIGNED_STORE64(void *p, uint64 v) { memcpy(p, &v, sizeof v); }

#endif

// The UNALIGNED_LOADW and UNALIGNED_STOREW macros load and store values
// of type uword_t.
#ifdef _LP64
#define UNALIGNED_LOADW(_p) UNALIGNED_LOAD64(_p)
#define UNALIGNED_STOREW(_p, _val) UNALIGNED_STORE64(_p, _val)
#else
#define UNALIGNED_LOADW(_p) UNALIGNED_LOAD32(_p)
#define UNALIGNED_STOREW(_p, _val) UNALIGNED_STORE32(_p, _val)
#endif

inline void UnalignedCopy16(const void *src, void *dst) {
  UNALIGNED_STORE16(dst, UNALIGNED_LOAD16(src));
}

inline void UnalignedCopy32(const void *src, void *dst) {
  UNALIGNED_STORE32(dst, UNALIGNED_LOAD32(src));
}

inline void UnalignedCopy64(const void *src, void *dst) {
  if (sizeof(void *) == 8) {
    UNALIGNED_STORE64(dst, UNALIGNED_LOAD64(src));
  } else {
    const char *src_char = reinterpret_cast<const char *>(src);
    char *dst_char = reinterpret_cast<char *>(dst);

    UNALIGNED_STORE32(dst_char, UNALIGNED_LOAD32(src_char));
    UNALIGNED_STORE32(dst_char + 4, UNALIGNED_LOAD32(src_char + 4));
  }
}

#endif  // defined(__cplusplus), end of unaligned API

#endif  // S2_BASE_PORT_H_
