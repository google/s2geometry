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
// Various portability macros, type definitions, and inline functions
// This file is used for both C and C++!
//
// These are weird things we need to do to get this compiling on
// random systems (and on SWIG).
//

#ifndef S2_THIRD_PARTY_ABSL_BASE_PORT_H_
#define S2_THIRD_PARTY_ABSL_BASE_PORT_H_

#include <climits>         // So we can set the bounds of our types
#include <cstring>         // for memcpy()
#include <cstdlib>         // for free()

#include "s2/third_party/absl/base/config.h"

#if defined(__cplusplus)
#include <cstddef>
#if defined(_STLPORT_VERSION)
#error "STLPort is not supported."
#endif
#endif

// Enforce C++11 as the minimum.  Note that Visual Studio has not
// advanced __cplusplus despite being good enough for our purposes, so
// so we exempt it from the check.
#if defined(__cplusplus) && !defined(_MSC_VER) && !defined(SWIG)
#if __cplusplus < 201103L
#error "C++ versions less than C++11 are not supported."
#endif
#endif

#if defined(OS_CYGWIN) || defined(__CYGWIN__)
#error "Cygwin is not supported."
#endif

#if defined(__APPLE__)
// traditionally defined __APPLE__ themselves via other build systems, since mac
// TODO(user): Remove this when all toolchains make the proper defines.
#include <TargetConditionals.h>
#if defined(TARGET_OS_IPHONE) && TARGET_OS_IPHONE
#ifndef OS_IOS
#define OS_IOS 1
#endif
#define SUPPRESS_MOBILE_IOS_BASE_PORT_H
#endif  // defined(TARGET_OS_IPHONE) && TARGET_OS_IPHONE
#endif  // defined(__APPLE__)

#if defined(OS_CYGWIN) || defined(__ANDROID__)
#include <malloc.h>         // for memalign()
#elif defined(_MSC_VER)
#include <cstdio>          // declare snprintf/vsnprintf before overriding
#endif

#include "s2/third_party/absl/base/integral_types.h"

// We support gcc 4.7 and later.
#if defined(__GNUC__) && !defined(__clang__)
#if __GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 7)
#error "This package requires gcc 4.7 or higher"
#endif
#endif

// We support MSVC++ 14.0 update 2 and later.
#if defined(_MSC_FULL_VER) && _MSC_FULL_VER < 190023918
#error "This package requires Visual Studio 2015 Update 2 or higher"
#endif

// We support Apple Xcode clang 4.2.1 (version 421.11.65) and later.
// This corresponds to Apple Xcode version 4.5.
#if defined(__apple_build_version__) && __apple_build_version__ < 4211165
#error "This package requires __apple_build_version__ of 4211165 or higher"
#endif

// Must happens before inttypes.h inclusion */
#if defined(__APPLE__)
/* From MacOSX's inttypes.h:
 * "C++ implementations should define these macros only when
 *  __STDC_FORMAT_MACROS is defined before <inttypes.h> is included." */
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif  /* __STDC_FORMAT_MACROS */
#endif  /* __APPLE__ */

/* Default for most OSes */
/* We use SIGPWR since that seems unlikely to be used for other reasons. */
#define GOOGLE_OBSCURE_SIGNAL  SIGPWR

#if defined OS_LINUX || defined OS_CYGWIN || defined OS_ANDROID || \
    defined(__ANDROID__)
// _BIG_ENDIAN
#include <endian.h>
#endif

#if defined OS_LINUX || defined OS_CYGWIN

// GLIBC-related macros.
#include <features.h>

#ifndef __GLIBC_PREREQ
#define __GLIBC_PREREQ(a, b) 0  // not a GLIBC system
#endif

// The uint mess:
// mysql.h sets _GNU_SOURCE which sets __USE_MISC in <features.h>
// sys/types.h typedefs uint if __USE_MISC
// mysql typedefs uint if HAVE_UINT not set
// The following typedef is carefully considered, and should not cause
//  any clashes
#if !defined(__USE_MISC)
#if !defined(HAVE_UINT)
#define HAVE_UINT 1
typedef unsigned int uint;
#endif
#if !defined(HAVE_USHORT)
#define HAVE_USHORT 1
typedef unsigned short ushort;
#endif
#if !defined(HAVE_ULONG)
#define HAVE_ULONG 1
typedef unsigned long ulong;
#endif
#endif

#if !defined(HAVE_TLS) && \
    (defined(GOOGLE_LIBCXX) || defined(_GLIBCXX_HAVE_TLS)) && \
    (defined(ARCH_K8) || defined(ARCH_POWERPC64) || defined(ARCH_PPC) || \
     defined(ARCH_ARM))
#define HAVE_TLS 1
#endif

#elif defined OS_FREEBSD

// _BIG_ENDIAN
#include <machine/endian.h>

#elif defined(__APPLE__) || defined(OS_IOS)

// BIG_ENDIAN
#include <machine/endian.h>  // NOLINT(build/include)
/* Let's try and follow the Linux convention */
#define __BYTE_ORDER  BYTE_ORDER
#define __LITTLE_ENDIAN LITTLE_ENDIAN
#define __BIG_ENDIAN BIG_ENDIAN

#endif

// The following guarantees declaration of the byte swap functions, and
// defines __BYTE_ORDER for MSVC
#ifdef _MSC_VER
#include <cstdlib>  // NOLINT(build/include)
#define __BYTE_ORDER __LITTLE_ENDIAN
#define bswap_16(x) _byteswap_ushort(x)
#define bswap_32(x) _byteswap_ulong(x)
#define bswap_64(x) _byteswap_uint64(x)

#elif defined(__APPLE__) || defined(OS_IOS)
// Mac OS X / Darwin features
#include <libkern/OSByteOrder.h>
#define bswap_16(x) OSSwapInt16(x)
#define bswap_32(x) OSSwapInt32(x)
#define bswap_64(x) OSSwapInt64(x)

#elif defined(__GLIBC__) || defined(__CYGWIN__)
#include <byteswap.h>  // IWYU pragma: export

#else

static inline uint16 bswap_16(uint16 x) {
  return (uint16)(((x & 0xFF) << 8) | ((x & 0xFF00) >> 8));  // NOLINT
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
  return (((x & GG_ULONGLONG(0xFF)) << 56) |
          ((x & GG_ULONGLONG(0xFF00)) << 40) |
          ((x & GG_ULONGLONG(0xFF0000)) << 24) |
          ((x & GG_ULONGLONG(0xFF000000)) << 8) |
          ((x & GG_ULONGLONG(0xFF00000000)) >> 8) |
          ((x & GG_ULONGLONG(0xFF0000000000)) >> 24) |
          ((x & GG_ULONGLONG(0xFF000000000000)) >> 40) |
          ((x & GG_ULONGLONG(0xFF00000000000000)) >> 56));
}
#define bswap_64(x) bswap_64(x)

#endif


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

#else

#if defined(__LITTLE_ENDIAN__)
#define IS_LITTLE_ENDIAN
#elif defined(__BIG_ENDIAN__)
#define IS_BIG_ENDIAN
#endif

// there is also PDP endian ...

#endif  // __BYTE_ORDER

// Define the OS's path separator
#ifdef __cplusplus  // C won't merge duplicate const variables at link time
// Some headers provide a macro for this (GCC's system.h), remove it so that we
// can use our own.
#undef PATH_SEPARATOR
#if defined(OS_WINDOWS)
const char PATH_SEPARATOR = '\\';
#else
const char PATH_SEPARATOR = '/';
#endif
#endif

// Windows has O_BINARY as a flag to open() (like "b" for fopen).
// Linux doesn't need make this distinction.
#if defined OS_LINUX && !defined O_BINARY
#define O_BINARY 0
#endif

#ifdef _MSC_VER
// doesn't have uid_t
typedef int uid_t;
#endif

// Mac OS X / Darwin and iOS features

#if defined(__APPLE__) || defined(OS_IOS)

// For mmap, Linux defines both MAP_ANONYMOUS and MAP_ANON and says MAP_ANON is
// deprecated. In Darwin, MAP_ANON is all there is.
#if !defined MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

// Linux has this in <sys/cdefs.h>
#define __ptr_t void *

// Linux has this in <linux/errno.h>
#define EXFULL      ENOMEM  // not really that great a translation...

// Labeled sections are not supported on Darwin/iOS.
#define HAVE_ATTRIBUTE_SECTION 0

// Darwin doesn't have strnlen. No comment.
inline size_t strnlen(const char *s, size_t maxlen) {
  const char* end = (const char *)memchr(s, '\0', maxlen);
  if (end)
    return end - s;
  return maxlen;
}

// Doesn't exist on OSX.
#define MSG_NOSIGNAL 0

// No SIGPWR on MacOSX.  SIGINFO seems suitably obscure.
#undef GOOGLE_OBSCURE_SIGNAL
#define GOOGLE_OBSCURE_SIGNAL  SIGINFO

#elif defined(OS_CYGWIN)  // Cygwin-specific behavior.

#if defined(__CYGWIN32__)
#define __WORDSIZE 32
#else
// It's probably possible to support 64-bit, but the #defines will need checked.
#error "Cygwin is currently only 32-bit."
#endif

// No signalling on Windows.
#undef GOOGLE_OBSCURE_SIGNAL
#define GOOGLE_OBSCURE_SIGNAL 0

struct stack_t {
  void* ss_sp;
  int ss_flags;
  size_t ss_size;
};
inline int sigaltstack(stack_t* ss, stack_t* oss) { return 0; }

#define PTHREAD_STACK_MIN 0  // Not provided by cygwin

// Scans memory for a character.
// memrchr is used in a few places, but it's linux-specific.
inline void* memrchr(const void* bytes, int find_char, size_t len) {
  const unsigned char* cursor =
      reinterpret_cast<const unsigned char*>(bytes) + len - 1;
  unsigned char actual_char = find_char;
  for (; cursor >= bytes; --cursor) {
    if (*cursor == actual_char) {
      return const_cast<void*>(reinterpret_cast<const void*>(cursor));
    }
  }
  return nullptr;
}

#endif

// Klocwork static analysis tool's C/C++ complier kwcc
#if defined(__KLOCWORK__)
#define STATIC_ANALYSIS
#endif // __KLOCWORK__

// GCC-specific features

#if (defined(__GNUC__) || defined(__APPLE__) || defined(OS_IOS)) && \
    !defined(SWIG)

//
// Tell the compiler to do printf format string checking if the
// compiler supports it; see the 'format' attribute in
// <http://gcc.gnu.org/onlinedocs/gcc-4.3.0/gcc/Function-Attributes.html>.
//
// N.B.: As the GCC manual states, "[s]ince non-static C++ methods
// have an implicit 'this' argument, the arguments of such methods
// should be counted from two, not one."
//
#define PRINTF_ATTRIBUTE(string_index, first_to_check) \
    __attribute__((__format__ (__printf__, string_index, first_to_check)))
#define SCANF_ATTRIBUTE(string_index, first_to_check) \
    __attribute__((__format__ (__scanf__, string_index, first_to_check)))

// Cache line alignment
#if defined(__i386__) || defined(__x86_64__)
#define CACHELINE_SIZE 64
#elif defined(__powerpc64__)
#define CACHELINE_SIZE 128
#elif defined(__aarch64__)
// We would need to read special register ctr_el0 to find out L1 dcache size.
// This value is a good estimate based on a real aarch64 machine.
#define CACHELINE_SIZE 64
#elif defined(__arm__)
// Cache line sizes for ARM: These values are not strictly correct since
// cache line sizes depend on implementations, not architectures.  There
// are even implementations with cache line sizes configurable at boot
// time.
#if defined(__ARM_ARCH_5T__)
#define CACHELINE_SIZE 32
#elif defined(__ARM_ARCH_7A__)
#define CACHELINE_SIZE 64
#endif
#endif

#ifndef CACHELINE_SIZE
// A reasonable default guess.  Note that overestimates tend to waste more
// space, while underestimates tend to waste more time.
#define CACHELINE_SIZE 64
#endif

// On some compilers, expands to __attribute__((aligned(CACHELINE_SIZE))).
// For compilers where this is not known to work, expands to nothing.
//
// No further guarantees are made here.  The result of applying the macro
// to variables and types is always implementation defined.
//
// WARNING: It is easy to use this attribute incorrectly, even to the point
// of causing bugs that are difficult to diagnose, crash, etc.  It does not
// guarantee that objects are aligned to a cache line.
//
// Recommendations:
//
// 1) Consult compiler documentation; this comment is not kept in sync as
//    toolchains evolve.
// 2) Verify your use has the intended effect. This often requires inspecting
//    the generated machine code.
// 3) Prefer applying this attribute to individual variables.  Avoid
//    applying it to types.  This tends to localize the effect.
#define CACHELINE_ALIGNED __attribute__((aligned(CACHELINE_SIZE)))

// Prevent the compiler from complaining about or optimizing away variables
// that appear unused
#undef ATTRIBUTE_UNUSED
#define ATTRIBUTE_UNUSED __attribute__ ((__unused__))

// For functions we want to force inline or not inline.
// Introduced in gcc 3.1.
#define ATTRIBUTE_ALWAYS_INLINE  __attribute__ ((always_inline))
#define HAVE_ATTRIBUTE_ALWAYS_INLINE 1
#define ATTRIBUTE_NOINLINE __attribute__ ((noinline))
#define HAVE_ATTRIBUTE_NOINLINE 1

// Prevent the compiler from optimizing away stack frames for functions which
// end in a call to another function.
#ifdef __clang__
#define ATTRIBUTE_NO_TAIL_CALL __attribute__((disable_tail_calls))
#else
#define ATTRIBUTE_NO_TAIL_CALL \
  __attribute__((optimize("no-optimize-sibling-calls")))
#endif
#define HAVE_ATTRIBUTE_NO_TAIL_CALL 1

// For weak functions
#undef ATTRIBUTE_WEAK
#define ATTRIBUTE_WEAK __attribute__ ((weak))
#define HAVE_ATTRIBUTE_WEAK 1

// Tell the compiler to use "initial-exec" mode for a thread-local variable.
// See http://people.redhat.com/drepper/tls.pdf for the gory details.
#define ATTRIBUTE_INITIAL_EXEC __attribute__ ((tls_model ("initial-exec")))

// Tell the compiler either that a particular function parameter
// should be a non-null pointer, or that all pointer arguments should
// be non-null.
//
// Note: As the GCC manual states, "[s]ince non-static C++ methods
// have an implicit 'this' argument, the arguments of such methods
// should be counted from two, not one."
//
// Args are indexed starting at 1.
// For non-static class member functions, the implicit "this" argument
// is arg 1, and the first explicit argument is arg 2.
// For static class member functions, there is no implicit "this", and
// the first explicit argument is arg 1.
//
//   /* arg_a cannot be null, but arg_b can */
//   void Function(void* arg_a, void* arg_b) ATTRIBUTE_NONNULL(1);
//
//   class C {
//     /* arg_a cannot be null, but arg_b can */
//     void Method(void* arg_a, void* arg_b) ATTRIBUTE_NONNULL(2);
//
//     /* arg_a cannot be null, but arg_b can */
//     static void StaticMethod(void* arg_a, void* arg_b) ATTRIBUTE_NONNULL(1);
//   };
//
// If no arguments are provided, then all pointer arguments should be non-null.
//
//  /* No pointer arguments may be null. */
//  void Function(void* arg_a, void* arg_b, int arg_c) ATTRIBUTE_NONNULL();
//
// NOTE: The GCC nonnull attribute actually accepts a list of arguments, but
// ATTRIBUTE_NONNULL does not.
#define ATTRIBUTE_NONNULL(arg_index) __attribute__((nonnull(arg_index)))

//
// Tell the compiler that a given function never returns
//
#define ATTRIBUTE_NORETURN __attribute__((noreturn))

// Tell AddressSanitizer (or other memory testing tools) to ignore a given
// function. Useful for cases when a function reads random locations on stack,
// calls _exit from a cloned subprocess, deliberately accesses buffer
// out of bounds or does other scary things with memory.
#ifdef ADDRESS_SANITIZER
#define ATTRIBUTE_NO_SANITIZE_ADDRESS \
    __attribute__((no_sanitize_address))
#else
#define ATTRIBUTE_NO_SANITIZE_ADDRESS
#endif

// Tell MemorySanitizer to relax the handling of a given function. All "Use of
// uninitialized value" warnings from such functions will be suppressed, and all
// values loaded from memory will be considered fully initialized.
// This is similar to the ADDRESS_SANITIZER attribute above, but deals with
// initializedness rather than addressability issues.
#ifdef MEMORY_SANITIZER
#define ATTRIBUTE_NO_SANITIZE_MEMORY \
    __attribute__((no_sanitize_memory))
#else
#define ATTRIBUTE_NO_SANITIZE_MEMORY
#endif

// Tell ThreadSanitizer to not instrument a given function.
// If you are adding this attribute, please cc dynamic-tools@ on the cl.
#ifdef THREAD_SANITIZER
#define ATTRIBUTE_NO_SANITIZE_THREAD __attribute__((no_sanitize_thread))
#else
#define ATTRIBUTE_NO_SANITIZE_THREAD
#endif

// Tell UndefinedSanitizer to ignore a given function. Useful for cases
// where certain behavior (eg. devision by zero) is being used intentionally.
#ifdef UNDEFINED_BEHAVIOR_SANITIZER
#define ATTRIBUTE_NO_SANITIZE_UNDEFINED \
  __attribute__((no_sanitize("undefined")))
#else
#define ATTRIBUTE_NO_SANITIZE_UNDEFINED
#endif

// Tell ControlFlowIntegrity sanitizer to not instrument a given function.
#ifdef CONTROL_FLOW_INTEGRITY
#define ATTRIBUTE_NO_SANITIZE_CFI __attribute__((no_sanitize("cfi")))
#else
#define ATTRIBUTE_NO_SANITIZE_CFI
#endif

#ifndef HAVE_ATTRIBUTE_SECTION  // may have been pre-set to 0, e.g. for Darwin
#define HAVE_ATTRIBUTE_SECTION 1
#endif

#if HAVE_ATTRIBUTE_SECTION  // define section support for the case of GCC

//
// Tell the compiler/linker to put a given function into a section and define
// "__start_ ## name" and "__stop_ ## name" symbols to bracket the section.
// This functionality is supported by GNU linker.
// Any function with ATTRIBUTE_SECTION must not be inlined, or it will
// be placed into whatever section its caller is placed into.
//
#ifndef ATTRIBUTE_SECTION
#define ATTRIBUTE_SECTION(name) \
  __attribute__ ((section (#name))) __attribute__ ((noinline))
#endif

//
// Weak section declaration to be used as a global declaration
// for ATTRIBUTE_SECTION_START|STOP(name) to compile and link
// even without functions with ATTRIBUTE_SECTION(name).
// DEFINE_ATTRIBUTE_SECTION should be in the exactly one file; it's
// a no-op on ELF but not on Mach-O.
//
#ifndef DECLARE_ATTRIBUTE_SECTION_VARS
#define DECLARE_ATTRIBUTE_SECTION_VARS(name) \
  extern char __start_##name[] ATTRIBUTE_WEAK; \
  extern char __stop_##name[] ATTRIBUTE_WEAK
#endif
#ifndef DEFINE_ATTRIBUTE_SECTION_VARS
#define INIT_ATTRIBUTE_SECTION_VARS(name)
#define DEFINE_ATTRIBUTE_SECTION_VARS(name)
#endif

//
// Return void* pointers to start/end of a section of code with
// functions having ATTRIBUTE_SECTION(name).
// Returns 0 if no such functions exits.
// One must DECLARE_ATTRIBUTE_SECTION_VARS(name) for this to compile and link.
//
#define ATTRIBUTE_SECTION_START(name) (reinterpret_cast<void*>(__start_##name))
#define ATTRIBUTE_SECTION_STOP(name) (reinterpret_cast<void*>(__stop_##name))

#endif  // HAVE_ATTRIBUTE_SECTION

// Support for aligning the stack on 32-bit x86.

#if defined(__i386__) && (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 2))
#define ATTRIBUTE_STACK_ALIGN_FOR_OLD_LIBC __attribute__((force_align_arg_pointer))
#define REQUIRE_STACK_ALIGN_TRAMPOLINE (0)
#elif defined(__i386__) || defined(__x86_64__)
#define REQUIRE_STACK_ALIGN_TRAMPOLINE (1)
#define ATTRIBUTE_STACK_ALIGN_FOR_OLD_LIBC
#else
#define REQUIRE_STACK_ALIGN_TRAMPOLINE (0)
#define ATTRIBUTE_STACK_ALIGN_FOR_OLD_LIBC
#endif

// Tell the compiler to warn about unused return values for functions declared
// with this macro. The macro must appear as the very first part of a function
// declaration or definition:
//
//   MUST_USE_RESULT Sprocket* AllocateSprocket();
//
// This placement has the broadest compatibility with GCC, Clang, and MSVC, with
// both defs and decls, and with GCC-style attributes, MSVC declspec, and C++11
// attributes. Note: past advice was to place the macro after the argument list.
#if defined(SWIG)
#define MUST_USE_RESULT
#elif __GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4)
#define MUST_USE_RESULT __attribute__ ((warn_unused_result))
#else
#define MUST_USE_RESULT
#endif

//
// Prevent the compiler from padding a structure to natural alignment
//
#if __GNUC__ && !defined(SWIG)
#define ATTRIBUTE_PACKED __attribute__((__packed__))
#else
#define ATTRIBUTE_PACKED
#endif

#ifdef __cplusplus
#if defined(__GNUC__) || defined(__llvm__)
// Defined behavior on some of the uarchs:
// PREFETCH_HINT_T0:
//   prefetch to all levels of the hierarchy (except on p4: prefetch to L2)
// PREFETCH_HINT_NTA:
//   p4: fetch to L2, but limit to 1 way (out of the 8 ways)
//   core: skip L2, go directly to L1
//   k8 rev E and later: skip L2, can go to either of the 2-ways in L1
enum PrefetchHint {
  PREFETCH_HINT_T0 = 3,  // More temporal locality
  PREFETCH_HINT_T1 = 2,
  PREFETCH_HINT_T2 = 1,  // Less temporal locality
  PREFETCH_HINT_NTA = 0  // No temporal locality
};
#else
// prefetch is a no-op for this target. Feel free to add more sections above.
#endif

// The default behavior of prefetch is to speculatively load for read only. This
// is safe for all currently supported platforms. However, prefetch for store
// may have problems depending on the target platform (x86, PPC, arm). Check
// with the platforms team (platforms-servers@) before introducing any changes
// to this function to identify potential impact on current and future servers.
extern inline void prefetch(const void *x, int hint) {
#if defined(__llvm__)
  // In the gcc version of prefetch(), hint is only a constant _after_ inlining
  // (assumed to have been successful).  llvm views things differently, and
  // checks constant-ness _before_ inlining.  This leads to compilation errors
  // with using the other version of this code with llvm.
  //
  // One way round this is to use a switch statement to explicitly match
  // prefetch hint enumerations, and invoke __builtin_prefetch for each valid
  // value.  llvm's optimization removes the switch and unused case statements
  // after inlining, so that this boils down in the end to the same as for gcc;
  // that is, a single inlined prefetchX instruction.
  //
  // Note that this version of prefetch() cannot verify constant-ness of hint.
  // If client code calls prefetch() with a variable value for hint, it will
  // receive the full expansion of the switch below, perhaps also not inlined.
  // This should however not be a problem in the general case of well behaved
  // caller code that uses the supplied prefetch hint enumerations.
  switch (hint) {
    case PREFETCH_HINT_T0:
      __builtin_prefetch(x, 0, PREFETCH_HINT_T0);
      break;
    case PREFETCH_HINT_T1:
      __builtin_prefetch(x, 0, PREFETCH_HINT_T1);
      break;
    case PREFETCH_HINT_T2:
      __builtin_prefetch(x, 0, PREFETCH_HINT_T2);
      break;
    case PREFETCH_HINT_NTA:
      __builtin_prefetch(x, 0, PREFETCH_HINT_NTA);
      break;
    default:
      __builtin_prefetch(x);
      break;
  }
#elif defined(__GNUC__)
 #if !defined(ARCH_PIII) || defined(__SSE__)
  if (__builtin_constant_p(hint)) {
    __builtin_prefetch(x, 0, hint);
  } else {
    // Defaults to PREFETCH_HINT_T0
    __builtin_prefetch(x);
  }
#else
  // We want a __builtin_prefetch, but we build with the default -march=i386
  // where __builtin_prefetch quietly turns into nothing.
  // Once we crank up to -march=pentium3 or higher the __SSE__
  // clause above will kick in with the builtin.
  // -- mec 2006-06-06
  if (hint == PREFETCH_HINT_NTA)
    __asm__ __volatile__("prefetchnta (%0)" : : "r"(x));
 #endif
#else
  // You get no effect.  Feel free to add more sections above.
#endif
}

// prefetch intrinsic (bring data to L1 without polluting L2 cache)
extern inline void prefetch(const void *x) {
  return prefetch(x, 0);
}
#endif  // ifdef __cplusplus

//
// GCC can be told that a certain branch is not likely to be taken (for
// instance, a CHECK failure), and use that information in static analysis.
// Giving it this information can help it optimize for the common case in
// the absence of better information (ie. -fprofile-arcs).
//
#if defined(__GNUC__)
#define PREDICT_FALSE(x) (__builtin_expect(x, 0))
#define PREDICT_TRUE(x) (__builtin_expect(!!(x), 1))
#else
#define PREDICT_FALSE(x) x
#define PREDICT_TRUE(x) x
#endif

//
// Tell GCC that a function is hot or cold. GCC can use this information to
// improve static analysis, i.e. a conditional branch to a cold function
// is likely to be not-taken.
// This annotation is used for function declarations, e.g.:
//   int foo() ATTRIBUTE_HOT;
//
#if __GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 3)
#define ATTRIBUTE_HOT __attribute__ ((hot))
#define ATTRIBUTE_COLD __attribute__ ((cold))
#else
#define ATTRIBUTE_HOT
#define ATTRIBUTE_COLD
#endif

#define FTELLO ftello
#define FSEEKO fseeko

#else   // not GCC

#define PRINTF_ATTRIBUTE(string_index, first_to_check)
#define SCANF_ATTRIBUTE(string_index, first_to_check)
#define CACHELINE_SIZE 64
#define CACHELINE_ALIGNED
#define ATTRIBUTE_UNUSED
#define ATTRIBUTE_ALWAYS_INLINE
#define ATTRIBUTE_NOINLINE
#define ATTRIBUTE_NO_TAIL_CALL
#define HAVE_ATTRIBUTE_NO_TAIL_CALL 0
#define ATTRIBUTE_HOT
#define ATTRIBUTE_COLD
#define ATTRIBUTE_WEAK
#define HAVE_ATTRIBUTE_WEAK 0
#define ATTRIBUTE_INITIAL_EXEC
#define ATTRIBUTE_NONNULL(...)
#if defined(_MSC_VER)
#define ATTRIBUTE_NORETURN __declspec(noreturn)
#else
#define ATTRIBUTE_NORETURN
#endif
#define ATTRIBUTE_NO_SANITIZE_ADDRESS
#define ATTRIBUTE_NO_SANITIZE_MEMORY
#if !defined(SWIG) || SWIG_VERSION >= 0x020000
#define HAVE_ATTRIBUTE_SECTION 0
#endif
#define ATTRIBUTE_PACKED
#define ATTRIBUTE_STACK_ALIGN_FOR_OLD_LIBC
#define REQUIRE_STACK_ALIGN_TRAMPOLINE (0)
#define MUST_USE_RESULT
#if defined(__cplusplus)
extern inline void prefetch(const void*) {}
#endif
#define PREDICT_FALSE(x) x
#define PREDICT_TRUE(x) x

// These should be redefined appropriately if better alternatives to
// ftell/fseek exist in the compiler
#define FTELLO ftell
#define FSEEKO fseek

#endif  // GCC

#if defined(__cplusplus) &&                                               \
    (((defined(__GNUC__) || defined(__APPLE__) || defined(OS_IOS) || \
       defined(__NVCC__)) &&                                              \
      !defined(SWIG)) ||                                                  \
     ((__GNUC__ >= 3 || defined(__clang__)) && defined(__ANDROID__)))

inline void *aligned_malloc(size_t size, int minimum_alignment) {
#if defined(__ANDROID__) || defined(OS_ANDROID) || defined(OS_CYGWIN)
  return memalign(minimum_alignment, size);
#else  // !__ANDROID__ && !OS_ANDROID && !OS_CYGWIN
  void *ptr = nullptr;
  // posix_memalign requires that the requested alignment be at least
  // sizeof(void*). In this case, fall back on malloc which should return memory
  // aligned to at least the size of a pointer.
  const int required_alignment = sizeof(void*);
  if (minimum_alignment < required_alignment)
    return malloc(size);
  if (posix_memalign(&ptr, minimum_alignment, size) != 0)
    return nullptr;
  else
    return ptr;
#endif
}

inline void aligned_free(void *aligned_memory) {
  free(aligned_memory);
}

#endif  // aligned_malloc

//
// Provides a char array with the exact same alignment as another type. The
// first parameter must be a complete type, the second parameter is how many
// of that type to provide space for.
//
//   ALIGNED_CHAR_ARRAY(struct stat, 16) storage_;
//
#if defined(__cplusplus)
#undef ALIGNED_CHAR_ARRAY
// Because MSVC and older GCCs require that the argument to their alignment
// construct to be a literal constant integer, we use a template instantiated
// at all the possible powers of two.
#ifndef SWIG
template<int alignment, int size> struct AlignType { };
template<int size> struct AlignType<0, size> { typedef char result[size]; };
#if defined(_MSC_VER)
#define BASE_PORT_H_ALIGN_ATTRIBUTE(X) __declspec(align(X))
#define BASE_PORT_H_ALIGN_OF(T) __alignof(T)
#elif defined(__GNUC__) || defined(__INTEL_COMPILER)
#define BASE_PORT_H_ALIGN_ATTRIBUTE(X) __attribute__((aligned(X)))
#define BASE_PORT_H_ALIGN_OF(T) __alignof__(T)
#endif

#if defined(BASE_PORT_H_ALIGN_ATTRIBUTE)

#define BASE_PORT_H_ALIGNTYPE_TEMPLATE(X) \
  template<int size> struct AlignType<X, size> { \
    typedef BASE_PORT_H_ALIGN_ATTRIBUTE(X) char result[size]; \
  }

BASE_PORT_H_ALIGNTYPE_TEMPLATE(1);
BASE_PORT_H_ALIGNTYPE_TEMPLATE(2);
BASE_PORT_H_ALIGNTYPE_TEMPLATE(4);
BASE_PORT_H_ALIGNTYPE_TEMPLATE(8);
BASE_PORT_H_ALIGNTYPE_TEMPLATE(16);
BASE_PORT_H_ALIGNTYPE_TEMPLATE(32);
BASE_PORT_H_ALIGNTYPE_TEMPLATE(64);
BASE_PORT_H_ALIGNTYPE_TEMPLATE(128);
BASE_PORT_H_ALIGNTYPE_TEMPLATE(256);
BASE_PORT_H_ALIGNTYPE_TEMPLATE(512);
BASE_PORT_H_ALIGNTYPE_TEMPLATE(1024);
BASE_PORT_H_ALIGNTYPE_TEMPLATE(2048);
BASE_PORT_H_ALIGNTYPE_TEMPLATE(4096);
BASE_PORT_H_ALIGNTYPE_TEMPLATE(8192);
// Any larger and MSVC++ will complain.

#define ALIGNED_CHAR_ARRAY(T, Size) \
  typename AlignType<BASE_PORT_H_ALIGN_OF(T), sizeof(T) * Size>::result

#undef BASE_PORT_H_ALIGNTYPE_TEMPLATE
#undef BASE_PORT_H_ALIGN_ATTRIBUTE

#else  // defined(BASE_PORT_H_ALIGN_ATTRIBUTE)
#define ALIGNED_CHAR_ARRAY \
  you_must_define_ALIGNED_CHAR_ARRAY_for_your_compiler_in_base_port_h
#endif  // defined(BASE_PORT_H_ALIGN_ATTRIBUTE)

#else  // !SWIG

// SWIG can't represent alignment and doesn't care about alignment on data
// members (it works fine without it).
template<typename Size>
struct AlignType { typedef char result[Size]; };
#define ALIGNED_CHAR_ARRAY(T, Size) AlignType<Size * sizeof(T)>::result

// Enough to parse with SWIG, will never be used by running code.
#define BASE_PORT_H_ALIGN_OF(Type) 16

#endif // !SWIG
#else  // __cplusplus
#define ALIGNED_CHAR_ARRAY ALIGNED_CHAR_ARRAY_is_not_available_without_Cplusplus
#endif // __cplusplus

#if !HAVE_ATTRIBUTE_SECTION  // provide dummy definitions

#define ATTRIBUTE_SECTION(name)
#define INIT_ATTRIBUTE_SECTION_VARS(name)
#define DEFINE_ATTRIBUTE_SECTION_VARS(name)
#define DECLARE_ATTRIBUTE_SECTION_VARS(name)
#define ATTRIBUTE_SECTION_START(name) (reinterpret_cast<void*>(0))
#define ATTRIBUTE_SECTION_STOP(name) (reinterpret_cast<void*>(0))

#endif  // !HAVE_ATTRIBUTE_SECTION


#ifdef _MSC_VER     /* if Visual C++ */

// This compiler flag can be easily overlooked on MSVC.
// _CHAR_UNSIGNED gets set with the /J flag.
#ifndef _CHAR_UNSIGNED
#error chars must be unsigned!  Use the /J flag on the compiler command line.
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

// Don't warn about unused local variables.
//
// extension to silence particular instances of this warning.  There's no way
// to define ATTRIBUTE_UNUSED to quiet particular instances of this warning in
// VC++, so we disable it globally.  Currently, there aren't many false
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

#include <winsock2.h>
#include <cassert>
#include <windows.h>
#include <process.h>  // _getpid()
#undef ERROR

#include <cmath>  // for HUGE_VAL

#ifndef HUGE_VALF
#define HUGE_VALF (static_cast<float>(HUGE_VAL))
#endif

#ifdef __cplusplus
// Define a minimal set of things typically available in the global
// namespace in Google code.  ::string is handled elsewhere, and uniformly
// for all targets.
#include <functional>
using std::hash;
#endif  // __cplusplus

// VC++ doesn't understand "uint"
#ifndef HAVE_UINT
#define HAVE_UINT 1
typedef unsigned int uint;
#endif

// VC++ doesn't understand "ssize_t"
// <windows.h> from above includes <BaseTsd.h> and <BaseTsd.h> defines SSIZE_T
#ifndef HAVE_SSIZET
#define HAVE_SSIZET 1
typedef SSIZE_T ssize_t;
#endif

#define strtoq   _strtoi64
#define strtouq  _strtoui64
#define strtoll  _strtoi64
#define strtoull _strtoui64
#define atoll    _atoi64


// You say tomato, I say atotom
#define PATH_MAX MAX_PATH

// MSVC requires to know how code will be linked in order to compile it.
// This information can be provided via a .def file or __declspec() annotations.
// The following macro can be set on the compiler command line by those wishing
// to use __declspec.
#ifndef BASE_PORT_MSVC_DLL_MACRO
#define BASE_PORT_MSVC_DLL_MACRO
#endif

// Wrap Microsoft _snprintf/_vsnprintf calls so they nul-terminate on buffer
// overflow.
#define vsnprintf base_port_MSVC_vsnprintf
BASE_PORT_MSVC_DLL_MACRO
    int base_port_MSVC_vsnprintf(char *str, size_t size,
                                 const char *format, va_list ap);
#define snprintf base_port_MSVC_snprintf
BASE_PORT_MSVC_DLL_MACRO
    int base_port_MSVC_snprintf(char *str, size_t size, const char *fmt, ...);

// You say tomato, I say _tomato
#define strcasecmp _stricmp
#define strncasecmp _strnicmp
#define strdup _strdup
#define tempnam _tempnam
#define chdir  _chdir
#define getcwd _getcwd
#define putenv  _putenv
#if _MSC_VER >= 1900  // Only needed for VS2015+
#define getpid _getpid
#define timezone _timezone
#define tzname _tzname
#endif

// You say tomato, I say toma
inline int random() { return rand(); }
inline void srandom(unsigned int seed) { srand(seed); }

// You say juxtapose, I say transpose
#define bcopy(s, d, n) memcpy(d, s, n)

inline void *aligned_malloc(size_t size, int minimum_alignment) {
  return _aligned_malloc(size, minimum_alignment);
}

inline void aligned_free(void *aligned_memory) {
  _aligned_free(aligned_memory);
}

// ----- BEGIN VC++ STUBS & FAKE DEFINITIONS ---------------------------------

typedef void (*sig_t)(int);

// This actually belongs in errno.h but there's a name conflict in errno
// on WinNT. They (and a ton more) are also found in Winsock2.h, but
// if'd out under NT. We need this subset at minimum.
#define EXFULL      ENOMEM  // not really that great a translation...

//
// Really from <string.h>
//

inline void bzero(void *s, int n) {
  memset(s, 0, n);
}

// From glob.h
#define __ptr_t   void *

// Defined all over the place.
typedef int pid_t;

// From stat.h
typedef unsigned int mode_t;

// u_int16_t, int16_t don't exist in MSVC
typedef unsigned short u_int16_t;
typedef short int16_t;

// ----- END VC++ STUBS & FAKE DEFINITIONS ----------------------------------

#endif  // _MSC_VER

#ifdef __cplusplus
#ifdef STL_MSVC  // not always the same as _MSC_VER
#include "s2/third_party/absl/base/internal/port_hash.inc"
#else
struct PortableHashBase { };
#endif
#endif  // __cplusplus

#if defined(OS_WINDOWS) || defined(__APPLE__) || defined(OS_IOS)
// gethostbyname() *is* thread-safe for Windows native threads. It is also
// safe on Mac OS X and iOS, where it uses thread-local storage, even though the
// manpages claim otherwise. For details, see
// http://lists.apple.com/archives/Darwin-dev/2006/May/msg00008.html
#else
// gethostbyname() is not thread-safe.  So disallow its use.  People
// should either use the HostLookup::Lookup*() methods, or gethostbyname_r()
#define gethostbyname gethostbyname_is_not_thread_safe_DO_NOT_USE
#endif

// Define the namespace for pre-C++11 functors for hash_map and hash_set.
// This is not the namespace for C++11 functors (that namespace is "std").
//
// We used to require that the build tool or Makefile provide this definition.
// Now we usually get it from testing target macros. If the testing target
// macros are different from an external definition, you will get a build
// error.
//
// TODO(user): always get HASH_NAMESPACE from testing target macros.

#if defined(__GNUC__) && defined(GOOGLE_GLIBCXX_VERSION)
// Crosstool v17 or later.
#define HASH_NAMESPACE __gnu_cxx
#elif defined(_MSC_VER)
// MSVC.
// http://msdn.microsoft.com/en-us/library/6x7w9f6z(v=vs.100).aspx
#define HASH_NAMESPACE stdext
#elif defined(__APPLE__)
// Xcode.
#define HASH_NAMESPACE __gnu_cxx
#elif defined(__GNUC__)
// Some other version of gcc.
#define HASH_NAMESPACE __gnu_cxx
#else
// HASH_NAMESPACE defined externally.
// TODO(user): make this an error. Do not use external value of HASH_NAMESPACE.
#endif

#ifndef HASH_NAMESPACE
// TODO(user): try to delete this.
// I think gcc 2.95.3 was the last toolchain to use this.
#define HASH_NAMESPACE_DECLARATION_START
#define HASH_NAMESPACE_DECLARATION_END
#else
#define HASH_NAMESPACE_DECLARATION_START  namespace HASH_NAMESPACE {
#define HASH_NAMESPACE_DECLARATION_END    }
#endif

// Our STL-like classes use __STD.
#if defined(__GNUC__) || defined(__APPLE__) || defined(OS_IOS) || \
    defined(_MSC_VER)
#define __STD std
#endif

#if defined __GNUC__
#define STREAM_SET(s, bit) (s).setstate(std::ios_base::bit)
#define STREAM_SETF(s, flag) (s).setf(std::ios_base::flag)
#else
#define STREAM_SET(s, bit) (s).set(std::ios::bit)
#define STREAM_SETF(s, flag) (s).setf(std::ios::flag)
#endif

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

// The unaligned API is C++ only.  The declarations use C++ features
// (namespaces, inline) which are absent or incompatible in C.
#if defined(__cplusplus)

#if defined(ADDRESS_SANITIZER) || defined(THREAD_SANITIZER) ||\
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

#elif defined(ARCH_PIII) || defined(ARCH_ATHLON) || \
     defined(ARCH_K8) || defined(ARCH_PPC)

// x86 and x86-64 can perform unaligned loads/stores directly;
// modern PowerPC hardware can also do unaligned integer loads and stores;
// but note: the FPU still sends unaligned loads and stores to a trap handler!

#define UNALIGNED_LOAD16(_p) (*reinterpret_cast<const uint16 *>(_p))
#define UNALIGNED_LOAD32(_p) (*reinterpret_cast<const uint32 *>(_p))
#define UNALIGNED_LOAD64(_p) (*reinterpret_cast<const uint64 *>(_p))

#define UNALIGNED_STORE16(_p, _val) (*reinterpret_cast<uint16 *>(_p) = (_val))
#define UNALIGNED_STORE32(_p, _val) (*reinterpret_cast<uint32 *>(_p) = (_val))
#define UNALIGNED_STORE64(_p, _val) (*reinterpret_cast<uint64 *>(_p) = (_val))

#elif defined(__arm__) && \
      !defined(__ARM_ARCH_5__) && \
      !defined(__ARM_ARCH_5T__) && \
      !defined(__ARM_ARCH_5TE__) && \
      !defined(__ARM_ARCH_5TEJ__) && \
      !defined(__ARM_ARCH_6__) && \
      !defined(__ARM_ARCH_6J__) && \
      !defined(__ARM_ARCH_6K__) && \
      !defined(__ARM_ARCH_6Z__) && \
      !defined(__ARM_ARCH_6ZK__) && \
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
} ATTRIBUTE_PACKED;

struct Unaligned32Struct {
  uint32 value;
  uint8 dummy;  // To make the size non-power-of-two.
} ATTRIBUTE_PACKED;

}  // namespace internal
}  // namespace base

#define UNALIGNED_LOAD16(_p) \
    ((reinterpret_cast<const ::base::internal::Unaligned16Struct *>(_p))->value)
#define UNALIGNED_LOAD32(_p) \
    ((reinterpret_cast<const ::base::internal::Unaligned32Struct *>(_p))->value)

#define UNALIGNED_STORE16(_p, _val) \
    ((reinterpret_cast< ::base::internal::Unaligned16Struct *>(_p))->value = \
         (_val))
#define UNALIGNED_STORE32(_p, _val) \
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

inline void UNALIGNED_STORE64(void *p, uint64 v) {
  memcpy(p, &v, sizeof v);
}

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

inline void UNALIGNED_STORE16(void *p, uint16 v) {
  memcpy(p, &v, sizeof v);
}

inline void UNALIGNED_STORE32(void *p, uint32 v) {
  memcpy(p, &v, sizeof v);
}

inline void UNALIGNED_STORE64(void *p, uint64 v) {
  memcpy(p, &v, sizeof v);
}

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

// printf macros for size_t, in the style of inttypes.h
#if defined(_LP64) || defined(OS_IOS)
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

#define GPRIuPTHREAD "lu"
#define GPRIxPTHREAD "lx"
#ifdef OS_CYGWIN
#define PRINTABLE_PTHREAD(pthreadt) reinterpret_cast<uintptr_t>(pthreadt)
#else
#define PRINTABLE_PTHREAD(pthreadt) pthreadt
#endif

#define SIZEOF_MEMBER(t, f)   sizeof(((t*) 4096)->f)

#define OFFSETOF_MEMBER(t, f)         \
  (reinterpret_cast<char*>(           \
     &reinterpret_cast<t*>(16)->f) -  \
   reinterpret_cast<char*>(16))

#ifdef PTHREADS_REDHAT_WIN32
#include <iosfwd>     // NOLINT(build/include)
#include <pthread.h>  // NOLINT(build/include)
// pthread_t is not a simple integer or pointer on Win32
std::ostream& operator << (std::ostream& out, const pthread_t& thread_id);
#endif

// GXX_EXPERIMENTAL_CXX0X is defined by gcc and clang up to at least
// gcc-4.7 and clang-3.1 (2011-12-13).  __cplusplus was defined to 1
// in gcc before 4.7 (Crosstool 16) and clang before 3.1, but is
// defined according to the language version in effect thereafter.
// Microsoft Visual Studio 14 (2015) sets __cplusplus==199711 despite
// reasonably good C++11 support, so we set LANG_CXX for it and
// newer versions (_MSC_VER >= 1900).
#if (defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L || \
     (defined(_MSC_VER) && _MSC_VER >= 1900))
// Define this to 1 if the code is compiled in C++11 mode; leave it
// undefined otherwise.  Do NOT define it to 0 -- that causes
// '#ifdef LANG_CXX11' to behave differently from '#if LANG_CXX11'.
#define LANG_CXX11 1
#endif

// On some platforms, a "function pointer" points to a function descriptor
// rather than directly to the function itself.  Use FUNC_PTR_TO_CHAR_PTR(func)
// to get a char-pointer to the first instruction of the function func.
#if (defined(__powerpc__) && !(_CALL_ELF > 1)) || \
    defined(__ia64)
// use opd section for function descriptors on these platforms, the function
// address is the first word of the descriptor
enum { kPlatformUsesOPDSections = 1 };
#define FUNC_PTR_TO_CHAR_PTR(func) (reinterpret_cast<char**>(func)[0])
#else
enum { kPlatformUsesOPDSections = 0 };
#define FUNC_PTR_TO_CHAR_PTR(func)  (reinterpret_cast<char *>(func))
#endif

// Private implementation detail: __has_extension is useful to implement
// static_assert, and defining it for all toolchains avoids an extra level of
// nesting of #if/#ifdef/#ifndef.
#ifndef __has_extension
#define __has_extension(x) 0  // MSVC 10's preprocessor can't handle 'false'.
#endif

#ifdef __cplusplus
// CompileAssert<T> is deprecated.  Use static_assert instead.
template <bool>
struct CompileAssert {
};
#endif  // __cplusplus

// Some platforms have a ::string class that is different from ::std::string
// (although the interface is the same, of course).  On other platforms,
// ::string is the same as ::std::string.
#if defined(__cplusplus) && !defined(SWIG)
#include <string>
#ifndef HAS_GLOBAL_STRING
using std::basic_string;
using std::string;
// TODO(user): using std::wstring?
#endif  // HAS_GLOBAL_STRING
#endif  // SWIG, __cplusplus

#ifdef __cplusplus
namespace base {
// We support C++14's sized deallocation for all C++ builds,
// though for other toolchains, we fall back to using delete.
inline void sized_delete(void* ptr, size_t size) {
#ifdef GOOGLE_HAVE_SIZED_DELETE
  ::operator delete(ptr, size);
#else
  (void)size;
  ::operator delete(ptr);
#endif  // GOOGLE_HAVE_SIZED_DELETE
}
}  // namespace base
#endif  // __cplusplus

// This sanity check can be removed when all references to
// LANG_CXX11 is removed from the code base.
#if defined(__cplusplus) && !defined(LANG_CXX11) && !defined(SWIG)
#error "LANG_CXX11 is required."
#endif

// A variable declaration annotated with the ABSL_CONST_INIT attribute will
// not compile (on supported platforms) unless the variable has a constant
// initializer. This is useful for variables with static and thread storage
// duration, because it guarantees that they will not suffer from the so-called
// "static init order fiasco".
//
// Sample usage:
//
// ABSL_CONST_INIT static MyType my_var = MakeMyType(...);
//
// Note that this attribute is redundant if the variable is declared constexpr.
#if defined(__cplusplus) && defined(__has_cpp_attribute)
// NOTE: requiring __cplusplus above should not be necessary, but
// works around https://bugs.llvm.org/show_bug.cgi?id=23435.
#if __has_cpp_attribute(clang::require_constant_initialization)
// NOLINTNEXTLINE(whitespace/braces) (b/36288871)
#define ABSL_CONST_INIT [[clang::require_constant_initialization]]
#else
#define ABSL_CONST_INIT
#endif  // __has_cpp_attribute(clang::require_constant_initialization)
#else
#define ABSL_CONST_INIT
#endif  // defined(__has_cpp_attribute)

#endif  // S2_THIRD_PARTY_ABSL_BASE_PORT_H_
