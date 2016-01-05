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
//
// These are weird things we need to do to get this compiling on
// random systems (and on SWIG).
//

#ifndef S2GEOMETRY_BASE_PORT_H_
#define S2GEOMETRY_BASE_PORT_H_

#include <climits>         // So we can set the bounds of our types
#include <cstring>         // for memcpy()
#include <cstdlib>         // for free()

#if defined(OS_CYGWIN)
#error "Cygwin is not supported."
#endif

#if defined(__CYGWIN__)
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

#if defined(__APPLE__) || defined(OS_IOS)
#include <unistd.h>         // for getpagesize() on mac
#elif defined(OS_CYGWIN) || defined(__ANDROID__)
#include <malloc.h>         // for memalign()
#elif defined(_MSC_VER)
#include <cstdio>          // declare snprintf/vsnprintf before overriding
#endif

#include "base/integral_types.h"

// We support gcc 4.6 and later.
#if defined(__GNUC__) && !defined(__clang__)
#if __GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 6)
#error "This package requires gcc 4.6 or higher"
#endif
#endif

// We support MSVC++ 12.0 and later.
#if defined(_MSC_VER) && _MSC_VER < 1800
#error "This package requires _MSC_VER of 1800 or higher"
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

#if defined(__cplusplus)
#include <cstddef>              // For _GLIBCXX macros
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
  return static_cast<uint16>(((x & 0xFF) << 8) | ((x & 0xFF00) >> 8));
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

// Mach-O supports sections (albeit with small names), but doesn't have
// vars at the beginning and end.  Instead you should call the function
// getsectdata("__DATA", name, &size).
#define HAVE_ATTRIBUTE_SECTION 1

// Any function with ATTRIBUTE_SECTION must not be inlined, or it will
// be placed into whatever section its caller is placed into.
#define ATTRIBUTE_SECTION(name) \
  __attribute__ ((section ("__DATA, " #name))) __attribute__ ((noinline))

#define ENUM_DYLD_BOOL  // so that we don't pollute the global namespace
extern "C" {
  #include <mach-o/getsect.h>
  #include <mach-o/dyld.h>
}
class AssignAttributeStartEnd {
 public:
  AssignAttributeStartEnd(const char* name, char** pstart, char** pend) {
    // Find out what dynamic library name is defined in
    for (int i = _dyld_image_count() - 1; i >= 0; --i) {
      const mach_header* hdr = _dyld_get_image_header(i);
      uint32_t len;
      *pstart = getsectdatafromheader(hdr, "__DATA", name, &len);
      if (*pstart) {   // nullptr if not defined in this dynamic library
        *pstart += _dyld_get_image_vmaddr_slide(i);   // correct for reloc
        *pend = *pstart + len;
        return;
      }
    }
    // If we get here, not defined in a dll at all.  See if defined statically.
    unsigned long len;    // don't ask me why this type isn't uint32_t too...
    *pstart = getsectdata("__DATA", name, &len);
    *pend = *pstart + len;
  }
};

// 1) DEFINE_ATTRIBUTE_SECTION_VARS: must be called once per unique
//    name.  You want to make sure this is executed before any
//    DECLARE_ATTRIBUTE_SECTION_VARS; the easiest way is to put them
//    in the same .cc file.  Put this call at the global level.
// 2) INIT_ATTRIBUTE_SECTION_VARS: you can scatter calls to this in
//    multiple places to help ensure execution before any
//    DECLARE_ATTRIBUTE_SECTION_VARS.  You must have at least one
//    DEFINE, but you can have many INITs.  Put each in its own scope.
// 3) DECLARE_ATTRIBUTE_SECTION_VARS: must be called before using
//    ATTRIBUTE_SECTION_START or ATTRIBUTE_SECTION_STOP on a name.
//    Put this call at the global level.
#define DECLARE_ATTRIBUTE_SECTION_VARS(name) \
  extern char* __start_##name; \
  extern char* __stop_##name;

#define INIT_ATTRIBUTE_SECTION_VARS(name)               \
  DECLARE_ATTRIBUTE_SECTION_VARS(name);                 \
  static const AssignAttributeStartEnd __assign_##name( \
    #name, &__start_##name, &__stop_##name)

#define DEFINE_ATTRIBUTE_SECTION_VARS(name)             \
  char* __start_##name, *__stop_##name;                 \
  INIT_ATTRIBUTE_SECTION_VARS(name)

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
// We would need to read special regiter ctr_el0 to find out L1 dcache size.
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

#define CACHELINE_ALIGNED __attribute__((aligned(CACHELINE_SIZE)))

//
// Prevent the compiler from complaining about or optimizing away variables
// that appear unused
#undef ATTRIBUTE_UNUSED
#define ATTRIBUTE_UNUSED __attribute__ ((unused))

//
// For functions we want to force inline or not inline.
// Introduced in gcc 3.1.
#define ATTRIBUTE_ALWAYS_INLINE  __attribute__ ((always_inline))
#define HAVE_ATTRIBUTE_ALWAYS_INLINE 1
#define ATTRIBUTE_NOINLINE __attribute__ ((noinline))
#define HAVE_ATTRIBUTE_NOINLINE 1

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
//   /* arg_a cannot be nullptr, but arg_b can */
//   void Function(void* arg_a, void* arg_b) ATTRIBUTE_NONNULL(1);
//
//   class C {
//     /* arg_a cannot be nullptr, but arg_b can */
//     void Method(void* arg_a, void* arg_b) ATTRIBUTE_NONNULL(2);
//
//     /* arg_a cannot be nullptr, but arg_b can */
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


//
// Tell the compiler to warn about unused return values for functions declared
// with this macro.  The macro should be used on function declarations
// following the argument list:
//
//   Sprocket* AllocateSprocket() MUST_USE_RESULT;
//
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

#ifdef __cplusplus
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
#define ATTRIBUTE_HOT
#define ATTRIBUTE_COLD
#define ATTRIBUTE_WEAK
#define HAVE_ATTRIBUTE_WEAK 0
#define ATTRIBUTE_INITIAL_EXEC
#define ATTRIBUTE_NONNULL(arg_index)
#define ATTRIBUTE_NORETURN
#define ATTRIBUTE_NO_SANITIZE_ADDRESS
#define ATTRIBUTE_NO_SANITIZE_MEMORY
#define HAVE_ATTRIBUTE_SECTION 0
#define ATTRIBUTE_PACKED
#define ATTRIBUTE_STACK_ALIGN_FOR_OLD_LIBC
#define REQUIRE_STACK_ALIGN_TRAMPOLINE (0)
#define MUST_USE_RESULT
extern inline void prefetch(const void*) {}
#define PREDICT_FALSE(x) x
#define PREDICT_TRUE(x) x

// These should be redefined appropriately if better alternatives to
// ftell/fseek exist in the compiler
#define FTELLO ftell
#define FSEEKO fseek

#endif  // GCC

#if ((defined(__GNUC__) || defined(__APPLE__) || defined(OS_IOS) ||  \
      defined(__NVCC__)) && !defined(SWIG)) ||                            \
    ((__GNUC__ >= 3 || defined(__clang__)) && defined(__ANDROID__))

#if !defined(__cplusplus) && !defined(__APPLE__) && !defined(OS_IOS) && \
    !defined(OS_CYGWIN)
// stdlib.h only declares this in C++, not in C, so we declare it here.
// Also make sure to avoid declaring it on platforms which don't support it.
extern int posix_memalign(void **memptr, size_t alignment, size_t size);
#endif

inline void *aligned_malloc(size_t size, int minimum_alignment) {
#if defined(__APPLE__)
  // mac lacks memalign(), posix_memalign(), however, according to
  // http://stackoverflow.com/questions/196329/osx-lacks-memalign
  // mac allocs are already 16-byte aligned.
  if (minimum_alignment <= 16)
    return malloc(size);
  // next, try to return page-aligned memory. perhaps overkill
  if (minimum_alignment <= getpagesize())
    return valloc(size);
  // give up
  return nullptr;
#elif defined(__ANDROID__) || defined(OS_ANDROID) || defined(OS_CYGWIN)
  return memalign(minimum_alignment, size);
#else  // !__ANDROID__ && !OS_ANDROID && !__APPLE__ && !OS_CYGWIN
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

#endif
// #if ((defined(__GNUC__) || defined(__APPLE__) || defined(OS_IOS) ||
// defined(__NVCC__)) && !defined(SWIG)) ||
// ((__GNUC__ >= 3 || defined(__clang__)) && defined(__ANDROID__))

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
#elif defined(__GNUC__)
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

#endif // !SWIG
#else  // __cpluscplus
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

// MSVC is a little hyper-active in its warnings
// Signed vs. unsigned comparison is ok.
#pragma warning(disable : 4018 )
// We know casting from a long to a char may lose data
#pragma warning(disable : 4244 )
// Don't need performance warnings about converting ints to bools
#pragma warning(disable : 4800 )
// Integral constant overflow is apparently ok too
// for example:
//  short k;  int n;
//  k = k + n;
#pragma warning(disable : 4307 )
// It's ok to use this* in constructor
// Example:
//  class C {
//   Container cont_;
//   C() : cont_(this) { ...
#pragma warning(disable : 4355 )
// Truncating from double to float is ok
#pragma warning(disable : 4305 )

#include <winsock2.h>
#include <cassert>
#include <windows.h>
#undef ERROR

#include <cfloat>  // for nextafter functionality on windows
#include <cmath>  // for HUGE_VAL

#ifndef HUGE_VALF
#define HUGE_VALF (static_cast<float>(HUGE_VAL))
#endif

namespace std {}  // Avoid error if we didn't see std.
using namespace std;

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
#define nextafter _nextafter
#define strdup _strdup
#define tempnam _tempnam
#define chdir  _chdir
#define getcwd _getcwd
#define putenv  _putenv


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

// See http://en.wikipedia.org/wiki/IEEE_754 for details of
// floating point format.

inline int fpclassify_double(double x) {
  const int float_point_class =_fpclass(x);
  int c99_class;
  switch  (float_point_class) {
  case _FPCLASS_SNAN:  // Signaling NaN
  case _FPCLASS_QNAN:  // Quiet NaN
    c99_class = FP_NAN;
    break;
  case _FPCLASS_NZ:  // Negative zero ( -0)
  case _FPCLASS_PZ:  // Positive 0 (+0)
    c99_class = FP_ZERO;
    break;
  case _FPCLASS_NINF:  // Negative infinity ( -INF)
  case _FPCLASS_PINF:  // Positive infinity (+INF)
    c99_class = FP_INFINITE;
    break;
  case _FPCLASS_ND:  // Negative denormalized
  case _FPCLASS_PD:  // Positive denormalized
    c99_class = FP_SUBNORMAL;
    break;
  case _FPCLASS_NN:  // Negative normalized non-zero
  case _FPCLASS_PN:  // Positive normalized non-zero
    c99_class = FP_NORMAL;
    break;
  default:
    c99_class = FP_NAN;  // Should never happen
    break;
  }
  return  c99_class;
}

// This function handle the special subnormal case for float; it will
// become a normal number while casting to double.
// bit_cast is avoided to simplify dependency and to create a code that is
// easy to deploy in C code
inline int fpclassify_float(float x) {
  uint32 bitwise_representation;
  memcpy(&bitwise_representation, &x, 4);
  if ((bitwise_representation & 0x7f800000) == 0 &&
      (bitwise_representation & 0x007fffff) != 0)
    return FP_SUBNORMAL;
  return fpclassify_double(x);
}
//
// This define takes care of the denormalized float; the casting to
// double make it a normal number
#define fpclassify(x) ((sizeof(x) == sizeof(float)) ? fpclassify_float(x) : fpclassify_double(x))

#define isnan _isnan

inline int isinf(double x) {
  const int float_point_class =_fpclass(x);
  if (float_point_class == _FPCLASS_PINF) return 1;
  if (float_point_class == _FPCLASS_NINF) return -1;
  return 0;
}

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

#ifdef STL_MSVC  // not always the same as _MSC_VER
#include "base/port_hash.h"
#else
struct PortableHashBase { };
#endif

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
#elif defined(__GNUC__) && (defined(STLPORT) || defined(USE_STD_HASH))
// A version of gcc with stlport.
#define HASH_NAMESPACE std
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
#define STREAM_SET(s, bit) (s).setstate(ios_base::bit)
#define STREAM_SETF(s, flag) (s).setf(ios_base::flag)
#else
#define STREAM_SET(s, bit) (s).set(ios::bit)
#define STREAM_SETF(s, flag) (s).setf(ios::flag)
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

#ifdef __cplusplus
extern "C" {
#endif  // __cplusplus
uint16_t __sanitizer_unaligned_load16(const void *p);
uint32_t __sanitizer_unaligned_load32(const void *p);
uint64_t __sanitizer_unaligned_load64(const void *p);
void __sanitizer_unaligned_store16(void *p, uint16_t v);
void __sanitizer_unaligned_store32(void *p, uint32_t v);
void __sanitizer_unaligned_store64(void *p, uint64_t v);
#ifdef __cplusplus
}  // extern "C"
#endif  // __cplusplus

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

// NOTE(user): These are only exported to C++ because the macros they depend on
// use C++-only syntax. This #ifdef can be removed if/when the macros are fixed.

#if defined(__cplusplus)

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

#endif  // defined(__cpluscplus)

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
// defined according to the language version in effect thereafter.  I
// believe MSVC will also define __cplusplus according to the language
// version, but haven't checked that. Stlport is used by many Android projects
// and does not have full C++11 STL support.
#if (defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L) && \
    !defined(STLPORT)
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
// We support C++11's static_assert(expression, message) for all C++
// builds, though for some pre-C++11 toolchains we fall back to using
// GG_PRIVATE_STATIC_ASSERT, which has two limitations: (1) the
// expression argument will need to be parenthesized if it would
// otherwise contain commas outside of parentheses, and (2) the
// message is ignored (though the compiler error will likely mention
// "static_assert_failed" and point to the line with the failing assertion).

// Something else (perhaps libc++) may have provided its own definition of
// static_assert.
#ifndef static_assert
#if LANG_CXX11 || __has_extension(cxx_static_assert) || defined(_MSC_VER)
// There's a native implementation of static_assert, no need to define our own.
#elif __has_extension(c_static_assert)
// C11's _Static_assert is available, and makes a great static_assert.
#define static_assert _Static_assert
#else
// Fall back on our home-grown implementation, with its limitations.
#define static_assert GG_PRIVATE_STATIC_ASSERT
#endif
#endif

// CompileAssert is an implementation detail of COMPILE_ASSERT and
// GG_PRIVATE_STATIC_ASSERT.
template <bool>
struct CompileAssert {
};

// GG_PRIVATE_STATIC_ASSERT: A poor man's static_assert.  This doesn't handle
// condition expressions that contain unparenthesized top-level commas;
// write GG_PRIVATE_STATIC_ASSERT((expr), "comment") when needed.
#define GG_PRIVATE_CAT_IMMEDIATE(a, b) a ## b
#define GG_PRIVATE_CAT(a, b) GG_PRIVATE_CAT_IMMEDIATE(a, b)
#define GG_PRIVATE_STATIC_ASSERT(expr, ignored) \
  typedef CompileAssert<(static_cast<bool>(expr))> \
  GG_PRIVATE_CAT(static_assert_failed_at_line, __LINE__)[bool(expr) ? 1 : -1] \
  ATTRIBUTE_UNUSED

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

#endif  // S2GEOMETRY_BASE_PORT_H_
