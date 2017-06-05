// Copyright 2017 Google Inc. All Rights Reserved.
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

// Defines preprocessor macros describing the presence of "features" available.
// This facilitates writing portable code by parameterizing the compilation
// based on the presence or lack of a feature.
//
// We define a feature as some interface we wish to program to: for example,
// some library function or system call.
//
// For example, suppose a programmer wants to write a program that uses the
// 'mmap' system call. Then one might write:
//
// #include "third_party/absl/base/config.h"
//
// #ifdef GOOGLE_HAVE_MMAP
// #include "sys/mman.h"
// #endif  // GOOGLE_HAVE_MMAP
//
// ...
// #ifdef GOOGLE_HAVE_MMAP
// void *ptr = mmap(...);
// ...
// #endif  // GOOGLE_HAVE_MMAP
//
// As a special note, using feature macros from config.h to determine whether
// to include a particular header requires violating the style guide's required
// ordering for headers: this is permitted.


#ifndef S2_THIRD_PARTY_ABSL_BASE_CONFIG_H_
#define S2_THIRD_PARTY_ABSL_BASE_CONFIG_H_

// Included for the __GLIBC__ macro (or similar macros on other systems).
#include <climits>

#ifdef __cplusplus
// Included for __GLIBCXX__, _LIBCPP_VERSION
#include <cstddef>
#endif  // __cplusplus

// If we're using glibc, make sure we meet a minimum version requirement
// before we proceed much further.
//
// We have chosen glibc 2.12 as the minimum as it was tagged for release
// in May, 2010 and includes some functionality used in Google software
// (for instance pthread_setname_np):
// https://sourceware.org/ml/libc-alpha/2010-05/msg00000.html
#ifdef __GLIBC_PREREQ
#if !__GLIBC_PREREQ(2, 12)
#error "Minimum required version of glibc is 2.12."
#endif
#endif

// ABSL_HAVE_BUILTIN is a function-like feature checking macro.
// It's a wrapper around __has_builtin, which is defined by only clang now.
// It evaluates to 1 if the builtin is supported or 0 if not.
// Define it to avoid an extra level of #ifdef __has_builtin check.
// http://releases.llvm.org/3.3/tools/clang/docs/LanguageExtensions.html
#ifdef __has_builtin
#define ABSL_HAVE_BUILTIN(x) __has_builtin(x)
#else
#define ABSL_HAVE_BUILTIN(x) 0
#endif

// GOOGLE_HAVE_SIZED_DELETE is defined when C++14's sized deallocation
// operators are available.
#if (defined(__clang__) && defined(__cpp_sized_deallocation)) || \
    defined(__GXX_DELETE_WITH_SIZE__)
#define GOOGLE_HAVE_SIZED_DELETE
#endif

// GOOGLE_HAVE_STD_IS_TRIVIALLY_DESTRUCTIBLE is defined when
// std::is_trivially_destructible<T> is supported.
//
// All supported compilers using libc++ have it, as does gcc >= 4.8
// using libstdc++, as does Visual Studio.
// https://gcc.gnu.org/onlinedocs/gcc-4.8.1/libstdc++/manual/manual/status.html#status.iso.2011
// is the first version where std::is_trivially_destructible no longer
// appeared as missing in the Type properties row.
#if defined(_LIBCPP_VERSION) ||                                          \
    (!defined(__clang__) && defined(__GNUC__) && defined(__GLIBCXX__) && \
     (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 8))) ||        \
    defined(_MSC_VER)
#define GOOGLE_HAVE_STD_IS_TRIVIALLY_DESTRUCTIBLE
#endif

// GOOGLE_HAVE_STD_IS_TRIVIALLY_CONSTRUCTIBLE is defined when
// std::is_trivially_default_constructible<T> and
// std::is_trivially_copy_constructible<T> are supported.
//
// GOOGLE_HAVE_STD_IS_TRIVIALLY_ASSIGNABLE is defined when
// std::is_trivially_copy_assignable<T> is supported.
//
// Clang with libc++ supports it, as does gcc >= 5.1 with either
// libc++ or libstdc++, as does Visual Studio.
// https://gcc.gnu.org/gcc-5/changes.html lists as new
// "std::is_trivially_constructible, std::is_trivially_assignable
// etc."
#if (defined(__clang__) && defined(_LIBCPP_VERSION)) ||          \
    (!defined(__clang__) && defined(__GNUC__) &&                 \
     (__GNUC__ > 5 || (__GNUC__ == 5 && __GNUC_MINOR__ >= 1)) && \
     (defined(_LIBCPP_VERSION) || defined(__GLIBCXX__))) ||      \
    defined(_MSC_VER)
#define GOOGLE_HAVE_STD_IS_TRIVIALLY_CONSTRUCTIBLE
#define GOOGLE_HAVE_STD_IS_TRIVIALLY_ASSIGNABLE
#endif

// ABSL_HAVE_THREAD_LOCAL is defined when C++11's thread_local is available.
// Clang implements thread_local keyword but Xcode did not support the
// implementation until Xcode 8.
#if !defined(__apple_build_version__) || __apple_build_version__ >= 8000042
#define ABSL_HAVE_THREAD_LOCAL
#endif

// ABSL_HAVE_INTRINSIC_INT128 is defined when the implementation provides the
// 128 bit integral type: __int128.
//
// __SIZEOF_INT128__ is defined by Clang and GCC when __int128 is supported.
// Clang on ppc64 and aarch64 are exceptions where __int128 exists but has a
// sporadic compiler crashing bug. Nvidia's nvcc also defines __GNUC__ and
// __SIZEOF_INT128__ but not all versions that do this support __int128. Support
// has been tested for versions >= 7.
#if (defined(__clang__) && defined(__SIZEOF_INT128__) &&                 \
     !defined(__ppc64__) && !defined(__aarch64__)) ||                    \
    (defined(__CUDACC__) && defined(__SIZEOF_INT128__) &&                \
     __CUDACC_VER__ >= 70000) ||                                         \
    (!defined(__clang__) && !defined(__CUDACC__) && defined(__GNUC__) && \
     defined(__SIZEOF_INT128__))
#define ABSL_HAVE_INTRINSIC_INT128
#endif

// Operating system-specific features.
//
// Currently supported operating systems and associated preprocessor
// symbols:
//
//   Linux and Linux-derived           __linux__
//   Android                           __ANDROID__ (implies __linux__)
//   Linux (non-Android)               __linux__ && !__ANDROID__
//   Darwin (Mac OS X and iOS)         __APPLE__ && __MACH__
//   Akaros (http://akaros.org)        __ros__
//   Windows                           _WIN32
//   NaCL                              __native_client__
//   AsmJS                             __asmjs__
//   Fuschia                           __Fuchsia__
//
// Note that since Android defines both __ANDROID__ and __linux__, one
// may probe for either Linux or Android by simply testing for __linux__.
//

// GOOGLE_HAVE_FADVISE is defined when the system provides the posix_fadvise(2)
// system call as defined in POSIX.1-2001.
#if defined(__linux__) || defined(__ros__)
#define GOOGLE_HAVE_FADVISE 1
#else
#undef GOOGLE_HAVE_FADVISE
#endif

// GOOGLE_HAVE_FORK is defined when the system provides the fork(2) function
// to create a new process.  fork() has existed on every version of Unix since
// 1970, and many other systems as well.
#if defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) || \
    defined(__ros__) || defined(__native_client__)
#define GOOGLE_HAVE_FORK 1
#else
#undef GOOGLE_HAVE_FORK
#endif
//
// GOOGLE_HAVE_GETPAGESIZE is defined when the system has a getpagesize(2)
// implementation.  Note: getpagesize(2) was removed in POSIX.1-2001.  New
// code should use `sysconf(_SC_PAGESIZE)` instead.
#if defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) || \
    defined(__ros__) || defined(__native_client__)
#define GOOGLE_HAVE_GETPAGESIZE 1
#else
#undef GOOGLE_HAVE_GETPAGESIZE
#endif

// GOOGLE_HAVE_MLOCK is defined when the system has an mlock(2) implementation
// as defined in POSIX.1-2001.
#if defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) || \
    defined(__ros__)
#define GOOGLE_HAVE_MLOCK 1
#else
#undef GOOGLE_HAVE_MLOCK
#endif

// GOOGLE_HAVE_MLOCKALL is defined when the system has an mlockall(2)
// implementation as defined in POSIX.1-2001.  See also mlock(2) and
// GOOGLE_HAVE_MLOCK.
#if (defined(__linux__) && !defined(__ANDROID__)) || \
    (defined(__APPLE__) && defined(__MACH__))
#define GOOGLE_HAVE_MLOCKALL 1
#else
#undef GOOGLE_HAVE_MLOCKALL
#endif

// GOOGLE_HAVE_MMAP is defined when the system has an mmap(2) implementation
// as defined in POSIX.1-2001.
#if defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) ||      \
    defined(__ros__) || defined(__native_client__) || defined(__asmjs__) || \
    defined(__Fuchsia__)
#define GOOGLE_HAVE_MMAP 1
#else
#undef GOOGLE_HAVE_MMAP
#endif

// GOOGLE_HAVE_POSIX_MEMALIGN is defined when the system provides the
// posix_memalign(3) function to allocate memory aligned on a boundary,
// as defined in POSIX.1-2001.
#if (defined(__linux__) && !defined(__ANDROID__)) ||                 \
    (defined(__APPLE__) && defined(__MACH__)) || defined(__ros__) || \
    defined(__native_client__)
#define GOOGLE_HAVE_POSIX_MEMALIGN 1
#else
#undef GOOGLE_HAVE_POSIX_MEMALIGN
#endif

// GOOGLE_HAVE_POSIX_SIGNAL_STACK is defined on systems that provide
// support for separate signals stacks via the sigaltstack(2) call,
// as defined by POSIX.1-2008.  Note that "sigaltstack" looks like a
// typo, but is not: it is "Sig Alt Stack" not "signal stack".
#if defined(__linux__) || (defined(__APPLE__) && defined(__MACH__))
#define GOOGLE_HAVE_POSIX_SIGNAL_ALT_STACK 1
#else
#undef GOOGLE_HAVE_POSIX_SIGNAL_ALT_STACK
#endif

// GOOGLE_HAVE_POSIX_SPAWN is defined when the system provides the
// posix_spawn(3) call, as defined by the POSIX advanced realtime support
// supplement and version 3 of the Single UNIX Specification (SUSv3).
//
// See also Austin T. Clements et al, "The Scalable Commutativity
// Rule: Designing Scalable Software for Multicore Processors"
// (https://people.csail.mit.edu/nickolai/papers/clements-sc.pdf).
#if defined(__linux__) || (defined(__APPLE__) && defined(__MACH__))
#define GOOGLE_HAVE_POSIX_SPAWN 1
#else
#undef GOOGLE_HAVE_POSIX_SPAWN
#endif

// GOOGLE_HAVE_POSIX_TIMER is defined on systems that provide the
// timer_create(2) family of functions as specified in POSIX.1-2001.
#if defined(__linux__) || (defined(__APPLE__) && defined(__MACH__))
#define GOOGLE_HAVE_POSIX_TIMER 1
#else
#undef GOOGLE_HAVE_POSIX_TIMER
#endif

// GOOGLE_HAS_PTHREAD_GETSCHEDPARAM is defined when the system implements the
// pthread_(get|set)schedparam(3) functions as defined in POSIX.1-2001.
#if defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) || \
    defined(__ros__)
// TODO(user): This should probably be called GOOGLE_HAVE_PTHREAD_SCHEDPARAM.
// It may not be grammatically correct, but it would be consistent with other
// symbols.
#define GOOGLE_HAS_PTHREAD_GETSCHEDPARAM 1
#else
#undef GOOGLE_HAS_PTHREAD_GETSCHEDPARAM
#endif

// GOOGLE_HAVE_PTHREAD_SETNAME_NP is defined when the system provides the
// pthread_setname_np(3) function as defined by the behavior of the
// implementation in glibc.  Note that this is a non-standard but common
// extension to the pthreads interface.
#if (defined(__linux__) && !defined(__ANDROID__)) || \
    (defined(__APPLE__) && defined(__MACH__)) || \
    defined(__native_client__) || \
    (defined(__ANDROID__) && defined(__ANDROID_API__) && __ANDROID_API__ >= 10)
#define GOOGLE_HAVE_PTHREAD_SETNAME_NP 1
#else
#undef GOOGLE_HAVE_PTHREAD_SETNAME_NP
#endif

// GOOGLE_HAVE_SCHED_GETCPU is defined when the system implements
// sched_getcpu(3) as by glibc and it's imitators.
#if defined(__linux__) || defined(__ros__)
#define GOOGLE_HAVE_SCHED_GETCPU 1
#else
#undef GOOGLE_HAVE_SCHED_GETCPU
#endif

// GOOGLE_HAVE_SCHED_YIELD is defined when the system implements
// sched_yield(2) as defined in POSIX.1-2001.
#if defined(__linux__) || defined(__ros__) || defined(__native_client__)
#define GOOGLE_HAVE_SCHED_YIELD 1
#else
#undef GOOGLE_HAVE_SCHED_YIELD
#endif

// GOOGLE_HAVE_SEMAPHORE_H is defined when the system supports the <semaphore.h>
// header and sem_open(3) family of functions as standardized in POSIX.1-2001.
//
// Note: While Apple does have <semaphore.h> for both iOS and macOS, it is
// explicity deprecated and will cause build failures if enabled for those
// systems.  We side-step the issue by not defining it here for Apple platforms.
#if defined(__linux__) || defined(__ros__)
#define GOOGLE_HAVE_SEMAPHORE_H 1
#else
#undef GOOGLE_HAVE_SEMAPHORE_H
#endif

// GOOGLE_HAVE_SETGID is defined on systems that provide the setgid(2) system
// call as defined by POSIX.1-2001.
#if defined(__linux__) || (defined(__APPLE__) && defined(__MACH__))
#define GOOGLE_HAVE_SETGID 1
#else
#undef GOOGLE_HAVE_SETGID
#endif

// GOOGLE_HAVE_SETUID is defined on systems that provide the setuid(2) system
// call as defined by POSIX.1-2001.
#if defined(__linux__) || (defined(__APPLE__) && defined(__MACH__))
#define GOOGLE_HAVE_SETUID 1
#else
#undef GOOGLE_HAVE_SETUID
#endif

// GOOGLE_HAVE_SIGINFO_T is defined when the implementation provides the
// siginfo_t for the sigaction(2) interface, as standardized in POSIX.1-2001.
#if defined(__linux__) || (defined(__APPLE__) && defined(__MACH__)) || \
    defined(__ros__)
#define GOOGLE_HAVE_SIGINFO_T 1
#else
#undef GOOGLE_HAVE_SIGINFO_T
#endif

// ABI-specific features.
#if defined(__x86_64__)
// x86_64 System V/Intel C++ ABI-specific definitions.
#elif defined(_M_X64)
// x86_64 Win32 ABI-specific definitions.
#elif defined(__aarch64__)
// ARMv8 aarch64 EABI specific definitions.
#elif defined(__arm__)
// ARM EABI specific definitions.
#elif defined(__ppc64__)
// PowerPC64 ELF specific definitions.
#elif defined(__ppc__)
// PowerPC ELF specific definitions.
#endif

// Library-specific features.
#if defined(__GOOGLE_GRTE_VERSION__)
// feature tests for Google's GRTE
#define GOOGLE_HAVE_ALARM 1
#elif defined(__GLIBC__)
// feature test for glibc
#define GOOGLE_HAVE_ALARM 1
#elif defined(_MSC_VER)
// feature tests for Microsoft's library
#undef GOOGLE_HAVE_ALARM
#elif defined(__native_client__)
#undef GOOGLE_HAVE_ALARM
#else
// other standard libraries
#define GOOGLE_HAVE_ALARM 1
#endif

// POSIX library support (not part of the C/C++ standard).
#if defined(_POSIX_C_SOURCE) && _POSIX_C_SOURCE >= 2
#define GOOGLE_HAVE_FNMATCH 1
#endif

#if defined(_STLPORT_VERSION)
#error "STLPort is not supported."
#endif

#endif  // S2_THIRD_PARTY_ABSL_BASE_CONFIG_H_
