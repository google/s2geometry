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
// #include "s2/third_party/absl/base/config.h"
//
// #ifdef GOOGLE_HAVE_MMAP
// #include "s2/sys/mman.h"
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

// Compiler-specific features.
//
// Note: compilers such as clang and ICC define the symbol __GNUC__.
// Most of these implement GCC-compatible frontends, so if the __GNUC__
// macro is defined we bring in a generic GCC-compatibility layer first
// and compiler-specific includes later.
#if defined(__GNUC__)

#define GOOGLE_GCC_VERSION \
  (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)

#if defined(__clang__)

// clang specific feature tests
#if defined(__cpp_sized_deallocation)
#define GOOGLE_HAVE_SIZED_DELETE
#else
#undef GOOGLE_HAVE_SIZED_DELETE
#endif

// clang only claims to be compatible with libstdc++ 4.2, which doesn't
// implement is_trivially_xxx traits
#ifndef __GLIBCXX__
#define GOOGLE_HAVE_STD_IS_TRIVIALLY_DESTRUCTIBLE
#define GOOGLE_HAVE_STD_IS_TRIVIALLY_CONSTRUCTIBLE
#define GOOGLE_HAVE_STD_IS_TRIVIALLY_ASSIGNABLE
#endif

#else  // not clang, so assume gcc

// gcc specific feature tests
#if defined(__GXX_DELETE_WITH_SIZE__)
#define GOOGLE_HAVE_SIZED_DELETE
#else
#undef GOOGLE_HAVE_SIZED_DELETE
#endif  // __GXX_DELETE_WITH_SIZE__

#if defined(__GLIBCXX__)  // libstdc++
// using gcc and libstdc++
// std::is_trivially_destructible is implemented since libstdc++ 4.8.
#if GOOGLE_GCC_VERSION >= 40800
#define GOOGLE_HAVE_STD_IS_TRIVIALLY_DESTRUCTIBLE
#endif
// std::is_trivially_constructible and std::is_trivially_assignable are
// implemented since libstdc++ 5.1.
#if GOOGLE_GCC_VERSION >= 50100
#define GOOGLE_HAVE_STD_IS_TRIVIALLY_CONSTRUCTIBLE
#define GOOGLE_HAVE_STD_IS_TRIVIALLY_ASSIGNABLE
#endif

#elif defined(_LIBCPP_VERSION)  // libc++
// using gcc and libc++
// libc++ fully implements is_trivially_destructible if GCC version >= 4.3.
#define GOOGLE_HAVE_STD_IS_TRIVIALLY_DESTRUCTIBLE
// libc++ fully implements is_trivially_constructible and
// is_trivially_assignable if GCC version >= 5.1.
#if GOOGLE_GCC_VERSION >= 50100
#define GOOGLE_HAVE_STD_IS_TRIVIALLY_CONSTRUCTIBLE
#define GOOGLE_HAVE_STD_IS_TRIVIALLY_ASSIGNABLE
#endif

#endif  // libstdc++ or libc++
#endif  // gcc or clang

#elif defined(_MSC_VER)

#undef GOOGLE_HAVE_SIZED_DELETE

// is_trivially_xxx are supported since Visual Studio 2012.
#define GOOGLE_HAVE_STD_IS_TRIVIALLY_DESTRUCTIBLE
#define GOOGLE_HAVE_STD_IS_TRIVIALLY_CONSTRUCTIBLE
#define GOOGLE_HAVE_STD_IS_TRIVIALLY_ASSIGNABLE

#endif  // _MSC_VER

// Operating system-specific features.
#if defined(__linux__)
#define GOOGLE_HAVE_GETPAGESIZE 1
#define GOOGLE_HAVE_MLOCK 1
#define GOOGLE_HAVE_MMAP 1
#define GOOGLE_HAS_PTHREAD_GETSCHEDPARAM 1
#define GOOGLE_HAVE_SCHED_GETCPU 1
#define GOOGLE_HAVE_SCHED_YIELD 1
#define GOOGLE_HAVE_SEMAPHORE_H 1
#define GOOGLE_HAVE_SIGINFO_T 1
#elif defined(__APPLE__) && defined(__MACH__)
#define GOOGLE_HAVE_GETPAGESIZE 1
#define GOOGLE_HAVE_MLOCK 1
#define GOOGLE_HAVE_MMAP 1
#define GOOGLE_HAS_PTHREAD_GETSCHEDPARAM 1
#undef GOOGLE_HAVE_SCHED_GETCPU
#define GOOGLE_HAVE_SCHED_YIELD 1
#undef GOOGLE_HAVE_SEMAPHORE_H
#define GOOGLE_HAVE_SIGINFO_T 1
#elif defined(__ros__)
#define GOOGLE_HAVE_GETPAGESIZE 1
#define GOOGLE_HAVE_MLOCK 1
#define GOOGLE_HAVE_MMAP 1
#define GOOGLE_HAS_PTHREAD_GETSCHEDPARAM 1
#define GOOGLE_HAVE_SCHED_GETCPU 1
#define GOOGLE_HAVE_SCHED_YIELD 1
#define GOOGLE_HAVE_SEMAPHORE_H 1
#define GOOGLE_HAVE_SIGINFO_T 1
#elif defined(_WIN32)
#undef GOOGLE_HAVE_GETPAGESIZE
#undef GOOGLE_HAVE_MLOCK
#undef GOOGLE_HAVE_MMAP
#undef GOOGLE_HAS_PTHREAD_GETSCHEDPARAM
#undef GOOGLE_HAVE_SCHED_GETCPU
#undef GOOGLE_HAVE_SCHED_YIELD
#undef GOOGLE_HAVE_SEMAPHORE_H
#undef GOOGLE_HAVE_SIGINFO_T
#elif defined(__native_client__)
#define GOOGLE_HAVE_GETPAGESIZE 1
#undef GOOGLE_HAVE_MLOCK
#define GOOGLE_HAVE_MMAP 1
#undef GOOGLE_HAS_PTHREAD_GETSCHEDPARAM
#undef GOOGLE_HAVE_SCHED_GETCPU
#define GOOGLE_HAVE_SCHED_YIELD 1
#undef GOOGLE_HAVE_SEMAPHORE_H
#undef GOOGLE_HAVE_SIGINFO_T
#elif defined(__asmjs__)
#undef GOOGLE_HAVE_GETPAGESIZE
#undef GOOGLE_HAVE_MLOCK
#undef GOOGLE_HAVE_MMAP
#undef GOOGLE_HAS_PTHREAD_GETSCHEDPARAM
#undef GOOGLE_HAVE_SCHED_GETCPU
#undef GOOGLE_HAVE_SCHED_YIELD
#undef GOOGLE_HAVE_SEMAPHORE_H
#undef GOOGLE_HAVE_SIGINFO_T
#elif defined(__Fuchsia__)
#undef GOOGLE_HAVE_GETPAGESIZE
#undef GOOGLE_HAVE_MLOCK
#define GOOGLE_HAVE_MMAP 1
#undef GOOGLE_HAS_PTHREAD_GETSCHEDPARAM
#undef GOOGLE_HAVE_SCHED_GETCPU
#undef GOOGLE_HAVE_SCHED_YIELD
#undef GOOGLE_HAVE_SEMAPHORE_H
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
