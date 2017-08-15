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

// Various portability macros, type definitions, and inline functions
// This file is used for both C and C++!
//
// These are weird things we need to do to get this compiling on
// random systems (and on SWIG).
//
// This files is structured into the following high-level categories:
// - Platform checks (OS, Compiler, C++, Library)
// - Utility macros
// - Utility functions
// - Type alias
// - Predefined system/language macros
// - Predefined system/language functions
// - Obsolete
//
// NOTE FOR GOOGLERS:
//
// IWYU pragma: private, include "base/port.h"

#ifndef S2_THIRD_PARTY_ABSL_BASE_PORT_H_
#define S2_THIRD_PARTY_ABSL_BASE_PORT_H_

#include <cassert>
#include <climits>         // So we can set the bounds of our types
#include <cstdlib>         // for free()
#include <cstring>         // for memcpy()

#ifdef _MSC_VER
#include <intrin.h>
#endif

#include "s2/third_party/absl/base/attributes.h"
#include "s2/third_party/absl/base/config.h"
#include "s2/third_party/absl/base/integral_types.h"
#include "s2/third_party/absl/base/optimization.h"

#ifdef SWIG
%include "third_party/absl/base/attributes.h"
#endif

// -----------------------------------------------------------------------------
// Operating System Check
// -----------------------------------------------------------------------------

#if defined(__CYGWIN__)
#error "Cygwin is not supported."
#endif

// -----------------------------------------------------------------------------
// Compiler Check
// -----------------------------------------------------------------------------

// We support MSVC++ 14.0 update 2 and later.
// This minimum will go up.
#if defined(_MSC_FULL_VER) && _MSC_FULL_VER < 190023918
#error "This package requires Visual Studio 2015 Update 2 or higher"
#endif

// We support gcc 4.7 and later.
// This minimum will go up.
// Do not test clang. As of crosstool v18, clang identifies as gcc 4.2.
#if defined(__GNUC__) && !defined(__clang__)
#if __GNUC__ < 4 || (__GNUC__ == 4 && __GNUC_MINOR__ < 7)
#error "This package requires gcc 4.7 or higher"
#endif
#endif

// We support Apple Xcode clang 4.2.1 (version 421.11.65) and later.
// This corresponds to Apple Xcode version 4.5.
// This minimum will go up.
#if defined(__apple_build_version__) && __apple_build_version__ < 4211165
#error "This package requires __apple_build_version__ of 4211165 or higher"
#endif

// -----------------------------------------------------------------------------
// C++ Version Check
// -----------------------------------------------------------------------------

// Enforce C++11 as the minimum.  Note that Visual Studio has not
// advanced __cplusplus despite being good enough for our purposes, so
// so we exempt it from the check.
#if defined(__cplusplus) && !defined(_MSC_VER) && !defined(SWIG)
#if __cplusplus < 201103L
#error "C++ versions less than C++11 are not supported."
#endif
#endif

// -----------------------------------------------------------------------------
// C++ Standard Library Check
// -----------------------------------------------------------------------------

#if defined(__cplusplus)
#include <cstddef>
#if defined(_STLPORT_VERSION)
#error "STLPort is not supported."
#endif
#endif

// -----------------------------------------------------------------------------
// Utility Macros
// -----------------------------------------------------------------------------

// ABSL_FUNC_PTR_TO_CHAR_PTR
// On some platforms, a "function pointer" points to a function descriptor
// rather than directly to the function itself.
// Use ABSL_FUNC_PTR_TO_CHAR_PTR(func) to get a char-pointer to the first
// instruction of the function func.
#if defined(__cplusplus)
#if (defined(__powerpc__) && !(_CALL_ELF > 1)) || defined(__ia64)
// use opd section for function descriptors on these platforms, the function
// address is the first word of the descriptor
namespace absl {
enum { kPlatformUsesOPDSections = 1 };
}  // namespace absl
#define ABSL_FUNC_PTR_TO_CHAR_PTR(func) (reinterpret_cast<char **>(func)[0])
#else  // not PPC or IA64
namespace absl {
enum { kPlatformUsesOPDSections = 0 };
}  // namespace absl
#define ABSL_FUNC_PTR_TO_CHAR_PTR(func) (reinterpret_cast<char *>(func))
#endif  // PPC or IA64
#endif  // __cplusplus

// -----------------------------------------------------------------------------
// Utility Functions
// -----------------------------------------------------------------------------

#if !defined(SWIG)
#if defined(__cplusplus)
namespace absl {
constexpr char PathSeparator() {
#ifdef _WIN32
  return '\\';
#else
  return '/';
#endif
}
}  // namespace absl
#endif  // __cplusplus
#endif

// -----------------------------------------------------------------------------
// Type Alias
// -----------------------------------------------------------------------------

#ifdef _MSC_VER
// uid_t
// MSVC doesn't have uid_t
typedef int uid_t;

// pid_t
// Defined all over the place.
typedef int pid_t;
#endif  // _MSC_VER

// -----------------------------------------------------------------------------
// Predefined System/Language Macros
// -----------------------------------------------------------------------------

// MAP_ANONYMOUS
#if defined(__APPLE__)
// For mmap, Linux defines both MAP_ANONYMOUS and MAP_ANON and says MAP_ANON is
// deprecated. In Darwin, MAP_ANON is all there is.
#if !defined MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif  // !MAP_ANONYMOUS
#endif  // __APPLE__

// PATH_MAX
// You say tomato, I say atotom
#ifdef _MSC_VER
#define PATH_MAX MAX_PATH
#endif

// -----------------------------------------------------------------------------
// Obsolete (to be removed)
// -----------------------------------------------------------------------------

// HAS_GLOBAL_STRING
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

// NOTE: These live in Abseil purely as a short-term layering workaround to
// resolve a dependency chain between util/hash/hash, absl/strings, and //base:
// in order for //base to depend on absl/strings, the includes of hash need
// to be in absl, not //base.  string_view defines hashes.
//
// -----------------------------------------------------------------------------
// HASH_NAMESPACE, HASH_NAMESPACE_DECLARATION_START/END
// -----------------------------------------------------------------------------

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
// TODO(user): make this an error. Do not use external value of
// HASH_NAMESPACE.
#endif

#ifndef HASH_NAMESPACE
// TODO(user): try to delete this.
// I think gcc 2.95.3 was the last toolchain to use this.
#define HASH_NAMESPACE_DECLARATION_START
#define HASH_NAMESPACE_DECLARATION_END
#else
#define HASH_NAMESPACE_DECLARATION_START namespace HASH_NAMESPACE {
#define HASH_NAMESPACE_DECLARATION_END }
#endif

#endif  // S2_THIRD_PARTY_ABSL_BASE_PORT_H_
