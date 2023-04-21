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
// s2geometry. It is structed into the following high-level categories:
// - Endianness
// - Performance optimization (alignment)

#include <cstring>

#include "s2/base/integral_types.h"

// -----------------------------------------------------------------------------
// Endianness
// -----------------------------------------------------------------------------

// IS_LITTLE_ENDIAN, IS_BIG_ENDIAN
#if defined(__linux__) || defined(__ANDROID__)
#include <endian.h>

#elif defined(__APPLE__)

// BIG_ENDIAN
#include <machine/endian.h>  // NOLINT(build/include)

/* Let's try and follow the Linux convention */
#define __BYTE_ORDER BYTE_ORDER
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

// -----------------------------------------------------------------------------
// Performance Optimization
// -----------------------------------------------------------------------------

// Alignment

// Unaligned APIs

// Portable handling of unaligned loads, stores, and copies. These are simply
// constant-length memcpy calls.
//
// TODO(user): These APIs are forked in Abseil, see
// "third_party/absl/base/internal/unaligned_access.h".
//
// The unaligned API is C++ only.  The declarations use C++ features
// (namespaces, inline) which are absent or incompatible in C.
#if defined(__cplusplus)

namespace base {

// Can't use ATTRIBUTE_NO_SANITIZE_MEMORY because this file is included before
// attributes.h is.
#ifdef __has_attribute
#if __has_attribute(no_sanitize_memory)
#define NO_SANITIZE_MEMORY __attribute__((no_sanitize_memory))
#endif  // __has_attribute(no_sanitize_memory)
#endif  // defined __has_attribute

#ifndef NO_SANITIZE_MEMORY
#define NO_SANITIZE_MEMORY /**/
#endif

template <typename T>
T NO_SANITIZE_MEMORY UnalignedLoad(const void *p) {
  T t;
  memcpy(&t, p, sizeof t);
  return t;
}

#undef NO_SANITIZE_MEMORY

template <typename T>
void UnalignedStore(void *p, T t) {
  memcpy(p, &t, sizeof t);
}
}  // namespace base

inline uint16 UNALIGNED_LOAD16(const void *p) {
  return base::UnalignedLoad<uint16>(p);
}

inline uint32 UNALIGNED_LOAD32(const void *p) {
  return base::UnalignedLoad<uint32>(p);
}

inline uint64 UNALIGNED_LOAD64(const void *p) {
  return base::UnalignedLoad<uint64>(p);
}

inline void UNALIGNED_STORE16(void *p, uint16 v) {
  base::UnalignedStore(p, v);
}

inline void UNALIGNED_STORE32(void *p, uint32 v) {
  base::UnalignedStore(p, v);
}

inline void UNALIGNED_STORE64(void *p, uint64 v) {
  base::UnalignedStore(p, v);
}

#endif  // defined(__cplusplus), end of unaligned API

#endif  // S2_BASE_PORT_H_
