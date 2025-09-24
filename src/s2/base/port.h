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
// - Utility functions
// - Performance optimization (alignment)

#include <cstdint>
#include <cstring>
#include <new>

// -----------------------------------------------------------------------------
// Utility Functions
// -----------------------------------------------------------------------------

// C++14 sized deallocation
namespace base {
inline void sized_delete(void *ptr, size_t size) {
  ::operator delete(ptr, size);
}

inline void sized_delete_array(void *ptr, size_t size) {
  ::operator delete[](ptr, size);
}
}  // namespace base

// -----------------------------------------------------------------------------
// Performance Optimization
// -----------------------------------------------------------------------------

// Alignment

// Unaligned APIs

// Portable handling of unaligned loads, stores, and copies. These are simply
// constant-length memcpy calls.

namespace base {
template <typename T>
T UnalignedLoad(const void *p) {
  T t;
  memcpy(&t, p, sizeof t);
  return t;
}

template <typename T>
void UnalignedStore(void *p, T t) {
  memcpy(p, &t, sizeof t);
}
}  // namespace base

inline uint16_t UNALIGNED_LOAD16(const void *p) {
  return base::UnalignedLoad<uint16_t>(p);
}

inline uint32_t UNALIGNED_LOAD32(const void *p) {
  return base::UnalignedLoad<uint32_t>(p);
}

inline uint64_t UNALIGNED_LOAD64(const void *p) {
  return base::UnalignedLoad<uint64_t>(p);
}

inline void UNALIGNED_STORE16(void *p, uint16_t v) {
  base::UnalignedStore(p, v);
}

inline void UNALIGNED_STORE32(void *p, uint32_t v) {
  base::UnalignedStore(p, v);
}

inline void UNALIGNED_STORE64(void *p, uint64_t v) {
  base::UnalignedStore(p, v);
}

#endif  // S2_BASE_PORT_H_
