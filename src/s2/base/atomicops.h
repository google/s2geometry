// Copyright 2016 Google Inc. All Rights Reserved.
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

#ifndef S2_BASE_ATOMICOPS_H_
#define S2_BASE_ATOMICOPS_H_


#include <atomic>

typedef std::atomic<int32_t> Atomic32;

namespace base {
namespace subtle {

inline int32_t Barrier_AtomicIncrement(Atomic32* atomic, int32_t increment) {
  return *atomic += increment;
}

inline int32_t NoBarrier_Load(Atomic32 const* atomic) {
  return atomic->load(std::memory_order_relaxed);
}

inline void NoBarrier_Store(Atomic32* atomic, int32_t value) {
  atomic->store(value, std::memory_order_relaxed);
}

inline int32_t Acquire_Load(Atomic32 const* atomic) {
  return atomic->load(std::memory_order_acquire);
}

inline void Release_Store(Atomic32* atomic, int32_t value) {
  atomic->store(value, std::memory_order_release);
}

}  // namespace subtle
}  // namespace base

#endif  // S2_BASE_ATOMICOPS_H_
