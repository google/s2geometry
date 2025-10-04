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

#ifndef S2_BASE_SPINLOCK_H_
#define S2_BASE_SPINLOCK_H_

#include <atomic>

#include "absl/base/thread_annotations.h"

class ABSL_LOCKABLE ABSL_ATTRIBUTE_WARN_UNUSED SpinLock {
 public:
  SpinLock() = default;
  ~SpinLock() = default;
  SpinLock(SpinLock const&) = delete;
  SpinLock& operator=(SpinLock const&) = delete;

  inline void lock() ABSL_EXCLUSIVE_LOCK_FUNCTION() {
    while (locked_.exchange(true, std::memory_order_acquire)) {
      // Spin.
      continue;
    }
  }

  inline void unlock() ABSL_UNLOCK_FUNCTION() {
    locked_.store(false, std::memory_order_release);
  }

  [[nodiscard]] inline bool IsHeld() const {
    return locked_.load(std::memory_order_relaxed);
  }

 private:
  std::atomic_bool locked_{false};
};

class [[nodiscard]] ABSL_SCOPED_LOCKABLE SpinLockHolder {
 public:
  inline explicit SpinLockHolder(SpinLock& l) ABSL_EXCLUSIVE_LOCK_FUNCTION(l)
      : lock_(l) {
    l.lock();
  }
  inline ~SpinLockHolder() ABSL_UNLOCK_FUNCTION() { lock_.unlock(); }

  SpinLockHolder(const SpinLockHolder&) = delete;
  SpinLockHolder& operator=(const SpinLockHolder&) = delete;

 private:
  SpinLock& lock_;
};

#endif  // S2_BASE_SPINLOCK_H_
