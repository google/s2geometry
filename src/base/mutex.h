// Copyright 2015 Google Inc. All Rights Reserved.
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
#ifndef UTIL_GEOMETRY_OPENSOURCE_BASE_MUTEX_H_
#define UTIL_GEOMETRY_OPENSOURCE_BASE_MUTEX_H_

#include <condition_variable>
#include <mutex>

class Mutex {
 public:
  inline Mutex() : lock_(mutex_, std::defer_lock) {}
  ~Mutex() = default;
  Mutex(Mutex const&) = delete;
  Mutex& operator=(Mutex const&) = delete;

  inline void Lock() { lock_.lock(); }
  inline void Unlock() { lock_.unlock(); }

 private:
  std::mutex mutex_;
  // Use unique_lock because we need to interact with condition_variables.
  std::unique_lock<std::mutex> lock_;

  friend class CondVar;
};

class CondVar {
 public:
  CondVar() = default;
  ~CondVar() = default;
  CondVar(CondVar const&) = delete;
  CondVar& operator=(CondVar const&) = delete;

  inline void Wait(Mutex* mu) { cond_var_.wait(mu->lock_); };
  inline void Signal() { cond_var_.notify_one(); }
  inline void SignalAll() { cond_var_.notify_all(); }

 private:
  std::condition_variable cond_var_;
};

#endif  // UTIL_GEOMETRY_OPENSOURCE_BASE_MUTEX_H_