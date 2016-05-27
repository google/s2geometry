// Copyright 2006 Google Inc. All Rights Reserved.
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

// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// ---
//
//
// ManualConstructor statically-allocates space in which to store some
// object, but does not initialize it.  You can then call the constructor
// and destructor for the object yourself as you see fit.  This is useful
// for memory management optimizations, where you want to initialize and
// destroy an object multiple times but only allocate it once.
//
// (When I say ManualConstructor statically allocates space, I mean that
// the ManualConstructor object itself is forced to be the right size.)
//
// For example usage, check out util/gtl/small_map.h.
//

#ifndef S2_UTIL_GTL_MANUAL_CONSTRUCTOR_H_
#define S2_UTIL_GTL_MANUAL_CONSTRUCTOR_H_

#include <cstddef>
#include <new>
#include <utility>

#include "s2/base/port.h"

namespace google {

template <typename Type>
class ManualConstructor {
 public:
  // No constructor or destructor because one of the most useful uses of
  // this class is as part of a union, and members of a union cannot have
  // constructors or destructors.  And, anyway, the whole point of this
  // class is to bypass these.

  // Support users creating arrays of ManualConstructor<>s.  This ensures that
  // the array itself has the correct alignment.
  static void* operator new[](size_t size) {
    return aligned_malloc(size, BASE_PORT_H_ALIGN_OF(Type));
  }
  static void operator delete[](void* mem) {
    aligned_free(mem);
  }

  inline Type* get() {
    return reinterpret_cast<Type*>(space_);
  }
  inline const Type* get() const  {
    return reinterpret_cast<const Type*>(space_);
  }

  inline Type* operator->() { return get(); }
  inline const Type* operator->() const { return get(); }

  inline Type& operator*() { return *get(); }
  inline const Type& operator*() const { return *get(); }

  inline void Init() {
    new(space_) Type;
  }

  // Init() constructs the Type instance using the given arguments
  // (which are forwarded to Type's constructor). In C++11, Init() can
  // take any number of arguments of any type, and forwards them perfectly.
  // On pre-C++11 platforms, it can take up to 11 arguments, and may not be
  // able to forward certain kinds of arguments.
  //
  // Note that Init() with no arguments performs default-initialization,
  // not zero-initialization (i.e it behaves the same as "new Type;", not
  // "new Type();"), so it will leave non-class types uninitialized.
  //
  // Rvalue references and std::forward are used by arbiter permission;
#ifdef LANG_CXX11
  template <typename ... Ts>
  inline void Init(Ts&&... args) {  // NOLINT
    new(space_) Type(std::forward<Ts>(args)...);  // NOLINT
  }
#else  // !defined(LANG_CXX11)
  template <typename T1>
  inline void Init(const T1& p1) {
    new(space_) Type(p1);
  }

  template <typename T1, typename T2>
  inline void Init(const T1& p1, const T2& p2) {
    new(space_) Type(p1, p2);
  }

  template <typename T1, typename T2, typename T3>
  inline void Init(const T1& p1, const T2& p2, const T3& p3) {
    new(space_) Type(p1, p2, p3);
  }

  template <typename T1, typename T2, typename T3, typename T4>
  inline void Init(const T1& p1, const T2& p2, const T3& p3, const T4& p4) {
    new(space_) Type(p1, p2, p3, p4);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  inline void Init(const T1& p1, const T2& p2, const T3& p3, const T4& p4,
                   const T5& p5) {
    new(space_) Type(p1, p2, p3, p4, p5);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5,
            typename T6>
  inline void Init(const T1& p1, const T2& p2, const T3& p3, const T4& p4,
                   const T5& p5, const T6& p6) {
    new(space_) Type(p1, p2, p3, p4, p5, p6);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5,
            typename T6, typename T7>
  inline void Init(const T1& p1, const T2& p2, const T3& p3, const T4& p4,
                   const T5& p5, const T6& p6, const T7& p7) {
    new(space_) Type(p1, p2, p3, p4, p5, p6, p7);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5,
            typename T6, typename T7, typename T8>
  inline void Init(const T1& p1, const T2& p2, const T3& p3, const T4& p4,
                   const T5& p5, const T6& p6, const T7& p7, const T8& p8) {
    new(space_) Type(p1, p2, p3, p4, p5, p6, p7, p8);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5,
            typename T6, typename T7, typename T8, typename T9>
  inline void Init(const T1& p1, const T2& p2, const T3& p3, const T4& p4,
                   const T5& p5, const T6& p6, const T7& p7, const T8& p8,
                   const T9& p9) {
    new(space_) Type(p1, p2, p3, p4, p5, p6, p7, p8, p9);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5,
            typename T6, typename T7, typename T8, typename T9, typename T10>
  inline void Init(const T1& p1, const T2& p2, const T3& p3, const T4& p4,
                   const T5& p5, const T6& p6, const T7& p7, const T8& p8,
                   const T9& p9, const T10& p10) {
    new(space_) Type(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10);
  }

  template <typename T1, typename T2, typename T3, typename T4, typename T5,
            typename T6, typename T7, typename T8, typename T9, typename T10,
            typename T11>
  inline void Init(const T1& p1, const T2& p2, const T3& p3, const T4& p4,
                   const T5& p5, const T6& p6, const T7& p7, const T8& p8,
                   const T9& p9, const T10& p10, const T11& p11) {
    new(space_) Type(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11);
  }
#endif  // LANG_CXX11

  inline void Destroy() {
    get()->~Type();
  }

 private:
  ALIGNED_CHAR_ARRAY(Type, 1) space_;
};

}

#endif  // S2_UTIL_GTL_MANUAL_CONSTRUCTOR_H_
