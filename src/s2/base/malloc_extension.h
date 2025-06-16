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

#ifndef S2_BASE_MALLOC_EXTENSION_H_
#define S2_BASE_MALLOC_EXTENSION_H_

// Some s2geometry functions use TCMalloc extensions.  Here, we provide
// default implementations so s2geometry can be used with other allocators.
//
// See
// https://github.com/google/tcmalloc/blob/master/tcmalloc/malloc_extension.h

#include <cstddef>

extern "C" {

// `nallocx` performs the same size computation that `malloc` would perform.
// The default weak implementation just returns the `size` argument.
// See https://jemalloc.net/jemalloc.3.html
size_t nallocx(size_t size, int flags) noexcept;

// `__size_returning_new` returns a pointer to the allocated memory along
// with the (possibly rounded up) actual allocation size used.  The default
// weak implementation just returns the requested size.
// This is a TCMalloc prototype of P0901R5 "Size feedback in operator new".
// https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2019/p0901r5.html
struct __sized_ptr_t {
  void* p;
  size_t n;
};
__sized_ptr_t __size_returning_new(size_t size);

}  // extern "C"

#endif  // S2_BASE_MALLOC_EXTENSION_H_
