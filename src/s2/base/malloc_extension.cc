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

#include "s2/base/malloc_extension.h"

#include <cstddef>
#include <new>

#include "absl/base/attributes.h"

ABSL_ATTRIBUTE_WEAK ABSL_ATTRIBUTE_NOINLINE size_t nallocx(size_t size,
                                                           int) noexcept {
  return size;
}

ABSL_ATTRIBUTE_WEAK ABSL_ATTRIBUTE_NOINLINE __sized_ptr_t
__size_returning_new(size_t size) {
  return {::operator new(size), size};
}
