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

#ifndef S2_UTIL_RANDOM_SHARED_BIT_GEN_H_
#define S2_UTIL_RANDOM_SHARED_BIT_GEN_H_

#include "absl/base/no_destructor.h"
#include "absl/random/random.h"

namespace util_random {

// A URBG that is implemented using a thread_local bit generator.
class SharedBitGen {
 public:
  SharedBitGen() = default;
  // SharedBitGen is move-only.
  SharedBitGen(SharedBitGen&&) noexcept = default;
  SharedBitGen& operator=(SharedBitGen&&) noexcept = default;
  SharedBitGen(const SharedBitGen&) = delete;
  SharedBitGen& operator=(const SharedBitGen&) = delete;

  using result_type = typename absl::BitGen::result_type;

  static constexpr result_type(min)() { return (absl::BitGen::min)(); }
  static constexpr result_type(max)() { return (absl::BitGen::max)(); }

  result_type operator()() {
    static thread_local absl::NoDestructor<absl::BitGen> bitgen;
    return (*bitgen)();
  }
};

}  // namespace util_random

#endif  // S2_UTIL_RANDOM_SHARED_BIT_GEN_H_
