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

#ifndef S2_UTIL_BITS_BITS_H_
#define S2_UTIL_BITS_BITS_H_

#include <cstdint>

#include "absl/base/optimization.h"
#include "absl/numeric/bits.h"

// Use namespace because this used to be a static class.
namespace Bits {

inline int Log2FloorNonZero(uint32_t n) {
  ABSL_ASSUME(n != 0);
  return absl::bit_width(n) - 1;
}

inline int Log2FloorNonZero64(uint64_t n) {
  ABSL_ASSUME(n != 0);
  return absl::bit_width(n) - 1;
}

inline int Log2Ceiling(uint32_t n) {
  int floor = absl::bit_width(n) - 1;
  if ((n & (n - 1)) == 0) {  // zero or a power of two
    return floor;
  } else {
    return floor + 1;
  }
}

}  // namespace Bits

#endif  // S2_UTIL_BITS_BITS_H_
