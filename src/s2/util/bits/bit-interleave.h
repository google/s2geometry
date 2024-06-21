// Copyright 2009 Google Inc. All Rights Reserved.
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

// Author: jyrki@google.com (Jyrki Alakuijala)
//
// Interleaving bits quickly by table lookup.

#ifndef S2_UTIL_BITS_BIT_INTERLEAVE_H_
#define S2_UTIL_BITS_BIT_INTERLEAVE_H_

#include <cstdint>

#include "s2/base/types.h"

namespace util_bits {

// These functions interleave the given arguments into the return value.
//
// The 0-bit in val0 will be the 0-bit in the return value.
// The 0-bit in val1 will be the 1-bit in the return value.
// The 1-bit of val0 will be the 2-bit in the return value, and so on.
uint16_t InterleaveUint8(uint8_t val0, uint8_t val1);
uint32_t InterleaveUint16(uint16_t val0, uint16_t val1);
uint64_t InterleaveUint32(uint32_t val0, uint32_t val1);

// These functions will decode the interleaved values.
void DeinterleaveUint8(uint16_t code, uint8_t *val0, uint8_t *val1);
void DeinterleaveUint16(uint32_t code, uint16_t *val0, uint16_t *val1);
void DeinterleaveUint32(uint64_t code, uint32_t *val0, uint32_t *val1);

// These functions interleave three arguments into the return value.
// The 0-bit in val0 will be the 0-bit in the return value.
// The 0-bit in val1 will be the 1-bit in the return value.
// The 0-bit in val2 will be the 2-bit in the return value.
// The 1-bit of val0 will be the 3-bit in the return value, and so on.
uint32_t InterleaveUint8(uint8_t val0, uint8_t val1, uint8_t val2);

// These functions will decode the interleaved values.
void DeinterleaveUint8(uint32_t code, uint8_t *val0, uint8_t *val1,
                       uint8_t *val2);

}  // namespace util_bits

#endif  // S2_UTIL_BITS_BIT_INTERLEAVE_H_
