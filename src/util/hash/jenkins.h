// Copyright 2011 Google Inc. All Rights Reserved.
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

//
// The core Jenkins Lookup2-based hashing routines. These are legacy hashing
// routines and should be avoided in new code. Their implementations are dated
// and cannot be changed due to values being recorded and breaking if not
// preserved. New code which explicitly desires this property should use the
// consistent hashing libraries. New code which does not explicitly desire this
// behavior should use the generic hashing routines in hash.h.

#ifndef S2GEOMETRY_UTIL_HASH_JENKINS_H_
#define S2GEOMETRY_UTIL_HASH_JENKINS_H_

#include <stdlib.h>
#include "base/integral_types.h"

static const uint32 MIX32 = 0x12b9b0a1UL;           // pi; an arbitrary number
static const uint64 MIX64 = GG_ULONGLONG(0x2b992ddfa23249d6);  // more of pi

const uint32 kIllegalHash32 = static_cast<uint32>(0xffffffffU);
const uint16 kIllegalHash16 = static_cast<uint16>(0xffffU);

// ----------------------------------------------------------------------
// Hash32StringWithSeed()
// Hash64StringWithSeed()
// Hash32NumWithSeed()
// Hash64NumWithSeed()
//   These are Bob Jenkins' hash functions, one for 32 bit numbers
//   and one for 64 bit numbers.  Each takes a string as input and
//   a start seed.  Hashing the same string with two different seeds
//   should give two independent hash values.
//      The *Num*() functions just do a single mix, in order to
//   convert the given number into something *random*.
//
// Note that these methods may return any value for the given size, while
// the corresponding HashToXX() methods avoids certain reserved values.
// ----------------------------------------------------------------------
// These hash functions are no longer recommended.
uint32 Hash32StringWithSeed(const char *s, size_t len, uint32 c);
uint64 Hash64StringWithSeed(const char *s, size_t len, uint64 c);

#endif  // S2GEOMETRY_UTIL_HASH_JENKINS_H_
