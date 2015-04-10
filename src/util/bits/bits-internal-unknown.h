// Copyright 2007 Google Inc. All Rights Reserved.
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

#ifndef UTIL_BITS_BITS_INTERNAL_UNKNOWN_H__
#define UTIL_BITS_BITS_INTERNAL_UNKNOWN_H__

inline int Bits::Log2Floor(uint32 n) {
  return Bits::Log2Floor_Portable(n);
}

inline int Bits::Log2FloorNonZero(uint32 n) {
  return Bits::Log2FloorNonZero_Portable(n);
}

inline int Bits::FindLSBSetNonZero(uint32 n) {
  return Bits::FindLSBSetNonZero_Portable(n);
}

inline int Bits::Log2Floor64(uint64 n) {
  return Bits::Log2Floor64_Portable(n);
}

inline int Bits::Log2FloorNonZero64(uint64 n) {
  return Bits::Log2FloorNonZero64_Portable(n);
}

inline int Bits::FindLSBSetNonZero64(uint64 n) {
  return Bits::FindLSBSetNonZero64_Portable(n);
}

#endif  // UTIL_BITS_BITS_INTERNAL_UNKNOWN_H__