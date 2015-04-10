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
// This is the legacy unified hash library implementation. Its components are
// being split up into smaller, dedicated libraries. What remains here are
// things still being migrated.
//
// To find the implementation of the core Bob Jenkins lookup2 hash, look in
// jenkins.cc.

#include "util/hash/hash.h"

static const uint32 kFingerprintSeed0 = 0xabc;
static const uint32 kFingerprintSeed1 = 0xdef;

uint64 FingerprintReferenceImplementation(const char *s, size_t len) {
  uint32 hi = Hash32StringWithSeed(s, len, kFingerprintSeed0);
  uint32 lo = Hash32StringWithSeed(s, len, kFingerprintSeed1);
  return CombineFingerprintHalves(hi, lo);
}

uint64 FingerprintInterleavedImplementation(const char *s, size_t len) {
  uint32 hi = Hash32StringWithSeed(s, len, kFingerprintSeed0);
  uint32 lo = Hash32StringWithSeed(s, len, kFingerprintSeed1);
  return CombineFingerprintHalves(hi, lo);
}

#if defined(__GNUC__) && !defined(__QNX__)
HASH_NAMESPACE_DECLARATION_START
template class hash_set<string>;
template class hash_map<string, string>;
HASH_NAMESPACE_DECLARATION_END
#endif
