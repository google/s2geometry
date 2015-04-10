// Copyright 1999 Google Inc. All Rights Reserved.
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
//
//
// This file contains routines for hashing and fingerprinting.
//
// A hash function takes an arbitrary input bitstring (string, char*,
// number) and turns it into a hash value (a fixed-size number) such
// that unequal input values have a high likelihood of generating
// unequal hash values.  A fingerprint is a hash whose design is
// biased towards avoiding hash collisions, possibly at the expense of
// other characteristics such as execution speed.
//
// In general, if you are only using the hash values inside a single
// executable -- you're not writing the values to disk, and you don't
// depend on another instance of your program, running on another
// machine, generating the same hash values as you -- you want to use
// a HASH.  Otherwise, you want to use a FINGERPRINT.
//
// RECOMMENDED HASH FOR STRINGS:    GoodFastHash
//
// It is a functor, so you can use it like this:
//     hash_map<string, xxx, GoodFastHash<string> >
//     hash_set<char *, GoodFastHash<char*> >
//
// RECOMMENDED HASH FOR NUMBERS:    hash<>
//
// Note that this is likely the identity hash, so if your
// numbers are "non-random" (especially in the low bits), another
// choice is better.  You can use it like this:
//     hash_map<int, xxx>
//     hash_set<uint64>
//
// RECOMMENDED HASH FOR POINTERS:    hash<>
//
// This is a lightly mixed hash. An identity function can be a poor choice
// because pointers can have low entropy in the least significant bits, and
// those bits are important for efficiency in some uses, e.g., dense_hash_map<>.
//
// RECOMMENDED HASH FOR std::pair<T, U>: GoodFastHash<std::pair<T, U>>.
//
// This is significantly faster than hash<std::pair<T, U>> and calls
// GoodFastHash<T> and GoodFastHash<U> (instead of hash<T> and hash<U>), where
// they are available.
//
// RECOMMENDED HASH FOR std::tuple: hash<std::tuple<T, U, ...>>
//
// RECOMMENDED HASH FOR std::array: hash<std::array<T, N>>
//
// RECOMMENDED HASH FOR STRUCTS: ::util_hash::Hash
//
// ::util_hash::Hash(x, y, z) hashes x, y, and z (using GoodFastHash<> if it's
// available; otherwise using hash<>) and then combines the hashes.
//
// RECOMMENDED FINGERPRINT:
//
// For string input, use Fingerprint2011 from fingerprint2011.h. Do *not* use
// Fingerprint in new code; it has problems with excess collisions (see below).
//
// For integer input, Fingerprint is still recommended and has no known
// collision problems.
//
// OTHER HASHES AND FINGERPRINTS:
//
//
// The wiki page also has good advice for when to use a fingerprint vs
// a hash.
//
//
// Note: if your file declares hash_map<string, ...> or
// hash_set<string>, it will use the default hash function,
// hash<string>.  This is not a great choice.  Always provide an
// explicit functor, such as GoodFastHash, as a template argument.
// (Either way, you will need to #include this file to get the
// necessary definition.)
//
// Some of the hash functions below are documented to be fixed
// forever; the rest (whether they're documented as so or not) may
// change over time.  If you require a hash function that does not
// change over time, you should have unittests enforcing this
// property.  We already have several such functions; see
// hash_unittest.cc for the details and unittests.

#ifndef UTIL_HASH_HASH_H_
#define UTIL_HASH_HASH_H_

#include <stddef.h>
#include <stdint.h>     // for uintptr_t
#include <string.h>
#include <algorithm>
#include <ext/hash_map>
using __gnu_cxx::hash;
using __gnu_cxx::hash_map;     // hacky way to make sure we import standard hash<> fns
#include <ext/hash_set>
using __gnu_cxx::hash;
using __gnu_cxx::hash_set;
#include <limits>
#include <string>
#include <utility>

#include "base/casts.h"
#include "base/int128.h"
#include "base/integral_types.h"
#include "base/macros.h"
#include "base/port.h"
#include "base/type_traits.h"
#include "util/hash/builtin_type_hash.h"
#include "util/hash/city.h"
#include "util/hash/hash128to64.h"
#include "util/hash/jenkins.h"
#include "util/hash/jenkins_lookup2.h"
#include "util/hash/string_hash.h"

// TODO(jyrki): This is one of the legacy hashes, kept here only
// temporarily for making it easier to remove the physical dependency
// to legacy_hash.h
inline uint32 HashTo32(const char *s, size_t slen) {
  uint32 retval = Hash32StringWithSeed(s, slen, MIX32);
  return retval == kIllegalHash32 ? static_cast<uint32>(retval-1) : retval;
}

HASH_NAMESPACE_DECLARATION_START

// STLport and MSVC 10.0 above already define these.
#if !defined(_STLP_LONG_LONG) && !(defined(_MSC_VER) && _MSC_VER >= 1600)

#if defined(_MSC_VER)
// MSVC's stl implementation with _MSC_VER less than 1600 doesn't have
// this hash struct. STLport already defines this.
template <typename T>
struct hash {
  size_t operator()(const T& t) const;
};
#endif  // defined(_MSC_VER)

template<> struct hash<int64> {
  size_t operator()(int64 x) const { return static_cast<size_t>(x); }
};

template<> struct hash<uint64> {
  size_t operator()(uint64 x) const { return static_cast<size_t>(x); }
};

#endif  // !defined(_STLP_LONG_LONG) && !(defined(_MSC_VER) && _MSC_VER >= 1600)

template<> struct hash<bool> {
  size_t operator()(bool x) const { return static_cast<size_t>(x); }
};

HASH_NAMESPACE_DECLARATION_END


// ----------------------------------------------------------------------
// Fingerprint()
//   When used for string input, this is not recommended for new code.
//   Instead, use Fingerprint2011(), a higher-quality and faster hash function.
//   (See fingerprint2011.h.) The functions below that take integer input are
//   still recommended.
//
//   Fingerprinting a string (or char*) will never return 0 or 1,
//   in case you want a couple of special values.  However,
//   fingerprinting a numeric type may produce 0 or 1.
// ----------------------------------------------------------------------
extern uint64 FingerprintReferenceImplementation(const char *s, size_t len);
extern uint64 FingerprintInterleavedImplementation(const char *s, size_t len);
inline uint64 Fingerprint(const char *s, size_t len) {
  if (sizeof(s) == 8) {  // 64-bit systems have 8-byte pointers.
    // The better choice when we have a decent number of registers.
    return FingerprintInterleavedImplementation(s, len);
  } else {
    return FingerprintReferenceImplementation(s, len);
  }
}

// Routine that combines together the hi/lo part of a fingerprint
// and changes the result appropriately to avoid returning 0/1.
inline uint64 CombineFingerprintHalves(uint64 hi, uint32 lo) {
  uint64 result = (hi << 32) | lo;
  // (result >> 1) is here the same as (result > 1), but slightly faster.
  if (PREDICT_TRUE(result >> 1)) {
    return result;  // Not 0 or 1, return as is.
  }
  return result ^ GG_ULONGLONG(0x130f9bef94a0a928);
}

inline uint64 Fingerprint(const string& s) {
  return Fingerprint(s.data(), s.size());
}
inline uint64 Hash64StringWithSeed(const string& s, uint64 c) {
  return Hash64StringWithSeed(s.data(), s.size(), c);
}
inline uint64 Fingerprint(schar c) {
  return Hash64NumWithSeed(static_cast<uint64>(c), MIX64);
}
inline uint64 Fingerprint(char c) {
  return Hash64NumWithSeed(static_cast<uint64>(c), MIX64);
}
inline uint64 Fingerprint(uint16 c) {
  return Hash64NumWithSeed(static_cast<uint64>(c), MIX64);
}
inline uint64 Fingerprint(int16 c) {
  return Hash64NumWithSeed(static_cast<uint64>(c), MIX64);
}
inline uint64 Fingerprint(uint32 c) {
  return Hash64NumWithSeed(static_cast<uint64>(c), MIX64);
}
inline uint64 Fingerprint(int32 c) {
  return Hash64NumWithSeed(static_cast<uint64>(c), MIX64);
}
inline uint64 Fingerprint(uint64 c) {
  return Hash64NumWithSeed(static_cast<uint64>(c), MIX64);
}
inline uint64 Fingerprint(int64 c) {
  return Hash64NumWithSeed(static_cast<uint64>(c), MIX64);
}

// This concatenates two 64-bit fingerprints. It is a convenience function to
// get a fingerprint for a combination of already fingerprinted components.
// It assumes that each input is already a good fingerprint itself.
// Note that this is legacy code and new code should use its replacement
// FingerprintCat2011().
//
// Note that in general it's impossible to construct Fingerprint(str)
// from the fingerprints of substrings of str.  One shouldn't expect
// FingerprintCat(Fingerprint(x), Fingerprint(y)) to indicate
// anything about Fingerprint(StrCat(x, y)).
inline uint64 FingerprintCat(uint64 fp1, uint64 fp2) {
  return Hash64NumWithSeed(fp1, fp2);
}

HASH_NAMESPACE_DECLARATION_START

// This intended to be a "good" hash function.  It may change from time to time.
template<> struct hash<uint128> {
  size_t operator()(const uint128& x) const {
    if (sizeof(&x) == 8) {  // 64-bit systems have 8-byte pointers.
      return Hash128to64(x);
    } else {
      uint32 a = static_cast<uint32>(Uint128Low64(x)) +
          static_cast<uint32>(0x9e3779b9UL);
      uint32 b = static_cast<uint32>(Uint128Low64(x) >> 32) +
          static_cast<uint32>(0x9e3779b9UL);
      uint32 c = static_cast<uint32>(Uint128High64(x)) + MIX32;
      mix(a, b, c);
      a += static_cast<uint32>(Uint128High64(x) >> 32);
      mix(a, b, c);
      return c;
    }
  }
  // Less than operator for MSVC use.
  bool operator()(const uint128& a, const uint128& b) const {
    return a < b;
  }
  static const size_t bucket_size = 4;  // These are required by MSVC
  static const size_t min_buckets = 8;  // 4 and 8 are defaults.
};

// Avoid collision with definition in port_hash.h (via port.h).
#ifndef HAVE_DEFINED_HASH_FOR_POINTERS
#define HAVE_DEFINED_HASH_FOR_POINTERS
// Hash pointers as if they were int's, but bring more entropy to
// the lower bits.
template<class T> struct hash<T*> {
  size_t operator()(T *x) const {
    size_t k = reinterpret_cast<size_t>(x);
    return k + (k >> 6);
  }
};
#endif  // HAVE_DEFINED_HASH_FOR_POINTERS

#if defined(__GNUC__)
#if !defined(STLPORT)
// Use our nice hash function for strings
template<class _CharT, class _Traits, class _Alloc>
struct hash<basic_string<_CharT, _Traits, _Alloc> > {
  size_t operator()(const basic_string<_CharT, _Traits, _Alloc>& k) const {
    return Hash32StringWithSeed(k.data(), k.size(), MIX32);
  }
};
#endif  // !defined(STLPORT)

// they don't define a hash for const string at all
template<> struct hash<const string> {
  size_t operator()(const string& k) const {
    return HashStringThoroughly(k.data(), k.size());
  }
};
#endif  // defined(__GNUC__)

// MSVC's STL requires an ever-so slightly different decl
#if defined(STL_MSVC)
template<> struct hash<char const*> {
  size_t operator()(char const* const k) const {
    return HashTo32(k, strlen(k));
  }
  // Less than operator:
  bool operator()(char const* const a, char const* const b) const {
    return strcmp(a, b) < 0;
  }
  static const size_t bucket_size = 4;  // These are required by MSVC
  static const size_t min_buckets = 8;  // 4 and 8 are defaults.
};

// MSVC 10.0 and above have already defined this.
#if !defined(_MSC_VER) || _MSC_VER < 1600
template<> struct hash<string> {
  size_t operator()(const string& k) const {
    return HashTo32(k.data(), k.size());
  }
  // Less than operator:
  bool operator()(const string& a, const string& b) const {
    return a < b;
  }
  static const size_t bucket_size = 4;  // These are required by MSVC
  static const size_t min_buckets = 8;  // 4 and 8 are defaults.
};
#endif  // !defined(_MSC_VER) || _MSC_VER < 1600

#endif  // defined(STL_MSVC)

// Hasher for STL pairs. Requires hashers for both members to be defined.
// Prefer GoodFastHash, particularly if speed is important.
template<class First, class Second>
struct hash<pair<First, Second> > {
  size_t operator()(const pair<First, Second>& p) const {
    size_t h1 = hash<First>()(p.first);
    size_t h2 = hash<Second>()(p.second);
    // The decision below is at compile time
    return (sizeof(h1) <= sizeof(uint32)) ?
            Hash32NumWithSeed(h1, h2)
            : Hash64NumWithSeed(h1, h2);
  }
  // Less than operator for MSVC.
  bool operator()(const pair<First, Second>& a,
                  const pair<First, Second>& b) const {
    return a < b;
  }
  static const size_t bucket_size = 4;  // These are required by MSVC
  static const size_t min_buckets = 8;  // 4 and 8 are defaults.
};

HASH_NAMESPACE_DECLARATION_END

// If you want an excellent string hash function, and you don't mind if it
// might change when you sync and recompile, please use GoodFastHash<>.
// For most applications, GoodFastHash<> is a good choice, better than
// hash<string> or hash<char*> or similar.  GoodFastHash<> can change
// from time to time and may differ across platforms, and we'll strive
// to keep improving it.
//
// By the way, when deleting the contents of a hash_set of pointers, it is
// unsafe to delete *iterator because the hash function may be called on
// the next iterator advance.  Use STLDeleteContainerPointers().

template<class X> struct GoodFastHash;

// This intended to be a "good" hash function.  It may change from time to time.
template<> struct GoodFastHash<char*> {
  size_t operator()(const char* s) const {
    return HashStringThoroughly(s, strlen(s));
  }
  // Less than operator for MSVC.
  bool operator()(const char* a, const char* b) const {
    return strcmp(a, b) < 0;
  }
  static const size_t bucket_size = 4;  // These are required by MSVC
  static const size_t min_buckets = 8;  // 4 and 8 are defaults.
};

// This intended to be a "good" hash function.  It may change from time to time.
template<> struct GoodFastHash<const char*> {
  size_t operator()(const char* s) const {
    return HashStringThoroughly(s, strlen(s));
  }
  // Less than operator for MSVC.
  bool operator()(const char* a, const char* b) const {
    return strcmp(a, b) < 0;
  }
  static const size_t bucket_size = 4;  // These are required by MSVC
  static const size_t min_buckets = 8;  // 4 and 8 are defaults.
};

// This intended to be a "good" hash function.  It may change from time to time.
template<class _CharT, class _Traits, class _Alloc>
struct GoodFastHash<basic_string<_CharT, _Traits, _Alloc> > {
  size_t operator()(const basic_string<_CharT, _Traits, _Alloc>& k) const {
    return HashStringThoroughly(k.data(), k.size() * sizeof(k[0]));
  }
  // Less than operator for MSVC.
  bool operator()(const basic_string<_CharT, _Traits, _Alloc>& a,
                  const basic_string<_CharT, _Traits, _Alloc>& b) const {
    return a < b;
  }
  static const size_t bucket_size = 4;  // These are required by MSVC
  static const size_t min_buckets = 8;  // 4 and 8 are defaults.
};

// This intended to be a "good" hash function.  It may change from time to time.
template<class _CharT, class _Traits, class _Alloc>
struct GoodFastHash<const basic_string<_CharT, _Traits, _Alloc> > {
  size_t operator()(const basic_string<_CharT, _Traits, _Alloc>& k) const {
    return HashStringThoroughly(k.data(), k.size() * sizeof(k[0]));
  }
  // Less than operator for MSVC.
  bool operator()(const basic_string<_CharT, _Traits, _Alloc>& a,
                  const basic_string<_CharT, _Traits, _Alloc>& b) const {
    return a < b;
  }
  static const size_t bucket_size = 4;  // These are required by MSVC
  static const size_t min_buckets = 8;  // 4 and 8 are defaults.
};

// This intended to be a "good" hash function; it is much faster than
// hash<std::pair<T, U>>.  It may change from time to time.  Requires
// GoodFastHash<> or hash<> to be defined for T and U.
template<class T, class U>
struct GoodFastHash<std::pair<T, U> > {
  size_t operator()(const std::pair<T, U>& k) const {
    size_t h1 = HashPart(k.first, /* dummy */ 0);
    size_t h2 = HashPart(k.second, /* dummy */ 0);

    // Mix the hashes together.  Multiplicative hashing mixes the high-order
    // bits better than the low-order bits, and rotating moves the high-order
    // bits down to the low end, where they matter more for most hashtable
    // implementations.
    static const size_t kMul = static_cast<size_t>(0xc6a4a7935bd1e995ULL);
    if (base::is_integral<T>::value) {
      // We want to avoid GoodFastHash({x, y}) == 0 for common values of {x, y}.
      // hash<X> is the identity function for integral types X, so without this,
      // GoodFastHash({0, 0}) would be 0.
      h1 += 109;
    }
    h1 = h1 * kMul;
    h1 = (h1 << 21) | (h1 >> (std::numeric_limits<size_t>::digits - 21));
    return h1 + h2;
  }

  // Less than operator for MSVC.
  bool operator()(const std::pair<T, U>& a, const std::pair<T, U>& b) const {
    return a < b;
  }
  static const size_t bucket_size = 4;  // These are required by MSVC
  static const size_t min_buckets = 8;  // 4 and 8 are defaults.

 private:
  // Use GoodFastHash<V> if it exists and has an operator(); otherwise, use
  // hash<V>.  This relies on decltype, which requires C++11.
#ifdef LANG_CXX11
  template<class V, decltype(GoodFastHash<V>()(std::declval<V>())) = 0>
  size_t HashPart(const V& v, int) const {
    return GoodFastHash<V>()(v);
  }
#endif

  template<class V>
  size_t HashPart(const V& v, ...) const {
    return HASH_NAMESPACE::hash<V>()(v);
  }
};

#ifdef LANG_CXX11
#include <array>
#include <tuple>

HASH_NAMESPACE_DECLARATION_START

namespace hash_internal {

// Lightly hash two hash codes together. When used repetitively to mix more
// than two values, the new values should be in the first argument.
inline size_t Mix(size_t new_hash, size_t accu) {
  static const size_t kMul = static_cast<size_t>(0xc6a4a7935bd1e995ULL);
  // Multiplicative hashing will mix bits better in the msb end ...
  accu *= kMul;
  // ... and rotating will move the better mixed msb-bits to lsb-bits.
  return ((accu << 21) |
          (accu >> (std::numeric_limits<size_t>::digits - 21))) +
      new_hash;
}

// HashPart uses operator() from GoodFastHash<V>, preferably, or hash<V>.
template <class V, decltype(GoodFastHash<V>()(std::declval<V>())) = 0>
inline size_t HashPart(const V& v, int) {return GoodFastHash<V>()(v); }

template <class V>
inline size_t HashPart(const V& v, ...) { return hash<V>()(v); }

// Hash a single value using either GoodFastHash<> or hash<>.
template <typename V> static size_t InternalHash(const V& v) {
  return HashPart(v, /* dummy */ 0);
}

template <int N> struct Helper {
  // Hash fields 0 ... N-1 of u.
  template <typename U>
  static size_t TupleHash(const U& u, size_t seed) {
    static_assert(N >= 1, "N == 0 is special cased, and N < 0 is silly.");
    constexpr int n1 = N - 1;
    return Helper<n1>::TupleHash(u, Mix(InternalHash(std::get<n1>(u)), seed));
  }
};

template <> struct Helper<0> {
  template <typename U>
  static size_t TupleHash(const U&, size_t seed) { return seed; }
};

}  // namespace hash_internal

// Hash functions for tuples.  These are intended to be "good" hash functions.
// They may change from time to time.  GoodFastHash<> or hash<> must be defined
// for the tuple elements.
template <typename... T>
struct hash<std::tuple<T...>> {
 public:
  size_t operator()(const std::tuple<T...>& t) const {
    const int kSeed = 113;
    return hash_internal::Helper<sizeof...(T)>::
        template TupleHash<std::tuple<T...>>(t, kSeed);
  }
  // Less than operator for MSVC.
  bool operator()(const std::tuple<T...>& a,
                  const std::tuple<T...>& b) const {
    return a < b;
  }
  static const size_t bucket_size = 4;  // These are required by MSVC
  static const size_t min_buckets = 8;  // 4 and 8 are defaults.
};

// Hash functions for std::array.  These are intended to be "good" hash
// functions.  They may change from time to time.  GoodFastHash<> or hash<> must
// be defined for T.
template <typename T, std::size_t N>
struct hash<std::array<T, N>> {
 public:
  size_t operator()(const std::array<T, N>& t) const {
    const int kSeed = 71;
    return hash_internal::Helper<N>::
        template TupleHash<std::array<T, N>>(t, kSeed);
  }
  // Less than operator for MSVC.
  bool operator()(const std::array<T, N>& a,
                  const std::array<T, N>& b) const {
    return a < b;
  }
  static const size_t bucket_size = 4;  // These are required by MSVC
  static const size_t min_buckets = 8;  // 4 and 8 are defaults.
};

HASH_NAMESPACE_DECLARATION_END

namespace util_hash {
// ::util_hash::Hash(a, b, c, ...) hashes a, b, c, and so on (using
// GoodFastHash<>, if it's available for the given type, otherwise using
// hash<>), and then combines the individual hashes.  This is intended to be a
// pretty good hash function, which may change from time to time.  (Its quality
// mostly depends on the quality of GoodFastHash<> and/or hash<>.)
//
// In the somewhat unusual case of nested calls to Hash(), it is best if
// the new values should appear first in the arguments list.  For example:
//
//  size_t Hash(int x, int y, vector<T> v, vector<T> w) {
//    auto combine = [](size_t h, const T& elem) {
//      return util_hash::Hash(elem, h);  // Note that elem is the first arg.
//    };
//    size_t vh =
//        std::accumulate(v.begin(), v.end(), static_cast<size_t>(0), combine);
//    size_t wh =
//        std::accumulate(w.begin(), w.end(), static_cast<size_t>(0), combine);
//    // Note that x and y come before vh and wh.
//    return util_hash::Hash(x, y, vh, wh);
//  }
//
// A stronger (and slower) way to combine multiple hash codes together is to
// use hash<uint128>.  The order of args in hash<uint128> doesn't matter.  For
// example:
//
//  size_t Hash(T x, U y) {
//    return hash<uint128>()(uint128(util_hash::Hash(x), util_hash::Hash(y)));
//  }

inline size_t Hash() {
  return 113;
}

template <typename First, typename... T>
size_t Hash(const First& f, const T&... t) {
  return HASH_NAMESPACE::hash_internal::Mix(
      HASH_NAMESPACE::hash_internal::InternalHash(f), Hash(t...));
}
}  // namespace util_hash

#endif  // LANG_CXX11

// Extern template declarations.
//
// gcc only for now.  msvc and others: this technique is likely to work with
// your compiler too.  changelists welcome.
//
// This technique is limited to template specializations whose hash key
// functions are declared in this file.

#if defined(__GNUC__)
HASH_NAMESPACE_DECLARATION_START
extern template class hash_set<string>;
extern template class hash_map<string, string>;
HASH_NAMESPACE_DECLARATION_END
#endif  // defined(__GNUC__)

#endif  // UTIL_HASH_HASH_H_