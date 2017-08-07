//
// Copyright (C) 2005 and onwards Google, Inc.
//
// Author: Jeffrey Rennie
//
// Portable hash_sets and hash_maps.
//
#ifndef THIRD_PARTY_ABSL_BASE_INTERNAL_PORT_HASH_INC_
#define THIRD_PARTY_ABSL_BASE_INTERNAL_PORT_HASH_INC_

#include <cassert>
#include <hash_set>
#include <hash_map>
#include <functional>

#if !defined(_MSC_VER) || !defined(STL_MSVC)
#error port_hash.h should only be included when using microsoft compiler and STL.
#endif

// GNU-based hash containers take 4 template args:
//   KeyType, HashFunc, EqualFunc, Allocator
// MSVC-based hash containers take 3 template args:
//   KeyType, HashCompare, Allocator
// where HashCompare has two operator overloads, a single-arg overload
// for hash, and a 2-arg overload for less-than.

// All the code in this file is dedicated to reconciling these differences.
// See http://www.google.com/search?q=hash_compare for more info.
#define NATIVE_HASH_NAMESPACE stdext

struct PortableHashBase {
  // These two public members are required by msvc.  4 and 8 are the
  // default values.
  static const size_t bucket_size = 4;
  static const size_t min_buckets = 8;
};

namespace msvchash {
template <typename Key>
struct hash;

// These are missing from MSVC.
template<> struct hash<int> {
  size_t operator()(int n) const {
    return static_cast<size_t>(n);
  }
};

template<> struct hash<unsigned int> {
  size_t operator()(unsigned int n) const {
    return static_cast<size_t>(n);
  }
};

// If the 3rd template parameter of the GNU interface (KeyEqual) is
// omitted, then we know that it's using the == operator, so we can
// safely use the < operator.
//
// If the third parameter is specified, then we get a compile time
// error, and we know we have to go back and add some #ifdefs.
template <typename Key, typename Hash>
struct HashAndLessOperator : PortableHashBase {
  bool operator()(const Key& a, const Key& b) const {
    return a < b;
  }
  size_t operator()(const Key& key) const {
    return hasher_(key);
  }
  Hash hasher_;
};

template <class Key,
          class Hash>
class hash_set : public NATIVE_HASH_NAMESPACE
      ::hash_set<Key, HashAndLessOperator<Key, Hash> > {
 public:
  hash_set() {}
  explicit hash_set(int buckets) {}
  typedef equal_to<Key> key_equal;
  size_type bucket_count() const {
    return size() / bucket_size;
  }
};

template <class Key, class Val,
          class Hash>
class hash_map : public NATIVE_HASH_NAMESPACE
      ::hash_map<Key, Val, HashAndLessOperator<Key, Hash> > {
 public:
  hash_map() {}
  explicit hash_map(int buckets) {}
  typedef equal_to<Key> key_equal;
  size_type bucket_count() const {
    return size() / bucket_size;
  }
};

template <class Key,
          class Hash>
class hash_multiset : public NATIVE_HASH_NAMESPACE
      ::hash_multiset<Key, HashAndLessOperator<Key, Hash> > {
 public:
  hash_multiset() {}
  explicit hash_multiset(int buckets) {}
  typedef equal_to<Key> key_equal;
  size_type bucket_count() const {
    return size() / bucket_size;
  }
};

template <class Key, class Val,
          class Hash>
class hash_multimap : public NATIVE_HASH_NAMESPACE
      ::hash_multimap<Key, Val, HashAndLessOperator<Key, Hash> > {
 public:
  hash_multimap() {}
  explicit hash_multimap(int buckets) {}
  typedef equal_to<Key> key_equal;
  size_type bucket_count() const {
    return size() / bucket_size;
  }
};

}  // end namespace msvchash

#endif  // THIRD_PARTY_ABSL_BASE_INTERNAL_PORT_HASH_INC_
