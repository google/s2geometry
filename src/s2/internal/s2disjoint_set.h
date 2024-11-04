// Copyright 2022 Google Inc. All Rights Reserved.
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


#ifndef S2_INTERNAL_S2DISJOINT_SET_H_
#define S2_INTERNAL_S2DISJOINT_SET_H_

#include <optional>

#include "absl/container/flat_hash_map.h"
#include "absl/log/absl_check.h"

namespace internal {

// A disjoint set (AKA union-find set, AKA merge-find set) is a data structure
// that stores a partition of a set into disjoint subsets.  It allows us to
// efficiently add, merge and find a representative member of a set.
//
// In practice, disjoint sets are used for connected-component like algorithms
// where we want to build up relationships between subsets and query the final
// connected pieces.
//
// This implementation uses both path compression and union-by-size, so both the
// Union and FindRoot operations have O(a(N)) amortized complexity, where a(n)
// is the inverse Ackermann function.  For any value that we'll ever care about,
// a(n) is less than 5, which is, for all practical purposes, constant.
template <typename T>
class DisjointSet {
 public:
  // Adds a new element to the set.  The root of the element is itself, that is,
  // it's disjoint from all the other elements of the set.
  //
  // If the element is already in the set, no changes are made and false is
  // returned, otherwise true.
  //
  // Immediately after calling Add(val), FindRoot(val) == val.
  bool Add(const T& val);

  // Finds the root of the given element. If the element isn't in the set,
  // the returned optional has no value, otherwise it will contain the root.
  std::optional<T> FindRoot(const T& val);

  // Performs union of two subsets.  Makes the parent of a's subset equal to the
  // parent of b's subset.  If either a or b isn't in the set, returns false and
  // does not modify the set.
  bool Union(const T& a, const T& b);

  // Returns the total number of elements in the set.
  int Size() const { return nodes_.size(); }

  // Clears the set so that Size() == 0.
  void Clear() { nodes_.clear(); }

  // Reserve space for at least num elements.
  void Reserve(int num) { nodes_.reserve(num); }

 private:
  struct Node {
    explicit Node(const T& value) : value(value) {}
    T value;
    int size = 1;
  };

  // Implementation for find root, expects that val is in the set.  Recursively
  // seeks to the root of the set and sets each element's parent to the root as
  // it goes, returning a reference to the root node.
  Node& FindRootImpl(const T& val);

  absl::flat_hash_map<T, Node> nodes_;
};

template <typename T>
bool DisjointSet<T>::Add(const T& val) {
  return nodes_.emplace(val, Node(val)).second;
}

template <typename T>
std::optional<T> DisjointSet<T>::FindRoot(const T& val) {
  // Val might not be in the table at all, so use find for first lookup.
  auto iter = nodes_.find(val);
  if (iter == nodes_.end()) {
    return {};
  }

  Node& node = iter->second;
  if (node.value == val) {
    return val;
  }
  return node.value = FindRootImpl(node.value).value;
}

template <typename T>
typename DisjointSet<T>::Node& DisjointSet<T>::FindRootImpl(const T& val) {
  auto iter = nodes_.find(val);
  ABSL_DCHECK(iter != nodes_.end());

  Node& node = iter->second;
  if (node.value == val) {
    return node;
  }

  // It can be shown that a disjoint set of size N has an upper bound on its
  // height of floor(log(N)).  So even a billion element set will only recurse
  // 30 levels in the worst case here.
  Node& root = FindRootImpl(node.value);
  node.value = root.value;
  return root;
}

template <typename T>
bool DisjointSet<T>::Union(const T& a, const T& b) {
  auto iter_a = nodes_.find(a);
  if (iter_a == nodes_.end()) {
    return false;
  }

  auto iter_b = nodes_.find(b);
  if (iter_b == nodes_.end()) {
    return false;
  }

  Node& root_a = FindRootImpl(iter_a->second.value);
  Node& root_b = FindRootImpl(iter_b->second.value);

  if (root_a.value != root_b.value) {
    if (root_a.size < root_b.size) {
      root_a.value = root_b.value;
      root_b.size += root_a.size;
    } else {
      root_b.value = root_a.value;
      root_a.size += root_b.size;
    }
  }
  return true;
}

}  // namespace internal

#endif  // S2_INTERNAL_S2DISJOINT_SET_H_
