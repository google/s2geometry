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

#ifndef S2_UTIL_GTL_UNALIGNED_H_
#define S2_UTIL_GTL_UNALIGNED_H_

#include <array>
#include <cstddef>
#include <cstring>
#include <type_traits>
#include <utility>

#include "absl/base/casts.h"
#include "s2/util/gtl/type_traits.h"
#include "s2/util/gtl/unaligned_internal.h"

namespace gtl {

// Loads the byte representation (in host endianness) from the address pointed
// by 'data'. [data, data + sizeof(T)) address space must be valid. It does not
// matter whether 'data' is properly aligned for T or not.
//
// Unchecked runtime constraint:
//   * [data, data + sizeof(T)) address space must be valid.
//
// Compile-time constraints:
//   * T must be bit-castable (absl::bit_cast).
//   * Must specify template argument T.
//   * Type of 'data' must be one of char*, unsigned char*, std::byte* or void*,
//     or a const-qualified version of one of those types.
template <typename T>
T UnalignedLoad(const char* data) {
  static_assert(internal_unaligned::IsBitCastableTo<T>);
  std::array<char, sizeof(T)> tmp;
  std::memcpy(&tmp, data, sizeof(T));
  return absl::bit_cast<T>(tmp);
}
template <typename T>
T UnalignedLoad(const unsigned char* data) {
  return UnalignedLoad<T>(reinterpret_cast<const char*>(data));
}
template <typename T>
T UnalignedLoad(const std::byte* data) {
  return UnalignedLoad<T>(reinterpret_cast<const char*>(data));
}
template <typename T, int&... ExplicitArgumentBarrier, typename U,
          typename = std::enable_if_t<std::is_same_v<U, void>>>
T UnalignedLoad(const U* data) {
  return UnalignedLoad<T>(reinterpret_cast<const char*>(data));
}

// Stores the byte representation (in host endianness) of 'value' to the address
// pointed by 'data'. It does not matter whether 'data' is properly aligned for
// T or not.
//
// Unchecked runtime constraint:
//   * [data, data + sizeof(T)) address space must be valid.
//
// Compile-time constraints:
//   * T must be bit-castable (absl::bit_cast).
//   * Must specify template argument T.
//   * Type of 'data' must be one of char*, unsigned char*, std::byte* or void*.
template <typename T>
void UnalignedStore(const type_identity_t<T>& value, char* data) {
  static_assert(internal_unaligned::IsBitCastableTo<T>);
  std::memcpy(data, &value, sizeof(T));
}
template <typename T>
void UnalignedStore(const type_identity_t<T>& value, unsigned char* data) {
  UnalignedStore<T>(value, reinterpret_cast<char*>(data));
}
template <typename T>
void UnalignedStore(const type_identity_t<T>& value, std::byte* data) {
  UnalignedStore<T>(value, reinterpret_cast<char*>(data));
}
template <typename T, int&... ExplicitArgumentBarrier, typename U,
          typename = std::enable_if_t<std::is_same_v<U, void>>>
void UnalignedStore(const type_identity_t<T>& value, U* data) {
  UnalignedStore<T>(value, reinterpret_cast<char*>(data));
}

// Unaligned<T> is a wrapper that allows a T value to be stored at an unaligned
// memory location, and then safely accessed with Load() and Store() methods.
// When used for struct fields this provides for packed data structures.
//
// Example:
//
//   struct Unpadded {
//     gtl::Unaligned<uint64> a;
//     gtl::Unaligned<uint16> b;
//     char c;
//   };
//   static_assert(sizeof(Unpadded) == 8 + 2 + 1);
//   static_assert(alignof(Unpadded) == 1);
//   static_assert(sizeof(std::array<Unpadded, 3>) == 11 * 3);
//
// Reading fields:
//
//   uint64 ReadFields(const Unpadded& u) {
//     return u.a.Load() + u.b.Load() + u.c;
//   }
//
// Writing fields:
//
//   void WriteFields(uint64 a, uint16 b, char c, Unpadded& u) {
//     u.a.Store(a);
//     u.b.Store(b);
//     u.c = c;
//   }
//
// Writing the struct into a char buffer:
//
//   void WriteToBuffer(const Unpadded& u, char* data) {
//     gtl::UnalignedStore<Unpadded>(u, data);
//   }
//
// Reading the struct from a char buffer:
//
//   Unpadded ReadFromBuffer(const char* data) {
//     return gtl::UnalignedLoad<Unpadded>(data);
//   }
//
// Constraints:
//   * T must be bit-castable (absl::bit_cast).
//   * T must not be a pointer.
//
// Guarantees:
//   * alignof(Unaligned<T>) == 1 for all T.
//   * sizeof(Unaligned<T>) == sizeof(T).
//   * Unaligned<T> is trivial, trivially default constructible, default
//     constructible iff T is trivial, trivially default constructible, default
//     constructible, respectively.
//   * Unaligned<T> is trivially copy constructible, copy assignable,
//     move constructible, move assignable, destructible, copyable.
//   * Unaligned<T> is bit-castable (absl::bit_cast).
template <typename T>
class Unaligned final : public internal_unaligned::UnalignedBase<T> {
  // TODO(user): File a bug and add it here. go/gtl-unaligned-and-leak-checker.
  static_assert(!std::is_pointer_v<T>);
  static_assert(internal_unaligned::IsBitCastableTo<T>);

 public:
  // Default constructor is default, custom or deleted based on T.
  Unaligned() = default;

  // Trivial copy and assign.
  Unaligned(const Unaligned&) = default;
  Unaligned& operator=(const Unaligned&) = default;

  // Explicit 'from T' constructor.
  explicit Unaligned(const T& value)
      : internal_unaligned::UnalignedBase<T>(std::in_place) {
    Store(value);
  }

  // Trivial destructor.
  ~Unaligned() = default;

  T Load() const { return absl::bit_cast<T>(*this); }
  void Store(const T& value) { *this = absl::bit_cast<Unaligned>(value); }

  // `gtl::Unaligned<T>` can be used as a hash container key if `T` is hashable.
  template <typename H>
  friend H AbslHashValue(H h, const Unaligned& u) {
    return H::combine(std::move(h), u.Load());
  }

  friend bool operator==(const Unaligned& u1, const Unaligned& u2) {
    return u1.Load() == u2.Load();
  }
  friend bool operator!=(const Unaligned& u1, const Unaligned& u2) {
    return u1.Load() != u2.Load();
  }

 private:
  char data_[sizeof(T)];
};

// CTAD deduction guide for making Unaligned<T> objects. E.g. with the Unpadded
// struct above we can write:
//
//   Unpadded MakeUnpadded(uint64 a, uint16 b, char c) {
//     return {.a = gtl::Unaligned(a), .b = gtl::Unaligned(b), .c = c};
//   }
template <class T>
Unaligned(T) -> Unaligned<T>;

}  // namespace gtl

#endif  // S2_UTIL_GTL_UNALIGNED_H_
