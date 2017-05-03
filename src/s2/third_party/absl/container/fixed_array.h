// Copyright 2005 Google Inc. All Rights Reserved.
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


#ifndef S2_THIRD_PARTY_ABSL_CONTAINER_FIXED_ARRAY_H_
#define S2_THIRD_PARTY_ABSL_CONTAINER_FIXED_ARRAY_H_

#include <cstddef>
#include <algorithm>
#include <cassert>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <memory>
#include <new>
#include <type_traits>

#include "s2/third_party/absl/base/macros.h"
#include "s2/third_party/absl/base/port.h"
#include "s2/third_party/dynamic_annotations/dynamic_annotations.h"

// A FixedArray<T> represents a non-resizable array of T where the
// length of the array does not need to be a compile time constant.
//
// FixedArray allocates small arrays inline, and large arrays on
// the heap.  It is a good replacement for non-standard and deprecated
// uses of alloca() and variable length arrays (a GCC extension).
//
// FixedArray keeps performance fast for small arrays, because it
// avoids heap operations.  It also helps reduce the chances of
// accidentally overflowing your stack if large input is passed to
// your function.
//
// Also, FixedArray is useful for writing portable code.  Not all
// compilers support arrays of dynamic size.

// Most users should not specify an inline_elements argument and let
// FixedArray<> automatically determine the number of elements
// to store inline based on sizeof(T).
//
// If inline_elements is specified, the FixedArray<> implementation
// will store arrays of length <= inline_elements inline.
//
// Note that a FixedArray constructed with a size_type argument will
// default-initialize its values. This matches the behavior of regular
// arrays, but not of vector. This means trivially constructible types
// (e.g. int, int[4], double) will be left uninitialized, and the others
// will be default-constructed.
//
// Note that if a heap allocation is needed, it is done with global
// ::operator new[]() and ::operator delete[](), not with any class-scope
// overrides in T.  As FixedArray is emulating stack-based variable-length
// arrays, this should not be surprising.
namespace absl {
namespace fixed_array_internal {
const std::size_t kInlineElementsDefault = static_cast<std::size_t>(-1);
}  // namespace fixed_array_internal

template <typename T, std::size_t inline_elements =
                          fixed_array_internal::kInlineElementsDefault>
class FixedArray {
 public:
  // For playing nicely with stl:
  using value_type = T;
  using iterator = T*;
  using const_iterator = const T*;
  using reference = T&;
  using const_reference = const T&;
  using pointer = T*;
  using const_pointer = const T*;
  using difference_type = ptrdiff_t;
  using size_type = size_t;

  // Creates an array object that can store "n" elements.
  // Note that trivially constructible elements will be uninitialized.
  // Please see class-level documentation above.
  explicit FixedArray(size_type n) : rep_(MakeRep(n)) { }

  // Creates an array initialized with the elements from the input
  // range. The size will always be "std::distance(first, last)".
  // REQUIRES: Iter must be a forward_iterator or better.
  template <typename Iter>
  FixedArray(Iter first, Iter last) : rep_(MakeRep(first, last)) { }

  // Create the array from an initializer_list.
  FixedArray(std::initializer_list<T> init_list)
      : FixedArray(init_list.begin(), init_list.end()) {}

  ~FixedArray() {
    CleanUpRep(&rep_);
  }

  FixedArray(const FixedArray&) = delete;
  void operator=(const FixedArray&) = delete;

  // Returns the length of the array.
  size_type size() const { return rep_.size(); }

  // Returns the largest possible value of std::distance(begin(), end()) for a
  // FixedArray<T>.  This is equivalent to the most possible addressable bytes
  // over the number of bytes taken by T.
  constexpr size_type max_size() const {
    return std::numeric_limits<difference_type>::max() / sizeof(value_type);
  }

  // Returns whether or not the array is empty.
  bool empty() const { return size() == 0; }

  // Returns the memory size of the array in bytes.
  size_t memsize() const { return size() * sizeof(value_type); }

  // Returns a pointer to the underlying element array.
  const_pointer data() const { return AsValue(rep_.begin()); }
  pointer data() { return AsValue(rep_.begin()); }

  // An older name for the more standard-friendly .data().
  const_pointer get() const { return data(); }
  pointer get() { return data(); }

  // REQUIRES: 0 <= i < size()
  // Returns a reference to the "i"th element.
  reference operator[](size_type i) {
    assert(i < size());
    return get()[i];
  }

  // REQUIRES: 0 <= i < size()
  // Returns a reference to the "i"th element.
  const_reference operator[](size_type i) const {
    assert(i < size());
    return get()[i];
  }

  iterator begin() { return get(); }
  iterator end() { return get() + size(); }

  const_iterator begin() const { return get(); }
  const_iterator end() const { return get() + size(); }

 private:
  // ----------------------------------------
  // HolderTraits:
  // Wrapper to hold elements of type T for the case where T is an array type.
  // If 'T' is an array type, HolderTraits::type is a struct with a 'T v;'.
  // Otherwise, HolderTraits::type is simply 'T'.
  //
  // Maintainer's Note: The simpler solution would be to simply wrap T in a
  // struct whether it's an array or not: 'struct Holder { T v; };', but
  // that causes some paranoid diagnostics to misfire about uses of get(),
  // believing that 'get()' (aka '&rep_.begin().v') is a pointer to a single
  // element, rather than the packed array that it really is.
  // e.g.:
  //
  //     FixedArray<char> buf(1);
  //     sprintf(buf.get(), "foo");
  //
  //     error: call to int __builtin___sprintf_chk(etc...)
  //     will always overflow destination buffer [-Werror]
  //
  class HolderTraits {
    template <typename U>
    struct SelectImpl {
      using type = U;
      static pointer AsValue(type* p) { return p; }
    };

    // Partial specialization for elements of array type.
    template <typename U, size_t N>
    struct SelectImpl<U[N]> {
      struct Holder { U v[N]; };
      using type = Holder;
      static pointer AsValue(type* p) { return &p->v; }
    };
    using Impl = SelectImpl<value_type>;

   public:
    using type = typename Impl::type;

    static pointer AsValue(type *p) { return Impl::AsValue(p); }

    // TODO(user): fix the type aliasing violation
    // this assertion hints at.
    static_assert(sizeof(type) == sizeof(value_type),
                  "Holder must be same size as value_type");
  };

  using Holder = typename HolderTraits::type;
  static pointer AsValue(Holder *p) { return HolderTraits::AsValue(p); }

  // ----------------------------------------
  // InlineSpace:
  // Allocate some space, not an array of elements of type T, so that we can
  // skip calling the T constructors and destructors for space we never use.
  // How many elements should we store inline?
  //   a. If not specified, use a default of 256 bytes (256 bytes
  //      seems small enough to not cause stack overflow or unnecessary
  //      stack pollution, while still allowing stack allocation for
  //      reasonably long character arrays.
  //   b. Never use 0 length arrays (not ISO C++)
  //
  class InlineSpace {
    using Buffer =
        typename std::aligned_storage<sizeof(Holder), alignof(Holder)>::type;
    static const size_type kDefaultBytes = 256;

    template <std::size_t N, typename Ignored>
    struct Impl {
      static const size_type kSize = N;
      Buffer* get() { return space_; }
      // Annotate left_redzone_, right_redzone_, and the unusable portion of
      // space_ with ANNOTATE_CONTIGUOUS_CONTAINER.
      void Annotate(size_t n, bool creating) const {
#ifdef ADDRESS_SANITIZER
        if (!n) return;
        const void* bot = &left_redzone_;
        const void* beg = &space_[0];
        const void* end = &space_[n];
        const void* top = &right_redzone_ + 1;
        // args: (beg, end, old_mid, new_mid)
        if (creating) {
          ANNOTATE_CONTIGUOUS_CONTAINER(beg, top, top, end);
          ANNOTATE_CONTIGUOUS_CONTAINER(bot, beg, beg, bot);
        } else {
          ANNOTATE_CONTIGUOUS_CONTAINER(beg, top, end, top);
          ANNOTATE_CONTIGUOUS_CONTAINER(bot, beg, bot, beg);
        }
#endif  // ADDRESS_SANITIZER
      }

     private:
      static_assert(kSize > 0, "kSize must be positive");
      ADDRESS_SANITIZER_REDZONE(left_redzone_);
      Buffer space_[kSize];
      ADDRESS_SANITIZER_REDZONE(right_redzone_);
    };

    // specialize for 0-element case: no 'space_' array.
    template <typename Ignored>
    struct Impl<0, Ignored> {
      static const size_type kSize = 0;
      Buffer* get() {
        static Buffer buffer;
        return &buffer;
      }
      void Annotate(size_t n, bool creating) const {}
    };

    // specialize for kInlineElementsDefault case.
    // Use up to kDefaultBytes.
    template <typename Ignored>
    struct Impl<fixed_array_internal::kInlineElementsDefault, Ignored>
        : Impl<kDefaultBytes / sizeof(value_type), Ignored> {};

    using ImplType = Impl<inline_elements, void>;

    ImplType space_;

   public:
    static const size_type kSize = ImplType::kSize;

    Holder* get() { return reinterpret_cast<Holder*>(&space_.get()[0]); }
    void AnnotateConstruct(size_t n) const { space_.Annotate(n, true); }
    void AnnotateDestruct(size_t n) const { space_.Annotate(n, false); }
  };

  static const size_type kInlineElements = InlineSpace::kSize;

  Holder* inline_space() { return inline_space_.get(); }

  // ----------------------------------------
  // Rep:
  // A const Rep object holds FixedArray's size and data pointer.
  //
  class Rep {
   public:
    Rep(size_type n, Holder* p) : n_(n), p_(p) { }
    Holder* begin() const { return p_; }
    Holder* end() const { return p_ + n_; }
    size_type size() const { return n_; }
   private:
    size_type n_;
    Holder* p_;
  };

  void CleanUpRep(const Rep* rep) {
    // Destruction must be in reverse order.
    // Loop optimizes to nothing for trivially destructible T.
    for (Holder* p = rep->end(); p != rep->begin(); )
      (--p)->~Holder();
    if (IsAllocated(rep->size())) {
      Deallocate(rep->begin());
    } else {
      inline_space_.AnnotateDestruct(rep->size());
    }
  }

  Rep MakeRep(size_type n) {
    Holder *pa = IsAllocated(n) ? Allocate(n) : inline_space();
    // Loop optimizes to nothing for trivially constructible T.
    for (Holder *p = pa; p != pa + n; ++p)
      // Note: no parens: default init only.
      // Also note '::' to avoid Holder class placement new operator.
      ::new(p) Holder;
    if (!IsAllocated(n)) {
      inline_space_.AnnotateConstruct(n);
    }
    return Rep(n, pa);
  }

  template <typename Iter>
  Rep MakeRep(Iter first, Iter last, std::forward_iterator_tag) {
    size_type n = std::distance(first, last);
    Holder *pa = IsAllocated(n) ? Allocate(n) : inline_space();
    if (!IsAllocated(n)) {
      inline_space_.AnnotateConstruct(n);
    }
    std::uninitialized_copy(first, last, AsValue(pa));
    return Rep(n, pa);
  }

  template <typename Iter>
  Rep MakeRep(Iter first, Iter last) {
    using IterTraits = std::iterator_traits<Iter>;
    return MakeRep(first, last, typename IterTraits::iterator_category());
  }

  Holder* Allocate(size_type n) {
    return static_cast<Holder*>(::operator new[](n * sizeof(Holder)));
  }

  void Deallocate(Holder* p) {
    return ::operator delete[](p);
  }

  bool IsAllocated(size_type n) const { return n > kInlineElements; }

  // ----------------------------------------
  // Data members
  //
  Rep const rep_;
  InlineSpace inline_space_;
};

}  // namespace absl

// Temporary aliases while moving into the absl namespace.
// TODO(user): Delete temporary aliases after namespace update.
namespace gtl {
namespace fixed_array {
using absl::fixed_array_internal::kInlineElementsDefault;
}
}
using absl::FixedArray;

#endif  // S2_THIRD_PARTY_ABSL_CONTAINER_FIXED_ARRAY_H_
