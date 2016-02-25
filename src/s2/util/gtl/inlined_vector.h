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

//
// An InlinedVector<T,N,A> is like a vector<T,A>, except that storage
// for sequences of length <= N are provided inline without requiring
// any heap allocation.  Typically N is very small (e.g., 4) so that
// sequences that are expected to be short do not require allocations.
//
// Only some of the vector<> operations are currently implemented.
// Other operations may be added as needed to facilitate migrating
// code that uses vector<> to InlinedVector<>.
//
// NOTE: If you want an inlined version to replace use of a
// vector<bool>, consider using util::bitmap::InlinedBitVector<NBITS>
// in util/bitmap/inlined_bitvector.h
//

#ifndef S2_UTIL_GTL_INLINED_VECTOR_H_
#define S2_UTIL_GTL_INLINED_VECTOR_H_

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <iterator>
#include <memory>

#include <glog/logging.h>
#include "s2/base/macros.h"
#include "s2/base/port.h"
#include <type_traits>
#include "s2/util/gtl/gtl_namespace.h"
#include "s2/util/gtl/manual_constructor.h"

// Must come after "base/port.h"
#ifdef LANG_CXX11
#include <initializer_list>  // NOLINT(build/include_order)
#endif  // LANG_CXX11


GTL_NAMESPACE_BEGIN

template <typename T, int N, typename A = std::allocator<T> >
class InlinedVector {
 public:
  typedef A allocator_type;
  typedef typename allocator_type::value_type value_type;
  typedef typename allocator_type::pointer pointer;
  typedef typename allocator_type::const_pointer const_pointer;
  typedef typename allocator_type::reference reference;
  typedef typename allocator_type::const_reference const_reference;
  typedef typename allocator_type::size_type size_type;
  typedef typename allocator_type::difference_type difference_type;
  typedef pointer iterator;
  typedef const_pointer const_iterator;
  typedef std::reverse_iterator<iterator> reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

  InlinedVector() : allocator_and_tag_(allocator_type()) {}

  explicit InlinedVector(const allocator_type& alloc)
      : allocator_and_tag_(alloc) {}

  // Create a vector with n copies of value_type().
  explicit InlinedVector(size_type n) : allocator_and_tag_(allocator_type()) {
    InitAssign(n);
  }

  // Create a vector with n copies of elem
  InlinedVector(size_type n, const value_type& elem,
                const allocator_type& alloc = allocator_type())
      : allocator_and_tag_(alloc) {
    InitAssign(n, elem);
  }

  // Create and initialize with the elements [range_start .. range_end).
  // The unused enable_if argument restricts this constructor so that it is
  // elided when value_type is an integral type.  This prevents ambiguous
  // interpretation between a call to this constructor with two integral
  // arguments and a call to the preceding (n, elem) constructor.
  template <typename InputIterator>
  InlinedVector(InputIterator range_start, InputIterator range_end,
                const allocator_type& alloc = allocator_type(),
                typename std::enable_if<
                    !std::is_integral<InputIterator>::value>::type* = nullptr)
      : allocator_and_tag_(alloc) {
    AppendRange(range_start, range_end);
  }

#ifdef LANG_CXX11
  InlinedVector(std::initializer_list<value_type> init,
                const allocator_type& alloc = allocator_type())
      : allocator_and_tag_(alloc) {
    AppendRange(init.begin(), init.end());
  }
#endif  // LANG_CXX11

  InlinedVector(const InlinedVector& v);

  ~InlinedVector() { clear(); }

  InlinedVector& operator=(const InlinedVector& v) {
    // Optimized to avoid reallocation.
    // Prefer reassignment to copy construction for elements.
    if (size() < v.size()) {  // grow
      reserve(v.size());
      std::copy(v.begin(), v.begin() + size(), begin());
      std::copy(v.begin() + size(), v.end(), std::back_inserter(*this));
    } else {  // maybe shrink
      erase(begin() + v.size(), end());
      std::copy(v.begin(), v.end(), begin());
    }
    return *this;
  }

  template <class InputIterator>
  void assign(InputIterator first, InputIterator last,
              typename std::enable_if<
                  !std::is_integral<InputIterator>::value>::type* = nullptr) {
    AssignRange(first, last);
  }

#ifdef LANG_CXX11
  void assign(std::initializer_list<value_type> init) {
    AssignRange(init.begin(), init.end());
  }
#endif  // LANG_CXX11

  void assign(size_type n, const value_type& elem) {
    if (n <= size()) {  // Possibly shrink
      std::fill_n(begin(), n, elem);
      erase(begin() + n, end());
      return;
    }
    // Grow
    reserve(n);
    std::fill_n(begin(), size(), elem);
    if (allocated()) {
      UninitializedFillAllocated(allocated_space() + size(),
                                 allocated_space() + n, elem);
    } else {
      UninitializedFillInlined(inlined_space() + size(),
                               inlined_space() + n, elem);
    }
    set_size_internal(n);
  }

  size_type size() const {
    return allocated() ? allocation().size() : tag().size();
  }

  bool empty() const { return (size() == 0); }

  // Return number of elements that can be stored in vector
  // without requiring a reallocation of underlying memory
  size_type capacity() const {
    return allocated() ? allocation().capacity() : N;
  }

  // Return a pointer to the underlying array.
  // Only result[0,size()-1] are defined.
  const_pointer data() const {
    return allocated() ? allocated_space() : inlined_space();
  }
  pointer data() {
    return allocated() ? allocated_space() : inlined_space();
  }

  // An older name for the more standard-friendly .data().
  const_pointer array() const { return data(); }
  pointer mutable_array() { return data(); }

  // Remove all elements
  void clear() {
    size_type s = size();
    if (allocated()) {
      DestroyAllocated(allocated_space(), allocated_space() + s);
      allocation().Dealloc(allocator());
    } else {
      DestroyInlined(inlined_space(), inlined_space() + s);
    }
    tag() = Tag();
  }

  // Return the ith element
  // REQUIRES: 0 <= i < size()
  const value_type& at(size_type i) const {
    DCHECK_LT(i, size());
    return array()[i];
  }
  const value_type& operator[](size_type i) const {
    DCHECK_LT(i, size());
    return array()[i];
  }

  // Return a non-const reference to the ith element
  // REQUIRES: 0 <= i < size()
  value_type& at(size_type i) {
    DCHECK_LT(i, size());
    return mutable_array()[i];
  }
  value_type& operator[](size_type i) {
    DCHECK_LT(i, size());
    return mutable_array()[i];
  }

  value_type& back() {
    DCHECK(!empty());
    return at(size() - 1);
  }

  const value_type& back() const {
    DCHECK(!empty());
    return at(size() - 1);
  }

  value_type& front() {
    DCHECK(!empty());
    return at(0);
  }

  const value_type& front() const {
    DCHECK(!empty());
    return at(0);
  }

  // Append t to the vector.
  // Increases size() by one.
  // Amortized complexity: O(1)
  // Worst-case complexity: O(size())
  void push_back(const value_type& t) {
    size_type s = size();
    DCHECK_LE(s, capacity());
    if (s == capacity()) {
      return GrowAndPushBack(t);
    }
    DCHECK_LT(s, capacity());

    if (allocated()) {
      ConstructAllocated(allocated_space() + s, t);
    } else {
      ConstructInlined(inlined_space() + s, t);
    }

    set_size_internal(s + 1);
  }

  void pop_back() {
    DCHECK(!empty());
    size_type s = size();
    if (allocated()) {
      DestroyAllocated(allocated_space() + s - 1, allocated_space() + s);
    } else {
      DestroyInlined(inlined_space() + s - 1, inlined_space() + s);
    }
    set_size_internal(s - 1);
  }

  // Resizes the vector to contain "n" elements.
  // If "n" is smaller than the initial size, extra elements are destroyed.
  // If "n" is larger than the initial size, enough copies of "elem"
  // are appended to increase the size to "n". If "elem" is omitted,
  // new elements are value-initialized.
  void resize(size_type n);
  void resize(size_type n, const value_type& elem);

  iterator begin() { return mutable_array(); }
  const_iterator begin() const { return array(); }
  const_iterator cbegin() const { return begin(); }

  iterator end() { return mutable_array() + size(); }
  const_iterator end() const { return array() + size(); }
  const_iterator cend() const { return end(); }

  reverse_iterator rbegin() { return reverse_iterator(end()); }
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }
  const_reverse_iterator crbegin() const { return rbegin(); }

  reverse_iterator rend() { return reverse_iterator(begin()); }
  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }
  const_reverse_iterator crend() const { return rend(); }

  iterator insert(iterator pos, const value_type& v);

  iterator erase(iterator pos) {
    DCHECK_LT(pos, end());
    DCHECK_GE(pos, begin());
    std::copy(pos + 1, end(), pos);
    pop_back();
    return pos;
  }

  iterator erase(iterator first, iterator last);

  // Enlarges the underlying representation so it can hold at least
  // "n" elements without reallocation.
  // Does not change size() or the actual contents of the vector.
  void reserve(size_type n) {
    if (n > capacity()) {
      // Make room for new elements
      EnlargeBy(n - size());
    }
  }

  // Swap the contents of *this with other.
  // REQUIRES: value_type is swappable and copyable.
  void swap(InlinedVector& other);

  allocator_type get_allocator() const { return allocator(); }

 private:
  static_assert(N > 0, "inlined vector with nonpositive size");

  // TODO(user): Some Android NDK builds falsely claim C++11 support.
  // http://test/OCL:59047547:BASE:60708067:1391139009376:52f2bf25
  // This whole class could be replaced with:
  // using AllocatorTraits = std::allocator_traits<allocator_type>;
  struct AllocatorTraits {
    typedef typename allocator_type::value_type value_type;
    typedef typename allocator_type::pointer pointer;
    typedef typename allocator_type::size_type size_type;

    static void construct(allocator_type& a,  // NOLINT(runtime/references)
                          pointer p) {
      // Tricky: do we support non-copyable types, or support allocators
      // that do special things with construct()? Non-copyable types are
      // needed today, so they are more important. When we sort out the
      // Android NDK C++11 problem, we will be able to use the proper
      // std::allocator_traits<A>::construct(p, ...).
      //
      // a.construct(p, value_type());
      new (p) value_type();
    }
    static void construct(allocator_type& a,  // NOLINT(runtime/references)
                          pointer p, const value_type& t) {
      a.construct(p, t);
    }
    static void destroy(allocator_type& a,  // NOLINT(runtime/references)
                        pointer p) {
      a.destroy(p);
    }
    static pointer allocate(allocator_type& a,  // NOLINT(runtime/references)
                            size_type n) {
      return a.allocate(n);
    }
    static void deallocate(allocator_type& a,  // NOLINT(runtime/references)
                           pointer p, size_type n) {
      a.deallocate(p, n);
    }
  };

  // If the vector is inlined, holds the size of the vector.
  // If the vector is allocated, holds the special value kAllocated,
  // and the size is stored in the vector's Allocation.
  class Tag {
   public:
    Tag() : size_(0) {}
    size_type size() const { return size_; }
    void set_size(size_type n) { size_ = n; }
    bool allocated() const { return size_ == kAllocated; }
    void set_allocated() { size_ = kAllocated; }
   private:
    static const size_type kAllocated = -1;
    size_type size_;
  };

  // Derives from allocator_type to use the empty base class optimization.
  // If the allocator_type is stateless, we can 'store'
  // our instance of it for free.
  class AllocatorAndTag : private allocator_type {
   public:
    explicit AllocatorAndTag(const allocator_type& a, Tag t = Tag())
        : allocator_type(a), tag_(t) {
    }
    Tag& tag() { return tag_; }
    const Tag& tag() const { return tag_; }
    allocator_type& allocator() { return *this; }
    const allocator_type& allocator() const { return *this; }
   private:
    Tag tag_;
  };

  class Allocation {
   public:
    Allocation(allocator_type& a,  // NOLINT(runtime/references)
               size_type capacity)
        : size_(0),
          capacity_(capacity),
          buffer_(AllocatorTraits::allocate(a, capacity_)) {}

    void Dealloc(allocator_type& a) {  // NOLINT(runtime/references)
      AllocatorTraits::deallocate(a, buffer(), capacity());
    }

    size_type size() const { return size_; }
    void set_size(size_type s) { size_ = s; }
    size_type capacity() const { return capacity_; }
    const value_type* buffer() const { return buffer_; }
    value_type* buffer() { return buffer_; }

   private:
    size_type size_;
    size_type capacity_;
    value_type* buffer_;
  };

  const Tag& tag() const { return allocator_and_tag_.tag(); }
  Tag& tag() { return allocator_and_tag_.tag(); }

  Allocation& allocation() {
    return *rep_.allocation_storage.allocation.get();
  }
  const Allocation& allocation() const {
    return *rep_.allocation_storage.allocation.get();
  }
  void init_allocation(const Allocation& allocation) {
    rep_.allocation_storage.allocation.Init(allocation);
  }

  value_type* inlined_space() {
    return rep_.inlined_storage.inlined[0].get();
  }
  const value_type* inlined_space() const {
    return rep_.inlined_storage.inlined[0].get();
  }

  value_type* allocated_space() {
    return allocation().buffer();
  }
  const value_type* allocated_space() const {
    return allocation().buffer();
  }

  const allocator_type& allocator() const {
    return allocator_and_tag_.allocator();
  }
  allocator_type& allocator() {
    return allocator_and_tag_.allocator();
  }

  bool allocated() const { return tag().allocated(); }
  void set_allocated() { return tag().set_allocated(); }

  void set_size_internal(size_type n) {
    if (allocated()) {
      allocation().set_size(n);
    } else {
      tag().set_size(n);
    }
  }

  // Enlarge the underlying representation so we can store size_ + delta elems.
  // The size is not changed, and any newly added memory is not initialized.
  void EnlargeBy(size_type delta);

  void ResetAllocation(Allocation new_allocation) {
    if (allocated()) {
      DestroyAllocated(allocated_space(), allocated_space() + size());
      DCHECK_EQ(begin(), allocated_space());
      allocation().Dealloc(allocator());
      allocation() = new_allocation;
    } else {
      DestroyInlined(inlined_space(), inlined_space() + size());
      init_allocation(new_allocation);  // bug: only init once
      set_allocated();
    }
  }

  void GrowAndPushBack(const value_type& t) {
    DCHECK_EQ(size(), capacity());
    const size_type s = size();

    Allocation new_allocation(allocator(), 2 * capacity());
    new_allocation.set_size(s + 1);

    UninitializedCopyAllocated(array(), array() + s, new_allocation.buffer());
    ConstructAllocated(new_allocation.buffer() + s, t);

    ResetAllocation(new_allocation);
  }

  void InitAssign(size_type n);
  void InitAssign(size_type n, const value_type& t);

  void ConstructInlined(pointer p) {
    new(p) value_type();
  }

  void ConstructInlined(pointer p, const value_type& t) {
    new(p) value_type(t);
  }

  void ConstructAllocated(pointer p) {
    AllocatorTraits::construct(allocator(), p);
  }
  void ConstructAllocated(pointer p, const value_type& t) {
    AllocatorTraits::construct(allocator(), p, t);
  }

  template <typename Iter>
  void UninitializedCopyInlined(Iter src, Iter src_last, value_type* dst) {
    std::uninitialized_copy(src, src_last, dst);
  }

  template <typename Iter>
  void UninitializedCopyAllocated(Iter src, Iter src_last, value_type* dst) {
    for (; src != src_last; ++dst, ++src) ConstructAllocated(dst, *src);
  }

  void UninitializedFillInlined(value_type* dst, value_type* dst_last) {
    for (; dst != dst_last; ++dst) ConstructInlined(dst);
  }
  void UninitializedFillInlined(value_type* dst, value_type* dst_last,
                                const value_type& t) {
    std::uninitialized_fill(dst, dst_last, t);
  }

  void UninitializedFillAllocated(value_type* dst, value_type* dst_last) {
    for (; dst != dst_last; ++dst) ConstructAllocated(dst);
  }
  void UninitializedFillAllocated(value_type* dst, value_type* dst_last,
                                  const value_type& t) {
    for (; dst != dst_last; ++dst) ConstructAllocated(dst, t);
  }

  // Destroy [ptr, ptr_last) in place.
  void DestroyInlined(value_type* ptr, value_type* ptr_last);
  void DestroyAllocated(value_type* ptr, value_type* ptr_last);

  template <typename Iter>
  void AppendRange(Iter first, Iter last, std::input_iterator_tag) {
    std::copy(first, last, std::back_inserter(*this));
  }

  // Faster path for forward iterators.
  template <typename Iter>
  void AppendRange(Iter first, Iter last, std::forward_iterator_tag);

  template <typename Iter>
  void AppendRange(Iter first, Iter last) {
    typedef typename std::iterator_traits<Iter>::iterator_category IterTag;
    AppendRange(first, last, IterTag());
  }

  template <typename Iter>
  void AssignRange(Iter first, Iter last, std::input_iterator_tag);

  // Faster path for forward iterators.
  template <typename Iter>
  void AssignRange(Iter first, Iter last, std::forward_iterator_tag);

  template <typename Iter>
  void AssignRange(Iter first, Iter last) {
    typedef typename std::iterator_traits<Iter>::iterator_category IterTag;
    AssignRange(first, last, IterTag());
  }

  AllocatorAndTag allocator_and_tag_;

  // Either the inlined or allocated representation
  union Rep {
    // Use struct to perform indirection that solves a bizarre compilation
    // error on Visual Studio (all known versions).
    struct {
      google::ManualConstructor<value_type>
          inlined[N];
    } inlined_storage;
    struct {
      google::ManualConstructor<Allocation>
          allocation;
    } allocation_storage;
  } rep_;
};

template <typename T, int N, typename A>
const typename InlinedVector<T, N, A>::size_type
    InlinedVector<T, N, A>::Tag::kAllocated;

template <typename T, int N, typename A>
void swap(InlinedVector<T, N, A>& a, InlinedVector<T, N, A>& b) {
  a.swap(b);
}

template <typename T, int N, typename A>
bool operator==(const InlinedVector<T, N, A>& a,
                const InlinedVector<T, N, A>& b) {
  return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin());
}

template <typename T, int N, typename A>
bool operator!=(const InlinedVector<T, N, A>& a,
                const InlinedVector<T, N, A>& b) {
  return !(a == b);
}

template <typename T, int N, typename A>
bool operator<(const InlinedVector<T, N, A>& a,
               const InlinedVector<T, N, A>& b) {
  return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}

template <typename T, int N, typename A>
bool operator>(const InlinedVector<T, N, A>& a,
               const InlinedVector<T, N, A>& b) {
  return b < a;
}

template <typename T, int N, typename A>
bool operator<=(const InlinedVector<T, N, A>& a,
                const InlinedVector<T, N, A>& b) {
  return !(b < a);
}

template <typename T, int N, typename A>
bool operator>=(const InlinedVector<T, N, A>& a,
                const InlinedVector<T, N, A>& b) {
  return !(a < b);
}


// ========================================
// Implementation

template <typename T, int N, typename A>
InlinedVector<T, N, A>::InlinedVector(const InlinedVector& v)
    : allocator_and_tag_(v.allocator()) {
  reserve(v.size());
  if (allocated()) {
    UninitializedCopyAllocated(v.begin(), v.end(), allocated_space());
  } else {
    UninitializedCopyInlined(v.begin(), v.end(), inlined_space());
  }
  set_size_internal(v.size());
}

template <typename T, int N, typename A>
void InlinedVector<T, N, A>::InitAssign(size_type n, const value_type& t) {
  if (n > static_cast<size_type>(N)) {
    Allocation new_allocation(allocator(), n);
    init_allocation(new_allocation);
    set_allocated();
    UninitializedFillAllocated(allocated_space(),
                               allocated_space() + n, t);
  } else {
    UninitializedFillInlined(inlined_space(),
                             inlined_space() + n, t);
  }
  set_size_internal(n);
}

template <typename T, int N, typename A>
void InlinedVector<T, N, A>::InitAssign(size_type n) {
  if (n > static_cast<size_type>(N)) {
    Allocation new_allocation(allocator(), n);
    init_allocation(new_allocation);
    set_allocated();
    UninitializedFillAllocated(allocated_space(),
                               allocated_space() + n);
  } else {
    UninitializedFillInlined(inlined_space(),
                             inlined_space() + n);
  }
  set_size_internal(n);
}

template <typename T, int N, typename A>
void InlinedVector<T, N, A>::resize(size_type n) {
  size_type s = size();
  if (n < s) {
    erase(begin() + n, end());
    return;
  }
  reserve(n);
  DCHECK_GE(capacity(), n);

  // Fill new space with elements constructed in-place.
  if (allocated()) {
    UninitializedFillAllocated(allocated_space() + s,
                               allocated_space() + n);
  } else {
    UninitializedFillInlined(inlined_space() + s,
                             inlined_space() + n);
  }
  set_size_internal(n);
}

template <typename T, int N, typename A>
void InlinedVector<T, N, A>::resize(size_type n, const value_type& elem) {
  size_type s = size();
  if (n < s) {
    erase(begin() + n, end());
    return;
  }
  reserve(n);
  DCHECK_GE(capacity(), n);

  // Fill new space with copies of 'elem'.
  if (allocated()) {
    UninitializedFillAllocated(allocated_space() + s,
                               allocated_space() + n,
                               elem);
  } else {
    UninitializedFillInlined(inlined_space() + s,
                             inlined_space() + n,
                             elem);
  }
  set_size_internal(n);
}

template <typename T, int N, typename A>
typename InlinedVector<T, N, A>::iterator
InlinedVector<T, N, A>::insert(iterator pos, const value_type& v) {
  DCHECK_GE(pos, begin());
  DCHECK_LE(pos, end());
  if (pos == end()) {
    push_back(v);
    return end() - 1;
  }
  size_type s = size();
  size_type idx = std::distance(begin(), pos);
  if (s == capacity()) {
    EnlargeBy(1);
  }
  CHECK_LT(s, capacity());
  pos = begin() + idx;  // Reset 'pos' into a post-enlarge iterator.

  if (allocated()) {
    ConstructAllocated(allocated_space() + s,
                       *(allocated_space() + s - 1));
    std::copy_backward(pos,
                       allocated_space() + s - 1,
                       allocated_space() + s);
  } else {
    ConstructInlined(inlined_space() + s,
                     *(inlined_space() + s - 1));
    std::copy_backward(pos,
                       inlined_space() + s - 1,
                       inlined_space() + s);
  }

  *pos = v;

  set_size_internal(s + 1);
  return pos;
}

template <typename T, int N, typename A>
typename InlinedVector<T, N, A>::iterator
InlinedVector<T, N, A>::erase(iterator first, iterator last) {
  DCHECK_LE(begin(), first);
  DCHECK_LE(first, last);
  DCHECK_LE(last, end());

  size_type s = size();
  ptrdiff_t erase_gap = std::distance(first, last);

  if (allocated()) {
    std::copy(last, allocated_space() + s, first);
    DestroyAllocated(allocated_space() + s - erase_gap,
                     allocated_space() + s);
  } else {
    std::copy(last, inlined_space() + s, first);
    DestroyInlined(inlined_space() + s - erase_gap,
                   inlined_space() + s);
  }

  set_size_internal(size() - erase_gap);

  return first;
}

template <typename T, int N, typename A>
void InlinedVector<T, N, A>::swap(InlinedVector& other) {
  using std::swap;  // Augment ADL with std::swap.
  if (&other == this) {
    return;
  }
  if (allocated() && other.allocated()) {
    // Both out of line, so just swap the tag, allocation, and allocator.
    swap(tag(), other.tag());
    swap(allocation(), other.allocation());
    swap(allocator(), other.allocator());
    return;
  }
  if (!allocated() && !other.allocated()) {
    // Both inlined: swap up to smaller size, then move remaining elements.
    InlinedVector* a = this;
    InlinedVector* b = &other;
    if (size() < other.size()) {
      swap(a, b);
    }

    const size_type a_size = a->size();
    const size_type b_size = b->size();
    DCHECK_GE(a_size, b_size);
    // 'a' is larger. Swap the elements up to the smaller array size.
    std::swap_ranges(a->inlined_space(),
                     a->inlined_space() + b_size,
                     b->inlined_space());

    // Move the remaining elements: A[b_size,a_size) -> B[b_size,a_size)
    b->UninitializedCopyInlined(a->inlined_space() + b_size,
                                a->inlined_space() + a_size,
                                b->inlined_space() + b_size);
    a->DestroyInlined(a->inlined_space() + b_size,
                      a->inlined_space() + a_size);

    swap(a->tag(), b->tag());
    swap(a->allocator(), b->allocator());
    DCHECK_EQ(b->size(), a_size);
    DCHECK_EQ(a->size(), b_size);
    return;
  }
  // One is out of line, one is inline.
  // We first move the elements from the inlined vector into the
  // inlined space in the other vector.  We then put the other vector's
  // pointer/capacity into the originally inlined vector and swap
  // the tags.
  InlinedVector* a = this;
  InlinedVector* b = &other;
  if (a->allocated()) {
    swap(a, b);
  }
  DCHECK(!a->allocated());
  DCHECK(b->allocated());
  const size_type a_size = a->size();
  const size_type b_size = b->size();

  // Made Local copies of size(), don't need tag() accurate anymore
  swap(a->tag(), b->tag());

  // Copy b_allocation out before b's union gets clobbered by inline_space.
  Allocation b_allocation = b->allocation();

  b->UninitializedCopyInlined(a->inlined_space(),
                              a->inlined_space() + a_size,
                              b->inlined_space());
  a->DestroyInlined(a->inlined_space(),
                    a->inlined_space() + a_size);

  a->allocation() = b_allocation;

  if (a->allocator() != b->allocator()) {
    swap(a->allocator(), b->allocator());
  }

  DCHECK_EQ(b->size(), a_size);
  DCHECK_EQ(a->size(), b_size);
}

template <typename T, int N, typename A>
void InlinedVector<T, N, A>::EnlargeBy(size_type delta) {
  const size_type s = size();
  DCHECK_LE(s, capacity());

  size_type target = std::max(static_cast<size_type>(N), s + delta);

  // Compute new capacity by repeatedly doubling current capacity
  // TODO(user): Check and avoid overflow?
  size_type new_capacity = capacity();
  while (new_capacity < target) {
    new_capacity <<= 1;
  }

  Allocation new_allocation(allocator(), new_capacity);
  new_allocation.set_size(s);

  UninitializedCopyAllocated(array(), array() + s, new_allocation.buffer());

  ResetAllocation(new_allocation);
}

template <typename T, int N, typename A>
void InlinedVector<T, N, A>::DestroyInlined(value_type* ptr,
                                            value_type* ptr_last) {
  for (value_type* p = ptr; p != ptr_last; ++p) {
    p->~value_type();
  }

  // Overwrite unused memory with 0xab so we can catch uninitialized usage.
  // Cast to void* to tell the compiler that we don't care that we might be
  // scribbling on a vtable pointer.
#ifndef NDEBUG
  if (ptr != ptr_last) {
    memset(reinterpret_cast<void*>(ptr), 0xab,
           sizeof(*ptr) * (ptr_last - ptr));
  }
#endif
}

template <typename T, int N, typename A>
void InlinedVector<T, N, A>::DestroyAllocated(value_type* ptr,
                                              value_type* ptr_last) {
  for (value_type* p = ptr; p != ptr_last; ++p) {
    AllocatorTraits::destroy(allocator(), p);
  }

  // Overwrite unused memory with 0xab so we can catch uninitialized usage.
  // Cast to void* to tell the compiler that we don't care that we might be
  // scribbling on a vtable pointer.
#ifndef NDEBUG
  if (ptr != ptr_last) {
    memset(reinterpret_cast<void*>(ptr), 0xab,
           sizeof(*ptr) * (ptr_last - ptr));
  }
#endif
}

template <typename T, int N, typename A>
template <typename Iter>
void InlinedVector<T, N, A>::AppendRange(Iter first, Iter last,
                                         std::forward_iterator_tag) {
  typedef typename std::iterator_traits<Iter>::difference_type Length;
  Length length = std::distance(first, last);
  reserve(size() + length);
  if (allocated()) {
    UninitializedCopyAllocated(first, last, allocated_space() + size());
  } else {
    UninitializedCopyInlined(first, last, inlined_space() + size());
  }
  set_size_internal(size() + length);
}

template <typename T, int N, typename A>
template <typename Iter>
void InlinedVector<T, N, A>::AssignRange(Iter first, Iter last,
                                         std::input_iterator_tag) {
  // Optimized to avoid reallocation.
  // Prefer reassignment to copy construction for elements.
  iterator out = begin();
  for ( ; first != last && out != end(); ++first, ++out)
    *out = *first;
  erase(out, end());
  std::copy(first, last, std::back_inserter(*this));
}

template <typename T, int N, typename A>
template <typename Iter>
void InlinedVector<T, N, A>::AssignRange(Iter first, Iter last,
                                         std::forward_iterator_tag) {
  typedef typename std::iterator_traits<Iter>::difference_type Length;
  Length length = std::distance(first, last);
  // Prefer reassignment to copy construction for elements.
  if (length <= size()) {
    erase(std::copy(first, last, begin()), end());
    return;
  }
  reserve(length);
  iterator out = begin();
  for (; out != end(); ++first, ++out) *out = *first;
  if (allocated())
    UninitializedCopyAllocated(first, last, out);
  else
    UninitializedCopyInlined(first, last, out);
  set_size_internal(length);
}

GTL_NAMESPACE_END

#endif  // S2_UTIL_GTL_INLINED_VECTOR_H_
