// Copyright 2012 Google Inc. All Rights Reserved.
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
// DO NOT INCLUDE THIS FILE. It contains implementation details of
// scoped_ptr and related classes; its contents, as well as the file
// itself, should be considered implementation details of scoped_ptr.h

#ifndef BASE_SCOPED_PTR_INTERNALS_H_
#define BASE_SCOPED_PTR_INTERNALS_H_

#include <algorithm>  // for std::swap

namespace base {
namespace internal {


// Minimal implementation of the core logic of scoped_ptr, suitable for
// reuse in both scoped_ptr and its specialization.
template <class C, class D>
class scoped_ptr_impl {
 public:
  explicit scoped_ptr_impl(C* p) : data_(p) { }

  ~scoped_ptr_impl() {
    if (data_.ptr != NULL) {
      (static_cast<D&>(data_))(data_.ptr);
    }
  }

  void reset(C* p) {
    if (data_.ptr != NULL) {
      // Note that this can lead to undefined behavior and memory leaks
      // in the unlikely but possible case that get_deleter()(get())
      // indirectly deletes this. The fix is to reset ptr_ before deleting
      // its old value, but first we need to clean up the code that relies
      // on the current sequencing.
      (static_cast<D&>(data_))(data_.ptr);
    }
    data_.ptr = p;
  }

  C* get() const { return data_.ptr; }

  void swap(scoped_ptr_impl& p2) {
    // Standard swap idiom: 'using std::swap' ensures that std::swap is
    // present in the overload set, but we call swap unqualified so that
    // any more-specific overloads can be used, if available.
    using std::swap;
    swap(static_cast<D&>(data_), static_cast<D&>(p2.data_));
    swap(data_.ptr, p2.data_.ptr);
  }

  C* release() {
    C* retVal = data_.ptr;
    data_.ptr = NULL;
    return retVal;
  }

 private:
  // Use the empty base class optimization to allow us to have a D member,
  // while avoiding any space overhead for it when D is an empty class.
  // See e.g. http://www.cantrip.org/emptyopt.html for a good discussion of
  // this technique.
  struct Data : public D {
    explicit Data(C* ptr_in) : ptr(ptr_in) {}

    C* ptr;
  };

  Data data_;

  // Disallow copy and assignment.
  scoped_ptr_impl(const scoped_ptr_impl&);
  scoped_ptr_impl& operator=(const scoped_ptr_impl&);
};

}  // namespace internal
}  // namespace base

#endif  // BASE_SCOPED_PTR_INTERNALS_H_