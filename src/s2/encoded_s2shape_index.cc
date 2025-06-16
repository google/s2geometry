// Copyright 2018 Google Inc. All Rights Reserved.
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

// Author: ericv@google.com (Eric Veach)

#include "s2/encoded_s2shape_index.h"

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include "s2/base/casts.h"
#include "absl/base/optimization.h"
#include "absl/numeric/bits.h"
#include "s2/util/coding/coder.h"
#include "s2/util/coding/varint.h"
#include "s2/encoded_s2cell_id_vector.h"
#include "s2/encoded_string_vector.h"
#include "s2/mutable_s2shape_index.h"
#include "s2/s2cell_id.h"
#include "s2/s2point.h"
#include "s2/s2shape.h"
#include "s2/s2shape_index.h"

using std::make_unique;
using std::unique_ptr;
using std::vector;

S2Shape* EncodedS2ShapeIndex::GetShape(int id) const {
  // This method is called when a shape has not been decoded yet.
  unique_ptr<S2Shape> shape = (*shape_factory_)[id];
  S2Shape* expected = kUndecodedShape();
  if (shapes_[id].compare_exchange_strong(expected, shape.get(),
                                          std::memory_order_acq_rel)) {
    return shape.release();  // Ownership has been transferred to shapes_.
  }
  return expected;  // Another thread updated shapes_[id] first.
}

const S2ShapeIndexCell* EncodedS2ShapeIndex::GetCell(int i) const {
  // memory_order_release ensures that no reads or writes in the current
  // thread can be reordered after this store, and all writes in the current
  // thread are visible to other threads that acquire the same atomic
  // variable.
  //
  // memory_order_acquire ensures that no reads or writes in the current
  // thread can be reordered before this load, and all writes in other threads
  // that release the same atomic variable are visible in this thread.
  //
  // We use this to implement lock-free synchronization on the read path as
  // follows:
  //
  //  1. cells_decoded(i) is updated using acquire/release semantics
  //  2. cells_[i] is written before cells_decoded(i)
  //  3. cells_[i] is read after cells_decoded(i)
  //
  // Note that we do still use a lock for the write path to ensure that
  // cells_[i] and cell_decoded(i) are updated together atomically.
  if (cell_decoded(i)) return cells_[i];

  // Decode the cell before acquiring the spinlock in order to minimize the
  // time that the lock is held.
  auto cell = make_unique<S2ShapeIndexCell>();
  Decoder decoder = encoded_cells_.GetDecoder(i);
  if (!cell->Decode(num_shape_ids(), &decoder)) {
    return nullptr;
  }
  // Recheck cell_decoded(i) once we hold the lock in case another thread
  // has decoded this cell in the meantime.
  SpinLockHolder l(&cells_lock_);
  if (cell_decoded(i)) return cells_[i];

  // Update the cell, setting cells_[i] before cell_decoded(i).
  cells_[i] = cell.get();
  set_cell_decoded(i);
  if (cell_cache_.size() < static_cast<size_t>(max_cell_cache_size())) {
    cell_cache_.push_back(i);
  }
  return cell.release();  // Ownership has been transferred to cells_.
}

EncodedS2ShapeIndex::EncodedS2ShapeIndex() = default;

EncodedS2ShapeIndex::~EncodedS2ShapeIndex() {
  // Although Minimize() does slightly more than required for destruction
  // (i.e., it resets vector elements to their default values), this does not
  // affect benchmark times.
  Minimize();
}

bool EncodedS2ShapeIndex::Init(Decoder* decoder,
                               const ShapeFactory& shape_factory) {
  Minimize();
  uint64_t max_edges_version;
  if (!decoder->get_varint64(&max_edges_version)) return false;
  version_ = max_edges_version & 3;
  if (version_ != MutableS2ShapeIndex::kCurrentEncodingVersionNumber) {
    return false;
  }
  options_.set_max_edges_per_cell(max_edges_version >> 2);

  // AtomicShape is a subtype of std::atomic<S2Shape*> that changes the
  // default constructor value to kUndecodedShape().  This saves the effort of
  // initializing all the elements twice.
  shapes_ = vector<AtomicShape>(shape_factory.size());
  shape_factory_ = shape_factory.Clone();
  if (!cell_ids_.Init(decoder)) return false;

  // The cells_ elements are *uninitialized memory*.  Instead we have bit
  // vector (cells_decoded_) to indicate which elements of cells_ are valid.
  // This reduces constructor times by more than a factor of 50, since rather
  // than needing to initialize one 64-bit pointer per cell to zero, we only
  // need to initialize one bit per cell to zero.
  //
  // For very large S2ShapeIndexes the internal memset() call to initialize
  // cells_decoded_ still takes about 1.3 microseconds per million cells
  // (assuming an optimized implementation that writes 32 bytes per cycle),
  // but this seems reasonable relative to other likely costs (I/O, etc).
  //
  // NOTE(ericv): DO NOT use make_unique<> here!  make_unique<> allocates memory
  // using "new T[n]()", which initializes all elements of the array.  This
  // slows down some benchmarks by over 100x.
  //
  // cells_ = make_unique<S2ShapeIndexCell*>[](cell_ids_.size());
  // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  //                                NO NO NO
  cells_.reset(new S2ShapeIndexCell*[cell_ids_.size()]);
  cells_decoded_ = vector<std::atomic<uint64_t>>((cell_ids_.size() + 63) >> 6);

  return encoded_cells_.Init(decoder);
}

void EncodedS2ShapeIndex::Encode(Encoder* encoder) const {
  // Re-encode the max edges and version number.
  encoder->Ensure(Varint::kMax64);
  uint64_t max_edges = options_.max_edges_per_cell();
  encoder->put_varint64(max_edges << 2 | version_);

  // And copy the encoded cell ids and cells.
  cell_ids_.Encode(encoder);
  encoded_cells_.Encode(encoder);
}

void EncodedS2ShapeIndex::Minimize() {
  if (cells_ == nullptr) return;  // Not initialized yet.

  for (auto& atomic_shape : shapes_) {
    S2Shape* shape = atomic_shape.load(std::memory_order_relaxed);
    if (shape != kUndecodedShape() && shape != nullptr) {
      atomic_shape.store(kUndecodedShape(), std::memory_order_relaxed);
      delete shape;
    }
  }
  if (cell_cache_.size() < static_cast<size_t>(max_cell_cache_size())) {
    // When only a tiny fraction of the cells are decoded, we keep track of
    // those cells in cell_cache_ to avoid the cost of scanning the
    // cells_decoded_ vector.  (The cost is only about 1 cycle per 64 cells,
    // but for a huge polygon with 1 million cells that's still 16000 cycles.)
    for (int pos : cell_cache_) {
      cells_decoded_[pos >> 6].store(0, std::memory_order_relaxed);
      delete cells_[pos];
    }
  } else {
    // Scan the cells_decoded_ vector looking for cells that must be deleted.
    for (int i = cells_decoded_.size(); --i >= 0;) {
      uint64_t bits = cells_decoded_[i].load(std::memory_order_relaxed);
      if (bits == 0) continue;
      do {
        ABSL_ASSUME(bits != 0);
        int offset = absl::countr_zero(bits);
        delete cells_[(i << 6) + offset];
        bits &= bits - 1;
      } while (bits != 0);
      cells_decoded_[i].store(0, std::memory_order_relaxed);
    }
  }
  cell_cache_.clear();
}

size_t EncodedS2ShapeIndex::SpaceUsed() const {
  // TODO(ericv): Add SpaceUsed() method to S2Shape base class,and include
  // memory owned by the allocated S2Shapes (here and in S2ShapeIndex).
  size_t size = sizeof(*this);
  size += shapes_.capacity() * sizeof(std::atomic<S2Shape*>);
  size += cell_ids_.size() * sizeof(std::atomic<S2ShapeIndexCell*>);  // cells_
  size += cells_decoded_.capacity() * sizeof(std::atomic<uint64_t>);
  size += cell_cache_.capacity() * sizeof(int);
  return size;
}
