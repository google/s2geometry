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

#ifndef S2_S2DENSITY_TREE_INTERNAL_H_
#define S2_S2DENSITY_TREE_INTERNAL_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <string>

#include "absl/log/absl_check.h"
#include "absl/strings/string_view.h"
#include "s2/util/coding/varint.h"
#include "s2/s2cell_id.h"

// ReversibleBytes is a simple wrapper around a string of bytes that we can
// append to and get the reversed output of.
class ReversibleBytes {
 public:
  ReversibleBytes() = default;

  void Clear() { bytes_.clear(); }

  void AppendBytes(absl::string_view bytes) {
    // Can't use `append(bytes)` since `absl::string_view` might not be
    // `std::string_view` -- `std::string` won't know about
    // `absl::string_view`.
    bytes_.append(bytes.data(), bytes.size());
  }
  void AppendVarint64(uint64_t v) { Varint::Append64(&bytes_, v); }

  void ReverseFrom(size_t start) {
    std::reverse(bytes_.begin() + start, bytes_.end());
  }

  size_t size() const { return bytes_.size(); }

  std::string Reversed() const {
    return std::string(bytes_.rbegin(), bytes_.rend());
  }

 private:
  std::string bytes_;
};

// ReversedCellEncoder records the offsets of each in a sequence of S2Density
// Cells written by the caller in reverse order, and then writes a mask followed
// by the offset of each encoded cell in reverse order.  This allows
// unpredictable varints to be written at the start of each parent without
// having a more complicated scheme than simply reversing all the bytes after
// encoding.
class ReversedCellEncoder {
 public:
  explicit ReversedCellEncoder(ReversibleBytes* output)
      : output_(output), start_(output_->size()) {}

  // Records the size of an encoded child cell.  The caller should encode
  // child cells in reverse order.
  void Next() {
    ABSL_DCHECK_LT(size_, lengths_.size());
    lengths_[size_++] = output_->size() - start_;
    start_ = output_->size();
  }

  // Writes the mask and lengths in normal order and then reverses them, where
  // the last length is not included because the position of whatever comes
  // after can be inferred.
  void Finish(uint64_t v) {
    output_->AppendVarint64(v);

    // Since the child cells were encoded in reverse order, if there are n
    // children present then this loop encodes child cell lengths [0..n-2].
    for (int i = size_ - 1; i > 0; --i) {
      output_->AppendVarint64(lengths_[i]);
    }

    // The child cells are already encoded in reverse order and with reversed
    // bytes, so all that's needed to encode this cell in reverse is to reverse
    // the header we just encoded.
    output_->ReverseFrom(start_);
    start_ = output_->size();
  }

 private:
  ReversibleBytes* output_;
  std::array<uint64_t, S2CellId::kNumFaces> lengths_;
  size_t size_ = 0;
  size_t start_;
};

#endif  // S2_S2DENSITY_TREE_INTERNAL_H_
