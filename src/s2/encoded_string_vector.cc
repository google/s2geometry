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

#include "s2/encoded_string_vector.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "s2/util/coding/coder.h"
#include "s2/encoded_uint_vector.h"

using absl::MakeSpan;
using absl::Span;
using absl::string_view;
using std::string;
using std::vector;

namespace s2coding {

StringVectorEncoder::StringVectorEncoder() = default;

void StringVectorEncoder::Encode(Encoder* encoder) {
  offsets_.push_back(data_.length());
  // We don't encode the first element of "offsets_", which is always zero.
  EncodeUintVector<uint64_t>(
      MakeSpan(offsets_.data() + 1, offsets_.data() + offsets_.size()),
      encoder);
  encoder->Ensure(data_.length());
  encoder->putn(data_.base(), data_.length());
}

void StringVectorEncoder::Encode(Span<const string> v, Encoder* encoder) {
  StringVectorEncoder string_vector;
  for (const auto& str : v) string_vector.Add(str);
  string_vector.Encode(encoder);
}

bool EncodedStringVector::Init(Decoder* decoder) {
  if (!offsets_.Init(decoder)) return false;
  data_ = decoder->skip(0);
  uint64_t length = 0;
  for (int i = 0, n = offsets_.size(); i < n; ++i) {
    // Strings are packed sequentially, offsets_[i] representing end of item i.
    // String length is diff of two sequential offsets, thus the offsets should
    // not decrease, otherwise we might read memory beyond the end of data.
    if (offsets_[i] < length) return false;
    length = offsets_[i];
  }
  if (decoder->avail() < length) return false;
  decoder->skip(length);
  return true;
}

vector<string_view> EncodedStringVector::Decode() const {
  size_t n = size();
  vector<string_view> result(n);
  for (size_t i = 0; i < n; ++i) {
    result[i] = (*this)[i];
  }
  return result;
}

// The encoding must be identical to StringVectorEncoder::Encode().
void EncodedStringVector::Encode(Encoder* encoder) const {
  offsets_.Encode(encoder);

  if (offsets_.size() > 0) {
    const uint64_t length = offsets_[offsets_.size() - 1];
    encoder->Ensure(length);
    encoder->putn(data_, length);
  }
}

}  // namespace s2coding
