// Copyright 2000 Google Inc. All Rights Reserved.
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

#include "s2/util/coding/coder.h"

#include <cassert>
#include <algorithm>

#include "s2/third_party/absl/base/integral_types.h"
#include <glog/logging.h>

// An initialization value used when we are allowed to
unsigned char Encoder::kEmptyBuffer = 0;

Encoder::Encoder()
  : underlying_buffer_(&kEmptyBuffer) {
}

Encoder::~Encoder() {
  if (underlying_buffer_ != &kEmptyBuffer) {
    delete[] underlying_buffer_;
  }
}

int Encoder::varint32_length(uint32 v) {
  return Varint::Length32(v);
}

int Encoder::varint64_length(uint64 v) {
  return Varint::Length64(v);
}

void Encoder::EnsureSlowPath(size_t N) {
  CHECK(ensure_allowed());
  assert(avail() < N);
  assert(length() == 0 || orig_ == underlying_buffer_);

  // Double buffer size, but make sure we always have at least N extra bytes
  const size_t current_len = length();
  const size_t new_capacity = std::max(current_len + N, 2 * current_len);

  unsigned char* new_buffer = new unsigned char[new_capacity];
  memcpy(new_buffer, underlying_buffer_, current_len);
  if (underlying_buffer_ != &kEmptyBuffer) {
    delete[] underlying_buffer_;
  }
  underlying_buffer_ = new_buffer;

  orig_ = new_buffer;
  limit_ = new_buffer + new_capacity;
  buf_ = orig_ + current_len;
  CHECK(avail() >= N);
}

namespace {
// Copies N bytes from *src to *dst then advances both pointers by N bytes.
// Template parameter N specifies the number of bytes to copy. Passing
// constant size results in optimized code from memcpy for the size.
template <size_t N>
void CopyAndAdvance(const uint8** src, uint8** dst) {
  memcpy(*dst, *src, N);
  *dst += N;
  *src += N;
}
}  // namespace

// Tries a fast path if both the decoder and the encoder have enough room for
// max varint64 (10 bytes). With enough room, we don't need boundary checks at
// every iterations. Also, memcpy with known size is faster than copying a byte
// at a time (e.g. one movq vs. eight movb's).
//
// If either the decoder or the encoder doesn't have enough room, it falls back
// to previous example where copy and boundary check happen at every byte.
bool Encoder::PutVarint64FromDecoderLessCommonSizes(Decoder* dec) {
  const unsigned char* dec_ptr = dec->buf_;
  const unsigned char* dec_limit = dec->limit_;

  // Check once if both the encoder and the decoder have enough room for
  // maximum varint64 (kVarintMax64) instead of checking at every bytes.
  if (ABSL_PREDICT_TRUE(dec_ptr <= dec_limit - kVarintMax64 &&
                        buf_ <= limit_ - kVarintMax64)) {
    if (dec_ptr[2] < 128) {
      CopyAndAdvance<3>(&dec->buf_, &buf_);
    } else if (dec_ptr[3] < 128) {
      CopyAndAdvance<4>(&dec->buf_, &buf_);
    } else if (dec_ptr[4] < 128) {
      CopyAndAdvance<5>(&dec->buf_, &buf_);
    } else if (dec_ptr[5] < 128) {
      CopyAndAdvance<6>(&dec->buf_, &buf_);
    } else if (dec_ptr[6] < 128) {
      CopyAndAdvance<7>(&dec->buf_, &buf_);
    } else if (dec_ptr[7] < 128) {
      CopyAndAdvance<8>(&dec->buf_, &buf_);
    } else if (dec_ptr[8] < 128) {
      CopyAndAdvance<9>(&dec->buf_, &buf_);
    } else if (dec_ptr[9] < 2) {
      // 10th byte stores at most 1 bit for varint64.
      CopyAndAdvance<10>(&dec->buf_, &buf_);
    } else {
      return false;
    }
    return true;
  }

  unsigned char c;
  unsigned char* enc_ptr = buf_;

  // The loop executes at most (kVarintMax64 - 1) iterations because either the
  // decoder or the encoder has less availability than kVarintMax64. We must be
  // careful about the cost of moving any computation out of the loop.
  // Xref cl/133546957 for more details of various implementations we explored.
  do {
    if (dec_ptr >= dec_limit)
      return false;
    if (enc_ptr >= limit_)
      return false;
    c = *dec_ptr;
    *enc_ptr = c;
    ++dec_ptr;
    ++enc_ptr;
  } while (c >= 128);

  dec->buf_ = dec_ptr;
  buf_ = enc_ptr;
  return true;
}

void Encoder::RemoveLast(size_t N) {
  CHECK(length() >= N);
  buf_ -= N;
}

void Encoder::Resize(size_t N) {
  CHECK(length() >= N);
  buf_ = orig_ + N;
  assert(length() == N);
}
