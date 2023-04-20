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


#ifndef S2_S2CODER_H_
#define S2_S2CODER_H_

#include <type_traits>

#include "s2/base/integral_types.h"
#include "absl/status/status.h"
#include "s2/util/coding/coder.h"
#include "s2/s2error.h"

// A general purpose encoding/decoding interface for S2 data types; This
// interface abstracts away differences in how classes marshal themselves to
// provide a simpler, more consistent API. Each value is serialized into a
// compact internally-restrained format, such that users of the coder do not
// need to independently include the length of variable-length values.
//
// Concrete types used with this API are required to be move-constructible and
// default-constructible.
//
// Note that, in general, encoded forms don't store type information, leaving
// correct tracking of types to the end user, although polymorphism is an
// internal feature of some S2Coders.

namespace s2coding {

// Controls whether to optimize for speed or size when encoding shapes.  (Note
// that encoding is always lossless, and that compact encodings are currently
// only possible when points have been snapped to S2CellId centers.)
enum class CodingHint : uint8 { FAST, COMPACT };

// S2Coder interface.
template <class T>
class S2Coder {
 public:
  virtual ~S2Coder() = default;

  // Encodes the value of T into the given encoder, including its own length
  // information if not implicit from the type.
  virtual void Encode(Encoder&, const T&) const = 0;

  // Decodes the next value of T from the decoder and leaves it positioned at
  // the first byte after the encoded representation of T.
  virtual bool Decode(Decoder&, T&, S2Error&) const = 0;

  absl::Status Decode(Decoder& decoder, T& v) const {
    S2Error error;
    v.Init(&decoder, error);
    return ToStatus(error);
  }

  // Disallow assignment to avoid slicing through references.
  S2Coder& operator=(const S2Coder&) = delete;
  S2Coder& operator=(S2Coder&&) = delete;
};

// Generic coder for any class with an Encode method that does not take a coding
// hint.
template <typename T>
class S2BasicCoder : public S2Coder<T> {
 public:
  void Encode(Encoder& encoder, const T& v) const override {
    v.Encode(&encoder);
  }

  bool Decode(Decoder& decoder, T& v, S2Error& error) const override {
    return v.Init(&decoder, error);
  }

  using S2Coder<T>::Decode;
};

// Generic coder for any class with an Encode method that takes a coding hint.
template <typename T>
class S2HintCoder : public S2Coder<T> {
 public:
  explicit S2HintCoder(CodingHint hint = CodingHint::FAST) : hint_(hint) {}

  void Encode(Encoder& encoder, const T& v) const override {
    v.Encode(&encoder, hint_);
  }

  bool Decode(Decoder& decoder, T& v, S2Error& error) const override {
    return v.Init(&decoder, error);
  }

  using S2Coder<T>::Decode;

 private:
  CodingHint hint_;
};

namespace internal {

// Coders for now legacy Decode style decoding.  Decode methods don't populate
// an S2Error so new types should instead favor the Init style of decoding.

// Generic coder for any class with an Encode method that does not take a coding
// hint, and a Decode method for decoding the type.
template <typename T>
class S2LegacyCoder : public S2Coder<T> {
 public:
  void Encode(Encoder& encoder, const T& v) const override {
    v.Encode(&encoder);
  }

  bool Decode(Decoder& decoder, T& v, S2Error& error) const override {
    if (!v.Decode(&decoder)) {
      error.Init(S2Error::DATA_LOSS, "Unknown decoding error");
      return false;
    }
    return true;
  }
};

// Generic coder for any class with an Encode method that takes a coding hint,
// and a Decode method for decoding the type.
template <typename T>
class S2LegacyHintCoder : public S2Coder<T> {
 public:
  explicit S2LegacyHintCoder(CodingHint hint = CodingHint::FAST)
      : hint_(hint) {}

  void Encode(Encoder& encoder, const T& v) const override {
    v.Encode(&encoder, hint_);
  }

  bool Decode(Decoder& decoder, T& v, S2Error& error) const override {
    if (!v.Decode(&decoder)) {
      error.Init(S2Error::DATA_LOSS, "Unknown decoding error");
      return false;
    }
    return true;
  }

 private:
  CodingHint hint_;
};

}  // namespace internal

}  // namespace s2coding

#endif  // S2_S2CODER_H_
