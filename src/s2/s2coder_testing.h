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


#ifndef S2_S2CODER_TESTING_H_
#define S2_S2CODER_TESTING_H_

#include "s2/util/coding/coder.h"
#include "s2/s2coder.h"
#include "s2/s2error.h"

namespace s2coding {

// Encodes a type T into a new Encoder instance.
template <typename T>
Encoder EncodeToEncoder(const S2Coder<T>& coder, const T& shape) {
  Encoder encoder;
  coder.Encode(encoder, shape);
  return encoder;
}

// Decodes a type T encoded into an Encoder instance.
template <typename T>
bool DecodeFromEncoded(const S2Coder<T>& coder, Encoder& encoder, T& value,
                       S2Error& error) {
  Decoder decoder(encoder.base(), encoder.length());
  return coder.Decode(decoder, value, error);
}

// Encodes then decodes a type T via S2Coder<T>, returning decoded object.
template <typename T>
T RoundTrip(const S2Coder<T>& coder, const T& shape, S2Error& error) {
  Encoder encoder = EncodeToEncoder(coder, shape);
  T val;
  S2_CHECK(DecodeFromEncoded(coder, encoder, val, error));
  return val;
}

}  // namespace s2coding

#endif  // S2_S2CODER_TESTING_H_
