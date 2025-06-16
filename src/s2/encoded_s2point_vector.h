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

#ifndef S2_ENCODED_S2POINT_VECTOR_H_
#define S2_ENCODED_S2POINT_VECTOR_H_

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "absl/base/attributes.h"
#include "absl/base/optimization.h"
#include "absl/log/absl_log.h"
#include "absl/types/span.h"
#include "s2/util/coding/coder.h"
#include "s2/encoded_string_vector.h"
#include "s2/encoded_uint_vector.h"
#include "s2/s2coder.h"
#include "s2/s2error.h"
#include "s2/s2point.h"

namespace s2coding {

// Encodes a vector of S2Points in a format that can later be decoded as an
// EncodedS2PointVector.
//
// REQUIRES: "encoder" uses the default constructor, so that its buffer
//           can be enlarged as necessary by calling Ensure(int).
void EncodeS2PointVector(absl::Span<const S2Point> points, CodingHint hint,
                         Encoder* encoder);

// This class represents an encoded vector of S2Points.  Values are decoded
// only when they are accessed.  This allows for very fast initialization and
// no additional memory use beyond the encoded data.  The encoded data is not
// owned by this class; typically it points into a large contiguous buffer
// that contains other encoded data as well.
//
// This is one of several helper classes that allow complex data structures to
// be initialized from an encoded format in constant time and then decoded on
// demand.  This can be a big performance advantage when only a small part of
// the data structure is actually used.
class EncodedS2PointVector {
 public:
  // Constructs an uninitialized object; requires Init() to be called.
  EncodedS2PointVector() {}

  // Initializes the EncodedS2PointVector.
  //
  // Encoded types do not necessarily read all of their data when calling
  // Init(), instead only data necessary to decode further is read, which means
  // that encoded types may in general log ERROR errors when reading from an
  // untrusted byte stream using unchecked accessors (e.g. operator[]).
  //
  // Checked accessors such as Decode(S2Error&) will instead set an error
  // condition, at the cost of potentially more overhead when decoding.
  //
  // REQUIRES: The Decoder data buffer must outlive this object.
  ABSL_MUST_USE_RESULT bool Init(Decoder* decoder);

  // Initializes the EncodedS2PointVector as above.  Setting S2Error if
  // initialization fails.  A true return status does -not- mean that the byte
  // stream is without errors.  The individual points must still be checked with
  // a checked accessor such as Decode(S2Error&).
  ABSL_MUST_USE_RESULT bool Init(Decoder* decoder, S2Error& error);

  // Returns the size of the original vector.
  size_t size() const;

  // Returns the element at the given index.
  S2Point operator[](int i) const;

  // Decodes and returns the entire original vector.
  std::vector<S2Point> Decode() const;

  // Decodes the entire original vector and stores it in the given span.
  bool Decode(absl::Span<S2Point> points) const;

  // Decodes and returns the entire original vector, checking for errors as it
  // goes.  If an error occurs, a truncated vector is returned.
  std::vector<S2Point> Decode(S2Error& error) const;

  // Decodes the entire original vector and stores it in the given span.  Checks
  // for errors as it goes.  If an error occurs returns false and the span may
  // be incompletely written.
  bool Decode(absl::Span<S2Point> points, S2Error& error) const;

  // Copy the encoded data to the encoder. This allows for "reserialization" of
  // encoded shapes created through lazy decoding.
  void Encode(Encoder* encoder) const;

  // TODO(ericv): Consider adding a method that returns an adjacent pair of
  // points.  This would save some decoding overhead.

 private:
  friend void EncodeS2PointVector(absl::Span<const S2Point>, CodingHint,
                                  Encoder*);
  friend void EncodeS2PointVectorFast(absl::Span<const S2Point>, Encoder*);
  friend void EncodeS2PointVectorCompact(absl::Span<const S2Point>, Encoder*);

  // Decodes and returns the point at the given index.  If any errors occur when
  // decoding, the error is set and a default constructed point returned.
  S2Point At(int i, S2Error& error) const;

  bool InitUncompressedFormat(Decoder* decoder);
  bool InitCellIdsFormat(Decoder* decoder);
  S2Point DecodeCellIdsFormat(int i, S2Error* error) const;

  // We use a tagged union to represent multiple formats, as opposed to an
  // abstract base class or templating.  This represents the best compromise
  // between performance, space, and convenience.  Note that the overhead of
  // checking the tag is trivial and will typically be branch-predicted
  // perfectly.
  //
  // TODO(ericv): Once additional formats have been implemented, consider
  // using std::variant<> instead.  It's unclear whether this would have
  // better or worse performance than the current approach.
  enum Format : uint8_t {
    UNCOMPRESSED = 0,
    CELL_IDS = 1,
  };
  Format format_;
  uint32_t size_;
  union {
    struct {
      const S2Point* points;
    } uncompressed_;
    struct {
      EncodedStringVector blocks;
      uint64_t base;
      uint8_t level;
      bool have_exceptions;

      // TODO(ericv): Use std::atomic_flag to cache the last point decoded in
      // a thread-safe way.  This reduces benchmark times for actual polygon
      // operations (e.g. S2ClosestEdgeQuery) by about 15%.
    } cell_ids_;
  };
};


//////////////////   Implementation details follow   ////////////////////


inline size_t EncodedS2PointVector::size() const {
  return size_;
}

inline S2Point EncodedS2PointVector::At(int i, S2Error& error) const {
  switch (format_) {
    case Format::UNCOMPRESSED:
      return uncompressed_.points[i];

    case Format::CELL_IDS:
      return DecodeCellIdsFormat(i, &error);

    default:
      error = S2Error::DataLoss("Unrecognized format");
      return {};
  }

  ABSL_UNREACHABLE();
}

inline S2Point EncodedS2PointVector::operator[](int i) const {
  // Check error conditions in debug mode, elide them otherwise.
#ifndef NDEBUG
  S2Error error;
  S2Point point = At(i, error);
  if (!error.ok()) {
    ABSL_LOG(ERROR) << error.message();
  }
  return point;
#else
  switch (format_) {
    case Format::UNCOMPRESSED:
      return uncompressed_.points[i];

    case Format::CELL_IDS:
      return DecodeCellIdsFormat(i, nullptr);

    default:
      return {};
  }

  ABSL_UNREACHABLE();
#endif
}

}  // namespace s2coding

#endif  // S2_ENCODED_S2POINT_VECTOR_H_
