// Copyright 2013 Google Inc. All Rights Reserved.
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

#include "s2/s2lax_polygon_shape.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <memory>
#include <utility>
#include <vector>

#include "absl/log/absl_check.h"
#include "absl/memory/memory.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"

#include "s2/encoded_s2point_vector.h"
#include "s2/encoded_uint_vector.h"
#include "s2/s2coder.h"
#include "s2/s2error.h"
#include "s2/s2loop.h"
#include "s2/s2point.h"
#include "s2/s2polygon.h"
#include "s2/s2shape.h"
#include "s2/s2shapeutil_get_reference_point.h"
#include "s2/util/coding/coder.h"
#include "s2/util/coding/varint.h"

using absl::MakeSpan;
using absl::Span;
using std::make_unique;
using std::unique_ptr;
using std::vector;
using ChainPosition = S2Shape::ChainPosition;

namespace {
template <typename T>
unique_ptr<T> make_unique_for_overwrite(size_t n) {
  // We only need to support this one variant.
  static_assert(std::is_array<T>::value, "T must be an array type");
  return unique_ptr<T>(new typename absl::remove_extent_t<T>[n]);
}
}  // namespace


// When adding a new encoding, be aware that old binaries will not be able
// to decode it.
static const unsigned char kCurrentEncodingVersionNumber = 1;

S2LaxPolygonShape::S2LaxPolygonShape(
    absl::Span<const S2LaxPolygonShape::Loop> loops) {
  Init(loops);
}

S2LaxPolygonShape::S2LaxPolygonShape(Span<const Span<const S2Point>> loops) {
  Init(loops);
}

S2LaxPolygonShape::S2LaxPolygonShape(const S2Polygon& polygon) {
  Init(polygon);
}

S2LaxPolygonShape::S2LaxPolygonShape(S2LaxPolygonShape&& b) noexcept
    : num_loops_(std::exchange(b.num_loops_, 0)),
      prev_loop_(b.prev_loop_.exchange(0, std::memory_order_relaxed)),
      num_vertices_(std::exchange(b.num_vertices_, 0)),
      vertices_(std::move(b.vertices_)),
      loop_starts_(std::move(b.loop_starts_)) {}

S2LaxPolygonShape& S2LaxPolygonShape::operator=(
    S2LaxPolygonShape&& b) noexcept {
  using std::memory_order_relaxed;

  num_loops_ = std::exchange(b.num_loops_, 0);
  prev_loop_.store(b.prev_loop_.exchange(0, memory_order_relaxed),
                   memory_order_relaxed);
  num_vertices_ = std::exchange(b.num_vertices_, 0);
  vertices_ = std::move(b.vertices_);
  loop_starts_ = std::move(b.loop_starts_);
  return *this;
}

void S2LaxPolygonShape::Init(absl::Span<const S2LaxPolygonShape::Loop> loops) {
  vector<Span<const S2Point>> spans;
  spans.reserve(loops.size());
  for (const S2LaxPolygonShape::Loop& loop : loops) {
    spans.emplace_back(loop);
  }
  Init(spans);
}

void S2LaxPolygonShape::Init(const S2Polygon& polygon) {
  vector<Span<const S2Point>> spans;
  for (int i = 0; i < polygon.num_loops(); ++i) {
    const S2Loop* loop = polygon.loop(i);
    if (loop->is_full()) {
      spans.emplace_back();  // Empty span.
    } else {
      spans.emplace_back(&loop->vertex(0), loop->num_vertices());
    }
  }
  Init(spans);

  // S2Polygon and S2LaxPolygonShape holes are oriented oppositely, so we need
  // to reverse the orientation of any loops representing holes.
  for (int i = 0; i < polygon.num_loops(); ++i) {
    if (polygon.loop(i)->is_hole()) {
      S2Point* v0 = &vertices_[loop_starts_[i]];
      std::reverse(v0, v0 + num_loop_vertices(i));
    }
  }
}

void S2LaxPolygonShape::Init(Span<const Span<const S2Point>> loops) {
  num_loops_ = loops.size();
  if (num_loops_ == 0) {
    num_vertices_ = 0;
  } else if (num_loops_ == 1) {
    num_vertices_ = loops[0].size();
    // TODO(ericv): Use std::allocator to obtain uninitialized memory instead.
    // This would avoid default-constructing all the elements before we
    // overwrite them, and it would also save 8 bytes of memory allocation
    // since "new T[]" stores its own copy of the array size.
    //
    // Note that even absl::make_unique_for_overwrite<> and c++20's
    // std::make_unique_for_overwrite<T[]> default-construct all elements when
    // T is a class type.
    vertices_ = make_unique<S2Point[]>(num_vertices_);
    std::copy(loops[0].begin(), loops[0].end(), vertices_.get());
  } else {
    // Don't use make_unique<> here in order to avoid zero initialization.
    loop_starts_ = make_unique_for_overwrite<uint32_t[]>(num_loops_ + 1);
    num_vertices_ = 0;
    for (int i = 0; i < num_loops_; ++i) {
      loop_starts_[i] = num_vertices_;
      num_vertices_ += loops[i].size();
    }
    loop_starts_[num_loops_] = num_vertices_;
    vertices_ = make_unique<S2Point[]>(num_vertices_);  // TODO(see above)
    for (int i = 0; i < num_loops_; ++i) {
      std::copy(loops[i].begin(), loops[i].end(),
                vertices_.get() + loop_starts_[i]);
    }
  }
}

int S2LaxPolygonShape::num_loop_vertices(int i) const {
  ABSL_DCHECK_LT(i, num_loops());
  if (num_loops() == 1) {
    return num_vertices_;
  } else {
    return loop_starts_[i + 1] - loop_starts_[i];
  }
}

const S2Point& S2LaxPolygonShape::loop_vertex(int i, int j) const {
  ABSL_DCHECK_LT(i, num_loops());
  ABSL_DCHECK_LT(j, num_loop_vertices(i));
  if (i == 0) {
    return vertices_[j];
  } else {
    return vertices_[loop_starts_[i] + j];
  }
}

void S2LaxPolygonShape::Encode(Encoder* encoder,
                               s2coding::CodingHint hint) const {
  encoder->Ensure(1 + Varint::kMax32);
  encoder->put8(kCurrentEncodingVersionNumber);
  encoder->put_varint32(num_loops());
  s2coding::EncodeS2PointVector(MakeSpan(vertices_.get(), num_vertices()),
                                hint, encoder);
  if (num_loops() > 1) {
    s2coding::EncodeUintVector<uint32_t>(
        MakeSpan(loop_starts_.get(), num_loops() + 1), encoder);
  }
}

bool S2LaxPolygonShape::Init(Decoder* decoder, S2Error* error) {
  const auto Error = [&error](absl::string_view message) {
    if (error != nullptr) {
      *error = S2Error::DataLoss(message);
    }
    return false;
  };

  if (decoder->avail() < 1) {
    return Error("Insufficient data in decoder");
  }

  uint8_t version = decoder->get8();
  if (version != kCurrentEncodingVersionNumber) {
    return Error("Bad version number in byte string");
  }

  uint32_t num_loops;
  if (!decoder->get_varint32(&num_loops)) {
    return Error("Failed to decode number of loops");
  }
  // `loop_starts_` is indexed by an integer, and contains an extra entry,
  // so we limit num_loops to INT32_MAX - 1.
  if (num_loops > std::numeric_limits<int32_t>::max() - 1u) {
    return Error("Number of loops too large");
  }

  num_loops_ = num_loops;
  s2coding::EncodedS2PointVector vertices;
  if (!vertices.Init(decoder)) {
    return Error("Failed to decode vertices");
  }

  if (num_loops_ == 0) {
    num_vertices_ = 0;
  } else {
    num_vertices_ = vertices.size();
    vertices_ = make_unique<S2Point[]>(num_vertices_);  // TODO(see above)

    // Load the polygon vertices from the encoded s2point vector.
    if (error == nullptr) {
      vertices.Decode(absl::MakeSpan(vertices_.get(), num_vertices_));
    } else {
      vertices.Decode(absl::MakeSpan(vertices_.get(), num_vertices_), *error);
      if (!error->ok()) {
        return false;
      }
    }

    if (num_loops_ > 1) {
      s2coding::EncodedUintVector<uint32_t> loop_starts;
      if (!loop_starts.Init(decoder) || loop_starts.size() != num_loops_ + 1) {
        return Error("Failed to decode loop offsets");
      }

      loop_starts_ = make_unique_for_overwrite<uint32_t[]>(loop_starts.size());
      for (size_t i = 0; i < loop_starts.size(); ++i) {
        loop_starts_[i] = loop_starts[i];
      }
    }
  }

  return true;
}

S2Shape::Edge S2LaxPolygonShape::edge(int e) const {
  // Method names are fully specified to enable inlining.
  ChainPosition pos = S2LaxPolygonShape::chain_position(e);
  return S2LaxPolygonShape::chain_edge(pos.chain_id, pos.offset);
}

S2Shape::ReferencePoint S2LaxPolygonShape::GetReferencePoint() const {
  return s2shapeutil::GetReferencePoint(*this);
}

S2Shape::Chain S2LaxPolygonShape::chain(int i) const {
  ABSL_DCHECK_LT(i, num_loops());
  if (num_loops() == 1) {
    return Chain(0, num_vertices_);
  } else {
    int start = loop_starts_[i];
    return Chain(start, loop_starts_[i + 1] - start);
  }
}

EncodedS2LaxPolygonShape::EncodedS2LaxPolygonShape(
    EncodedS2LaxPolygonShape&& b) noexcept
    : num_loops_(std::exchange(b.num_loops_, 0)),
      prev_loop_(b.prev_loop_.exchange(0, std::memory_order_relaxed)),
      vertices_(std::move(b.vertices_)),
      loop_starts_(std::move(b.loop_starts_)) {}

EncodedS2LaxPolygonShape& EncodedS2LaxPolygonShape::operator=(
    EncodedS2LaxPolygonShape&& b) noexcept {
  num_loops_ = std::exchange(b.num_loops_, 0);
  prev_loop_.store(b.prev_loop_.exchange(0, std::memory_order_relaxed),
                   std::memory_order_relaxed);
  vertices_ = std::move(b.vertices_);
  loop_starts_ = std::move(b.loop_starts_);
  return *this;
}

bool EncodedS2LaxPolygonShape::Init(Decoder* decoder) {
  if (decoder->avail() < 1) return false;
  uint8_t version = decoder->get8();
  if (version != kCurrentEncodingVersionNumber) return false;

  uint32_t num_loops;
  if (!decoder->get_varint32(&num_loops)) return false;
  // `loop_starts_` is indexed by an integer, and contains an extra entry,
  // so we limit num_loops to INT32_MAX - 1.
  if (num_loops > std::numeric_limits<int32_t>::max() - 1u) return false;
  num_loops_ = num_loops;

  if (!vertices_.Init(decoder)) return false;

  if (num_loops_ > 1) {
    if (!loop_starts_.Init(decoder)) return false;
  }
  return true;
}

// The encoding must be identical to S2LaxPolygonShape::Encode().
void EncodedS2LaxPolygonShape::Encode(Encoder* encoder,
                                      s2coding::CodingHint) const {
  encoder->Ensure(1 + Varint::kMax32);
  encoder->put8(kCurrentEncodingVersionNumber);
  encoder->put_varint32(num_loops_);
  vertices_.Encode(encoder);
  if (num_loops_ > 1) {
    loop_starts_.Encode(encoder);
  }
}

int EncodedS2LaxPolygonShape::num_vertices() const {
  if (num_loops() <= 1) {
    return vertices_.size();
  } else {
    return loop_starts_[num_loops()];
  }
}

int EncodedS2LaxPolygonShape::num_loop_vertices(int i) const {
  ABSL_DCHECK_LT(i, num_loops());
  if (num_loops() == 1) {
    return vertices_.size();
  } else {
    return loop_starts_[i + 1] - loop_starts_[i];
  }
}

S2Point EncodedS2LaxPolygonShape::loop_vertex(int i, int j) const {
  ABSL_DCHECK_LT(i, num_loops());
  ABSL_DCHECK_LT(j, num_loop_vertices(i));
  if (num_loops() == 1) {
    return vertices_[j];
  } else {
    return vertices_[loop_starts_[i] + j];
  }
}

S2Shape::Edge EncodedS2LaxPolygonShape::edge(int e) const {
  ABSL_DCHECK_LT(e, num_edges());
  size_t e1 = e + 1;
  if (num_loops() == 1) {
    if (e1 == vertices_.size()) { e1 = 0; }
    return Edge(vertices_[e], vertices_[e1]);
  } else {
    // Method names are fully specified to enable inlining.
    ChainPosition pos = EncodedS2LaxPolygonShape::chain_position(e);
    return EncodedS2LaxPolygonShape::chain_edge(pos.chain_id, pos.offset);
  }
}

S2Shape::ReferencePoint EncodedS2LaxPolygonShape::GetReferencePoint() const {
  return s2shapeutil::GetReferencePoint(*this);
}

S2Shape::Chain EncodedS2LaxPolygonShape::chain(int i) const {
  ABSL_DCHECK_LT(i, num_loops());
  if (num_loops() == 1) {
    return Chain(0, vertices_.size());
  } else {
    int start = loop_starts_[i];
    return Chain(start, loop_starts_[i + 1] - start);
  }
}
