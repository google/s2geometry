// Copyright 2005 Google Inc. All Rights Reserved.
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

#ifndef S2_S2POINT_REGION_H_
#define S2_S2POINT_REGION_H_

#include <vector>

#include "absl/base/macros.h"
#include "absl/log/absl_check.h"

#include "s2/_fp_contract_off.h"  // IWYU pragma: keep
#include "s2/s1angle.h"
#include "s2/s2point.h"
#include "s2/s2pointutil.h"
#include "s2/s2region.h"
#include "s2/util/coding/coder.h"

class S2Cap;
class S2Cell;
class S2LatLngRect;

// An S2PointRegion is a region that contains a single point.  It is more
// expensive than the raw S2Point type and is useful mainly for completeness.
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.
class S2PointRegion final : public S2Region {
 public:
  // Create a region containing the given point, which must be unit length.
  explicit S2PointRegion(const S2Point& point);

  ~S2PointRegion() override;

  const S2Point& point() const { return point_; }

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  S2PointRegion* Clone() const override;
  S2Cap GetCapBound() const override;
  S2LatLngRect GetRectBound() const override;
  void GetCellUnionBound(std::vector<S2CellId>* cell_ids) const override;
  bool Contains(const S2Cell& cell) const override { return false; }
  bool MayIntersect(const S2Cell& cell) const override;
  bool Contains(const S2Point& p) const override { return (point_ == p); }

  // Appends a serialized representation of the S2Point to "encoder".
  //
  // REQUIRES: "encoder" uses the default constructor, so that its buffer
  //           can be enlarged as necessary by calling Ensure(int).
  void Encode(Encoder* encoder) const;

  // Decodes an S2Point encoded with Encode().  Returns true on success.
  // (Returns false if the encoded point is not unit length.)
  bool Decode(Decoder* decoder);

 private:
  S2Point point_;
};

inline S2PointRegion::S2PointRegion(const S2Point& point) : point_(point) {
  ABSL_DCHECK(S2::IsUnitLength(point));
}

#endif  // S2_S2POINT_REGION_H_
