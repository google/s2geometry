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

#ifndef S2_S2POINTREGION_H_
#define S2_S2POINTREGION_H_

#include <glog/logging.h>
#include "s2/base/macros.h"
#include "s2/fpcontractoff.h"
#include "s2/s1angle.h"
#include "s2/s2pointutil.h"
#include "s2/s2region.h"

class Decoder;
class Encoder;
class S2Cap;
class S2Cell;
class S2LatLngRect;

// An S2PointRegion is a region that contains a single point.  It is more
// expensive than the raw S2Point type and is useful mainly for completeness.
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator.
class S2PointRegion : public S2Region {
 public:
  // Create a region containing the given point, which must be unit length.
  explicit S2PointRegion(S2Point const& point);

  ~S2PointRegion() override;

  S2Point const& point() const { return point_; }

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  S2PointRegion* Clone() const override;
  S2Cap GetCapBound() const override;
  S2LatLngRect GetRectBound() const override;
  bool Contains(S2Cell const& cell) const override { return false; }
  bool MayIntersect(S2Cell const& cell) const override;
  bool VirtualContainsPoint(S2Point const& p) const override {
    return Contains(p);
  }
  bool Contains(S2Point const& p) const { return (point_ == p); }
  void Encode(Encoder* const encoder) const override;
  // Ensures the decoded point has unit length.
  bool Decode(Decoder* const decoder) override;

 private:
  S2Point point_;
};

inline S2PointRegion::S2PointRegion(S2Point const& point) : point_(point) {
  DCHECK(S2::IsUnitLength(point));
}

#endif  // S2_S2POINTREGION_H_
