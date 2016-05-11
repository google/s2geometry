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

#include "s2/s2regionunion.h"

#include <memory>
#include <vector>

#include <gtest/gtest.h>
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cellid.h"
#include "s2/s2latlng.h"
#include "s2/s2latlngrect.h"
#include "s2/s2pointregion.h"
#include "s2/s2regioncoverer.h"
#include "s2/util/gtl/ptr_util.h"

using std::unique_ptr;
using std::vector;

namespace {

TEST(S2RegionUnionTest, Basic) {
  vector<S2Region*> regions;
  S2RegionUnion ru_empty(&regions);
  EXPECT_EQ(0, ru_empty.num_regions());
  EXPECT_EQ(S2Cap::Empty(), ru_empty.GetCapBound());
  EXPECT_EQ(S2LatLngRect::Empty(), ru_empty.GetRectBound());
  unique_ptr<S2Region> empty_clone(ru_empty.Clone());

  regions.push_back(new S2PointRegion(S2LatLng::FromDegrees(35, 40)
                                      .ToPoint()));
  regions.push_back(new S2PointRegion(S2LatLng::FromDegrees(-35, -40)
                                      .ToPoint()));

  // Check that Clone() returns a deep copy.
  auto two_points_orig = gtl::MakeUnique<S2RegionUnion>(&regions);
  EXPECT_TRUE(regions.empty());

  unique_ptr<S2RegionUnion> two_points(two_points_orig->Clone());
  two_points_orig.reset();
  EXPECT_EQ(S2LatLngRect(S2LatLng::FromDegrees(-35, -40),
                         S2LatLng::FromDegrees(35, 40)),
            two_points->GetRectBound());

  S2Cell face0 = S2Cell::FromFace(0);
  EXPECT_TRUE(two_points->MayIntersect(face0));
  EXPECT_FALSE(two_points->Contains(face0));

  EXPECT_TRUE(two_points->Contains(S2LatLng::FromDegrees(35, 40).ToPoint()));
  EXPECT_TRUE(two_points->Contains(S2LatLng::FromDegrees(-35, -40).ToPoint()));
  EXPECT_FALSE(two_points->Contains(S2LatLng::FromDegrees(0, 0).ToPoint()));

  // Check that we can Add() another region.
  unique_ptr<S2RegionUnion> three_points(two_points->Clone());
  EXPECT_FALSE(three_points->Contains(S2LatLng::FromDegrees(10, 10).ToPoint()));
  three_points->Add(new S2PointRegion(S2LatLng::FromDegrees(10, 10)
                                          .ToPoint()));
  EXPECT_TRUE(three_points->Contains(S2LatLng::FromDegrees(10, 10).ToPoint()));

  S2RegionCoverer coverer;
  coverer.set_max_cells(1);
  vector<S2CellId> covering;
  coverer.GetCovering(*two_points, &covering);
  EXPECT_EQ(1, covering.size());
  EXPECT_EQ(face0.id(), covering[0]);
}

}  // namespace
