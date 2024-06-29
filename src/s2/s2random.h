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

#ifndef S2_S2RANDOM_H_
#define S2_S2RANDOM_H_

#include "absl/random/bit_gen_ref.h"
#include "s2/s2cap.h"
#include "s2/s2cell_id.h"
#include "s2/s2latlng_rect.h"
#include "s2/s2point.h"
#include "s2/util/math/matrix3x3.h"

// This file defines functions for generating random points, cell ids, etc.
namespace s2random {

// Returns a random number in the closed interval `[lo, hi]`, whose log is
// uniformly distributed.  Note this is a continuous distribution, whereas
// `absl::LogUniform` is (as of 2024-05) a discrete distribution.
// REQUIRES: 0 < lo < hi < infinity
double LogUniform(absl::BitGenRef bitgen, double lo, double hi);

// Returns a random number skewed towards smaller values.  First a `base` is
// picked uniformly from range `[0,max_log]` and then `base` random bits are
// returned.  The effect is to pick a number in the range `[0,2^max_log-1]` with
// bias towards smaller numbers.
// REQUIRES: max_log > 0
// REQUIRES: max_log < 32
int SkewedInt(absl::BitGenRef bitgen, int max_log);

// Returns a random unit-length vector.
S2Point Point(absl::BitGenRef bitgen);

// Returns a right-handed coordinate frame (three orthonormal vectors).
void Frame(absl::BitGenRef bitgen, S2Point& x, S2Point& y, S2Point& z);
Matrix3x3_d Frame(absl::BitGenRef bitgen);

// Given a unit-length z-axis, computes x- and y-axes such that `(x,y,z)` is a
// right-handed coordinate frame (three orthonormal vectors).
void FrameAt(absl::BitGenRef bitgen, const S2Point& z, S2Point& x, S2Point& y);
Matrix3x3_d FrameAt(absl::BitGenRef bitgen, const S2Point& z);

// Returns a cap with a random axis such that the log of its area is
// uniformly distributed between the *logs* of the two given values.
// (The log of the cap angle is also approximately uniformly distributed.)
// REQUIRES: min_area > 0
// REQUIRES: max_area <= 4π
// REQUIRES: min_area <= max_area
S2Cap Cap(absl::BitGenRef bitgen, double min_area, double max_area);

// Returns a point chosen uniformly at random (with respect to area)
// from the given cap.
// REQUIRES: !cap.empty()
S2Point SamplePoint(absl::BitGenRef bitgen, const S2Cap& cap);

// Returns a point chosen uniformly at random (with respect to area on the
// sphere) from the given latitude-longitude rectangle.
// REQUIRES: !rect.empty()
S2Point SamplePoint(absl::BitGenRef bitgen, const S2LatLngRect& rect);

// Returns a random cell id at the given level or at a randomly chosen
// level.  The distribution is uniform over the space of cell ids,
// but only approximately uniform over the surface of the sphere.
S2CellId CellId(absl::BitGenRef bitgen, int level);
S2CellId CellId(absl::BitGenRef bitgen);

}  // namespace s2random

#endif  // S2_S2RANDOM_H_
