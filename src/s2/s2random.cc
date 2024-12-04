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

#include "s2/s2random.h"

#include <cmath>
#include <cstdint>

#include "absl/log/absl_check.h"
#include "absl/random/bit_gen_ref.h"
#include "absl/random/distributions.h"
#include "s2/s2cap.h"
#include "s2/s2cell_id.h"
#include "s2/s2edge_crossings.h"
#include "s2/s2latlng_rect.h"
#include "s2/s2point.h"
#include "s2/util/math/matrix3x3.h"

namespace s2random {

double LogUniform(absl::BitGenRef bitgen, double lo, double hi) {
  ABSL_DCHECK_LT(0, lo);
  ABSL_DCHECK_LT(lo, hi);
  ABSL_DCHECK_LT(hi, std::numeric_limits<double>::infinity());
  return std::exp2(absl::Uniform(bitgen, std::log2(lo), std::log2(hi)));
}

int SkewedInt(absl::BitGenRef bitgen, int max_log) {
  ABSL_DCHECK_GE(max_log, 0);
  ABSL_DCHECK_LT(max_log, sizeof(int) * 8);
  int base = absl::Uniform(absl::IntervalClosedClosed, bitgen, 0, max_log);
  return absl::Uniform(bitgen, 0, 1 << base);
}

S2Point Point(absl::BitGenRef bitgen) {
  // The order of evaluation of function arguments is unspecified,
  // so we may not just call S2Point with three RandDouble-based args.
  // Use temporaries to induce sequence points between calls.
  double x = absl::Uniform(absl::IntervalClosedClosed, bitgen, -1.0, 1.0);
  double y = absl::Uniform(absl::IntervalClosedClosed, bitgen, -1.0, 1.0);
  double z = absl::Uniform(absl::IntervalClosedClosed, bitgen, -1.0, 1.0);
  return S2Point(x, y, z).Normalize();
}

void Frame(absl::BitGenRef bitgen, S2Point& x, S2Point& y, S2Point& z) {
  z = Point(bitgen);
  FrameAt(bitgen, z, x, y);
}

Matrix3x3_d Frame(absl::BitGenRef bitgen) {
  S2Point p = Point(bitgen);
  return FrameAt(bitgen, p);
}

void FrameAt(absl::BitGenRef bitgen, const S2Point& z,  //
             S2Point& x, S2Point& y) {
  x = S2::RobustCrossProd(z, Point(bitgen)).Normalize();
  y = S2::RobustCrossProd(z, x).Normalize();
}

Matrix3x3_d FrameAt(absl::BitGenRef bitgen, const S2Point& z) {
  S2Point x, y;
  FrameAt(bitgen, z, x, y);
  return Matrix3x3_d::FromCols(x, y, z);
}

S2CellId CellId(absl::BitGenRef bitgen, int level) {
  int face = absl::Uniform(bitgen, 0, S2CellId::kNumFaces);
  uint64_t pos =
      absl::Uniform<uint64_t>(bitgen, 0, uint64_t{1} << S2CellId::kPosBits);
  return S2CellId::FromFacePosLevel(face, pos, level);
}

S2CellId CellId(absl::BitGenRef bitgen) {
  int level =
      absl::Uniform(absl::IntervalClosedClosed, bitgen, 0, S2CellId::kMaxLevel);
  return CellId(bitgen, level);
}

S2Cap Cap(absl::BitGenRef bitgen, double min_area, double max_area) {
  // min_area cannot be 0 with the formula we're using.
  // XXX: Probably need a bigger value for numerical reasons?
  ABSL_DCHECK_GE(min_area, 0);
  ABSL_DCHECK_LE(min_area, max_area);
  ABSL_DCHECK_LE(max_area, 4 * M_PI);
  double exponent = absl::Uniform(absl::IntervalClosedClosed, bitgen, 0.0, 1.0);
  double cap_area = max_area * pow(min_area / max_area, exponent);
  ABSL_DCHECK_GE(cap_area, min_area);
  ABSL_DCHECK_LE(cap_area, max_area);

  // The surface area of a cap is 2*Pi times its height.
  return S2Cap::FromCenterArea(Point(bitgen), cap_area);
}

S2Point SamplePoint(absl::BitGenRef bitgen, const S2Cap& cap) {
  ABSL_DCHECK(!cap.is_empty());
  // We consider the cap axis to be the "z" axis.  We choose two other axes to
  // complete the coordinate frame.

  Matrix3x3_d m;
  S2::GetFrame(cap.center(), &m);

  // The surface area of a spherical cap is directly proportional to its
  // height.  First we choose a random height, and then we choose a random
  // point along the circle at that height.

  double h =
      absl::Uniform(absl::IntervalClosedClosed, bitgen, 0.0, cap.height());
  double theta = absl::Uniform(bitgen, 0.0, 2 * M_PI);
  double r = sqrt(h * (2 - h));  // Radius of circle.

  // The result should already be very close to unit-length, but we might as
  // well make it accurate as possible.
  return S2::FromFrame(m, S2Point(cos(theta) * r, sin(theta) * r, 1 - h))
      .Normalize();
}

S2Point SamplePoint(absl::BitGenRef bitgen, const S2LatLngRect& rect) {
  ABSL_DCHECK(!rect.is_empty());
  // First choose a latitude uniformly with respect to area on the sphere.
  double sin_lo = sin(rect.lat().lo());
  double sin_hi = sin(rect.lat().hi());
  double lat =
      asin(absl::Uniform(absl::IntervalClosedClosed, bitgen, sin_lo, sin_hi));

  // Now choose longitude uniformly within the given range.
  // XXX: Is closed-closed right for a full lng span?
  double lng = absl::Uniform(absl::IntervalClosedClosed, bitgen,
                             rect.lng().lo(), rect.lng().hi());
  return S2LatLng::FromRadians(lat, lng).Normalized().ToPoint();
}

}  // namespace s2random
