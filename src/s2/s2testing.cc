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

#include "s2/s2testing.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <ios>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "absl/log/absl_check.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"

#include "s2/base/commandlineflags.h"
#include "s2/base/types.h"
#include "s2/r1interval.h"
#include "s2/r2.h"
#include "s2/s1angle.h"
#include "s2/s1interval.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2cell_union.h"
#include "s2/s2edge_distances.h"
#include "s2/s2latlng.h"
#include "s2/s2latlng_rect.h"
#include "s2/s2lax_polygon_shape.h"
#include "s2/s2lax_polyline_shape.h"
#include "s2/s2loop.h"
#include "s2/s2point.h"
#include "s2/s2pointutil.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"
#include "s2/s2region.h"
#include "s2/s2shape_index.h"
#include "s2/s2text_format.h"
#include "s2/util/math/matrix3x3.h"

using absl::string_view;
using std::make_unique;
using std::string;
using std::unique_ptr;
using std::vector;

S2_DEFINE_int32(s2_random_seed, 1,
             "Seed value that can be passed to S2Testing::rnd.Reset()");

const double S2Testing::kEarthRadiusKm = 6371.01;

S2Testing::Random::Random() {
  // Unfortunately we can't use FLAGS_s2_random_seed here, because the default
  // S2Testing::Random instance is initialized before command-line flags have
  // been parsed.
  srandom(1);
}

void S2Testing::Random::Reset(int seed) {
  srandom(seed);
}

// Return a 64-bit unsigned integer whose lowest "num_bits" are random, and
// whose other bits are zero.
inline uint64 GetBits(int num_bits) {
  ABSL_DCHECK_GE(num_bits, 0);
  ABSL_DCHECK_LE(num_bits, 64);

  // This code uses random(), which returns an integer in the range
  // from 0 to (2^31)-1 inclusive (i.e. all of the lower 31 bits are
  // in play, and the 32nd bit and higher are 0) regardless of whether
  // its return type (long) is larger than 32 bits.  See
  //
  // www.gnu.org/software/libc/manual/html_node/BSD-Random.html#BSD-Random
  //
  // Note that at some point the manual page in linux claimed that the range
  // is 0 to RAND_MAX as defined in stdlib.h.  RAND_MAX however is part only
  // of the ISO rand() interface.  At least as of glibc-2.21, rand() is
  // simply an alias for random().  On other systems, rand() may differ,
  // but random() should always adhere to the behavior specified in BSD.
  static const int RAND_BITS = 31;

  uint64 result = 0;
  for (int bits = 0; bits < num_bits; bits += RAND_BITS) {
    result = (result << RAND_BITS) + random();
  }
  if (num_bits < 64) {  // Not legal to shift by full bitwidth of type
    result &= ((1ULL << num_bits) - 1);
  }
  return result;
}

uint64 S2Testing::Random::Rand64() { return GetBits(64); }

uint32 S2Testing::Random::Rand32() { return GetBits(32); }

double S2Testing::Random::RandDouble() {
  const int NUM_BITS = 53;
  return ldexp(GetBits(NUM_BITS), -NUM_BITS);
}

int32 S2Testing::Random::Uniform(int32 n) {
  ABSL_DCHECK_GT(n, 0);
  return static_cast<uint32>(RandDouble() * n);
}

double S2Testing::Random::UniformDouble(double min, double limit) {
  ABSL_DCHECK_LT(min, limit);
  return min + RandDouble() * (limit - min);
}

bool S2Testing::Random::OneIn(int32 n) { return Uniform(n) == 0; }

int32 S2Testing::Random::Skewed(int max_log) {
  ABSL_DCHECK_GE(max_log, 0);
  ABSL_DCHECK_LE(max_log, 31);
  int32 base = Uniform(max_log + 1);
  return GetBits(31) & ((1U << base) - 1);
}

S2Testing::Random S2Testing::rnd;

void S2Testing::AppendLoopVertices(const S2Loop& loop,
                                   vector<S2Point>* vertices) {
  int n = loop.num_vertices();
  const S2Point* base = &loop.vertex(0);
  ABSL_DCHECK_EQ(&loop.vertex(n - 1), base + n - 1);
  vertices->insert(vertices->end(), base, base + n);
}

vector<S2Point> S2Testing::MakeRegularPoints(const S2Point& center,
                                             S1Angle radius,
                                             int num_vertices) {
  unique_ptr<S2Loop> loop(
      S2Loop::MakeRegularLoop(center, radius, num_vertices));
  vector<S2Point> points;
  points.reserve(loop->num_vertices());
  for (int i = 0; i < loop->num_vertices(); i++) {
    points.push_back(loop->vertex(i));
  }
  return points;
}

S1Angle S2Testing::MetersToAngle(double meters) {
  return KmToAngle(0.001 * meters);
}

S1Angle S2Testing::KmToAngle(double km) {
  return S1Angle::Radians(km / kEarthRadiusKm);
}

double S2Testing::AreaToMeters2(double steradians) {
  return 1e6 * AreaToKm2(steradians);
}

double S2Testing::AreaToKm2(double steradians) {
  return steradians * kEarthRadiusKm * kEarthRadiusKm;
}

S2Point S2Testing::RandomPoint() {
  // The order of evaluation of function arguments is unspecified,
  // so we may not just call S2Point with three RandDouble-based args.
  // Use temporaries to induce sequence points between calls.
  double x = rnd.UniformDouble(-1, 1);
  double y = rnd.UniformDouble(-1, 1);
  double z = rnd.UniformDouble(-1, 1);
  return S2Point(x, y, z).Normalize();
}

void S2Testing::GetRandomFrame(S2Point* x, S2Point* y, S2Point* z) {
  *z = RandomPoint();
  GetRandomFrameAt(*z, x, y);
}

Matrix3x3_d S2Testing::GetRandomFrame() {
  return GetRandomFrameAt(RandomPoint());
}

void S2Testing::GetRandomFrameAt(const S2Point& z, S2Point* x, S2Point *y) {
  *x = z.CrossProd(RandomPoint()).Normalize();
  *y = z.CrossProd(*x).Normalize();
}

Matrix3x3_d S2Testing::GetRandomFrameAt(const S2Point& z) {
  S2Point x, y;
  GetRandomFrameAt(z, &x, &y);
  return Matrix3x3_d::FromCols(x, y, z);
}

S2CellId S2Testing::GetRandomCellId(int level) {
  int face = rnd.Uniform(S2CellId::kNumFaces);
  uint64 pos = rnd.Rand64() & ((1ULL << S2CellId::kPosBits) - 1);
  return S2CellId::FromFacePosLevel(face, pos, level);
}

S2CellId S2Testing::GetRandomCellId() {
  return GetRandomCellId(rnd.Uniform(S2CellId::kMaxLevel + 1));
}

S2Cap S2Testing::GetRandomCap(double min_area, double max_area) {
  double cap_area = max_area * pow(min_area / max_area, rnd.RandDouble());
  ABSL_DCHECK_GE(cap_area, min_area);
  ABSL_DCHECK_LE(cap_area, max_area);

  // The surface area of a cap is 2*Pi times its height.
  return S2Cap::FromCenterArea(RandomPoint(), cap_area);
}

void S2Testing::ConcentricLoopsPolygon(const S2Point& center,
                                       int num_loops,
                                       int num_vertices_per_loop,
                                       S2Polygon* polygon) {
  Matrix3x3_d m;
  S2::GetFrame(center, &m);
  vector<unique_ptr<S2Loop>> loops;
  for (int li = 0; li < num_loops; ++li) {
    vector<S2Point> vertices;
    double radius = 0.005 * (li + 1) / num_loops;
    double radian_step = 2 * M_PI / num_vertices_per_loop;
    for (int vi = 0; vi < num_vertices_per_loop; ++vi) {
      double angle = vi * radian_step;
      S2Point p(radius * cos(angle), radius * sin(angle), 1);
      vertices.push_back(S2::FromFrame(m, p.Normalize()));
    }
    loops.push_back(make_unique<S2Loop>(vertices));
  }
  polygon->InitNested(std::move(loops));
}

S2Point S2Testing::SamplePoint(const S2Cap& cap) {
  // We consider the cap axis to be the "z" axis.  We choose two other axes to
  // complete the coordinate frame.

  Matrix3x3_d m;
  S2::GetFrame(cap.center(), &m);

  // The surface area of a spherical cap is directly proportional to its
  // height.  First we choose a random height, and then we choose a random
  // point along the circle at that height.

  double h = rnd.RandDouble() * cap.height();
  double theta = 2 * M_PI * rnd.RandDouble();
  double r = sqrt(h * (2 - h));  // Radius of circle.

  // The result should already be very close to unit-length, but we might as
  // well make it accurate as possible.
  return S2::FromFrame(m, S2Point(cos(theta) * r, sin(theta) * r, 1 - h))
         .Normalize();
}

S2Point S2Testing::SamplePoint(const S2LatLngRect& rect) {
  // First choose a latitude uniformly with respect to area on the sphere.
  double sin_lo = sin(rect.lat().lo());
  double sin_hi = sin(rect.lat().hi());
  double lat = asin(rnd.UniformDouble(sin_lo, sin_hi));

  // Now choose longitude uniformly within the given range.
  double lng = rect.lng().lo() + rnd.RandDouble() * rect.lng().GetLength();
  return S2LatLng::FromRadians(lat, lng).Normalized().ToPoint();
}

void S2Testing::CheckCovering(const S2Region& region,
                              const S2CellUnion& covering,
                              bool check_tight, S2CellId id) {
  if (!id.is_valid()) {
    for (int face = 0; face < 6; ++face) {
      CheckCovering(region, covering, check_tight, S2CellId::FromFace(face));
    }
    return;
  }

  if (!region.MayIntersect(S2Cell(id))) {
    // If region does not intersect id, then neither should the covering.
    if (check_tight) ABSL_CHECK(!covering.Intersects(id));

  } else if (!covering.Contains(id)) {
    // The region may intersect id, but we can't assert that the covering
    // intersects id because we may discover that the region does not actually
    // intersect upon further subdivision.  (MayIntersect is not exact.)
    ABSL_CHECK(!region.Contains(S2Cell(id)));
    ABSL_CHECK(!id.is_leaf());
    S2CellId end = id.child_end();
    S2CellId child;
    for (child = id.child_begin(); child != end; child = child.next()) {
      CheckCovering(region, covering, check_tight, child);
    }
  }
}

