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

#include "s2testing.h"

#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <sys/resource.h>   // for rusage, RUSAGE_SELF
#include <sys/time.h>
#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include <gflags/gflags.h>
#include "base/integral_types.h"
#include <glog/logging.h>
#include "base/stringprintf.h"
#include "strings/serialize.h"
#include "strings/split.h"
#include "r1interval.h"
#include "s1angle.h"
#include "s1interval.h"
#include "s2cap.h"
#include "s2cell.h"
#include "s2cellunion.h"
#include "s2latlng.h"
#include "s2latlngrect.h"
#include "s2loop.h"
#include "s2polygon.h"
#include "s2polyline.h"
#include "s2region.h"
#include "s2textformat.h"
#include "util/gtl/stl_util.h"
#include "util/math/matrix3x3.h"

using std::max;
using std::pair;
using std::unique_ptr;
using std::vector;

DEFINE_int32(s2_random_seed, 1, "Initial random seed for S2Testing::rnd");

double const S2Testing::kEarthRadiusKm = 6371.01;

S2Testing::Random::Random() {
  Reset(FLAGS_s2_random_seed);
}

void S2Testing::Random::Reset(int seed) {
  srandom(seed);
}

// Return a 64-bit unsigned integer whose lowest "num_bits" are random, and
// whose other bits are zero.
inline uint64 GetBits(int num_bits) {
  DCHECK_GE(num_bits, 0);
  DCHECK_LE(num_bits, 64);

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
  static int const RAND_BITS = 31;

  uint64 result = 0;
  for (int bits = 0; bits < num_bits; bits += RAND_BITS) {
    result = (result << RAND_BITS) + random();
  }
  if (num_bits < 64) {  // Not legal to shift by full bitwidth of type
    result &= ((1ULL << num_bits) - 1);
  }
  return result;
}

uint64 S2Testing::Random::Rand64() {
  return GetBits(64);
}

uint32 S2Testing::Random::Rand32() {
  return GetBits(32);
}

double S2Testing::Random::RandDouble() {
  int const NUM_BITS = 53;
  return ldexp(GetBits(NUM_BITS), -NUM_BITS);
}

int32 S2Testing::Random::Uniform(int32 n) {
  return static_cast<uint32>(RandDouble() * n);
}

double S2Testing::Random::UniformDouble(double min, double limit) {
  return min + RandDouble() * (limit - min);
}

bool S2Testing::Random::OneIn(int32 n) {
  return Uniform(n) == 0;
}

int32 S2Testing::Random::Skewed(int max_log) {
  DCHECK_GE(max_log, 0);
  DCHECK_LE(max_log, 31);
  int32 base = Uniform(max_log + 1);
  return GetBits(31) & ((1U << base) - 1);
}

S2Testing::Random S2Testing::rnd;

void S2Testing::AppendLoopVertices(S2Loop const& loop,
                                   vector<S2Point>* vertices) {
  int n = loop.num_vertices();
  S2Point const* base = &loop.vertex(0);
  DCHECK_EQ(&loop.vertex(n - 1), base + n - 1);
  vertices->insert(vertices->end(), base, base + n);
}

vector<S2Point> S2Testing::MakeRegularPoints(S2Point const& center,
                                             S1Angle radius,
                                             int num_vertices) {
  unique_ptr<S2Loop> loop(MakeRegularLoop(center, radius, num_vertices));
  vector<S2Point> points;
  for (int i = 0; i < loop.get()->num_vertices(); i++) {
    points.push_back(loop.get()->vertex(i));
  }
  return points;
}

S2Loop* S2Testing::MakeRegularLoop(S2Point const& center,
                                   S1Angle radius,
                                   int num_vertices) {
  return S2Loop::MakeRegularLoop(center, radius, num_vertices);
}

S2Loop* S2Testing::MakeRegularLoop(Matrix3x3_d const& frame,
                                   S1Angle radius, int num_vertices) {
  // We construct the loop in the given frame coordinates, with the center at
  // (0, 0, 1).  For a loop of radius "r", the loop vertices have the form
  // (x, y, z) where x^2 + y^2 = sin(r) and z = cos(r).  The distance on the
  // sphere (arc length) from each vertex to the center is acos(cos(r)) = r.
  double z = cos(radius.radians());
  double r = sin(radius.radians());
  double radian_step = 2 * M_PI / num_vertices;
  vector<S2Point> vertices;
  for (int i = 0; i < num_vertices; ++i) {
    double angle = i * radian_step;
    S2Point p(r * cos(angle), r * sin(angle), z);
    vertices.push_back(S2::FromFrame(frame, p).Normalize());
  }
  return new S2Loop(vertices);
}

S1Angle S2Testing::KmToAngle(double km) {
  return S1Angle::Radians(km / kEarthRadiusKm);
}

S1Angle S2Testing::MetersToAngle(double meters) {
  return KmToAngle(0.001 * meters);
}

void DumpLoop(S2Loop const* loop) {
  // Only for calling from a debugger.
  std::cout << s2textformat::ToString(loop) << "\n";
}

void DumpPolyline(S2Polyline const* polyline) {
  // Only for calling from a debugger.
  std::cout << s2textformat::ToString(polyline) << "\n";
}

void DumpPolygon(S2Polygon const* polygon) {
  // Only for calling from a debugger.
  std::cout << s2textformat::ToString(polygon) << "\n";
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

void S2Testing::GetRandomFrame(Vector3_d* x, Vector3_d* y, Vector3_d* z) {
  *z = RandomPoint();
  GetRandomFrameAt(*z, x, y);
}

Matrix3x3_d S2Testing::GetRandomFrame() {
  return GetRandomFrameAt(RandomPoint());
}

void S2Testing::GetRandomFrameAt(S2Point const& z, S2Point* x, S2Point *y) {
  *x = z.CrossProd(RandomPoint()).Normalize();
  *y = z.CrossProd(*x).Normalize();
}

Matrix3x3_d S2Testing::GetRandomFrameAt(S2Point const& z) {
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
  DCHECK_GE(cap_area, min_area);
  DCHECK_LE(cap_area, max_area);

  // The surface area of a cap is 2*Pi times its height.
  return S2Cap::FromCenterArea(RandomPoint(), cap_area);
}

void S2Testing::ConcentricLoopsPolygon(S2Point const& center,
                                       int num_loops,
                                       int num_vertices_per_loop,
                                       S2Polygon* polygon) {
  Matrix3x3_d m;
  S2::GetFrame(center, &m);
  vector<S2Loop*> loops;
  for (int li = 0; li < num_loops; ++li) {
    vector<S2Point> vertices;
    double radius = 0.005 * (li + 1) / num_loops;
    double radian_step = 2 * M_PI / num_vertices_per_loop;
    for (int vi = 0; vi < num_vertices_per_loop; ++vi) {
      double angle = vi * radian_step;
      S2Point p(radius * cos(angle), radius * sin(angle), 1);
      vertices.push_back(S2::FromFrame(m, p.Normalize()));
    }
    loops.push_back(new S2Loop(vertices));
  }
  polygon->Init(&loops);
}

S2Point S2Testing::SamplePoint(S2Cap const& cap) {
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

S2Point S2Testing::SamplePoint(S2LatLngRect const& rect) {
  // First choose a latitude uniformly with respect to area on the sphere.
  double sin_lo = sin(rect.lat().lo());
  double sin_hi = sin(rect.lat().hi());
  double lat = asin(rnd.UniformDouble(sin_lo, sin_hi));

  // Now choose longitude uniformly within the given range.
  double lng = rect.lng().lo() + rnd.RandDouble() * rect.lng().GetLength();
  return S2LatLng::FromRadians(lat, lng).Normalized().ToPoint();
}

void S2Testing::CheckCovering(S2Region const& region,
                              S2CellUnion const& covering,
                              bool check_tight, S2CellId id) {
  if (!id.is_valid()) {
    for (int face = 0; face < 6; ++face) {
      CheckCovering(region, covering, check_tight, S2CellId::FromFace(face));
    }
    return;
  }

  if (!region.MayIntersect(S2Cell(id))) {
    // If region does not intersect id, then neither should the covering.
    if (check_tight) CHECK(!covering.Intersects(id));

  } else if (!covering.Contains(id)) {
    // The region may intersect id, but we can't assert that the covering
    // intersects id because we may discover that the region does not actually
    // intersect upon further subdivision.  (MayIntersect is not exact.)
    CHECK(!region.Contains(S2Cell(id)));
    CHECK(!id.is_leaf());
    S2CellId end = id.child_end();
    S2CellId child;
    for (child = id.child_begin(); child != end; child = child.next()) {
      CheckCovering(region, covering, check_tight, child);
    }
  }
}

double S2Testing::GetCpuTime() {
  struct rusage ru;
  CHECK_EQ(getrusage(RUSAGE_SELF, &ru), 0);
  return ru.ru_utime.tv_sec + ru.ru_utime.tv_usec / 1e6;
}

void S2Testing::DeleteLoops(vector<S2Loop*>* loops) {
  // TODO(user): Update all callers and delete this function.
  STLDeleteElements(loops);
}


S2Testing::Fractal::Fractal()
    : max_level_(-1), min_level_arg_(-1), min_level_(-1),
      dimension_(log(4)/log(3)), /* standard Koch curve */
      edge_fraction_(0), offset_fraction_(0) {
  ComputeOffsets();
}

void S2Testing::Fractal::set_max_level(int max_level) {
  DCHECK_GE(max_level, 0);
  max_level_ = max_level;
  ComputeMinLevel();
}

void S2Testing::Fractal::set_min_level(int min_level_arg) {
  DCHECK_GE(min_level_arg, -1);
  min_level_arg_ = min_level_arg;
  ComputeMinLevel();
}

void S2Testing::Fractal::ComputeMinLevel() {
  if (min_level_arg_ >= 0 && min_level_arg_ <= max_level_) {
    min_level_ = min_level_arg_;
  } else {
    min_level_ = max_level_;
  }
}

void S2Testing::Fractal::set_fractal_dimension(double dimension) {
  DCHECK_GE(dimension, 1.0);
  DCHECK_LT(dimension, 2.0);
  dimension_ = dimension;
  ComputeOffsets();
}

void S2Testing::Fractal::ComputeOffsets() {
  edge_fraction_ = pow(4.0, -1.0 / dimension_);
  offset_fraction_ = sqrt(edge_fraction_ - 0.25);
}

void S2Testing::Fractal::SetLevelForApproxMinEdges(int min_edges) {
  // Map values in the range [3*(4**n)/2, 3*(4**n)*2) to level n.
  set_min_level(round(0.5 * log2(min_edges / 3)));
}

void S2Testing::Fractal::SetLevelForApproxMaxEdges(int max_edges) {
  // Map values in the range [3*(4**n)/2, 3*(4**n)*2) to level n.
  set_max_level(round(0.5 * log2(max_edges / 3)));
}

double S2Testing::Fractal::min_radius_factor() const {
  // The minimum radius is attained at one of the vertices created by the
  // first subdivision step as long as the dimension is not too small (at
  // least kMinDimensionForMinRadiusAtLevel1, see below).  Otherwise we fall
  // back on the incircle radius of the original triangle, which is always a
  // lower bound (and is attained when dimension = 1).
  //
  // The value below was obtained by letting AE be an original triangle edge,
  // letting ABCDE be the corresponding polyline after one subdivision step,
  // and then letting BC be tangent to the inscribed circle at the center of
  // the fractal O.  This gives rise to a pair of similar triangles whose edge
  // length ratios can be used to solve for the corresponding "edge fraction".
  // This method is slightly conservative because it is computed using planar
  // rather than spherical geometry.  The value below is equal to
  // -log(4)/log((2 + cbrt(2) - cbrt(4))/6).
  double const kMinDimensionForMinRadiusAtLevel1 = 1.0852230903040407;
  if (dimension_ >= kMinDimensionForMinRadiusAtLevel1) {
    return sqrt(1 + 3 * edge_fraction_ * (edge_fraction_ - 1));
  }
  return 0.5;
}

double S2Testing::Fractal::max_radius_factor() const {
  // The maximum radius is always attained at either an original triangle
  // vertex or at a middle vertex from the first subdivision step.
  return max(1.0, offset_fraction_ * sqrt(3) + 0.5);
}

void S2Testing::Fractal::GetR2Vertices(vector<R2Point>* vertices) const {
  // The Koch "snowflake" consists of three Koch curves whose initial edges
  // form an equilateral triangle.
  R2Point v0(1.0, 0.0);
  R2Point v1(-0.5, sqrt(3)/2);
  R2Point v2(-0.5, -sqrt(3)/2);
  GetR2VerticesHelper(v0, v1, 0, vertices);
  GetR2VerticesHelper(v1, v2, 0, vertices);
  GetR2VerticesHelper(v2, v0, 0, vertices);
}

// Given the two endpoints (v0,v4) of an edge, recursively subdivide the edge
// to the desired level, and insert all vertices of the resulting curve up to
// but not including the endpoint "v4".
void S2Testing::Fractal::GetR2VerticesHelper(R2Point const& v0,
                                             R2Point const& v4, int level,
                                             vector<R2Point>* vertices) const {
  if (level >= min_level_ && S2Testing::rnd.OneIn(max_level_ - level + 1)) {
    // Stop subdivision at this level.
    vertices->push_back(v0);
    return;
  }
  // Otherwise compute the intermediate vertices v1, v2, and v3.
  Vector2_d dir = v4 - v0;
  R2Point v1 = v0 + edge_fraction_ * dir;
  R2Point v2 = 0.5 * (v0 + v4) - offset_fraction_ * dir.Ortho();
  R2Point v3 = v4 - edge_fraction_ * dir;

  // And recurse on the four sub-edges.
  GetR2VerticesHelper(v0, v1, level+1, vertices);
  GetR2VerticesHelper(v1, v2, level+1, vertices);
  GetR2VerticesHelper(v2, v3, level+1, vertices);
  GetR2VerticesHelper(v3, v4, level+1, vertices);
}

S2Loop* S2Testing::Fractal::MakeLoop(Matrix3x3_d const& frame,
                                     S1Angle nominal_radius) const {
  vector<R2Point> r2vertices;
  GetR2Vertices(&r2vertices);
  vector<S2Point> vertices;
  for (int i = 0; i < r2vertices.size(); ++i) {
    // Convert each vertex to polar coordinates.
    R2Point const& v = r2vertices[i];
    double theta = atan2(v[1], v[0]);
    double radius = nominal_radius.radians() * v.Norm();

    // See the comments in MakeRegularLoop.
    double z = cos(radius);
    double r = sin(radius);
    S2Point p(r * cos(theta), r * sin(theta), z);
    vertices.push_back(S2::FromFrame(frame, p).Normalize());
  }
  return new S2Loop(vertices);
}
