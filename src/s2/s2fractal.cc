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

#include "s2/s2fractal.h"

#include <memory>
#include <vector>

#include "absl/log/absl_check.h"
#include "absl/random/distributions.h"
#include "s2/r2.h"
#include "s2/s1angle.h"

using std::make_unique;
using std::max;
using std::unique_ptr;
using std::vector;

S2Fractal::S2Fractal(absl::BitGenRef bitgen) : bitgen_(bitgen) {
  ComputeOffsets();
}

void S2Fractal::set_max_level(int max_level) {
  ABSL_DCHECK_GE(max_level, 0);
  max_level_ = max_level;
  ComputeMinLevel();
}

void S2Fractal::set_min_level(int min_level_arg) {
  ABSL_DCHECK_GE(min_level_arg, -1);
  min_level_arg_ = min_level_arg;
  ComputeMinLevel();
}

void S2Fractal::ComputeMinLevel() {
  if (min_level_arg_ >= 0 && min_level_arg_ <= max_level_) {
    min_level_ = min_level_arg_;
  } else {
    min_level_ = max_level_;
  }
}

void S2Fractal::set_fractal_dimension(double dimension) {
  ABSL_DCHECK_GE(dimension, 1.0);
  ABSL_DCHECK_LT(dimension, 2.0);
  dimension_ = dimension;
  ComputeOffsets();
}

void S2Fractal::ComputeOffsets() {
  edge_fraction_ = pow(4.0, -1.0 / dimension_);
  offset_fraction_ = sqrt(edge_fraction_ - 0.25);
}

void S2Fractal::SetLevelForApproxMinEdges(int min_edges) {
  // Map values in the range [3*(4**n)/2, 3*(4**n)*2) to level n.
  set_min_level(round(0.5 * log2(min_edges / 3)));
}

void S2Fractal::SetLevelForApproxMaxEdges(int max_edges) {
  // Map values in the range [3*(4**n)/2, 3*(4**n)*2) to level n.
  set_max_level(round(0.5 * log2(max_edges / 3)));
}

double S2Fractal::min_radius_factor() const {
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
  const double kMinDimensionForMinRadiusAtLevel1 = 1.0852230903040407;
  if (dimension_ >= kMinDimensionForMinRadiusAtLevel1) {
    return sqrt(1 + 3 * edge_fraction_ * (edge_fraction_ - 1));
  }
  return 0.5;
}

double S2Fractal::max_radius_factor() const {
  // The maximum radius is always attained at either an original triangle
  // vertex or at a middle vertex from the first subdivision step.
  return max(1.0, offset_fraction_ * sqrt(3) + 0.5);
}

vector<R2Point> S2Fractal::GetR2Vertices() {
  // The Koch "snowflake" consists of three Koch curves whose initial edges
  // form an equilateral triangle.
  vector<R2Point> vertices;
  R2Point v0(1.0, 0.0);
  R2Point v1(-0.5, sqrt(3) / 2);
  R2Point v2(-0.5, -sqrt(3) / 2);
  GetR2VerticesHelper(v0, v1, 0, vertices);
  GetR2VerticesHelper(v1, v2, 0, vertices);
  GetR2VerticesHelper(v2, v0, 0, vertices);
  return vertices;
}

// Given the two endpoints (v0,v4) of an edge, recursively subdivide the edge
// to the desired level, and insert all vertices of the resulting curve up to
// but not including the endpoint "v4".
void S2Fractal::GetR2VerticesHelper(const R2Point& v0, const R2Point& v4,
                                    int level, vector<R2Point>& vertices) {
  const int levels_remaining = max_level_ - level + 1;
  if (level >= min_level_ &&
      absl::Bernoulli(bitgen_, 1.0 / levels_remaining)) {
    // Stop subdivision at this level.
    vertices.push_back(v0);
    return;
  }
  // Otherwise compute the intermediate vertices v1, v2, and v3.
  Vector2_d dir = v4 - v0;
  R2Point v1 = v0 + edge_fraction_ * dir;
  R2Point v2 = 0.5 * (v0 + v4) - offset_fraction_ * dir.Ortho();
  R2Point v3 = v4 - edge_fraction_ * dir;

  // And recurse on the four sub-edges.
  GetR2VerticesHelper(v0, v1, level + 1, vertices);
  GetR2VerticesHelper(v1, v2, level + 1, vertices);
  GetR2VerticesHelper(v2, v3, level + 1, vertices);
  GetR2VerticesHelper(v3, v4, level + 1, vertices);
}

unique_ptr<S2Loop> S2Fractal::MakeLoop(const Matrix3x3_d& frame,
                                       S1Angle nominal_radius) {
  vector<R2Point> r2vertices = GetR2Vertices();
  vector<S2Point> vertices;
  vertices.reserve(r2vertices.size());
  double r = nominal_radius.radians();
  for (const R2Point& v : r2vertices) {
    S2Point p(v[0] * r, v[1] * r, 1);
    vertices.push_back(S2::FromFrame(frame, p).Normalize());
  }
  return make_unique<S2Loop>(vertices);
}
