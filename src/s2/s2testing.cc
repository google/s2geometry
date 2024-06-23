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
#include <cstdint>
#include <cstdlib>
#include <ios>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "absl/flags/flag.h"
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
                "Seed value that can be used in benchmarks.");

const double S2Testing::kEarthRadiusKm = 6371.01;

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

std::seed_seq S2Testing::MakeTaggedSeedSeq(string_view name,
                                           std::ostream& strm) {
  const int32_t seed = absl::GetFlag(FLAGS_s2_random_seed);
  const string seed_str = absl::StrCat(name, seed);
  strm << "Seeding " << name << " with " << seed;
  return std::seed_seq(seed_str.begin(), seed_str.end());
}



