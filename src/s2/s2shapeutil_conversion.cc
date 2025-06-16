// Copyright 2021 Google Inc. All Rights Reserved.
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

#include "s2/s2shapeutil_conversion.h"

#include <memory>
#include <utility>
#include <vector>

#include "absl/log/absl_check.h"
#include "s2/s2loop.h"
#include "s2/s2point.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"
#include "s2/s2shape.h"
#include "s2/s2shape_measures.h"

namespace s2shapeutil {

using std::make_unique;
using std::unique_ptr;
using std::vector;

vector<S2Point> ShapeToS2Points(const S2Shape& multipoint) {
  ABSL_DCHECK_EQ(multipoint.dimension(), 0);
  vector<S2Point> points;
  points.reserve(multipoint.num_edges());
  for (int i = 0; i < multipoint.num_edges(); ++i) {
    points.push_back(multipoint.edge(i).v0);
  }
  return points;
}

unique_ptr<S2Polyline> ShapeToS2Polyline(const S2Shape& line) {
  ABSL_DCHECK_EQ(line.dimension(), 1);
  ABSL_DCHECK_EQ(line.num_chains(), 1);
  vector<S2Point> vertices;
  S2::GetChainVertices(line, 0, &vertices);
  return make_unique<S2Polyline>(std::move(vertices));
}

unique_ptr<S2Polygon> ShapeToS2Polygon(const S2Shape& poly) {
  if (poly.is_full()) {
    return make_unique<S2Polygon>(make_unique<S2Loop>(S2Loop::kFull()));
  }
  ABSL_DCHECK_EQ(poly.dimension(), 2);
  vector<unique_ptr<S2Loop>> loops;
  vector<S2Point> vertices;
  for (int i = 0; i < poly.num_chains(); ++i) {
    S2::GetChainVertices(poly, i, &vertices);
    loops.push_back(make_unique<S2Loop>(vertices));
  }
  auto output_poly = make_unique<S2Polygon>();
  if (loops.size() == 1) {
    output_poly->Init(std::move(loops[0]));
  } else {
    output_poly->InitOriented(std::move(loops));
  }
  return output_poly;
}

}  // namespace s2shapeutil
