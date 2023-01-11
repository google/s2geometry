// Copyright 2018 Google Inc. All Rights Reserved.
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

#include "s2/s2polyline_measures.h"

#include "absl/types/span.h"
#include "s2/s1angle.h"
#include "s2/s2centroids.h"
#include "s2/s2point.h"
#include "s2/s2point_span.h"

namespace S2 {

S1Angle GetLength(S2PointSpan polyline) {
  S1Angle length;
  for (size_t i = 1; i < polyline.size(); ++i) {
    length += S1Angle(polyline[i - 1], polyline[i]);
  }
  return length;
}

S2Point GetCentroid(S2PointSpan polyline) {
  S2Point centroid;
  for (size_t i = 1; i < polyline.size(); ++i) {
    centroid += S2::TrueCentroid(polyline[i - 1], polyline[i]);
  }
  return centroid;
}

}  // namespace S2
