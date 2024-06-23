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

#ifndef S2_S2TESTING_H_
#define S2_S2TESTING_H_

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "absl/strings/string_view.h"

#include "s2/_fp_contract_off.h"  // IWYU pragma: keep
#include "s2/base/commandlineflags.h"
#include "s2/base/commandlineflags_declare.h"
#include "s2/r2.h"
#include "s2/s1angle.h"
#include "s2/s1chord_angle.h"
#include "s2/s2cell_id.h"
#include "s2/s2point.h"
#include "s2/s2region.h"
#include "s2/util/math/matrix3x3.h"

class S1Angle;
class S2Cap;
class S2CellUnion;
class S2LatLng;
class S2LatLngRect;
class S2Loop;
class S2Polygon;
class S2Polyline;
class S2Region;

// For benchmarks, this seed can be set on the command line, then combined
// with other data to vary the seed.  Typical usage is:
// ```
// const std::string seed_str =
//     absl::StrCat(__func__, absl::GetFlag(FLAGS_s2_random_seed);
// const std::seed_seq seed(seed_str.begin(), seed_str.end());
// std::mt19937_64 bitgen(seed);
// // Use `bitgen`.
// ```
S2_DECLARE_int32(s2_random_seed);

// This class defines various static functions that are useful for writing
// unit tests.
class S2Testing {
 public:
  // Returns a vector of points shaped as a regular polygon with
  // num_vertices vertices, all on a circle of the specified angular
  // radius around the center.  The radius is the actual distance from
  // the center to the circle along the sphere.
  //
  // If you want to construct a regular polygon, try this:
  //   S2Polygon polygon(S2Loop::MakeRegularLoop(center, radius, num_vertices));
  static std::vector<S2Point> MakeRegularPoints(const S2Point& center,
                                           S1Angle radius,
                                           int num_vertices);

  // Append the vertices of "loop" to "vertices".
  static void AppendLoopVertices(const S2Loop& loop,
                                 std::vector<S2Point>* vertices);

  // Convert a distance on the Earth's surface to an angle.
  // Do not use these methods in non-testing code; use s2earth.h instead.
  static S1Angle MetersToAngle(double meters);
  static S1Angle KmToAngle(double km);

  // Convert an area in steradians (as returned by the S2 area methods) to
  // square meters or square kilometers.
  static double AreaToMeters2(double steradians);
  static double AreaToKm2(double steradians);

  // The Earth's mean radius in kilometers (according to NASA).
  static const double kEarthRadiusKm;

  // Return a polygon with the specified center, number of concentric loops
  // and vertices per loop.
  static void ConcentricLoopsPolygon(const S2Point& center,
                                     int num_loops,
                                     int num_vertices_per_loop,
                                     S2Polygon* polygon);

  // Checks that "covering" completely covers the given region.  If
  // "check_tight" is true, also checks that it does not contain any cells
  // that do not intersect the given region.  ("id" is only used internally.)
  static void CheckCovering(const S2Region& region,
                            const S2CellUnion& covering,
                            bool check_tight,
                            S2CellId id = S2CellId());

  // Combines `name` with `FLAGS_s2_random_seed` to make a seed sequence.
  // TODO(user): Remove this when `S2Testing::MakeTaggedSeedSeq` is released.
  static std::seed_seq MakeTaggedSeedSeq(absl::string_view name,
                                         std::ostream& strm);


 private:
  // Contains static methods
  S2Testing() = delete;
  S2Testing(const S2Testing&) = delete;
  void operator=(const S2Testing&) = delete;
};

// Compare two sets of "closest" items, where "expected" is computed via brute
// force (i.e., considering every possible candidate) and "actual" is computed
// using a spatial data structure.  Here "max_size" is a bound on the maximum
// number of items, "max_distance" is a limit on the distance to any item, and
// "max_error" is the maximum error allowed when selecting which items are
// closest (see S2ClosestEdgeQuery::Options::max_error).
template <typename Id, typename Distance>
bool CheckDistanceResults(
    const std::vector<std::pair<Distance, Id>>& expected,
    const std::vector<std::pair<Distance, Id>>& actual,
    int max_size, Distance max_distance, typename Distance::Delta max_error);


//////////////////// Implementation Details Follow ////////////////////////


namespace S2 {
namespace internal {

// Check that result set "x" contains all the expected results from "y", and
// does not include any duplicate results.
template <typename Id, typename Distance>
bool CheckResultSet(const std::vector<std::pair<Distance, Id>>& x,
                    const std::vector<std::pair<Distance, Id>>& y,
                    int max_size, Distance max_distance,
                    typename Distance::Delta max_error,
                    typename Distance::Delta max_pruning_error,
                    absl::string_view label) {
  using Result = std::pair<Distance, Id>;
  // Results should be sorted by distance, but not necessarily then by Id.
  EXPECT_TRUE(std::is_sorted(x.begin(), x.end(),
                             [](const Result& x, const Result& y) {
                               return x.first < y.first;
                             }));

  // Result set X should contain all the items from Y whose distance is less
  // than "limit" computed below.
  Distance limit = Distance::Zero();
  if (x.size() < max_size) {
    // Result set X was not limited by "max_size", so it should contain all
    // the items up to "max_distance", except that a few items right near the
    // distance limit may be missed because the distance measurements used for
    // pruning S2Cells are not conservative.
    if (max_distance == Distance::Infinity()) {
      limit = max_distance;
    } else {
      limit = max_distance - max_pruning_error;
    }
  } else if (!x.empty()) {
    // Result set X contains only the closest "max_size" items, to within a
    // tolerance of "max_error + max_pruning_error".
    limit = (x.back().first - max_error) - max_pruning_error;
  }

  bool result = true;
  for (const auto& yp : y) {
    // Note that this test also catches duplicate values.
    int count = std::count_if(x.begin(), x.end(), [&yp](const Result& xp) {
        return xp.second == yp.second;
      });
    if (yp.first < limit && count != 1) {
      result = false;
      std::cout << (count > 1 ? "Duplicate" : label) << " distance = "
                << S1ChordAngle(yp.first) << ", id = " << yp.second
                << std::endl;
    }
  }

  return result;
}

}  // namespace internal
}  // namespace S2

template <typename Id, typename Distance>
bool CheckDistanceResults(
    const std::vector<std::pair<Distance, Id>>& expected,
    const std::vector<std::pair<Distance, Id>>& actual,
    int max_size, Distance max_distance, typename Distance::Delta max_error) {
  // This is a conservative bound on the error in computing the distance from
  // the target geometry to an S2Cell.  Such errors can cause candidates to be
  // pruned from the result set even though they may be slightly closer.
  static const typename Distance::Delta kMaxPruningError(
      S1ChordAngle::Radians(1e-15));
  // Use `&` instead of `&&` to evaluate both sides and cast to int to avoid
  // `bitwise-instead-of-logical` warning.
  return (static_cast<int>(S2::internal::CheckResultSet(
              actual, expected, max_size, max_distance, max_error,
              kMaxPruningError, "Missing")) & /*not &&*/
          static_cast<int>(S2::internal::CheckResultSet(
              expected, actual, max_size, max_distance, max_error,
              Distance::Delta::Zero(), "Extra")));
}

#endif  // S2_S2TESTING_H_
