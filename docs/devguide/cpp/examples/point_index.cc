// Copyright 2017 Google Inc. All Rights Reserved.
// Author: ericv@google.com (Eric Veach)
//
// This example shows how to build and query an in-memory index of points
// using S2PointIndex.

#include <cinttypes>
#include <cstdint>
#include <cstdio>

#include "base/commandlineflags.h"
#include "base/init_google.h"  // MOE:strip_line
#include "geostore/geomodel/public/s2earth.h"
#include "third_party/absl/flags/flag.h"
#include "util/geometry/s1chord_angle.h"
#include "util/geometry/s2closest_point_query.h"
#include "util/geometry/s2point_index.h"
#include "util/geometry/s2testing.h"

DEFINE_int32(num_index_points, 10000, "Number of points to index");
DEFINE_int32(num_queries, 10000, "Number of queries");
DEFINE_double(query_radius_km, 100, "Query radius in kilometers");

// MOE:begin_strip
using geostore::S2Earth;

// MOE:end_strip
int main(int argc, char **argv) {
  // MOE:begin_strip
  // TODO(jrosenstock,b/209396162): We should call
  // gflags::ParseCommandLineFlags here, but that would break WITH_FLAGS=off.
  InitGoogle(argv[0], &argc, &argv, true);
  // MOE:end_strip
  // Build an index containing random points anywhere on the Earth.
  S2PointIndex<int> index;
  for (int i = 0; i < absl::GetFlag(FLAGS_num_index_points); ++i) {
    index.Add(S2Testing::RandomPoint(), i);
  }

  // Create a query to search within the given radius of a target point.
  S2ClosestPointQuery<int> query(&index);
  query.mutable_options()->set_max_distance(S1Angle::Radians(
      S2Earth::KmToRadians(absl::GetFlag(FLAGS_query_radius_km))));

  // Repeatedly choose a random target point, and count how many index points
  // are within the given radius of that point.
  int64_t num_found = 0;
  for (int i = 0; i < absl::GetFlag(FLAGS_num_queries); ++i) {
    S2ClosestPointQuery<int>::PointTarget target(S2Testing::RandomPoint());
    num_found += query.FindClosestPoints(&target).size();
  }

  std::printf("Found %" PRId64 " points in %d queries\n", num_found,
              absl::GetFlag(FLAGS_num_queries));
  return 0;
}
