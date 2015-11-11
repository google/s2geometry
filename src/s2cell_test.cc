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

#include "s2cell.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <map>
#include <utility>
#include <vector>

#include <glog/logging.h>

#include "base/macros.h"
#include "base/stringprintf.h"
#include <gtest/gtest.h>
#include "r2.h"
#include "r2rect.h"
#include "s1angle.h"
#include "s1chordangle.h"
#include "s1interval.h"
#include "s2.h"
#include "s2cap.h"
#include "s2edgequery.h"
#include "s2edgeutil.h"
#include "s2latlng.h"
#include "s2latlngrect.h"
#include "s2loop.h"
#include "s2shapeindex.h"
#include "s2testing.h"

using std::map;
using std::max;
using std::min;
using std::vector;

TEST(S2Cell, TestFaces) {
  map<S2Point, int> edge_counts;
  map<S2Point, int> vertex_counts;
  for (int face = 0; face < 6; ++face) {
    S2CellId id = S2CellId::FromFace(face);
    S2Cell cell(id);
    EXPECT_EQ(id, cell.id());
    EXPECT_EQ(face, cell.face());
    EXPECT_EQ(0, cell.level());
    // Top-level faces have alternating orientations to get RHS coordinates.
    EXPECT_EQ(face & S2::kSwapMask, cell.orientation());
    EXPECT_FALSE(cell.is_leaf());
    for (int k = 0; k < 4; ++k) {
      edge_counts[cell.GetEdgeRaw(k)] += 1;
      vertex_counts[cell.GetVertexRaw(k)] += 1;
      EXPECT_DOUBLE_EQ(0.0, cell.GetVertexRaw(k).DotProd(cell.GetEdgeRaw(k)));
      EXPECT_DOUBLE_EQ(0.0,
                       cell.GetVertexRaw((k+1)&3).DotProd(cell.GetEdgeRaw(k)));
      EXPECT_DOUBLE_EQ(1.0,
                       cell.GetVertexRaw(k)
                       .CrossProd(cell.GetVertexRaw((k+1)&3))
                       .Normalize().DotProd(cell.GetEdge(k)));
    }
  }
  // Check that edges have multiplicity 2 and vertices have multiplicity 3.
  for (map<S2Point, int>::iterator i = edge_counts.begin();
       i != edge_counts.end(); ++i) {
    EXPECT_EQ(2, i->second);
  }
  for (map<S2Point, int>::iterator i = vertex_counts.begin();
       i != vertex_counts.end(); ++i) {
    EXPECT_EQ(3, i->second);
  }
}

struct LevelStats {
  double count;
  double min_area, max_area, avg_area;
  double min_width, max_width, avg_width;
  double min_edge, max_edge, avg_edge, max_edge_aspect;
  double min_diag, max_diag, avg_diag, max_diag_aspect;
  double min_angle_span, max_angle_span, avg_angle_span;
  double min_approx_ratio, max_approx_ratio;
  LevelStats()
    : count(0), min_area(100), max_area(0), avg_area(0),
      min_width(100), max_width(0), avg_width(0),
      min_edge(100), max_edge(0), avg_edge(0), max_edge_aspect(0),
      min_diag(100), max_diag(0), avg_diag(0), max_diag_aspect(0),
      min_angle_span(100), max_angle_span(0), avg_angle_span(0),
      min_approx_ratio(100), max_approx_ratio(0) {}
};
static vector<LevelStats> level_stats(S2CellId::kMaxLevel+1);

static void GatherStats(S2Cell const& cell) {
  LevelStats* s = &level_stats[cell.level()];
  double exact_area = cell.ExactArea();
  double approx_area = cell.ApproxArea();
  double min_edge = 100, max_edge = 0, avg_edge = 0;
  double min_diag = 100, max_diag = 0;
  double min_width = 100, max_width = 0;
  double min_angle_span = 100, max_angle_span = 0;
  for (int i = 0; i < 4; ++i) {
    double edge = cell.GetVertexRaw(i).Angle(cell.GetVertexRaw((i+1)&3));
    min_edge = min(edge, min_edge);
    max_edge = max(edge, max_edge);
    avg_edge += 0.25 * edge;
    S2Point mid = cell.GetVertexRaw(i) + cell.GetVertexRaw((i+1)&3);
    double width = M_PI_2 - mid.Angle(cell.GetEdgeRaw(i^2));
    min_width = min(width, min_width);
    max_width = max(width, max_width);
    if (i < 2) {
      double diag = cell.GetVertexRaw(i).Angle(cell.GetVertexRaw(i^2));
      min_diag = min(diag, min_diag);
      max_diag = max(diag, max_diag);
      double angle_span = cell.GetEdgeRaw(i).Angle(-cell.GetEdgeRaw(i^2));
      min_angle_span = min(angle_span, min_angle_span);
      max_angle_span = max(angle_span, max_angle_span);
    }
  }
  s->count += 1;
  s->min_area = min(exact_area, s->min_area);
  s->max_area = max(exact_area, s->max_area);
  s->avg_area += exact_area;
  s->min_width = min(min_width, s->min_width);
  s->max_width = max(max_width, s->max_width);
  s->avg_width += 0.5 * (min_width + max_width);
  s->min_edge = min(min_edge, s->min_edge);
  s->max_edge = max(max_edge, s->max_edge);
  s->avg_edge += avg_edge;
  s->max_edge_aspect = max(max_edge / min_edge, s->max_edge_aspect);
  s->min_diag = min(min_diag, s->min_diag);
  s->max_diag = max(max_diag, s->max_diag);
  s->avg_diag += 0.5 * (min_diag + max_diag);
  s->max_diag_aspect = max(max_diag / min_diag, s->max_diag_aspect);
  s->min_angle_span = min(min_angle_span, s->min_angle_span);
  s->max_angle_span = max(max_angle_span, s->max_angle_span);
  s->avg_angle_span += 0.5 * (min_angle_span + max_angle_span);
  double approx_ratio = approx_area / exact_area;
  s->min_approx_ratio = min(approx_ratio, s->min_approx_ratio);
  s->max_approx_ratio = max(approx_ratio, s->max_approx_ratio);
}

static void TestSubdivide(S2Cell const& cell) {
  GatherStats(cell);
  if (cell.is_leaf()) return;

  S2Cell children[4];
  CHECK(cell.Subdivide(children));
  S2CellId child_id = cell.id().child_begin();
  double exact_area = 0;
  double approx_area = 0;
  double average_area = 0;
  for (int i = 0; i < 4; ++i, child_id = child_id.next()) {
    exact_area += children[i].ExactArea();
    approx_area += children[i].ApproxArea();
    average_area += children[i].AverageArea();

    // Check that the child geometry is consistent with its cell ID.
    EXPECT_EQ(child_id, children[i].id());
    EXPECT_TRUE(S2::ApproxEquals(children[i].GetCenter(), child_id.ToPoint()));
    S2Cell direct(child_id);
    EXPECT_EQ(direct.face(), children[i].face());
    EXPECT_EQ(direct.level(), children[i].level());
    EXPECT_EQ(direct.orientation(), children[i].orientation());
    EXPECT_EQ(direct.GetCenterRaw(), children[i].GetCenterRaw());
    for (int k = 0; k < 4; ++k) {
      EXPECT_EQ(direct.GetVertexRaw(k), children[i].GetVertexRaw(k));
      EXPECT_EQ(direct.GetEdgeRaw(k), children[i].GetEdgeRaw(k));
    }

    // Test Contains() and MayIntersect().
    EXPECT_TRUE(cell.Contains(children[i]));
    EXPECT_TRUE(cell.MayIntersect(children[i]));
    EXPECT_FALSE(children[i].Contains(cell));
    EXPECT_TRUE(cell.Contains(children[i].GetCenterRaw()));
    EXPECT_TRUE(cell.VirtualContainsPoint(children[i].GetCenterRaw()));
    for (int j = 0; j < 4; ++j) {
      EXPECT_TRUE(cell.Contains(children[i].GetVertexRaw(j)));
      if (j != i) {
        EXPECT_FALSE(children[i].Contains(children[j].GetCenterRaw()));
        EXPECT_FALSE(
            children[i].VirtualContainsPoint(children[j].GetCenterRaw()));
        EXPECT_FALSE(children[i].MayIntersect(children[j]));
      }
    }

    // Test GetCapBound and GetRectBound.
    S2Cap parent_cap = cell.GetCapBound();
    S2LatLngRect parent_rect = cell.GetRectBound();
    if (cell.Contains(S2Point(0, 0, 1)) || cell.Contains(S2Point(0, 0, -1))) {
      EXPECT_TRUE(parent_rect.lng().is_full());
    }
    S2Cap child_cap = children[i].GetCapBound();
    S2LatLngRect child_rect = children[i].GetRectBound();
    EXPECT_TRUE(child_cap.Contains(children[i].GetCenter()));
    EXPECT_TRUE(child_rect.Contains(children[i].GetCenterRaw()));
    EXPECT_TRUE(parent_cap.Contains(children[i].GetCenter()));
    EXPECT_TRUE(parent_rect.Contains(children[i].GetCenterRaw()));
    for (int j = 0; j < 4; ++j) {
      EXPECT_TRUE(child_cap.Contains(children[i].GetVertex(j)));
      EXPECT_TRUE(child_rect.Contains(children[i].GetVertex(j)));
      EXPECT_TRUE(child_rect.Contains(children[i].GetVertexRaw(j)));
      EXPECT_TRUE(parent_cap.Contains(children[i].GetVertex(j)));
      EXPECT_TRUE(parent_rect.Contains(children[i].GetVertex(j)));
      EXPECT_TRUE(parent_rect.Contains(children[i].GetVertexRaw(j)));
      if (j != i) {
        // The bounding caps and rectangles should be tight enough so that
        // they exclude at least two vertices of each adjacent cell.
        int cap_count = 0;
        int rect_count = 0;
        for (int k = 0; k < 4; ++k) {
          if (child_cap.Contains(children[j].GetVertex(k)))
            ++cap_count;
          if (child_rect.Contains(children[j].GetVertexRaw(k)))
            ++rect_count;
        }
        EXPECT_LE(cap_count, 2);
        if (child_rect.lat_lo().radians() > -M_PI_2 &&
            child_rect.lat_hi().radians() < M_PI_2) {
          // Bounding rectangles may be too large at the poles because the
          // pole itself has an arbitrary fixed longitude.
          EXPECT_LE(rect_count, 2);
        }
      }
    }

    // Check all children for the first few levels, and then sample randomly.
    // We also always subdivide the cells containing a few chosen points so
    // that we have a better chance of sampling the minimum and maximum metric
    // values.  kMaxSizeUV is the absolute value of the u- and v-coordinate
    // where the cell size at a given level is maximal.
    double const kMaxSizeUV = 0.3964182625366691;
    R2Point const special_uv[] = {
      R2Point(DBL_EPSILON, DBL_EPSILON),  // Face center
      R2Point(DBL_EPSILON, 1),            // Edge midpoint
      R2Point(1, 1),                      // Face corner
      R2Point(kMaxSizeUV, kMaxSizeUV),    // Largest cell area
      R2Point(DBL_EPSILON, kMaxSizeUV),   // Longest edge/diagonal
    };
    bool force_subdivide = false;
    for (int k = 0; k < ARRAYSIZE(special_uv); ++k) {
      if (children[i].GetBoundUV().Contains(special_uv[k]))
        force_subdivide = true;
    }
    if (force_subdivide ||
        cell.level() < (google::DEBUG_MODE ? 5 : 6) ||
        S2Testing::rnd.OneIn(google::DEBUG_MODE ? 5 : 4)) {
      TestSubdivide(children[i]);
    }
  }

  // Check sum of child areas equals parent area.
  //
  // For ExactArea(), the best relative error we can expect is about 1e-6
  // because the precision of the unit vector coordinates is only about 1e-15
  // and the edge length of a leaf cell is about 1e-9.
  //
  // For ApproxArea(), the areas are accurate to within a few percent.
  //
  // For AverageArea(), the areas themselves are not very accurate, but
  // the average area of a parent is exactly 4 times the area of a child.

  EXPECT_LE(fabs(log(exact_area / cell.ExactArea())), fabs(log(1 + 1e-6)));
  EXPECT_LE(fabs(log(approx_area / cell.ApproxArea())), fabs(log(1.03)));
  EXPECT_LE(fabs(log(average_area / cell.AverageArea())), fabs(log(1 + 1e-15)));
}

template <int dim>
static void CheckMinMaxAvg(
    char const* label, int level, double count, double abs_error,
    double min_value, double max_value, double avg_value,
    S2::Metric<dim> const& min_metric,
    S2::Metric<dim> const& max_metric,
    S2::Metric<dim> const& avg_metric) {

  // All metrics are minimums, maximums, or averages of differential
  // quantities, and therefore will not be exact for cells at any finite
  // level.  The differential minimum is always a lower bound, and the maximum
  // is always an upper bound, but these minimums and maximums may not be
  // achieved for two different reasons.  First, the cells at each level are
  // sampled and we may miss the most extreme examples.  Second, the actual
  // metric for a cell is obtained by integrating the differential quantity,
  // which is not constant across the cell.  Therefore cells at low levels
  // (bigger cells) have smaller variations.
  //
  // The "tolerance" below is an attempt to model both of these effects.
  // At low levels, error is dominated by the variation of differential
  // quantities across the cells, while at high levels error is dominated by
  // the effects of random sampling.
  double tolerance = (max_metric.GetValue(level) - min_metric.GetValue(level)) /
                     sqrt(min(count, 0.5 * double(1 << level)));
  if (tolerance == 0) tolerance = abs_error;

  double min_error = min_value - min_metric.GetValue(level);
  double max_error = max_metric.GetValue(level) - max_value;
  double avg_error = fabs(avg_metric.GetValue(level) - avg_value);
  printf("%-10s (%6.0f samples, tolerance %8.3g) - min (%9.3g : %9.3g) "
         "max (%9.3g : %9.3g), avg (%9.3g : %9.3g)\n",
         label, count, tolerance,
         min_error / min_value, min_error / tolerance,
         max_error / max_value, max_error / tolerance,
         avg_error / avg_value, avg_error / tolerance);

  EXPECT_LE(min_metric.GetValue(level), min_value + abs_error);
  EXPECT_GE(min_metric.GetValue(level), min_value - tolerance);
  EXPECT_LE(max_metric.GetValue(level), max_value + tolerance);
  EXPECT_GE(max_metric.GetValue(level), max_value - abs_error);
  EXPECT_NEAR(avg_metric.GetValue(level), avg_value, 10 * tolerance);
}

TEST(S2Cell, TestSubdivide) {
  // Only test a sample of faces to reduce the runtime.
  TestSubdivide(S2Cell::FromFace(0));
  TestSubdivide(S2Cell::FromFace(3));
  TestSubdivide(S2Cell::FromFace(5));

  // The maximum edge *ratio* is the ratio of the longest edge of any cell to
  // the shortest edge of any cell at the same level (and similarly for the
  // maximum diagonal ratio).
  //
  // The maximum edge *aspect* is the maximum ratio of the longest edge of a
  // cell to the shortest edge of that same cell (and similarly for the
  // maximum diagonal aspect).

  printf("Level    Area      Edge          Diag          Approx       Average\n");
  printf("        Ratio  Ratio Aspect  Ratio Aspect    Min    Max    Min    Max\n");
  for (int i = 0; i <= S2CellId::kMaxLevel; ++i) {
    LevelStats* s = &level_stats[i];
    if (s->count > 0) {
      s->avg_area /= s->count;
      s->avg_width /= s->count;
      s->avg_edge /= s->count;
      s->avg_diag /= s->count;
      s->avg_angle_span /= s->count;
    }
    printf("%5d  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
           i, s->max_area / s->min_area,
           s->max_edge / s->min_edge, s->max_edge_aspect,
           s->max_diag / s->min_diag, s->max_diag_aspect,
           s->min_approx_ratio, s->max_approx_ratio,
           S2Cell::AverageArea(i) / s->max_area,
           S2Cell::AverageArea(i) / s->min_area);
  }

  // Now check the validity of the S2 length and area metrics.
  for (int i = 0; i <= S2CellId::kMaxLevel; ++i) {
    LevelStats const* s = &level_stats[i];
    if (s->count == 0) continue;

    printf("Level %2d - metric (error/actual : error/tolerance)\n", i);

    // The various length calculations are only accurate to 1e-15 or so,
    // so we need to allow for this amount of discrepancy with the theoretical
    // minimums and maximums.  The area calculation is accurate to about 1e-15
    // times the cell width.
    CheckMinMaxAvg("area", i, s->count, 1e-15 * s->min_width,
                   s->min_area, s->max_area, s->avg_area,
                   S2::kMinArea, S2::kMaxArea, S2::kAvgArea);
    CheckMinMaxAvg("width", i, s->count, 1e-15,
                   s->min_width, s->max_width, s->avg_width,
                   S2::kMinWidth, S2::kMaxWidth, S2::kAvgWidth);
    CheckMinMaxAvg("edge", i, s->count, 1e-15,
                   s->min_edge, s->max_edge, s->avg_edge,
                   S2::kMinEdge, S2::kMaxEdge, S2::kAvgEdge);
    CheckMinMaxAvg("diagonal", i, s->count, 1e-15,
                   s->min_diag, s->max_diag, s->avg_diag,
                   S2::kMinDiag, S2::kMaxDiag, S2::kAvgDiag);
    CheckMinMaxAvg("angle span", i, s->count, 1e-15,
                   s->min_angle_span, s->max_angle_span, s->avg_angle_span,
                   S2::kMinAngleSpan, S2::kMaxAngleSpan, S2::kAvgAngleSpan);

    // The aspect ratio calculations are ratios of lengths and are therefore
    // less accurate at higher subdivision levels.
    EXPECT_LE(s->max_edge_aspect, S2::kMaxEdgeAspect + 1e-15 * (1 << i));
    EXPECT_LE(s->max_diag_aspect, S2::kMaxDiagAspect + 1e-15 * (1 << i));
  }
}

static int const kMaxLevel = google::DEBUG_MODE ? 6 : 11;

static void ExpandChildren1(S2Cell const& cell) {
  S2Cell children[4];
  CHECK(cell.Subdivide(children));
  if (children[0].level() < kMaxLevel) {
    for (int pos = 0; pos < 4; ++pos) {
      ExpandChildren1(children[pos]);
    }
  }
}

static void ExpandChildren2(S2Cell const& cell) {
  S2CellId id = cell.id().child_begin();
  for (int pos = 0; pos < 4; ++pos, id = id.next()) {
    S2Cell child(id);
    if (child.level() < kMaxLevel) ExpandChildren2(child);
  }
}

TEST(S2Cell, TestPerformance) {
  double subdivide_start = S2Testing::GetCpuTime();
  ExpandChildren1(S2Cell::FromFace(0));
  double subdivide_time = S2Testing::GetCpuTime() - subdivide_start;
  fprintf(stderr, "Subdivide: %.3f seconds\n", subdivide_time);

  double constructor_start = S2Testing::GetCpuTime();
  ExpandChildren2(S2Cell::FromFace(0));
  double constructor_time = S2Testing::GetCpuTime() - constructor_start;
  fprintf(stderr, "Constructor: %.3f seconds\n", constructor_time);
}

TEST(S2Cell, CellVsLoopRectBound) {
  // This test verifies that the S2Cell and S2Loop bounds contain each other
  // to within their maximum errors.
  //
  // The S2Cell and S2Loop calculations for the latitude of a vertex can differ
  // by up to 2 * DBL_EPSILON, therefore the S2Cell bound should never exceed
  // the S2Loop bound by more than this (the reverse is not true, because the
  // S2Loop code sometimes thinks that the maximum occurs along an edge).
  // Similarly, the longitude bounds can differ by up to 4 * DBL_EPSILON since
  // the S2Cell bound has an error of 2 * DBL_EPSILON and then expands by this
  // amount, while the S2Loop bound does no expansion at all.

  // Possible additional S2Cell error compared to S2Loop error:
  static S2LatLng kCellError = S2LatLng::FromRadians(2 * DBL_EPSILON,
                                                     4 * DBL_EPSILON);
  // Possible additional S2Loop error compared to S2Cell error:
  static S2LatLng kLoopError = S2EdgeUtil::RectBounder::MaxErrorForTests();

  for (int iter = 0; iter < 1000; ++iter) {
    S2Cell cell(S2Testing::GetRandomCellId());
    S2Loop loop(cell);
    S2LatLngRect cell_bound = cell.GetRectBound();
    S2LatLngRect loop_bound = loop.GetRectBound();
    EXPECT_TRUE(loop_bound.Expanded(kCellError).Contains(cell_bound));
    EXPECT_TRUE(cell_bound.Expanded(kLoopError).Contains(loop_bound));
  }
}

TEST(S2Cell, RectBoundIsLargeEnough) {
  // Construct many points that are nearly on an S2Cell edge, and verify that
  // whenever the cell contains a point P then its bound contains S2LatLng(P).

  for (int iter = 0; iter < 1000; /* advanced in loop below */) {
    S2Cell cell(S2Testing::GetRandomCellId());
    int i1 = S2Testing::rnd.Uniform(4);
    int i2 = (i1 + 1) & 3;
    S2Point v1 = cell.GetVertex(i1);
    S2Point v2 = S2Testing::SamplePoint(
        S2Cap(cell.GetVertex(i2), S1Angle::Radians(1e-15)));
    S2Point p = S2EdgeUtil::Interpolate(S2Testing::rnd.RandDouble(), v1, v2);
    if (S2Loop(cell).Contains(p)) {
      EXPECT_TRUE(cell.GetRectBound().Contains(S2LatLng(p)));
      ++iter;
    }
  }
}

static S1ChordAngle GetDistanceToPointBruteForce(S2Cell const& cell,
                                                 S2Point const& target) {
  if (cell.Contains(target)) return S1ChordAngle::Zero();
  S1ChordAngle min_distance = S1ChordAngle::Infinity();
  for (int i = 0; i < 4; ++i) {
    S2EdgeUtil::UpdateMinDistance(target, cell.GetVertex(i),
                                  cell.GetVertex((i + 1) & 3), &min_distance);
  }
  return min_distance;
}

TEST(S2Cell, GetDistanceToPoint) {
  S2Testing::rnd.Reset(FLAGS_s2_random_seed);
  for (int iter = 0; iter < 1000; ++iter) {
    SCOPED_TRACE(StringPrintf("Iteration %d", iter));
    S2Cell cell(S2Testing::GetRandomCellId());
    S2Point target = S2Testing::RandomPoint();
    S1Angle expected = GetDistanceToPointBruteForce(cell, target).ToAngle();
    S1Angle actual = cell.GetDistance(target).ToAngle();
    // The error has a peak near Pi/2 for edge distance, and another peak near
    // Pi for vertex distance.
    EXPECT_NEAR(expected.radians(), actual.radians(), 1e-12);
    if (expected.radians() <= M_PI / 3) {
      EXPECT_NEAR(expected.radians(), actual.radians(), 1e-15);
    }
  }
}

static void ChooseEdgeNearCell(S2Cell const& cell, S2Point* a, S2Point* b) {
  S2Cap cap = cell.GetCapBound();
  if (S2Testing::rnd.OneIn(5)) {
    // Choose a point anywhere on the sphere.
    *a = S2Testing::RandomPoint();
  } else {
    // Choose a point inside or somewhere near the cell.
    *a = S2Testing::SamplePoint(S2Cap(cap.center(), 1.5 * cap.GetRadius()));
  }
  // Now choose a maximum edge length ranging from very short to very long
  // relative to the cell size, and choose the other endpoint.
  double max_length = min(100 * pow(1e-4, S2Testing::rnd.RandDouble()) *
                          cap.GetRadius().radians(), M_PI_2);
  *b = S2Testing::SamplePoint(S2Cap(*a, S1Angle::Radians(max_length)));
}

static S1ChordAngle GetDistanceToEdgeBruteForce(
    S2Cell const& cell, S2Point const& a, S2Point const& b) {
  if (cell.Contains(a) || cell.Contains(b)) {
    return S1ChordAngle::Zero();
  }
  S2Loop loop(cell);
  S2ShapeIndex index;
  index.Add(new S2Loop::Shape(&loop));
  S2EdgeQuery query(index);
  vector<int> edges;
  if (query.GetCrossings(a, b, index.shape(0), &edges)) {
    return S1ChordAngle::Zero();
  }
  S1ChordAngle min_dist = S1ChordAngle::Infinity();
  for (int i = 0; i < 4; ++i) {
    S2Point v0 = cell.GetVertex(i);
    S2Point v1 = cell.GetVertex((i + 1) & 3);
    S2EdgeUtil::UpdateMinDistance(a, v0, v1, &min_dist);
    S2EdgeUtil::UpdateMinDistance(b, v0, v1, &min_dist);
    S2EdgeUtil::UpdateMinDistance(v0, a, b, &min_dist);
  }
  return min_dist;
}

TEST(S2Cell, GetDistanceToEdge) {
  S2Testing::rnd.Reset(FLAGS_s2_random_seed);
  for (int iter = 0; iter < 1000; ++iter) {
    SCOPED_TRACE(StringPrintf("Iteration %d", iter));
    S2Cell cell(S2Testing::GetRandomCellId());
    S2Point a, b;
    ChooseEdgeNearCell(cell, &a, &b);
    S1Angle expected = GetDistanceToEdgeBruteForce(cell, a, b).ToAngle();
    S1Angle actual = cell.GetDistanceToEdge(a, b).ToAngle();
    // The error has a peak near Pi/2 for edge distance, and another peak near
    // Pi for vertex distance.
    EXPECT_NEAR(expected.radians(), actual.radians(), 1e-12);
    if (expected.radians() <= M_PI / 3) {
      EXPECT_NEAR(expected.radians(), actual.radians(), 1e-15);
    }
  }
}
