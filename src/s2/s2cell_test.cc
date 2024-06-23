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

#include "s2/s2cell.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <string>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include "absl/container/flat_hash_map.h"
#include "absl/flags/flag.h"
#include "absl/log/absl_check.h"
#include "absl/log/log_streamer.h"
#include "absl/random/bit_gen_ref.h"
#include "absl/random/random.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_format.h"
#include "absl/strings/string_view.h"

#include "s2/base/log_severity.h"
#include "s2/r2.h"
#include "s2/r2rect.h"
#include "s2/s1angle.h"
#include "s2/s1chord_angle.h"
#include "s2/s1interval.h"
#include "s2/s2cap.h"
#include "s2/s2cell_id.h"
#include "s2/s2edge_crossings.h"
#include "s2/s2edge_distances.h"
#include "s2/s2latlng.h"
#include "s2/s2latlng_rect.h"
#include "s2/s2latlng_rect_bounder.h"
#include "s2/s2loop.h"
#include "s2/s2metrics.h"
#include "s2/s2point.h"
#include "s2/s2pointutil.h"
#include "s2/s2random.h"
#include "s2/s2testing.h"
#include "s2/s2text_format.h"
#include "s2/util/coding/coder.h"

using absl::flat_hash_map;
using absl::StrCat;
using absl::string_view;
using S2::internal::kSwapMask;
using s2textformat::MakePointOrDie;
using std::fabs;
using std::max;
using std::min;
using std::pow;
using std::vector;
using ::testing::Eq;

// Reflects the center point of a cell across the given boundary.
S2Point ReflectCenter(const S2Cell& cell, int k) {
  Vector3_d normal = cell.GetEdgeRaw(k);
  return Matrix3x3_d::Householder(normal) * cell.GetCenter();
}

TEST(S2Cell, TestFaces) {
  flat_hash_map<S2Point, int> edge_counts;
  flat_hash_map<S2Point, int> vertex_counts;
  for (int face = 0; face < 6; ++face) {
    S2CellId id = S2CellId::FromFace(face);
    S2Cell cell(id);
    EXPECT_EQ(id, cell.id());
    EXPECT_EQ(face, cell.face());
    EXPECT_EQ(0, cell.level());
    // Top-level faces have alternating orientations to get RHS coordinates.
    EXPECT_EQ(face & kSwapMask, cell.orientation());
    EXPECT_FALSE(cell.is_leaf());
    for (int k = 0; k < 4; ++k) {
      edge_counts[cell.GetEdgeRaw(k)] += 1;
      vertex_counts[cell.GetVertexRaw(k)] += 1;
      EXPECT_DOUBLE_EQ(0.0, cell.GetVertexRaw(k).DotProd(cell.GetEdgeRaw(k)));
      EXPECT_DOUBLE_EQ(0.0,
                       cell.GetVertexRaw(k + 1).DotProd(cell.GetEdgeRaw(k)));
      EXPECT_DOUBLE_EQ(1.0,
                       cell.GetVertexRaw(k).CrossProd(cell.GetVertexRaw(k + 1)).
                       Normalize().DotProd(cell.GetEdge(k)));
    }
  }
  // Check that edges have multiplicity 2 and vertices have multiplicity 3.
  for (const auto& p : edge_counts) {
    EXPECT_EQ(2, p.second);
  }
  for (const auto& p : vertex_counts) {
    EXPECT_EQ(3, p.second);
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

static void GatherStats(const S2Cell& cell) {
  LevelStats* s = &level_stats[cell.level()];
  double exact_area = cell.ExactArea();
  double approx_area = cell.ApproxArea();
  double min_edge = 100, max_edge = 0, avg_edge = 0;
  double min_diag = 100, max_diag = 0;
  double min_width = 100, max_width = 0;
  double min_angle_span = 100, max_angle_span = 0;
  for (int i = 0; i < 4; ++i) {
    double edge = cell.GetVertexRaw(i).Angle(cell.GetVertexRaw(i + 1));
    min_edge = min(edge, min_edge);
    max_edge = max(edge, max_edge);
    avg_edge += 0.25 * edge;
    S2Point mid = cell.GetVertexRaw(i) + cell.GetVertexRaw(i + 1);
    double width = M_PI_2 - mid.Angle(cell.GetEdgeRaw(i + 2));
    min_width = min(width, min_width);
    max_width = max(width, max_width);
    if (i < 2) {
      double diag = cell.GetVertexRaw(i).Angle(cell.GetVertexRaw(i + 2));
      min_diag = min(diag, min_diag);
      max_diag = max(diag, max_diag);
      double angle_span = cell.GetEdgeRaw(i).Angle(-cell.GetEdgeRaw(i + 2));
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

static void TestSubdivide(absl::BitGenRef bitgen, const S2Cell& cell) {
  GatherStats(cell);
  if (cell.is_leaf()) return;

  S2Cell children[4];
  ABSL_CHECK(cell.Subdivide(children));
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
    for (int j = 0; j < 4; ++j) {
      EXPECT_TRUE(cell.Contains(children[i].GetVertexRaw(j)));
      if (j != i) {
        EXPECT_FALSE(children[i].Contains(children[j].GetCenterRaw()));
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
    const double kMaxSizeUV = 0.3964182625366691;
    const R2Point special_uv[] = {
      R2Point(DBL_EPSILON, DBL_EPSILON),  // Face center
      R2Point(DBL_EPSILON, 1),            // Edge midpoint
      R2Point(1, 1),                      // Face corner
      R2Point(kMaxSizeUV, kMaxSizeUV),    // Largest cell area
      R2Point(DBL_EPSILON, kMaxSizeUV),   // Longest edge/diagonal
    };
    bool force_subdivide = false;
    for (const R2Point& uv : special_uv) {
      if (children[i].GetBoundUV().Contains(uv))
        force_subdivide = true;
    }
    if (force_subdivide || cell.level() < (S2_DEBUG_MODE ? 5 : 6) ||
        absl::Bernoulli(bitgen, S2_DEBUG_MODE ? 0.2 : 0.25)) {
      TestSubdivide(bitgen, children[i]);
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
static void CheckMinMaxAvg(string_view label, int level, double count,
                           double abs_error, double min_value, double max_value,
                           double avg_value, const S2::Metric<dim>& min_metric,
                           const S2::Metric<dim>& max_metric,
                           const S2::Metric<dim>& avg_metric) {
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
  absl::PrintF(
      "%-10s (%6.0f samples, tolerance %8.3g) - min %9.4g (%9.3g : %9.3g) "
      "max %9.4g (%9.3g : %9.3g), avg %9.4g (%9.3g : %9.3g)\n",
      label, count, tolerance, min_value, min_error / min_value,
      min_error / tolerance, max_value, max_error / max_value,
      max_error / tolerance, avg_value, avg_error / avg_value,
      avg_error / tolerance);

  EXPECT_LE(min_metric.GetValue(level), min_value + abs_error);
  EXPECT_GE(min_metric.GetValue(level), min_value - tolerance);
  EXPECT_LE(max_metric.GetValue(level), max_value + tolerance);
  EXPECT_GE(max_metric.GetValue(level), max_value - abs_error);
  EXPECT_NEAR(avg_metric.GetValue(level), avg_value, 10 * tolerance);
}

TEST(S2Cell, TestSubdivide) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "TEST_SUBDIVIDE",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  // Only test a sample of faces to reduce the runtime.
  TestSubdivide(bitgen, S2Cell::FromFace(0));
  TestSubdivide(bitgen, S2Cell::FromFace(3));
  TestSubdivide(bitgen, S2Cell::FromFace(5));

  // This table is useful in evaluating the quality of the various S2
  // projections.
  //
  // The maximum edge *ratio* is the ratio of the longest edge of any cell to
  // the shortest edge of any cell at the same level (and similarly for the
  // maximum diagonal ratio).
  //
  // The maximum edge *aspect* is the maximum ratio of the longest edge of a
  // cell to the shortest edge of that same cell (and similarly for the
  // maximum diagonal aspect).
  absl::PrintF(
      "Ratio:  (Max value for any cell) / (Min value for any cell)\n"
      "Aspect: (Max value / min value) for any cell\n"
      "                   Edge          Diag       Approx Area/    Avg Area/\n"
      "         Area     Length        Length       Exact Area    Exact Area\n"
      "Level   Ratio  Ratio Aspect  Ratio Aspect    Min    Max    Min    Max\n"
      "--------------------------------------------------------------------\n");
  for (int i = 0; i <= S2CellId::kMaxLevel; ++i) {
    LevelStats* s = &level_stats[i];
    if (s->count > 0) {
      s->avg_area /= s->count;
      s->avg_width /= s->count;
      s->avg_edge /= s->count;
      s->avg_diag /= s->count;
      s->avg_angle_span /= s->count;
    }
    absl::PrintF("%5d  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
                 i, s->max_area / s->min_area, s->max_edge / s->min_edge,
                 s->max_edge_aspect, s->max_diag / s->min_diag,
                 s->max_diag_aspect, s->min_approx_ratio, s->max_approx_ratio,
                 S2Cell::AverageArea(i) / s->max_area,
                 S2Cell::AverageArea(i) / s->min_area);
  }

  // Now check the validity of the S2 length and area metrics.
  for (int i = 0; i <= S2CellId::kMaxLevel; ++i) {
    const LevelStats* s = &level_stats[i];
    if (s->count == 0) continue;

    absl::PrintF("Level %2d - metric value (error/actual : error/tolerance)\n",
                 i);

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
  static S2LatLng kLoopError = S2LatLngRectBounder::MaxErrorForTests();

  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "CELL_VS_LOOP_RECT_BOUND",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int iter = 0; iter < 1000; ++iter) {
    S2Cell cell(s2random::CellId(bitgen));
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
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "RECT_BOUND_IS_LARGE_ENOUGH",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  for (int iter = 0; iter < 1000; /* advanced in loop below */) {
    S2Cell cell(s2random::CellId(bitgen));
    int i = absl::Uniform(bitgen, 0, 4);
    S2Point v1 = cell.GetVertex(i);
    S2Point v2 = s2random::SamplePoint(
        bitgen, S2Cap(cell.GetVertex(i + 1), S1Angle::Radians(1e-15)));
    S2Point p = S2::Interpolate(v1, v2, absl::Uniform(bitgen, 0.0, 1.0));
    if (S2Loop(cell).Contains(p)) {
      EXPECT_TRUE(cell.GetRectBound().Contains(S2LatLng(p)));
      ++iter;
    }
  }
}

TEST(S2Cell, ConsistentWithS2CellIdFromPoint) {
  // Construct many points that are nearly on an S2Cell edge, and verify that
  // S2Cell(S2CellId(p)).Contains(p) is always true.
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "CONSISTENT_WITH_S2CELL_ID_FROM_POINT",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  for (int iter = 0; iter < 1000; ++iter) {
    S2Cell cell(s2random::CellId(bitgen));
    int i = absl::Uniform(bitgen, 0, 4);
    S2Point v1 = cell.GetVertex(i);
    S2Point v2 = s2random::SamplePoint(
        bitgen, S2Cap(cell.GetVertex(i + 1), S1Angle::Radians(1e-15)));
    S2Point p = S2::Interpolate(v1, v2, absl::Uniform(bitgen, 0.0, 1.0));
    EXPECT_TRUE(S2Cell(S2CellId(p)).Contains(p));
  }
}

TEST(S2CellId, AmbiguousContainsPoint) {
  // This tests a case where S2CellId returns the "wrong" cell for a point
  // that is very close to the cell edge. (ConsistentWithS2CellIdFromPoint
  // generates more examples like this.)
  //
  // The S2Point below should have x = 0, but conversion from latlng to
  // (x,y,z) gives x = 6.1e-17.  When xyz is converted to uv, this gives u =
  // -6.1e-17.  However when converting to st, which is centered at 0.5 rather
  // than 0, the low precision bits of u are lost and we wind up with s = 0.5.
  // S2CellId(const S2Point&) then chooses an arbitrary neighboring cell.
  //
  // This tests that S2Cell::Contains() expands the cell bounds sufficiently
  // so that the returned cell is still considered to contain "p".
  S2Point p = S2LatLng::FromDegrees(-2, 90).ToPoint();
  S2CellId cell_id = S2CellId(p).parent(1);
  S2Cell cell(cell_id);
  EXPECT_TRUE(cell.Contains(p));
}

static S1ChordAngle GetDistanceToPointBruteForce(const S2Cell& cell,
                                                 const S2Point& target) {
  S1ChordAngle min_distance = S1ChordAngle::Infinity();
  for (int i = 0; i < 4; ++i) {
    S2::UpdateMinDistance(target, cell.GetVertex(i),
                                  cell.GetVertex(i + 1), &min_distance);
  }
  return min_distance;
}

static S1ChordAngle GetMaxDistanceToPointBruteForce(const S2Cell& cell,
                                                    const S2Point& target) {
  if (cell.Contains(-target)) {
    return S1ChordAngle::Straight();
  }
  S1ChordAngle max_distance = S1ChordAngle::Negative();
  for (int i = 0; i < 4; ++i) {
    S2::UpdateMaxDistance(target, cell.GetVertex(i),
                                  cell.GetVertex(i + 1), &max_distance);
  }
  return max_distance;
}

TEST(S2Cell, GetDistanceToPoint) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "GET_DISTANCE_TO_POINT",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int iter = 0; iter < 1000; ++iter) {
    SCOPED_TRACE(StrCat("Iteration ", iter));
    S2Cell cell(s2random::CellId(bitgen));
    S2Point target = s2random::Point(bitgen);
    S1Angle expected_to_boundary =
        GetDistanceToPointBruteForce(cell, target).ToAngle();
    S1Angle expected_to_interior =
        cell.Contains(target) ? S1Angle::Zero() : expected_to_boundary;
    S1Angle expected_max =
        GetMaxDistanceToPointBruteForce(cell, target).ToAngle();
    S1Angle actual_to_boundary = cell.GetBoundaryDistance(target).ToAngle();
    S1Angle actual_to_interior = cell.GetDistance(target).ToAngle();
    S1Angle actual_max = cell.GetMaxDistance(target).ToAngle();
    // The error has a peak near Pi/2 for edge distance, and another peak near
    // Pi for vertex distance.
    EXPECT_NEAR(expected_to_boundary.radians(),
                actual_to_boundary.radians(), 1e-12);
    EXPECT_NEAR(expected_to_interior.radians(),
                actual_to_interior.radians(), 1e-12);
    EXPECT_NEAR(expected_max.radians(),
                actual_max.radians(), 1e-12);
    if (expected_to_boundary.radians() <= M_PI / 3) {
      EXPECT_NEAR(expected_to_boundary.radians(),
                  actual_to_boundary.radians(), 1e-15);
      EXPECT_NEAR(expected_to_interior.radians(),
                  actual_to_interior.radians(), 1e-15);
    }
    if (expected_max.radians() <= M_PI / 3) {
      EXPECT_NEAR(expected_max.radians(), actual_max.radians(), 1e-15);
    }
  }
}

static void ChooseEdgeNearCell(absl::BitGenRef bitgen, const S2Cell& cell,
                               S2Point* a, S2Point* b) {
  S2Cap cap = cell.GetCapBound();
  if (absl::Bernoulli(bitgen, 0.2)) {
    // Choose a point anywhere on the sphere.
    *a = s2random::Point(bitgen);
  } else {
    // Choose a point inside or somewhere near the cell.
    *a = s2random::SamplePoint(bitgen,
                               S2Cap(cap.center(), 1.5 * cap.GetRadius()));
  }
  // Now choose a maximum edge length ranging from very short to very long
  // relative to the cell size, and choose the other endpoint.
  double max_length =
      min(s2random::LogUniform(bitgen, 1e-2, 1e2) * cap.GetRadius().radians(),
          M_PI_2);
  *b = s2random::SamplePoint(bitgen, S2Cap(*a, S1Angle::Radians(max_length)));

  if (absl::Bernoulli(bitgen, 0.05)) {
    // Occasionally replace edge with antipodal edge.
    *a = -*a;
    *b = -*b;
  }
}

static S1ChordAngle GetDistanceToEdgeBruteForce(
    const S2Cell& cell, const S2Point& a, const S2Point& b) {
  if (cell.Contains(a) || cell.Contains(b)) {
    return S1ChordAngle::Zero();
  }

  S1ChordAngle min_dist = S1ChordAngle::Infinity();
  for (int i = 0; i < 4; ++i) {
    S2Point v0 = cell.GetVertex(i);
    S2Point v1 = cell.GetVertex(i + 1);
    // If the edge crosses through the cell, max distance is 0.
    if (S2::CrossingSign(a, b, v0, v1) >= 0) {
      return S1ChordAngle::Zero();
    }
    S2::UpdateMinDistance(a, v0, v1, &min_dist);
    S2::UpdateMinDistance(b, v0, v1, &min_dist);
    S2::UpdateMinDistance(v0, a, b, &min_dist);
  }
  return min_dist;
}

static S1ChordAngle GetMaxDistanceToEdgeBruteForce(
    const S2Cell& cell, const S2Point& a, const S2Point& b) {
  // If any antipodal endpoint is within the cell, the max distance is Pi.
  if (cell.Contains(-a) || cell.Contains(-b)) {
    return S1ChordAngle::Straight();
  }

  S1ChordAngle max_dist = S1ChordAngle::Negative();
  for (int i = 0; i < 4; ++i) {
    S2Point v0 = cell.GetVertex(i);
    S2Point v1 = cell.GetVertex(i + 1);
    // If the antipodal edge crosses through the cell, max distance is Pi.
    if (S2::CrossingSign(-a, -b, v0, v1) >= 0) {
      return S1ChordAngle::Straight();
    }
    S2::UpdateMaxDistance(a, v0, v1, &max_dist);
    S2::UpdateMaxDistance(b, v0, v1, &max_dist);
    S2::UpdateMaxDistance(v0, a, b, &max_dist);
  }
  return max_dist;
}

TEST(S2Cell, GetDistanceToEdge) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "GET_DISTANCE_TO_EDGE",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int iter = 0; iter < 1000; ++iter) {
    SCOPED_TRACE(StrCat("Iteration ", iter));
    S2Cell cell(s2random::CellId(bitgen));
    S2Point a, b;
    ChooseEdgeNearCell(bitgen, cell, &a, &b);
    S1Angle expected_min = GetDistanceToEdgeBruteForce(cell, a, b).ToAngle();
    S1Angle expected_max =
        GetMaxDistanceToEdgeBruteForce(cell, a, b).ToAngle();
    S1Angle actual_min = cell.GetDistance(a, b).ToAngle();
    S1Angle actual_max = cell.GetMaxDistance(a, b).ToAngle();
    // The error has a peak near Pi/2 for edge distance, and another peak near
    // Pi for vertex distance.
    if (expected_min.radians() > M_PI/2) {
      // Max error for S1ChordAngle as it approaches Pi is about 3e-8.
      EXPECT_NEAR(expected_min.radians(), actual_min.radians(), 3e-8);
    } else if (expected_min.radians() <= M_PI / 3) {
      EXPECT_NEAR(expected_min.radians(), actual_min.radians(), 1e-15);
    } else {
      EXPECT_NEAR(expected_min.radians(), actual_min.radians(), 1e-12);
    }

    EXPECT_NEAR(expected_max.radians(), actual_max.radians(), 1e-12);
    if (expected_max.radians() <= M_PI / 3) {
      EXPECT_NEAR(expected_max.radians(), actual_max.radians(), 1e-15);
    }
  }
}

TEST(S2Cell, GetMaxDistanceToEdge) {
  // Test an edge for which its antipode crosses the cell. Validates both the
  // standard and brute force implementations for this case.
  S2Cell cell = S2Cell::FromFacePosLevel(0, 0, 20);
  S2Point a = -S2::Interpolate(cell.GetCenter(), cell.GetVertex(0), 2.0);
  S2Point b = -S2::Interpolate(cell.GetCenter(), cell.GetVertex(2), 2.0);

  S1ChordAngle actual = cell.GetMaxDistance(a, b);
  S1ChordAngle expected = GetMaxDistanceToEdgeBruteForce(cell, a, b);

  EXPECT_NEAR(expected.radians(), S1ChordAngle::Straight().radians(), 1e-15);
  EXPECT_NEAR(actual.radians(), S1ChordAngle::Straight().radians(), 1e-15);
}

TEST(S2Cell, GetMaxDistanceToCellAntipodal) {
  S2Point p = MakePointOrDie("0:0");
  S2Cell cell(p);
  S2Cell antipodal_cell(-p);
  S1ChordAngle dist = cell.GetMaxDistance(antipodal_cell);
  EXPECT_EQ(S1ChordAngle::Straight(), dist);
}

TEST(S2Cell, GetMaxDistanceToCell) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "GET_MAX_DISTANCE_TO_CELL",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < 1000; i++) {
    S2Cell cell(s2random::CellId(bitgen));
    S2Cell test_cell(s2random::CellId(bitgen));
    S2CellId antipodal_leaf_id(-test_cell.GetCenter());
    S2Cell antipodal_test_cell(antipodal_leaf_id.parent(test_cell.level()));

    S1ChordAngle dist_from_min = S1ChordAngle::Straight() -
        cell.GetDistance(antipodal_test_cell);
    S1ChordAngle dist_from_max = cell.GetMaxDistance(test_cell);
    EXPECT_NEAR(dist_from_min.radians(), dist_from_max.radians(), 1e-8);
  }
}

TEST(S2Cell, EncodeDecode) {
  S2Cell orig_cell(S2LatLng::FromDegrees(40.7406264, -74.0029963));
  Encoder encoder;
  orig_cell.Encode(&encoder);

  S2Cell decoded_cell(S2LatLng::FromDegrees(51.494987, -0.146585));
  Decoder decoder(encoder.base(), encoder.length());
  ASSERT_TRUE(decoded_cell.Decode(&decoder));

  EXPECT_EQ(orig_cell, decoded_cell);
  EXPECT_EQ(orig_cell.face(), decoded_cell.face());
  EXPECT_EQ(orig_cell.level(), decoded_cell.level());
  EXPECT_EQ(orig_cell.orientation(), decoded_cell.orientation());
  EXPECT_EQ(orig_cell.id(), decoded_cell.id());
  EXPECT_EQ(orig_cell.GetBoundUV(), decoded_cell.GetBoundUV());
}

TEST(S2Cell, GetUVCoordOfEdge) {
  // Four cells on face 0 with two boundaries each on 0/0.
  S2Cell kCell0[4] = {
      S2Cell(S2CellId::FromToken("0f")), S2Cell(S2CellId::FromToken("05")),
      S2Cell(S2CellId::FromToken("1b")), S2Cell(S2CellId::FromToken("11"))};

  // And four cells on face 4 which is rotated w.r.t face 0.
  S2Cell kCell4[4] = {
      S2Cell(S2CellId::FromToken("8f")), S2Cell(S2CellId::FromToken("85")),
      S2Cell(S2CellId::FromToken("9b")), S2Cell(S2CellId::FromToken("91"))};

  for (int k = 0; k < 4; ++k) {
    EXPECT_THAT(kCell0[k].GetUVCoordOfEdge(k + 0), Eq(0));
    EXPECT_THAT(kCell0[k].GetUVCoordOfEdge(k + 1), Eq(0));
    EXPECT_THAT(kCell4[k].GetUVCoordOfEdge(k + 0), Eq(0));
    EXPECT_THAT(kCell4[k].GetUVCoordOfEdge(k + 1), Eq(0));
  }
}

TEST(S2Cell, GetSizeIJAgreesWithCellId) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "GET_SIZE_IJ_AGREES_WITH_CELL_ID",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < 100; ++i) {
    S2CellId id = s2random::CellId(bitgen);
    S2Cell cell(id);
    EXPECT_EQ(cell.GetSizeIJ(), id.GetSizeIJ(id.level()));
  }
}

TEST(S2Cell, GetIJCoordOfEdge) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "GET_IJ_COORD_OF_EDGE",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < 100; ++i) {
    S2CellId id = s2random::CellId(bitgen);
    S2Cell cell(id);

    // Look up the canonical IJ coordinates of the cell boundary.
    int ij[2];
    int orientation;
    id.ToFaceIJOrientation(ij, ij + 1, &orientation);

    int ij_size = cell.GetSizeIJ();
    R2Rect ij_bounds;
    for (int k = 0; k < 2; ++k) {
      int ij_lo = ij[k] & -ij_size;
      ij_bounds[k][0] = ij_lo;
      ij_bounds[k][1] = ij_lo + ij_size;
    }

    // Check that each boundary coordinate is correct.
    for (int k = 0; k < 4; ++k) {
      EXPECT_THAT(cell.GetIJCoordOfEdge(k),
                  Eq(ij_bounds.GetVertex(k)[(k + 1) % 2]));
    }
  }
}

