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

#include "s2/s2random.h"

#include <cmath>
#include <cstdint>

#include "s2/base/commandlineflags.h"
#include <gtest/gtest.h>
#include "absl/flags/flag.h"
#include "absl/log/log_streamer.h"
#include "absl/random/random.h"
#include "s2/s2cap.h"
#include "s2/s2cell_id.h"
#include "s2/s2latlng_rect.h"
#include "s2/s2point.h"
#include "s2/s2pointutil.h"
#include "s2/s2testing.h"
#include "s2/util/math/matrix3x3.h"

// The default value should be chosen so the test runs in a few seconds.
S2_DEFINE_int32(num_samples, 100'000,
                "Number of random samples to use in tests.");

namespace {

TEST(LogUniform, InRange) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "LOG_UNIFORM_IN_RANGE",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    // 1e-30 to 1e30 is in the range of values we typically use.
    const double log10_lo = absl::Uniform(bitgen, -30.0, 30.0);
    const double log10_hi = absl::Uniform(bitgen, log10_lo, 30.0);

    const double lo = std::pow(10.0, log10_lo);
    const double hi = std::pow(10.0, log10_hi);

    const double v = s2random::LogUniform(bitgen, lo, hi);
    EXPECT_LE(lo, v);
    EXPECT_LE(v, hi);
  }
}

TEST(SkewedInt, InRange) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "SKEWED_INT_IN_RANGE",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const int max_log = absl::Uniform(bitgen, 1, 32);
    const int v = s2random::SkewedInt(bitgen, max_log);
    EXPECT_GE(v, 0);
    EXPECT_LT(v, uint32_t{1} << max_log);
  }
}

TEST(RandomPoint, UnitLength) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "RANDOM_POINT_VALID",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const S2Point p = s2random::Point(bitgen);
    EXPECT_TRUE(S2::IsUnitLength(p));
  }
}

// These thresholds have been experimentally determined to make the tests
// pass.  If you know how to do the numerical analysis to prove a bound,
// please update them!
//
// Maximum absolute value for dot product of orthogonal vectors from
// `Frame()`/`FrameAt()`.
constexpr double kFrameOrthoEps = 1e-12;
// Maximum absolute value for squared norm in right-handed test for
// `Frame()`/`FrameAt()`.
constexpr double kFrameRightHandedEps2 = 1e-25;

TEST(FrameVectors, UnitLength) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRAME_VECTORS_UNIT_LENGTH",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    S2Point x, y, z;
    s2random::Frame(bitgen, x, y, z);
    EXPECT_TRUE(S2::IsUnitLength(x));
    EXPECT_TRUE(S2::IsUnitLength(y));
    EXPECT_TRUE(S2::IsUnitLength(z));
  }
}

TEST(FrameVectors, Orthogonal) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRAME_VECTORS_ORTHOGONAL",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    S2Point x, y, z;
    s2random::Frame(bitgen, x, y, z);
    EXPECT_LE(std::abs(x.DotProd(y)), kFrameOrthoEps);
    EXPECT_LE(std::abs(x.DotProd(z)), kFrameOrthoEps);
    EXPECT_LE(std::abs(y.DotProd(z)), kFrameOrthoEps);
  }
}

TEST(FrameVectors, RightHanded) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRAME_VECTORS_RIGHT_HANDED",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    S2Point x, y, z;
    s2random::Frame(bitgen, x, y, z);
    // We know these vectors are orthogonal, so stability is not an issue.
    // Plain `CrossProd` can be used.
    const S2Point x_cross_y = x.CrossProd(y);
    const S2Point y_cross_z = y.CrossProd(z);
    const S2Point z_cross_x = z.CrossProd(x);
    EXPECT_LE((z - x_cross_y).Norm2(), kFrameRightHandedEps2);
    EXPECT_LE((x - y_cross_z).Norm2(), kFrameRightHandedEps2);
    EXPECT_LE((y - z_cross_x).Norm2(), kFrameRightHandedEps2);
  }
}

// `Frame(BitGenRef, Matrix3x3_d)` and `FrameAt()` tests are based on the
// `Frame(BitGenRef, S2Point, S2Point, S2Point)` tests, so omit comments for
// brevity.

TEST(FrameMatrix, UnitLength) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRAME_MATRIX_UNIT_LENGTH",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const Matrix3x3_d frame = s2random::Frame(bitgen);
    EXPECT_TRUE(S2::IsUnitLength(frame.Col(0)));
    EXPECT_TRUE(S2::IsUnitLength(frame.Col(1)));
    EXPECT_TRUE(S2::IsUnitLength(frame.Col(2)));
  }
}

TEST(FrameMatrix, Orthogonal) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRAME_MATRIX_ORTHOGONAL",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const Matrix3x3_d frame = s2random::Frame(bitgen);
    EXPECT_LE(std::abs(frame.Col(0).DotProd(frame.Col(1))), kFrameOrthoEps);
    EXPECT_LE(std::abs(frame.Col(0).DotProd(frame.Col(2))), kFrameOrthoEps);
    EXPECT_LE(std::abs(frame.Col(1).DotProd(frame.Col(2))), kFrameOrthoEps);
  }
}

TEST(FrameMatrix, RightHanded) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRAME_MATRIX_RIGHT_HANDED",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const Matrix3x3_d frame = s2random::Frame(bitgen);
    const S2Point x_cross_y = frame.Col(0).CrossProd(frame.Col(1));
    const S2Point y_cross_z = frame.Col(1).CrossProd(frame.Col(2));
    const S2Point z_cross_x = frame.Col(2).CrossProd(frame.Col(0));
    EXPECT_LE((frame.Col(2) - x_cross_y).Norm2(), kFrameRightHandedEps2);
    EXPECT_LE((frame.Col(0) - y_cross_z).Norm2(), kFrameRightHandedEps2);
    EXPECT_LE((frame.Col(1) - z_cross_x).Norm2(), kFrameRightHandedEps2);
  }
}

TEST(FrameAtVectors, UnitLength) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRAME_AT_VECTORS_UNIT_LENGTH",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const S2Point z = s2random::Point(bitgen);
    S2Point x, y;
    s2random::FrameAt(bitgen, z, x, y);
    EXPECT_TRUE(S2::IsUnitLength(x));
    EXPECT_TRUE(S2::IsUnitLength(y));
  }
}

TEST(FrameAtVectors, Orthogonal) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRAME_AT_VECTORS_ORTHOGONAL",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const S2Point z = s2random::Point(bitgen);
    S2Point x, y;
    s2random::FrameAt(bitgen, z, x, y);
    EXPECT_LE(std::abs(x.DotProd(y)), kFrameOrthoEps);
    EXPECT_LE(std::abs(x.DotProd(z)), kFrameOrthoEps);
    EXPECT_LE(std::abs(y.DotProd(z)), kFrameOrthoEps);
  }
}

TEST(FrameAtVectors, RightHanded) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRAME_AT_VECTORS_RIGHT_HANDED",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const S2Point z = s2random::Point(bitgen);
    S2Point x, y;
    s2random::FrameAt(bitgen, z, x, y);
    const S2Point x_cross_y = x.CrossProd(y);
    const S2Point y_cross_z = y.CrossProd(z);
    const S2Point z_cross_x = z.CrossProd(x);
    EXPECT_LE((z - x_cross_y).Norm2(), kFrameRightHandedEps2);
    EXPECT_LE((x - y_cross_z).Norm2(), kFrameRightHandedEps2);
    EXPECT_LE((y - z_cross_x).Norm2(), kFrameRightHandedEps2);
  }
}

TEST(FrameAtMatrix, ZIsCopy) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRAME_AT_MATRIX_Z_IS_COPY",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const S2Point z = s2random::Point(bitgen);
    const Matrix3x3_d frame = s2random::FrameAt(bitgen, z);
    // Last column is a copy of the input.  The other columns are tested below.
    EXPECT_EQ(z, frame.Col(2));
  }
}

TEST(FrameAtMatrix, UnitLength) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRAME_AT_MATRIX_UNIT_LENGTH",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const S2Point z = s2random::Point(bitgen);
    const Matrix3x3_d frame = s2random::FrameAt(bitgen, z);
    EXPECT_TRUE(S2::IsUnitLength(frame.Col(0)));
    EXPECT_TRUE(S2::IsUnitLength(frame.Col(1)));
    EXPECT_TRUE(S2::IsUnitLength(frame.Col(2)));
  }
}

TEST(FrameAtMatrix, Orthogonal) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRAME_AT_MATRIX_ORTHOGONAL",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const S2Point z = s2random::Point(bitgen);
    const Matrix3x3_d frame = s2random::FrameAt(bitgen, z);
    EXPECT_LE(std::abs(frame.Col(0).DotProd(frame.Col(1))), kFrameOrthoEps);
    EXPECT_LE(std::abs(frame.Col(0).DotProd(frame.Col(2))), kFrameOrthoEps);
    EXPECT_LE(std::abs(frame.Col(1).DotProd(frame.Col(2))), kFrameOrthoEps);
  }
}

TEST(FrameAtMatrix, RightHanded) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "FRAME_AT_MATRIX_RIGHT_HANDED",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const S2Point z = s2random::Point(bitgen);
    const Matrix3x3_d frame = s2random::FrameAt(bitgen, z);
    const S2Point x_cross_y = frame.Col(0).CrossProd(frame.Col(1));
    const S2Point y_cross_z = frame.Col(1).CrossProd(frame.Col(2));
    const S2Point z_cross_x = frame.Col(2).CrossProd(frame.Col(0));
    EXPECT_LE((frame.Col(2) - x_cross_y).Norm2(), kFrameRightHandedEps2);
    EXPECT_LE((frame.Col(0) - y_cross_z).Norm2(), kFrameRightHandedEps2);
    EXPECT_LE((frame.Col(1) - z_cross_x).Norm2(), kFrameRightHandedEps2);
  }
}

TEST(RandomCap, AreaInRange) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "RANDOM_CAP_AREA_IN_RANGE",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const double min_area =
        absl::Uniform(absl::IntervalOpenClosed, bitgen, 1e-16, 4 * M_PI);
    const double max_area =
        absl::Uniform(absl::IntervalOpenClosed, bitgen, min_area, 4 * M_PI);
    const S2Cap cap = s2random::Cap(bitgen, min_area, max_area);
    EXPECT_GE(cap.GetArea(), min_area);
    EXPECT_LE(cap.GetArea(), max_area);
  }
}

TEST(SamplePoint, InsideCap) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "SAMPLE_POINT_INSIDE_CAP",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    // We must avoid empty caps, so set a small min_area.
    const S2Cap cap =
        s2random::Cap(bitgen, /*min_area=*/1e-16, /*max_area=*/4 * M_PI);
    const S2Point p = s2random::SamplePoint(bitgen, cap);
    EXPECT_TRUE(cap.Contains(p)) << "cap: " << cap << " point: " << p;
  }
}

TEST(SamplePoint, InsideRect) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "SAMPLE_POINT_INSIDE_RECT",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  // There is no s2random::Rect(), so just use a fixed rect.
  const S2LatLngRect rect(S2LatLng::FromDegrees(-10, 15),
                          S2LatLng::FromDegrees(30, 50));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const S2Point p = s2random::SamplePoint(bitgen, rect);
    EXPECT_TRUE(rect.Contains(p)) << "rect: " << rect << " point: " << p;
  }
}

TEST(CellId, AtSpecifiedLevel) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "CELL_ID_AT_SPECIFIED_LEVEL",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const int level = absl::Uniform(absl::IntervalClosedClosed, bitgen, 0,
                                    S2CellId::kMaxLevel);
    const S2CellId cell_id = s2random::CellId(bitgen, level);
    EXPECT_TRUE(cell_id.is_valid());
    EXPECT_EQ(cell_id.level(), level);
  }
}

TEST(CellId, Valid) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "CELL_ID_VALID",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int i = 0; i < absl::GetFlag(FLAGS_num_samples); ++i) {
    const S2CellId cell_id = s2random::CellId(bitgen);
    EXPECT_TRUE(cell_id.is_valid());
  }
}

}  // namespace
