// Copyright 2011 Google Inc. All Rights Reserved.
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


#include "s2pointcompression.h"

#include <stddef.h>
#include <string>
#include <vector>

#include <gflags/gflags.h>
#include <glog/logging.h>
#include <gtest/gtest.h>
#include "util/coding/coder.h"
#include "s1angle.h"
#include "s2.h"
#include "s2cellid.h"
#include "s2testing.h"
#include "s2textformat.h"
#include "util/gtl/fixedarray.h"

using std::vector;

DEFINE_int32(s2pointcompression_bm_level, 30,
             "Level to encode at for benchmarks.");
DEFINE_double(s2pointcompression_bm_radius_km, 1000.0,
              "Radius to use for loop for benchmarks.");

namespace {

S2Point SnapPointToLevel(S2Point const& point, int level) {
  return S2CellId::FromPoint(point).parent(level).ToPoint();
}

vector<S2Point> SnapPointsToLevel(vector<S2Point> const& points,
                                  int level) {
  vector<S2Point> snapped_points(points.size());
  for (int i = 0; i < points.size(); ++i) {
    snapped_points[i] = SnapPointToLevel(points[i], level);
  }
  return snapped_points;
}

// Make a regular loop around the corner of faces 0, 1, and 2 with the
// specified radius in meters (on the earth) and number of vertices.
vector<S2Point> MakeRegularPoints(int num_vertices,
                                  double radius_km,
                                  int level) {
  S2Point const center = S2Point(1.0, 1.0, 1.0).Normalize();
  S1Angle const radius_angle = S2Testing::KmToAngle(radius_km);

  vector<S2Point> const unsnapped_points =
      S2Testing::MakeRegularPoints(center, radius_angle, num_vertices);

  return SnapPointsToLevel(unsnapped_points, level);
}

void MakeXYZFaceSiTiPoints(S2Point const* points, int num_points,
                           S2XYZFaceSiTi* result) {
  for (int i = 0; i < num_points; ++i) {
    result[i].xyz = points[i];
    result[i].cell_level = S2::XYZtoFaceSiTi(points[i], &result[i].face,
                                             &result[i].si, &result[i].ti);
  }
}

bool PointsEqual(S2Point const* a, int num_a, S2Point const* b, int num_b) {
  if (num_a != num_b) return false;
  for (int i = 0; i < num_a; ++i) {
    if (a[i] != b[i]) {
      return false;
    }
  }
  return true;
}

class S2PointCompressionTest : public ::testing::Test {
 protected:
  virtual void SetUp() {
    loop_4_ = MakeRegularPoints(4, 0.1, S2::kMaxCellLevel);

    S2Point const center = S2Point(1.0, 1.0, 1.0).Normalize();
    S1Angle const radius = S2Testing::KmToAngle(0.1);
    loop_4_unsnapped_ = S2Testing::MakeRegularPoints(center, radius, 4);

    // Radius is 100m, so points are about 141 meters apart.
    // Snapping to level 14 will move them by < 47m.
    loop_4_level_14_ = MakeRegularPoints(4, 0.1, 14);

    loop_100_ = MakeRegularPoints(100, 0.1, S2::kMaxCellLevel);

    loop_100_unsnapped_ = S2Testing::MakeRegularPoints(center, radius, 100);

    loop_100_mixed_15_ = S2Testing::MakeRegularPoints(center, radius, 100);
    for (int i = 0; i < 15; ++i) {
      loop_100_mixed_15_[3 * i] = SnapPointToLevel(loop_100_mixed_15_[3 * i],
                                                   S2::kMaxCellLevel);
    }

    loop_100_mixed_25_ = S2Testing::MakeRegularPoints(center, radius, 100);
    for (int i = 0; i < 25; ++i) {
      loop_100_mixed_25_[4 * i] = SnapPointToLevel(loop_100_mixed_25_[4 * i],
                                                   S2::kMaxCellLevel);
    }

    // Circumference is 628m, so points are about 6 meters apart.
    // Snapping to level 22 will move them by < 2m.
    loop_100_level_22_ = MakeRegularPoints(100, 0.1, 22);

    vector<S2Point> multi_face_points(6);
    multi_face_points[0] = S2::FaceUVtoXYZ(0, -0.5, 0.5).Normalize();
    multi_face_points[1] = S2::FaceUVtoXYZ(1, -0.5, 0.5).Normalize();
    multi_face_points[2] = S2::FaceUVtoXYZ(1, 0.5, -0.5).Normalize();
    multi_face_points[3] = S2::FaceUVtoXYZ(2, -0.5, 0.5).Normalize();
    multi_face_points[4] = S2::FaceUVtoXYZ(2, 0.5, -0.5).Normalize();
    multi_face_points[5] = S2::FaceUVtoXYZ(2, 0.5, 0.5).Normalize();
    loop_multi_face_ = SnapPointsToLevel(multi_face_points, S2::kMaxCellLevel);

    vector<S2Point> line_points(100);
    for (int i = 0; i < line_points.size(); ++i) {
      double const s = 0.01 + 0.005 * i;
      double const t = 0.01 + 0.009 * i;
      double const u = S2::STtoUV(s);
      double const v = S2::STtoUV(t);
      line_points[i] = S2::FaceUVtoXYZ(0, u, v).Normalize();
    }
    line_ = SnapPointsToLevel(line_points, S2::kMaxCellLevel);
  }

  void Encode(S2Point const* points, int num_points, int level) {
    FixedArray<S2XYZFaceSiTi> pts(num_points);
    MakeXYZFaceSiTiPoints(points, num_points, pts.get());
    S2EncodePointsCompressed(pts.get(), num_points, level, &encoder_);
  }

  void Decode(int num_points, int level, S2Point* points) {
    decoder_.reset(encoder_.base(), encoder_.length());
    EXPECT_TRUE(S2DecodePointsCompressed(&decoder_, num_points, level, points));
  }

  void Roundtrip(const vector<S2Point>& loop, int level) {
    Encode(&loop[0], loop.size(), level);
    vector<S2Point> points(loop.size());
    Decode(loop.size(), level, &points[0]);

    EXPECT_TRUE(PointsEqual(&loop[0], loop.size(), &points[0], points.size()))
        << "Decoded points\n" << s2textformat::ToString(points)
        << "\ndo not match original points\n"
        << s2textformat::ToString(loop);
  }

  Encoder encoder_;
  Decoder decoder_;

  // Four vertex loop near the corner of faces 0, 1, and 2.
  vector<S2Point> loop_4_;
  // Four vertex loop near the corner of faces 0, 1, and 2;
  // unsnapped.
  vector<S2Point> loop_4_unsnapped_;
  // Four vertex loop near the corner of faces 0, 1, and 2;
  // snapped to level 14.
  vector<S2Point> loop_4_level_14_;
  // 100 vertex loop near the corner of faces 0, 1, and 2.
  vector<S2Point> loop_100_;
  // 100 vertex loop near the corner of faces 0, 1, and 2;
  // unsnapped.
  vector<S2Point> loop_100_unsnapped_;
  // 100 vertex loop near the corner of faces 0, 1, and 2;
  // 15 points snapped to kMakCellLevel, the others not snapped.
  vector<S2Point> loop_100_mixed_15_;
  // 100 vertex loop near the corner of faces 0, 1, and 2;
  // 25 points snapped to kMakCellLevel, the others not snapped.
  vector<S2Point> loop_100_mixed_25_;
  // 100 vertex loop near the corner of faces 0, 1, and 2;
  // snapped to level 22.
  vector<S2Point> loop_100_level_22_;
  // A loop with two vertices on each of three faces.
  vector<S2Point> loop_multi_face_;
  // A straight line of 100 vertices on face 0 that should compress well.
  vector<S2Point> line_;
};

TEST_F(S2PointCompressionTest, RoundtripsEmpty) {
  // Just check this doesn't crash.
  Encode(NULL, 0, S2::kMaxCellLevel);
  Decode(0, S2::kMaxCellLevel, NULL);
}

TEST_F(S2PointCompressionTest, RoundtripsFourVertexLoop) {
  Roundtrip(loop_4_, S2::kMaxCellLevel);
}

TEST_F(S2PointCompressionTest, RoundtripsFourVertexLoopUnsnapped) {
  Roundtrip(loop_4_unsnapped_, S2::kMaxCellLevel);
}

TEST_F(S2PointCompressionTest, FourVertexLoopSize) {
  Encode(&loop_4_[0], loop_4_.size(), S2::kMaxCellLevel);
  // It would take 32 bytes uncompressed.
  EXPECT_EQ(39, encoder_.length());
}

TEST_F(S2PointCompressionTest, RoundtripsFourVertexLevel14Loop) {
  int const level = 14;
  Roundtrip(loop_4_level_14_, level);
}

TEST_F(S2PointCompressionTest, FourVertexLevel14LoopSize) {
  int const level = 14;
  Encode(&loop_4_level_14_[0], loop_4_level_14_.size(), level);
  // It would take 4 bytes per vertex without compression.
  EXPECT_EQ(23, encoder_.length());
}

TEST_F(S2PointCompressionTest, Roundtrips100VertexLoop) {
  Roundtrip(loop_100_, S2::kMaxCellLevel);
}

TEST_F(S2PointCompressionTest, Roundtrips100VertexLoopUnsnapped) {
  Roundtrip(loop_100_unsnapped_, S2::kMaxCellLevel);
}

TEST_F(S2PointCompressionTest, Roundtrips100VertexLoopMixed15) {
  Roundtrip(loop_100_mixed_15_, S2::kMaxCellLevel);
  EXPECT_EQ(2381, encoder_.length());
}

TEST_F(S2PointCompressionTest, Roundtrips100VertexLoopMixed25) {
  Roundtrip(loop_100_mixed_25_, S2::kMaxCellLevel);
  EXPECT_EQ(2131, encoder_.length());
}

TEST_F(S2PointCompressionTest, OneHundredVertexLoopSize) {
  Encode(&loop_100_[0], loop_100_.size(), S2::kMaxCellLevel);
  EXPECT_EQ(257, encoder_.length());
}

TEST_F(S2PointCompressionTest, OneHundredVertexLoopUnsnappedSize) {
  Encode(&loop_100_unsnapped_[0], loop_100_unsnapped_.size(),
         S2::kMaxCellLevel);
  EXPECT_EQ(2756, encoder_.length());
}

TEST_F(S2PointCompressionTest, Roundtrips100VertexLevel22Loop) {
  int const level = 22;
  Roundtrip(loop_100_level_22_, level);
}

TEST_F(S2PointCompressionTest, OneHundredVertexLoopLevel22Size) {
  Encode(&loop_100_level_22_[0], loop_100_level_22_.size(), 22);
  EXPECT_EQ(148, encoder_.length());
}

TEST_F(S2PointCompressionTest, MultiFaceLoop) {
  Roundtrip(loop_multi_face_, S2::kMaxCellLevel);
}

TEST_F(S2PointCompressionTest, StraightLineCompressesWell) {
  Roundtrip(line_, S2::kMaxCellLevel);
  // About 1 byte / vertex.
  EXPECT_EQ(line_.size() + 17, encoder_.length());
}


}  // namespace
