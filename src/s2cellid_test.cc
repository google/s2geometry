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

#include "s2cellid.h"

#include <math.h>
#include <stdio.h>
#include <algorithm>
#include <ext/hash_map>
using __gnu_cxx::hash;
using __gnu_cxx::hash_map;
#include <iosfwd>
#include <iostream>
#include <vector>

#include <gflags/gflags.h>
#include <glog/logging.h>

#include "base/macros.h"
#include "gtest/gtest.h"
#include "s2.h"
#include "s2latlng.h"
#include "s2testing.h"


using std::min;
using std::vector;

DEFINE_int32(iters, 20000000,
             "Number of iterations for timing tests with optimized build");

DEFINE_int32(htm_level, 29, "Maximum HTM level to use");
DEFINE_int32(build_level, 5, "HTM build level to use");

static S2CellId GetCellId(double lat_degrees, double lng_degrees) {
  S2CellId id = S2CellId::FromLatLng(S2LatLng::FromDegrees(lat_degrees,
                                                           lng_degrees));
  LOG(INFO) << std::hex << id.id();
  return id;
}

TEST(S2CellId, DefaultConstructor) {
  S2CellId id;
  EXPECT_EQ(id.id(), 0);
  EXPECT_FALSE(id.is_valid());
}

TEST(S2CellId, S2CellIdHasher) {
  EXPECT_EQ(S2CellIdHasher()(GetCellId(0, 90)),
            S2CellIdHasher()(GetCellId(0, 90)));
}

TEST(S2CellId, FaceDefinitions) {
  EXPECT_EQ(GetCellId(0, 0).face(), 0);
  EXPECT_EQ(GetCellId(0, 90).face(), 1);
  EXPECT_EQ(GetCellId(90, 0).face(), 2);
  EXPECT_EQ(GetCellId(0, 180).face(), 3);
  EXPECT_EQ(GetCellId(0, -90).face(), 4);
  EXPECT_EQ(GetCellId(-90, 0).face(), 5);
}

TEST(S2CellId, FromFace) {
  for (int face = 0; face < 6; ++face) {
    EXPECT_EQ(S2CellId::FromFacePosLevel(face, 0, 0), S2CellId::FromFace(face));
  }
}

TEST(S2CellId, ParentChildRelationships) {
  S2CellId id = S2CellId::FromFacePosLevel(3, 0x12345678,
                                           S2CellId::kMaxLevel - 4);
  EXPECT_TRUE(id.is_valid());
  EXPECT_EQ(id.face(), 3);
  EXPECT_EQ(id.pos(), 0x12345700);
  EXPECT_EQ(id.level(), S2CellId::kMaxLevel - 4);
  EXPECT_FALSE(id.is_leaf());

  EXPECT_EQ(id.child_begin(id.level() + 2).pos(), 0x12345610);
  EXPECT_EQ(id.child_begin().pos(), 0x12345640);
  EXPECT_EQ(id.parent().pos(), 0x12345400);
  EXPECT_EQ(id.parent(id.level() - 2).pos(), 0x12345000);

  // Check ordering of children relative to parents.
  EXPECT_LT(id.child_begin(), id);
  EXPECT_GT(id.child_end(), id);
  EXPECT_EQ(id.child_begin().next().next().next().next(), id.child_end());
  EXPECT_EQ(id.child_begin(S2CellId::kMaxLevel), id.range_min());
  EXPECT_EQ(id.child_end(S2CellId::kMaxLevel), id.range_max().next());

  // Check that cells are represented by the position of their center
  // along the Hilbert curve.
  EXPECT_EQ(id.range_min().id() + id.range_max().id(), 2 * id.id());
}

TEST(S2CellId, CenterSiTi) {
  S2CellId id = S2CellId::FromFacePosLevel(3, 0x12345678,
                                           S2CellId::kMaxLevel);
  // Check that the (si, ti) coordinates of the center end in a
  // 1 followed by (30 - level) 0s.
  int si, ti;

  // Leaf level, 30.
  id.GetCenterSiTi(&si, &ti);
  CHECK_EQ(1 << 0, si & 1);
  CHECK_EQ(1 << 0, ti & 1);

  // Level 29.
  id.parent(S2CellId::kMaxLevel - 1).GetCenterSiTi(&si, &ti);
  CHECK_EQ(1 << 1, si & 3);
  CHECK_EQ(1 << 1, ti & 3);

  // Level 28.
  id.parent(S2CellId::kMaxLevel - 2).GetCenterSiTi(&si, &ti);
  CHECK_EQ(1 << 2, si & 7);
  CHECK_EQ(1 << 2, ti & 7);

  // Level 20.
  id.parent(S2CellId::kMaxLevel - 10).GetCenterSiTi(&si, &ti);
  CHECK_EQ(1 << 10, si & ((1 << 11) - 1));
  CHECK_EQ(1 << 10, ti & ((1 << 11) - 1));

  // Level 10.
  id.parent(S2CellId::kMaxLevel - 20).GetCenterSiTi(&si, &ti);
  CHECK_EQ(1 << 20, si & ((1 << 21) - 1));
  CHECK_EQ(1 << 20, ti & ((1 << 21) - 1));

  // Level 0.
  id.parent(0).GetCenterSiTi(&si, &ti);
  CHECK_EQ(1 << 30, si & ((1U << 31) - 1));
  CHECK_EQ(1 << 30, ti & ((1U << 31) - 1));
}

TEST(S2CellId, Wrapping) {
  // Check wrapping from beginning of Hilbert curve to end and vice versa.
  EXPECT_EQ(S2CellId::Begin(0).prev_wrap(), S2CellId::End(0).prev());

  EXPECT_EQ(S2CellId::Begin(S2CellId::kMaxLevel).prev_wrap(),
            S2CellId::FromFacePosLevel(
                5, ~static_cast<uint64>(0) >> S2CellId::kFaceBits,
                S2CellId::kMaxLevel));
  EXPECT_EQ(S2CellId::Begin(S2CellId::kMaxLevel).advance_wrap(-1),
            S2CellId::FromFacePosLevel(
                5, ~static_cast<uint64>(0) >> S2CellId::kFaceBits,
                S2CellId::kMaxLevel));

  EXPECT_EQ(S2CellId::End(4).prev().next_wrap(), S2CellId::Begin(4));
  EXPECT_EQ(S2CellId::End(4).advance(-1).advance_wrap(1), S2CellId::Begin(4));

  EXPECT_EQ(S2CellId::End(S2CellId::kMaxLevel).prev().next_wrap(),
            S2CellId::FromFacePosLevel(0, 0, S2CellId::kMaxLevel));
  EXPECT_EQ(S2CellId::End(S2CellId::kMaxLevel).advance(-1).advance_wrap(1),
            S2CellId::FromFacePosLevel(0, 0, S2CellId::kMaxLevel));
}

TEST(S2CellId, Advance) {
  S2CellId id = S2CellId::FromFacePosLevel(3, 0x12345678,
                                           S2CellId::kMaxLevel - 4);
  // Check basic properties of advance().
  EXPECT_EQ(S2CellId::Begin(0).advance(7), S2CellId::End(0));
  EXPECT_EQ(S2CellId::Begin(0).advance(12), S2CellId::End(0));
  EXPECT_EQ(S2CellId::End(0).advance(-7), S2CellId::Begin(0));
  EXPECT_EQ(S2CellId::End(0).advance(-12000000), S2CellId::Begin(0));
  int num_level_5_cells = 6 << (2 * 5);
  EXPECT_EQ(S2CellId::Begin(5).advance(500),
            S2CellId::End(5).advance(500 - num_level_5_cells));
  EXPECT_EQ(id.child_begin(S2CellId::kMaxLevel).advance(256),
            id.next().child_begin(S2CellId::kMaxLevel));
  EXPECT_EQ(S2CellId::FromFacePosLevel(1, 0, S2CellId::kMaxLevel)
            .advance(static_cast<int64>(4) << (2 * S2CellId::kMaxLevel)),
            S2CellId::FromFacePosLevel(5, 0, S2CellId::kMaxLevel));

  // Check basic properties of advance_wrap().
  EXPECT_EQ(S2CellId::Begin(0).advance_wrap(7), S2CellId::FromFace(1));
  EXPECT_EQ(S2CellId::Begin(0).advance_wrap(12), S2CellId::Begin(0));
  EXPECT_EQ(S2CellId::FromFace(5).advance_wrap(-7), S2CellId::FromFace(4));
  EXPECT_EQ(S2CellId::Begin(0).advance_wrap(-12000000), S2CellId::Begin(0));
  EXPECT_EQ(S2CellId::Begin(5).advance_wrap(6644),
            S2CellId::Begin(5).advance_wrap(-11788));
  EXPECT_EQ(id.child_begin(S2CellId::kMaxLevel).advance_wrap(256),
            id.next().child_begin(S2CellId::kMaxLevel));
  EXPECT_EQ(S2CellId::FromFacePosLevel(5, 0, S2CellId::kMaxLevel)
            .advance_wrap(static_cast<int64>(2) << (2 * S2CellId::kMaxLevel)),
            S2CellId::FromFacePosLevel(1, 0, S2CellId::kMaxLevel));
}

TEST(S2CellId, MaximumTile) {
  // This method is tested more thoroughly in s2cellunion_test.cc.
  for (int iter = 0; iter < 1000; ++iter) {
    S2CellId id = S2Testing::GetRandomCellId(10);

    // Check that "limit" is returned for tiles at or beyond "limit".
    EXPECT_EQ(id, id.maximum_tile(id));
    EXPECT_EQ(id, id.child(0).maximum_tile(id));
    EXPECT_EQ(id, id.child(1).maximum_tile(id));
    EXPECT_EQ(id, id.next().maximum_tile(id));
    EXPECT_EQ(id.child(0), id.maximum_tile(id.child(0)));

    // Check that the tile size is increased when possible.
    EXPECT_EQ(id, id.child(0).maximum_tile(id.next()));
    EXPECT_EQ(id, id.child(0).maximum_tile(id.next().child(0)));
    EXPECT_EQ(id, id.child(0).maximum_tile(id.next().child(1).child(0)));
    EXPECT_EQ(id, id.child(0).child(0).maximum_tile(id.next()));
    EXPECT_EQ(id, id.child(0).child(0).child(0).maximum_tile(id.next()));

    // Check that the tile size is decreased when necessary.
    EXPECT_EQ(id.child(0), id.maximum_tile(id.child(0).next()));
    EXPECT_EQ(id.child(0), id.maximum_tile(id.child(0).next().child(0)));
    EXPECT_EQ(id.child(0), id.maximum_tile(id.child(0).next().child(1)));
    EXPECT_EQ(id.child(0).child(0),
              id.maximum_tile(id.child(0).child(0).next()));
    EXPECT_EQ(id.child(0).child(0).child(0),
              id.maximum_tile(id.child(0).child(0).child(0).next()));

    // Check that the tile size is otherwise unchanged.
    EXPECT_EQ(id, id.maximum_tile(id.next()));
    EXPECT_EQ(id, id.maximum_tile(id.next().child(0)));
    EXPECT_EQ(id, id.maximum_tile(id.next().child(1).child(0)));
  }
}

TEST(S2CellId, GetCommonAncestorLevel) {
  // Two identical cell ids.
  EXPECT_EQ(0, S2CellId::FromFace(0).
            GetCommonAncestorLevel(S2CellId::FromFace(0)));
  EXPECT_EQ(30, S2CellId::FromFace(0).child_begin(30).
            GetCommonAncestorLevel(S2CellId::FromFace(0).child_begin(30)));

  // One cell id is a descendant of the other.
  EXPECT_EQ(0, S2CellId::FromFace(0).child_begin(30).
            GetCommonAncestorLevel(S2CellId::FromFace(0)));
  EXPECT_EQ(0, S2CellId::FromFace(5).
            GetCommonAncestorLevel(S2CellId::FromFace(5).child_end(30).prev()));

  // Two cells that have no common ancestor.
  EXPECT_EQ(-1, S2CellId::FromFace(0).
            GetCommonAncestorLevel(S2CellId::FromFace(5)));
  EXPECT_EQ(-1, S2CellId::FromFace(2).child_begin(30).
            GetCommonAncestorLevel(S2CellId::FromFace(3).child_end(20)));

  // Two cells that have a common ancestor distinct from both of them.
  EXPECT_EQ(8, S2CellId::FromFace(5).child_begin(9).next().child_begin(15).
            GetCommonAncestorLevel(
                S2CellId::FromFace(5).child_begin(9).child_begin(20)));
  EXPECT_EQ(1, S2CellId::FromFace(0).child_begin(2).child_begin(30).
            GetCommonAncestorLevel(
                S2CellId::FromFace(0).child_begin(2).next().child_begin(5)));
}

TEST(S2CellId, Inverses) {
  // Check the conversion of random leaf cells to S2LatLngs and back.
  for (int i = 0; i < 200000; ++i) {
    S2CellId id = S2Testing::GetRandomCellId(S2CellId::kMaxLevel);
    EXPECT_TRUE(id.is_leaf());
    EXPECT_EQ(id.level(), S2CellId::kMaxLevel);
    S2LatLng center = id.ToLatLng();
    EXPECT_EQ(S2CellId::FromLatLng(center).id(), id.id());
  }
}

TEST(S2CellId, Tokens) {
  // Test random cell ids at all levels.
  for (int i = 0; i < 10000; ++i) {
    S2CellId id = S2Testing::GetRandomCellId();
    string token = id.ToToken();
    EXPECT_LE(token.size(), 16);
    EXPECT_EQ(S2CellId::FromToken(token), id);
    EXPECT_EQ(S2CellId::FromToken(token.data(), token.size()), id);
  }
  // Check that invalid cell ids can be encoded.
  string token = S2CellId::None().ToToken();
  EXPECT_EQ(S2CellId::FromToken(token), S2CellId::None());
  EXPECT_EQ(S2CellId::FromToken(token.data(), token.size()), S2CellId::None());

  // Check that supplying tokens with non-alphanumeric characters
  // returns S2CellId::None().
  EXPECT_EQ(S2CellId::FromToken("876b e99"), S2CellId::None());
  EXPECT_EQ(S2CellId::FromToken("876bee99\n"), S2CellId::None());
  EXPECT_EQ(S2CellId::FromToken("876[ee99"), S2CellId::None());
  EXPECT_EQ(S2CellId::FromToken(" 876bee99"), S2CellId::None());
}


static const int kMaxExpandLevel = 3;

static void ExpandCell(S2CellId parent, vector<S2CellId>* cells,
                       hash_map<S2CellId, S2CellId>* parent_map) {
  cells->push_back(parent);
  if (parent.level() == kMaxExpandLevel) return;
  int i, j, orientation;
  int face = parent.ToFaceIJOrientation(&i, &j, &orientation);
  EXPECT_EQ(face, parent.face());

  S2CellId child = parent.child_begin();
  for (int pos = 0; child != parent.child_end(); child = child.next(), ++pos) {
    (*parent_map)[child] = parent;
    // Do some basic checks on the children.
    EXPECT_EQ(parent.child(pos), child);
    EXPECT_EQ(pos, child.child_position());
    // Test child_position(level) on all the child's ancestors.
    for (S2CellId ancestor = child; ancestor.level() >= 1;
         ancestor = (*parent_map)[ancestor]) {
      EXPECT_EQ(ancestor.child_position(),
                child.child_position(ancestor.level()));
    }
    EXPECT_EQ(pos, child.child_position(child.level()));
    EXPECT_EQ(child.level(), parent.level() + 1);
    EXPECT_FALSE(child.is_leaf());
    int child_orientation;
    EXPECT_EQ(child.ToFaceIJOrientation(&i, &j, &child_orientation), face);
    EXPECT_EQ(child_orientation, orientation ^ S2::kPosToOrientation[pos]);
    ExpandCell(child, cells, parent_map);
  }
}

TEST(S2CellId, Containment) {
  // Test contains() and intersects().
  hash_map<S2CellId, S2CellId> parent_map;
  vector<S2CellId> cells;
  for (int face = 0; face < 6; ++face) {
    ExpandCell(S2CellId::FromFace(face), &cells, &parent_map);
  }
  for (int i = 0; i < cells.size(); ++i) {
    for (int j = 0; j < cells.size(); ++j) {
      bool contained = true;
      for (S2CellId id = cells[j]; id != cells[i]; id = parent_map[id]) {
        if (parent_map.find(id) == parent_map.end()) {
          contained = false;
          break;
        }
      }
      EXPECT_EQ(cells[i].contains(cells[j]), contained);
      EXPECT_EQ(cells[j] >= cells[i].range_min() &&
                cells[j] <= cells[i].range_max(), contained);
      EXPECT_EQ(cells[i].intersects(cells[j]),
                cells[i].contains(cells[j]) || cells[j].contains(cells[i]));
    }
  }
}

static int const kMaxWalkLevel = 8;

TEST(S2CellId, Continuity) {
  // Make sure that sequentially increasing cell ids form a continuous
  // path over the surface of the sphere, i.e. there are no
  // discontinuous jumps from one region to another.

  double max_dist = S2::kMaxEdge.GetValue(kMaxWalkLevel);
  S2CellId end = S2CellId::End(kMaxWalkLevel);
  S2CellId id = S2CellId::Begin(kMaxWalkLevel);
  for (; id != end; id = id.next()) {
    EXPECT_LE(id.ToPointRaw().Angle(id.next_wrap().ToPointRaw()), max_dist);
    EXPECT_EQ(id.advance_wrap(1), id.next_wrap());
    EXPECT_EQ(id.next_wrap().advance_wrap(-1), id);

    // Check that the ToPointRaw() returns the center of each cell
    // in (s,t) coordinates.
    double u, v;
    S2::XYZtoFaceUV(id.ToPointRaw(), &u, &v);
    static double const kCellSize = 1.0 / (1 << kMaxWalkLevel);
    EXPECT_NEAR(remainder(S2::UVtoST(u), 0.5 * kCellSize), 0.0, 1e-15);
    EXPECT_NEAR(remainder(S2::UVtoST(v), 0.5 * kCellSize), 0.0, 1e-15);
  }
}

TEST(S2CellId, Coverage) {
  // Make sure that random points on the sphere can be represented to the
  // expected level of accuracy, which in the worst case is sqrt(2/3) times
  // the maximum arc length between the points on the sphere associated with
  // adjacent values of "i" or "j".  (It is sqrt(2/3) rather than 1/2 because
  // the cells at the corners of each face are stretched -- they have 60 and
  // 120 degree angles.)

  double max_dist = 0.5 * S2::kMaxDiag.GetValue(S2CellId::kMaxLevel);
  for (int i = 0; i < 1000000; ++i) {
    S2Point p = S2Testing::RandomPoint();
    S2Point q = S2CellId::FromPoint(p).ToPointRaw();
    EXPECT_LE(p.Angle(q), max_dist);
  }
}

static void TestAllNeighbors(S2CellId id, int level) {
  DCHECK_GE(level, id.level());
  DCHECK_LT(level, S2CellId::kMaxLevel);

  // We compute AppendAllNeighbors, and then add in all the children of "id"
  // at the given level.  We then compare this against the result of finding
  // all the vertex neighbors of all the vertices of children of "id" at the
  // given level.  These should give the same result.
  vector<S2CellId> all, expected;
  id.AppendAllNeighbors(level, &all);
  S2CellId end = id.child_end(level + 1);
  for (S2CellId c = id.child_begin(level + 1); c != end; c = c.next()) {
    all.push_back(c.parent());
    c.AppendVertexNeighbors(level, &expected);
  }
  // Sort the results and eliminate duplicates.
  std::sort(all.begin(), all.end());
  std::sort(expected.begin(), expected.end());
  all.erase(std::unique(all.begin(), all.end()), all.end());
  expected.erase(std::unique(expected.begin(), expected.end()), expected.end());
  EXPECT_EQ(expected, all);
}

TEST(S2CellId, Neighbors) {
  // Check the edge neighbors of face 1.
  static int out_faces[] = { 5, 3, 2, 0 };
  S2CellId face_nbrs[4];
  S2CellId::FromFace(1).GetEdgeNeighbors(face_nbrs);
  for (int i = 0; i < 4; ++i) {
    EXPECT_TRUE(face_nbrs[i].is_face());
    EXPECT_EQ(face_nbrs[i].face(), out_faces[i]);
  }

  // Check the edge neighbors of the corner cells at all levels.  This case is
  // trickier because it requires projecting onto adjacent faces.
  static int const kMaxIJ = S2CellId::kMaxSize - 1;
  for (int level = 1; level <= S2CellId::kMaxLevel; ++level) {
    S2CellId id = S2CellId::FromFaceIJ(1, 0, 0).parent(level);
    S2CellId nbrs[4];
    id.GetEdgeNeighbors(nbrs);
    // These neighbors were determined manually using the face and axis
    // relationships defined in s2.cc.
    int size_ij = S2CellId::GetSizeIJ(level);
    EXPECT_EQ(S2CellId::FromFaceIJ(5, kMaxIJ, kMaxIJ).parent(level), nbrs[0]);
    EXPECT_EQ(S2CellId::FromFaceIJ(1, size_ij, 0).parent(level), nbrs[1]);
    EXPECT_EQ(S2CellId::FromFaceIJ(1, 0, size_ij).parent(level), nbrs[2]);
    EXPECT_EQ(S2CellId::FromFaceIJ(0, kMaxIJ, 0).parent(level), nbrs[3]);
  }

  // Check the vertex neighbors of the center of face 2 at level 5.
  vector<S2CellId> nbrs;
  S2CellId::FromPoint(S2Point(0, 0, 1)).AppendVertexNeighbors(5, &nbrs);
  std::sort(nbrs.begin(), nbrs.end());
  for (int i = 0; i < 4; ++i) {
    EXPECT_EQ(nbrs[i], S2CellId::FromFaceIJ(
                 2, (1 << 29) - (i < 2), (1 << 29) - (i == 0 || i == 3))
             .parent(5));
  }
  nbrs.clear();

  // Check the vertex neighbors of the corner of faces 0, 4, and 5.
  S2CellId id = S2CellId::FromFacePosLevel(0, 0, S2CellId::kMaxLevel);
  id.AppendVertexNeighbors(0, &nbrs);
  std::sort(nbrs.begin(), nbrs.end());
  EXPECT_EQ(nbrs.size(), 3);
  EXPECT_EQ(nbrs[0], S2CellId::FromFace(0));
  EXPECT_EQ(nbrs[1], S2CellId::FromFace(4));
  EXPECT_EQ(nbrs[2], S2CellId::FromFace(5));

  // Check that AppendAllNeighbors produces results that are consistent
  // with AppendVertexNeighbors for a bunch of random cells.
  for (int i = 0; i < 1000; ++i) {
    S2CellId id = S2Testing::GetRandomCellId();
    if (id.is_leaf()) id = id.parent();

    // TestAllNeighbors computes approximately 2**(2*(diff+1)) cell ids,
    // so it's not reasonable to use large values of "diff".
    int max_diff = min(6, S2CellId::kMaxLevel - id.level() - 1);
    int level = id.level() + S2Testing::rnd.Uniform(max_diff);
    TestAllNeighbors(id, level);
  }
}

TEST(S2CellId, OutputOperator) {
  S2CellId cell(0xbb04000000000000ULL);
  std::ostringstream s;
  s << cell;
  EXPECT_EQ("5/31200", s.str());
}

TEST(S2CellId, ToPointBenchmark) {
  // This "test" is really a benchmark, so skip it unless we're optimized.
  if (google::DEBUG_MODE) return;

  // Test speed of conversions from points to leaf cells.
  double control_start = S2Testing::GetCpuTime();
  S2CellId begin = S2CellId::Begin(S2CellId::kMaxLevel);
  S2CellId end = S2CellId::End(S2CellId::kMaxLevel);
  uint64 delta = (end.id() - begin.id()) / FLAGS_iters;
  delta &= ~static_cast<uint64>(1);  // Make sure all ids are leaf cells.

  S2CellId id = begin;
  double sum = 0;
  for (int i = FLAGS_iters; i > 0; --i) {
    sum += static_cast<double>(id.id());
    id = S2CellId(id.id() + delta);
  }
  double control_time = S2Testing::GetCpuTime() - control_start;
  printf("\tControl:    %8.3f usecs\n", 1e6 * control_time / FLAGS_iters);
  EXPECT_NE(sum, 0);  // Don't let the loop get optimized away.

  double test_start = S2Testing::GetCpuTime();
  sum = 0;
  id = begin;
  for (int i = FLAGS_iters; i > 0; --i) {
    sum += id.ToPointRaw()[0];
    id = S2CellId(id.id() + delta);
  }
  double test_time = S2Testing::GetCpuTime() - test_start - control_time;
  printf("\tToPointRaw: %8.3f usecs\n", 1e6 * test_time / FLAGS_iters);
  EXPECT_NE(sum, 0);  // Don't let the loop get optimized away.

  test_start = S2Testing::GetCpuTime();
  sum = 0;
  id = begin;
  for (int i = FLAGS_iters; i > 0; --i) {
    sum += id.ToPoint()[0];
    id = S2CellId(id.id() + delta);
  }
  test_time = S2Testing::GetCpuTime() - test_start - control_time;
  printf("\tToPoint:    %8.3f usecs\n", 1e6 * test_time / FLAGS_iters);
  EXPECT_NE(sum, 0);  // Don't let the loop get optimized away.
}

TEST(S2CellId, FromPointBenchmark) {
  // This "test" is really a benchmark, so skip it unless we're optimized.
  if (google::DEBUG_MODE) return;

  // The sample points follow a spiral curve that completes one revolution
  // around the z-axis every 1/dt samples.  The z-coordinate increases
  // from -4 to +4 over FLAGS_iters samples.

  S2Point start(1, 0, -4);
  double dz = (-2 * start.z()) / FLAGS_iters;
  double dt = 1.37482937133e-4;

  // Test speed of conversions from leaf cells to points.
  double control_start = S2Testing::GetCpuTime();
  uint64 isum = 0;
  S2Point p = start;
  for (int i = FLAGS_iters; i > 0; --i) {
    // Cheap rotation around the z-axis (spirals inward slightly
    // each revolution).
    p += S2Point(-dt * p.y(), dt * p.x(), dz);
    isum += MathUtil::FastIntRound(p[0] + p[1] + p[2]);
  }
  double control_time = S2Testing::GetCpuTime() - control_start;
  printf("\tControl:    %8.3f usecs\n", 1e6 * control_time / FLAGS_iters);
  EXPECT_NE(isum, 0);  // Don't let the loop get optimized away.

  double test_start = S2Testing::GetCpuTime();
  isum = 0;
  p = start;
  for (int i = FLAGS_iters; i > 0; --i) {
    p += S2Point(-dt * p.y(), dt * p.x(), dz);
    isum += S2CellId::FromPoint(p).id();
  }
  double test_time = S2Testing::GetCpuTime() - test_start - control_time;
  printf("\tFromPoint:  %8.3f usecs\n", 1e6 * test_time / FLAGS_iters);
  EXPECT_NE(isum, 0);  // Don't let the loop get optimized away.
}
