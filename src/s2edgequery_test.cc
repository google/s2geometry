// Copyright 2013 Google Inc. All Rights Reserved.
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

#include "s2edgequery.h"

#include <utility>  // pair<>
#include <vector>
#include "base/stringprintf.h"
#include "gtest/gtest.h"
#include "s1angle.h"
#include "s2cap.h"
#include "s2cell.h"
#include "s2cellid.h"
#include "s2edgeutil.h"
#include "s2shapeutil.h"
#include "s2testing.h"
#include "util/gtl/algorithm.h"

using std::pair;
using std::vector;

namespace {

typedef pair<S2Point, S2Point> TestEdge;

S2Point PerturbAtDistance(S1Angle distance,
                          S2Point const& a0, S2Point const& b0) {
  S2Point x = S2EdgeUtil::InterpolateAtDistance(distance, a0, b0);
  if (S2Testing::rnd.OneIn(2)) {
    for (int i = 0; i < 3; ++i) {
      x[i] = nextafter(x[i], S2Testing::rnd.OneIn(2) ? 1 : -1);
    }
    x = x.Normalize();
  }
  return x;
}

// Generate sub-edges of some given edge (a0,b0).  The length of the sub-edges
// is distributed exponentially over a large range, and the endpoints may be
// slightly perturbed to one side of (a0,b0) or the other.
void GetPerturbedSubEdges(S2Point a0, S2Point b0, int count,
                          vector<TestEdge>* edges) {
  edges->clear();
  a0 = a0.Normalize();
  b0 = b0.Normalize();
  S1Angle length0(a0, b0);
  for (int i = 0; i < count; ++i) {
    S1Angle length = length0 * pow(1e-15, S2Testing::rnd.RandDouble());
    S1Angle offset = (length0 - length) * S2Testing::rnd.RandDouble();
    edges->push_back(
        std::make_pair(PerturbAtDistance(offset, a0, b0),
                       PerturbAtDistance(offset + length, a0, b0)));
  }
}

// Generate edges whose center is randomly chosen from the given S2Cap, and
// whose length is randomly chosen up to "max_length".
void GetCapEdges(S2Cap const& center_cap, S1Angle max_length,
                 int count, vector<TestEdge>* edges) {
  edges->clear();
  for (int i = 0; i < count; ++i) {
    S2Point center = S2Testing::SamplePoint(center_cap);
    S2Cap edge_cap(center, 0.5 * max_length);
    S2Point p1 = S2Testing::SamplePoint(edge_cap);
    // Compute p1 reflected through "center", and normalize for good measure.
    S2Point p2 = (2 * p1.DotProd(center) * center - p1).Normalize();
    edges->push_back(std::make_pair(p1, p2));
  }
}

void TestAllCrossings(vector<TestEdge> const& edges) {
  S2EdgeVectorShape* shape = new S2EdgeVectorShape;
  for (int i = 0; i < edges.size(); ++i) {
    shape->Add(edges[i].first, edges[i].second);
  }
  // Force more subdivision than usual to make the test more challenging.
  S2ShapeIndexOptions options;
  options.set_max_edges_per_cell(1);
  S2ShapeIndex index(options);
  index.Insert(shape);  // Takes ownership
  // To check that candidates are being filtered reasonably, we count the
  // total number of candidates that the total number of edge pairs that
  // either intersect or are very close to intersecting.
  int num_candidates = 0, num_nearby_pairs = 0;
  for (int i = 0; i < edges.size(); ++i) {
    SCOPED_TRACE(StringPrintf("Iteration %d", i));
    S2Point const& a = edges[i].first;
    S2Point const& b = edges[i].second;
    vector<int> candidates;
    S2EdgeQuery query(index);
    query.GetCandidates(a, b, shape, &candidates);

    // Verify that the second version of GetCandidates returns the same result.
    S2EdgeQuery::EdgeMap edge_map;
    query.GetCandidates(a, b, &edge_map);
    EXPECT_EQ(1, edge_map.size());
    EXPECT_EQ(shape, edge_map.begin()->first);
    EXPECT_EQ(candidates, edge_map.begin()->second);
    EXPECT_TRUE(!candidates.empty());

    // Now check the actual candidates.
    EXPECT_TRUE(util::gtl::is_sorted(candidates.begin(), candidates.end()));
    EXPECT_GE(candidates.front(), 0);
    EXPECT_LT(candidates.back(), shape->num_edges());
    num_candidates += candidates.size();
    string missing_candidates;
    for (int i = 0; i < shape->num_edges(); ++i) {
      S2Point const *c, *d;
      shape->GetEdge(i, &c, &d);
      if (*c == a || *c == b || *d == a || *d == b ||
          S2EdgeUtil::RobustCrossing(a, b, *c, *d) > 0) {
        ++num_nearby_pairs;
        if (!std::binary_search(candidates.begin(), candidates.end(), i)) {
          StringAppendF(&missing_candidates, " %d", i);
        }
      } else {
        double const kMaxDist = S2::kMaxDiag.GetValue(S2::kMaxCellLevel);
        if (S2EdgeUtil::GetDistance(a, *c, *d).radians() < kMaxDist ||
            S2EdgeUtil::GetDistance(b, *c, *d).radians() < kMaxDist ||
            S2EdgeUtil::GetDistance(*c, a, b).radians() < kMaxDist ||
            S2EdgeUtil::GetDistance(*d, a, b).radians() < kMaxDist) {
          ++num_nearby_pairs;
        }
      }
    }
    EXPECT_TRUE(missing_candidates.empty()) << missing_candidates;
  }
  // There is nothing magical about this particular ratio; this check exists
  // to catch changes that dramatically increase the number of candidates.
  EXPECT_LE(num_candidates, 3 * num_nearby_pairs);
}

// Test edges that lie in the plane of one of the S2 cube edges.  Such edges
// may lie on the boundary between two cube faces, or pass through a cube
// vertex, or follow a 45 diagonal across a cube face toward its center.
//
// This test is sufficient to demonstrate that padding the cell boundaries
// is necessary for correctness.  (It fails if S2ShapeIndex::kCellPadding is
// set to zero.)
TEST(GetCrossingCandidates, PerturbedCubeEdges) {
  S2Testing::Random* rnd = &S2Testing::rnd;
  vector<TestEdge> edges;
  for (int iter = 0; iter < 10; ++iter) {
    int face = rnd->Uniform(6);
    double scale = pow(1e-15, rnd->RandDouble());
    R2Point uv(2 * rnd->Uniform(2) - 1, 2 * rnd->Uniform(2) - 1);  // vertex
    S2Point a0 = S2::FaceUVtoXYZ(face, scale * uv);
    S2Point b0 = a0 - 2 * S2::GetNorm(face);
    // TODO(ericv): This test is currently slow because *every* crossing test
    // needs to invoke S2::ExpensiveCCW().
    GetPerturbedSubEdges(a0, b0, 30, &edges);
    TestAllCrossings(edges);
  }
}

// Test edges that lie in the plane of one of the S2 cube face axes.  These
// edges are special because one coordinate is zero, and they lie on the
// boundaries between the immediate child cells of the cube face.
TEST(GetCrossingCandidates, PerturbedCubeFaceAxes) {
  S2Testing::Random* rnd = &S2Testing::rnd;
  vector<TestEdge> edges;
  for (int iter = 0; iter < 5; ++iter) {
    int face = rnd->Uniform(6);
    double scale = pow(1e-15, rnd->RandDouble());
    S2Point axis = S2::GetUVWAxis(face, rnd->Uniform(2));
    S2Point a0 = scale * axis + S2::GetNorm(face);
    S2Point b0 = scale * axis - S2::GetNorm(face);
    GetPerturbedSubEdges(a0, b0, 30, &edges);
    TestAllCrossings(edges);
  }
}

TEST(GetCrossingCandidates, CapEdgesNearCubeVertex) {
  // Test a random collection of edges near the S2 cube vertex where the
  // Hilbert curve starts and ends.
  vector<TestEdge> edges;
  GetCapEdges(S2Cap(S2Point(-1, -1, 1).Normalize(), S1Angle::Radians(1e-3)),
              S1Angle::Radians(1e-4), 1000, &edges);
  TestAllCrossings(edges);
}

TEST(GetCrossingCandidates, DegenerateEdgeOnCellVertexIsItsOwnCandidate) {
  for (int i = 0; i < 100; ++i) {
    vector<TestEdge> edges;
    S2Cell cell(S2Testing::GetRandomCellId());
    edges.push_back(std::make_pair(cell.GetVertex(0), cell.GetVertex(0)));
    TestAllCrossings(edges);
  }
}

TEST(GetCrossingCandidates, CollinearEdgesOnCellBoundaries) {
  const int kNumEdgeIntervals = 8;  // 9*8/2 = 36 edges
  for (int level = 0; level <= S2CellId::kMaxLevel; ++level) {
    S2Cell cell(S2Testing::GetRandomCellId(level));
    int v1 = S2Testing::rnd.Uniform(4);
    int v2 = (v1 + 1) & 3;
    S2Point p1 = cell.GetVertexRaw(v1);
    S2Point p2 = cell.GetVertexRaw(v2);
    S2Point delta = (p2 - p1) / kNumEdgeIntervals;
    vector<TestEdge> edges;
    for (int i = 0; i <= kNumEdgeIntervals; ++i) {
      for (int j = 0; j < i; ++j) {
        edges.push_back(std::make_pair((p1 + i * delta).Normalize(),
                                       (p1 + j * delta).Normalize()));
      }
    }
    TestAllCrossings(edges);
  }
}

}  // namespace