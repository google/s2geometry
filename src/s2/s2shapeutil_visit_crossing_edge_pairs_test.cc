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

#include "s2/s2shapeutil_visit_crossing_edge_pairs.h"

#include <algorithm>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include "absl/log/absl_log.h"
#include "absl/strings/string_view.h"

#include "s2/mutable_s2shape_index.h"
#include "s2/s2crossing_edge_query.h"
#include "s2/s2debug.h"
#include "s2/s2edge_crossings.h"
#include "s2/s2edge_vector_shape.h"
#include "s2/s2error.h"
#include "s2/s2latlng.h"
#include "s2/s2loop.h"
#include "s2/s2point.h"
#include "s2/s2polygon.h"
#include "s2/s2shape.h"
#include "s2/s2shape_index.h"
#include "s2/s2shapeutil_edge_iterator.h"
#include "s2/s2shapeutil_shape_edge.h"
#include "s2/s2shapeutil_shape_edge_id.h"
#include "s2/s2text_format.h"

using absl::string_view;
using std::make_unique;
using std::string;
using std::unique_ptr;
using std::vector;

namespace s2shapeutil {

// A set of edge pairs within an S2ShapeIndex.
using EdgePairVector = vector<std::pair<ShapeEdgeId, ShapeEdgeId>>;

// Get crossings in one index.
EdgePairVector GetCrossings(const S2ShapeIndex& index, CrossingType type) {
  EdgePairVector edge_pairs;
  VisitCrossingEdgePairs(
      index, type, [&edge_pairs](const ShapeEdge& a, const ShapeEdge& b, bool) {
        edge_pairs.push_back(std::make_pair(a.id(), b.id()));
        return true;  // Continue visiting.
      });
  if (edge_pairs.size() > 1) {
    std::sort(edge_pairs.begin(), edge_pairs.end());
    edge_pairs.erase(std::unique(edge_pairs.begin(), edge_pairs.end()),
                     edge_pairs.end());
  }
  return edge_pairs;
}

// Get crossings between two indexes.
EdgePairVector GetCrossings(const S2ShapeIndex& indexA,
                            const S2ShapeIndex& indexB, CrossingType type) {
  EdgePairVector edge_pairs;
  VisitCrossingEdgePairs(
      indexA, indexB, type,
      [&edge_pairs](const ShapeEdge& a, const ShapeEdge& b, bool) {
        edge_pairs.push_back(std::make_pair(a.id(), b.id()));
        return true;  // Continue visiting.
      });
  if (edge_pairs.size() > 1) {
    std::sort(edge_pairs.begin(), edge_pairs.end());
    edge_pairs.erase(std::unique(edge_pairs.begin(), edge_pairs.end()),
                     edge_pairs.end());
  }
  return edge_pairs;
}

// Brute force crossings in one index.
EdgePairVector GetCrossingEdgePairsBruteForce(const S2ShapeIndex& index,
                                              CrossingType type) {
  EdgePairVector result;
  int min_sign = (type == CrossingType::ALL) ? 0 : 1;
  for (EdgeIterator a_iter(&index); !a_iter.Done(); a_iter.Next()) {
    auto a = a_iter.edge();
    EdgeIterator b_iter = a_iter;
    for (b_iter.Next(); !b_iter.Done(); b_iter.Next()) {
      auto b = b_iter.edge();
      if (S2::CrossingSign(a.v0, a.v1, b.v0, b.v1) >= min_sign) {
        result.push_back(
            std::make_pair(a_iter.shape_edge_id(), b_iter.shape_edge_id()));
      }
    }
  }
  return result;
}

// Brute force crossings between two indexes.
EdgePairVector GetCrossingEdgePairsBruteForce(const S2ShapeIndex& indexA,
                                              const S2ShapeIndex& indexB,
                                              CrossingType type) {
  EdgePairVector result;
  int min_sign = (type == CrossingType::ALL) ? 0 : 1;
  for (EdgeIterator a_iter(&indexA); !a_iter.Done(); a_iter.Next()) {
    auto a = a_iter.edge();
    for (EdgeIterator b_iter(&indexB); !b_iter.Done(); b_iter.Next()) {
      auto b = b_iter.edge();
      if (S2::CrossingSign(a.v0, a.v1, b.v0, b.v1) >= min_sign) {
        result.push_back(
            std::make_pair(a_iter.shape_edge_id(), b_iter.shape_edge_id()));
      }
    }
  }
  return result;
}

std::ostream& operator<<(std::ostream& os,
                         const std::pair<ShapeEdgeId, ShapeEdgeId>& pair) {
  return os << "(" << pair.first << "," << pair.second << ")";
}


// Compares an 'expected' and an 'actual' set of edge crossings.
void TestCrossingEdgePairs(const EdgePairVector& expected,
                           const EdgePairVector& actual) {
  if (actual != expected) {
    ADD_FAILURE() << "Unexpected edge pairs; see details below."
                  << "\nExpected number of edge pairs: " << expected.size()
                  << "\nActual number of edge pairs: " << actual.size();
    for (const auto& edge_pair : expected) {
      if (std::count(actual.begin(), actual.end(), edge_pair) != 1) {
        std::cout << "Missing value: " << edge_pair << std::endl;
      }
    }
    for (const auto& edge_pair : actual) {
      if (std::count(expected.begin(), expected.end(), edge_pair) != 1) {
        std::cout << "Extra value: " << edge_pair << std::endl;
      }
    }
  }
}

// Compare brute force to optimized algorithms for finding crossings within one
// index. Checks that both algorithms find the expected number of crossings.
void TestGetCrossingEdgePairs(const S2ShapeIndex& index,
                              CrossingType type, int expected_crossing_count) {
  EdgePairVector expected = GetCrossingEdgePairsBruteForce(index, type);
  EdgePairVector actual = GetCrossings(index, type);
  EXPECT_EQ(expected_crossing_count, expected.size());
  EXPECT_EQ(expected_crossing_count, actual.size());
  TestCrossingEdgePairs(expected, actual);
}

// Compare brute force to optimized algorithms for finding crossings between two
// indexes. Checks that both algorithms find the expected number of crossings.
void TestGetCrossingEdgePairs(const S2ShapeIndex& indexA,
                              const S2ShapeIndex& indexB, CrossingType type,
                             int expected_crossing_count) {
  EdgePairVector expected =
      GetCrossingEdgePairsBruteForce(indexA, indexB, type);
  EdgePairVector actual = GetCrossings(indexA, indexB, type);
  EXPECT_EQ(expected_crossing_count, expected.size());
  EXPECT_EQ(expected_crossing_count, actual.size());
  TestCrossingEdgePairs(expected, actual);
}

TEST(GetCrossingEdgePairs, NoIntersectionsOneIndex) {
  MutableS2ShapeIndex index;
  TestGetCrossingEdgePairs(index, CrossingType::ALL, 0);
  TestGetCrossingEdgePairs(index, CrossingType::INTERIOR, 0);
}

TEST(GetCrossingEdgePairs, NoIntersectionsTwoIndexes) {
  MutableS2ShapeIndex indexA;
  MutableS2ShapeIndex indexB;
  TestGetCrossingEdgePairs(indexA, indexB, CrossingType::ALL, 0);
  TestGetCrossingEdgePairs(indexA, indexB, CrossingType::INTERIOR, 0);
}

TEST(GetCrossingEdgePairs, EdgeGridOneIndex) {
  constexpr int kGridSize = 10;
  double epsilon = 1e-10;
  // There are 11 horizontal and 11 vertical lines. The expected number of
  // interior crossings is 9x9, plus 9 "touching" intersections along each of
  // the left, right, and bottom edges. "epsilon" is used to make the interior
  // lines slightly longer so the "touches" actually cross, otherwise 3 of the
  // 27 touches are not considered intersecting.
  // However, the vertical lines do not reach the top line as it curves on the
  // surface of the sphere: despite "epsilon" those 9 are not even very close
  // to intersecting. Thus 9 * 12 = 108 interior and four more at the corners
  // when CrossingType::ALL is used.
  MutableS2ShapeIndex index;
  auto shape = make_unique<S2EdgeVectorShape>();
  for (int i = 0; i <= kGridSize; ++i) {
    double e = (i == 0 || i == kGridSize) ? 0 : epsilon;
    shape->Add(S2LatLng::FromDegrees(-e, i).ToPoint(),
               S2LatLng::FromDegrees(kGridSize + e, i).ToPoint());
    shape->Add(S2LatLng::FromDegrees(i, -e).ToPoint(),
               S2LatLng::FromDegrees(i, kGridSize + e).ToPoint());
  }
  index.Add(std::move(shape));
  TestGetCrossingEdgePairs(index, CrossingType::ALL, 112);
  TestGetCrossingEdgePairs(index, CrossingType::INTERIOR, 108);
}

TEST(GetCrossingEdgePairs, EdgeGridTwoIndexes) {
  constexpr int kGridSize = 10;
  double epsilon = 1e-10;

  MutableS2ShapeIndex indexA;
  MutableS2ShapeIndex indexB;
  auto shapeA = make_unique<S2EdgeVectorShape>();
  auto shapeB = make_unique<S2EdgeVectorShape>();
  for (int i = 0; i <= kGridSize; ++i) {
    double e = (i == 0 || i == kGridSize) ? 0 : epsilon;
    shapeA->Add(S2LatLng::FromDegrees(-e, i).ToPoint(),
                S2LatLng::FromDegrees(kGridSize + e, i).ToPoint());
    shapeB->Add(S2LatLng::FromDegrees(i, -e).ToPoint(),
                S2LatLng::FromDegrees(i, kGridSize + e).ToPoint());
  }
  indexA.Add(std::move(shapeA));
  indexB.Add(std::move(shapeB));
  // See comments on the previous test regarding the number of crossings.
  TestGetCrossingEdgePairs(indexA, indexB, CrossingType::ALL, 112);
  TestGetCrossingEdgePairs(indexA, indexB, CrossingType::INTERIOR, 108);
}

// Return true if any loop crosses any other loop (including vertex crossings
// and duplicate edges), or any loop has a self-intersection (including
// duplicate vertices).
static bool HasSelfIntersection(const MutableS2ShapeIndex& index) {
  S2Error error;
  if (s2shapeutil::FindSelfIntersection(index, &error)) {
    ABSL_VLOG(1) << error;
    return true;
  }
  return false;
}

// This function recursively verifies that HasCrossing returns the given
// result for all possible cyclic permutations of the loop vertices for the
// given set of loops.
void TestHasCrossingPermutations(vector<unique_ptr<S2Loop>>* loops, int i,
                                 bool has_crossing) {
  if (i == loops->size()) {
    MutableS2ShapeIndex index;
    S2Polygon polygon(std::move(*loops), S2Debug::DISABLE);
    index.Add(make_unique<S2Polygon::Shape>(&polygon));
    EXPECT_EQ(has_crossing, HasSelfIntersection(index));
    *loops = polygon.Release();
  } else {
    unique_ptr<S2Loop> orig_loop = std::move((*loops)[i]);
    for (int j = 0; j < orig_loop->num_vertices(); ++j) {
      vector<S2Point> vertices;
      for (int k = 0; k < orig_loop->num_vertices(); ++k) {
        vertices.push_back(orig_loop->vertex(j + k));
      }
      (*loops)[i] = make_unique<S2Loop>(vertices, S2Debug::DISABLE);
      TestHasCrossingPermutations(loops, i+1, has_crossing);
    }
    (*loops)[i] = std::move(orig_loop);
  }
}

// Given a string representing a polygon, and a boolean indicating whether this
// polygon has any self-intersections or loop crossings, verify that
// HasSelfIntersection returns the expected result for all possible cyclic
// permutations of the loop vertices.
void TestHasCrossing(string_view polygon_str, bool has_crossing) {
  // Set S2Debug::DISABLE to allow invalid polygons.
  unique_ptr<S2Polygon> polygon =
      s2textformat::MakePolygonOrDie(polygon_str, S2Debug::DISABLE);
  vector<unique_ptr<S2Loop>> loops = polygon->Release();
  TestHasCrossingPermutations(&loops, 0, has_crossing);
}

TEST(FindSelfIntersection, Basic) {
  // Coordinates are (lat,lng), which can be visualized as (y,x).
  TestHasCrossing("0:0, 0:1, 0:2, 1:2, 1:1, 1:0", false);
  TestHasCrossing("0:0, 0:1, 0:2, 1:2, 0:1, 1:0", true);  // duplicate vertex
  TestHasCrossing("0:0, 0:1, 1:0, 1:1", true);  // edge crossing
  TestHasCrossing("0:0, 1:1, 0:1; 0:0, 1:1, 1:0", true);  // duplicate edge
  TestHasCrossing("0:0, 1:1, 0:1; 1:1, 0:0, 1:0", true);  // reversed edge
  TestHasCrossing("0:0, 0:2, 2:2, 2:0; 1:1, 0:2, 3:1, 2:0",
                  true);  // vertex crossing
}

}  // namespace s2shapeutil
