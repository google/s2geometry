// Copyright 2016 Google Inc. All Rights Reserved.
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

#include "s2/s2builderutil_layers.h"

#include <algorithm>
#include <map>
#include <memory>
#include <set>
#include <string>
#include "s2/third_party/absl/base/integral_types.h"
#include <gtest/gtest.h>
#include "s2/third_party/absl/memory/memory.h"
#include "s2/third_party/absl/strings/str_join.h"
#include "s2/s2builderutil_snap_functions.h"
#include "s2/s2debug.h"
#include "s2/s2textformat.h"
#include "s2/util/gtl/stl_util.h"

using absl::MakeUnique;
using s2builderutil::S2PolygonLayer;
using s2builderutil::S2PolylineLayer;
using s2builderutil::S2PolylineVectorLayer;
using s2builderutil::S2PointVectorLayer;
using s2textformat::MakePolyline;
using s2textformat::MakePoint;
using std::map;
using std::set;
using std::unique_ptr;
using std::vector;
using EdgeType = S2Builder::EdgeType;
using PolylineType = S2PolylineVectorLayer::Options::PolylineType;

namespace {

void TestS2Polygon(vector<char const*> const& input_strs,
                   char const* expected_str, EdgeType edge_type) {
  SCOPED_TRACE(edge_type == EdgeType::DIRECTED ? "DIRECTED" : "UNDIRECTED");
  S2Builder builder((S2Builder::Options()));
  S2Polygon output;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(
      &output, S2PolygonLayer::Options(edge_type)));
  for (auto input_str : input_strs) {
    unique_ptr<S2Polygon> input(s2textformat::MakeVerbatimPolygon(input_str));
    builder.AddPolygon(*input);
  }
  S2Error error;
  ASSERT_TRUE(builder.Build(&error));
  // The input strings in tests may not be in normalized form, so we build an
  // S2Polygon and convert it back to a string.
  unique_ptr<S2Polygon> expected(s2textformat::MakePolygon(expected_str));
  EXPECT_EQ(s2textformat::ToString(*expected),
            s2textformat::ToString(output));
}

void TestS2Polygon(vector<char const*> const& input_strs,
                   char const* expected_str) {
  TestS2Polygon(input_strs, expected_str, EdgeType::DIRECTED);
  TestS2Polygon(input_strs, expected_str, EdgeType::UNDIRECTED);
}

void TestS2PolygonUnchanged(char const* input_str) {
  TestS2Polygon(vector<char const*>{input_str}, input_str);
}

TEST(S2PolygonLayer, NoLoops) {
  TestS2PolygonUnchanged("");
}

TEST(S2PolygonLayer, SmallLoop) {
  TestS2PolygonUnchanged("0:0, 0:1, 1:1");
}

TEST(S2PolygonLayer, ThreeLoops) {
  // The second two loops are nested.
  TestS2PolygonUnchanged("0:1, 1:1, 0:0; "
                         "3:3, 3:6, 6:6, 6:3; "
                         "4:4, 4:5, 5:5, 5:4");
}

TEST(S2PolygonLayer, PartialLoop) {
  S2Builder builder((S2Builder::Options()));
  S2Polygon output;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output));
  builder.AddPolyline(*MakePolyline("0:1, 2:3, 4:5"));
  S2Error error;
  EXPECT_FALSE(builder.Build(&error));
  EXPECT_EQ(S2Error::BUILDER_EDGES_DO_NOT_FORM_LOOPS, error.code());
  EXPECT_TRUE(output.is_empty());
}

TEST(S2PolygonLayer, InvalidPolygon) {
  S2Builder builder((S2Builder::Options()));
  S2Polygon output;
  S2PolygonLayer::Options options;
  options.set_validate(true);
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output, options));
  builder.AddPolyline(*MakePolyline("0:0, 0:10, 10:0, 10:10, 0:0"));
  S2Error error;
  EXPECT_FALSE(builder.Build(&error));
  EXPECT_EQ(S2Error::LOOP_SELF_INTERSECTION, error.code());
}

TEST(S2PolygonLayer, DuplicateInputEdges) {
  // Check that S2PolygonLayer can assemble polygons even when there are
  // duplicate edges (after sibling pairs are removed).
  S2Builder builder((S2Builder::Options()));
  S2Polygon output;
  S2PolygonLayer::Options options;
  options.set_validate(true);
  builder.StartLayer(MakeUnique<S2PolygonLayer>(&output, options));
  builder.AddPolyline(*MakePolyline("0:0, 0:2, 2:2, 1:1, 0:2, 2:2, 2:0, 0:0"));
  S2Error error;
  EXPECT_FALSE(builder.Build(&error));
  EXPECT_EQ(S2Error::POLYGON_LOOPS_SHARE_EDGE, error.code());
  ASSERT_EQ(2, output.num_loops());
  unique_ptr<S2Loop> loop0(s2textformat::MakeLoop("0:0, 0:2, 2:2, 2:0"));
  unique_ptr<S2Loop> loop1(s2textformat::MakeLoop("0:2, 2:2, 1:1"));
  EXPECT_TRUE(loop0->Equals(output.loop(0)));
  EXPECT_TRUE(loop1->Equals(output.loop(1)));
}

// Since we don't expect to have any crossing edges, the key for each edge is
// simply the sum of its endpoints.  This key has the advantage of being
// unchanged when the endpoints of an edge are swapped.
using EdgeLabelMap = map<S2Point, set<int32>>;

void AddPolylineWithLabels(S2Polyline const& polyline, EdgeType edge_type,
                           int32 label_begin, S2Builder* builder,
                           EdgeLabelMap *edge_label_map) {
  for (int i = 0; i + 1 < polyline.num_vertices(); ++i) {
    int32 label = label_begin + i;
    builder->set_label(label);
    // With undirected edges, reverse the direction of every other input edge.
    int dir = edge_type == EdgeType::DIRECTED ? 1 : (i & 1);
    builder->AddEdge(polyline.vertex(i + (1 - dir)), polyline.vertex(i + dir));
    S2Point key = polyline.vertex(i) + polyline.vertex(i + 1);
    (*edge_label_map)[key].insert(label);
  }
}

static void TestEdgeLabels(EdgeType edge_type) {
  S2Builder builder((S2Builder::Options()));
  S2Polygon output;
  S2PolygonLayer::LabelSetIds label_set_ids;
  IdSetLexicon label_set_lexicon;
  builder.StartLayer(MakeUnique<S2PolygonLayer>(
      &output, &label_set_ids, &label_set_lexicon,
      S2PolygonLayer::Options(edge_type)));

  // We use a polygon consisting of 3 loops.  The loops are reordered and
  // some of the loops are inverted during S2Polygon construction.
  EdgeLabelMap edge_label_map;
  AddPolylineWithLabels(*MakePolyline("0:0, 9:1, 1:9, 0:0, 2:8, 8:2, "
                                      "0:0, 0:10, 10:10, 10:0, 0:0"),
                        edge_type, 0, &builder, &edge_label_map);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error));
  vector<int> expected_loop_sizes = { 4, 3, 3 };
  ASSERT_EQ(expected_loop_sizes.size(), label_set_ids.size());
  for (int i = 0; i < expected_loop_sizes.size(); ++i) {
    ASSERT_EQ(expected_loop_sizes[i], label_set_ids[i].size());
    for (int j = 0; j < label_set_ids[i].size(); ++j) {
      S2Point key = output.loop(i)->vertex(j) + output.loop(i)->vertex(j + 1);
      set<int32> const& expected_labels = edge_label_map[key];
      ASSERT_EQ(expected_labels.size(),
                label_set_lexicon.id_set(label_set_ids[i][j]).size());
      EXPECT_TRUE(std::equal(
          expected_labels.begin(), expected_labels.end(),
          label_set_lexicon.id_set(label_set_ids[i][j]).begin()));
    }
  }
}

TEST(S2PolygonLayer, DirectedEdgeLabels) {
  TestEdgeLabels(EdgeType::DIRECTED);
}

TEST(S2PolygonLayer, UndirectedEdgeLabels) {
  TestEdgeLabels(EdgeType::UNDIRECTED);
}

TEST(S2PolygonLayer, ThreeLoopsIntoOne) {
  // Three loops (two shells and one hole) that combine into one.
  TestS2Polygon(
      {"10:0, 0:0, 0:10, 5:10, 10:10, 10:5",
       "0:10, 0:15, 5:15, 5:10",
       "10:10, 5:10, 5:5, 10:5"},
      "10:5, 10:0, 0:0, 0:10, 0:15, 5:15, 5:10, 5:5");
}

TEST(S2PolygonLayer, TrianglePyramid) {
  // A big CCW triangle containing 3 CW triangular holes.  The whole thing
  // looks like a pyramid of nine triangles.  The output consists of 6
  // positive triangles with no holes.
  TestS2Polygon(
      {"0:0, 0:2, 0:4, 0:6, 1:5, 2:4, 3:3, 2:2, 1:1",
       "0:2, 1:1, 1:3",
       "0:4, 1:3, 1:5",
       "1:3, 2:2, 2:4"},
      "0:4, 0:6, 1:5; 2:4, 3:3, 2:2; 2:2, 1:1, 1:3; "
      "1:1, 0:0, 0:2; 1:3, 0:2, 0:4; 1:3, 1:5, 2:4");
}

TEST(S2PolygonLayer, ComplexNesting) {
  // A complex set of nested polygons, with the loops in random order and the
  // vertices in random cyclic order within each loop.  This test checks that
  // the order (after S2Polygon::InitNested is called) is preserved exactly,
  // whether directed or undirected edges are used.
  TestS2PolygonUnchanged(
      "47:15, 47:5, 5:5, 5:15; "
      "35:12, 35:7, 27:7, 27:12; "
      "1:50, 50:50, 50:1, 1:1; "
      "42:22, 10:22, 10:25, 42:25; "
      "47:30, 47:17, 5:17, 5:30; "
      "7:27, 45:27, 45:20, 7:20; "
      "37:7, 37:12, 45:12, 45:7; "
      "47:47, 47:32, 5:32, 5:47; "
      "50:60, 50:55, 1:55, 1:60; "
      "25:7, 17:7, 17:12, 25:12; "
      "7:7, 7:12, 15:12, 15:7");
}

TEST(S2PolygonLayer, FiveLoopsTouchingAtOneCommonPoint) {
  // Five nested loops that touch at one common point.
  TestS2PolygonUnchanged("0:0, 0:10, 10:10, 10:0; "
                         "0:0, 1:9, 9:9, 9:1; "
                         "0:0, 2:8, 8:8, 8:2; "
                         "0:0, 3:7, 7:7, 7:3; "
                         "0:0, 4:6, 6:6, 6:4");
}

TEST(S2PolygonLayer, FourNestedDiamondsTouchingAtTwoPointsPerPair) {
  // Four diamonds nested inside each other, where each diamond shares two
  // vertices with the diamond inside it and shares its other two vertices
  // with the diamond that contains it.  The resulting shape looks vaguely
  // like an eye made out of chevrons.
  TestS2Polygon(
      {"0:10, -10:0, 0:-10, 10:0",
       "0:-20, -10:0, 0:20, 10:0",
       "0:-10, -5:0, 0:10, 5:0",
       "0:5, -5:0, 0:-5, 5:0"},
      "10:0, 0:10, -10:0, 0:20; "
      "0:-20, -10:0, 0:-10, 10:0; "
      "5:0, 0:-10, -5:0, 0:-5; "
      "0:5, -5:0, 0:10, 5:0");
}

TEST(S2PolygonLayer, SevenDiamondsTouchingAtOnePointPerPair) {
  // Seven diamonds nested within each other touching at one
  // point between each nested pair.
  TestS2PolygonUnchanged("0:-70, -70:0, 0:70, 70:0; "
                         "0:-70, -60:0, 0:60, 60:0; "
                         "0:-50, -60:0, 0:50, 50:0; "
                         "0:-40, -40:0, 0:50, 40:0; "
                         "0:-30, -30:0, 0:30, 40:0; "
                         "0:-20, -20:0, 0:30, 20:0; "
                         "0:-10, -20:0, 0:10, 10:0");
}

void TestS2Polyline(
    vector<char const*> const& input_strs,
    char const* expected_str, EdgeType edge_type,
    S2Builder::Options const& options = S2Builder::Options()) {
  SCOPED_TRACE(edge_type == EdgeType::DIRECTED ? "DIRECTED" : "UNDIRECTED");
  S2Builder builder(options);
  S2Polyline output;
  builder.StartLayer(MakeUnique<S2PolylineLayer>(
      &output, S2PolylineLayer::Options(edge_type)));
  for (auto input_str : input_strs) {
    builder.AddPolyline(*MakePolyline(input_str));
  }
  S2Error error;
  ASSERT_TRUE(builder.Build(&error));
  EXPECT_EQ(expected_str, s2textformat::ToString(output));
}

// Convenience function that tests both directed and undirected edges.
void TestS2Polyline(
    vector<char const*> const& input_strs, char const* expected_str,
    S2Builder::Options const& options = S2Builder::Options()) {
  TestS2Polyline(input_strs, expected_str, EdgeType::DIRECTED, options);
  TestS2Polyline(input_strs, expected_str, EdgeType::UNDIRECTED, options);
}

void TestS2PolylineUnchanged(char const* input_str) {
  TestS2Polyline(vector<char const*>{input_str}, input_str);
}

TEST(S2PolylineLayer, NoEdges) {
  TestS2Polyline({}, "");
}

TEST(S2PolylineLayer, OneEdge) {
  // Even with undirected edges, S2PolylineLayer prefers to reconstruct edges
  // in their original direction.
  TestS2PolylineUnchanged("3:4, 1:1");
  TestS2PolylineUnchanged("1:1, 3:4");
}

TEST(S2PolylineLayer, StraightLineWithBacktracking) {
  TestS2PolylineUnchanged("0:0, 1:0, 2:0, 3:0, 2:0, 1:0, 2:0, 3:0, 4:0");
}

TEST(S2PolylineLayer, EarlyWalkTerminationWithEndLoop1) {
  // Test that the "early walk termination" code (which is needed by
  // S2PolylineVectorLayer in order to implement idempotency) does not create
  // two polylines when it is possible to assemble the edges into one.
  //
  // This example tests a code path where the early walk termination code
  // should not be triggered at all (but was at one point due to a bug).
  S2Builder::Options options;
  options.set_snap_function(s2builderutil::IntLatLngSnapFunction(2));
  TestS2Polyline({"0:0, 0:2, 0:1"}, "0:0, 0:1, 0:2, 0:1", options);
}

TEST(S2PolylineLayer, EarlyWalkTerminationWithEndLoop2) {
  // This tests a different code path where the walk is terminated early
  // (yield a polyline with one edge), and then the walk is "maximimzed" by
  // appending a two-edge loop to the end.
  TestS2Polyline({"0:0, 0:1", "0:2, 0:1", "0:1, 0:2"},
                 "0:0, 0:1, 0:2, 0:1");
}

TEST(S2PolylineLayer, SimpleLoop) {
  TestS2PolylineUnchanged("0:0, 0:5, 5:5, 5:0, 0:0");
}

TEST(S2PolylineLayer, ManyLoops) {
  // This polyline consists of many overlapping loops that keep returning to
  // the same starting vertex (2:2).  This tests whether the implementation is
  // able to assemble the polyline in the original order.
  TestS2PolylineUnchanged(
      "0:0, 2:2, 2:4, 2:2, 2:4, 4:4, 4:2, 2:2, 4:4, 4:2, 2:2, 2:0, 2:2, "
      "2:0, 4:0, 2:2, 4:2, 2:2, 0:2, 0:4, 2:2, 0:4, 0:2, 2:2, 0:4, 2:2, "
      "0:2, 2:2, 0:0, 0:2, 2:2, 0:0");
}

TEST(S2PolylineLayer, UnorderedLoops) {
  // This test consists of 5 squares that touch diagonally, similar to the 5
  // white squares of a 3x3 chessboard.  The edges of these squares need to be
  // reordered to assemble them into a single unbroken polyline.
  TestS2Polyline({
      "3:3, 3:2, 2:2, 2:3, 3:3",
      "1:0, 0:0, 0:1, 1:1, 1:0",
      "3:1, 3:0, 2:0, 2:1, 3:1",
      "1:3, 1:2, 0:2, 0:1, 1:3",
      "1:1, 1:2, 2:2, 2:1, 1:1",  // Central square
      },
    "3:3, 3:2, 2:2, 2:1, 3:1, 3:0, 2:0, 2:1, 1:1, 1:0, 0:0, "
    "0:1, 1:1, 1:2, 0:2, 0:1, 1:3, 1:2, 2:2, 2:3, 3:3");
}

TEST(S2PolylineLayer, SplitEdges) {
  // Test reconstruction of a polyline where two edges have been split into
  // many pieces by crossing edges.  This example is particularly challenging
  // because (1) the edges form a loop, and (2) the first and last edges are
  // identical (but reversed).  This is designed to test the heuristics that
  // attempt to find the first edge of the input polyline.
  S2Builder::Options options;
  options.set_split_crossing_edges(true);
  options.set_snap_function(s2builderutil::IntLatLngSnapFunction(7));
  TestS2Polyline(
      {"0:10, 0:0, 1:0, -1:2, 1:4, -1:6, 1:8, -1:10, -5:0, 0:0, 0:10"},
      "0:10, 0:9, 0:7, 0:5, 0:3, 0:1, 0:0, 1:0, 0:1, -1:2, 0:3, 1:4, 0:5, "
      "-1:6, 0:7, 1:8, 0:9, -1:10, -5:0, 0:0, 0:1, 0:3, 0:5, 0:7, 0:9, 0:10",
      options);
}

TEST(S2PolylineLayer, SimpleEdgeLabels) {
  S2Builder builder((S2Builder::Options()));
  S2Polyline output;
  S2PolylineLayer::LabelSetIds label_set_ids;
  IdSetLexicon label_set_lexicon;
  builder.StartLayer(MakeUnique<S2PolylineLayer>(
      &output, &label_set_ids, &label_set_lexicon,
      S2PolylineLayer::Options(EdgeType::UNDIRECTED)));
  builder.set_label(5);
  builder.AddPolyline(*MakePolyline("0:0, 0:1, 0:2"));
  builder.push_label(7);
  builder.AddPolyline(*MakePolyline("0:3, 0:2"));
  builder.clear_labels();
  builder.AddPolyline(*MakePolyline("0:3, 0:4, 0:5"));
  builder.set_label(11);
  builder.AddPolyline(*MakePolyline("0:6, 0:5"));
  S2Error error;
  ASSERT_TRUE(builder.Build(&error));
  vector<vector<int32>> expected = {{5}, {5}, {5, 7}, {}, {}, {11}};
  ASSERT_EQ(expected.size(), label_set_ids.size());
  for (int i = 0; i < expected.size(); ++i) {
    ASSERT_EQ(expected[i].size(),
              label_set_lexicon.id_set(label_set_ids[i]).size());
    int j = 0;
    for (int32 label : label_set_lexicon.id_set(label_set_ids[i])) {
      EXPECT_EQ(expected[i][j++], label);
    }
  }
}

TEST(S2PolylineLayer, InvalidPolyline) {
  S2Builder builder((S2Builder::Options()));
  S2Polyline output;
  S2PolylineLayer::Options options;
  options.set_validate(true);
  builder.StartLayer(MakeUnique<S2PolylineLayer>(&output, options));
  vector<S2Point> vertices;
  vertices.push_back(S2Point(1, 0, 0));
  vertices.push_back(S2Point(-1, 0, 0));
  S2Polyline input(vertices, S2Debug::DISABLE);
  builder.AddPolyline(input);
  S2Error error;
  EXPECT_FALSE(builder.Build(&error));
  EXPECT_EQ(S2Error::ANTIPODAL_VERTICES, error.code());
}

void TestS2PolylineVector(
    vector<char const*> const& input_strs,
    vector<char const*> const& expected_strs,
    EdgeType edge_type,
    S2PolylineVectorLayer::Options layer_options =  // by value
    S2PolylineVectorLayer::Options(),
    S2Builder::Options const& builder_options = S2Builder::Options()) {
  layer_options.set_edge_type(edge_type);
  SCOPED_TRACE(edge_type == EdgeType::DIRECTED ? "DIRECTED" : "UNDIRECTED");
  S2Builder builder(builder_options);
  vector<unique_ptr<S2Polyline>> output;
  builder.StartLayer(MakeUnique<S2PolylineVectorLayer>(&output, layer_options));
  for (auto input_str : input_strs) {
    builder.AddPolyline(*MakePolyline(input_str));
  }
  S2Error error;
  ASSERT_TRUE(builder.Build(&error));
  vector<string> output_strs;
  for (auto const& polyline : output) {
    output_strs.push_back(s2textformat::ToString(*polyline));
  }
  EXPECT_EQ(strings::Join(expected_strs, "; "),
            strings::Join(output_strs, "; "));
}

// Convenience function that tests both directed and undirected edges.
void TestS2PolylineVector(
    vector<char const*> const& input_strs,
    vector<char const*> const& expected_strs,
    S2PolylineVectorLayer::Options const& layer_options =
    S2PolylineVectorLayer::Options(),
    S2Builder::Options const& builder_options = S2Builder::Options()) {
  TestS2PolylineVector(input_strs, expected_strs, EdgeType::DIRECTED,
                       layer_options, builder_options);
  TestS2PolylineVector(input_strs, expected_strs, EdgeType::UNDIRECTED,
                       layer_options, builder_options);
}

void TestS2PolylineVectorUnchanged(vector<char const*> const& input_strs) {
  TestS2PolylineVector(input_strs, input_strs);
}

TEST(S2PolylineVectorLayer, NoEdges) {
  TestS2PolylineVectorUnchanged({});
}

TEST(S2PolylineVectorLayer, TwoPolylines) {
  TestS2PolylineVectorUnchanged({"0:0, 1:1, 2:2", "4:4, 3:3"});
}

TEST(S2PolylineVectorLayer, JoiningPolylines) {
  // Check that polylines are joined together when possible, even if they were
  // not adjacent in the input.  For undirected edges, the polyline direction
  // should be chosen such that the first edge of the polyline was added to
  // S2Builder before the last edge of the polyline.
  TestS2PolylineVector({"1:1, 2:2", "3:3, 2:2", "0:0, 1:1"},
                       {"3:3, 2:2", "0:0, 1:1, 2:2"}, EdgeType::DIRECTED);
  TestS2PolylineVector({"1:1, 2:2", "3:3, 2:2", "0:0, 1:1"},
                       {"3:3, 2:2, 1:1, 0:0"}, EdgeType::UNDIRECTED);
}

TEST(S2PolylineVectorLayer, SegmentNetwork) {
  // Test a complex network of polylines that meet at shared vertices.
  TestS2PolylineVectorUnchanged({
      "0:0, 1:1, 2:2",
      "2:2, 2:3, 2:4",
      "2:4, 3:4, 4:4",
      "2:2, 3:2, 4:2",
      "4:2, 4:3, 4:4",
      "1:0, 2:2",
      "0:1, 2:2",
      "5:4, 4:4",
      "4:5, 4:4",
      "2:4, 2:5, 1:5, 1:4, 2:4",
      "4:2, 6:1, 5:0",  // Two nested loops
      "4:2, 7:0, 6:-1",
      "11:1, 11:0, 10:0, 10:1, 11:1"  // Isolated loop
    });
}

TEST(S2PolylineVectorLayer, MultipleIntersectingWalks) {
  // This checks idempotency for directed edges in the case of several
  // polylines that share edges (and that even share loops).  The test
  // happens to pass for undirected edges as well.
  S2PolylineVectorLayer::Options layer_options;
  layer_options.set_polyline_type(PolylineType::WALK);
  vector<char const*> input = {
    "5:5, 5:6, 6:5, 5:5, 5:4, 5:3",
    "4:4, 5:5, 6:5, 5:6, 5:5, 5:6, 6:5, 5:5, 4:5",
    "3:5, 5:5, 5:6, 6:5, 5:5, 5:6, 6:6, 7:7",
  };
  TestS2PolylineVector(input, input, layer_options);
}

TEST(S2PolylineVectorLayer, EarlyWalkTermination) {
  // This checks idempotency for cases where earlier polylines in the input
  // happen to terminate in the middle of later polylines.  This requires
  // building non-maximal polylines.
  S2PolylineVectorLayer::Options layer_options;
  layer_options.set_polyline_type(PolylineType::WALK);
  vector<char const*> input = {
    "0:1, 1:1",
    "1:0, 1:1, 1:2",
    "0:2, 1:2, 2:2",
    "2:1, 2:2, 2:3"
  };
  TestS2PolylineVector(input, input, layer_options);
}

TEST(S2PolylineVectorLayer, InputEdgeStartsMultipleLoops) {
  // A single input edge is split into several segments by removing portions
  // of it, and then each of those segments becomes one edge of a loop.
  S2PolylineVectorLayer::Options layer_options;
  layer_options.set_polyline_type(PolylineType::WALK);
  layer_options.set_sibling_pairs(
      S2PolylineVectorLayer::Options::SiblingPairs::DISCARD);
  S2Builder::Options builder_options;
  builder_options.set_snap_function(s2builderutil::IntLatLngSnapFunction(7));
  vector<char const*> input = {
    "0:10, 0:0",
    "0:6, 1:6, 1:7, 0:7, 0:8",
    "0:8, 1:8, 1:9, 0:9, 0:10",
    "0:2, 1:2, 1:3, 0:3, 0:4",
    "0:0, 1:0, 1:1, 0:1, 0:2",
    "0:4, 1:4, 1:5, 0:5, 0:6",
  };
  vector<char const*> expected = {
    "0:1, 0:0, 1:0, 1:1, 0:1",
    "0:3, 0:2, 1:2, 1:3, 0:3",
    "0:5, 0:4, 1:4, 1:5, 0:5",
    "0:7, 0:6, 1:6, 1:7, 0:7",
    "0:9, 0:8, 1:8, 1:9, 0:9",
  };
  TestS2PolylineVector(input, expected, layer_options, builder_options);
}

TEST(S2PolylineVectorLayer, SimpleEdgeLabels) {
  S2Builder builder((S2Builder::Options()));
  vector<unique_ptr<S2Polyline>> output;
  S2PolylineVectorLayer::LabelSetIds label_set_ids;
  IdSetLexicon label_set_lexicon;
  S2PolylineVectorLayer::Options layer_options;
  layer_options.set_edge_type(EdgeType::UNDIRECTED);
  layer_options.set_duplicate_edges(
      S2PolylineVectorLayer::Options::DuplicateEdges::MERGE);
  builder.StartLayer(MakeUnique<S2PolylineVectorLayer>(
      &output, &label_set_ids, &label_set_lexicon, layer_options));
  builder.set_label(1);
  builder.AddPolyline(*MakePolyline("0:0, 0:1, 0:2"));
  builder.set_label(2);
  builder.AddPolyline(*MakePolyline("0:3, 0:2, 0:1"));
  builder.clear_labels();
  builder.AddPolyline(*MakePolyline("0:4, 0:5"));
  S2Error error;
  ASSERT_TRUE(builder.Build(&error));
  vector<vector<vector<int32>>> expected = {{{1}, {1, 2}, {2}}, {{}}};
  ASSERT_EQ(expected.size(), label_set_ids.size());
  for (int i = 0; i < expected.size(); ++i) {
    ASSERT_EQ(expected[i].size(), label_set_ids[i].size());
    for (int j = 0; j < expected[i].size(); ++j) {
      ASSERT_EQ(expected[i][j].size(),
                label_set_lexicon.id_set(label_set_ids[i][j]).size());
      int k = 0;
      for (int32 label : label_set_lexicon.id_set(label_set_ids[i][j])) {
        EXPECT_EQ(expected[i][j][k++], label);
      }
    }
  }
}

void VerifyS2PointVectorLayerResults(
    S2PointVectorLayer::LabelSetIds const& label_set_ids,
    IdSetLexicon const& label_set_lexicon, vector<S2Point> const& output,
    string const& str_expected_points,
    vector<vector<int32>> const& expected_labels) {
  vector<S2Point> expected_points =
      s2textformat::ParsePoints(str_expected_points);

  ASSERT_EQ(expected_labels.size(), label_set_ids.size());
  for (int i = 0; i < output.size(); ++i) {
    EXPECT_EQ(expected_points[i], output[i]);
    ASSERT_EQ(expected_labels[i].size(),
              label_set_lexicon.id_set(label_set_ids[i]).size());
    int k = 0;
    for (int32 label : label_set_lexicon.id_set(label_set_ids[i])) {
      EXPECT_EQ(expected_labels[i][k++], label);
    }
  }
}

void AddPoint(S2Point p, S2Builder* builder) { builder->AddEdge(p, p); }

TEST(S2PointVectorLayer, MergeDuplicates) {
  S2Builder builder((S2Builder::Options()));
  std::vector<S2Point> output;
  IdSetLexicon label_set_lexicon;
  S2PointVectorLayer::LabelSetIds label_set_ids;
  builder.StartLayer(MakeUnique<S2PointVectorLayer>(
      &output, &label_set_ids, &label_set_lexicon,
      S2PointVectorLayer::Options(
          S2Builder::GraphOptions::DuplicateEdges::MERGE)));

  builder.set_label(1);
  AddPoint(MakePoint("0:1"), &builder);
  AddPoint(MakePoint("0:2"), &builder);
  builder.set_label(2);
  AddPoint(MakePoint("0:1"), &builder);
  AddPoint(MakePoint("0:4"), &builder);
  AddPoint(MakePoint("0:5"), &builder);
  builder.clear_labels();
  AddPoint(MakePoint("0:5"), &builder);
  AddPoint(MakePoint("0:6"), &builder);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error));

  vector<vector<int32>> expected_labels = {{1, 2}, {1}, {2}, {2}, {}};
  string expected_points = "0:1, 0:2, 0:4, 0:5, 0:6";

  VerifyS2PointVectorLayerResults(label_set_ids, label_set_lexicon, output,
                                  expected_points, expected_labels);
}

TEST(S2PointVectorLayer, KeepDuplicates) {
  S2Builder builder((S2Builder::Options()));
  std::vector<S2Point> output;
  IdSetLexicon label_set_lexicon;
  S2PointVectorLayer::LabelSetIds label_set_ids;
  builder.StartLayer(MakeUnique<S2PointVectorLayer>(
      &output, &label_set_ids, &label_set_lexicon,
      S2PointVectorLayer::Options(
          S2Builder::GraphOptions::DuplicateEdges::KEEP)));

  builder.set_label(1);
  AddPoint(MakePoint("0:1"), &builder);
  AddPoint(MakePoint("0:2"), &builder);
  builder.set_label(2);
  AddPoint(MakePoint("0:1"), &builder);
  AddPoint(MakePoint("0:4"), &builder);
  AddPoint(MakePoint("0:5"), &builder);
  builder.clear_labels();
  AddPoint(MakePoint("0:5"), &builder);
  AddPoint(MakePoint("0:6"), &builder);
  S2Error error;
  ASSERT_TRUE(builder.Build(&error));

  vector<vector<int32>> expected_labels = {{1}, {2}, {1}, {2}, {2}, {}, {}};
  string expected_points = "0:1, 0:1, 0:2, 0:4, 0:5, 0:5, 0:6";

  VerifyS2PointVectorLayerResults(label_set_ids, label_set_lexicon, output,
                                  expected_points, expected_labels);
}

TEST(S2PointVectorLayer, Error) {
  S2Builder builder((S2Builder::Options()));
  std::vector<S2Point> output;
  builder.StartLayer(MakeUnique<S2PointVectorLayer>(
      &output, S2PointVectorLayer::Options(
                   S2Builder::GraphOptions::DuplicateEdges::KEEP)));

  AddPoint(MakePoint("0:1"), &builder);
  builder.AddEdge(MakePoint("0:3"), MakePoint("0:4"));
  AddPoint(MakePoint("0:5"), &builder);
  S2Error error;
  EXPECT_FALSE(builder.Build(&error));
  EXPECT_EQ(error.code(), S2Error::INVALID_ARGUMENT);
  EXPECT_EQ(error.text(), "Found non-degenerate edges");

  EXPECT_EQ(2, output.size());
  EXPECT_EQ(MakePoint("0:1"), output[0]);
  EXPECT_EQ(MakePoint("0:5"), output[1]);
}

#if 0
// Sketch of a test that converts a road network into a polygon mesh.
TEST(LaxPolygonVectorLayer, RoadNetwork) {
  S2Builder::Options builder_options;
  builder_options.set_snap_function(s2builderutil::IntLatLngSnapFunction(7));
  builder_options.set_split_crossing_edges(true);
  S2Builder builder(builder_options);
  LaxPolygonVectorLayer::Options layer_options;
  layer_options.set_degenerate_boundaries(
      LaxPolygonVectorLayer::Options::DegenerateBoundaries::KEEP);
  LaxPolygonVector polygons;
  LaxPolygonVectorLayer::LabelSetIds label_set_ids;
  IdSetLexicon label_set_lexicon;
  builder.StartLayer(MakeUnique<LaxPolygonVectorLayer>(
      &polygons, &label_set_ids, &label_set_lexicon, layer_options));
  ValueLexicon<FeatureId> feature_id_lexicon;
  for (auto const& feature : features) {
    builder.set_label(feature_id_lexicon.Add(segment.feature_id()));
    unique_ptr<S2Polyline> polyline(ToS2Polyline(feature.polyline(0)));
    builder.AddPolyline(*polyline);
  }
  S2Error error;
  ASSERT_TRUE(builder.Build(&error));
  // Should be able to walk over the polygons and find the corresponding sets
  // of labels.
}
#endif

}  // namespace
