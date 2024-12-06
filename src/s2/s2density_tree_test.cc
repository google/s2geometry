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

#include "s2/s2density_tree.h"

#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gtest/gtest.h>
#include "absl/container/btree_map.h"
#include "absl/container/flat_hash_map.h"
#include "absl/log/absl_check.h"
#include "absl/log/absl_log.h"
#include "absl/log/absl_check.h"
#include "absl/log/log_streamer.h"
#include "absl/meta/type_traits.h"
#include "absl/random/bit_gen_ref.h"
#include "absl/random/random.h"
#include "absl/types/span.h"
#include "s2/util/coding/coder.h"
#include "s2/mutable_s2shape_index.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2cell_union.h"
#include "s2/s2coder_testing.h"
#include "s2/s2coords.h"
#include "s2/s2density_tree_internal.h"
#include "s2/s2error.h"
#include "s2/s2latlng.h"
#include "s2/s2lax_polygon_shape.h"
#include "s2/s2lax_polyline_shape.h"
#include "s2/s2loop.h"
#include "s2/s2point.h"
#include "s2/s2point_vector_shape.h"
#include "s2/s2random.h"
#include "s2/s2region_coverer.h"
#include "s2/s2shape.h"
#include "s2/s2shape_index_region.h"
#include "s2/s2testing.h"
#include "s2/s2text_format.h"
#include "s2/s2wrapped_shape.h"

using ::std::make_unique;
using ::std::unique_ptr;
using ::std::vector;
using ::testing::Eq;

absl::btree_map<S2CellId, int64_t> SumToRoot(
    const absl::btree_map<S2CellId, int64_t>& bases) {
  absl::btree_map<S2CellId, int64_t> sum;

  for (const auto& [cell, weight] : bases) {
    for (int level = 0; level <= cell.level(); ++level) {
      if (auto iter = sum.find(cell.parent(level)); iter != sum.end()) {
        iter->second += weight;
      } else {
        sum.insert({cell.parent(level), weight});
      }
    }
  }

  return sum;
}

// Returns a vector of S2CellIds from a density tree.  If only_leaves is true,
// only returns leaf cells.
std::vector<S2CellId> TreeCells(const S2DensityTree& tree,
                                bool only_leaves = false) {
  std::vector<S2CellId> cell_ids;

  S2Error error;
  tree.VisitCells(
      [&](S2CellId id, const S2DensityTree::Cell& cell) {
        if (!only_leaves || !cell.has_children()) {
          cell_ids.emplace_back(id);
        }
        return S2DensityTree::VisitAction::ENTER_CELL;
      },
      &error);
  EXPECT_TRUE(error.ok());

  return cell_ids;
};

bool ExpectDecodedTreesEqual(const absl::btree_map<S2CellId, int64_t>& got,
                             const absl::btree_map<S2CellId, int64_t>& wanted) {
  if (wanted.size() != got.size()) {
    ABSL_LOG(WARNING) << "size mismatch: got " << got.size() << ", wanted "
                      << wanted.size();
    return false;
  }

  for (const auto& [k, v] : wanted) {
    if (auto iter = got.find(k); iter != got.end()) {
      if (iter->second != v) {
        // Item found, but value does not match.
        ABSL_LOG(WARNING) << "value mismatch for cell=" << k.ToString()
                          << ": got " << iter->second << ", wanted " << v;
        return false;
      }
    } else {
      // Item not found.
      ABSL_LOG(WARNING) << "key mismatch: wanted cell=" << k.ToString()
                        << ", but not found";
      return false;
    }
  }

  return true;
}

TEST(ReversibleBytes, ClearWorks) {
  ReversibleBytes out;
  out.AppendBytes("1");
  EXPECT_THAT(out.size(), Eq(1));
  out.Clear();
  EXPECT_THAT(out.size(), Eq(0));
}

TEST(ReversedLengthsWriter, Codec) {
  // Create an output with some data in it before we get started.
  ReversibleBytes out;
  out.AppendBytes("1");

  // Write blocks of length 0, 2, 0, 1.
  ReversedCellEncoder encoder(&out);
  encoder.Next();
  out.AppendBytes("2");
  out.AppendBytes("3");
  encoder.Next();
  encoder.Next();
  out.AppendBytes("4");
  encoder.Next();

  // Write the offsets we have collected, and a mask.  The children have no
  // meaning here, and we already know which children are present, so the mask
  // value doesn't matter.
  int mask = 123;
  encoder.Finish(mask);

  // Read from a reversed copy.  We leave off the last block, so we expect
  // lengths 1, 0, 2.
  std::string in = out.Reversed();
  Decoder decoder(in.data(), in.size());

  uint64_t v;
  ASSERT_TRUE(decoder.get_varint64(&v));
  EXPECT_EQ(v, mask);

  ASSERT_TRUE(decoder.get_varint64(&v));
  EXPECT_EQ(v, 1);

  ASSERT_TRUE(decoder.get_varint64(&v));
  EXPECT_EQ(v, 0);

  ASSERT_TRUE(decoder.get_varint64(&v));
  EXPECT_EQ(v, 2);

  ASSERT_EQ(decoder.get8(), '4');
  ASSERT_EQ(decoder.get8(), '3');
  ASSERT_EQ(decoder.get8(), '2');

  // At this point we're outside the lengths.  We wrote a 1 before we started,
  // so make sure that's all we have left.
  ASSERT_EQ(decoder.get8(), '1');
  EXPECT_EQ(0, decoder.avail());
}

class TreeEncoderTest : public ::testing::Test {
 protected:
  void Put(S2CellId cell, int64_t weight) { encoder_.Put(cell, weight); }

  S2DensityTree BuildTree() {
    S2DensityTree tree;
    encoder_.Build(&tree);
    return tree;
  }

  absl::btree_map<S2CellId, int64_t> BuildAndDecode() {
    S2Error error;
    auto decoded = BuildTree().Decode(&error);
    EXPECT_EQ(error.code(), S2Error::OK);
    return decoded;
  }

  // Returns a new encoder.
  S2DensityTree::TreeEncoder Encoder() const { return {}; }

  void Clear() { encoder_.Clear(); }

 private:
  S2DensityTree::TreeEncoder encoder_;
};

TEST_F(TreeEncoderTest, EncodeEmpty) { EXPECT_TRUE(BuildAndDecode().empty()); }

TEST_F(TreeEncoderTest, EncodeOneFace) {
  Put(S2CellId::FromFace(3), 17);

  auto entries = BuildAndDecode();
  ASSERT_EQ(1, entries.size());

  auto iter = entries.find(S2CellId::FromFace(3));
  ASSERT_NE(iter, entries.end());
  EXPECT_EQ(iter->second, 17);
}

TEST_F(TreeEncoderTest, EncodeOneLeaf) {
  auto expected = SumToRoot(absl::btree_map<S2CellId, int64_t>{
      {S2CellId(S2Point(0, 1, 0)), 123},
  });

  for (const auto& [cell, weight] : expected) {
    Put(cell, weight);
  }

  auto actual = BuildAndDecode();
  EXPECT_TRUE(ExpectDecodedTreesEqual(actual, expected));
}

TEST_F(TreeEncoderTest, EncodeOneBranch) {
  S2CellId split = S2CellId::FromFaceIJ(1, 1 << 10, 2 << 10).parent(10);
  auto expected = SumToRoot(absl::btree_map<S2CellId, int64_t>{
      {split.child_begin(20), 1},
      {split.child_end(20), 17},
  });

  for (const auto& [cell, weight] : expected) {
    Put(cell, weight);
  }

  auto actual = BuildAndDecode();
  EXPECT_TRUE(ExpectDecodedTreesEqual(actual, expected));
}

TEST_F(TreeEncoderTest, EncodeEachFace) {
  absl::btree_map<S2CellId, int64_t> expected;
  for (int i = 0; i < S2CellId::kNumFaces; ++i) {
    expected.insert({S2CellId::FromFace(i), 10 + i});
  }

  for (const auto& [cell, weight] : expected) {
    Put(cell, weight);
  }

  auto actual = BuildAndDecode();
  EXPECT_TRUE(ExpectDecodedTreesEqual(actual, expected));
}

TEST_F(TreeEncoderTest, EncodeRandomBranches) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "ENCODE_RANDOM_BRANCHES",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));
  for (int64_t weight = 1; weight < 1000; ++weight) {
    Clear();

    absl::btree_map<S2CellId, int64_t> expected;
    for (int j = 0; j < 50; ++j) {
      expected.insert({s2random::CellId(bitgen), weight});
    }
    expected = SumToRoot(expected);
    for (const auto& [cell, weight] : expected) {
      Put(cell, weight);
    }
    auto actual = BuildAndDecode();
    EXPECT_TRUE(ExpectDecodedTreesEqual(actual, expected));
  }
}

TEST(S2DensityTreeTest, LimitsToMaxWeight) {
  MutableS2ShapeIndex index;
  index.Add(make_unique<S2PointVectorShape>(
      vector<S2Point>{S2Point(1, 2, 3).Normalize()}));
  index.Add(make_unique<S2PointVectorShape>(
      vector<S2Point>{S2Point(1, 4, 9).Normalize()}));
  index.Add(make_unique<S2PointVectorShape>(
      vector<S2Point>{S2Point(1, 6, 10).Normalize()}));

  S2Error error;
  S2DensityTree tree;
  ASSERT_TRUE(tree.InitToShapeDensity(
      index, [](const S2Shape&) { return S2DensityTree::kMaxWeight; }, 10000,
      30, &error));

  auto parsed = tree.Decode(&error);
  ASSERT_EQ(error.code(), S2Error::OK);
  for (const auto& [cell, weight] : parsed) {
    EXPECT_EQ(S2DensityTree::kMaxWeight, weight);
  }
}

TEST(S2DensityTreeTest, VisitorCancellation) {
  MutableS2ShapeIndex index;
  index.Add(make_unique<S2PointVectorShape>(
      vector<S2Point>{S2Point(1, 2, 3).Normalize()}));

  S2Error error;
  S2DensityTree tree;
  ASSERT_TRUE(tree.InitToVertexDensity(index, 10000, 30, &error))
      << error.message();

  ASSERT_FALSE(tree.VisitCells(
      [](S2CellId, const S2DensityTree::Cell&) {
        return S2DensityTree::VisitAction::STOP;
      },
      &error));
  EXPECT_TRUE(error.ok());
}

TEST(S2DensityTreeTest, VisitUninitializedTree) {
  S2DensityTree tree;

  int cell_count = 0;
  S2Error error;
  tree.VisitCells(
      [&](S2CellId, const S2DensityTree::Cell&) {
        ++cell_count;
        return S2DensityTree::VisitAction::ENTER_CELL;
      },
      &error);

  EXPECT_EQ(0, cell_count);
  EXPECT_TRUE(error.ok()) << error;
}

TEST(S2DensityTreeTest, Encode) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "ENCODE",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  MutableS2ShapeIndex index;
  for (int i = 0; i < 10; ++i) {
    index.Add(make_unique<S2PointVectorShape>(
        vector<S2Point>{s2random::Point(bitgen)}));
  }

  S2DensityTree tree;
  S2Error error;
  tree.InitToVertexDensity(index, 10'000, 20, &error);
  ASSERT_TRUE(error.ok()) << error;

  // Encode once using a buffer managed by the Encoder.
  Encoder encoder;
  tree.Encode(&encoder);
  std::string encoded(encoder.base(), encoder.length());

  // Encode again using our own buffer
  std::string buffer;
  buffer.resize(tree.EncodedSize());
  encoder.reset(buffer.data(), buffer.size());
  tree.Encode(&encoder);

  // The two encoding methods must yield the same result (size and contents).
  EXPECT_EQ(buffer, encoded);
}

TEST(S2DensityTreeTest, S2CoderWorks) {
  auto index = s2textformat::MakeIndexOrDie("0:0 | 1:1 | 2:2 | 3:3 | 4:4 # #");

  S2DensityTree tree;
  S2Error error;
  tree.InitToVertexDensity(*index, 10'000, 20, &error);
  ASSERT_TRUE(error.ok()) << error;

  error = S2Error::Ok();
  auto decoded = s2coding::RoundTrip(S2DensityTree::Coder(), tree, error);

  // Fully decode density trees and compare to check the decoding.
  error = S2Error::Ok();
  absl::btree_map<S2CellId, int64_t> tree_a = tree.Decode(&error);
  EXPECT_TRUE(error.ok());

  error = S2Error::Ok();
  absl::btree_map<S2CellId, int64_t> tree_b = decoded.Decode(&error);
  EXPECT_TRUE(error.ok());

  EXPECT_EQ(tree_a, tree_b);
}

TEST(S2DensityTreeTest, S2CoderWorks_UninitializedTree) {
  S2DensityTree tree;
  S2Error error;

  auto decoded = s2coding::RoundTrip(S2DensityTree::Coder(), tree, error);

  absl::btree_map<S2CellId, int64_t> cell_map = tree.Decode(&error);
  EXPECT_TRUE(error.ok());

  EXPECT_THAT(cell_map, testing::IsEmpty());
}

TEST(S2DensityTreeTest, InitToFeatureDensity) {
  MutableS2ShapeIndex index;

  using Feature = std::string;

  absl::flat_hash_map<const S2Shape*, Feature*> features;
  absl::flat_hash_map<Feature*, int64_t> weights;

  S2Point p = S2LatLng::FromDegrees(5, 5).ToPoint();
  S2Point q = S2LatLng::FromDegrees(-5, 5).ToPoint();

  // Define a feature with 2 shapes.
  Feature two_shapes{"TwoShapes"};
  weights.insert({&two_shapes, 1});
  {
    auto shape_id =
        index.Add(make_unique<S2PointVectorShape>(vector<S2Point>{p}));
    features.insert({index.shape(shape_id), &two_shapes});
  }
  {
    auto shape_id =
        index.Add(make_unique<S2PointVectorShape>(vector<S2Point>{q}));
    features.insert({index.shape(shape_id), &two_shapes});
  }

  // Define another feature with 1 shape.
  Feature one_shapes{"OneShapes"};
  weights.insert({&one_shapes, 5});
  {
    auto shape_id =
        index.Add(make_unique<S2PointVectorShape>(vector<S2Point>{p}));
    features.insert({index.shape(shape_id), &one_shapes});
  }

  S2DensityTree tree;
  S2Error error;
  tree.InitToFeatureDensity<std::string>(
      index,
      [&](const S2Shape& shape) -> const std::string* {
        const auto iter = features.find(&shape);
        return iter != features.end() ? iter->second : nullptr;
      },
      [&](const Feature& feature) { return weights.at(&feature); }, 100, 1,
      &error);
  ASSERT_TRUE(error.ok()) << error;

  auto parsed = tree.Decode(&error);
  ASSERT_TRUE(error.ok()) << error;

  // TwoShapes isn't double counted
  EXPECT_THAT(parsed, testing::UnorderedElementsAre(
                          testing::Pair(S2CellId(p).parent(0), 6),
                          testing::Pair(S2CellId(p).parent(1), 6),
                          testing::Pair(S2CellId(q).parent(1), 1)));
}

TEST(S2DensityTreeTest, CanNormalizeTree) {
  constexpr int kNumPoints = 1000;

  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "CAN_NORMALIZE_TREE",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  MutableS2ShapeIndex index;
  for (int i = 0; i < kNumPoints; ++i) {
    index.Add(make_unique<S2PointVectorShape>(
        vector<S2Point>{s2random::Point(bitgen)}));
  }

  // Checks that each cell's weight is the sum of its children.
  const auto ExpectNormalized = [](const S2DensityTree& tree) {
    S2DensityTree::DecodedPath path(&tree);
    S2Error error;

    tree.VisitCells(
        [&](S2CellId id, const S2DensityTree::Cell& cell) {
          // Sum up child weights.
          if (cell.has_children()) {
            int64_t child_weight = 0;
            for (int i = 0; i < 4; ++i) {
              const S2DensityTree::Cell* child =
                  path.GetCell(id.child(i), &error);
              EXPECT_TRUE(error.ok());
              child_weight += child->weight();
            }

            EXPECT_TRUE(cell.weight() == child_weight ||
                        cell.weight() + 1 == child_weight);
          }
          return S2DensityTree::VisitAction::ENTER_CELL;
        },
        &error);
  };

  S2DensityTree tree;
  S2Error error;
  tree.InitToVertexDensity(index, 10'000, 20, &error);
  ASSERT_TRUE(error.ok()) << error;
  ASSERT_GT(TreeCells(tree).size(), kNumPoints);

  // Normalize the tree.  The cells shouldn't change but the weights should be
  // normalized now.
  S2DensityTree normalized = tree.Normalize(&error);
  EXPECT_TRUE(error.ok());
  EXPECT_THAT(TreeCells(tree), Eq(TreeCells(normalized)));
  ExpectNormalized(normalized);
}

TEST_F(TreeEncoderTest, NormalizeBalances) {
  // Tests that children with more weight than their parent are linearly
  // rebalanced.
  const S2CellId kFace0 = S2CellId::FromFace(0);

  S2DensityTree tree, expected;
  {
    const absl::btree_map<S2CellId, int64_t> leaves = {
        {kFace0, 3}, {kFace0.child(0), 2}, {kFace0.child(1), 4}};

    auto encoder = Encoder();
    for (const auto& weighted_cell : SumToRoot(leaves)) {
      encoder.Put(weighted_cell.first, weighted_cell.second);
    }
    encoder.Build(&tree);
  }

  {
    const absl::btree_map<S2CellId, int64_t> leaves = {
        {kFace0, 9}, {kFace0.child(0), 3}, {kFace0.child(1), 6}};

    auto encoder = Encoder();
    for (const auto& weighted_cell : SumToRoot(leaves)) {
      encoder.Put(weighted_cell.first, weighted_cell.second);
    }
    encoder.Build(&expected);
  }

  S2Error error;
  tree = tree.Normalize(&error);
  ASSERT_TRUE(error.ok());

  EXPECT_EQ(TreeCells(expected), TreeCells(tree));
}

TEST_F(TreeEncoderTest, NormalizeDoesNotAffectDisjointPaths) {
  // Tests that 3 disjoint paths are unaffected by normalize.
  const S2CellId kFace0 = S2CellId::FromFace(0);

  const absl::btree_map<S2CellId, int64_t> leaves = {
      {kFace0.child(0), 1},
      {kFace0.child(1).child(2), 1},
      {kFace0.child(2), 1}};

  auto encoder = Encoder();
  for (const auto& weighted_cell : SumToRoot(leaves)) {
    encoder.Put(weighted_cell.first, weighted_cell.second);
  }

  S2DensityTree tree;
  encoder.Build(&tree);

  S2Error error;
  S2DensityTree normalized = tree.Normalize(&error);
  ASSERT_TRUE(error.ok());

  EXPECT_EQ(TreeCells(tree), TreeCells(normalized));
}

TEST_F(TreeEncoderTest, NormalizeDoesNotOverflow) {
  // Tests that perfectly divided large weights normalize without overflow.a
  const S2CellId kFace0 = S2CellId::FromFace(0);

  const int64_t kMax32 = std::numeric_limits<int32_t>::max();
  const int64_t kMax64 = std::numeric_limits<int64_t>::max();

  const absl::btree_map<S2CellId, int64_t> leaves = {
      {kFace0.child(1).child(2), kMax32},
      {kFace0.child(1).child(3), kMax64 - kMax32 - 1},
      {kFace0.child(2), 1}};

  auto encoder = Encoder();
  for (const auto& weighted_cell : SumToRoot(leaves)) {
    encoder.Put(weighted_cell.first, weighted_cell.second);
  }

  S2DensityTree tree;
  encoder.Build(&tree);

  S2Error error;
  S2DensityTree normalized = tree.Normalize(&error);
  ASSERT_TRUE(error.ok());

  EXPECT_EQ(TreeCells(tree), TreeCells(normalized));
}

TEST(S2DensityTreeTest, LeavesReturnsLeavesOfTree) {
  constexpr int kNumPoints = 1000;

  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "LEAVES_ARE_SUBSET",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  MutableS2ShapeIndex index;
  for (int i = 0; i < kNumPoints; ++i) {
    index.Add(make_unique<S2PointVectorShape>(
        vector<S2Point>{s2random::Point(bitgen)}));
  }

  S2DensityTree tree;
  S2Error error;
  tree.InitToVertexDensity(index, 10'000, 20, &error);
  ASSERT_TRUE(error.ok()) << error;
  ASSERT_GT(TreeCells(tree).size(), kNumPoints);

  S2CellUnion leaves = tree.Leaves(&error);
  EXPECT_TRUE(error.ok());
  EXPECT_THAT(leaves.cell_ids(), Eq(TreeCells(tree, true)));
}

class DecodedPathTest : public TreeEncoderTest {};

TEST_F(DecodedPathTest, DecoderScalesWeightsBasedOnParent) {
  absl::btree_map<S2CellId, int64_t> base;

  // Create a tree where 4 leaves all share the same weight as their parent.
  S2CellId parent = S2CellId::FromFacePosLevel(0, 0, 5);
  base.insert({parent, 100});
  for (const auto& weighted_cell : SumToRoot(base)) {
    Put(weighted_cell.first, weighted_cell.second);
  }
  for (int i = 0; i < 4; ++i) {
    Put(parent.child(i), 100);
  }
  auto tree = BuildTree();

  S2DensityTree::DecodedPath decoder(&tree);
  S2Error error;

  // Children all have the same weight normalized as 25% of the parent's weight.
  for (int i = 0; i < 4; ++i) {
    const int64_t normal_weight =
        tree.GetNormalCellWeight(parent.child(i), &decoder, &error);
    ASSERT_TRUE(error.ok()) << error;
    EXPECT_EQ(25, normal_weight);

    const int64_t weight =
        tree.GetCellWeight(parent.child(i), &decoder, &error);
    ASSERT_TRUE(error.ok()) << error;
    EXPECT_EQ(100, weight);
  }
}

TEST_F(DecodedPathTest, DecodesPathsCorrectly) {
  const S2CellId kFace0 = S2CellId::FromFace(1);
  const S2CellId kFace2 = S2CellId::FromFace(2);

  S2CellId kCell22 = kFace2.child(2);

  const absl::btree_map<S2CellId, int64_t> base = {{kCell22.child(2), 100},
                                                   {kCell22.child(3), 120}};

  for (const auto& weighted_cell : SumToRoot(base)) {
    Put(weighted_cell.first, weighted_cell.second);
  }

  // Create cell decoder wrapped around a density tree.
  S2DensityTree tree = BuildTree();
  S2DensityTree::DecodedPath decoder(&tree);
  S2Error error;

  // Asserts that the given cell id has the given weight.
  const auto ExpectWeight = [&](S2CellId id, int weight) {
    S2Error error;
    EXPECT_EQ(decoder.GetCell(id, &error)->weight(), weight);
    EXPECT_TRUE(error.ok());
  };

  // Returns a random descendent of the given cell id.
  const auto RandomDescendant = [](absl::BitGenRef bitgen, S2CellId id) {
    ABSL_CHECK_LT(id.level(), S2::kMaxCellLevel);
    int cnt = absl::Uniform(bitgen, 0, S2::kMaxCellLevel - (id.level() + 1));
    for (int i = 0; i < cnt; ++i) {
      id = id.child(absl::Uniform(bitgen, 0, 4));
    }
    return id;
  };

  // Face ids not in the tree should return an empty cell.
  for (int face = 0; face < 6; ++face) {
    if (face == 2) {
      continue;
    }
    ExpectWeight(S2CellId::FromFace(face), 0);
  }

  // But face two is in the tree so we should get it's summed children's weight.
  ExpectWeight(S2CellId::FromFace(2), 220);

  // Children of a face not in the tree should be empty.
  ExpectWeight(kFace0.child(0), 0);
  ExpectWeight(kFace0.child(1), 0);

  // Level 1 children of a face in the tree.
  ExpectWeight(kFace2.child(2), 220);
  ExpectWeight(kFace2.child(3), 0);

  // Cell 22 child 2 and 3 are in the tree, its other children are not.
  ExpectWeight(kCell22.child(0), 0);
  ExpectWeight(kCell22.child(1), 0);
  ExpectWeight(kCell22.child(2), 100);
  ExpectWeight(kCell22.child(3), 120);

  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "DECODES_PATH_CORRECTLY",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  // Cells not in the tree should resolve to their deepest ancestral leaf cell.
  // Therefore, a random descendant of a non-leaf cell should return zero.
  for (int iter = 0; iter < 100; ++iter) {
    ExpectWeight(RandomDescendant(bitgen, kFace2.child(3)), 0);
  }

  // But random descendants of a leaf cell should resolve to that cell.
  for (int iter = 0; iter < 100; ++iter) {
    ExpectWeight(RandomDescendant(bitgen, kCell22.child(2)), 100);
    ExpectWeight(RandomDescendant(bitgen, kCell22.child(3)), 120);
  }
}

class GetPartitioningTest : public TreeEncoderTest {};

TEST_F(GetPartitioningTest, RemovesPointlessSplits) {
  absl::btree_map<S2CellId, int64_t> base{
      {S2CellId::FromFacePosLevel(0, 0, 4), 20},
  };
  for (const auto& weighted_cell : SumToRoot(base)) {
    Put(weighted_cell.first, weighted_cell.second);
  }

  // Add children with all the same weight as the parent.
  for (int i = 0; i < 4; ++i) {
    Put(S2CellId::FromFacePosLevel(0, 0, 4).child(i), 20);
  }

  auto tree = BuildTree();

  S2Error error;
  vector<S2CellUnion> partitioning = tree.GetPartitioning(100, &error);
  ASSERT_TRUE(error.ok());

  for (const auto& covering : partitioning) {
    for (const auto cell : covering) {
      // Level 5 children has been removed.
      EXPECT_EQ(4, cell.level());
    }
  }
}

TEST_F(GetPartitioningTest, ReplacesChildrenWithParent) {
  absl::btree_map<S2CellId, int64_t> base{
      {S2CellId::FromFacePosLevel(0, 0, 4), 20},
      {S2CellId::FromFacePosLevel(1, 0, 4), 40},
  };
  for (const auto& weighted_cell : SumToRoot(base)) {
    Put(weighted_cell.first, weighted_cell.second);
  }

  // On face 0, add an additional level of children. We expect these will be
  // removed because the parent isn't too large.
  for (int i = 0; i < 4; ++i) {
    Put(S2CellId::FromFacePosLevel(0, 0, 4).child(i), 18);
  }

  // On face 1, add an additional level of children. We expect these will NOT be
  // merged because the parent is too large.
  for (int i = 0; i < 4; ++i) {
    Put(S2CellId::FromFacePosLevel(1, 0, 4).child(i), 18);
  }

  auto tree = BuildTree();

  S2Error error;
  vector<S2CellUnion> partitioning = tree.GetPartitioning(100, &error);
  ASSERT_TRUE(error.ok());

  for (const auto& covering : partitioning) {
    for (const auto cell : covering) {
      if (cell.face() == 0) {
        // Cells are on level 4 based on merging the children into the parent.
        EXPECT_EQ(cell.level(), 4);
      } else if (cell.face() == 1) {
        // Cells are on level 5 based on NOT merging the children into the
        // parent.
        EXPECT_EQ(cell.level(), 5);
      } else {
        ADD_FAILURE() << "Unexpected cell face in partitioning: "
                      << cell.face();
      }
    }
  }
}

TEST_F(GetPartitioningTest, OversizeCells) {
  absl::btree_map<S2CellId, int64_t> base;
  for (int i = 0; i < S2CellId::kNumFaces; ++i) {
    base.insert({S2CellId::FromFacePosLevel(i, 0, 10), 1000});
  }

  for (const auto& weighted_cell : SumToRoot(base)) {
    Put(weighted_cell.first, weighted_cell.second);
  }
  auto tree = BuildTree();

  S2Error error;
  vector<S2CellUnion> partitioning = tree.GetPartitioning(10, &error);

  // Each element of the partitioning is a single oversized cell.
  EXPECT_EQ(partitioning.size(), S2CellId::kNumFaces);
  for (const auto& cover : partitioning) {
    EXPECT_EQ(cover.num_cells(), 1);
  }
}

class SumDensityTreesTest : public ::testing::TestWithParam<bool> {
 public:
  SumDensityTreesTest() = default;

  void Insert(S2DensityTree::TreeEncoder& encoder, S2CellId cell) {
    if (auto iter = weights_.find(cell); iter != weights_.end()) {
      encoder.Put(cell, iter->second);
      if (!cell.is_leaf()) {
        for (int i = 0; i < 4; ++i) {
          Insert(encoder, cell.child(i));
        }
      }
    }
  }

  S2DensityTree BuildTree(
      int64_t approx_size_bytes, int max_level,
      S2DensityTree::BreadthFirstTreeBuilder::CellWeightFunction sum,
      S2Error* error) {
    S2DensityTree tree;
    S2DensityTree::TreeEncoder encoder;
    S2DensityTree::BreadthFirstTreeBuilder(approx_size_bytes, max_level,
                                           encoder)
        .Build(sum, &tree, error);
    return tree;
  }

  void CheckSum(const absl::btree_map<S2CellId, int64_t>& expected,
                absl::Span<const S2CellId> roots, int max_level = 30) {
    S2Error error;
    vector<S2DensityTree> trees;
    for (S2CellId root : roots) {
      S2DensityTree::TreeEncoder encoder;

      // Add weights under the given roots
      Insert(encoder, root);

      // And then if this isn't a leaf cell, insert the root weight into all
      // cells above it.
      int64_t weight = weights_.find(root)->second;
      while (root.level() > 0) {
        root = root.parent();
        encoder.Put(root, weight);
      }
      encoder.Build(&trees.emplace_back());

      auto parsed = trees.back().Decode(&error);
    }

    vector<const S2DensityTree*> tree_ptrs;
    for (auto& tree : trees) {
      tree_ptrs.push_back(&tree);
    }

    S2DensityTree sum_tree;

    if (GetParam()) {
      sum_tree.InitToSumDensity(tree_ptrs, 1'000, max_level, &error);
    } else {
      sum_tree.InitToSumDensity(tree_ptrs, max_level, &error);
    }

    auto sum_parsed = sum_tree.Decode(&error);

    EXPECT_TRUE(ExpectDecodedTreesEqual(sum_parsed, expected));
  }

 protected:
  // clang-format off
  absl::btree_map<S2CellId, int64_t> weights_{
    {S2CellId::FromFace(1), 3},
      {S2CellId::FromFace(1).child(1), 1},
      {S2CellId::FromFace(1).child(2), 2},
      {S2CellId::FromFacePosLevel(1, 0, 30), 4},
    {S2CellId::FromFace(2), 4},
    {S2CellId::FromFace(3), 2},
      {S2CellId::FromFace(3).child(0), 2},
      {S2CellId::FromFacePosLevel(3, 0, 30), 2},
  };
  // clang-format on
};

TEST_P(SumDensityTreesTest, SumEmpty) { CheckSum({}, {}); }

TEST_P(SumDensityTreesTest, SumOne) {
  CheckSum(
      {
          {S2CellId::FromFace(1), 3},
          {S2CellId::FromFace(1).child(1), 1},
          {S2CellId::FromFace(1).child(2), 2},
      },
      {
          S2CellId::FromFace(1),
      });
}

TEST_P(SumDensityTreesTest, SumNested) {
  CheckSum(
      {
          // 3, 1, 2 + 1, 1, 0 = 4, 2, 2
          {S2CellId::FromFace(1), 4},
          {S2CellId::FromFace(1).child(1), 2},
          {S2CellId::FromFace(1).child(2), 2},
      },
      {
          S2CellId::FromFace(1),
          S2CellId::FromFace(1).child(1),
      });
}

TEST_P(SumDensityTreesTest, SumDisjoint) {
  CheckSum(
      {
          {S2CellId::FromFace(2), 4},
          {S2CellId::FromFace(3), 2},
          {S2CellId::FromFace(3).child(0), 2},
      },
      {
          S2CellId::FromFace(2),
          S2CellId::FromFace(3),
      });
}

TEST_P(SumDensityTreesTest, SumLeaves) {
  CheckSum(SumToRoot({
               {S2CellId::FromFacePosLevel(1, 0, 30), 4},
               {S2CellId::FromFacePosLevel(3, 0, 30), 2},
           }),
           {
               S2CellId::FromFacePosLevel(1, 0, 30),
               S2CellId::FromFacePosLevel(3, 0, 30),
           });
}

TEST_P(SumDensityTreesTest, SumLeavesLevelLimited) {
  CheckSum(SumToRoot({
               {S2CellId::FromFacePosLevel(1, 0, 20), 4},
               {S2CellId::FromFacePosLevel(3, 0, 20), 2},
           }),
           {
               S2CellId::FromFacePosLevel(1, 0, 30),
               S2CellId::FromFacePosLevel(3, 0, 30),
           },
           20);
}

TEST_P(SumDensityTreesTest, SumMaxLevel) {
  const S2CellId cell = S2CellId::FromFace(5).child(2).child(1).child(0);

  S2Error error;
  for (int max_level = 0; max_level <= cell.level(); ++max_level) {
    auto tree = BuildTree(
        10'000, max_level,
        [&](S2CellId cell_id, S2Error* error) {
          return cell_id.intersects(cell);
        },
        &error);
    ASSERT_EQ(error.code(), S2Error::OK);

    auto actual = tree.Decode(&error);
    ASSERT_EQ(error.code(), S2Error::OK);

    EXPECT_TRUE(ExpectDecodedTreesEqual(
        actual, SumToRoot({{cell.parent(max_level), 1}})));
  }
}

TEST_P(SumDensityTreesTest, SumEmptyAndNonEmpty) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "SUM_EMPTY_AND_NON_EMPTY",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  S2DensityTree empty_tree;
  S2DensityTree tree;

  MutableS2ShapeIndex index;
  index.Add(make_unique<S2PointVectorShape>(
      vector<S2Point>{s2random::Point(bitgen)}));

  S2Error error;
  tree.InitToVertexDensity(index, 1'000, 10, &error);
  ASSERT_TRUE(error.ok()) << error;

  vector<const S2DensityTree*> trees{&tree, &empty_tree};

  S2DensityTree sum_tree;
  sum_tree.InitToSumDensity(trees, 1'000, 10, &error);
  ASSERT_TRUE(error.ok()) << error;

  auto decoded = tree.Decode(&error);
  ASSERT_TRUE(error.ok()) << error;

  // The sum tree should be equivalent to the original non-empty tree.
  EXPECT_FALSE(decoded.empty());
  EXPECT_EQ(decoded, sum_tree.Decode(&error));
}

INSTANTIATE_TEST_SUITE_P(SumTrees, SumDensityTreesTest, testing::Bool(),
                         [](const testing::TestParamInfo<bool>& info) {
                           return info.param ? "WithSizeLimit"
                                             : "WithoutSizeLimit";
                         });

class CoveringsTest : public ::testing::Test {
 protected:
  CoveringsTest() = default;

  // Verifies the given shape/weight pairs' density measures against the
  // intersects/contains results from S2ShapeIndexRegion, which has the same
  // semantics.
  void CheckCoverings(
      vector<std::pair<unique_ptr<S2Shape>, int64_t>> shape_weights) {
    absl::btree_map<const S2Shape*, int64_t> weightmap;
    MutableS2ShapeIndex index;
    for (auto& [shape, weight] : shape_weights) {
      weightmap.insert({shape.get(), weight});
      index.Add(std::move(shape));
    }

    S2ShapeIndexRegion<MutableS2ShapeIndex> region(&index);
    S2RegionCoverer::Options options;
    options.set_max_cells(64);
    S2RegionCoverer coverer(options);
    S2CellUnion cover = coverer.GetCovering(region);

    // Verify cover cells always intersect.
    S2DensityTree::IndexCellWeightFunction measure(
        &index, [&](const S2Shape& shape) -> int64_t {
          if (auto iter = weightmap.find(&shape); iter != weightmap.end()) {
            return iter->second;
          }
          return 0;
        });
    S2Error error;
    for (const S2CellId cell : cover) {
      int64_t weight = Weight(weightmap, cell);
      if (region.Contains(S2Cell(cell))) {
        weight = -weight;
      }
      EXPECT_EQ(weight, measure.WeighCell(cell, &error));
      EXPECT_EQ(error.code(), S2Error::OK);
    }

    // Verify the complement does not intersect.
    S2CellUnion complement(
        coverer.GetCovering(S2Loop(S2Loop::kFull()).GetCapBound()).cell_ids());
    for (const S2CellId cell : complement.Difference(
             coverer.GetCovering(MakeS2ShapeIndexRegion(&index)))) {
      // ...unless the cell is on the border between cover and its complement,
      // since both S2ShapeIndexRegion and S2DensityTree use a small amount of
      // outward padding.
      int64_t expected = cover.Intersects(cell) ? Weight(weightmap, cell) : 0;
      EXPECT_EQ(expected, measure.WeighCell(cell, &error));
      EXPECT_EQ(error.code(), S2Error::OK);
    }
  }

 private:
  int64_t Weight(const absl::btree_map<const S2Shape*, int64_t>& weights,
                 const S2CellId cell) const {
    int64_t sum = 0;
    for (const auto [shape, weight] : weights) {
      MutableS2ShapeIndex index;
      index.Add(make_unique<S2WrappedShape>(shape));
      if (MakeS2ShapeIndexRegion(&index).MayIntersect(S2Cell(cell))) {
        sum += weight;
      }
    }
    return sum;
  }
};

TEST_F(CoveringsTest, ShapeIndexEmpty) {
  vector<std::pair<unique_ptr<S2Shape>, int64_t>> weights;
  weights.emplace_back(
      make_unique<S2Loop::OwningShape>(make_unique<S2Loop>(S2Loop::kEmpty())),
      1);

  CheckCoverings(std::move(weights));
}

TEST_F(CoveringsTest, ShapeIndexPoint) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "SHAPE_INDEX_POINT",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  vector<std::pair<unique_ptr<S2Shape>, int64_t>> weights;
  weights.emplace_back(
      make_unique<S2PointVectorShape>(vector<S2Point>{s2random::Point(bitgen)}),
      1);

  CheckCoverings(std::move(weights));
}

TEST_F(CoveringsTest, ShapeIndexLine) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "SHAPE_INDEX_LINE",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  vector<std::pair<unique_ptr<S2Shape>, int64_t>> weights;
  weights.emplace_back(
      make_unique<S2LaxPolylineShape>(S2Testing::MakeRegularPoints(
          s2random::Point(bitgen), S2Testing::KmToAngle(1), 3)),
      1);

  CheckCoverings(std::move(weights));
}

TEST_F(CoveringsTest, ShapeIndexPolygon) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "SHAPE_INDEX_POLYGON",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  vector<std::pair<unique_ptr<S2Shape>, int64_t>> weights;
  weights.emplace_back(
      make_unique<S2LaxPolygonShape>(
          vector<vector<S2Point>>{S2Testing::MakeRegularPoints(
              s2random::Point(bitgen), S2Testing::KmToAngle(1), 5)}),
      1);

  CheckCoverings(std::move(weights));
}

TEST_F(CoveringsTest, ShapeIndexMultiple) {
  absl::BitGen bitgen(S2Testing::MakeTaggedSeedSeq(
      "SHAPE_INDEX_MULTIPLE",
      absl::LogInfoStreamer(__FILE__, __LINE__).stream()));

  vector<std::pair<unique_ptr<S2Shape>, int64_t>> weights;
  weights.emplace_back(
      make_unique<S2PointVectorShape>(vector<S2Point>{s2random::Point(bitgen)}),
      1);
  weights.emplace_back(
      make_unique<S2LaxPolygonShape>(
          vector<vector<S2Point>>{S2Testing::MakeRegularPoints(
              s2random::Point(bitgen), S2Testing::KmToAngle(1), 5)}),
      2);
  weights.emplace_back(
      make_unique<S2LaxPolylineShape>(S2Testing::MakeRegularPoints(
          s2random::Point(bitgen), S2Testing::KmToAngle(1), 3)),
      3);

  CheckCoverings(std::move(weights));
}
