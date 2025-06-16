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

#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <functional>
#include <iterator>
#include <limits>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "absl/base/nullability.h"
#include "absl/container/btree_map.h"
#include "absl/container/btree_set.h"
#include "absl/container/flat_hash_map.h"
#include "absl/log/absl_check.h"
#include "absl/log/absl_log.h"
#include "absl/numeric/int128.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "s2/util/coding/coder.h"
#include "s2/util/coding/varint.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2cell_union.h"
#include "s2/s2density_tree_internal.h"
#include "s2/s2error.h"
#include "s2/s2shape.h"
#include "s2/s2shape_index.h"
#include "s2/s2shape_index_region.h"
#include "s2/util/math/mathutil.h"

using absl::btree_map;
using absl::btree_set;
using absl::string_view;
using std::string;
using std::vector;

static constexpr string_view kVersion = "S2DensityTree0";
static constexpr int kNumChildrenPerCell = 4;

// A Node associates an S2CellId and its decoded state within an
// S2DensityTree.
class Node {
 public:
  // Constructs a new Node from the given cell_id, cell and decoder.  'cell'
  // may be nullptr, indicating that this Node does not exist in the
  // S2DensityTree.
  constexpr Node(S2CellId cell_id, const S2DensityTree::Cell* cell,
                 S2DensityTree::DecodedPath* cell_path)
      : cell_id_(cell_id), cell_path_(cell_path) {
    if (cell != nullptr) {
      this->cell_ = *cell;
    }
  }

  // Returns true if the Node is valid and represents a cell within the
  // S2DensityTree.
  bool valid() const {
    return cell_id_.is_valid() && cell_.has_value() && cell_->weight() > 0;
  }

  // Returns the positive weight of this Node, or 0 if the
  int64_t weight() const {
    if (!valid()) {
      return 0;
    }

    return cell_->weight();
  }

  // Returns the Node representing the parent of this Node.
  Node parent() const {
    if (!valid() || cell_id_.is_face()) {
      return Node::Invalid();
    }

    S2Error error;
    const S2DensityTree::Cell* parent_cell =
        cell_path_->GetCell(cell_id_.parent(), &error);
    ABSL_DCHECK(error.ok()) << error;
    return Node(cell_id_.parent(), parent_cell, cell_path_);
  }

  // Returns the Node representing the child of this Node at the given index.
  Node child(int child_index) const {
    if (!valid() || cell_id_.is_leaf() ||
        cell_->child_offset(child_index) < 0) {
      return Invalid();
    }

    S2Error error;
    const S2DensityTree::Cell* child_cell =
        cell_path_->GetCell(cell_id_.child(child_index), &error);
    ABSL_DCHECK(error.ok()) << error;
    return Node(cell_id_.child(child_index), child_cell, cell_path_);
  }

  // Returns true if this Node has more than 1 child with a positive weight.
  bool HasMultipleWeightedChildren() const {
    bool found = false;
    for (int i = 0; i < kNumChildrenPerCell; ++i) {
      if (!child(i).valid()) {
        continue;
      }

      if (found) {
        return true;
      } else {
        found = true;
      }
    }

    return false;
  }

  // Returns true if this Node has multiple children that all have the same
  // weight as the node.
  bool AllChildrenHaveSameWeight() const {
    return HasMultipleWeightedChildren() &&
           VisitWeightedChildren([this](const Node& child_node) {
             return child_node.weight() == weight();
           });
  }

  // Visit all children of this Node present in the tree, calling 'predicate' on
  // each. If predicate returns false, visitation stops and this method returns
  // false.
  bool VisitWeightedChildren(std::function<bool(const Node&)> predicate) const {
    for (int i = 0; i < kNumChildrenPerCell; ++i) {
      const Node child_node = child(i);
      if (child_node.valid() && !predicate(child_node)) {
        return false;
      }
    }

    return true;
  }

  // A convenience method around calling S2DensityTree::GetNormalCellWeight.
  int64_t GetNormalCellWeight() const {
    if (!valid()) {
      return 0;
    }

    S2Error error;
    const S2DensityTree& tree = cell_path_->tree();
    const int64_t normal_weight =
        tree.GetNormalCellWeight(cell_id_, *cell_, cell_path_, &error);
    ABSL_DCHECK(error.ok()) << error;
    return normal_weight;
  }

  // Returns the S2CellId associated with this Node.
  S2CellId cell_id() const { return cell_id_; }

  bool operator<(const Node& other) const { return cell_id_ < other.cell_id_; }

  // Because we use S2CellId as our basis for comparison, we follow suit with
  // S2CellId itself and use a linear search instead of a binary search in
  // btree sets and maps.
  using absl_btree_prefer_linear_node_search = std::true_type;

 private:
  // A sentinel value representing an invalid Node.
  static constexpr Node Invalid() {
    return Node(S2CellId::Sentinel(), nullptr, nullptr);
  }

  S2CellId cell_id_;
  std::optional<S2DensityTree::Cell> cell_;
  S2DensityTree::DecodedPath* cell_path_;
};

// S2DensityTree /////////////////////////////////////////

bool S2DensityTree::InitToShapeDensity(const S2ShapeIndex& index,
                                       const ShapeWeightFunction& weight_fn,
                                       int64_t approximate_size_bytes,
                                       int max_level, S2Error* error) {
  ABSL_DCHECK(error != nullptr) << "error must be non-nullptr";
  *error = S2Error::Ok();

  IndexCellWeightFunction index_cell_weight_fn(&index, weight_fn);

  TreeEncoder encoder;
  BreadthFirstTreeBuilder builder(approximate_size_bytes, max_level, encoder);
  return builder.Build(
      [&](const S2CellId cell_id, S2Error* error) {
        return index_cell_weight_fn.WeighCell(cell_id, error);
      },
      this, error);
}

bool S2DensityTree::InitToVertexDensity(const S2ShapeIndex& index,
                                        int64_t approximate_size_bytes,
                                        int max_level, S2Error* error) {
  return InitToShapeDensity(
      index,
      [&](const S2Shape& shape) {
        switch (shape.dimension()) {
          case 0:
            return shape.num_chains();
          case 1:
            return shape.num_chains() + shape.num_edges();
          case 2:
            return shape.num_edges();
        }
        ABSL_LOG(ERROR) << "unexpected shape with " << shape.dimension()
                         << " dimensions";
        return 0;
      },
      approximate_size_bytes, max_level, error);
}

bool S2DensityTree::InitToSumDensity(vector<const S2DensityTree*>& trees,
                                     int64_t approximate_size_bytes,
                                     int max_level, S2Error* error) {
  ABSL_DCHECK(error != nullptr) << "error must be non-nullptr";
  *error = S2Error::Ok();

  vector<DecodedPath> cell_paths;
  cell_paths.reserve(trees.size());
  for (const auto* tree : trees) {
    cell_paths.emplace_back(tree);
  }

  TreeEncoder encoder;
  BreadthFirstTreeBuilder builder(approximate_size_bytes, max_level, encoder);
  return builder.Build(
      [&](S2CellId cell_id, S2Error* error) -> int64_t {
        int64_t sum = 0;
        bool contained = true;

        for (auto& cell_path : cell_paths) {
          const Cell* cell = cell_path.GetCell(cell_id, error);
          if (!error->ok()) {
            return 0;
          }

          sum += cell->weight();
          contained &= !cell->has_children();

          sum = std::min(sum, kMaxWeight);
        }

        return contained ? -sum : sum;
      },
      this, error);
}

bool S2DensityTree::InitToSumDensity(vector<const S2DensityTree*>& trees,
                                     int max_level, S2Error* error) {
  ABSL_DCHECK(error != nullptr) << "error must be non-nullptr";
  *error = S2Error::Ok();

  TreeEncoder encoder;

  for (const auto* tree : trees) {
    tree->VisitCells(
        [&](S2CellId cell_id, const Cell& cell) {
          if (cell_id.level() > max_level) {
            return VisitAction::SKIP_CELL;
          }

          encoder.Put(cell_id, cell.weight());
          return VisitAction::ENTER_CELL;
        },
        error);

    if (!error->ok()) {
      return false;
    }
  }

  encoder.Build(this);
  return true;
}

bool S2DensityTree::VisitCells(const CellVisitor& visitor_fn,
                               S2Error* error) const {
  ABSL_DCHECK(error != nullptr) << "error must be non-nullptr";
  *error = S2Error::Ok();

  for (int face = 0; face < decoded_faces_.size(); ++face) {
    const int64_t offset = decoded_faces_[face];
    if (offset < 0) {
      continue;
    }

    if (!VisitRecursive(visitor_fn, S2CellId::FromFace(face), offset, error)) {
      return false;
    }
  }

  return true;
}

bool S2DensityTree::VisitRecursive(const CellVisitor& visitor_fn,
                                   S2CellId cell_id, int64_t position,
                                   S2Error* error) const {
  Decoder decoder(encoded_.data(), encoded_.size());
  decoder.skip(position);

  Cell cell;
  if (!cell.Decode(decoder, error)) {
    return false;
  }

  // Visit this node.
  switch (visitor_fn(cell_id, cell)) {
    case VisitAction::ENTER_CELL:
      break;
    case VisitAction::SKIP_CELL:
      return true;
    case VisitAction::STOP:
      return false;
  }

  // Visit the children.
  for (int i = 0; i < 4; ++i) {
    const int64_t offset = cell.child_offset(i);
    if (offset < 0) {
      continue;
    }

    if (!VisitRecursive(visitor_fn, cell_id.child(i), offset, error)) {
      return false;
    }
  }

  return true;
}

int64_t S2DensityTree::GetCellWeight(const S2CellId cell_id,
                                     DecodedPath* cell_path,
                                     S2Error* error) const {
  ABSL_DCHECK(error != nullptr) << "error must be non-nullptr";
  *error = S2Error::Ok();

  const Cell* cell = cell_path->GetCell(cell_id, error);
  if (!error->ok()) {
    return 0;
  }

  return cell->weight();
}

int64_t S2DensityTree::GetNormalCellWeight(const S2CellId cell_id,
                                           DecodedPath* cell_path,
                                           S2Error* error) const {
  ABSL_DCHECK(error != nullptr) << "error must be non-nullptr";
  *error = S2Error::Ok();
  ABSL_DCHECK(cell_path != nullptr) << "decoder must be non-nullptr";

  const Cell* cell = cell_path->GetCell(cell_id, error);
  if (!error->ok()) {
    return 0;
  }

  if (cell->weight() == 0) {
    // Cell is not present in this tree.
    return 0;
  }

  return GetNormalCellWeight(cell_id, *cell, cell_path, error);
}

int64_t S2DensityTree::GetNormalCellWeight(const S2CellId cell_id,
                                           const Cell& cell,
                                           DecodedPath* cell_path,
                                           S2Error* error) const {
  double scale = 1.0;
  Node node(cell_id, &cell, cell_path);

  while (node.parent().valid()) {
    const int64_t weight = node.weight();
    node = node.parent();
    int64_t sum = 0;

    for (int i = 0; i < kNumChildrenPerCell; ++i) {
      sum += node.child(i).weight();
    }

    scale *= static_cast<double>(weight) / sum;
  }

  return MathUtil::Round<int64_t>(scale * node.weight());
}

vector<S2CellUnion> S2DensityTree::GetPartitioning(int64_t max_weight,
                                                   S2Error* error) const {
  ABSL_DCHECK(error != nullptr) << "error must be non-nullptr";
  *error = S2Error::Ok();

  // Sample cells at 1/16th of the desired size. This yields a more efficient
  // bin-packing of our S2CellUnions compared to sampling cells closer to the
  // max_weight.
  const int64_t target_weight = max_weight / 16;

  DecodedPath decoder(this);
  btree_set<Node> candidates;

  // Collect an initial set of cells which are either less than target_weight,
  // or single cells with no children (which may be greater than
  // target_weight).
  VisitCells(
      [&](const S2CellId cell_id, const S2DensityTree::Cell& cell) {
        if (!error->ok()) {
          return VisitAction::STOP;
        }

        if (cell.weight() > target_weight && cell.has_children()) {
          return VisitAction::ENTER_CELL;
        }

        candidates.insert(Node(cell_id, &cell, &decoder));
        return VisitAction::SKIP_CELL;
      },
      error);
  if (!error->ok()) {
    return {};
  }

  // Revise the initial set of candidates by looking for optimizations.
  btree_set<Node> nodes;
  for (Node node : candidates) {
    // Skip cells contained by previously computed nodes.
    if (!nodes.empty() &&
        std::prev(nodes.end())->cell_id().intersects(node.cell_id())) {
      continue;
    }

    // While the current node is a pointless split, move up the tree.
    while (node.parent().valid() && node.parent().AllChildrenHaveSameWeight()) {
      node = node.parent();

      // Moving up the tree may cover previously-added nodes, so remove them
      // as we go.
      while (!nodes.empty() &&
             std::prev(nodes.end())->cell_id().intersects(node.cell_id())) {
        nodes.erase(std::prev(nodes.end()));
      }
    }

    // Add the node, potentially after selecting one higher in the tree.
    nodes.insert(node);

    // If all children have been added, replace with the parent if not too
    // much larger.
    for (Node parent = node.parent();
         parent.valid() && parent.GetNormalCellWeight() < max_weight / 4 &&
         parent.HasMultipleWeightedChildren() &&
         parent.VisitWeightedChildren(
             [&](const Node& node) { return nodes.contains(node); });
         parent = parent.parent()) {
      parent.VisitWeightedChildren([&](const Node& n) {
        nodes.erase(n);
        return true;
      });
      nodes.insert(parent);
    }
  }

  vector<S2CellUnion> partitioning;
  vector<S2CellId> cover;

  int64_t current_weight = 0;
  for (const Node& node : nodes) {
    const int64_t normal_weight = node.GetNormalCellWeight();

    if (!cover.empty() && current_weight + normal_weight >= max_weight) {
      partitioning.push_back(S2CellUnion::FromVerbatim(std::move(cover)));
      cover.clear();
      current_weight = 0;
    }

    cover.push_back(node.cell_id());
    current_weight += normal_weight;
  }

  if (!cover.empty()) {
    partitioning.push_back(S2CellUnion::FromVerbatim(std::move(cover)));
  }

  return partitioning;
}

btree_map<S2CellId, int64_t> S2DensityTree::Decode(S2Error* error) const {
  btree_map<S2CellId, int64_t> weights;

  VisitCells(
      [&](const S2CellId cell_id, const Cell& node) {
        weights.insert({cell_id, node.weight()});
        return VisitAction::ENTER_CELL;
      },
      error);

  return weights;
}

bool S2DensityTree::Init(Decoder* decoder, S2Error& error) {
  encoded_ = {decoder->skip(0), decoder->avail()};
  if (encoded_.empty()) {
    // An uninitialized tree is still valid.
    return true;
  }

  return DecodeHeader(decoder, &decoded_faces_, &error);
}

void S2DensityTree::Encode(Encoder* encoder) const {
  // Allow interop with both pre-allocated and self-managed Encoders.
  if (encoder->ensure_allowed()) {
    encoder->Ensure(encoded_.size());
  }
  ABSL_CHECK_GE(encoder->avail(), encoded_.size());

  encoder->putn(encoded_.data(), encoded_.size());
}

// IndexCellWeightFunction ///////////////////////////////////

int64_t S2DensityTree::IndexCellWeightFunction::WeighCell(
    const S2CellId cell_id, S2Error*) {
  int64_t sum = 0;
  bool all_contained = true;

  index_region_.VisitIntersectingShapes(
      S2Cell(cell_id), [&](const S2Shape* shape, bool contains_target) {
        const int64_t weight = weight_fn_(*shape);
        ABSL_DCHECK_GE(weight, 0);
        ABSL_DCHECK_LE(weight, kMaxWeight);
        sum += weight;
        all_contained &= contains_target;
        return true;
      });

  sum = std::min(sum, kMaxWeight);
  return all_contained ? -sum : sum;
}

// BreadthFirstTreeBuilder ///////////////////////////////////

bool S2DensityTree::BreadthFirstTreeBuilder::Build(
    const CellWeightFunction& weight_fn, S2DensityTree* tree,
    S2Error* error) const {
  vector<std::pair<S2CellId, S2CellId>> ranges{{
      S2CellId::Begin(S2CellId::kMaxLevel),
      S2CellId::End(S2CellId::kMaxLevel),
  }};
  vector<std::pair<S2CellId, S2CellId>> next_level_ranges;

  for (int level = 0, size_estimate_bytes = 0;
       !ranges.empty() && level <= max_level_ &&
       size_estimate_bytes < approximate_size_bytes_;
       ++level) {
    S2CellId last_range_end = S2CellId::Sentinel();

    for (auto& range : ranges) {
      for (S2CellId cell_id = range.first.parent(level); cell_id < range.second;
           cell_id = cell_id.next()) {
        // Get the weight and skip this cell_id unless it's larger than 0.
        int64_t weight = weight_fn(cell_id, error);
        if (!error->ok()) {
          return false;
        }

        if (weight == 0) {
          // Skip disjoint cells.
          continue;
        } else if (weight < 0) {
          // Get the absolute weight and skip searching the children.
          weight = -weight;
        } else {
          // Add this hilbert range to the ranges to scan at the next level.
          const S2CellId begin = cell_id.range_min();
          const S2CellId end = cell_id.range_max().next();
          if (begin == last_range_end) {
            // Extend the existing range.
            next_level_ranges.back().second = end;
          } else {
            // Add a new range.
            next_level_ranges.push_back({begin, end});
          }
          last_range_end = end;
        }

        // Save the weight for repacking later and estimate the size it will
        // consume.
        ABSL_DCHECK_LE(weight, kMaxWeight)
            << "CellIdWeightFn produced weight greater than kMaxWeight: "
            << weight;
        encoder_.Put(cell_id, std::min(weight, kMaxWeight));
        size_estimate_bytes += TreeEncoder::EstimateSize(weight);
      }
    }

    ranges = std::move(next_level_ranges);
    // Ensure that the moved-from vector is reset to a known state.
    next_level_ranges.clear();
  }

  encoder_.Build(tree);
  return true;
}

// Cell ///////////////////////////////////////////

bool S2DensityTree::Cell::Decode(Decoder& decoder, S2Error* error) {
  uint64_t bits;
  if (!decoder.get_varint64(&bits)) {
    *error = S2Error::Internal(absl::StrCat(
        "Failed to decode bits for Cell at position ", decoder.pos()));
    return false;
  }

  weight_ = bits >> S2DensityTree::kChildMaskBits;

  std::bitset<kChildMaskBits> child_mask(bits);

  // Ensure the offsets are default-initialized to guard against Cell reuse,
  // which might be done by a DecodedPath cache.
  offsets_.fill(-1);

  // If this is a leaf without children, we can return early (skipping a bunch
  // of branches below).
  if (child_mask.to_ulong() == 0) {
    return true;
  }

  // Read the cumulative offsets for each child in the child mask.
  int64_t offset = 0;
  for (int i = 0; i < offsets_.size(); ++i) {
    if (child_mask[i]) {
      // Child is set.  Set the current offset and discard the bit.
      offsets_[i] = offset;
      child_mask.reset(i);
      if (child_mask.to_ulong() > 0) {
        // There's another child, so increment the cumulative offset by
        // the offset of the next child.
        uint64_t v;
        if (!decoder.get_varint64(&v)) {
          *error = S2Error::Internal(absl::StrCat(
              "Failed to decode child offset at position ", decoder.pos()));
          return false;
        }
        offset += v;
      }
    }
    // Else child is unset. We leave the offset at the initial value of -1.
  }

  // Increment the offsets by the decoder position at the end of the offsets,
  // which is what they are relative to.
  for (auto& pos : offsets_) {
    if (pos >= 0) {
      pos += decoder.pos();
    }
  }

  return true;
}

void S2DensityTree::Cell::Clear() {
  offsets_ = kNoChildren;
  weight_ = 0;
}

bool S2DensityTree::Cell::DecodeAt(const S2DensityTree* tree, uint64_t pos,
                                   S2Error* error) {
  Decoder decoder(tree->encoded_.data(), tree->encoded_.size());
  decoder.skip(pos);
  return Decode(decoder, error);
}

// TreeEncoder ///////////////////////////////////////////
//
// The encoded format is fairly straightforward, starting with a header:
//   * Version: Marks this buffer as an S2DensityTree.
//   * Face mask (varint): 6 bits indicating whether a face is encoded in this
//     buffer.
//   * Face lengths (varint64 array): length of each encoded face for which
//     the face mask bit is set, except the last face length which isn't
//     needed to compute its start position.
//
// Each encoded face is a sequence of encoded cells.  Each encoded cell is:
//   * Masked Weight (varint64): the high 60 bits contains the unsigned weight
//     of the current cell, and the low 4 bits contains the child mask
//     indicating which offsets follow.
//   * Lengths (varint64 array): length of each encoded child cell for which
//     the child mask bit is set, except the last child length which isn't
//     needed to compute its start position.
//
// The face and child masks use bit 0 for child 0, bit 1 for child 1, etc.
// These indicate the child is present, and that the weight is greater than 0.

void S2DensityTree::TreeEncoder::Put(const S2CellId cell_id, int64_t weight) {
  weights_[cell_id] += weight;
}

void S2DensityTree::TreeEncoder::Build(S2DensityTree* tree) {
  ReversibleBytes output;
  EncodeTreeReversed(&output);
  output.AppendBytes(kVersion);
  output.ReverseFrom(output.size() - kVersion.size());

  string bytes = output.Reversed();
  Decoder decoder(bytes.data(), bytes.size());

  S2Error error;
  bool ok = tree->Init(&decoder, error);
  ABSL_DCHECK(ok) << error.message();
}

void S2DensityTree::TreeEncoder::EncodeTreeReversed(ReversibleBytes* output) {
  ReversedCellEncoder lengths(output);
  uint64_t mask = 0;

  for (int face = 5; face >= 0; --face) {
    const S2CellId face_cell = S2CellId::FromFace(face);

    if (auto iter = weights_.find(face_cell); iter != weights_.end()) {
      const int64_t weight = iter->second;

      EncodeSubtreeReversed(face_cell, weight, output);
      lengths.Next();
      mask |= 1 << face;
    }
  }

  lengths.Finish(mask);
}

void S2DensityTree::TreeEncoder::EncodeSubtreeReversed(
    const S2CellId cell_id, int64_t weight, ReversibleBytes* output) {
  ReversedCellEncoder lengths(output);
  uint64_t mask = 0;

  if (!cell_id.is_leaf()) {
    for (int i = 3; i >= 0; --i) {
      const S2CellId child = cell_id.child(i);

      if (const auto iter = weights_.find(child); iter != weights_.end()) {
        const int64_t weight = iter->second;

        EncodeSubtreeReversed(child, weight, output);
        lengths.Next();
        mask |= 1 << i;
      }
    }
  }

  lengths.Finish((weight << S2DensityTree::kChildMaskBits) | mask);
}

// static
int S2DensityTree::TreeEncoder::EstimateSize(int64_t weight) {
  const int weight_size =
      Varint::Length64(weight << S2DensityTree::kChildMaskBits | 0xF);
  return weight_size + 2 * Varint::Length64(weight_size);
}

void S2DensityTree::TreeEncoder::Clear() {
  weights_.erase(weights_.begin(), weights_.end());
}

// static
bool S2DensityTree::DecodeHeader(Decoder* decoder, DecodedFaces* decoded_faces,
                                 S2Error* error) {
  // Verify the version string.  Check the available length first because
  // Decoder::getn() does not check bounds for us.
  if (decoder->avail() < kVersion.size()) {
    *error = S2Error::InvalidArgument(
        "Not enough bytes to decode magic value for S2DensityTree");
    return false;
  }

  char magic[kVersion.size()];
  decoder->getn(magic, sizeof(magic));
  if (string_view(magic, sizeof(magic)) != kVersion) {
    *error = S2Error::InvalidArgument("Bad magic value for S2DensityTree");
    return false;
  }

  // Read the lengths of each cube face except the last.

  uint64_t bits;
  if (!decoder->get_varint64(&bits)) {
    *error = S2Error::InvalidArgument(absl::StrCat(
        "Failed to decode face mask at position %d", decoder->pos()));
    return false;
  }

  const std::bitset<S2CellId::kNumFaces> face_mask(bits);

  // A list of pairings of <face, offset>.  This is a temporary representation
  // needed because the offsets are relative to the decoder positions *after*
  // all of the faces have been decoded.
  vector<std::pair<int, uint64_t>> coded_faces;
  coded_faces.reserve(S2CellId::kNumFaces);

  int64_t offset = 0;
  for (int i = 0; i < face_mask.size(); ++i) {
    if (face_mask.test(i)) {
      coded_faces.push_back({i, offset});

      // There are (face_mask.count()-1) lengths to read.
      if (coded_faces.size() < face_mask.count()) {
        uint64_t v;
        if (!decoder->get_varint64(&v)) {
          ABSL_LOG(ERROR) << "Failed to decode offset at position "
                          << decoder->pos();
          return false;
        }
        offset += v;
      }
    }
  }

  // Convert the face/offset pairs into decoded faces.
  for (const auto& coded_face : coded_faces) {
    const int face = coded_face.first;
    const int64_t offset = coded_face.second;

    decoded_faces->at(face) = decoder->pos() + offset;
  }

  return true;
}

const S2DensityTree::Cell* S2DensityTree::DecodedPath::GetCell(
    const S2CellId cell_id, S2Error* error) {
  // If the new cell is on a different face, load the new face and reset the
  // cell stack before proceeding.
  if (last_.face() != cell_id.face()) {
    last_ = cell_id.parent(0);
    LoadFace(cell_id.face(), error);
    if (!error->ok()) {
      return nullptr;
    }
  }

  // Load cells up to and including the requested level.
  return LoadCell(cell_id, error);
}

void S2DensityTree::DecodedPath::LoadFace(int face, S2Error* error) {
  S2DensityTree::Cell* cell = &stack_[0];

  const int64_t offset = tree_->decoded_faces_[face];
  if (offset < 0) {
    cell->Clear();
  } else {
    cell->DecodeAt(tree_, offset, error);
  }
}

const S2DensityTree::Cell* S2DensityTree::DecodedPath::LoadCell(
    S2CellId cell_id, S2Error* error) {
  const int start_level = last_.GetCommonAncestorLevel(cell_id);
  const int cell_level = cell_id.level();

  // Loading starts at the most recently loaded level.
  S2DensityTree::Cell* cell = &stack_[start_level];
  S2DensityTree::Cell* prev = nullptr;

  int level = start_level + 1;
  for (; level <= cell_level; ++level) {
    // Find the next child position, if available.
    const int64_t offset = cell->child_offset(cell_id.child_position(level));

    // Examine the frame associated with the child.
    prev = cell;
    cell = &stack_[level];

    if (offset < 0) {
      // The cell isn't in the tree. If its ancestor is a leaf node, we can
      // return that, otherwise we'll return an empty cell.
      if (!prev->has_children()) {
        cell = prev;
      } else {
        cell->Clear();
      }
      break;
    } else if (!cell->DecodeAt(tree_, offset, error)) {
      return nullptr;
    }
  }

  last_ = cell_id.parent(level - 1);
  return cell;
}

S2DensityTree S2DensityTree::Normalize(absl::Nonnull<S2Error*> error) const {
  *error = S2Error::Ok();

  DecodedPath path(this);

  absl::flat_hash_map<S2CellId, int64_t> weights;
  VisitCells(
      [&](S2CellId id, const Cell& cell) {
        absl::int128 weight = cell.weight();
        if (!id.is_face()) {
          const S2CellId parent = id.parent();

          absl::int128 sibling_weight = 0;
          for (int i = 0; i < 4; ++i) {
            const Cell* sibling = path.GetCell(parent.child(i), error);
            if (!error->ok()) {
              return VisitAction::STOP;
            }
            sibling_weight += sibling->weight();
          }

          weight = (weight * weights[parent] - 1) / sibling_weight + 1;
        }

        weights[id] = static_cast<int64_t>(weight);
        return VisitAction::ENTER_CELL;
      },
      error);

  TreeEncoder encoder;
  for (const auto& entry : weights) {
    encoder.Put(entry.first, entry.second);
  }

  S2DensityTree tree;
  encoder.Build(&tree);
  return tree;
}

S2CellUnion S2DensityTree::Leaves(absl::Nonnull<S2Error*> error) const {
  std::vector<S2CellId> leaves;

  VisitCells(
      [&](S2CellId cell_id, const Cell& cell) {
        if (cell.has_children()) {
          return VisitAction::ENTER_CELL;
        }
        leaves.emplace_back(cell_id);
        return VisitAction::SKIP_CELL;
      },
      error);

  return S2CellUnion::FromVerbatim(std::move(leaves));
}
