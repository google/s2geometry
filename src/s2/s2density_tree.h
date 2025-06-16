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

#ifndef S2_S2DENSITY_TREE_H_
#define S2_S2DENSITY_TREE_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "absl/base/nullability.h"
#include "absl/container/btree_map.h"
#include "absl/container/flat_hash_set.h"
#include "absl/log/absl_check.h"
#include "s2/util/coding/coder.h"
#include "s2/s2cell_id.h"
#include "s2/s2cell_union.h"
#include "s2/s2coder.h"
#include "s2/s2density_tree_internal.h"
#include "s2/s2error.h"
#include "s2/s2shape.h"
#include "s2/s2shape_index.h"
#include "s2/s2shape_index_region.h"

// A density tree is a kind of spatial histogram, computed over a collection of
// shapes in an S2ShapeIndex.  It's most often used to do spatial clustering
// into equal-size shards, where spatial datasets are well-known for their skew,
// which would defeat a regular tiling of the surface.  But with detailed
// spatial density, we can easily divide the surface into equal-weight shards.
//
// A density tree is more formally a map from S2CellIds to weights with the
// property that if any cell is present in the tree, then all of its ancestors
// are present in the tree.  The weight of each cell is the sum of the weights
// of the shapes that intersect that cell.  (Note that the weight of a parent
// cell is generally not the sum of the weights of its children.)
//
// Weights have no units.  Their meaning is defined by the operation that
// computes them.  For example, the weight of a shape could be the number of
// vertices in each shape.  Generally a client will compute shape density for a
// whole set of shapes at once by putting them in an S2ShapeIndex and providing
// a ShapeWeightFunction to map a shape to the weight of that shape.  If there
// are multiple indices to weigh, intermediate trees can be combined via
// S2DensityTree::InitToSumDensity.
//
// As a limitation, if a shape overlaps several cells at a given level,
// S2DensityTree does not keep track of the fact that the weights in those cells
// are correlated.  So while it is capable of accurately measuring the weight in
// individual cells, it can very significantly overcount the weight in a
// collection of cells such as the S2CellUnions returned by GetPartitioning.
//
// Density trees are fast and cheap to initialize from a Decoder, but the
// entries are parsed lazily.  If lookups will be visited more often than there
// are entries, it may be faster to 'Decode' into a <cell, weight> map that
// has faster lookups, but note the in-memory size is several times larger.
//
// S2DensityTrees can be simplified to keep their size below a threshold.
// Initializers such as InitToShapeDensity accept approximate byte size and cell
// level limits.
//
// As a more complete example:
//
//   MutableS2ShapeIndex index;
//   index.Add(...);
//
//   S2DensityTree tree;
//   S2Error error;
//   ABSL_DCHECK(tree.InitToVertexDensity(index, 10'000, 15, &error))
//       << error.message();
//
// An example of summing trees:
//
//   std::vector<S2DensityTree> trees;
//   // Eliding tree initialization as demonstrated above.
//
//   std::vector<const S2DensityTree*> tree_ptrs;
//   for (const auto& tree : trees) {
//     tree_ptrs.push_back(&tree);
//   }
//
//   S2DensityTree sum_tree;
//   S2Error error;
//   ABSL_DCHECK(sum_tree.InitToSumDensity(tree_ptrs, 100'000, 15, &error))
//       << error.message();
//
// Example of getting equal sized shards from a given tree:
//
//   S2DensityTree tree;
//   // Eliding tree initialization as demonstrated above.
//
//   // Create shards of no more than 100MiB each.
//   S2Error error;
//   auto partitions = tree.GetPartitioning(100 * 2 << 10, &error);
//   ABSL_DCHECK(error.ok()) << error.message();

class S2DensityTree {
 public:
  static constexpr int kChildMaskBits = 4;
  static constexpr int64_t kMaxWeight =
      std::numeric_limits<int64_t>::max() >> kChildMaskBits;

  typedef s2coding::S2BasicCoder<S2DensityTree> Coder;

  // Default constructor. Creates an uninitialized tree without any stored
  // weights. Use one of the below Init* methods to populate the tree.
  S2DensityTree() = default;

  ~S2DensityTree() = default;

  // S2DensityTrees can be freely moved and copied.
  S2DensityTree(const S2DensityTree&) = default;
  S2DensityTree& operator=(const S2DensityTree&) = default;
  S2DensityTree(S2DensityTree&&) = default;
  S2DensityTree& operator=(S2DensityTree&&) = default;

  // Type definition for a function that returns the weight of S2Shapes (see
  // class-level comment for the meaning of weights).  It is an error to return
  // a weight outside of the closed range [0, kMaxWeight].
  using ShapeWeightFunction = std::function<int64_t(const S2Shape&)>;

  // Initialize the S2DensityTree with the given index with 'weight_fn'
  // providing the weights for each shape.  The 'approximate_size_bytes' and
  // 'max_level' act to limit the size of the tree.  Note that the byte size
  // limit will never interrupt the addition of a level.  Only full levels will
  // ever be added to the tree.  This means that the byte size limit may be
  // exceeded by a factor of 4x in the worst case, where deeper trees typically
  // have lower error margins than shallower trees in practice.
  //
  // If this function returns false, the 'error' object can be expected to
  // contain a reason the tree could not be initialized.
  bool InitToShapeDensity(const S2ShapeIndex& index,
                          const ShapeWeightFunction& weight_fn,
                          int64_t approximate_size_bytes, int max_level,
                          S2Error* error);

  // A wrapper around InitToShapeDensity which uses the number of vertices in
  // each shape to calculate weights.
  bool InitToVertexDensity(const S2ShapeIndex& index,
                           int64_t approximate_size_bytes, int max_level,
                           S2Error* error);

  // Type definition for a function that returns a pointer to the associated T
  // for a given S2Shape. This function is allowed to return nullptr if the user
  // wishes not to create an association to a given S2Shape.
  template <typename T>
  using FeatureLookupFunction = std::function<const T*(const S2Shape&)>;

  // Type definition for a function that returns the weight of a T. Weights
  // must be non-negative integers.
  template <typename T>
  using FeatureWeightFunction = std::function<int64_t(const T&)>;

  // Similar to InitToShapeDensity but allows the caller to relate arbitrary
  // objects of type T (referred to here as features) back to shapes in the
  // index.  'feature_lookup_fn' provides the relation from a given S2Shape
  // within the 'index' to the user provided feature and 'feature_weight_fn'
  // returns the weight of a given feature.
  //
  // Each feature is guaranteed to be weighed only once per intersection with a
  // cell (e.g. in the case of a feature having multiple shapes).  However, this
  // guarantee is provided via pointer identity, and it is up to the caller to
  // ensure that for a given FeatureLookupFunction, the same input S2Shape
  // always returns a pointer to the same underlying object.
  template <typename T>
  bool InitToFeatureDensity(const S2ShapeIndex& index,
                            const FeatureLookupFunction<T>& feature_lookup_fn,
                            const FeatureWeightFunction<T>& feature_weight_fn,
                            int64_t approximate_size_bytes, int max_level,
                            S2Error* error);

  // Returns a new S2DensityTree that contains the combined weights across the
  // cells the given input trees.  The new tree will be held to the constraints
  // of 'approximate_size_bytes' and 'max_level' which may result in a lossy
  // operation if any of the input trees are more detailed than the sum
  // operation allows.
  bool InitToSumDensity(std::vector<const S2DensityTree*>& trees,
                        int64_t approximate_size_bytes, int max_level,
                        S2Error* error);

  // Same as above, but sums trees without regard for a maximum encoded size.
  // This enables the use of a substantially faster summing algorithm.
  bool InitToSumDensity(std::vector<const S2DensityTree*>& trees, int max_level,
                        S2Error* error);

  // Initialize the S2DensityTree from the given decoder.
  bool Init(Decoder* decoder, S2Error& error);

  // Write the encoded form of the S2DensityTree to the given encoder.
  void Encode(Encoder* encoder) const;

  // Returns the number of bytes this tree occupies in its encoded form.
  size_t EncodedSize() const { return encoded_.size(); }

  class Cell;

  enum class VisitAction : int8_t {
    ENTER_CELL,  // Continue visitation with the children of this node.
    SKIP_CELL,   // Continue visitation but skip the children of this node.
    STOP,        // Abort visitation.  VisitCells will return false.
  };

  // Perform a depth-first traversal of the cells within the tree.  When this
  // function stops early, either due to a decoding error or 'visitor_fn'
  // returning STOP, this function will return false, and the 'error' object
  // will describe the reason.
  using CellVisitor = std::function<VisitAction(S2CellId, const Cell&)>;
  bool VisitCells(const CellVisitor& visitor_fn, S2Error* error) const;

  // Returns a partitioning of the sphere into S2CellUnions such that the total
  // weight of the cells in each union is no more than 'max_weight'.  Cells in
  // each union tend to be geographically close to each other, following along
  // the Hilbert curve to fill each union.
  //
  // Each partition's S2CellUnion will not be normalized. This ensures that the
  // only cells that are returned are cells which are represented in the
  // S2DensityTree.
  //
  // CAVEAT: If the density tree contains single cells with no children, whose
  // weight exceeds 'max_weight', each such cell will be an element of the
  // partitioning. This can also occur as a result of the density tree lacking
  // sufficient detail (e.g. because it was truncated due to a size or cell
  // level limit).
  std::vector<S2CellUnion> GetPartitioning(int64_t max_weight,
                                           S2Error* error) const;

  // Returns a fully-decoded map of this tree.  This is only useful if far more
  // lookups will be done than there are entries in the map, and the lookups are
  // sufficiently random that a DecodedPath is not sufficient.  In that case
  // fully decoding all the cells may be faster than lazily-decoding cells on
  // each lookup.
  absl::btree_map<S2CellId, int64_t> Decode(S2Error* error) const;

  class DecodedPath;

  // Returns the normalized weight of the given cell_id, which is much more
  // accurate if the resulting weights will be added together.  If a given cell
  // is not present in the density tree, a weight of zero is returned (this is
  // not an error).  Error will be set if the requested cell_id could not be
  // decoded.
  int64_t GetNormalCellWeight(S2CellId cell_id, DecodedPath* cell_path,
                              S2Error* error) const;

  // Returns the weight of the given cell_id.  If a given cell is not present in
  // the density tree, a weight of zero is returned (this is not an error).
  // Error will be set if the requested cell_id could not be decoded.
  int64_t GetCellWeight(S2CellId cell_id, DecodedPath* cell_path,
                        S2Error* error) const;

  // Density trees map cells to weights that intersect that cell, but larger
  // features may intersect many density cells, so the deeper into the tree we
  // look, the more a feature will contribute undetected duplicates of its
  // weight many times.
  //
  // So we offer building a normalized version of the tree where every node's
  // weight is scaled by (its parent's weight / the sum of weights of the node
  // and its siblings). This makes the weight of a parent equal to the sum of
  // its children.
  S2DensityTree Normalize(absl::Nonnull<S2Error*> error) const;

  // Returns an S2CellUnion containing the leaves of this tree.  The cell union
  // is not necessarily normalized.
  S2CellUnion Leaves(absl::Nonnull<S2Error*> error) const;

  // The decoded weight and offsets of encoded cells.
  class Cell {
   public:
    Cell() = default;

    // Returns the weight of the cell.
    int64_t weight() const { return weight_; }

    // Returns true if this Cell has any children that can be visited.
    bool has_children() const { return offsets_ != kNoChildren; }

   private:
    friend class S2DensityTree;
    friend class DecodedPath;
    friend class Node;

    static constexpr std::array<int64_t, 4> kNoChildren{-1, -1, -1, -1};

    // Reset the Cell to an uninitialized state.
    void Clear();

    // Loads this cell with the weight and child offsets at the given Decoder
    // position.
    bool Decode(Decoder& decoder, S2Error* error);

    // A convenience method for decoding cells when the bookkeeping of the
    // decoder state isn't necessary.
    bool DecodeAt(const S2DensityTree* tree, uint64_t pos, S2Error* error);

    // Return the offset of the given child index.  May return -1, indicating
    // that there is no such child.  No bounds checking is done, and the index
    // must be in the range [0,3].
    int64_t child_offset(int index) const { return offsets_[index]; }

    int64_t weight_ = 0;
    std::array<int64_t, 4> offsets_ = kNoChildren;
  };

  // Represents a decoded path of cells in the tree.  This only decodes as much
  // of the path from the root down to the cell given to 'weight' as is
  // different from the last call, which greatly accelerates calls to many
  // nearby cells.
  class DecodedPath {
   public:
    explicit DecodedPath(const S2DensityTree* tree) : tree_(tree) {}

    // Returns the decoded cell for the given S2CellId.  The weight should only
    // be used if error->ok() == true.  If a given cell_id is not present in the
    // density tree, a Cell object with a weight of zero is returned (this is
    // not an error).
    const S2DensityTree::Cell* GetCell(const S2CellId cell_id, S2Error* error);

    // Returns the tree passed during construction.
    const S2DensityTree& tree() const { return *tree_; }

   private:
    // Decodes a face cell and sets it as the first cell in the stack.  Error
    // will be set if the requested face cell could not be loaded.
    void LoadFace(int face, S2Error* error);

    // Decodes all parents of the given cell_id and the cell_id itself.  This
    // method assumes that the correct face cell has already been loaded via a
    // call to LoadFace.  Error will be set if the requested cell could not be
    // loaded. If the requested cell is not present in the tree, a Cell object
    // with a weight of zero is returned.
    const S2DensityTree::Cell* LoadCell(S2CellId cell_id, S2Error* error);

    const S2DensityTree* const tree_;
    // Valid cell levels range from {0 ... kMaxLevel}.
    std::vector<Cell> stack_{S2CellId::kMaxLevel + 1};
    S2CellId last_ = S2CellId::Sentinel();
  };

  // A collector of cell/weight pairs and an encoder of them.  Since the face
  // and cell encodings both require varint offsets to the children in the
  // header, we write everything backwards, and then reverse the array at the
  // end.
  class TreeEncoder {
   public:
    TreeEncoder() = default;

    // Inserts the given cell/weight pair into the current encoder.
    void Put(S2CellId cell, int64_t weight);

    // Encodes the current set of cell/weight pairs and returns the new density
    // tree.
    void Build(S2DensityTree* tree);

    // Returns the estimated encoded size of 'weight', where we can precompute
    // the weight size, but must guess that the offset to the weight has to skip
    // an average of two siblings with that weight.
    static int EstimateSize(int64_t weight);

    void Clear();

   private:
    void EncodeTreeReversed(ReversibleBytes* output);
    void EncodeSubtreeReversed(S2CellId cell, int64_t weight,
                               ReversibleBytes* output);

    absl::btree_map<S2CellId, int64_t> weights_;
  };

  // A builder of density trees that visits a 'weight_fn' function in
  // breadth-first order, starting with the top level face cells, and visits
  // every S2 cell where the weight_fn result is greater than 0, the cell
  // level is below the given max level, and the final tree size is estimated
  // to be at least the given approximate size.
  class BreadthFirstTreeBuilder {
   public:
    BreadthFirstTreeBuilder(int64_t approx_size_bytes, int max_level,
                            TreeEncoder& encoder)
        : approximate_size_bytes_(approx_size_bytes),
          max_level_(max_level),
          encoder_(encoder) {}

    using CellWeightFunction = std::function<int64_t(S2CellId, S2Error*)>;

    // Builds the density tree that forms from a breadth-first visitation of the
    // given weight_fn.  The behavior of the builder is steered by the weight
    // returned:
    //
    //   Positive values - The weight is assigned to the cell and the builder
    //     will visit the cell's children and weigh them.
    //
    //   Negative values - The absolute value of the weight is assigned to the
    //     cell and the builder will -not- visit the cell's children.
    //
    //   Zero - The builder will skip the cell and its children.
    bool Build(const CellWeightFunction& weight_fn, S2DensityTree* tree,
               S2Error* error) const;

   private:
    const int64_t approximate_size_bytes_;
    const int max_level_;
    TreeEncoder& encoder_;
  };

 private:
  // A reusable measurer of shape index density that computes the weight at
  // each cell the index intersects, as the sum of the weight_fn for each shape
  // in the index that intersects a cell.
  class IndexCellWeightFunction {
   public:
    IndexCellWeightFunction(const S2ShapeIndex* index,
                            ShapeWeightFunction weight_fn)
        : index_region_(index), weight_fn_(std::move(weight_fn)) {}

    // Returns the weight of the given 'cell_id', calculated by summing the
    // weight of the shapes intersecting with the cell.  This function can
    // return negative weights, which indicates that the cell is fully contained
    // by the index.
    int64_t WeighCell(S2CellId cell_id, S2Error* error);

   private:
    S2ShapeIndexRegion<S2ShapeIndex> index_region_;
    ShapeWeightFunction weight_fn_;
  };

  // A reusable measurer of feature density that computes the weight at each
  // feature the index intersects, according to the weight of each feature.
  template <typename T>
  class FeatureCellWeightFunction {
   public:
    FeatureCellWeightFunction(const S2ShapeIndex* index,
                              const FeatureLookupFunction<T>& feature_lookup_fn,
                              const FeatureWeightFunction<T>& feature_weight_fn)
        : index_region_(index),
          feature_lookup_fn_(feature_lookup_fn),
          feature_weight_fn_(feature_weight_fn) {}

    // Returns the weight of the given 'cell_id', calculated by summing the
    // weight of the features intersecting with the cell.  This function can
    // return negative weights, which indicates that the cell is fully contained
    // by the index.
    int64_t WeighCell(S2CellId cell_id, S2Error* error);

   private:
    S2ShapeIndexRegion<S2ShapeIndex> index_region_;
    const FeatureLookupFunction<T>& feature_lookup_fn_;
    const FeatureWeightFunction<T>& feature_weight_fn_;
  };

  friend class S2DensityClusterQueryTest;
  friend class CoveringsTest;
  friend class SumDensityTreesTest;
  friend class TreeEncoderTest;
  friend class Node;

  // DecodedFaces stores the offsets to each of the face cells encoded in the
  // tree. The index into the array indicates the face cell. A negative offset
  // indicates that the tree does not contain this face.
  using DecodedFaces = std::array<int64_t, S2CellId::kNumFaces>;

  // Loads face positions.  The cursor must be positioned at the start of the
  // version string.
  static bool DecodeHeader(Decoder* decoder, DecodedFaces* decoded_faces,
                           S2Error* error);

  // The recursive portion of the public Visit method.
  bool VisitRecursive(const CellVisitor& visitor_fn, S2CellId cell_id,
                      int64_t position, S2Error* error) const;

  // An overload of the public GetNormalWeight that provides its own decoded
  // Cell. It's assumed that the cell corresponds to the passed cell_id.
  int64_t GetNormalCellWeight(S2CellId cell_id, const Cell& cell,
                              DecodedPath* cell_path, S2Error* error) const;

  std::string encoded_;
  DecodedFaces decoded_faces_{-1, -1, -1, -1, -1, -1};
};

template <typename T>
bool S2DensityTree::InitToFeatureDensity(
    const S2ShapeIndex& index,
    const FeatureLookupFunction<T>& feature_lookup_fn,
    const FeatureWeightFunction<T>& feature_weight_fn,
    int64_t approximate_size_bytes, int max_level, S2Error* error) {
  ABSL_DCHECK(error != nullptr) << "error must be non-nullptr";
  *error = S2Error::Ok();

  FeatureCellWeightFunction<T> index_cell_weight_fn(&index, feature_lookup_fn,
                                                    feature_weight_fn);

  TreeEncoder encoder;
  BreadthFirstTreeBuilder builder(approximate_size_bytes, max_level, encoder);
  return builder.Build(
      [&](const S2CellId cell_id, S2Error* error) {
        return index_cell_weight_fn.WeighCell(cell_id, error);
      },
      this, error);
}

template <typename T>
int64_t S2DensityTree::FeatureCellWeightFunction<T>::WeighCell(
    const S2CellId cell_id, S2Error* error) {
  int64_t sum = 0;
  bool all_contained = true;

  absl::flat_hash_set<const T*> found;
  index_region_.VisitIntersectingShapes(
      S2Cell(cell_id), [&](const S2Shape* shape, bool contains_target) {
        const T* feature = feature_lookup_fn_(*shape);

        if (feature && found.insert(feature).second) {
          const int64_t weight = feature_weight_fn_(*feature);
          ABSL_DCHECK_GE(weight, 0);
          ABSL_DCHECK_LE(weight, kMaxWeight);
          sum += weight;
          all_contained &= contains_target;
        }

        return true;
      });

  sum = std::min(sum, kMaxWeight);
  return all_contained ? -sum : sum;
}

#endif  // S2_S2DENSITY_TREE_H_
