// Copyright 2018 Google Inc. All Rights Reserved.
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

#ifndef S2_ENCODED_S2SHAPE_INDEX_H_
#define S2_ENCODED_S2SHAPE_INDEX_H_

#include "s2/encoded_s2cell_id_vector.h"
#include "s2/encoded_string_vector.h"
#include "s2/mutable_s2shape_index.h"

class EncodedS2ShapeIndex final : public S2ShapeIndex {
 public:
  using Options = MutableS2ShapeIndex::Options;
  using ShapeFactory = S2ShapeIndex::ShapeFactory;

  // Creates an index that must be initialized by calling Init().
  EncodedS2ShapeIndex();

  ~EncodedS2ShapeIndex() override;

  // Initializes the EncodedS2ShapeIndex, returning true on success.
  //
  // This method does not decode the S2Shape objects in the index; this is
  // the responsibility of the client-provided function "shape_factory"
  // (see s2shapeutil_coding.h).  Example usage:
  //
  //   index.Init(decoder, s2shapeutil::LazyDecodeShapeFactory(decoder));
  //
  // Note that the encoded shape vector must *precede* the encoded S2ShapeIndex
  // in the Decoder's data buffer in this example.
  bool Init(Decoder* decoder, const ShapeFactory& shape_factory);

  const Options& options() const { return options_; }

  // The number of distinct shape ids in the index.  This equals the number of
  // shapes in the index provided that no shapes have ever been removed.
  // (Shape ids are not reused.)
  int num_shape_ids() const override { return shapes_.size(); }

  // Return a pointer to the shape with the given id, or nullptr if the shape
  // has been removed from the index.
  S2Shape* shape(int id) const override;

  // Minimizes memory usage by requesting that any data structures that can be
  // rebuilt should be discarded.  This method invalidates all iterators.
  //
  // Like all non-const methods, this method is not thread-safe.
  void Minimize() override;

  class Iterator final : public IteratorBase {
   public:
    // Default constructor; must be followed by a call to Init().
    Iterator();

    // Constructs an iterator positioned as specified.  By default iterators
    // are unpositioned, since this avoids an extra seek in this situation
    // where one of the seek methods (such as Locate) is immediately called.
    //
    // If you want to position the iterator at the beginning, e.g. in order to
    // loop through the entire index, do this instead:
    //
    //   for (EncodedS2ShapeIndex::Iterator it(&index, S2ShapeIndex::BEGIN);
    //        !it.done(); it.Next()) { ... }
    explicit Iterator(const EncodedS2ShapeIndex* index,
                      InitialPosition pos = UNPOSITIONED);

    // Initializes an iterator for the given EncodedS2ShapeIndex.
    void Init(const EncodedS2ShapeIndex* index,
              InitialPosition pos = UNPOSITIONED);

    // Inherited non-virtual methods:
    //   S2CellId id() const;
    //   const S2ShapeIndexCell& cell() const;
    //   bool done() const;
    //   S2Point center() const;

    // IteratorBase API:
    void Begin() override;
    void Finish() override;
    void Next() override;
    bool Prev() override;
    void Seek(S2CellId target) override;
    bool Locate(const S2Point& target) override;
    CellRelation Locate(S2CellId target) override;

   protected:
    const S2ShapeIndexCell* GetCell() const override;
    std::unique_ptr<IteratorBase> Clone() const override;
    void Copy(const IteratorBase& other) override;

   private:
    void Refresh();  // Updates the IteratorBase fields.
    const EncodedS2ShapeIndex* index_;
    int32 cell_pos_;  // Current position in the vector of index cells.
    int32 num_cells_;
  };

  // Returns the number of bytes currently occupied by the index (including any
  // unused space at the end of vectors, etc). It has the same thread safety
  // as the other "const" methods (see introduction).
  size_t SpaceUsed() const override;

 protected:
  std::unique_ptr<IteratorBase> NewIterator(InitialPosition pos) const override;

 private:
  friend class Iterator;

  // Returns a value indicating that a shape has not been decoded yet.
  inline static S2Shape* kUndecodedShape() {
    return reinterpret_cast<S2Shape*>(1);
  }

  // Like std::atomic<S2Shape*>, but defaults to kUndecodedShape().
  class AtomicShape : public std::atomic<S2Shape*> {
   public:
    AtomicShape() : std::atomic<S2Shape*>(kUndecodedShape()) {}
  };

  S2Shape* GetShape(int id) const;
  const S2ShapeIndexCell* GetCell(int i) const;
  bool cell_decoded(int i) const;
  bool test_and_set_cell_decoded(int i) const;
  int max_cell_cache_size() const;

  std::unique_ptr<ShapeFactory> shape_factory_;

  // The options specified for this index.
  Options options_;

  // A vector containing all shapes in the index.  Initially all shapes are
  // set to kUndecodedShape(); as shapes are decoded, they are added to the
  // vector using std::atomic::compare_exchange_strong.
  mutable std::vector<AtomicShape> shapes_;

  // A vector containing the S2CellIds of each cell in the index.
  s2coding::EncodedS2CellIdVector cell_ids_;

  // A vector containing the encoded contents of each cell in the index.
  s2coding::EncodedStringVector encoded_cells_;

  // A raw array containing the decoded contents of each cell in the index.
  // Initially all values are *uninitialized memory*.  The cells_decoded_
  // field below keeps track of which elements are present.
  mutable std::unique_ptr<std::atomic<S2ShapeIndexCell*>[]> cells_;

  // A bit vector indicating which elements of cells_ have been decoded.
  // All other elements of cells_ contain uninitialized (random) memory.
  mutable std::vector<std::atomic<uint64>> cells_decoded_;

  // In order to minimize destructor time when very few cells of a large
  // S2ShapeIndex are needed, we keep track of the indices of the first few
  // cells to be decoded.  This lets us avoid scanning the cells_decoded_
  // vector when the number of cells decoded is very small.
  mutable std::vector<int> cell_cache_;

  // Protects all updates to cells_ and cells_decoded_.
  mutable SpinLock cells_lock_;

  EncodedS2ShapeIndex(const EncodedS2ShapeIndex&) = delete;
  void operator=(const EncodedS2ShapeIndex&) = delete;
};


//////////////////   Implementation details follow   ////////////////////


inline EncodedS2ShapeIndex::Iterator::Iterator() : index_(nullptr) {
}

inline EncodedS2ShapeIndex::Iterator::Iterator(
    const EncodedS2ShapeIndex* index, InitialPosition pos) {
  Init(index, pos);
}

inline void EncodedS2ShapeIndex::Iterator::Init(
    const EncodedS2ShapeIndex* index, InitialPosition pos) {
  index_ = index;
  num_cells_ = index->cell_ids_.size();
  cell_pos_ = (pos == BEGIN) ? 0 : num_cells_;
  Refresh();
}

inline void EncodedS2ShapeIndex::Iterator::Refresh() {
  if (cell_pos_ == num_cells_) {
    set_finished();
  } else {
    // It's faster to initialize the cell to nullptr even if it has already
    // been decoded, since algorithms frequently don't need it (i.e., based on
    // the S2CellId they might not need to look at the cell contents).
    set_state(index_->cell_ids_[cell_pos_], nullptr);
  }
}

inline void EncodedS2ShapeIndex::Iterator::Begin() {
  cell_pos_ = 0;
  Refresh();
}

inline void EncodedS2ShapeIndex::Iterator::Finish() {
  cell_pos_ = num_cells_;
  Refresh();
}

inline void EncodedS2ShapeIndex::Iterator::Next() {
  S2_DCHECK(!done());
  ++cell_pos_;
  Refresh();
}

inline bool EncodedS2ShapeIndex::Iterator::Prev() {
  if (cell_pos_ == 0) return false;
  --cell_pos_;
  Refresh();
  return true;
}

inline void EncodedS2ShapeIndex::Iterator::Seek(S2CellId target) {
  cell_pos_ = index_->cell_ids_.lower_bound(target);
  Refresh();
}

inline std::unique_ptr<EncodedS2ShapeIndex::IteratorBase>
EncodedS2ShapeIndex::NewIterator(InitialPosition pos) const {
  return absl::make_unique<Iterator>(this, pos);
}

inline S2Shape* EncodedS2ShapeIndex::shape(int id) const {
  S2Shape* shape = shapes_[id].load(std::memory_order_relaxed);
  if (shape != kUndecodedShape()) return shape;
  return GetShape(id);
}

// Returns true if the given cell has been decoded yet.
inline bool EncodedS2ShapeIndex::cell_decoded(int i) const {
  uint64 group_bits = cells_decoded_[i >> 6].load(std::memory_order_relaxed);
  return (group_bits & (1ULL << (i & 63))) != 0;
}

// Marks the given cell as decoded and returns true if it was already marked.
inline bool EncodedS2ShapeIndex::test_and_set_cell_decoded(int i) const {
  std::atomic<uint64>* group = &cells_decoded_[i >> 6];
  uint64 group_bits = group->load(std::memory_order_relaxed);
  uint64 test_bit = 1ULL << (i & 63);
  group->store(group_bits | test_bit, std::memory_order_relaxed);
  return (group_bits & test_bit) != 0;
}

inline int EncodedS2ShapeIndex::max_cell_cache_size() const {
  // The cell cache is sized so that scanning decoded_cells_ in the destructor
  // costs about 30 cycles per decoded cell in the worst case.  (This overhead
  // is acceptable relative to the other costs of decoding each cell.)
  //
  // For example, if there are 65,536 cells then we won't need to scan
  // encoded_cells_ unless we decode at least (65536/2048) == 32 cells.  It
  // takes about 1 cycle per 64 cells to scan encoded_cells_, so that works
  // out to (65536/64) == 1024 cycles.  However this cost is amortized over
  // the 32 cells decoded, which works out to 32 cycles per cell.
  return cell_ids_.size() >> 11;
}

#endif  // S2_ENCODED_S2SHAPE_INDEX_H_
