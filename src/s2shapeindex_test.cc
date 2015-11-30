// Copyright 2012 Google Inc. All Rights Reserved.
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

#include "s2shapeindex.h"

#include <pthread.h>
#include <memory>
#include <string>
#include <vector>

#include <gflags/gflags.h>
#include "base/mutex.h"
#include <gtest/gtest.h>
#include "r2.h"
#include "r2rect.h"
#include "s1angle.h"
#include "s2.h"
#include "s2cap.h"
#include "s2cellid.h"
#include "s2cellunion.h"
#include "s2edgeutil.h"
#include "s2error.h"
#include "s2loop.h"
#include "s2polygon.h"
#include "s2shapeutil.h"
#include "s2testing.h"
#include "s2textformat.h"
#include "util/gtl/stl_util.h"

using s2shapeutil::S2EdgeVectorShape;
using std::unique_ptr;
using std::vector;


class S2ShapeIndexTest : public ::testing::Test {
 protected:
  // This test harness owns an S2ShapeIndex for convenience.  All S2Shapes in
  // the index are automatically deleted upon destruction.

  S2ShapeIndex index_;

  // Verify that that every cell of the index contains the correct edges, and
  // that no cells are missing from the index.  The running time of this
  // function is quadratic in the number of edges.
  void QuadraticValidate();

  // Given an edge and a cell id, determine whether or not the edge should be
  // present in that cell and verify that this matches "index_has_edge".
  void ValidateEdge(S2Point const& a, S2Point const& b,
                    S2CellId id, bool index_has_edge);

  // Given a shape and a cell id, determine whether or not the shape contains
  // the cell center and verify that this matches "index_contains_center".
  void ValidateInterior(S2Shape const* shape, S2CellId id,
                        bool index_contains_center);
};

void S2ShapeIndexTest::QuadraticValidate() {
  // Iterate through a sequence of nonoverlapping cell ids that cover the
  // sphere and include as a subset all the cell ids used in the index.  For
  // each cell id, verify that the expected set of edges is present.

  // "min_cellid" is the first S2CellId that has not been validated yet.
  S2CellId min_cellid = S2CellId::Begin(S2CellId::kMaxLevel);
  for (S2ShapeIndex::Iterator it(index_); ; it.Next()) {
    // Generate a list of S2CellIds ("skipped cells") that cover the gap
    // between the last cell we validated and the next cell in the index.
    S2CellUnion skipped;
    if (!it.Done()) {
      S2CellId cellid = it.id();
      EXPECT_GE(cellid, min_cellid);
      skipped.InitFromBeginEnd(min_cellid, cellid.range_min());
      min_cellid = cellid.range_max().next();
    } else {
      // Validate the empty cells beyond the last cell in the index.
      skipped.InitFromBeginEnd(min_cellid,
                               S2CellId::End(S2CellId::kMaxLevel));
    }
    // Iterate through all the shapes, simultaneously validating the current
    // index cell and all the skipped cells.
    int short_edges = 0;  // number of edges counted toward subdivision
    for (int id = 0; id < index_.num_shape_ids(); ++id) {
      S2Shape const* shape = index_.shape(id);
      S2ClippedShape const* clipped = NULL;
      if (!it.Done()) clipped = it.cell()->find_clipped(shape);

     // First check that contains_center() is set correctly.
      for (int j = 0; j < skipped.num_cells(); ++j) {
        ValidateInterior(shape, skipped.cell_id(j), false);
      }
      if (!it.Done()) {
        bool contains_center = clipped && clipped->contains_center();
        ValidateInterior(shape, it.id(), contains_center);
      }
      // Then check that the appropriate edges are present.
      for (int e = 0; e < shape->num_edges(); ++e) {
        S2Point const *a, *b;
        shape->GetEdge(e, &a, &b);
        for (int j = 0; j < skipped.num_cells(); ++j) {
          ValidateEdge(*a, *b, skipped.cell_id(j), false);
        }
        if (!it.Done()) {
          bool has_edge = clipped && clipped->ContainsEdge(e);
          ValidateEdge(*a, *b, it.id(), has_edge);
          if (has_edge && it.id().level() < index_.GetEdgeMaxLevel(*a, *b)) {
            ++short_edges;
          }
        }
      }
    }
    EXPECT_LE(short_edges, index_.options().max_edges_per_cell());
    if (it.Done()) break;
  }
}

// Verify that "index_has_edge" is true if and only if the edge AB intersects
// the given cell id.
void S2ShapeIndexTest::ValidateEdge(S2Point const& a, S2Point const& b,
                                    S2CellId id, bool index_has_edge) {
  // Expand or shrink the padding slightly to account for errors in the
  // function we use to test for intersection (IntersectsRect).
  double padding = S2ShapeIndex::kCellPadding;
  padding += (index_has_edge ? 1 : -1) * S2EdgeUtil::kIntersectsRectErrorUVDist;
  R2Rect bound = id.GetBoundUV().Expanded(padding);
  R2Point a_uv, b_uv;
  EXPECT_EQ(S2EdgeUtil::ClipToPaddedFace(a, b, id.face(), padding, &a_uv, &b_uv)
            && S2EdgeUtil::IntersectsRect(a_uv, b_uv, bound),
            index_has_edge);
}

void S2ShapeIndexTest::ValidateInterior(S2Shape const* shape, S2CellId id,
                                        bool index_contains_center) {
  if (!shape->has_interior()) return;
  S2Point a = S2::Origin(), b = id.ToPoint();
  S2EdgeUtil::EdgeCrosser crosser(&a, &b);
  bool contains_center = shape->contains_origin();
  for (int e = 0; e < shape->num_edges(); ++e) {
    S2Point const *c, *d;
    shape->GetEdge(e, &c, &d);
    contains_center ^= crosser.EdgeOrVertexCrossing(c, d);
  }
  EXPECT_EQ(contains_center, index_contains_center);
}

namespace {

void TestIteratorMethods(S2ShapeIndex const& index) {
  S2ShapeIndex::Iterator it(index);
  EXPECT_TRUE(it.AtBegin());
  it.Finish();
  EXPECT_TRUE(it.Done());
  vector<S2CellId> ids;
  S2CellId min_cellid = S2CellId::Begin(S2CellId::kMaxLevel);
  for (it.Reset(); !it.Done(); it.Next()) {
    S2CellId cellid = it.id();
    S2CellUnion skipped;
    skipped.InitFromBeginEnd(min_cellid, cellid.range_min());
    S2ShapeIndex::Iterator it2(index);
    for (int i = 0; i < skipped.num_cells(); ++i) {
      EXPECT_FALSE(it2.Locate(skipped.cell_id(i).ToPoint()));
      EXPECT_EQ(S2ShapeIndex::DISJOINT, it2.Locate(skipped.cell_id(i)));
      it2.Reset();
      it2.Seek(skipped.cell_id(i));
      EXPECT_EQ(cellid, it2.id());
    }
    if (!ids.empty()) {
      EXPECT_FALSE(it.AtBegin());
      it2 = it;
      it2.Prev();
      EXPECT_EQ(ids.back(), it2.id());
      it2.Next();
      EXPECT_EQ(cellid, it2.id());
      it2.Seek(ids.back());
      EXPECT_EQ(ids.back(), it2.id());
      it2.SeekForward(cellid);
      EXPECT_EQ(cellid, it2.id());
      it2.SeekForward(ids.back());
      EXPECT_EQ(cellid, it2.id());
    }
    it2.Reset();
    EXPECT_EQ(cellid.ToPoint(), it.center());
    EXPECT_TRUE(it2.Locate(it.center()));
    EXPECT_EQ(cellid, it2.id());
    it2.Reset();
    EXPECT_EQ(S2ShapeIndex::INDEXED, it2.Locate(cellid));
    EXPECT_EQ(cellid, it2.id());
    if (!cellid.is_face()) {
      it2.Reset();
      EXPECT_EQ(S2ShapeIndex::SUBDIVIDED, it2.Locate(cellid.parent()));
      EXPECT_LE(it2.id(), cellid);
      EXPECT_GE(it2.id(), cellid.parent().range_min());
    }
    if (!cellid.is_leaf()) {
      for (int i = 0; i < 4; ++i) {
        it2.Reset();
        EXPECT_EQ(S2ShapeIndex::INDEXED, it2.Locate(cellid.child(i)));
        EXPECT_EQ(cellid, it2.id());
      }
    }
    ids.push_back(cellid);
    min_cellid = cellid.range_max().next();
  }
}

TEST_F(S2ShapeIndexTest, NoEdges) {
  EXPECT_TRUE(S2ShapeIndex::Iterator(index_).Done());
  TestIteratorMethods(index_);
}

TEST_F(S2ShapeIndexTest, OneEdge) {
  index_.Add(new S2EdgeVectorShape(S2Point(1, 0, 0), S2Point(0, 1, 0)));
  QuadraticValidate();
  TestIteratorMethods(index_);
}

TEST_F(S2ShapeIndexTest, LoopsSpanningThreeFaces) {
  S2Polygon polygon;
  int const kNumEdges = 100;  // Validation is quadratic
  // Construct two loops consisting of kNumEdges vertices each, centered
  // around the cube vertex at the start of the Hilbert curve.
  S2Testing::ConcentricLoopsPolygon(S2Point(1, -1, -1).Normalize(), 2,
                                    kNumEdges, &polygon);
  vector<S2Loop*> loops;
  polygon.Release(&loops);
  for (int i = 0; i < loops.size(); ++i) {
    index_.Add(new S2Loop::Shape(loops[i]));
  }
  QuadraticValidate();
  TestIteratorMethods(index_);
  STLDeleteElements(&loops);
}


TEST_F(S2ShapeIndexTest, ManyIdenticalEdges) {
  int const kNumEdges = 100;  // Validation is quadratic
  S2Point a = S2Point(0.99, 0.99, 1).Normalize();
  S2Point b = S2Point(-0.99, -0.99, 1).Normalize();
  for (int i = 0; i < kNumEdges; ++i) {
    index_.Add(new S2EdgeVectorShape(a, b));
  }
  QuadraticValidate();
  TestIteratorMethods(index_);
  // Since all edges span the diagonal of a face, no subdivision should
  // have occurred (with the default index options).
  for (S2ShapeIndex::Iterator it(index_); !it.Done(); it.Next()) {
    EXPECT_EQ(0, it.id().level());
  }
}

TEST_F(S2ShapeIndexTest, ManyTinyEdges) {
  // This test adds many edges to a single leaf cell, to check that
  // subdivision stops when no further subdivision is possible.
  int const kNumEdges = 100;  // Validation is quadratic
  // Construct two points in the same leaf cell.
  S2Point a = S2CellId::FromPoint(S2Point(1, 0, 0)).ToPoint();
  S2Point b = (a + S2Point(0, 1e-12, 0)).Normalize();
  S2EdgeVectorShape* shape = new S2EdgeVectorShape;
  for (int i = 0; i < kNumEdges; ++i) {
    shape->Add(a, b);
  }
  index_.Add(shape);
  QuadraticValidate();
  // Check that there is exactly one index cell and that it is a leaf cell.
  S2ShapeIndex::Iterator it(index_);
  ASSERT_TRUE(!it.Done());
  EXPECT_TRUE(it.id().is_leaf());
  it.Next();
  EXPECT_TRUE(it.Done());
}


// Add the loops to the given index.
void AddLoops(vector<S2Loop*> const& loops, S2ShapeIndex* index) {
  for (int i = 0; i < loops.size(); ++i) {
    index->Add(new S2Loop::Shape(loops[i]));
  }
}

// Return true if any loop crosses any other loop (including vertex crossings
// and duplicate edges), or any loop has a self-intersection (including
// duplicate vertices).
static bool HasAnyCrossing(S2ShapeIndex const& index,
                           vector<S2Loop*> const& loops) {
  S2Error error;
  if (s2shapeutil::FindAnyCrossing(index, loops, &error)) {
    VLOG(1) << error;
    return true;
  }
  return false;
}

// This function recursively verifies that HasCrossing returns the given
// result for all possible cyclic permutations of the loop vertices for the
// given set of loops.
void TestHasCrossingPermutations(vector<S2Loop*>* loops, int i,
                                 bool has_crossing) {
  if (i == loops->size()) {
    S2ShapeIndex index;
    AddLoops(*loops, &index);
    EXPECT_EQ(has_crossing, HasAnyCrossing(index, *loops));
  } else {
    S2Loop* orig_loop = (*loops)[i];
    for (int j = 0; j < orig_loop->num_vertices(); ++j) {
      vector<S2Point> vertices;
      for (int k = 0; k < orig_loop->num_vertices(); ++k) {
        vertices.push_back(orig_loop->vertex(j + k));
      }
      unique_ptr<S2Loop> new_loop(new S2Loop(vertices));
      (*loops)[i] = new_loop.get();
      TestHasCrossingPermutations(loops, i+1, has_crossing);
    }
    (*loops)[i] = orig_loop;
  }
}

// Given a string reprsenting a polygon, and a boolean indicating whether this
// polygon has any self-intersections or loop crossings, verify that all
// HasAnyCrossing returns the expected result for all possible cyclic
// permutations of the loop vertices.
void TestHasCrossing(const string& polygon_str, bool has_crossing) {
  FLAGS_s2debug = false;  // Allow invalid polygons (restored by gUnit)
  unique_ptr<S2Polygon> polygon(s2textformat::MakePolygon(polygon_str));
  vector<S2Loop*> loops;
  polygon->Release(&loops);
  TestHasCrossingPermutations(&loops, 0, has_crossing);
  STLDeleteElements(&loops);
}

TEST_F(S2ShapeIndexTest, HasCrossing) {
  // Coordinates are (lat,lng), which can be visualized as (y,x).
  TestHasCrossing("0:0, 0:1, 0:2, 1:2, 1:1, 1:0", false);
  TestHasCrossing("0:0, 0:1, 0:2, 1:2, 0:1, 1:0", true);  // duplicate vertex
  TestHasCrossing("0:0, 0:1, 1:0, 1:1", true);  // edge crossing
  TestHasCrossing("0:0, 1:1, 0:1; 0:0, 1:1, 1:0", true);  // duplicate edge
  TestHasCrossing("0:0, 1:1, 0:1; 1:1, 0:0, 1:0", true);  // reversed edge
  TestHasCrossing("0:0, 0:2, 2:2, 2:0; 1:1, 0:2, 3:1, 2:0",
                  true);  // vertex crossing
}

// A test that repeatedly updates "index_" in one thread and attempts to
// concurrently read the index_ from several other threads.  When all threads
// have finished reading, the first thread makes another update.
//
// Note that we only test concurrent read access, since S2ShapeIndex requires
// all updates to be single-threaded and not concurrent with any reads.
class LazyUpdatesTest : public ::testing::Test {
 public:
  LazyUpdatesTest() : num_updates_(0), num_readers_left_(0) {
  }

  // The function executed by each reader thread.
  void ReaderThread();

 protected:
  S2ShapeIndex index_;
  // The following fields are guarded by lock_.
  Mutex lock_;
  int num_updates_;
  int num_readers_left_;

  // Signalled when a new update is ready to be processed.
  CondVar update_ready_;
  // Signalled when all readers have processed the latest update.
  CondVar all_readers_done_;
};

void LazyUpdatesTest::ReaderThread() {
  lock_.Lock();
  for (int last_update = 0; ; last_update = num_updates_) {
    while (num_updates_ == last_update) {
      update_ready_.Wait(&lock_);
    }
    if (num_updates_ < 0) break;

    // The index is built on demand the first time we attempt to use it.
    // We intentionally release the lock so that many threads have a chance
    // to access the S2ShapeIndex in parallel.
    lock_.Unlock();
    for (S2ShapeIndex::Iterator it(index_); !it.Done(); it.Next())
      continue;
    lock_.Lock();
    if (--num_readers_left_ == 0) {
      all_readers_done_.Signal();
    }
  }
  lock_.Unlock();
}

static void* StartReader(void* arg) {
  static_cast<LazyUpdatesTest*>(arg)->ReaderThread();
  return NULL;
}

TEST_F(LazyUpdatesTest, ConstMethodsThreadSafe) {
  // Ensure that lazy updates are thread-safe.  In other words, make sure that
  // nothing bad happens when multiple threads call "const" methods that
  // cause pending updates to be applied.

  // The number of readers should be large enough so that it is likely that
  // several readers will be running at once (with a multiple-core CPU).
  int const kNumReaders = 8;
  pthread_t readers[kNumReaders];
  for (int i = 0; i < kNumReaders; ++i) {
    CHECK_EQ(0, pthread_create(&readers[i], NULL, StartReader,
                               static_cast<void*>(this)));
  }
  lock_.Lock();
  int const kIters = 100;
  for (int iter = 0; iter < kIters; ++iter) {
    // Loop invariant: lock_ is held and num_readers_left_ == 0.
    DCHECK_EQ(0, num_readers_left_);
    // Since there are no readers, it is safe to modify the index.
    index_.Reset();
    int num_vertices = 4 * S2Testing::rnd.Skewed(10);  // Up to 4K vertices
    unique_ptr<S2Loop> loop(S2Testing::MakeRegularLoop(
        S2Testing::RandomPoint(), S2Testing::KmToAngle(5), num_vertices));
    index_.Add(new S2Loop::Shape(loop.get()));
    num_readers_left_ = kNumReaders;
    ++num_updates_;
    update_ready_.SignalAll();
    while (num_readers_left_ > 0) {
      all_readers_done_.Wait(&lock_);
    }
  }
  // Signal the readers to exit.
  num_updates_ = -1;
  update_ready_.SignalAll();
  lock_.Unlock();
  for (int i = 0; i < kNumReaders; ++i) {
    CHECK_EQ(0, pthread_join(readers[i], NULL));
  }
}

TEST(S2ShapeIndex, GetContainingShapes) {
  // Also tests Contains(S2Shape const*, S2Point).
  int const kNumVerticesPerLoop = 10;
  S1Angle const kMaxLoopRadius = S2Testing::KmToAngle(10);
  S2Cap const center_cap(S2Testing::RandomPoint(), kMaxLoopRadius);
  using LoopShape = s2shapeutil::S2LoopOwningShape;
  S2ShapeIndex index;
  for (int i = 0; i < 100; ++i) {
    S2Loop* loop = S2Testing::MakeRegularLoop(
        S2Testing::SamplePoint(center_cap),
        S2Testing::rnd.RandDouble() * kMaxLoopRadius,
        kNumVerticesPerLoop);
    index.Add(new LoopShape(loop));
  }
  for (int i = 0; i < 100; ++i) {
    S2Point p = S2Testing::SamplePoint(center_cap);
    vector<S2Shape const*> expected, actual;
    for (int j = 0; j < index.num_shape_ids(); ++j) {
      S2Shape const* shape = index.shape(j);
      S2Loop const* loop = static_cast<LoopShape const*>(shape)->loop();
      if (loop->Contains(p)) {
        EXPECT_TRUE(index.ShapeContains(shape, p));
        expected.push_back(shape);
      } else {
        EXPECT_FALSE(index.ShapeContains(shape, p));
      }
    }
    index.GetContainingShapes(p, &actual);
    EXPECT_EQ(expected, actual);
  }
}


}  // namespace
