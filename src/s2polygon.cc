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

#include "s2polygon.h"

#include <math.h>
#include <stddef.h>
#include <algorithm>
#include <ext/hash_map>
using __gnu_cxx::hash;
using __gnu_cxx::hash_map;
#include <set>
#include <utility>
#include <vector>

#include "base/atomicops.h"
#include "base/casts.h"
#include <gflags/gflags.h>
#include <glog/logging.h>
#include "util/coding/coder.h"
#include "s1angle.h"
#include "s1interval.h"
#include "s2.h"
#include "s2cap.h"
#include "s2cell.h"
#include "s2cellid.h"
#include "s2cellunion.h"
#include "s2closestedgequery.h"
#include "s2edgequery.h"
#include "s2edgeutil.h"
#include "s2error.h"
#include "s2latlng.h"
#include "s2latlngrect.h"
#include "s2loop.h"
#include "s2pointcompression.h"
#include "s2polygonbuilder.h"
#include "s2polyline.h"
#include "s2shapeindex.h"
#include "s2shapeutil.h"
#include "util/gtl/fixedarray.h"
#include "util/gtl/stl_util.h"

using std::max;
using std::min;
using std::pair;
using std::set;
using std::vector;

DEFINE_bool(
    s2polygon_lazy_indexing, true,
    "Build the S2ShapeIndex only when it is first needed.  This can save "
    "significant amounts of memory and time when geometry is constructed but "
    "never queried, for example when converting from one format to another.");

// When adding a new encoding, be aware that old binaries will not
// be able to decode it.
static const unsigned char kCurrentLosslessEncodingVersionNumber = 1;
static const unsigned char kCurrentCompressedEncodingVersionNumber = 4;

S2Polygon::S2Polygon()
    : has_holes_(false),
      s2debug_override_(ALLOW_S2DEBUG),
      num_vertices_(0),
      unindexed_contains_calls_(0) {
}

S2Polygon::S2Polygon(vector<S2Loop*>* loops)
    : s2debug_override_(ALLOW_S2DEBUG) {
  InitNested(loops);
}

S2Polygon::S2Polygon(vector<S2Loop*>* loops, S2debugOverride override)
    : s2debug_override_(override) {
  InitNested(loops);
}

S2Polygon::S2Polygon(S2Cell const& cell)
    : loops_(1, new S2Loop(cell)),
      s2debug_override_(ALLOW_S2DEBUG),
      unindexed_contains_calls_(0) {
  InitOneLoop();
}

void S2Polygon::set_s2debug_override(S2debugOverride override) {
  s2debug_override_ = override;
}

S2debugOverride S2Polygon::s2debug_override() const {
  return static_cast<S2debugOverride>(s2debug_override_);
}

void S2Polygon::Copy(S2Polygon const* src) {
  DCHECK_EQ(0, num_loops());
  for (int i = 0; i < src->num_loops(); ++i) {
    loops_.push_back(src->loop(i)->Clone());
  }
  has_holes_ = src->has_holes_;
  s2debug_override_ = src->s2debug_override_;
  num_vertices_ = src->num_vertices();
  base::subtle::NoBarrier_Store(&unindexed_contains_calls_, 0);
  bound_ = src->bound_;
  subregion_bound_ = src->subregion_bound_;
  InitIndex();  // TODO(ericv): Copy the index efficiently.
}

S2Polygon* S2Polygon::Clone() const {
  S2Polygon* result = new S2Polygon;
  result->Copy(this);
  return result;
}

void S2Polygon::Release(vector<S2Loop*>* loops) {
  if (loops != NULL) {
    loops->insert(loops->end(), loops_.begin(), loops_.end());
  }
  // Reset the polygon to be empty.
  loops_.clear();
  has_holes_ = false;
  num_vertices_ = 0;
  bound_ = S2LatLngRect::Empty();
  subregion_bound_ = S2LatLngRect::Empty();
  index_.Reset();
}

void S2Polygon::ClearLoops() {
  index_.Reset();
  for (int i = 0; i < loops_.size(); ++i) {
    delete loops_[i];
  }
  loops_.clear();
  base::subtle::NoBarrier_Store(&unindexed_contains_calls_, 0);
}

S2Polygon::~S2Polygon() {
  ClearLoops();
}

bool S2Polygon::IsValid() const {
  S2Error error;
  if (FindValidationError(&error)) {
    LOG_IF(ERROR, FLAGS_s2debug) << error;
    return false;
  }
  return true;
}

bool S2Polygon::FindValidationError(S2Error* error) const {
  for (int i = 0; i < num_loops(); ++i) {
    // Check for loop errors that don't require building an S2ShapeIndex.
    if (loop(i)->FindValidationErrorNoIndex(error)) {
      error->Init(error->code(),
                  "Loop %d: %s", i, error->text().c_str());
      return true;
    }
    // Check that no loop is empty, and that the full loop only appears in the
    // full polygon.
    if (loop(i)->is_empty()) {
      error->Init(S2Error::POLYGON_EMPTY_LOOP,
                  "Loop %d: empty loops are not allowed", i);
      return true;
    }
    if (loop(i)->is_full() && num_loops() > 1) {
      error->Init(S2Error::POLYGON_EXCESS_FULL_LOOP,
                  "Loop %d: full loop appears in non-full polygon", i);
      return true;
    }
  }
  // Finally, check for loop self-intersections and loop pairs that cross
  // (including duplicate edges and vertices).
  //
  // TODO(ericv): Also verify the nesting hierarchy.
  return s2shapeutil::FindAnyCrossing(index_, loops_, error);
}

void S2Polygon::InsertLoop(S2Loop* new_loop, S2Loop* parent,
                           LoopMap* loop_map) {
  vector<S2Loop*>* children = &(*loop_map)[parent];
  for (int i = 0; i < children->size(); ++i) {
    S2Loop* child = (*children)[i];
    if (child->ContainsNested(new_loop)) {
      InsertLoop(new_loop, child, loop_map);
      return;
    }
  }
  // No loop may contain the complement of another loop.  (Handling this case
  // is significantly more complicated.)
  if (FLAGS_s2debug && s2debug_override_ == ALLOW_S2DEBUG) {
    CHECK(parent == NULL || !new_loop->ContainsNested(parent));
  }

  // Some of the children of the parent loop may now be children of
  // the new loop.
  vector<S2Loop*>* new_children = &(*loop_map)[new_loop];
  for (int i = 0; i < children->size();) {
    S2Loop* child = (*children)[i];
    if (new_loop->ContainsNested(child)) {
      new_children->push_back(child);
      children->erase(children->begin() + i);
    } else {
      ++i;
    }
  }
  children->push_back(new_loop);
}

void S2Polygon::InitLoop(S2Loop* loop, int depth, LoopMap* loop_map) {
  if (loop) {
    loop->set_depth(depth);
    loops_.push_back(loop);
  }
  vector<S2Loop*> const& children = (*loop_map)[loop];
  for (int i = 0; i < children.size(); ++i) {
    InitLoop(children[i], depth + 1, loop_map);
  }
}

bool S2Polygon::ContainsChild(S2Loop* a, S2Loop* b, LoopMap const& loop_map) {
  // This function is just used to verify that the loop map was
  // constructed correctly.

  if (a == b) return true;
  vector<S2Loop*> const& children = loop_map.find(a)->second;
  for (int i = 0; i < children.size(); ++i) {
    if (ContainsChild(children[i], b, loop_map)) return true;
  }
  return false;
}

void S2Polygon::InitIndex() {
  DCHECK_EQ(0, index_.num_shape_ids());
  for (int i = 0; i < num_loops(); ++i) {
    index_.Insert(new S2Loop::Shape(loop(i)));
  }
  if (!FLAGS_s2polygon_lazy_indexing) {
    index_.ForceApplyUpdates();  // Force index construction now.
  }
  if (FLAGS_s2debug && s2debug_override_ == ALLOW_S2DEBUG) {
    // Note that FLAGS_s2debug is false in optimized builds (by default).
    CHECK(IsValid());
  }
}

void S2Polygon::InitNested(vector<S2Loop*>* loops) {
  ClearLoops();
  loops_.swap(*loops);

  if (num_loops() == 1) {
    InitOneLoop();
    return;
  }

  LoopMap loop_map;
  for (int i = 0; i < num_loops(); ++i) {
    InsertLoop(loop(i), NULL, &loop_map);
  }
  // Reorder the loops in depth-first traversal order.
  loops_.clear();
  InitLoop(NULL, -1, &loop_map);

  if (FLAGS_s2debug && s2debug_override_ == ALLOW_S2DEBUG) {
    // Check that the LoopMap is correct (this is fairly cheap).
    for (int i = 0; i < num_loops(); ++i) {
      for (int j = 0; j < num_loops(); ++j) {
        if (i == j) continue;
        CHECK_EQ(ContainsChild(loop(i), loop(j), loop_map),
                 loop(i)->ContainsNested(loop(j)));
      }
    }
  }

  // Compute has_holes_, num_vertices_, bound_, subregion_bound_.
  InitLoopProperties();
}

void S2Polygon::InitOneLoop() {
  DCHECK_EQ(1, num_loops());
  S2Loop* loop = loops_[0];
  loop->set_depth(0);
  has_holes_ = false;
  num_vertices_ = loop->num_vertices();
  bound_ = loop->GetRectBound();
  subregion_bound_ = S2EdgeUtil::RectBounder::ExpandForSubregions(bound_);
  InitIndex();
}

void S2Polygon::InitOriented(vector<S2Loop*>* loops) {
  // Here is the algorithm:
  //
  // 1. Remember which of the given loops contain S2::Origin().
  //
  // 2. Invert loops as necessary to ensure that they are nestable (i.e., no
  //    loop contains the complement of any other loop).  This may result in a
  //    set of loops corresponding to the complement of the given polygon, but
  //    we will fix that problem later.
  //
  //    We make the loops nestable by first normalizing all the loops (i.e.,
  //    inverting any loops whose turning angle is negative).  This handles
  //    all loops except those whose turning angle is very close to zero
  //    (within the maximum error tolerance).  Any such loops are inverted if
  //    and only if they contain S2::Origin().  (In theory this step is only
  //    necessary if there are at least two such loops.)  The resulting set of
  //    loops is guaranteed to be nestable.
  //
  // 3. Build the polygon.  This yields either the desired polygon or its
  //    complement.
  //
  // 4. If there is at least one loop, we find a loop L that is adjacent to
  //    S2::Origin() (where "adjacent" means that there exists a path
  //    connecting S2::Origin() to some vertex of L such that the path does
  //    not cross anythat does not cross any loop).  There may be a single
  //    such adjacent loop, or there may be several (in which case they should
  //    all have the same contains_origin() value).  We choose L to be the
  //    loop containing the origin whose depth is greatest, or loop(0) (a
  //    top-level shell) if no such loop exists.
  //
  // 5. If (L originally contained origin) != (polygon contains origin), we
  //    invert the polygon.  This is done by inverting a top-level shell whose
  //    turning angle is minimal and then fixing the nesting hierarchy.  Note
  //    that because we normalized all the loops initially, this step is only
  //    necessary if the polygon requires at least one non-normalized loop to
  //    represent it.

  set<S2Loop const*> contained_origin;
  for (int i = 0; i < loops->size(); ++i) {
    S2Loop* loop = (*loops)[i];
    if (loop->contains_origin()) {
      contained_origin.insert(loop);
    }
    double angle = loop->GetTurningAngle();
    if (fabs(angle) > loop->GetTurningAngleMaxError()) {
      // Normalize the loop.
      if (angle < 0) loop->Invert();
    } else {
      // Ensure that the loop does not contain the origin.
      if (loop->contains_origin()) loop->Invert();
    }
  }
  InitNested(loops);
  if (num_loops() > 0) {
    S2Loop* origin_loop = loop(0);
    bool polygon_contains_origin = false;
    for (int i = 0; i < num_loops(); ++i) {
      if (loop(i)->contains_origin()) {
        polygon_contains_origin ^= true;
        origin_loop = loop(i);
      }
    }
    if (contained_origin.count(origin_loop) != polygon_contains_origin) {
      Invert();
    }
  }
  if (FLAGS_s2debug && s2debug_override_ == ALLOW_S2DEBUG) {
    // Verify that the original loops had consistent shell/hole orientations.
    // Each original loop L should have been inverted if and only if it now
    // represents a hole.
    for (int i = 0; i < loops_.size(); ++i) {
      DCHECK_EQ(contained_origin.count(loop(i)) != loop(i)->contains_origin(),
                loop(i)->is_hole());
    }
  }
}

void S2Polygon::InitLoopProperties() {
  has_holes_ = false;
  num_vertices_ = 0;
  bound_ = S2LatLngRect::Empty();
  for (int i = 0; i < num_loops(); ++i) {
    if (loop(i)->is_hole()) {
      has_holes_ = true;
    } else {
      bound_ = bound_.Union(loop(i)->GetRectBound());
    }
    num_vertices_ += loop(i)->num_vertices();
  }
  subregion_bound_ = S2EdgeUtil::RectBounder::ExpandForSubregions(bound_);
  InitIndex();
}

int S2Polygon::GetParent(int k) const {
  int depth = loop(k)->depth();
  if (depth == 0) return -1;  // Optimization.
  while (--k >= 0 && loop(k)->depth() >= depth) continue;
  return k;
}

int S2Polygon::GetLastDescendant(int k) const {
  if (k < 0) return num_loops() - 1;
  int depth = loop(k)->depth();
  while (++k < num_loops() && loop(k)->depth() > depth) continue;
  return k - 1;
}

double S2Polygon::GetArea() const {
  double area = 0;
  for (int i = 0; i < num_loops(); ++i) {
    area += loop(i)->sign() * loop(i)->GetArea();
  }
  return area;
}

S2Point S2Polygon::GetCentroid() const {
  S2Point centroid;
  for (int i = 0; i < num_loops(); ++i) {
    centroid += loop(i)->sign() * loop(i)->GetCentroid();
  }
  return centroid;
}

int S2Polygon::GetSnapLevel() const {
  int snap_level = -1;
  for (int i = 0; i < loops_.size(); ++i) {
    S2Loop const* child = loop(i);
    for (int j = 0; j < child->num_vertices(); ++j) {
      int face;
      unsigned int si, ti;
      int level = S2::XYZtoFaceSiTi(child->vertex(j), &face, &si, &ti);
      if (level < 0) return level;  // Vertex is not a cell center.
      if (level != snap_level) {
        if (snap_level < 0) {
          snap_level = level;  // First vertex.
        } else {
          return -1;  // Vertices at more than one cell level.
        }
      }
    }
  }
  return snap_level;
}

S1Angle S2Polygon::GetDistance(S2Point const& x) const {
  if (Contains(x)) return S1Angle::Zero();
  S2ClosestEdgeQuery query(index_);
  return query.GetDistance(x);
}

S1Angle S2Polygon::GetDistanceToBoundary(S2Point const& x) const {
  S2ClosestEdgeQuery query(index_);
  return query.GetDistance(x);
}

S2Point S2Polygon::Project(S2Point const& x) const {
  if (Contains(x)) return x;
  S2ClosestEdgeQuery query(index_);
  return query.Project(x);
}

S2Point S2Polygon::ProjectToBoundary(S2Point const& x) const {
  S2ClosestEdgeQuery query(index_);
  return query.Project(x);
}

// Return +1 if this polygon (A) contains the boundary of B, -1 if A excludes
// the boundary of B, and 0 if the boundaries of A and B cross.
int S2Polygon::CompareBoundary(S2Loop const* b) const {
  int result = -1;
  for (int i = 0; i < num_loops() && result != 0; ++i) {
    // If B crosses any loop of A, the result is 0.  Otherwise the result
    // changes sign each time B is contained by a loop of A.
    result *= -loop(i)->CompareBoundary(b);
  }
  return result;
}

// Return true if this polygon (A) contains the entire boundary of B.
bool S2Polygon::ContainsBoundary(S2Polygon const* b) const {
  for (int j = 0; j < b->num_loops(); ++j) {
    if (CompareBoundary(b->loop(j)) <= 0) return false;
  }
  return true;
}

// Return true if this polygon (A) excludes the entire boundary of B.
bool S2Polygon::ExcludesBoundary(S2Polygon const* b) const {
  for (int j = 0; j < b->num_loops(); ++j) {
    if (CompareBoundary(b->loop(j)) >= 0) return false;
  }
  return true;
}

// Given a polygon A and a loop B whose boundaries do not cross, return true
// if A contains the boundary of B.  Shared edges are handled according to the
// rule described in S2Loop::ContainsNonCrossingBoundary().
bool S2Polygon::ContainsNonCrossingBoundary(S2Loop const* b,
                                            bool reverse_b) const {
  bool inside = false;
  for (int i = 0; i < num_loops(); ++i) {
    inside ^= loop(i)->ContainsNonCrossingBoundary(b, reverse_b);
  }
  return inside;
}

// Given two polygons A and B such that the boundary of A does not cross any
// loop of B, return true if A excludes all shell boundaries of B.
bool S2Polygon::ExcludesNonCrossingShells(S2Polygon const* b) const {
  for (int j = 0; j < b->num_loops(); ++j) {
    if (b->loop(j)->is_hole()) continue;
    if (ContainsNonCrossingBoundary(b->loop(j), false /*reverse_b*/))
      return false;
  }
  return true;
}

// Given two polygons A and B such that the boundary of A does not cross any
// loop of B, return true if A excludes all shell boundaries of the complement
// of B.
bool S2Polygon::ExcludesNonCrossingComplementShells(S2Polygon const* b)
    const {
  // Special case to handle the complement of the empty or full polygons.
  if (b->is_empty()) return !is_full();
  if (b->is_full()) return true;

  // Otherwise the complement of B may be obtained by inverting loop(0) and
  // then swapping the shell/hole status of all other loops.  This implies
  // that the shells of the complement consist of loop 0 plus all the holes of
  // the original polygon.
  for (int j = 0; j < b->num_loops(); ++j) {
    if (j > 0 && !b->loop(j)->is_hole()) continue;

    // The interior of the complement is to the right of loop 0, and to the
    // left of the loops that were originally holes.
    if (ContainsNonCrossingBoundary(b->loop(j), j == 0 /*reverse_b*/))
      return false;
  }
  return true;
}

bool S2Polygon::AnyLoopContains(S2Loop const* b) const {
  // Return true if any loop contains the given loop.
  for (int i = 0; i < num_loops(); ++i) {
    if (loop(i)->Contains(b)) return true;
  }
  return false;
}

bool S2Polygon::AnyLoopIntersects(S2Loop const* b) const {
  // Return true if any loop intersects the given loop.
  for (int i = 0; i < num_loops(); ++i) {
    if (loop(i)->Intersects(b)) return true;
  }
  return false;
}

bool S2Polygon::Contains(S2Polygon const* b) const {
  // If both polygons have one loop, use the more efficient S2Loop method.
  // Note that S2Loop::Contains does its own bounding rectangle check.
  if (num_loops() == 1 && b->num_loops() == 1) {
    return loop(0)->Contains(b->loop(0));
  }

  // Otherwise if neither polygon has holes, we can still use the more
  // efficient S2Loop::Contains method (rather than CompareBoundary),
  // but it's worthwhile to do our own bounds check first.
  if (!subregion_bound_.Contains(b->bound_)) {
    // Even though Bound(A) does not contain Bound(B), it is still possible
    // that A contains B.  This can only happen when union of the two bounds
    // spans all longitudes.  For example, suppose that B consists of two
    // shells with a longitude gap between them, while A consists of one shell
    // that surrounds both shells of B but goes the other way around the
    // sphere (so that it does not intersect the longitude gap).
    if (!bound_.lng().Union(b->bound_.lng()).is_full()) return false;
  }
  if (!has_holes_ && !b->has_holes_) {
    for (int j = 0; j < b->num_loops(); ++j) {
      if (!AnyLoopContains(b->loop(j))) return false;
    }
    return true;
  }

  // Polygon A contains B iff B does not intersect the complement of A.  From
  // the intersection algorithm below, this means that the complement of A
  // must exclude the entire boundary of B, and B must exclude all shell
  // boundaries of the complement of A.  (It can be shown that B must then
  // exclude the entire boundary of the complement of A.)  The first call
  // below returns false if the boundaries cross, therefore the second call
  // does not need to check for any crossing edges (which makes it cheaper).
  return ContainsBoundary(b) && b->ExcludesNonCrossingComplementShells(this);
}

bool S2Polygon::Intersects(S2Polygon const* b) const {
  // If both polygons have one loop, use the more efficient S2Loop method.
  // Note that S2Loop::Intersects does its own bounding rectangle check.
  if (num_loops() == 1 && b->num_loops() == 1) {
    return loop(0)->Intersects(b->loop(0));
  }

  // Otherwise if neither polygon has holes, we can still use the more
  // efficient S2Loop::Intersects method.  The polygons intersect if and
  // only if some pair of loop regions intersect.
  if (!bound_.Intersects(b->bound_)) return false;
  if (!has_holes_ && !b->has_holes_) {
    for (int j = 0; j < b->num_loops(); ++j) {
      if (AnyLoopIntersects(b->loop(j))) return true;
    }
    return false;
  }

  // Polygon A is disjoint from B if A excludes the entire boundary of B and B
  // excludes all shell boundaries of A.  (It can be shown that B must then
  // exclude the entire boundary of A.)  The first call below returns false if
  // the boundaries cross, therefore the second call does not need to check
  // for crossing edges.
  return !ExcludesBoundary(b) || !b->ExcludesNonCrossingShells(this);
}

S2Cap S2Polygon::GetCapBound() const {
  return bound_.GetCapBound();
}

bool S2Polygon::Contains(S2Cell const& target) const {
  S2ShapeIndex::Iterator it(index_);
  S2ShapeIndex::CellRelation relation = it.Locate(target.id());

  // If "target" is disjoint from all index cells, it is not contained.
  // Similarly, if "target" is subdivided into one or more index cells then it
  // is not contained, since index cells are subdivided only if they (nearly)
  // intersect a sufficient number of edges.  (But note that if "target" itself
  // is an index cell then it may be contained, since it could be a cell with
  // no indexed edges in the polygon interior.)
  if (relation != S2ShapeIndex::INDEXED) return false;

  // Otherwise check if any edges intersect "target".  At this point, the
  // iterator is guaranteed to point to an index cell containing "target".
  if (BoundaryApproxIntersects(it, target)) return false;

  // Otherwise check if the polygon contains the center of "target".
  return Contains(it, target.GetCenter());
}

bool S2Polygon::ApproxContains(S2Polygon const* b,
                               S1Angle vertex_merge_radius) const {
  S2Polygon difference;
  difference.InitToDifferenceSloppy(b, this, vertex_merge_radius);
  return difference.num_loops() == 0;
}

bool S2Polygon::MayIntersect(S2Cell const& target) const {
  S2ShapeIndex::Iterator it(index_);
  S2ShapeIndex::CellRelation relation = it.Locate(target.id());

  // If "target" does not overlap any index cell, there is no intersection.
  if (relation == S2ShapeIndex::DISJOINT) return false;

  // If "target" properly contains an index cell then there is an intersection
  // to within the S2ShapeIndex error bound (see Contains).  But if "target"
  // itself is an index cell then it may not be contained, since it could be a
  // cell with no indexed edges that is contained by both a shell and a hole,
  // and is therefore in the polygon exterior.
  if (relation == S2ShapeIndex::SUBDIVIDED) return true;

  // Otherwise check if any edges intersect "target".
  if (BoundaryApproxIntersects(it, target)) return true;

  // Otherwise check if the polygon contains the center of "target".
  return Contains(it, target.GetCenter());
}

bool S2Polygon::BoundaryApproxIntersects(S2ShapeIndex::Iterator const& it,
                                         S2Cell const& target) const {
  DCHECK(it.id().contains(target.id()));
  S2ShapeIndexCell const* cell = it.cell();
  for (int a = 0; a < cell->num_shapes(); ++a) {
    S2ClippedShape const& a_clipped = cell->clipped(a);
    int a_num_clipped = a_clipped.num_edges();
    if (a_num_clipped == 0) continue;

    // We can save some work if "target" is the index cell itself (given that
    // there is at least one indexed edge).
    if (it.id() == target.id()) return true;

    // Otherwise check whether any of the edges intersect "target".
    static double const kMaxError = (S2EdgeUtil::kFaceClipErrorUVCoord +
                                     S2EdgeUtil::kIntersectsRectErrorUVDist);
    R2Rect bound = target.GetBoundUV().Expanded(kMaxError);
    S2Loop const* a_loop = loop(a_clipped.shape_id());
    for (int i = 0; i < a_num_clipped; ++i) {
      int ai = a_clipped.edge(i);
      R2Point v0, v1;
      if (S2EdgeUtil::ClipToPaddedFace(a_loop->vertex(ai), a_loop->vertex(ai+1),
                                       target.face(), kMaxError, &v0, &v1) &&
          S2EdgeUtil::IntersectsRect(v0, v1, bound)) {
        return true;
      }
    }
  }
  return false;
}

bool S2Polygon::VirtualContainsPoint(S2Point const& p) const {
  return Contains(p);  // The same as Contains() below, just virtual.
}

bool S2Polygon::Contains(S2Point const& p) const {
  // NOTE(ericv): A bounds check slows down this function by about 50%.  It is
  // worthwhile only when it might allow us to delay building the index.
  if (!index_.is_fresh() && !bound_.Contains(p)) return false;

  // For small polygons it is faster to just check all the crossings.
  // Otherwise we keep track of the number of calls to Contains() and only
  // build the index once enough calls have been made so that we think it is
  // worth the effort.  See S2Loop::Contains(S2Point) for detailed comments.
  static int const kMaxBruteForceVertices = 32;
  static int const kMaxUnindexedContainsCalls = 20;
  if (num_vertices() <= kMaxBruteForceVertices ||
      (!index_.is_fresh() &&
       base::subtle::Barrier_AtomicIncrement(&unindexed_contains_calls_, 1) !=
       kMaxUnindexedContainsCalls)) {
    bool inside = false;
    for (int i = 0; i < num_loops(); ++i) {
      // Use brute force to avoid building the loop's S2ShapeIndex.
      inside ^= loop(i)->BruteForceContains(p);
    }
    return inside;
  }
  // Otherwise we look up the S2ShapeIndex cell containing this point.
  S2ShapeIndex::Iterator it(index_);
  if (!it.Locate(p)) return false;
  return Contains(it, p);
}

bool S2Polygon::Contains(S2ShapeIndex::Iterator const& it,
                         S2Point const& p) const {
  // Test containment by drawing a line segment from the cell center to the
  // given point and counting edge crossings.
  S2ShapeIndexCell const* cell = it.cell();
  bool inside = false;
  bool crosser_initialized = false;
  S2EdgeUtil::EdgeCrosser crosser;
  S2Point center;
  for (int a = 0; a < cell->num_shapes(); ++a) {
    S2ClippedShape const& a_clipped = cell->clipped(a);
    inside ^= a_clipped.contains_center();
    int a_num_clipped = a_clipped.num_edges();
    if (a_num_clipped > 0) {
      if (!crosser_initialized) {
        center = it.center();
        crosser.Init(&center, &p);
        crosser_initialized = true;
      }
      int ai_prev = -2;
      S2Loop const* a_loop = loop(a_clipped.shape_id());
      for (int i = 0; i < a_num_clipped; ++i) {
        int ai = a_clipped.edge(i);
        if (ai != ai_prev + 1) crosser.RestartAt(&a_loop->vertex(ai));
        ai_prev = ai;
        inside ^= crosser.EdgeOrVertexCrossing(&a_loop->vertex(ai+1));
      }
    }
  }
  return inside;
}

void S2Polygon::Encode(Encoder* const encoder) const {
  if (num_vertices_ == 0) {
    EncodeCompressed(encoder, NULL, S2::kMaxCellLevel);
    return;
  }
  // Converts all the polygon vertices to S2XYZFaceSiTi format.
  FixedArray<S2XYZFaceSiTi> all_vertices(num_vertices_);
  S2XYZFaceSiTi* current_loop_vertices = all_vertices.get();
  for (int i = 0; i < loops_.size(); ++i) {
    loops_[i]->GetXYZFaceSiTiVertices(current_loop_vertices);
    current_loop_vertices += loops_[i]->num_vertices();
  }
  // Computes an histogram of the cell levels at which the vertices are snapped.
  // (histogram[0] is the number of unsnapped vertices, histogram[i] the number
  // of vertices snapped at level i-1).
  int histogram[S2::kMaxCellLevel + 2];
  for (int i = 0; i <= S2::kMaxCellLevel + 1; ++i) {
    histogram[i] = 0;
  }
  for (int i = 0; i < num_vertices_; ++i) {
    histogram[all_vertices[i].cell_level + 1] += 1;
  }
  // Computes the level at which most of the vertices are snapped.
  int snap_level = 0;
  int num_snapped = histogram[snap_level + 1];
  for (int i = 1; i <= S2::kMaxCellLevel; ++i) {
    if (histogram[i + 1] > num_snapped) {
      snap_level = i;
      num_snapped = histogram[i + 1];
    }
  }
  // Choose an encoding format based on the number of unsnapped vertices and a
  // rough estimate of the encoded sizes.

  // The compressed encoding requires approximately 4 bytes per vertex plus
  // "exact_point_size" for each unsnapped vertex (encoded as an S2Point plus
  // the index at which it is located).
  int exact_point_size = sizeof(S2Point) + 2;
  int num_unsnapped = num_vertices_ - num_snapped;
  int compressed_size = 4 * num_vertices_ + exact_point_size * num_unsnapped;
  int lossless_size = sizeof(S2Point) * num_vertices_;
  if (compressed_size < lossless_size) {
    EncodeCompressed(encoder, all_vertices.get(), snap_level);
  } else {
    EncodeLossless(encoder);
  }
}

void S2Polygon::EncodeLossless(Encoder* const encoder) const {
  encoder->Ensure(10);  // Sufficient
  encoder->put8(kCurrentLosslessEncodingVersionNumber);
  // This code used to write "owns_loops_", so write "true" for compatibility.
  encoder->put8(true);
  encoder->put8(has_holes_);
  encoder->put32(loops_.size());
  DCHECK_GE(encoder->avail(), 0);

  for (int i = 0; i < num_loops(); ++i) {
    loop(i)->Encode(encoder);
  }
  bound_.Encode(encoder);
}

bool S2Polygon::Decode(Decoder* const decoder) {
  if (decoder->avail() < sizeof(unsigned char)) return false;
  unsigned char version = decoder->get8();
  switch (version) {
    case kCurrentLosslessEncodingVersionNumber:
      return DecodeInternal(decoder, false);
    case kCurrentCompressedEncodingVersionNumber:
      return DecodeCompressed(decoder);
  }
  return false;
}

bool S2Polygon::DecodeWithinScope(Decoder* const decoder) {
  if (decoder->avail() < sizeof(unsigned char)) return false;
  unsigned char version = decoder->get8();
  switch (version) {
    case kCurrentLosslessEncodingVersionNumber:
      return DecodeInternal(decoder, true);
    case kCurrentCompressedEncodingVersionNumber:
      return DecodeCompressed(decoder);
  }
  return false;
}

bool S2Polygon::DecodeInternal(Decoder* const decoder, bool within_scope) {
  if (decoder->avail() < 2 * sizeof(uint8) + sizeof(uint32)) return false;
  ClearLoops();
  decoder->get8();  // Ignore irrelevant serialized owns_loops_ value.
  has_holes_ = decoder->get8();
  int num_loops = decoder->get32();
  loops_.reserve(num_loops);
  num_vertices_ = 0;
  for (int i = 0; i < num_loops; ++i) {
    loops_.push_back(new S2Loop);
    loops_.back()->set_s2debug_override(s2debug_override());
    if (within_scope) {
      if (!loops_.back()->DecodeWithinScope(decoder)) return false;
    } else {
      if (!loops_.back()->Decode(decoder)) return false;
    }
    num_vertices_ += loops_.back()->num_vertices();
  }
  if (!bound_.Decode(decoder)) return false;
  subregion_bound_ = S2EdgeUtil::RectBounder::ExpandForSubregions(bound_);
  InitIndex();
  return true;
}

// IntersectionSet represents a set of intersections with a particular edge.
// Each intersection is a pair (t, P) where "t" is the interpolation parameter
// along the edge (a fraction in the range [0,1]), and P is the actual
// intersection point.
typedef vector<pair<double, S2Point> > IntersectionSet;

// EdgeClipper find all the intersections of a given edge with the edges
// contained in an S2ShapeIndex.  It is used to implement polygon operations
// such as intersection and union.
class EdgeClipper {
 public:
  // Initialize an EdgeClipper for the given S2ShapeIndex.  If the query edge
  // is the same as an index edge (a "shared edge"), then the edge will be
  // included in the output if and only if "add_shared_edges" is true.  The
  // "reverse_edges" function allows the edges of any index shape to be
  // reversed before this test is performed (this is used to reverse the loop
  // orientation of "holes" in certain algorithms).
  EdgeClipper(S2ShapeIndex const& index, bool add_shared_edges,
              bool (*reverse_edges)(S2Shape const* shape));

  // Find all intersections of the edge (a0, a1) with edges in the index and
  // append them to "intersections".  (The result is unsorted.)
  void ClipEdge(S2Point const& a0, S2Point const& a1,
                IntersectionSet* intersections);

 private:
  // Given two edges A and B such that RobustCrossing(A, B) >= 0, determine if
  // they intersect and add any intersection point to "intersections_".
  // "b_shape" is the S2Shape containing edge B, and "crossing" is the result
  // of RobustCrossing(A, B).
  void AddIntersection(S2Point const& a0, S2Point const& a1,
                       S2Point const& b0, S2Point const& b1,
                       S2Shape const* b_shape, int crossing);

  bool const add_shared_edges_;
  bool (* const reverse_edges_)(S2Shape const* shape);

  // Temporary variables used while processing a query, declared here to
  // reduce memory allocation and argument-passing overhead.
  S2EdgeQuery query_;
  S2EdgeQuery::EdgeMap edge_map_;
  IntersectionSet* intersections_;
};

EdgeClipper::EdgeClipper(S2ShapeIndex const& index, bool add_shared_edges,
                         bool (*reverse_edges)(S2Shape const* shape))
    : add_shared_edges_(add_shared_edges),
      reverse_edges_(reverse_edges) {
  query_.Init(index);
}

void EdgeClipper::AddIntersection(S2Point const& a0, S2Point const& a1,
                                  S2Point const& b0, S2Point const& b1,
                                  S2Shape const* b_shape, int crossing) {
  DCHECK_GE(crossing, 0);
  if (crossing > 0) {
    // There is a proper edge crossing.
    S2Point x = S2EdgeUtil::GetIntersection(a0, a1, b0, b1);
    double t = S2EdgeUtil::GetDistanceFraction(x, a0, a1);
    intersections_->push_back(std::make_pair(t, x));

  } else if (S2EdgeUtil::VertexCrossing(a0, a1, b0, b1)) {
    // There is a crossing at one of the vertices.  The basic rule is simple:
    // if a0 equals one of the "b" vertices, the crossing occurs at t=0;
    // otherwise, it occurs at t=1.
    //
    // This has the effect that when two reversed edges exist (i.e., a0==b1
    // and a1==b0) and each edge is clipped against the other, neither one is
    // included in the output (which is correct).  However when two shared
    // edges exist (i.e., a0==b0 and a1==b1), both are included in the output
    // (which is incorrect).  The "add_shared_edges" flag gives explicit
    // control over whether these shared edges are included in the output; if
    // it is false, shared edges are excluded by changing their intersection
    // parameter from 0 to 1.  This allows exactly one copy of shared edges to
    // be preserved, by calling ClipBoundary() twice with opposite values of
    // the "add_shared_edges" flag.
    double t = (a0 == b0 || a0 == b1) ? 0 : 1;
    if (!add_shared_edges_ && a1 == (reverse_edges_(b_shape) ? b0 : b1)) {
      t = 1;  // Excludes (a0,a1) from the output.
    }
    intersections_->push_back(std::make_pair(t, t == 0 ? a0 : a1));
  }
}

// Find all points where the polygon B intersects the edge (a0,a1),
// and add the corresponding parameter values (in the range [0,1]) to
// "intersections".
void EdgeClipper::ClipEdge(S2Point const& a0, S2Point const& a1,
                           IntersectionSet* intersections) {
  if (!query_.GetCandidates(a0, a1, &edge_map_)) return;

  // Iterate through the candidate loops, and then the candidate edges within
  // each loop.
  intersections_ = intersections;
  S2EdgeUtil::EdgeCrosser crosser(&a0, &a1);
  for (S2EdgeQuery::EdgeMap::const_iterator it = edge_map_.begin();
       it != edge_map_.end(); ++it) {
    S2Shape const* b_shape = it->first;
    vector<int> const& b_candidates = it->second;
    int n = b_candidates.size();
    S2Point const *b0, *b1 = NULL;
    for (int j = 0; j < n; ++j) {
      S2Point const* b1_prev = b1;
      b_shape->GetEdge(b_candidates[j], &b0, &b1);
      if (b0 != b1_prev) crosser.RestartAt(b0);
      int crossing = crosser.RobustCrossing(b1);
      if (crossing < 0) continue;
      AddIntersection(a0, a1, *b0, *b1, b_shape, crossing);
    }
  }
}

static bool ReverseHoles(S2Shape const* shape) {
  return down_cast<S2Loop::Shape const*>(shape)->loop()->is_hole();
}

// Clip the boundary of A to the interior of B, and add the resulting edges to
// "builder".  Shells are directed CCW and holes are directed clockwise.  If
// "reverse_a" is true, these directions are reversed in polygon A.  If
// "invert_b" is true, the boundary of A is clipped to the exterior rather
// than the interior of B.  If "add_shared_edges" is true, then the output
// will include any edges that are shared between A and B (both edges must be
// in the same direction after any edge reversals are taken into account).
void S2Polygon::ClipBoundary(S2Polygon const* a, bool reverse_a,
                             S2Polygon const* b, bool invert_b,
                             bool add_shared_edges, S2PolygonBuilder* builder) {
  EdgeClipper clipper(b->index_, add_shared_edges, ReverseHoles);
  IntersectionSet intersections;
  for (int i = 0; i < a->num_loops(); ++i) {
    S2Loop const* a_loop = a->loop(i);
    int n = a_loop->num_vertices();
    int dir = (a_loop->is_hole() ^ reverse_a) ? -1 : 1;
    bool inside = b->Contains(a_loop->vertex(0)) ^ invert_b;
    for (int j = (dir > 0) ? 0 : n; n > 0; --n, j += dir) {
      S2Point const& a0 = a_loop->vertex(j);
      S2Point const& a1 = a_loop->vertex(j + dir);
      clipper.ClipEdge(a0, a1, &intersections);
      if (inside) intersections.push_back(std::make_pair(0, a0));
      inside = (intersections.size() & 1);
      DCHECK_EQ((b->Contains(a1) ^ invert_b), inside);
      if (inside) intersections.push_back(std::make_pair(1, a1));
      std::sort(intersections.begin(), intersections.end());
      for (int k = 0; k < intersections.size(); k += 2) {
        if (intersections[k] == intersections[k+1]) continue;
        builder->AddEdge(intersections[k].second, intersections[k+1].second);
      }
      intersections.clear();  // Reuse to avoid repeated memory allocation
    }
  }
}

void S2Polygon::Invert() {
  // Inverting any one loop will invert the polygon.  The best loop to invert
  // is the one whose area is largest, since this yields the smallest area
  // after inversion.  The loop with the largest area is always at depth 0.
  // The descendents of this loop all have their depth reduced by 1, while the
  // former siblings of this loop all have their depth increased by 1.

  // The empty and full polygons are handled specially.
  if (is_empty()) {
    loops_.push_back(new S2Loop(S2Loop::kFull()));
  } else if (is_full()) {
    ClearLoops();
  } else {
    // Find the loop whose area is largest (i.e., whose turning angle is
    // smallest), minimizing calls to GetTurningAngle().  In particular, for
    // polygons with a single shell at level 0 there is not need to call
    // GetTurningAngle() at all.  (This method is relatively expensive.)
    int best = 0;
    double const kNone = 10.0;  // Flag that means "not computed yet"
    double best_angle = kNone;
    for (int i = 1; i < num_loops(); ++i) {
      if (loop(i)->depth() == 0) {
        // We defer computing the turning angle of loop 0 until we discover
        // that the polygon has another top-level shell.
        if (best_angle == kNone) best_angle = loop(best)->GetTurningAngle();
        double angle = loop(i)->GetTurningAngle();
        if (angle < best_angle) {
          best = i;
          best_angle = angle;
        }
      }
    }
    // Build the new loops vector, starting with the inverted loop.
    loop(best)->Invert();
    vector<S2Loop*> new_loops;
    new_loops.reserve(num_loops());
    new_loops.push_back(loop(best));
    // Add the former siblings of this loop as descendants.
    int last_best = GetLastDescendant(best);
    for (int i = 0; i < num_loops(); ++i) {
      if (i < best || i > last_best) {
        loop(i)->set_depth(loop(i)->depth() + 1);
        new_loops.push_back(loop(i));
      }
    }
    // Add the former children of this loop as siblings.
    for (int i = 0; i < num_loops(); ++i) {
      if (i > best && i <= last_best) {
        loop(i)->set_depth(loop(i)->depth() - 1);
        new_loops.push_back(loop(i));
      }
    }
    loops_.swap(new_loops);
    DCHECK_EQ(new_loops.size(), num_loops());
  }
  index_.Reset();
  InitLoopProperties();
}

void S2Polygon::InitToComplement(S2Polygon const* a) {
  Copy(a);
  Invert();
}

void S2Polygon::InitToIntersection(S2Polygon const* a, S2Polygon const* b) {
  InitToIntersectionSloppy(a, b, S2EdgeUtil::kIntersectionTolerance);
}

void S2Polygon::InitToIntersectionSloppy(S2Polygon const* a, S2Polygon const* b,
                                         S1Angle vertex_merge_radius) {
  if (!a->bound_.Intersects(b->bound_)) return;

  // We want the boundary of A clipped to the interior of B,
  // plus the boundary of B clipped to the interior of A,
  // plus one copy of any directed edges that are in both boundaries.

  S2PolygonBuilderOptions options(S2PolygonBuilderOptions::DIRECTED_XOR());
  options.set_vertex_merge_radius(vertex_merge_radius);
  S2PolygonBuilder builder(options);
  ClipBoundary(a, false, b, false, true, &builder);
  ClipBoundary(b, false, a, false, false, &builder);
  if (!builder.AssemblePolygon(this, NULL)) {
    LOG(DFATAL) << "Bad directed edges in InitToIntersection";
  }
  // If the result had a non-empty boundary then we are done.  Unfortunately,
  // if the boundary is empty then there are two possible results: the empty
  // polygon or the full polygon.  This choice would be trivial to resolve
  // except for the existence of "vertex_merge_radius" and also numerical
  // errors when computing edge intersection points.  In particular:
  //
  //  - The intersection of two non-full polygons may be full.  For example,
  //    one or both polygons may have tiny cracks that are eliminated due to
  //    vertex merging/edge splicing.
  //
  //  - The intersection of two polygons that both contain S2::Origin() (or
  //    any other point) may be empty.  For example, both polygons may have
  //    tiny shells that surround the common point but that are eliminated.
  //
  //  - Even before any vertex merging/edge splicing, the computed boundary
  //    edges are not useful in distinguishing almost-full polygons from
  //    almost-empty due to numerical errors in computing edge intersections.
  //    Such errors can reverse the orientation of narrow cracks or slivers.
  //
  // So instead we fall back to heuristics.  Essentially we compute the
  // minimum and maximum intersection area based on the areas of the two input
  // polygons.  If only one of {0, 4*Pi} is possible then we return that
  // result.  If neither is possible (before vertex merging, etc) then we
  // return the one that is closest to being possible.  (It never true that
  // both are possible.)
  if (num_loops() == 0) {
    // We know that both polygons are non-empty due to the initial bounds
    // check.  By far the most common case is that the intersection is empty,
    // so we want to make that case fast.  The intersection area satisfies:
    //
    //   max(0, A + B - 4*Pi) <= Intersection(A, B) <= min(A, B)
    //
    // where A, B refer to a polygon and/or its area.  Note that if either A
    // or B is at most 2*Pi, the result must be "empty".  We can use the
    // bounding rectangle areas as upper bounds on the polygon areas.
    if (a->bound_.Area() <= 2 * M_PI || b->bound_.Area() <= 2 * M_PI) return;
    double a_area = a->GetArea(), b_area = b->GetArea();
    double min_area = max(0.0, a_area + b_area - 4 * M_PI);
    double max_area = min(a_area, b_area);
    if (min_area > 4 * M_PI - max_area)
      Invert();
  }
}

void S2Polygon::InitToUnion(S2Polygon const* a, S2Polygon const* b) {
  InitToUnionSloppy(a, b, S2EdgeUtil::kIntersectionTolerance);
}

void S2Polygon::InitToUnionSloppy(S2Polygon const* a, S2Polygon const* b,
                                  S1Angle vertex_merge_radius) {
  // We want the boundary of A clipped to the exterior of B,
  // plus the boundary of B clipped to the exterior of A,
  // plus one copy of any directed edges that are in both boundaries.

  S2PolygonBuilderOptions options(S2PolygonBuilderOptions::DIRECTED_XOR());
  options.set_vertex_merge_radius(vertex_merge_radius);
  S2PolygonBuilder builder(options);
  ClipBoundary(a, false, b, true, true, &builder);
  ClipBoundary(b, false, a, true, false, &builder);
  if (!builder.AssemblePolygon(this, NULL)) {
    LOG(DFATAL) << "Bad directed edges";
  }
  if (num_loops() == 0) {
    // See comments in InitToIntersectionSloppy().  In this case, the union
    // area satisfies:
    //
    //   max(A, B) <= Union(A, B) <= min(4*Pi, A + B)
    //
    // The most common case is that neither input polygon is empty, but the
    // union is empty due to vertex merging/simplification.
    if (a->bound_.Area() + b->bound_.Area() <= 2 * M_PI) return;
    double a_area = a->GetArea(), b_area = b->GetArea();
    double min_area = max(a_area, b_area);
    double max_area = min(4 * M_PI, a_area + b_area);
    if (min_area > 4 * M_PI - max_area)
      Invert();
  }
}

void S2Polygon::InitToDifference(S2Polygon const* a, S2Polygon const* b) {
  InitToDifferenceSloppy(a, b, S2EdgeUtil::kIntersectionTolerance);
}

void S2Polygon::InitToDifferenceSloppy(S2Polygon const* a, S2Polygon const* b,
                                       S1Angle vertex_merge_radius) {
  // We want the boundary of A clipped to the exterior of B,
  // plus the reversed boundary of B clipped to the interior of A,
  // plus one copy of any edge in A that is also a reverse edge in B.

  S2PolygonBuilderOptions options(S2PolygonBuilderOptions::DIRECTED_XOR());
  options.set_vertex_merge_radius(vertex_merge_radius);
  S2PolygonBuilder builder(options);
  ClipBoundary(a, false, b, true, true, &builder);
  ClipBoundary(b, true, a, false, false, &builder);
  if (!builder.AssemblePolygon(this, NULL)) {
    LOG(DFATAL) << "Bad directed edges in InitToDifference";
  }
  if (num_loops() == 0) {
    // See comments in InitToIntersectionSloppy().  In this case, the
    // difference area satisfies:
    //
    //   max(0, A - B) <= Difference(A, B) <= min(A, 4*Pi - B)
    //
    // By far the most common case is that result is empty.
    if (a->bound_.Area() <= 2 * M_PI || b->bound_.Area() >= 2 * M_PI) return;
    double a_area = a->GetArea(), b_area = b->GetArea();
    double min_area = max(0.0, a_area - b_area);
    double max_area = min(a_area, 4 * M_PI - b_area);
    if (min_area > 4 * M_PI - max_area)
      Invert();
  }
}

// An abstract class used to specify which vertices must be kept when
// simplifying a loop in SimplifyLoopAsPolyline below.
class S2VertexFilter {
 public:
  S2VertexFilter() { }
  virtual ~S2VertexFilter() { }

  // Returns true if the given vertex must be kept in the simplified loop
  // (returning false does not mean that the vertex *must* be removed, but that
  // it *can* be removed if it satisfies the simplification criteria).
  virtual bool ShouldKeepVertex(const S2Point& vertex) const = 0;

 private:
  DISALLOW_COPY_AND_ASSIGN(S2VertexFilter);
};

// An S2VertexFilter subclass to keep the vertices that are close to the
// boundary of a given S2 cell.
namespace {
class CellBoundaryVertexFilter : public S2VertexFilter {
 public:
  explicit CellBoundaryVertexFilter(const S2Cell& cell) {
    // Computes and stores the cell vertices and their dot products, to optimize
    // computations ShouldKeepVertex.
    for (int i = 0; i < 4; ++i) {
      cell_vertices_[i] = cell.GetVertex(i);
    }
  }

  virtual bool ShouldKeepVertex(const S2Point& vertex) const {
    const S1Angle kTolerance = S1Angle::Radians(1e-15);
    for (int i = 0; i < 4; ++i) {
      if (S2EdgeUtil::GetDistance(vertex,
                                  cell_vertices_[i],
                                  cell_vertices_[(i + 1) % 4]) < kTolerance) {
        return true;
      }
    }
    return false;
  }

 private:
  S2Point cell_vertices_[4];
};
}  // namespace

// Takes a loop and simplifies it.  This may return a self-intersecting
// polyline.  Always keeps the first vertex from the loop, as well as the
// vertices for which should_keep returns true. This function does not take
// ownership of should_keep, which can be NULL.
vector<S2Point>* SimplifyLoopAsPolyline(S2Loop const* loop, S1Angle tolerance,
                                        const S2VertexFilter* should_keep) {
  vector<S2Point> points(loop->num_vertices() + 1);
  // Add the first vertex at the beginning and at the end.
  for (int i = 0; i <= loop->num_vertices(); ++i) {
    points[i] = loop->vertex(i);
  }
  S2Polyline line(points);
  vector<int> indices;
  line.SubsampleVertices(tolerance, &indices);
  if (indices.size() <= 2) return NULL;
  if (should_keep != NULL) {
    vector<int> to_keep;
    // Remember the indices of the vertices that we must keep.
    for (int i = 0; i <= loop->num_vertices(); ++i) {
      if (should_keep->ShouldKeepVertex(loop->vertex(i))) {
        to_keep.push_back(i);
      }
    }
    if (!to_keep.empty()) {
      vector<int> new_indices;
      STLSetUnion(indices, to_keep, &new_indices);
      swap(indices, new_indices);
    }
  }
  // Add them all except the last: it is the same as the first.
  vector<S2Point>* simplified_line = new vector<S2Point>(indices.size() - 1);
  VLOG(4) << "Now simplified to: ";
  for (int i = 0; i + 1 < indices.size(); ++i) {
    (*simplified_line)[i] = line.vertex(indices[i]);
    VLOG(4) << S2LatLng(line.vertex(indices[i]));
  }
  return simplified_line;
}

static bool ReverseNone(S2Shape const* shape) {
  return false;
}

// Takes a set of possibly intersecting edges, stored in an S2ShapeIndex.
// Breaks the edges into small pieces so that there is no intersection
// anymore, and adds all these edges to the builder.
void BreakEdgesAndAddToBuilder(S2ShapeIndex const& index,
                               S2PolygonBuilder* builder) {
  // If there are self intersections, we add the pieces separately.
  // add_shared_edges ("true" below) can be false or true: it makes no
  // difference due to the way we call ClipEdge.
  EdgeClipper clipper(index, true, ReverseNone);
  IntersectionSet intersections;
  for (int id = 0; id < index.num_shape_ids(); ++id) {
    S2Shape const* shape = index.shape(id);
    int num_edges = shape->num_edges();
    for (int e = 0; e < num_edges; ++e) {
      S2Point const *from, *to;
      shape->GetEdge(e, &from, &to);
      intersections.push_back(std::make_pair(0, *from));
      clipper.ClipEdge(*from, *to, &intersections);
      intersections.push_back(std::make_pair(1, *to));
      std::sort(intersections.begin(), intersections.end());
      for (int k = 0; k + 1 < intersections.size(); ++k) {
        if (intersections[k] == intersections[k+1]) continue;
        builder->AddEdge(intersections[k].second, intersections[k+1].second);
      }
      intersections.clear();
    }
  }
}

namespace {
// A loop represented as a vector of S2Points.  As with S2Loop, the last
// vertex is implicitly joined to the first.
class PointVectorLoopShape : public S2Shape {
 public:
  explicit PointVectorLoopShape(vector<S2Point> const* vertices)
      : num_vertices_(vertices->size()),
        vertices_(&*vertices->begin()) {
  }
  int num_edges() const { return num_vertices_; }
  void GetEdge(int i, S2Point const** a, S2Point const** b) const {
    *a = &vertices_[i];
    *b = &vertices_[i+1 == num_vertices_ ? 0 : i+1];
  }
  bool has_interior() const { return false; }  // Not needed
  bool contains_origin() const { return false; }
  void Release() const { delete this; }
 protected:
  int num_vertices_;
  S2Point const* vertices_;
};
}  // namespace

// Simplifies the polygon.   The algorithm is straightforward and naive:
//   1. Simplify each loop by removing points while staying in the
//      tolerance zone.  This uses S2Polyline::SubsampleVertices(),
//      which is not guaranteed to be optimal in terms of number of
//      points.
//   2. Break any edge in pieces such that no piece intersects any
//      other.
//   3. Use the polygon builder to regenerate the full polygon.
//   4. If should_keep is not NULL, the vertices for which it returns true are
//      kept in the simplified polygon.
void S2Polygon::InitToSimplifiedInternal(S2Polygon const* a,
                                         S1Angle tolerance,
                                         bool snap_to_cell_centers,
                                         S2VertexFilter const* should_keep) {
  S2PolygonBuilderOptions builder_options =
      S2PolygonBuilderOptions::UNDIRECTED_XOR();
  builder_options.set_validate(false);
  if (should_keep != NULL) {
    // If there is a vertex filter, then we want to do as little vertex
    // merging as possible so that the vertices we want to keep don't move.
    // But on the other hand, when we break intersecting edges into pieces
    // there is some error in the intersection point.  S2PolygonBuilder needs
    // to be able to move vertices by up to this amount in order to produce
    // valid output.
    builder_options.set_vertex_merge_radius(S2EdgeUtil::kIntersectionTolerance);
  } else {
    // Ideally, we would want to set the vertex_merge_radius of the
    // builder roughly to tolerance (and in fact forego the edge
    // splitting step).  Alas, if we do that, we are liable to the
    // 'chain effect', where vertices are merged with closeby vertices
    // and so on, so that a vertex can move by an arbitrary distance.
    // So we remain conservative:
    builder_options.set_vertex_merge_radius(tolerance * 0.10);
    builder_options.set_snap_to_cell_centers(snap_to_cell_centers);
  }
  S2PolygonBuilder builder(builder_options);

  // Simplify each loop separately and add to the edge index.
  S2ShapeIndex index;
  vector<vector<S2Point>*> simplified_loops;
  for (int i = 0; i < a->num_loops(); ++i) {
    vector<S2Point>* simpler = SimplifyLoopAsPolyline(a->loop(i), tolerance,
        should_keep);
    if (NULL == simpler) continue;
    simplified_loops.push_back(simpler);
    index.Insert(new PointVectorLoopShape(simpler));
  }
  if (index.num_shape_ids() > 0) {
    BreakEdgesAndAddToBuilder(index, &builder);
    if (!builder.AssemblePolygon(this, NULL)) {
      LOG(DFATAL) << "Bad edges in InitToSimplified.";
    }
    // If there are no loops, check whether the result should be the full
    // polygon rather than the empty one.  (See InitToIntersectionSloppy.)
    if (num_loops() == 0) {
      if (a->bound_.Area() > 2 * M_PI && a->GetArea() > 2 * M_PI) Invert();
    }
  }
  for (int i = 0; i < simplified_loops.size(); ++i) {
    delete simplified_loops[i];
  }
  simplified_loops.clear();
}

void S2Polygon::InitToSimplified(S2Polygon const* a, S1Angle tolerance,
                                 bool snap_to_cell_centers) {
  InitToSimplifiedInternal(a, tolerance, snap_to_cell_centers, NULL);
}

void S2Polygon::InitToSimplifiedInCell(S2Polygon const* a, S2Cell const& cell,
                                       S1Angle tolerance) {
  CellBoundaryVertexFilter cell_filter(cell);
  InitToSimplifiedInternal(a, tolerance, false, &cell_filter);
}

void S2Polygon::InitToSnapped(S2Polygon const* a, int snap_level) {
  // TODO(user): Remove (tolerance * 0.1) from InitToSimplified
  //   and use that instead.
  S2PolygonBuilderOptions options;
  // Ensure that there will be no two vertices within the max leaf cell
  // diagonal of each other, therefore no two vertices in the same leaf cell,
  // and that no vertex will cross an edge after the points have been snapped
  // to the centers of leaf cells.  Add 1e-15 to the tolerance so we don't
  // set a tighter than leaf cell level because of numerical inaccuracy.
  options.SetRobustnessRadius(
      S1Angle::Radians(S2::kMaxDiag.GetValue(snap_level) / 2 + 1e-15));
  options.set_snap_to_cell_centers(true);

  S2PolygonBuilder polygon_builder(options);
  polygon_builder.AddPolygon(a);
  if (!polygon_builder.AssemblePolygon(this, NULL)) {
    LOG(DFATAL) << "AssemblePolygon failed in BuildSnappedPolygon";
  }
  // If there are no loops, check whether the result should be the full
  // polygon rather than the empty one.  (See InitToIntersectionSloppy.)
  if (num_loops() == 0) {
    if (a->bound_.Area() > 2 * M_PI && a->GetArea() > 2 * M_PI) Invert();
  }
}

void S2Polygon::InternalClipPolyline(bool invert,
                                     S2Polyline const* a,
                                     vector<S2Polyline*> *out,
                                     S1Angle merge_radius) const {
  // Clip the polyline A to the interior of this polygon.
  // The resulting polyline(s) will be appended to 'out'.
  // If invert is true, we clip A to the exterior of this polygon instead.
  // Vertices will be dropped such that adjacent vertices will not
  // be closer than 'merge_radius'.
  //
  // We do the intersection/subtraction by walking the polyline edges.
  // For each edge, we compute all intersections with the polygon boundary
  // and sort them in increasing order of distance along that edge.
  // We then divide the intersection points into pairs, and output a
  // clipped polyline segment for each one.
  // We keep track of whether we're inside or outside of the polygon at
  // all times to decide which segments to output.

  DCHECK(out->empty());

  EdgeClipper clipper(index_, true, ReverseNone);
  IntersectionSet intersections;
  vector<S2Point> vertices;
  int n = a->num_vertices();
  bool inside = Contains(a->vertex(0)) ^ invert;
  for (int j = 0; j < n-1; j++) {
    S2Point const& a0 = a->vertex(j);
    S2Point const& a1 = a->vertex(j + 1);
    clipper.ClipEdge(a0, a1, &intersections);
    if (inside) intersections.push_back(std::make_pair(0, a0));
    inside = (intersections.size() & 1);
    DCHECK_EQ((Contains(a1) ^ invert), inside);
    if (inside) intersections.push_back(std::make_pair(1, a1));
    std::sort(intersections.begin(), intersections.end());
    // At this point we have a sorted array of vertex pairs representing
    // the edge(s) obtained after clipping (a0,a1) against the polygon.
    for (int k = 0; k < intersections.size(); k += 2) {
      if (intersections[k] == intersections[k+1]) continue;
      S2Point const& v0 = intersections[k].second;
      S2Point const& v1 = intersections[k+1].second;

      // If the gap from the previous vertex to this one is large
      // enough, start a new polyline.
      if (!vertices.empty() && S1Angle(vertices.back(), v0) > merge_radius) {
        out->push_back(new S2Polyline(vertices));
        vertices.clear();
      }
      // Append this segment to the current polyline, ignoring any
      // vertices that are too close to the previous vertex.
      if (vertices.empty()) vertices.push_back(v0);
      if (S1Angle(vertices.back(), v1) > merge_radius) {
        vertices.push_back(v1);
      }
    }
    intersections.clear();
  }
  if (!vertices.empty()) {
    out->push_back(new S2Polyline(vertices));
  }
}

void S2Polygon::IntersectWithPolyline(
    S2Polyline const* a,
    vector<S2Polyline*> *out) const {
  IntersectWithPolylineSloppy(a, out, S2EdgeUtil::kIntersectionTolerance);
}

void S2Polygon::IntersectWithPolylineSloppy(
    S2Polyline const* a,
    vector<S2Polyline*> *out,
    S1Angle vertex_merge_radius) const {
  InternalClipPolyline(false, a, out, vertex_merge_radius);
}

void S2Polygon::SubtractFromPolyline(S2Polyline const* a,
                                     vector<S2Polyline*> *out) const {
  SubtractFromPolylineSloppy(a, out, S2EdgeUtil::kIntersectionTolerance);
}

void S2Polygon::SubtractFromPolylineSloppy(
    S2Polyline const* a,
    vector<S2Polyline*> *out,
    S1Angle vertex_merge_radius) const {
  InternalClipPolyline(true, a, out, vertex_merge_radius);
}


S2Polygon* S2Polygon::DestructiveUnion(vector<S2Polygon*>* polygons) {
  return DestructiveUnionSloppy(polygons, S2EdgeUtil::kIntersectionTolerance);
}

S2Polygon* S2Polygon::DestructiveUnionSloppy(vector<S2Polygon*>* polygons,
                                             S1Angle vertex_merge_radius) {
  // Effectively create a priority queue of polygons in order of number of
  // vertices.  Repeatedly union the two smallest polygons and add the result
  // to the queue until we have a single polygon to return.
  typedef std::multimap<int, S2Polygon*> QueueType;
  QueueType queue;  // Map from # of vertices to polygon.
  for (int i = 0; i < polygons->size(); ++i)
    queue.insert(
        std::make_pair((*polygons)[i]->num_vertices(), (*polygons)[i]));
  polygons->clear();

  while (queue.size() > 1) {
    // Pop two simplest polygons from queue.
    QueueType::iterator smallest_it = queue.begin();
    int a_size = smallest_it->first;
    S2Polygon* a_polygon = smallest_it->second;
    queue.erase(smallest_it);
    smallest_it = queue.begin();
    int b_size = smallest_it->first;
    S2Polygon* b_polygon = smallest_it->second;
    queue.erase(smallest_it);

    // Union and add result back to queue.
    S2Polygon* union_polygon = new S2Polygon();
    union_polygon->InitToUnionSloppy(a_polygon, b_polygon, vertex_merge_radius);
    delete a_polygon;
    delete b_polygon;
    queue.insert(std::make_pair(a_size + b_size, union_polygon));
    // We assume that the number of vertices in the union polygon is the
    // sum of the number of vertices in the original polygons, which is not
    // always true, but will almost always be a decent approximation, and
    // faster than recomputing.
  }

  if (queue.empty())
    return new S2Polygon();
  else
    return queue.begin()->second;
}

void S2Polygon::InitToCellUnionBorder(S2CellUnion const& cells) {
  // Use a polygon builder to union the cells in the union.  Due to rounding
  // errors, we can't do an exact union - when a small cell is adjacent to a
  // larger cell, the shared edges can fail to line up exactly.  Two cell edges
  // cannot come closer then kMinWidth, so if we have the polygon builder merge
  // edges within half that distance, we should always merge shared edges
  // without merging different edges.
  S2PolygonBuilderOptions options(S2PolygonBuilderOptions::DIRECTED_XOR());
  double min_cell_angle = S2::kMinWidth.GetValue(S2CellId::kMaxLevel);
  options.set_vertex_merge_radius(S1Angle::Radians(min_cell_angle / 2));
  S2PolygonBuilder builder(options);
  for (int i = 0; i < cells.num_cells(); ++i) {
    S2Loop cell_loop(S2Cell(cells.cell_id(i)));
    builder.AddLoop(&cell_loop);
  }
  if (!builder.AssemblePolygon(this, NULL)) {
    LOG(DFATAL) << "AssemblePolygon failed in InitToCellUnionBorder";
  }
  // If there are no loops, check whether the result should be the full
  // polygon rather than the empty one.  There are only two ways that this
  // should happen: either the cell union is empty, or it consists of all six
  // faces.
  if (num_loops() == 0) {
    if (cells.num_cells() == 0) return;
    DCHECK_EQ(static_cast<uint64>(6) << (2 * S2CellId::kMaxLevel),
              cells.LeafCellsCovered());
    Invert();
  }
}

bool S2Polygon::IsNormalized() const {
  // TODO(ericv): The condition tested here is insufficient.  The correct
  // condition is that each *connected component* of child loops can share at
  // most one vertex with their parent loop.  Example: suppose loop A has
  // children B, C, D, and the following pairs are connected: AB, BC, CD, DA.
  // Then the polygon is not normalized.
  set<S2Point> vertices;
  S2Loop const* last_parent = NULL;
  for (int i = 0; i < num_loops(); ++i) {
    S2Loop const* child = loop(i);
    if (child->depth() == 0) continue;
    S2Loop const* parent = loop(GetParent(i));
    if (parent != last_parent) {
      vertices.clear();
      for (int j = 0; j < parent->num_vertices(); ++j) {
        vertices.insert(parent->vertex(j));
      }
      last_parent = parent;
    }
    int count = 0;
    for (int j = 0; j < child->num_vertices(); ++j) {
      if (vertices.count(child->vertex(j)) > 0) ++count;
    }
    if (count > 1) return false;
  }
  return true;
}

bool S2Polygon::BoundaryEquals(S2Polygon const* b) const {
  if (num_loops() != b->num_loops()) return false;

  for (int i = 0; i < num_loops(); ++i) {
    S2Loop const* a_loop = loop(i);
    bool success = false;
    for (int j = 0; j < num_loops(); ++j) {
      S2Loop const* b_loop = b->loop(j);
      if ((b_loop->depth() == a_loop->depth()) &&
          b_loop->BoundaryEquals(a_loop)) {
        success = true;
        break;
      }
    }
    if (!success) return false;
  }
  return true;
}

bool S2Polygon::BoundaryApproxEquals(S2Polygon const* b,
                                     double max_error) const {
  if (num_loops() != b->num_loops()) return false;

  // For now, we assume that there is at most one candidate match for each
  // loop.  (So far this method is just used for testing.)

  for (int i = 0; i < num_loops(); ++i) {
    S2Loop const* a_loop = loop(i);
    bool success = false;
    for (int j = 0; j < num_loops(); ++j) {
      S2Loop const* b_loop = b->loop(j);
      if (b_loop->depth() == a_loop->depth() &&
          b_loop->BoundaryApproxEquals(a_loop, max_error)) {
        success = true;
        break;
      }
    }
    if (!success) return false;
  }
  return true;
}

bool S2Polygon::BoundaryNear(S2Polygon const* b, double max_error) const {
  if (num_loops() != b->num_loops()) return false;

  // For now, we assume that there is at most one candidate match for each
  // loop.  (So far this method is just used for testing.)

  for (int i = 0; i < num_loops(); ++i) {
    S2Loop const* a_loop = loop(i);
    bool success = false;
    for (int j = 0; j < num_loops(); ++j) {
      S2Loop const* b_loop = b->loop(j);
      if (b_loop->depth() == a_loop->depth() &&
          b_loop->BoundaryNear(a_loop, max_error)) {
        success = true;
        break;
      }
    }
    if (!success) return false;
  }
  return true;
}

void S2Polygon::EncodeCompressed(Encoder* encoder,
                                 S2XYZFaceSiTi const* all_vertices,
                                 int snap_level) const {
  CHECK_GE(snap_level, 0);
  // Sufficient for what we write. Typically enough for a 4 vertex polygon.
  encoder->Ensure(40);
  encoder->put8(kCurrentCompressedEncodingVersionNumber);
  encoder->put8(snap_level);
  encoder->put_varint32(num_loops());
  DCHECK_GE(encoder->avail(), 0);
  S2XYZFaceSiTi const* current_loop_vertices = all_vertices;
  for (int i = 0; i < num_loops(); ++i) {
    loops_[i]->EncodeCompressed(encoder, current_loop_vertices, snap_level);
    current_loop_vertices += loops_[i]->num_vertices();
  }
  // Do not write the bound, num_vertices, or has_holes_ as they can be
  // cheaply recomputed by DecodeCompressed.  Microbenchmarks show the
  // speed difference is inconsequential.
}

bool S2Polygon::DecodeCompressed(Decoder* decoder) {
  if (decoder->avail() < sizeof(uint8)) return false;
  ClearLoops();
  int snap_level = decoder->get8();
  uint32 num_loops;
  if (!decoder->get_varint32(&num_loops)) {
    return false;
  }
  loops_.reserve(num_loops);
  for (int i = 0; i < num_loops; ++i) {
    S2Loop* loop = new S2Loop;
    loop->set_s2debug_override(s2debug_override());
    loops_.push_back(loop);
    if (!loop->DecodeCompressed(decoder, snap_level)) {
      return false;
    }
  }
  InitLoopProperties();
  return true;
}
