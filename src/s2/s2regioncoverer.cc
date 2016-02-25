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

#include "s2/s2regioncoverer.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <functional>
#include <queue>
#include <unordered_set>
#include <vector>

#include <glog/logging.h>
#include "s2/s1angle.h"
#include "s2/s2.h"
#include "s2/s2cap.h"
#include "s2/s2cellunion.h"
#include "s2/s2region.h"
#include <algorithm>

using std::max;
using std::min;
using std::unordered_set;
using std::vector;

// Define storage for header file constants (the values are not needed here).
int const S2RegionCoverer::kDefaultMaxCells;

S2RegionCoverer::S2RegionCoverer() :
  min_level_(0),
  max_level_(S2CellId::kMaxLevel),
  level_mod_(1),
  max_cells_(kDefaultMaxCells),
  region_(nullptr) {
}

S2RegionCoverer::~S2RegionCoverer() {
  // Prevent inline destructor bloat by providing a definition.
}

void S2RegionCoverer::set_min_level(int min_level) {
  DCHECK_GE(min_level, 0);
  DCHECK_LE(min_level, S2CellId::kMaxLevel);
  min_level_ = max(0, min(S2CellId::kMaxLevel, min_level));
}

void S2RegionCoverer::set_max_level(int max_level) {
  DCHECK_GE(max_level, 0);
  DCHECK_LE(max_level, S2CellId::kMaxLevel);
  max_level_ = max(0, min(S2CellId::kMaxLevel, max_level));
}

void S2RegionCoverer::set_fixed_level(int level) {
  set_min_level(level);
  set_max_level(level);
}

void S2RegionCoverer::set_level_mod(int level_mod) {
  DCHECK_GE(level_mod, 1);
  DCHECK_LE(level_mod, 3);
  level_mod_ = max(1, min(3, level_mod));
}

void S2RegionCoverer::set_max_cells(int max_cells) {
  max_cells_ = max_cells;
}

S2RegionCoverer::Candidate* S2RegionCoverer::NewCandidate(S2Cell const& cell) {
  if (!region_->MayIntersect(cell)) return nullptr;

  bool is_terminal = false;
  if (cell.level() >= min_level_) {
    if (interior_covering_) {
      if (region_->Contains(cell)) {
        is_terminal = true;
      } else if (cell.level() + level_mod_ > max_level_) {
        return nullptr;
      }
    } else {
      if (cell.level() + level_mod_ > max_level_ || region_->Contains(cell)) {
        is_terminal = true;
      }
    }
  }
  size_t children_size = 0;
  if (!is_terminal) {
    children_size = sizeof(Candidate*) << max_children_shift();
  }
  Candidate* candidate =
      static_cast<Candidate*>(malloc(sizeof(Candidate) + children_size));
  candidate->cell = cell;
  candidate->is_terminal = is_terminal;
  candidate->num_children = 0;
  if (!is_terminal) {
    std::fill_n(&candidate->children[0], 1 << max_children_shift(),
                implicit_cast<Candidate*>(nullptr));
  }
  ++candidates_created_counter_;
  return candidate;
}

void S2RegionCoverer::DeleteCandidate(Candidate* candidate,
                                      bool delete_children) {
  if (delete_children) {
    for (int i = 0; i < candidate->num_children; ++i)
      DeleteCandidate(candidate->children[i], true);
  }
  free(candidate);
}

int S2RegionCoverer::ExpandChildren(Candidate* candidate,
                                    S2Cell const& cell, int num_levels) {
  num_levels--;
  S2Cell child_cells[4];
  cell.Subdivide(child_cells);
  int num_terminals = 0;
  for (int i = 0; i < 4; ++i) {
    if (num_levels > 0) {
      if (region_->MayIntersect(child_cells[i])) {
        num_terminals += ExpandChildren(candidate, child_cells[i], num_levels);
      }
      continue;
    }
    Candidate* child = NewCandidate(child_cells[i]);
    if (child) {
      candidate->children[candidate->num_children++] = child;
      if (child->is_terminal) ++num_terminals;
    }
  }
  return num_terminals;
}

void S2RegionCoverer::AddCandidate(Candidate* candidate) {
  if (candidate == nullptr) return;

  if (candidate->is_terminal) {
    result_.push_back(candidate->cell.id());
    DeleteCandidate(candidate, true);
    return;
  }
  DCHECK_EQ(0, candidate->num_children);

  // Expand one level at a time until we hit min_level_ to ensure that
  // we don't skip over it.
  int num_levels = (candidate->cell.level() < min_level_) ? 1 : level_mod_;
  int num_terminals = ExpandChildren(candidate, candidate->cell, num_levels);

  if (candidate->num_children == 0) {
    DeleteCandidate(candidate, false);

  } else if (!interior_covering_ &&
             num_terminals == 1 << max_children_shift() &&
             candidate->cell.level() >= min_level_) {
    // Optimization: add the parent cell rather than all of its children.
    // We can't do this for interior coverings, since the children just
    // intersect the region, but may not be contained by it - we need to
    // subdivide them further.
    candidate->is_terminal = true;
    AddCandidate(candidate);

  } else {
    // We negate the priority so that smaller absolute priorities are returned
    // first.  The heuristic is designed to refine the largest cells first,
    // since those are where we have the largest potential gain.  Among cells
    // of the same size, we prefer the cells with the fewest children.
    // Finally, among cells with equal numbers of children we prefer those
    // with the smallest number of children that cannot be refined further.
    int priority = -((((candidate->cell.level() << max_children_shift())
                       + candidate->num_children) << max_children_shift())
                     + num_terminals);
    pq_.push(std::make_pair(priority, candidate));
    VLOG(2) << "Push: " << candidate->cell.id() << " (" << priority << ") ";
  }
}

inline int S2RegionCoverer::AdjustLevel(int level) const {
  if (level_mod() > 1 && level > min_level()) {
    level -= (level - min_level()) % level_mod();
  }
  return level;
}

void S2RegionCoverer::AdjustCellLevels(vector<S2CellId>* cells) const {
  DCHECK(std::is_sorted(cells->begin(), cells->end()));
  if (level_mod() == 1) return;

  int out = 0;
  for (S2CellId id : *cells) {
    int level = id.level();
    int new_level = AdjustLevel(level);
    if (new_level != level) id = id.parent(new_level);
    if (out > 0 && (*cells)[out-1].contains(id)) continue;
    while (out > 0 && id.contains((*cells)[out-1])) --out;
    (*cells)[out++] = id;
  }
  cells->resize(out);
}

void S2RegionCoverer::GetInitialCandidates() {
  // Optimization: start with a small (usually 4 cell) covering of the
  // region's bounding cap.
  S2RegionCoverer tmp_coverer;
  tmp_coverer.set_max_cells(min(4, max_cells()));
  tmp_coverer.set_max_level(max_level());
  vector<S2CellId> cells;
  tmp_coverer.GetFastCovering(region_->GetCapBound(), &cells);
  AdjustCellLevels(&cells);
  for (S2CellId cell_id : cells) {
    AddCandidate(NewCandidate(S2Cell(cell_id)));
  }
}

void S2RegionCoverer::GetCoveringInternal(S2Region const& region) {
  // Strategy: Start with the 6 faces of the cube.  Discard any
  // that do not intersect the shape.  Then repeatedly choose the
  // largest cell that intersects the shape and subdivide it.
  //
  // result_ contains the cells that will be part of the output, while pq_
  // contains cells that we may still subdivide further.  Cells that are
  // entirely contained within the region are immediately added to the output,
  // while cells that do not intersect the region are immediately discarded.
  // Therefore pq_ only contains cells that partially intersect the region.
  // Candidates are prioritized first according to cell size (larger cells
  // first), then by the number of intersecting children they have (fewest
  // children first), and then by the number of fully contained children
  // (fewest children first).

  DCHECK(pq_.empty());
  DCHECK(result_.empty());
  region_ = &region;
  candidates_created_counter_ = 0;

  GetInitialCandidates();
  while (!pq_.empty() &&
         (!interior_covering_ || result_.size() < max_cells_)) {
    Candidate* candidate = pq_.top().second;
    pq_.pop();
    VLOG(2) << "Pop: " << candidate->cell.id();
    // For interior covering we keep subdividing no matter how many children
    // candidate has. If we reach max_cells_ before expanding all children,
    // we will just use some of them.
    // For exterior covering we cannot do this, because result has to cover the
    // whole region, so all children have to be used.
    // candidate->num_children == 1 case takes care of the situation when we
    // already have more then max_cells_ in relust (min_level is too high).
    // Subdividing of the candidate with one child does no harm in this case.
    if (interior_covering_ ||
        candidate->cell.level() < min_level_ ||
        candidate->num_children == 1 ||
        result_.size() + pq_.size() + candidate->num_children <= max_cells_) {
      // Expand this candidate into its children.
      for (int i = 0; i < candidate->num_children; ++i) {
        if (interior_covering_ && result_.size() >= max_cells_) {
          DeleteCandidate(candidate->children[i], true);
        } else {
          AddCandidate(candidate->children[i]);
        }
      }
      DeleteCandidate(candidate, false);
    } else {
      candidate->is_terminal = true;
      AddCandidate(candidate);
    }
  }
  VLOG(2) << "Created " << result_.size() << " cells, " <<
      candidates_created_counter_ << " candidates created, " <<
      pq_.size() << " left";
  while (!pq_.empty()) {
    DeleteCandidate(pq_.top().second, true);
    pq_.pop();
  }
  region_ = nullptr;
}

void S2RegionCoverer::GetCovering(S2Region const& region,
                                  vector<S2CellId>* covering) {
  // Rather than just returning the raw list of cell ids generated by
  // GetCoveringInternal(), we construct a cell union and then denormalize it.
  // This has the effect of replacing four child cells with their parent
  // whenever this does not violate the covering parameters specified
  // (min_level, level_mod, etc).  This strategy significantly reduces the
  // number of cells returned in many cases, and it is cheap compared to
  // computing the covering in the first place.

  S2CellUnion tmp;
  GetCellUnion(region, &tmp);
  tmp.Denormalize(min_level(), level_mod(), covering);
}

void S2RegionCoverer::GetInteriorCovering(S2Region const& region,
                                          vector<S2CellId>* interior) {
  S2CellUnion tmp;
  GetInteriorCellUnion(region, &tmp);
  tmp.Denormalize(min_level(), level_mod(), interior);
}

void S2RegionCoverer::GetCellUnion(S2Region const& region,
                                   S2CellUnion* covering) {
  interior_covering_ = false;
  GetCoveringInternal(region);
  covering->InitSwap(&result_);
}

void S2RegionCoverer::GetInteriorCellUnion(S2Region const& region,
                                           S2CellUnion* interior) {
  interior_covering_ = true;
  GetCoveringInternal(region);
  interior->InitSwap(&result_);
}

void S2RegionCoverer::GetFastCovering(S2Cap const& cap,
                                      vector<S2CellId>* covering) {
  GetRawFastCovering(cap, max_cells(), covering);
  NormalizeCovering(covering);
}

void S2RegionCoverer::GetRawFastCovering(S2Cap const& cap,
                                         int max_cells_hint,
                                         vector<S2CellId>* covering) {
  // TODO(ericv): The covering could be made quite a bit tighter by mapping
  // the cap to a rectangle in (i,j)-space and finding a covering for that.
  covering->clear();

  // Find the maximum level such that the cap contains at most one cell vertex
  // and such that S2CellId::AppendVertexNeighbors() can be called.
  int level = S2::kMinWidth.GetMaxLevel(2 * cap.GetRadius().radians());
  level = min(level, S2CellId::kMaxLevel - 1);

  // Don't bother trying to optimize the level == 0 case, since more than
  // four face cells may be required.
  if (level == 0) {
    covering->reserve(6);
    for (int face = 0; face < 6; ++face) {
      covering->push_back(S2CellId::FromFace(face));
    }
  } else {
    // The covering consists of the 4 cells at the given level that share the
    // cell vertex that is closest to the cap center.
    covering->reserve(4);
    S2CellId id = S2CellId::FromPoint(cap.center());
    id.AppendVertexNeighbors(level, covering);
  }
}

void S2RegionCoverer::NormalizeCovering(vector<S2CellId>* covering) {
  // This method makes no attempt to be optimal.  In particular, if
  // min_level() > 0 or level_mod() > 1 then it may return more than the
  // desired number of cells even when this isn't necessary.
  //
  // Note that when the covering parameters have their default values, almost
  // all of the code in this function is skipped.

  // If any cells are too small, or don't satisfy level_mod(), then replace
  // them with ancestors.
  if (max_level() < S2CellId::kMaxLevel || level_mod() > 1) {
    for (int i = 0; i < covering->size(); ++i) {
      S2CellId id = (*covering)[i];
      int level = id.level();
      int new_level = AdjustLevel(min(level, max_level()));
      if (new_level != level) {
        (*covering)[i] = id.parent(new_level);
      }
    }
  }
  // Sort the cells and simplify them.
  S2CellUnion::Normalize(covering);

  // If there are still too many cells, then repeatedly replace two adjacent
  // cells in S2CellId order by their lowest common ancestor.
  while (covering->size() > max_cells()) {
    int best_index = -1, best_level = -1;
    for (int i = 0; i + 1 < covering->size(); ++i) {
      int level = (*covering)[i].GetCommonAncestorLevel((*covering)[i+1]);
      level = AdjustLevel(level);
      if (level > best_level) {
        best_level = level;
        best_index = i;
      }
    }
    if (best_level < min_level()) break;
    (*covering)[best_index] = (*covering)[best_index].parent(best_level);
    S2CellUnion::Normalize(covering);
  }
  // Make sure that the covering satisfies min_level() and level_mod(),
  // possibly at the expense of satisfying max_cells().
  if (min_level() > 0 || level_mod() > 1) {
    S2CellUnion result;
    result.InitRawSwap(covering);
    result.Denormalize(min_level(), level_mod(), covering);
  }
}

void S2RegionCoverer::FloodFill(S2Region const& region, S2CellId start,
                                vector<S2CellId>* output) {
  unordered_set<S2CellId, S2CellIdHash> all;
  vector<S2CellId> frontier;
  output->clear();
  all.insert(start);
  frontier.push_back(start);
  while (!frontier.empty()) {
    S2CellId id = frontier.back();
    frontier.pop_back();
    if (!region.MayIntersect(S2Cell(id))) continue;
    output->push_back(id);

    S2CellId neighbors[4];
    id.GetEdgeNeighbors(neighbors);
    for (int edge = 0; edge < 4; ++edge) {
      S2CellId nbr = neighbors[edge];
      if (all.insert(nbr).second) {
        frontier.push_back(nbr);
      }
    }
  }
}

void S2RegionCoverer::GetSimpleCovering(
    S2Region const& region, S2Point const& start,
    int level, vector<S2CellId>* output) {
  return FloodFill(region, S2CellId::FromPoint(start).parent(level), output);
}
