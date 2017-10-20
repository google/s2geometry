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
#include "s2/s2cap.h"
#include "s2/s2cellunion.h"
#include "s2/s2metrics.h"
#include "s2/s2region.h"

using std::is_sorted;
using std::max;
using std::min;
using std::unordered_set;
using std::vector;

S2RegionCoverer::S2RegionCoverer(S2RegionCoverer::Options const& options) :
  options_(options) {
}

// Defaulted in the implementation to prevent inline bloat.
S2RegionCoverer::S2RegionCoverer() = default;
S2RegionCoverer::~S2RegionCoverer() = default;
S2RegionCoverer::S2RegionCoverer(S2RegionCoverer&&) = default;
S2RegionCoverer& S2RegionCoverer::operator=(S2RegionCoverer&&) = default;

void S2RegionCoverer::Options::set_min_level(int min_level) {
  DCHECK_GE(min_level, 0);
  DCHECK_LE(min_level, S2CellId::kMaxLevel);
  min_level_ = max(0, min(S2CellId::kMaxLevel, min_level));
}

void S2RegionCoverer::Options::set_max_level(int max_level) {
  DCHECK_GE(max_level, 0);
  DCHECK_LE(max_level, S2CellId::kMaxLevel);
  max_level_ = max(0, min(S2CellId::kMaxLevel, max_level));
}

void S2RegionCoverer::Options::set_fixed_level(int level) {
  set_min_level(level);
  set_max_level(level);
}

void S2RegionCoverer::Options::set_level_mod(int level_mod) {
  DCHECK_GE(level_mod, 1);
  DCHECK_LE(level_mod, 3);
  level_mod_ = max(1, min(3, level_mod));
}

void S2RegionCoverer::Options::set_max_cells(int max_cells) {
  max_cells_ = max_cells;
}

S2RegionCoverer::Candidate* S2RegionCoverer::NewCandidate(const S2Cell& cell) {
  if (!region_->MayIntersect(cell)) return nullptr;

  bool is_terminal = false;
  if (cell.level() >= options_.min_level()) {
    if (interior_covering_) {
      if (region_->Contains(cell)) {
        is_terminal = true;
      } else if (cell.level() + options_.level_mod() > options_.max_level()) {
        return nullptr;
      }
    } else {
      if (cell.level() + options_.level_mod() > options_.max_level() ||
          region_->Contains(cell)) {
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
                                    const S2Cell& cell, int num_levels) {
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

  // Expand one level at a time until we hit min_level() to ensure that we
  // don't skip over it.
  int num_levels = ((candidate->cell.level() < options_.min_level()) ?
                    1 : options_.level_mod());
  int num_terminals = ExpandChildren(candidate, candidate->cell, num_levels);

  if (candidate->num_children == 0) {
    DeleteCandidate(candidate, false);

  } else if (!interior_covering_ &&
             num_terminals == 1 << max_children_shift() &&
             candidate->cell.level() >= options_.min_level()) {
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
  if (options_.level_mod() > 1 && level > options_.min_level()) {
    level -= (level - options_.min_level()) % options_.level_mod();
  }
  return level;
}

void S2RegionCoverer::AdjustCellLevels(vector<S2CellId>* cells) const {
  DCHECK(is_sorted(cells->begin(), cells->end()));
  if (options_.level_mod() == 1) return;

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
  tmp_coverer.mutable_options()->set_max_cells(min(4, options_.max_cells()));
  tmp_coverer.mutable_options()->set_max_level(options_.max_level());
  vector<S2CellId> cells;
  tmp_coverer.GetFastCovering(*region_, &cells);
  AdjustCellLevels(&cells);
  for (S2CellId cell_id : cells) {
    AddCandidate(NewCandidate(S2Cell(cell_id)));
  }
}

void S2RegionCoverer::GetCoveringInternal(const S2Region& region) {
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
         (!interior_covering_ || result_.size() < options_.max_cells())) {
    Candidate* candidate = pq_.top().second;
    pq_.pop();
    VLOG(2) << "Pop: " << candidate->cell.id();
    // For interior coverings we keep subdividing no matter how many children
    // the candidate has.  If we reach max_cells() before expanding all
    // children, we will just use some of them.  For exterior coverings we
    // cannot do this, because the result has to cover the whole region, so
    // all children have to be used.  The (candidate->num_children == 1) case
    // takes care of the situation when we already have more than max_cells()
    // in results (min_level is too high).  Subdividing the candidate with one
    // child does no harm in this case.
    if (interior_covering_ ||
        candidate->cell.level() < options_.min_level() ||
        candidate->num_children == 1 ||
        (result_.size() + pq_.size() + candidate->num_children <=
         options_.max_cells())) {
      // Expand this candidate into its children.
      for (int i = 0; i < candidate->num_children; ++i) {
        if (interior_covering_ && result_.size() >= options_.max_cells()) {
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

  // Rather than just returning the raw list of cell ids, we construct a cell
  // union and then denormalize it.  This has the effect of replacing four
  // child cells with their parent whenever this does not violate the covering
  // parameters specified (min_level, level_mod, etc).  This significantly
  // reduces the number of cells returned in many cases, and it is cheap
  // compared to computing the covering in the first place.
  S2CellUnion::Normalize(&result_);
  if (options_.min_level() > 0 || options_.level_mod() > 1) {
    auto result_copy = result_;
    S2CellUnion::Denormalize(result_copy, options_.min_level(),
                             options_.level_mod(), &result_);
  }
}

void S2RegionCoverer::GetCovering(const S2Region& region,
                                  vector<S2CellId>* covering) {
  interior_covering_ = false;
  GetCoveringInternal(region);
  *covering = std::move(result_);
}

void S2RegionCoverer::GetInteriorCovering(const S2Region& region,
                                          vector<S2CellId>* interior) {
  interior_covering_ = true;
  GetCoveringInternal(region);
  *interior = std::move(result_);
}

S2CellUnion S2RegionCoverer::GetCovering(const S2Region& region) {
  interior_covering_ = false;
  GetCoveringInternal(region);
  return S2CellUnion::FromVerbatim(std::move(result_));
}

S2CellUnion S2RegionCoverer::GetInteriorCovering(const S2Region& region) {
  interior_covering_ = true;
  GetCoveringInternal(region);
  return S2CellUnion::FromVerbatim(std::move(result_));
}

void S2RegionCoverer::GetFastCovering(const S2Region& region,
                                      vector<S2CellId>* covering) {
  region.GetCellUnionBound(covering);
  NormalizeCovering(covering);
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
  if (options_.max_level() < S2CellId::kMaxLevel || options_.level_mod() > 1) {
    for (int i = 0; i < covering->size(); ++i) {
      S2CellId id = (*covering)[i];
      int level = id.level();
      int new_level = AdjustLevel(min(level, options_.max_level()));
      if (new_level != level) {
        (*covering)[i] = id.parent(new_level);
      }
    }
  }
  // Sort the cells and simplify them.
  S2CellUnion::Normalize(covering);

  // If there are still too many cells, then repeatedly replace two adjacent
  // cells in S2CellId order by their lowest common ancestor.
  while (covering->size() > options_.max_cells()) {
    int best_index = -1, best_level = -1;
    for (int i = 0; i + 1 < covering->size(); ++i) {
      int level = (*covering)[i].GetCommonAncestorLevel((*covering)[i+1]);
      level = AdjustLevel(level);
      if (level > best_level) {
        best_level = level;
        best_index = i;
      }
    }
    if (best_level < options_.min_level()) break;
    (*covering)[best_index] = (*covering)[best_index].parent(best_level);
    S2CellUnion::Normalize(covering);
  }
  // Make sure that the covering satisfies min_level() and level_mod(),
  // possibly at the expense of satisfying max_cells().
  if (options_.min_level() > 0 || options_.level_mod() > 1) {
    S2CellUnion::Denormalize(*covering, options_.min_level(),
                             options_.level_mod(), &result_);
    *covering = result_;
  }
}

void S2RegionCoverer::FloodFill(const S2Region& region, S2CellId start,
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
    const S2Region& region, const S2Point& start,
    int level, vector<S2CellId>* output) {
  return FloodFill(region, S2CellId(start).parent(level), output);
}
