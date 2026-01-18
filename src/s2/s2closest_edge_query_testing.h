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

// Common code for testing furthest/closest edge/point queries.

#ifndef S2_S2CLOSEST_EDGE_QUERY_TESTING_H_
#define S2_S2CLOSEST_EDGE_QUERY_TESTING_H_

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "absl/base/attributes.h"
#include "absl/log/absl_check.h"
#include "absl/random/bit_gen_ref.h"
#include "absl/random/random.h"
#include "s2/_fp_contract_off.h"  // IWYU pragma: keep
#include "s2/mutable_s2shape_index.h"
#include "s2/s1angle.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2edge_distances.h"
#include "s2/s2edge_vector_shape.h"
#include "s2/s2fractal.h"
#include "s2/s2loop.h"
#include "s2/s2metrics.h"
#include "s2/s2point.h"
#include "s2/s2point_vector_shape.h"
#include "s2/s2random.h"
#include "s2/s2shape.h"
#include "s2/s2shapeutil_count_edges.h"
#include "s2/s2shapeutil_shape_edge_id.h"
#include "s2/s2testing.h"

namespace s2testing {

// An abstract class that adds edges to a MutableS2ShapeIndex for benchmarking.
class ShapeIndexFactory {
 public:
  virtual ~ShapeIndexFactory() = default;

  // Requests that approximately "num_edges" edges located within the given
  // S2Cap bound should be added to "index".
  virtual void AddEdges(const S2Cap& index_cap, int num_edges,
                        MutableS2ShapeIndex* index) const = 0;
};

// Generates a regular loop that approximately fills the given S2Cap.
//
// Regular loops are nearly the worst case for distance calculations, since
// many edges are nearly equidistant from any query point that is not
// immediately adjacent to the loop.
class RegularLoopShapeIndexFactory : public ShapeIndexFactory {
 public:
  RegularLoopShapeIndexFactory() = default;
  // We present the same interface as the other factories, even though we
  // do not use the BitGenRef here.
  explicit RegularLoopShapeIndexFactory(
      absl::BitGenRef unused_bitgen ABSL_ATTRIBUTE_LIFETIME_BOUND) {}

  void AddEdges(const S2Cap& index_cap, int num_edges,
                MutableS2ShapeIndex* index) const override {
    index->Add(std::make_unique<S2Loop::OwningShape>(S2Loop::MakeRegularLoop(
        index_cap.center(), index_cap.GetRadius(), num_edges)));
  }
};

// Generates a fractal loop that approximately fills the given S2Cap.
class FractalLoopShapeIndexFactory : public ShapeIndexFactory {
 public:
  explicit FractalLoopShapeIndexFactory(
      absl::BitGenRef bitgen ABSL_ATTRIBUTE_LIFETIME_BOUND)
      : bitgen_(bitgen) {}

  void AddEdges(const S2Cap& index_cap, int num_edges,
                MutableS2ShapeIndex* index) const override {
    S2Fractal fractal(bitgen_);
    fractal.SetLevelForApproxMaxEdges(num_edges);
    index->Add(std::make_unique<S2Loop::OwningShape>(
        fractal.MakeLoop(s2random::FrameAt(bitgen_, index_cap.center()),
                         index_cap.GetRadius())));
  }

 private:
  absl::BitGenRef bitgen_;
};

// Generates a cloud of points that approximately fills the given S2Cap.
class PointCloudShapeIndexFactory : public ShapeIndexFactory {
 public:
  explicit PointCloudShapeIndexFactory(
      absl::BitGenRef bitgen ABSL_ATTRIBUTE_LIFETIME_BOUND)
      : bitgen_(bitgen) {}

  void AddEdges(const S2Cap& index_cap, int num_edges,
                MutableS2ShapeIndex* index) const override {
    std::vector<S2Point> points;
    points.reserve(num_edges);
    for (int i = 0; i < num_edges; ++i) {
      points.push_back(s2random::SamplePoint(bitgen_, index_cap));
    }
    index->Add(std::make_unique<S2PointVectorShape>(std::move(points)));
  }

 private:
  absl::BitGenRef bitgen_;
};

enum class TargetType { POINT, EDGE, CELL, INDEX };

static S2Shape::Edge SampleEdge(absl::BitGenRef bitgen,
                                const MutableS2ShapeIndex& index) {
  int e = absl::Uniform(bitgen, 0, s2shapeutil::CountEdges(index));
  for (int s = 0; ; ++s) {
    const S2Shape* shape = index.shape(s);
    if (shape == nullptr) continue;
    if (e < shape->num_edges()) {
      return shape->edge(e);
    }
    e -= shape->num_edges();
  }
}

static S2CellId SampleCell(absl::BitGenRef bitgen,
                           const MutableS2ShapeIndex& index) {
  int num_cells = 0;
  MutableS2ShapeIndex::Iterator iter(&index);
  for (; !iter.done(); iter.Next()) ++num_cells;
  iter.Begin();
  for (int i = absl::Uniform(bitgen, 0, num_cells); --i >= 0; iter.Next())
    continue;
  return iter.id();
}

static S2Point SampleBoundary(absl::BitGenRef bitgen, const S2Cap& cap) {
  return S2::GetPointOnLine(cap.center(), s2random::Point(bitgen),
                            cap.GetRadius());
}

static S1Angle FractionToRadius(absl::BitGenRef bitgen, double fraction,
                                double radius_km) {
  if (fraction < 0) {
    fraction = -absl::Uniform(bitgen, 0.0, 1.0) * fraction;
  }
  return fraction * S2Testing::KmToAngle(radius_km);
}

// Generates and adds geometry to a MutableS2ShapeIndex for use in an edge
// query, via its one public method GenerateEdgeQueryWithTargets.
//
// Approximately "num_index_edges" will be generated by "factory" and
// inserted.  The geometry is generated within an S2Cap of the radius specified
// by "radius_km" (the "index radius").  Parameters with "fraction" in their
// names are expressed as a fraction of this radius.
//
// Also generates a set of target geometries for the query, based on the
// "target_type" and the input parameters.  If "target_type" is INDEX, then:
//   (i) the target will have approximately "num_target_edges" edges.
//   (ii) "include_interiors" will be set on the target index.
//
//   - If "choose_target_from_index" is true, then the target will be chosen
//     from the geometry in the index itself, otherwise it will be chosen
//     randomly according to the parameters below:
//
//   - If target_radius_fraction > 0, the target radius will be approximately
//     the given fraction of the index radius; if target_radius_fraction < 0,
//     it will be chosen randomly up to corresponding positive fraction.
//
//   - If center_separation_fraction > 0, then the centers of index and target
//     bounding caps will be separated by the given fraction of the index
//     radius; if center_separation_fraction < 0, they will be separated by up
//     to the corresponding positive fraction.
//
//   - The "bitgen" is used as a source of randomness.  It must outlive
//     this instance.
//
class S2ClosestEdgeQueryBenchmarkHarness {
 public:
  S2ClosestEdgeQueryBenchmarkHarness(
      const ShapeIndexFactory& factory, int num_index_edges,
      bool include_interiors, TargetType target_type, int num_target_edges,
      bool choose_target_from_index, double radius_km,
      double target_radius_fraction, double center_separation_fraction,
      absl::BitGenRef bitgen ABSL_ATTRIBUTE_LIFETIME_BOUND)
      : factory_(factory),
        num_index_edges_(num_index_edges),
        include_interiors_(include_interiors),
        target_type_(target_type),
        num_target_edges_(num_target_edges),
        choose_target_from_index_(choose_target_from_index),
        radius_km_(radius_km),
        target_radius_fraction_(target_radius_fraction),
        center_separation_fraction_(center_separation_fraction),
        bitgen_(bitgen) {}

  template <class EdgeQueryType>
  void GenerateEdgeQueryWithTargets(
      EdgeQueryType* query,
      MutableS2ShapeIndex* query_index,
      std::vector<std::unique_ptr<typename EdgeQueryType::Target>>* targets,
      std::vector<std::unique_ptr<MutableS2ShapeIndex>>* target_indexes) {
    ABSL_CHECK_EQ(query_index, &query->index());

    // To save time, we generate at most this many distinct targets per index.
    static const int kMaxTargetsPerIndex = 100;

    query_index->Clear();
    S2Cap index_cap(s2random::Point(bitgen_), S2Testing::KmToAngle(radius_km_));
    factory_.AddEdges(index_cap, num_index_edges_, query_index);
    query_index->ForceBuild();
    query->ReInit();
    targets->clear();
    target_indexes->clear();

    int num_targets = kMaxTargetsPerIndex;
    if (target_type_ == TargetType::INDEX) {
      // Limit the total number of target edges to reduce the benchmark
      // running times.
      num_targets = std::min(num_targets, 500000 / num_target_edges_);
    }
    for (int i = 0; i < num_targets; ++i) {
      S1Angle target_dist =
          FractionToRadius(bitgen_, center_separation_fraction_, radius_km_);
      S2Point boundary_point =
          SampleBoundary(bitgen_, S2Cap(index_cap.center(), target_dist));
      S2Cap target_cap(
          boundary_point,
          FractionToRadius(bitgen_, target_radius_fraction_, radius_km_));

      if (target_type_ == TargetType::POINT) {
        S2Point v0;
        if (choose_target_from_index_) {
          v0 = SampleEdge(bitgen_, *query_index).v0;
        } else {
          v0 = target_cap.center();
        }
        targets->push_back(
            std::make_unique<typename EdgeQueryType::PointTarget>(v0));

      } else if (target_type_ == TargetType::EDGE) {
        S2Point v0, v1;
        if (choose_target_from_index_) {
          auto edge = SampleEdge(bitgen_, *query_index);
          v0 = edge.v0;
          v1 = edge.v1;
        } else {
          v0 = SampleBoundary(bitgen_, target_cap);
          v1 = SampleBoundary(bitgen_, target_cap);
        }
        targets->push_back(
            std::make_unique<typename EdgeQueryType::EdgeTarget>(v0, v1));

      } else if (target_type_ == TargetType::CELL) {
        S2CellId cellid;
        if (choose_target_from_index_) {
          cellid = SampleCell(bitgen_, *query_index);
        } else {
          cellid = S2CellId(target_cap.center()).parent(
              S2::kMaxDiag.GetClosestLevel(target_cap.GetRadius().radians()));
        }
        targets->push_back(std::make_unique<typename EdgeQueryType::CellTarget>(
            S2Cell(cellid)));

      } else {
        ABSL_DCHECK(target_type_ == TargetType::INDEX);
        auto target_index = std::make_unique<MutableS2ShapeIndex>();
        if (choose_target_from_index_) {
          auto shape = std::make_unique<S2EdgeVectorShape>();
          for (int i = 0; i < num_target_edges_; ++i) {
            auto edge = SampleEdge(bitgen_, *query_index);
            shape->Add(edge.v0, edge.v1);
          }
          target_index->Add(std::move(shape));
        } else {
          factory_.AddEdges(target_cap, num_target_edges_, &*target_index);
        }
        target_index->ForceBuild();
        auto target =
            std::make_unique<typename EdgeQueryType::ShapeIndexTarget>(
                &*target_index);
        target->set_include_interiors(include_interiors_);
        targets->push_back(std::move(target));
        target_indexes->push_back(std::move(target_index));
      }
    }
  }

 private:
  const ShapeIndexFactory& factory_;
  int num_index_edges_;
  bool include_interiors_;
  TargetType target_type_;
  int num_target_edges_;
  bool choose_target_from_index_;
  double radius_km_;
  double target_radius_fraction_;
  double center_separation_fraction_;
  absl::BitGenRef bitgen_;
};
}  // namespace s2testing
#endif  // S2_S2CLOSEST_EDGE_QUERY_TESTING_H_
