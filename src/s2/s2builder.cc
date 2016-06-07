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
//
// The algorithm is based on the idea of choosing a set of sites and computing
// their "limited radius Voronoi diagram", which is obtained by intersecting
// each Voronoi region with a disc of fixed radius (the "snap radius")
// centered around the corresponding site.
//
// For each input edge, we then determine the sequence of Voronoi regions
// crossed by that edge, and snap the edge to the corresponding sequence of
// sites.  (In other words, each input each is replaced by an edge chain.)
//
// The sites are chosen by starting with the set of input vertices, optionally
// snapping them to discrete point set (such as S2CellId centers or lat/lng E7
// coordinates), and then choosing a subset such that no two sites are closer
// than the given "snap_radius".  Note that the sites do not need to be spaced
// regularly -- their positions are completely arbitrary.
//
// Rather than computing the full limited radius Voronoi diagram, instead we
// compute on demand the sequence of Voronoi regions intersected by each edge.
// We do this by first finding all the sites that are within "snap_radius" of
// the edge, sorting them by distance from the edge origin, and then using an
// incremental algorithm.
//
// We implement the minimum edge-vertex separation property by snapping all
// the input edges and checking whether any site (the "site to avoid") would
// then be too close to the edge.  If so we add another site (the "separation
// site") along the input edge, positioned so that the new snapped edge is
// guaranteed to be far enough away from the site to avoid.  We then find all
// the input edges that are within "snap_radius" of the new site, and resnap
// those edges.  (It is very rare that we need to add a separation site, even
// when sites are densely spaced.)
//
// Idempotency is implemented by explicitly checking whether the input
// geometry already meets the output criteria.  This is not as sophisticated
// as Stable Snap Rounding (Hershberger); I have worked out the details and it
// is possible to adapt those ideas here, but it would make the implementation
// significantly more complex.
//
// The only way that different output layers interact is in the choice of
// Voronoi sites:
//
//  - Vertices from all layers contribute to the initial selection of sites.
//  - Edges in any layer that pass too close to a site can cause new sites to
//    be added (which affects snapping in all layers).
//  - Simplification can be thought of as removing sites.  A site can be
//    removed only if the snapped edges stay within the error bounds of the
//    corresponding input edges in all layers.
//
// Otherwise all layers are processed independently.  For example, sibling
// edge pairs can only cancel each other within a single layer (if desired).

#include "s2/s2builder.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include "s2/base/casts.h"
#include <glog/logging.h>
#include "s2/id_set_lexicon.h"
#include "s2/s1angle.h"
#include "s2/s1chordangle.h"
#include "s2/s2.h"
#include "s2/s2builder_graph.h"
#include "s2/s2builder_layer.h"
#include "s2/s2builderutil_snap_functions.h"
#include "s2/s2closestedgequery.h"
#include "s2/s2closestpointquery.h"
#include "s2/s2edgeutil.h"
#include "s2/s2error.h"
#include "s2/s2loop.h"
#include "s2/s2pointindex.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"
#include "s2/s2shapeindex.h"
#include "s2/s2shapeutil.h"
#include "s2/s2textformat.h"

using std::make_pair;
using std::max;
using std::pair;
using std::unique_ptr;
using std::vector;

// Internal flag intended to be set from within a debugger.
static bool s2builder_verbose = false;

S1Angle S2Builder::SnapFunction::max_edge_deviation() const {
  // We want max_edge_deviation() to be large enough compared to snap_radius()
  // such that edge splitting is rare.
  //
  // Using spherical trigonometry, if the endpoints of an edge of length L
  // move by at most a distance R, the center of the edge moves by at most
  // asin(sin(R) / cos(L / 2)).  Thus the (max_edge_deviation / snap_radius)
  // ratio increases with both the snap radius R and the edge length L.
  //
  // We arbitrarily limit the edge deviation to be at most 10% more than the
  // snap radius.  With the maximum allowed snap radius of 70 degrees, this
  // means that edges up to 30.6 degrees long are never split.  For smaller
  // snap radii, edges up to 49 degrees long are never split.  (Edges of any
  // length are not split unless their endpoints move far enough so that the
  // actual edge deviation exceeds the limit; in practice, splitting is rare
  // even with long edges.)  Note that it is always possible to split edges
  // when max_edge_deviation() is exceeded; see MaybeAddExtraSites().
  DCHECK_LE(snap_radius(), kMaxSnapRadius());
  double const kMaxEdgeDeviationRatio = 1.1;
  return kMaxEdgeDeviationRatio * snap_radius();
}

S2Builder::Options::Options()
    : snap_function_(unique_ptr<SnapFunction>(
          new s2builderutil::IdentitySnapFunction(S1Angle::Zero()))),
      split_crossing_edges_(false),
      simplify_edge_chains_(false),
      idempotent_(true) {
}

S2Builder::Options::Options(SnapFunction const& snap_function)
    : snap_function_(snap_function.Clone()),
      split_crossing_edges_(false),
      simplify_edge_chains_(false),
      idempotent_(true) {
}

S2Builder::Options::Options(Options const& options)
    :  snap_function_(options.snap_function_->Clone()),
       split_crossing_edges_(options.split_crossing_edges_),
       simplify_edge_chains_(options.simplify_edge_chains_),
       idempotent_(options.idempotent_) {
}

S2Builder::Options& S2Builder::Options::operator=(Options const& options) {
  snap_function_ = options.snap_function_->Clone();
  split_crossing_edges_ = options.split_crossing_edges_;
  simplify_edge_chains_ = options.simplify_edge_chains_;
  idempotent_ = options.idempotent_;
  return *this;
}

static S1ChordAngle GetMaxPointToPointDist(S1Angle a) {
  S1ChordAngle ca(a);
  return ca.PlusError(ca.GetS1AngleConstructorMaxError() +
                      ca.GetS2PointConstructorMaxError());
}

static S1ChordAngle GetMinPointToPointDist(S1Angle a) {
  S1ChordAngle ca(a);
  return ca.PlusError(-(ca.GetS1AngleConstructorMaxError() +
                        ca.GetS2PointConstructorMaxError()));
}

static S1ChordAngle GetMaxPointToEdgeDist(S1Angle a) {
  S1ChordAngle ca(a);
  return ca.PlusError(ca.GetS1AngleConstructorMaxError() +
                      S2EdgeUtil::GetUpdateMinDistanceMaxError(ca));
}

static S1ChordAngle GetMinPointToEdgeDist(S1Angle a) {
  S1ChordAngle ca(a);
  return ca.PlusError(-(ca.GetS1AngleConstructorMaxError() +
                        S2EdgeUtil::GetUpdateMinDistanceMaxError(ca)));
}

S2Builder::S2Builder(Options const& options) {
  Init(options);
}

void S2Builder::Init(Options const& options) {
  options_ = options;
  SnapFunction const& snap_function = options.snap_function();
  snap_radius_ = snap_function.snap_radius();
  edge_snap_radius_ = snap_radius_;
  DCHECK_LE(snap_radius_, SnapFunction::kMaxSnapRadius());

  if (options.split_crossing_edges()) {
    // Increase the snap radius for edges so that both edges are snapped to
    // the edge intersection location.  The problem is that the intersection
    // vertex itself may have been snapped to another vertex that is up to
    // snap_radius_ away.  Since the intersection vertex was not exactly at
    // the intersection point, we need to increase the snap radius slightly so
    // that both edges are guaranteed to be snapped to its new position.
    edge_snap_radius_ += S2EdgeUtil::kIntersectionError;
  }
  snapping_requested_ = (edge_snap_radius_ > S1Angle::Zero());

  // Increase the snap radius to account for errors in converting
  // snap_radius() to an S1ChordAngle and in measuring the S1ChordAngle
  // between two points.
  //
  // This does *not* account for errors in the implementation of
  // SnapPoint().  For example, S2CellId::ToPoint() has a maximum error of
  // 1.5 * DBL_EPSILON.  These errors are handled by having the SnapFunction
  // increase the snap_radius() before passing it to S2Builder.
  site_snap_radius_ca_ = GetMaxPointToPointDist(snap_radius_);

  // Convert the edge snap radius to an S1ChordAngle, and increase it to
  // account for errors.
  edge_snap_radius_ca_ = GetMaxPointToEdgeDist(edge_snap_radius_);

  // Compute the maximum distance that a vertex can be separated from an
  // edge while still affecting how that edge is snapped.
  max_edge_deviation_ = snap_function.max_edge_deviation();
  edge_site_query_radius_ = (max_edge_deviation_ +
                             snap_function.min_edge_vertex_separation());

  // Compute the maximum edge length such that even if both endpoints move by
  // the maximum distance allowed (i.e., snap_radius), the center of the edge
  // will still move by less than max_edge_deviation().  This saves us a lot
  // of work since then we don't need to check the actual deviation.
  min_edge_length_to_split_ca_ = GetMinPointToPointDist(S1Angle::Radians(
      2 * acos(sin(snap_radius_) / sin(max_edge_deviation_))));

  // If the condition below is violated, then AddExtraSites() needs to be
  // modified to check that snapped edges pass on the same side of each "site
  // to avoid" as the input edge.  Currently it doesn't need to do this
  // because the condition below guarantees that if the snapped edge passes on
  // the wrong side of the site then it is also too close, which will cause a
  // separation site to be added.
  //
  // Currently max_edge_deviation() is at most 1.1 * snap_radius(), whereas
  // min_edge_vertex_separation() is at least 0.219 * snap_radius() (based on
  // S2CellIdSnapFunction, which is currently the worst case).
  DCHECK_LE(snap_function.max_edge_deviation(),
            snap_function.snap_radius() +
            snap_function.min_edge_vertex_separation());

  // Adjust the other snap function parameters to account for errors.  Note
  // that min_site_separation_ca_ is adjusted *downwards* since it is used to
  // test whether input vertices too close together, and we want this test to
  // be conservative (in order to implement idempotency).  Similarly,
  // min_edge_site_separation_ca_ is used to test whether we need to add
  // "separation sites" so that edges do not approach non-incident vertices
  // too closely, and we only want to do this when we are sure that the edge
  // and vertex are definitely too close.
  min_site_separation_ = snap_function.min_vertex_separation();
  min_site_separation_ca_ = GetMinPointToPointDist(min_site_separation_);
  min_edge_site_separation_ = snap_function.min_edge_vertex_separation();
  min_edge_site_separation_ca_ = GetMinPointToEdgeDist(
      min_edge_site_separation_);

  // Compute the maximum possible distance between two sites whose Voronoi
  // regions touch.  (The maximum radius of each Voronoi region is
  // edge_snap_radius_.)  Then increase this bound to account for errors.
  max_adjacent_site_separation_ca_ =
      GetMaxPointToPointDist(2 * edge_snap_radius_);

  // Finally, we also precompute sin^2(edge_snap_radius), which is simply the
  // squared distance between a vertex and an edge measured perpendicular to
  // the plane containing the edge, and increase this value by the maximum
  // error in the calculation to compare this distance against the bound.
  double d = sin(edge_snap_radius_);
  edge_snap_radius_sin2_ = d * d;
  edge_snap_radius_sin2_ += ((9.5 * d + 2.5 + 2 * sqrt(3)) * d +
                             9 * DBL_EPSILON) * DBL_EPSILON;

  // Initialize the current label set.
  label_set_id_ = label_set_lexicon_.AddEmptySet();
  label_set_modified_ = false;

  // If snapping was requested, we try to determine whether the input geometry
  // already meets the output requirements.  This is necessary for
  // idempotency, and can also save work.  If we discover any reason that the
  // input geometry needs to be modified, snapping_needed_ is set to true.
  snapping_needed_ = false;
}

void S2Builder::clear_labels() {
  label_set_.clear();
  label_set_modified_ = true;
}

void S2Builder::push_label(Label label) {
  DCHECK_GE(label, 0);
  label_set_.push_back(label);
  label_set_modified_ = true;
}

void S2Builder::pop_label() {
  label_set_.pop_back();
  label_set_modified_ = true;
}

void S2Builder::set_label(Label label) {
  DCHECK_GE(label, 0);
  label_set_.resize(1);
  label_set_[0] = label;
  label_set_modified_ = true;
}

void S2Builder::StartLayer(unique_ptr<Layer> layer) {
  layer_begins_.push_back(input_edges_.size());
  layer_options_.push_back(layer->graph_options());
  layers_.push_back(std::move(layer));
}

// Input vertices are stored in a vector, with some removal of duplicates.
// Edges are represented as (VertexId, VertexId) pairs.  All edges are stored
// in a single vector; each layer corresponds to a contiguous range.

S2Builder::InputVertexId S2Builder::AddVertex(S2Point const& v) {
  // Remove duplicate vertices that follow the pattern AB, BC, CD.  If we want
  // to do anything more sophisticated, either use a ValueLexicon, or sort the
  // vertices once they have all been added, remove duplicates, and update the
  // edges.
  if (input_vertices_.empty() || v != input_vertices_.back()) {
    input_vertices_.push_back(v);
  }
  return input_vertices_.size() - 1;
}

void S2Builder::AddEdge(S2Point const& v0, S2Point const& v1) {
  DCHECK(!layers_.empty()) << "Call StartLayer before adding any edges";

  if (v0 == v1 && (layer_options_.back().degenerate_edges() ==
                   GraphOptions::DegenerateEdges::DISCARD)) {
    return;
  }
  InputVertexId j0 = AddVertex(v0);
  InputVertexId j1 = AddVertex(v1);
  input_edges_.push_back(InputEdge(j0, j1));

  // If there are any labels, then attach them to this input edge.
  if (label_set_modified_) {
    if (label_set_ids_.empty()) {
      // Populate the missing entries with empty label sets.
      label_set_ids_.assign(input_edges_.size() - 1, label_set_id_);
    }
    label_set_id_ = label_set_lexicon_.Add(label_set_);
    label_set_ids_.push_back(label_set_id_);
    label_set_modified_ = false;
  } else if (!label_set_ids_.empty()) {
    label_set_ids_.push_back(label_set_id_);
  }
}

void S2Builder::AddPolyline(S2Polyline const& polyline) {
  const int n = polyline.num_vertices();
  for (int i = 1; i < n; ++i) {
    AddEdge(polyline.vertex(i - 1), polyline.vertex(i));
  }
}

void S2Builder::AddLoop(S2Loop const& loop) {
  // Ignore loops that do not have a boundary.
  if (loop.is_empty_or_full()) return;

  const int n = loop.num_vertices();
  if (loop.sign() > 0 || (layer_options_.back().edge_type() ==
                          EdgeType::UNDIRECTED)) {
    for (int i = 0; i < n; ++i) {
      AddEdge(loop.vertex(i), loop.vertex(i + 1));
    }
  } else {
    // In order to preserve the cyclic order of the loop vertices, we add the
    // edge from vertex n-1 to vertex n-2 first.  This is because these edges
    // will be assembed into a clockwise loop, which will eventually be
    // normalized by S2Polygon by calling S2Loop::Invert().  S2Loop::Invert()
    // reverses the order of the vertices, so if we want to end up with the
    // original vertex order (0, 1, ..., n-1) then we need to build a
    // clockwise loop with vertex order (n-1, n-2, ..., 0).  This is done by
    // adding the edge (n-1, n-2) first, and then ensuring that Build()
    // assembles loops starting from edges in the order they were added.
    for (int i = 2 * n - 1; i >= n; --i) {
      AddEdge(loop.vertex(i), loop.vertex(i - 1));
    }
  }
}

void S2Builder::AddPolygon(S2Polygon const& polygon) {
  for (int i = 0; i < polygon.num_loops(); ++i) {
    AddLoop(*polygon.loop(i));
  }
}

void S2Builder::ForceVertex(S2Point const& vertex) {
  sites_.push_back(vertex);
}

// VertexIdLoopShape is just like LoopShape, except that vertices are
// specified as indices into a vertex vector.  This representation can be more
// compact when many loops are arranged in a mesh structure.
class VertexIdEdgeVectorShape : public S2Shape {
 public:
  // Requires that "edges" is constant for the lifetime of this object.
  VertexIdEdgeVectorShape(vector<pair<int32, int32>> const& edges,
                          vector<S2Point> const& vertices)
      : edges_(edges), vertices_(vertices) {
  }

  S2Point const& vertex0(int e) const { return vertex(edges_[e].first); }
  S2Point const& vertex1(int e) const { return vertex(edges_[e].second); }

  // S2Shape interface:
  int num_edges() const { return edges_.size(); }
  void GetEdge(int e, S2Point const** v0, S2Point const** v1) const;
  bool has_interior() const { return false; }
  bool contains_origin() const { return false; }

 private:
  S2Point const& vertex(int i) const { return vertices_[i]; }

  vector<std::pair<int32, int32>> const& edges_;
  vector<S2Point> const& vertices_;
};

void VertexIdEdgeVectorShape::GetEdge(int e, S2Point const** v0,
                                      S2Point const** v1) const {
  auto const& edge = edges_[e];
  *v0 = &vertices_[edge.first];
  *v1 = &vertices_[edge.second];
}

bool S2Builder::Build(S2Error* error) {
  // CHECK rather than DCHECK because this is friendlier than crashing on the
  // "error->ok()" call below.  It would be easy to allow (error == nullptr)
  // by declaring a local "tmp_error", but it seems better to make clients
  // think about error handling.
  CHECK(error != nullptr);
  error_ = error;

  // Mark the end of the last layer.
  layer_begins_.push_back(input_edges_.size());

  // See the algorithm overview at the top of this file.
  if (snapping_requested_ && !options_.idempotent()) {
    snapping_needed_ = true;
  }
  ChooseSites();
  BuildLayers();
  Reset();
  return error->ok();
}

void S2Builder::Reset() {
  input_vertices_.clear();
  input_edges_.clear();
  layers_.clear();
  layer_options_.clear();
  layer_begins_.clear();
  label_set_ids_.clear();
  label_set_lexicon_.Clear();
  label_set_.clear();
  label_set_modified_ = false;
  sites_.clear();
  edge_sites_.clear();
  snapping_needed_ = false;
}

void S2Builder::ChooseSites() {
  if (input_vertices_.empty()) return;

  S2ShapeIndex input_edge_index;
  input_edge_index.Add(new VertexIdEdgeVectorShape(input_edges_,
                                                   input_vertices_));
  if (options_.split_crossing_edges()) {
    AddEdgeCrossings(input_edge_index);
  }
  if (snapping_requested_) {
    S2PointIndex<SiteId> site_index;
    AddForcedSites(&site_index);
    // "rejected_vertex_index" contains the input vertices that were not
    // selected as sites.  It is only needed to implement idempotency; it is
    // not maintained or used once we determine that the input needs to be
    // snapped (snapping_needed_ is set to true).
    S2PointIndex<InputVertexId> rejected_vertex_index;
    ChooseInitialSites(&site_index, &rejected_vertex_index);
    CollectSiteEdges(site_index, rejected_vertex_index);
  }
  if (snapping_needed_) {
    AddExtraSites(input_edge_index);
  } else {
    CopyInputEdges();
  }
}

void S2Builder::CopyInputEdges() {
  // Sort the input vertices, discard duplicates, and update the input edges
  // to refer to the pruned vertex list.  (We sort in the same order used by
  // ChooseInitialSites() to avoid inconsistencies in tests.)
  vector<InputVertexKey> sorted = SortInputVertices();
  vector<InputVertexId> vmap(input_vertices_.size());
  sites_.clear();
  sites_.reserve(input_vertices_.size());
  for (int in = 0; in < sorted.size(); ) {
    S2Point const& site = input_vertices_[sorted[in].second];
    vmap[sorted[in].second] = sites_.size();
    while (++in < sorted.size() && input_vertices_[sorted[in].second] == site) {
      vmap[sorted[in].second] = sites_.size();
    }
    sites_.push_back(site);
  }
  input_vertices_ = sites_;
  for (int i = 0; i < input_edges_.size(); ++i) {
    InputEdge& e = input_edges_[i];
    e.first = vmap[e.first];
    e.second = vmap[e.second];
  }
}

vector<S2Builder::InputVertexKey> S2Builder::SortInputVertices() {
  // Sort all the candidate sites (vertices) first by S2CellId and then by
  // S2Point value.  This improves the spatial locality of queries when
  // choosing the Voronoi sites.  It also will also allow the S2PointIndex
  // approach to be replaced by a divide-and-conquer algorithm if necessary.
  //
  // TODO(ericv): If the user is snapping to discrete points that are much
  // more closely spaced than snap_radius() (or not snapping to discrete
  // points at all), this scheme tends to bias the Voronoi site locations
  // towards points that are earlier on the S2CellId Hilbert curve.  A
  // reasonable compromise would be to sort by S2CellId at a coarse level
  // (down to O(snap_radius) in size), and then sort by fingerprint, or
  // possibly by distance to the cell center.
  vector<InputVertexKey> keys;
  keys.reserve(input_vertices_.size());
  for (InputVertexId i = 0; i < input_vertices_.size(); ++i) {
    keys.push_back(InputVertexKey(S2CellId::FromPoint(input_vertices_[i]), i));
  }
  std::sort(keys.begin(), keys.end(),
            [this](InputVertexKey const& a, InputVertexKey const& b) {
      if (a.first < b.first) return true;
      if (b.first < a.first) return false;
      return input_vertices_[a.second] < input_vertices_[b.second];
    });
  return keys;
}

// Check all edge pairs for crossings, and add the corresponding intersection
// points to input_vertices_.  (The intersection points will be snapped and
// merged with the other vertices during site selection.)
void S2Builder::AddEdgeCrossings(S2ShapeIndex const& input_edge_index) {
  s2shapeutil::EdgePairList edge_pairs;
  s2shapeutil::GetCrossingEdgePairs(
      input_edge_index, s2shapeutil::CrossingType::INTERIOR, &edge_pairs);
  auto shape = down_cast<VertexIdEdgeVectorShape*>(input_edge_index.shape(0));
  for (auto const& edge_pair : edge_pairs) {
    snapping_needed_ = true;
    InputEdgeId e0 = edge_pair.first.edge_id();
    InputEdgeId e1 = edge_pair.second.edge_id();
    AddVertex(S2EdgeUtil::GetIntersection(
        shape->vertex0(e0), shape->vertex1(e0),
        shape->vertex0(e1), shape->vertex1(e1)));
  }
}

void S2Builder::AddForcedSites(S2PointIndex<SiteId>* site_index) {
  // Sort the forced sites and remove duplicates.
  std::sort(sites_.begin(), sites_.end());
  sites_.erase(std::unique(sites_.begin(), sites_.end()), sites_.end());
  // Add the forced sites to the index.
  for (SiteId id = 0; id < sites_.size(); ++id) {
    site_index->Add(sites_[id], id);
  }
}

void S2Builder::ChooseInitialSites(
    S2PointIndex<SiteId>* site_index,
    S2PointIndex<InputVertexId>* rejected_vertex_index) {
  // Sort the input vertices in a determinstic order.
  vector<InputVertexKey> sorted = SortInputVertices();

  // Consider the vertices in order.  If an existing site is already close
  // enough, then nothing more needs to be done.  Otherwise we snap the vertex
  // and add it to the list of sites.
  S2ClosestPointQuery<SiteId> site_query(*site_index);
  site_query.set_max_distance(snap_radius_);
  S2ClosestPointQuery<InputVertexId> vertex_query(*rejected_vertex_index);
  vertex_query.set_max_distance(min_site_separation_);
  for (int i = 0; i < sorted.size(); ++i) {
    S2Point const& vertex = input_vertices_[sorted[i].second];
    site_query.FindClosestPoint(vertex);
    if (site_query.num_points() == 0) {
      S2Point site = SnapSite(vertex);
      if (site != vertex) snapping_needed_ = true;
      site_index->Add(site, sites_.size());
      sites_.push_back(site);
      site_query.Reset();
    } else if (!snapping_needed_) {
      // Check whether this vertex is too close to any other input vertex
      // (including input vertices that were not selected as Voronoi sites).
      if (site_query.distance_ca(0) < min_site_separation_ca_) {
        if (site_query.point(0) == vertex) continue;
        snapping_needed_ = true;
      } else {
        vertex_query.FindClosestPoint(vertex);
        if (vertex_query.num_points() > 0 &&
            vertex_query.distance_ca(0) < min_site_separation_ca_) {
          if (vertex_query.point(0) == vertex) continue;
          snapping_needed_ = true;
        } else {
          rejected_vertex_index->Add(vertex, sorted[i].second);
          vertex_query.Reset();
        }
      }
    }
  }
}

S2Point S2Builder::SnapSite(S2Point const& point) const {
  if (!snapping_requested_) {
    return point;
  }
  S2Point site = options_.snap_function().SnapPoint(point);
  S1ChordAngle dist_moved(site, point);
  if (dist_moved > site_snap_radius_ca_) {
    error_->Init(S2Error::BUILDER_SNAP_RADIUS_TOO_SMALL,
                 "Snap function moved vertex (%.15g, %.15g, %.15g) "
                 "by %.15g, which is more than the specified snap "
                 "radius of %.15g", point.x(), point.y(), point.z(),
                 dist_moved.ToAngle().radians(), snap_radius_.radians());
  }
  return site;
}

// For each edge, find all sites within min_edge_site_query_radius_ and store
// them in edge_sites_.  Also, to implement idempotency this method also
// checks whether the input vertices and edges may already satisfy the output
// criteria.  If any problems are found then snapping_needed_ is set to true.
void S2Builder::CollectSiteEdges(
    S2PointIndex<SiteId> const& site_index,
    S2PointIndex<InputVertexId> const& rejected_vertex_index) {
  edge_sites_.resize(input_edges_.size());
  S2ClosestPointQuery<SiteId> site_query(site_index);
  site_query.set_max_distance(edge_site_query_radius_);
  S2ClosestPointQuery<InputVertexId> vertex_query(rejected_vertex_index);
  vertex_query.set_max_distance(min_edge_site_separation_);
  for (InputEdgeId e = 0; e < input_edges_.size(); ++e) {
    InputEdge const& edge = input_edges_[e];
    S2Point const& v0 = input_vertices_[edge.first];
    S2Point const& v1 = input_vertices_[edge.second];
    site_query.FindClosestPointsToEdge(v0, v1);
    auto* sites = &edge_sites_[e];
    sites->reserve(site_query.num_points());
    for (int j = 0; j < site_query.num_points(); ++j) {
      sites->push_back(site_query.data(j));
      if (!snapping_needed_ &&
          site_query.distance_ca(j) < min_edge_site_separation_ca_ &&
          site_query.point(j) != v0 && site_query.point(j) != v1) {
        snapping_needed_ = true;
      }
    }
    if (!snapping_needed_) {
      vertex_query.FindClosestPointsToEdge(v0, v1);
      for (int j = 0; j < vertex_query.num_points(); ++j) {
        if (vertex_query.point(j) != v0 && vertex_query.point(j) != v1) {
          snapping_needed_ = true;
        }
      }
    }
    SortSitesByDistance(v0, sites);
  }
}

static bool IsCloserTo(S2Point const& x, S2Point const& a, S2Point const& b) {
  S1ChordAngle xa(x, a), xb(x, b);
  if (xa < xb) return true;
  if (xb < xa) return false;
  return a < b;
}

void S2Builder::SortSitesByDistance(S2Point const& x,
                                    compact_array<SiteId>* sites) const {
  // Sort sites in increasing order of distance to X.
  std::sort(sites->begin(), sites->end(),
            [&x, this](SiteId i, SiteId j) {
      return IsCloserTo(x, sites_[i], sites_[j]);
    });
}

// There are two situatons where we need to add extra Voronoi sites in order
// to ensure that the snapped edges meet the output requirements:
//
//  (1) If a snapped edge deviates from its input edge by more than
//      max_edge_deviation(), we add a new site on the input edge near the
//      middle of the snapped edge.  This causes the snapped edge to split
//      into two pieces, so that it follows the input edge more closely.
//
//  (2) If a snapped edge is closer than min_edge_vertex_separation() to any
//      nearby site (the "site to avoid"), then we add a new site (the
//      "separation site") on the input edge near the site to avoid.  This
//      causes the snapped edge to follow the input edge more closely and is
//      guaranteed to increase the separation to the required distance.
//
// We check these conditions by snapping all the input edges to a chain of
// Voronoi sites and then testing each edge in the chain.  If a site needs to
// be added, we mark all nearby edges for re-snapping.
void S2Builder::AddExtraSites(S2ShapeIndex const& input_edge_index) {
  vector<SiteId> chain;  // Temporary
  vector<InputEdgeId> snap_queue;
  for (InputEdgeId max_e = 0; max_e < input_edges_.size(); ++max_e) {
    snap_queue.push_back(max_e);
    while (!snap_queue.empty()) {
      InputEdgeId e = snap_queue.back();
      snap_queue.pop_back();
      SnapEdge(e, &chain);
      // We could save the snapped chain here in a snapped_chains_ vector, to
      // avoid resnapping it in AddSnappedEdges() below, however currently
      // SnapEdge only accounts for less than 5% of the runtime.
      MaybeAddExtraSites(e, max_e, chain, input_edge_index, &snap_queue);
    }
  }
}

void S2Builder::MaybeAddExtraSites(InputEdgeId edge_id,
                                   InputEdgeId max_edge_id,
                                   vector<SiteId> const& chain,
                                   S2ShapeIndex const& input_edge_index,
                                   vector<InputEdgeId>* snap_queue) {
  // The snapped chain is always a *subsequence* of the nearby sites
  // (edge_sites_), so we walk through the two arrays in parallel looking for
  // sites that weren't snapped.  We also keep track of the current snapped
  // edge, since it is the only edge that can be too close.
  int i = 0;
  for (SiteId id : edge_sites_[edge_id]) {
    if (id == chain[i]) {
      if (++i == chain.size()) break;
      // Check whether this snapped edge deviates too far from its original
      // position.  If so, we split the edge by adding an extra site.
      S2Point const& v0 = sites_[chain[i - 1]];
      S2Point const& v1 = sites_[chain[i]];
      if (S1ChordAngle(v0, v1) < min_edge_length_to_split_ca_) continue;

      InputEdge const& edge = input_edges_[edge_id];
      S2Point const& a0 = input_vertices_[edge.first];
      S2Point const& a1 = input_vertices_[edge.second];
      if (!S2EdgeUtil::IsEdgeBNearEdgeA(a0, a1, v0, v1, max_edge_deviation_)) {
        // Add a new site on the input edge, positioned so that it splits the
        // snapped edge into two approximately equal pieces.  Then we find all
        // the edges near the new site (including this one) and add them to
        // the snap queue.
        //
        // Note that with large snap radii, it is possible that the snapped
        // edge wraps around the sphere the "wrong way".  To handle this we
        // find the preferred split location by projecting both endpoints onto
        // the input edge and taking their midpoint.
        S2Point mid = (S2EdgeUtil::GetClosestPoint(v0, a0, a1) +
                       S2EdgeUtil::GetClosestPoint(v1, a0, a1)).Normalize();
        S2Point new_site = GetSeparationSite(mid, v0, v1, edge_id);
        AddExtraSite(new_site, max_edge_id, input_edge_index, snap_queue);
        return;
      }
    } else if (i > 0) {
      // Check whether this "site to avoid" is closer to the snapped edge than
      // min_edge_vertex_separation().  We only need to consider the edge
      // interior, because all sites are separated by at least snap_radius().
      // Similarly, *only* this edge can be too close because its vertices
      // must span the point where "site_to_avoid" projects onto the input
      // edge XY.  Both of these claims rely on the fact that all sites are
      // separated by at least snap_radius().  If ForceVertex() is used to
      // force sites to be closer than this, nothing bad happens except that
      // we may not be able to ensure sufficient separation.
      S2Point const& site_to_avoid = sites_[id];
      S2Point const& v0 = sites_[chain[i - 1]];
      S2Point const& v1 = sites_[chain[i]];
      if (S2EdgeUtil::IsInteriorDistanceLess(
              site_to_avoid, v0, v1, min_edge_site_separation_ca_)) {
        // A snapped edge can only approach a site too closely when there are
        // no sites near the input edge near that point.  We fix that by
        // adding a new site along the input edge (a "separation site"), then
        // we find all the edges near the new site (including this one) and
        // add them to the snap queue.
        S2Point new_site = GetSeparationSite(site_to_avoid, v0, v1, edge_id);
        DCHECK_NE(site_to_avoid, new_site);
        AddExtraSite(new_site, max_edge_id, input_edge_index, snap_queue);
        return;
      }
    }
  }
}

// Adds a new site, then updates "edge_sites"_ for all edges near the new site
// and adds them to "snap_queue" for resnapping (unless their edge id exceeds
// "max_edge_id", since those edges have not been snapped the first time yet).
void S2Builder::AddExtraSite(S2Point const& new_site,
                             InputEdgeId max_edge_id,
                             S2ShapeIndex const& input_edge_index,
                             vector<InputEdgeId>* snap_queue) {
  SiteId new_site_id = sites_.size();
  sites_.push_back(new_site);
  S2ClosestEdgeQuery query(input_edge_index);
  query.set_max_distance(edge_site_query_radius_);
  query.FindClosestEdges(new_site);
  for (int k = 0; k < query.num_edges(); ++k) {
    InputEdgeId e = query.edge_id(k);
    auto* site_ids = &edge_sites_[e];
    site_ids->push_back(new_site_id);
    SortSitesByDistance(input_vertices_[input_edges_[e].first], site_ids);
    if (e <= max_edge_id) snap_queue->push_back(e);
  }
}

S2Point S2Builder::GetSeparationSite(S2Point const& site_to_avoid,
                                     S2Point const& v0, S2Point const& v1,
                                     InputEdgeId input_edge_id) const {
  // Given an edge XY and a site S, consider the part of XY that intersects
  // the coverage disc of S.  We call this the "coverage interval" of S.  The
  // SnapFunction implementations guarantee that the only way that a snapped
  // edge can be closer than min_edge_vertex_separation() to a non-snapped
  // site (i.e., site_to_avoid) if is there is a gap in the coverage of XY
  // near this site.  We can fix this problem simply by adding a new site to
  // fill this gap, located as closely as possible to the site to avoid.
  //
  // To calculate the coverage gap, we look at the two snapped sites on
  // either side of site_to_avoid, and find the endpoints of their coverage
  // intervals.  The we place a new site in the gap, located as closely as
  // possible to the site to avoid.  Note that the new site may move when it
  // is snapped by the snap_function, but it is guaranteed not to move by
  // more than snap_radius and therefore its coverage interval will still
  // intersect the gap.
  InputEdge const& edge = input_edges_[input_edge_id];
  S2Point const& x = input_vertices_[edge.first];
  S2Point const& y = input_vertices_[edge.second];
  Vector3_d xy_dir = y - x;
  S2Point n = S2::RobustCrossProd(x, y);
  S2Point new_site = S2EdgeUtil::GetClosestPoint(site_to_avoid, x, y, n);
  S2Point gap_min = GetCoverageEndpoint(v0, x, y, n);
  S2Point gap_max = GetCoverageEndpoint(v1, y, x, -n);
  if ((new_site - gap_min).DotProd(xy_dir) < 0) {
    new_site = gap_min;
  } else if ((gap_max - new_site).DotProd(xy_dir) < 0) {
    new_site = gap_max;
  }
  new_site = SnapSite(new_site);
  DCHECK_NE(v0, new_site);
  DCHECK_NE(v1, new_site);
  return new_site;
}

// Given a site P and an edge XY with normal N, intersect XY with the disc of
// radius snap_radius() around P, and return the intersection point that is
// further along the edge XY toward Y.
S2Point S2Builder::GetCoverageEndpoint(S2Point const& p, S2Point const& x,
                                       S2Point const& y, S2Point const& n)
    const {
  // Consider the plane perpendicular to P that cuts off a spherical cap of
  // radius snap_radius().  This plane intersects the plane through the edge
  // XY (perpendicular to N) along a line, and that line intersects the unit
  // sphere at two points Q and R, and we want to return the point R that is
  // further along the edge XY toward Y.
  //
  // Let M be the midpoint of QR.  This is the point along QR that is closest
  // to P.  We can now express R as the sum of two perpendicular vectors OM
  // and MR in the plane XY.  Vector MR is in the direction N x P, while
  // vector OM is in the direction (N x P) x N, where N = X x Y.
  //
  // The length of OM can be found using the Pythagorean theorem on triangle
  // OPM, and the length of MR can be found using the Pythagorean theorem on
  // triangle OMR.
  //
  // In the calculations below, we save some work by scaling all the vectors
  // by n.CrossProd(p).Norm2(), and normalizing at the end.
  double n2 = n.Norm2();
  double nDp = n.DotProd(p);
  S2Point nXp = n.CrossProd(p);
  S2Point nXpXn = n2 * p - nDp * n;
  Vector3_d om = sqrt(1 - edge_snap_radius_sin2_) * nXpXn;
  double mr2 = edge_snap_radius_sin2_ * n2 - nDp * nDp;

  // MR is constructed so that it points toward Y (rather than X).
  Vector3_d mr = sqrt(max(0.0, mr2)) * nXp;
  return (om + mr).Normalize();
}

void S2Builder::SnapEdge(InputEdgeId e, vector<SiteId>* chain) {
  chain->clear();
  InputEdge const& edge = input_edges_[e];
  if (!snapping_needed_) {
    chain->push_back(edge.first);
    chain->push_back(edge.second);
    return;
  }

  // Optimization: if there is only one nearby site, return.
  // Optimization: if there are exactly two nearby sites, and one is close
  // enough to each vertex, then return.

  // Compute the edge normal.
  S2Point const& x = input_vertices_[edge.first];
  S2Point const& y = input_vertices_[edge.second];
  S2Point n = S2::RobustCrossProd(x, y);  // Not normalized.

  // Now iterate through the sites.  We keep track of the sequence of sites
  // that are visited.
  auto const& candidates = edge_sites_[e];
  for (int i = 0; i < candidates.size(); ++i) {
    S2Point const& c = sites_[candidates[i]];
    // Skip any sites that are too far away.  (There will be some of these,
    // because we also keep track of "sites to avoid".)  Note that some sites
    // may be close enough to the line containing the edge, but not to the
    // edge itself, so we can just use the dot product with the edge normal.
    if (!S2EdgeUtil::IsDistanceLess(c, x, y, edge_snap_radius_ca_)) continue;

    // Check whether the new site C excludes the previous site B.  If so,
    // repeat with the previous site, and so on.
    bool add_site_c = true;
    for (; !chain->empty(); chain->pop_back()) {
      S2Point b = sites_[chain->back()];

      // First, check whether B and C are so far apart that their clipped
      // Voronoi regions can't intersect.
      S1ChordAngle bc(b, c);
      if (bc >= max_adjacent_site_separation_ca_) break;

      // Otherwise, we want to check whether site C prevents the Voronoi
      // region of B from intersecting XY, or vice versa.  This can determined
      // by computing the "coverage interval" (the segment of XY intersected
      // by the coverage disc of radius snap_radius) for each site.  If the
      // coverage interval of one site contains the coverage interval of the
      // other, then the contained site can be excluded.

      ExclusionResult result = GetVoronoiSiteExclusion(b, c, x, y, n);
      if (result == FIRST_EXCLUDED) continue;  // Site B is excluded by C.
      if (result == SECOND_EXCLUDED) {
        add_site_c = false;  // Site C is excluded by B.
        break;
      }
      // Otherwise check whether the previous site A is close enough to B and
      // C that it might further clip the Voronoi region of B.
      if (chain->size() < 2) break;
      S2Point a = sites_[chain->end()[-2]];
      S1ChordAngle ac(a, c);
      if (ac >= max_adjacent_site_separation_ca_) break;

      // Optional: if AB > max_adjacent_site_separation_ca_ then keep B.
      // Optional: if d(B, XY) < 0.5 * min(AB, BC) then keep B.

      // Compute the Voronoi vertex of A, B, C.
      int xyb = S2::Sign(x, y, b);
      if (S2::Sign(a, b, c) == xyb) {
        break;  // The Voronoi vertex is on the same side as B but further away.
      }
      S2Point norm_ab = S2::RobustCrossProd(a, b);
      S2Point norm_bc = S2::RobustCrossProd(b, c);
      S2Point mid_ab = a + b;
      S2Point mid_bc = b + c;
      S2Point edge_ab = norm_ab.CrossProd(mid_ab);
      S2Point edge_bc = norm_bc.CrossProd(mid_bc);
      S2Point v = edge_ab.CrossProd(edge_bc).Normalize();
      if (v.DotProd(b) < 0) v = -v;
      if (S1ChordAngle(b, v) > edge_snap_radius_ca_ ||
          S2::Sign(x, y, v) != xyb) {
        break;  // The Voronoi vertex is on the opposite side of XY from B.
      }
      // Otherwise B is excluded by A and C combined.
    }
    if (add_site_c) {
      chain->push_back(candidates[i]);
    }
  }
  if (s2builder_verbose) {
    std::cout << "(" << edge.first << "," << edge.second << "): ";
    for (SiteId id : *chain) {
      std::cout << id << " ";
    }
    std::cout << "\n";
  }
}

// Given two sites A and B and an edge XY with normal N, determine whether one
// site excludes the other site's limited radius Voronoi region from
// intersecting the edge XY.
S2Builder::ExclusionResult S2Builder::GetVoronoiSiteExclusion(
    S2Point const& a, S2Point const& b, S2Point const& x, S2Point const& y,
    S2Point const& n) const {
  // The following test ensures that that when two sites A and B are exactly
  // symmetric around XY (mirror images), that (1) we only snap to one of
  // them, and (2) this choice doesn't depend on the direction of XY.  This is
  // true because IsCloserTo() resolves ties by comparing "a" and "b"
  // lexicographically, so the smaller of "a" and "b" always wins.
  DCHECK(IsCloserTo(x, a, b));
  if (IsCloserTo(y, a, b)) {
    return SECOND_EXCLUDED;  // Site A is closer to every point on XY.
  }

  // Given an edge XY and a site S, consider the part of XY that intersects
  // the coverage disc of S.  We call this the "coverage interval" of S.  It
  // is essentially an interval along the great circle through XY, and can be
  // represented as the point at the center of the interval and an angle that
  // measures the semi-width or "radius" of the interval.
  //
  // To test whether site A excludes site B along the input edge XY, we test
  // whether the coverage interval of A contains the coverage interval of B.
  // Let "ra" and "rb" be the radii (semi-widths) of the two intervals, and
  // let "d" be the angle between their center points.  Then "a" contains "b"
  // if (rb + d <= ra), and "b" contains "a" if (rb - d >= ra).
  //
  // It is not efficient to work with the actual angles, so instead the code
  // below uses the sines of the angles, e.g. sin(rb + d) <= sin(ra).  This
  // works fine most of the time, but for the first condition (rb + d <= ra),
  // we also need to check that rb + d < 90 degrees.
  //
  // The actual calculation is based on the following.  Let A1 and B1 be the
  // points A and B scaled by sqrt(1 - snap_radius**2) [these points are
  // inside the sphere].  The planes perpendicular to OA1 and OA2 cut off two
  // discs of radius "snap_radius" around A and B.  Now consider the two lines
  // where these planes intersect the plane containing the input edge XY, and
  // let A2 and B2 be the points on these lines that are closest to A and B.
  // The coverage intervals of A and B can be represented as two intervals
  // along these lines, centered around A2 and B2.  Let P1 and P2 be the
  // endpoints of the coverage interval for A, and let Q1 and Q2 be the
  // endpoints of the coverage interval for B.
  //
  // To check whether B's interval is contained by A's interval, we test
  // whether both endpoints of B's interval (Q1 and Q2) are contained by A's
  // interval.  We could test whether Qi.DotProd(A2) > A2.Norm2(), where Qi is
  // the endpoint whose dot product with A2 is smaller.  However it is better
  // numerically to test whether Q1.DotProd(P1 - A2) < (P1 - A2).Norm2().
  // Essentially, the first test is equivalent to using the cosines of the
  // angles involved, while the second test is equivalent to using the sines
  // of the angles.
  //
  // The code below is more complicated than this because we don't actually
  // construct the points A1, A2, P1, Q1, and so on.  Instead we compute the
  // various lengths directly.  Also, we split up many of the calculations
  // into perpendicular components, e.g. Q1.DotProd(A2) is computed as
  // Q1.DotProd(B2) + Q1.DotProd(Q1 - B2).  Finally, everything is scaled up
  // by a factor of aXn.Norm() * bXn.Norm2() to avoid divisions.
  //
  // Potential optimization: Some of these calculations are also done in
  // SnapEdge and could be passed in here.
  double n2 = n.Norm2();
  double aDn = a.DotProd(n);
  double bDn = b.DotProd(n);
  S2Point aXn = a.CrossProd(n);
  S2Point bXn = b.CrossProd(n);
  double edge_snap_radius_cos = sqrt(1 - edge_snap_radius_sin2_);
  double a_radius = sqrt(max(0.0, edge_snap_radius_sin2_ * n2 - aDn * aDn));
  double b_radius = sqrt(max(0.0, edge_snap_radius_sin2_ * n2 - bDn * bDn));

  // Even though the sites are sorted so that A is closer to X and B is closer
  // to Y, it is possible for the quantity below (equivalent to (AxB).(XxY))
  // to be negative if the edge XY is very long and the sites A and B are past
  // the endpoints, so that AB wraps around the sphere in the opposite
  // direction from XY.
  double aXbDn = std::fabs(a.DotProd(bXn));
  double aXnDbXn = aXn.DotProd(bXn);
  double sin_d_cos_rb = edge_snap_radius_cos * n2 * aXbDn;
  double sin_rb_cos_d = b_radius * aXnDbXn;
  double sin_ra = a_radius * bXn.Norm2();
  if (sin_rb_cos_d - sin_d_cos_rb >= sin_ra) return FIRST_EXCLUDED;
  if (sin_rb_cos_d + sin_d_cos_rb <= sin_ra &&
      edge_snap_radius_cos * aXnDbXn > b_radius * aXbDn) {
    return SECOND_EXCLUDED;
  }
  return NEITHER_EXCLUDED;
}

void S2Builder::BuildLayers() {
  // Each output edge has an "input edge id set id" (an int32) representing
  // the set of input edge ids that were snapped to this edge.  The actual
  // InputEdgeIds can be retrieved using "input_edge_id_set_lexicon".
  vector<vector<Edge>> layer_edges;
  vector<vector<InputEdgeIdSetId>> layer_input_edge_ids;
  IdSetLexicon input_edge_id_set_lexicon;
  BuildLayerEdges(&layer_edges, &layer_input_edge_ids,
                  &input_edge_id_set_lexicon);

  // If there are a large number of layers, then we build a minimal subset of
  // vertices for each layer.  This ensures that layer types that iterate over
  // vertices will run in time proportional to the size of that layer rather
  // than the size of all layers combined.
  vector<vector<S2Point>> layer_vertices;
  if (layers_.size() > 10) {
    vector<Graph::VertexId> filter_tmp;  // Temporary used by FilterVertices.
    layer_vertices.resize(layers_.size());
    for (int i = 0; i < layers_.size(); ++i) {
      layer_vertices[i] = Graph::FilterVertices(sites_, &layer_edges[i],
                                                &filter_tmp);
    }
    vector<S2Point>().swap(sites_);  // Release memory
  }
  for (int i = 0; i < layers_.size(); ++i) {
    GraphOptions const& options = layer_options_[i];
    vector<Edge>& edges = layer_edges[i];
    vector<InputEdgeIdSetId>& input_edge_ids = layer_input_edge_ids[i];
    vector<S2Point>& vertices = (layer_vertices.empty() ? sites_
                                                        : layer_vertices[i]);
    Graph graph(options, vertices, edges, input_edge_ids,
                input_edge_id_set_lexicon, label_set_ids_,
                label_set_lexicon_);
    layers_[i]->Build(graph, error_);

    // Clear data as we go along to save space.
    vector<Edge>().swap(edges);
    vector<InputEdgeIdSetId>().swap(input_edge_ids);
    if (&vertices != &sites_) vector<S2Point>().swap(vertices);
  }
}

static void DumpEdges(vector<S2Builder::Graph::Edge> const& edges,
                      vector<S2Point> const& vertices) {
  for (auto const& e : edges) {
    vector<S2Point> v;
    v.push_back(vertices[e.first]);
    v.push_back(vertices[e.second]);
    std::cout << "S2Polyline: " << s2textformat::ToString(v)
              << "(" << e.first << "," << e.second << ")\n";
  }
}

void S2Builder::BuildLayerEdges(
    vector<vector<Edge>>* layer_edges,
    vector<vector<InputEdgeIdSetId>>* layer_input_edge_ids,
    IdSetLexicon* input_edge_id_set_lexicon) {
  layer_edges->resize(layers_.size());
  layer_input_edge_ids->resize(layers_.size());
  for (int i = 0; i < layers_.size(); ++i) {
    GraphOptions* options = &layer_options_[i];
    vector<Edge>& edges = (*layer_edges)[i];
    vector<InputEdgeIdSetId>& input_edge_ids = (*layer_input_edge_ids)[i];
    AddSnappedEdges(layer_begins_[i], layer_begins_[i+1], *options,
                    &edges, &input_edge_ids, input_edge_id_set_lexicon);
    if (s2builder_verbose) DumpEdges(edges, sites_);
    // The errors generated by ProcessEdges are really warnings, so we simply
    // record them and continue.
    Graph::ProcessEdges(options, &edges, &input_edge_ids,
                        input_edge_id_set_lexicon, error_);
  }
  if (snapping_needed_ && options().simplify_edge_chains()) {
    SimplifyEdgeChains(layer_edges, layer_input_edge_ids,
                       input_edge_id_set_lexicon);
  }
  // At this point we have no further need for the input geometry or nearby
  // site data, so we clear those fields to save space.
  vector<S2Point>().swap(input_vertices_);
  vector<InputEdge>().swap(input_edges_);
  vector<compact_array<SiteId>>().swap(edge_sites_);
}

void S2Builder::AddSnappedEdges(
    InputEdgeId begin, InputEdgeId end, GraphOptions const& options,
    vector<Edge>* edges, vector<InputEdgeIdSetId>* input_edge_ids,
    IdSetLexicon* input_edge_id_set_lexicon) {
  bool keep_degenerate_edges = (options.degenerate_edges() ==
                                GraphOptions::DegenerateEdges::KEEP);
  vector<SiteId> chain;
  for (InputEdgeId e = begin; e < end; ++e) {
    InputEdgeIdSetId id = input_edge_id_set_lexicon->AddSingleton(e);
    SnapEdge(e, &chain);
    if (chain.size() == 1) {
      if (!keep_degenerate_edges) continue;
      AddSnappedEdge(options, chain[0], chain[0], id,
                     edges, input_edge_ids);
    } else {
      for (int i = 1; i < chain.size(); ++i) {
        AddSnappedEdge(options, chain[i-1], chain[i], id,
                     edges, input_edge_ids);
      }
    }
  }
}

inline void S2Builder::AddSnappedEdge(
    GraphOptions const& options, SiteId src, SiteId dst, InputEdgeIdSetId id,
    vector<Edge>* edges, vector<InputEdgeIdSetId>* input_edge_ids) {
  edges->push_back(Edge(src, dst));
  input_edge_ids->push_back(id);
  if (options.edge_type() == EdgeType::UNDIRECTED) {
    edges->push_back(Edge(dst, src));
    input_edge_ids->push_back(id);
  }
}

void S2Builder::SimplifyEdgeChains(
    vector<vector<Edge>>* layer_edges,
    vector<vector<InputEdgeIdSetId>>* layer_input_edge_ids,
    IdSetLexicon* input_edge_id_set_lexicon) {
  LOG(FATAL) << "simplify_edge_chains() not implemented";
}
