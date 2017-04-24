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

#include "s2/s2builderutil_layers.h"

#include <algorithm>
#include <array>
#include <memory>
#include "s2/third_party/absl/memory/memory.h"
#include "s2/s2builder_graph.h"
#include "s2/s2debug.h"
#include "s2/s2polygon.h"

using std::make_pair;
using std::pair;
using std::unique_ptr;
using std::vector;

using EdgeType = S2Builder::EdgeType;
using Graph = S2Builder::Graph;
using GraphOptions = S2Builder::GraphOptions;
using Label = S2Builder::Label;

using DegenerateEdges = GraphOptions::DegenerateEdges;
using DuplicateEdges = GraphOptions::DuplicateEdges;
using SiblingPairs = GraphOptions::SiblingPairs;

using EdgeId = Graph::EdgeId;
using InputEdgeId = Graph::InputEdgeId;
using LoopType = Graph::LoopType;
using PolylineType = Graph::PolylineType;
using VertexId = Graph::VertexId;

namespace s2builderutil {

S2PolygonLayer::S2PolygonLayer(S2Polygon* polygon, Options const& options) {
  Init(polygon, nullptr, nullptr, options);
}

S2PolygonLayer::S2PolygonLayer(
    S2Polygon* polygon, LabelSetIds* label_set_ids,
    IdSetLexicon* label_set_lexicon, Options const& options) {
  Init(polygon, label_set_ids, label_set_lexicon, options);
}

void S2PolygonLayer::Init(
    S2Polygon* polygon, LabelSetIds* label_set_ids,
    IdSetLexicon* label_set_lexicon, Options const& options) {
  DCHECK_EQ(label_set_ids == nullptr, label_set_lexicon == nullptr);
  polygon_ = polygon;
  label_set_ids_ = label_set_ids;
  label_set_lexicon_ = label_set_lexicon;
  options_ = options;

  if (options_.validate()) {
    polygon_->set_s2debug_override(S2Debug::DISABLE);
  }
}

GraphOptions S2PolygonLayer::graph_options() const {
  // Prevent degenerate edges and sibling edge pairs.  There should not be any
  // duplicate edges if the input is valid, but if there are then we keep them
  // since this tends to produce more comprehensible errors.
  GraphOptions graph_options;
  graph_options.set_edge_type(options_.edge_type());
  graph_options.set_degenerate_edges(DegenerateEdges::DISCARD);
  graph_options.set_duplicate_edges(DuplicateEdges::KEEP);
  graph_options.set_sibling_pairs(SiblingPairs::DISCARD);
  return graph_options;
}

void S2PolygonLayer::AppendS2Loops(Graph const& g,
                                   vector<Graph::EdgeLoop> const& edge_loops,
                                   vector<unique_ptr<S2Loop>>* loops) const {
  vector<S2Point> vertices;
  for (auto const& edge_loop : edge_loops) {
    vertices.reserve(edge_loop.size());
    for (auto edge_id : edge_loop) {
      vertices.push_back(g.vertex(g.edge(edge_id).first));
    }
    loops->push_back(
        absl::MakeUnique<S2Loop>(vertices, polygon_->s2debug_override()));
    vertices.clear();
  }
}

void S2PolygonLayer::AppendEdgeLabels(
    Graph const& g,
    vector<Graph::EdgeLoop> const& edge_loops) {
  if (!label_set_ids_) return;

  vector<Label> labels;  // Temporary storage for labels.
  Graph::LabelFetcher fetcher(g, options_.edge_type());
  for (auto const& edge_loop : edge_loops) {
    vector<LabelSetId> loop_label_set_ids;
    loop_label_set_ids.reserve(edge_loop.size());
    for (auto edge_id : edge_loop) {
      fetcher.Fetch(edge_id, &labels);
      loop_label_set_ids.push_back(label_set_lexicon_->Add(labels));
    }
    label_set_ids_->push_back(std::move(loop_label_set_ids));
  }
}

void S2PolygonLayer::InitLoopMap(vector<unique_ptr<S2Loop>> const& loops,
                                 LoopMap* loop_map) const {
  if (!label_set_ids_) return;
  for (auto const& loop : loops) {
    (*loop_map)[&*loop] =
        make_pair(&loop - &loops[0], loop->contains_origin());
  }
}

void S2PolygonLayer::ReorderEdgeLabels(LoopMap const& loop_map) {
  if (!label_set_ids_) return;
  LabelSetIds new_ids(label_set_ids_->size());
  for (int i = 0; i < polygon_->num_loops(); ++i) {
    S2Loop* loop = polygon_->loop(i);
    pair<int, bool> const& old = loop_map.find(loop)->second;
    new_ids[i].swap((*label_set_ids_)[old.first]);
    if (loop->contains_origin() != old.second) {
      // S2Loop::Invert() reverses the order of the vertices, which leaves
      // the last edge unchanged.  For example, the loop ABCD (with edges
      // AB, BC, CD, DA) becomes the loop DCBA (with edges DC, CB, BA, AD).
      std::reverse(new_ids[i].begin(), new_ids[i].end() - 1);
    }
  }
  label_set_ids_->swap(new_ids);
}

void S2PolygonLayer::Build(Graph const& g, S2Error* error) {
  if (label_set_ids_) label_set_ids_->clear();

  // It's tricky to compute the edge labels for S2Polygons because the
  // S2Polygon::Init methods can reorder and/or invert the loops.  We handle
  // this by remembering the original vector index of each loop and whether or
  // not the loop contained S2::Origin().  By comparing this with the final
  // S2Polygon loops we can fix up the edge labels appropriately.
  LoopMap loop_map;
  if (g.options().edge_type() == EdgeType::DIRECTED) {
    vector<Graph::EdgeLoop> edge_loops;
    if (!g.GetDirectedLoops(LoopType::SIMPLE, &edge_loops, error)) {
      return;
    }
    vector<unique_ptr<S2Loop>> loops;
    AppendS2Loops(g, edge_loops, &loops);
    AppendEdgeLabels(g, edge_loops);
    vector<Graph::EdgeLoop>().swap(edge_loops);  // Release memory
    InitLoopMap(loops, &loop_map);
    polygon_->InitOriented(std::move(loops));
  } else {
    vector<Graph::UndirectedComponent> components;
    if (!g.GetUndirectedComponents(LoopType::SIMPLE, &components, error)) {
      return;
    }
    // It doesn't really matter which complement of each component we use,
    // since S2Polygon::InitNested() automatically normalizes all the loops.
    // The only reason to prefer one over the other is that when there are
    // multiple loops that touch, only one of the two complements matches the
    // structure of the input loops.  GetUndirectedComponents() tries to
    // ensure that this is always component 0.
    vector<unique_ptr<S2Loop>> loops;
    for (auto const& component : components) {
      AppendS2Loops(g, component[0], &loops);
      AppendEdgeLabels(g, component[0]);
    }
    vector<Graph::UndirectedComponent>().swap(components);  // Release memory
    InitLoopMap(loops, &loop_map);
    for (auto const& loop : loops) loop->Normalize();
    polygon_->InitNested(std::move(loops));
  }
  ReorderEdgeLabels(loop_map);
  if (options_.validate()) {
    polygon_->FindValidationError(error);
  }
}

S2PolylineLayer::S2PolylineLayer(S2Polyline* polyline,
                                 S2PolylineLayer::Options const& options) {
  Init(polyline, nullptr, nullptr, options);
}

S2PolylineLayer::S2PolylineLayer(
    S2Polyline* polyline, LabelSetIds* label_set_ids,
    IdSetLexicon* label_set_lexicon, Options const& options) {
  Init(polyline, label_set_ids, label_set_lexicon, options);
}

void S2PolylineLayer::Init(S2Polyline* polyline, LabelSetIds* label_set_ids,
                           IdSetLexicon* label_set_lexicon,
                           Options const& options) {
  DCHECK_EQ(label_set_ids == nullptr, label_set_lexicon == nullptr);
  polyline_ = polyline;
  label_set_ids_ = label_set_ids;
  label_set_lexicon_ = label_set_lexicon;
  options_ = options;

  if (options_.validate()) {
    polyline_->set_s2debug_override(S2Debug::DISABLE);
  }
}

GraphOptions S2PolylineLayer::graph_options() const {
  // Remove edges that collapse to a single vertex, but keep duplicate and
  // sibling edges, since merging duplicates or discarding siblings can make
  // it impossible to assemble the edges into a single polyline.
  GraphOptions graph_options;
  graph_options.set_edge_type(options_.edge_type());
  graph_options.set_degenerate_edges(DegenerateEdges::DISCARD);
  graph_options.set_duplicate_edges(DuplicateEdges::KEEP);
  graph_options.set_sibling_pairs(SiblingPairs::KEEP);
  return graph_options;
}

void S2PolylineLayer::Build(Graph const& g, S2Error* error) {
  if (g.num_edges() == 0) {
    polyline_->Init(vector<S2Point>{});
    return;
  }
  vector<Graph::EdgePolyline> edge_polylines =
      g.GetPolylines(PolylineType::WALK);
  if (edge_polylines.size() != 1) {
    error->Init(S2Error::BUILDER_EDGES_DO_NOT_FORM_POLYLINE,
                "Input edges cannot be assembled into polyline");
    return;
  }
  Graph::EdgePolyline const& edge_polyline = edge_polylines[0];
  vector<S2Point> vertices;  // Temporary storage for vertices.
  vertices.reserve(edge_polyline.size());
  vertices.push_back(g.vertex(g.edge(edge_polyline[0]).first));
  for (EdgeId e : edge_polyline) {
    vertices.push_back(g.vertex(g.edge(e).second));
  }
  if (label_set_ids_) {
    Graph::LabelFetcher fetcher(g, options_.edge_type());
    vector<Label> labels;  // Temporary storage for labels.
    label_set_ids_->reserve(edge_polyline.size());
    for (EdgeId e : edge_polyline) {
      fetcher.Fetch(e, &labels);
      label_set_ids_->push_back(label_set_lexicon_->Add(labels));
    }
  }
  polyline_->Init(vertices);
  if (options_.validate()) {
    polyline_->FindValidationError(error);
  }
}

S2PolylineVectorLayer::S2PolylineVectorLayer(
    vector<unique_ptr<S2Polyline>>* polylines,
    S2PolylineVectorLayer::Options const& options) {
  Init(polylines, nullptr, nullptr, options);
}

S2PolylineVectorLayer::S2PolylineVectorLayer(
    vector<unique_ptr<S2Polyline>>* polylines, LabelSetIds* label_set_ids,
    IdSetLexicon* label_set_lexicon, Options const& options) {
  Init(polylines, label_set_ids, label_set_lexicon, options);
}

void S2PolylineVectorLayer::Init(vector<unique_ptr<S2Polyline>>* polylines,
                                 LabelSetIds* label_set_ids,
                                 IdSetLexicon* label_set_lexicon,
                                 Options const& options) {
  DCHECK_EQ(label_set_ids == nullptr, label_set_lexicon == nullptr);
  polylines_ = polylines;
  label_set_ids_ = label_set_ids;
  label_set_lexicon_ = label_set_lexicon;
  options_ = options;
}

GraphOptions S2PolylineVectorLayer::graph_options() const {
  GraphOptions graph_options;
  graph_options.set_edge_type(options_.edge_type());
  graph_options.set_degenerate_edges(DegenerateEdges::DISCARD);
  graph_options.set_duplicate_edges(options_.duplicate_edges());
  graph_options.set_sibling_pairs(options_.sibling_pairs());
  return graph_options;
}

void S2PolylineVectorLayer::Build(Graph const& g, S2Error* error) {
  vector<Graph::EdgePolyline> edge_polylines = g.GetPolylines(
      options_.polyline_type());
  polylines_->reserve(edge_polylines.size());
  if (label_set_ids_) label_set_ids_->reserve(edge_polylines.size());
  vector<S2Point> vertices;  // Temporary storage for vertices.
  vector<Label> labels;  // Temporary storage for labels.
  for (auto const& edge_polyline : edge_polylines) {
    vertices.push_back(g.vertex(g.edge(edge_polyline[0]).first));
    for (EdgeId e : edge_polyline) {
      vertices.push_back(g.vertex(g.edge(e).second));
    }
    S2Polyline* polyline = new S2Polyline(vertices,
                                          options_.s2debug_override());
    vertices.clear();
    if (options_.validate()) {
      polyline->FindValidationError(error);
    }
    polylines_->emplace_back(polyline);
    if (label_set_ids_) {
      Graph::LabelFetcher fetcher(g, options_.edge_type());
      vector<LabelSetId> polyline_labels;
      polyline_labels.reserve(edge_polyline.size());
      for (EdgeId e : edge_polyline) {
        fetcher.Fetch(e, &labels);
        polyline_labels.push_back(label_set_lexicon_->Add(labels));
      }
      label_set_ids_->push_back(std::move(polyline_labels));
    }
  }
}

}  // namespace s2builderutil
