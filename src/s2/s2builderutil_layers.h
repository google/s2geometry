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
// A collection of standard Layer types.  There are essentially four types of
// output (where types in [square brackets] are not implemented yet):
//
//   Single polyline output:
//         S2PolylineLayer
//   Multiple polyline output:
//         S2PolylineVectorLayer
//   Single polygon output:
//         S2PolygonLayer
//         [LaxPolygonLayer]
//   Polygon mesh output:
//         [LaxPolygonVectorLayer]
//         [S2PolygonMeshLayer]
//
// Notice that there are two supported output types for polygons: S2Polygon
// and LaxPolygon. Use S2Polygon if you need the full range of operations
// that S2Polygon implements.  Use LaxPolygon if you want to represent
// polygons with zero-area degenerate regions, or if you need a type that has
// low memory overhead and fast initialization.  However, be aware that to
// convert from a LaxPolygon to an S2Polygon you will need to use S2Builder
// again.
//
// Similarly, there are two supported output formats for polygon meshes:
// LaxPolygonVector and S2PolygonMesh.  Use S2PolygonMesh if you need to be
// able to determine which polygons are adjacent to each edge or vertex;
// otherwise use LaxPolygonVector, which uses less memory and is faster to
// construct.

#ifndef S2_S2BUILDERUTIL_LAYERS_H_
#define S2_S2BUILDERUTIL_LAYERS_H_

#include <memory>
#include <utility>
#include <vector>
#include <glog/logging.h>
#include "s2/util/btree/btree_map.h"
#include "s2/id_set_lexicon.h"
#include "s2/s2builder.h"
#include "s2/s2builder_graph.h"
#include "s2/s2builder_layer.h"
#include "s2/s2debug.h"
#include "s2/s2error.h"
#include "s2/s2loop.h"
#include "s2/s2point.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"
#include "s2/s2shapeutil.h"

namespace s2builderutil {

// A layer type that assembles edges (directed or undirected) into an
// S2Polygon.  Returns an error if the edges cannot be assembled into loops.
//
// If the input edges are directed, they must be oriented such that the
// polygon interior is to the left of all edges.  Directed edges are always
// preferred (see S2Builder::EdgeType).
//
// Before the edges are assembled into loops, "sibling pairs" consisting of an
// edge and its reverse edge are automatically removed.  Such edge pairs
// represent zero-area degenerate regions, which S2Polygon does not allow.
// (If you need to build polygons with degeneracies, use LaxPolygonLayer
// instead.)
//
// S2PolygonLayer is implemented such that if the input to S2Builder is a
// polygon and is not modified, then the output has the same cyclic ordering
// of loop vertices and the same loop ordering as the input polygon.
//
// CAVEAT: Because polygons are constructed from their boundaries, this method
// cannot distinguish between the empty and full polygons.  An empty boundary
// always yields an empty polygon.  If the result should sometimes be the full
// polygon, such logic must be implemented outside of this class (and will
// need to consider factors other than the polygon's boundary).
class S2PolygonLayer : public S2Builder::Layer {
 public:
  class Options {
   public:
    // Constructor that uses the default options (listed below).
    Options();

    // Constructor that specifies the edge type.
    explicit Options(S2Builder::EdgeType edge_type);

    // Indicates whether the input edges provided to S2Builder are directed or
    // undirected.  Directed edges should be used whenever possible (see
    // S2Builder::EdgeType for details).
    //
    // If the input edges are directed, they should be oriented so that the
    // polygon interior is to the left of all edges.  This means that for a
    // polygon with holes, the outer loops ("shells") should be directed
    // counter-clockwise while the inner loops ("holes") should be directed
    // clockwise.  Note that S2Builder::AddPolygon() does this automatically.
    //
    // DEFAULT: S2Builder::EdgeType::DIRECTED
    S2Builder::EdgeType edge_type() const;
    void set_edge_type(S2Builder::EdgeType edge_type);

    // If true, calls FindValidationError() on the output polygon.  If any
    // error is found, it will be returned by S2Builder::Build().
    //
    // Note that this option calls set_s2debug_override(S2Debug::DISABLE) in
    // order to turn off the default error checking in debug builds.
    //
    // DEFAULT: false
    bool validate() const;
    void set_validate(bool validate);

   private:
    S2Builder::EdgeType edge_type_;
    bool validate_;
  };

  // Specifies that a polygon should be constructed using the given options.
  explicit S2PolygonLayer(S2Polygon* polygon,
                          Options const& options = Options());

  // Specifies that a polygon should be constructed using the given options,
  // and that any labels attached to the input edges should be returned in
  // "label_set_ids" and "label_set_lexicion".
  //
  // The labels associated with the edge "polygon.loop(i).vertex({j, j+1})"
  // can be retrieved as follows:
  //
  //   for (int32 label : label_set_lexicon.id_set(label_set_ids[i][j])) {...}
  using LabelSetIds = std::vector<std::vector<LabelSetId>>;
  S2PolygonLayer(S2Polygon* polygon, LabelSetIds* label_set_ids,
                 IdSetLexicon* label_set_lexicon,
                 Options const& options = Options());

  // Layer interface:
  GraphOptions graph_options() const override;
  void Build(Graph const& g, S2Error* error) override;

 private:
  void Init(S2Polygon* polygon, LabelSetIds* label_set_ids,
            IdSetLexicon* label_set_lexicon, Options const& options);
  void AppendS2Loops(Graph const& g,
                     std::vector<Graph::EdgeLoop> const& edge_loops,
                     std::vector<std::unique_ptr<S2Loop>>* loops) const;
  void AppendEdgeLabels(Graph const& g,
                        std::vector<Graph::EdgeLoop> const& edge_loops);
  using LoopMap = util::btree::btree_map<S2Loop*, std::pair<int, bool>>;
  void InitLoopMap(std::vector<std::unique_ptr<S2Loop>> const& loops,
                   LoopMap* loop_map) const;
  void ReorderEdgeLabels(LoopMap const& loop_map);

  S2Polygon* polygon_;
  LabelSetIds* label_set_ids_;
  IdSetLexicon* label_set_lexicon_;
  Options options_;
};

// Like S2PolygonLayer, but adds the polygon to an S2ShapeIndex.
class IndexedS2PolygonLayer : public S2Builder::Layer {
 public:
  using Options = S2PolygonLayer::Options;
  explicit IndexedS2PolygonLayer(S2ShapeIndex* index,
                                 Options const& options = Options())
      : index_(index), polygon_(new S2Polygon),
        layer_(polygon_.get(), options) {}

  GraphOptions graph_options() const override {
    return layer_.graph_options();
  }

  void Build(Graph const& g, S2Error* error) override {
    layer_.Build(g, error);
    if (error->ok()) {
      index_->Add(new S2Polygon::OwningShape(std::move(polygon_)));
    }
  }

 private:
  S2ShapeIndex* index_;
  std::unique_ptr<S2Polygon> polygon_;
  S2PolygonLayer layer_;
};

// a layer type that assembles edges (directed or undirected) into an
// S2Polyline.  Returns an error if the edges cannot be assembled into a
// single unbroken polyline.
//
// Duplicate edges are handled correctly (e.g., if a polyline backtracks on
// itself, or loops around and retraces some of its previous edges.)  The
// implementation attempts to preserve the order of directed input edges
// whenever possible, so that if the input is a polyline and it is not
// modified by S2Builder, then the output will be the same polyline (even if
// the polyline backtracks on itself or forms a loop).  With undirected edges,
// there are no such guarantees; for example, even if the input consists of a
// single undirected edge, then either directed edge may be returned.
//
// S2PolylineLayer does not support options such as discarding sibling pairs
// or merging duplicate edges because these options can split the polyline
// into several pieces.  Use S2PolylineVectorLayer if you need these features.
class S2PolylineLayer : public S2Builder::Layer {
 public:
  class Options {
   public:
    // Constructor that uses the default options (listed below).
    Options();

    // Constructor that specifies the edge type.
    explicit Options(S2Builder::EdgeType edge_type);

    // Indicates whether the input edges provided to S2Builder are directed or
    // undirected.  Directed edges should be used whenever possible to avoid
    // ambiguity.
    //
    // DEFAULT: S2Builder::EdgeType::DIRECTED
    S2Builder::EdgeType edge_type() const;
    void set_edge_type(S2Builder::EdgeType edge_type);

    // If true, calls FindValidationError() on the output polyline.  If any
    // error is found, it will be returned by S2Builder::Build().
    //
    // Note that this option calls set_s2debug_override(S2Debug::DISABLE) in
    // order to turn off the default error checking in debug builds.
    //
    // DEFAULT: false
    bool validate() const;
    void set_validate(bool validate);

   private:
    S2Builder::EdgeType edge_type_;
    bool validate_;
  };

  // Specifies that a polyline should be constructed using the given options.
  explicit S2PolylineLayer(S2Polyline* polyline,
                           Options const& options = Options());

  // Specifies that a polyline should be constructed using the given options,
  // and that any labels attached to the input edges should be returned in
  // "label_set_ids" and "label_set_lexicion".
  //
  // The labels associated with the edge "polyline.vertex({j, j+1})" can be
  // retrieved as follows:
  //
  //   for (int32 label : label_set_lexicon.id_set(label_set_ids[j])) {...}
  using LabelSetIds = std::vector<LabelSetId>;
  S2PolylineLayer(S2Polyline* polyline, LabelSetIds* label_set_ids,
                 IdSetLexicon* label_set_lexicon,
                 Options const& options = Options());

  // Layer interface:
  GraphOptions graph_options() const override;
  void Build(Graph const& g, S2Error* error) override;

 private:
  void Init(S2Polyline* polyline, LabelSetIds* label_set_ids,
            IdSetLexicon* label_set_lexicon, Options const& options);

  S2Polyline* polyline_;
  LabelSetIds* label_set_ids_;
  IdSetLexicon* label_set_lexicon_;
  Options options_;
};

// Like S2PolylineLayer, but adds the polyline to an S2ShapeIndex.
class IndexedS2PolylineLayer : public S2Builder::Layer {
 public:
  using Options = S2PolylineLayer::Options;
  explicit IndexedS2PolylineLayer(S2ShapeIndex* index,
                                 Options const& options = Options())
      : index_(index), polyline_(new S2Polyline),
        layer_(polyline_.get(), options) {}

  GraphOptions graph_options() const override {
    return layer_.graph_options();
  }

  void Build(Graph const& g, S2Error* error) override {
    layer_.Build(g, error);
    if (error->ok()) {
      index_->Add(new S2Polyline::OwningShape(std::move(polyline_)));
    }
  }

 private:
  S2ShapeIndex* index_;
  std::unique_ptr<S2Polyline> polyline_;
  S2PolylineLayer layer_;
};

// A layer type that assembles edges (directed or undirected) into multiple
// S2Polylines.  Returns an error if S2Builder found any problem with the
// input edges; this layer type does not generate any errors of its own.
//
// Duplicate edges are handled correctly (e.g., if a polyline backtracks on
// itself, or loops around and retraces some of its previous edges.)  The
// implementation attempts to preserve the order of the input edges whenever
// possible, so that if the input is a polyline and it is not modified by
// S2Builder, then the output will be the same polyline even if the polyline
// forms a loop.  However, note that this is not guaranteed when undirected
// edges are used: for example, if the input consists of a single undirected
// edge, then either directed edge may be returned.
class S2PolylineVectorLayer : public S2Builder::Layer {
 public:
  class Options {
   public:
    // Constructor that uses the default options (listed below).
    Options();

    // Constructor that specifies the edge type.
    explicit Options(S2Builder::EdgeType edge_type);

    // Indicates whether the input edges provided to S2Builder are directed or
    // undirected.
    //
    // Directed edges should be used whenever possible to avoid ambiguity.
    // The implementation attempts to preserve the structure of directed input
    // edges whenever possible, so that if the input is a vector of disjoint
    // polylines and none of them need to be modified then the output will be
    // the same polylines in the same order.  With undirected edges, there are
    // no such guarantees.
    //
    // DEFAULT: S2Builder::EdgeType::DIRECTED
    S2Builder::EdgeType edge_type() const;
    void set_edge_type(S2Builder::EdgeType edge_type);

    // Indicates whether polylines should be "paths" (which don't allow
    // duplicate vertices, except possibly the first and last vertex) or
    // "walks" (which allow duplicate vertices and edges).
    //
    // If your input consists of polylines, and you want to split them into
    // separate pieces whenever they self-intersect or cross each other, then
    // use PolylineType::PATH (and probably use split_crossing_edges()).  If
    // you don't mind if your polylines backtrack or contain loops, then use
    // PolylineType::WALK.
    //
    // DEFAULT: PolylineType::PATH
    using PolylineType = S2Builder::Graph::PolylineType;
    PolylineType polyline_type() const;
    void set_polyline_type(PolylineType polyline_type);

    // Indicates whether duplicate edges in the input should be kept (KEEP) or
    // merged together (MERGE).  Note you can use edge labels to determine
    // which input edges were merged into a given output edge.
    //
    // DEFAULT: DuplicateEdges::KEEP
    using DuplicateEdges = GraphOptions::DuplicateEdges;
    DuplicateEdges duplicate_edges() const;
    void set_duplicate_edges(DuplicateEdges duplicate_edges);

    // Indicates whether sibling edge pairs (i.e., pairs consisting of an edge
    // and its reverse edge) should be kept (KEEP) or discarded (DISCARD).
    // For example, if a polyline backtracks on itself, the DISCARD option
    // would cause this section of the polyline to be removed.  Note that this
    // option may cause a single polyline to split into several pieces (e.g.,
    // if a polyline has a "lollipop" shape).
    //
    // REQUIRES: sibling_pairs == { DISCARD, KEEP }
    //           (the CREATE and REQUIRE options are not allowed)
    //
    // DEFAULT: SiblingPairs::KEEP
    using SiblingPairs = GraphOptions::SiblingPairs;
    SiblingPairs sibling_pairs() const;
    void set_sibling_pairs(SiblingPairs sibling_pairs);

    // If true, calls FindValidationError() on each output polyline.  If any
    // error is found, it will be returned by S2Builder::Build().
    //
    // Note that this option calls set_s2debug_override(S2Debug::DISABLE) in
    // order to turn off the default error checking in debug builds.
    //
    // DEFAULT: false
    bool validate() const;
    void set_validate(bool validate);

    // This method can turn off the automatic validity checks triggered by the
    // --s2debug flag (which is on by default in debug builds).  The main
    // reason to do this is if your code already does its own error checking,
    // or if you need to work with invalid geometry for some reason.
    //
    // In any case, polylines have very few restrictions so they are unlikely
    // to have errors.  Errors include vertices that aren't unit length (which
    // can only happen if they are present in the input data), or adjacent
    // vertices that are at antipodal points on the sphere (unlikely with real
    // data).  The other possible error is adjacent identical vertices, but
    // this can't happen because S2Builder does not generate such polylines.
    //
    // DEFAULT: S2Debug::ALLOW
    S2Debug s2debug_override() const;
    void set_s2debug_override(S2Debug override);

   private:
    S2Builder::EdgeType edge_type_;
    PolylineType polyline_type_;
    DuplicateEdges duplicate_edges_;
    SiblingPairs sibling_pairs_;
    bool validate_;
    S2Debug s2debug_override_;
  };

  // Specifies that a vector of polylines should be constructed using the
  // given options.
  explicit S2PolylineVectorLayer(
      std::vector<std::unique_ptr<S2Polyline>>* polylines,
      Options const& options = Options());

  // Specifies that a vector of polylines should be constructed using the
  // given options, and that any labels attached to the input edges should be
  // returned in "label_set_ids" and "label_set_lexicion".
  //
  // The labels associated with the edge "polyline[i].vertex({j, j+1})" can be
  // retrieved as follows:
  //
  //   for (int32 label : label_set_lexicon.id_set(label_set_ids[i][j])) {...}
  using LabelSetIds = std::vector<std::vector<LabelSetId>>;
  S2PolylineVectorLayer(std::vector<std::unique_ptr<S2Polyline>>* polylines,
                        LabelSetIds* label_set_ids,
                        IdSetLexicon* label_set_lexicon,
                        Options const& options = Options());

  // Layer interface:
  GraphOptions graph_options() const override;
  void Build(Graph const& g, S2Error* error) override;

 private:
  void Init(std::vector<std::unique_ptr<S2Polyline>>* polylines,
            LabelSetIds* label_set_ids, IdSetLexicon* label_set_lexicon,
            Options const& options);

  std::vector<std::unique_ptr<S2Polyline>>* polylines_;
  LabelSetIds* label_set_ids_;
  IdSetLexicon* label_set_lexicon_;
  Options options_;
};

// Like S2PolylineVectorLayer, but adds the polylines to an S2ShapeIndex.
class IndexedS2PolylineVectorLayer : public S2Builder::Layer {
 public:
  using Options = S2PolylineVectorLayer::Options;
  explicit IndexedS2PolylineVectorLayer(S2ShapeIndex* index,
                                        Options const& options = Options())
      : index_(index), layer_(&polylines_, options) {}

  GraphOptions graph_options() const override {
    return layer_.graph_options();
  }

  void Build(Graph const& g, S2Error* error) override {
    layer_.Build(g, error);
    if (error->ok()) {
      for (auto& polyline : polylines_) {
        index_->Add(new S2Polyline::OwningShape(std::move(polyline)));
      }
    }
  }

 private:
  S2ShapeIndex* index_;
  std::vector<std::unique_ptr<S2Polyline>> polylines_;
  S2PolylineVectorLayer layer_;
};

// A layer type that collects degenerate edges as points.
// This layer expects all edges to be degenerate. In case of finding
// non-degenerate edges it sets S2Error but it still generates the
// output with degenerate edges.
class S2PointVectorLayer : public S2Builder::Layer {
 public:
  class Options {
   public:
    using DuplicateEdges = GraphOptions::DuplicateEdges;
    Options();
    explicit Options(DuplicateEdges duplicate_edges);

    // DEFAULT: DuplicateEdges::MERGE
    DuplicateEdges duplicate_edges() const;
    void set_duplicate_edges(DuplicateEdges duplicate_edges);

   private:
    DuplicateEdges duplicate_edges_;
  };

  explicit S2PointVectorLayer(std::vector<S2Point>* points,
                              Options const& options = Options());

  using LabelSetIds = std::vector<LabelSetId>;
  S2PointVectorLayer(std::vector<S2Point>* points, LabelSetIds* label_set_ids,
                     IdSetLexicon* label_set_lexicon,
                     Options const& options = Options());

  // Layer interface:
  GraphOptions graph_options() const override;
  void Build(Graph const& g, S2Error* error) override;

 private:
  std::vector<S2Point>* points_;
  LabelSetIds* label_set_ids_;
  IdSetLexicon* label_set_lexicon_;
  Options options_;
};

// Like S2PointVectorLayer, but adds the points to an S2ShapeIndex.
class IndexedS2PointVectorLayer : public S2Builder::Layer {
 public:
  using Options = S2PointVectorLayer::Options;
  explicit IndexedS2PointVectorLayer(S2ShapeIndex* index,
                                     Options const& options = Options())
      : index_(index), layer_(&points_, options) {}

  GraphOptions graph_options() const override {
    return layer_.graph_options();
  }

  void Build(Graph const& g, S2Error* error) override {
    layer_.Build(g, error);
    if (error->ok()) {
      index_->Add(new s2shapeutil::PointVectorShape(&points_));
    }
  }

 private:
  S2ShapeIndex* index_;
  std::vector<S2Point> points_;
  S2PointVectorLayer layer_;
};

#if 0
// A layer type that assembles edges into a polygon mesh.  A polygon mesh
// represents a subdivision of the sphere into faces, where each face is
// bounded by one or more edge loops.  Depending on the data you are starting
// with, you may want to turn on one or more of the following options:
//
//  - If the input edges can intersect, then use split_crossing_edges().
//
//  - If you want each input edge to be mapped into two edges (one in each
//    direction) then use EdgeType::UNDIRECTED.  For example, if the input
//    consists of polylines then you would use undirected edges, whereas if
//    the input consisted of a collection of loops that share boundaries
//    (e.g. the US states) you would use directed edges.  (Some edges may
//    exist in only one direction, in which case an edge will automatically be
//    created in the other direction.  You can distinguish such edges by using
//    labels, since automatically created edges will not have any labels.)
class S2MeshLayer;

// Similar to S2PolygonMeshLayer, except that the polygons are returned as a
// vector of s2shapeutil::LaxPolygon rather than an S2Mesh.  LaxPolygon is
// similar to S2Polygon except that it has more relaxed rules about duplicate
// vertices and edges (which is why it is being used here).
//
// This has the same effect as building an S2Mesh and extracting the faces,
// but it is faster and uses less memory.
using LaxPolygonVector = std::vector<std::unique_ptr<s2shapeutil::LaxPolygon>>;
class LaxPolygonVectorLayer;
#endif


//////////////////   Implementation details follow   ////////////////////


inline S2PolygonLayer::Options::Options()
    : edge_type_(S2Builder::EdgeType::DIRECTED), validate_(false) {
}

inline S2PolygonLayer::Options::Options(S2Builder::EdgeType edge_type)
    : edge_type_(edge_type), validate_(false) {
}

inline S2Builder::EdgeType S2PolygonLayer::Options::edge_type() const {
  return edge_type_;
}

inline void S2PolygonLayer::Options::set_edge_type(
    S2Builder::EdgeType edge_type) {
  edge_type_ = edge_type;
}

inline bool S2PolygonLayer::Options::validate() const {
  return validate_;
}

inline void S2PolygonLayer::Options::set_validate(bool validate) {
  validate_ = validate;
}

inline S2PolylineLayer::Options::Options()
    : edge_type_(S2Builder::EdgeType::DIRECTED), validate_(false) {
}

inline S2PolylineLayer::Options::Options(S2Builder::EdgeType edge_type)
    : edge_type_(edge_type), validate_(false) {
}

inline S2Builder::EdgeType S2PolylineLayer::Options::edge_type() const {
  return edge_type_;
}

inline void S2PolylineLayer::Options::set_edge_type(
    S2Builder::EdgeType edge_type) {
  edge_type_ = edge_type;
}

inline bool S2PolylineLayer::Options::validate() const {
  return validate_;
}

inline void S2PolylineLayer::Options::set_validate(bool validate) {
  validate_ = validate;
}

inline S2PolylineVectorLayer::Options::Options()
    : edge_type_(S2Builder::EdgeType::DIRECTED),
      polyline_type_(PolylineType::PATH),
      duplicate_edges_(DuplicateEdges::KEEP),
      sibling_pairs_(SiblingPairs::KEEP),
      validate_(false),
      s2debug_override_(S2Debug::ALLOW) {
}

inline S2PolylineVectorLayer::Options::Options(S2Builder::EdgeType edge_type)
    : edge_type_(edge_type),
      polyline_type_(PolylineType::PATH),
      duplicate_edges_(DuplicateEdges::KEEP),
      sibling_pairs_(SiblingPairs::KEEP),
      validate_(false),
      s2debug_override_(S2Debug::ALLOW) {
}

inline S2Builder::EdgeType S2PolylineVectorLayer::Options::edge_type() const {
  return edge_type_;
}

inline void S2PolylineVectorLayer::Options::set_edge_type(
    S2Builder::EdgeType edge_type) {
  edge_type_ = edge_type;
}

inline S2PolylineVectorLayer::Options::PolylineType
S2PolylineVectorLayer::Options::polyline_type() const {
  return polyline_type_;
}

inline void S2PolylineVectorLayer::Options::set_polyline_type(
    PolylineType polyline_type) {
  polyline_type_ = polyline_type;
}

inline S2PolylineVectorLayer::Options::DuplicateEdges
S2PolylineVectorLayer::Options::duplicate_edges() const {
  return duplicate_edges_;
}

inline void S2PolylineVectorLayer::Options::set_duplicate_edges(
    DuplicateEdges duplicate_edges) {
  duplicate_edges_ = duplicate_edges;
}

inline S2PolylineVectorLayer::Options::SiblingPairs
S2PolylineVectorLayer::Options::sibling_pairs() const {
  return sibling_pairs_;
}

inline void S2PolylineVectorLayer::Options::set_sibling_pairs(
    SiblingPairs sibling_pairs) {
  DCHECK(sibling_pairs == SiblingPairs::KEEP ||
         sibling_pairs == SiblingPairs::DISCARD);
  sibling_pairs_ = sibling_pairs;
}

inline bool S2PolylineVectorLayer::Options::validate() const {
  return validate_;
}

inline void S2PolylineVectorLayer::Options::set_validate(bool validate) {
  validate_ = validate;
  set_s2debug_override(S2Debug::DISABLE);
}

inline S2Debug S2PolylineVectorLayer::Options::s2debug_override() const {
  return s2debug_override_;
}

inline void S2PolylineVectorLayer::Options::set_s2debug_override(
    S2Debug s2debug_override) {
  s2debug_override_ = s2debug_override;
}

inline S2PointVectorLayer::Options::Options()
    : duplicate_edges_(DuplicateEdges::MERGE) {}

inline S2PointVectorLayer::Options::Options(DuplicateEdges duplicate_edges)
    : duplicate_edges_(duplicate_edges) {}

inline S2Builder::GraphOptions::DuplicateEdges
S2PointVectorLayer::Options::duplicate_edges() const {
  return duplicate_edges_;
}

inline void S2PointVectorLayer::Options::set_duplicate_edges(
    S2Builder::GraphOptions::DuplicateEdges duplicate_edges) {
  duplicate_edges_ = duplicate_edges;
}

}  // namespace s2builderutil

#endif  // S2_S2BUILDERUTIL_LAYERS_H_
