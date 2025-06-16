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

#ifndef S2_S2SHAPE_H_
#define S2_S2SHAPE_H_

#include <cstdint>
#include <iterator>

#include "absl/log/absl_log.h"

#include "s2/base/types.h"
#include "s2/util/coding/coder.h"
#include "s2/s2coder.h"
#include "s2/s2point.h"
#include "s2/s2pointutil.h"
#include "s2/util/coding/coder.h"

// The purpose of S2Shape is to represent polygonal geometry in a flexible
// way.  It is organized as a collection of edges that optionally defines an
// interior.  All geometry represented by an S2Shape must have the same
// dimension, which means that an S2Shape can represent either a set of
// points, a set of polylines, or a set of polygons.
//
// S2Shape is defined as an abstract base class in order to give clients
// control over the underlying data representation.  Sometimes an S2Shape does
// not have any data of its own, but instead "wraps" some other class.  There
// are various useful subtypes defined in *_shape.h, and some S2 classes also
// have a nested "Shape" class (e.g., S2Polygon::Shape).  It is easy for
// clients to implement their own subtypes, since the interface is minimal.
//
// S2Shape operations are typically defined on S2ShapeIndex objects rather
// than individual shapes.  An S2ShapeIndex is simply a collection of
// S2Shapes, possibly of different dimensions (e.g. 10 points and 3 polygons),
// organized into a data structure for efficient edge access.
//
// The edges of an S2Shape are identified by a contiguous range of "edge ids"
// starting at 0.  The edges are further subdivided into "chains", where each
// chain consists of a sequence of edges connected end-to-end (a polyline).
// For example, an S2Shape representing two polylines AB and CDE would have
// three edges (AB, CD, DE) grouped into two chains: (AB) and (CD, DE).
// Similarly, an S2Shape representing 5 points would have 5 chains consisting
// of one edge each.
//
// S2Shape has methods that allow edges to be accessed either using the global
// numbering (edge id) or within a particular chain.  The global numbering is
// sufficient for most purposes, but the chain representation is useful for
// certain algorithms such as intersection (see S2BooleanOperation).
class S2Shape {
 public:
  // An edge, consisting of two vertices "v0" and "v1".  Zero-length edges are
  // allowed, and can be used to represent points.
  struct Edge {
    S2Point v0, v1;
    constexpr Edge() = default;
    constexpr Edge(const S2Point& _v0, const S2Point& _v1) : v0(_v0), v1(_v1) {}

    // Returns the edge with the vertices reversed.
    Edge Reversed() const { return {v1, v0}; }

    // Returns true if the edge is degenerate.
    bool IsDegenerate() const { return v0 == v1; }

    // Returns true if point equals v1, indicating this edge is arriving.
    bool Incoming(const S2Point& point) const { return v1 == point; }

    // Returns true if point equals v0, indicating this edge is leaving.
    bool Outgoing(const S2Point& point) const { return v0 == point; }

    // Returns true if point is one of the vertices of this edge.
    bool IncidentOn(const S2Point& point) const {
      return Incoming(point) || Outgoing(point);
    }

    // TODO(ericv): Define all 6 comparisons.
    friend bool operator==(const Edge& x, const Edge& y) {
      return x.v0 == y.v0 && x.v1 == y.v1;
    }

    friend bool operator!=(const Edge& x, const Edge& y) { return !(x == y); }

    friend bool operator<(const Edge& x, const Edge& y) {
      return x.v0 < y.v0 || (x.v0 == y.v0 && x.v1 < y.v1);
    }

    template <typename H>
    friend H AbslHashValue(H h, const Edge& e) {
      return H::combine(std::move(h), e.v0, e.v1);
    }
  };

  // A range of edge ids corresponding to a chain of zero or more connected
  // edges, specified as a (start, length) pair.  The chain is defined to
  // consist of edge ids {start, start + 1, ..., start + length - 1}.
  struct Chain {
    int32_t start, length;
    Chain() = default;
    constexpr Chain(int32_t _start, int32_t _length)
        : start(_start), length(_length) {}

    friend bool operator==(const Chain& x, const Chain& y) {
      return x.start == y.start && x.length == y.length;
    }
  };

  // The position of an edge within a given edge chain, specified as a
  // (chain_id, offset) pair.  Chains are numbered sequentially starting from
  // zero, and offsets are measured from the start of each chain.
  struct ChainPosition {
    int32_t chain_id, offset;
    ChainPosition() = default;
    ChainPosition(int32_t _chain_id, int32_t _offset)
        : chain_id(_chain_id), offset(_offset) {}

    friend bool operator==(const ChainPosition& x, const ChainPosition& y) {
       return x.chain_id == y.chain_id && x.offset == y.offset;
    }
  };

  // A ReferencePoint consists of a point P and a boolean indicating whether P
  // is contained by a particular shape.
  struct ReferencePoint {
    S2Point point;
    bool contained;
    ReferencePoint() = default;
    ReferencePoint(S2Point _point, bool _contained)
        : point(_point), contained(_contained) {}

    // Returns a ReferencePoint with the given "contained" value and a default
    // "point".  It should be used when all points or no points are contained.
    static ReferencePoint Contained(bool _contained) {
      return ReferencePoint(S2::Origin(), _contained);
    }

    friend bool operator==(const ReferencePoint& x, const ReferencePoint& y) {
      return x.point == y.point && x.contained == y.contained;
    }
  };

  // A 32-bit tag that can be used to identify the type of an encoded S2Shape.
  // All encodable types have a non-zero type tag.  The tag associated with a
  // given shape type can be accessed as Shape::kTypeTag, while the tag
  // associated with a given object can be accessed as shape.type_tag().
  //
  // Type tags in the range 0..8191 are reserved for use by the S2 library.
  using TypeTag = uint32_t;

  // Indicates that a given S2Shape type cannot be encoded.
  static constexpr TypeTag kNoTypeTag = 0;

  // The following constant should be updated whenever new types are added.
  static constexpr TypeTag kNextAvailableTypeTag = 6;

  // The minimum allowable tag for user-defined S2Shape types.
  static constexpr TypeTag kMinUserTypeTag = 8192;

  S2Shape() = default;
  virtual ~S2Shape() = default;

  // Returns the number of edges in this shape, or points, if the shape's
  // dimension is 0.  Edges or points have ids ranging from 0 to
  // num_edges() - 1.
  virtual int num_edges() const = 0;

  // Returns the edge or point for the given edge id.  Points are represented as
  // degenerate edges, with equal endpoints, but not all degenerate edges are
  // points.
  //
  // REQUIRES: 0 <= id < num_edges()
  virtual Edge edge(int edge_id) const = 0;

  // Returns the dimension of the geometry represented by this shape.
  //
  //  0 - Point geometry.  Each point is represented as a degenerate edge.
  //
  //  1 - Polyline geometry.  Polyline edges may be degenerate.  A shape may
  //      represent any number of polylines.  Polylines edges may intersect.
  //
  //  2 - Polygon geometry.  Edges should be oriented such that the polygon
  //      interior is always on the left.  In theory the edges may be returned
  //      in any order, but typically the edges are organized as a collection
  //      of edge chains where each chain represents one polygon loop.
  //      Polygons may have degeneracies (e.g., degenerate edges or sibling
  //      pairs consisting of an edge and its corresponding reversed edge).
  //      A polygon loop may also be full (containing all points on the
  //      sphere); by convention this is represented as a chain with no edges.
  //      (See S2LaxPolygonShape for details.)
  //
  // Note that this method allows degenerate geometry of different dimensions
  // to be distinguished, e.g. it allows a point to be distinguished from a
  // polyline or polygon that has been simplified to a single point.
  virtual int dimension() const = 0;

  // Returns true if the shape contains no points.  (Note that the full
  // polygon is represented as a chain with zero edges.)
  bool is_empty() const {
    return num_edges() == 0 && (dimension() < 2 || num_chains() == 0);
  }
  // Returns true if the shape contains all points on the sphere and has no
  // edges.
  bool is_full() const {
    return num_edges() == 0 && dimension() == 2 && num_chains() > 0;
  }

  // Returns an arbitrary point P along with a boolean indicating whether P is
  // contained by the shape.  (The boolean value must be false for shapes that
  // do not have an interior.)
  //
  // This ReferencePoint may then be used to compute the containment of other
  // points by counting edge crossings.
  virtual ReferencePoint GetReferencePoint() const = 0;

  // Returns the number of contiguous edge chains in the shape.  For example,
  // a shape whose edges are [AB, BC, CD, AE, EF] would consist of two chains
  // (AB,BC,CD and AE,EF).  Every chain is assigned a "chain id" numbered
  // sequentially starting from zero.
  //
  // Note that it is always acceptable to implement this method by returning
  // num_edges() (i.e. every chain consists of a single edge), but this may
  // reduce the efficiency of some algorithms.
  virtual int num_chains() const = 0;

  // Returns the range of edge ids corresponding to the given edge chain.  The
  // edge chains must form contiguous, non-overlapping ranges that cover the
  // entire range of edge ids.  This is spelled out more formally below:
  //
  // REQUIRES: 0 <= i < num_chains()
  // REQUIRES: chain(i).length >= 0, for all i
  // REQUIRES: chain(0).start == 0
  // REQUIRES: chain(i).start + chain(i).length == chain(i+1).start,
  //           for i < num_chains() - 1
  // REQUIRES: chain(i).start + chain(i).length == num_edges(),
  //           for i == num_chains() - 1
  virtual Chain chain(int chain_id) const = 0;

  // Returns the edge at offset "offset" within edge chain "chain_id".
  // Equivalent to "shape.edge(shape.chain(chain_id).start + offset)"
  // but may be more efficient.
  virtual Edge chain_edge(int chain_id, int offset) const = 0;

  // Finds the chain containing the given edge, and returns the position of
  // that edge as a (chain_id, offset) pair.
  //
  // REQUIRES: shape.chain(pos.chain_id).start + pos.offset == edge_id
  // REQUIRES: shape.chain(pos.chain_id + 1).start > edge_id
  //
  // where     pos == shape.chain_position(edge_id).
  virtual ChainPosition chain_position(int edge_id) const = 0;

  // Returns an integer that can be used to identify the type of an encoded
  // S2Shape (see TypeTag above).
  virtual TypeTag type_tag() const { return kNoTypeTag; }

  // Appends an encoded representation of the S2Shape to "encoder".  Note that
  // the encoding does *not* include the type_tag(), so the tag will need to
  // be encoded separately if the shape type will be unknown at decoding time
  // (see s2shapeutil::EncodeTaggedShapes() and related functions).
  //
  // The encoded representation should satisfy the following:
  //
  //  - It should include a version number, so that the encoding may be changed
  //    or improved in the future.
  //
  //  - If "hint" is CodingHint::FAST, the encoding should be optimized for
  //    decoding performance.  If "hint" is CodingHint::COMPACT, the encoding
  //    should be optimized for space.
  //
  // REQUIRES: "encoder" uses the default constructor, so that its buffer
  //           can be enlarged as necessary by calling Ensure(int).
  virtual void Encode(Encoder* encoder, s2coding::CodingHint hint) const {
    ABSL_DLOG(FATAL) << "Encoding not implemented for this S2Shape type";
  }

  // Virtual methods that return pointers of your choice.  These methods are
  // intended to help with the problem of attaching additional data to S2Shape
  // objects.  For example, you could return a pointer to a source object, or
  // a pointer to a bundle of additional data allocated with the S2Shape.
  // Because this method exists in all S2Shapes, you can override it in each
  // type of shape you have and call it without knowing the concrete subtype.
  // For example, if you have polyline and polygon shapes, you can do this:
  //
  //   class MyPolyline : public S2Polyline::Shape {
  //    public:
  //     virtual void* mutable_user_data() { return &my_data_; }
  //    private:
  //     MyData my_data_;
  //   };
  //   class MyPolygon : public S2Polygon::Shape {
  //    public:
  //     virtual void* mutable_user_data() { return &my_data_; }
  //    private:
  //     MyData my_data_;
  //   };
  //   ...
  //   S2Shape* shape = index.shape(id);
  //   const MyData* data = static_cast<const MyData*>(shape->user_data());
  //
  // This is not the only way to map from an S2Shape back to your source
  // data.  Other reasonable techniques include:
  //
  //  - If all of your shapes are the same type, then you can create your own
  //    subclass of some existing S2Shape type (such as S2Polyline::Shape) and
  //    add your own methods and fields.  You can access this data by
  //    downcasting the S2Shape pointers returned by S2ShapeIndex methods.
  virtual const void* user_data() const { return nullptr; }
  virtual void* mutable_user_data() { return nullptr; }

  // ChainVertexIterator allows the use of iterator syntax for accessing
  // vertices of a shape's chain, e.g.:
  //
  //   S2Shape::Chain chain = shape->chain(chain_id);
  //   for (const S2Point& p : shape.vertices(chain)) { ... }
  // or
  //   int chain_id = ...
  //   for (const S2Point& p : shape.vertices(chain_id)) { ... }
  //
  // ChainVertexIterator supports dereference, increment and equality/inequality
  // operators.
  //
  // A ChainVertexIterator is valid as long as the shape it has been created
  // from exists. It doesn't own the shape it points to, so destroying it has no
  // effect on the shape.
  class ChainVertexIterator {
   public:
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;
    using pointer = const S2Point*;
    using reference = const S2Point&;
    using value_type = S2Point;

    ChainVertexIterator() = default;

    // Creates iterator pointing to the chain vertex with given index. The index
    // is zero-based vertex offset from the beginning of the chain.
    ChainVertexIterator(const S2Shape* shape, const Chain& chain,
                        int vertex_index);

    // Dereference operators.
    reference operator*() const;
    pointer operator->() const;

    // Prefix and postfix increment operators.
    ChainVertexIterator& operator++();
    ChainVertexIterator operator++(int);

    // Equality operator overload.
    bool operator==(ChainVertexIterator it) const;

    // Inquality operator overload.
    bool operator!=(ChainVertexIterator it) const;

   private:
    const S2Shape* shape_ = nullptr;
    Chain chain_;

    // Offset of the current vertex in the chain.
    int vertex_index_ = 0;

    // Cached copy of the edge used to access the current vertex.
    Edge edge_;

    // Offset of the current edge in the chain.
    int edge_offset_ = -1;

    // Offset of the current vertex in the current edge (0 for v0, 1 for v1).
    int edge_vertex_ = 0;

    // Updates the cached edge data if necessary, to make sure that the current
    // edge contains the current vertex.
    void UpdateCurrentEdge();
  };

  // ChainVertexRange provides access to the iterable vertex range of the given
  // chain:
  //
  //   S2Shape::ChainVertexRange vertices(shape, chain);
  //   for (const S2Point& p : vertices) { ... }
  //
  // A ChainVertexRange is valid as long as the shape it has been created from
  // exists. It doesn't own the shape it points to, so destroying it has no
  // effect on the shape.
  class ChainVertexRange {
   public:
    ChainVertexRange() = default;

    // Creates the vertex range for the given chain of the shape.
    ChainVertexRange(const S2Shape* shape, const Chain& chain);

    // Returns the vertex iterator pointing to the first vertex of the chain.
    ChainVertexIterator begin() const;

    // Returns the vertex iterator pointing to the end of the chain.
    ChainVertexIterator end() const;

    // Returns the number of vertices in the chain.
    // Note: in case of polyline shapes (shape->dimension() == 1) the number of
    // vertices differs from the number of edges by 1.
    int num_vertices() const;

   private:
    const S2Shape* shape_ = nullptr;
    Chain chain_;
  };

  // Returns the iterable vertex range for given chain.
  ChainVertexRange vertices(const Chain& chain) const;

  // Returns the iterable vertex range for given chain id.
  ChainVertexRange vertices(int chain_id) const;

  // ChainIterator allows iterating over the chains of the shape using
  // range-based for loops:
  //
  //   for (const S2Shape::Chain& chain : shape) { ... }
  //
  // A ChainIterator is valid as long as the shape it has been created from
  // exists. It doesn't own the shape it points to, so destroying it has no
  // effect on the shape.
  class ChainIterator {
   public:
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::forward_iterator_tag;
    using pointer = Chain*;
    using reference = Chain;
    using value_type = Chain;

    ChainIterator() = default;

    // Creates the iterator pointing to the shapes' chain with given chain id.
    ChainIterator(const S2Shape* shape, int chain_id)
        : shape_(shape), chain_id_(chain_id) {}

    // Dereference operator.
    reference operator*() const;

    // Prefix and postfix increment operators.
    ChainIterator& operator++();
    ChainIterator operator++(int);

    // REQUIRES: "it" and *this must reference the same S2Shape.
    bool operator==(ChainIterator it) const;

    // REQUIRES: "it" and *this must reference the same S2Shape.
    bool operator!=(ChainIterator it) const;

   private:
    const S2Shape* shape_ = nullptr;
    int chain_id_ = 0;
  };

  // ChainRange provides access to the iterable range of chains of the given
  // shape:
  //
  //   S2Shape::ChainRange chains(shape);
  //   for (const S2Shape::Chain& chain : chains) { ... }
  //
  // or:
  //
  //   for (const S2Shape::Chain& chain : shape->chains()) { ... }
  //
  // A ChainRange is valid as long as the shape it has been created from exists.
  // It doesn't own the shape it points to, so destroying it has no effect on
  // the shape.
  class ChainRange {
   public:
    ChainRange() = default;

    // Creates an instance of ChainRange for given shape.
    explicit ChainRange(const S2Shape* shape) : shape_(shape) {}

    // Returns the chain iterator pointing to the first chain of the shape.
    ChainIterator begin() const;

    // Returns the chain iterator pointing to the end of the chain range.
    ChainIterator end() const;

   private:
    const S2Shape* shape_ = nullptr;
  };

  // Returns the chain range of the shape.
  ChainRange chains() const;

 protected:
  // S2Shape has some state used by the S2ShapeIndex classes.  If we want to be
  // able to move or copy derived classes, we need to have those operations
  // available on S2Shape as well.  Having them defined however presents the
  // risk of accidental slicing as in:
  //
  //     S2Shape& a = S2LaxPolygon(...);
  //     S2Shape& b = S2LaxPolygon(...);
  //     a = b; // <-- Does not call S2LaxPolygon's copy assignment.
  //
  // So we make these protected to allow Derived classes to decide on their
  // own move/copy semantics without exposing them to broader use.
  S2Shape(S2Shape&&) = default;
  S2Shape(const S2Shape&) = default;
  S2Shape& operator=(const S2Shape&) = default;
  S2Shape& operator=(S2Shape&&) = default;

 private:
  friend class EncodedS2ShapeIndex;
  friend class MutableS2ShapeIndex;
};

//////////////////   Implementation details follow   ////////////////////

inline S2Shape::ChainIterator::reference S2Shape::ChainIterator::operator*()
    const {
  return shape_->chain(chain_id_);
}

inline S2Shape::ChainIterator& S2Shape::ChainIterator::operator++() {
  ++chain_id_;
  return *this;
}

inline S2Shape::ChainIterator S2Shape::ChainIterator::operator++(int) {
  return ChainIterator(shape_, chain_id_++);
}

inline bool S2Shape::ChainIterator::operator==(ChainIterator it) const {
  return chain_id_ == it.chain_id_ && shape_ == it.shape_;
}

inline bool S2Shape::ChainIterator::operator!=(ChainIterator it) const {
  return !(*this == it);
}

inline S2Shape::ChainRange S2Shape::chains() const { return ChainRange(this); }

inline S2Shape::ChainIterator S2Shape::ChainRange::begin() const {
  return ChainIterator(shape_, 0);
}

inline S2Shape::ChainIterator S2Shape::ChainRange::end() const {
  return ChainIterator(shape_, shape_->num_chains());
}

inline S2Shape::ChainVertexIterator::ChainVertexIterator(const S2Shape* shape,
                                                         const Chain& chain,
                                                         int vertex_index)
    : shape_(shape), chain_(chain), vertex_index_(vertex_index) {
  UpdateCurrentEdge();
}

inline S2Shape::ChainVertexIterator::reference
S2Shape::ChainVertexIterator::operator*() const {
  return edge_vertex_ ? edge_.v1 : edge_.v0;
}

inline S2Shape::ChainVertexIterator::pointer
S2Shape::ChainVertexIterator::operator->() const {
  return edge_vertex_ ? &edge_.v1 : &edge_.v0;
}

inline S2Shape::ChainVertexIterator&
S2Shape::ChainVertexIterator::operator++() {
  ++vertex_index_;
  UpdateCurrentEdge();
  return *this;
}

inline S2Shape::ChainVertexIterator S2Shape::ChainVertexIterator::operator++(
    int) {
  auto result = *this;
  ++vertex_index_;
  UpdateCurrentEdge();
  return result;
}

inline bool S2Shape::ChainVertexIterator::operator==(
    S2Shape::ChainVertexIterator it) const {
  return vertex_index_ == it.vertex_index_ && chain_.start == it.chain_.start &&
         shape_ == it.shape_;
}

inline bool S2Shape::ChainVertexIterator::operator!=(
    S2Shape::ChainVertexIterator it) const {
  return !(*this == it);
}

// This function makes sure that cached copy of the edge contains the vertex
// at the current vertex index. The cached copy is updated only if needed.
// Except for the last edge of even-sized chains, only the edges with even
// offsets are cached.
//
// We want to cache the current edge for two reasons:
//  (1) use only even edges to cut the edge generation in half;
//  (2) be able to return pointer/reference to S2Point in ChainVertexIterator.
inline void S2Shape::ChainVertexIterator::UpdateCurrentEdge() {
  // No need to update in case of the iterator pointing beyond the end of the
  // range.
  if (vertex_index_ > chain_.length) {
    return;
  }

  // Compute the new edge offset.
  int edge_offset = 0;

  // For point shapes (dim == 0) the edge offset is always 0.
  if (shape_->dimension() > 0) {
    // An edge with even offset (0, 2, ...) can provide access to two vertices
    // at once.
    edge_offset = vertex_index_ % 2 ? vertex_index_ - 1 : vertex_index_;
    if (edge_offset < chain_.length) {
      // Vertices with odd indices are accessed via v1, and those with even
      // indices via v0 of the even edge.
      edge_vertex_ = vertex_index_ % 2;
    } else {
      // A corner case: the last vertex of a chain with even number of edges
      // can be accessed only via the last edge (which has an odd offset).
      edge_offset = chain_.length - 1;
      edge_vertex_ = 1;
    }
  }
  // Update the cached copy of the edge only if needed.
  if (edge_offset_ != edge_offset) {
    edge_offset_ = edge_offset;
    edge_ = shape_->edge(chain_.start + edge_offset_);
  }
}

inline S2Shape::ChainVertexRange::ChainVertexRange(const S2Shape* shape,
                                                   const Chain& chain)
    : shape_(shape), chain_(chain) {}

inline S2Shape::ChainVertexIterator S2Shape::ChainVertexRange::begin() const {
  return ChainVertexIterator(shape_, chain_, 0);
}

inline S2Shape::ChainVertexIterator S2Shape::ChainVertexRange::end() const {
  return ChainVertexIterator(shape_, chain_, num_vertices());
}

inline int S2Shape::ChainVertexRange::num_vertices() const {
  // The number of vertices for point shapes (dimension = 0) is always 1, and it
  // is equal to the chain length, which is also = 1 for the point shapes.
  //
  // For the polygon shapes (dimension = 2) the number of vertices is equal to
  // the number of edges which is given by the chain length. For example, a
  // rectangle has 4 edges and 4 vertices.
  //
  // For polyline shape (dimension = 1) the number of vertices equals the number
  // of edges + 1 (one vertex at the start of each edge, plus one exta vertex
  // at the end of the last edge), hence chain_.length + 1.
  return shape_->dimension() == 1 ? chain_.length + 1 : chain_.length;
}

inline S2Shape::ChainVertexRange S2Shape::vertices(
    const S2Shape::Chain& chain) const {
  return ChainVertexRange(this, chain);
}

inline S2Shape::ChainVertexRange S2Shape::vertices(int chain_id) const {
  return vertices(chain(chain_id));
}

#endif  // S2_S2SHAPE_H_
