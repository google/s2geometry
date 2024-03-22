// Copyright 2021 Google Inc. All Rights Reserved.
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


#ifndef S2_S2CHAIN_INTERPOLATION_QUERY_H_
#define S2_S2CHAIN_INTERPOLATION_QUERY_H_

#include <memory>
#include <vector>

#include "s2/s1angle.h"
#include "s2/s2point.h"
#include "s2/s2shape.h"

// S2ChainInterpolationQuery is a helper class for querying points on S2Shape's
// edges by spherical distance.  The distance is computed cumulatively along the
// edges contained in the shape, using the order in which the edges are stored
// by the S2Shape object.
//
// If a particular chain id is specified at the query initialization, then the
// distance values are computed along that single chain, which allows per-chain
// interpolation.  If no chain is specified, then the interpolated point as a
// function of distance is discontinuous at chain boundaries.  Using multiple
// chains can be used in such algorithms as quasi-random sampling along the
// total span of a multiline.
//
// Once the query object is initialized, the complexity of each subsequent query
// is O( log(number of edges) ).  The complexity of the initialization and the
// memory footprint of the query object are both O(number of edges).
class S2ChainInterpolationQuery {
 public:
  S2ChainInterpolationQuery() = default;

  // Creates an instance of S2ChainInterpolationQuery with the given shape
  // and chain id.  See Init() method below for the usage of chain_id.
  explicit S2ChainInterpolationQuery(const S2Shape* shape, int chain_id = -1);

  // Destructor.
  virtual ~S2ChainInterpolationQuery() = default;

  // Initializes the instance of S2ChainInterpolationQuery with the given shape
  // and chain id.
  //
  // If a non-negative chain_id is supplied by the caller, then only the edges
  // belonging to the chain with that chain_id are used, otherwise the edges
  // from all chains are used in the order they are stored in the shape object.
  void Init(const S2Shape* shape, int chain_id = -1);

  // Gets the total length of the chain(s), which corresponds to the distance at
  // the end vertex of the last edge of the chain(s).  Returns zero length for
  // shapes containing no edges.
  S1Angle GetLength() const;

  // Returns the cumulative length along the edges being interpolated up to the
  // end of the given edge ID.  Uses the same edge ID numbering as the
  // S2ChainInterpolationQuery::Result.  Returns S1Angle::Infinity() if the edge
  // ID does not lie within the set of edges being interpolated. Returns
  // S1Angle::Zero() if the S2ChainInterpolationQuery is empty.
  S1Angle GetLengthAtEdgeEnd(int edge_id) const;

  // Result is the container class for the results of S2ChainInterpolationQuery
  // operations.  When valid, it contains the actual interpolated S2Point, the
  // id of the edge to which the point belongs, and the actual distance to the
  // point along the chain(s).
  class Result {
   public:
    // Creates an invalid result object.
    Result() : valid_(false) {}

    // Creates a valid result from the point, edge id and distance.
    Result(const S2Point& point, int edge_id, const S1Angle& distance)
        : valid_(true), point_(point), edge_id_(edge_id), distance_(distance) {}

    // Tells whether the result is valid. An invalid result can be returned
    // only by a query that is either uninitialized, or initialized with a shape
    // containing no edges.
    bool is_valid() const { return valid_; }

    // Returns the point which is the result of a query.  The point is defined
    // iff the Result object is valid.
    S2Point point() const { return point_; }

    // Returns the index of the edge on which the point is located.  The edge id
    // is defined iff the Result object is valid.
    int edge_id() const { return edge_id_; }

    // Returns the actual distance corresponding to the resulting point.  It may
    // differ from the input distance if the latter exceeds the maximum
    // parameter value (total length) of the shape's chain(s), in which case the
    // resulting distance is snapped to total length.  The distance value is
    // valid iff the Result object is valid.
    S1Angle distance() const { return distance_; }

   private:
    bool valid_;
    S2Point point_;
    int edge_id_;
    S1Angle distance_;
  };

  // Computes the S2Point located at the given distance along the edges from the
  // first vertex of the first edge.  The resulting S2Point is an approximation
  // of the true position with the error bound by kGetPointOnLineError.  Also
  // computes the edge id and the actual distance corresponding to the resulting
  // point.
  //
  // This method returns a valid Result iff the query has been initialized with
  // at least one edge.
  //
  // If the input distance exceeds the total length, then the resulting point is
  // the end vertex of the last edge, and the resulting distance is set to the
  // total length.
  //
  // If there are one or more degenerate (zero-length) edges corresponding to
  // the given distance, then the resulting point is located on the first of
  // these edges.
  Result AtDistance(const S1Angle& distance) const;

  // Similar to the above function, but takes the normalized fraction of the
  // distance as input, with fraction = 0 corresponding to the beginning of the
  // shape or chain and fraction = 1 to the end.  Forwards the call to
  // AtDistance().  A small precision loss may occur due to converting the
  // fraction to distance by multiplying it by the total length.
  Result AtFraction(double fraction) const;

  // Returns the vector of points that is a slice of the chain from
  // begin_fraction to end_fraction. If begin_fraction is greater than
  // end_fraction, then the points are returned in reverse order.
  //
  // For example, Slice(0,1) returns the entire chain, Slice(0, 0.5) returns the
  // first half of the chain, and Slice(1, 0.5) returns the second half of the
  // chain in reverse.
  //
  // The endpoints of the slice are interpolated (except when coinciding with an
  // existing vertex of the chain), and all the internal points are copied from
  // the chain as is.
  //
  // If the query is either uninitialized, or initialized with a shape
  // containing no edges, then an empty vector is returned.
  std::vector<S2Point> Slice(double begin_fraction, double end_fraction) const;

  // Appends the chain slice from begin_fraction to end_fraction to the given
  // slice. If begin_fraction is greater than end_fraction, then the points are
  // appended in reverse order. If the query is either uninitialized, or
  // initialized with a shape containing no edges, then no points are appended.
  void AddSlice(double begin_fraction, double end_fraction,
                std::vector<S2Point>& slice) const;

 private:
  const S2Shape* shape_;
  std::vector<S1Angle> cumulative_values_;
  int first_edge_id_;
  int last_edge_id_;
};

#endif  // S2_S2CHAIN_INTERPOLATION_QUERY_H_
