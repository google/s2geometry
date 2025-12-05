---
title: Interpolation along Chains
---

The purpose of this document is to provide a description of
[`S2ChainInterpolationQuery`](https://source.corp.google.com/piper///depot/google3/util/geometry/s2chain_interpolation_query.h)
and to illustrate its usage with code snippets.

## Overview

[`S2ChainInterpolationQuery`](https://source.corp.google.com/piper///depot/google3/util/geometry/s2chain_interpolation_query.h)
is a tool for interpolating point locations along chains. For example, it can
compute the position that is 50% along a polyline, or generate all the
positions on a polygon boundary that are 1e-4 radians apart. The interpolation
is based on spherical distance that is computed cumulatively along the edges of
the chain. It works with geometry classes that are derived from
[`S2Shape`](https://source.corp.google.com/piper///depot/google3/util/geometry/s2shape.h)
and represent polylines, polygons, or sets thereof.

An `S2Shape` is organized as a collection of edges, which are identified by a
contiguous range of integer edge ids, starting at 0. The edges are grouped into
chains, where each chain consists of a sequence of edges connected end-to-end
and thus forming a polyline, or a closed polygon. For example, a shape
representing a single open polyline with _N_ vertices would consist of just one
chain containing _N - 1_ edges. A shape representing an island with a lake in it
would consist of two chains, one chain for the outer contour of the island and
the other for the contour of the lake.

If the input shape contains multiple chains, the caller can specify the index of
the chain to be used, or let the query use all chains contained in the shape.
When multiple chains are used, they are assumed to be arranged end to end in the
distance space, in the same order as they are stored in the `S2Shape`, see
[Working with Multiple Chains](https://g3doc.corp.google.com/devguide/s2chain_interpolation_query.md#working-with-multiple-chains)
section for more details.

`S2ChainInterpolationQuery` can be used for:

  * Finding the point that is located at a given distance from the first vertex
   along the chain
  * Finding the point that is located at a given normalized distance, that is a
   fraction of the polyline's total length.
  * Obtaining the edge ids of the resulting points

## Initializing an `S2ChainInterpolationQuery`

An instance of `S2ChainInterpolationQuery` can be initialized as follows:

```
const S2Shape* shape = ...;
int chain_id = ...;

S2ChainInterpolationQuery query(shape, chain_id);
```

where `chain_id` is the non-negative id of the shape's chain on which the points
are going to be queried (0 <= `chain_id` < shape->num_chains()).

If `chain_id` is negative (or altogether omitted) then all shape's chains are
going to be used in the order of their respective ids, as if there were a single
polyline. Note that in this case the resulting point may not be a continuous
function of distance, since the beginning of a chain may not coincide with the
end of the previous chain.

When the input shape is known to contain only one chain, the `chain_id`
parameter can be either set to 0 or omitted.

## Using the Query Methods

Once an instance of `S2ChainInterpolationQuery` is initialized, the total length
of a chain (or chains) can be accessed via its `GetLength()` method:

```
S1Angle total_length = query.GetLength();
```

Note that if the query is uninitialized, or initialized with an empty shape
(i.e. a shape that contains no chains), the returned total length is zero
(`S1Angle::Zero()`.

To find the point at a given spherical distance _d_ along the chain(s), use the
`AtDistance()` method. This method also computes the id of the edge on which the
resulting point is located, as well as the actual parametric distance of the
point. The resulting point, edge id and distance can be obtained via the
respective accessor methods of `S2ChainInterpolationQuery::Result`:

```
S1Angle d = ...;
S2ChainInterpolationQuery::Result result = query.AtDistance(d);

if (result.is_valid()) {}
  S2Point point_at_distance = result.point();
  int edge_id = result.edge_id();
  S1Angle actual_distance = result.distance();
}
```

If the input distance _d_ is less than or equal to the chain's total length,
then the actual distance returned by `result.distance()` equals _d_. If,
however, _d_ is greater than the chain's total length, then the resulting point
is snapped to the last vertex of the chain(s) and `result.distance()` returns
the total length of the chain.

For example, consider the polyline _abc_ formed by three vertices _A_, _B_ and
_C_ with the respective latitudes/longitudes of (0, 0), (0, 1) and (0, 2). It
has two edges, _AB_ and _BC_, each of them having the spherical length of 1
degree. For this polyline we have:

```
// Vertices A, B and C
auto vertices = s2textformat::ParsePointsOrDie("0:0, 0:1, 0:2"); 

S2LaxPolylineShape shape(vertices);

// The above shape is known to contain just one chain (A,B,C), therefore one can
// omit the chain_id parameter in the query initialization.
S2ChainInterpolationQuery query(&shape);

// The total length is 2 degrees = sum of lengths of edges (A,B) and (B,C).
S1Angle total_length = query.GetLength();

// Find the point at the spherical distance of 1.5 degrees.
S2ChainInterpolationQuery::Result result = query.AtDistance(S1Angle::Degrees(1.5)));
S2Point point = result.point();
S1Angle distance = result.distance();
int edge_id = result.edge_id();
```

The resulting point is located half-way between vertices *B* and *C*, the
resulting distance is the same as the input distance of 1.5 degrees and the edge
id is 1, corresponding to the second edge of the chain. Alternatively, calling
`query.AtDistance(2.5)` yields the resulting point coinciding with *C*, a
resulting distance equal to the chain's total length (2 degrees), and an edge id
of 1, since the results of queries with distances exceeding the chain's length
are snapped to the last vertex of the chain.

To find the point at a given normalized distance (i.e. at a fraction of the
total length) _f_, where 0 <= _f_ <= 1, use the `AtFraction()` method:

```
double f = ...;  

auto result = query.AtFraction(f);
S2Point point_at_distance = result.point();
int edge_id = result.edge_id();
S1Angle actual_distance = result.distance();
```

The actual distance returned by `AtFraction()` is always within the interval of
[0, _total_length_]. If passed an input fraction value greater than 1, then
the resulting point is snapped to the last vertex of the chain(s) and the
distance returned by `AtFraction()` equals the total length of the chain.

Calling `query.AtFraction(0.75)` for the polyline shape in the above
example would yield the resulting point half-way between _B_ and _C_, resulting
distance of 1.5 and the edge id = 1 (the last edge of the chain). And
`query.AtFraction(1.5)` gives the resulting point coinciding with the
last vertex _C_, resulting length = 2 (the total length) and the resulting edge
id = 1 (the edge corresponding to the last vertex).

Both `AtDistance()` and `AtFraction()` return valid results for all input
distance values if and only if the query has been initialized with a shape
containing at least one valid edge.

Let us now consider an
[`S2LaxPolygonShape`](https://source.corp.google.com/piper///depot/google3/util/geometry/s2lax_polygon_shape.h)
containing two chains, each being 1 degree long and containing two edges.

If the query for this shape is initialized with a `chain_id` of 0:

```
S2ChainInterpolationQuery query(&shape, 0);
```

then only the first chain (chain_id = 0) is considered and the other chain is
ignored. In this case `query.GetLength()` returns the value corresponding to the
angle of 1 degree, which is the length of the single chain.

Now, if the query is initialized with a `chain_id` that is negative (or
omitted), both chains are considered as if they were part of the same polyline,
and `query.GetLength()` returns an angle of 2 degrees, which is the length of
both chains combined:

```
S2ChainInterpolationQuery query(&shape);
S1Angle total_length = query.GetLength();  // total_length.degrees() == 2
```

In this case calling `query.AtDistance(1.5)` or `query.AtFraction(0.75)` would
yield a point halfway down the second chain, and the resulting an edge id = 1
(the second edge of the shape corresponding to the first and the only edge of
the second chain).

## Working with Multiple Chains

An important caveat when considering multiple chain is that the resulting point
is not a continuous function of the input distance, since the last vertex of a
chain may not coincide with the first vertex of the next chain. In the above
case, calling `query.AtFraction( 0.5 )` yields the resulting point = last vertex
of the first chain, while `query.AtFraction( 0.5 + epsilon )` yields the first
vertex of the second chain, which may be located far apart from the last vertex
of the previous chain.

A typical use case for using queries with multiple chains is when one needs to
generate an *approximately* uniform sampling of all chains contained in a shape.
In these cases one can just call `query.AtFraction(t)` with *t* ranging from 0.0
to 1.0 with a certain step, depending on the desired sampling density. Note that
such sampling is only approximately uniform, since the chain lengths won't
necessarily be multiples of the step size.

## Computational Complexity

*   `S2ChainInterpolationQuery` initialization has the computational complexity
    of *O(N)*, where *N* is the total number of edges.
*   `AtDistance()` and `AtFraction()` calls have the complexity of *O( log(N)
    )*.
*   `GetLength()` call has the constant complexity *O(1)*.
*   The memory footprint of the query is *O(N)*.

## Numerical Accuracy

The `S2Point` returned by `result.point()` accessor is an approximation of the
true position with the error bound by
[`kGetPointOnLineError`](https://source.corp.google.com/piper///depot/google3/util/geometry/s2edge_distances.h;l=148?q=kGetPointOnLineError).

The computation error of the total length is bound by 3.25 * DBL_EPSILON * _N_,
where _N_ is the number of edges. An additional loss of precision may occur due
to floating point multiplication of `fraction * total_length` during
`AtFraction()` call.
