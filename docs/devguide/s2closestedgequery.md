---
title: Finding Nearby Edges
---

## Overview

`S2ClosestEdgeQuery` is a tool for finding the closest edge(s) between two
geometries.  By using the appropriate options, this class can answer questions
such as:

 - Find the minimum distance between two geometries *A* and *B*.

 - Find all edges of geometry *A* that are within a distance *d* of geometry
   *B*.

 - Find the *k* edges of geometry *A* that are closest to a given point *P*.

The input geometries may consist of any number of points, polylines,
and polygons (collectively known as *shapes*).  Shapes do not need to
be disjoint; they may overlap or intersect arbitrarily.  The
implementation is designed to be fast for both simple and complex
geometries (e.g., one edge vs. 100 million edges).

## The Index

`S2ClosestEdgeQuery` operates on two geometries, known as the *index*
and the *target*.  Each query returns the subset of edges in the index
that have a specified relationship to the target (e.g., the closest
edge to the target, or all edges within a distance *d* of the target).

The indexed geometry is provided as an `S2ShapeIndex` object, which is
simply a collection of points, polylines, and/or polygons (possibly
intersecting, as mentioned above).  Here is a simple example showing
how to build an index containing several polylines:

    vector<S2Polyline*> polylines = ...;
    S2ShapeIndex index;
    for (S2Polyline* polyline : polylines) {
      index.Add(new S2Polyline::Shape(polyline));
    }

An `S2ClosestEdgeQuery` object can then be constructed like so:

    S2ClosestEdgeQuery query(&index);

Note that `index` is passed as a const pointer rather than a reference
to reflect the fact that the index must persist for the lifetime of the
query object.

## The Target

The *target* represents the geometry that distances are measured to.
It is specified as an object of type `S2ClosestEdgeQuery::Target`.
There are `Target` subtypes for measuring the distance to a point, an
edge, an `S2Cell`, or an entire `S2ShapeIndex` (which allows finding
closest edges or measuring distances between two arbitrary collections
of geometry).

For example, the following code tests whether any of the polylines in
the example above is closer than `dist_limit` to a given point
`test_point`:

    S2Point test_point = ...;
    S1ChordAngle dist_limit = ...;
    S2ClosestEdgeQuery::PointTarget target(test_point);
    if (query.IsDistanceLess(&target, dist_limit)) { ... }

Note that `target` is passed as a pointer; this is because `Target`
objects may contain temporary data and can therefore only be used by
one query at a time.  (The target and query object are intended for use
only within a single thread, as opposed to the `S2ShapeIndex` object
which can be used in many different threads simultaneously.)

## Query Methods

The fastest and simplest query method is `IsDistanceLess`, which
returns true if the distance to the target is less than a given limit:

    bool IsDistanceLess(Target* target, S1ChordAngle limit);

If you want to test whether the distance is less than or equal to "limit"
instead, you can use `IsDistanceLessOrEqual` instead:

    if (query.IsDistanceLessOrEqual(&target, limit)) { ... }

To compute the actual distance, you can use the `GetDistance` method:

    S1ChordAngle GetDistance(Target* target);

However note that if you only need to compare the distance against some
threshold value, it is much faster to call `IsDistanceLess` instead.

The return type of this method (`S1ChordAngle`) is an efficient
representation of the distance between two points on the unit sphere.
It can easily be converted to an `S1Angle` if desired (which measures
angles in radians or degrees), and from there to an actual distance on
the Earth spheroid.  (See also [Modeling Accuracy](#modeling-accuracy).)

The most general query method is `FindClosestEdges`:

    std::vector<Result> FindClosestEdges(Target* target);

This method finds the closest edge(s) to the target that satisfy the
query options (described in the following section).  Each edge is
returned in the form of a `Result` object that identifies the edge and
specifies the distance to that edge:

    class Result {
      Distance distance() const;  // The distance from the target to this edge.
      int32 shape_id() const;     // Identifies an indexed shape.
      int32 edge_id() const;      // Identifies an edge within the shape.
    };

The `Distance` type is simply a thin wrapper around `S1ChordAngle` that
provides some additional methods needed by the `S2ClosestEdgeQuery`
implementation.  The other two fields identify an edge in the
`S2ShapeIndex` that was passed to the `S2ClosestEdgeQuery`
constructor.  The query object has two convenience methods that may be
useful for processing result edges:

    // Returns the endpoints of the given result edge.
    S2Shape::Edge GetEdge(const Result& result) const;

    // Returns the point on given result edge that is closest to "point".
    S2Point Project(const S2Point& point, const Result& result) const;

For example, to retrieve the endpoints of each edge and also find the
closest point on each result edge to a given target point, you could do
this:

    for (const auto& result : query.FindClosestEdges(&target)) {
      S2Shape::Edge edge = query.GetEdge(result);
      S2Point closest_point = query.Project(point, result);
    }

There is also a convenience method for the common case where only one
edge (the closest edge is needed).  This method returns a single
`Result` object rather than a vector:

    Result FindClosestEdge(Target* target);

## Options

`S2ClosestEdgeQuery` supports a number of options that control which
edges are returned by each query.  These options can be specified as a
second argument to the constructor:

    S2ClosestEdgeQuery::Options options;
    S2ClosestEdgeQuery query(&index, options);

Unlike most classes in the library, `S2ClosestEdgeQuery` also allows
options to be modified between calls to query methods via the
`mutable_options()` accessor method.  For example:

    S1ChordAngle distance1 = query.GetDistance(&target1);
    query.mutable_options()->set_max_distance(dist_limit);
    S1ChordAngle distance2 = query.GetDistance(&target2);

This can be significantly more efficient compared to changing the
options by constructing a new `S2ClosestEdgeQuery` object, because the
query object contains a significant amount of internal state that is
reused from one query to the next.

The most important options are the following:

    // Specifies that at most "max_results" edges should be returned.
    //
    // REQUIRES: max_results >= 1
    // DEFAULT: numeric_limits<int>::max()
    int max_results() const;
    void set_max_results(int max_results);

    // Specifies that only edges whose distance to the target is less than
    // "max_distance" should be returned.
    //
    // Note that edges whose distance is exactly equal to "max_distance" are
    // not returned.  Normally this doesn't matter, because distances are not
    // computed exactly in the first place, but if such edges are needed then
    // see set_inclusive_max_distance() below.
    //
    // DEFAULT: Distance::Infinity()
    void set_max_distance(S1ChordAngle max_distance);

    // Like set_max_distance(), except that edges whose distance is exactly
    // equal to "max_distance" are also returned.  Equivalent to calling
    // set_max_distance(max_distance.Successor()).
    void set_inclusive_max_distance(S1ChordAngle max_distance);

Note that you will always want to set one of these options before
calling `FindClosestEdges`, because otherwise all the edges in the
entire `S2ShapeIndex` will be returned.  (This is not necessary
when calling `FindClosestEdge`, `GetDistance`, or `IsDistanceLess`,
because these methods all implicitly limit max_results() to 1.)

Another important option is `include_interiors()`, which determines
whether distances are measured to the boundary and interior of polygons
or only to their boundaries:

    // Specifies that polygon interiors should be included when measuring
    // distances.
    //
    // DEFAULT: true
    bool include_interiors() const;
    void set_include_interiors(bool include_interiors);

By default this option is true, which means for example that a point
inside a polygon has a distance of zero.  Note that in this situation
there is no "closest edge" to return (since the minimum distance is
attained at a point interior to the polygon), and therefore the
`FindClosestEdge(s)` methods indicate this in the `Result` object(s) by
setting `edge_id` to -1.  (Such results can be identified using the
`is_interior()` method.)  For example:

    for (const auto& result : query.FindClosestEdges(&target)) {
      if (result.is_interior()) {
        ProcessInterior(result.shape_id());
      } else {
        ProcessEdge(result.shape_id, result.edge_id);
      }
    }

The following option can be useful when performance is critical:

    // Specifies that edges up to max_error() further away than the true
    // closest edges may be substituted in the result set, as long as such
    // edges satisfy all the remaining search criteria (such as max_distance).
    //
    // DEFAULT: Distance::Zero()
    Distance max_error() const;
    void set_max_error(Distance max_error);

This option give the algorithm permission to stop the search early when
the best possible improvement in distance drops below `max_error()`.
For example, this option is used internally by the `IsDistanceLess`
predicate to stop the search as soon as it finds any edge whose
distance is less than the given limit, rather than continuing to search
for an edge that is even closer (which would be wasted effort).

Note that this options only has an effect if `max_results()` is also
specified; otherwise all edges closer than `max_distance()` will always
be returned.

## Target Types

The following target types exist in addition to the `PointTarget` type
described earlier.

`EdgeTarget` computes the closest distance to a spherical geodesic edge.  For
example, to find the closest indexed edge to a target edge `ab`, you could do
this:

    S2ClosestEdgeQuery::EdgeTarget target(a, b);
    auto result = query.FindClosestEdge(&target);

`CellTarget` computes the distance to an `S2Cell`, including the interior of the
cell.  This target type is used extensively internally to implement other
algorithms.  For example, to test whether the distance to `cell` is less than
`limit`:

    S2ClosestEdgeQuery::CellTarget target(cell);
    if (query.IsDistanceLess(&target, limit)) { ... }

Finally, `ShapeIndexTarget` computes the distance to an `S2ShapeIndex`
(containing an arbitrary collection of points, polylines, and polygons, possibly
overlapping).  For example, to calculate the distance between two
`S2ShapeIndexes`, you can do this:

    S2ClosestEdgeQuery query(&index1);
    S2ClosestEdgeQuery::ShapeIndexTarget target(&index2);
    S1ChordAngle distance = query.GetDistance(&target);

(As always, if you only want to compare the distance against a threshold value
then it is much more efficient to call `IsDistanceLess` instead.)

The `ShapeIndexTarget` type also supports its own `Options`.  In particular, if
you want to measure the distance to the boundaries of polygons in the target
index (ignoring the polygon interiors), you can call

    target.set_include_interiors(false);

before executing the distance query.  Note that the `include_interiors()` option
on the query and the target are independent, so for example you can measure
distance from the boundary of one polygon to the boundary and interior of
another.

## Modeling Accuracy

This class models the Earth as a sphere, which can lead to distance errors of up
to 0.56% compared with modeling the Earth as an ellipsoid.  While this is
perfectly adequate for many purposes (e.g., determining nearby restaurants),
some applications need higher accuracy.

This can be achieved by recomputing distances using using an external geodesy
library such as `geographiclib`.  For exact results, all distances would need to
be computed this way (and in fact there is a templatized base class
`S2ClosestEdgeQueryBase` that allows this), but in most cases it is sufficient
to use spherical distances to determine which pair of points is closest, and
then recalculate the distance accurately for that pair of points only.  (This is
because the spherical distance error varies relatively slowly over the Earth's
surface.)  For example, if the target is a point then you could do the
following:

    S2ClosestEdgeQuery::PointTarget target(target_point);
    auto result = query.FindClosestEdge(&target);
    double meters;
    if (result.is_empty()) {
      // The index is empty.
      meters = std::numeric_limits<double>::infinity();
    } else if (result.is_interior()) {
      // The target point is inside an indexed shape.
      meters = 0.0;
    } else {
      // Find the point on the result edge that is closest to the target.
      S2Point closest_point = query.Project(target_point, result);
      meters = GeoidDistance(target_point, closest_point);
    }

The same technique can be used for more complex targets (`ShapeIndexTarget`) by
executing a second query to find the closest pair of edges between the two
geometry collections, finding the points on that edge pair that achieve the
minimum spherical distance (see `S2::GetEdgePairClosestPoints`), and then
calculating the ellipsoidal distance between those points:

    S2ClosestEdgeQuery query1(&index1);
    S2ClosestEdgeQuery::ShapeIndexTarget target2(&index2);
    auto result1 = query1.FindClosestEdge(&target2);
    double meters;
    ... like the code above ...
    } else {
      // Get the edge from index1 (edge1) that is closest to index2.
      S2Shape::Edge edge1 = query1.GetEdge(result1);

      // Now find the edge from index2 (edge2) that is closest to edge1.
      S2ClosestEdgeQuery query2(&index2);
      S2ClosestEdgeQuery::EdgeTarget target1(edge1.v0, edge1.v1);
      auto result2 = query2.FindClosestEdge(&target1);
      S2Shape::Edge edge2 = query2.GetEdge(result2);

      // Find the closest point pair on edge1 and edge2.
      auto closest = S2::GetEdgePairClosestPoints(edge1.v0, edge1.v1,
                                                  edge2.v0, edge2.v1);
      meters = GeoidDistance(closest.first, closest.second);
    }

## Numerical Accuracy

The spherical geodesic distance calculation used by this class is fast and
accurate, but not exact.  Expressing the errors in terms of distances on the
Earth's surface, for distances up to 10,000 km the error is less than 6
nanometers.  However for points that are nearly antipodal, the error can be as
large as 50 cm (dropping off to less than 1 micrometer for points that are 50 km
away from being antipodal).

If you would like to measure distances conservatively by including these error
bounds, you can use the following option:

    options.set_conservative_max_distance(limit);

This will automatically increase `limit` by the maximum error in the distance
calculation, ensuring that all edges whose true, exact distance is less than or
equal to `limit` are returned (along with some edges whose true distance is
slightly greater).

Algorithms that need to do exact distance comparisons can use this option to
compute a set of candidate edges that can then be filtered further using exact
predicates (such as `s2pred::CompareEdgeDistance`).

There is also a conservative version of the `IsDistanceLessOrEqual` predicate,

    bool IsConservativeDistanceLessOrEqual(Target* target, S1ChordAngle limit);

which automatically increases the given limit by the maximum error in the
distance calculation.
