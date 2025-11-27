---
title: Overview
---

S2 is a library for spherical geometry that aims to have the same
robustness, flexibility, and performance as the very best planar
geometry libraries.

## Spherical Geometry

Let's break down the elements of this goal.  First, why spherical
geometry?  (By the way, the name "S2" is derived from the mathematical
notation for the unit sphere, *S²*.)

Traditional cartography is based on *map projections*, which are simply
functions that map points on the Earth's surface to points on a planar
map.  Map projections create distortions due to the fact that the shape
of the Earth is not very close to the shape of a plane.  For example,
the well-known Mercator projection is discontinuous along the 180 degree
meridian, has large scale distortions at high latitudes, and cannot
represent the north and south poles at all.  Other projections make
different compromises, but no planar projection does a good job of
representing the entire surface of the Earth.

S2 approaches this problem by working exclusively with *spherical
projections*.  As the name implies, spherical projections map points on
the Earth's surface to a perfect mathematical sphere.  Such mappings
still create some distortion, of course, because the Earth is not quite
spherical--but as it turns out, the Earth is much closer to being a
sphere than a plane.  With spherical projections, it is possible to
approximate the entire Earth's surface with a maximum distortion of
0.56%.  Perhaps more importantly, spherical projections preserve the
correct topology of the Earth -- there are no singularities or
discontinuities to deal with.

Why not project onto an ellipsoid?  (The Earth isn't quite ellipsoidal
either, but it is even closer to being an ellipsoid than a sphere.)
The answer relates to the other goals stated above, namely performance
and robustness.  Ellipsoidal operations are still orders of magnitude
slower than the corresponding operations on a sphere.  Furthermore,
robust geometric algorithms require the implementation of exact
geometric predicates that are not subject to numerical errors.  While
this is fairly straightforward for planar geometry, and somewhat harder
for spherical geometry, it is not known how to implement all of the
necessary predicates for ellipsoidal geometry.

## Robustness

Which brings up the question, what do we mean by "robust"?

In the S2 library, the core operations are designed to be 100% robust.
This means that each operation makes strict mathematical guarantees
about its output, and is implemented in such a way that it meets those
guarantees for all possible valid inputs.  For example, if you compute
the intersection of two polygons, not only is the output guaranteed to
be topologically correct (up to the creation of degeneracies), but it
is also guaranteed that the boundary of the output stays within a
user-specified tolerance of true, mathematically exact result.

Robustness is very important when building higher-level algorithms,
since unexpected results from low-level operations can be very
difficult to handle.  S2 achieves this goal using a combination of
techniques from computational geometry, including *conservative error
bounds*, *exact geometric predicates*, and *snap rounding*.

S2 attempts to be precise both in terms of mathematical definitions
(e.g., whether regions include their boundaries, and how degeneracies
are handled) and numerical accuracy (e.g., minimizing cancellation
error).

## Flexibility

S2 is organized as a toolkit that gives clients as much control as
possible.  For example, S2 makes it easy to implement your own geometry
subtypes that store data in a format of your choice (the `S2Shape`
interface).  Similarly, most operations are designed to give clients
fine control over the semantics, e.g. whether polygon boundaries are
considered to be open or closed (or semi-open).  S2 has APIs at several
different levels to accommodate various uses of the library.

## Performance

S2 is designed to have good performance on large geographic datasets.
Most operations are accelerated using an in-memory edge index data
structure (`S2ShapeIndex`).  For example if you have a million
polygons, finding the polygon(s) that contain a given point typically
takes a few hundred nanoseconds.  Similarly it is fast to find objects
that are near each other, such as finding all the places of business
near a given road, or all the roads near a given location.

Many operations on indexed data also have output-sensitive running
times.  For example, if you intersect the Pacific Ocean with a small
rectangle near Seattle, the running time essentially depends on the
complexity of the Pacific Ocean near Seattle, not the complexity of the
Pacific Ocean overall.

## Scope

The S2 library provides the following:

*   Representations of angles, intervals, latitude-longitude points,
    unit vectors, and so on, and various operations on these types.

*   Geometric shapes over the unit sphere, such as spherical caps
    ("discs"), latitude-longitude rectangles, polylines, and polygons.

*   Robust constructive operations (e.g., union) and boolean predicates
    (e.g., containment) for arbitrary collections of points, polylines,
    and polygons.

*   Fast in-memory indexing of collections of points, polylines, and
    polygons.

*   Algorithms for measuring distances and finding nearby objects.

*   Robust algorithms for snapping and simplifying geometry (with
    accuracy and topology guarantees).

*   A collection of efficient yet exact mathematical predicates for
    testing relationships among geometric objects.

*   Support for spatial indexing, including the ability to approximate regions
    as collections of discrete "S2 cells".  This feature makes it easy to
    build large distributed spatial indexes.

*   Extensive testing on Google's vast collection of geographic data.

*   Flexible Apache 2.0 license.

On the other hand, the following are outside the scope of S2:

*   Planar geometry.  (There are many fine existing planar geometry
    libraries to choose from.)

*   Conversions to/from common GIS formats.  (To read such formats, use
    an external library such as
    [OGR](http://gdal.org/1.11/ogr/){:target="_blank"}.)

## Language Support

The [reference implementation of the S2
library](https://github.com/google/s2geometry) is written in C++.  Versions in
other languages are ported from the C++ code, and may not have the same
robustness, performance, or features as the C++ version.  As of November 2017,
the "official" ports are as follows:

*   [S2 library in Go](https://github.com/golang/geo).  This project is being
    actively developed and aims to have full equivalence with the C++ version.

*   [S2 library in Java](https://github.com/google/s2-geometry-library-java).
    Be aware that this port is based on the 2011 release of the S2 library, and
    does not have the same robustness, performance, or features as the current
    C++ implementation.

*   [S2 library in Python](https://github.com/google/s2geometry/tree/master/src/python).
    A subset of the library has been ported to Python using SWIG, and is
    included in the `python` directory of the C++ repository.  The supported
    feature set includes polylines, polygons, discs, rectangles, and `S2Cell`
    coverings.

## Layered Design

The S2 library is structured as a toolkit with various layers.  At the
lowest level, it provides a set of operations on points, edges, and simple
shapes such as rectangles and discs.  For example, this includes classes
such as `S2Point`, `S2CellId`, `S1Angle`, `S1Interval`, `S2LatLng`,
`S2Region`.  It also includes the exact predicates defined in
`s2predicates.h` and `s2edge_crossings.h`, utility functions such as those
in `s2edge_distances.h`, and classes related to the `S2Cell` hierarchy (a
hierarchical decomposition of the sphere) such as `S2CellId`, `S2Region`,
and `S2RegionCoverer`.

The next layer provides flexible support for polygonal geometry, including
point sets and polylines.  It is built around the following core classes:

*   `S2Builder`: a robust tool for constructing, and simplifying
    geometry, built on the framework of snap rounding.

*   `S2ShapeIndex`: an indexed collection of points, polylines, and/or
    polygons, possibly overlapping.

*   `S2BooleanOperation`: evaluates boolean operations and predicates for
    arbitrary collections of polygonal geometry.

*   `S2ClosestEdgeQuery`: measures distances and finds closest edges for
    arbitrary collections of polygonal geometry.

These are only the most important classes, but what they have in common
is that (1) they work with arbitrary collections of points, polylines,
and polygons; (2) they allow clients to control the underlying geometry
representation (via the `S2Shape` interface), and (3) they give clients
fine control over the semantics of operations (e.g., how to handle
degeneracies, whether polygons are open or closed, whether and how the
output should be snapped, etc).

The final layer consists of convenience classes, such as `S2Polygon`,
`S2Loop`, and `S2Polyline`.  These types have wide interfaces and
convenience methods that implement many of the operations above (e.g.,
intersection, distances), however they don't offer the same flexibility
in terms of data representation and control over semantics (e.g.,
polygon boundaries are always semi-open).

For example, you can measure the distance from a point to a polygon as
follows:

    S2Point p = ...;
    S2Polygon a = ...;
    S1Angle dist = a.GetDistance(p);

whereas with the `S2ShapeIndex` layer you would write

    S2Point p = ...;
    S2Polygon a = ...;
    S2ShapeIndex index;
    index.Add(absl::make_unique<S2Polygon::Shape>(&a));
    S2ClosestEdgeQuery query(&index);
    S1ChordAngle dist = query.GetDistance(p);

The first case is simpler, but the second case supports many more
options (e.g. measuring the distance to a collection of geometry, or
finding the objects within a certain distance).

## Design Choices

This section summarizes some of the major design choices made by the library.

### Spherical Geodesic Edges

In S2, all edges are "spherical geodesics", i.e. shortest paths on the
sphere.  For example, the edge between two points 100 meters apart on
opposite side of the North pole goes directly through the North pole.
In contrast, the same edge in the standard "plate carr&eacute;e"
projection (i.e., raw latitude/longitude coordinates) would follow a
line of constant latitude, which in this case happens to be a
semicircular path that always stays exactly 50 meters away from the
North pole.  (Lines of latitude except for the equator are never
shortest paths.)

With geodesic edges there are no special cases near the poles or the
180 degree meridian; for example, if two points 10km apart are
separated by the 180 degree meridian, the edge between them simply
crosses the meridian rather than going all the way around the other
side of the Earth.

### Orientation Matters

In S2, polygon loops contain the region on their left.  This differs
from planar libraries that might define loops as being "clockwise" or
"counter-clockwise", or that accept loops of either orientation.  None
of these concepts really applies to polygons drawn on a sphere; for
example, consider a loop that goes around the equator.  Is it clockwise
or counter-clockwise?  (For those familiar with the term, the concept
of "winding number" does not exist on the sphere.)

For these reasons, S2 follows the "interior is on the left" rule.  (This
is the same as the counter-clockwise rule for small loops, but differs
for very large loops that enclose more than half of the sphere.)

### Unit Vectors vs. Latitude/Longitude Pairs

In S2, points are represented internally as unit-length vectors (points
on the surface of a three-dimensional unit sphere) as opposed to
traditional (latitude, longitude) pairs.  This is for two reasons:

*   Unit vectors are much more efficient when working with geodesic
    edges.  Using (latitude, longitude) pairs would require constantly
    evaluating trigonometric functions (sin, cos, etc), which is slow
    even on modern CPU architectures.

*   Unit vectors allow exact geometric predicates to be evaluated
    efficiently.  To do this with (latitude, longitude) pairs, you
    would essentially need to convert them to the unit vector
    representation first.

For precision and efficiency, coordinates are representated as
double-precision numbers.  Robustness is achieved through the use of
the techniques mentioned previously: conservative error bounds (which
measure the error in certain calculations), exact geometric predicates
(which determine the true mathematical relationship among geometric
objects), and snap rounding (a technique that allows rounding results
to finite precision without creating topological errors).

### Classes vs. Functions

Most algorithms in S2 are implemented as classes rather than functions.
So for example, rather than having a `GetDistance` function, there is
an `S2ClosestEdgeQuery` class with a `GetDistance` method.  This design
has two advantages:

*   It allows caching and reuse of data structures.  For example, if an
    algorithm needs to make thousands of distance queries, it can
    declare one `S2ClosestEdgeQuery` object and call its methods
    thousands of times.  This provides the opportunity to avoid
    recomputing data that can be used from one query to the next.

*   It provides a convenient way to specify options (via a nested
    `Options` class), and to avoid respecifying those options many
    times when repeated operations need to be performed.

### Geocentric vs. Geodetic Coordinates

Geocentric and geodetic coordinates are two methods for defining a
(latitude, longitude) coordinate system over the Earth's surface.  The
question sometimes arises: which system does the S2 library use?

The answer is: neither (or both).  Points in S2 are modeled as lying on
a perfect sphere (just as points in a traditional planar geometry
library are modeled as lying on a perfect plane).  How points are
mapped onto the sphere is outside the scope of S2; it works equally
well with geocentric or geodetic coordinates.  (Actual geographic
datasets use geodetic coordinates exclusively.)

## Authors

The S2 library was written primarily by Eric Veach. Other significant
contributors include Jesse Rosenstock, Eric Engle (Java port lead), Robert
Snedegar (Go port lead), Julien Basch, and Tom Manshreck
