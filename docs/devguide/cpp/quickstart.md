---
title: Quick Start
---

Suppose that we have a few terabytes (or petabytes) of geographic data, and we
want to be able to query it efficiently.  This brief tutorial shows how the S2
library is useful for solving this problem.

## Prerequisites

Before going through this Quickstart, first makes sure your development
environment satisfies our
[platform requirements](/about/platforms) and
follow the installation instructions in the
[S2 Install Guide](/about/install).

## S2PointIndex

To keep things simple, let's suppose that our input data consists only of points
(rather than polylines, polygons, etc), and furthermore let's suppose that we
have only a few million points, so that they all fit in memory.  (We will look
at the more general case below.)

The easiest way to build an in-memory point index is to use an `S2PointIndex`,
as follows:

    // Build an index containing random points anywhere on the Earth.
    S2PointIndex<int> index;
    for (int i = 0; i < FLAGS_num_index_points; ++i) {
      S2Point point = S2Testing::RandomPoint();
      index.Add(point, i);
    }

This example uses `S2Testing::RandomPoint()` to generate points that are
randomly distributed over the Earth. Note that in the S2 library, [points are
represented as unit-length
vectors](/about/overview#unit-vectors-vs-latitudelongitude-pairs)
(the [`S2Point`](/devguide/basic_types#s2point) class),
rather than as latitude-longitude pairs. You can convert latitude-longitude
pairs to `S2Point` values using the
[`S2LatLng`](/devguide/basic_types#s2latlng) class, like
this:

    S2Point point(S2LatLng::FromDegrees(lat_degrees, lng_degrees));

Also notice that each point in the index has some auxiliary data associated
with it (in this case, the integer `i`).  You can use this feature to map the
points back to other data structures.

## S2ClosestPointQuery

Now that we have our index, how do we query it?  For example, suppose that we
want to find all the points within 100km of a target point.  This can be done
using `S2ClosestPointQuery`, which provides methods that find the closest
point(s) to a given target object (such as a point, polyline, or polygon).
You can either find the *k* closest points to the target, or all points within
a given radius of the target, or both.

In this example, we will find all points within 100km of a point target.
First we construct the `S2ClosestPointQuery` object:

    S2ClosestPointQuery<int> query(&index);
    query.mutable_options()->set_max_distance(
        S1Angle::Radians(S2Earth::KmToRadians(FLAGS_query_radius_km)));

Notice that we need to convert the radius in kilometers to an `S1Angle`.  This
is because S2 models the Earth as *unit sphere* (a sphere with radius 1),
which means that the distance between two points is the same as the angle
between those points (in radians) measured from the sphere's center.  S2
provides helper function for converting distances to angles, such as the
`S2Earth::KmToRadians` function used above.  (Note that because the Earth is
not quite spherical, distances computed in this way can have errors of up to
0.56%.  You can use a geodesy library to compute distances more accurately,
but spherical distances are perfectly adequate for many purposes, such as
finding nearby restaurants.)

Finally, we can find the points within the given radius of our target point
like this:

    S2ClosestPointQuery<int>::PointTarget target(S2Testing::RandomPoint());
    auto results = query.FindClosestPoints(&target);

This returns a vector of `Result` objects, where each result includes the
following fields:

    S2Point point() const;          // The indexed point.
    Data data() const;              // The auxiliary data for this point.
    S1ChordAngle distance() const;  // The distance to this point.

In our case `data()` is the integer value that we passed to `index.Add()`
above.

Putting this all together, our program for indexing and querying points looks
like this (see `examples/point_index.cc`):

```c++
  // Build an index containing random points anywhere on the Earth.
  S2PointIndex<int> index;
  for (int i = 0; i < FLAGS_num_index_points; ++i) {
    index.Add(S2Testing::RandomPoint(), i);
  }

  // Create a query to search within the given radius of a target point.
  S2ClosestPointQuery<int> query(&index);
  query.mutable_options()->set_max_distance(
      S1Angle::Radians(S2Earth::KmToRadians(FLAGS_query_radius_km)));

  // Repeatedly choose a random target point, and count how many index points
  // are within the given radius of that point.
  int64 num_found = 0;
  for (int i = 0; i < FLAGS_num_queries; ++i) {
    S2ClosestPointQuery<int>::PointTarget target(S2Testing::RandomPoint());
    num_found += query.FindClosestPoints(&target).size();
  }

  printf("Found %lld points in %d queries\n", num_found, FLAGS_num_queries);
```

## Indexing for Search

Let's suppose again that we have a lot of geographic data, and we want to
query it efficiently along with other attributes (such as names and
descriptions) using an *information retrieval system*.  Such systems index a
collection of *documents*, where each document consists of text and other
attributes.  Each document is converted into a collection of *index terms*
(e.g., significant words or phrases), which are gathered into an *inverted
index* that maps each term to a list of documents where that term occurs.
(For a comprehensive introduction to information retrieval, see [Introduction
to Information Retrieval](https://nlp.stanford.edu/IR-book/){:target="_blank"}
by Christopher D. Manning et al., Cambridge University Press, 2008.)

Our goal is not to solve the information retrieval problem in this tutorial,
but rather just to show how to add spatial data to an existing system. S2
defines a special class called `S2RegionTermIndexer` that is designed to deal
with the problem of converting spatial data into index terms, which can then
be indexed along with the other document information. (In fact the method
described here can be used with any key-value storage system.)

To keep things simple, let's assume that each document consists of a single
point (with no other information).  Each document is assigned a "document id"
that can be used to retrieve it later:

    std::vector<S2Point> documents;
    for (int docid = 0; docid < FLAGS_num_documents; ++docid) {
      documents.push_back(S2Testing::RandomPoint());
    }

As a substitute for the "inverted index" of an information retrieval system,
we will use a hash map.  The key of the hash map is an index term, and the
value is the set of document ids where this term is present:

    std::unordered_map<string, std::vector<int>> index;

We can now use `S2RegionTermIndexer` to index these documents as follows:

    // Create an indexer with appropriate options.
    S2RegionTermIndexer::Options options;
    S2RegionTermIndexer indexer(options);

    // Add the documents to the index.
    static const char kPrefix[] = "s2:";
    for (int docid = 0; docid < documents.size(); ++docid) {
      S2Point index_region = documents[docid];
      for (const auto& term : indexer.GetIndexTerms(index_region, kPrefix)) {
        index[term].push_back(docid);
      }
    }

For each document, the indexer returns a set of terms that should be
associated with the document for later retrieval.  Each term is simply an
alphanumeric string.  To distinguish these spatial indexing terms from other
index terms (such as words), the interface allows a prefix to be added (in
this case "s2:").  For example, an index term might look like "s2:34f14a9".

## Querying for Search

Let's now examine how to query our inverted index.  Again, for simplicity
let's assume that we want to retrieve all documents that are within a given
radius of a test point.  As before, the first step is to convert the query
radius to an `S1Angle`:

    S1Angle radius =
        S1Angle::Radians(S2Earth::KmToRadians(FLAGS_query_radius_km));

Our query region is a disc centered around a target point.  On the sphere, a
disc-shaped region such as this one is called a *spherical cap*
([`S2Cap`](/devguide/basic_types#s2cap)):

    S2Cap query_region(S2Testing::RandomPoint(), radius);

Now we convert the query region to a set of *query terms*.  The query terms
are constructed such that if the query region intersects the document region,
then the query terms and index terms are guaranteed to have at least one value
in common.  This means that to find the set of intersecting documents, we
simply need to look up the documents associated with each query term and
compute their union.  (Actual information retrieval systems do something more
sophisticated, but that doesn't concern us here.)

    std::set<int> candidates;
    for (const auto& term : indexer.GetQueryTerms(query_region, kPrefix)) {
      candidates.insert(index[term].begin(), index[term].end());
    }

"candidates" now contains all documents that intersect the query region, along
with some documents that nearly intersect it.  At this point, an information
retrieval system would "score" documents by examining their content more
closely.  In our case, for example, we might want to filter the candidates by
retrieving the original "document" and checking the distance to the target
more precisely. We can do this as follows:

    std::vector<int> result;
    for (int docid : candidates) {
      if (!query_region.Contains(documents[docid])) continue;
      result.push_back(docid);
    }

Now the `results` vector contains exactly the set of input documents that
intersect the query region.

## Optimization for Points

One final detail: when a spatial index contains only points, like this one, it
turns out that fewer query terms are needed than when the index also contains
non-point regions (such as polylines or polygons).  You can take advantage of
this optimization by modifying the code above as follows:

    S2RegionTermIndexer::Options options;
    options.set_index_contains_points_only(true);
    S2RegionTermIndexer indexer(options);

This reduces the number of query terms by approximately a factor of two, but
can only be used when an index contains only points.

The complete example program looks like this:

```c++
  // Create a set of "documents" to be indexed.  Each document consists of a
  // single point.
  std::vector<S2Point> documents;
  for (int docid = 0; docid < FLAGS_num_documents; ++docid) {
    documents.push_back(S2Testing::RandomPoint());
  }

  // We use a hash map as our inverted index.  The key is an index term, and
  // the value is the set of "document ids" where this index term is present.
  std::unordered_map<string, std::vector<int>> index;

  // Create an indexer suitable for an index that contains points only.
  S2RegionTermIndexer::Options options;
  options.set_index_contains_points_only(true);
  S2RegionTermIndexer indexer(options);
  static const char kPrefix[] = "s2:";

  // Add the documents to the index.
  for (int docid = 0; docid < documents.size(); ++docid) {
    S2Point index_region = documents[docid];
    for (const auto& term : indexer.GetIndexTerms(index_region, kPrefix)) {
      index[term].push_back(docid);
    }
  }

  // Convert the query radius to an angle representation.
  S1Angle radius =
      S1Angle::Radians(S2Earth::KmToRadians(FLAGS_query_radius_km));

  // Count the number of documents (points) found in all queries.
  int64 num_found = 0;
  for (int i = 0; i < FLAGS_num_queries; ++i) {
    // Choose a random center for query.
    S2Cap query_region(S2Testing::RandomPoint(), radius);

    // Convert the query region to a set of terms, and compute the union of
    // the document ids associated with those terms.
    std::set<int> candidates;
    for (const auto& term : indexer.GetQueryTerms(query_region, kPrefix)) {
      candidates.insert(index[term].begin(), index[term].end());
    }

    // "candidates" now contains all documents that intersect the query
    // region, along with some documents that nearly intersect it.  We can
    // prune the results by retrieving the original "document" and checking
    // the distance more precisely.
    std::vector<int> result;
    for (int docid : candidates) {
      if (!query_region.Contains(documents[docid])) continue;
      result.push_back(docid);
    }
    // Now do something with the results (in this example we just count them).
    num_found += result.size();
  }
  printf("Found %lld points in %d queries\n", num_found, FLAGS_num_queries);
```

## Other S2Region Types

As its name implies, `S2RegionTermIndexer` supports indexing and querying any
type of `S2Region`.  In addition to points and discs (`S2Cap`), other useful
`S2Region` types include:

*   [`S2LatLngRect`](/devguide/basic_types#s2latlngrect):
    a rectangle in latitude-longitude coordinates.
*   [`S2Polyline`](/devguide/basic_types#s2polyline): a
    polyline.
*   [`S2Polygon`](/devguide/basic_types#s2polygon): a
    polygon (can have holes and multiple shells).
*   [`S2CellUnion`](/devguide/s2cell_hierarchy#s2cellunion):
    a region approximated as a collection of `S2CellIds`.
*   `S2ShapeIndexRegion` - an arbitrary collection of points,
    polylines, and polygons.
*   `S2ShapeIndexBufferedRegion` - like the above, but expanded by a
    given radius.
*   `S2RegionUnion` - the union of arbitrary other regions.
*   `S2RegionIntersection` - the intersection of arbitrary other regions.

For example, you could use `S2RegionTermIndexer` to index a set of polylines,
and then query which polylines intersect a given polygon.

## Comparing Indexing Classes

This tutorial gave a brief introduction to two useful indexing classes
(`S2PointIndex` and `S2RegionTermIndexer`).  Here is a brief comparison of
their features:

`S2PointIndex`:

*   Represents a collection of points in memory.
*   Supports incremental updates.
*   Can use `S2ClosestPointQuery` to find the points near a given target:
    *   Target types include points, polylines, polygons, and collections.
    *   Can also restrict results to a given `S2Region`.

`S2RegionTermIndexer`:

*   A tool for adding spatial data to an information retrieval system.
*   Can add one or more `S2Region` objects to each document.
*   Queries return all documents that intersect a given `S2Region`.
*   Requires an external system for storing and retrieving index terms.

Another useful indexing class not described here is
[`S2ShapeIndex`](/devguide/s2shapeindex), which indexes an
arbitrary collection of points, polylines and polygons (collectively known as
*shapes*). `S2ShapeIndex` is the most generally useful of these classes for
working with geometry in memory.  It compares to the classes above as follows:

`S2ShapeIndex`:

*   Represents a collection of points, polylines, and polygons in memory.
*   Supports incremental updates.
*   Useful query classes include:
    *   `S2ContainsPointQuery`: returns the shape(s) that contain a given
        point.
    *   [`S2ClosestEdgeQuery`](/devguide/s2closestedgequery):
        returns the closest edges to a given point, polyline, polygon, or
        geometry collection.
    *   `S2CrossingEdgeQuery`: returns the edge(s) that cross a given edge.
    *   `S2BooleanOperation`: computes boolean operations such as union, and
        boolean predicates such as containment.
    *   `S2ShapeIndexRegion`: allows geometry collections to be approximated.
    *   `S2ShapeIndexBufferedRegion`: computes approximations that have been
        expanded by a given radius.
