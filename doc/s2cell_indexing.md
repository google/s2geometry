[TOC]

# Spatial Indexing

The [Quick Start](quickstart.md) document gives an introduction to spatial
indexing using the S2 library.  The most important classes for this purposes are
`S2PointIndex`, [`S2ShapeIndex`](s2shapeindex.md), and `S2RegionTermIndexer`.
All of these classes use the `S2Cell` hierarchy to accelerate

This
For indexing geometry in memory, the best approach 
There are many approaches to using the `S2Cell` hierarchy for spatial indexing,
depending on what you want to index (points vs. regions), what operations your
indexing data structure support (range queries vs. key lookups only), and what
sort of queries you want to make (point location vs. region queries). Let's go
over a few examples.

## Indexing Points

First, suppose you have a bunch of points in memory, and you want to repeatedly
find the set of points enclosed by various disc-shaped regions (e.g. all points
within 5km of the north pole). Do the following:

*   Convert each point to a corresponding leaf cell id (`S2CellId::FromPoint`).
*   Put all the cell ids in a vector and sort them.
*   Construct an `S2Cap` corresponding to the query region
    (`S2Cap::FromAxisAngle`).
*   Find a set of cell ids that cover the cap (an `S2CellUnion`) by calling
    `S2RegionCoverer::GetCellUnion`.
*   For each cell "x" in the covering, do a range query on the vector of points
    to find all the leaf cell ids that are contained by "x".
*   If you only want the points that are exactly contained by the disc, rather
    than just a superset of them, test each point against the cap.

Here is some example code:

```c++
vector<S2CellId> index;
for (const S2Point& point : points) {
  index.push_back(S2CellId(point));
}
std::sort(index.begin(), index.end());
S2RegionCoverer coverer;
coverer.mutable_options()->set_max_cells(12);  // Default is 8.
for (S2CellId id : coverer.GetCellUnion(S2Cap(center, radius))) {
  int lo = lower_bound(index.begin(), index.end(), id.range_min());
  int hi = upper_bound(index.begin(), index.end(), id.range_max());
  for (int j = lo; j < hi; ++j) {
    if (cap.Contains(index[j].ToPoint())) {
      LOG(INFO) << "Cap contains point " << index[j].ToPoint();
    }
  }
}
```

Some observations:

*   If you want to add or delete points from the index, use a `set<S2CellId>`
    rather than a vector. Sets have built-in `lower_bound()` and `upper_bound()`
    methods. If you have additional data that you want to associate with each
    cell id, use a `map<S2CellId, MyData*>` or `multimap<S2CellId, MyData*>`.

*   If the points don't all fit in memory, suitable choices for the index
    include `SSTable` and `BigTable`, both of which support efficient range
    queries.

*   The code above can be made slightly more efficient by skipping the
    upper_bound() call and instead advancing sequentially until the id exceeds
    range_max(). To use this approach, you should insert an S2CellId::Sentinel()
    in the index so that you don't also need to check for the end of the vector.
    For example:

```c++
index.push_back(S2CellId::Sentinel());
...
for (int i = 0; i < cells.size(); ++i) {
 int j = lower_bound(index.begin(), index.end(), cells[i].id().range_min());
 int max_id = cells[i].id().range_max();
 for (; index[j] <= max_id; ++j) {
   ...
 }
}
```

## Indexing Regions

Suppose that the objects to be indexed are regions rather than points. The
easiest approach is convert each region to a collection of cell ids, and create
one index entry for each cell id. This increases the index size by a small
constant factor that depends on how accurate an approximation is desired.

To execute a query, convert the query region to a collection of cell ids and do
a range query for each one just as before. However, since index entries now
consist of cells of all sizes, an extra step is required. For each cell id in
the query region, you also need to compute all the ancestors of those cell ids
and do key lookups on those as well:

```c++
typedef multimap<S2CellId, MyRegion*> MyIndex;
MyIndex index;
// Build index similar to before, except that each MyRegion* can
// have several index entries (one for each cell in its covering).
// Build query region and compute its covering as before.

// Add a sentinel to the index so that we don't need to worry about
// hitting index.end().
index[S2CellId::Sentinel()] = NULL;

// New query code:
set<MyRegion*> result;
set<S2CellId> ancestors;
for (int i = 0; i < cells.size(); ++i) {
  S2CellId id = cells[i].id();
  MyIndex::const_iterator j = index.lower_bound(id.range_min());
  S2CellId max_id = id.range_max();
  for (; j->first <= max_id; ++j) result.insert(j->second);
  while (!id.is_face()) {
    id = id.parent();
    if (ancestors.count(id) > 0) break;
    j = index.lower_bound(id);
    // Iterates over entries in index (a multi-map) that match the id.
    for (; j->first == id; ++j) result.insert(j->second);
    ancestors.insert(id);
  }
}
```

This code computes the set of original regions whose coverings intersect the
covering of the query region. The regions can be filtered further if a more
exact intersection test is desired.

## Indexes Without Range Queries

Indexing systems that do not support efficient range queries need a different
alternative for indexing regions, which is to expand the index by inserting all
the ancestors of each cell. To make this efficient, the ancestor cells are
labelled differently than the cells that actually belong to the coverings. In
other words, there are two separate terms (tokens) for each cell id: one for
cells that belong to a covering ("covering" terms), and another for ancestors of
these cells ("ancestor" terms).

To execute a query, we simply compute a covering for the query region, and then
look up all the "ancestor" terms for the cells in this covering, plus all the
"covering" terms for these cells and their ancestors. This effectively divides
the query into two parts: one to find all of the regions whose covering contains
a cell that is descendant of a cell in the query region, and another to find all
the regions whose covering contains a cell that is a ancestor of a cell in the
query region. (Cells that are identical to a cell in the query region are
handled as part of the second group, since we only inserted "covering" terms for
these cells into the index.)

This technique does not require any range queries, and can efficiently compute
overlaps between query and index regions of any size. The total number of
lookups per query is about 4-10 for the base covering (depending on how accurate
an approximation is desired) plus about 10-20 for the ancestors of these cells
(assuming the query region diameter is at least 10m).

If it is important to reduce the number of lookups, here are a few ways to do
this:

*   Set a minimum subdivision level for cell ids. For example, cells at
    subdivision level 10 happen to be about 10km wide. If the index only
    contains cells that are 10km or smaller (i.e. level 10 or higher), then
    ancestors at levels less than 10 can be skipped during the query process
    (eliminating 10 lookups). To use this technique, just call `GetCellUnion`
    normally, and then iterate through the children of any cells that are too
    big as outlined below. Of course, very large regions such as Canada will
    need a lot of cells.

    ```c++
    S2CellId id = cells[i].id();
    int level = max(id.level(), kMinIndexLevel);
    S2CellId end = id.child_end(level);
    for (S2CellId c = id.child_begin(level); c != end; c = c.next()) {
      index[c] = region;
    }
    ```

*   Set a maximum subdivision level for cell ids. For example, cells at
    subdivision level 20 happen to be about 10m wide. If query regions never use
    cells smaller than this, then they will never need to look up more than 20
    levels of ancestors (or if the minimum level is 10, then only 10 levels of
    ancestors). To use this technique, just call
    `S2RegionCoverer::set_max_level()`.

*   Skip some levels of the `S2Cell` hierarchy. For example, if only even levels
    are used, then the hierarchy effectively has a fan-out of 16 rather than 4,
    and the number of ancestors of each cell is halved. This can be implemented
    by modifying the output of `GetCellUnion` similar to the example above.

*   You can reduce the number of query terms by creating both "ancestor" and
    "covering" terms at index time for the cells in each covering. Then queries
    only need to include "ancestor" terms for the cells in the query region,
    rather than both "ancestor" and "covering" terms.

Note that all of these techniques can also be used with indexes that support
range queries (maps, `SSTable`, `BigTable`). In that case, the minimum
subdivision level only needs to be enforced when indexing regions, and the
maximum subdivision level only needs to be used when querying regions.

## Unique Indexes

Suppose that furthermore, we only want to have one index entry for each region
in our index. For example, we might have a collection of features stored in a
`BigTable` and want to look them up and edit them in place.

The problem in this case is that features cannot in general be covered with one
cell. (Consider a small circle at the intersection of three face cells, or a
large region like the southern hemisphere.) The solution is to have a two-level
index with keys of the form `(radius_bucket, cell_id)`. When indexing a feature,
we first compute a bounding spherical cap, which effectively gives us a center
and a radius. The center point is converted to a leaf cell id, and the radius is
discretized into one of a small number of "buckets". For example:

```c++
// Make the smallest bucket have a radius of 100m (converted to radians),
// and use a ratio of 4 between buckets.  Buckets are numbered from 0
// to kMaxBucket (which equals 10 for this choice of parameters).

double kMinRadius = 100 / S2Earth::RadiusMeters();
double kLogRatio = log(4.0);

pair<int, S2CellId> GetKey(const S2Region& region) {
  S2Cap cap = region.GetCapBound();
  // lrint(x + 0.5) is similar to int(x + 1) but much faster.
  int bucket = lrint(0.5 + log(cap.angle().radians() / kMinRadius) / kLogRatio);
  S2CellId cell_id = S2CellId::FromPoint(cap.axis());
  return make_pair(bucket, cell_id);
}
```

To index a set of regions, we sort them by their `(bucket, cell_id)` key (e.g.
by inserting them in a `BigTable`). Then given a query, we can find all the
regions that intersect the query region by doing one lookup per radius bucket.
To find the matching regions in a given bucket, we first expand the query region
on all sides by a distance equal to the bucket radius (i.e. the maximum radius
of regions that were inserted into that bucket). We then do a simple point
location query to find all of the cell ids at that level that lie within the
expanded query region. This is equivalent to finding all the bounding caps at
that level that intersect the original, non-expanded query region. Here is some
example code:

```c++
typedef btree_map<pair<int, S2CellId>, MyRegion*> MyIndex;
static MyIndex index;
static const int kNumBuckets = GetKey(S2Cap::Full()).first + 1;

void FindMatches(const S2LatLngRect& query, vector<MyRegion *> *result) {
  S2RegionCoverer coverer;
  for (int bucket = 0; bucket < kNumBuckets; ++bucket) {
    double radius = min(2 * M_PI, kMinRadius * exp(kLogRatio * bucket));
    S2CellUnion covering = coverer.GetCellUnion(
        query.ExpandedByDistance(S1Angle::Radians(radius)));
    for (S2CellId id : covering) {
      MyIndex::const_iterator
        lo = index.lower_bound(make_pair(bucket, id.range_min())),
        hi = index.lower_bound(make_pair(bucket, id.range_max()));
      for (; lo != hi; ++lo) result->insert(lo->second);
    }
  }
}
```

A few notes about this type of indexing:

*   Since all regions are converted to discs before indexing, long skinny
    regions are not handled efficiently. It is not a problem if there are
    relatively few such regions, but this technique is not a good choice if your
    data is dominated by long linear features.

*   If your regions need to be indexed in more than one way, e.g. if you have a
    unique identifier so that regions can refer to other regions, then having
    only one index entry per region is not much of an advantage. For example, a
    Maps Wiki might have one table in a `BigTable` that stores each feature that
    has been edited, keyed by its unique identifier, and a secondary table that
    is used only for spatial indexing. There is no particular reason that the
    secondary table needs to have only one entry per region.
