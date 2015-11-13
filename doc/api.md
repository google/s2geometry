S2 Geometry Library Overview
=========================

This is a package for manipulating geometric shapes.  It is rather
heavily weighted towards spherical geometry because currently it is
mainly used for geographic data.  Other types of geometric primitives
will be added as necessary.

Currently the package consists of:

* Basic representations of angles, latitude-longitude points, unit 3-vectors,
  and conversions among them.

* Various shapes over the unit sphere, such as spherical caps,
  latitude-longitude rectangles, and polygons.  These
  are collectively known as "regions".

* A hierarchical decomposition of the sphere into regions called "cells".
  The hierarchy starts with the six faces of a projected cube and
  recursively subdivides them in a quadtree-like fashion.

* The ability to approximate arbitrary regions as a collection of
  cells.  This is useful for building inverted indexes that allow
  queries over arbitrarily shaped regions.

The implementations attempt to be precise both in terms of mathematical
definitions (e.g. whether regions include their boundaries,
representations of empty and full regions) and numerical accuracy
(e.g. avoiding cancellation error).

Note that the intent of this package is to represent geometry as a
mathematical abstraction.  For example, although the unit sphere is
obviously a useful approximation for the Earth's surface, functions that
are specifically related to geography should be put elsewhere
(e.g. easting/northing conversions, the Earth's diameter, etc).

Basic Types
===========

Angles
------

The `S1Angle` class represents a one-dimensional angle (as opposed to a
2D solid angle).  It has methods for converting angles to or from
radians, degrees, and the E6/E7 representations (i.e. degrees multiplied
by 1e6/1e7 and rounded to the nearest integer).

    class S1Angle {
     public:
      inline static S1Angle Radians(double radians);
      inline static S1Angle Degrees(double degrees);
      inline static S1Angle E6(long e6);
      inline static S1Angle E7(long e7);

      double radians() const;
      double degrees() const;
      long e6() const;
      long e7() const;
    };

Points
------

The `S2Point` class represents a point on the unit sphere as a 3D
vector.  Usually points are normalized to be unit length, but some
methods do not require this.  The `S2Point` class is simply a synonym
for the `Vector_3d` class from `util/math/vector3-inl.h`, which defines
overloaded operators that make it convenient to write arithmetic
expressions (e.g. `x*p1 + (1-x)*p2`).

Some of its more useful methods include:

    class S2Point /*Vector3_d*/ {
     public:
      S2Point(double x, double y, double z);
      double x(), y(), z();                 // Named component accessors
      double& operator[](int i);            // Return component i (0, 1, 2)
      bool operator==(S2Point const& v);    // Equality testing
      bool operator!=(S2Point const& v);    // Inequality testing
      S2Point operator+=(S2Point const& v); // Add another vector
      S2Point operator-=(S2Point const& v); // Subtract another vector
      S2Point operator*=(double k);         // Multiply by a scalar
      S2Point operator/=(double k);         // Divide by a scalar
      S2Point operator+(S2Point const& v);  // Add two vectors
      S2Point operator-(S2Point const& v);  // Subtract two vectors
      S2Point operator*(double k);          // Multiply by a scalar
      S2Point operator/(double k);          // Divide by a scalar
      double DotProd(S2Point const& v);     // Dot product
      S2Point CrossProd(S2Point const& v);  // Cross product
      double Norm2();                       // Squared L2 norm
      double Norm();                        // L2 norm
      S2Point Normalize();                  // Return a normalized **copy**
      double Angle(S2Point const&v);        // Angle between two vectors (radians)
    };

See `util/math/vector3-inl.h` for details and additional methods.

Latitudes and Longitudes
------------------------

The `S2LatLng` class represents a point on the unit sphere as a pair of
latitude-longitude coordinates.

    class S2LatLng : public Vector2_d {
     public:
      inline S2LatLng(S1Angle const& lat, S1Angle const& lng);
      // Basic constructor.  The latitude and longitude must be within
      // the ranges allowed by is_valid() below.

      inline static S2LatLng FromRadians(double lat_radians, double lng_radians);
      inline static S2LatLng FromDegrees(double lat_degrees, double lng_degrees);
      inline static S2LatLng FromE6(long lat_e6, long lng_e6);
      inline static S2LatLng FromE7(long lat_e7, long lng_e7);
      // Convenience functions -- shorter than calling S1Angle::Radians(), etc.

      explicit S2LatLng(S2Point const& p);
      // Convert a direction vector (not necessarily unit length) to an S2LatLng.

      S1Angle lat() const;
      S1Angle lng() const;
      // Accessor methods.

      inline bool is_valid() const;
      // Return true if the latitude is between -90 and 90 degrees inclusive
      // and the longitude is between -180 and 180 degrees inclusive.

      S2Point ToPoint() const;
      // Convert an S2LatLng to the equivalent unit-length vector (S2Point).
    };

Since `S2LatLng` inherits from `Vector2_d`, there are also operators for
arithmetic expressions, etc, as discussed above for `S2Point`.  Here is
an example code snippet:

      S2LatLng ll = S2LatLng::FromE6(lat32, lng32);
      double d = ll.lat().degrees();
      CHECK_EQ(S2LatLng::FromDegrees(90, 45), S2Point(0, 0, 1));

Regions
=======

An `S2Region` represents a two-dimensional region over the unit sphere.
It is an abstract interface with various concrete subtypes.

The main purpose of this interface is to allow complex regions to be
approximated as simpler regions.  So rather than having a wide variety
of virtual methods that are implemented by all subtypes, the interface
is restricted to methods that are useful for computing approximations.
Here they are:

    class S2Region {
     public:
      virtual S2Cap GetCapBound() const;
      // Return a bounding spherical cap.

      virtual S2LatLngRect GetRectBound() const;
      // Return a bounding latitude-longitude rectangle.

      virtual bool Contains(S2Cell const& cell) const;
      // If this method returns true, the region completely contains the given
      // cell.  Otherwise, either the region does not contain the cell or the
      // containment relationship could not be determined.

      virtual bool MayIntersect(S2Cell const& cell) const = 0;
      // If this method returns false, the region does not intersect the given
      // cell.  Otherwise, either region intersects the cell, or the intersection
      // relationship could not be determined.
    };

In addition, all `S2Region` subtypes implement a
`Contains(S2Point const& p)` method that returns true if the given point
`p` is contained by the region.  The point is generally required to be unit
length, although some subtypes may relax this restriction.

Latitude-Longitude Rectangles
-----------------------------

An `S2LatLngRect` is a type of `S2Region` that represents a rectangle in
latitude-longitude space.  It is capable of representing the empty and
full rectangles as well as single points.  It has an `AddPoint` method
that makes it easy to construct a bounding rectangle for a set of points,
including point sets that span the 180 degree meridian.  Here are its
methods:

    class S2LatLngRect : public S2Region {
     public:
      inline S2LatLngRect(S2LatLng const& lo, S2LatLng const& hi);
      // Construct a rectangle from minimum and maximum latitudes and longitudes.
      // If lo.lng() > hi.lng(), the rectangle spans the 180 degree longitude line.

      inline S2LatLngRect(R1Interval lat, S1Interval const& lng);
      // Construct a rectangle from latitude and longitude intervals.

      S2LatLng lo() const;
      S2LatLng hi() const;
      R1Interval const& lat() const;
      S1Interval const& lng() const;
      S1Angle lat_lo() const;
      S1Angle lat_hi() const;
      S1Angle lng_lo() const;
      S1Angle lng_hi() const;
      // Accessor methods.

      static S2LatLngRect Empty();
      static S2LatLngRect Full();
      // The canonical empty and full rectangles.

      inline bool is_valid() const;
      // Return true if the rectangle is valid, which essentially just means
      // that the latitude bounds do not exceed Pi/2 in absolute value and
      // the longitude bounds do not exceed Pi in absolute value.

      inline bool is_empty() const;
      // Return true if the rectangle is empty, i.e. it contains no points at all.

      inline bool is_full() const;
      // Return true if the rectangle is full, i.e. it contains all points.

      bool is_inverted() const { return lng_.is_inverted(); }
      // Return true if lng_.lo() > lng_.hi(), i.e. the rectangle crosses
      // the 180 degree latitude line.

      S2LatLng GetVertex(int k) const;
      // Return the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order.

      S2LatLng GetCenter() const;
      // Return the center of the rectangle in latitude-longitude space
      // (in general this is not the center of the region on the sphere).

      S2LatLng GetSize() const;
      // Return the width and height of this rectangle in latitude-longitude space.

      bool Contains(S2LatLng const& ll) const;
      // More efficient version of Contains() that accepts a S2LatLng rather than
      // an S2Point.

      bool InteriorContains(S2Point const& p) const;
      // Return true if and only if the given point is contained in the interior
      // of the region (i.e. the region excluding its boundary).

      bool InteriorContains(S2LatLng const& ll) const;
      // More efficient version of InteriorContains() that accepts a S2LatLng
      // rather than an S2Point.

      bool Contains(S2LatLngRect const& other) const;
      // Return true if and only if the rectangle contains the given other
      // rectangle.

      bool InteriorContains(S2LatLngRect const& other) const;
      // Return true if and only if the interior of this rectangle contains all
      // points of the given other rectangle (including its boundary).

      bool Intersects(S2LatLngRect const& other) const;
      // Return true if this rectangle and the given other rectangle have any
      // points in common.

      bool InteriorIntersects(S2LatLngRect const& other) const;
      // Return true if and only if the interior of this rectangle intersects
      // any point (including the boundary) of the given other rectangle.

      void AddPoint(S2Point const& p);
      void AddPoint(S2LatLng const& ll);
      // Increase the size of the bounding rectangle to include the given point.
      // The rectangle is expanded by the minimum amount possible.

      S2LatLngRect Union(S2LatLngRect const& other) const;
      // Return the smallest rectangle containing the union of this rectangle and
      // the given rectangle.

      S2LatLngRect Intersection(S2LatLngRect const& other) const;
      // Return the smallest rectangle containing the intersection of this
      // rectangle and the given rectangle.  Note that the region of intersection
      // may consist of two disjoint rectangles, in which case a single rectangle
      // spanning both of them is returned.

      S2LatLngRect ExpandedByDistance(S1Angle const& angle) const;
      // Return a rectangle that contains the convolution of this rectangle with
      // a cap of the given angle.  This expands the rectangle by a fixed
      // distance (as opposed to growing the rectangle in latitude-longitude
      // space).  The returned rectangle includes all points whose minimum
      // distance to the original rectangle is at most the given angle.

      inline bool operator==(S2LatLngRect const& other) const;
      // Return true if two rectangles contains the same set of points.

      ////////////////////////////////////////////////////////////////////////
      // S2Region interface (see s2region.h for details):

      virtual S2Cap GetCapBound() const;
      virtual S2LatLngRect GetRectBound() const;
      virtual bool Contains(S2Cell const& cell) const;
      virtual bool MayIntersect(S2Cell const& cell) const;
      bool Contains(S2Point const& p) const;
      // The point 'p' is not required to be unit length.
    };

Spherical Caps
--------------

An `S2Cap` represents a spherical cap, i.e. a portion of a sphere cut off
by a plane.  The cap is defined by its axis and height.  This
representation has good numerical accuracy for very small caps, unlike the
(axis, min-distance-from-origin) representation, and is also efficient for
containment tests, unlike the (axis, angle) representation.

Here are the methods available:

    class S2Cap : public S2Region {
     public:
      // Caps may be constructed from either an axis and a height, or an axis and
      // an angle.  To avoid ambiguity, there are no public constructors except
      // the default constructor.

      inline static S2Cap FromAxisHeight(S2Point const& axis, double height);
      // Create a cap given its axis and the cap height, i.e. the maximum
      // projected distance along the cap axis from the cap center.
      // 'axis' should be a unit-length vector.

      inline static S2Cap FromAxisAngle(S2Point const& axis,
                                        S1Angle const& angle);
      // Create a cap given its axis and the cap opening angle, i.e. maximum
      // angle between the axis and a point on the cap.  'axis' should be a
      // unit-length vector, and 'angle' should be between 0 and 180 degrees.

      static S2Cap Empty();
      // Return an empty cap, i.e. a cap that contains no points.

      static S2Cap Full();
      // Return a full cap, i.e. a cap that contains all points.

      S2Point const axis() const;
      double height() const;
      // Accessor methods.

      S1Angle angle() const;
      // Return the cap opening angle in radians, or a negative number for
      // empty caps.

      bool is_valid() const;
      // We allow negative heights (to represent empty caps) but not heights
      // greater than 2.

      bool is_empty() const;
      // Return true if the cap is empty, i.e. it contains no points.

      bool is_full() const;
      // Return true if the cap is full, i.e. it contains all points.

      S2Cap Complement() const;
      // Return the complement of the interior of the cap.  A cap and its
      // complement have the same boundary but do not share any interior points.
      // The complement operator is not a bijection, since the complement of a
      // singleton cap (containing a single point) is the same as the complement
      // of an empty cap.

      bool Contains(S2Cap const& other) const;
      // Return true if and only if this cap contains the given other cap
      // (in a set containment sense, e.g. every cap contains the empty cap).

      bool InteriorIntersects(S2Cap const& other) const;
      // Return true if and only if the interiors of the two caps intersect.
      // (Therefore returns false if either cap is empty or a singleton.)

      bool InteriorContains(S2Point const& p) const;
      // Return true if and only if the given point is contained in the interior
      // of the region (i.e. the region excluding its boundary).

      void AddPoint(S2Point const& p);
      // Increase the cap height if necessary to include the given point.
      // If the cap is empty the axis is set to the given point, but otherwise
      // it is left unchanged.

      ////////////////////////////////////////////////////////////////////////
      // S2Region interface (see s2region.h for details):

      virtual S2Cap GetCapBound() const;
      virtual S2LatLngRect GetRectBound() const;
      virtual bool Contains(S2Cell const& cell) const;
      virtual bool MayIntersect(S2Cell const& cell) const;
      bool Contains(S2Point const& p) const;
    };

S2Cell Hierarchy
================

Note: Octavian Procopiuc (tavi@google.com) created a nice
<a href="https://docs.google.com/presentation/d/1Hl4KapfAENAOf4gv-pSngKwvS_jwNVHRPZTTDzXXn6Q/edit?usp=sharing">
slide deck</a> that discusses the S2 cell hierarchy and coordinate systems.

The library provides methods for subdividing the sphere into a hierarchical
collection of "cells".  The cells have a space-filling curve structure
that makes them useful for spatial indexing.  There are methods for
approximating arbitrary regions as a collection of cells.

The `S2Cell` structure is defined as follows.  There are six top-level
"face cells", obtained by projecting the six faces of a cube -- (face, u, v) coordinates below --
onto the unit sphere -- (x, y, z) coordinates below. We call this cube the uv-cube.

Each face is then subdivided recursively into four cells
in a quadtree-like fashion. On the uv-cube, a cell is a rectangle whose edges are aligned with the
sides of its face.  On the sphere, it is a spherical quadrilateral bounded by
four geodesics (great circle segments).

There are a total of 30
levels of subdivision defined (i.e. 6 * 4<sup>30</sup> leaf cells),
which gives a resolution of about 1cm
everywhere on a sphere the size of the earth. Details on the cell areas at
each level appear on the [S2 Statistics page](cell_statistics.md).  Each cell
is uniquely identified by a 64-bit _cell id_.

Coordinate Systems
------------------

In order for cells to be roughly the same size on the sphere, they are not
the same size on the uv-cube.  We introduce another cube -- the st-cube.  On this cube.
the cells are perfectly square and divided through their center.  The st-cube is
projected on the uv-cube so that cells at the periphery of an st-face are larger than
cells at their center.  In a 2-d projection, it looks like this:

<a href="xyz_to_uv_to_st.jpg">
<img src="xyz_to_uv_to_st.jpg" alt="xyz_to_uv_to_st.jpg"  width="320" height="240"  />
</a>

Note that a cell on the st-cube is not bounded by straight lines.

There are a few more coordinate systems worth introducing:

* (id)<br> Cell id.  An `S2CellId` is a 64-bit encoding of a face and a
  Hilbert curve position on that face, as discussed below.  The Hilbert
  curve position implicitly encodes both the position of a cell and its
  subdivision level.

* (face, i, j)<br> Leaf-cell coordinates.  "i" and "j" are integers in
  the range [0,(2**30)-1] that identify a particular leaf cell on a
  given face.  The (i, j) coordinate system is right-handed on every
  face, and the faces are oriented such that Hilbert curves connect

* (face, s, t)<br> Cell-space coordinates.  "s" and "t" are real numbers
  in the range [0,1] that identify a point on the given face.  For
  example, the point (s, t) = (0.5, 0.5) corresponds to the center of the
  top-level face cell.  This point is also a vertex of exactly four
  cells at each subdivision level greater than zero.

* (face, si, ti)<br> Discrete cell-space coordinates.  These are
  obtained by multiplying "s" and "t" by 2**31 and rounding to the
  nearest integer.  Discrete coordinates lie in the range
  [0,2**31].  This coordinate system can represent the edge and
  center positions of all cells with no loss of precision (including
  non-leaf cells).

* (face, u, v)<br> Cube-space coordinates.  To make the cells at each
  level more uniform in size after they are projected onto the sphere,
  we apply a nonlinear transformation of the form u=f(s), v=f(t), where
  f(s) = (4 s^2 - 1) / 3 if s >= 1/2 and (1 - 4 (1 - s)^2) / 3 if s < 1/2.
  The (u, v) coordinates after this transformation give the actual
  coordinates of a point on the cube face (modulo some 90 degree
  rotations) before it is projected onto the unit sphere.

* (x, y, z)<br> Direction vector (`S2Point`).  Direction vectors are not
  necessarily unit length, and are often chosen to be points on the
  biunit cube [-1,+1]x[-1,+1]x[-1,+1].  They can be be normalized
  to obtain the corresponding point on the unit sphere.

* (lat, lng)<br> Latitude and longitude (`S2LatLng`).  Latitudes must be
  between -90 and 90 degrees inclusive, and longitudes must be between
  -180 and 180 degrees inclusive.

Note that the (i, j), (s, t), (si, ti), and (u, v) coordinate systems are
right-handed on all six faces.

Conversion between coordinate systems can lead to precision problems.

Cell Ids
--------

Cell ids are a convenient representation for both points and regions
on the unit sphere.  Points are generally represented as leaf
cells, while regions are represented as collections of cells
at any level.

An `S2CellId` is a 64-bit unsigned integer that uniquely identifies a
cell in the S2 cell decomposition.  It has the following format:

```
   id = [face][face_pos]
```

where "face" is a 3-bit number (range 0..5) encoding the cube face,
and "face_pos" is a 61-bit number encoding the position of the center
of this cell along a space-filling curve over this face (see below).

In particular, the id of a cell at level "k" consists of a 3-bit face
number, followed by k pairs of bits that recursively select one of the
four children of each cell.  The remainder consists of a single 1-bit
followed by zeros (this is the Hilbert curve position corresponding
to the center of the cell as discussed above).  For example:

    01010000...0   The top-level cell of face 2.  (The Hilbert curve
                   position 0.10* corresponds to the center of the cell.)
    00110100...0   Subcell 10 of the top-level cell of face 1.

Cell ids have the following convenient properties:

* The level of a cell id can be determined by looking at the position
  of its lowest-numbered 1-bit.  For a cell at level k, the
  lowest-numbered 1-bit is at position 2*(kMaxLevel-k).

* The id of a parent cell is at the midpoint of the range of ids
  spanned by its children (or by its descendants at any level).

S2CellId to st Coordinates
--------------------------

Each cell is a spherical quadrilateral bounded by four geodesics (great
circle segments).  The level 0 subdivision (the top level) consists of six
face cells.  The children of each face cell are numbered according
to their traversal order along a Hilbert curve.

The canonical Hilbert curve on a cell visits its children in the
following (i,j) order: (0,0), (0,1), (1,1), (1,0).  Graphically, the
traversal order looks like an inverted U:

             (0,1) ^----> (1,1)
                   |    |
                   |    |
             (0,0) x    v (1,0)

The traversal order in each subcell is obtained by applying the
following transformations to the ordering in its parent cell:

          (0,0):   (i,j) -> (j,i)       [axes swapped]
          (0,1):   (i,j) -> (i,j)
          (1,1):   (i,j) -> (i,j)
          (1,0):   (i,j) -> (1-j,1-i)   [axes swapped and inverted]

For example, the traversal order in the four subcells of cell (0,0) is
(0,0), (1,0), (1,1), (0,1) (the i- and j-axes have been swapped).  After
one level of expansion according to these rules, the curve looks like this:

              ^-->  ^-->
              |  |  |  |
              ^  v-->  v (*)
              |        |
              <--^  <--v
                 |  |
              x-->  v-->

Each leaf cell can be assigned a label according to its position in the
traversal order.  For example, the traversal above has 16 positions
ranging from 0000 to 1111 in binary notation.  The cell marked (*) has
the label 1011.  Notice that we can obtain the traversal position of the
parent of any cell (within its own traversal order) by stripping the
last two bits from its label.

We typically think of the Hilbert curve position as a real number
between 0 and 1, by prefixing "0." to the binary number above.
For example, position 0.1000 in binary (0.5 in decimal) is at
the beginning of curve in subcell 10 of the root cell, which
turns out to be the middle of the top-level face.  Note that
it is true in general that the Hilbert curve position at the
middle of any cell is the average of the positions at the points
where the curve enters and exits that cell.

A particular cell can be uniquely identified by the Hilbert
curve position at the center of that cell.  Note that no two
cells have the same center position.  For example, the Hilbert
curve position at the center of the top-level face cell is
0.10* in binary notation (where 0* denotes infinitely many
zeroes).  The position at the center of subcell (0,0) of the
top-level face cell is 0.0010*.  Note that a Hilbert curve
position specifies both the position **and** the subdivision
level of a particular cell.

[Images of the Earth projected onto S2 cells can be found
here.](earthcube/earthcube.md)

`S2CellId` class
----------------

The `S2CellId` class is a thin wrapper over a 64-bit cell id that
provides methods for navigating the cell hierarchy (finding parents,
children, containment tests, etc).  Since leaf cells are often used
to represent points on the unit sphere, the `S2CellId` class also
provides methods for converting directly to and from an `S2Point`.
Here are its methods:

    class S2CellId {
     public:
      static int const kFaceBits = 3;
      static int const kNumFaces = 6;
      static int const kMaxLevel = 30;  // Valid levels: 0..kMaxLevel
      static int const kPosBits = 2 * kMaxLevel + 1;
      // Although only 60 bits are needed to represent the index of a leaf
      // cell, we need an extra bit in order to represent the position of
      // the center of the leaf cell along the Hilbert curve.

      inline explicit S2CellId(uint64 id) : id_(id) {}

      inline S2CellId() : id_(0) {}
      inline static S2CellId None() { return S2CellId(); }
      // The default constructor returns an invalid cell id.

      inline static S2CellId Sentinel() { return S2CellId(~uint64(0)); }
      // Returns an invalid cell id guaranteed to be larger than any
      // valid cell id.  Useful for creating indexes.

      static S2CellId FromFacePosLevel(int face, uint64 pos, int level);
      // Return a cell give its face (range 0..5), 61-bit Hilbert curve position
      // within that face, and level (range 0..kMaxLevel).  The given position
      // will be modified to correspond to the Hilbert curve position at the
      // center of the returned cell.  This is a static function rather than a
      // constructor in order to give names to the arguments.

      static S2CellId FromPoint(S2Point const& p);
      // Return the leaf cell containing the given point (a direction
      // vector, not necessarily unit length).

      static S2CellId FromLatLng(S2LatLng const& ll);
      // Return the leaf cell containing the given S2LatLng.

      S2Point ToPoint() const { return ToPointRaw().Normalize(); }
      S2Point ToPointRaw() const;
      // Return the direction vector corresponding to the center of the given
      // cell.  The vector returned by ToPointRaw is not necessarily unit length.

      S2LatLng ToLatLng() const;
      // Return the S2LatLng corresponding to the center of the given cell.

      inline uint64 id() const { return id_; }
      // The 64-bit unique identifier for this cell.

      inline bool is_valid() const;
      // Return true if id() represents a valid cell.

      inline int face() const;
      // Which cube face this cell belongs to, in the range 0..5.

      inline uint64 pos() const;
      // The position of the cell center along the Hilbert curve over this face,
      // in the range 0..(2**kPosBits-1).

      int level() const;
      // Return the subdivision level of the cell (range 0..kMaxLevel).

      inline bool is_leaf() const;
      // Return true if this is a leaf cell (more efficient than checking
      // whether level() == kMaxLevel).

      inline bool is_face() const;
      // Return true if this is a top-level face cell (more efficient than
      // checking whether level() == 0).

      inline S2CellId range_min() const;
      inline S2CellId range_max() const;
      // Methods that return the range of cell ids that are contained
      // within this cell (including itself).  The range is **inclusive**
      // (i.e. test using >= and <=) and the return values of both
      // methods are valid leaf cell ids.
      //
      // These methods should not be used for iteration.  If you want to
      // iterate through all the leaf cells, call child_begin(kMaxLevel) and
      // child_end(kMaxLevel) instead.
      //
      // It would in fact be error-prone to define a range_end() method,
      // because (range_max().id() + 1) is not always a valid cell id, and the
      // iterator would need to be tested using "<" rather that the usual "!=".

      inline bool contains(S2CellId const& other) const;
      // Return true if the given cell is contained within this one.

      inline S2CellId parent() const;
      inline S2CellId parent(int level) const;
      // Return the cell at the previous level or at the given level.

      inline S2CellId child_begin() const;
      inline S2CellId child_begin(int level) const;
      inline S2CellId child_end() const;
      inline S2CellId child_end(int level) const;
      // Iterator-style methods for traversing the immediate children of a cell
      // or all of the children at a given level.  Note that the end value is
      // exclusive, just like standard STL iterators.  You should iterate
      // using code like this:
      //
      //   for(S2CellId c = id.child_begin(); c != id.child_end(); c = c.next())
      //     ...
      //
      // The convention for advancing the iterator is "c = c.next()" rather
      // than "++c" to avoid possible confusion with incrementing the
      // underlying 64-bit cell id.

      inline S2CellId next() const;
      inline S2CellId prev() const;
      // Return the next/previous cell at the same level along the Hilbert curve.
      // Works correctly when advancing from one face to the next, but
      // does **not** wrap around from the last face to the first or vice versa.

      static inline S2CellId GetBegin(int level);
      static inline S2CellId GetEnd(int level);
      // Iterator-style methods for traversing all the cells along the Hilbert
      // curve at a given level (across all 6 faces of the cube).  Note that
      // the end value is exclusive (just like standard STL iterators).

      /////////////////////////////////////////////////////////////////////
      // "Internal" Methods
      //
      // These methods are used within the S2 library and are also exposed
      // for testing purposes.

      int ToFaceIJOrientation(int* pi, int* pj, int* orientation) const;
      // Return the (face, i, j) coordinates for the leaf cell corresponding to
      // this cell id.  Since cells are represented by the Hilbert curve position
      // at the center of the cell, the returned (i,j) for non-leaf cells will be
      // a leaf cell adjacent to the cell center.  If "orientation" is non-NULL,
      // also return the Hilbert curve orientation for the current cell.

      int64 num_leaf_cells() const;
      // Return the number of leaf cells contained in this cell.

      static S2CellId FromFaceIJ(int face, int i, int j);
      // Return a leaf cell given its cube face (range 0..5) and
      // i- and j-coordinates (see s2.h).
    };

Cell Regions
------------

An `S2Cell` is an `S2Region` object that represents a cell.  Unlike `S2CellId`,
it views a cell as a representing a spherical quadrilateral rather than a point,
and it supports efficient containment and intersection tests.  However, it is also
a more expensive representation (currently 48 bytes rather than 8).

Here are its methods:

    class S2Cell : public S2Region {
     public:
      S2Cell() {}
      // The default constructor is required in order to use freelists.
      // It should not be used otherwise.

      explicit S2Cell(S2CellId const& id);
      // An S2Cell always corresponds to a particular S2CellId.  The other
      // constructors are just convenience methods.

      static S2Cell FromFacePosLevel(int face, uint64 pos, int level);
      explicit S2Cell(S2Point const& p);
      explicit S2Cell(S2LatLng const& ll);
      // Convenience methods.

      static S2Cell* New(S2Cell const& cell);
      static void Delete(S2Cell* cell);
      // We provide a global, thread-safe freelist for clients that need to
      // allocate and delete lots of cells efficiently.

      inline S2CellId id() const;
      inline int face() const;
      inline int level() const;
      inline int orientation() const;
      inline bool is_leaf() const;

      S2Point GetVertex(int k) const { return GetVertexRaw(k).Normalize(); }
      S2Point GetVertexRaw(int k) const;
      // Return the k-th vertex of the cell (k = 0,1,2,3).  Vertices are returned
      // in CCW order.  The points returned by GetVertexRaw are not necessarily
      // unit length.

      S2Point GetEdge(int k) const { return GetEdgeRaw(k).Normalize(); }
      S2Point GetEdgeRaw(int k) const;
      // Return the inward-facing normal of the great circle passing through
      // the edge from vertex k to vertex k+1 (mod 4).  The normals returned
      // by GetEdgeRaw are not necessarily unit length.

      bool Subdivide(S2Cell children[4]) const;
      // If this is not a leaf cell, set children[0..3] to the four children of
      // this cell (in traversal order) and return true.  Otherwise returns false.
      // This method is equivalent to the following:
      //
      // for (pos=0, id=child_begin(); id != child_end(); id = id.next(), ++pos)
      //   children[i] = S2Cell(id);
      //
      // except that it is about 2.5 times faster.

      S2Point GetCenter() const { return GetCenterRaw().Normalize(); }
      S2Point GetCenterRaw() const;
      // Return the direction vector corresponding to the center in (s,t)-space of
      // the given cell.  This is the point at which the cell is divided into four
      // subcells; it is not necessarily the centroid of the cell in (u,v)-space
      // or (x,y,z)-space.  The point returned by GetCenterRaw is not necessarily
      // unit length.

      static double AverageArea(int level);
      // Return the average area for cells at the given level.

      double AverageArea() const { return AverageArea(level_); }
      // Return the average area of cells at this level.  This is accurate to
      // within a factor of 1.7 (for S2_QUADRATIC_PROJECTION) and is extremely
      // cheap to compute.

      double ApproxArea() const;
      // Return the approximate area of this cell.  This method is accurate to
      // within 3% percent for all cell sizes and accurate to within 0.1% for
      // cells at level 5 or higher (i.e. 300km square or smaller).  It is
      // moderately cheap to compute.

      double ExactArea() const;
      // Return the area of this cell as accurately as possible.  This method is
      // more expensive but it is accurate to 6 digits of precision even for leaf
      // cells (whose area is approximately 1e-18).

      ////////////////////////////////////////////////////////////////////////
      // S2Region interface (see s2region.h for details):

      virtual S2Cap GetCapBound() const;
      virtual S2LatLngRect GetRectBound() const;
      virtual bool Contains(S2Cell const& cell) const;
      virtual bool MayIntersect(S2Cell const& cell) const;
      bool Contains(S2Point const& p) const;
      // The point 'p' is not required to be unit length.
    };

Cell Unions
-----------

An `S2CellUnion` is a region consisting of cells of various sizes.  Typically
a cell union is used to approximate some other shape.  There is a tradeoff
between the accuracy of the approximation and how many cells are used.
Unlike polygons, cells have a fixed hierarchical structure.  This makes
them more suitable for spatial indexing.

Here is the interface:

    class S2CellUnion : public S2Region {
     public:
      S2CellUnion();

      void Init(vector<S2Cell> const& cells);
      void Init(vector<S2CellId> const& ids);
      void Init(vector<uint64> const& ids);
      // Populates a cell union with the given cells, S2CellIds, or 64-bit cell
      // ids.  May be called repeatedly to change the contents of the cell union.

      int num_cells() const;
      S2CellId const& cell_id(int i) const;
      vector<S2CellId> const& cell_ids() const;
      // There are convenience methods as well as direct access to the
      // underlying vector.

      ////////////////////////////////////////////////////////////////////////
      // S2Region interface (see s2region.h for details):

      virtual S2Cap GetCapBound() const;
      virtual S2LatLngRect GetRectBound() const;
      virtual bool Contains(S2Cell const& cell) const;
      virtual bool MayIntersect(S2Cell const& cell) const;
      bool Contains(S2Point const& p) const;
    };

Approximating Regions
---------------------

An `S2RegionCoverer` is a class that allows arbitrary regions to be
approximated as unions of cells (`S2CellUnion`).  This is useful for
implementing various sorts of search and precomputation operations.

Typical usage:

    S2RegionCoverer coverer;
    coverer.set_max_cells(5);
    S2Cap cap = S2Cap::FromAxisAngle(...);
    S2CellUnion covering;
    coverer.GetCellUnion(cap, &covering);

This yields a cell union of at most 5 cells that is guaranteed to cover the
given cap (a disc-shaped region on the sphere).

The approximation algorithm is not optimal but does a pretty good job in
practice.  The output does not always use the maximum number of cells
allowed, both because this would not always yield a better approximation,
and because `max_cells()` is a limit on how much work is done exploring the
possible covering as well as a limit on the final output size.

The default `max_cells()` value of 8 is a reasonable point on the
precision/performance tradeoff curve for caps, but you may want to use
higher values for odd-shaped regions, such as long skinny lat-lng rects.
[Here are some examples of approximations.](coverer/coverer.md)

Here is the interface to `S2RegionCoverer`:

    class S2RegionCoverer {
     public:
      static int const kDefaultMaxCells = 8;
      static int const kDefaultMaxLevel = S2CellId::kMaxLevel;
      // By default, the covering uses at most 8 cells at any level.  This gives
      // a reasonable tradeoff between the number of cells used and the accuracy
      // of the approximation.
      //
      // Here are the median and worst case values for the area ratio when
      // approximating 100,000 discs of random size:
      //
      //   max_cells:        3      4     5     6     8    12    20   100   1000
      //   median ratio:  5.31   3.31  2.72  2.33  1.98  1.66  1.42  1.11   1.01
      //   worst case:     7e5  15.83  9.84  5.10  4.03  2.74  1.94  1.19   1.03
      //
      // Note that if max_cells() is set to less than 4, the area of the covering
      // may be arbitrarily large compared to the area of the original region,
      // even if the original region is convex (e.g. an S2Cap or S2LatLngRect).
      // With max_cells() == 4, the maximum ratio of the area of the covering
      // to the area of an S2Cap is about 16 (this is a pathological worst case).

      S2RegionCoverer();
      ~S2RegionCoverer();
      int max_cells() const;
      int max_level() const;
      void set_max_cells(int max_cells);
      void set_max_level(int max_level);
      // Sets the maximum desired number of cells in the approximation and the
      // maximum subdivision level of the cells to be used.  max_cells should
      // not be set to less than 4 if you want a reasonable approximation
      // (see note above).
      //
      // Note that for any setting of max_cells, up to 6 cells may be returned if
      // that is the minimum number of cells required (i.e. if the region
      // intersects all six face cells).  Up to 3 cells may be returned even for
      // very tiny convex regions if they happen to be located at the intersection
      // of three cube faces.

      void GetCovering(S2Region const& region, vector<S2CellId>* covering);
      void GetInteriorCovering(S2Region const& region, vector<S2CellId>* interior);
      // Return a vector of cell ids that covers (GetCovering) or is contained
      // within (GetInteriorCovering) the given region and satisfies the various
      // restrictions specified above.

      void GetCellUnion(S2Region const& region, S2CellUnion* covering);
      void GetInteriorCellUnion(S2Region const& region, S2CellUnion* interior);
      // Return a normalized cell union that covers (GetCellUnion) or is contained
      // within (GetInteriorCellUnion) the given region and satisfies the
      // restrictions **EXCEPT** for min_level() and level_mod().  These criteria
      // cannot be satisfied using a cell union because cell unions are
      // automatically normalized by replacing four child cells with their parent
      // whenever possible.  (Note that the list of cell ids passed to the cell
      // union constructor does in fact satisfy all the given restrictions.)

    };

Spatial Indexing
----------------

There are several reasonable approaches to using cells for spatial
indexing, depending on what you want to index (points vs. regions),
what operations your indexing data structure support (range queries
vs. key lookups only), and what sort of queries you want to make
(point location vs. region queries).  Let's go over a few examples.

### Indexing Points

First, suppose you have a bunch of points in memory, and you want
to repeatedly find the set of points enclosed by various disc-shaped
regions (e.g. all points within 5km of the north pole).  Do the
following:

* Convert each point to a corresponding leaf cell id
  (`S2CellId::FromPoint`).
* Put all the cell ids in a vector and sort them.
* Construct an `S2Cap` corresponding to the query region
  (`S2Cap::FromAxisAngle`).
* Find a set of cell ids that cover the cap (an `S2CellUnion`)
   by calling `S2RegionCoverer::GetCellUnion`.
* For each cell "x" in the covering, do a range query on the vector
  of points to find all the leaf cell ids that are contained by "x".
* If you only want the points that are exactly contained by the disc,
  rather than just a superset of them, test each point against the cap.

Here is some example code:

    vector<S2CellId> index;
    for (int i = 0; i < points.size(); ++i) {
      index.push_back(S2CellId::FromPoint(points[i]));
    }
    sort(index.begin(), index.end());
    S2Cap cap = S2Cap::FromAxisAngle(center, radius);
    S2RegionCoverer coverer;
    coverer.set_max_cells(12);  // Default is 8.
    S2CellUnion covering;
    coverer.GetCellUnion(cap, &covering);
    vector<S2CellId> const& cell_ids = covering.cell_ids();
    for (int i = 0; i < cell_ids.size(); ++i) {
      int lo = lower_bound(index.begin(), index.end(), cell_ids[i].range_min());
      int hi = upper_bound(index.begin(), index.end(), cell_ids[i].range_max());
      for (int j = lo; j < hi; ++j) {
        if (cap.Contains(index[j].ToPoint())) {
          LOG(INFO) << "Cap contains point " << index[j].ToPoint();
        }
      }
    }

Some observations:

* If you want to add or delete points from the index, use a
  `set<S2CellId>` rather than a vector.  Sets have built-in
  `lower_bound()` and `upper_bound()` methods.  If you have
  additional data that you want to associate with each cell id,
  use a `map<S2CellId, MyData*>` or `multimap<S2CellId, MyData*>`.

* If the points don't all fit in memory, suitable choices for the
  index include `SSTable` and `BigTable`, both of which support
  efficient range queries.

* The code above can be made slightly more efficient by skipping
  the upper_bound() call and instead advancing sequentially until
  the id exceeds range_max().  To use this approach, you should
  insert an S2CellId::Sentinel() in the index so that you don't
  also need to check for the end of the vector.  For example:

```
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

### Indexing Regions

Suppose that the objects to be indexed are regions rather than points.
The easiest approach is convert each region to a collection of cell ids,
and create one index entry for each cell id.  This increases the index
size by a small constant factor that depends on how accurate an
approximation is desired.

To execute a query, convert the query region to a collection of cell ids
and do a range query for each one just as before.  However, since index
entries now consist of cells of all sizes, an extra step is required.
For each cell id in the query region, you also need to compute all the
ancestors of those cell ids and do key lookups on those as well:

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

This code computes the set of original regions whose coverings intersect
the covering of the query region.  The regions can be filtered further
if a more exact intersection test is desired.

See s2edgeindex.cc for a fully worked-out example with edges
(where a lot of optimizations can be done).

### Indexes Without Range Queries

Indexing systems that do not support efficient range queries need a different
alternative for indexing regions, which is to expand the index by
inserting all the ancestors of each cell.  To make this efficient, the
ancestor cells are labelled differently than the cells that actually
belong to the coverings.  In other words, there are two separate terms
(tokens) for each cell id: one for cells that belong to a covering
("covering" terms), and another for ancestors of these cells
("ancestor" terms).

To execute a query, we simply compute a covering for the query region,
and then look up all the "ancestor" terms for the cells in this covering,
plus all the "covering" terms for these cells and their ancestors.
This effectively divides the query into two parts: one to find all of
the regions whose covering contains a cell that is descendant of a cell
in the query region, and another to find all the regions whose covering
contains a cell that is a ancestor of a cell in the query region.
(Cells that are identical to a cell in the query region are handled as
part of the second group, since we only inserted "covering" terms for
these cells into the index.)

This technique does not require any range queries, and can efficiently
compute overlaps between query and index regions of any size.
The total number of lookups per query is about 4-10 for
the base covering (depending on how accurate an approximation is
desired) plus about 10-20 for the ancestors of these cells (assuming the
query region diameter is at least 10m).

If it is important to reduce the number of lookups, here are a few ways
to do this:

* Set a minimum subdivision level for cell ids.  For example, cells
  at subdivision level 10 happen to be about 10km wide.  If the index
  only contains cells that are 10km or smaller (i.e. level 10 or
  higher), then ancestors at levels less than 10 can be skipped during
  the query process (eliminating 10 lookups).  To use this technique,
  just call `GetCellUnion` normally, and then iterate through the
  children of any cells that are too big as outlined below.  Of
  course, very large regions such as Canada will need a lot of cells.

  ```
   S2CellId id = cells[i].id();
   int level = max(id.level(), kMinIndexLevel);
   S2CellId end = id.child_end(level);
   for (S2CellId c = id.child_begin(level); c != end; c = c.next()) {
     index[c] = region;
   }
  ```

* Set a maximum subdivision level for cell ids.  For example, cells
  at subdivision level 20 happen to be about 10m wide.  If query
  regions never use cells smaller than this, then they will never need
  to look up more than 20 levels of ancestors (or if the minimum level
  is 10, then only 10 levels of ancestors).  To use this technique,
  just call `S2RegionCoverer::set_max_level()`.

* Skip some levels of the `S2Cell` hierarchy.  For example, if only
  even levels are used, then the hierarchy effectively has a fan-out
  of 16 rather than 4, and the number of ancestors of each cell is
  halved.  This can be implemented by modifying the output of
  `GetCellUnion` similar to the example above.

* You can reduce the number of query terms by creating both "ancestor"
  and "covering" terms at index time for the cells in each covering.
  Then queries only need to include "ancestor" terms for the cells in
  the query region, rather than both "ancestor" and "covering" terms.

Note that all of these techniques can also be used with
indexes that support range queries (maps, `SSTable`, `BigTable`).
In that case, the minimum subdivision level only needs to be enforced
when indexing regions, and the maximum subdivision level only needs to
be used when querying regions.

### Unique Indexes

Suppose that furthermore, we only want to have one index entry for each
region in our index.  For example, we might have a collection of
features stored in a `BigTable` and want to look them up and edit
them in place.

The problem in this case is that features cannot in general be covered
with one cell.  (Consider a small circle at the intersection of three
face cells, or a large region like the southern hemisphere.)  The
solution is to have a two-level index with keys of the form
`(radius_bucket, cell_id)`.  When indexing a feature, we first compute a
bounding spherical cap, which effectively gives us a center and a
radius.  The center point is converted to a leaf cell id, and the radius
is discretized into one of a small number of "buckets".  For example:

```
    // Make the smallest bucket have a radius of 100m (converted to radians),
    // and use a ratio of 4 between buckets.  Buckets are numbered from 0
    // to kMaxBucket (which equals 10 for this choice of parameters).

    double kMinRadius = 100 / S2Earth::RadiusMeters();
    double kLogRatio = log(4.0);

    pair<int, S2CellId> GetKey(S2Region const& region) {
      S2Cap cap = region.GetCapBound();
      // lrint(x + 0.5) is similar to int(x + 1) but much faster.
      int bucket = lrint(0.5 + log(cap.angle().radians() / kMinRadius) / kLogRatio);
      S2CellId cell_id = S2CellId::FromPoint(cap.axis());
      return make_pair(bucket, cell_id);
    }
```

To index a set of regions, we sort them by their `(bucket, cell_id)` key
(e.g. by inserting them in a `BigTable`).  Then given a query, we can
find all the regions that intersect the query region by doing one lookup
per radius bucket.  To find the matching regions in a given bucket, we
first expand the query region on all sides by a distance equal to the
bucket radius (i.e. the maximum radius of regions that were inserted
into that bucket).  We then do a simple point location query to find all
of the cell ids at that level that lie within the expanded query
region.  This is equivalent to finding all the bounding caps at that
level that intersect the original, non-expanded query region.  Here is
some example code:

    typedef map<pair<int, S2CellId>, MyRegion*> MyIndex;
    static MyIndex index;
    static int const kNumBuckets = GetKey(S2Cap::Full()).first + 1;

    void FindMatches(S2LatLngRect const& query, vector<MyRegion *> *result) {
      S2RegionCoverer coverer;
      S2CellUnion covering;
      for (int bucket = 0; bucket < kNumBuckets; ++bucket) {
        double radius = min(2 * M_PI, kMinRadius * exp(kLogRatio * bucket));
        coverer.GetCellUnion(query.ExpandedByDistance(S1Angle::Radians(radius)),
                             &covering);
        for (int j = 0; j < covering.num_cells(); ++j) {
          S2CellId id = covering.cell_id(j);
          MyIndex::const_iterator
            lo = index.lower_bound(make_pair(bucket, id.range_min())),
            hi = index.lower_bound(make_pair(bucket, id.range_max()));
          for (; lo != hi; ++lo) result->insert(lo->second);
        }
      }
    }

A few notes about this type of indexing:

* Since all regions are converted to discs before indexing, long
  skinny regions are not handled efficiently.  It is not a problem
  if there are relatively few such regions, but this technique is
  not a good choice if your data is dominated by long linear features.

* If your regions need to be indexed in more than one way, e.g. if
  you have a unique identifier so that regions can refer to other
  regions, then having only one index entry per region is not much
  of an advantage.  For example, a Maps Wiki might have one table
  in a `BigTable` that stores each feature that has been edited,
  keyed by its unique identifier, and a secondary table that is
  used only for spatial indexing.  There is no particular reason
  that the secondary table needs to have only one entry per region.

Appendix: Alternatives Considered
=================================

The S2 subdivision scheme has a number of advantages over the Hierarchical
Triangular Mesh (http://skyserver.org/HTM) framework currently used by
local search and others:

* Converting from a cell id to a point or vice versa is about
  100 times faster than the corresponding HTM operation.  For example,
  a unit vector can be converted to an S2CellId in about
  0.15 microseconds, while converting a unit vector to an HTM triangle
  id takes about 25 microseconds.  Similarly, converting an S2CellId to
  a vector takes about 0.04 microseconds, while the same HTM operation
  takes about 19 microseconds.  (If full resolution is not needed, the
  HTM times can be improved by using fewer levels -- e.g. by using
  15 levels rather than 30, conversions are about twice as fast, but
  the maximum resolution is reduced from 1cm to about 300 meters.
  The S2 library is still about 100 times faster, though.)

* The S2 library has no storage requirements.  In contrast, the HTM
  library precomputes the upper levels of the mesh during initialization
  to get even the performance numbers mentioned above.  Computing the
  first 6 levels (the default) takes about 2MB of memory and 5
  milliseconds.  It is not practical to precompute much more than this
  because memory requirements increase by a factor of 4 for each
  additional level.

* A linear scan of a set of S2CellIds follows a space-filling curve,
  which maximizes locality of reference.  For example, if each cell id
  is looked up in a spatial index of some sort, a space-filling curve
  minimizes the number of times that we switch from one index element to
  the next.  In theory this should help cache performance (including
  non-local caches such as BigTable).  A linear scan through HTM
  triangles also has fairly good locality of reference due to the
  hierarchical structure, but it does not define a space-filling curve
  (due to the triangle ordering chosen by its inventors) and therefore
  the locality of reference is not quite as good.  (At any given level,
  the HTM path is about 2.2 times longer than the corresonding S2CellId
  path.)

Another frequently mentioned alternative is HEALPix
[http://www.eso.org/science/healpix](http://www.eso.org/science/healpix).
The main advantage of HEALPix is that it is suitable for calculations
involving spherical harmonics.  This isn't relevant to any of our
current applications, and the scheme otherwise has several disadvantages
from our point of view (e.g. cell boundaries are not geodesics, base
structure is more complicated).

The COBE quadrilateralized spherical cube projection
(http://lambda.gsfc.nasa.gov/product/cobe/skymap_info_new.cfm)
also starts by projecting the six faces of a cube onto the unit
sphere, and subdivides each face in a quadtree fashion.  However,
it does not use a space-filling curve scheme for labelling the
cells, the cell edges are not geodesics, and it uses a much more
complicated projection scheme designed to minimize distortion.
