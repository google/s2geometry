# S2 Geometry Library

[TOC]

## Overview

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

## Basic Types

### S1Angle

The `S1Angle` class represents a one-dimensional angle (as opposed to a
2D solid angle).  It has methods for converting angles to or from
radians, degrees, and the E5/E6/E7 representations (i.e. degrees multiplied
by 1e5/1e6/1e7 and rounded to the nearest integer).

```c++
class S1Angle {
 public:
  // These methods construct S1Angle objects from their measure in radians
  // or degrees.
  static S1Angle Radians(double radians);
  static S1Angle Degrees(double degrees);
  static S1Angle E5(int32 e5);
  static S1Angle E6(int32 e6);
  static S1Angle E7(int32 e7);

  // Convenience functions -- to use when args have been fixed32s in protos.
  //
  // The arguments are static_cast into int32, so very large unsigned values
  // are treated as negative numbers.
  static S1Angle UnsignedE6(uint32 e6);
  static S1Angle UnsignedE7(uint32 e7);

  // The default constructor yields a zero angle.  This is useful for STL
  // containers and class methods with output arguments.
  S1Angle();

  // Return an angle larger than any finite angle.
  static S1Angle Infinity();

  // A explicit shorthand for the default constructor.
  static S1Angle Zero();

  // Return the angle between two points, which is also equal to the distance
  // between these points on the unit sphere.  The points do not need to be
  // normalized.
  S1Angle(S2Point const& x, S2Point const& y);

  // Like the constructor above, but return the angle (i.e., distance)
  // between two S2LatLng points.
  S1Angle(S2LatLng const& x, S2LatLng const& y);

  double radians() const;
  double degrees() const;

  int32 e5() const;
  int32 e6() const;
  int32 e7() const;

  // Return the absolute value of an angle.
  S1Angle abs() const;

  // Comparison operators.
  friend bool operator==(S1Angle x, S1Angle y);
  friend bool operator!=(S1Angle x, S1Angle y);
  friend bool operator<(S1Angle x, S1Angle y);
  friend bool operator>(S1Angle x, S1Angle y);
  friend bool operator<=(S1Angle x, S1Angle y);
  friend bool operator>=(S1Angle x, S1Angle y);

  // Simple arithmetic operators for manipulating S1Angles.
  friend S1Angle operator-(S1Angle a);
  friend S1Angle operator+(S1Angle a, S1Angle b);
  friend S1Angle operator-(S1Angle a, S1Angle b);
  friend S1Angle operator*(double m, S1Angle a);
  friend S1Angle operator*(S1Angle a, double m);
  friend S1Angle operator/(S1Angle a, double m);
  friend double operator/(S1Angle a, S1Angle b);
  S1Angle& operator+=(S1Angle a);
  S1Angle& operator-=(S1Angle a);
  S1Angle& operator*=(double m);
  S1Angle& operator/=(double m);

  // Trigonmetric functions (not necessary but slightly more convenient).
  friend double sin(S1Angle a);
  friend double cos(S1Angle a);
  friend double tan(S1Angle a);

  // Return the angle normalized to the range (-180, 180] degrees.
  S1Angle Normalized() const;

  // Normalize this angle to the range (-180, 180] degrees.
  void Normalize();
};
```

### S2Point

The `S2Point` class represents a point on the unit sphere as a 3D
vector.  Usually points are normalized to be unit length, but some
methods do not require this.  The `S2Point` class is simply a synonym
for the `Vector_3d` class from `util/math/vector3-inl.h`, which defines
overloaded operators that make it convenient to write arithmetic
expressions (e.g. `x*p1 + (1-x)*p2`).

Some of its more useful methods include:

```c++
class S2Point /* Vector3_d */ {
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
```

Note that the `S2Point` type is technically defined as follows:

```c++
typedef Vector3<double> Vector3_d;
typedef Vector3_d S2Point;
```

See `util/math/vector3-inl.h` for details and additional methods.

### S2LatLng

The `S2LatLng` class represents a point on the unit sphere as a pair of
latitude-longitude coordinates.

```c++
class S2LatLng {
 public:
  // Constructor.  The latitude and longitude are allowed to be outside
  // the is_valid() range.  However, note that most methods that accept
  // S2LatLngs expect them to be normalized (see Normalized() below).
  S2LatLng(S1Angle lat, S1Angle lng);

  // The default constructor sets the latitude and longitude to zero.  This is
  // mainly useful when declaring arrays, STL containers, etc.
  S2LatLng();

  // Convert a direction vector (not necessarily unit length) to an S2LatLng.
  explicit S2LatLng(S2Point const& p);

  // Returns an S2LatLng for which is_valid() will return false.
  static S2LatLng Invalid();

  // Convenience functions -- shorter than calling S1Angle::Radians(), etc.
  static S2LatLng FromRadians(double lat_radians, double lng_radians);
  static S2LatLng FromDegrees(double lat_degrees, double lng_degrees);
  static S2LatLng FromE5(int32 lat_e5, int32 lng_e5);
  static S2LatLng FromE6(int32 lat_e6, int32 lng_e6);
  static S2LatLng FromE7(int32 lat_e7, int32 lng_e7);

  // Convenience functions -- to use when args have been fixed32s in protos.
  //
  // The arguments are static_cast into int32, so very large unsigned values
  // are treated as negative numbers.
  static S2LatLng FromUnsignedE6(uint32 lat_e6, uint32 lng_e6);
  static S2LatLng FromUnsignedE7(uint32 lat_e7, uint32 lng_e7);

  // Methods to compute the latitude and longitude of a point separately.
  static S1Angle Latitude(S2Point const& p);
  static S1Angle Longitude(S2Point const& p);

  // Accessor methods.
  S1Angle lat() const;
  S1Angle lng() const;
  R2Point const& coords() const;

  // Return true if the latitude is between -90 and 90 degrees inclusive
  // and the longitude is between -180 and 180 degrees inclusive.
  bool is_valid() const;

  // Clamps the latitude to the range [-90, 90] degrees, and adds or subtracts
  // a multiple of 360 degrees to the longitude if necessary to reduce it to
  // the range [-180, 180].
  S2LatLng Normalized() const;

  // Convert a normalized S2LatLng to the equivalent unit-length vector.
  S2Point ToPoint() const;

  // Return the distance (measured along the surface of the sphere) to the
  // given S2LatLng.  This is mathematically equivalent to:
  //
  //   S1Angle(ToPoint(), o.ToPoint())
  //
  // but this implementation is slightly more efficient.  Both S2LatLngs
  // must be normalized.
  S1Angle GetDistance(S2LatLng const& o) const;

  // Simple arithmetic operations for manipulating latitude-longitude pairs.
  // The results are not normalized (see Normalized()).
  friend S2LatLng operator+(S2LatLng const& a, S2LatLng const& b);
  friend S2LatLng operator-(S2LatLng const& a, S2LatLng const& b);
  friend S2LatLng operator*(double m, S2LatLng const& a);
  friend S2LatLng operator*(S2LatLng const& a, double m);

  bool operator==(S2LatLng const& o) const;
  bool operator!=(S2LatLng const& o) const;
  bool operator<(S2LatLng const& o) const;
  bool operator>(S2LatLng const& o) const;
  bool operator<=(S2LatLng const& o) const;
  bool operator>=(S2LatLng const& o) const;

  bool ApproxEquals(S2LatLng const& o, double max_error = 1e-15) const;
};
```

### S2Region

An `S2Region` represents a two-dimensional region over the unit sphere.
It is an abstract interface with various concrete subtypes.

The main purpose of this interface is to allow complex regions to be
approximated as simpler regions.  So rather than having a wide variety
of virtual methods that are implemented by all subtypes, the interface
is restricted to methods that are useful for computing approximations.
Here they are:

```c++
class S2Region {
 public:
  virtual ~S2Region();

  // Return a deep copy of this region.  If you want to narrow the result to a
  // specific known region type, use down_cast<T*> from casts.h.
  // Subtypes return pointers to that subtype from their Clone() methods.
  virtual S2Region* Clone() const = 0;

  // Return a bounding spherical cap. This is not guaranteed to be exact.
  virtual S2Cap GetCapBound() const = 0;

  // Return a bounding latitude-longitude rectangle that contains the region.
  // The bounds are not guaranteed to be tight.
  virtual S2LatLngRect GetRectBound() const = 0;

  // If this method returns true, the region completely contains the given
  // cell.  Otherwise, either the region does not contain the cell or the
  // containment relationship could not be determined.
  virtual bool Contains(S2Cell const& cell) const = 0;

  // If this method returns false, the region does not intersect the given
  // cell.  Otherwise, either region intersects the cell, or the intersection
  // relationship could not be determined.
  virtual bool MayIntersect(S2Cell const& cell) const = 0;

  // Return true if and only if the given point is contained by the region.
  // The point 'p' is generally required to be unit length, although some
  // subtypes may relax this restriction.
  //
  // NOTE: If you will be calling this function on one specific subtype only,
  // or if performance is a consideration, please use the non-virtual
  // method Contains(S2Point const& p) declared below!
  virtual bool VirtualContainsPoint(S2Point const& p) const = 0;

  // Use encoder to generate a serialized representation of this region.
  // Assumes that encoder can be enlarged using calls to Ensure(int).
  //
  // The representation chosen is left up to the sub-classes but it should
  // satisfy the following constraints:
  // - It should encode a version number.
  // - It should be deserializable using the corresponding Decode method.
  // - Performance, not space, should be the chief consideration. Encode() and
  //   Decode() should be implemented such that the combination is equivalent
  //   to calling Clone().
  virtual void Encode(Encoder* const encoder) const = 0;

  // Reconstruct a region from the serialized representation generated by
  // Encode(). Note that since this method is virtual, it requires that a
  // Region object of the appropriate concrete type has already been
  // constructed. It is not possible to decode regions of unknown type.
  //
  // Whenever the Decode method is changed to deal with new serialized
  // representations, it should be done so in a manner that allows for
  // older versions to be decoded i.e. the version number in the serialized
  // representation should be used to decide how to decode the data.
  //
  // Returns true on success.
  virtual bool Decode(Decoder* const decoder) = 0;

  // Provide the same functionality as Decode, except that decoded regions are
  // allowed to point directly into the Decoder's memory buffer rather than
  // copying the data.  This method can be much faster for regions that have
  // a lot of data (such as polygons), but the decoded region is only valid
  // within the scope (lifetime) of the Decoder's memory buffer.
  // Default implementation just calls Decode.
  virtual bool DecodeWithinScope(Decoder* const decoder);
};
```

In addition, all `S2Region` subtypes implement a
`Contains(S2Point const& p)` method that returns true if the given point
`p` is contained by the region.  The point is generally required to be unit
length, although some subtypes may relax this restriction.


### S2LatLngRect

An `S2LatLngRect` is a type of `S2Region` that represents a rectangle in
latitude-longitude space.  It is capable of representing the empty and
full rectangles as well as single points.  It has an `AddPoint` method
that makes it easy to construct a bounding rectangle for a set of points,
including point sets that span the 180 degree meridian.  Here are its
methods:

```c++
class S2LatLngRect : public S2Region {
 public:
  // Construct a rectangle from minimum and maximum latitudes and longitudes.
  // If lo.lng() > hi.lng(), the rectangle spans the 180 degree longitude
  // line. Both points must be normalized, with lo.lat() <= hi.lat().
  // The rectangle contains all the points p such that 'lo' <= p <= 'hi',
  // where '<=' is defined in the obvious way.
  S2LatLngRect(S2LatLng const& lo, S2LatLng const& hi);

  // Construct a rectangle from latitude and longitude intervals.  The two
  // intervals must either be both empty or both non-empty, and the latitude
  // interval must not extend outside [-90, +90] degrees.
  // Note that both intervals (and hence the rectangle) are closed.
  S2LatLngRect(R1Interval const& lat, S1Interval const& lng);

  // The default constructor creates an empty S2LatLngRect.
  S2LatLngRect();

  // Construct a rectangle of the given size centered around the given point.
  // "center" needs to be normalized, but "size" does not.  The latitude
  // interval of the result is clamped to [-90,90] degrees, and the longitude
  // interval of the result is Full() if and only if the longitude size is
  // 360 degrees or more.  Examples of clamping (in degrees):
  //
  //   center=(80,170),  size=(40,60)   -> lat=[60,90],   lng=[140,-160]
  //   center=(10,40),   size=(210,400) -> lat=[-90,90],  lng=[-180,180]
  //   center=(-90,180), size=(20,50)   -> lat=[-90,-80], lng=[155,-155]
  static S2LatLngRect FromCenterSize(S2LatLng const& center,
                                     S2LatLng const& size);

  // Construct a rectangle containing a single (normalized) point.
  static S2LatLngRect FromPoint(S2LatLng const& p);

  // Construct the minimal bounding rectangle containing the two given
  // normalized points.  This is equivalent to starting with an empty
  // rectangle and calling AddPoint() twice.  Note that it is different than
  // the S2LatLngRect(lo, hi) constructor, where the first point is always
  // used as the lower-left corner of the resulting rectangle.
  static S2LatLngRect FromPointPair(S2LatLng const& p1, S2LatLng const& p2);

  // Accessor methods.
  S1Angle lat_lo() const;
  S1Angle lat_hi() const;
  S1Angle lng_lo() const;
  S1Angle lng_hi() const;
  R1Interval const& lat() const;
  S1Interval const& lng() const;
  R1Interval *mutable_lat();
  S1Interval *mutable_lng();
  S2LatLng lo() const;
  S2LatLng hi() const;

  // The canonical empty and full rectangles, as derived from the Empty
  // and Full R1 and S1 Intervals.
  // Empty: lat_lo=1, lat_hi=0, lng_lo=Pi, lng_hi=-Pi (radians)
  static S2LatLngRect Empty();
  // Full: lat_lo=-Pi/2, lat_hi=Pi/2, lng_lo=-Pi, lng_hi=Pi (radians)
  static S2LatLngRect Full();

  // The full allowable range of latitudes and longitudes.
  static R1Interval FullLat();
  static S1Interval FullLng();

  // Return true if the rectangle is valid, which essentially just means
  // that the latitude bounds do not exceed Pi/2 in absolute value and
  // the longitude bounds do not exceed Pi in absolute value.  Also, if
  // either the latitude or longitude bound is empty then both must be.
  bool is_valid() const;

  // Return true if the rectangle is empty, i.e. it contains no points at all.
  bool is_empty() const;

  // Return true if the rectangle is full, i.e. it contains all points.
  bool is_full() const;

  // Return true if the rectangle is a point, i.e. lo() == hi()
  bool is_point() const;

  // Return true if lng_.lo() > lng_.hi(), i.e. the rectangle crosses
  // the 180 degree longitude line.
  bool is_inverted() const { return lng_.is_inverted(); }

  // Return the k-th vertex of the rectangle (k = 0,1,2,3) in CCW order (lower
  // left, lower right, upper right, upper left).
  S2LatLng GetVertex(int k) const;

  // Return the center of the rectangle in latitude-longitude space
  // (in general this is not the center of the region on the sphere).
  S2LatLng GetCenter() const;

  // Return the width and height of this rectangle in latitude-longitude
  // space.  Empty rectangles have a negative width and height.
  S2LatLng GetSize() const;

  // Returns the surface area of this rectangle on the unit sphere.
  double Area() const;

  // Return the true centroid of the rectangle multiplied by its surface area
  // (see s2.h for details on centroids).  The result is not unit length, so
  // you may want to normalize it.  Note that in general the centroid is
  // *not* at the center of the rectangle, and in fact it may not even be
  // contained by the rectangle.  (It is the "center of mass" of the rectangle
  // viewed as subset of the unit sphere, i.e. it is the point in space about
  // which this curved shape would rotate.)
  //
  // The reason for multiplying the result by the rectangle area is to make it
  // easier to compute the centroid of more complicated shapes.  The centroid
  // of a union of disjoint regions can be computed simply by adding their
  // GetCentroid() results.
  S2Point GetCentroid() const;

  // More efficient version of Contains() that accepts a S2LatLng rather than
  // an S2Point.  The argument must be normalized.
  bool Contains(S2LatLng const& ll) const;

  // Return true if and only if the given point is contained in the interior
  // of the region (i.e. the region excluding its boundary).  The point 'p'
  // does not need to be normalized.
  bool InteriorContains(S2Point const& p) const;

  // More efficient version of InteriorContains() that accepts a S2LatLng
  // rather than an S2Point.  The argument must be normalized.
  bool InteriorContains(S2LatLng const& ll) const;

  // Return true if and only if the rectangle contains the given other
  // rectangle.
  bool Contains(S2LatLngRect const& other) const;

  // Return true if and only if the interior of this rectangle contains all
  // points of the given other rectangle (including its boundary).
  bool InteriorContains(S2LatLngRect const& other) const;

  // Return true if this rectangle and the given other rectangle have any
  // points in common.
  bool Intersects(S2LatLngRect const& other) const;

  // Returns true if this rectangle intersects the given cell.  (This is an
  // exact test and may be fairly expensive, see also MayIntersect below.)
  bool Intersects(S2Cell const& cell) const;

  // Return true if and only if the interior of this rectangle intersects
  // any point (including the boundary) of the given other rectangle.
  bool InteriorIntersects(S2LatLngRect const& other) const;

  // Increase the size of the bounding rectangle to include the given point.
  // The rectangle is expanded by the minimum amount possible.  The S2LatLng
  // argument must be normalized.
  void AddPoint(S2Point const& p);
  void AddPoint(S2LatLng const& ll);

  // Return a rectangle that has been expanded by margin.lat() on each side in
  // the latitude direction, and by margin.lng() on each side in the longitude
  // direction.  If either margin is negative, then shrink the rectangle on
  // the corresponding sides instead.  The resulting rectangle may be empty.
  //
  // As noted above, the latitude-longitude space has the topology of a
  // cylinder.  Longitudes "wrap around" at +/-180 degrees, while latitudes
  // are clamped to range [-90, 90].  This means that any expansion (positive
  // or negative) of the full longitude range remains full (since the
  // "rectangle" is actually a continuous band around the cylinder), while
  // expansion of the full latitude range remains full only if the margin is
  // positive.
  //
  // If either the latitude or longitude interval becomes empty after
  // expansion by a negative margin, the result is empty.
  //
  // Note that if an expanded rectangle contains a pole, it may not contain
  // all possible lat/lng representations of that pole (see header above).
  // Use the PolarClosure() method if you do not want this behavior.
  //
  // If you are trying to grow a rectangle by a certain *distance* on the
  // sphere (e.g. 5km), use the ExpandedByDistance() method instead.
  S2LatLngRect Expanded(S2LatLng const& margin) const;

  // If the rectangle does not include either pole, return it unmodified.
  // Otherwise expand the longitude range to Full() so that the rectangle
  // contains all possible representations of the contained pole(s).
  S2LatLngRect PolarClosure() const;

  // Return the smallest rectangle containing the union of this rectangle and
  // the given rectangle.
  S2LatLngRect Union(S2LatLngRect const& other) const;

  // Return the smallest rectangle containing the intersection of this
  // rectangle and the given rectangle.  Note that the region of intersection
  // may consist of two disjoint rectangles, in which case a single rectangle
  // spanning both of them is returned.
  S2LatLngRect Intersection(S2LatLngRect const& other) const;

  // Expand this rectangle so that it contains all points within the given
  // distance of the boundary, and return the smallest such rectangle.  If the
  // distance is negative, then instead shrink this rectangle so that it
  // excludes all points within the given absolute distance of the boundary,
  // and return the largest such rectangle.
  //
  // Unlike Expanded(), this method treats the rectangle as a set of points on
  // the sphere, and measures distances on the sphere.  For example, you can
  // use this method to find a rectangle that contains all points within 5km
  // of a given rectangle.  Because this method uses the topology of the
  // sphere, note the following:
  //
  //  - The full and empty rectangles have no boundary on the sphere.  Any
  //    expansion (positive or negative) of these rectangles leaves them
  //    unchanged.
  //
  //  - Any rectangle that covers the full longitude range does not have an
  //    east or west boundary, therefore no expansion (positive or negative)
  //    will occur in that direction.
  //
  //  - Any rectangle that covers the full longitude range and also includes
  //    a pole will not be expanded or contracted at that pole, because it
  //    does not have a boundary there.
  //
  //  - If a rectangle is within the given distance of a pole, the result will
  //    include the full longitude range (because all longitudes are present
  //    at the poles).
  //
  // Expansion and contraction are defined such that they are inverses whenver
  // possible, i.e.
  //
  //   rect.ExpandedByDistance(x).ExpandedByDistance(-x) == rect
  //
  // (approximately), so long as the first operation does not cause a
  // rectangle boundary to disappear (i.e., the longitude range newly becomes
  // full or empty, or the latitude range expands to include a pole).
  S2LatLngRect ExpandedByDistance(S1Angle distance) const;

  // Returns the minimum distance (measured along the surface of the sphere) to
  // the given S2LatLngRect. Both S2LatLngRects must be non-empty.
  S1Angle GetDistance(S2LatLngRect const& other) const;

  // Returns the minimum distance (measured along the surface of the sphere)
  // from a given point to the rectangle (both its boundary and its interior).
  // The latlng must be valid.
  S1Angle GetDistance(S2LatLng const& p) const;

  // Returns the (directed or undirected) Hausdorff distance (measured along the
  // surface of the sphere) to the given S2LatLngRect. The directed Hausdorff
  // distance from rectangle A to rectangle B is given by
  //     h(A, B) = max_{p in A} min_{q in B} d(p, q).
  // The Hausdorff distance between rectangle A and rectangle B is given by
  //     H(A, B) = max{h(A, B), h(B, A)}.
  S1Angle GetDirectedHausdorffDistance(S2LatLngRect const& other) const;
  S1Angle GetHausdorffDistance(S2LatLngRect const& other) const;

  // Return true if two rectangles contains the same set of points.
  bool operator==(S2LatLngRect const& other) const;

  // Return the opposite of what operator == returns.
  bool operator!=(S2LatLngRect const& other) const;

  // Return true if the latitude and longitude intervals of the two rectangles
  // are the same up to the given tolerance (see r1interval.h and s1interval.h
  // for details).
  bool ApproxEquals(S2LatLngRect const& other, double max_error = 1e-15) const;

  // ApproxEquals() with separate tolerances for latitude and longitude.
  bool ApproxEquals(S2LatLngRect const& other, S2LatLng const& max_error) const;

  // Return true if the edge AB intersects the given edge of constant
  // longitude.
  static bool IntersectsLngEdge(S2Point const& a, S2Point const& b,
                                R1Interval const& lat, double lng);

  // Return true if the edge AB intersects the given edge of constant
  // latitude.  Requires the vectors to have unit length.
  static bool IntersectsLatEdge(S2Point const& a, S2Point const& b,
                                double lat, S1Interval const& lng);

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  virtual S2LatLngRect* Clone() const;
  virtual S2Cap GetCapBound() const;
  virtual S2LatLngRect GetRectBound() const;
  virtual bool Contains(S2Cell const& cell) const;
  virtual bool VirtualContainsPoint(S2Point const& p);

  // This test is cheap but is NOT exact.  Use Intersects() if you want a more
  // accurate and more expensive test.  Note that when this method is used by
  // an S2RegionCoverer, the accuracy isn't all that important since if a cell
  // may intersect the region then it is subdivided, and the accuracy of this
  // method goes up as the cells get smaller.
  virtual bool MayIntersect(S2Cell const& cell) const;

  // The point 'p' does not need to be normalized.
  bool Contains(S2Point const& p) const;

  virtual void Encode(Encoder* const encoder) const;
  virtual bool Decode(Decoder* const decoder);
};
```

### S2Cap

An `S2Cap` represents a spherical cap, i.e. a portion of a sphere cut off
by a plane.  The cap is defined by its axis and height.  This
representation has good numerical accuracy for very small caps, unlike the
(axis, min-distance-from-origin) representation, and is also efficient for
containment tests, unlike the (axis, angle) representation.

Here are the methods available:

```c++
class S2Cap : public S2Region {
 public:
  // The default constructor returns an empty S2Cap.
  S2Cap();

  // Construct a cap with the given center and radius.  A negative radius
  // yields an empty cap; a radius of 180 degrees or more yields a full cap
  // (containing the entire sphere).  "center" should be unit length.
  S2Cap(S2Point const& center, S1Angle radius);

  // Convenience function that creates a cap containing a single point.  This
  // method is more efficient that the S2Cap(center, radius) constructor.
  static S2Cap FromPoint(S2Point const& center);

  // Return a cap with the given center and height (see comments above).  A
  // negative height yields an empty cap; a height of 2 or more yields a full
  // cap.  "center" should be unit length.
  static S2Cap FromCenterHeight(S2Point const& center, double height);

  // Return a cap with the given center and surface area.  Note that the area
  // can also be interpreted as the solid angle subtended by the cap (because
  // the sphere has unit radius).  A negative area yields an empty cap; an
  // area of 4*Pi or more yields a full cap.  "center" should be unit length.
  static S2Cap FromCenterArea(S2Point const& center, double area);

  // Return an empty cap, i.e. a cap that contains no points.
  static S2Cap Empty();

  // Return a full cap, i.e. a cap that contains all points.
  static S2Cap Full();

  ~S2Cap();

  // Accessor methods.
  S2Point const& center() const;
  double height() const;

  // Return the cap radius.  (This method is relatively expensive since only
  // the cap height is stored, and may yield a slightly different result than
  // the value passed to the (center, radius) constructor.)
  S1Angle GetRadius() const;

  // Return the area of the cap.
  double GetArea() const;

  // Return the true centroid of the cap multiplied by its surface area
  // (see s2.h for details on centroids). The result lies on the ray from
  // the origin through the cap's center, but it is not unit length. Note
  // that if you just want the "surface centroid", i.e. the normalized result,
  // then it is much simpler just to call center().
  //
  // The reason for multiplying the result by the cap area is to make it
  // easier to compute the centroid of more complicated shapes.  The centroid
  // of a union of disjoint regions can be computed simply by adding their
  // GetCentroid() results. Caveat: for caps that contain a single point
  // (i.e., zero radius), this method always returns the origin (0, 0, 0).
  // This is because shapes with no area don't affect the centroid of a
  // union whose total area is positive.
  S2Point GetCentroid() const;

  // We allow negative heights (to represent empty caps) but heights are
  // normalized so that they do not exceed 2.
  bool is_valid() const;

  // Return true if the cap is empty, i.e. it contains no points.
  bool is_empty() const;

  // Return true if the cap is full, i.e. it contains all points.
  bool is_full() const;

  // Return the complement of the interior of the cap.  A cap and its
  // complement have the same boundary but do not share any interior points.
  // The complement operator is not a bijection because the complement of a
  // singleton cap (containing a single point) is the same as the complement
  // of an empty cap.
  S2Cap Complement() const;

  // Return true if and only if this cap contains the given other cap
  // (in a set containment sense, e.g. every cap contains the empty cap).
  bool Contains(S2Cap const& other) const;

  // Return true if and only if this cap intersects the given other cap,
  // i.e. whether they have any points in common.
  bool Intersects(S2Cap const& other) const;

  // Return true if and only if the interior of this cap intersects the
  // given other cap.  (This relationship is not symmetric, since only
  // the interior of this cap is used.)
  bool InteriorIntersects(S2Cap const& other) const;

  // Return true if and only if the given point is contained in the interior
  // of the cap (i.e. the cap excluding its boundary).  "p" should be be a
  // unit-length vector.
  bool InteriorContains(S2Point const& p) const;

  // Increase the cap height if necessary to include the given point.  If the
  // cap is empty then the center is set to the given point, but otherwise the
  // center is not changed.  "p" should be a unit-length vector.
  void AddPoint(S2Point const& p);

  // Increase the cap height if necessary to include "other".  If the current
  // cap is empty it is set to the given other cap.
  void AddCap(S2Cap const& other);

  // Return a cap that contains all points within a given distance of this
  // cap.  Note that any expansion of the empty cap is still empty.
  S2Cap Expanded(S1Angle distance) const;

  // Return the cap height corresponding to the given radius.  This method can
  // be used to efficiently construct many caps of the same radius (since
  // converting the radius to a height is relatively expensive).  Negative
  // radii are mapped to negative heights (i.e., empty caps), while radii of
  // of 180 degrees or more are mapped to a height of 2 (i.e., a full cap).
  static double RadiusToHeight(S1Angle radius);

  // Return the smallest cap which encloses this cap and "other".
  S2Cap Union(S2Cap const& other) const;

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  virtual S2Cap* Clone() const;
  virtual S2Cap GetCapBound() const;
  virtual S2LatLngRect GetRectBound() const;
  virtual bool Contains(S2Cell const& cell) const;
  virtual bool MayIntersect(S2Cell const& cell) const;
  virtual bool VirtualContainsPoint(S2Point const& p);

  // The point "p" should be a unit-length vector.
  bool Contains(S2Point const& p) const;

  // Encoding/decoding not implemented.
  virtual void Encode(Encoder* const encoder) const;
  virtual bool Decode(Decoder* const decoder);

  ///////////////////////////////////////////////////////////////////////
  // The following static methods are convenience functions for assertions
  // and testing purposes only.

  // Return true if two caps are identical.
  bool operator==(S2Cap const& other) const;

  // Return true if the cap center and height differ by at most "max_error"
  // from the given cap "other".
  bool ApproxEquals(S2Cap const& other, double max_error = 1e-14) const;
};
```

## Cell Hierarchy

The library provides methods for subdividing the sphere into a hierarchical
collection of "cells".  The cells have a space-filling curve structure
that makes them useful for spatial indexing.  There are methods for
approximating arbitrary regions as a collection of cells.

The S2 cell structure is defined as follows.  There are six top-level
*face cells*, obtained by projecting the six faces of a cube -- (face, u, v)
coordinates below -- onto the unit sphere -- (x, y, z) coordinates below.
We call this cube the uv-cube.

Each face is then subdivided recursively into four cells
in a quadtree-like fashion. On the uv-cube, a cell is a rectangle whose edges
are aligned with the sides of its face.
On the sphere, it is a spherical quadrilateral bounded by four geodesics
(great circle segments).

There are a total of 30
levels of subdivision defined (i.e. 6 * 4<sup>30</sup> leaf cells),
which gives a resolution of about 1cm
everywhere on a sphere the size of the earth. Details on the cell areas at
each level appear on the [S2 Cell Statistics page](cell_statistics.md).
Each cell is uniquely identified by a 64-bit *cell id*.

### Coordinate Systems
------------------

In order for cells to be roughly the same size on the sphere, they are not
the same size on the uv-cube.  We introduce another cube -- the st-cube.
On this cube the cells are perfectly square and divided through their center.
The st-cube is projected on the uv-cube so that cells at the periphery of an
st-face are larger than cells at their center.
In a 2-d projection, it looks like this:

<a href="xyz_to_uv_to_st.jpg">
<img src="xyz_to_uv_to_st.jpg" alt="xyz_to_uv_to_st.jpg"  width="320" height="240"  />
</a>

Note that a cell on the st-cube is not bounded by straight lines.

There are a few more coordinate systems worth introducing:

* (id)<br> Cell id.  A 64-bit encoding of a face and a
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

### Cell Ids

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

### Cell Id to ST Coordinates

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

### `S2CellId` Class

The `S2CellId` class is a thin wrapper over a 64-bit cell id that
provides methods for navigating the cell hierarchy (finding parents,
children, containment tests, etc).  Since leaf cells are often used
to represent points on the unit sphere, the `S2CellId` class also
provides methods for converting directly to and from an `S2Point`.
Here are its methods:

```c++
class S2CellId {
 public:
  static int const kFaceBits = 3;
  static int const kNumFaces = 6;
  static int const kMaxLevel = 30;  // Valid levels: 0..kMaxLevel
  static int const kMaxSize = 1 << kMaxLevel;

  // Although only 60 bits are needed to represent the index of a leaf
  // cell, we need an extra bit in order to represent the position of
  // the center of the leaf cell along the Hilbert curve.
  static int const kPosBits = 2 * kMaxLevel + 1;

  explicit S2CellId(uint64 id) : id_(id) {}

  // The default constructor returns an invalid cell id.
  S2CellId();
  static S2CellId None();

  // Return an invalid cell id guaranteed to be larger than any
  // valid cell id.  Useful for creating indexes.
  static S2CellId Sentinel();

  // Return the cell corresponding to a given S2 cube face.
  static S2CellId FromFace(int face);

  // Return a cell given its face (range 0..5), Hilbert curve position within
  // that face (an unsigned integer with S2CellId::kPosBits bits), and level
  // (range 0..kMaxLevel).  The given position will be modified to correspond
  // to the Hilbert curve position at the center of the returned cell.  This
  // is a static function rather than a constructor in order to indicate what
  // the arguments represent.
  static S2CellId FromFacePosLevel(int face, uint64 pos, int level);

  // Return a leaf cell containing the given point "p".  Usually there is
  // exactly one such cell, but for points along the edge of a cell, any
  // adjacent cell may be (deterministically) chosen.  This is because
  // S2CellIds are considered to be closed sets.  The returned cell will
  // always contain the given point, i.e.
  //
  //   S2Cell(S2CellId::FromPoint(p)).Contains(p)
  //
  // is always true.  The point "p" does not need to be normalized.
  static S2CellId FromPoint(S2Point const& p);

  // Return the leaf cell containing the given normalized S2LatLng.
  static S2CellId FromLatLng(S2LatLng const& ll);

  // Return the direction vector corresponding to the center of the given
  // cell.  The vector returned by ToPointRaw is not necessarily unit length.
  // This method returns the same result as S2Cell::GetCenter().
  S2Point ToPoint() const;
  S2Point ToPointRaw() const;

  // Return the center of the cell in (s,t) coordinates (see s2.h).
  R2Point GetCenterST() const;

  // Return the edge length of this cell in (s,t)-space.
  double GetSizeST() const;

  // Return the edge length in (s,t)-space of cells at the given level.
  static double GetSizeST(int level);

  // Return the bounds of this cell in (s,t)-space.
  R2Rect GetBoundST() const;

  // Return the center of the cell in (u,v) coordinates (see s2.h).  Note that
  // the center of the cell is defined as the point at which it is recursively
  // subdivided into four children; in general, it is not at the midpoint of
  // the (u,v) rectangle covered by the cell.
  R2Point GetCenterUV() const;

  // Return the bounds of this cell in (u,v)-space.
  R2Rect GetBoundUV() const;

  // Return the (face, si, ti) coordinates of the center of the cell.  Note
  // that although (si,ti) coordinates span the range [0,2**31] in general,
  // the cell center coordinates are always in the range [1,2**31-1] and
  // therefore can be represented using a signed 32-bit integer.
  int GetCenterSiTi(int* psi, int* pti) const;

  // Return the S2LatLng corresponding to the center of the given cell.
  S2LatLng ToLatLng() const;

  // The 64-bit unique identifier for this cell.
  uint64 id() const;

  // Return true if id() represents a valid cell.
  bool is_valid() const;

  // Which cube face this cell belongs to, in the range 0..5.
  int face() const;

  // The position of the cell center along the Hilbert curve over this face,
  // in the range 0..(2**kPosBits-1).
  uint64 pos() const;

  // Return the subdivision level of the cell (range 0..kMaxLevel).
  int level() const;

  // Return the edge length of this cell in (i,j)-space.
  int GetSizeIJ() const;

  // Like the above, but return the size of cells at the given level.
  static int GetSizeIJ(int level);

  // Return true if this is a leaf cell (more efficient than checking
  // whether level() == kMaxLevel).
  bool is_leaf() const;

  // Return true if this is a top-level face cell (more efficient than
  // checking whether level() == 0).
  bool is_face() const;

  // Return the child position (0..3) of this cell within its parent.
  // REQUIRES: level() >= 1.
  int child_position() const;

  // Return the child position (0..3) of this cell's ancestor at the given
  // level within its parent.  For example, child_position(1) returns the
  // position of this cell's level-1 ancestor within its top-level face cell.
  // REQUIRES: 1 <= level <= this->level().
  int child_position(int level) const;

  // Methods that return the range of cell ids that are contained
  // within this cell (including itself).  The range is *inclusive*
  // (i.e. test using >= and <=) and the return values of both
  // methods are valid leaf cell ids.
  //
  // These methods should not be used for iteration.  If you want to
  // iterate through all the leaf cells, call child_begin(kMaxLevel) and
  // child_end(kMaxLevel) instead.  Also see maximum_tile(), which can be used
  // to iterate through cell ranges using cells at different levels.
  //
  // It would in fact be error-prone to define a range_end() method, because
  // this method would need to return (range_max().id() + 1) which is not
  // always a valid cell id.  This also means that iterators would need to be
  // tested using "<" rather that the usual "!=".
  S2CellId range_min() const;
  S2CellId range_max() const;

  // Return true if the given cell is contained within this one.
  bool contains(S2CellId other) const;

  // Return true if the given cell intersects this one.
  bool intersects(S2CellId other) const;

  // Return the cell at the previous level or at the given level (which must
  // be less than or equal to the current level).
  S2CellId parent() const;
  S2CellId parent(int level) const;

  // Return the immediate child of this cell at the given traversal order
  // position (in the range 0 to 3).  This cell must not be a leaf cell.
  S2CellId child(int position) const;

  // Iterator-style methods for traversing the immediate children of a cell or
  // all of the children at a given level (greater than or equal to the current
  // level).  Note that the end value is exclusive, just like standard STL
  // iterators, and may not even be a valid cell id.  You should iterate using
  // code like this:
  //
  //   for(S2CellId c = id.child_begin(); c != id.child_end(); c = c.next())
  //     ...
  //
  // The convention for advancing the iterator is "c = c.next()" rather
  // than "++c" to avoid possible confusion with incrementing the
  // underlying 64-bit cell id.
  S2CellId child_begin() const;
  S2CellId child_begin(int level) const;
  S2CellId child_end() const;
  S2CellId child_end(int level) const;

  // Return the next/previous cell at the same level along the Hilbert curve.
  // Works correctly when advancing from one face to the next, but
  // does *not* wrap around from the last face to the first or vice versa.
  S2CellId next() const;
  S2CellId prev() const;

  // This method advances or retreats the indicated number of steps along the
  // Hilbert curve at the current level, and returns the new position.  The
  // position is never advanced past End() or before Begin().
  S2CellId advance(int64 steps) const;

  // Like next() and prev(), but these methods wrap around from the last face
  // to the first and vice versa.  They should *not* be used for iteration in
  // conjunction with child_begin(), child_end(), Begin(), or End().  The
  // input must be a valid cell id.
  S2CellId next_wrap() const;
  S2CellId prev_wrap() const;

  // This method advances or retreats the indicated number of steps along the
  // Hilbert curve at the current level, and returns the new position.  The
  // position wraps between the first and last faces as necessary.  The input
  // must be a valid cell id.
  S2CellId advance_wrap(int64 steps) const;

  // Return the largest cell with the same range_min() and such that
  // range_max() < limit.range_min().  Returns "limit" if no such cell exists.
  // This method can be used to generate a small set of S2CellIds that covers
  // a given range (a "tiling").  This example shows how to generate a tiling
  // for a semi-open range of leaf cells [start, limit):
  //
  //   for (S2CellId id = start.maximum_tile(limit);
  //        id != limit; id = id.next().maximum_tile(limit)) { ... }
  //
  // Note that in general the cells in the tiling will be of different sizes;
  // they gradually get larger (near the middle of the range) and then
  // gradually get smaller (as "limit" is approached).
  S2CellId maximum_tile(S2CellId limit) const;

  // Return the level of the "lowest common ancestor" of this cell and
  // "other".  Note that because of the way that cell levels are numbered,
  // this is actually the *highest* level of any shared ancestor.  Return -1
  // if the two cells do not have any common ancestor (i.e., they are from
  // different faces).
  int GetCommonAncestorLevel(S2CellId other) const;

  // Iterator-style methods for traversing all the cells along the Hilbert
  // curve at a given level (across all 6 faces of the cube).  Note that the
  // end value is exclusive (just like standard STL iterators), and is not a
  // valid cell id.
  static S2CellId Begin(int level);
  static S2CellId End(int level);

  // Methods to encode and decode cell ids to compact text strings suitable
  // for display or indexing.  Cells at lower levels (i.e. larger cells) are
  // encoded into fewer characters.  The maximum token length is 16.
  //
  // ToToken() returns a string by value for convenience; the compiler
  // does this without intermediate copying in most cases.
  //
  // These methods guarantee that FromToken(ToToken(x)) == x even when
  // "x" is an invalid cell id.  All tokens are alphanumeric strings.
  // FromToken() returns S2CellId::None() for malformed inputs.
  string ToToken() const;
  static S2CellId FromToken(const char* token, size_t length);
  static S2CellId FromToken(string const& token);

  // Creates a debug human readable string. Used for << and available for direct
  // usage as well.
  string ToString() const;

  // Return the four cells that are adjacent across the cell's four edges.
  // Neighbors are returned in the order defined by S2Cell::GetEdge.  All
  // neighbors are guaranteed to be distinct.
  void GetEdgeNeighbors(S2CellId neighbors[4]) const;

  // Return the neighbors of closest vertex to this cell at the given level,
  // by appending them to "output".  Normally there are four neighbors, but
  // the closest vertex may only have three neighbors if it is one of the 8
  // cube vertices.
  //
  // Requires: level < this->level(), so that we can determine which vertex is
  // closest (in particular, level == kMaxLevel is not allowed).
  void AppendVertexNeighbors(int level, std::vector<S2CellId>* output) const;

  // Append all neighbors of this cell at the given level to "output".  Two
  // cells X and Y are neighbors if their boundaries intersect but their
  // interiors do not.  In particular, two cells that intersect at a single
  // point are neighbors.
  //
  // Requires: nbr_level >= this->level().  Note that for cells adjacent to a
  // face vertex, the same neighbor may be appended more than once.
  void AppendAllNeighbors(int nbr_level, std::vector<S2CellId>* output) const;

  /////////////////////////////////////////////////////////////////////
  // Low-level methods.

  // Return a leaf cell given its cube face (range 0..5) and
  // i- and j-coordinates (see s2.h).
  static S2CellId FromFaceIJ(int face, int i, int j);

  // Return the (face, i, j) coordinates for the leaf cell corresponding to
  // this cell id.  Since cells are represented by the Hilbert curve position
  // at the center of the cell, the returned (i,j) for non-leaf cells will be
  // a leaf cell adjacent to the cell center.  If "orientation" is non-NULL,
  // also return the Hilbert curve orientation for the current cell.
  int ToFaceIJOrientation(int* pi, int* pj, int* orientation) const;

  // Return the lowest-numbered bit that is on for this cell id, which is
  // equal to (uint64(1) << (2 * (kMaxLevel - level))).  So for example,
  // a.lsb() <= b.lsb() if and only if a.level() >= b.level(), but the
  // first test is more efficient.
  uint64 lsb() const;

  // Return the lowest-numbered bit that is on for cells at the given level.
  static uint64 lsb_for_level(int level);

  // Return the bound in (u,v)-space for the cell at the given level containing
  // the leaf cell with the given (i,j)-coordinates.
  static R2Rect IJLevelToBoundUV(int ij[2], int level);
};
```

### `S2Cell` Class

An `S2Cell` is an `S2Region` object that represents a cell.  Unlike `S2CellId`,
it views a cell as a representing a spherical quadrilateral rather than a point,
and it supports efficient containment and intersection tests.  However, it is
also a more expensive representation (currently 48 bytes rather than 8).

Here are its methods:

```c++
class S2Cell : public S2Region {
 public:
  // The default constructor is required in order to use freelists.
  // Cells should otherwise always be constructed explicitly.
  S2Cell();

  // An S2Cell always corresponds to a particular S2CellId.  The other
  // constructors are just convenience methods.
  explicit S2Cell(S2CellId id);

  // Return the cell corresponding to the given S2 cube face.
  static S2Cell FromFace(int face);

  // Return a cell given its face (range 0..5), Hilbert curve position within
  // that face (an unsigned integer with S2CellId::kPosBits bits), and level
  // (range 0..kMaxLevel).  The given position will be modified to correspond
  // to the Hilbert curve position at the center of the returned cell.  This
  // is a static function rather than a constructor in order to indicate what
  // the arguments represent.
  static S2Cell FromFacePosLevel(int face, uint64 pos, int level);

  // Convenience methods.  The S2LatLng must be normalized.
  explicit S2Cell(S2Point const& p);
  explicit S2Cell(S2LatLng const& ll);

  S2CellId id() const;
  int face() const;
  int level() const;
  int orientation() const;
  bool is_leaf() const;

  // These are equivalent to the S2CellId methods, but have a more efficient
  // implementation since the level has been precomputed.
  int GetSizeIJ() const;
  double GetSizeST() const;

  // Return the k-th vertex of the cell (k = 0,1,2,3).  Vertices are returned
  // in CCW order (lower left, lower right, upper right, upper left in the UV
  // plane).  The points returned by GetVertexRaw are not normalized.
  S2Point GetVertex(int k) const;
  S2Point GetVertexRaw(int k) const;

  // Return the inward-facing normal of the great circle passing through
  // the edge from vertex k to vertex k+1 (mod 4).  The normals returned
  // by GetEdgeRaw are not necessarily unit length.
  S2Point GetEdge(int k) const;
  S2Point GetEdgeRaw(int k) const;

  // If this is not a leaf cell, set children[0..3] to the four children of
  // this cell (in traversal order) and return true.  Otherwise returns false.
  // This method is equivalent to the following:
  //
  // for (pos=0, id=child_begin(); id != child_end(); id = id.next(), ++pos)
  //   children[pos] = S2Cell(id);
  //
  // except that it is more than two times faster.
  bool Subdivide(S2Cell children[4]) const;

  // Return the direction vector corresponding to the center in (s,t)-space of
  // the given cell.  This is the point at which the cell is divided into four
  // subcells; it is not necessarily the centroid of the cell in (u,v)-space
  // or (x,y,z)-space.  The point returned by GetCenterRaw is not necessarily
  // unit length.
  S2Point GetCenter() const;
  S2Point GetCenterRaw() const;

  // Return the average area for cells at the given level.
  static double AverageArea(int level);

  // Return the average area of cells at this level.  This is accurate to
  // within a factor of 1.7 (for S2_QUADRATIC_PROJECTION) and is extremely
  // cheap to compute.
  double AverageArea() const;

  // Return the approximate area of this cell.  This method is accurate to
  // within 3% percent for all cell sizes and accurate to within 0.1% for
  // cells at level 5 or higher (i.e. squares 350km to a side or smaller
  // on the Earth's surface).  It is moderately cheap to compute.
  double ApproxArea() const;

  // Return the area of this cell as accurately as possible.  This method is
  // more expensive but it is accurate to 6 digits of precision even for leaf
  // cells (whose area is approximately 1e-18).
  double ExactArea() const;

  // Return the bounds of this cell in (u,v)-space.
  R2Rect GetBoundUV() const;

  // Return the distance from the cell to the given point.  Returns zero if
  // the point is inside the cell.
  S1ChordAngle GetDistance(S2Point const& target) const;

  // Return the minimum distance from the cell to the given edge AB.  Returns
  // zero if the edge intersects the cell interior.
  S1ChordAngle GetDistanceToEdge(S2Point const& a, S2Point const& b) const;

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  virtual S2Cell* Clone() const;
  virtual S2Cap GetCapBound() const;
  virtual S2LatLngRect GetRectBound() const;
  virtual bool Contains(S2Cell const& cell) const;
  virtual bool MayIntersect(S2Cell const& cell) const;
  virtual bool VirtualContainsPoint(S2Point const& p) const {
    return Contains(p);  // The same as Contains() below, just virtual.
  }

  // Return true if the cell contains the given point "p".  Note that unlike
  // S2Loop/S2Polygon, S2Cells are considered to be closed sets.  This means
  // that points along an S2Cell edge (or at a vertex) belong to the adjacent
  // cell(s) as well.
  //
  // If instead you want every point to be contained by exactly one S2Cell,
  // you will need to convert the S2Cells to S2Loops (which implement point
  // containment this way).
  //
  // The point "p" does not need to be normalized.
  bool Contains(S2Point const& p) const;

  // Encoding/decoding not implemented;
  virtual void Encode(Encoder* const encoder) const;
  virtual bool Decode(Decoder* const decoder);
```

### Cell Unions

An `S2CellUnion` is an `S2Region` consisting of cells of various sizes.
A cell union is typically used to approximate some other shape.
There is a tradeoff between the accuracy of the approximation and how many
cells are used. Unlike polygons, cells have a fixed hierarchical structure.
This makes them more suitable for spatial indexing.

Here is the interface:

```c++
class S2CellUnion : public S2Region {
 public:
  // The default constructor does nothing.  The cell union cannot be used
  // until one of the Init() methods is called.
  S2CellUnion();

  // Populates a cell union with the given S2CellIds or 64-bit cells ids, and
  // then calls Normalize().  The InitSwap() version takes ownership of the
  // vector data without copying and clears the given vector.  These methods
  // may be called multiple times.
  void Init(vector<S2CellId> const& cell_ids);
  void Init(vector<uint64> const& cell_ids);
  void InitSwap(std::vector<S2CellId>* cell_ids);

  // Like Init(), but does not call Normalize().  The cell union *must* be
  // normalized before doing any calculations with it, so it is the caller's
  // responsibility to make sure that the input is normalized.  This method is
  // useful when converting cell unions to another representation and back.
  // These methods may be called multiple times.
  void InitRaw(std::vector<S2CellId> const& cell_ids);
  void InitRaw(std::vector<uint64> const& cell_ids);
  void InitRawSwap(std::vector<S2CellId>* cell_ids);

  // Gives ownership of the vector data to the client without copying, and
  // clears the content of the cell union.  The original data in cell_ids
  // is lost if there was any.  This is the opposite of InitRawSwap().
  void Detach(std::vector<S2CellId>* cell_ids);

  // Convenience methods for accessing the individual cell ids.
  int num_cells() const;
  S2CellId const& cell_id(int i) const;

  // Direct access to the underlying vector for STL algorithms.
  std::vector<S2CellId> const& cell_ids() const;

  // Normalizes the cell union by discarding cells that are contained by other
  // cells, replacing groups of 4 child cells by their parent cell whenever
  // possible, and sorting all the cell ids in increasing order.  Returns true
  // if the number of cells was reduced.
  //
  // If InitRaw() was used then this method *must* be called before doing any
  // calculations on the cell union, such as Intersects() or Contains().
  bool Normalize();

  // Replaces "output" with an expanded version of the cell union where any
  // cells whose level is less than "min_level" or where (level - min_level)
  // is not a multiple of "level_mod" are replaced by their children, until
  // either both of these conditions are satisfied or the maximum level is
  // reached.
  //
  // This method allows a covering generated by S2RegionCoverer using
  // min_level() or level_mod() constraints to be stored as a normalized cell
  // union (which allows various geometric computations to be done) and then
  // converted back to the original list of cell ids that satisfies the
  // desired constraints.
  void Denormalize(int min_level, int level_mod,
                   std::vector<S2CellId>* output) const;

  // If there are more than "excess" elements of the cell_ids() vector that
  // are allocated but unused, reallocate the array to eliminate the excess
  // space.  This reduces memory usage when many cell unions need to be held
  // in memory at once.
  void Pack(int excess = 0);

  // Return true if the cell union contains the given cell id.  Containment is
  // defined with respect to regions, e.g. a cell contains its 4 children.
  // This is a fast operation (logarithmic in the size of the cell union).
  bool Contains(S2CellId id) const;

  // Return true if the cell union intersects the given cell id.
  // This is a fast operation (logarithmic in the size of the cell union).
  bool Intersects(S2CellId id) const;

  // Return true if this cell union contain/intersects the given other cell
  // union.
  bool Contains(S2CellUnion const* y) const;
  bool Intersects(S2CellUnion const* y) const;

  // Initialize this cell union to the union, intersection, or
  // difference (x - y) of the two given cell unions.
  // Requires: x != this and y != this.
  void GetUnion(S2CellUnion const* x, S2CellUnion const* y);
  void GetIntersection(S2CellUnion const* x, S2CellUnion const* y);
  void GetDifference(S2CellUnion const* x, S2CellUnion const* y);

  // Specialized version of GetIntersection() that gets the intersection of a
  // cell union with the given cell id.  This can be useful for "splitting" a
  // cell union into chunks.
  void GetIntersection(S2CellUnion const* x, S2CellId id);

  // Expands the cell union by adding a "rim" of cells on expand_level
  // around the union boundary.
  //
  // For each cell c in the union, we add all cells at level
  // expand_level that abut c.  There are typically eight of those
  // (four edge-abutting and four sharing a vertex).  However, if c is
  // finer than expand_level, we add all cells abutting
  // c.parent(expand_level) as well as c.parent(expand_level) itself,
  // as an expand_level cell rarely abuts a smaller cell.
  //
  // Note that the size of the output is exponential in
  // "expand_level".  For example, if expand_level == 20 and the input
  // has a cell at level 10, there will be on the order of 4000
  // adjacent cells in the output.  For most applications the
  // Expand(min_radius, max_level_diff) method below is easier to use.
  void Expand(int expand_level);

  // Expand the cell union such that it contains all points whose distance to
  // the cell union is at most "min_radius", but do not use cells that are
  // more than "max_level_diff" levels higher than the largest cell in the
  // input.  The second parameter controls the tradeoff between accuracy and
  // output size when a large region is being expanded by a small amount
  // (e.g. expanding Canada by 1km).  For example, if max_level_diff == 4 the
  // region will always be expanded by approximately 1/16 the width of its
  // largest cell.  Note that in the worst case, the number of cells in the
  // output can be up to 4 * (1 + 2 ** max_level_diff) times larger than the
  // number of cells in the input.
  void Expand(S1Angle min_radius, int max_level_diff);

  // Create a cell union that corresponds to a continuous range of cell ids.
  // The output is a normalized collection of cell ids that covers the leaf
  // cells between "min_id" and "max_id" inclusive.
  // REQUIRES: min_id.is_leaf(), max_id.is_leaf(), min_id <= max_id.
  void InitFromRange(S2CellId min_id, S2CellId max_id);

  // Like InitFromRange(), except that the union covers the range of leaf
  // cells from "begin" (inclusive) to "end" (exclusive), as with Python
  // ranges or STL iterator ranges.  If (begin == end) the result is empty.
  // REQUIRES: begin.is_leaf(), end.is_leaf(), begin <= end.
  void InitFromBeginEnd(S2CellId begin, S2CellId end);

  // The number of leaf cells covered by the union.
  // This will be no more than 6*2^60 for the whole sphere.
  uint64 LeafCellsCovered() const;

  // Approximate this cell union's area by summing the average area of
  // each contained cell's average area, using the AverageArea method
  // from the S2Cell class.
  // This is equivalent to the number of leaves covered, multiplied by
  // the average area of a leaf.
  // Note that AverageArea does not take into account distortion of cell, and
  // thus may be off by up to a factor of 1.7.
  // NOTE: Since this is proportional to LeafCellsCovered(), it is
  // always better to use the other function if all you care about is
  // the relative average area between objects.
  double AverageBasedArea() const;

  // Calculate this cell union's area by summing the approximate area for each
  // contained cell, using the ApproxArea method from the S2Cell class.
  double ApproxArea() const;

  // Calculate this cell union's area by summing the exact area for each
  // contained cell, using the Exact method from the S2Cell class.
  double ExactArea() const;

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  virtual S2CellUnion* Clone() const;
  virtual S2Cap GetCapBound() const;
  virtual S2LatLngRect GetRectBound() const;

  // This is a fast operation (logarithmic in the size of the cell union).
  virtual bool Contains(S2Cell const& cell) const;

  // This is a fast operation (logarithmic in the size of the cell union).
  virtual bool MayIntersect(S2Cell const& cell) const;

  virtual bool VirtualContainsPoint(S2Point const& p);

  // Encoding/decoding not implemented.
  virtual void Encode(Encoder* const encoder);
  virtual bool Decode(Decoder* const decoder);

  // The point 'p' does not need to be normalized.
  // This is a fast operation (logarithmic in the size of the cell union).
  bool Contains(S2Point const& p) const;

  ////////////////////////////////////////////////////////////////////////
  // Static methods intended for high-performance clients that prefer to
  // manage their own storage.

  // Like Normalize() above, but works directly with a vector of S2CellIds.
  // Equivalent to the following:
  //    S2CellUnion cell_union;
  //    cell_union.Init(*cell_ids);
  //    cell_union.Detach(cell_ids);
  static bool Normalize(std::vector<S2CellId>* cell_ids);

  // Like GetIntersection() above, but works directly with vectors of S2CellIds.
  // Equivalent to the following:
  //    S2CellUnion x_union, y_union, result_union;
  //    x_union.Init(x);
  //    y_union.Init(y);
  //    result_union.GetIntersection(&x_union, &y_union);
  //    result_union.Detach(out);
  // except that this method has slightly more relaxed normalization
  // requirements: the input vectors may contain groups of 4 child cells that
  // all have the same parent.  (In a normalized S2CellUnion, such groups are
  // always replaced by the parent cell.)
  static void GetIntersection(std::vector<S2CellId> const& x,
                              std::vector<S2CellId> const& y,
                              std::vector<S2CellId>* out);
};
```

### Approximating Regions

An `S2RegionCoverer` is a class that allows arbitrary regions to be
approximated as unions of cells (`S2CellUnion`).  This is useful for
implementing various sorts of search and precomputation operations.

Typical usage:

    S2RegionCoverer coverer;
    coverer.set_max_cells(5);
    S2Cap cap = S2Cap::FromAxisAngle(...);
    S2CellUnion covering;
    coverer.GetCellUnion(cap, &covering);

The result is a cell union of at most 5 cells that is guaranteed to cover the
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

Here is the interface of `S2RegionCoverer`:

```c++
class S2RegionCoverer {
 public:
  // By default, the covering uses at most 8 cells at any level.  This gives
  // a reasonable tradeoff between the number of cells used and the accuracy
  // of the approximation (see table below).
  static int const kDefaultMaxCells = 8;

  S2RegionCoverer();
  ~S2RegionCoverer();

  // Set the minimum and maximum cell level to be used.  The default is to use
  // all cell levels.  Requires: max_level() >= min_level().
  //
  // To find the cell level corresponding to a given physical distance, use
  // the S2Cell metrics defined in s2.h.  For example, to find the cell
  // level that corresponds to an average edge length of 10km, use:
  //
  //     int level = S2::kAvgEdge.GetClosestLevel(
  //                 geostore::S2Earth::KmToRadians(length_km));
  //
  // Note: min_level() takes priority over max_cells(), i.e. cells below the
  // given level will never be used even if this causes a large number of
  // cells to be returned.
  void set_min_level(int min_level);
  void set_max_level(int max_level);
  int min_level() const;
  int max_level() const;

  // Convenience function that sets both the maximum and minimum cell levels.
  // Note that since min_level() takes priority over max_cells(), an arbitrary
  // number of cells may be returned even if max_cells() is small.
  void set_fixed_level(int level);

  // If specified, then only cells where (level - min_level) is a multiple of
  // "level_mod" will be used (default 1).  This effectively allows the
  // branching factor of the S2CellId hierarchy to be increased.  Currently
  // the only parameter values allowed are 1, 2, or 3, corresponding to
  // branching factors of 4, 16, and 64 respectively.
  void set_level_mod(int level_mod);
  int level_mod() const;

  // Sets the maximum desired number of cells in the approximation (defaults
  // to kDefaultMaxCells).  Note the following:
  //
  //  - For any setting of max_cells(), up to 6 cells may be returned if that
  //    is the minimum number of cells required (e.g. if the region intersects
  //    all six face cells).  Up to 3 cells may be returned even for very tiny
  //    convex regions if they happen to be located at the intersection of
  //    three cube faces.
  //
  //  - For any setting of max_cells(), an arbitrary number of cells may be
  //    returned if min_level() is too high for the region being approximated.
  //
  //  - If max_cells() is less than 4, the area of the covering may be
  //    arbitrarily large compared to the area of the original region even if
  //    the region is convex (e.g. an S2Cap or S2LatLngRect).
  //
  // Accuracy is measured by dividing the area of the covering by the area of
  // the original region.  The following table shows the median and worst case
  // values for this area ratio on a test case consisting of 100,000 spherical
  // caps of random size (generated using s2regioncoverer_unittest):
  //
  //   max_cells:        3      4     5     6     8    12    20   100   1000
  //   median ratio:  5.33   3.32  2.73  2.34  1.98  1.66  1.42  1.11   1.01
  //   worst case:  215518  14.41  9.72  5.26  3.91  2.75  1.92  1.20   1.02
  void set_max_cells(int max_cells);
  int max_cells() const;

  // Return a vector of cell ids that covers (GetCovering) or is contained
  // within (GetInteriorCovering) the given region and satisfies the various
  // restrictions specified above.
  void GetCovering(S2Region const& region, std::vector<S2CellId>* covering);
  void GetInteriorCovering(S2Region const& region,
                           std::vector<S2CellId>* interior);

  // Return a normalized cell union that covers (GetCellUnion) or is contained
  // within (GetInteriorCellUnion) the given region and satisfies the
  // restrictions *EXCEPT* for min_level() and level_mod().  These criteria
  // cannot be satisfied using a cell union because cell unions are
  // automatically normalized by replacing four child cells with their parent
  // whenever possible.  (Note that the list of cell ids passed to the cell
  // union constructor does in fact satisfy all the given restrictions.)
  void GetCellUnion(S2Region const& region, S2CellUnion* covering);
  void GetInteriorCellUnion(S2Region const& region, S2CellUnion* interior);

  // Like GetCovering(), except that this method is much faster and the
  // coverings are not as tight.  All of the usual parameters are respected
  // (max_cells, min_level, max_level, and level_mod), except that the
  // implementation makes no attempt to take advantage of large values of
  // max_cells().  (A small number of cells will always be returned.)
  //
  // This function is useful as a starting point for algorithms that
  // recursively subdivide cells.
  void GetFastCovering(S2Cap const& cap, std::vector<S2CellId>* covering);

  // Given a connected region and a starting point, return a set of cells at
  // the given level that cover the region.
  //
  // Note that this method is *not* faster than the regular GetCovering()
  // method for most region types, such as S2Cap or S2Polygon, and in fact it
  // can be much slower when the output consists of a large number of cells.
  // Currently it can be faster at generating coverings of long narrow regions
  // such as polylines, but this may change in the future, in which case this
  // method will most likely be removed.
  static void GetSimpleCovering(S2Region const& region, S2Point const& start,
                                int level, std::vector<S2CellId>* output);
};
```

## Spatial Indexing

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

```c++
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
```

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

  ```c++
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

```c++
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

```c++
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
```

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

## Appendix: Alternatives Considered

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
