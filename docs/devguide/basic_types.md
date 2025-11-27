---
title: Basic Types
---



## S1Angle

The `S1Angle` class represents a one-dimensional angle (as opposed to a 2D solid
angle). It has methods for converting angles to or from radians, degrees, and
the E5/E6/E7 representations (i.e. degrees multiplied by 1e5/1e6/1e7 and rounded
to the nearest integer).

```c++
class S1Angle {
 public:
  // These methods construct S1Angle objects from their measure in radians
  // or degrees.
  static constexpr S1Angle Radians(double radians);
  static constexpr S1Angle Degrees(double degrees);
  static constexpr S1Angle E5(int32 e5);
  static constexpr S1Angle E6(int32 e6);
  static constexpr S1Angle E7(int32 e7);

  // The default constructor yields a zero angle.  This is useful for STL
  // containers and class methods with output arguments.
  constexpr S1Angle();

  // Return an angle larger than any finite angle.
  static constexpr S1Angle Infinity();

  // A explicit shorthand for the default constructor.
  static constexpr S1Angle Zero();

  // Return the angle between two points, which is also equal to the distance
  // between these points on the unit sphere.  The points do not need to be
  // normalized.
  S1Angle(const S2Point& x, const S2Point& y);

  // Like the constructor above, but return the angle (i.e., distance)
  // between two S2LatLng points.
  S1Angle(const S2LatLng& x, const S2LatLng& y);

  constexpr double radians() const;
  constexpr double degrees() const;

  int32 e5() const;
  int32 e6() const;
  int32 e7() const;

  // Return the absolute value of an angle.
  S1Angle abs() const;

  // Return the angle normalized to the range (-180, 180] degrees.
  S1Angle Normalized() const;

  // Normalize this angle to the range (-180, 180] degrees.
  void Normalize();
};
```

See `s1angle.h` for additional methods, including comparison and arithmetic
operators. For example, if `x` and `y` are angles, then you can write

    if (sin(0.5 * x) > (x + y) / (x - y)) { ... }

## S2Point

The `S2Point` class represents a point on the unit sphere as a 3D vector.
Usually points are normalized to be unit length, but some methods do not require
this. The `S2Point` class is simply a synonym for the `Vector_3d` class from
`util/math/vector.h`, which defines overloaded operators that make it convenient
to write arithmetic expressions (e.g. `x*p1 + (1-x)*p2`).

Some of its more useful methods include:

```c++
class S2Point /* Vector3_d */ {
 public:
  S2Point(double x, double y, double z);
  double x(), y(), z();                 // Named component accessors
  double& operator[](int i);            // Return component i (0, 1, 2)
  bool operator==(const S2Point& v);    // Equality testing
  bool operator!=(const S2Point& v);    // Inequality testing
  S2Point operator+=(const S2Point& v); // Add another vector
  S2Point operator-=(const S2Point& v); // Subtract another vector
  S2Point operator*=(double k);         // Multiply by a scalar
  S2Point operator/=(double k);         // Divide by a scalar
  S2Point operator+(const S2Point& v);  // Add two vectors
  S2Point operator-(const S2Point& v);  // Subtract two vectors
  S2Point operator*(double k);          // Multiply by a scalar
  S2Point operator/(double k);          // Divide by a scalar
  double DotProd(const S2Point& v);     // Dot product
  S2Point CrossProd(const S2Point& v);  // Cross product
  double Norm2();                       // Squared L2 norm
  double Norm();                        // L2 norm
  S2Point Normalize();                  // Return a normalized **copy**
  double Angle(const S2Point&v);        // Angle between two vectors (radians)
};
```

Note that the `S2Point` type is technically defined as follows:

    typedef Vector3<double> Vector3_d;
    typedef Vector3_d S2Point;

See `util/math/vector.h` for details and additional methods.

## S2LatLng

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
  explicit S2LatLng(const S2Point& p);

  // Returns an S2LatLng for which is_valid() will return false.
  static S2LatLng Invalid();

  // Convenience functions -- shorter than calling S1Angle::Radians(), etc.
  static S2LatLng FromRadians(double lat_radians, double lng_radians);
  static S2LatLng FromDegrees(double lat_degrees, double lng_degrees);
  static S2LatLng FromE5(int32 lat_e5, int32 lng_e5);
  static S2LatLng FromE6(int32 lat_e6, int32 lng_e6);
  static S2LatLng FromE7(int32 lat_e7, int32 lng_e7);

  // Accessor methods.
  S1Angle lat() const;
  S1Angle lng() const;
  const R2Point& coords() const;

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
  S1Angle GetDistance(const S2LatLng& o) const;
};
```

See `s2latlng.h` for additional methods, including comparison and arithmetic
operators.

## S2Region

An `S2Region` represents a two-dimensional region over the unit sphere. It is an
abstract interface with various concrete subtypes, such as discs, rectangles,
polylines, polygons, geometry collections, buffered shapes, etc.

The main purpose of this interface is to allow complex regions to be
approximated as simpler regions. So rather than having a wide variety of virtual
methods that are implemented by all subtypes, the interface is restricted to
methods that are useful for computing approximations. Here they are:

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
  virtual bool Contains(const S2Cell& cell) const = 0;

  // If this method returns false, the region does not intersect the given
  // cell.  Otherwise, either region intersects the cell, or the intersection
  // relationship could not be determined.
  virtual bool MayIntersect(const S2Cell& cell) const = 0;

  // Return true if and only if the given point is contained by the region.
  // The point 'p' is generally required to be unit length, although some
  // subtypes may relax this restriction.
  virtual bool Contains(const S2Point& p) const = 0;

  //////////////////////////////////////////////////////////////////////////
  // Many S2Region subtypes also define the following non-virtual methods.

  // Appends a serialized representation of the region to "encoder".
  void Encode(Encoder* const encoder) const;

  // Decodes an S2Region encoded with Encode().  Note that this method
  // requires that an S2Region object of the appropriate concrete type has
  // already been constructed.  It is not possible to decode regions of
  // unknown type.
  bool Decode(Decoder* const decoder);
};
```

See `s2region.h` for the full interface.

## S2LatLngRect

An `S2LatLngRect` is a type of `S2Region` that represents a rectangle in
latitude-longitude space. It is capable of representing the empty and full
rectangles as well as single points. It has an `AddPoint` method that makes it
easy to construct a bounding rectangle for a set of points, including point sets
that span the 180 degree meridian. Here are its methods:

```c++
class S2LatLngRect final : public S2Region {
 public:
  // Construct a rectangle from minimum and maximum latitudes and longitudes.
  // If lo.lng() > hi.lng(), the rectangle spans the 180 degree longitude
  // line. Both points must be normalized, with lo.lat() <= hi.lat().
  // The rectangle contains all the points p such that 'lo' <= p <= 'hi',
  // where '<=' is defined in the obvious way.
  S2LatLngRect(const S2LatLng& lo, const S2LatLng& hi);

  // Construct a rectangle from latitude and longitude intervals.  The two
  // intervals must either be both empty or both non-empty, and the latitude
  // interval must not extend outside [-90, +90] degrees.
  S2LatLngRect(const R1Interval& lat, const S1Interval& lng);

  // The default constructor creates an empty S2LatLngRect.
  S2LatLngRect();

  // Construct a rectangle of the given size centered around the given point.
  static S2LatLngRect FromCenterSize(const S2LatLng& center,
                                     const S2LatLng& size);

  // Construct a rectangle containing a single (normalized) point.
  static S2LatLngRect FromPoint(const S2LatLng& p);

  // Construct the minimal bounding rectangle containing the two given
  // normalized points.  The points do not need to be ordered.
  static S2LatLngRect FromPointPair(const S2LatLng& p1, const S2LatLng& p2);

  // Accessor methods.
  S1Angle lat_lo() const;
  S1Angle lat_hi() const;
  S1Angle lng_lo() const;
  S1Angle lng_hi() const;
  const R1Interval& lat() const;
  const S1Interval& lng() const;
  R1Interval *mutable_lat();
  S1Interval *mutable_lng();
  S2LatLng lo() const;
  S2LatLng hi() const;

  // The canonical empty and full rectangles.
  static S2LatLngRect Empty();
  static S2LatLngRect Full();

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

  // Return true if the rectangle crosses the 180 degree longitude line.
  bool is_inverted() const;

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

  // Return the true centroid of the rectangle multiplied by its surface area.
  // The true centroid is the "center of mass" of the rectangle viewed as
  // subset of the unit sphere, i.e. it is the point in space about which this
  // curved shape would rotate.)
  //
  // The reason for multiplying the result by the rectangle area is to make it
  // easier to compute the centroid of more complicated shapes.
  S2Point GetCentroid() const;

  // More efficient version of Contains() that accepts a S2LatLng rather than
  // an S2Point.
  bool Contains(const S2LatLng& ll) const;

  // Return true if and only if the given point is contained in the interior
  // of the region (i.e. the region excluding its boundary).
  bool InteriorContains(const S2Point& p) const;

  // More efficient version of InteriorContains() that accepts a S2LatLng
  // rather than an S2Point.
  bool InteriorContains(const S2LatLng& ll) const;

  // Return true if and only if the rectangle contains the given other
  // rectangle.
  bool Contains(const S2LatLngRect& other) const;

  // Return true if and only if the interior of this rectangle contains all
  // points of the given other rectangle (including its boundary).
  bool InteriorContains(const S2LatLngRect& other) const;

  // Return true if this rectangle and the given other rectangle have any
  // points in common.
  bool Intersects(const S2LatLngRect& other) const;

  // Returns true if this rectangle intersects the given cell.  (This is an
  // exact test and may be fairly expensive, see also MayIntersect below.)
  bool Intersects(const S2Cell& cell) const;

  // Return true if and only if the interior of this rectangle intersects
  // any point (including the boundary) of the given other rectangle.
  bool InteriorIntersects(const S2LatLngRect& other) const;

  // Increase the size of the bounding rectangle to include the given point.
  // The rectangle is expanded by the minimum amount possible.
  void AddPoint(const S2Point& p);
  void AddPoint(const S2LatLng& ll);

  // Return a rectangle that has been expanded by margin.lat() on each side in
  // the latitude direction, and by margin.lng() on each side in the longitude
  // direction.  If either margin is negative, then shrink the rectangle on
  // the corresponding sides instead.  The resulting rectangle may be empty.
  //
  // If you are trying to grow a rectangle by a certain *distance* on the
  // sphere (e.g. 5km), use the ExpandedByDistance() method instead.
  S2LatLngRect Expanded(const S2LatLng& margin) const;

  // If the rectangle does not include either pole, return it unmodified.
  // Otherwise expand the longitude range to Full() so that the rectangle
  // contains all possible representations of the contained pole(s).
  S2LatLngRect PolarClosure() const;

  // Return the smallest rectangle containing the union of this rectangle and
  // the given rectangle.
  S2LatLngRect Union(const S2LatLngRect& other) const;

  // Return the smallest rectangle containing the intersection of this
  // rectangle and the given rectangle.  Note that the region of intersection
  // may consist of two disjoint rectangles, in which case a single rectangle
  // spanning both of them is returned.
  S2LatLngRect Intersection(const S2LatLngRect& other) const;

  // Expand this rectangle so that it contains all points within the given
  // distance of the boundary, and return the smallest such rectangle.  If the
  // distance is negative, then instead shrink this rectangle so that it
  // excludes all points within the given absolute distance of the boundary,
  // and return the largest such rectangle.
  //
  // Unlike Expanded(), this method treats the rectangle as a set of points on
  // the sphere, and measures distances on the sphere.  For example, you can
  // use this method to find a rectangle that contains all points within 5km
  // of a given rectangle.
  S2LatLngRect ExpandedByDistance(S1Angle distance) const;

  // Returns the minimum distance (measured along the surface of the sphere) to
  // the given S2LatLngRect. Both S2LatLngRects must be non-empty.
  S1Angle GetDistance(const S2LatLngRect& other) const;

  // Returns the minimum distance (measured along the surface of the sphere)
  // from a given point to the rectangle (both its boundary and its interior).
  S1Angle GetDistance(const S2LatLng& p) const;

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  S2LatLngRect* Clone() const override;
  S2Cap GetCapBound() const override;
  S2LatLngRect GetRectBound() const override;
  // The point 'p' does not need to be normalized.
  bool Contains(const S2Cell& cell) const override;

  // This test is cheap but is NOT exact.  Use Intersects() if you want a more
  // accurate and more expensive test.
  bool MayIntersect(const S2Cell& cell) const override;
};
```

See `s2latlngrect.h` for additional methods.

## S1ChordAngle

The S2 library actually defines two angle representations, `S1Angle` and
`S1ChordAngle`.

`S1ChordAngle` is a specialized angle type that represents the distance
between two points on the sphere.  (Note that the distance between two
points can be expressed as the angle between those points measured from the
sphere's center.)  Its representation makes it very efficient for computing
and comparing distances, but unlike `S1Angle` it is only capable of
representing angles between 0 and 180 degrees.  (Internally, `S1ChordAngle`
computes the squared distance between two points through the interior of
the sphere, i.e. the squared length of the chord between those points).

`S1ChordAngle` is the preferred representation for distances internally,
because it is much faster to compute than `S1Angle` and also supports the
exact predicates needed for robust geometric algorithms.  `S1ChordAngle`
supports many of the same methods as `S1Angle`, and it is also easy to
convert back and forth as required.

## S2Cap

An `S2Cap` represents a spherical cap, i.e. a portion of a sphere cut off by a
plane. This is the equivalent of a *closed disc* in planar geometry (i.e., a
circle together with its interior), so if you are looking for a way to represent
a circle or disc then you are probably looking for an `S2Cap`.

Here are the methods available:

```c++
class S2Cap final : public S2Region {
 public:
  // The default constructor returns an empty S2Cap.
  S2Cap();

  // Construct a cap with the given center and radius.  A negative radius
  // yields an empty cap; a radius of 180 degrees or more yields a full cap
  // (containing the entire sphere).
  S2Cap(const S2Point& center, S1Angle radius);

  // Convenience function that creates a cap containing a single point.  This
  // method is more efficient that the S2Cap(center, radius) constructor.
  static S2Cap FromPoint(const S2Point& center);

  // Return a cap with the given center and height.  (The height of a cap is
  // the distance from the cap center to the cutoff plane; see comments
  // above.)  A negative height yields an empty cap; a height of 2 or more
  // yields a full cap.
  static S2Cap FromCenterHeight(const S2Point& center, double height);

  // Return a cap with the given center and surface area.  Note that the area
  // can also be interpreted as the solid angle subtended by the cap (because
  // the sphere has unit radius).  A negative area yields an empty cap; an
  // area of 4*Pi or more yields a full cap.
  static S2Cap FromCenterArea(const S2Point& center, double area);

  // Return an empty cap, i.e. a cap that contains no points.
  static S2Cap Empty();

  // Return a full cap, i.e. a cap that contains all points.
  static S2Cap Full();

  // Accessor methods.
  const S2Point& center() const;
  S1ChordAngle radius() const;
  double height() const;

  // Return the cap radius.  (This method is relatively expensive since only
  // the cap height is stored internally, and may yield a slightly different
  // result than the value passed to the (center, radius) constructor.)
  S1Angle GetRadius() const;

  // Return the area of the cap.
  double GetArea() const;

  // Return the true centroid of the cap multiplied by its surface area.  Note
  // that if you just want the "surface centroid", i.e. the normalized result,
  // then it is much simpler just to call center().
  //
  // The reason for multiplying the result by the cap area is to make it
  // easier to compute the centroid of more complicated shapes.
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
  S2Cap Complement() const;

  // Return true if and only if this cap contains the given other cap
  // (in a set containment sense, e.g. every cap contains the empty cap).
  bool Contains(const S2Cap& other) const;

  // Return true if and only if this cap intersects the given other cap,
  // i.e. whether they have any points in common.
  bool Intersects(const S2Cap& other) const;

  // Return true if and only if the interior of this cap intersects the
  // given other cap.  (This relationship is not symmetric, since only
  // the interior of this cap is used.)
  bool InteriorIntersects(const S2Cap& other) const;

  // Return true if and only if the given point is contained in the interior
  // of the cap (i.e. the cap excluding its boundary).
  bool InteriorContains(const S2Point& p) const;

  // Increase the cap height if necessary to include the given point.  If the
  // cap is empty then the center is set to the given point, but otherwise the
  // center is not changed.
  void AddPoint(const S2Point& p);

  // Increase the cap height if necessary to include "other".  If the current
  // cap is empty it is set to the given other cap.
  void AddCap(const S2Cap& other);

  // Return a cap that contains all points within a given distance of this
  // cap.  Note that any expansion of the empty cap is still empty.
  S2Cap Expanded(S1Angle distance) const;

  // Return the smallest cap which encloses this cap and "other".
  S2Cap Union(const S2Cap& other) const;

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  S2Cap* Clone() const override;
  S2Cap GetCapBound() const override;
  S2LatLngRect GetRectBound() const override;
  // The point "p" should be a unit-length vector.
  bool Contains(const S2Cell& cell) const override;
  bool MayIntersect(const S2Cell& cell) const override;
};
```

See `s2cap.h` for additional methods.

## S2Polyline

An `S2Polyline` represents a sequence of zero or more vertices connected by
straight edges (geodesics). Edges of length 0 and 180 degrees are not allowed,
i.e. adjacent vertices should not be identical or antipodal.

```c++
class S2Polyline final : public S2Region {
 public:
  // Creates an empty S2Polyline that should be initialized by calling Init()
  // or Decode().
  S2Polyline();

  // Convenience constructors that call Init() with the given vertices.
  explicit S2Polyline(const std::vector<S2Point>& vertices);
  explicit S2Polyline(const std::vector<S2LatLng>& vertices);

  ~S2Polyline();

  // Initialize a polyline that connects the given vertices. Empty polylines are
  // allowed.  Adjacent vertices should not be identical or antipodal.  All
  // vertices should be unit length.
  void Init(const std::vector<S2Point>& vertices);

  // Convenience initialization function that accepts latitude-longitude
  // coordinates rather than S2Points.
  void Init(const std::vector<S2LatLng>& vertices);

  // Return true if the given vertices form a valid polyline.
  bool IsValid() const;

  // Returns true if this is *not* a valid polyline and sets "error"
  // appropriately.  Otherwise returns false and leaves "error" unchanged.
  bool FindValidationError(S2Error* error) const;

  // Accessor methods.
  int num_vertices() const;
  const S2Point& vertex(int k) const;

  // Return the length of the polyline.
  S1Angle GetLength() const;

  // Return the true centroid of the polyline multiplied by the length of the
  // polyline.  Prescaling by the polyline length makes it easy to compute the
  // centroid of several polylines (by simply adding up their centroids).
  S2Point GetCentroid() const;

  // Return the point whose distance from vertex 0 along the polyline is the
  // given fraction of the polyline's total length.  Fractions less than zero
  // or greater than one are clamped.  The return value is unit length.
  S2Point Interpolate(double fraction) const;

  // Like Interpolate(), but also return the index of the next polyline
  // vertex after the interpolated point P.  This allows the caller to easily
  // construct a given suffix of the polyline by concatenating P with the
  // polyline vertices starting at "next_vertex".
  S2Point GetSuffix(double fraction, int* next_vertex) const;

  // The inverse operation of GetSuffix/Interpolate.  Given a point on the
  // polyline, returns the ratio of the distance to the point from the
  // beginning of the polyline over the length of the polyline.  The return
  // value is always betwen 0 and 1 inclusive.  See GetSuffix() for the
  // meaning of "next_vertex".
  double UnInterpolate(const S2Point& point, int next_vertex) const;

  // Given a point, returns a point on the polyline that is closest to the given
  // point.  See GetSuffix() for the meaning of "next_vertex".
  S2Point Project(const S2Point& point, int* next_vertex) const;

  // Returns true if the point given is on the right hand side of the polyline,
  // using a naive definition of "right-hand-sideness" where the point is on
  // the RHS of the polyline iff the point is on the RHS of the line segment in
  // the polyline which it is closest to.
  bool IsOnRight(const S2Point& point) const;

  // Return true if this polyline intersects the given polyline. If the
  // polylines share a vertex they are considered to be intersecting.
  bool Intersects(const S2Polyline* line) const;

  // Reverse the order of the polyline vertices.
  void Reverse();

  // Return a subsequence of vertex indices such that the polyline connecting
  // these vertices is never further than "tolerance" from the original
  // polyline.
  void SubsampleVertices(S1Angle tolerance, std::vector<int>* indices) const;

  // Return true if "covered" is within "max_error" of a contiguous subpath of
  // this polyline over its entire length.  Specifically, this method returns
  // true if this polyline has parameterization a:[0,1] -> S^2, "covered" has
  // parameterization b:[0,1] -> S^2, and there is a non-decreasing function
  // f:[0,1] -> [0,1] such that distance(a(f(t)), b(t)) <= max_error for all t.
  //
  // You can think of this as testing whether it is possible to drive a car
  // along "covered" and a car along some subpath of this polyline such that no
  // car ever goes backward, and the cars are always within "max_error" of each
  // other.
  bool NearlyCoversPolyline(const S2Polyline& covered,
                            S1Angle max_error) const;

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  S2Polyline* Clone() const override;
  S2Cap GetCapBound() const override;
  S2LatLngRect GetRectBound() const override;
  bool Contains(const S2Cell& cell) const override;
  bool MayIntersect(const S2Cell& cell) const override;
};
```

## S2Loop

An `S2Loop` represents a simple spherical polygon. It consists of a single chain
of vertices where the first vertex is implicitly connected to the last. All
loops are defined to have a CCW orientation, i.e. the interior of the loop is on
the left side of the edges. This implies that a clockwise loop enclosing a small
area is interpreted to be a CCW loop enclosing a very large area.

Loops are not allowed to have any duplicate vertices (whether adjacent or not),
and non-adjacent edges are not allowed to intersect. Loops must have at least 3
vertices (except for the "empty" and "full" loops discussed below). These
restrictions make it possible to implement exact polygon-polygon containment and
intersection tests very efficiently. See `S2Builder` if your data does not meet
these requirements.

There are two special loops: the "empty" loop contains no points, while the
"full" loop contains all points. These loops do not have any edges, but to
preserve the invariant that every loop can be represented as a vertex chain,
they are defined as having exactly one vertex each (see kEmpty and kFull).

Point containment of loops is defined such that if the sphere is subdivided into
faces (loops), every point is contained by exactly one face. This implies that
loops do not necessarily contain their vertices.

```c++
class S2Loop final : public S2Region {
 public:
  // Default constructor.  The loop must be initialized by calling Init() or
  // Decode() before it is used.
  S2Loop();

  // Convenience constructor that calls Init() with the given vertices.
  explicit S2Loop(const std::vector<S2Point>& vertices);

  // Construct a loop corresponding to the given cell.
  //
  // Note that the loop and cell *do not* contain exactly the same set of
  // points, because S2Loop and S2Cell have slightly different definitions of
  // point containment.  For example, an S2Cell vertex is contained by all
  // four neighboring S2Cells, but it is contained by exactly one of four
  // S2Loops constructed from those cells.
  explicit S2Loop(const S2Cell& cell);

  ~S2Loop();

  // Initialize a loop with given vertices.  The last vertex is implicitly
  // connected to the first.  All points should be unit length.  Loops must
  // have at least 3 vertices (except for the "empty" and "full" loops, see
  // kEmpty and kFull).  This method may be called multiple times.
  void Init(const std::vector<S2Point>& vertices);

  // A special vertex chain of length 1 that creates an empty loop (i.e., a
  // loop with no edges that contains no points).  Example usage:
  //    S2Loop empty(S2Loop::kEmpty());
  static std::vector<S2Point> kEmpty();

  // A special vertex chain of length 1 that creates a full loop (i.e., a loop
  // with no edges that contains all points).  See kEmpty() for details.
  static std::vector<S2Point> kFull();

  // Returns true if this is a valid loop.
  bool IsValid() const;

  // Returns true if this is *not* a valid loop and sets "error"
  // appropriately.  Otherwise returns false and leaves "error" unchanged.
  bool FindValidationError(S2Error* error) const;

  int num_vertices() const;

  // For convenience, we make two entire copies of the vertex list available:
  // vertex(n..2*n-1) is mapped to vertex(0..n-1), where n == num_vertices().
  const S2Point& vertex(int i) const;

  // Like vertex(), but this method returns vertices in reverse order if the
  // loop represents a polygon hole. This ensures that the interior of the
  // polygon is always to the left of the vertex chain.
  const S2Point& oriented_vertex(int i) const;

  // Return true if this is the special "empty" loop that contains no points.
  bool is_empty() const;

  // Return true if this is the special "full" loop that contains all points.
  bool is_full() const;

  // Return true if this loop is either "empty" or "full".
  bool is_empty_or_full() const;

  // The depth of a loop is defined as its nesting level within its containing
  // polygon.  "Outer shell" loops have depth 0, holes within those loops have
  // depth 1, shells within those holes have depth 2, etc.  This field is only
  // used by the S2Polygon implementation.
  int depth() const;
  void set_depth(int depth);

  // Return true if this loop represents a hole in its containing polygon.
  bool is_hole() const;

  // The sign of a loop is -1 if the loop represents a hole in its containing
  // polygon, and +1 otherwise.
  int sign() const;

  // Return true if the loop area is at most 2*Pi.
  bool IsNormalized() const;

  // Invert the loop if necessary so that the area enclosed by the loop is at
  // most 2*Pi.
  void Normalize();

  // Reverse the order of the loop vertices, effectively complementing
  // the region represented by the loop.
  void Invert();

  // Return the area of the loop interior, i.e. the region on the left side of
  // the loop.  The return value is between 0 and 4*Pi.
  double GetArea() const;

  // Return the true centroid of the loop multiplied by the area of the loop.
  // We prescale by the loop area for two reasons: (1) it is cheaper to
  // compute this way, and (2) it makes it easier to compute the centroid of
  // more complicated shapes (by splitting them into disjoint regions and
  // adding their centroids).
  S2Point GetCentroid() const;

  // Return the sum of the turning angles at each vertex.  The return value is
  // positive if the loop is counter-clockwise, negative if the loop is
  // clockwise, and zero if the loop is a great circle.
  //
  // This quantity is also called the "geodesic curvature" of the loop.
  double GetCurvature() const;

  // Return the maximum error in GetCurvature().  The return value is not
  // constant; it depends on the loop.
  double GetCurvatureMaxError() const;

  // Return the distance from the given point to the loop interior.  If the
  // loop is empty, return S1Angle::Infinity().
  S1Angle GetDistance(const S2Point& x) const;

  // Return the distance from the given point to the loop boundary.  If the
  // loop is empty or full, return S1Angle::Infinity() (since the loop has no
  // boundary).
  S1Angle GetDistanceToBoundary(const S2Point& x) const;

  // If the given point is contained by the loop, return it.  Otherwise return
  // the closest point on the loop boundary.  If the loop is empty, return the
  // input argument.
  S2Point Project(const S2Point& x) const;

  // Return the closest point on the loop boundary to the given point.  If the
  // loop is empty or full, return the input argument (since the loop has no
  // boundary).
  S2Point ProjectToBoundary(const S2Point& x) const;

  // Return true if the region contained by this loop is a superset of the
  // region contained by the given other loop.
  bool Contains(const S2Loop* b) const;

  // Return true if the region contained by this loop intersects the region
  // contained by the given other loop.
  bool Intersects(const S2Loop* b) const;

  // Return true if two loops have the same vertices in the same linear order
  // (i.e., cyclic rotations are not allowed).
  bool Equals(const S2Loop* b) const;

  // Return true if two loops have the same boundary.  This is true if and
  // only if the loops have the same vertices in the same cyclic order (i.e.,
  // the vertices may be cyclically rotated).  The empty and full loops are
  // considered to have different boundaries.
  bool BoundaryEquals(const S2Loop* b) const;

  // Return true if two loops have the same boundary except for vertex
  // perturbations.  More precisely, the vertices in the two loops must be in
  // the same cyclic order, and corresponding vertex pairs must be separated
  // by no more than "max_error".
  bool BoundaryApproxEquals(const S2Loop& b,
                            S1Angle max_error = S1Angle::Radians(1e-15)) const;

  // Return true if the two loop boundaries are within "max_error" of each
  // other along their entire lengths.  The two loops may have different
  // numbers of vertices.  More precisely, this method returns true if the two
  // loops have parameterizations a:[0,1] -> S^2, b:[0,1] -> S^2 such that
  // distance(a(t), b(t)) <= max_error for all t.  You can think of this as
  // testing whether it is possible to drive two cars all the way around the
  // two loops such that no car ever goes backward and the cars are always
  // within "max_error" of each other.
  bool BoundaryNear(const S2Loop& b,
                    S1Angle max_error = S1Angle::Radians(1e-15)) const;

  // This method computes the oriented surface integral of some quantity f(x)
  // over the loop interior, given a function f_tri(A,B,C) that returns the
  // corresponding integral over the spherical triangle ABC.
  template <class T>
  T GetSurfaceIntegral(T f_tri(const S2Point&, const S2Point&, const S2Point&))
      const;

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  S2Loop* Clone() const override;

  // GetRectBound() returns essentially tight results, while GetCapBound()
  // might have a lot of extra padding.  Both bounds are conservative in that
  // if the loop contains a point P, then the bound contains P also.
  S2Cap GetCapBound() const override;
  S2LatLngRect GetRectBound() const override;
  bool Contains(const S2Cell& cell) const override;
  bool MayIntersect(const S2Cell& cell) const override;
  // Return true if the loop contains the given point.  Point containment is
  // defined such that if the sphere is subdivided into faces (loops), every
  // point is contained by exactly one face.  This implies that loops do not
  // necessarily contain their vertices.
  bool Contains(const S2Point& p) const override;

  // Generally clients should not use S2Loop::Encode().  Instead they should
  // encode an S2Polygon, which unlike this method supports (lossless)
  // compression.
  void Encode(Encoder* const encoder) const override;

  // Decode a loop encoded with Encode() or EncodeCompressed().  These methods
  // may be called with loops that have already been initialized.
  bool Decode(Decoder* const decoder) override;
  bool DecodeWithinScope(Decoder* const decoder) override;
};
```

## S2Polygon

An `S2Polygon` is an `S2Region` object that represents a polygon. A polygon is
defined by zero or more loops; recall that the interior of a loop is defined to
be its left-hand side (see `S2Loop`). There are two different conventions for
creating an `S2Polygon`:

*   `InitNested()` expects the input loops to be nested hierarchically. The
    polygon interior then consists of the set of points contained by an odd
    number of loops. So for example, a circular region with a hole in it would
    be defined as two CCW loops, with one loop containing the other. The loops
    can be provided in any order.

    When the orientation of the input loops is unknown, the nesting requirement
    is typically met by calling `S2Loop::Normalize()` on each loop (which
    inverts the loop if necessary so that it encloses at most half the sphere).
    But in fact any set of loops can be used as long as (1) there is no pair of
    loops that cross, and (2) there is no pair of loops whose union is the
    entire sphere.

*   `InitOriented()` expects the input loops to be oriented such that the
    polygon interior is on the left-hand side of every loop. So for example, a
    circular region with a hole in it would be defined using a CCW outer loop
    and a CW inner loop. The loop orientations must all be consistent; for
    example, it is not valid to have one CCW loop nested inside another CCW
    loop, because the region between the two loops is on the left-hand side of
    one loop and the right-hand side of the other.

Most clients will not call these methods directly; instead they should use
`S2Builder`, which has better support for dealing with imperfect data.

When the polygon is initialized, the given loops are automatically converted
into a canonical form consisting of "shells" and "holes". Shells and holes are
both oriented CCW, and are nested hierarchically. The loops are reordered to
correspond to a preorder traversal of the nesting hierarchy; `InitOriented()`
may also invert some loops.

Polygons may represent any region of the sphere with a polygonal boundary,
including the entire sphere (known as the "full" polygon). The full polygon
consists of a single full loop (see `S2Loop`), whereas the empty polygon has no
loops at all.

Polygons have the following restrictions:

*   Loops may not cross, i.e. the boundary of a loop may not intersect both the
    interior and exterior of any other loop.

*   Loops may not share edges, i.e. if a loop contains an edge AB, then no other
    loop may contain AB or BA.

*   Loops may share vertices, however no vertex may appear twice in a single
    loop (see `S2Loop`).

*   No loop may be empty. The full loop may appear only in the full polygon.


```c++
class S2Polygon final : public S2Region {
 public:
  // The default constructor creates an empty polygon.  It can be made
  // non-empty by calling Init(), Decode(), etc.
  S2Polygon();

  // Convenience constructor that calls InitNested() with the given loops.
  // Takes ownership of the loops and clears the vector.
  explicit S2Polygon(std::vector<std::unique_ptr<S2Loop>> loops);

  // Convenience constructor that creates a polygon with a single loop
  // corresponding to the given cell.
  explicit S2Polygon(const S2Cell& cell);

  // Convenience constructor that calls Init(unique_ptr<S2Loop>).
  explicit S2Polygon(std::unique_ptr<S2Loop> loop);

  // Destroys the polygon and frees its loops.
  ~S2Polygon();

  // Create a polygon from a set of hierarchically nested loops.  The polygon
  // interior consists of the points contained by an odd number of loops.
  // (Recall that a loop contains the set of points on its left-hand side.)
  //
  // This method takes ownership of the given loops and clears the given
  // vector.  It then figures out the loop nesting hierarchy and assigns every
  // loop a depth.  Shells have even depths, and holes have odd depths.  Note
  // that the loops are reordered so the hierarchy can be traversed more
  // easily (see GetParent(), GetLastDescendant(), and S2Loop::depth()).
  //
  // This method may be called more than once, in which case any existing
  // loops are deleted before being replaced by the input loops.
  void InitNested(std::vector<std::unique_ptr<S2Loop>> loops);

  // Like InitNested(), but expects loops to be oriented such that the polygon
  // interior is on the left-hand side of all loops.  This implies that shells
  // and holes should have opposite orientations in the input to this method.
  void InitOriented(std::vector<std::unique_ptr<S2Loop>> loops);

  // Initialize a polygon from a single loop.  Takes ownership of the loop.
  void Init(std::unique_ptr<S2Loop> loop);

  // Releases ownership of and returns the loops of this polygon, and resets
  // the polygon to be empty.
  std::vector<std::unique_ptr<S2Loop>> Release();

  // Makes a deep copy of the given source polygon.  The destination polygon
  // will be cleared if necessary.
  void Copy(const S2Polygon* src);

  // Returns true if this is a valid polygon (including checking whether all
  // the loops are themselves valid).
  bool IsValid() const;

  // Returns true if this is *not* a valid polygon and sets "error"
  // appropriately.  Otherwise returns false and leaves "error" unchanged.
  bool FindValidationError(S2Error* error) const;

  // Return true if this is the empty polygon (consisting of no loops).
  bool is_empty() const;

  // Return true if this is the full polygon (consisting of a single loop that
  // encompasses the entire sphere).
  bool is_full() const;

  // Return the number of loops in this polygon.
  int num_loops() const;

  // Total number of vertices in all loops.
  int num_vertices() const;

  // Return the loop at the given index.  Note that during initialization, the
  // given loops are reordered according to a preorder traversal of the loop
  // nesting hierarchy.  This implies that every loop is immediately followed
  // by its descendants.  This hierarchy can be traversed using the methods
  // GetParent(), GetLastDescendant(), and S2Loop::depth().
  const S2Loop* loop(int k) const;
  S2Loop* loop(int k);

  // Return the index of the parent of loop k, or -1 if it has no parent.
  int GetParent(int k) const;

  // Return the index of the last loop that is contained within loop k.
  int GetLastDescendant(int k) const;

  // Return the area of the polygon interior, i.e. the region on the left side
  // of an odd number of loops.  The return value is between 0 and 4*Pi.
  double GetArea() const;

  // Return the true centroid of the polygon multiplied by the area of the
  // polygon.  We prescale by the polygon area for two reasons: (1) it is
  // cheaper to compute this way, and (2) it makes it easier to compute the
  // centroid of more complicated shapes (by splitting them into disjoint
  // regions and adding their centroids).
  S2Point GetCentroid() const;

  // If all of the polygon's vertices happen to be the centers of S2Cells at
  // some level, then return that level, otherwise return -1.
  int GetSnapLevel() const;

  // Return the distance from the given point to the polygon interior.  If the
  // polygon is empty, return S1Angle::Infinity().
  S1Angle GetDistance(const S2Point& x) const;

  // Return the distance from the given point to the polygon boundary.  If the
  // polygon is empty or full, return S1Angle::Infinity() (since the polygon
  // has no boundary).
  S1Angle GetDistanceToBoundary(const S2Point& x) const;

  // Return the overlap fractions between two polygons, i.e. the ratios of the
  // area of intersection to the area of each polygon.
  static std::pair<double, double> GetOverlapFractions(const S2Polygon* a,
                                                       const S2Polygon* b);

  // If the given point is contained by the polygon, return it.  Otherwise
  // return the closest point on the polygon boundary.  If the polygon is
  // empty, return the input argument.
  S2Point Project(const S2Point& x) const;

  // Return the closest point on the polygon boundary to the given point.  If
  // the polygon is empty or full, return the input argument (since the
  // polygon has no boundary).
  S2Point ProjectToBoundary(const S2Point& x) const;

  // Return true if this polygon contains the given other polygon, i.e.
  // if polygon A contains all points contained by polygon B.
  bool Contains(const S2Polygon* b) const;

  // Returns true if this polgyon (A) approximately contains the given other
  // polygon (B). This is true if it is possible to move the vertices of B
  // no further than "tolerance" such that A contains the modified B.
  //
  // For example, the empty polygon will contain any polygon whose maximum
  // width is no more than "tolerance".
  bool ApproxContains(const S2Polygon* b, S1Angle tolerance) const;

  // Return true if this polygon intersects the given other polygon, i.e.
  // if there is a point that is contained by both polygons.
  bool Intersects(const S2Polygon* b) const;

  // Returns true if this polgyon (A) and the given polygon (B) are
  // approximately disjoint.  This is true if it is possible to ensure that A
  // and B do not intersect by moving their vertices no further than
  // "tolerance".
  //
  // For example, any polygon is approximately disjoint from a polygon whose
  // maximum width is no more than "tolerance".
  bool ApproxDisjoint(const S2Polygon& b, S1Angle tolerance) const;

  // Initialize this polygon to the intersection, union, difference (A - B),
  // or symmetric difference (XOR) of the given two polygons.
  //
  // "snap_function" allows you to specify a minimum spacing between output
  // vertices, and/or that the vertices should be snapped to a discrete set of
  // points (e.g. S2CellId centers or E7 lat/lng coordinates).  Any snap
  // function can be used, including the IdentitySnapFunction with a
  // snap_radius of zero (which preserves the input vertices exactly).
  void InitToIntersection(const S2Polygon* a, const S2Polygon* b);
  void InitToIntersection(const S2Polygon& a, const S2Polygon& b,
                          const S2Builder::SnapFunction& snap_function);

  void InitToUnion(const S2Polygon* a, const S2Polygon* b);
  void InitToUnion(const S2Polygon& a, const S2Polygon& b,
                   const S2Builder::SnapFunction& snap_function);

  void InitToDifference(const S2Polygon* a, const S2Polygon* b);
  void InitToDifference(const S2Polygon& a, const S2Polygon& b,
                        const S2Builder::SnapFunction& snap_function);

  void InitToSymmetricDifference(const S2Polygon* a, const S2Polygon* b);
  void InitToSymmetricDifference(
      const S2Polygon& a, const S2Polygon& b,
      const S2Builder::SnapFunction& snap_function);

  // Convenience functions that use the IdentitySnapFunction with the given
  // snap radius.
  void InitToApproxIntersection(const S2Polygon* a, const S2Polygon* b,
                                S1Angle vertex_merge_radius);
  void InitToApproxUnion(const S2Polygon* a, const S2Polygon* b,
                         S1Angle vertex_merge_radius);
  void InitToApproxDifference(const S2Polygon* a, const S2Polygon* b,
                              S1Angle vertex_merge_radius);
  void InitToApproxSymmetricDifference(const S2Polygon* a, const S2Polygon* b,
                                       S1Angle vertex_merge_radius);

  // Initializes this polygon according to the given "snap_function" and
  // reduces the number of vertices if possible, while ensuring that no vertex
  // moves further than snap_function.snap_radius().
  //
  // Simplification works by replacing nearly straight chains of short edges
  // with longer edges, in a way that preserves the topology of the input
  // polygon up to the creation of degeneracies.  This means that loops or
  // portions of loops may become degenerate, in which case they are removed.
  void InitToSimplified(const S2Polygon& a,
                        const S2Builder::SnapFunction& snap_function);

  // Use S2Builder to build this polygon by assembling the edges of a
  // given polygon after snapping its vertices to the center of leaf cells at
  // the given "snap_level".  The default snap level corresponds to a
  // tolerance of approximately 1.5cm on the surface of the Earth.
  // Such a polygon can be efficiently compressed when serialized.
  void InitToSnapped(const S2Polygon* polygon,
                     int snap_level = S2CellId::kMaxLevel);

  // Initialize this polygon to the complement of the given polygon.
  void InitToComplement(const S2Polygon* a);

  // Invert the polygon (replace it by its complement).
  void Invert();

  // Intersect this polygon with the polyline "in" and append the resulting
  // zero or more polylines to "out" (which must be empty).  The polylines
  // are appended in the order they would be encountered by traversing "in"
  // from beginning to end.  Note that the output may include polylines with
  // only one vertex, but there will not be any zero-vertex polylines.
  std::vector<std::unique_ptr<S2Polyline>>
  IntersectWithPolyline(const S2Polyline& in) const;

  // Similar to IntersectWithPolyline(), except that vertices will be
  // dropped as necessary to ensure that all adjacent vertices in the
  // sequence obtained by concatenating the output polylines will be
  // farther than "vertex_merge_radius" apart.
  std::vector<std::unique_ptr<S2Polyline>>
  ApproxIntersectWithPolyline(const S2Polyline& in,
                              S1Angle vertex_merge_radius) const;

  // Like IntersectWithPolyline, but subtracts this polygon from
  // the given polyline.
  std::vector<std::unique_ptr<S2Polyline>>
  SubtractFromPolyline(const S2Polyline& in) const;

  // Like ApproxIntersectWithPolyline, but subtracts this polygon
  // from the given polyline.
  std::vector<std::unique_ptr<S2Polyline>>
  ApproxSubtractFromPolyline(const S2Polyline& in,
                             S1Angle vertex_merge_radius) const;

  // Return a polygon which is the union of the given polygons.
  // Clears the vector and deletes the polygons.
  static std::unique_ptr<S2Polygon>
  DestructiveUnion(std::vector<std::unique_ptr<S2Polygon>> polygons);
  static std::unique_ptr<S2Polygon>
  DestructiveApproxUnion(std::vector<std::unique_ptr<S2Polygon>> polygons,
                         S1Angle vertex_merge_radius);

  // Initialize this polygon to the outline of the given cell union.
  // In principle this polygon should exactly contain the cell union and
  // this polygon's inverse should not intersect the cell union, but rounding
  // issues may cause this not to be the case.
  void InitToCellUnionBorder(const S2CellUnion& cells);

  // Return true if every loop of this polygon shares at most one vertex with
  // its parent loop.  Every polygon has a unique normalized form.  Normalized
  // polygons are useful for testing since it is easy to compare whether two
  // polygons represent the same region.
  bool IsNormalized() const;

  // Return true if two polygons have exactly the same loops.  The loops must
  // appear in the same order, and corresponding loops must have the same
  // linear vertex ordering (i.e., cyclic rotations are not allowed).
  bool Equals(const S2Polygon* b) const;

  // Return true if two polygons have the same boundary.  More precisely, this
  // method requires that both polygons have loops with the same cyclic vertex
  // order and the same nesting hierarchy.  (This implies that vertices may be
  // cyclically rotated between corresponding loops, and the loop ordering may
  // be different between the two polygons as long as the nesting hierarchy is
  // the same.)
  bool BoundaryEquals(const S2Polygon* b) const;

  // Return true if two polygons have the same boundary except for vertex
  // perturbations.  Both polygons must have loops with the same cyclic vertex
  // order and the same nesting hierarchy, but the vertex locations are
  // allowed to differ by up to "max_error".
  bool BoundaryApproxEquals(const S2Polygon& b,
                            S1Angle max_error = S1Angle::Radians(1e-15)) const;

  // Return true if two polygons have boundaries that are within "max_error"
  // of each other along their entire lengths.  More precisely, there must be
  // a bijection between the two sets of loops such that for each pair of
  // loops, "a_loop->BoundaryNear(b_loop)" is true.
  bool BoundaryNear(const S2Polygon& b,
                    S1Angle max_error = S1Angle::Radians(1e-15)) const;

  ////////////////////////////////////////////////////////////////////////
  // S2Region interface (see s2region.h for details):

  // GetRectBound() returns essentially tight results, while GetCapBound()
  // might have a lot of extra padding.  Both bounds are conservative in that
  // if the loop contains a point P, then the bound contains P also.
  S2Polygon* Clone() const override;
  S2Cap GetCapBound() const override;
  S2LatLngRect GetRectBound() const override;

  bool Contains(const S2Cell& cell) const override;
  bool MayIntersect(const S2Cell& cell) const override;
  // Return true if the polygon contains the given point.  Point containment
  // is defined such that if the sphere is subdivided into polygons, every
  // point is contained by exactly one polygons.  This implies that polygons
  // do not necessarily contain their vertices.
  bool Contains(const S2Point& p) const override;

  // Encode the polygon with about 4 bytes per vertex, assuming the vertices
  // have all been snapped to the centers of S2Cells at a given level
  // (typically with InitToSnapped).  The other vertices are stored using 24
  // bytes.  Decoding a polygon encoded this way always returns the original
  // polygon, without any loss of precision.
  void Encode(Encoder* const encoder) const override;

  // Decode a polygon encoded with Encode().
  bool Decode(Decoder* const decoder) override;
};
```
