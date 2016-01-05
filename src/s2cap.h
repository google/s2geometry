// Copyright 2005 Google Inc. All Rights Reserved.
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

#ifndef S2GEOMETRY_S2CAP_H_
#define S2GEOMETRY_S2CAP_H_

#include <cmath>
#include <algorithm>
#include <iosfwd>

#include <glog/logging.h>
#include "fpcontractoff.h"
#include "s1angle.h"
#include "s2.h"
#include "s2region.h"

class Decoder;
class Encoder;
class S2Cell;
class S2LatLngRect;

// S2Cap represents a disc-shaped region defined by a center and radius.
// Technically this shape is called a "spherical cap" (rather than disc)
// because it is not planar; the cap represents a portion of the sphere that
// has been cut off by a plane.  The boundary of the cap is the circle defined
// by the intersection of the sphere and the plane.  For containment purposes,
// the cap is a closed set, i.e. it contains its boundary.
//
// For the most part, you can use a spherical cap wherever you would use a
// disc in planar geometry.  The radius of the cap is measured along the
// surface of the sphere (rather than the straight-line distance through the
// interior).  Thus a cap of radius Pi/2 is a hemisphere, and a cap of radius
// Pi covers the entire sphere.
//
// Internally, the cap is represented by its center and "height".  The height
// is simply the distance from the center point to the cutoff plane.  This
// representation is much more efficient for containment tests than the
// (center, radius) representation.  There is also support for "empty" and
// "full" caps, which contain no points and all points respectively.
//
// This class is intended to be copied by value as desired.  It uses
// the default copy constructor and assignment operator, however it is
// not a "plain old datatype" (POD) because it has virtual functions.
class S2Cap : public S2Region {
 public:
  // The default constructor returns an empty S2Cap.
  S2Cap() : center_(1, 0, 0), height_(-1) {}

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
  static S2Cap Empty() { return S2Cap(); }

  // Return a full cap, i.e. a cap that contains all points.
  static S2Cap Full() { return S2Cap(S2Point(1, 0, 0), 2); }

  ~S2Cap() {}

  // Accessor methods.
  S2Point const& center() const { return center_; }
  double height() const { return height_; }

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
  bool is_valid() const { return S2::IsUnitLength(center_) && height_ <= 2; }

  // Return true if the cap is empty, i.e. it contains no points.
  bool is_empty() const { return height_ < 0; }

  // Return true if the cap is full, i.e. it contains all points.
  bool is_full() const { return height_ == 2; }

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
  virtual bool VirtualContainsPoint(S2Point const& p) const {
    return Contains(p);  // The same as Contains() below, just virtual.
  }

  // The point "p" should be a unit-length vector.
  bool Contains(S2Point const& p) const;

  virtual void Encode(Encoder* const encoder) const {
    LOG(FATAL) << "Unimplemented";
  }
  virtual bool Decode(Decoder* const decoder) { return false; }

  ///////////////////////////////////////////////////////////////////////
  // The following static methods are convenience functions for assertions
  // and testing purposes only.

  // Return true if two caps are identical.
  bool operator==(S2Cap const& other) const;

  // Return true if the cap center and height differ by at most "max_error"
  // from the given cap "other".
  bool ApproxEquals(S2Cap const& other, double max_error = 1e-14) const;

 private:
  // Here are some useful relationships between the cap height (h), the cap
  // radius (r), the maximum chord length from the cap's center (d), and the
  // radius of cap's base (a).
  //
  //     h = 1 - cos(r)
  //       = 2 * sin^2(r/2)
  //   d^2 = 2 * h
  //       = a^2 + h^2

  S2Cap(S2Point const& center, double height)
    : center_(center), height_(height) {
    DCHECK(is_valid());
  }

  // Return true if the cap intersects "cell", given that the cap does contain
  // any of the cell vertices (supplied in "vertices", an array of length 4).
  bool Intersects(S2Cell const& cell, S2Point const* vertices) const;

  S2Point center_;
  double height_;
};

std::ostream& operator<<(std::ostream& os, S2Cap const& cap);


//////////////////   Implementation details follow   ////////////////////


inline S2Cap::S2Cap(S2Point const& center, S1Angle radius)
    : center_(center), height_(RadiusToHeight(radius)) {
  DCHECK(is_valid());
}

inline S2Cap S2Cap::FromPoint(S2Point const& center) {
  return S2Cap(center, 0);
}

inline S2Cap S2Cap::FromCenterHeight(S2Point const& center, double height) {
  return S2Cap(center, height);
}

inline S2Cap S2Cap::FromCenterArea(S2Point const& center, double area) {
  return S2Cap(center, area / (2 * M_PI));
}

#endif  // S2GEOMETRY_S2CAP_H_
