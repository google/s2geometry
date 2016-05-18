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

#include "s2/s2cap.h"

#include <iosfwd>

#include <glog/logging.h>

#include "s2/base/integral_types.h"
#include "s2/r1interval.h"
#include "s2/s1interval.h"
#include "s2/s2.h"
#include "s2/s2cell.h"
#include "s2/s2edgeutil.h"
#include "s2/s2latlng.h"
#include "s2/s2latlngrect.h"
#include "s2/util/math/vector.h"

using std::max;

// Multiply a positive number by this constant to ensure that the result
// of a floating point operation is at least as large as the true
// infinite-precision result.
static double const kRoundUp = 1.0 + 1.0 / (GG_LONGLONG(1) << 52);

/*static*/ double S2Cap::RadiusToHeight(S1Angle radius) {
  if (radius.radians() < 0) return -1;     // empty
  if (radius.radians() >= M_PI) return 2;  // full

  // The height of the cap can be computed as 1 - cos(r), but this isn't very
  // accurate for angles close to zero (where cos(r) is almost 1).  The
  // formula below has much better precision.
  double d = sin(0.5 * radius.radians());
  return 2 * d * d;
}

S1Angle S2Cap::GetRadius() const {
  // This could also be computed as acos(1 - height_), but the following
  // formula is much more accurate when the cap height is small.  It
  // follows from the relationship h = 1 - cos(r) = 2 sin^2(r/2).
  if (is_empty()) return S1Angle::Radians(-1);
  return S1Angle::Radians(2 * asin(sqrt(0.5 * height_)));
}

double S2Cap::GetArea() const {
  return 2 * M_PI * max(0.0, height_);
}

S2Point S2Cap::GetCentroid() const {
  // From symmetry, the centroid of the cap must be somewhere on the line
  // from the origin to the center of the cap on the surface of the sphere.
  // When a sphere is divided into slices of constant thickness by a set of
  // parallel planes, all slices have the same surface area. This implies
  // that the radial component of the centroid is simply the midpoint of the
  // range of radial distances spanned by the cap. That is easily computed
  // from the cap height.

  if (is_empty()) return S2Point();
  double r = 1.0 - height_ / 2.0;
  return center_ * r * GetArea();
}

S2Cap S2Cap::Complement() const {
  // The complement of a full cap is an empty cap, not a singleton.
  // Also make sure that the complement of an empty cap has height 2.
  double height = is_full() ? -1 : 2 - max(height_, 0.0);
  return S2Cap::FromCenterHeight(-center_, height);
}

bool S2Cap::Contains(S2Cap const& other) const {
  if (is_full() || other.is_empty()) return true;
  return GetRadius() >= S1Angle(center_, other.center_) + other.GetRadius();
}

bool S2Cap::Intersects(S2Cap const& other) const {
  if (is_empty() || other.is_empty()) return false;

  return GetRadius() + other.GetRadius() >= S1Angle(center_, other.center_);
}

bool S2Cap::InteriorIntersects(S2Cap const& other) const {
  // Make sure this cap has an interior and the other cap is non-empty.
  if (height_ <= 0 || other.is_empty()) return false;

  return GetRadius() + other.GetRadius() > S1Angle(center_, other.center_);
}

void S2Cap::AddPoint(S2Point const& p) {
  // Compute the squared chord length, then convert it into a height.
  DCHECK(S2::IsUnitLength(p));
  if (is_empty()) {
    center_ = p;
    height_ = 0;
  } else {
    // To make sure that the resulting cap actually includes this point,
    // we need to round up the distance calculation.  That is, after
    // calling cap.AddPoint(p), cap.Contains(p) should be true.
    double dist2 = (center_ - p).Norm2();
    height_ = max(height_, kRoundUp * 0.5 * dist2);
  }
}

void S2Cap::AddCap(S2Cap const& other) {
  if (is_empty()) {
    *this = other;
  } else {
    // See comments for AddPoint().  This could be optimized by doing the
    // calculation in terms of cap heights rather than cap opening angles.
    S1Angle radius = S1Angle(center_, other.center_) + other.GetRadius();
    height_ = max(height_, kRoundUp * RadiusToHeight(radius));
  }
}

S2Cap S2Cap::Expanded(S1Angle distance) const {
  DCHECK_GE(distance.radians(), 0);
  if (is_empty()) return Empty();
  return S2Cap(center_, GetRadius() + distance);
}

S2Cap S2Cap::Union(S2Cap const& other) const {
  if (height() < other.height()) {
    return other.Union(*this);
  }
  if (is_full() || other.is_empty()) {
    return *this;
  }
  S1Angle this_radius = GetRadius();
  S1Angle other_radius = other.GetRadius();
  S1Angle distance(center(), other.center());
  if (this_radius >= distance + other_radius) {
    return *this;
  } else {
    S1Angle result_radius = 0.5 * (distance + this_radius + other_radius);
    S2Point result_center = S2EdgeUtil::InterpolateAtDistance(
        0.5 * (distance - this_radius + other_radius),
        center(),
        other.center());
    return S2Cap(result_center, result_radius);
  }
}

S2Cap* S2Cap::Clone() const {
  return new S2Cap(*this);
}

S2Cap S2Cap::GetCapBound() const {
  return *this;
}

S2LatLngRect S2Cap::GetRectBound() const {
  if (is_empty()) return S2LatLngRect::Empty();

  // Convert the center to a (lat,lng) pair, and compute the cap angle.
  S2LatLng center_ll(center_);
  double cap_angle = GetRadius().radians();

  bool all_longitudes = false;
  double lat[2], lng[2];
  lng[0] = -M_PI;
  lng[1] = M_PI;

  // Check whether cap includes the south pole.
  lat[0] = center_ll.lat().radians() - cap_angle;
  if (lat[0] <= -M_PI_2) {
    lat[0] = -M_PI_2;
    all_longitudes = true;
  }
  // Check whether cap includes the north pole.
  lat[1] = center_ll.lat().radians() + cap_angle;
  if (lat[1] >= M_PI_2) {
    lat[1] = M_PI_2;
    all_longitudes = true;
  }
  if (!all_longitudes) {
    // Compute the range of longitudes covered by the cap.  We use the law
    // of sines for spherical triangles.  Consider the triangle ABC where
    // A is the north pole, B is the center of the cap, and C is the point
    // of tangency between the cap boundary and a line of longitude.  Then
    // C is a right angle, and letting a,b,c denote the sides opposite A,B,C,
    // we have sin(a)/sin(A) = sin(c)/sin(C), or sin(A) = sin(a)/sin(c).
    // Here "a" is the cap angle, and "c" is the colatitude (90 degrees
    // minus the latitude).  This formula also works for negative latitudes.
    //
    // The formula for sin(a) follows from the relationship h = 1 - cos(a).

    double sin_a = sqrt(height_ * (2 - height_));
    double sin_c = cos(center_ll.lat().radians());
    if (sin_a <= sin_c) {
      double angle_A = asin(sin_a / sin_c);
      lng[0] = remainder(center_ll.lng().radians() - angle_A, 2 * M_PI);
      lng[1] = remainder(center_ll.lng().radians() + angle_A, 2 * M_PI);
    }
  }
  return S2LatLngRect(R1Interval(lat[0], lat[1]),
                      S1Interval(lng[0], lng[1]));
}

bool S2Cap::Intersects(S2Cell const& cell, S2Point const* vertices) const {
  // Return true if this cap intersects any point of 'cell' excluding its
  // vertices (which are assumed to already have been checked).

  // If the cap is a hemisphere or larger, the cell and the complement of the
  // cap are both convex.  Therefore since no vertex of the cell is contained,
  // no other interior point of the cell is contained either.
  if (height_ >= 1) return false;

  // We need to check for empty caps due to the center check just below.
  if (is_empty()) return false;

  // Optimization: return true if the cell contains the cap center.  (This
  // allows half of the edge checks below to be skipped.)
  if (cell.Contains(center_)) return true;

  // At this point we know that the cell does not contain the cap center,
  // and the cap does not contain any cell vertex.  The only way that they
  // can intersect is if the cap intersects the interior of some edge.

  double sin2_angle = height_ * (2 - height_);  // sin^2(cap_angle)
  for (int k = 0; k < 4; ++k) {
    S2Point edge = cell.GetEdgeRaw(k);
    double dot = center_.DotProd(edge);
    if (dot > 0) {
      // The center is in the interior half-space defined by the edge.  We don't
      // need to consider these edges, since if the cap intersects this edge
      // then it also intersects the edge on the opposite side of the cell
      // (because we know the center is not contained with the cell).
      continue;
    }
    // The Norm2() factor is necessary because "edge" is not normalized.
    if (dot * dot > sin2_angle * edge.Norm2()) {
      return false;  // Entire cap is on the exterior side of this edge.
    }
    // Otherwise, the great circle containing this edge intersects
    // the interior of the cap.  We just need to check whether the point
    // of closest approach occurs between the two edge endpoints.
    Vector3_d dir = edge.CrossProd(center_);
    if (dir.DotProd(vertices[k]) < 0 && dir.DotProd(vertices[(k+1)&3]) > 0)
      return true;
  }
  return false;
}

bool S2Cap::Contains(S2Cell const& cell) const {
  // If the cap does not contain all cell vertices, return false.
  // We check the vertices before taking the Complement() because we can't
  // accurately represent the complement of a very small cap (a height
  // of 2-epsilon is rounded off to 2).
  S2Point vertices[4];
  for (int k = 0; k < 4; ++k) {
    vertices[k] = cell.GetVertex(k);
    if (!Contains(vertices[k])) return false;
  }
  // Otherwise, return true if the complement of the cap does not intersect
  // the cell.  (This test is slightly conservative, because technically we
  // want Complement().InteriorIntersects() here.)
  return !Complement().Intersects(cell, vertices);
}

bool S2Cap::MayIntersect(S2Cell const& cell) const {
  // If the cap contains any cell vertex, return true.
  S2Point vertices[4];
  for (int k = 0; k < 4; ++k) {
    vertices[k] = cell.GetVertex(k);
    if (Contains(vertices[k])) return true;
  }
  return Intersects(cell, vertices);
}

bool S2Cap::Contains(S2Point const& p) const {
  DCHECK(S2::IsUnitLength(p));
  return (center_ - p).Norm2() <= 2 * height_;
}

bool S2Cap::InteriorContains(S2Point const& p) const {
  DCHECK(S2::IsUnitLength(p));
  return is_full() || (center_ - p).Norm2() < 2 * height_;
}

bool S2Cap::operator==(S2Cap const& other) const {
  return (center_ == other.center_ && height_ == other.height_) ||
         (is_empty() && other.is_empty()) ||
         (is_full() && other.is_full());
}

bool S2Cap::ApproxEquals(S2Cap const& other, double max_error) const {
  return (S2::ApproxEquals(center_, other.center_, max_error) &&
          fabs(height_ - other.height_) <= max_error) ||
         (is_empty() && other.height_ <= max_error) ||
         (other.is_empty() && height_ <= max_error) ||
         (is_full() && other.height_ >= 2 - max_error) ||
         (other.is_full() && height_ >= 2 - max_error);
}

std::ostream& operator<<(std::ostream& os, S2Cap const& cap) {
  return os << "[Center=" << cap.center()
            << ", Radius=" << cap.GetRadius() << "]";
}
