#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <cmath>
#include <sstream>
#include <vector>

#include "absl/hash/hash.h"
#include "s2/r1interval.h"
#include "s2/s1angle.h"
#include "s2/s1interval.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2latlng.h"
#include "s2/s2latlng_rect.h"
#include "s2/s2point.h"

namespace py = pybind11;

void MaybeThrowInvalidLatInterval(const R1Interval& lat) {
  if (std::fabs(lat.lo()) > M_PI_2 || std::fabs(lat.hi()) > M_PI_2) {
    throw py::value_error("lat interval must be within [-pi/2, pi/2]");
  }
}

void MaybeThrowEmptyMismatch(const R1Interval& lat, const S1Interval& lng) {
  if (lat.is_empty() != lng.is_empty()) {
    throw py::value_error(
        "lat and lng intervals must both be empty or both non-empty");
  }
}

void bind_s2latlng_rect(py::module& m) {
  py::class_<S2LatLngRect>(m, "S2LatLngRect",
      "A closed latitude-longitude rectangle on the sphere.\n\n"
      "Capable of representing the empty and full rectangles as well as\n"
      "single points. Uses cylindrical topology: longitudes wrap at +/-180\n"
      "degrees while latitudes are clamped to [-90, 90]. See\n"
      "s2/s2latlng_rect.h for comprehensive documentation.")

      // Constructors
      .def(py::init<>(),
           "Construct an empty S2LatLngRect.")
      .def(py::init<const S2LatLng&, const S2LatLng&>(),
           py::arg("lo"), py::arg("hi"),
           "Construct a rectangle from lo and hi corner points.\n\n"
           "If lo.lng() > hi.lng() the rectangle spans the 180-degree line.\n"
           "Both points must be normalized with lo.lat() <= hi.lat().")
      .def(py::init([](const R1Interval& lat, const S1Interval& lng) {
               MaybeThrowInvalidLatInterval(lat);
               MaybeThrowEmptyMismatch(lat, lng);
               return S2LatLngRect(lat, lng);
           }),
           py::arg("lat"), py::arg("lng"),
           "Construct a rectangle from latitude and longitude intervals.\n\n"
           "Both intervals must be empty or both non-empty. The latitude\n"
           "interval must lie within [-pi/2, pi/2] radians.\n\n"
           "Raises ValueError if either condition is violated.")

      // Factory methods
      .def_static("from_center_size", &S2LatLngRect::FromCenterSize,
           py::arg("center"), py::arg("size"),
           "Construct a rectangle of the given size centered around center.\n\n"
           "center must be normalized; size does not need to be. The latitude\n"
           "interval is clamped to [-90, 90] degrees, and the longitude\n"
           "interval is Full() if the longitude size is 360 degrees or more.")
      .def_static("from_point", &S2LatLngRect::FromPoint,
           py::arg("p"),
           "Construct a rectangle containing a single normalized point.")
      .def_static("from_point_pair", &S2LatLngRect::FromPointPair,
           py::arg("p1"), py::arg("p2"),
           "Construct the minimal bounding rectangle containing two points.\n\n"
           "Equivalent to starting with an empty rectangle and calling\n"
           "add_point() twice. Both points must be normalized.")
      .def_static("empty", &S2LatLngRect::Empty,
           "Return an empty rectangle (contains no points).")
      .def_static("full", &S2LatLngRect::Full,
           "Return the full rectangle (contains all points).")
      .def_static("full_lat", &S2LatLngRect::FullLat,
           "Return the full allowable latitude range as an R1Interval.")
      .def_static("full_lng", &S2LatLngRect::FullLng,
           "Return the full allowable longitude range as an S1Interval.")
      .def_static("from_latlngs", [](const std::vector<S2LatLng>& latlngs) {
               // Replaces the mutable AddPoint pattern: the Python interface is
               // immutable, so callers pass a list and get back a new rectangle.
               if (latlngs.empty()) return S2LatLngRect::Empty();
               S2LatLngRect result = S2LatLngRect::FromPoint(latlngs[0]);
               for (size_t i = 1; i < latlngs.size(); ++i) result.AddPoint(latlngs[i]);
               return result;
           }, py::arg("latlngs"),
           "Return the smallest rectangle containing all given S2LatLng points.\n\n"
           "Returns the empty rectangle if the list is empty.")

      // Properties
      .def_property_readonly("lat", &S2LatLngRect::lat,
           "The latitude interval as an R1Interval (in radians).")
      .def_property_readonly("lng", &S2LatLngRect::lng,
           "The longitude interval as an S1Interval (in radians).")
      .def_property_readonly("lat_lo", &S2LatLngRect::lat_lo,
           "The minimum latitude as an S1Angle.")
      .def_property_readonly("lat_hi", &S2LatLngRect::lat_hi,
           "The maximum latitude as an S1Angle.")
      .def_property_readonly("lng_lo", &S2LatLngRect::lng_lo,
           "The minimum longitude as an S1Angle.")
      .def_property_readonly("lng_hi", &S2LatLngRect::lng_hi,
           "The maximum longitude as an S1Angle.")
      .def_property_readonly("lo", &S2LatLngRect::lo,
           "The lower-left corner as an S2LatLng.")
      .def_property_readonly("hi", &S2LatLngRect::hi,
           "The upper-right corner as an S2LatLng.")

      // Predicates
      .def("is_empty", &S2LatLngRect::is_empty,
           "Return true if the rectangle contains no points.")
      .def("is_full", &S2LatLngRect::is_full,
           "Return true if the rectangle contains all points.")
      .def("is_point", &S2LatLngRect::is_point,
           "Return true if the rectangle is a single point (lo == hi).")
      .def("is_inverted", &S2LatLngRect::is_inverted,
           "Return true if the longitude interval crosses the 180-degree line.")

      // Geometric operations
      .def("vertex", &S2LatLngRect::GetVertex, py::arg("k"),
           "Return the k-th vertex in CCW order (k = 0,1,2,3).\n\n"
           "Order: lower-left, lower-right, upper-right, upper-left.\n"
           "The argument is reduced modulo 4.")
      .def("center", &S2LatLngRect::GetCenter,
           "Return the center in latitude-longitude space.\n\n"
           "Note: this is not the center of the region on the sphere.")
      .def("size", &S2LatLngRect::GetSize,
           "Return the width and height of this rectangle as an S2LatLng.\n\n"
           "Empty rectangles have a negative width and height.")
      .def("area", &S2LatLngRect::Area,
           "Return the surface area of this rectangle on the unit sphere.")
      .def("centroid", &S2LatLngRect::GetCentroid,
           "Return the true centroid multiplied by its surface area.\n\n"
           "The result is not unit length. The centroid may not be contained\n"
           "by the rectangle.")
      .def("contains", py::overload_cast<const S2LatLngRect&>(
               &S2LatLngRect::Contains, py::const_),
           py::arg("other"),
           "Return true if this rectangle contains the given rectangle.")
      .def("contains_latlng", py::overload_cast<const S2LatLng&>(
               &S2LatLngRect::Contains, py::const_),
           py::arg("ll"),
           "Return true if this rectangle contains the given S2LatLng.\n\n"
           "The argument must be normalized.")
      .def("contains_point", py::overload_cast<const S2Point&>(
               &S2LatLngRect::Contains, py::const_),
           py::arg("p"),
           "Return true if this rectangle contains the given S2Point.\n\n"
           "The point does not need to be normalized.")
      .def("contains_cell", py::overload_cast<const S2Cell&>(
               &S2LatLngRect::Contains, py::const_),
           py::arg("cell"),
           "Return true if this rectangle contains the given cell.")
      .def("interior_contains_point", py::overload_cast<const S2Point&>(
               &S2LatLngRect::InteriorContains, py::const_),
           py::arg("p"),
           "Return true if the interior contains the given point.\n\n"
           "The point does not need to be normalized.")
      .def("interior_contains_latlng", py::overload_cast<const S2LatLng&>(
               &S2LatLngRect::InteriorContains, py::const_),
           py::arg("ll"),
           "Return true if the interior contains the given S2LatLng.\n\n"
           "The argument must be normalized.")
      .def("interior_contains", py::overload_cast<const S2LatLngRect&>(
               &S2LatLngRect::InteriorContains, py::const_),
           py::arg("other"),
           "Return true if the interior of this rectangle contains all\n"
           "points of the given rectangle (including its boundary).")
      .def("intersects", py::overload_cast<const S2LatLngRect&>(
               &S2LatLngRect::Intersects, py::const_),
           py::arg("other"),
           "Return true if this rectangle and other have any points in common.")
      .def("intersects_cell", py::overload_cast<const S2Cell&>(
               &S2LatLngRect::Intersects, py::const_),
           py::arg("cell"),
           "Return true if this rectangle intersects the given cell.\n\n"
           "This is an exact test and may be expensive. Use may_intersect()\n"
           "for a cheaper approximate test.")
      .def("interior_intersects", &S2LatLngRect::InteriorIntersects,
           py::arg("other"),
           "Return true if the interior of this rectangle intersects\n"
           "any point (including the boundary) of the given rectangle.")
      .def("boundary_intersects", &S2LatLngRect::BoundaryIntersects,
           py::arg("v0"), py::arg("v1"),
           "Return true if the boundary of this rectangle intersects\n"
           "the given geodesic edge (v0, v1).")
      .def("may_intersect", &S2LatLngRect::MayIntersect, py::arg("cell"),
           "Return true if this rectangle may intersect the given cell.\n\n"
           "Cheap but not exact; use intersects_cell() for an exact test.")
      .def("expanded", &S2LatLngRect::Expanded, py::arg("margin"),
           "Return this rectangle expanded by margin.lat() on each latitude\n"
           "side and margin.lng() on each longitude side.\n\n"
           "A negative margin shrinks the rectangle. Uses lat/lng space\n"
           "topology (cylindrical). Use expanded_by_distance() to expand by\n"
           "a spherical distance.")
      .def("polar_closure", &S2LatLngRect::PolarClosure,
           "Return this rectangle with longitude range expanded to Full()\n"
           "if it contains a pole, so that all representations of any\n"
           "contained pole are included.")
      .def("union", &S2LatLngRect::Union, py::arg("other"),
           "Return the smallest rectangle containing the union of this\n"
           "rectangle and other.")
      .def("intersection", &S2LatLngRect::Intersection, py::arg("other"),
           "Return the smallest rectangle containing the intersection of\n"
           "this rectangle and other.\n\n"
           "If the intersection consists of two disjoint rectangles, returns\n"
           "a single rectangle spanning both.")
      .def("expanded_by_distance", &S2LatLngRect::ExpandedByDistance,
           py::arg("distance"),
           "Return this rectangle expanded so that it contains all points\n"
           "within the given spherical distance of its boundary.\n\n"
           "A negative distance shrinks the rectangle instead. Unlike\n"
           "expanded(), this method measures distances on the sphere.")
      .def("distance", py::overload_cast<const S2LatLngRect&>(
               &S2LatLngRect::GetDistance, py::const_),
           py::arg("other"),
           "Return the minimum spherical distance to the given rectangle.\n\n"
           "Both rectangles must be non-empty.")
      .def("distance_latlng", py::overload_cast<const S2LatLng&>(
               &S2LatLngRect::GetDistance, py::const_),
           py::arg("p"),
           "Return the minimum spherical distance from this rectangle to\n"
           "the given point. The point must be valid.")
      .def("directed_hausdorff_distance",
           py::overload_cast<const S2LatLngRect&>(
               &S2LatLngRect::GetDirectedHausdorffDistance, py::const_),
           py::arg("other"),
           "Return the directed Hausdorff distance from this rectangle to other.\n\n"
           "h(A,B) = max_{p in A} min_{q in B} d(p,q).")
      .def("hausdorff_distance", &S2LatLngRect::GetHausdorffDistance,
           py::arg("other"),
           "Return the Hausdorff distance between this rectangle and other.\n\n"
           "H(A,B) = max(h(A,B), h(B,A)).")
      .def("cap_bound", &S2LatLngRect::GetCapBound,
           "Return the smallest cap containing this rectangle.")
      .def("rect_bound", &S2LatLngRect::GetRectBound,
           "Return a bounding rectangle for this rectangle (returns self).")
      .def("cell_union_bound", [](const S2LatLngRect& self) {
               // GetCellUnionBound uses an output parameter; return by value instead.
               std::vector<S2CellId> cell_ids;
               self.GetCellUnionBound(&cell_ids);
               return cell_ids;
           },
           "Return a list of S2CellIds whose union covers this rectangle.")

      // Operators
      .def(py::self == py::self, "Return true if rectangles are equal.")
      .def(py::self != py::self, "Return true if rectangles are not equal.")
      .def("approx_equals", py::overload_cast<const S2LatLngRect&, S1Angle>(
               &S2LatLngRect::ApproxEquals, py::const_),
           py::arg("other"),
           py::arg("max_error") = S1Angle::Radians(1e-15),
           "Return true if this rectangle is approximately equal to other.\n\n"
           "Uses a uniform tolerance for both lat and lng.")
      .def("__hash__", [](const S2LatLngRect& self) {
           return absl::HashOf(self);
      })

      // String representation
      .def("__repr__", [](const S2LatLngRect& self) {
           std::ostringstream oss;
           oss << "S2LatLngRect(" << self << ")";
           return oss.str();
      })
      .def("__str__", [](const S2LatLngRect& self) {
           std::ostringstream oss;
           oss << self;
           return oss.str();
      });
}
