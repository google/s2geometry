#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <sstream>

#include "s2/r2rect.h"

namespace py = pybind11;

namespace {

void MaybeThrowInvalidRect(const R2Point& lo, const R2Point& hi) {
  // Both x and y intervals must be either both empty or both non-empty.
  bool x_empty = lo.x() > hi.x();
  bool y_empty = lo.y() > hi.y();
  if (x_empty != y_empty) {
    throw py::value_error(
        "Invalid R2Rect: x and y intervals must be both empty or both "
        "non-empty");
  }
}

}  // namespace

void bind_r2rect(py::module& m) {
  py::class_<R2Rect>(m, "R2Rect")
      // Constructors
      .def(py::init<>(), "Default constructor creates an empty rectangle")
      .def(py::init([](const R2Point& lo, const R2Point& hi) {
             MaybeThrowInvalidRect(lo, hi);
             return R2Rect(lo, hi);
           }),
           py::arg("lo"), py::arg("hi"),
           "Constructor from lower-left and upper-right points.\n\n"
           "Raises ValueError if one axis interval is empty and the other "
           "is not.")
      .def(py::init([](const R1Interval& x, const R1Interval& y) {
             if (x.is_empty() != y.is_empty()) {
               throw py::value_error(
                   "Invalid R2Rect: x and y intervals must be both empty or "
                   "both non-empty");
             }
             return R2Rect(x, y);
           }),
           py::arg("x"), py::arg("y"),
           "Constructor from x and y intervals.\n\n"
           "Both intervals must be either both empty or both non-empty.\n"
           "Raises ValueError if one is empty and the other is not.")

      // Static factory methods
      .def_static("empty", &R2Rect::Empty,
                  "Returns the canonical empty rectangle")
      .def_static("from_center_size",
                  [](const R2Point& center, const R2Point& size) {
                    if (size.x() < 0 || size.y() < 0) {
                      throw py::value_error(
                          "Both components of size must be non-negative");
                    }
                    return R2Rect::FromCenterSize(center, size);
                  },
                  py::arg("center"), py::arg("size"),
                  "Construct a rectangle from a center point and size.\n\n"
                  "Both components of size must be non-negative.\n"
                  "Raises ValueError if any size component is negative.")
      .def_static("from_point", &R2Rect::FromPoint, py::arg("p"),
                  "Constructs a rectangle containing a single point")
      .def_static("from_point_pair", &R2Rect::FromPointPair,
                  py::arg("p1"), py::arg("p2"),
                  "Constructs the minimal bounding rectangle containing two "
                  "points")

      // Properties
      .def_property_readonly("lo", &R2Rect::lo,
                             "Lower-left corner of the rectangle")
      .def_property_readonly("hi", &R2Rect::hi,
                             "Upper-right corner of the rectangle")
      .def_property_readonly("x", &R2Rect::x, "The x-interval")
      .def_property_readonly("y", &R2Rect::y, "The y-interval")

      // Predicates
      .def("is_valid", &R2Rect::is_valid,
           "Return true if the rectangle is valid")
      .def("is_empty", &R2Rect::is_empty,
           "Return true if the rectangle is empty")

      // Geometric operations
      .def("vertex",
           py::overload_cast<int>(&R2Rect::GetVertex, py::const_),
           py::arg("k"),
           "Return the k-th vertex of the rectangle (k = 0,1,2,3) in CCW "
           "order.\n\n"
           "Vertex 0 is in the lower-left corner. The argument is reduced "
           "modulo 4.")
      .def("vertex_ij",
           [](const R2Rect& self, int i, int j) {
             return self.GetVertex(i, j);
           },
           py::arg("i"), py::arg("j"),
           "Return the vertex in direction i along x-axis (0=left, 1=right) "
           "and direction j along y-axis (0=down, 1=up)")
      .def("center", &R2Rect::GetCenter,
           "Return the center of the rectangle")
      .def("size", &R2Rect::GetSize,
           "Return the width and height of the rectangle.\n\n"
           "Empty rectangles have a negative width and height.")
      .def("contains", py::overload_cast<const R2Point&>(
               &R2Rect::Contains, py::const_),
           py::arg("p"),
           "Return true if the rectangle contains the given point")
      .def("interior_contains", py::overload_cast<const R2Point&>(
               &R2Rect::InteriorContains, py::const_),
           py::arg("p"),
           "Return true if the interior of the rectangle contains the given "
           "point")
      .def("contains", py::overload_cast<const R2Rect&>(
               &R2Rect::Contains, py::const_),
           py::arg("other"),
           "Return true if the rectangle contains the given other rectangle")
      .def("interior_contains", py::overload_cast<const R2Rect&>(
               &R2Rect::InteriorContains, py::const_),
           py::arg("other"),
           "Return true if the interior of this rectangle contains all points "
           "of 'other'")
      .def("intersects", &R2Rect::Intersects, py::arg("other"),
           "Return true if this rectangle and 'other' have any points in "
           "common")
      .def("interior_intersects", &R2Rect::InteriorIntersects,
           py::arg("other"),
           "Return true if the interior of this rectangle intersects any "
           "point of 'other'")
      .def("add_point", &R2Rect::AddPoint, py::arg("p"),
           "Expand the rectangle to include the given point")
      .def("add_rect", &R2Rect::AddRect, py::arg("other"),
           "Expand the rectangle to include the given other rectangle")
      .def("project", [](const R2Rect& self, const R2Point& p) {
               if (self.is_empty()) {
                 throw py::value_error(
                     "Cannot project onto an empty rectangle");
               }
               return self.Project(p);
             }, py::arg("p"),
           "Return the closest point in the rectangle to 'p'.\n\n"
           "The rectangle must be non-empty.")
      .def("expanded",
           py::overload_cast<const R2Point&>(&R2Rect::Expanded, py::const_),
           py::arg("margin"),
           "Return a rectangle expanded on each side by margin.\n\n"
           "margin.x() is applied on each side in the x-direction and "
           "margin.y() in the y-direction. If a margin component is negative, "
           "shrink the corresponding sides instead.")
      .def("union", &R2Rect::Union, py::arg("other"),
           "Return the smallest rectangle containing this rectangle and "
           "'other'")
      .def("intersection", &R2Rect::Intersection, py::arg("other"),
           "Return the smallest rectangle containing the intersection of "
           "this rectangle and 'other'")
      // Note: default value must match C++ signature in r2rect.h
      .def("approx_equals", &R2Rect::ApproxEquals,
           py::arg("other"), py::arg("max_error") = 1e-15,
           "Return true if approximately equal to 'other'")

      // Operators
      .def(py::self == py::self, "Return true if two rectangles are equal")
      .def(py::self != py::self, "Return true if two rectangles are not equal")

      // String representation
      .def("__repr__", [](const R2Rect& r) {
        std::ostringstream oss;
        oss << "R2Rect(" << r << ")";
        return oss.str();
      })
      .def("__str__", [](const R2Rect& r) {
        std::ostringstream oss;
        oss << r;
        return oss.str();
      });
}
