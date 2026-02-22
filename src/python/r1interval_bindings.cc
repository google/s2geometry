#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <sstream>

#include "s2/r1interval.h"

namespace py = pybind11;

void bind_r1interval(py::module& m) {
  py::class_<R1Interval>(m, "R1Interval")
      // Constructors
      .def(py::init<>(), "Default constructor creates an empty interval")
      .def(py::init<double, double>(),
           py::arg("lo"), py::arg("hi"),
           "Constructor that accepts the endpoints of the interval.\n\n"
           "If lo > hi, the interval is empty.")

      // Static factory methods
      .def_static("empty", &R1Interval::Empty, "Returns the empty interval")
      .def_static("from_point", &R1Interval::FromPoint, py::arg("p"),
                  "Constructs an interval containing a single point")
      .def_static("from_point_pair", &R1Interval::FromPointPair,
                  py::arg("p1"), py::arg("p2"),
                  "Constructs the minimal interval containing two points")

      // Properties
      .def_property_readonly("lo", &R1Interval::lo, "Lower bound")
      .def_property_readonly("hi", &R1Interval::hi, "Upper bound")
      .def("bounds", [](const R1Interval& self) {
        return py::make_tuple(self.lo(), self.hi());
      }, "Return bounds as a tuple (lo, hi)")

      // Predicates
      .def("is_empty", &R1Interval::is_empty,
           "Return true if the interval is empty, i.e. it contains no points")

      // Geometric operations
      .def("center", &R1Interval::GetCenter,
           "Return the center of the interval.\n\n"
           "For empty intervals, the result is arbitrary.")
      .def("length", &R1Interval::GetLength,
           "Return the length of the interval.\n\n"
           "The length of an empty interval is negative.")
      .def("contains", py::overload_cast<double>(
               &R1Interval::Contains, py::const_),
           py::arg("p"),
           "Return true if the interval contains the point 'p'")
      .def("interior_contains", py::overload_cast<double>(
               &R1Interval::InteriorContains, py::const_),
           py::arg("p"),
           "Return true if the interior of the interval contains the point 'p'")
      .def("contains", py::overload_cast<const R1Interval&>(
               &R1Interval::Contains, py::const_),
           py::arg("other"),
           "Return true if the interval contains the given interval")
      .def("interior_contains", py::overload_cast<const R1Interval&>(
               &R1Interval::InteriorContains, py::const_),
           py::arg("other"),
           "Return true if the interior of this interval contains the entire "
           "interval 'other'")
      .def("intersects", &R1Interval::Intersects, py::arg("other"),
           "Return true if the two intervals contain any points in common")
      .def("interior_intersects", &R1Interval::InteriorIntersects,
           py::arg("other"),
           "Return true if the interior of this interval intersects any point "
           "of the given interval")
      .def("add_point", &R1Interval::AddPoint, py::arg("p"),
           "Expand the interval to include the given point 'p'")
      .def("add_interval", &R1Interval::AddInterval, py::arg("other"),
           "Expand the interval to include the given interval")
      .def("project", [](const R1Interval& self, double p) {
               if (self.is_empty()) {
                 throw py::value_error("Cannot project onto an empty interval");
               }
               return self.Project(p);
             }, py::arg("p"),
           "Return the closest point in the interval to 'p'.\n\n"
           "The interval must be non-empty.")
      .def("expanded", &R1Interval::Expanded, py::arg("margin"),
           "Return interval expanded on each side by 'margin'.\n\n"
           "If 'margin' is negative, shrink the interval instead. The "
           "resulting interval may be empty. Any expansion of an empty "
           "interval remains empty.")
      .def("union", &R1Interval::Union, py::arg("other"),
           "Return the smallest interval containing this interval and 'other'")
      .def("intersection", &R1Interval::Intersection, py::arg("other"),
           "Return the intersection of this interval with 'other'")
      .def("directed_hausdorff_distance",
           &R1Interval::GetDirectedHausdorffDistance,
           py::arg("other"),
           "Return the directed Hausdorff distance to 'other'")
      // Note: default value must match C++ signature in r1interval.h
      .def("approx_equals", &R1Interval::ApproxEquals,
           py::arg("other"), py::arg("max_error") = 1e-15,
           "Return true if approximately equal to 'other'.\n\n"
           "The empty interval is considered to be positioned arbitrarily, "
           "thus any interval with length <= 2*max_error matches the empty "
           "interval.")

      // Operators
      .def(py::self == py::self, "Return true if two intervals contain the same set of points")
      .def(py::self != py::self, "Return true if two intervals do not contain the same set of points")

      // String representation
      .def("__repr__", [](const R1Interval& i) {
        std::ostringstream oss;
        oss << "R1Interval(" << i << ")";
        return oss.str();
      })
      .def("__str__", [](const R1Interval& i) {
        std::ostringstream oss;
        oss << i;
        return oss.str();
      });
}
