#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "s2/s1interval.h"

namespace py = pybind11;

void bind_s1interval(py::module& m) {
  py::class_<S1Interval>(m, "S1Interval")
      // Constructors
      .def(py::init<>(), "Default constructor (empty interval)")
      .def(py::init<double, double>(),
           py::arg("lo"), py::arg("hi"),
           "Construct interval from lo and hi bounds (in radians)")

      // Static factory methods
      .def_static("empty", &S1Interval::Empty, "Return an empty interval")
      .def_static("full", &S1Interval::Full, "Return a full interval")
      .def_static("from_point", &S1Interval::FromPoint, py::arg("p"),
                  "Construct interval containing a single point")
      .def_static("from_point_pair", &S1Interval::FromPointPair,
                  py::arg("p1"), py::arg("p2"),
                  "Construct minimal interval containing two points")

      // Properties
      .def_property("lo", &S1Interval::lo, &S1Interval::set_lo, "Lower bound")
      .def_property("hi", &S1Interval::hi, &S1Interval::set_hi, "Upper bound")
      .def("bounds", [](const S1Interval& self) {
        return py::make_tuple(self.lo(), self.hi());
      }, "Return bounds as a tuple (lo, hi)")

      // Predicates
      .def("is_valid", &S1Interval::is_valid, "Check if interval is valid")
      .def("is_full", &S1Interval::is_full, "Check if interval is full")
      .def("is_empty", &S1Interval::is_empty, "Check if interval is empty")
      .def("is_inverted", &S1Interval::is_inverted,
           "Check if interval is inverted (lo > hi)")

      // Geometric operations
      .def("get_center", &S1Interval::GetCenter, "Return center of interval")
      .def("get_length", &S1Interval::GetLength, "Return length of interval")
      .def("get_complement_center", &S1Interval::GetComplementCenter,
           "Return center of complement")
      .def("contains", py::overload_cast<double>(&S1Interval::Contains, py::const_),
           py::arg("p"), "Check if interval contains a point")
      .def("interior_contains", py::overload_cast<double>(
               &S1Interval::InteriorContains, py::const_),
           py::arg("p"), "Check if interval's interior contains a point")
      .def("contains", py::overload_cast<const S1Interval&>(
               &S1Interval::Contains, py::const_),
           py::arg("other"), "Check if interval contains another interval")
      .def("interior_contains", py::overload_cast<const S1Interval&>(
               &S1Interval::InteriorContains, py::const_),
           py::arg("other"),
           "Check if interval's interior contains another interval")
      .def("intersects", &S1Interval::Intersects, py::arg("other"),
           "Check if interval intersects another")
      .def("interior_intersects", &S1Interval::InteriorIntersects,
           py::arg("other"),
           "Check if interval's interior intersects another")
      .def("add_point", &S1Interval::AddPoint, py::arg("p"),
           "Expand interval to include a point")
      .def("project", &S1Interval::Project, py::arg("p"),
           "Return closest point in interval")
      .def("expanded", &S1Interval::Expanded, py::arg("margin"),
           "Return expanded interval")
      .def("union", &S1Interval::Union, py::arg("other"),
           "Return union with another interval")
      .def("intersection", &S1Interval::Intersection, py::arg("other"),
           "Return intersection with another interval")
      .def("complement", &S1Interval::Complement,
           "Return complement of this interval")
      .def("get_directed_hausdorff_distance", 
           &S1Interval::GetDirectedHausdorffDistance,
           py::arg("other"),
           "Return directed Hausdorff distance to another interval")
      .def("approx_equals", &S1Interval::ApproxEquals,
           py::arg("other"), py::arg("max_error") = 1e-15,
           "Check if approximately equal")

      // Operators
      .def(py::self == py::self, "Check equality")
      .def(py::self != py::self, "Check inequality")

      // String representation
      .def("__repr__", [](const S1Interval& i) {
        return "S1Interval(" + std::to_string(i.lo()) + ", " +
               std::to_string(i.hi()) + ")";
      })
      .def("__str__", [](const S1Interval& i) {
        if (i.is_empty()) return std::string("[∅]");
        if (i.is_full()) return std::string("[0, 2π)");
        return "[" + std::to_string(i.lo()) + ", " +
               std::to_string(i.hi()) + "]";
      });
}
