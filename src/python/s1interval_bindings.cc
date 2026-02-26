#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <sstream>

#include "s2/s1interval.h"

namespace py = pybind11;

namespace {

void MaybeThrowInvalidPoint(double p) {
  if (!S1Interval::IsValidPoint(p)) throw py::value_error("Invalid S1 point: " + std::to_string(p));
}

}  // namespace

void bind_s1interval(py::module& m) {
  py::class_<S1Interval>(m, "S1Interval")
      // Constructors
      .def(py::init<>(), "Default constructor creates an empty interval")
      .def(py::init([](double lo, double hi) {
             MaybeThrowInvalidPoint(lo);
             MaybeThrowInvalidPoint(hi);
             return S1Interval(lo, hi);
           }),
           py::arg("lo"), py::arg("hi"),
           "Constructor that accepts the endpoints of the interval.\n\n"
           "Both endpoints must be in the range -Pi to Pi inclusive.\n"
           "Raises ValueError if either bound is outside that range.")

      // Static factory methods
      .def_static("empty", &S1Interval::Empty, "Returns the empty interval")
      .def_static("full", &S1Interval::Full, "Returns the full interval")
      .def_static("from_point", &S1Interval::FromPoint, py::arg("p"),
                  "Constructs an interval containing a single point")
      .def_static("from_point_pair",
                  [](double p1, double p2) {
                    MaybeThrowInvalidPoint(p1);
                    MaybeThrowInvalidPoint(p2);
                    return S1Interval::FromPointPair(p1, p2);
                  },
                  py::arg("p1"), py::arg("p2"),
                  "Constructs the minimal interval containing two points")

      // Properties
      .def_property_readonly("lo", &S1Interval::lo, "Lower bound")
      .def_property_readonly("hi", &S1Interval::hi, "Upper bound")
      .def("bounds", [](const S1Interval& self) {
        return py::make_tuple(self.lo(), self.hi());
      }, "Return bounds as a tuple (lo, hi)")

      // Predicates
      .def("is_full", &S1Interval::is_full,
           "Return true if the interval contains all points on the unit circle")
      .def("is_empty", &S1Interval::is_empty,
           "Return true if the interval is empty, i.e. it contains no points")
      .def("is_inverted", &S1Interval::is_inverted,
           "Return true if lo() > hi(). (This is true for empty intervals.)")

      // Geometric operations
      .def("center", &S1Interval::GetCenter,
           "Return the midpoint of the interval.\n\n"
           "For full and empty intervals, the result is arbitrary.")
      .def("length", &S1Interval::GetLength,
           "Return the length of the interval.\n\n"
           "The length of an empty interval is negative.")
      .def("complement_center", &S1Interval::GetComplementCenter,
           "Return the midpoint of the complement of the interval.\n\n"
           "For full and empty intervals, the result is arbitrary. For a\n"
           "singleton interval, the result is its antipodal point on S1.")
      .def("contains", [](const S1Interval& self, double p) {
               MaybeThrowInvalidPoint(p);
               return self.Contains(p);
             }, py::arg("p"),
           "Return true if the interval (which is closed) contains the point 'p'")
      .def("interior_contains", [](const S1Interval& self, double p) {
               MaybeThrowInvalidPoint(p);
               return self.InteriorContains(p);
             }, py::arg("p"),
           "Return true if the interior of the interval contains the point 'p'")
      .def("contains", py::overload_cast<const S1Interval&>(
               &S1Interval::Contains, py::const_),
           py::arg("other"),
           "Return true if the interval contains the given interval 'y'")
      .def("interior_contains", py::overload_cast<const S1Interval&>(
               &S1Interval::InteriorContains, py::const_),
           py::arg("other"),
           "Return true if the interior of this interval contains the entire interval 'y'")
      .def("intersects", &S1Interval::Intersects, py::arg("other"),
           "Return true if the two intervals contain any points in common")
      .def("interior_intersects", &S1Interval::InteriorIntersects,
           py::arg("other"),
           "Return true if the interior of this interval contains any point of 'y'")
      .def("add_point", [](S1Interval& self, double p) {
               MaybeThrowInvalidPoint(p);
               self.AddPoint(p);
             }, py::arg("p"),
           "Expand the interval to contain the given point 'p'.\n\n"
           "The point should be an angle in the range [-Pi, Pi].")
      .def("project", [](const S1Interval& self, double p) {
               if (self.is_empty()) throw py::value_error("Invalid S1Interval");
               MaybeThrowInvalidPoint(p);
               return self.Project(p);
             }, py::arg("p"),
           "Return the closest point in the interval to 'p'.\n\n"
           "The interval must be non-empty.")
      .def("expanded", &S1Interval::Expanded, py::arg("margin"),
           "Return interval expanded on each side by 'margin' (radians).\n\n"
           "If 'margin' is negative, shrink the interval instead. The resulting\n"
           "interval may be empty or full. Any expansion of a full interval remains\n"
           "full, and any expansion of an empty interval remains empty.")
      .def("union", &S1Interval::Union, py::arg("other"),
           "Return the smallest interval containing this interval and 'y'")
      .def("intersection", &S1Interval::Intersection, py::arg("other"),
           "Return the smallest interval containing the intersection with 'y'.\n\n"
           "Note that the region of intersection may consist of two disjoint intervals.")
      .def("complement", &S1Interval::Complement,
           "Return the complement of the interior of the interval")
      .def("directed_hausdorff_distance", 
           &S1Interval::GetDirectedHausdorffDistance,
           py::arg("other"),
           "Return the directed Hausdorff distance to 'y'")
      // Note: default value must match C++ signature in s1interval.h
      .def("approx_equals", &S1Interval::ApproxEquals,
           py::arg("other"), py::arg("max_error") = 1e-15,
           "Return true if approximately equal to 'y'.\n\n"
           "Two intervals are approximately equal if each endpoint can be moved\n"
           "by at most 'max_error' (radians) to match the other interval.")

      // Operators
      .def(py::self == py::self, "Return true if two intervals contain the same set of points")
      .def(py::self != py::self, "Return true if two intervals do not contain the same set of points")

      // String representation
      .def("__repr__", [](const S1Interval& i) {
        std::ostringstream oss;
        oss << "S1Interval(" << i << ")";
        return oss.str();
      })
      .def("__str__", [](const S1Interval& i) {
        std::ostringstream oss;
        oss << i;
        return oss.str();
      });
}
