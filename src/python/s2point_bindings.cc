#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <sstream>

#include "s2/s2point.h"

namespace py = pybind11;

void bind_s2point(py::module& m) {
  py::class_<S2Point>(m, "S2Point")
      // Constructors
      .def(py::init<>(), "Default constructor creates a zero vector")
      .def(py::init<double, double, double>(),
           py::arg("x"), py::arg("y"), py::arg("z"),
           "Constructor that accepts x, y, z coordinates")
      .def(py::init([](py::tuple t) {
        if (t.size() != 3) {
          throw py::value_error("Tuple must have exactly 3 elements");
        }
        return S2Point(t[0].cast<double>(), t[1].cast<double>(), t[2].cast<double>());
      }), py::arg("coords"), "Constructor that accepts a tuple of (x, y, z) coordinates")

      // Properties
      .def_property_readonly("x", py::overload_cast<>(&S2Point::x, py::const_),
                             "The x coordinate")
      .def_property_readonly("y", py::overload_cast<>(&S2Point::y, py::const_),
                             "The y coordinate")
      .def_property_readonly("z", py::overload_cast<>(&S2Point::z, py::const_),
                             "The z coordinate")
      .def("data", [](const S2Point& self) {
        return py::make_tuple(self.x(), self.y(), self.z());
      }, "Return coordinates as a tuple (x, y, z)")

      // Vector operations
      .def("norm", &S2Point::Norm, "Return the Euclidean norm (length)")
      .def("norm2", &S2Point::Norm2, "Return the squared Euclidean norm (dot product with itself)")
      .def("normalize", &S2Point::Normalize,
           "Return a normalized copy of the vector.\n\n"
           "Returns a unit vector if the norm is nonzero.")
      .def("dot_prod", [](const S2Point& self, const S2Point& other) {
        return self.DotProd(other);
      }, py::arg("other"), "Return the dot product with another point")
      .def("cross_prod", [](const S2Point& self, const S2Point& other) -> S2Point {
        return self.CrossProd(other);
      }, py::arg("other"), "Return the cross product with another point")
      .def("angle", [](const S2Point& self, const S2Point& other) {
        return self.Angle(other);
      }, py::arg("other"),
           "Return the angle between this and another point (radians).\n\n"
           "The result is in the range [0, pi]. If either vector is zero-length,\n"
           "or nearly zero-length, the result will be zero.")

      // Operators
      .def(py::self + py::self, "Add two points (vector addition)")
      .def(py::self - py::self, "Subtract two points (vector subtraction)")
      .def(py::self * double(), "Multiply by scalar")
      .def("__rmul__", [](const S2Point& self, double v) -> S2Point {
        return self * v;
      }, "Multiply by scalar (reversed operands)")
      .def(py::self / double(), "Divide by scalar")
      .def(-py::self, "Negate point")
      .def(py::self == py::self, "Return true if points are exactly equal")
      .def(py::self != py::self, "Return true if points are not exactly equal")
      .def(py::self += py::self, "In-place addition")
      .def(py::self -= py::self, "In-place subtraction")
      .def("__imul__", [](S2Point& self, double v) -> S2Point& {
        return self *= v;
      }, py::arg("v"), "In-place multiplication by scalar")
      .def("__itruediv__", [](S2Point& self, double v) -> S2Point& {
        return self /= v;
      }, py::arg("v"), "In-place division by scalar")

      // String representation
      // __repr__ prefixes class name, __str__ delegates to C++ operator<<
      .def("__repr__", [](const S2Point& p) {
        std::ostringstream oss;
        oss << "S2Point(" << p << ")";
        return oss.str();
      })
      .def("__str__", [](const S2Point& p) {
        std::ostringstream oss;
        oss << p;
        return oss.str();
      });
}
