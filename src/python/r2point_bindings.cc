#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <sstream>

#include "s2/r2.h"

namespace py = pybind11;

void bind_r2point(py::module& m) {
  py::class_<R2Point>(m, "R2Point")
      // Constructors
      .def(py::init<>(), "Default constructor creates a zero vector")
      .def(py::init<double, double>(),
           py::arg("x"), py::arg("y"),
           "Constructor that accepts x, y coordinates")
      .def(py::init([](py::tuple t) {
        if (t.size() != 2) {
          throw py::value_error("Tuple must have exactly 2 elements");
        }
        return R2Point(t[0].cast<double>(), t[1].cast<double>());
      }), py::arg("coords"), "Constructor that accepts a tuple of (x, y) coordinates")

      // Properties
      .def_property_readonly("x", py::overload_cast<>(&R2Point::x, py::const_),
                             "The x coordinate")
      .def_property_readonly("y", py::overload_cast<>(&R2Point::y, py::const_),
                             "The y coordinate")
      .def("data", [](const R2Point& self) {
        return py::make_tuple(self.x(), self.y());
      }, "Return coordinates as a tuple (x, y)")

      // Vector operations
      .def("norm", &R2Point::Norm, "Return the Euclidean norm (length)")
      .def("norm2", &R2Point::Norm2, "Return the squared Euclidean norm")
      .def("normalize", &R2Point::Normalize,
           "Return a normalized copy of the vector.\n\n"
           "Returns a unit vector if the norm is nonzero.")
      .def("dot_prod", [](const R2Point& self, const R2Point& other) {
        return self.DotProd(other);
      }, py::arg("other"), "Return the dot product with another point")
      .def("cross_prod", [](const R2Point& self, const R2Point& other) {
        return self.CrossProd(other);
      }, py::arg("other"),
           "Return the cross product with another point (scalar for 2D)")
      .def("angle", [](const R2Point& self, const R2Point& other) {
        return self.Angle(other);
      }, py::arg("other"),
           "Return the angle between this and another point (radians).\n\n"
           "The result is in the range [-pi, pi].")
      .def("ortho", &R2Point::Ortho,
           "Return a perpendicular vector (rotated 90 degrees CCW)")
      .def("fabs", &R2Point::Fabs,
           "Return a vector with component-wise absolute values")

      // Operators
      .def(py::self + py::self, "Add two points (vector addition)")
      .def(py::self - py::self, "Subtract two points (vector subtraction)")
      .def(py::self * double(), "Multiply by scalar")
      .def("__rmul__", [](const R2Point& self, double v) -> R2Point {
        return self * v;
      }, "Multiply by scalar (reversed operands)")
      .def(py::self / double(), "Divide by scalar")
      .def(-py::self, "Negate point")
      .def(py::self == py::self, "Return true if points are exactly equal")
      .def(py::self != py::self, "Return true if points are not exactly equal")
      .def(py::self += py::self, "In-place addition")
      .def(py::self -= py::self, "In-place subtraction")
      .def("__imul__", [](R2Point& self, double v) -> R2Point& {
        return self *= v;
      }, py::arg("v"), "In-place multiplication by scalar")
      .def("__itruediv__", [](R2Point& self, double v) -> R2Point& {
        return self /= v;
      }, py::arg("v"), "In-place division by scalar")

      // String representation
      .def("__repr__", [](const R2Point& p) {
        std::ostringstream oss;
        oss << "R2Point(" << p << ")";
        return oss.str();
      })
      .def("__str__", [](const R2Point& p) {
        std::ostringstream oss;
        oss << p;
        return oss.str();
      });
}
