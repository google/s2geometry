#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "s2/s2point.h"

namespace py = pybind11;

void bind_s2point(py::module& m) {
  py::class_<S2Point>(m, "S2Point")
      // Constructors
      .def(py::init<>(), "Default constructor")
      .def(py::init<double, double, double>(),
           py::arg("x"), py::arg("y"), py::arg("z"),
           "Construct S2Point from x, y, z coordinates")

      // Accessors
      .def_property_readonly("x", py::overload_cast<>(&S2Point::x, py::const_),
                             "x coordinate")
      .def_property_readonly("y", py::overload_cast<>(&S2Point::y, py::const_),
                             "y coordinate")
      .def_property_readonly("z", py::overload_cast<>(&S2Point::z, py::const_),
                             "z coordinate")

      // Vector operations
      .def("norm", &S2Point::Norm, "Return the Euclidean norm (length)")
      .def("norm2", &S2Point::Norm2, "Return the squared Euclidean norm")
      .def("normalize", &S2Point::Normalize, "Return a normalized copy")
      .def("dot_prod", [](const S2Point& self, const S2Point& other) {
        return self.DotProd(other);
      }, py::arg("other"), "Return the dot product with another point")
      .def("cross_prod", [](const S2Point& self, const S2Point& other) -> S2Point {
        return self.CrossProd(other);
      }, py::arg("other"), "Return the cross product with another point")
      .def("angle", &S2Point::Angle, py::arg("other"),
           "Return the angle to another point")

      // Operators
      .def(py::self + py::self, "Add two points")
      .def(py::self - py::self, "Subtract two points")
      .def(py::self * double(), "Multiply by scalar")
      .def("__rmul__", [](const S2Point& self, double v) -> S2Point {
        return self * v;
      }, "Multiply by scalar")
      .def(py::self / double(), "Divide by scalar")
      .def(-py::self, "Negate point")
      .def(py::self == py::self, "Check equality")
      .def(py::self != py::self, "Check inequality")
      .def(py::self += py::self, "In-place addition")
      .def(py::self -= py::self, "In-place subtraction")
      .def("__imul__", [](S2Point& self, double v) -> S2Point& {
        return self *= v;
      }, py::arg("v"), "In-place multiplication")
      .def("__itruediv__", [](S2Point& self, double v) -> S2Point& {
        return self /= v;
      }, py::arg("v"), "In-place division")

      // String representation
      .def("__repr__", [](const S2Point& p) {
        return "S2Point(" + std::to_string(p.x()) + ", " +
               std::to_string(p.y()) + ", " + std::to_string(p.z()) + ")";
      })
      .def("__str__", [](const S2Point& p) {
        return "(" + std::to_string(p.x()) + ", " +
               std::to_string(p.y()) + ", " + std::to_string(p.z()) + ")";
      });
}
