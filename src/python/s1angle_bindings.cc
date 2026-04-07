#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <sstream>
#include <string>

#include "s2/s1angle.h"
#include "s2/s2point.h"

namespace py = pybind11;

namespace {

void MaybeThrowNotNormalized(const S1Angle& angle) {
  if (!angle.IsNormalizedAngle()) {
    throw py::value_error("Angle " + std::to_string(angle.degrees()) +
                          " degrees is not in the normalized range (-180, 180]");
  }
}

}  // namespace

void bind_s1angle(py::module& m) {
  py::class_<S1Angle>(m, "S1Angle",
      "Represents a one-dimensional angle.\n\n"
      "The internal representation is a double-precision value in radians,\n"
      "so conversion to and from radians is exact. Conversions between\n"
      "degrees and radians are not always exact due to floating-point\n"
      "arithmetic (e.g. from_degrees(60).degrees() != 60). Use e5/e6/e7\n"
      "representations for exact discrete comparisons.\n\n"
      "See s2/s1angle.h for comprehensive documentation including exact\n"
      "conversion guarantees and edge cases.")
      // Constructors
      .def(py::init<>(), "Default constructor creates a zero angle")
      .def(py::init<const S2Point&, const S2Point&>(),
           py::arg("x"), py::arg("y"),
           "Construct the angle between two points.\n\n"
           "This is also equal to the distance between the points on the\n"
           "unit sphere. The points do not need to be normalized.")

      // Factory methods
      .def_static("from_radians", &S1Angle::Radians, py::arg("radians"),
                  "Construct an angle from its measure in radians.\n\n"
                  "This conversion is exact.")
      .def_static("from_degrees", &S1Angle::Degrees, py::arg("degrees"),
                  "Construct an angle from its measure in degrees.\n\n"
                  "Note: the round-trip from_degrees(x).degrees() is not\n"
                  "always exact. For example, from_degrees(60).degrees() != 60.")
      .def_static("from_e5", &S1Angle::E5, py::arg("e5"),
                  "Construct an angle from its E5 representation.\n\n"
                  "E5 is degrees multiplied by 1e5 and rounded to the\n"
                  "nearest integer.\n\n"
                  "Note: E5 does not share the exact conversion guarantees\n"
                  "of E6/E7. Avoid testing E5 values for exact equality\n"
                  "with other formats.")
      .def_static("from_e6", &S1Angle::E6, py::arg("e6"),
                  "Construct an angle from its E6 representation.\n\n"
                  "E6 is degrees multiplied by 1e6 and rounded to the\n"
                  "nearest integer.\n\n"
                  "For any integer n: from_degrees(n) == from_e6(1000000 * n).")
      .def_static("from_e7", &S1Angle::E7, py::arg("e7"),
                  "Construct an angle from its E7 representation.\n\n"
                  "E7 is degrees multiplied by 1e7 and rounded to the\n"
                  "nearest integer.\n\n"
                  "For any integer n: from_degrees(n) == from_e7(10000000 * n).")
      .def_static("zero", &S1Angle::Zero, "Return a zero angle")
      .def_static("infinity", &S1Angle::Infinity,
                  "Return an angle larger than any finite angle")

      // Geometric operations
      .def("radians", &S1Angle::radians,
           "Return the angle in radians.\n\n"
           "This is the internal representation, so the conversion is exact.")
      .def("degrees", &S1Angle::degrees,
           "Return the angle in degrees.\n\n"
           "Note: from_degrees(x).degrees() is not always exactly x due to\n"
           "the intermediate conversion to radians.")
      .def("e5", [](const S1Angle& self) {
               MaybeThrowNotNormalized(self);
               return self.e5();
           },
           "Return the E5 representation (degrees * 1e5, rounded).\n\n"
           "The angle must be in the normalized range (-180, 180] degrees.\n"
           "Raises ValueError if out of range.")
      .def("e6", [](const S1Angle& self) {
               MaybeThrowNotNormalized(self);
               return self.e6();
           },
           "Return the E6 representation (degrees * 1e6, rounded).\n\n"
           "The angle must be in the normalized range (-180, 180] degrees.\n"
           "Raises ValueError if out of range.")
      .def("e7", [](const S1Angle& self) {
               MaybeThrowNotNormalized(self);
               return self.e7();
           },
           "Return the E7 representation (degrees * 1e7, rounded).\n\n"
           "The angle must be in the normalized range (-180, 180] degrees.\n"
           "Raises ValueError if out of range.")
      .def("abs", &S1Angle::abs, "Return the absolute value of the angle")
      .def("normalized", &S1Angle::Normalized,
           "Return the angle normalized to the range (-180, 180] degrees")
      .def("sin_cos", [](const S1Angle& self) {
               auto sc = self.SinCos();
               return py::make_tuple(sc.sin, sc.cos);
           },
           "Return (sin, cos) of the angle as a tuple.\n\n"
           "This may be more efficient than calling sin and cos separately.")

      // Operators
      .def(py::self == py::self, "Return true if angles are exactly equal")
      .def(py::self != py::self, "Return true if angles are not exactly equal")
      .def(py::self < py::self, "Return true if this angle is less than other")
      .def(py::self > py::self,
           "Return true if this angle is greater than other")
      .def(py::self <= py::self,
           "Return true if this angle is less than or equal to other")
      .def(py::self >= py::self,
           "Return true if this angle is greater than or equal to other")
      .def(-py::self, "Negate angle")
      .def(py::self + py::self, "Add two angles")
      .def(py::self - py::self, "Subtract two angles")
      .def(py::self * double(), "Multiply angle by scalar")
      .def("__rmul__", [](const S1Angle& self, double m) {
        return m * self;
      }, "Multiply angle by scalar (reversed operands)")
      .def(py::self / double(), "Divide angle by scalar")
      .def("__truediv__", [](const S1Angle& a, const S1Angle& b) -> double {
        return a / b;
      }, py::arg("other"), "Divide two angles, returning a scalar ratio")
      .def(py::self += py::self, "In-place addition")
      .def(py::self -= py::self, "In-place subtraction")
      .def("__imul__", [](S1Angle& self, double m) -> S1Angle& {
        return self *= m;
      }, py::arg("m"), "In-place multiplication by scalar")
      .def("__itruediv__", [](S1Angle& self, double m) -> S1Angle& {
        return self /= m;
      }, py::arg("m"), "In-place division by scalar")

      // String representation
      .def("__repr__", [](S1Angle a) {
        std::ostringstream oss;
        oss << "S1Angle(" << a << ")";
        return oss.str();
      })
      .def("__str__", [](S1Angle a) {
        std::ostringstream oss;
        oss << a;
        return oss.str();
      });

  // Module-level functions
  m.def("sin", [](S1Angle a) { return sin(a); }, py::arg("angle"),
        "Return the sine of an S1Angle");
  m.def("cos", [](S1Angle a) { return cos(a); }, py::arg("angle"),
        "Return the cosine of an S1Angle");
  m.def("tan", [](S1Angle a) { return tan(a); }, py::arg("angle"),
        "Return the tangent of an S1Angle");
}
