#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <cmath>
#include <sstream>

#include "absl/strings/str_cat.h"
#include "s2/s1angle.h"
#include "s2/s1chord_angle.h"
#include "s2/s2point.h"
#include "s2/s2pointutil.h"

namespace py = pybind11;

namespace {

void MaybeThrowNotUnitLength(const S2Point& p, const char* name) {
  if (!S2::IsUnitLength(p)) {
    throw py::value_error(
        absl::StrCat(name, " must be a unit-length vector (norm=",
                     p.Norm(), ")"));
  }
}

void MaybeThrowSpecialForArithmetic(const S1ChordAngle& a, const char* name) {
  if (a.is_special()) {
    throw py::value_error(
        absl::StrCat(name, " must not be Negative() or Infinity() "
                     "for arithmetic"));
  }
}

}  // namespace

void bind_s1chord_angle(py::module& m) {
  py::class_<S1ChordAngle>(m, "S1ChordAngle",
      "Represents the angle subtended by a chord on the unit sphere.\n\n"
      "S1ChordAngle can represent angles between 0 and Pi radians. It is\n"
      "more efficient than S1Angle for computing and comparing distances,\n"
      "but loses some accuracy as the angle approaches Pi radians.\n\n"
      "See s2/s1chord_angle.h for comprehensive documentation, including\n"
      "accuracy analysis and guidance on when to prefer S1Angle.")
      // Constructors
      .def(py::init<>(), "Default constructor creates a zero chord angle")
      .def(py::init([](const S2Point& x, const S2Point& y) {
               MaybeThrowNotUnitLength(x, "x");
               MaybeThrowNotUnitLength(y, "y");
               return S1ChordAngle(x, y);
           }),
           py::arg("x"), py::arg("y"),
           "Construct the chord angle between two unit-length points.\n\n"
           "Raises ValueError if either point is not unit-length.")
      .def(py::init<S1Angle>(), py::arg("angle"),
           "Construct from an S1Angle.\n\n"
           "Angles outside [0, Pi] are mapped as follows:\n"
           "  Infinity() -> Infinity()\n"
           "  negative   -> Negative()\n"
           "  > Pi       -> Straight()\n"
           "This conversion is relatively expensive; prefer to convert at\n"
           "the boundaries of your algorithm.")

      // Factory methods
      .def_static("from_radians", &S1ChordAngle::Radians, py::arg("radians"),
                  "Construct a chord angle from an angle in radians")
      .def_static("from_degrees", &S1ChordAngle::Degrees, py::arg("degrees"),
                  "Construct a chord angle from an angle in degrees")
      .def_static("from_e5", &S1ChordAngle::E5, py::arg("e5"),
                  "Construct a chord angle from the E5 representation")
      .def_static("from_e6", &S1ChordAngle::E6, py::arg("e6"),
                  "Construct a chord angle from the E6 representation")
      .def_static("from_e7", &S1ChordAngle::E7, py::arg("e7"),
                  "Construct a chord angle from the E7 representation")
      .def_static("zero", &S1ChordAngle::Zero, "Return the zero chord angle")
      .def_static("right", &S1ChordAngle::Right,
                  "Return a 90-degree chord angle")
      .def_static("straight", &S1ChordAngle::Straight,
                  "Return a 180-degree chord angle (the maximum finite value)")
      .def_static("infinity", &S1ChordAngle::Infinity,
                  "Return a chord angle larger than any finite chord angle")
      .def_static("negative", &S1ChordAngle::Negative,
                  "Return a chord angle smaller than Zero()")

      // Properties
      .def_property_readonly("radians", &S1ChordAngle::radians,
           "The angle in radians.\n\n"
           "Note: this performs a trigonometric conversion and should be\n"
           "avoided in inner loops.")
      .def_property_readonly("degrees", &S1ChordAngle::degrees,
           "The angle in degrees.\n\n"
           "Note: this performs a trigonometric conversion and should be\n"
           "avoided in inner loops.")
      .def_property_readonly("e5", &S1ChordAngle::e5,
                             "The E5 representation (degrees * 1e5, rounded)")
      .def_property_readonly("e6", &S1ChordAngle::e6,
                             "The E6 representation (degrees * 1e6, rounded)")
      .def_property_readonly("e7", &S1ChordAngle::e7,
                             "The E7 representation (degrees * 1e7, rounded)")

      // Predicates
      .def("is_zero", &S1ChordAngle::is_zero,
           "Return true if this is exactly zero")
      .def("is_negative", &S1ChordAngle::is_negative,
           "Return true if this is less than zero (e.g. Negative())")
      .def("is_infinity", &S1ChordAngle::is_infinity,
           "Return true if this is the Infinity() sentinel")
      .def("is_special", &S1ChordAngle::is_special,
           "Return true if this is Negative() or Infinity()")
      .def("is_valid", &S1ChordAngle::is_valid,
           "Return true if the internal representation is valid.\n\n"
           "Negative() and Infinity() are both considered valid.")

      // Geometric operations
      .def("to_angle", &S1ChordAngle::ToAngle,
           "Convert to an S1Angle.\n\n"
           "Infinity() converts to S1Angle::Infinity(); Negative() converts\n"
           "to an unspecified negative S1Angle. Uses trigonometric functions\n"
           "and should be avoided in inner loops.")
      .def("sin", [](const S1ChordAngle& self) { return sin(self); },
           "Return the sine of the chord angle.\n\n"
           "More accurate and efficient than converting to S1Angle first.")
      .def("cos", [](const S1ChordAngle& self) { return cos(self); },
           "Return the cosine of the chord angle.\n\n"
           "More accurate and efficient than converting to S1Angle first.")
      .def("tan", [](const S1ChordAngle& self) { return tan(self); },
           "Return the tangent of the chord angle.\n\n"
           "More accurate and efficient than converting to S1Angle first.")

      // Operators
      .def(py::self == py::self, "Return true if chord angles are equal")
      .def(py::self != py::self, "Return true if chord angles are not equal")
      .def(py::self < py::self,
           "Return true if this is less than other (by length2)")
      .def(py::self > py::self,
           "Return true if this is greater than other (by length2)")
      .def(py::self <= py::self,
           "Return true if this is less than or equal to other")
      .def(py::self >= py::self,
           "Return true if this is greater than or equal to other")
      .def("__add__", [](const S1ChordAngle& a, const S1ChordAngle& b) {
               MaybeThrowSpecialForArithmetic(a, "left operand");
               MaybeThrowSpecialForArithmetic(b, "right operand");
               return a + b;
           }, py::is_operator(),
           "Add two chord angles, clamping the result to [0, Pi].\n\n"
           "Raises ValueError if either operand is Negative() or Infinity().")
      .def("__sub__", [](const S1ChordAngle& a, const S1ChordAngle& b) {
               MaybeThrowSpecialForArithmetic(a, "left operand");
               MaybeThrowSpecialForArithmetic(b, "right operand");
               return a - b;
           }, py::is_operator(),
           "Subtract two chord angles, clamping the result to [0, Pi].\n\n"
           "Raises ValueError if either operand is Negative() or Infinity().")
      .def("__iadd__", [](S1ChordAngle& self, const S1ChordAngle& other) {
               MaybeThrowSpecialForArithmetic(self, "left operand");
               MaybeThrowSpecialForArithmetic(other, "right operand");
               self += other;
               return self;
           }, py::is_operator(), "In-place addition")
      .def("__isub__", [](S1ChordAngle& self, const S1ChordAngle& other) {
               MaybeThrowSpecialForArithmetic(self, "left operand");
               MaybeThrowSpecialForArithmetic(other, "right operand");
               self -= other;
               return self;
           }, py::is_operator(), "In-place subtraction")
      .def("__hash__", [](const S1ChordAngle& self) {
        return std::hash<double>()(self.length2());
      })

      // String representation
      .def("__repr__", [](const S1ChordAngle& a) {
        std::ostringstream oss;
        oss << "S1ChordAngle(" << a << ")";
        return oss.str();
      })
      .def("__str__", [](const S1ChordAngle& a) {
        std::ostringstream oss;
        oss << a;
        return oss.str();
      });
}
