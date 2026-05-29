#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include <cmath>
#include <sstream>

#include "absl/strings/str_cat.h"
#include "s2/s1angle.h"
#include "s2/s2latlng.h"
#include "s2/s2point.h"

namespace py = pybind11;

namespace {

void MaybeThrowNotValid(const S2LatLng& ll) {
  if (!ll.is_valid()) {
    throw py::value_error(absl::StrCat(
        "Invalid S2LatLng: (", ll.lat().degrees(), ", ", ll.lng().degrees(),
        ") (latitude must be in [-90, 90], longitude in [-180, 180])"));
  }
}

void MaybeThrowNotFinite(double m) {
  if (!std::isfinite(m)) {
    throw py::value_error(absl::StrCat("Scalar must be finite, got ", m));
  }
}

}  // namespace

void bind_s2latlng(py::module& m) {
  py::class_<S2LatLng>(m, "S2LatLng",
      "Represents a point on the unit sphere as a latitude-longitude pair.\n\n"
      "All S2LatLng values are guaranteed to be valid: latitude in [-90, 90]\n"
      "and longitude in [-180, 180] degrees. Constructors and factory methods\n"
      "raise ValueError for out-of-range inputs. Arithmetic operators\n"
      "automatically normalize the result.\n\n"
      "Note: the native C++ implementation does not normalize arithmetic\n"
      "results and allows invalid states.\n\n"
      "See s2/s2latlng.h for comprehensive documentation.")
      // Constructors
      .def(py::init<>(), "Default constructor sets latitude and longitude to zero")
      .def(py::init([](S1Angle lat, S1Angle lng) {
               S2LatLng ll(lat, lng);
               MaybeThrowNotValid(ll);
               return ll;
           }),
           py::arg("lat"), py::arg("lng"),
           "Construct from latitude and longitude angles.\n\n"
           "Raises ValueError if latitude is not in [-90, 90] or\n"
           "longitude is not in [-180, 180] degrees.")
      .def(py::init([](const S2Point& p) {
               S2LatLng ll(p);
               MaybeThrowNotValid(ll);
               return ll;
           }),
           py::arg("point"),
           "Convert a direction vector to an S2LatLng.\n\n"
           "The vector does not need to be unit length.\n"
           "Raises ValueError if the point has non-finite coordinates.")

      // Factory methods
      .def_static("from_radians", [](double lat_radians, double lng_radians) {
               S2LatLng ll = S2LatLng::FromRadians(lat_radians, lng_radians);
               MaybeThrowNotValid(ll);
               return ll;
           },
           py::arg("lat_radians"), py::arg("lng_radians"),
           "Construct from latitude and longitude in radians.\n\n"
           "Raises ValueError if out of range.")
      .def_static("normalized_from_radians",
           [](double lat_radians, double lng_radians) {
               MaybeThrowNotFinite(lat_radians);
               MaybeThrowNotFinite(lng_radians);
               return S2LatLng::FromRadians(lat_radians, lng_radians)
                   .Normalized();
           },
           py::arg("lat_radians"), py::arg("lng_radians"),
           "Construct from latitude and longitude in radians, accepting any\n"
           "finite values. Latitude is clamped to [-Pi/2, Pi/2] and longitude\n"
           "is wrapped to [-Pi, Pi].\n\n"
           "Raises ValueError if either value is non-finite.")
      .def_static("from_degrees", [](double lat_degrees, double lng_degrees) {
               S2LatLng ll = S2LatLng::FromDegrees(lat_degrees, lng_degrees);
               MaybeThrowNotValid(ll);
               return ll;
           },
           py::arg("lat_degrees"), py::arg("lng_degrees"),
           "Construct from latitude and longitude in degrees.\n\n"
           "Raises ValueError if out of range.")
      .def_static("from_e5", [](int32_t lat_e5, int32_t lng_e5) {
               S2LatLng ll = S2LatLng::FromE5(lat_e5, lng_e5);
               MaybeThrowNotValid(ll);
               return ll;
           },
           py::arg("lat_e5"), py::arg("lng_e5"),
           "Construct from latitude and longitude in E5 format.\n\n"
           "Raises ValueError if out of range.")
      .def_static("from_e6", [](int32_t lat_e6, int32_t lng_e6) {
               S2LatLng ll = S2LatLng::FromE6(lat_e6, lng_e6);
               MaybeThrowNotValid(ll);
               return ll;
           },
           py::arg("lat_e6"), py::arg("lng_e6"),
           "Construct from latitude and longitude in E6 format.\n\n"
           "Raises ValueError if out of range.")
      .def_static("from_e7", [](int32_t lat_e7, int32_t lng_e7) {
               S2LatLng ll = S2LatLng::FromE7(lat_e7, lng_e7);
               MaybeThrowNotValid(ll);
               return ll;
           },
           py::arg("lat_e7"), py::arg("lng_e7"),
           "Construct from latitude and longitude in E7 format.\n\n"
           "Raises ValueError if out of range.")
      .def_static("latitude", &S2LatLng::Latitude, py::arg("point"),
                  "Return the latitude of a direction vector")
      .def_static("longitude", &S2LatLng::Longitude, py::arg("point"),
                  "Return the longitude of a direction vector")

      // Properties
      .def_property_readonly("lat", &S2LatLng::lat, "The latitude as an S1Angle")
      .def_property_readonly("lng", &S2LatLng::lng,
                             "The longitude as an S1Angle")
      .def_property_readonly("coords", &S2LatLng::coords,
                             "The (lat, lng) coordinates in radians as an R2Point")

      // Geometric operations
      .def("to_point", &S2LatLng::ToPoint,
           "Convert to the equivalent unit-length S2Point")
      .def("get_distance", &S2LatLng::GetDistance,
           py::arg("other"),
           "Return the surface distance to another S2LatLng.\n\n"
           "Uses the Haversine formula.")
      // Note: default value must match C++ signature in s2latlng.h
      .def("approx_equals", [](const S2LatLng& self, const S2LatLng& o,
                                S1Angle max_error) {
               return self.ApproxEquals(o, max_error);
           },
           py::arg("other"),
           py::arg("max_error") = S1Angle::Radians(1e-15),
           "Return true if approximately equal to 'other'.\n\n"
           "Note: this compares coordinates in rectangular lat/lng space.\n"
           "For points near the poles, consider using get_distance() instead.")

      // Operators
      .def(py::self == py::self, "Return true if exactly equal")
      .def(py::self != py::self, "Return true if not exactly equal")
      .def(py::self < py::self,
           "Compare by latitude first, then longitude")
      .def(py::self > py::self,
           "Compare by latitude first, then longitude")
      .def(py::self <= py::self,
           "Compare by latitude first, then longitude")
      .def(py::self >= py::self,
           "Compare by latitude first, then longitude")
      .def("__add__", [](const S2LatLng& a, const S2LatLng& b) {
        return (a + b).Normalized();
      }, py::arg("other"),
           "Add two S2LatLngs.\n\n"
           "The result is automatically normalized.\n"
           "Note: the native C++ implementation does not normalize.")
      .def("__sub__", [](const S2LatLng& a, const S2LatLng& b) {
        return (a - b).Normalized();
      }, py::arg("other"),
           "Subtract two S2LatLngs.\n\n"
           "The result is automatically normalized.\n"
           "Note: the native C++ implementation does not normalize.")
      .def("__mul__", [](const S2LatLng& a, double m) {
        MaybeThrowNotFinite(m);
        return (a * m).Normalized();
      }, py::arg("m"),
           "Multiply by scalar.\n\n"
           "The result is automatically normalized.\n"
           "Note: the native C++ implementation does not normalize.")
      .def("__rmul__", [](const S2LatLng& a, double m) {
        MaybeThrowNotFinite(m);
        return (m * a).Normalized();
      }, py::arg("m"),
           "Multiply by scalar (reversed operands).\n\n"
           "The result is automatically normalized.\n"
           "Note: the native C++ implementation does not normalize.")

      // String representation
      .def("to_string_in_degrees", &S2LatLng::ToStringInDegrees,
           "Return \"lat,lng\" in degrees, e.g. \"37.794000,-122.395000\"")
      .def("__repr__", [](const S2LatLng& ll) {
        std::ostringstream oss;
        oss << "S2LatLng(" << ll << ")";
        return oss.str();
      })
      .def("__str__", [](const S2LatLng& ll) {
        std::ostringstream oss;
        oss << ll;
        return oss.str();
      });
}
