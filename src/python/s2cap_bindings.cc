#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <sstream>
#include <vector>

#include "absl/hash/hash.h"

#include "s2/s1angle.h"
#include "s2/s1chord_angle.h"
#include "s2/s2cap.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2point.h"

namespace py = pybind11;

void bind_s2cap(py::module& m) {
  py::class_<S2Cap>(m, "S2Cap",
      "A disc-shaped region defined by a center and radius on the sphere.\n\n"
      "The cap represents a portion of the unit sphere cut off by a plane.\n"
      "The boundary is a circle; the cap is a closed set (contains its\n"
      "boundary). A cap of radius Pi/2 is a hemisphere; Pi covers the full\n"
      "sphere. See s2/s2cap.h for comprehensive documentation.")

      // Constructors
      .def(py::init<>(),
           "Construct an empty S2Cap.")
      .def(py::init([](const S2Point& center, S1Angle radius) {
               return S2Cap(center.Normalize(), radius);
           }),
           py::arg("center"), py::arg("radius"),
           "Construct a cap with the given center and radius.\n\n"
           "center is normalized if not already unit length. Negative radius\n"
           "produces the empty cap; radius >= 180 degrees produces the full cap.")
      .def(py::init([](const S2Point& center, S1ChordAngle radius) {
               return S2Cap(center.Normalize(), radius);
           }),
           py::arg("center"), py::arg("radius"),
           "Construct a cap with the given center and radius as S1ChordAngle.\n\n"
           "center is normalized if not already unit length.")

      // Factory methods
      .def_static("from_point", [](const S2Point& center) {
               return S2Cap::FromPoint(center.Normalize());
           }, py::arg("center"),
           "Return a cap containing a single point.\n\n"
           "center is normalized if not already unit length.")
      .def_static("from_center_height", [](const S2Point& center, double height) {
               return S2Cap::FromCenterHeight(center.Normalize(), height);
           }, py::arg("center"), py::arg("height"),
           "Return a cap with the given center and height.\n\n"
           "Height is the distance from the center point to the cutoff plane.\n"
           "center is normalized if not already unit length.\n"
           "A negative height yields an empty cap; height >= 2 yields a full cap.")
      .def_static("from_center_area", [](const S2Point& center, double area) {
               return S2Cap::FromCenterArea(center.Normalize(), area);
           }, py::arg("center"), py::arg("area"),
           "Return a cap with the given center and surface area in steradians.\n\n"
           "The area also equals the solid angle subtended by the cap.\n"
           "center is normalized if not already unit length.\n"
           "A negative area yields an empty cap; area >= 4*Pi yields a full cap.")
      .def_static("empty", &S2Cap::Empty,
           "Return an empty cap (contains no points).")
      .def_static("full", &S2Cap::Full,
           "Return a full cap (contains all points).")
      .def_static("from_points", [](const std::vector<S2Point>& points) {
               if (points.empty()) return S2Cap::Empty();
               S2Cap result = S2Cap::FromPoint(points[0]);
               for (size_t i = 1; i < points.size(); ++i) result.AddPoint(points[i]);
               return result;
           }, py::arg("points"),
           "Return the smallest cap containing all given points.\n\n"
           "Returns the empty cap if the list is empty.")
      .def_static("from_caps", [](const std::vector<S2Cap>& caps) {
               if (caps.empty()) return S2Cap::Empty();
               S2Cap result = caps[0];
               for (size_t i = 1; i < caps.size(); ++i) result.AddCap(caps[i]);
               return result;
           }, py::arg("caps"),
           "Return the smallest cap containing all given caps.\n\n"
           "Returns the empty cap if the list is empty.")

      // Properties
      .def_property_readonly("center", &S2Cap::center,
           "The center of the cap as a unit-length S2Point.")
      .def_property_readonly("radius", &S2Cap::radius,
           "The radius as an S1ChordAngle.")
      .def_property_readonly("height", &S2Cap::height,
           "The height of the cap (distance from center point to cutoff plane).")

      // Geometric operations
      .def("radius_angle", &S2Cap::GetRadius,
           "Return the radius as an S1Angle.\n\n"
           "Requires a trigonometric operation; may differ slightly from the\n"
           "value passed to the S1Angle constructor.")
      .def("area", &S2Cap::GetArea,
           "Return the surface area of the cap in steradians.")
      .def("centroid", &S2Cap::GetCentroid,
           "Return the true centroid of the cap multiplied by its surface area.\n\n"
           "The result lies on the ray from the origin through the cap's center.\n"
           "For zero-radius caps, always returns the origin (0, 0, 0).")

      // Predicates
      .def("is_empty", &S2Cap::is_empty,
           "Return true if the cap contains no points.")
      .def("is_full", &S2Cap::is_full,
           "Return true if the cap contains all points.")

      .def("complement", &S2Cap::Complement,
           "Return the complement of the interior of the cap.\n\n"
           "Same boundary as this cap but no shared interior points.\n"
           "Note: complement of a singleton equals complement of an empty cap.")
      .def("expanded", &S2Cap::Expanded, py::arg("distance"),
           "Return a cap containing all points within distance of this cap.\n\n"
           "Any expansion of an empty cap is still empty.")
      .def("union", &S2Cap::Union, py::arg("other"),
           "Return the smallest cap enclosing this cap and other.")

      // Containment / intersection
      .def("contains", py::overload_cast<const S2Cap&>(
               &S2Cap::Contains, py::const_),
           py::arg("other"),
           "Return true if this cap contains the given cap.")
      .def("contains_point", py::overload_cast<const S2Point&>(
               &S2Cap::Contains, py::const_),
           py::arg("point"),
           "Return true if this cap contains the given point.\n\n"
           "point should be unit length.")
      .def("contains_cell", py::overload_cast<const S2Cell&>(
               &S2Cap::Contains, py::const_),
           py::arg("cell"),
           "Return true if this cap contains the given cell.")
      .def("intersects", py::overload_cast<const S2Cap&>(
               &S2Cap::Intersects, py::const_),
           py::arg("other"),
           "Return true if this cap intersects the given cap.")
      .def("interior_intersects", &S2Cap::InteriorIntersects, py::arg("other"),
           "Return true if the interior of this cap intersects other.\n\n"
           "This relationship is not symmetric: only the interior of this cap\n"
           "is tested, not the interior of other.")
      .def("interior_contains_point", &S2Cap::InteriorContains,
           py::arg("point"),
           "Return true if the interior of this cap contains the given point.\n\n"
           "point should be unit length.")
      .def("may_intersect", &S2Cap::MayIntersect, py::arg("cell"),
           "Return true if this cap may intersect the given cell.")

      .def("cap_bound", &S2Cap::GetCapBound,
           "Return a bounding cap for this cap (returns self).")
      // get_rect_bound() is deferred until S2LatLngRect is bound.
      .def("cell_union_bound", [](const S2Cap& self) {
               std::vector<S2CellId> cell_ids;
               self.GetCellUnionBound(&cell_ids);
               return cell_ids;
           },
           "Return a list of S2CellIds whose union covers this cap.")

      // Operators
      .def(py::self == py::self, "Return true if caps are identical.")
      .def(py::self != py::self, "Return true if caps are not identical.")
      .def("approx_equals", &S2Cap::ApproxEquals,
           py::arg("other"),
           py::arg("max_error") = S1Angle::Radians(1e-14),
           "Return true if this cap is approximately equal to other.\n\n"
           "Checks that the angle between centers and the difference in chord\n"
           "radii are both within max_error radians.")
      .def("__hash__", [](const S2Cap& self) {
           return absl::HashOf(self.center(), self.radius().length2());
      })

      // String representation
      .def("__repr__", [](const S2Cap& self) {
           std::ostringstream oss;
           oss << "S2Cap(" << self << ")";
           return oss.str();
      })
      .def("__str__", [](const S2Cap& self) {
           std::ostringstream oss;
           oss << self;
           return oss.str();
      });
}
