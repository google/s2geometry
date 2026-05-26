#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <sstream>
#include <vector>

#include "absl/strings/str_cat.h"
#include "s2/s1chord_angle.h"
#include "s2/s2cell.h"
#include "s2/s2cell_id.h"
#include "s2/s2latlng.h"
#include "s2/s2point.h"

namespace py = pybind11;

namespace {

void MaybeThrowFaceOutOfRange(int face) {
  if (face < 0 || face >= S2CellId::kNumFaces) {
    throw py::value_error(
        absl::StrCat("Face ", face, " out of range [0, ",
                     S2CellId::kNumFaces - 1, "]"));
  }
}

void MaybeThrowLevelOutOfRange(int level, int min, int max) {
  if (level < min || level > max) {
    throw py::value_error(
        absl::StrCat("Level ", level, " out of range [", min, ", ", max, "]"));
  }
}

void MaybeThrowPositionOutOfRange(uint64_t pos) {
  if (pos > S2CellId::kMaxPosition) {
    throw py::value_error(
        absl::StrCat("pos ", pos, " out of range [0, ", S2CellId::kMaxPosition,
                     "]"));
  }
}

}  // namespace

void bind_s2cell(py::module& m) {
  auto cls = py::class_<S2Cell>(m, "S2Cell",
      "An S2Region representing a single cell on the sphere.\n\n"
      "Unlike S2CellId (which is just a 64-bit identifier), S2Cell carries\n"
      "precomputed state that allows efficient containment and intersection\n"
      "tests. See s2/s2cell.h for comprehensive documentation.")
      // Constructors
      .def(py::init([](S2CellId id) {
               // No validity check needed: all S2CellId objects reachable from
               // Python are guaranteed valid by the S2CellId bindings.
               return S2Cell(id);
           }),
           py::arg("cell_id"),
           "Construct the cell corresponding to the given S2CellId.")
      .def(py::init([](const S2Point& p) {
               return S2Cell(p);
           }),
           py::arg("point"),
           "Construct a leaf cell containing the given point.\n\n"
           "The point does not need to be normalized.")
      .def(py::init([](const S2LatLng& ll) {
               return S2Cell(ll);
           }),
           py::arg("latlng"),
           "Construct a leaf cell containing the given S2LatLng.")

      // Factory methods
      .def_static("from_face", [](int face) {
               MaybeThrowFaceOutOfRange(face);
               return S2Cell::FromFace(face);
           }, py::arg("face"),
           "Return the cell corresponding to the given S2 cube face (0..5).\n\n"
           "Raises ValueError if face is out of range.")
      .def_static("from_face_pos_level", [](int face, uint64_t pos, int level) {
               MaybeThrowFaceOutOfRange(face);
               MaybeThrowPositionOutOfRange(pos);
               MaybeThrowLevelOutOfRange(level, 0, S2CellId::kMaxLevel);
               return S2Cell::FromFacePosLevel(face, pos, level);
           },
           py::arg("face"), py::arg("pos"), py::arg("level"),
           "Return a cell given its face, Hilbert curve position, and level.\n\n"
           "Raises ValueError if face, pos, or level is out of range.")

      // Properties
      .def_property_readonly("id", &S2Cell::id,
                             "The S2CellId this cell corresponds to")
      .def_property_readonly("face", &S2Cell::face,
                             "Which cube face this cell belongs to (0..5)")
      .def_property_readonly("level", &S2Cell::level,
                             "The subdivision level (0..kMaxLevel)")
      .def_property_readonly("orientation", &S2Cell::orientation,
                             "The Hilbert curve orientation of this cell.\n\n"
                             "A bitmask: bit 0 (kSwapMask) indicates swapped\n"
                             "axes; bit 1 (kInvertMask) indicates 180-degree\n"
                             "rotation. Values are in [0, 3].")

      // Predicates
      .def("is_leaf", &S2Cell::is_leaf,
           "Return true if this is a leaf cell (level == kMaxLevel)")

      // Geometric operations
      .def("get_size_ij", &S2Cell::GetSizeIJ,
           "Return the edge length of this cell in (i,j)-space")
      .def("get_size_st", &S2Cell::GetSizeST,
           "Return the edge length of this cell in (s,t)-space")
      .def("vertex", &S2Cell::GetVertex, py::arg("k"),
           "Return the k-th vertex of the cell (k = 0,1,2,3) in CCW order.\n\n"
           "Lower-left, lower-right, upper-right, upper-left in the UV plane.\n"
           "The argument is reduced modulo 4 to the range [0..3].\n"
           "The returned point is normalized.")
      .def("edge", &S2Cell::GetEdge, py::arg("k"),
           "Return the normalized inward-facing normal of the great circle\n"
           "passing through the edge from vertex k to vertex k+1 (mod 4).\n\n"
           "The argument is reduced modulo 4 to the range [0..3].")
      .def("uv_coord_of_edge", &S2Cell::GetUVCoordOfEdge, py::arg("k"),
           "Return either U or V for the given edge, whichever is constant\n"
           "along it.\n\n"
           "Boundaries 0 and 2 return V; boundaries 1 and 3 return U.\n"
           "The argument is reduced modulo 4 to the range [0..3].")
      .def("ij_coord_of_edge", &S2Cell::GetIJCoordOfEdge, py::arg("k"),
           "Return either I or J for the given edge, whichever is constant\n"
           "along it.\n\n"
           "Boundaries 0 and 2 return J; boundaries 1 and 3 return I.\n"
           "The argument is reduced modulo 4 to the range [0..3].")
      .def("center", &S2Cell::GetCenter,
           "Return the center of the cell as a normalized S2Point")
      .def_static("average_area_for_level", [](int level) {
               MaybeThrowLevelOutOfRange(level, 0, S2CellId::kMaxLevel);
               return S2Cell::AverageArea(level);
           }, py::arg("level"),
           "Return the average area of cells at the given level,\n"
           "in steradians.\n\n"
           "Raises ValueError if level is out of range.")
      .def("approx_area", &S2Cell::ApproxArea,
           "Return the approximate area of this cell in steradians.\n\n"
           "Accurate to within 3% for all cell sizes and within 0.1% for\n"
           "cells at level 5 or higher.")
      .def("exact_area", &S2Cell::ExactArea,
           "Return the area of this cell as accurately as possible,\n"
           "in steradians.\n\n"
           "More expensive than approx_area but accurate to 6 digits\n"
           "even for leaf cells.")
      .def("get_bound_uv", &S2Cell::GetBoundUV,
           "Return the bounds of this cell in (u,v)-space")
      .def("get_distance", py::overload_cast<const S2Point&>(
               &S2Cell::GetDistance, py::const_),
           py::arg("point"),
           "Return the distance from this cell to the given point.\n\n"
           "Returns zero if the point is inside the cell.")
      .def("get_boundary_distance", &S2Cell::GetBoundaryDistance,
           py::arg("point"),
           "Return the distance from the cell boundary to the given point.")
      .def("get_max_distance", py::overload_cast<const S2Point&>(
               &S2Cell::GetMaxDistance, py::const_),
           py::arg("point"),
           "Return the maximum distance from this cell to the given point.")
      .def("get_distance_to_edge",
           py::overload_cast<const S2Point&, const S2Point&>(
               &S2Cell::GetDistance, py::const_),
           py::arg("a"), py::arg("b"),
           "Return the minimum distance from this cell to the edge AB.\n\n"
           "Returns zero if the edge intersects the cell interior.")
      .def("get_max_distance_to_edge",
           py::overload_cast<const S2Point&, const S2Point&>(
               &S2Cell::GetMaxDistance, py::const_),
           py::arg("a"), py::arg("b"),
           "Return the maximum distance from this cell to the edge AB.")
      .def("get_distance_to_cell",
           py::overload_cast<const S2Cell&>(&S2Cell::GetDistance, py::const_),
           py::arg("cell"),
           "Return the distance from this cell to the given cell.\n\n"
           "Returns zero if one cell contains the other.")
      .def("get_max_distance_to_cell",
           py::overload_cast<const S2Cell&>(&S2Cell::GetMaxDistance, py::const_),
           py::arg("cell"),
           "Return the maximum distance from this cell to the given cell.")
      .def("get_cell_union_bound", [](const S2Cell& self) {
               std::vector<S2CellId> cell_ids;
               self.GetCellUnionBound(&cell_ids);
               return cell_ids;
           },
           "Return a list of S2CellIds whose union covers this cell.\n\n"
           "For a single S2Cell, this always returns a list containing\n"
           "just this cell's id.")
      .def("contains", py::overload_cast<const S2Cell&>(
               &S2Cell::Contains, py::const_),
           py::arg("cell"),
           "Return true if this cell contains the given cell")
      .def("contains_point", py::overload_cast<const S2Point&>(
               &S2Cell::Contains, py::const_),
           py::arg("point"),
           "Return true if this cell contains the given point.\n\n"
           "S2Cells are closed sets: points along an edge or vertex\n"
           "belong to the adjacent cell(s) as well.\n"
           "The point does not need to be normalized.")
      .def("may_intersect", &S2Cell::MayIntersect, py::arg("cell"),
           "Return true if this cell may intersect the given cell")

      // Traversal
      .def("subdivide", [](const S2Cell& self) {
               if (self.is_leaf()) {
                 throw py::value_error("Leaf cell has no children");
               }
               S2Cell children[4];
               self.Subdivide(children);
               return py::make_tuple(children[0], children[1],
                                     children[2], children[3]);
           },
           "Return the four children of this cell as a tuple.\n\n"
           "Raises ValueError if this is a leaf cell.")

      // Operators
      .def(py::self == py::self, "Return true if cells are equal")
      .def(py::self != py::self, "Return true if cells are not equal")
      .def(py::self < py::self, "Compare cells by their cell id")
      .def(py::self > py::self, "Compare cells by their cell id")
      .def(py::self <= py::self, "Compare cells by their cell id")
      .def(py::self >= py::self, "Compare cells by their cell id")
      .def("__hash__", [](const S2Cell& self) {
        return absl::Hash<S2CellId>()(self.id());
      })

      // String representation
      .def("__repr__", [](const S2Cell& cell) {
        std::ostringstream oss;
        oss << "S2Cell(" << cell.id() << ")";
        return oss.str();
      })
      .def("__str__", [](const S2Cell& cell) {
        std::ostringstream oss;
        oss << cell.id();
        return oss.str();
      });

  // Edge boundary constants (from S2Cell::Boundary enum).
  cls.attr("BOTTOM_EDGE") = static_cast<int>(S2Cell::kBottomEdge);
  cls.attr("RIGHT_EDGE")  = static_cast<int>(S2Cell::kRightEdge);
  cls.attr("TOP_EDGE")    = static_cast<int>(S2Cell::kTopEdge);
  cls.attr("LEFT_EDGE")   = static_cast<int>(S2Cell::kLeftEdge);

  // TODO: The following S2Cell methods are not yet bound because they depend
  // on types that have not been bound yet:
  //   - get_cap_bound()   -> S2Cap
  //   - get_rect_bound()  -> S2LatLngRect
}
