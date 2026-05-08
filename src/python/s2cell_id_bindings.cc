#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <sstream>
#include <vector>

#include "absl/strings/str_cat.h"
#include "s2/s2cell_id.h"
#include "s2/s2latlng.h"

namespace py = pybind11;

namespace {

void MaybeThrowNotValid(S2CellId id) {
  if (!id.is_valid()) {
    throw py::value_error(absl::StrCat("Invalid S2CellId: ", id.ToString()));
  }
}

void MaybeThrowLevelOutOfRange(int level, int min, int max) {
  if (level < min || level > max) {
    throw py::value_error(
        absl::StrCat("Level ", level, " out of range [", min, ", ", max, "]"));
  }
}

void MaybeThrowFaceOutOfRange(int face) {
  if (face < 0 || face >= S2CellId::kNumFaces) {
    throw py::value_error(
        absl::StrCat("Face ", face, " out of range [0, ",
                     S2CellId::kNumFaces - 1, "]"));
  }
}

void MaybeThrowIfLeaf(S2CellId id) {
  if (id.is_leaf()) {
    throw py::value_error("Function invalid for leaf cells");
  }
}

void MaybeThrowIfFace(S2CellId id) {
  if (id.is_face()) {
    throw py::value_error("Function invalid for face cells");
  }
}

void MaybeThrowCellIdOutOfRange(const char* name, int value, int min, int max) {
  if (value < min || value > max) {
    throw py::value_error(
        absl::StrCat(name, " ", value, " out of range [", min, ", ", max, "]"));
  }
}

void MaybeThrowChildPositionOutOfRange(int position) {
  if (position < 0 || position > 3) {
    throw py::value_error(
        absl::StrCat("Child position ", position, " out of range [0, 3]"));
  }
}

}  // namespace

void bind_s2cell_id(py::module& m) {
  py::class_<S2CellId>(m, "S2CellId",
      "A 64-bit unsigned integer that uniquely identifies a cell in the\n"
      "S2 cell decomposition.\n\n"
      "See s2/s2cell_id.h for comprehensive documentation.")
      // Constructors
      .def(py::init([](uint64_t id) {
               S2CellId cell(id);
               MaybeThrowNotValid(cell);
               return cell;
           }),
           py::arg("id"),
           "Construct from a 64-bit cell id.\n\n"
           "Raises ValueError if the id is not valid.")
      .def(py::init([](const S2Point& p) {
               return S2CellId(p);
           }),
           py::arg("point"),
           "Construct a leaf cell containing the given point.\n\n"
           "The point does not need to be normalized.")
      .def(py::init([](const S2LatLng& ll) {
               return S2CellId(ll);
           }),
           py::arg("latlng"),
           "Construct a leaf cell containing the given S2LatLng.")

      // Constants
      .def_readonly_static("MAX_LEVEL", &S2CellId::kMaxLevel,
           "Maximum cell subdivision level (30)")
      .def_readonly_static("NUM_FACES", &S2CellId::kNumFaces,
           "Number of cube faces (6)")

      // Factory methods
      .def_static("from_face", [](int face) {
               MaybeThrowFaceOutOfRange(face);
               return S2CellId::FromFace(face);
           }, py::arg("face"),
           "Return the cell corresponding to a given S2 cube face (0..5).\n\n"
           "Raises ValueError if face is out of range.")
      .def_static("from_face_pos_level", [](int face, uint64_t pos, int level) {
               MaybeThrowFaceOutOfRange(face);
               MaybeThrowLevelOutOfRange(level, 0, S2CellId::kMaxLevel);
               return S2CellId::FromFacePosLevel(face, pos, level);
           },
           py::arg("face"), py::arg("pos"), py::arg("level"),
           "Return a cell given its face, Hilbert curve position, and level.\n\n"
           "Raises ValueError if face or level is out of range.")
      .def_static("from_token", [](const std::string& token) {
               S2CellId cell = S2CellId::FromToken(token);
               if (cell == S2CellId::None()) {
                 throw py::value_error(
                     absl::StrCat("Invalid S2CellId token: '", token, "'"));
               }
               return cell;
           }, py::arg("token"),
           "Return a cell from its token string.\n\n"
           "Raises ValueError if the token is malformed.")
      .def_static("from_face_ij", [](int face, int i, int j) {
               MaybeThrowFaceOutOfRange(face);
               MaybeThrowCellIdOutOfRange("i", i, 0, S2CellId::kMaxSize - 1);
               MaybeThrowCellIdOutOfRange("j", j, 0, S2CellId::kMaxSize - 1);
               return S2CellId::FromFaceIJ(face, i, j);
           },
           py::arg("face"), py::arg("i"), py::arg("j"),
           "Return a leaf cell given its cube face and (i, j) coordinates")

      // Properties
      .def_property_readonly("id", &S2CellId::id,
                             "The 64-bit unique identifier for this cell")

      // Predicates
      .def("is_leaf", &S2CellId::is_leaf,
           "Return true if this is a leaf cell (level == kMaxLevel)")
      .def("is_face", &S2CellId::is_face,
           "Return true if this is a top-level face cell (level == 0)")

      // Geometric operations
      .def_property_readonly("face", &S2CellId::face,
           "Which cube face this cell belongs to (0..5)")
      .def_property_readonly("pos", &S2CellId::pos,
           "The position along the Hilbert curve over this face")
      .def_property_readonly("level", &S2CellId::level,
           "The subdivision level (0..kMaxLevel)")
      .def("get_size_ij",
           py::overload_cast<>(&S2CellId::GetSizeIJ, py::const_),
           "Return the edge length of this cell in (i,j)-space")
      .def_static("get_size_ij_for_level", [](int level) {
               MaybeThrowLevelOutOfRange(level, 0, S2CellId::kMaxLevel);
               return S2CellId::GetSizeIJ(level);
           }, py::arg("level"),
           "Return the edge length in (i,j)-space of cells at the given level")
      .def("get_size_st",
           py::overload_cast<>(&S2CellId::GetSizeST, py::const_),
           "Return the edge length of this cell in (s,t)-space")
      .def_static("get_size_st_for_level", [](int level) {
               MaybeThrowLevelOutOfRange(level, 0, S2CellId::kMaxLevel);
               return S2CellId::GetSizeST(level);
           }, py::arg("level"),
           "Return the edge length in (s,t)-space of cells at the given level")
      .def("to_point", &S2CellId::ToPoint,
           "Return the center of the cell as a normalized S2Point")
      .def("to_lat_lng", &S2CellId::ToLatLng,
           "Return the S2LatLng corresponding to the center of the cell")
      .def("get_center_st", &S2CellId::GetCenterST,
           "Return the center of the cell in (s,t)-space")
      .def("get_bound_st", &S2CellId::GetBoundST,
           "Return the bounds of this cell in (s,t)-space")
      .def("get_center_uv", &S2CellId::GetCenterUV,
           "Return the center of the cell in (u,v)-space")
      .def("get_bound_uv", &S2CellId::GetBoundUV,
           "Return the bounds of this cell in (u,v)-space")
      .def("child_position", [](S2CellId self) {
               MaybeThrowLevelOutOfRange(self.level(), 1, S2CellId::kMaxLevel);
               return self.child_position();
           },
           "Return the child position (0..3) within this cell's parent.\n\n"
           "Raises ValueError if level < 1.")
      .def("child_position_at_level", [](S2CellId self, int level) {
               MaybeThrowLevelOutOfRange(level, 1, self.level());
               return self.child_position(level);
           }, py::arg("level"),
           "Return the child position (0..3) of this cell's ancestor at\n"
           "the given level within its parent.\n\n"
           "Raises ValueError if level is out of range [1, self.level()].")
      .def("to_token", &S2CellId::ToToken,
           "Return a compact token string for this cell.\n\n"
           "Tokens preserve ordering and the round-trip\n"
           "from_token(cell.to_token()) == cell is guaranteed.")
      .def("range_min", &S2CellId::range_min,
           "Return the minimum cell id contained within this cell")
      .def("range_max", &S2CellId::range_max,
           "Return the maximum cell id contained within this cell")
      .def("contains", &S2CellId::contains, py::arg("other"),
           "Return true if the given cell is contained within this one")
      .def("intersects", &S2CellId::intersects, py::arg("other"),
           "Return true if the given cell intersects this one")
      .def("get_common_ancestor_level", &S2CellId::GetCommonAncestorLevel,
           py::arg("other"),
           "Return the level of the lowest common ancestor.\n\n"
           "Returns -1 if the cells are from different faces.")

      // Traversal
      .def("parent", [](S2CellId self) {
               MaybeThrowIfFace(self);
               return self.parent();
           },
           "Return the cell at the previous level.\n\n"
           "Raises ValueError if this is already a face cell.")
      .def("parent_at_level", [](S2CellId self, int level) {
               MaybeThrowLevelOutOfRange(level, 0, self.level());
               return self.parent(level);
           }, py::arg("level"),
           "Return the cell at the given level.\n\n"
           "Raises ValueError if level is out of range [0, self.level()].")
      .def("child", [](S2CellId self, int position) {
               MaybeThrowIfLeaf(self);
               MaybeThrowChildPositionOutOfRange(position);
               return self.child(position);
           }, py::arg("position"),
           "Return the immediate child at the given position (0..3).\n\n"
           "Raises ValueError if this is a leaf cell or position is out of range.")
      .def("get_edge_neighbors", [](S2CellId self) {
               S2CellId neighbors[4];
               self.GetEdgeNeighbors(neighbors);
               return py::make_tuple(neighbors[0], neighbors[1],
                                     neighbors[2], neighbors[3]);
           },
           "Return the four cells adjacent across this cell's edges")
      .def("get_vertex_neighbors", [](S2CellId self, int level) {
               MaybeThrowIfFace(self);
               MaybeThrowLevelOutOfRange(level, 0, self.level() - 1);
               std::vector<S2CellId> output;
               self.AppendVertexNeighbors(level, &output);
               return output;
           }, py::arg("level"),
           "Return the neighbors of the closest vertex at the given level.\n\n"
           "Normally returns 4 neighbors, but may return 3 for cube vertices.\n"
           "Raises ValueError if level >= self.level().")
      .def("get_all_neighbors", [](S2CellId self, int level) {
               MaybeThrowLevelOutOfRange(level, self.level(),
                                         S2CellId::kMaxLevel);
               std::vector<S2CellId> output;
               self.AppendAllNeighbors(level, &output);
               return output;
           }, py::arg("level"),
           "Return all neighbors of this cell at the given level.\n\n"
           "Raises ValueError if level < self.level().")

      // Operators
      .def(py::self == py::self, "Return true if cell ids are equal")
      .def(py::self != py::self, "Return true if cell ids are not equal")
      .def(py::self < py::self, "Compare cell ids by their numeric value")
      .def(py::self > py::self, "Compare cell ids by their numeric value")
      .def(py::self <= py::self, "Compare cell ids by their numeric value")
      .def(py::self >= py::self, "Compare cell ids by their numeric value")
      .def("__hash__", [](S2CellId self) {
        return absl::Hash<S2CellId>()(self);
      })

      // String representation
      .def("__repr__", [](S2CellId id) {
        std::ostringstream oss;
        oss << "S2CellId(" << id << ")";
        return oss.str();
      })
      .def("__str__", [](S2CellId id) {
        std::ostringstream oss;
        oss << id;
        return oss.str();
      });
}
