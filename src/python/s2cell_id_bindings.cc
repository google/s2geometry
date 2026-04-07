#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include "s2/s2cell_id.h"
#include "s2/s2latlng.h"

namespace py = pybind11;

namespace {

void MaybeThrowNotValid(S2CellId id) {
  if (!id.is_valid()) {
    throw py::value_error("Invalid S2CellId: " + id.ToString());
  }
}

void MaybeThrowLevelOutOfRange(int level, int min, int max) {
  if (level < min || level > max) {
    throw py::value_error("Level " + std::to_string(level) +
                          " out of range [" + std::to_string(min) +
                          ", " + std::to_string(max) + "]");
  }
}

void MaybeThrowFaceOutOfRange(int face) {
  if (face < 0 || face >= S2CellId::kNumFaces) {
    throw py::value_error("Face " + std::to_string(face) +
                          " out of range [0, " +
                          std::to_string(S2CellId::kNumFaces - 1) + "]");
  }
}

void MaybeThrowIfLeaf(S2CellId id) {
  if (id.is_leaf()) {
    throw py::value_error("Leaf cell has no children");
  }
}

void MaybeThrowIfFace(S2CellId id) {
  if (id.is_face()) {
    throw py::value_error("Face cell has no parent");
  }
}

void MaybeThrowChildPositionOutOfRange(int position) {
  if (position < 0 || position > 3) {
    throw py::value_error("Child position " + std::to_string(position) +
                          " out of range [0, 3]");
  }
}

// A range of S2CellIds at the same level, supporting len, indexing, iteration,
// and reverse iteration. The range is [begin, end) where end is exclusive.
struct S2CellIdRange {
  S2CellId begin;
  S2CellId end;

  int64_t size() const {
    return end.distance_from_begin() - begin.distance_from_begin();
  }

  S2CellId item(int64_t i) const {
    int64_t len = size();
    if (i < 0) i += len;
    if (i < 0 || i >= len) {
      throw py::index_error("index " + std::to_string(i) + " out of range");
    }
    return begin.advance(i);
  }

  bool contains(S2CellId cell) const {
    return cell >= begin && cell < end &&
           cell.level() == begin.level();
  }
};

// Forward iterator for S2CellIdRange.
struct S2CellIdForwardIter {
  S2CellId cur;
  S2CellId end;

  S2CellId next() {
    if (cur == end) throw py::stop_iteration();
    S2CellId result = cur;
    cur = cur.next();
    return result;
  }
};

// Reverse iterator for S2CellIdRange.
struct S2CellIdReverseIter {
  S2CellId cur;
  S2CellId begin;
  bool done;

  S2CellId next() {
    if (done) throw py::stop_iteration();
    S2CellId result = cur;
    if (cur == begin) {
      done = true;
    } else {
      cur = cur.prev();
    }
    return result;
  }
};

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
      .def_property_readonly_static("kMaxLevel",
           [](py::object) { return S2CellId::kMaxLevel; },
           "Maximum cell subdivision level (30)")
      .def_property_readonly_static("kNumFaces",
           [](py::object) { return S2CellId::kNumFaces; },
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
                 throw py::value_error("Invalid S2CellId token: '" + token + "'");
               }
               return cell;
           }, py::arg("token"),
           "Return a cell from its token string.\n\n"
           "Raises ValueError if the token is malformed.")
      .def_static("from_debug_string", [](const std::string& str) {
               S2CellId cell = S2CellId::FromDebugString(str);
               if (cell == S2CellId::None()) {
                 throw py::value_error(
                     "Invalid S2CellId debug string: '" + str + "'");
               }
               return cell;
           }, py::arg("str"),
           "Return a cell from its debug string (e.g. \"3/02\").\n\n"
           "Raises ValueError if the string is malformed.")
      .def_static("from_face_ij", [](int face, int i, int j) {
               MaybeThrowFaceOutOfRange(face);
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
      .def("face", &S2CellId::face,
           "Return which cube face this cell belongs to (0..5)")
      .def("pos", &S2CellId::pos,
           "Return the position along the Hilbert curve over this face")
      .def("level", &S2CellId::level,
           "Return the subdivision level (0..kMaxLevel)")
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
      .def("to_point_raw", &S2CellId::ToPointRaw,
           "Return the center of the cell as an S2Point (not normalized)")
      .def("to_lat_lng", &S2CellId::ToLatLng,
           "Return the S2LatLng corresponding to the center of the cell")
      .def("get_center_st", &S2CellId::GetCenterST,
           "Return the center of the cell in (s,t) coordinates")
      .def("get_bound_st", &S2CellId::GetBoundST,
           "Return the bounds of this cell in (s,t)-space")
      .def("get_center_uv", &S2CellId::GetCenterUV,
           "Return the center of the cell in (u,v) coordinates")
      .def("get_bound_uv", &S2CellId::GetBoundUV,
           "Return the bounds of this cell in (u,v)-space")
      .def_static("expanded_by_distance_uv", &S2CellId::ExpandedByDistanceUV,
           py::arg("uv"), py::arg("distance"),
           "Expand a (u,v) rectangle to contain all points within the\n"
           "given distance of the boundary")
      .def("get_center_si_ti", [](S2CellId self) {
               int si, ti;
               int face = self.GetCenterSiTi(&si, &ti);
               return py::make_tuple(face, si, ti);
           },
           "Return the (face, si, ti) coordinates of the cell center")
      .def("to_face_ij_orientation", [](S2CellId self) {
               int i, j, orientation;
               int face = self.ToFaceIJOrientation(&i, &j, &orientation);
               return py::make_tuple(face, i, j, orientation);
           },
           "Return (face, i, j, orientation) for this cell")
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
      .def("to_string", &S2CellId::ToString,
           "Return the debug string (e.g. \"3/02\")")
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
      .def("children", [](S2CellId self, py::object level_obj) {
               MaybeThrowIfLeaf(self);
               S2CellId begin, end;
               if (level_obj.is_none()) {
                 begin = self.child_begin();
                 end = self.child_end();
               } else {
                 int level = level_obj.cast<int>();
                 MaybeThrowLevelOutOfRange(level, self.level(),
                                           S2CellId::kMaxLevel);
                 begin = self.child_begin(level);
                 end = self.child_end(level);
               }
               return S2CellIdRange{begin, end};
           }, py::arg("level") = py::none(),
           "Return a range over the children of this cell.\n\n"
           "With no argument, returns the 4 immediate children.\n"
           "With a level argument, returns all descendants at that level.\n"
           "The range supports len(), indexing, iteration, and reversed().\n"
           "Raises ValueError if this is a leaf cell or level is out of range.")
      .def_static("cells", [](int level) {
               MaybeThrowLevelOutOfRange(level, 0, S2CellId::kMaxLevel);
               return S2CellIdRange{S2CellId::Begin(level),
                                    S2CellId::End(level)};
           }, py::arg("level"),
           "Return a range over all cells at the given level across all 6 faces.\n\n"
           "The range supports len(), indexing, iteration, and reversed().\n"
           "Warning: the number of cells grows as 6 * 4^level.")
      .def("__iter__", [](S2CellId self) {
               return S2CellIdForwardIter{self, S2CellId::End(self.level())};
           },
           "Iterate along the Hilbert curve at this cell's level,\n"
           "starting from this cell to the end of the level.")
      .def("__reversed__", [](S2CellId self) {
               return S2CellIdReverseIter{self, S2CellId::Begin(self.level()),
                                          false};
           },
           "Iterate in reverse along the Hilbert curve at this cell's level,\n"
           "starting from this cell back to the beginning of the level.")
      .def("get_edge_neighbors", [](S2CellId self) {
               S2CellId neighbors[4];
               self.GetEdgeNeighbors(neighbors);
               return py::make_tuple(neighbors[0], neighbors[1],
                                     neighbors[2], neighbors[3]);
           },
           "Return the four cells adjacent across this cell's edges")
      .def("get_vertex_neighbors", [](S2CellId self, int level) {
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

  py::class_<S2CellIdRange>(m, "S2CellIdRange")
      .def("__len__", &S2CellIdRange::size)
      .def("__getitem__", &S2CellIdRange::item, py::arg("index"))
      .def("__contains__", &S2CellIdRange::contains, py::arg("cell"))
      .def("__iter__", [](const S2CellIdRange& self) {
        return S2CellIdForwardIter{self.begin, self.end};
      })
      .def("__reversed__", [](const S2CellIdRange& self) {
        if (self.size() == 0) {
          return S2CellIdReverseIter{self.begin, self.begin, true};
        }
        return S2CellIdReverseIter{self.end.prev(), self.begin, false};
      });

  py::class_<S2CellIdForwardIter>(m, "S2CellIdForwardIter")
      .def("__iter__", [](S2CellIdForwardIter& self) -> S2CellIdForwardIter& {
        return self;
      })
      .def("__next__", &S2CellIdForwardIter::next);

  py::class_<S2CellIdReverseIter>(m, "S2CellIdReverseIter")
      .def("__iter__", [](S2CellIdReverseIter& self) -> S2CellIdReverseIter& {
        return self;
      })
      .def("__next__", &S2CellIdReverseIter::next);
}
