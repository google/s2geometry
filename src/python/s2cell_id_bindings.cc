#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include <cstdint>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "absl/strings/str_cat.h"
#include "s2/s2cell_id.h"
#include "s2/s2cellid_range.h"
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

void MaybeThrowPositionOutOfRange(uint64_t pos) {
  if (pos > S2CellId::kMaxPosition) {
    throw py::value_error(
        absl::StrCat("pos ", pos, " out of range [0, ", S2CellId::kMaxPosition,
                     "]"));
  }
}

// Raise OverflowError if n exceeds Py_ssize_t; otherwise return it cast.
// On 32-bit Python, cells() ranges above level ~15 exceed ssize_t max.
Py_ssize_t LenOrThrow(int64_t n) {
  if (n > std::numeric_limits<Py_ssize_t>::max()) {
    throw std::overflow_error(absl::StrCat(
        "S2CellIdRange size ", n, " exceeds Py_ssize_t max (",
        std::numeric_limits<Py_ssize_t>::max(),
        ") on this platform. On 32-bit Python, cells() ranges above "
        "level ~15 exceed ssize_t. Use .size() for the full int64 count."));
  }
  return static_cast<Py_ssize_t>(n);
}

// Normalize a Python-style index (negative counts from end) against a range
// of length n and return the 0-based index.  Raises IndexError if the index
// is still out of [0, n) after normalization.
// https://docs.python.org/3/reference/datamodel.html#object.__getitem__
int64_t IndexOrThrow(int64_t i, int64_t n) {
  if (i < 0) i += n;
  if (i < 0 || i >= n) throw py::index_error("index out of range");
  return i;
}

// Raise ValueError if the slice step is anything other than 1 (or None).
// S2CellIdRange slices must be contiguous; reversed() covers step=-1.
void StepOrThrow(py::slice s) {
  py::object step_obj = s.attr("step");
  if (step_obj.is_none()) return;
  int64_t step = step_obj.cast<int64_t>();
  if (step != 1) {
    throw py::value_error(absl::StrCat(
        "S2CellIdRange only supports slices with step 1 (got step ",
        step, ")"));
  }
}

// Resolve a Python slice against a range of length n and return the clamped
// [start, stop) bounds.  Slice indices are silently clamped per:
// https://docs.python.org/3/reference/datamodel.html#sequences
std::pair<int64_t, int64_t> ClampedSlice(int64_t n, py::slice s) {
  auto clamp = [n](int64_t i) -> int64_t {
    if (i < 0) i += n;
    if (i < 0) return int64_t{0};
    if (i > n) return n;
    return i;
  };
  py::object start_obj = s.attr("start");
  py::object stop_obj  = s.attr("stop");
  int64_t start = start_obj.is_none() ? int64_t{0} : clamp(start_obj.cast<int64_t>());
  int64_t stop  = stop_obj.is_none()  ? n           : clamp(stop_obj.cast<int64_t>());
  if (stop < start) stop = start;
  return {start, stop};
}

// Dereference an optional S2CellId from an iterator's next(), or raise
// StopIteration if the iterator is exhausted.
S2CellId NextOrStop(std::optional<S2CellId> result) {
  if (!result) throw py::stop_iteration();
  return *result;
}

}  // namespace

void bind_s2cell_id(py::module& m) {
  auto cls = py::class_<S2CellId>(m, "S2CellId",
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

      // Factory methods
      .def_static("from_face", [](int face) {
               MaybeThrowFaceOutOfRange(face);
               return S2CellId::FromFace(face);
           }, py::arg("face"),
           "Return the cell corresponding to a given S2 cube face (0..5).\n\n"
           "Raises ValueError if face is out of range.")
      .def_static("from_face_pos_level", [](int face, uint64_t pos, int level) {
               MaybeThrowFaceOutOfRange(face);
               MaybeThrowPositionOutOfRange(pos);
               MaybeThrowLevelOutOfRange(level, 0, S2CellId::kMaxLevel);
               return S2CellId::FromFacePosLevel(face, pos, level);
           },
           py::arg("face"), py::arg("pos"), py::arg("level"),
           "Return a cell given its face, Hilbert curve position, and level.\n\n"
           "Raises ValueError if face, pos, or level is out of range.")
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
      .def("size_ij",
           py::overload_cast<>(&S2CellId::GetSizeIJ, py::const_),
           "Return the edge length of this cell in (i,j)-space")
      .def_static("size_ij_for_level", [](int level) {
               MaybeThrowLevelOutOfRange(level, 0, S2CellId::kMaxLevel);
               return S2CellId::GetSizeIJ(level);
           }, py::arg("level"),
           "Return the edge length in (i,j)-space of cells at the given level")
      .def("size_st",
           py::overload_cast<>(&S2CellId::GetSizeST, py::const_),
           "Return the edge length of this cell in (s,t)-space")
      .def_static("size_st_for_level", [](int level) {
               MaybeThrowLevelOutOfRange(level, 0, S2CellId::kMaxLevel);
               return S2CellId::GetSizeST(level);
           }, py::arg("level"),
           "Return the edge length in (s,t)-space of cells at the given level")
      .def("to_point", &S2CellId::ToPoint,
           "Return the center of the cell as a normalized S2Point")
      .def("to_lat_lng", &S2CellId::ToLatLng,
           "Return the S2LatLng corresponding to the center of the cell")
      .def("center_st", &S2CellId::GetCenterST,
           "Return the center of the cell in (s,t)-space")
      .def("bound_st", &S2CellId::GetBoundST,
           "Return the bounds of this cell in (s,t)-space")
      .def("center_uv", &S2CellId::GetCenterUV,
           "Return the center of the cell in (u,v)-space")
      .def("bound_uv", &S2CellId::GetBoundUV,
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
      .def("common_ancestor_level", &S2CellId::GetCommonAncestorLevel,
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
           "The range supports len(), indexing, slicing, iteration, and\n"
           "reversed().\n"
           "Raises ValueError if this is a leaf cell or level is out of range.")
      .def_static("cells", [](int level) {
               MaybeThrowLevelOutOfRange(level, 0, S2CellId::kMaxLevel);
               return S2CellIdRange{S2CellId::Begin(level),
                                    S2CellId::End(level)};
           }, py::arg("level"),
           "Return a range over all cells at the given level across all 6 faces.\n\n"
           "The range supports len(), indexing, slicing, iteration, and\n"
           "reversed().\n"
           "Warning: the number of cells grows as 6 * 4^level.")
      .def("edge_neighbors", [](S2CellId self) {
               S2CellId neighbors[4];
               self.GetEdgeNeighbors(neighbors);
               return py::make_tuple(neighbors[0], neighbors[1],
                                     neighbors[2], neighbors[3]);
           },
           "Return the four cells adjacent across this cell's edges")
      .def("vertex_neighbors", [](S2CellId self, int level) {
               MaybeThrowIfFace(self);
               MaybeThrowLevelOutOfRange(level, 0, self.level() - 1);
               std::vector<S2CellId> output;
               self.AppendVertexNeighbors(level, &output);
               return output;
           }, py::arg("level"),
           "Return the neighbors of the closest vertex at the given level.\n\n"
           "Normally returns 4 neighbors, but may return 3 for cube vertices.\n"
           "Raises ValueError if level >= self.level().")
      .def("all_neighbors", [](S2CellId self, int level) {
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

  cls.attr("MAX_LEVEL")    = S2CellId::kMaxLevel;
  cls.attr("MAX_POSITION") = S2CellId::kMaxPosition;
  cls.attr("NUM_FACES")    = S2CellId::kNumFaces;

  py::class_<S2CellIdRange>(m, "S2CellIdRange",
      "A range of S2CellIds at the same level along the Hilbert curve.\n\n"
      "Supports len(), indexing, slicing, iteration, reversed(), and 'in'.\n"
      "Returned from S2CellId.children() and S2CellId.cells().\n\n"
      "Slicing: only step=1 is supported. Use reversed() for step=-1.\n"
      "Other step values raise ValueError.")
      .def("__len__",      [](const S2CellIdRange& self) {
        return LenOrThrow(self.size());
      })
      .def("size", &S2CellIdRange::size,
           "Return the number of cells in the range as a 64-bit integer.\n\n"
           "Equivalent to len() on 64-bit platforms; use this method when\n"
           "the range may exceed Py_ssize_t (e.g. very large cells() ranges\n"
           "on 32-bit Python).")
      .def("__getitem__",  [](const S2CellIdRange& self, int64_t i) {
        return self.at(IndexOrThrow(i, self.size()));
      }, py::arg("index"))
      .def("__getitem__",  [](const S2CellIdRange& self, py::slice s) {
        StepOrThrow(s);
        auto [start, stop] = ClampedSlice(self.size(), s);
        return self.slice(start, stop);
      }, py::arg("slice"))
      .def("__contains__", &S2CellIdRange::contains, py::arg("cell"))
      .def("__iter__", [](const S2CellIdRange& self) {
        return S2CellIdForwardIterator{self.begin, self.end};
      })
      .def("__reversed__", [](const S2CellIdRange& self) {
        if (self.size() == 0) {
          return S2CellIdReverseIterator{self.begin, self.begin, true};
        }
        return S2CellIdReverseIterator{self.end.prev(), self.begin, false};
      });

  py::class_<S2CellIdForwardIterator>(m, "S2CellIdForwardIterator", py::module_local())
      .def("__iter__", [](S2CellIdForwardIterator& self) -> S2CellIdForwardIterator& {
        return self;
      })
      .def("__next__", [](S2CellIdForwardIterator& self) {
        return NextOrStop(self.next());
      });

  py::class_<S2CellIdReverseIterator>(m, "S2CellIdReverseIterator", py::module_local())
      .def("__iter__", [](S2CellIdReverseIterator& self) -> S2CellIdReverseIterator& {
        return self;
      })
      .def("__next__", [](S2CellIdReverseIterator& self) {
        return NextOrStop(self.next());
      });
}
