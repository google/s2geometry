// Copyright Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS-IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//

#include "s2/s2text_format.h"

#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "absl/container/inlined_vector.h"
#include "absl/log/absl_check.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/ascii.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_format.h"
#include "absl/strings/str_split.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "s2/mutable_s2shape_index.h"
#include "s2/s1angle.h"
#include "s2/s2cell_id.h"
#include "s2/s2cell_union.h"
#include "s2/s2debug.h"
#include "s2/s2latlng.h"
#include "s2/s2latlng_rect.h"
#include "s2/s2lax_polygon_shape.h"
#include "s2/s2lax_polyline_shape.h"
#include "s2/s2loop.h"
#include "s2/s2point.h"
#include "s2/s2point_vector_shape.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"
#include "s2/s2shape.h"
#include "s2/s2shape_index.h"
#include "s2/util/task/status_macros.h"

using absl::Span;
using absl::string_view;
using std::make_unique;
using std::string;
using std::unique_ptr;
using std::vector;

namespace s2textformat {

namespace {
static vector<string_view> SplitString(string_view str, char separator) {
  vector<string_view> result =
      absl::StrSplit(str, separator, absl::SkipWhitespace());
  for (string_view& e : result) {
    e = absl::StripAsciiWhitespace(e);
  }
  return result;
}
}  // namespace

absl::StatusOr<vector<S2LatLng>> ParseLatLngs(string_view str) {
  vector<S2LatLng> latlngs;
  for (const string_view lat_lng_str :
       absl::StrSplit(str, ',', absl::SkipEmpty())) {
    const absl::InlinedVector<string_view, 2> lat_lng =
        absl::StrSplit(lat_lng_str, ':');
    if (lat_lng.size() != 2) {
      return absl::InvalidArgumentError(absl::StrCat(
          "Invalid latlng format (expected lat:lng): ", lat_lng_str));
    }
    double lat, lng;
    if (!absl::SimpleAtod(lat_lng[0], &lat)) {
      return absl::InvalidArgumentError(
          absl::StrCat("Failed to parse latitude: ", lat_lng[0]));
    }
    if (!absl::SimpleAtod(lat_lng[1], &lng)) {
      return absl::InvalidArgumentError(
          absl::StrCat("Failed to parse longitude: ", lat_lng[1]));
    }
    S2LatLng latlng = S2LatLng::FromDegrees(lat, lng);
    if (!latlng.is_valid()) {
      return absl::InvalidArgumentError(
          absl::StrCat("Invalid latlng coordinates: ", lat_lng_str));
    }
    latlngs.push_back(latlng);
  }
  return latlngs;
}

absl::StatusOr<vector<S2Point>> ParsePoints(string_view str) {
  ASSIGN_OR_RETURN(vector<S2LatLng> latlngs, ParseLatLngs(str));
  vector<S2Point> points;
  points.reserve(latlngs.size());
  for (const S2LatLng& latlng : latlngs) {
    points.push_back(latlng.ToPoint());
  }
  return points;
}

absl::StatusOr<S2Point> MakePoint(string_view str) {
  ASSIGN_OR_RETURN(vector<S2Point> points, ParsePoints(str));
  if (points.size() != 1) {
    return absl::InvalidArgumentError(
        absl::StrCat("Expected 1 point, got ", points.size()));
  }
  return points[0];
}

absl::StatusOr<S2LatLng> MakeLatLng(string_view str) {
  ASSIGN_OR_RETURN(vector<S2LatLng> latlngs, ParseLatLngs(str));
  if (latlngs.size() != 1) {
    return absl::InvalidArgumentError(
        absl::StrCat("Expected 1 latlng, got ", latlngs.size()));
  }
  return latlngs[0];
}

absl::StatusOr<S2LatLngRect> MakeLatLngRect(string_view str) {
  ASSIGN_OR_RETURN(vector<S2LatLng> latlngs, ParseLatLngs(str));
  if (latlngs.empty()) {
    return absl::InvalidArgumentError(
        "No points found to construct S2LatLngRect");
  }
  S2LatLngRect rect = S2LatLngRect::FromPoint(latlngs[0]);
  for (size_t i = 1; i < latlngs.size(); ++i) {
    rect.AddPoint(latlngs[i]);
  }
  return rect;
}

absl::StatusOr<S2CellId> MakeCellId(string_view str) {
  S2CellId cell_id = S2CellId::FromDebugString(str);
  if (cell_id == S2CellId::None()) {
    return absl::InvalidArgumentError(absl::StrCat("Invalid S2CellId: ", str));
  }
  return cell_id;
}

absl::StatusOr<S2CellUnion> MakeCellUnion(string_view str) {
  vector<S2CellId> cell_ids;
  for (const string_view cell_str : SplitString(str, ',')) {
    ASSIGN_OR_RETURN(S2CellId cell_id, MakeCellId(cell_str));
    cell_ids.push_back(cell_id);
  }
  return S2CellUnion(std::move(cell_ids));
}

absl::StatusOr<unique_ptr<S2Loop>> MakeLoop(string_view str,
                                            S2Debug debug_override) {
  if (str == "empty") {
    return make_unique<S2Loop>(S2Loop::kEmpty());
  }
  if (str == "full") {
    return make_unique<S2Loop>(S2Loop::kFull());
  }
  ASSIGN_OR_RETURN(vector<S2Point> vertices, ParsePoints(str));
  return make_unique<S2Loop>(vertices, debug_override);
}

absl::StatusOr<unique_ptr<S2Polyline>> MakePolyline(string_view str,
                                                    S2Debug debug_override) {
  ASSIGN_OR_RETURN(vector<S2Point> vertices, ParsePoints(str));
  return make_unique<S2Polyline>(vertices, debug_override);
}

absl::StatusOr<unique_ptr<S2LaxPolylineShape>> MakeLaxPolyline(
    string_view str) {
  ASSIGN_OR_RETURN(vector<S2Point> vertices, ParsePoints(str));
  return make_unique<S2LaxPolylineShape>(vertices);
}

static absl::StatusOr<unique_ptr<S2Polygon>> InternalMakePolygon(
    string_view str, S2Debug debug_override, bool normalize_loops) {
  if (str == "empty") str = "";
  vector<string_view> loop_strs = SplitString(str, ';');
  vector<unique_ptr<S2Loop>> loops;
  for (const string_view loop_str : loop_strs) {
    ASSIGN_OR_RETURN(unique_ptr<S2Loop> loop,
                     MakeLoop(loop_str, debug_override));
    // Don't normalize loops that were explicitly specified as "full".
    if (normalize_loops && !loop->is_full()) loop->Normalize();
    loops.push_back(std::move(loop));
  }
  return make_unique<S2Polygon>(std::move(loops), debug_override);
}

absl::StatusOr<unique_ptr<S2Polygon>> MakePolygon(string_view str,
                                                  S2Debug debug_override) {
  return InternalMakePolygon(str, debug_override, /*normalize_loops=*/true);
}

absl::StatusOr<unique_ptr<S2Polygon>> MakeVerbatimPolygon(string_view str) {
  return InternalMakePolygon(str, S2Debug::ALLOW, /*normalize_loops=*/false);
}

absl::StatusOr<unique_ptr<S2LaxPolygonShape>> MakeLaxPolygon(string_view str) {
  vector<string_view> loop_strs = SplitString(str, ';');
  vector<vector<S2Point>> loops;
  for (const string_view loop_str : loop_strs) {
    if (loop_str == "full") {
      loops.push_back(vector<S2Point>());
    } else if (loop_str != "empty") {
      ASSIGN_OR_RETURN(vector<S2Point> points, ParsePoints(loop_str));
      loops.push_back(std::move(points));
    }
  }
  return make_unique<S2LaxPolygonShape>(std::move(loops));
}

absl::StatusOr<unique_ptr<MutableS2ShapeIndex>> MakeIndex(string_view str) {
  absl::InlinedVector<string_view, 3> strs = absl::StrSplit(str, '#');
  if (strs.size() != 3) {
    return absl::InvalidArgumentError(
        absl::StrCat("Must contain two # characters: ", str));
  }
  auto index = make_unique<MutableS2ShapeIndex>();
  vector<S2Point> points;
  for (const string_view point_str : SplitString(strs[0], '|')) {
    ASSIGN_OR_RETURN(S2Point point, MakePoint(point_str));
    points.push_back(point);
  }
  if (!points.empty()) {
    index->Add(make_unique<S2PointVectorShape>(std::move(points)));
  }
  for (const string_view line_str : SplitString(strs[1], '|')) {
    ASSIGN_OR_RETURN(unique_ptr<S2LaxPolylineShape> lax_polyline,
                     MakeLaxPolyline(line_str));
    index->Add(std::move(lax_polyline));
  }
  for (const string_view polygon_str : SplitString(strs[2], '|')) {
    ASSIGN_OR_RETURN(unique_ptr<S2LaxPolygonShape> lax_polygon,
                     MakeLaxPolygon(polygon_str));
    index->Add(std::move(lax_polygon));
  }
  return index;
}

static void AppendVertex(const S2LatLng& ll, string* out,
                         bool roundtrip_precision = false) {
  if (roundtrip_precision) {
    absl::StrAppendFormat(out, "%.17g:%.17g", ll.lat().degrees(),
                          ll.lng().degrees());
  } else {
    absl::StrAppendFormat(out, "%.15g:%.15g", ll.lat().degrees(),
                          ll.lng().degrees());
  }
}

static void AppendVertex(const S2Point& p, string* out,
                         bool roundtrip_precision = false) {
  S2LatLng ll(p);
  return AppendVertex(ll, out, roundtrip_precision);
}

static void AppendVertices(const S2Point* v, int n, string* out) {
  for (int i = 0; i < n; ++i) {
    if (i > 0) *out += ", ";
    AppendVertex(v[i], out);
  }
}

string ToString(const S2Point& point) {
  string out;
  AppendVertex(point, &out);
  return out;
}

string ToString(const S2LatLng& latlng) {
  string out;
  AppendVertex(latlng, &out);
  return out;
}

string ToString(const S2LatLngRect& rect) {
  string out;
  AppendVertex(rect.lo(), &out);
  out += ", ";
  AppendVertex(rect.hi(), &out);
  return out;
}

string ToString(const S2CellId cell_id) {
  return cell_id.ToString();
}

string ToString(const S2CellUnion& cell_union) {
  string out;
  for (S2CellId cell_id : cell_union) {
    if (!out.empty()) out += ", ";
    out += cell_id.ToString();
  }
  return out;
}

string ToString(const S2Loop& loop) {
  if (loop.is_empty()) {
    return "empty";
  } else if (loop.is_full()) {
    return "full";
  }
  string out;
  if (loop.num_vertices() > 0) {
    AppendVertices(&loop.vertex(0), loop.num_vertices(), &out);
  }
  return out;
}

string ToString(const S2Polyline& polyline) {
  string out;
  if (polyline.num_vertices() > 0) {
    AppendVertices(&polyline.vertex(0), polyline.num_vertices(), &out);
  }
  return out;
}

string ToString(const S2Polygon& polygon, string_view loop_separator) {
  if (polygon.is_empty()) {
    return "empty";
  } else if (polygon.is_full()) {
    return "full";
  }
  string out;
  for (int i = 0; i < polygon.num_loops(); ++i) {
    if (i > 0) absl::StrAppend(&out, loop_separator);
    const S2Loop& loop = *polygon.loop(i);
    AppendVertices(&loop.vertex(0), loop.num_vertices(), &out);
  }
  return out;
}

string ToString(Span<const S2Point> points) {
  string out;
  AppendVertices(points.data(), points.size(), &out);
  return out;
}

string ToString(Span<const S2LatLng> latlngs) {
  string out;
  for (size_t i = 0; i < latlngs.size(); ++i) {
    if (i > 0) out += ", ";
    AppendVertex(latlngs[i], &out);
  }
  return out;
}

string ToString(const S2Shape& shape) {
  // Polygon chains are separated by a ; instead of |.
  const char* separator = shape.dimension() == 2 ? "; " : " | ";

  string out;
  if (shape.dimension() == 1) out += "# ";
  if (shape.dimension() == 2) out += "## ";

  int nchain = 0;
  for (const auto& chain : shape.chains()) {
    if (nchain++ > 0) {
      out += separator;
    }

    int nvertex = 0;
    for (const S2Point& vertex : shape.vertices(chain)) {
      if (nvertex++ > 0) {
        out += ", ";
      }
      AppendVertex(vertex, &out);
    }
  }

  if (shape.dimension() == 1) out += " #";
  if (shape.dimension() == 0) out += " ##";
  return out;
}

string ToString(const S2LaxPolylineShape& polyline) {
  string out;
  if (polyline.num_vertices() > 0) {
    AppendVertices(&polyline.vertex(0), polyline.num_vertices(), &out);
  }
  return out;
}

string ToString(const S2LaxPolygonShape& polygon, string_view loop_separator) {
  string out;
  for (int i = 0; i < polygon.num_loops(); ++i) {
    if (i > 0) absl::StrAppend(&out, loop_separator);
    int n = polygon.num_loop_vertices(i);
    if (n == 0) {
      out += "full";
    } else {
      AppendVertices(&polygon.loop_vertex(i, 0), n, &out);
    }
  }
  return out;
}

string ToString(const S2ShapeIndex& index, bool roundtrip_precision) {
  string out;
  for (int dim = 0; dim < 3; ++dim) {
    if (dim > 0) out += "#";
    int count = 0;
    for (const S2Shape* shape : index) {
      if (shape == nullptr || shape->dimension() != dim) continue;
      out += (count > 0) ? " | " : (dim > 0) ? " " : "";
      for (int i = 0; i < shape->num_chains(); ++i, ++count) {
        if (i > 0) out += (dim == 2) ? "; " : " | ";
        S2Shape::Chain chain = shape->chain(i);
        if (chain.length == 0) {
          ABSL_DCHECK_EQ(dim, 2);
          out += "full";
        } else {
          AppendVertex(shape->edge(chain.start).v0, &out, roundtrip_precision);
        }
        int limit = chain.start + chain.length;
        if (dim != 1) --limit;
        for (int e = chain.start; e < limit; ++e) {
          out += ", ";
          AppendVertex(shape->edge(e).v1, &out, roundtrip_precision);
        }
      }
    }
    // Example output: "# #", "0:0 # #", "# # 0:0, 0:1, 1:0"
    if (dim == 1 || (dim == 0 && count > 0)) out += " ";
  }
  return out;
}

}  // namespace s2textformat
