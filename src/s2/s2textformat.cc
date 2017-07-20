// Copyright 2017 Google Inc. All Rights Reserved.
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

#include "s2/s2textformat.h"

#include <string>
#include <vector>

#include <glog/logging.h>
#include "s2/base/stringprintf.h"
#include "s2/strings/serialize.h"
#include "s2/third_party/absl/memory/memory.h"
#include "s2/third_party/absl/strings/str_split.h"
#include "s2/third_party/absl/strings/string_view.h"
#include "s2/third_party/absl/strings/strip.h"
#include "s2/s2latlng.h"
#include "s2/s2loop.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"
#include "s2/s2shapeutil.h"

using absl::MakeUnique;
using absl::string_view;
using std::pair;
using std::unique_ptr;
using std::vector;

namespace s2textformat {

static vector<string_view> SplitString(string_view str, char separator) {
  vector<string_view> result =
      strings::Split(str, separator, strings::SkipWhitespace());
  strings::StripWhitespaceInCollection(&result);
  return result;
}

static double ParseDouble(const string& str) {
  char* end_ptr = nullptr;
  double value = strtod(str.c_str(), &end_ptr);
  CHECK(end_ptr && *end_ptr == 0) << ": str == \"" << str << "\"";
  return value;
}

vector<S2LatLng> ParseLatLngs(string_view str) {
  vector<pair<string, string>> ps;
  CHECK(DictionaryParse(str, &ps)) << ": str == \"" << str << "\"";
  vector<S2LatLng> latlngs;
  for (auto const& p : ps) {
    latlngs.push_back(S2LatLng::FromDegrees(ParseDouble(p.first),
                                            ParseDouble(p.second)));
  }
  return latlngs;
}

vector<S2Point> ParsePoints(string_view str) {
  vector<S2LatLng> latlngs = ParseLatLngs(str);
  vector<S2Point> vertices;
  for (auto const& latlng : latlngs) {
    vertices.push_back(latlng.ToPoint());
  }
  return vertices;
}

S2Point MakePoint(string_view str) {
  vector<S2Point> vertices = ParsePoints(str);
  CHECK_EQ(vertices.size(), 1);
  return vertices[0];
}

S2LatLngRect MakeLatLngRect(string_view str) {
  vector<S2LatLng> latlngs = ParseLatLngs(str);
  CHECK_GT(latlngs.size(), 0);
  S2LatLngRect rect = S2LatLngRect::FromPoint(latlngs[0]);
  for (int i = 1; i < latlngs.size(); ++i) {
    rect.AddPoint(latlngs[i]);
  }
  return rect;
}

unique_ptr<S2Loop> MakeLoop(string_view str) {
  if (str == "empty") return MakeUnique<S2Loop>(S2Loop::kEmpty());
  if (str == "full") return MakeUnique<S2Loop>(S2Loop::kFull());
  vector<S2Point> vertices = ParsePoints(str);
  return MakeUnique<S2Loop>(vertices);
}

unique_ptr<S2Polyline> MakePolyline(string_view str) {
  vector<S2Point> vertices = ParsePoints(str);
  return MakeUnique<S2Polyline>(vertices);
}

unique_ptr<s2shapeutil::LaxPolyline> MakeLaxPolyline(string_view str) {
  auto vertices = ParsePoints(str);
  return MakeUnique<s2shapeutil::LaxPolyline>(vertices);
}

static unique_ptr<S2Polygon> InternalMakePolygon(string_view str,
                                                 bool normalize_loops) {
  vector<string_view> loop_strs = SplitString(str, ';');
  vector<unique_ptr<S2Loop>> loops;
  for (auto const& loop_str : loop_strs) {
    auto loop = MakeLoop(loop_str);
    if (normalize_loops) loop->Normalize();
    loops.push_back(std::move(loop));
  }
  return MakeUnique<S2Polygon>(std::move(loops));
}

unique_ptr<S2Polygon> MakePolygon(string_view str) {
  return InternalMakePolygon(str, true);
}

unique_ptr<S2Polygon> MakeVerbatimPolygon(string_view str) {
  return InternalMakePolygon(str, false);
}

unique_ptr<s2shapeutil::LaxPolygon> MakeLaxPolygon(string_view str) {
  vector<string_view> loop_strs = SplitString(str, ';');
  vector<vector<S2Point>> loops;
  for (auto const& loop_str : loop_strs) {
    if (loop_str == "full") {
      loops.push_back(vector<S2Point>());
    } else {
      loops.push_back(ParsePoints(loop_str));
    }
  }
  return MakeUnique<s2shapeutil::LaxPolygon>(loops);
}

unique_ptr<S2ShapeIndex> MakeIndex(string_view str) {
  vector<string_view> strs = strings::Split(str, '#');
  DCHECK_EQ(3, strs.size()) << "Must contain two # characters: " << str;

  auto index = MakeUnique<S2ShapeIndex>();
  vector<S2Point> points;
  for (auto const& point_str : SplitString(strs[0], '|')) {
    points.push_back(MakePoint(point_str));
  }
  if (!points.empty()) {
    index->Add(MakeUnique<s2shapeutil::PointVectorShape>(&points));
  }
  for (auto const& line_str : SplitString(strs[1], '|')) {
    index->Add(MakeLaxPolyline(line_str));
  }
  for (auto const& polygon_str : SplitString(strs[2], '|')) {
    index->Add(MakeLaxPolygon(polygon_str));
  }
  return index;
}

static void AppendVertex(S2Point const& p, string* out) {
  S2LatLng ll(p);
  StringAppendF(out, "%.15g:%.15g", ll.lat().degrees(), ll.lng().degrees());
}

static void AppendVertices(S2Point const* v, int n, string* out) {
  for (int i = 0; i < n; ++i) {
    if (i > 0) *out += ", ";
    AppendVertex(v[i], out);
  }
}

string ToString(S2Point const& point) {
  string out;
  AppendVertex(point, &out);
  return out;
}

string ToString(S2LatLngRect const& rect) {
  string out;
  AppendVertex(rect.lo().ToPoint(), &out);
  out += ", ";
  AppendVertex(rect.hi().ToPoint(), &out);
  return out;
}

string ToString(S2Loop const& loop) {
  string out;
  if (loop.num_vertices() > 0) {
    AppendVertices(&loop.vertex(0), loop.num_vertices(), &out);
  }
  return out;
}

string ToString(S2Polyline const& polyline) {
  string out;
  if (polyline.num_vertices() > 0) {
    AppendVertices(&polyline.vertex(0), polyline.num_vertices(), &out);
  }
  return out;
}

string ToString(S2Polygon const& polygon) {
  string out;
  for (int i = 0; i < polygon.num_loops(); ++i) {
    if (i > 0) out += ";\n";
    S2Loop const& loop = *polygon.loop(i);
    AppendVertices(&loop.vertex(0), loop.num_vertices(), &out);
  }
  return out;
}

string ToString(vector<S2Point> const& points) {
  string out;
  AppendVertices(points.data(), points.size(), &out);
  return out;
}

string ToString(vector<S2LatLng> const& latlngs) {
  string out;
  for (int i = 0; i < latlngs.size(); ++i) {
    if (i > 0) out += ", ";
    AppendVertex(latlngs[i].ToPoint(), &out);
  }
  return out;
}

string ToString(s2shapeutil::LaxPolyline const& polyline) {
  string out;
  if (polyline.num_vertices() > 0) {
    AppendVertices(&polyline.vertex(0), polyline.num_vertices(), &out);
  }
  return out;
}

string ToString(s2shapeutil::LaxPolygon const& polygon) {
  string out;
  for (int i = 0; i < polygon.num_loops(); ++i) {
    if (i > 0) out += ";\n";
    int n = polygon.num_loop_vertices(i);
    if (n > 0) AppendVertices(&polygon.loop_vertex(i, 0), n, &out);
  }
  return out;
}

string ToString(S2ShapeIndex const& index) {
  string out;
  for (int dim = 0; dim < 3; ++dim) {
    if (dim > 0) out += "#";
    int count = 0;
    for (int s = 0; s < index.num_shape_ids(); ++s) {
      S2Shape* shape = index.shape(s);
      if (shape == nullptr || shape->dimension() != dim) continue;
      out += (count > 0) ? " | " : (dim > 0) ? " " : "";
      for (int i = 0; i < shape->num_chains(); ++i, ++count) {
        if (i > 0) out += (dim == 2) ? "; " : " | ";
        S2Shape::Chain chain = shape->chain(i);
        AppendVertex(shape->edge(chain.start).v0, &out);
        int limit = chain.start + chain.length;
        if (dim != 1) --limit;
        for (int e = chain.start; e < limit; ++e) {
          out += ", ";
          AppendVertex(shape->edge(e).v1, &out);
        }
      }
    }
    // Example output: "# #", "0:0 # #", "# # 0:0, 0:1, 1:0"
    if (dim == 1 || (dim == 0 && count > 0)) out += " ";
  }
  return out;
}

}  // namespace s2textformat
