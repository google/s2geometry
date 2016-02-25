// Copyright 2016 Google Inc. All Rights Reserved.
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
#include "s2/strings/split.h"
#include "s2/s2.h"
#include "s2/s2latlng.h"
#include "s2/s2loop.h"
#include "s2/s2polygon.h"
#include "s2/s2polyline.h"

using std::pair;
using std::vector;

namespace s2textformat {

static double ParseDouble(const string& str) {
  char* end_ptr = nullptr;
  double value = strtod(str.c_str(), &end_ptr);
  CHECK(end_ptr && *end_ptr == 0) << ": str == \"" << str << "\"";
  return value;
}

void ParseLatLngs(string const& str, vector<S2LatLng>* latlngs) {
  vector<pair<string, string>> ps;
  CHECK(DictionaryParse(str, &ps)) << ": str == \"" << str << "\"";
  latlngs->clear();
  for (auto const& p : ps) {
    latlngs->push_back(S2LatLng::FromDegrees(ParseDouble(p.first),
                                             ParseDouble(p.second)));
  }
}

void ParsePoints(string const& str, vector<S2Point>* vertices) {
  vector<S2LatLng> latlngs;
  ParseLatLngs(str, &latlngs);
  vertices->clear();
  for (auto const& latlng : latlngs) {
    vertices->push_back(latlng.ToPoint());
  }
}

S2Point MakePoint(string const& str) {
  vector<S2Point> vertices;
  ParsePoints(str, &vertices);
  CHECK_EQ(vertices.size(), 1);
  return vertices[0];
}

S2LatLngRect MakeLatLngRect(string const& str) {
  vector<S2LatLng> latlngs;
  ParseLatLngs(str, &latlngs);
  CHECK_GT(latlngs.size(), 0);
  S2LatLngRect rect = S2LatLngRect::FromPoint(latlngs[0]);
  for (int i = 1; i < latlngs.size(); ++i) {
    rect.AddPoint(latlngs[i]);
  }
  return rect;
}

S2Loop* MakeLoop(string const& str) {
  if (str == "empty") return new S2Loop(S2Loop::kEmpty());
  if (str == "full") return new S2Loop(S2Loop::kFull());
  vector<S2Point> vertices;
  ParsePoints(str, &vertices);
  return new S2Loop(vertices);
}

S2Polyline* MakePolyline(string const& str) {
  vector<S2Point> vertices;
  ParsePoints(str, &vertices);
  return new S2Polyline(vertices);
}

static S2Polygon* InternalMakePolygon(string const& str,
                                      bool normalize_loops) {
  vector<string> loop_strs = strings::Split(str, ';', strings::SkipEmpty());
  vector<S2Loop*> loops;
  for (auto const& loop_str : loop_strs) {
    S2Loop* loop = MakeLoop(loop_str);
    if (normalize_loops) loop->Normalize();
    loops.push_back(loop);
  }
  return new S2Polygon(&loops);  // Takes ownership.
}

S2Polygon* MakePolygon(string const& str) {
  return InternalMakePolygon(str, true);
}

S2Polygon* MakeVerbatimPolygon(string const& str) {
  return InternalMakePolygon(str, false);
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

string ToString(S2Loop const* loop) {
  string out;
  AppendVertices(&loop->vertex(0), loop->num_vertices(), &out);
  return out;
}

string ToString(S2Polyline const* polyline) {
  string out;
  AppendVertices(&polyline->vertex(0), polyline->num_vertices(), &out);
  return out;
}

string ToString(S2Polygon const* polygon) {
  string out;
  for (int i = 0; i < polygon->num_loops(); ++i) {
    if (i > 0) out += ";\n";
    S2Loop const* loop = polygon->loop(i);
    AppendVertices(&loop->vertex(0), loop->num_vertices(), &out);
  }
  return out;
}

string ToString(vector<S2Point> const& points) {
  string out;
  AppendVertices(&points[0], points.size(), &out);
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

}  // namespace s2textformat
