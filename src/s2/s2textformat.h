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

#ifndef S2_S2TEXTFORMAT_H_
#define S2_S2TEXTFORMAT_H_

// s2textformat contains a collection of functions for converting
// geometry to and from a human-readable format.  It is mainly
// intended for testing and debugging.  Be aware that the
// human-readable format is *not* designed to preserve the full
// precision of the original object, so it should not be used
// for data storage.

#include <memory>
#include <string>
#include <vector>

#include "s2/third_party/absl/strings/string_view.h"
#include "s2/s2latlngrect.h"

class S2LatLng;
class S2Loop;
class S2Polygon;
class S2Polyline;
class S2ShapeIndex;
class S2ShapeIndexBase;
namespace s2shapeutil { class LaxPolygon; }
namespace s2shapeutil { class LaxPolyline; }
namespace s2shapeutil { class LaxPolyline; }

namespace s2textformat {

// Returns an S2Point corresponding to the given a latitude-longitude
// coordinate in degrees.  Example of the input format:
//     "-20:150"
S2Point MakePoint(absl::string_view str);

// Parses a string of one or more latitude-longitude coordinates in degrees,
// and return the corresponding vector of S2LatLng points.
// Examples of the input format:
//     ""                            // no points
//     "-20:150"                     // one point
//     "-20:150, -20:151, -19:150"   // three points
std::vector<S2LatLng> ParseLatLngs(absl::string_view str);

// Parses a string in the same format as ParseLatLngs, and return the
// corresponding vector of S2Point values.
std::vector<S2Point> ParsePoints(absl::string_view str);

// Given a string in the same format as ParseLatLngs, returns the minimal
// bounding S2LatLngRect that contains the coordinates.
S2LatLngRect MakeLatLngRect(absl::string_view str);

// Given a string of latitude-longitude coordinates in degrees,
// returns a newly allocated loop.  Example of the input format:
//     "-20:150, 10:-120, 0.123:-170.652"
// The strings "empty" or "full" create an empty or full loop respectively.
std::unique_ptr<S2Loop> MakeLoop(absl::string_view str);

// Similar to MakeLoop(), but returns an S2Polyline rather than an S2Loop.
std::unique_ptr<S2Polyline> MakePolyline(absl::string_view str);

// Like MakePolyline, but returns an s2shapeutil::LaxPolyline instead.
std::unique_ptr<s2shapeutil::LaxPolyline> MakeLaxPolyline(
    absl::string_view str);

// Given a sequence of loops separated by semicolons, returns a newly
// allocated polygon.  Loops are automatically normalized by inverting them
// if necessary so that they enclose at most half of the unit sphere.
// (Historically this was once a requirement of polygon loops.  It also
// hides the problem that if the user thinks of the coordinates as X:Y
// rather than LAT:LNG, it yields a loop with the opposite orientation.)
//
// Examples of the input format:
//     "10:20, 90:0, 20:30"                                  // one loop
//     "10:20, 90:0, 20:30; 5.5:6.5, -90:-180, -15.2:20.3"   // two loops
//     ""       // the empty polygon (consisting of no loops)
//     "full"   // the full polygon (consisting of one full loop)
//     "empty"  // **INVALID** (a polygon consisting of one empty loop)
std::unique_ptr<S2Polygon> MakePolygon(absl::string_view str);

// Like MakePolygon(), except that it does not normalize loops (i.e., it
// gives you exactly what you asked for).
std::unique_ptr<S2Polygon> MakeVerbatimPolygon(absl::string_view str);

// Parses a string in the same format as MakePolygon, except that loops must
// be oriented so that the interior of the loop is always on the left, and
// polygons with degeneracies are supported.  As with MakePolygon, "full"
// denotes the full polygon and "empty" is not allowed (instead, create a
// polygon with no loops).
std::unique_ptr<s2shapeutil::LaxPolygon> MakeLaxPolygon(absl::string_view str);

// Returns an S2ShapeIndex containing the points, polylines, and loops (in the
// form of a single polygon) described by the following format:
//
//   point1|point2|... # line1|line2|... # polygon1|polygon2|...
//
// Examples:
//   1:2 | 2:3 # #                     // Two points
//   # 0:0, 1:1, 2:2 | 3:3, 4:4 #      // Two polylines
//   # # 0:0, 0:3, 3:0; 1:1, 2:1, 1:2  // Two nested loops (one polygon)
//   5:5 # 6:6, 7:7 # 0:0, 0:1, 1:0    // One of each
//
// Loops should be directed so that the region's interior is on the left.
// Loops can be degenerate (they do not need to meet S2Loop requirements).
std::unique_ptr<S2ShapeIndex> MakeIndex(absl::string_view str);

// Convert a point, lat-lng rect, loop, polyline, or polygon to the string
// format above.
string ToString(S2Point const& point);
string ToString(S2LatLngRect const& rect);
string ToString(S2Loop const& loop);
string ToString(S2Polyline const& polyline);
string ToString(S2Polygon const& polygon);
string ToString(std::vector<S2Point> const& points);
string ToString(std::vector<S2LatLng> const& points);
string ToString(s2shapeutil::LaxPolyline const& polyline);
string ToString(s2shapeutil::LaxPolygon const& polygon);

// Convert the contents of an S2ShapeIndex to the format above.  The index may
// contain S2Shapes of any type.  Shapes are reordered if necessary so that
// all point geometry (shapes of dimension 0) are first, followed by all
// polyline geometry, followed by all polygon geometry.
string ToString(S2ShapeIndexBase const& index);

}  // namespace s2textformat

#endif  // S2_S2TEXTFORMAT_H_
